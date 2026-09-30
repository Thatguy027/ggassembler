"""``/api/level3/*`` - multigene assembly.

Touches `levels/level3_multigene` and `core` only.
"""

from __future__ import annotations

import io
from typing import Any

from fastapi import APIRouter, Request, Response
from fastapi.responses import PlainTextResponse
from pydantic import BaseModel, Field

from ..core import protocol
from ..core.library import Library
from ..levels import level3_multigene as level3

router = APIRouter(prefix="/api/level3", tags=["level3"])


def get_library(request: Request) -> Library:
    return request.app.state.library


class DesignRequest(BaseModel):
    backbone: str | None = None
    transcription_units: list[str] = Field(default_factory=list)
    name: str = "multigene"

    def to_design(self) -> level3.MultigeneDesign:
        return level3.MultigeneDesign(
            backbone=self.backbone,
            transcription_units=list(self.transcription_units),
            name=self.name,
        )


def _payload(result: level3.MultigeneResult, library: Library | None = None) -> dict[str, Any]:
    assembly = result.assembly
    return {
        "ok": result.ok,
        "name": result.design.name,
        "length": assembly.length,
        "enzyme": assembly.enzyme,
        "integration": result.integration,
        "chain": result.chain,
        # what is inside each transcription unit, so the screen can show the
        # promoter driving the gene rather than only the chain of cassettes
        "units": _unit_inputs(library, result.design, result) if library else [],
        "suggestions": result.suggestions,
        "scars": result.scars,
        "parts": [
            {
                "name": p.name,
                "source_name": p.source_name,
                "component": p.component,
                "label": p.label,
                "start": p.start,
                "end": p.end,
                "length": p.length,
                "left_overhang": p.left_overhang,
                "right_overhang": p.right_overhang,
                "part_type": p.part_type,
            }
            for p in assembly.parts
        ],
        "junctions": [
            {
                "overhang": j.overhang,
                "upstream": j.upstream,
                "downstream": j.downstream,
                "position": j.position,
            }
            for j in assembly.junctions
        ],
        "check_primers": [
            {"name": p.name, "sequence": p.sequence, "tm": p.tm, "position": p.position}
            for p in result.check_primers
        ],
        "issues": [
            {"level": i.level, "code": i.code, "message": i.message} for i in assembly.issues
        ],
        "counts": {
            "errors": len(assembly.errors),
            "warnings": len(assembly.warnings),
            "units": result.design.count,
        },
    }


def _contributed_length(entry: Any) -> int:
    """How much of this plasmid ends up in the construct, not how big it is.

    The two differ by the whole backbone, and using the wrong one is why the
    diagram's columns did not line up with the bar beneath them - a 5,253 bp
    plasmid contributing a 2,999 bp cassette was drawn nearly twice as wide as
    the segment it builds.
    """
    span = entry.cassette_span
    if not span or not entry.length:
        return entry.length
    return (span[1] - span[0]) % entry.length or entry.length


def _unit_inputs(
    library: Library, design: level3.MultigeneDesign, result: level3.MultigeneResult
) -> list[dict[str, Any]]:
    """One column per fragment of the construct, in the order the bar draws them.

    Mirrors `result.parts` exactly - backbone included, with no parts of its
    own - so the diagram above and the construct below share one set of
    proportions rather than two that happen to look similar.
    """
    by_name = {name: name for name in design.transcription_units if name}
    out = []
    for position, part in enumerate(result.assembly.parts):
        entry = library.get(part.source_name)
        is_unit = part.source_name in by_name
        out.append({
            "role": f"TU{position}" if is_unit else "backbone",
            "is_unit": is_unit,
            "name": part.source_name,
            "display": entry.display if entry else part.source_name,
            "length": part.length,
            "left_overhang": part.left_overhang,
            "right_overhang": part.right_overhang,
            "parts": level3.unit_contents(library, entry) if (entry and is_unit) else [],
        })
    # number the units 1..n, ignoring the backbone's place in the bar
    unit_number = 0
    for column in out:
        if column["is_unit"]:
            unit_number += 1
            column["role"] = f"TU{unit_number}"
    return out


@router.get("/options")
def options(request: Request) -> dict[str, Any]:
    """The cassettes and destination backbones in the library."""
    library = get_library(request)
    payload = level3.options(library)
    payload["enzyme"] = library.scheme.multigene_enzyme.name
    payload["linearizer"] = library.scheme.linearizer.name
    payload["connectors"] = level3.connector_names(library)
    return payload


@router.get("/default")
def default(request: Request) -> dict[str, Any]:
    """A worked multigene assembly to open the screen on, from this library.

    Level 2 opens on a finished cassette; opening Level 3 on two red errors
    instead reads as a broken page rather than an empty one.
    """
    library = get_library(request)
    design = level3.default_design(library)
    payload = _payload(level3.build(library, design), library)
    payload["backbone"] = design.backbone
    payload["transcription_units"] = design.transcription_units
    payload["name"] = design.name
    return payload


class UnitRequest(BaseModel):
    promoter: str = ""
    cds: str = ""
    terminator: str = ""
    cds_b: str = ""
    terminator_b: str = ""
    name: str = ""

    def to_spec(self) -> level3.UnitSpec:
        return level3.UnitSpec(
            promoter=self.promoter, cds=self.cds, terminator=self.terminator,
            cds_b=self.cds_b, terminator_b=self.terminator_b, name=self.name,
        )


class DesignRequestBody(BaseModel):
    """A multigene construct described by what it should express."""

    units: list[UnitRequest] = Field(default_factory=list)
    backbone: str | None = None
    shared: dict[str, str] = Field(default_factory=dict)
    name: str = "pMultigene"


def _design(library: Library, body: DesignRequestBody) -> level3.DesignReport:
    return level3.design(
        library,
        [u.to_spec() for u in body.units],
        backbone=body.backbone,
        shared=body.shared,
        name=body.name,
    )


@router.post("/design")
def design(request: Request, body: DesignRequestBody) -> dict[str, Any]:
    """The eight-part plasmids needed to build a multigene construct.

    The design direction: say what you want expressed, get the list of things
    to build first, with the connectors already worked out so they chain.
    """
    library = get_library(request)
    report = _design(library, body)
    return {
        "ok": report.ok,
        "name": report.name,
        "backbone": report.backbone,
        "multigene_length": report.multigene_length,
        "cassettes": [
            {
                "name": c.name,
                "parts": c.parts,
                "order": [s for s in level3.SLOT_ORDER if s in c.parts],
                "left_overhang": c.left,
                "right_overhang": c.right,
                "left_label": c.left_label,
                "right_label": c.right_label,
                "length": c.length,
                "ok": c.ok,
                "issues": [
                    {"level": i.level, "code": i.code, "message": i.message}
                    for i in c.issues if i.level != "info"
                ],
            }
            for c in report.cassettes
        ],
        "issues": [
            {"level": i.level, "code": i.code, "message": i.message}
            for i in report.issues if i.level != "info"
        ],
        "report": level3.design_report(report, library),
    }


@router.post("/design.zip")
def design_zip(request: Request, body: DesignRequestBody) -> Response:
    """Every plasmid in the design as a .gb, plus the build order that indexes them."""
    import io
    import zipfile

    from Bio import SeqIO

    library = get_library(request)
    report = _design(library, body)

    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", zipfile.ZIP_DEFLATED) as archive:
        archive.writestr(f"{report.name}-build-order.txt",
                         level3.design_report(report, library))
        for index, plan in enumerate(report.cassettes, start=1):
            if plan.record is None:
                continue
            handle = io.StringIO()
            plan.record.annotations.setdefault("molecule_type", "DNA")
            SeqIO.write(plan.record, handle, "genbank")
            archive.writestr(f"{index:02d}_{plan.name}.gb", handle.getvalue())
        if report.multigene is not None:
            handle = io.StringIO()
            report.multigene.annotations.setdefault("molecule_type", "DNA")
            SeqIO.write(report.multigene, handle, "genbank")
            archive.writestr(f"{len(report.cassettes) + 1:02d}_{report.name}.gb",
                             handle.getvalue())

    return Response(
        buffer.getvalue(),
        media_type="application/zip",
        headers={"Content-Disposition": f'attachment; filename="{report.name}-design.zip"'},
    )


@router.post("/assemble")
def assemble(request: Request, body: DesignRequest) -> dict[str, Any]:
    library = get_library(request)
    return _payload(level3.build(library, body.to_design()), library)


def _protocol(library: Library, design: level3.MultigeneDesign) -> dict[str, Any]:
    """The reaction that joins the cassettes into the backbone.

    Every multigene vector is a dropout whose enzyme sites leave on the
    fragment being replaced, so this is the reaction the short programme was
    written for - and the one place in the app where such a plasmid can
    actually be in the tube. Level 2 never sees one: a reversed-site plasmid is
    not a part, so no slot offers it.
    """
    backbone = library.get(design.backbone) if design.backbone else None
    units = [library.get(name) for name in design.transcription_units if name]

    pieces = [
        (e.display, e.length, e.conc_ng_ul, e.path) for e in units if e is not None
    ]
    rx = protocol.reaction(
        name=design.name,
        enzyme=library.scheme.multigene_enzyme,
        parts=pieces,
        destination=(
            (backbone.display, backbone.length, backbone.conc_ng_ul, backbone.path)
            if backbone is not None
            else None
        ),
        reversed_dropout=bool(backbone is not None and backbone.call.reversed_sites),
        selection=(backbone.ecoli_marker or "") if backbone is not None else "",
    )
    return {
        "ok": rx.ok,
        "name": rx.name,
        "enzyme": rx.enzyme,
        "total_ul": rx.total_ul,
        "fmol_each": rx.fmol_each,
        "selection": rx.selection,
        "issues": rx.issues,
        "notes": rx.notes,
        "missing": rx.missing,
        "reversed_dropout": rx.reversed_dropout,
        "components": [
            {
                "name": c.name, "kind": c.kind, "path": c.path, "length": c.length,
                "conc_ng_ul": c.conc_ng_ul, "fmol": c.fmol, "ng": c.ng,
                "volume_ul": c.volume_ul, "measured": c.measured, "note": c.note,
            }
            for c in rx.components
        ],
        "steps": [{"label": s.label, "detail": s.detail} for s in rx.steps],
        "text": protocol.as_text(rx),
    }


@router.post("/protocol")
def protocol_route(request: Request, body: DesignRequest) -> dict[str, Any]:
    """How to set the multigene reaction up: what to pipette, and the cycling."""
    return _protocol(get_library(request), body.to_design())


@router.post("/protocol.txt", response_class=PlainTextResponse)
def protocol_text(request: Request, body: DesignRequest) -> PlainTextResponse:
    payload = _protocol(get_library(request), body.to_design())
    return PlainTextResponse(
        payload["text"],
        headers={"Content-Disposition": f'attachment; filename="{body.name}-protocol.txt"'},
    )


@router.post("/export", response_class=PlainTextResponse)
def export(request: Request, body: DesignRequest) -> PlainTextResponse:
    library = get_library(request)
    result = level3.build(library, body.to_design())
    if result.assembly.product is None:
        messages = "; ".join(i.message for i in result.assembly.errors)
        return PlainTextResponse(f"cannot assemble: {messages}", status_code=422)

    from Bio import SeqIO

    handle = io.StringIO()
    result.assembly.product.annotations.setdefault("molecule_type", "DNA")
    SeqIO.write(result.assembly.product, handle, "genbank")

    design = body.to_design()
    library.log_build(
        name=design.name,
        level="level3",
        parts=[
            e for e in (
                library.get(n)
                for n in [design.backbone, *design.transcription_units] if n
            ) if e
        ],
        length=result.assembly.length,
        issues=sorted({i.code for i in result.issues if i.level != "info"}),
    )
    # kept for the Sequencing screen to check clones against later
    library.remember_expected(design.name, result.assembly.product, level="level3")
    return PlainTextResponse(
        handle.getvalue(),
        headers={"Content-Disposition": f'attachment; filename="{body.name}.gb"'},
    )


@router.post("/save")
def save(request: Request, body: DesignRequest) -> dict[str, Any]:
    library = get_library(request)
    result = level3.build(library, body.to_design())
    path = level3.save(library, result, body.name)
    if path is None:
        return {"ok": False, "issues": _payload(result)["issues"]}
    return {"ok": True, "path": path, "name": body.name}
