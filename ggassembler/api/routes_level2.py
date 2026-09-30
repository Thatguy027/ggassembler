"""``/api/level2/*`` - cassette assembly.

Touches `levels/level2_cassette` and `core` only. Changing Level 1 or Level 3
cannot reach this file, and this file cannot reach them.
"""

from __future__ import annotations

import io
from typing import Any

from fastapi import APIRouter, HTTPException, Request
from fastapi.responses import PlainTextResponse
from pydantic import BaseModel, Field

from ..core.assembly import AssemblyResult
from ..core.library import Library
from ..core import protocol
from ..core.parttypes import badge_for, color_for
from ..core.seqio import write_genbank
from ..levels import level2_cassette as level2

router = APIRouter(prefix="/api/level2", tags=["level2"])


def get_library(request: Request) -> Library:
    return request.app.state.library


class DesignRequest(BaseModel):
    """The Level 2 screen's state, as the browser holds it."""

    selections: dict[str, str] = Field(default_factory=dict)
    split_3: bool = False
    split_4: bool = False
    split_8: bool = False
    composite_left: bool = False
    composite_right: bool = False
    name: str = "cassette"

    def to_design(self) -> level2.CassetteDesign:
        return level2.CassetteDesign(
            selections=dict(self.selections),
            split_3=self.split_3,
            split_4=self.split_4,
            split_8=self.split_8,
            composite_left=self.composite_left,
            composite_right=self.composite_right,
            name=self.name,
        )


def _slot_payload(library: Library, design: level2.CassetteDesign) -> list[dict[str, Any]]:
    """Each panel: its badge, its overhangs, and the parts that fit it."""
    out = []
    for slot in level2.slots(design, library.scheme):
        options = level2.options(library, slot)
        badge_bg, badge_fg = badge_for(slot.key, library.scheme)
        out.append(
            {
                "key": slot.key,
                "label": slot.label,
                "column": slot.column,
                "description": slot.description,
                "five_prime": slot.five_prime,
                "three_prime": slot.three_prime,
                "color": color_for(slot.key, library.scheme),
                "badge_bg": badge_bg,
                "badge_fg": badge_fg,
                "selected": design.selections.get(slot.key),
                "match_count": len(options),
                "options": [
                    {
                        "name": e.name,
                        "display": e.display,
                        "component": e.component,
                        "aliases": e.aliases,
                        "path": e.path,
                        "length": e.length,
                        "part_type": e.call.part_type,
                        "connector_overhang": e.connector_overhang,
                        "internal_sites": e.sites.part_enzyme_internal,
                        # what is annotated inside the fragment this plasmid
                        # contributes, so the picker can match on it without a
                        # round trip: `tdh` has to find pYTK009
                        "labels": [
                            f.label
                            for f in library.labels_in_part(e)
                        ],
                        "component_source": e.component_source,
                    }
                    for e in options
                ],
            }
        )
    return out


def _result_payload(result: AssemblyResult, design: level2.CassetteDesign) -> dict[str, Any]:
    return {
        "ok": result.ok,
        "name": design.name,
        "length": result.length,
        "enzyme": result.enzyme,
        "parts": [
            {
                "name": p.name,
                "source_name": p.source_name,
                "part_type": p.part_type,
                "component": p.component,
                "label": p.label,
                "start": p.start,
                "end": p.end,
                "length": p.length,
                "left_overhang": p.left_overhang,
                "right_overhang": p.right_overhang,
                "color": p.color,
            }
            for p in result.parts
        ],
        "junctions": [
            {
                "overhang": j.overhang,
                "upstream": j.upstream,
                "downstream": j.downstream,
                "position": j.position,
            }
            for j in result.junctions
        ],
        "issues": level2.validation_strip(result),
        "counts": {
            "errors": len(result.errors),
            "warnings": len(result.warnings),
            "junctions": len(result.junctions),
        },
    }


@router.get("/default")
def default(request: Request) -> dict[str, Any]:
    """A worked assembly to open the screen on, built from this library."""
    library = get_library(request)
    design = level2.default_design(library)
    result = level2.build(library, design)
    payload = _result_payload(result, design)
    payload["slots"] = _slot_payload(library, design)
    payload["is_integration"] = design.is_integration
    payload["selections"] = design.selections
    payload["name"] = design.name
    return payload


def _protocol_payload(library: Library, design: level2.CassetteDesign) -> dict[str, Any]:
    """The reaction for the current design, priced in fmol rather than ng.

    The pieces are the parts themselves; the destination is whichever of them
    supplies the backbone, which at Level 2 is the type 8 position.
    """
    result = level2.build(library, design)
    if not result.ok:
        return {"ok": False, "issues": [i.message for i in result.errors]}

    selected = [
        (key, library.get(name)) for key, name in design.selections.items() if name
    ]
    pieces = [
        (entry.display, entry.length, entry.conc_ng_ul, entry.path)
        for _, entry in selected
        if entry is not None
    ]
    rx = protocol.reaction(
        name=design.name,
        enzyme=library.scheme.part_enzyme,
        parts=pieces,
        selection=next(
            (f"{e.ecoli_marker} + green/white" for _, e in selected
             if e is not None and e.ecoli_marker),
            "",
        ),
    )
    return {
        "ok": rx.ok,
        "name": rx.name,
        "enzyme": rx.enzyme,
        "total_ul": rx.total_ul,
        "fmol_each": rx.fmol_each,
        "selection": rx.selection,
        "issues": rx.issues,
        "missing": rx.missing,
        "components": [
            {
                "name": c.name, "kind": c.kind, "path": c.path, "length": c.length,
                "conc_ng_ul": c.conc_ng_ul, "fmol": c.fmol, "ng": c.ng,
                "volume_ul": c.volume_ul, "note": c.note,
            }
            for c in rx.components
        ],
        "steps": [{"label": s.label, "detail": s.detail} for s in rx.steps],
        "text": protocol.as_text(rx),
    }


@router.post("/protocol")
def protocol_route(request: Request, body: DesignRequest) -> dict[str, Any]:
    """How to set the reaction up: what to pipette, and the cycling."""
    library = get_library(request)
    return _protocol_payload(library, body.to_design())


@router.post("/protocol.txt", response_class=PlainTextResponse)
def protocol_text(request: Request, body: DesignRequest) -> PlainTextResponse:
    library = get_library(request)
    payload = _protocol_payload(library, body.to_design())
    if "text" not in payload:
        return PlainTextResponse(
            "cannot assemble: " + "; ".join(payload.get("issues", [])), status_code=422
        )
    return PlainTextResponse(
        payload["text"],
        headers={"Content-Disposition": f'attachment; filename="{body.name}-protocol.txt"'},
    )


class DecomposeRequest(BaseModel):
    cassette: str
    name: str | None = None


@router.post("/decompose")
def decompose(request: Request, body: DecomposeRequest) -> dict[str, Any]:
    """Read a finished cassette back into the parts that built it.

    The screen's other direction: instead of picking eight parts, pick a
    construct you already have and get its panels filled in, ready to swap one.
    """
    library = get_library(request)
    entry = library.get(body.cassette)
    if entry is None:
        raise HTTPException(404, f"no plasmid named {body.cassette}")

    found = level2.decompose(library, entry)
    design = found.to_design(body.name or f"{entry.name}-v2")
    payload = _result_payload(level2.build(library, design), design)
    payload["slots"] = _slot_payload(library, design)
    payload["selections"] = design.selections
    payload["name"] = design.name
    payload["split_3"] = design.split_3
    payload["split_4"] = design.split_4
    payload["split_8"] = design.split_8
    payload["composite_left"] = design.composite_left
    payload["composite_right"] = design.composite_right
    payload["source"] = {
        "name": entry.name,
        "display": entry.display,
        "length": found.length,
        "complete": found.complete,
        "covered": found.covered,
        "matches": [
            {
                "name": m.name, "display": m.display, "part_type": m.part_type,
                "start": m.start, "end": m.end, "length": m.length,
                "component": m.component,
            }
            for m in found.matches
        ],
        "gaps": [{"start": s, "length": n} for s, n in found.gaps],
    }
    return payload


@router.post("/slots")
def slots(request: Request, body: DesignRequest) -> dict[str, Any]:
    """The panels for the current view modes, with their dropdown contents."""
    library = get_library(request)
    design = body.to_design()
    return {
        "slots": _slot_payload(library, design),
        "is_integration": design.is_integration,
        "enzyme": library.scheme.part_enzyme.name,
    }


@router.post("/assemble")
def assemble(request: Request, body: DesignRequest) -> dict[str, Any]:
    """Simulate the reaction and return everything the screen draws."""
    library = get_library(request)
    design = body.to_design()
    result = level2.build(library, design)
    payload = _result_payload(result, design)
    payload["slots"] = _slot_payload(library, design)
    payload["is_integration"] = design.is_integration
    return payload


@router.post("/export", response_class=PlainTextResponse)
def export(request: Request, body: DesignRequest) -> PlainTextResponse:
    """The assembled cassette as a GenBank file."""
    library = get_library(request)
    design = body.to_design()
    result = level2.build(library, design)
    if not result.ok or result.product is None:
        messages = "; ".join(i.message for i in result.errors)
        return PlainTextResponse(f"cannot assemble: {messages}", status_code=422)

    handle = io.StringIO()
    from Bio import SeqIO

    result.product.annotations.setdefault("molecule_type", "DNA")
    SeqIO.write(result.product, handle, "genbank")
    return PlainTextResponse(
        handle.getvalue(),
        headers={"Content-Disposition": f'attachment; filename="{design.name}.gb"'},
    )


@router.post("/save")
def save(request: Request, body: DesignRequest) -> dict[str, Any]:
    """Write the product into the library folder and re-index it."""
    library = get_library(request)
    design = body.to_design()
    result = level2.build(library, design)
    if not result.ok or result.product is None:
        return {"ok": False, "issues": level2.validation_strip(result)}

    path = library.folder / f"{design.name}.gb"
    write_genbank(result.product, path)
    library.scan()
    return {"ok": True, "path": str(path), "name": design.name}
