"""``/api/plate/*`` - the standing source plate, and the sweeps built from it.

Follows the same rule as every other screen's router: it touches `core` and its
own `levels/plate_batch`, and no other level's module.
"""

from __future__ import annotations

from typing import Any

import io
import zipfile

from Bio import SeqIO
from fastapi import APIRouter, HTTPException, Request, Response
from fastapi.responses import PlainTextResponse
from pydantic import BaseModel

from .. import __version__
from ..core.library import Library
from ..levels import plate_batch

router = APIRouter(prefix="/api/plate", tags=["plate"])


def get_library(request: Request) -> Library:
    return request.app.state.library


@router.get("/source")
def source(request: Request) -> dict[str, Any]:
    """The standing plate as it is, with everything wrong with it."""
    library = get_library(request)
    return plate_batch.summary(plate_batch.load(library), library)


class WellRequest(BaseModel):
    """One well of the standing plate. An empty plasmid clears it."""

    well: str
    plasmid: str = ""
    conc_ng_ul: float | None = None
    volume_ul: float | None = None
    note: str = ""


@router.post("/source/well")
def set_well(request: Request, body: WellRequest) -> dict[str, Any]:
    library = get_library(request)
    try:
        plate_batch.set_well(
            library,
            well=body.well,
            plasmid=body.plasmid,
            conc_ng_ul=body.conc_ng_ul,
            volume_ul=body.volume_ul,
            note=body.note,
        )
    except ValueError as error:
        raise HTTPException(422, str(error)) from None
    return plate_batch.summary(plate_batch.load(library), library)


class LayoutRequest(BaseModel):
    """The plate itself, rather than what is in it."""

    id: str = ""
    labware: str = ""
    dead_volume_ul: float | None = None


@router.post("/source/layout")
def set_layout(request: Request, body: LayoutRequest) -> dict[str, Any]:
    library = get_library(request)
    try:
        plate_batch.set_layout(
            library,
            plate_id=body.id,
            labware=body.labware,
            dead_volume_ul=body.dead_volume_ul,
        )
    except ValueError as error:
        raise HTTPException(422, str(error)) from None
    return plate_batch.summary(plate_batch.load(library), library)


@router.get("/candidates")
def candidates(request: Request, part_type: str | None = None) -> list[dict[str, Any]]:
    """Plasmids that could go in a well, for the assignment dropdown.

    Narrowed by part type when one is given, because the question being asked
    is almost always "which promoters do I have", not "show me 622 names".
    """
    library = get_library(request)
    entries = library.parts_of_type(part_type) if part_type else library.unique_entries()
    return [
        {
            "name": e.name,
            "display": e.display,
            "part_type": e.call.part_type,
            "component": e.component,
            "length": e.length,
            "conc_ng_ul": e.conc_ng_ul,
        }
        for e in sorted(entries, key=lambda e: e.name.lower())
    ]


@router.get("/wells")
def well_names() -> dict[str, Any]:
    """The geometry every grid on this screen is drawn against."""
    return {
        "rows": list(plate_batch.ROWS),
        "columns": list(plate_batch.COLUMNS),
        "row_major": plate_batch.wells(),
        "column_major": plate_batch.wells_column_major(),
        "default_labware": plate_batch.DEFAULT_LABWARE,
        "default_dead_volume_ul": plate_batch.DEFAULT_DEAD_VOLUME,
    }


# --------------------------------------------------------------------------- #
# the sweep
# --------------------------------------------------------------------------- #


class FactorRequest(BaseModel):
    position: str = ""
    candidates: list[str] = []


class SweepRequest(BaseModel):
    """A base cassette, and up to two positions varying across the plate."""

    base: dict[str, str] = {}
    row: FactorRequest = FactorRequest()
    column: FactorRequest = FactorRequest()
    split_3: bool = False
    split_4: bool = False
    split_8: bool = False
    composite_left: bool = False
    composite_right: bool = False
    name: str = "sweep"
    pattern: str = "{base}_{row}_{col}"

    def to_design(self) -> plate_batch.SweepDesign:
        return plate_batch.SweepDesign(
            base=dict(self.base),
            row=plate_batch.Factor(self.row.position, list(self.row.candidates)),
            column=plate_batch.Factor(self.column.position, list(self.column.candidates)),
            split_3=self.split_3, split_4=self.split_4, split_8=self.split_8,
            composite_left=self.composite_left, composite_right=self.composite_right,
            name=self.name, pattern=self.pattern,
        )


def _issues(issues) -> list[dict[str, str]]:
    return [{"level": i.level, "code": i.code, "message": i.message} for i in issues]


@router.get("/slots")
def slot_options(
    request: Request,
    split_3: bool = False,
    split_4: bool = False,
    split_8: bool = False,
    composite_left: bool = False,
    composite_right: bool = False,
) -> list[dict[str, Any]]:
    """The positions a sweep can vary, and what each one could hold."""
    library = get_library(request)
    keys = plate_batch.slot_keys(
        split_3=split_3, split_4=split_4, split_8=split_8,
        composite_left=composite_left, composite_right=composite_right,
    )
    out = []
    for slot in plate_batch.slots_for(keys, library.scheme):
        options = library.parts_with_overhangs(*slot.overhangs)
        out.append({
            "key": slot.key,
            "label": slot.label,
            "description": slot.description,
            "column": slot.column,
            "five_prime": slot.five_prime,
            "three_prime": slot.three_prime,
            "options": [
                {"name": e.name, "display": e.display, "component": e.component,
                 "length": e.length, "part_type": e.call.part_type}
                for e in options
            ],
        })
    return out


@router.post("/sweep")
def sweep(request: Request, body: SweepRequest) -> dict[str, Any]:
    """Assemble every well in silico, before anything is pipetted.

    Candidates that cannot sit in the position they were chosen for are
    reported on their own rather than as ninety-six identical assembly
    failures: that is a fact about the candidate, not about the plate.
    """
    library = get_library(request)
    design = body.to_design()

    bad = plate_batch.check_candidates(library, design)
    if any(i.level == "error" for i in bad):
        return {
            "ok": False,
            "blocked": True,
            "issues": _issues(bad),
            "wells": [], "axes": {"rows": [], "columns": []},
            "counts": {"ok": 0, "warn": 0, "err": 0, "total": 0},
        }

    built = plate_batch.build_plate(library, design)
    counts = {
        "ok": sum(1 for w in built if w.status == "ok"),
        "warn": sum(1 for w in built if w.status == "warn"),
        "err": sum(1 for w in built if w.status == "err"),
        "total": len(built),
    }
    usage = plate_batch.part_usage(built)

    return {
        "ok": counts["err"] == 0,
        "blocked": False,
        "issues": _issues(bad),
        "wells": [
            {
                "well": w.well,
                "row": w.row_index,
                "column": w.column_index,
                "name": w.name,
                "row_variant": w.row_variant,
                "column_variant": w.column_variant,
                "length": w.length,
                "status": w.status,
                "selections": w.selections,
                "issues": _issues(w.issues),
            }
            for w in built
        ],
        "axes": plate_batch.axis_summary(built, design),
        "counts": counts,
        "shared_positions": plate_batch.shared_positions(design),
        "usage": usage,
    }


# --------------------------------------------------------------------------- #
# the reaction, and what comes out of it
# --------------------------------------------------------------------------- #


class SetupRequest(BaseModel):
    """The bench's numbers. Defaults, never constants."""

    total_ul: float = 10.0
    part_ul: float = 1.0
    target_fmol: float = 20.0
    overage: float = 1.15
    buffer_ul: float = 1.0
    ligase_ul: float = 0.5
    enzyme_ul: float = 0.5

    def to_setup(self) -> plate_batch.ReactionSetup:
        return plate_batch.ReactionSetup(
            total_ul=self.total_ul, part_ul=self.part_ul,
            target_fmol=self.target_fmol, overage=self.overage,
            buffer_ul=self.buffer_ul, ligase_ul=self.ligase_ul,
            enzyme_ul=self.enzyme_ul,
        )


class RunRequest(SweepRequest):
    """A sweep, plus how to cost it."""

    setup: SetupRequest = SetupRequest()


def _planned(library, body: RunRequest):
    """Everything a costed run needs, or an HTTP error naming what stopped it."""
    design = body.to_design()
    bad = plate_batch.check_candidates(library, design)
    if any(i.level == "error" for i in bad):
        raise HTTPException(422, "; ".join(i.message for i in bad))

    built = plate_batch.build_plate(library, design)
    plate = plate_batch.load(library)
    plan = plate_batch.plan_reaction(library, design, built, plate, body.setup.to_setup())
    return design, built, plate, plan


@router.post("/reaction")
def reaction(request: Request, body: RunRequest) -> dict[str, Any]:
    """The master-mix recipe for the whole plate, and the per-well transfers."""
    library = get_library(request)
    design, built, plate, plan = _planned(library, body)

    errors = sum(1 for w in built if w.status == "err")
    return {
        "ok": plan.ok and errors == 0,
        "wells": plan.wells,
        "failed_wells": errors,
        "mix_per_well_ul": plan.mix_per_well_ul,
        "mix_total_ul": plan.mix_total_ul,
        "overage": plan.setup.overage,
        "mix": [
            {
                "name": c.name, "kind": c.kind,
                "per_well_ul": c.per_well_ul, "total_ul": c.total_ul,
                "conc_ng_ul": c.conc_ng_ul, "length": c.length,
                "fmol": c.fmol, "measured": c.measured, "note": c.note,
            }
            for c in plan.mix
        ],
        "transfers": [
            {"source": t.source, "destination": t.destination, "plasmid": t.plasmid,
             "volume_ul": t.volume_ul, "axis": t.axis}
            for t in plan.transfers
        ],
        # the count that says whether this is a protocol or a fantasy: 96 x 8
        # is 768 moves, and a sweep is meant to be a few dozen
        "aspirations": len({(t.source, t.axis) for t in plan.transfers}) + plan.wells // 8 + 1,
        "draw": [
            {"well": w, "volume_ul": round(v, 2)} for w, v in sorted(plan.draw.per_well.items())
        ],
        "issues": _issues(plan.issues),
    }


@router.post("/map.csv", response_class=PlainTextResponse)
def plate_map(request: Request, body: RunRequest) -> PlainTextResponse:
    """The file that gets printed and taped to the bench."""
    library = get_library(request)
    design, built, plate, _plan = _planned(library, body)
    csv_text = plate_batch.plate_map_csv(built, design, plate.id)
    return PlainTextResponse(
        csv_text,
        headers={"Content-Disposition": f'attachment; filename="{design.name}-plate.csv"'},
    )


@router.post("/products.gb", response_class=PlainTextResponse)
def products(request: Request, body: RunRequest) -> PlainTextResponse:
    """Every predicted plasmid, as one multi-record GenBank file."""
    library = get_library(request)
    design, built, _plate, _plan = _planned(library, body)

    made = plate_batch.products(library, design, built)
    if not made:
        raise HTTPException(422, "no well assembles, so there is nothing to write")

    handle = io.StringIO()
    for _well, record in made:
        record.annotations.setdefault("molecule_type", "DNA")
    SeqIO.write([record for _well, record in made], handle, "genbank")
    return PlainTextResponse(
        handle.getvalue(),
        headers={"Content-Disposition": f'attachment; filename="{design.name}-products.gb"'},
    )


@router.post("/products.zip")
def products_zip(request: Request, body: RunRequest) -> Response:
    """One GenBank per well, for anything that wants files rather than records."""
    library = get_library(request)
    design, built, plate, _plan = _planned(library, body)

    made = plate_batch.products(library, design, built)
    if not made:
        raise HTTPException(422, "no well assembles, so there is nothing to write")

    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", zipfile.ZIP_DEFLATED) as archive:
        for well, record in made:
            handle = io.StringIO()
            record.annotations.setdefault("molecule_type", "DNA")
            SeqIO.write(record, handle, "genbank")
            archive.writestr(f"{well.well}_{well.name}.gb", handle.getvalue())
        archive.writestr(f"{design.name}-plate.csv",
                         plate_batch.plate_map_csv(built, design, plate.id))
    buffer.seek(0)
    return Response(
        buffer.getvalue(),
        media_type="application/zip",
        headers={"Content-Disposition": f'attachment; filename="{design.name}-products.zip"'},
    )


@router.post("/commit")
def commit(request: Request, body: RunRequest) -> dict[str, Any]:
    """Log every well and take the run's volumes out of the source plate.

    Separate from generating a protocol on purpose. The source plate records a
    physical object, so decrementing it is a claim that the liquid actually
    moved - and that claim belongs to one explicit action rather than to
    whichever endpoint happened to be called last.
    """
    library = get_library(request)
    design, built, plate, plan = _planned(library, body)
    if not plan.ok:
        raise HTTPException(422, "; ".join(i.message for i in plan.issues) or "the run cannot be drawn")

    logged = plate_batch.log_builds(library, design, built, plate.id)
    plate_batch.apply_draw(library, plan.draw)
    return {
        "ok": True,
        "logged": logged,
        "plate": plate_batch.summary(plate_batch.load(library), library),
    }


@router.post("/protocol")
def protocol_preview(request: Request, body: RunRequest) -> dict[str, Any]:
    """The generated protocol, and what the simulator made of it.

    Simulation runs here rather than on download, because the point of it is to
    decide whether a download should be offered at all.
    """
    library = get_library(request)
    design, built, plate, plan = _planned(library, body)
    if not plan.ok:
        raise HTTPException(422, "; ".join(i.message for i in plan.issues))

    source = plate_batch.generate_protocol(
        library, design, built, plate, plan, version=__version__
    )
    checked = plate_batch.simulate(source)
    return {
        "source": source,
        "lines": source.count("\n") + 1,
        "transfers": len(plan.transfers),
        "simulation": {
            "ran": checked.ran,
            "ok": checked.ok,
            "verdict": checked.verdict,
            "detail": checked.detail,
            "output": checked.output,
        },
        # never offered on an unchecked protocol: "could not check" is not
        # "checked", and a guard that treats them alike is worse than none
        "downloadable": checked.ran and checked.ok,
    }


@router.post("/protocol.py", response_class=PlainTextResponse)
def protocol_file(request: Request, body: RunRequest) -> PlainTextResponse:
    library = get_library(request)
    design, built, plate, plan = _planned(library, body)
    if not plan.ok:
        raise HTTPException(422, "; ".join(i.message for i in plan.issues))

    source = plate_batch.generate_protocol(
        library, design, built, plate, plan, version=__version__
    )
    checked = plate_batch.simulate(source)
    if not (checked.ran and checked.ok):
        raise HTTPException(
            422,
            f"this protocol was {checked.verdict}"
            + (f": {checked.detail}" if checked.detail else ""),
        )
    return PlainTextResponse(
        source,
        headers={"Content-Disposition": f'attachment; filename="{design.name}-protocol.py"'},
    )
