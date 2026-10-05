"""``/api/level1/*`` - part-plasmid construction.

Touches `levels/level1_part` and `core` only.
"""

from __future__ import annotations

import io
from typing import Any

from fastapi import APIRouter, BackgroundTasks, Request
from fastapi.responses import PlainTextResponse
from pydantic import BaseModel, Field

from . import routes_sync
from ..core.library import Library
from ..core.parttypes import DESCRIPTIONS
from ..levels import level1_part as level1

router = APIRouter(prefix="/api/level1", tags=["level1"])


def get_library(request: Request) -> Library:
    return request.app.state.library


class ConventionBody(BaseModel):
    gly_ser_linker: bool = True
    strip_stop: bool = True
    stop_and_xhoi: bool = True
    strip_start: bool = True
    infer_from_sequence: bool = True


class DesignRequest(BaseModel):
    part_type: str = "3"
    template: str | None = None
    sequence: str | None = None
    start: int = 1
    end: int = 0
    entry_vector: str | None = None
    destination: str | None = None
    mode: str = "pcr"
    domesticate: bool = False
    name: str = "new_part"
    conventions: ConventionBody = Field(default_factory=ConventionBody)

    def to_request(self) -> level1.PartRequest:
        return level1.PartRequest(
            part_type=self.part_type,
            template=self.template,
            sequence=self.sequence,
            start=self.start,
            end=self.end,
            entry_vector=self.entry_vector,
            destination=self.destination or self.entry_vector,
            mode=self.mode,
            domesticate=self.domesticate,
            name=self.name,
            conventions=level1.Convention(**self.conventions.model_dump()),
        )


def _primer(primer: level1.Primer | None) -> dict[str, Any] | None:
    if primer is None:
        return None
    return {
        "name": primer.name,
        "sequence": primer.sequence,
        "length": primer.length,
        "tm": primer.tm,
        "annealing": primer.annealing,
        "segments": [{"role": s.role, "seq": s.seq} for s in primer.segments],
    }


def _payload(result: level1.Level1Design) -> dict[str, Any]:
    return {
        "ok": result.ok,
        "part_type": result.part_type,
        "five_prime": result.five_prime,
        "three_prime": result.three_prime,
        "body_length": len(result.body),
        "entry_vector": result.entry_vector,
        "destination": result.destination,
        "cloning_enzyme": result.cloning_enzyme,
        "destination_marker": result.destination_marker,
        "destination_universal": result.destination_universal,
        "accepts": list(result.accepts) if result.accepts else None,
        "validated_as": result.validated_as,
        "product_length": len(result.product.seq) if result.product else 0,
        "mode": result.mode,
        "fragments": [
            {
                "index": f.index,
                "start": f.start,
                "end": f.end,
                "length": f.length,
                "left_overhang": f.left_overhang,
                "right_overhang": f.right_overhang,
                "forward": _primer(f.forward),
                "reverse": _primer(f.reverse),
                "mutations": [
                    {
                        "position": m.position,
                        "from_base": m.from_base,
                        "to_base": m.to_base,
                        "enzyme": m.enzyme,
                        "silent": m.silent,
                        "gene": m.gene,
                        "codon_from": m.codon_from,
                        "codon_to": m.codon_to,
                        "description": m.description,
                    }
                    for m in f.mutations
                ],
            }
            for f in result.fragments
        ],
        "oligos": [_primer(o) for o in result.oligos],
        "gblock": result.gblock,
        "ligation": result.ligation,
        "binding": result.binding,
        "issues": [
            {"level": i.level, "code": i.code, "message": i.message} for i in result.issues
        ],
        "counts": {
            "errors": len(result.errors),
            "warnings": len(result.warnings),
            "fragments": len(result.fragments),
        },
    }


@router.get("/options")
def options(request: Request) -> dict[str, Any]:
    """Part types, entry vectors and templates the screen can offer."""
    library = get_library(request)
    scheme = library.scheme
    types = []
    for name in (*scheme.atoms, *(w for _, _, w in scheme.merges), "234", "678"):
        if name in {t["name"] for t in types}:
            continue
        types.append(
            {
                "name": name,
                "description": DESCRIPTIONS.get(name, "composite part"),
                "five_prime": None,
            }
        )
    from ..core.parttypes import type_overhangs

    for entry in types:
        pair = type_overhangs(entry["name"], scheme)
        entry["five_prime"], entry["three_prime"] = pair if pair else (None, None)
    types.sort(key=lambda t: (len(t["name"]), t["name"]))

    return {
        "types": types,
        "destinations": [
            {
                "name": d.name,
                "display": d.display,
                "component": d.component,
                "cloning_enzyme": d.cloning_enzyme,
                "accepts": list(d.accepts),
                "part_type": d.part_type,
                "universal": not d.supplies_part_sites,
                "kind": d.kind,
                "marker": d.marker,
                "length": d.backbone_length,
            }
            for d in level1.part_destinations(library, library.scheme)
        ],
        "entry_vectors": [
            {"name": e.name, "display": e.display, "component": e.component,
             "length": e.length, "marker": e.ecoli_marker}
            for e in library.entry_vectors()
        ],
        "templates": [
            {
                "name": e.name,
                "display": e.display,
                "component": e.component,
                "length": e.length,
                "part_type": e.call.part_type,
            }
            for e in library.unique_entries()
        ],
        "enzymes": {
            "part": scheme.part_enzyme.name,
            "multigene": scheme.multigene_enzyme.name,
            "linearizer": scheme.linearizer.name,
        },
        "entry_overhangs": list(scheme.entry_vector_overhangs),
    }


@router.post("/design")
def design(request: Request, body: DesignRequest) -> dict[str, Any]:
    library = get_library(request)
    return _payload(level1.design(library, body.to_request()))


@router.post("/export", response_class=PlainTextResponse)
def export(request: Request, body: DesignRequest) -> PlainTextResponse:
    """The predicted part plasmid as GenBank."""
    library = get_library(request)
    result = level1.design(library, body.to_request())
    if result.product is None:
        messages = "; ".join(i.message for i in result.errors) or "design incomplete"
        return PlainTextResponse(f"cannot build a part plasmid: {messages}", status_code=422)

    from Bio import SeqIO

    handle = io.StringIO()
    result.product.annotations.setdefault("molecule_type", "DNA")
    SeqIO.write(result.product, handle, "genbank")
    return PlainTextResponse(
        handle.getvalue(),
        headers={"Content-Disposition": f'attachment; filename="{body.name}.gb"'},
    )


@router.post("/save")
def save(request: Request, body: DesignRequest, tasks: BackgroundTasks) -> dict[str, Any]:
    """Write the predicted part plasmid into the library folder and re-index it."""
    library = get_library(request)
    result = level1.design(library, body.to_request())
    path = level1.save(library, result, body.name)
    if path is None:
        return {"ok": False, "issues": _payload(result)["issues"]}
    routes_sync.after_save(request, tasks, f"add {body.name}")
    return {"ok": True, "path": path, "name": body.name, "validated_as": result.validated_as}
