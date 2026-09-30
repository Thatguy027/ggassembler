"""``/api/level3/*`` - multigene assembly.

Touches `levels/level3_multigene` and `core` only.
"""

from __future__ import annotations

import io
from typing import Any

from fastapi import APIRouter, Request
from fastapi.responses import PlainTextResponse
from pydantic import BaseModel, Field

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


def _payload(result: level3.MultigeneResult) -> dict[str, Any]:
    assembly = result.assembly
    return {
        "ok": result.ok,
        "name": result.design.name,
        "length": assembly.length,
        "enzyme": assembly.enzyme,
        "integration": result.integration,
        "chain": result.chain,
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
                "color": p.color,
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
    payload = _payload(level3.build(library, design))
    payload["backbone"] = design.backbone
    payload["transcription_units"] = design.transcription_units
    payload["name"] = design.name
    return payload


@router.post("/assemble")
def assemble(request: Request, body: DesignRequest) -> dict[str, Any]:
    library = get_library(request)
    return _payload(level3.build(library, body.to_design()))


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
