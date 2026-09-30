"""Shared, read-only library endpoints (``/api/library/*``).

Every screen may read these; only the Library screen writes to them. Nothing
here knows what a level is, which is why all three levels can depend on it
without depending on each other.
"""

from __future__ import annotations

from typing import Any

import csv
import io

from fastapi import APIRouter, HTTPException, Request
from fastapi.responses import PlainTextResponse
from pydantic import BaseModel

from ..core.library import Library, PlasmidEntry, triage
from ..core.parttypes import badge_for, color_for, is_composite
from ..core.search import search as run_search

router = APIRouter(prefix="/api/library", tags=["library"])


def get_library(request: Request) -> Library:
    """The one library this server was started against."""
    return request.app.state.library


class OverrideRequest(BaseModel):
    path: str
    part_type: str | None = None
    reason: str = "set by hand"


class ConcentrationRequest(BaseModel):
    path: str
    conc_ng_ul: float | None = None


def compact(entry: PlasmidEntry, scheme_name: str = "YTK") -> dict[str, Any]:
    """The shape every dropdown and table row needs, and nothing more."""
    from ..core.parttypes import get_scheme

    scheme = get_scheme(entry.scheme or scheme_name)
    badge_bg, badge_fg = badge_for(entry.call.part_type, scheme)
    return {
        "path": entry.path,
        "name": entry.name,
        "display": entry.display,
        "component": entry.component,
        "component_source": entry.component_source,
        "conc_ng_ul": entry.conc_ng_ul,
        "aliases": entry.aliases,
        "length": entry.length,
        "part_type": entry.call.part_type,
        "confidence": entry.call.confidence,
        "reason": entry.call.reason,
        "source": entry.call.source,
        "conflict": entry.call.conflict,
        "five_prime": entry.call.five_prime,
        "three_prime": entry.call.three_prime,
        "reversed_sites": entry.call.reversed_sites,
        "roles": entry.roles,
        "connector_overhang": entry.connector_overhang,
        "level3_ready": entry.level3_ready,
        "cassette_overhangs": list(entry.cassette_overhangs) if entry.cassette_overhangs else None,
        "multigene_overhangs": entry.multigene_overhangs,
        "internal_multigene_positions": entry.internal_multigene_positions,
        "ecoli_marker": entry.ecoli_marker,
        "internal_sites": {
            "part_enzyme": entry.sites.part_enzyme_internal,
            "multigene_enzyme": entry.sites.multigene_enzyme_internal,
            "linearizer": entry.sites.linearizer_total,
        },
        "site_totals": {
            "part_enzyme": entry.sites.part_enzyme_total,
            "multigene_enzyme": entry.sites.multigene_enzyme_total,
            "linearizer": entry.sites.linearizer_total,
        },
        "enzymes": {
            "part": scheme.part_enzyme.name,
            "multigene": scheme.multigene_enzyme.name,
            "linearizer": scheme.linearizer.name,
        },
        "composite": is_composite(entry.call.part_type, scheme),
        "unrecognised": entry.call.part_type is None and not entry.is_assembled,
        "triage": triage(entry),
        "assumed_circular": entry.call.assumed_circular,
        "needs_attention": bool(
            entry.call.conflict
            or entry.sites.part_enzyme_internal
            or (entry.call.part_type is None and not entry.is_assembled)
            or (entry.call.part_type and entry.call.confidence != "digest")
            or entry.internal_multigene_positions
        ),
        "color": color_for(entry.call.part_type, scheme),
        "badge_bg": badge_bg,
        "badge_fg": badge_fg,
    }


@router.get("/summary")
def summary(request: Request) -> dict[str, Any]:
    """Counts for the header strip and the Library screen."""
    library = get_library(request)
    entries = library.unique_entries()
    return {
        "roots": [str(r) for r in library.roots],
        "base": str(library.base),
        "scheme": library.scheme.name,
        "enzymes": {
            "part": library.scheme.part_enzyme.name,
            "multigene": library.scheme.multigene_enzyme.name,
            "linearizer": library.scheme.linearizer.name,
        },
        "total": len(library.entries),
        "unique": len(entries),
        "merged": len(library.entries) - len(entries),
        "typed": sum(1 for e in entries if e.call.part_type),
        "by_digest": sum(1 for e in entries if e.call.confidence == "digest"),
        "parts": sum(1 for e in entries if e.is_part),
        "cassettes": len(library.cassettes()),
        "multigene_vectors": len(library.multigene_vectors()),
        "entry_vectors": len(library.entry_vectors()),
        "conflicts": sum(1 for e in entries if e.call.conflict),
        # composite *parts* you could build with; a dropout spans 2-3-4 too but
        # is consumed by the reaction, so counting it here would mislead
        "composite": sum(
            1
            for e in entries
            if is_composite(e.call.part_type, library.scheme) and not e.call.reversed_sites
        ),
        "unrecognised": sum(1 for e in entries if e.call.part_type is None and not e.is_assembled),
        "internal_sites": sum(1 for e in entries if e.sites.part_enzyme_internal),
        "internal_multigene": sum(1 for e in entries if e.internal_multigene_positions),
        "errors": library.errors,
        "duplicates": library.duplicates(),
    }


@router.get("/plasmids")
def plasmids(
    request: Request,
    five: str | None = None,
    three: str | None = None,
    part_type: str | None = None,
    role: str | None = None,
) -> list[dict[str, Any]]:
    """Every plasmid, optionally narrowed to the ones a slot can use."""
    library = get_library(request)
    if five and three:
        entries = library.parts_with_overhangs(five, three)
    elif part_type:
        entries = library.parts_of_type(part_type)
    else:
        entries = library.unique_entries()
    if role:
        entries = [e for e in entries if role in e.roles]
    return [compact(e) for e in entries]


@router.get("/search")
def search(
    request: Request,
    q: str = "",
    limit: int = 25,
    part_type: str | None = None,
    role: str | None = None,
    usable: str | None = None,
) -> dict[str, Any]:
    """Plasmids whose annotations match `q`, best first.

    Shared, so every screen can offer the same box; each screen decides for
    itself what to do with a hit, which is the only level-specific part.
    `usable` is how a screen says what it can work with: ``part`` for the
    cassette screen, ``multigene`` for the multigene one.
    """
    library = get_library(request)
    scheme = library.scheme
    hits, suppressed = run_search(
        library, q, limit=limit, part_type=part_type, role=role, usable=usable
    )
    out = []
    for hit in hits:
        badge_bg, badge_fg = badge_for(hit.part_type, scheme)
        out.append({
            "name": hit.name,
            "display": hit.display,
            "path": hit.path,
            "component": hit.component,
            "component_source": hit.component_source,
            "part_type": hit.part_type,
            "roles": hit.roles,
            "length": hit.length,
            "score": hit.score,
            "where": hit.where,
            "matched": hit.matched,
            "aliases": hit.aliases,
            "five_prime": hit.five_prime,
            "three_prime": hit.three_prime,
            "connector_overhang": hit.connector_overhang,
            "cassette_overhangs": hit.cassette_overhangs,
            "internal_sites": hit.internal_sites,
            "in_part": hit.in_part,
            "usable_as_part": hit.usable_as_part,
            "usable_as_unit": hit.usable_as_unit,
            "color": color_for(hit.part_type, scheme),
            "badge_bg": badge_bg,
            "badge_fg": badge_fg,
        })
    return {"query": q, "count": len(out), "hits": out, "suppressed": suppressed}


@router.get("/plasmids/{name}")
def plasmid(request: Request, name: str) -> dict[str, Any]:
    """One plasmid in full, features included."""
    library = get_library(request)
    entry = library.get(name)
    if entry is None:
        raise HTTPException(404, f"no plasmid named {name}")
    data = compact(entry)
    data["record_id"] = entry.record_id
    data["topology"] = entry.topology
    data["checksum"] = entry.checksum
    data["sites"] = {
        "part_enzyme_total": entry.sites.part_enzyme_total,
        "multigene_enzyme_total": entry.sites.multigene_enzyme_total,
        "linearizer_total": entry.sites.linearizer_total,
    }
    data["features"] = [
        {"type": f.type, "label": f.label, "start": f.start, "end": f.end, "strand": f.strand}
        for f in entry.features
    ]
    return data


@router.get("/connectors")
def connectors(request: Request) -> dict[str, list[str]]:
    """Connector overhang -> the parts carrying it, learned from the files."""
    return get_library(request).connector_overhangs()


class CanonicalRequest(BaseModel):
    name: str
    """The copy to promote, by file name, path or record id."""
    canonical: bool = True


@router.get("/duplicates")
def duplicate_groups(request: Request) -> list[dict[str, Any]]:
    """Each group of identical files and which copy currently names the group."""
    return get_library(request).duplicate_groups()


@router.post("/canonical")
def canonical(request: Request, body: CanonicalRequest) -> dict[str, Any]:
    """Choose which copy of an identical sequence speaks for the group.

    Takes the name the user actually sees, since a merged row shows names and
    aliases rather than paths.
    """
    library = get_library(request)
    target = next(
        (e for e in library.sorted_entries()
         if body.name in (e.name, e.path, e.record_id)),
        None,
    )
    if target is None:
        raise HTTPException(404, f"no plasmid named {body.name}")
    library.set_canonical(target.path, body.canonical)
    entry = next(
        (e for e in library.unique_entries() if e.checksum == target.checksum), None
    )
    return compact(entry) if entry else {"path": target.path}


class DesignRecord(BaseModel):
    level: str
    name: str
    design: dict[str, Any] = {}


@router.get("/designs")
def designs(request: Request, level: str | None = None) -> list[dict[str, Any]]:
    """Saved designs, newest first. Shared, because every level saves the same way."""
    return get_library(request).designs(level)


@router.post("/designs")
def save_design(request: Request, body: DesignRecord) -> dict[str, Any]:
    try:
        return get_library(request).save_design(body.level, body.name, body.design)
    except ValueError as error:
        raise HTTPException(422, str(error)) from None


@router.delete("/designs")
def delete_design(request: Request, level: str, name: str) -> dict[str, Any]:
    return {"deleted": get_library(request).delete_design(level, name)}


@router.get("/builds")
def builds(request: Request, limit: int = 100) -> list[dict[str, Any]]:
    """What has been exported, newest first, with each part's checksum."""
    return get_library(request).builds(limit)


@router.post("/rescan")
def rescan(request: Request, force: bool = False) -> dict[str, Any]:
    library = get_library(request)
    library.scan(force=force)
    return summary(request)


@router.post("/override")
def override(request: Request, body: OverrideRequest) -> dict[str, Any]:
    """Assign a type by hand. Manual assignments always beat detection."""
    library = get_library(request)
    if library.get(body.path) is None and body.path not in library.entries:
        raise HTTPException(404, f"no plasmid at {body.path}")
    library.set_override(body.path, body.part_type, body.reason)
    entry = library.entries.get(body.path)
    return compact(entry) if entry else {"path": body.path, "part_type": body.part_type}


@router.post("/concentration")
def concentration(request: Request, body: ConcentrationRequest) -> dict[str, Any]:
    """Record what a prep measured at, in ng/µL.

    This is the one number a reaction setup needs that no file can supply, so
    it is kept beside the index rather than in it.
    """
    library = get_library(request)
    if library.get(body.path) is None and body.path not in library.entries:
        raise HTTPException(404, f"no plasmid at {body.path}")
    try:
        library.set_concentration(body.path, body.conc_ng_ul)
    except ValueError as error:
        raise HTTPException(422, str(error)) from None
    entry = library.entries.get(body.path)
    return compact(entry) if entry else {"path": body.path, "conc_ng_ul": body.conc_ng_ul}


@router.get("/index.csv", response_class=PlainTextResponse)
def index_csv(request: Request) -> PlainTextResponse:
    """The whole index as a spreadsheet, for a lab notebook or a shared drive."""
    library = get_library(request)
    scheme = library.scheme
    buffer = io.StringIO()
    writer = csv.writer(buffer)
    writer.writerow([
        "path", "name", "length", "topology", "part_type", "confidence", "source",
        "five_prime", "three_prime", "reversed_sites", "roles", "connector_overhang",
        f"{scheme.part_enzyme.name}_total", f"{scheme.part_enzyme.name}_internal",
        f"{scheme.multigene_enzyme.name}_total", f"{scheme.linearizer.name}_total",
        "ecoli_marker", "evidence", "conflict",
    ])
    for entry in library.sorted_entries():
        writer.writerow([
            entry.path, entry.name, entry.length, entry.topology,
            entry.call.part_type or "", entry.call.confidence, entry.call.source,
            entry.call.five_prime or "", entry.call.three_prime or "",
            "yes" if entry.call.reversed_sites else "",
            " ".join(entry.roles), entry.connector_overhang or "",
            entry.sites.part_enzyme_total, entry.sites.part_enzyme_internal,
            entry.sites.multigene_enzyme_total, entry.sites.linearizer_total,
            entry.ecoli_marker or "", entry.call.reason, entry.call.conflict or "",
        ])
    return PlainTextResponse(
        buffer.getvalue(),
        media_type="text/csv",
        headers={"Content-Disposition": 'attachment; filename="ggasm-index.csv"'},
    )


@router.get("/types")
def types(request: Request) -> list[dict[str, object]]:
    """Every part type the scheme knows, for the filter and the assign dialog."""
    from ..core.parttypes import DESCRIPTIONS

    library = get_library(request)
    scheme = library.scheme
    counts: dict[str, int] = {}
    for entry in library.sorted_entries():
        if entry.call.part_type:
            counts[entry.call.part_type] = counts.get(entry.call.part_type, 0) + 1

    names = list(scheme.atoms) + [w for _, _, w in scheme.merges] + ["234", "678"]
    seen, out = set(), []
    for name in sorted(set(names) | set(counts), key=lambda n: (len(n), n)):
        if name in seen:
            continue
        seen.add(name)
        bg, fg = badge_for(name, scheme)
        out.append({
            "name": name,
            "description": DESCRIPTIONS.get(name, "composite part"),
            "count": counts.get(name, 0),
            "badge_bg": bg,
            "badge_fg": fg,
        })
    return out
