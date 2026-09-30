"""``/api/sequencing/*`` - check sequenced clones against what they should be.

Files arrive as text in the request body rather than as a multipart upload.
That is deliberate: multipart would add `python-multipart` to a tool whose
whole point is that it runs off a local folder with almost nothing installed,
and a browser can read a `.gb` and post its contents in three lines.

Like the Library screen, this is not an assembly level: it owns no `levelN`
endpoints and imports no level module.
"""

from __future__ import annotations

from typing import Any

from fastapi import APIRouter, HTTPException, Request
from pydantic import BaseModel

from ..core import align, seqio
from ..core.library import Library

router = APIRouter(prefix="/api/sequencing", tags=["sequencing"])

#: Past this the viewer's rows get heavy and the pairwise fallback gets slow.
#: A YTK multigene construct is ~15 kb, so this is well clear of real work.
MAX_LENGTH = 200_000

#: More clones than anyone picks off one plate in a single check.
MAX_CLONES = 48


def get_library(request: Request) -> Library:
    return request.app.state.library


class Upload(BaseModel):
    """One file the browser read, with the name it had on disk."""

    name: str = ""
    text: str = ""


class RunRequest(BaseModel):
    """What to check, and what to check it against."""

    reference: str = ""
    """A library plasmid's name, or a kept construct's name."""

    source: str = "library"
    """``library``, ``expected`` or ``upload``."""

    reference_file: Upload | None = None
    clones: list[Upload] = []


def _clone_name(upload: Upload, n: int) -> str:
    """The name to show for a clone: the file's, less the noise.

    Vendors name files after the tube, and the tube is what is written on the
    person's plate - so the file name is the most useful identifier there is,
    once the suffix and the vendor's own decoration are off it.
    """
    stem = (upload.name or f"clone {n}").rsplit("/", 1)[-1]
    for suffix in (".gb", ".gbk", ".genbank", ".fa", ".fasta", ".fna", ".seq", ".txt", ".dna"):
        if stem.lower().endswith(suffix):
            stem = stem[: -len(suffix)]
            break
    for tail in ("_consensus", ".consensus", "_assembly", "_plasmid"):
        if stem.lower().endswith(tail):
            stem = stem[: -len(tail)]
    return stem.strip() or f"clone {n}"


def _reference(library: Library, body: RunRequest) -> tuple[str, str, Any]:
    """The expected sequence: `(name, sequence, record or None)`."""
    if body.source == "upload":
        if not body.reference_file or not body.reference_file.text.strip():
            raise HTTPException(422, "no reference file was given")
        try:
            records = seqio.read_text_records(
                body.reference_file.text, body.reference_file.name or "reference"
            )
        except ValueError as error:
            raise HTTPException(422, str(error)) from None
        record = records[0]
        name = body.reference_file.name or record.id or "reference"
        return _clone_name(Upload(name=name), 0), seqio.sequence(record), record

    if not body.reference:
        raise HTTPException(422, "no reference was chosen")

    if body.source == "expected":
        record = library.expected_record(body.reference)
        if record is None:
            raise HTTPException(404, f"no kept construct called {body.reference}")
        return body.reference, seqio.sequence(record), record

    entry = library.get(body.reference)
    if entry is None:
        raise HTTPException(404, f"no plasmid called {body.reference}")
    record = library.record(entry)
    return entry.name, seqio.sequence(record), record


@router.get("/references")
def references(request: Request) -> dict[str, Any]:
    """Everything a clone can be checked against, in the two places they live."""
    library = get_library(request)
    return {
        "library": [
            {
                "name": e.name,
                "length": e.length,
                "part_type": e.call.part_type,
                "roles": list(e.roles),
                "component": e.component,
                "display": e.display,
            }
            # the deduplicated view: two files holding one sequence are one
            # thing to check a clone against, under the name you chose for it
            for e in sorted(library.unique_entries(), key=lambda e: e.name.lower())
        ],
        "expected": [
            {"name": d["name"], "length": d["length"]} for d in library.expected()
        ],
    }


def _difference(d: align.Difference) -> dict[str, Any]:
    return {
        "kind": d.kind,
        "start": d.start,
        "end": d.end,
        "length": d.length,
        "where": d.where,
        "expected": d.expected,
        "found": d.found,
        "regions": list(d.regions),
        "text": d.describe(),
    }


@router.post("/run")
def run(request: Request, body: RunRequest) -> dict[str, Any]:
    """Place every clone on the reference and report what differs.

    Each clone is compared independently and then merged into one set of
    columns, so a row on the screen means the same reference base on every
    clone - which is the only way two clones can be read against each other.
    """
    library = get_library(request)
    name, reference, record = _reference(library, body)

    if not reference:
        raise HTTPException(422, f"{name} has no sequence")
    if len(reference) > MAX_LENGTH:
        raise HTTPException(422, f"{name} is {len(reference):,} bp; the limit is {MAX_LENGTH:,}")
    if not body.clones:
        raise HTTPException(422, "no sequencing files were given")
    if len(body.clones) > MAX_CLONES:
        raise HTTPException(422, f"{len(body.clones)} files; the limit is {MAX_CLONES}")

    regions = align.regions_from_record(record) if record is not None else []

    results: list[align.CloneAlignment] = []
    rejected: list[dict[str, str]] = []
    for n, upload in enumerate(body.clones, start=1):
        label = _clone_name(upload, n)
        try:
            records = seqio.read_text_records(upload.text, label)
        except ValueError as error:
            # the parser names the file so its message stands alone; here the
            # name is already the row, so saying it twice is just noise
            why = str(error)
            rejected.append({"name": label, "why": why.removeprefix(f"{label}: ")})
            continue
        for index, parsed in enumerate(records):
            sequence = seqio.sequence(parsed)
            if not sequence:
                rejected.append({"name": label, "why": "the file holds no sequence"})
                continue
            if len(sequence) > MAX_LENGTH:
                rejected.append({"name": label, "why": f"{len(sequence):,} bp is too long"})
                continue
            # one file can hold several records; only then is the record's own
            # name needed to tell them apart
            title = label if len(records) == 1 else f"{label} · {parsed.id or index + 1}"
            results.append(align.compare(reference, sequence, name=title, regions=regions))

    placed = [r for r in results if r.placement.placed]
    merged = align.rows(reference, placed)

    return {
        "reference": {
            "name": name,
            "length": len(reference),
            "source": body.source,
            "regions": [
                {"start": r.start, "end": r.end, "label": r.label, "kind": r.kind}
                for r in regions
            ],
        },
        "clones": [
            {
                "name": r.name,
                "length": r.length,
                "placed": r.placement.placed,
                "strand": r.placement.strand,
                "offset": r.placement.offset,
                "votes": r.placement.votes,
                "identity": r.identity(reference),
                "clean": r.clean,
                "verdict": r.verdict(),
                "note": r.note,
                "changed_bases": r.changed_bases,
                "length_difference": r.length - len(reference),
                "differences": [_difference(d) for d in r.differences],
            }
            for r in results
        ],
        "alignment": {
            "columns": merged.columns,
            "reference": merged.reference,
            "rows": [{"name": n, "row": row} for n, row in merged.rows],
            "insertions": [list(pair) for pair in merged.insertions],
        },
        "rejected": rejected,
        "summary": _summary(results, rejected, len(reference)),
    }


def _summary(
    results: list[align.CloneAlignment], rejected: list[dict[str, str]], length: int
) -> dict[str, Any]:
    """The sentence at the top: how many clones are usable, and which.

    Named rather than counted. After a screen of eight clones the only thing
    anyone wants is which tube to grow up, and a count does not answer that.
    """
    clean = [r.name for r in results if r.clean]
    changed = [r.name for r in results if r.placement.placed and not r.clean]
    unplaced = [r.name for r in results if not r.placement.placed]

    if clean:
        headline = f"{', '.join(clean)} {'has' if len(clean) == 1 else 'have'} no mismatches"
    elif results and not unplaced:
        headline = "no clone matches the reference exactly"
    elif results:
        headline = "no clone is clean, and some do not match this plasmid at all"
    else:
        headline = "nothing could be read"

    return {
        "headline": headline,
        "clean": clean,
        "changed": changed,
        "unplaced": unplaced,
        "rejected": [r["name"] for r in rejected],
        "checked": len(results),
        "reference_length": length,
    }
