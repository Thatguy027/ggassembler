"""``/api/uniprot/*`` - find a protein, and turn it into a part.

The only part of this app that talks to the network. It touches `core` and
nothing else, and every route names a UniProt failure as a UniProt failure
rather than a server error, because being offline is an ordinary thing for a
tool that otherwise runs entirely off a local folder.
"""

from __future__ import annotations

from typing import Any

from fastapi import APIRouter, HTTPException, Request
from pydantic import BaseModel

from ..core import uniprot
from ..core.codons import (
    TABLE_NAMES,
    TABLES,
    BackTranslationError,
    back_translate,
    clean_protein,
    get_table,
    translate,
)
from ..core.enzymes import ENZYMES, get_enzyme
from ..core.library import Library

router = APIRouter(prefix="/api/uniprot", tags=["uniprot"])

def default_enzymes(library: Library) -> tuple:
    """What a synthetic part should avoid unless told otherwise.

    The scheme's own enzymes, the alternates and the linearizer - a part
    carrying any of them assembles today and cannot be cloned or integrated
    later, which is the expensive kind of mistake to find afterwards.
    """
    scheme = library.scheme
    return (
        scheme.part_enzyme,
        scheme.multigene_enzyme,
        scheme.linearizer,
        *scheme.alternate_cloning_enzymes,
    )


def _resolve(library: Library, names: list[str] | None) -> tuple:
    """The enzymes to avoid: the caller's list, or the scheme's if none given.

    An empty list is a choice, not an omission - it means avoid nothing - so
    `None` and `[]` are deliberately different.
    """
    if names is None:
        return default_enzymes(library)
    try:
        return tuple(get_enzyme(name) for name in names)
    except ValueError as error:
        raise HTTPException(422, str(error)) from None


def get_library(request: Request) -> Library:
    return request.app.state.library


def _entry(found: uniprot.Entry) -> dict[str, Any]:
    return {
        "accession": found.accession,
        "name": found.name,
        "protein": found.protein,
        "organism": found.organism,
        "genes": found.genes,
        "length": found.length,
        "reviewed": found.reviewed,
        "label": found.label,
    }


@router.get("/search")
def search(
    request: Request,
    q: str = "",
    organism: str | None = None,
    limit: int = 20,
) -> dict[str, Any]:
    """Proteins matching `q`, reviewed entries first."""
    try:
        found = uniprot.search(q, limit=limit, organism=organism or None)
    except uniprot.UniProtError as error:
        raise HTTPException(503, str(error)) from None
    return {"query": q, "count": len(found), "results": [_entry(e) for e in found]}


class PartRequest(BaseModel):
    """A protein to write as DNA - from UniProt, or pasted straight in."""

    accession: str = ""
    protein: str = ""
    stop: bool = True
    name: str = ""
    codon_table: str = "scerevisiae"
    avoid: list[str] | None = None
    """Enzyme names to keep out of the sequence. `None` means the scheme's
    own; `[]` means avoid nothing, which is a choice and not an omission."""


@router.get("/options")
def options(request: Request) -> dict[str, Any]:
    """The codon tables and the enzymes a caller can pick between."""
    library = get_library(request)
    return {
        "codon_tables": [
            {"key": key, "name": TABLE_NAMES.get(key, key)} for key in sorted(TABLES)
        ],
        "enzymes": [
            {"name": e.name, "site": e.site, "default": True}
            for e in default_enzymes(library)
        ],
        "other_enzymes": [
            {"name": e.name, "site": e.site, "default": False}
            for e in ENZYMES.values()
            if e.name not in {x.name for x in default_enzymes(library)}
        ],
    }


@router.post("/part")
def part(request: Request, body: PartRequest) -> dict[str, Any]:
    """Write a protein as a coding sequence that carries none of the kit's sites.

    Chosen codon by codon rather than domesticated afterwards: a synonymous
    change made to remove one site can spell another two codons along, whereas
    a codon that would spell a site is simply never used.
    """
    library = get_library(request)

    found = None
    protein = clean_protein(body.protein)
    if body.accession:
        try:
            found = uniprot.entry(body.accession)
        except uniprot.UniProtError as error:
            raise HTTPException(503, str(error)) from None
        protein = found.sequence
    if not protein:
        raise HTTPException(422, "give an accession or a protein sequence")

    avoid = _resolve(library, body.avoid)
    try:
        table = get_table(body.codon_table)
    except ValueError as error:
        raise HTTPException(422, str(error)) from None
    try:
        dna, compromised = back_translate(
            protein, avoid=avoid, stop="TAA" if body.stop else None, table=table
        )
    except BackTranslationError as error:
        raise HTTPException(422, str(error)) from None

    gc = sum(dna.count(base) for base in "GC")

    return {
        "name": body.name or (found.genes[0] if found and found.genes else "")
                or (found.accession if found else "part"),
        "accession": found.accession if found else "",
        "source": _entry(found) if found else None,
        "protein": protein,
        "protein_length": len(protein.rstrip("*")),
        "dna": dna,
        "length": len(dna),
        "avoided": [e.name for e in avoid],
        "codon_table": body.codon_table,
        "codon_table_name": TABLE_NAMES.get(body.codon_table, body.codon_table),
        "gc": round(100 * gc / len(dna), 1) if dna else 0.0,
        "compromised": compromised,
        "stop_added": bool(body.stop and not protein.endswith("*")),
        # the round trip is the only check that matters, so it is made here
        # rather than trusted: a part that does not spell its protein is not a
        # part, however clean its sites are
        "verified": translate(dna).rstrip("*") == protein.rstrip("*"),
    }
