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
from ..core.codons import BackTranslationError, back_translate, clean_protein, translate
from ..core.enzymes import get_enzyme
from ..core.library import Library

router = APIRouter(prefix="/api/uniprot", tags=["uniprot"])

#: Sites a synthetic part must not contain: the scheme's own enzymes, plus the
#: alternates and the linearizer, since a part carrying one cannot be cloned
#: or integrated later even though it assembles today.
def _avoid(library: Library) -> tuple:
    scheme = library.scheme
    return (
        scheme.part_enzyme,
        scheme.multigene_enzyme,
        scheme.linearizer,
        *scheme.alternate_cloning_enzymes,
    )


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

    avoid = _avoid(library)
    try:
        dna, compromised = back_translate(
            protein, avoid=avoid, stop="TAA" if body.stop else None
        )
    except BackTranslationError as error:
        raise HTTPException(422, str(error)) from None

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
        "compromised": compromised,
        "stop_added": bool(body.stop and not protein.endswith("*")),
        # the round trip is the only check that matters, so it is made here
        # rather than trusted: a part that does not spell its protein is not a
        # part, however clean its sites are
        "verified": translate(dna).rstrip("*") == protein.rstrip("*"),
    }
