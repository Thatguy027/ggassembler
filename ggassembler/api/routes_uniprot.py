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
    feasibility,
    get_table,
    translate,
)
from ..core.enzymes import ENZYMES, Enzyme, get_enzyme, IUPAC
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


def as_enzyme(text: str) -> Enzyme:
    """One thing to keep out of a sequence, named or spelled.

    A name when the app knows it, otherwise the site itself in IUPAC - which is
    the only way to ask for something outside the four-enzyme list without this
    module growing a copy of REBASE. The cut offsets are nominal: nothing here
    digests anything, it only refuses to spell the recognition sequence.
    """
    text = text.strip()
    if not text:
        raise ValueError("no enzyme or site given")
    try:
        return get_enzyme(text)
    except ValueError:
        pass

    site = text.upper().replace("U", "T")
    bad = sorted({c for c in site if c not in IUPAC})
    if bad:
        raise ValueError(
            f"{text!r} is not an enzyme this app knows, and not a site either "
            f"(found {', '.join(bad)})"
        )
    if len(site) < 4:
        raise ValueError(f"{text!r} is too short to be a useful site")
    return Enzyme(site, site, top_offset=0, bottom_offset=len(site))


def _resolve(library: Library, names: list[str] | None) -> tuple:
    """The enzymes to avoid: the caller's list, or the scheme's if none given.

    An empty list is a choice, not an omission - it means avoid nothing - so
    `None` and `[]` are deliberately different.
    """
    if names is None:
        return default_enzymes(library)
    try:
        return tuple(as_enzyme(name) for name in names)
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


#: What each target does about the stop codon, and why.
#:
#: This is the whole reason the screen has to ask. A Type 3 part *omits* its
#: stop - the Type 4 terminator supplies TAA straight after the ATCC overhang,
#: and omitting it is what makes read-through into a Type 4a fusion possible at
#: all. Handing Level 1 a sequence with a stop on it, for a part type whose
#: convention is to remove one, means the two screens immediately undo each
#: other's work.
PART_TARGETS: dict[str, dict[str, object]] = {
    "3": {
        "label": "Type 3 - coding sequence",
        "stop": False,
        "why": "a Type 3 part omits the stop; the Type 4 terminator supplies TAA "
               "right after the ATCC overhang",
    },
    "3a": {
        "label": "Type 3a - N-terminal half of a split CDS",
        "stop": False,
        "why": "3a reads through into 3b, so it cannot carry a stop",
    },
    "3b": {
        "label": "Type 3b - C-terminal half of a split CDS",
        "stop": False,
        "why": "like Type 3: the terminator part supplies the stop",
    },
    "4a": {
        "label": "Type 4a - C-terminal fusion",
        "stop": True,
        "why": "a 4a part ends the protein itself, with its own TAA and XhoI, "
               "rather than reading through its TGGC flank",
    },
    "": {
        "label": "A gene to order - no part type",
        "stop": True,
        "why": "a plain open reading frame, ending where you expect it to",
    },
}


class PartRequest(BaseModel):
    """A protein to write as DNA - from UniProt, or pasted straight in."""

    accession: str = ""
    protein: str = ""
    name: str = ""
    part_type: str = ""
    """Which YTK position this is for. The type decides the stop codon; see
    `PART_TARGETS`. Empty means a plain gene to order."""

    stop: bool | None = None
    """Overrides the part type's own answer. `None` - the normal case - means
    follow the convention rather than second-guess it."""

    codon_table: str = "scerevisiae"
    avoid: list[str] | None = None
    """Enzyme names, or IUPAC site strings, to keep out of the sequence. `None`
    means the scheme's own; `[]` means avoid nothing, which is a choice and not
    an omission."""


@router.get("/options")
def options(request: Request) -> dict[str, Any]:
    """The codon tables, part types and enzymes a caller can pick between."""
    library = get_library(request)
    return {
        "part_types": [
            {"key": key, "label": spec["label"], "stop": spec["stop"], "why": spec["why"]}
            for key, spec in PART_TARGETS.items()
        ],
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

    target = PART_TARGETS.get(body.part_type)
    if target is None:
        raise HTTPException(
            422,
            f"unknown part type: {body.part_type}. "
            f"Known: {', '.join(k or '(none)' for k in PART_TARGETS)}",
        )
    # The part type decides the stop unless the caller overrode it. Following
    # the convention is the point: a Type 3 sequence handed to Level 1 with a
    # stop on it is a sequence Level 1's own conventions will strip straight
    # back off, and the two screens disagreeing about the same part is worse
    # than either being wrong on its own.
    wants_stop = bool(target["stop"]) if body.stop is None else bool(body.stop)

    avoid = _resolve(library, body.avoid)
    try:
        table = get_table(body.codon_table)
    except ValueError as error:
        raise HTTPException(422, str(error)) from None
    try:
        dna, compromised = back_translate(
            protein, avoid=avoid, stop="TAA" if wants_stop else None, table=table
        )
    except BackTranslationError as error:
        raise HTTPException(422, str(error)) from None

    # A protein handed over with its own terminal stop, for a part type whose
    # convention is not to carry one, has to lose it - otherwise the "no stop"
    # answer depends on what the source record happened to include.
    if not wants_stop and translate(dna).endswith("*"):
        dna = dna[:-3]

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
        "part_type": body.part_type,
        "part_type_label": target["label"],
        "stop": wants_stop,
        "stop_reason": target["why"],
        "stop_followed_convention": body.stop is None,
        "stop_added": bool(wants_stop and not protein.endswith("*")),
        "stop_removed": bool(not wants_stop and protein.rstrip().endswith("*")),
        # will a vendor make it? none of this shows up in a translation check
        "feasibility": feasibility(dna),
        # the round trip is the only check that matters, so it is made here
        # rather than trusted: a part that does not spell its protein is not a
        # part, however clean its sites are
        "verified": translate(dna).rstrip("*") == protein.rstrip("*"),
    }
