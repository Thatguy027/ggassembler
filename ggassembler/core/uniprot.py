"""A thin read-only client for the UniProt REST API.

The rest of this app is offline by construction: it reads GenBank files off a
disk and never asks anything of the network. This module is the one exception,
and it is deliberately small because of that - search, fetch one entry, and
nothing else. Everything it returns is plain data that the rest of the app
treats like any other input.

The network call is injected rather than called directly, so the tests that
exercise the parsing never touch the internet and never depend on UniProt
being up or on its result for a given query staying the same. What the tests
cannot check is whether the live service still answers in this shape; that is
what `UniProtError` is for.

No new dependency: `urllib` from the standard library is enough for two GETs.
"""

from __future__ import annotations

import json
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass, field
from typing import Any, Callable

BASE = "https://rest.uniprot.org/uniprotkb"
TIMEOUT = 20.0

#: What a search asks for. Kept short: every field is one more thing to parse
#: and the screen shows five of them.
FIELDS = "accession,id,protein_name,gene_names,organism_name,length,reviewed"


class UniProtError(RuntimeError):
    """UniProt could not be reached, or did not answer in the expected shape."""


@dataclass
class Entry:
    """One UniProt record, reduced to what a part designer needs."""

    accession: str
    name: str
    protein: str
    organism: str
    genes: list[str] = field(default_factory=list)
    length: int = 0
    reviewed: bool = False
    sequence: str = ""

    @property
    def label(self) -> str:
        """`XKS1 - Xylulokinase`, the way a person refers to it."""
        gene = self.genes[0] if self.genes else ""
        return f"{gene} - {self.protein}" if gene else self.protein


def _get(url: str) -> Any:
    """One GET, returning parsed JSON, with every failure named the same way."""
    request = urllib.request.Request(url, headers={"Accept": "application/json"})
    try:
        with urllib.request.urlopen(request, timeout=TIMEOUT) as response:
            return json.load(response)
    except urllib.error.HTTPError as error:
        raise UniProtError(f"UniProt returned {error.code} for {url}") from None
    except urllib.error.URLError as error:
        raise UniProtError(f"could not reach UniProt: {error.reason}") from None
    except (TimeoutError, OSError) as error:
        raise UniProtError(f"could not reach UniProt: {error}") from None
    except json.JSONDecodeError:
        raise UniProtError("UniProt did not return JSON") from None


Fetcher = Callable[[str], Any]


def parse_entry(raw: dict[str, Any]) -> Entry:
    """One search result, as far as it can be trusted to be there.

    Written defensively on purpose: UniProt omits fields rather than nulling
    them, and a missing recommended name is normal for an unreviewed entry.
    """
    description = raw.get("proteinDescription") or {}
    named = description.get("recommendedName") or description.get("submissionNames") or {}
    if isinstance(named, list):
        named = named[0] if named else {}
    protein = ((named.get("fullName") or {}).get("value")) or ""

    genes = [
        (gene.get("geneName") or {}).get("value", "")
        for gene in raw.get("genes") or []
    ]
    return Entry(
        accession=raw.get("primaryAccession", ""),
        name=raw.get("uniProtkbId", ""),
        protein=protein or raw.get("primaryAccession", "unnamed"),
        organism=(raw.get("organism") or {}).get("scientificName", ""),
        genes=[g for g in genes if g],
        length=(raw.get("sequence") or {}).get("length", 0),
        reviewed=raw.get("entryType", "").startswith("UniProtKB reviewed"),
        sequence=(raw.get("sequence") or {}).get("value", ""),
    )


def search(
    query: str,
    limit: int = 20,
    organism: str | None = None,
    reviewed_first: bool = True,
    fetch: Fetcher = _get,
) -> list[Entry]:
    """Search UniProt, best first.

    `organism` takes a taxonomy id - 559292 is the *S. cerevisiae* reference
    strain - and narrows the query the way the database itself does rather than
    by filtering afterwards, so the result count means what it says.
    """
    query = query.strip()
    if not query:
        return []
    # kept before the organism clause is folded in, since the ranking below
    # compares it against gene names and accessions
    original = query
    if organism:
        query = f"({query}) AND organism_id:{organism}"

    url = f"{BASE}/search?" + urllib.parse.urlencode({
        "query": query,
        "fields": FIELDS,
        "format": "json",
        "size": str(max(1, min(limit, 100))),
    })
    entries = [parse_entry(raw) for raw in (fetch(url).get("results") or [])]
    if reviewed_first:
        entries.sort(key=lambda e: rank(e, original))
    return entries


def rank(entry: Entry, query: str) -> tuple:
    """Sort key: what was asked for first, then what has been curated.

    UniProt's own relevance put GRE3 above XKS1 for the query `XKS1`, which is
    defensible as text search and useless as an answer - someone typing a gene
    name has already decided which gene they want. So an exact match on the
    gene name or the accession outranks everything, then a gene name that
    merely starts with it, then reviewed entries, and only then whatever else
    came back.

    Length breaks the remaining ties, because among near-duplicates the shorter
    record is nearly always the reviewed canonical one rather than an isoform.
    """
    wanted = query.strip().lower()
    genes = [gene.lower() for gene in entry.genes]

    if wanted and (wanted == entry.accession.lower() or wanted in genes):
        exactness = 0
    elif wanted and any(gene.startswith(wanted) for gene in genes):
        exactness = 1
    elif wanted and entry.name.lower().startswith(wanted):
        exactness = 2
    else:
        exactness = 3

    return (exactness, not entry.reviewed, entry.length)


def entry(accession: str, fetch: Fetcher = _get) -> Entry:
    """One record in full, including the sequence a search does not return."""
    accession = accession.strip()
    if not accession:
        raise UniProtError("no accession given")
    raw = fetch(f"{BASE}/{urllib.parse.quote(accession)}?format=json")
    found = parse_entry(raw)
    if not found.sequence:
        raise UniProtError(f"{accession} came back without a sequence")
    return found
