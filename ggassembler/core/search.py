"""Find a plasmid by what is annotated inside it, not by what it is called.

A library of several hundred plasmids is unusable through a dropdown: you know
you want the TDH3 promoter, you do not know that it is filed as pYTK009. So the
question this module answers is "which plasmids have *TDH* in them", and the
answer has to distinguish two very different kinds of hit:

* the annotation sits **inside the fragment the plasmid contributes** - the
  part that would actually end up in your construct. That is the real hit.
* the annotation sits in the **backbone**, outside that fragment. Every plasmid
  in the kit carries AmpR and a ColE1 origin; matching those is noise, and
  ranking them alongside a genuine part would bury it.

So each hit records `where` it matched and is scored accordingly. Nothing here
renames anything or guesses: the text searched is the user's own labels.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

from .library import Library, PlasmidEntry, clean_label, features_within

#: Where a match landed, best first. `part` and `cassette` are the fragment the
#: plasmid contributes; `backbone` is everything else in the file.
PART, CASSETTE, COMPONENT, NAME, BACKBONE = (
    "part", "cassette", "component", "name", "backbone",
)

#: How much a match in each place is worth. A backbone hit is kept - you may
#: genuinely be looking for the plasmid carrying KanR - but it cannot outrank a
#: part whose own sequence is what you asked for.
WEIGHTS = {PART: 10.0, CASSETTE: 9.0, COMPONENT: 8.0, NAME: 6.0, BACKBONE: 1.0}

#: Annotations every plasmid carries, which match far too much to be useful on
#: their own. They still score, just far below anything specific.
_UBIQUITOUS = re.compile(
    r"^(ampr|kanr|cmr|specr|tetr|cole1|pmb1|ori|f1 ori|rep_origin|"
    r"bom|rop|t0|t1|lac ?operator|lacz|m13|bla)\b",
    re.IGNORECASE,
)

_SQUASH = re.compile(r"[^a-z0-9]+")


def _squash(text: str) -> str:
    return _SQUASH.sub("", text.lower())


def terms(query: str) -> list[str]:
    """The query split into the terms a hit must satisfy, all of them."""
    return [t for t in query.lower().split() if t]


def score_term(term: str, text: str) -> float:
    """How well one term matches one piece of text: 0 when it does not.

    An exact match beats a match at the start of a word, which beats a match
    buried mid-word, because ``TDH3`` typed in full should put the TDH3
    promoter above a plasmid whose note happens to contain ``ScTDH3prom-v2``.
    """
    if not term or not text:
        return 0.0
    lowered = text.lower()
    if lowered == term:
        return 3.0

    partial = 0.0
    start = 0
    while (i := lowered.find(term, start)) != -1:
        if i == 0 or not lowered[i - 1].isalnum():
            return 2.0
        # A camelCase break reads as a word boundary, which matters because
        # half this library is named for a species: searching `tdh` has to find
        # `ScTDH3 Promoter`, and `scarlet` has to find `mScarlet-I`.
        if text[i].isupper() and text[i - 1].islower():
            return 2.0
        partial = 1.0
        start = i + 1
    if partial:
        return partial

    # last resort: ignore punctuation, so "tdh3" finds "TDH-3" and "pTDH_3"
    squashed = _squash(text)
    if _squash(term) and _squash(term) in squashed:
        return 0.75
    return 0.0


@dataclass
class Hit:
    """One plasmid the query matched, and the evidence for it."""

    name: str
    display: str
    path: str
    component: str
    part_type: str | None
    roles: list[str]
    length: int
    score: float
    where: str
    """The best place the query matched: `part`, `cassette`, `component`,
    `name` or `backbone`."""
    matched: list[str] = field(default_factory=list)
    """The labels that matched, in the order they appear in the plasmid."""
    aliases: list[str] = field(default_factory=list)
    five_prime: str | None = None
    three_prime: str | None = None
    connector_overhang: str | None = None
    cassette_overhangs: list[str] | None = None
    internal_sites: int = 0
    component_source: str = ""
    in_part: bool = False
    """True when at least one matched label lies in the contributed fragment."""
    usable_as_part: bool = False
    """The part enzyme really does release a part from this plasmid."""
    usable_as_unit: bool = False
    """The multigene enzyme reads it as a finished cassette or a destination."""


def usable_as_part(entry: PlasmidEntry) -> bool:
    """Whether a part could actually be cut out of this plasmid and used.

    Being *called* a type is not enough. The fragment has to exist: a finished
    multigene construct still carries the labels of the parts that built it,
    and a dropout vector's fragment is the piece the reaction throws away.
    """
    return bool(entry.is_part and entry.part_span)


def usable_as_unit(entry: PlasmidEntry) -> bool:
    """Whether this is a finished cassette or a multigene destination vector."""
    return bool(entry.is_cassette or entry.is_multigene_vector)


#: What each screen can actually put to work, so one endpoint serves them all.
USES = {"part": usable_as_part, "multigene": usable_as_unit}


def _searchable(entry: PlasmidEntry) -> tuple[list[str], list[str]]:
    """The labels inside the contributed fragment, and those outside it.

    Labels are cleaned of pipeline boilerplate first, so a search for ``XYL2``
    reports the hit as ``XYL2`` rather than ``XYL2_amplicon_BsaI_largest``.
    """
    labels = [(f, clean_label(f.label)) for f in entry.features if f.label]
    span = entry.part_span or entry.cassette_span
    if not span:
        return [], [label for _, label in labels]
    inside = {id(f) for f in features_within(entry.features, span, entry.length)}
    within = [label for f, label in labels if id(f) in inside]
    without = [label for f, label in labels if id(f) not in inside]
    return within, without


def score_entry(entry: PlasmidEntry, wanted: list[str]) -> Hit | None:
    """Score one plasmid against every term, or None if a term never matched.

    Every term has to land somewhere - searching ``tdh terminator`` should not
    return the TDH3 promoter just because half the query hit.
    """
    if not wanted:
        return None

    within, without = _searchable(entry)
    fragment_kind = CASSETTE if (entry.cassette_span and not entry.part_span) else PART
    names = [entry.name, *entry.aliases]

    # a description read off the file name is worth no more than the name it
    # was read from, so it is scored as one
    component_place = NAME if entry.component_source == "filename" else COMPONENT

    fields: list[tuple[str, list[str]]] = [
        (fragment_kind, within),
        (component_place, [entry.component] if entry.component else []),
        (NAME, names),
        (BACKBONE, without),
    ]

    total = 0.0
    best_place, best_value = BACKBONE, -1.0
    by_place: dict[str, list[str]] = {}
    seen: set[str] = set()

    for term in wanted:
        term_best = 0.0
        term_place = BACKBONE
        for place, texts in fields:
            weight = WEIGHTS[place]
            for text in texts:
                raw = score_term(term, text)
                if not raw:
                    continue
                if place in (PART, CASSETTE, BACKBONE) and _UBIQUITOUS.match(text):
                    raw *= 0.2
                value = raw * weight
                if value > term_best:
                    term_best, term_place = value, place
                if text not in seen:
                    seen.add(text)
                    by_place.setdefault(place, []).append(text)
        if term_best == 0.0:
            return None  # this term matched nothing: the plasmid is not a hit
        total += term_best
        if term_best > best_value:
            best_value, best_place = term_best, term_place

    # Report only the evidence worth reading. A search for `tdh terminator`
    # matches the backbone's own CamR terminator in every plasmid in the kit;
    # listing that beside the real hit reads as though it were part of it.
    order = [best_place] + [p for p in (PART, CASSETTE, COMPONENT, NAME) if p != best_place]
    matched = [text for place in order for text in by_place.get(place, [])]
    if not matched:
        matched = by_place.get(BACKBONE, [])

    return Hit(
        name=entry.name,
        display=entry.display,
        path=entry.path,
        component=entry.component,
        part_type=entry.call.part_type if not entry.call.reversed_sites else None,
        roles=list(entry.roles),
        length=entry.length,
        score=round(total, 3),
        where=best_place,
        matched=matched[:6],
        aliases=list(entry.aliases),
        five_prime=entry.call.five_prime,
        three_prime=entry.call.three_prime,
        connector_overhang=entry.connector_overhang,
        cassette_overhangs=(
            list(entry.cassette_overhangs) if entry.cassette_overhangs else None
        ),
        internal_sites=entry.sites.part_enzyme_internal,
        component_source=entry.component_source,
        in_part=best_place in (PART, CASSETTE, COMPONENT),
        usable_as_part=usable_as_part(entry),
        usable_as_unit=usable_as_unit(entry),
    )


def search(
    library: Library,
    query: str,
    limit: int = 25,
    part_type: str | None = None,
    role: str | None = None,
    usable: str | None = None,
) -> tuple[list[Hit], int]:
    """Plasmids matching `query`, best first, and how many were held back.

    `part_type`, `role` and `usable` narrow the result to what a particular
    screen could take, so the same search serves all of them. `usable` is the
    one that matters in practice: a screen that offers a plasmid it cannot use
    is worse than one that offers nothing, so the count of suppressed matches
    comes back too and can be reported rather than hidden.
    """
    wanted = terms(query)
    if not wanted:
        return [], 0

    fits = USES.get(usable or "")
    hits: list[Hit] = []
    suppressed = 0
    for entry in library.unique_entries():
        if part_type and entry.call.part_type != part_type:
            continue
        if role and role not in entry.roles:
            continue
        hit = score_entry(entry, wanted)
        if hit is None:
            continue
        if fits and not fits(entry):
            suppressed += 1
            continue
        hits.append(hit)

    # score first; then a real part ahead of an untyped file; then the shorter
    # plasmid, which is the more specific one; then by name so it never wobbles
    hits.sort(key=lambda h: (-h.score, not h.in_part, h.part_type is None, h.length, h.name))
    return hits[:limit], suppressed
