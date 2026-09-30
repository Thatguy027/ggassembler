"""Level 2 - cassette (transcription unit) assembly.

Eight slots around the part circle, one part each, ligated by the scheme's part
enzyme into a cassette plasmid. The panel variations the UI offers - splitting
type 3 into 3a+3b, type 4 into 4a+4b, type 8 into 8a+8b, or collapsing 2-3-4
and 6-7-8 into single composite slots - are all **view modes over one slot
model**, not separate assembly paths. `slots()` returns whichever set of arcs
the current view asks for; everything downstream treats them identically.

This module owns no enzyme knowledge of its own and imports no other level.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from ..core import parttypes
from ..core.assembly import ERROR, INFO, WARNING, AssemblyResult, Issue, Piece, assemble
from ..core import seqio
from ..core.library import Library, PlasmidEntry
from ..core.parttypes import YTK, Scheme

#: The left column of the screen, then the right.
LEFT_POSITIONS = ("1", "2", "3", "4")
RIGHT_POSITIONS = ("5", "6", "7", "8")


@dataclass(frozen=True)
class Slot:
    """One dropdown on the Level 2 screen."""

    key: str
    """The part type this slot takes, e.g. ``3``, ``3a``, ``234``."""
    five_prime: str
    three_prime: str
    column: str
    """``left`` or ``right``."""

    @property
    def label(self) -> str:
        return self.key.replace("234", "2·4").replace("678", "6·8")

    @property
    def description(self) -> str:
        return parttypes.DESCRIPTIONS.get(self.key, "composite part")

    @property
    def overhangs(self) -> tuple[str, str]:
        return self.five_prime, self.three_prime


@dataclass
class CassetteDesign:
    """A Level 2 screen's state: which view modes are on, and what is chosen."""

    selections: dict[str, str] = field(default_factory=dict)
    """Slot key -> plasmid name."""
    split_3: bool = False
    """Show 3a + 3b instead of a whole type 3."""
    split_4: bool = False
    """Show 4a + 4b instead of a whole type 4."""
    split_8: bool = False
    """Show 8a + 8b - the chromosomal-integration configuration."""
    composite_left: bool = False
    """Collapse types 2, 3 and 4 into one 234 slot."""
    composite_right: bool = False
    """Collapse types 6, 7 and 8 into one 678 slot."""
    name: str = "cassette"

    @property
    def is_integration(self) -> bool:
        """8a + 8b means the construct integrates rather than replicating."""
        return self.split_8 and not self.composite_right


def slots(design: CassetteDesign, scheme: Scheme = YTK) -> list[Slot]:
    """The active slots, in circle order."""
    keys: list[str] = ["1"]
    if design.composite_left:
        keys.append("234")
    else:
        keys.append("2")
        keys.extend(("3a", "3b") if design.split_3 else ("3",))
        keys.extend(("4a", "4b") if design.split_4 else ("4",))
    keys.append("5")
    if design.composite_right:
        keys.append("678")
    else:
        keys.append("6")
        keys.append("7")
        keys.extend(("8a", "8b") if design.split_8 else ("8",))

    out = []
    for key in keys:
        pair = parttypes.type_overhangs(key, scheme)
        if pair is None:
            continue
        column = "left" if key in ("1", "2", "234") or key[0] in "34" else "right"
        out.append(Slot(key=key, five_prime=pair[0], three_prime=pair[1], column=column))
    return out


def options(library: Library, slot: Slot) -> list[PlasmidEntry]:
    """Library parts whose released overhangs match this slot exactly."""
    return library.parts_with_overhangs(*slot.overhangs)


def option_counts(library: Library, design: CassetteDesign) -> dict[str, int]:
    """How many parts each slot can offer - the count the panel shows."""
    return {slot.key: len(options(library, slot)) for slot in slots(design, library.scheme)}


def default_design(library: Library, name: str = "pCassette") -> CassetteDesign:
    """A complete, working cassette to open the screen with.

    An empty screen teaches nothing: eight empty dropdowns do not show what a
    finished construct looks like, or that the ring is clickable. So the page
    opens on a real assembly picked from this library.

    Nothing is hard-coded to a catalogue: each slot takes the simplest part that
    fits it - no internal sites, an annotated description, a short name, and for
    the two connector slots an overhang not already used, so the cassette can
    actually be placed in a multigene assembly afterwards.
    """
    design = CassetteDesign(name=name)
    used_connectors: set[str] = set()

    for slot in slots(design, library.scheme):
        candidates = options(library, slot)
        if not candidates:
            continue
        best = min(candidates, key=lambda e: _preference(e, slot, used_connectors))
        design.selections[slot.key] = best.name
        if best.connector_overhang:
            used_connectors.add(best.connector_overhang)

    return design


def _preference(entry: PlasmidEntry, slot: Slot, used_connectors: set[str]) -> tuple:
    """How good a default a part makes: clean, described, distinct, simple."""
    wants_connector = slot.key in ("1", "5")
    return (
        bool(entry.sites.part_enzyme_internal),
        bool(entry.internal_multigene_positions),
        wants_connector and not entry.connector_overhang,
        bool(entry.connector_overhang and entry.connector_overhang in used_connectors),
        not entry.component,
        len(entry.name),
        entry.name,
    )


def build(library: Library, design: CassetteDesign) -> AssemblyResult:
    """Resolve the design against the library and simulate the reaction."""
    scheme = library.scheme
    active = slots(design, scheme)
    pieces: list[Piece] = []
    issues: list[Issue] = []
    chosen: list[PlasmidEntry] = []

    for slot in active:
        name = design.selections.get(slot.key)
        if not name:
            issues.append(
                Issue(ERROR, "empty_slot",
                      f"slot {slot.label} ({slot.five_prime} -> {slot.three_prime}) is empty")
            )
            continue

        entry = library.get(name)
        if entry is None:
            issues.append(Issue(ERROR, "unknown_part", f"no plasmid named {name} in the library"))
            continue

        if entry.overhangs != slot.overhangs:
            issues.append(
                Issue(ERROR, "wrong_overhangs",
                      f"{entry.name} releases {entry.call.five_prime} -> "
                      f"{entry.call.three_prime}, but slot {slot.label} needs "
                      f"{slot.five_prime} -> {slot.three_prime}")
            )
            continue

        if entry.call.reversed_sites:
            issues.append(
                Issue(ERROR, "dropout_as_part",
                      f"{entry.name} is a dropout ({entry.call.part_type}); it is consumed "
                      f"by the reaction and never ends up in the product")
            )
            continue

        try:
            pieces.append(library.piece(entry))
        except (ValueError, OSError) as exc:
            issues.append(Issue(ERROR, "unreadable_part", f"{entry.name}: {exc}"))
            continue
        chosen.append(entry)

    if any(i.level == ERROR for i in issues):
        return AssemblyResult(False, None, issues=issues, enzyme=scheme.part_enzyme.name)

    result = assemble(pieces, scheme.part_enzyme, name=design.name)
    result.issues = issues + result.issues
    result.issues += _notes(design, chosen)
    result.ok = result.ok and not result.errors
    if not result.ok:
        result.product = None
    return result


def _notes(design: CassetteDesign, chosen: list[PlasmidEntry]) -> list[Issue]:
    """Advice that is not a defect: screening colour, integration, marker."""
    notes: list[Issue] = []

    dropout = _dropout_colour(chosen)
    if dropout:
        colour, owner = dropout
        notes.append(
            Issue(INFO, "screening",
                  f"{owner} carries a {colour} dropout: screen "
                  f"{colour}/white colonies at the next level")
        )

    if design.is_integration:
        notes.append(
            Issue(INFO, "integration",
                  "8a + 8b is the chromosomal-integration configuration: the type 7 slot "
                  "should hold a 3' genomic homology arm, and the finished plasmid must be "
                  "linearized with NotI before transformation")
        )

    notes += _connector_notes(chosen)

    # Only the part covering position 8 contributes an E. coli marker to the
    # product; the others leave theirs behind on their own part plasmids, so
    # comparing all eight would warn about something that cannot happen.
    markers = {
        e.ecoli_marker
        for e in chosen
        if e.ecoli_marker and "8" in parttypes.positions(e.call.part_type or "")
    }
    if len(markers) > 1:
        notes.append(
            Issue(WARNING, "mixed_markers",
                  f"more than one backbone part carries an E. coli marker "
                  f"({', '.join(sorted(markers))}); only one can select the product")
        )

    return notes


def _connector_notes(chosen: list[PlasmidEntry]) -> list[Issue]:
    """The two connectors must differ, or the cassette cannot be ordered later.

    A Level 3 assembly is ordered by matching each cassette's ConR overhang to
    the next cassette's ConL. A cassette whose own two connectors carry the same
    overhang - ConL1 with ConR1, say - has no place in that chain: it would have
    to follow itself.
    """
    connectors = [
        (e, e.connector_overhang)
        for e in chosen
        if e.connector_overhang and e.call.part_type in ("1", "5")
    ]
    if len(connectors) != 2:
        return []
    (left, left_overhang), (right, right_overhang) = connectors
    if left_overhang != right_overhang:
        return []
    return [
        Issue(WARNING, "same_connector",
              f"{left.name} and {right.name} both carry connector overhang "
              f"{left_overhang}; the finished cassette would have the same connector "
              f"at both ends and could not be placed in a multigene assembly")
    ]


def _dropout_colour(chosen: list[PlasmidEntry]) -> tuple[str, str] | None:
    """Which colour dropout, if any, the assembled backbone carries."""
    palette = (("green", ("gfp", "venus", "sfgfp", "egfp")), ("red", ("rfp", "mrfp", "mcherry")))
    for entry in chosen:
        haystack = " ".join(f.label for f in entry.features).lower()
        for colour, needles in palette:
            if any(needle in haystack for needle in needles):
                return colour, entry.name
    return None


def validation_strip(result: AssemblyResult) -> list[dict[str, str]]:
    """The live validation strip, as rows the UI can render without re-checking."""
    return [
        {"level": issue.level, "code": issue.code, "message": issue.message}
        for issue in result.issues
    ]


# --------------------------------------------------------------------------- #
# reading a finished cassette back into parts
# --------------------------------------------------------------------------- #


@dataclass
class Match:
    """One library part found inside an assembled cassette."""

    name: str
    display: str
    part_type: str
    start: int
    end: int
    length: int
    component: str = ""


@dataclass
class Decomposition:
    """What a finished cassette is made of, as far as the library can tell."""

    cassette: str
    length: int
    matches: list[Match] = field(default_factory=list)
    selections: dict[str, str] = field(default_factory=dict)
    split_3: bool = False
    split_4: bool = False
    split_8: bool = False
    composite_left: bool = False
    composite_right: bool = False
    gaps: list[tuple[int, int]] = field(default_factory=list)
    """Stretches no library part explains, as (start, length)."""

    @property
    def covered(self) -> int:
        return sum(m.length for m in self.matches)

    @property
    def complete(self) -> bool:
        return bool(self.matches) and not self.gaps

    def to_design(self, name: str) -> CassetteDesign:
        return CassetteDesign(
            selections=dict(self.selections),
            split_3=self.split_3,
            split_4=self.split_4,
            split_8=self.split_8,
            composite_left=self.composite_left,
            composite_right=self.composite_right,
            name=name,
        )


def decompose(library: Library, cassette: PlasmidEntry) -> Decomposition:
    """Read an assembled cassette back into the library parts that built it.

    The assembly consumed every BsaI site, so this cannot be a digest: the
    parts are gone as *sites* but still present as *sequence*. So each part's
    released fragment is looked for in the cassette directly, and the ones that
    abut end-to-end are chained into a tiling of the plasmid.

    Where the chain covers the whole circle the panels can be filled exactly.
    Where it does not, the gaps come back as gaps rather than being papered
    over - a part swapped in from outside the library is a real answer, and
    guessing at it would be worse than saying so.
    """
    out = Decomposition(cassette=cassette.name, length=cassette.length)
    sequence = seqio.sequence(library.record(cassette)).upper()
    if not sequence:
        return out

    # search the doubled sequence so a part spanning the origin is still found
    doubled = sequence + sequence
    n = len(sequence)

    starts: dict[int, list[tuple[PlasmidEntry, int]]] = {}
    for entry, fragment in library.part_fragments():
        if len(fragment) > n:
            continue
        at = doubled.find(fragment)
        while at != -1 and at < n:
            starts.setdefault(at, []).append((entry, (at + len(fragment)) % n))
            at = doubled.find(fragment, at + 1)

    if not starts:
        out.gaps = [(0, n)]
        return out

    chain = _tile(starts, n)
    if not chain:
        out.gaps = [(0, n)]
        return out

    for entry, start, end in chain:
        out.matches.append(
            Match(
                name=entry.name,
                display=entry.display,
                part_type=entry.call.part_type or "",
                start=start,
                end=end,
                length=(end - start) % n or n,
                component=entry.component,
            )
        )

    out.gaps = _gaps(chain, n)
    _fill_selections(out)
    return out


def _tile(
    starts: dict[int, list[tuple[PlasmidEntry, int]]], n: int
) -> list[tuple[PlasmidEntry, int, int]]:
    """The longest run of parts that abut end-to-end, closing the circle if it can.

    Tried from every match in turn, preferring atomic parts over composites so
    the panels come back as granular as the library allows: a 2-3-4 dropout and
    three separate parts can describe the same stretch, and three panels are
    more use than one.
    """
    def rank(pair: tuple[PlasmidEntry, int]) -> tuple:
        entry, _ = pair
        atoms = len(parttypes.positions(entry.call.part_type or "", YTK))
        return (atoms, entry.name)

    best: list[tuple[PlasmidEntry, int, int]] = []
    for origin in sorted(starts):
        used: set[str] = set()
        chain: list[tuple[PlasmidEntry, int, int]] = []
        at = origin
        while True:
            candidates = [c for c in sorted(starts.get(at, []), key=rank)
                          if c[0].name not in used]
            if not candidates:
                break
            entry, end = candidates[0]
            used.add(entry.name)
            chain.append((entry, at, end))
            at = end
            if at == origin:
                break  # the circle closed
        covered = sum((e - s) % n or n for _, s, e in chain)
        if covered > sum((e - s) % n or n for _, s, e in best):
            best = chain
        if covered == n:
            break
    return best


def _gaps(chain: list[tuple[PlasmidEntry, int, int]], n: int) -> list[tuple[int, int]]:
    """The stretches the chain leaves unexplained."""
    if not chain:
        return [(0, n)]
    covered = sum((e - s) % n or n for _, s, e in chain)
    if covered >= n:
        return []
    # the chain is contiguous, so there is exactly one gap: end of last to start
    return [(chain[-1][2], (chain[0][1] - chain[-1][2]) % n)]


#: A composite key and the atomic keys it replaces on screen. Only one of the
#: two can be shown at once, so only one may be selected.
_SUPERSEDES = {
    "3": ("3a", "3b"),
    "4": ("4a", "4b"),
    "8": ("8a", "8b"),
    "234": ("2", "3", "3a", "3b", "4", "4a", "4b"),
    "678": ("6", "7", "8", "8a", "8b"),
}


def _fill_selections(out: Decomposition) -> None:
    """Turn the matched parts into panel choices and the view flags they need.

    A partial chain can turn up two parts whose positions overlap - a whole
    type 8 from one stretch and an 8b from another - and the screen has no way
    to show both. The finer choice wins, because it is the one that says more,
    and the coarser key is dropped rather than left to contradict it.
    """
    for match in out.matches:
        if match.part_type:
            out.selections[match.part_type] = match.name

    for coarse, finer in _SUPERSEDES.items():
        if coarse in out.selections and any(k in out.selections for k in finer):
            del out.selections[coarse]

    keys = set(out.selections)
    out.split_3 = bool(keys & {"3a", "3b"})
    out.split_4 = bool(keys & {"4a", "4b"})
    out.split_8 = bool(keys & {"8a", "8b"})
    out.composite_left = "234" in keys
    out.composite_right = "678" in keys
