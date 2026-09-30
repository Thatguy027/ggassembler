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

from dataclasses import dataclass, field, replace
from typing import Sequence

from ..core import cassette as cassette_core
from ..core import decompose as decompose_core
from ..core import parttypes
from ..core.assembly import ERROR, INFO, WARNING, AssemblyResult, Issue, Piece, assemble
from ..core import seqio
from ..core.library import Library, PlasmidEntry
from ..core.parttypes import YTK, Scheme

#: The left column of the screen, then the right.
LEFT_POSITIONS = ("1", "2", "3", "4")
RIGHT_POSITIONS = ("5", "6", "7", "8")


#: One position of the part circle. Defined in `core.cassette` because the
#: Plate screen builds ninety-six cassettes and may not import this module.
Slot = cassette_core.Slot


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
    return cassette_core.slots_for(
        cassette_core.slot_keys(
            split_3=design.split_3,
            split_4=design.split_4,
            split_8=design.split_8,
            composite_left=design.composite_left,
            composite_right=design.composite_right,
        ),
        scheme,
    )


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
    """Resolve the design against the library and simulate the reaction.

    The resolution and the reaction live in `core.cassette`, shared with the
    Plate screen. What stays here is the advice - screening colour, marker,
    integration - which is about how this screen reads, not about the kit.
    """
    result, chosen = cassette_core.build(
        library, design.selections, slots(design, library.scheme), name=design.name
    )
    if result.product is None and not result.ok:
        return result

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
    notes += _multigene_notes(chosen)

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


def _multigene_notes(chosen: list[PlasmidEntry]) -> list[Issue]:
    """Warn when the cassette could never be cut out for a multigene assembly.

    The multigene enzyme releases a cassette by cutting inside the type 1 and
    type 5 parts at its ends. A part that has had that site domesticated away
    still assembles perfectly here - nothing at this level touches it - so this
    is the only moment to say so, before the plasmid is built and the failure
    turns up two steps later as a cassette with no ends.
    """
    blocked = [e for e in chosen if e.level3_ready is False]
    if not blocked:
        return []
    names = ", ".join(f"{e.name} (type {e.call.part_type})" for e in blocked)
    return [
        Issue(
            WARNING,
            "blocks_multigene",
            f"{names} carries no site for the multigene enzyme, so this cassette "
            f"assembles but can never be cut out of its plasmid for a Level 3 "
            f"assembly. Fine if this is the final construct.",
        )
    ]


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


#: Reading a construct back into parts is a core capability, not a Level 2 one -
#: Level 3 asks the same question of each transcription unit. The tiling lives
#: there; what stays here is the only part that *is* Level 2's: turning those
#: matches into panel choices and the view flags those panels need.
Match = decompose_core.Match
Unmatched = decompose_core.Unmatched


@dataclass
class Decomposition:
    """A tiling, expressed as something the Level 2 screen can load."""

    cassette: str
    length: int
    matches: list[Match] = field(default_factory=list)
    selections: dict[str, str] = field(default_factory=dict)
    split_3: bool = False
    split_4: bool = False
    split_8: bool = False
    composite_left: bool = False
    composite_right: bool = False
    unmatched: list[Unmatched] = field(default_factory=list)
    """Stretches no library part explains, with the ends they sit between."""

    @property
    def gaps(self) -> list[tuple[int, int]]:
        """The same stretches as ``(start, length)``, for callers that only count."""
        return [(u.start, u.length) for u in self.unmatched]

    @property
    def covered(self) -> int:
        return sum(m.length for m in self.matches)

    @property
    def complete(self) -> bool:
        return bool(self.matches) and not self.unmatched

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
    """Read an assembled cassette back into the parts that built it.

    The screen's other direction: instead of picking eight parts, pick a
    construct you already have and get its panels filled in, ready to swap one.
    """
    tiling = decompose_core.tile(library, cassette)
    out = Decomposition(
        cassette=cassette.name,
        length=tiling.length,
        matches=list(tiling.matches),
        unmatched=list(tiling.unmatched),
    )
    _fill_selections(out)
    return out


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


# --------------------------------------------------------------------------- #
# swapping one position across many candidates
# --------------------------------------------------------------------------- #


@dataclass
class SweepRow:
    """One candidate tried in one slot of an otherwise fixed design."""

    name: str
    display: str
    part_type: str | None
    part_length: int
    ok: bool
    length: int = 0
    issues: list[Issue] = field(default_factory=list)
    component: str = ""

    @property
    def errors(self) -> list[Issue]:
        return [i for i in self.issues if i.level == ERROR]

    @property
    def warnings(self) -> list[Issue]:
        return [i for i in self.issues if i.level == WARNING]


def sweep(
    library: Library,
    design: CassetteDesign,
    slot: str,
    candidates: Sequence[str],
) -> list[SweepRow]:
    """Build the same design once per candidate in one position.

    A promoter titration against a fixed coding sequence is the standard form
    of this: forty-seven type 2 parts, one construct each, everything else held
    still. Doing it by hand means forty-seven passes through the same screen,
    and the interesting part - which of them assemble cleanly and how long each
    comes out - is the same two facts every time.

    Every candidate is reported, including the ones that fail: a promoter that
    carries an internal site is a result, not an omission.
    """
    rows: list[SweepRow] = []
    for name in candidates:
        entry = library.get(name)
        if entry is None:
            rows.append(SweepRow(
                name=name, display=name, part_type=None, part_length=0, ok=False,
                issues=[Issue(ERROR, "unknown_part", f"no plasmid named {name}")],
            ))
            continue

        trial = replace(design, selections={**design.selections, slot: entry.name})
        result = build(library, trial)
        rows.append(SweepRow(
            name=entry.name,
            display=entry.display,
            part_type=entry.call.part_type,
            part_length=entry.length,
            ok=result.ok,
            length=result.length,
            issues=list(result.issues),
            component=entry.component,
        ))
    return rows


def sweep_products(
    library: Library,
    design: CassetteDesign,
    slot: str,
    candidates: Sequence[str],
    name_template: str = "{design}_{part}",
):
    """The same sweep, yielding ``(filename, record)`` for the ones that built.

    Kept apart from `sweep` because assembling the records is the expensive
    half and the screen only needs the table until someone asks to write.
    """
    for name in candidates:
        entry = library.get(name)
        if entry is None:
            continue
        trial = replace(design, selections={**design.selections, slot: entry.name})
        trial.name = name_template.format(design=design.name, part=entry.name)
        result = build(library, trial)
        if result.product is not None:
            yield f"{trial.name}.gb", result.product


def picklist(rows: list[SweepRow], design: CassetteDesign, slot: str) -> str:
    """The sweep as a tab-separated sheet, for a notebook or a plate map."""
    fixed = ", ".join(
        f"{key}={value}" for key, value in sorted(design.selections.items()) if key != slot
    )
    lines = [
        f"# {design.name}: position {slot} swept across {len(rows)} candidates",
        f"# held fixed: {fixed}" if fixed else "# nothing else selected",
        "",
        "\t".join(["part", "component", "part_bp", "construct_bp", "assembles",
                   "errors", "warnings", "notes"]),
    ]
    for row in rows:
        lines.append("\t".join([
            row.name,
            row.component,
            str(row.part_length),
            str(row.length) if row.ok else "",
            "yes" if row.ok else "no",
            str(len(row.errors)),
            str(len(row.warnings)),
            "; ".join(i.message for i in row.issues if i.level != INFO),
        ]))
    return "\n".join(lines) + "\n"
