"""Resolving eight chosen plasmids into a cassette reaction.

Lifted out of `levels/level2_cassette.py` when a second screen needed it. The
Plate screen builds ninety-six of these, and it may not import Level 2 - levels
never import one another - so the choice was to duplicate the rules or to move
them here. They are rules about what the *kit* permits, not about what either
screen looks like, so here is where they belong.

What lives here is the part that would be wrong to have two copies of: the slot
set for a given set of view modes, and the loop that turns a slot-to-plasmid
mapping into pieces, refusing the four ways it can be wrong. What stays in
Level 2 is everything about the screen - the advice it prints, the decomposition
path, the panel counts.

A drifting second copy of "a dropout is consumed by the reaction and never ends
up in the product" would not look like a bug. It would look like the Plate
screen being more permissive than the Cassette screen, which is worse.
"""

from __future__ import annotations

from dataclasses import dataclass

from . import parttypes
from .assembly import ERROR, AssemblyResult, Issue, Piece, assemble
from .library import Library, PlasmidEntry
from .parttypes import YTK, Scheme


@dataclass(frozen=True)
class Slot:
    """One position of the part circle, and the overhangs it must match."""

    key: str
    """The part type this slot takes, e.g. ``3``, ``3a``, ``234``."""
    five_prime: str
    three_prime: str
    column: str
    """``left`` or ``right`` - which side of the Level 2 screen it sits on."""

    @property
    def label(self) -> str:
        return self.key.replace("234", "2·4").replace("678", "6·8")

    @property
    def description(self) -> str:
        return parttypes.DESCRIPTIONS.get(self.key, "composite part")

    @property
    def overhangs(self) -> tuple[str, str]:
        return self.five_prime, self.three_prime


def slot_keys(
    split_3: bool = False,
    split_4: bool = False,
    split_8: bool = False,
    composite_left: bool = False,
    composite_right: bool = False,
) -> list[str]:
    """The part types a cassette needs, in circle order.

    The view modes - splitting 3, 4 or 8, collapsing 2-3-4 or 6-7-8 - are
    modes over one slot model rather than separate assembly paths. Everything
    downstream treats whatever comes back here identically.
    """
    keys: list[str] = ["1"]
    if composite_left:
        keys.append("234")
    else:
        keys.append("2")
        keys.extend(("3a", "3b") if split_3 else ("3",))
        keys.extend(("4a", "4b") if split_4 else ("4",))
    keys.append("5")
    if composite_right:
        keys.append("678")
    else:
        keys.append("6")
        keys.append("7")
        keys.extend(("8a", "8b") if split_8 else ("8",))
    return keys


def slots_for(keys: list[str], scheme: Scheme = YTK) -> list[Slot]:
    """Those keys as slots, dropping any the scheme has no overhangs for."""
    out = []
    for key in keys:
        pair = parttypes.type_overhangs(key, scheme)
        if pair is None:
            continue
        column = "left" if key in ("1", "2", "234") or key[0] in "34" else "right"
        out.append(Slot(key=key, five_prime=pair[0], three_prime=pair[1], column=column))
    return out


@dataclass
class Resolved:
    """What a set of choices amounts to, before the reaction is simulated."""

    pieces: list[Piece]
    chosen: list[PlasmidEntry]
    issues: list[Issue]

    @property
    def ok(self) -> bool:
        return not any(i.level == ERROR for i in self.issues)


def resolve(library: Library, selections: dict[str, str], slots: list[Slot]) -> Resolved:
    """Turn slot-to-plasmid choices into pieces, saying what is wrong.

    Four ways a choice can be wrong, and each is worth its own message:

    - the slot is empty
    - the name is not in the library at all
    - the plasmid is real but releases the wrong overhangs for this position
    - the plasmid is a dropout, so the reaction consumes it and it never ends
      up in the product - which looks like a part on every screen until the
      construct comes back the wrong length
    """
    pieces: list[Piece] = []
    chosen: list[PlasmidEntry] = []
    issues: list[Issue] = []

    for slot in slots:
        name = selections.get(slot.key)
        if not name:
            issues.append(Issue(
                ERROR, "empty_slot",
                f"slot {slot.label} ({slot.five_prime} -> {slot.three_prime}) is empty",
            ))
            continue

        entry = library.get(name)
        if entry is None:
            issues.append(Issue(
                ERROR, "unknown_part", f"no plasmid named {name} in the library"))
            continue

        if entry.overhangs != slot.overhangs:
            issues.append(Issue(
                ERROR, "wrong_overhangs",
                f"{entry.name} releases {entry.call.five_prime} -> "
                f"{entry.call.three_prime}, but slot {slot.label} needs "
                f"{slot.five_prime} -> {slot.three_prime}",
            ))
            continue

        if entry.call.reversed_sites:
            issues.append(Issue(
                ERROR, "dropout_as_part",
                f"{entry.name} is a dropout ({entry.call.part_type}); it is consumed "
                f"by the reaction and never ends up in the product",
            ))
            continue

        try:
            pieces.append(library.piece(entry))
        except (ValueError, OSError) as exc:
            issues.append(Issue(ERROR, "unreadable_part", f"{entry.name}: {exc}"))
            continue
        chosen.append(entry)

    return Resolved(pieces=pieces, chosen=chosen, issues=issues)


def build(
    library: Library,
    selections: dict[str, str],
    slots: list[Slot],
    name: str = "cassette",
) -> tuple[AssemblyResult, list[PlasmidEntry]]:
    """Resolve and simulate, returning the result and what went into it.

    The parts are handed back because every caller needs them for something
    the assembly does not carry - Level 2 for its advice, the Plate screen for
    the volumes it will pipette.
    """
    scheme = library.scheme
    found = resolve(library, selections, slots)

    if not found.ok:
        return (
            AssemblyResult(False, None, issues=found.issues, enzyme=scheme.part_enzyme.name),
            found.chosen,
        )

    result = assemble(found.pieces, scheme.part_enzyme, name=name)
    result.issues = found.issues + result.issues
    result.ok = result.ok and not result.errors
    if not result.ok:
        result.product = None
    return result, found.chosen
