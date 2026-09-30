"""The standing source plate: what is in each well, and how much is left.

A combinatorial sweep is only tractable at the bench because the parts live in
a plate that persists between runs. Nobody re-racks eight plasmids before every
build; they keep a normalized working plate and draw from it. So this is not a
description of one run - it is a description of a physical object on a shelf,
and the app's job is to know what is in it and how much of each is left.

Two things follow from that, and both are why this file exists at all.

**Concentration enters the data model here.** The fmol maths needs it and
nothing upstream has it: a GenBank file records a sequence, not how much of it
is in your freezer. The library can carry a measured concentration per plasmid;
this lets a well override it, because the same plasmid diluted into two wells
is two concentrations.

**Volume has to be tracked across runs, not within one.** A sweep draws from a
connector well ninety-six times. Running dry in column 9 is not an error the
app can report afterwards - the plate is already half-pipetted and the reaction
that matters has no ligase in it. So a run is checked against the remaining
volume before it is allowed to start, and the plate is decremented when it is
generated.

Imports `core` only, and no other screen's module.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any, Iterable

from ..core.assembly import ERROR, INFO, WARNING, Issue
from ..core.cassette import build as build_cassette
from ..core.cassette import slot_keys, slots_for
from ..core import protocol
from ..core.library import Library

#: The file the standing plate lives in, beside the index.
SOURCE_FILE = "source_plate.json"

#: The rows and columns of a standard 96-well plate.
ROWS = "ABCDEFGH"
COLUMNS = tuple(range(1, 13))

#: Opentrons' name for the plate this is written for. Stored rather than
#: assumed, because the protocol has to load the labware the person owns.
DEFAULT_LABWARE = "opentrons_96_wellplate_200ul_pcr_full_skirt"

#: What a pipette cannot reach. Below this the well is empty as far as a run is
#: concerned, whatever the number says.
DEFAULT_DEAD_VOLUME = 5.0

WELL_NAME = re.compile(r"^([A-H])(1[0-2]|[1-9])$")


def wells() -> list[str]:
    """Every well, row-major: A1 ... A12, B1 ... H12."""
    return [f"{row}{column}" for row in ROWS for column in COLUMNS]


def wells_column_major() -> list[str]:
    """Every well down the columns: A1, B1 ... H1, A2 ...

    The order a one-factor sweep fills, and the order an 8-channel pipette
    works in - one aspiration serves a whole column, so a design that fills
    column-major is a design that can be dispensed column-wise.
    """
    return [f"{row}{column}" for column in COLUMNS for row in ROWS]


def is_well(name: str) -> bool:
    return bool(WELL_NAME.match(name.strip().upper()))


def normalise_well(name: str) -> str:
    """``a01`` and ``A1`` are the same well; the app says ``A1``."""
    text = name.strip().upper()
    match = re.match(r"^([A-H])0*(\d{1,2})$", text)
    if not match or not 1 <= int(match.group(2)) <= 12:
        raise ValueError(f"{name!r} is not a well on a 96-well plate")
    return f"{match.group(1)}{int(match.group(2))}"


@dataclass
class SourceWell:
    """One well of the standing plate."""

    plasmid: str
    conc_ng_ul: float = 0.0
    volume_ul: float = 0.0
    note: str = ""

    def to_dict(self) -> dict[str, Any]:
        out: dict[str, Any] = {
            "plasmid": self.plasmid,
            "conc_ng_ul": round(self.conc_ng_ul, 3),
            "volume_ul": round(self.volume_ul, 3),
        }
        if self.note:
            out["note"] = self.note
        return out

    @classmethod
    def from_dict(cls, raw: dict[str, Any]) -> SourceWell:
        return cls(
            plasmid=str(raw.get("plasmid", "")).strip(),
            conc_ng_ul=float(raw.get("conc_ng_ul") or 0.0),
            volume_ul=float(raw.get("volume_ul") or 0.0),
            note=str(raw.get("note", "")),
        )


@dataclass
class SourcePlate:
    """The standing plate, as it is right now."""

    id: str = "SRC-01"
    labware: str = DEFAULT_LABWARE
    dead_volume_ul: float = DEFAULT_DEAD_VOLUME
    wells: dict[str, SourceWell] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return {
            "id": self.id,
            "labware": self.labware,
            "dead_volume_ul": self.dead_volume_ul,
            "wells": {w: self.wells[w].to_dict() for w in sorted(self.wells, key=_order)},
        }

    @classmethod
    def from_dict(cls, raw: dict[str, Any]) -> SourcePlate:
        plate = cls(
            id=str(raw.get("id") or "SRC-01"),
            labware=str(raw.get("labware") or DEFAULT_LABWARE),
            dead_volume_ul=float(raw.get("dead_volume_ul") or DEFAULT_DEAD_VOLUME),
        )
        for name, well in (raw.get("wells") or {}).items():
            try:
                plate.wells[normalise_well(name)] = SourceWell.from_dict(well)
            except (ValueError, TypeError, AttributeError):
                continue  # reported by `check`, not dropped silently on a load
        return plate

    # -- lookups ----------------------------------------------------------- #

    def find(self, plasmid: str) -> list[str]:
        """Every well holding this plasmid, in pipetting order."""
        return [
            w for w in wells_column_major()
            if w in self.wells and self.wells[w].plasmid == plasmid
        ]

    def plasmids(self) -> set[str]:
        return {w.plasmid for w in self.wells.values() if w.plasmid}

    def usable(self, well: str) -> float:
        """Volume a pipette can actually reach in this well."""
        entry = self.wells.get(well)
        return max(0.0, entry.volume_ul - self.dead_volume_ul) if entry else 0.0


def _order(well: str) -> tuple[int, str]:
    match = WELL_NAME.match(well)
    return (int(match.group(2)), match.group(1)) if match else (99, well)


# --------------------------------------------------------------------------- #
# reading and writing
# --------------------------------------------------------------------------- #


def path_for(library: Library):
    return library.cache_dir / SOURCE_FILE


def load(library: Library) -> SourcePlate:
    """The standing plate, or an empty one the first time."""
    try:
        raw = json.loads(path_for(library).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return SourcePlate()
    return SourcePlate.from_dict(raw)


def save(library: Library, plate: SourcePlate) -> SourcePlate:
    target = path_for(library)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(plate.to_dict(), indent=2, sort_keys=True), encoding="utf-8")
    return plate


def set_well(
    library: Library,
    well: str,
    plasmid: str,
    conc_ng_ul: float | None = None,
    volume_ul: float | None = None,
    note: str = "",
) -> SourcePlate:
    """Put a plasmid in a well, or clear it by passing an empty name.

    A concentration that is not given falls back to the library's own measured
    value, which is the number the Library screen records. A well may override
    it, because the same plasmid diluted into two wells is two concentrations
    and only the well knows which.
    """
    plate = load(library)
    name = normalise_well(well)

    if not plasmid.strip():
        plate.wells.pop(name, None)
        return save(library, plate)

    entry = library.get(plasmid)
    if entry is None:
        raise ValueError(f"no plasmid called {plasmid} in this library")

    # what was already in this well, but only if it is the same plasmid -
    # putting a different one in is a new well, not an edit of the old one
    previous = plate.wells.get(name)
    if previous is not None and previous.plasmid != entry.name:
        previous = None

    if conc_ng_ul is None:
        # the well's own number, else the library's measured one. The Library
        # screen records a concentration per plasmid; a well may override it,
        # because the same plasmid diluted into two wells is two numbers.
        conc_ng_ul = previous.conc_ng_ul if previous else 0.0
        if not conc_ng_ul:
            conc_ng_ul = entry.conc_ng_ul or 0.0

    plate.wells[name] = SourceWell(
        plasmid=entry.name,
        conc_ng_ul=float(conc_ng_ul),
        volume_ul=float(volume_ul if volume_ul is not None
                        else (previous.volume_ul if previous else 0.0)),
        note=note or (previous.note if previous else ""),
    )
    return save(library, plate)


def set_layout(library: Library, plate_id: str = "", labware: str = "",
               dead_volume_ul: float | None = None) -> SourcePlate:
    """The plate's own properties, without touching its wells."""
    plate = load(library)
    if plate_id.strip():
        plate.id = plate_id.strip()
    if labware.strip():
        plate.labware = labware.strip()
    if dead_volume_ul is not None:
        if dead_volume_ul < 0:
            raise ValueError("dead volume cannot be negative")
        plate.dead_volume_ul = float(dead_volume_ul)
    return save(library, plate)


# --------------------------------------------------------------------------- #
# is this plate usable?
# --------------------------------------------------------------------------- #


def check(plate: SourcePlate, library: Library) -> list[Issue]:
    """Everything wrong with the plate as it stands, worst first.

    A plate is edited by hand over weeks and the library is re-indexed under
    it, so a well can name a plasmid that has since been renamed or archived.
    Found on load rather than at the moment a protocol is generated, when the
    tips are already on.
    """
    issues: list[Issue] = []

    for well in sorted(plate.wells, key=_order):
        entry = plate.wells[well]
        if not entry.plasmid:
            issues.append(Issue(WARNING, "empty_well", f"{well} has no plasmid named"))
            continue
        if library.get(entry.plasmid) is None:
            issues.append(Issue(
                ERROR, "unknown_plasmid",
                f"{well} names {entry.plasmid}, which is not in this library - "
                f"it may have been renamed or archived since the plate was made",
            ))
        if entry.conc_ng_ul <= 0:
            issues.append(Issue(
                WARNING, "no_concentration",
                f"{well} ({entry.plasmid}) has no concentration, so no fmol can be "
                f"costed from it",
            ))
        if entry.volume_ul <= plate.dead_volume_ul:
            issues.append(Issue(
                WARNING, "at_dead_volume",
                f"{well} ({entry.plasmid}) holds {entry.volume_ul:g} µL, at or below "
                f"the {plate.dead_volume_ul:g} µL a pipette can reach",
            ))

    duplicates = _duplicates(plate)
    for plasmid, where in sorted(duplicates.items()):
        issues.append(Issue(
            INFO, "duplicated",
            f"{plasmid} is in {len(where)} wells ({', '.join(where)}) - a run will "
            f"draw from the first with enough volume",
        ))
    return issues


def _duplicates(plate: SourcePlate) -> dict[str, list[str]]:
    seen: dict[str, list[str]] = {}
    for well in wells_column_major():
        entry = plate.wells.get(well)
        if entry and entry.plasmid:
            seen.setdefault(entry.plasmid, []).append(well)
    return {name: where for name, where in seen.items() if len(where) > 1}


@dataclass
class Draw:
    """What one run takes out of the plate."""

    per_well: dict[str, float] = field(default_factory=dict)
    """How much is drawn from each well, in µL."""

    missing: list[str] = field(default_factory=list)
    """Plasmids the run needs that are not on the plate at all."""

    short: list[tuple[str, str, float, float]] = field(default_factory=list)
    """``(well, plasmid, needed, usable)`` where the well cannot cover the run."""

    @property
    def ok(self) -> bool:
        return not self.missing and not self.short


def plan_draw(plate: SourcePlate, needs: dict[str, float]) -> Draw:
    """Work out which wells a run draws from, and whether they can cover it.

    `needs` is µL per plasmid across the whole run. A plasmid in several wells
    is taken from the first that can cover the whole amount rather than split
    across two: splitting is possible at the bench and is a different protocol,
    and quietly generating one nobody asked for is worse than saying the well
    is short.
    """
    draw = Draw()
    for plasmid in sorted(needs):
        wanted = needs[plasmid]
        if wanted <= 0:
            continue
        candidates = plate.find(plasmid)
        if not candidates:
            draw.missing.append(plasmid)
            continue
        chosen = next((w for w in candidates if plate.usable(w) >= wanted), None)
        if chosen is None:
            best = max(candidates, key=plate.usable)
            draw.short.append((best, plasmid, wanted, plate.usable(best)))
            continue
        draw.per_well[chosen] = draw.per_well.get(chosen, 0.0) + wanted
    return draw


def apply_draw(library: Library, draw: Draw) -> SourcePlate:
    """Take the run's volumes out of the plate and write it back.

    Only ever called once a protocol has actually been generated. The plate is
    a record of a physical object, so decrementing it for a run that was never
    pipetted is how it stops being one.
    """
    if not draw.ok:
        raise ValueError("this run cannot be drawn from the plate as it stands")
    plate = load(library)
    for well, volume in draw.per_well.items():
        entry = plate.wells.get(well)
        if entry is None:
            raise ValueError(f"{well} is no longer on the plate")
        plate.wells[well] = replace(entry, volume_ul=max(0.0, entry.volume_ul - volume))
    return save(library, plate)


def summary(plate: SourcePlate, library: Library, needs: dict[str, float] | None = None
            ) -> dict[str, Any]:
    """Everything the screen draws, in one shape."""
    draw = plan_draw(plate, needs or {})
    issues = check(plate, library)
    return {
        "id": plate.id,
        "labware": plate.labware,
        "dead_volume_ul": plate.dead_volume_ul,
        "wells": [
            {
                "well": well,
                "plasmid": entry.plasmid,
                "conc_ng_ul": entry.conc_ng_ul,
                "volume_ul": entry.volume_ul,
                "usable_ul": plate.usable(well),
                "note": entry.note,
                "known": library.get(entry.plasmid) is not None,
                "part_type": _part_type(library, entry.plasmid),
                "draw_ul": round(draw.per_well.get(well, 0.0), 2),
            }
            for well, entry in sorted(plate.wells.items(), key=lambda kv: _order(kv[0]))
        ],
        "issues": [{"level": i.level, "code": i.code, "message": i.message} for i in issues],
        "missing": draw.missing,
        "short": [
            {"well": w, "plasmid": p, "needed_ul": round(n, 2), "usable_ul": round(u, 2)}
            for w, p, n, u in draw.short
        ],
        "ok": draw.ok and not any(i.level == ERROR for i in issues),
    }


def _part_type(library: Library, plasmid: str) -> str | None:
    entry = library.get(plasmid)
    return entry.call.part_type if entry else None


def needs_from(counts: Iterable[tuple[str, float]]) -> dict[str, float]:
    """Fold ``(plasmid, µL)`` pairs into one total per plasmid."""
    out: dict[str, float] = {}
    for plasmid, volume in counts:
        out[plasmid] = out.get(plasmid, 0.0) + volume
    return out


# --------------------------------------------------------------------------- #
# the sweep
# --------------------------------------------------------------------------- #
#
# One base cassette with the backbone and connectors fixed, and one or two
# positions varying across the plate: 8 variants down the rows, 12 across the
# columns.
#
# The layout is the factorial, which is the whole point. Row A is promoter 1;
# column 3 is CDS 3. Position already carries identity, so the well fill is
# free to carry something you cannot get any other way - whether that well
# assembles. And because a factor is shared along an axis, a bad part lights up
# an entire row or column: eight identical warnings are one part problem, and
# saying it once on the axis is worth more than saying it eight times.

@dataclass
class Factor:
    """One axis of the sweep: a position, and what varies in it."""

    position: str = ""
    """A slot key - ``2``, ``3``, ``4a`` and so on. Empty means no factor."""

    candidates: list[str] = field(default_factory=list)
    """Plasmid names, in the order they are laid out along the axis."""

    @property
    def active(self) -> bool:
        return bool(self.position and self.candidates)


@dataclass
class SweepDesign:
    """A base cassette, and up to two positions varying across a plate."""

    base: dict[str, str] = field(default_factory=dict)
    """Slot key -> plasmid, for every position that does not vary."""

    row: Factor = field(default_factory=Factor)
    column: Factor = field(default_factory=Factor)

    split_3: bool = False
    split_4: bool = False
    split_8: bool = False
    composite_left: bool = False
    composite_right: bool = False

    name: str = "sweep"
    pattern: str = "{base}_{row}_{col}"
    """Tokens: ``{base}``, ``{row}``, ``{col}``, ``{well}``."""

    def keys(self) -> list[str]:
        return slot_keys(
            split_3=self.split_3, split_4=self.split_4, split_8=self.split_8,
            composite_left=self.composite_left, composite_right=self.composite_right,
        )

    @property
    def factors(self) -> list[Factor]:
        return [f for f in (self.row, self.column) if f.active]


@dataclass
class PlateWell:
    """One well of the destination plate, assembled in silico."""

    well: str
    row_index: int
    column_index: int
    name: str
    selections: dict[str, str] = field(default_factory=dict)
    row_variant: str = ""
    column_variant: str = ""
    length: int = 0
    issues: list[Issue] = field(default_factory=list)

    @property
    def errors(self) -> list[Issue]:
        return [i for i in self.issues if i.level == ERROR]

    @property
    def warnings(self) -> list[Issue]:
        return [i for i in self.issues if i.level == WARNING]

    @property
    def status(self) -> str:
        if self.errors:
            return "err"
        return "warn" if self.warnings else "ok"


def _fits(library: Library, position: str, plasmid: str, keys: list[str],
          scheme=None) -> str:
    """Why this candidate cannot go in this position, or an empty string.

    Checked at selection time rather than at assembly. Ninety-six identical
    "wrong overhangs" errors is the same fact said ninety-six times, and it is
    a fact about the candidate you picked, not about the plate.
    """
    slot = next((s for s in slots_for(keys, library.scheme) if s.key == position), None)
    if slot is None:
        return f"{position} is not a position in this cassette"
    entry = library.get(plasmid)
    if entry is None:
        return f"{plasmid} is not in the library"
    if entry.call.reversed_sites:
        return f"{plasmid} is a dropout; the reaction consumes it"
    if entry.overhangs != slot.overhangs:
        return (f"{plasmid} releases {entry.call.five_prime} → {entry.call.three_prime}, "
                f"but position {slot.label} needs {slot.five_prime} → {slot.three_prime}")
    return ""


def check_candidates(library: Library, design: SweepDesign) -> list[Issue]:
    """Every candidate that cannot sit in the position it was chosen for."""
    keys = design.keys()
    issues: list[Issue] = []
    for axis, factor in (("row", design.row), ("column", design.column)):
        if not factor.position:
            continue
        limit = 8 if axis == "row" else 12
        if len(factor.candidates) > limit:
            issues.append(Issue(
                ERROR, "too_many",
                f"the {axis} axis holds {limit} variants; {len(factor.candidates)} were chosen",
            ))
        for plasmid in factor.candidates:
            why = _fits(library, factor.position, plasmid, keys)
            if why:
                issues.append(Issue(ERROR, "bad_candidate", f"{axis} axis: {why}"))
    if not design.factors:
        issues.append(Issue(ERROR, "no_factor", "nothing varies, so there is no sweep"))
    return issues


def well_name(design: SweepDesign, base: str, row: str, column: str, well: str) -> str:
    """The construct's name, from the pattern."""
    out = design.pattern or "{base}_{row}_{col}"
    for token, value in (("{base}", base), ("{row}", row), ("{col}", column),
                         ("{well}", well)):
        out = out.replace(token, value)
    return re.sub(r"_+", "_", out).strip("_") or well


def layout(design: SweepDesign) -> list[tuple[str, int, int, str, str]]:
    """Which variants land in which well: ``(well, row, column, row var, col var)``.

    Two factors make the obvious grid. One factor fills column-major, because
    that is the order an 8-channel works in - a design that fills down the
    columns is a design that can be dispensed a column at a time.
    """
    rows, columns = design.row, design.column
    out: list[tuple[str, int, int, str, str]] = []

    if rows.active and columns.active:
        for c, column_variant in enumerate(columns.candidates):
            for r, row_variant in enumerate(rows.candidates):
                out.append((f"{ROWS[r]}{c + 1}", r, c, row_variant, column_variant))
        return out

    only = rows if rows.active else columns
    for n, variant in enumerate(only.candidates[:96]):
        well = wells_column_major()[n]
        r, c = ROWS.index(well[0]), int(well[1:]) - 1
        out.append((well, r, c, variant if rows.active else "",
                    "" if rows.active else variant))
    return out


def build_plate(library: Library, design: SweepDesign) -> list[PlateWell]:
    """Assemble every well in silico.

    This is the thing the app can do that nobody does by hand: ninety-six
    assemblies, each checked for junction compatibility, repeated overhangs,
    internal sites and near-matched overhangs, before a single tip goes on.
    """
    keys = design.keys()
    slots = slots_for(keys, library.scheme)
    out: list[PlateWell] = []

    for well, r, c, row_variant, column_variant in layout(design):
        selections = dict(design.base)
        if row_variant:
            selections[design.row.position] = row_variant
        if column_variant:
            selections[design.column.position] = column_variant

        name = well_name(
            design, design.name,
            _short(library, row_variant), _short(library, column_variant), well,
        )
        result, _chosen = build_cassette(library, selections, slots, name=name)
        out.append(PlateWell(
            well=well, row_index=r, column_index=c, name=name,
            selections=selections,
            row_variant=row_variant, column_variant=column_variant,
            length=result.length if result.ok else 0,
            issues=list(result.issues),
        ))
    return out


def _short(library: Library, plasmid: str) -> str:
    """A name fit for a construct name: the component where there is one."""
    if not plasmid:
        return ""
    entry = library.get(plasmid)
    label = (entry.component if entry and entry.component else plasmid)
    return re.sub(r"[^A-Za-z0-9]+", "", label.split("·")[0])[:16] or plasmid


def axis_summary(built: list[PlateWell], design: SweepDesign) -> dict[str, Any]:
    """What is wrong with each row and column, rather than with each well.

    A factor is shared along an axis, so a part with an internal BsaI site
    fails every well in its column. Reporting that eight times is reporting the
    wrong thing: it is one part problem, and the axis is where it belongs.
    """
    rows: dict[int, list[PlateWell]] = {}
    columns: dict[int, list[PlateWell]] = {}
    for well in built:
        rows.setdefault(well.row_index, []).append(well)
        columns.setdefault(well.column_index, []).append(well)

    def describe(group: list[PlateWell], variant: str) -> dict[str, Any]:
        codes: dict[str, int] = {}
        for well in group:
            for issue in well.issues:
                if issue.level in (ERROR, WARNING):
                    codes[issue.code] = codes.get(issue.code, 0) + 1
        whole = [code for code, n in codes.items() if n == len(group) and len(group) > 1]
        return {
            "variant": variant,
            "wells": len(group),
            "err": sum(1 for w in group if w.status == "err"),
            "warn": sum(1 for w in group if w.status == "warn"),
            # a code that fired in *every* well on this axis is the part, not
            # the wells - that is the sentence worth printing
            "shared": sorted(whole),
        }

    return {
        "rows": [
            describe(rows[i], design.row.candidates[i] if i < len(design.row.candidates) else "")
            for i in sorted(rows)
        ],
        "columns": [
            describe(columns[i],
                     design.column.candidates[i] if i < len(design.column.candidates) else "")
            for i in sorted(columns)
        ],
    }


def part_usage(built: list[PlateWell]) -> dict[str, int]:
    """How many wells each plasmid goes into - the basis of every volume."""
    counts: dict[str, int] = {}
    for well in built:
        for plasmid in well.selections.values():
            counts[plasmid] = counts.get(plasmid, 0) + 1
    return counts


def shared_positions(design: SweepDesign) -> list[str]:
    """The positions identical in every well, which go into the master mix.

    The insight that makes this tractable at the bench. Ninety-six reactions
    times eight parts is 768 transfers, which is not a protocol anyone runs.
    But six of the eight positions are the same in every well, so they are
    mixed once and distributed with the 8-channel, and only the varying
    positions are cherry-picked.
    """
    varying = {f.position for f in design.factors}
    return [key for key in design.keys() if key not in varying]


# --------------------------------------------------------------------------- #
# the reaction
# --------------------------------------------------------------------------- #
#
# Ninety-six reactions times eight parts is 768 transfers, which nobody runs.
# The sweep's own structure is what makes it tractable: six of the eight
# positions are identical in every well, so they go into one master mix with
# the buffer, ligase and enzyme, are distributed once with the 8-channel, and
# only the varying positions are cherry-picked.
#
# So the maths here is in two halves - what goes into the mix once, and what is
# transferred per well - and those halves are the protocol's shape as well as
# its arithmetic. The fmol is `core.protocol`'s; nothing here recomputes it.


@dataclass
class ReactionSetup:
    """The numbers the bench decides, with defaults rather than constants."""

    total_ul: float = 10.0
    """Final volume of each well's reaction."""

    part_ul: float = 1.0
    """Volume transferred per cherry-picked part."""

    target_fmol: float = protocol.DEFAULT_FMOL
    """Equimolar target per part, used to check the transfer actually delivers it."""

    overage: float = 1.15
    """Extra master mix, for what the tips and the reservoir keep."""

    buffer_ul: float = 1.0
    ligase_ul: float = 0.5
    enzyme_ul: float = 0.5
    """T4 buffer, T4 ligase and the part enzyme, per reaction."""

    # How much mix each well gets is not a setting: it is everything that is
    # not cherry-picked, which depends on how many positions vary. `plan`
    # works it out and puts it on the plan.


@dataclass
class MixComponent:
    """One line of the master-mix recipe."""

    name: str
    kind: str
    per_well_ul: float
    total_ul: float
    conc_ng_ul: float = 0.0
    length: int = 0
    fmol: float = 0.0
    measured: bool = True
    note: str = ""


@dataclass
class Transfer:
    """One cherry-picked move: source well to destination well."""

    source: str
    destination: str
    plasmid: str
    volume_ul: float
    axis: str
    """``row`` or ``column`` - which factor this transfer serves."""


@dataclass
class ReactionPlan:
    """Everything needed to cost, check and then pipette a plate."""

    setup: ReactionSetup
    mix: list[MixComponent] = field(default_factory=list)
    transfers: list[Transfer] = field(default_factory=list)
    mix_per_well_ul: float = 0.0
    mix_total_ul: float = 0.0
    wells: int = 0
    draw: Draw = field(default_factory=Draw)
    issues: list[Issue] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return self.draw.ok and not any(i.level == ERROR for i in self.issues)


def plan_reaction(
    library: Library,
    design: SweepDesign,
    built: list[PlateWell],
    plate: SourcePlate,
    setup: ReactionSetup | None = None,
) -> ReactionPlan:
    """Cost the whole plate: what goes in the mix, and what is moved per well."""
    setup = setup or ReactionSetup()
    varying = [f.position for f in design.factors]
    shared = shared_positions(design)
    wells = len(built)
    plan = ReactionPlan(setup=setup, wells=wells)
    if not wells:
        plan.issues.append(Issue(ERROR, "no_wells", "this sweep has no wells"))
        return plan

    # -- the mix ---------------------------------------------------------- #
    for key in shared:
        plasmid = design.base.get(key)
        entry = library.get(plasmid) if plasmid else None
        if entry is None:
            plan.issues.append(Issue(
                ERROR, "mix_part_missing",
                f"position {key} is in every well but names no plasmid in the library",
            ))
            continue
        conc = _concentration(plate, entry)
        measured = conc > 0
        if not measured:
            conc = protocol.ASSUMED_CONC
        ng = protocol.ng_for(setup.target_fmol, entry.length)
        per_well = protocol.volume_for(ng, conc)
        plan.mix.append(MixComponent(
            name=entry.name, kind=f"part {key}",
            per_well_ul=round(per_well, 3),
            total_ul=round(per_well * wells * setup.overage, 2),
            conc_ng_ul=conc, length=entry.length,
            fmol=setup.target_fmol, measured=measured,
            note="" if measured else f"assumed {protocol.ASSUMED_CONC:g} ng/µL",
        ))
        if not measured:
            plan.issues.append(Issue(
                WARNING, "assumed_concentration",
                f"{entry.name} has no measured concentration; costed at "
                f"{protocol.ASSUMED_CONC:g} ng/µL",
            ))

    for name, kind, volume in (
        ("T4 DNA ligase buffer", "reagent", setup.buffer_ul),
        ("T4 DNA ligase", "reagent", setup.ligase_ul),
        (library.scheme.part_enzyme.name, "enzyme", setup.enzyme_ul),
    ):
        plan.mix.append(MixComponent(
            name=name, kind=kind, per_well_ul=volume,
            total_ul=round(volume * wells * setup.overage, 2),
        ))

    dna_and_reagents = sum(c.per_well_ul for c in plan.mix)
    picked = len(varying) * setup.part_ul
    water = setup.total_ul - dna_and_reagents - picked
    if water < -0.005:
        plan.issues.append(Issue(
            ERROR, "over_volume",
            f"the parts and reagents come to {dna_and_reagents + picked:.2f} µL, more than "
            f"the {setup.total_ul:g} µL reaction - raise the volume or lower the fmol target",
        ))
    else:
        plan.mix.append(MixComponent(
            name="Nuclease-free water", kind="reagent",
            per_well_ul=round(max(0.0, water), 3),
            total_ul=round(max(0.0, water) * wells * setup.overage, 2),
        ))

    plan.mix_per_well_ul = round(sum(c.per_well_ul for c in plan.mix), 3)
    plan.mix_total_ul = round(plan.mix_per_well_ul * wells * setup.overage, 2)

    # -- the cherry-picked parts ------------------------------------------ #
    #
    # The volume is counted for every well that needs the part, whether or not
    # the plate has it. A part with nowhere to come from produces no transfer,
    # and if it also produced no *demand* the run would look costed and would
    # be short one component in every well - which is the exact failure this
    # whole screen exists to catch, arrived at by the tool that was meant to
    # catch it.
    picked_needs: list[tuple[str, float]] = []
    for well in built:
        for axis, position in (("row", design.row.position), ("column", design.column.position)):
            if position not in varying:
                continue
            plasmid = well.selections.get(position)
            if not plasmid:
                continue
            picked_needs.append((plasmid, setup.part_ul))
            source = _source_well(plate, plasmid)
            if source is not None:
                plan.transfers.append(Transfer(
                    source=source, destination=well.well, plasmid=plasmid,
                    volume_ul=setup.part_ul, axis=axis,
                ))

    # -- what the run takes out of the plate ------------------------------- #
    needs = needs_from(
        picked_needs
        + [(c.name, c.total_ul) for c in plan.mix if c.kind.startswith("part")]
    )
    plan.draw = plan_draw(plate, needs)
    for plasmid in plan.draw.missing:
        plan.issues.append(Issue(
            ERROR, "not_on_plate",
            f"{plasmid} is not on the source plate - add it to a well before generating",
        ))
    for well, plasmid, needed, usable in plan.draw.short:
        plan.issues.append(Issue(
            ERROR, "well_short",
            f"{well} holds {usable:.1f} µL of {plasmid} above the dead volume, and this "
            f"run needs {needed:.1f} µL",
        ))

    # the fmol a cherry-picked transfer actually delivers, which is the number
    # the target is only a request for
    for position in varying:
        for plasmid in sorted({w.selections.get(position, "") for w in built} - {""}):
            entry = library.get(plasmid)
            if entry is None:
                continue
            conc = _concentration(plate, entry)
            if conc <= 0:
                continue
            delivered = protocol.fmol_for(conc * setup.part_ul, entry.length)
            if delivered < setup.target_fmol * 0.5 or delivered > setup.target_fmol * 2.0:
                plan.issues.append(Issue(
                    WARNING, "off_target_fmol",
                    f"{setup.part_ul:g} µL of {plasmid} delivers {delivered:.0f} fmol, "
                    f"against a {setup.target_fmol:g} fmol target - dilute it, or change "
                    f"the per-part volume",
                ))
    return plan


def _concentration(plate: SourcePlate, entry) -> float:
    """The well's number if it has one, else the library's measured value."""
    for well in plate.find(entry.name):
        if plate.wells[well].conc_ng_ul > 0:
            return plate.wells[well].conc_ng_ul
    return entry.conc_ng_ul or 0.0


def _source_well(plate: SourcePlate, plasmid: str) -> str | None:
    found = plate.find(plasmid)
    return found[0] if found else None


# --------------------------------------------------------------------------- #
# outputs
# --------------------------------------------------------------------------- #


def plate_map_csv(built: list[PlateWell], design: SweepDesign, plate_id: str = "") -> str:
    """Well, construct, every position, and the predicted length.

    The file that gets printed and taped to the bench, so it is written for a
    person reading a row rather than for a parser: one column per position, in
    circle order, named by the position rather than by what happens to be in it.
    """
    import csv
    import io

    keys = design.keys()
    out = io.StringIO()
    writer = csv.writer(out, lineterminator="\n")
    writer.writerow(
        ["plate", "well", "construct", *[f"part_{k}" for k in keys], "length_bp", "status"]
    )
    for well in built:
        writer.writerow([
            plate_id, well.well, well.name,
            *[well.selections.get(k, "") for k in keys],
            well.length or "", well.status,
        ])
    return out.getvalue()


def products(library: Library, design: SweepDesign, built: list[PlateWell]):
    """The predicted plasmid for every well that assembles, as records.

    Rebuilt rather than kept from validation: holding ninety-six SeqRecords in
    memory through a screen's lifetime to save a second of CPU is the wrong
    trade, and this way the file that is written is produced by the same path
    that said the well was clean.
    """
    slots = slots_for(design.keys(), library.scheme)
    out = []
    for well in built:
        if well.status == "err":
            continue
        result, _ = build_cassette(library, well.selections, slots, name=well.name)
        if result.product is not None:
            record = result.product
            record.id = well.name[:16] or well.well
            record.name = record.id
            record.description = f"{design.name} {well.well}"
            out.append((well, record))
    return out


def log_builds(library: Library, design: SweepDesign, built: list[PlateWell],
               plate_id: str) -> int:
    """One row per well, carrying the plate and well it came from.

    The point of the plate id and the well: a colony picked off H7 three weeks
    later is traceable back to what went into H7, which is the question the
    build log exists to answer and the one a per-plate record cannot.
    """
    written = 0
    for well in built:
        if well.status == "err":
            continue
        parts = [library.get(name) for name in well.selections.values()]
        library.log_build(
            name=f"{well.name} [{plate_id} {well.well}]",
            level="plate",
            parts=[e for e in parts if e is not None],
            length=well.length,
            issues=sorted({i.code for i in well.issues if i.level != INFO}),
        )
        written += 1
    return written


# --------------------------------------------------------------------------- #
# the protocol
# --------------------------------------------------------------------------- #
#
# A self-contained protocol.py per run means every plate ships code that has
# never been executed. So the file is deliberately two halves:
#
#   - a literal data block (PLATE, MIX, TRANSFERS) that changes every run and
#     is therefore diffable against the last one;
#   - a runner below it that is byte-identical every time, generated from
#     RUNNER, so the code path is the one that was tested.
#
# And it is simulated before anyone is offered the download. Deck plan and API
# names below are from the current Opentrons Flex documentation, not memory:
# the Thermocycler takes no slot argument because it occupies fixed slots, and
# `apiLevel` belongs in `requirements` rather than in `metadata`.

#: What the generated protocol asks the robot for. Every one of these is a load
#: name the caller may need to change for their own deck, so they are data in
#: the generated file rather than buried in the runner.
DECK = {
    "api_level": "2.20",
    "robot": "Flex",
    # the Thermocycler occupies fixed slots and is loaded without one
    "thermocycler": "thermocyclerModuleV2",
    "destination_labware": "opentrons_96_wellplate_200ul_pcr_full_skirt",
    "source_slot": "C1",
    "reservoir_slot": "C2",
    "reservoir_labware": "opentrons_24_tuberack_nest_1.5ml_snapcap",
    "reservoir_well": "A1",
    "tiprack_slots": ["B2", "B3"],
    "tiprack_labware": "opentrons_flex_96_tiprack_50ul",
    "trash_slot": "A3",
    "single_pipette": "flex_1channel_50",
    "multi_pipette": "flex_8channel_50",
    "single_mount": "left",
    "multi_mount": "right",
}

#: The thermocycling the kit calls for: 30 cycles of digest and ligate, then a
#: final digest, then heat inactivation.
def _cycling(enzyme: str) -> dict[str, Any]:
    return {
        "cycles": protocol.CYCLES,
        "digest_c": protocol.CUT_TEMPERATURE.get(enzyme, protocol.DEFAULT_CUT_TEMPERATURE),
        "ligate_c": 16,
        "digest_seconds": 300,
        "ligate_seconds": 300,
        "final_digest_c": 60,
        "final_digest_minutes": 10,
        "inactivate_c": 80,
        "inactivate_minutes": 10,
        "hold_c": 10,
        "lid_c": 105,
    }


RUNNER = '''
# --------------------------------------------------------------------------- #
# Everything below this line is the same in every protocol this app generates.
# It is not edited per run: the run is the data above. If you need to change
# how the plate is built, change it in the app so the next hundred plates get
# the same fix.
# --------------------------------------------------------------------------- #

def _groups(transfers):
    """One tip per plasmid, not per transfer.

    A tip may revisit wells holding the same DNA; it may never carry one part
    into another part's well. Cross-contamination here does not look like
    contamination three days later - it looks like a cloning failure.
    """
    out = {}
    for move in transfers:
        out.setdefault((move["source"], move["plasmid"]), []).append(move)
    return out


def run(protocol):
    tips = [
        protocol.load_labware(DECK["tiprack_labware"], slot)
        for slot in DECK["tiprack_slots"]
    ]
    protocol.load_trash_bin(DECK["trash_slot"])

    single = protocol.load_instrument(
        DECK["single_pipette"], DECK["single_mount"], tip_racks=tips)
    multi = protocol.load_instrument(
        DECK["multi_pipette"], DECK["multi_mount"], tip_racks=tips)

    thermocycler = protocol.load_module(DECK["thermocycler"])
    destination = thermocycler.load_labware(DECK["destination_labware"])
    source = protocol.load_labware(PLATE["source_labware"], DECK["source_slot"])
    reservoir = protocol.load_labware(DECK["reservoir_labware"], DECK["reservoir_slot"])
    mix = reservoir[DECK["reservoir_well"]]

    thermocycler.open_lid()

    # 1. master mix, one column at a time with the 8-channel. Six of the eight
    #    positions are identical in every well, so they are already in this
    #    tube: 96 reactions cost one dispense per column rather than 768.
    multi.pick_up_tip()
    for column in PLATE["columns"]:
        multi.transfer(
            PLATE["mix_per_well_ul"],
            mix,
            destination.columns_by_name()[str(column)][0],
            new_tip="never",
            mix_after=(1, PLATE["mix_per_well_ul"] / 2),
        )
    multi.drop_tip()

    # 2. the parts that vary, cherry-picked. One tip per plasmid, multi-
    #    dispensed to every well that gets it.
    for (source_well, plasmid), moves in _groups(TRANSFERS).items():
        single.pick_up_tip()
        for move in moves:
            single.transfer(
                move["volume_ul"],
                source[source_well],
                destination[move["destination"]],
                new_tip="never",
            )
        single.drop_tip()

    # 3. mix each well once the parts are in
    for well in PLATE["wells"]:
        single.pick_up_tip()
        single.mix(PLATE["mix_repeats"], PLATE["mix_volume_ul"], destination[well])
        single.drop_tip()

    # 4. seal by hand. The thermocycler lid is not a seal and a 10 uL reaction
    #    at 60 C for ten minutes will not survive being unsealed.
    protocol.pause(
        "Seal the plate, return it to the thermocycler, and resume."
    )

    # 5. cycle
    thermocycler.close_lid()
    thermocycler.set_lid_temperature(CYCLING["lid_c"])
    thermocycler.execute_profile(
        steps=[
            {"temperature": CYCLING["digest_c"], "hold_time_seconds": CYCLING["digest_seconds"]},
            {"temperature": CYCLING["ligate_c"], "hold_time_seconds": CYCLING["ligate_seconds"]},
        ],
        repetitions=CYCLING["cycles"],
        block_max_volume=PLATE["reaction_ul"],
    )
    thermocycler.set_block_temperature(
        CYCLING["final_digest_c"], hold_time_minutes=CYCLING["final_digest_minutes"],
        block_max_volume=PLATE["reaction_ul"])
    thermocycler.set_block_temperature(
        CYCLING["inactivate_c"], hold_time_minutes=CYCLING["inactivate_minutes"],
        block_max_volume=PLATE["reaction_ul"])
    thermocycler.set_block_temperature(CYCLING["hold_c"])
    thermocycler.deactivate_lid()
    thermocycler.open_lid()
'''


def source_hash(plate: SourcePlate) -> str:
    """A fingerprint of the plate this run was planned against.

    In the header so that a protocol found on a USB stick in three weeks can be
    checked against the plate as it is now. A plate that has been re-racked
    since makes every source well in the file a guess.
    """
    import hashlib

    return hashlib.sha1(
        json.dumps(plate.to_dict(), sort_keys=True).encode("utf-8")
    ).hexdigest()[:12]


def generate_protocol(
    library: Library,
    design: SweepDesign,
    built: list[PlateWell],
    plate: SourcePlate,
    plan: ReactionPlan,
    version: str = "",
) -> str:
    """The run as data, then the runner that is the same every time."""
    import datetime

    enzyme = library.scheme.part_enzyme.name
    used = sorted({w.well[1:] for w in built}, key=int)
    header = [
        f"# Generated by GG Assembler {version or 'dev'}",
        f"# {datetime.datetime.now().astimezone().isoformat(timespec='seconds')}",
        "#",
        f"# Source plate   {plate.id} ({plate.labware})",
        f"# Plate map hash {source_hash(plate)}",
        f"# Base cassette  " + ", ".join(
            f"{k}={v}" for k, v in sorted(design.base.items()) if k not in
            {f.position for f in design.factors}
        ),
        f"# Row factor     " + (
            f"{design.row.position}: {', '.join(design.row.candidates)}"
            if design.row.active else "none"),
        f"# Column factor  " + (
            f"{design.column.position}: {', '.join(design.column.candidates)}"
            if design.column.active else "none"),
        f"# Wells          {len(built)}",
        "#",
        "# Simulate before running:  opentrons_simulate this_file.py",
    ]

    data = {
        "PLATE": {
            "name": design.name,
            "description": (
                f"{len(built)}-well Golden Gate sweep, {enzyme}, "
                f"from source plate {plate.id}"
            ),
            "source_labware": plate.labware,
            "columns": [int(c) for c in used],
            "wells": [w.well for w in built],
            "mix_per_well_ul": plan.mix_per_well_ul,
            "mix_total_ul": plan.mix_total_ul,
            "reaction_ul": plan.setup.total_ul,
            "mix_volume_ul": round(plan.setup.total_ul / 2, 2),
            "mix_repeats": 3,
        },
        "DECK": DECK,
        "CYCLING": _cycling(enzyme),
        "MIX": [
            {"name": c.name, "kind": c.kind, "per_well_ul": c.per_well_ul,
             "total_ul": c.total_ul, "note": c.note}
            for c in plan.mix
        ],
        "TRANSFERS": [
            {"source": t.source, "destination": t.destination, "plasmid": t.plasmid,
             "volume_ul": t.volume_ul, "axis": t.axis}
            for t in plan.transfers
        ],
    }

    # `metadata` and `requirements` are read *statically* by the Opentrons
    # parser, before the file is ever executed: it walks the AST and refuses
    # anything that is not a literal - no subscripts, no calls. So they are
    # written out as literals here rather than assembled in the runner from
    # PLATE and DECK, which is what the first version did and which
    # `opentrons_simulate` rejected with "Could not read the contents of the
    # metadata dict". No documentation page says this; the simulator does.
    body = ["\n".join(header), ""]
    body.append("metadata = " + json.dumps({
        "protocolName": design.name,
        "author": "GG Assembler",
        "description": data["PLATE"]["description"],
    }, indent=4))
    body.append("")
    body.append("requirements = " + json.dumps({
        "robotType": DECK["robot"], "apiLevel": DECK["api_level"],
    }, indent=4))
    body.append("")
    for key in ("DECK", "PLATE", "CYCLING", "MIX", "TRANSFERS"):
        body.append(f"{key} = " + json.dumps(data[key], indent=4))
        body.append("")
    body.append(RUNNER.strip())
    return "\n".join(body) + "\n"


@dataclass
class Simulation:
    """What `opentrons_simulate` made of a generated protocol."""

    ran: bool
    """False when the simulator is not installed - which is not a pass."""

    ok: bool = False
    output: str = ""
    detail: str = ""

    @property
    def verdict(self) -> str:
        if not self.ran:
            return "not simulated"
        return "passed" if self.ok else "failed"


def simulate(source: str, timeout: float = 180.0) -> Simulation:
    """Run the generated file through `opentrons_simulate`.

    The download is not offered until this passes. A protocol that has never
    been executed is exactly the thing that fails at the bench with the tips
    already on, and the simulator is the only cheap way to find out.

    A missing simulator is reported as a missing simulator. Treating "could not
    check" as "checked" is how a guard becomes worse than no guard.
    """
    import shutil
    import subprocess
    import tempfile

    binary = shutil.which("opentrons_simulate")
    if binary is None:
        return Simulation(
            ran=False,
            detail="opentrons_simulate is not installed here, so this protocol has "
                   "not been checked. Install the opentrons package, or simulate it "
                   "yourself before running the plate.",
        )

    with tempfile.TemporaryDirectory() as folder:
        path = Path(folder) / "protocol.py"
        path.write_text(source, encoding="utf-8")
        try:
            done = subprocess.run(
                [binary, str(path)], capture_output=True, text=True, timeout=timeout
            )
        except subprocess.TimeoutExpired:
            return Simulation(ran=True, ok=False,
                              detail=f"the simulator did not finish within {timeout:g}s")

    output = (done.stdout or "") + (done.stderr or "")
    return Simulation(
        ran=True,
        ok=done.returncode == 0,
        output=output[-8000:],
        detail="" if done.returncode == 0 else _first_error(output),
    )


def _first_error(output: str) -> str:
    """The line worth putting on screen, out of a page of simulator output."""
    for line in reversed(output.splitlines()):
        if line.strip() and not line.startswith((" ", "\t")):
            return line.strip()[:300]
    return "the simulator rejected this protocol"
