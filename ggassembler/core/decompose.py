"""Read an assembled construct back into the library parts that built it.

This cannot be a digest. Assembly consumed every site the part enzyme would
have cut, so the parts are gone *as sites* while still being present *as
sequence*. So each library part's released fragment is looked for in the
construct directly, and the ones that abut end-to-end are chained into a tiling
of the circle.

It lives in `core` because the question is not level-specific: a Level 2 screen
asks it of a cassette to fill its eight panels, and a Level 3 screen asks it of
each transcription unit to show what is inside it. Neither level may import the
other, and neither should own an answer they both need.

Where the chain covers the whole circle, everything is accounted for. Where it
does not, the gap comes back *as* a gap, carrying the overhangs at both ends -
those name the position the missing part occupied even when the part itself is
nowhere on the shelf, and "1,750 bp of type 3 that isn't a part I know" is a
far more useful thing to be told than "1,750 bp unexplained".
"""

from __future__ import annotations

from dataclasses import dataclass, field

from . import parttypes, seqio
from .library import Library, PlasmidEntry
from .parttypes import YTK, Scheme


@dataclass
class Match:
    """One library part found inside an assembled construct."""

    name: str
    display: str
    part_type: str
    start: int
    end: int
    length: int
    component: str = ""


@dataclass
class Unmatched:
    """A stretch of an assembled construct that no library part accounts for.

    Reported rather than dropped, and reported with its *ends*: the four bases
    at each boundary are the overhangs the missing part was ligated by, and
    they say which position it occupies even though the part is unknown.
    """

    start: int
    end: int
    length: int
    left_overhang: str = ""
    right_overhang: str = ""
    part_type: str | None = None
    """The position those two overhangs delimit, when they delimit a known one."""


@dataclass
class Tiling:
    """What a construct is made of, as far as the library can tell."""

    source: str
    length: int
    matches: list[Match] = field(default_factory=list)
    unmatched: list[Unmatched] = field(default_factory=list)

    @property
    def covered(self) -> int:
        return sum(m.length for m in self.matches)

    @property
    def complete(self) -> bool:
        return bool(self.matches) and not self.unmatched

    def of_types(self, wanted: tuple[str, ...]) -> list[Match]:
        """The matches at the given positions, in the order they appear."""
        return [m for m in self.matches if m.part_type in wanted]


def tile(library: Library, entry: PlasmidEntry) -> Tiling:
    """Match every library part against `entry` and chain the ones that abut."""
    out = Tiling(source=entry.name, length=entry.length)
    sequence = seqio.sequence(library.record(entry)).upper()
    if not sequence:
        return out

    # search the doubled sequence so a part spanning the origin is still found
    doubled = sequence + sequence
    n = len(sequence)

    starts: dict[int, list[tuple[PlasmidEntry, int]]] = {}
    for part, fragment in library.part_fragments():
        if len(fragment) > n:
            continue
        at = doubled.find(fragment)
        while at != -1 and at < n:
            starts.setdefault(at, []).append((part, (at + len(fragment)) % n))
            at = doubled.find(fragment, at + 1)

    chain = _chain(starts, n) if starts else []
    if not chain:
        out.unmatched = [Unmatched(start=0, end=0, length=n)]
        return out

    for part, start, end in chain:
        out.matches.append(
            Match(
                name=part.name,
                display=part.display,
                part_type=part.call.part_type or "",
                start=start,
                end=end,
                length=(end - start) % n or n,
                component=part.component,
            )
        )
    out.unmatched = _unmatched(chain, sequence, library.scheme)
    return out


def _chain(
    starts: dict[int, list[tuple[PlasmidEntry, int]]], n: int
) -> list[tuple[PlasmidEntry, int, int]]:
    """The longest run of parts that abut end-to-end, closing the circle if it can.

    Tried from every match in turn, preferring atomic parts over composites, so
    the answer is as granular as the library allows: a 2-3-4 dropout and three
    separate parts describe the same stretch, and three are more use than one.
    """
    def rank(pair: tuple[PlasmidEntry, int]) -> tuple:
        part, _ = pair
        return (len(parttypes.positions(part.call.part_type or "", YTK)), part.name)

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
            part, end = candidates[0]
            used.add(part.name)
            chain.append((part, at, end))
            at = end
            if at == origin:
                break  # the circle closed
        covered = sum((e - s) % n or n for _, s, e in chain)
        if covered > sum((e - s) % n or n for _, s, e in best):
            best = chain
        if covered == n:
            break
    return best


def _unmatched(
    chain: list[tuple[PlasmidEntry, int, int]], sequence: str, scheme: Scheme
) -> list[Unmatched]:
    """The stretches the chain leaves unexplained, with the ends they sit between."""
    n = len(sequence)
    width = scheme.part_enzyme.overhang_length

    def overhang_at(position: int) -> str:
        return "".join(sequence[(position + i) % n] for i in range(width))

    if not chain:
        return [Unmatched(start=0, end=0, length=n)]
    if sum((e - s) % n or n for _, s, e in chain) >= n:
        return []

    # the chain is contiguous, so there is exactly one gap: the end of the last
    # match round to the start of the first
    start, end = chain[-1][2], chain[0][1]
    five, three = overhang_at(start), overhang_at(end)
    call = parttypes.span(five, three, scheme=scheme)
    return [
        Unmatched(
            start=start,
            end=end,
            length=(end - start) % n or n,
            left_overhang=five,
            right_overhang=three,
            part_type=call.name if call and call.legal else None,
        )
    ]
