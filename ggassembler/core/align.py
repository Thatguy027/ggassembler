"""Compare a sequencing result against the plasmid it was meant to be.

A vendor returns a consensus for each clone: the whole circle, starting at
whatever base the assembler happened to open it at, and as often as not on the
opposite strand. None of that is a difference, and a tool that reports it as
one is worse than useless - so the first job here is to work out *how the two
molecules correspond* before a single base is compared.

That happens in three steps.

**Place it.** Shared k-mers vote on a diagonal: for each one, `ref_pos -
clone_pos` modulo the reference length. The same molecule read from a different
start gives one overwhelming winner, on one strand. If neither strand produces
a clear diagonal, the clone is not this plasmid and is reported as such rather
than aligned into nonsense.

**Cut the circle inside an exact match.** Both sequences are then rotated to
begin at the same long exact match. This matters more than it looks: rotating
by the diagonal alone puts the seam at an arbitrary base, and if the estimate
is off by three the alignment opens with a three-base indel that is not real.
Starting both strings inside a shared k-mer cannot do that.

**Align between anchors.** Unique shared k-mers are chained into collinear
blocks, and only the stretches *between* them are handed to a real pairwise
alignment. A correct clone is nearly all anchor, so this is fast; a clone with
a dropout has one large window, which is exactly where the expensive algorithm
is worth paying for.

What comes out is per reference base: which base the clone had there, what it
inserted before it, and nothing else. Differences are read off that, grouped
into runs, and named by the feature they land in - because "a substitution at
2,341" is a coordinate, and "a substitution in the XI coding sequence" is the
thing you actually have to decide about.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, Sequence

from Bio import Align
from Bio.SeqRecord import SeqRecord

#: Anchor length. Long enough to be unique in a plasmid, short enough that a
#: consensus with the odd error still produces plenty of them.
ANCHOR_K = 20

#: Retried with this when the first pass finds too few anchors - a noisier
#: read, or a short one, still places at 12.
FALLBACK_K = 12

#: A diagonal needs this many votes before the clone is called placed. Two
#: unrelated plasmids sharing a backbone can produce a handful.
MIN_VOTES = 8

#: The largest between-anchor window handed to the pairwise aligner. Beyond
#: this the region is reported as unalignable rather than the app hanging on a
#: quadratic algorithm; a window this big means the clone is not the plasmid.
MAX_WINDOW = 6000

COMPLEMENT = str.maketrans("ACGTNRYSWKMBDHVacgtnryswkmbdhv", "TGCANYRSWMKVHDBtgcanyrswmkvhdb")

#: What the clone had at a reference base, when it had nothing. Three different
#: things, and telling them apart is the whole job: a gap is a deletion the
#: clone really carries, a dot is reference no read covered, and a query mark
#: is a stretch too divergent to align - which is a statement about the clone,
#: not a missing measurement.
GAP = "-"
UNCOVERED = "."
UNALIGNED = "?"


def revcomp(sequence: str) -> str:
    """The reverse complement, leaving unknown characters alone."""
    return sequence.translate(COMPLEMENT)[::-1]


# --------------------------------------------------------------------------- #
# what a difference is
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class Region:
    """A named stretch of the reference - a part, or an annotated feature.

    `end` is exclusive and may exceed the length, which is how a feature that
    crosses the origin of a circular plasmid is written.
    """

    start: int
    end: int
    label: str
    kind: str = ""

    def covers(self, start: int, end: int, length: int) -> bool:
        """Whether this region overlaps `[start, end)` on a circle of `length`."""
        end = max(end, start + 1)
        for shift in (0, length, -length):
            if self.start + shift < end and start < self.end + shift:
                return True
        return False


@dataclass(frozen=True)
class Difference:
    """One place a clone is not the plasmid it was meant to be."""

    kind: str
    """``substitution``, ``insertion``, ``deletion`` or ``unaligned``."""

    start: int
    """0-based reference position. For an insertion, the base it sits before."""

    expected: str
    found: str
    regions: tuple[str, ...] = ()

    @property
    def length(self) -> int:
        return len(self.expected) if self.kind != "insertion" else len(self.found)

    @property
    def end(self) -> int:
        return self.start + (len(self.expected) if self.kind != "insertion" else 0)

    @property
    def where(self) -> str:
        """``2,341`` or ``2,341-2,346``, 1-based the way a person reads a map."""
        if self.end <= self.start + 1:
            return f"{self.start + 1:,}"
        return f"{self.start + 1:,}-{self.end:,}"

    def describe(self, most: int = 3) -> str:
        """One line, in the terms the person has to make a decision in.

        A large deletion can cross a dozen features, and naming all twelve
        turns the one line that has to be readable into a paragraph. The full
        list is still on `regions` for anything that wants it.
        """
        where = self.where
        named = list(self.regions[:most])
        if len(self.regions) > most:
            named.append(f"{len(self.regions) - most} more")
        inside = f" in {', '.join(named)}" if named else ""
        if self.kind == "substitution":
            return f"{self.expected}→{self.found} at {where}{inside}"
        if self.kind == "deletion":
            return f"{self.length} bp deleted at {where}{inside}"
        if self.kind == "insertion":
            return f"{len(self.found)} bp inserted before {where}{inside}"
        # Name the limit. Without it this reads as a failure of the app rather
        # than a statement about the clone, and the number is the one thing
        # that says which it is: a stretch this long is not a clone with
        # mismatches in it, it is a different construct.
        return (
            f"{self.length:,} bp could not be aligned at {where}{inside}: the clone "
            f"diverges from the reference over more than {MAX_WINDOW:,} bp in one "
            f"stretch, which is past what this app will align base by base"
        )


@dataclass
class Placement:
    """How the clone's own coordinates sit on the reference."""

    strand: str = "+"
    offset: int = 0
    """The reference position the clone's first base corresponds to."""

    votes: int = 0
    placed: bool = True

    @property
    def flipped(self) -> bool:
        return self.strand == "-"


@dataclass
class CloneAlignment:
    """One sequenced clone, placed on the reference and read off base by base."""

    name: str
    length: int
    placement: Placement
    at: list[str] = field(default_factory=list)
    """One character per reference base: the clone's base, `-`, or `.`."""

    ins: dict[int, str] = field(default_factory=dict)
    """What the clone carries *before* each reference position."""

    differences: list[Difference] = field(default_factory=list)
    note: str = ""

    # -- counts ------------------------------------------------------------ #

    @property
    def covered(self) -> int:
        return sum(1 for c in self.at if c != UNCOVERED)

    @property
    def substitutions(self) -> int:
        return sum(1 for d in self.differences if d.kind == "substitution")

    @property
    def insertions(self) -> int:
        return sum(1 for d in self.differences if d.kind == "insertion")

    @property
    def deletions(self) -> int:
        return sum(1 for d in self.differences if d.kind == "deletion")

    @property
    def changed_bases(self) -> int:
        return sum(d.length for d in self.differences)

    @property
    def clean(self) -> bool:
        return self.placement.placed and not self.differences

    def identity(self, reference: str) -> float:
        """Matching bases over aligned columns, as a percentage."""
        columns = self.covered + sum(len(v) for v in self.ins.values())
        if not columns:
            return 0.0
        same = sum(
            1
            for i, c in enumerate(self.at)
            if c != UNCOVERED and i < len(reference) and c == reference[i]
        )
        return round(100.0 * same / columns, 3)

    def verdict(self) -> str:
        """The sentence that goes next to the clone's name."""
        if not self.placement.placed:
            return "does not match the reference"
        if not self.differences:
            return "no mismatches"
        n = len(self.differences)
        kinds = []
        if self.substitutions:
            kinds.append(f"{self.substitutions} substitution" + ("s" if self.substitutions > 1 else ""))
        if self.deletions:
            kinds.append(f"{self.deletions} deletion" + ("s" if self.deletions > 1 else ""))
        if self.insertions:
            kinds.append(f"{self.insertions} insertion" + ("s" if self.insertions > 1 else ""))
        if any(d.kind == "unaligned" for d in self.differences):
            kinds.append("an unalignable stretch")
        return f"{n} difference" + ("s" if n > 1 else "") + " - " + ", ".join(kinds)


# --------------------------------------------------------------------------- #
# placing the clone on the reference
# --------------------------------------------------------------------------- #


def _unique_kmers(sequence: str, k: int, circular: bool = False) -> dict[str, int]:
    """Every k-mer that occurs exactly once, mapped to where it starts.

    Repeats are dropped rather than resolved. A plasmid carries real repeats -
    two copies of the same terminator, a duplicated connector - and a k-mer
    inside one of them cannot say which copy it came from. The stretches
    between the unambiguous anchors get aligned properly instead of guessed at.
    """
    text = sequence + sequence[: k - 1] if circular else sequence
    seen: dict[str, int] = {}
    repeated: set[str] = set()
    for i in range(len(text) - k + 1):
        if circular and i >= len(sequence):
            break
        kmer = text[i : i + k]
        if kmer in seen:
            repeated.add(kmer)
        else:
            seen[kmer] = i
    for kmer in repeated:
        del seen[kmer]
    return seen


def place(reference: str, clone: str, k: int = ANCHOR_K) -> Placement:
    """Which strand the clone was read on, and where its first base sits.

    Every shared k-mer votes for a diagonal - the constant offset between the
    two coordinate systems. The right answer wins by orders of magnitude, so
    this is a count rather than a search.
    """
    length = len(reference)
    if not length or not clone:
        return Placement(placed=False)

    index = _unique_kmers(reference.upper(), k, circular=True)
    best = Placement(placed=False)

    for strand in ("+", "-"):
        text = (clone if strand == "+" else revcomp(clone)).upper()
        votes: dict[int, int] = {}
        for i in range(len(text) - k + 1):
            at = index.get(text[i : i + k])
            if at is None:
                continue
            diagonal = (at - i) % length
            votes[diagonal] = votes.get(diagonal, 0) + 1
        if not votes:
            continue
        diagonal, count = max(votes.items(), key=lambda kv: kv[1])
        if count > best.votes:
            best = Placement(strand=strand, offset=diagonal, votes=count, placed=True)

    if best.votes < MIN_VOTES:
        if k > FALLBACK_K:
            # a short clone, or a noisy one: fewer long k-mers survive intact
            return place(reference, clone, k=FALLBACK_K)
        return Placement(strand=best.strand, offset=best.offset, votes=best.votes, placed=False)
    return best


def _seam(reference: str, clone: str, placement: Placement, k: int) -> tuple[int, int] | None:
    """A shared k-mer on the winning diagonal, to open both circles at.

    Cutting a circle anywhere else risks inventing an indel at the join: the
    diagonal is an estimate, and an estimate three bases out opens the
    alignment with a three-base gap that no clone actually carries.
    """
    length = len(reference)
    index = _unique_kmers(reference.upper(), k, circular=True)
    text = (clone if placement.strand == "+" else revcomp(clone)).upper()
    for i in range(len(text) - k + 1):
        at = index.get(text[i : i + k])
        if at is not None and (at - i) % length == placement.offset:
            return at, i
    return None


# --------------------------------------------------------------------------- #
# aligning between anchors
# --------------------------------------------------------------------------- #


def _aligner() -> Align.PairwiseAligner:
    """Scoring for two copies of the same molecule, not two homologues.

    Gaps are made expensive to open and cheap to extend, so a real dropout
    comes back as one event rather than a scatter of small ones, and a
    substitution is never re-explained as an insertion beside a deletion.
    """
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2.0
    aligner.mismatch_score = -3.0
    aligner.open_gap_score = -8.0
    aligner.extend_gap_score = -0.5
    # end gaps are left on the normal gap scores: in global mode that is the
    # default, and naming them explicitly only picks a fight with whichever
    # Biopython release last renamed the properties
    return aligner


def _columns(target: str, query: str) -> list[tuple[str, str]]:
    """Column pairs for two stretches, using a real alignment only when needed."""
    if not target and not query:
        return []
    if target == query:
        return list(zip(target, query))
    if not target:
        return [(GAP, c) for c in query]
    if not query:
        return [(c, GAP) for c in target]
    if len(target) == len(query) and len(target) <= 4:
        # too short to align meaningfully; equal length means substitutions
        return list(zip(target, query))
    if len(target) > MAX_WINDOW or len(query) > MAX_WINDOW:
        raise _TooBig(len(target), len(query))

    alignment = _aligner().align(target, query)[0]
    return list(zip(str(alignment[0]), str(alignment[1])))


class _TooBig(Exception):
    """A between-anchor window no pairwise alignment should be asked to do."""

    def __init__(self, target: int, query: int) -> None:
        super().__init__(f"{target} x {query}")
        self.target = target
        self.query = query


def _chain(reference: str, clone: str, k: int) -> list[tuple[int, int, int]]:
    """Collinear exact anchors as (ref start, clone start, length), in order.

    Shared unique k-mers, kept only where they advance in both sequences, then
    merged where they overlap - which turns a run of k-mers stepping one base
    at a time into the single long exact match it describes.
    """
    index = _unique_kmers(reference, k)
    hits = [
        (index[clone[i : i + k]], i)
        for i in range(len(clone) - k + 1)
        if clone[i : i + k] in index
    ]

    chained: list[tuple[int, int, int]] = []
    for ref_at, clone_at in hits:
        if chained:
            r, c, n = chained[-1]
            if ref_at - r == clone_at - c and ref_at <= r + n:
                chained[-1] = (r, c, max(n, ref_at - r + k))
                continue
            if ref_at < r + n or clone_at < c + n:
                continue  # out of order: a repeat the uniqueness filter let through
        chained.append((ref_at, clone_at, k))
    return chained


def compare(
    reference: str,
    clone: str,
    name: str = "clone",
    regions: Sequence[Region] = (),
    k: int = ANCHOR_K,
) -> CloneAlignment:
    """Place one clone on the reference and read off every difference."""
    reference = reference.upper()
    placement = place(reference, clone, k=k)
    result = CloneAlignment(
        name=name,
        length=len(clone),
        placement=placement,
        at=[UNCOVERED] * len(reference),
    )
    if not placement.placed:
        result.note = (
            "no stretch of this sequence matches the reference on either strand"
        )
        return result

    seam = _seam(reference, clone, placement, k) or _seam(
        reference, clone, placement, FALLBACK_K
    )
    if seam is None:  # pragma: no cover - placement implies a seam exists
        result.placement = Placement(
            strand=placement.strand, offset=placement.offset, votes=placement.votes,
            placed=False,
        )
        result.note = "matched the reference but could not be opened at a shared point"
        return result

    ref_at, clone_at = seam
    turned = (clone if placement.strand == "+" else revcomp(clone)).upper()
    ref_rotated = reference[ref_at:] + reference[:ref_at]
    clone_rotated = turned[clone_at:] + turned[:clone_at]

    anchors = _chain(ref_rotated, clone_rotated, k)

    # walk the two strings, taking anchors free and aligning what falls between
    pairs: list[tuple[str, str]] = []
    r = c = 0
    for ref_start, clone_start, span in anchors:
        if ref_start < r or clone_start < c:
            continue
        try:
            pairs.extend(_columns(ref_rotated[r:ref_start], clone_rotated[c:clone_start]))
        except _TooBig:
            pairs.extend((UNALIGNED, UNALIGNED) for _ in range(ref_start - r))
        pairs.extend(zip(ref_rotated[ref_start : ref_start + span],
                         clone_rotated[clone_start : clone_start + span]))
        r, c = ref_start + span, clone_start + span
    try:
        pairs.extend(_columns(ref_rotated[r:], clone_rotated[c:]))
    except _TooBig:
        pairs.extend((UNALIGNED, UNALIGNED) for _ in range(len(ref_rotated) - r))

    _read_off(result, pairs, ref_at, len(reference))
    result.differences = _differences(result, reference, regions)
    return result


def _read_off(
    result: CloneAlignment, pairs: Iterable[tuple[str, str]], rotation: int, length: int
) -> None:
    """Turn alignment columns into one clone base per reference base.

    Everything downstream - the difference list, the viewer, the summary -
    reads this and only this, so the rotation is undone exactly once, here.
    """
    i = 0
    for ref_char, clone_char in pairs:
        position = (i + rotation) % length
        if ref_char == GAP:
            result.ins[position] = result.ins.get(position, "") + clone_char
            continue
        result.at[position] = clone_char
        i += 1


def _regions_at(regions: Sequence[Region], start: int, end: int, length: int) -> tuple[str, ...]:
    labels = [r.label for r in regions if r.covers(start, end, length)]
    seen: list[str] = []
    for label in labels:
        if label not in seen:
            seen.append(label)
    return tuple(seen)


def _differences(
    result: CloneAlignment, reference: str, regions: Sequence[Region]
) -> list[Difference]:
    """Group the per-base record into events, each named by where it lands.

    Consecutive changed bases are one event: six substitutions in a row is a
    six-base change to look at once, not six things to work through.
    """
    length = len(reference)
    out: list[Difference] = []

    def flush(kind: str, start: int, expected: str, found: str) -> None:
        out.append(
            Difference(
                kind=kind,
                start=start,
                expected=expected,
                found=found,
                regions=_regions_at(regions, start, start + max(len(expected), 1), length),
            )
        )

    run_kind = ""
    run_start = 0
    expected = found = ""

    for i in range(length):
        got = result.at[i]
        if got == UNCOVERED:
            kind = ""
        elif got == UNALIGNED:
            kind = "unaligned"
        elif got == GAP:
            kind = "deletion"
        elif got != reference[i]:
            kind = "substitution"
        else:
            kind = ""

        if kind != run_kind or (kind and i != run_start + len(expected)):
            if run_kind:
                flush(run_kind, run_start, expected, found)
            run_kind, run_start, expected, found = kind, i, "", ""
        if kind:
            expected += reference[i]
            found += got if got != GAP else ""

    if run_kind:
        flush(run_kind, run_start, expected, found)

    for position, inserted in sorted(result.ins.items()):
        out.append(
            Difference(
                kind="insertion",
                start=position,
                expected="",
                found=inserted,
                regions=_regions_at(regions, position, position + 1, length),
            )
        )

    out.sort(key=lambda d: (d.start, d.kind))
    return out


# --------------------------------------------------------------------------- #
# the viewer's rows
# --------------------------------------------------------------------------- #


@dataclass
class Rows:
    """One multiple alignment: the reference and every clone, same width."""

    reference: str
    rows: list[tuple[str, str]] = field(default_factory=list)
    insertions: list[tuple[int, int]] = field(default_factory=list)
    """``(reference position, width)`` for each column group that is an
    insertion, so a viewer can label columns without being sent an index."""

    @property
    def columns(self) -> int:
        return len(self.reference)


def rows(reference: str, clones: Sequence[CloneAlignment]) -> Rows:
    """Merge the per-clone records into one alignment every row shares.

    Each clone is placed against the reference independently, so they disagree
    about how many columns an insertion needs. The merge takes the widest and
    pads the rest, which is what makes the columns on screen mean the same
    thing on every row.
    """
    reference = reference.upper()
    width: dict[int, int] = {}
    for clone in clones:
        for position, inserted in clone.ins.items():
            width[position] = max(width.get(position, 0), len(inserted))

    ref_row: list[str] = []
    built: list[list[str]] = [[] for _ in clones]
    gaps: list[tuple[int, int]] = []

    for i in range(len(reference)):
        pad = width.get(i, 0)
        if pad:
            gaps.append((i, pad))
            ref_row.append(GAP * pad)
            for n, clone in enumerate(clones):
                built[n].append(clone.ins.get(i, "").ljust(pad, GAP))
        ref_row.append(reference[i])
        for n, clone in enumerate(clones):
            built[n].append(clone.at[i] if i < len(clone.at) else UNCOVERED)

    return Rows(
        reference="".join(ref_row),
        rows=[(clone.name, "".join(built[n])) for n, clone in enumerate(clones)],
        insertions=gaps,
    )


# --------------------------------------------------------------------------- #
# regions from a GenBank record
# --------------------------------------------------------------------------- #

#: Feature types worth naming a difference by. A `source` covers the whole
#: plasmid and would label every difference with it, which says nothing.
NAMED_FEATURES = {
    "CDS", "gene", "promoter", "terminator", "misc_feature", "rep_origin",
    "regulatory", "mobile_element", "protein_bind", "primer_bind", "RBS",
    "sig_peptide", "misc_RNA", "ncRNA", "tRNA", "rRNA", "oriT", "STS",
}

#: Where a label lives, in the order worth trying.
LABEL_KEYS = ("label", "gene", "product", "note", "standard_name", "locus_tag")


def regions_from_record(record: SeqRecord, limit: int = 400) -> list[Region]:
    """The annotated features of a reference, as regions a difference can name.

    Annotations are used here and nowhere else in this app: naming *where* a
    difference falls is a convenience, and being wrong about it costs a
    misleading label, not a wrong construct. What a plasmid *is* still comes
    from the digest.
    """
    out: list[Region] = []
    for feature in record.features[: limit * 4]:
        if feature.type not in NAMED_FEATURES:
            continue
        label = ""
        for key in LABEL_KEYS:
            values = feature.qualifiers.get(key)
            if values:
                label = str(values[0]).strip()
                break
        if not label:
            continue
        try:
            start = int(feature.location.start)
            end = int(feature.location.end)
        except (TypeError, ValueError):
            continue
        if end <= start:
            continue
        out.append(Region(start=start, end=end, label=label, kind=feature.type))
        if len(out) >= limit:
            break
    return out
