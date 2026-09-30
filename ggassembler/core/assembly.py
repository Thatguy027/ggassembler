"""Generic Golden Gate simulator: a set of fragments to one circular product.

One function, `assemble`, serves all three levels. It knows nothing about part
types, cassettes or connectors - only about overhangs. Give it fragments and
the enzyme, and it builds the junction graph, walks it, and either closes one
circle or tells you exactly which overhang stopped it.

Fragment convention
-------------------
A fragment's sequence **includes its left overhang and excludes its right**,
which is the convention `enzymes.digest` produces. That makes the product the
plain concatenation of the fragments in order, with every base accounted for
exactly once, and makes each junction's overhang the first four bases of the
downstream fragment.

A biological problem is never an exception. Everything wrong with a reaction
comes back as an `Issue` in the result, so a UI can show all of them at once
rather than the first one that raised.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from .enzymes import (
    Enzyme,
    find_sites,
    is_palindromic,
    shares_noncontiguous_three,
    shares_three_of_four,
)

ERROR, WARNING, INFO = "error", "warning", "info"


@dataclass
class Issue:
    """Something worth telling the user about a reaction."""

    level: str
    code: str
    message: str

    def __str__(self) -> str:
        return f"[{self.level}] {self.code}: {self.message}"


@dataclass
class Piece:
    """One fragment going into a reaction."""

    record: SeqRecord
    left_overhang: str
    right_overhang: str
    source_name: str
    part_type: str | None = None
    color: str | None = None
    component: str = ""
    """What this fragment is - ``ScTDH3 Promoter`` - as opposed to which plasmid
    it came out of. This is what the product gets labelled with."""

    @property
    def label(self) -> str:
        """What to call this fragment on a map or in a GenBank file."""
        return self.component or self.source_name

    @property
    def seq(self) -> str:
        return str(self.record.seq).upper()

    def __len__(self) -> int:
        return len(self.record.seq)


@dataclass
class PlacedPart:
    """Where a piece ended up in the product."""

    name: str
    source_name: str
    part_type: str | None
    start: int
    end: int
    length: int
    left_overhang: str
    right_overhang: str
    color: str | None = None
    component: str = ""

    @property
    def label(self) -> str:
        return self.component or self.source_name


@dataclass
class Junction:
    """Where two pieces were joined, and by which overhang."""

    overhang: str
    upstream: str
    downstream: str
    position: int


@dataclass
class AssemblyResult:
    """Everything a caller needs to show, export or diagnose a reaction."""

    ok: bool
    product: SeqRecord | None
    parts: list[PlacedPart] = field(default_factory=list)
    junctions: list[Junction] = field(default_factory=list)
    issues: list[Issue] = field(default_factory=list)
    enzyme: str = ""

    @property
    def errors(self) -> list[Issue]:
        return [i for i in self.issues if i.level == ERROR]

    @property
    def warnings(self) -> list[Issue]:
        return [i for i in self.issues if i.level == WARNING]

    @property
    def length(self) -> int:
        return len(self.product.seq) if self.product else 0

    def codes(self) -> set[str]:
        return {i.code for i in self.issues}


def assemble(
    pieces: list[Piece],
    enzyme: Enzyme,
    name: str = "assembly",
    circular: bool = True,
) -> AssemblyResult:
    """Ligate `pieces` into one circle, or explain why they do not close.

    The order pieces arrive in does not matter; the junction graph decides it.
    """
    issues: list[Issue] = []

    if not pieces:
        return AssemblyResult(False, None, issues=[Issue(ERROR, "no_parts", "no fragments given")],
                              enzyme=enzyme.name)

    issues += _check_overhangs(pieces, enzyme)
    issues += _check_internal_sites(pieces, enzyme)

    order, walk_issues = _walk(pieces)
    issues += walk_issues
    if order is None:
        return AssemblyResult(False, None, issues=issues, enzyme=enzyme.name)

    issues += _check_near_matches(order)

    product, parts, junctions = _build(order, enzyme, name, circular)
    ok = not any(i.level == ERROR for i in issues)
    return AssemblyResult(
        ok=ok,
        product=product if ok else None,
        parts=parts,
        junctions=junctions,
        issues=issues,
        enzyme=enzyme.name,
    )


# --------------------------------------------------------------------------- #
# validation
# --------------------------------------------------------------------------- #


def _check_overhangs(pieces: list[Piece], enzyme: Enzyme) -> list[Issue]:
    issues: list[Issue] = []
    width = enzyme.overhang_length

    for piece in pieces:
        for side, overhang in (("5'", piece.left_overhang), ("3'", piece.right_overhang)):
            if not overhang or len(overhang) != width:
                issues.append(
                    Issue(ERROR, "bad_overhang",
                          f"{piece.source_name}: {side} overhang {overhang!r} is not "
                          f"{width} nt, as {enzyme.name} would leave")
                )
            elif set(overhang) - set("ACGT"):
                issues.append(
                    Issue(ERROR, "bad_overhang",
                          f"{piece.source_name}: {side} overhang {overhang!r} is not plain DNA")
                )
            elif is_palindromic(overhang):
                issues.append(
                    Issue(ERROR, "palindromic_overhang",
                          f"{piece.source_name}: {side} overhang {overhang} is palindromic "
                          f"and would ligate to itself in either orientation")
                )

    for label, getter in (("donor", lambda p: p.right_overhang), ("acceptor", lambda p: p.left_overhang)):
        seen: dict[str, list[str]] = {}
        for piece in pieces:
            seen.setdefault(getter(piece), []).append(piece.source_name)
        for overhang, owners in seen.items():
            if len(owners) > 1:
                issues.append(
                    Issue(ERROR, "duplicate_overhang",
                          f"overhang {overhang} is used as a {label} by "
                          f"{len(owners)} fragments ({', '.join(owners)}); each overhang "
                          f"must appear exactly once in each direction")
                )

    donors = {p.right_overhang for p in pieces}
    acceptors = {p.left_overhang for p in pieces}
    for overhang in sorted(donors - acceptors):
        owner = next(p.source_name for p in pieces if p.right_overhang == overhang)
        issues.append(
            Issue(ERROR, "missing_overhang",
                  f"nothing accepts overhang {overhang}, left by {owner}'s 3' end")
        )
    for overhang in sorted(acceptors - donors):
        owner = next(p.source_name for p in pieces if p.left_overhang == overhang)
        issues.append(
            Issue(ERROR, "missing_overhang",
                  f"nothing donates overhang {overhang}, needed by {owner}'s 5' end")
        )

    return issues


def _check_internal_sites(pieces: list[Piece], enzyme: Enzyme) -> list[Issue]:
    issues = []
    for piece in pieces:
        sites = find_sites(piece.seq, enzyme, circular=False)
        if sites:
            issues.append(
                Issue(WARNING, "internal_site",
                      f"{piece.source_name} carries {len(sites)} internal {enzyme.name} "
                      f"site(s); it will be cut again and the product may not survive")
            )
    return issues


def _check_near_matches(order: list[Piece]) -> list[Issue]:
    """Overhang pairs close enough to misligate.

    The SI names these as the modes the YTK overhang set was chosen to avoid:
    a 3-of-4 match, and a 3-nt match that is not contiguous.
    """
    issues = []
    overhangs = [(p.left_overhang, p.source_name) for p in order]
    for i, (a, owner_a) in enumerate(overhangs):
        for b, owner_b in overhangs[i + 1:]:
            if not shares_three_of_four(a, b):
                continue
            if shares_noncontiguous_three(a, b):
                issues.append(
                    Issue(WARNING, "near_match_noncontiguous",
                          f"{a} and {b} share three bases either side of one mismatch "
                          f"({owner_a} / {owner_b}); they can misligate")
                )
            else:
                issues.append(
                    Issue(WARNING, "near_match_3of4",
                          f"{a} and {b} differ at one end position only "
                          f"({owner_a} / {owner_b}); they can misligate")
                )
    return issues


# --------------------------------------------------------------------------- #
# the junction graph
# --------------------------------------------------------------------------- #


def _walk(pieces: list[Piece]) -> tuple[list[Piece] | None, list[Issue]]:
    """Order the pieces by following overhangs, or say where the chain breaks."""
    by_left: dict[str, Piece] = {}
    for piece in pieces:
        by_left.setdefault(piece.left_overhang, piece)

    start = pieces[0]
    order = [start]
    seen = {id(start)}
    current = start

    while True:
        nxt = by_left.get(current.right_overhang)
        if nxt is None:
            return None, [
                Issue(ERROR, "open_end",
                      f"the chain stops at {current.source_name}: no fragment starts "
                      f"with {current.right_overhang}")
            ]
        if id(nxt) in seen:
            if nxt is start and len(order) == len(pieces):
                return order, []
            if nxt is start:
                unused = [p.source_name for p in pieces if id(p) not in seen]
                return None, [
                    Issue(ERROR, "multiple_circles",
                          f"{len(order)} of {len(pieces)} fragments close a circle on their "
                          f"own; left over: {', '.join(unused)}")
                ]
            return None, [
                Issue(ERROR, "branch", f"{nxt.source_name} is reached twice")
            ]
        order.append(nxt)
        seen.add(id(nxt))
        current = nxt


# --------------------------------------------------------------------------- #
# the product
# --------------------------------------------------------------------------- #


def _build(
    order: list[Piece], enzyme: Enzyme, name: str, circular: bool
) -> tuple[SeqRecord, list[PlacedPart], list[Junction]]:
    sequence = "".join(p.seq for p in order)
    product = SeqRecord(
        Seq(sequence),
        id=name,
        name=name[:16],
        description=f"{name}, assembled in silico with {enzyme.name}",
        annotations={
            "molecule_type": "DNA",
            "topology": "circular" if circular else "linear",
            "source": "synthetic DNA construct",
            "comment": (
                f"Assembled by GG Assembler with {enzyme.name} from: "
                + ", ".join(p.source_name for p in order)
            ),
        },
    )

    parts: list[PlacedPart] = []
    junctions: list[Junction] = []
    offset = 0
    for index, piece in enumerate(order):
        end = offset + len(piece)
        parts.append(
            PlacedPart(
                name=piece.record.id or piece.source_name,
                source_name=piece.source_name,
                part_type=piece.part_type,
                start=offset,
                end=end,
                length=len(piece),
                left_overhang=piece.left_overhang,
                right_overhang=piece.right_overhang,
                color=piece.color,
                component=piece.component,
            )
        )

        for feature in piece.record.features:
            product.features.append(_shift(feature, offset))

        # Label the span with what the DNA is, not which tube it came from; the
        # source plasmid stays in the note, so provenance is never lost.
        #
        # The note says "position 2", never "type 2": this file will be read back
        # by the detector, and a note reading "type 2 part from ..." would be
        # taken as an annotation claiming the *whole plasmid* is a type 2 part.
        kind = f"position {piece.part_type} " if piece.part_type else ""
        product.features.append(
            SeqFeature(
                SimpleLocation(offset, end, 1),
                type="misc_feature",
                qualifiers={
                    "label": [piece.label],
                    "note": [f"{kind}part from {piece.source_name}".strip()],
                },
            )
        )

        downstream = order[(index + 1) % len(order)]
        junctions.append(
            Junction(
                overhang=piece.right_overhang,
                upstream=piece.label,
                downstream=downstream.label,
                position=end % len(sequence) if sequence else 0,
            )
        )
        offset = end

    for junction in junctions:
        product.features.append(_junction_feature(junction, len(sequence)))

    return product, parts, junctions


def _shift(feature: SeqFeature, offset: int) -> SeqFeature:
    shifted = feature.location + offset if offset else feature.location
    return SeqFeature(shifted, type=feature.type,
                      qualifiers={k: list(v) for k, v in feature.qualifiers.items()})


def _junction_feature(junction: Junction, total: int) -> SeqFeature:
    """A 4-nt marker on the overhang itself, naming the two parts it joins."""
    start = junction.position
    end = start + len(junction.overhang)
    if end > total:  # the junction at the origin
        start, end = 0, len(junction.overhang)
    return SeqFeature(
        SimpleLocation(start, end, 1),
        type="misc_feature",
        qualifiers={
            "label": [f"junction {junction.overhang}"],
            "note": [f"{junction.upstream} -> {junction.downstream} via {junction.overhang}"],
        },
    )


def redigest(product: SeqRecord, enzyme: Enzyme) -> list[str]:
    """The overhangs a finished product would give back, for round-trip checks."""
    from .enzymes import digest

    return [f.left_overhang or "" for f in digest(str(product.seq), enzyme)]
