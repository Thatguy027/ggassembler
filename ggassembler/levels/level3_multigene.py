"""Level 3 - multigene assembly.

Finished cassettes are joined into one plasmid by the scheme's multigene enzyme
(BsmBI in YTK). What orders them is the connectors: each cassette was built with
a type 1 part at its left and a type 5 part at its right, and each of those
carries a connector-specific overhang. Cassette *i*'s right connector has to be
met by cassette *i+1*'s left, and the two ends of the chain have to meet the
destination vector's backbone.

Those overhangs are read from the files, never assumed - the same rule as
everywhere else in this package. So this module does no arithmetic on
connectors at all: it collects the fragments each plasmid releases and hands
them to the generic simulator, which either closes one circle or names the
overhang that stopped it.

This module imports no other level.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Callable, Sequence

from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqUtils import MeltingTemp

from ..core import decompose as core_decompose
from ..core import seqio
from ..core.assembly import ERROR, INFO, WARNING, AssemblyResult, Issue, Piece, assemble
from ..core.enzymes import digest
from ..core.library import Library, PlasmidEntry
from ..core.parttypes import YTK, Scheme, color_for

#: How much of each side of a junction the connector barcode scar covers.
SCAR = 20

#: Colony-PCR primers verifying each junction.
CHECK_PRIMER = 20

#: ConS, ConE, Con1..Con9, and the ConLS / ConRE spellings of the same things.
_CONNECTOR_NAME = re.compile(r"\bCon(?:L|R)?[A-Z0-9][\u2032']?\b", re.IGNORECASE)

#: Labels that describe the vector or the provenance of a part rather than the
#: cargo. A cassette's dropdown should say what gene it carries, not that it has
#: an origin or that it came from pYTK013.
_BACKBONE_NOISE = re.compile(
    r"^(amp|kan|cam|cm|spec|tet)r\b|ampr|cole1|cen6|ars4|\bori\b|"
    r"^promoter$|^terminator$|^con[lrs0-9]|scar|^pytk\d+$|\bfp$|\brp$|ribosome",
    re.IGNORECASE,
)

#: How each feature type is drawn on the annotation map, reusing the part palette.
FEATURE_COLORS = {
    "promoter": "var(--part-2)",
    "CDS": "var(--part-3)",
    "gene": "var(--part-3)",
    "terminator": "var(--part-4)",
    "rep_origin": "var(--part-7)",
    "protein_bind": "var(--part-5)",
    "primer": "var(--part-1)",
}
DEFAULT_FEATURE_COLOR = "var(--part-8)"


@dataclass
class MultigeneDesign:
    """A Level 3 screen's state: a backbone and an ordered list of cassettes."""

    backbone: str | None = None
    transcription_units: list[str] = field(default_factory=list)
    name: str = "multigene"

    @property
    def count(self) -> int:
        return len(self.transcription_units)


@dataclass
class CheckPrimer:
    name: str
    sequence: str
    tm: float
    position: int


@dataclass
class MultigeneResult:
    """Everything the Level 3 screen shows."""

    assembly: AssemblyResult
    design: MultigeneDesign
    chain: list[dict[str, str]] = field(default_factory=list)
    suggestions: list[dict[str, str]] = field(default_factory=list)
    scars: list[dict[str, object]] = field(default_factory=list)
    check_primers: list[CheckPrimer] = field(default_factory=list)
    integration: bool = False

    @property
    def ok(self) -> bool:
        return self.assembly.ok

    @property
    def issues(self) -> list[Issue]:
        return self.assembly.issues

    def codes(self) -> set[str]:
        return self.assembly.codes()


# --------------------------------------------------------------------------- #
# what the library offers
# --------------------------------------------------------------------------- #


def connector_names(library: Library) -> dict[str, dict[str, list[str]]]:
    """What each connector overhang is called, learned from the library.

    Two sources, both the user's own files. First, which type 1 and type 5 parts
    carry each overhang. Second - and much more useful - the ``Con1``, ``ConS``,
    ``ConE`` features that annotated cassettes already carry: for each one, the
    boundary of the released fragment it sits nearest names that overhang. A
    name is only used when the library agrees on it, so a stray label cannot
    rename a connector on its own.

    Nothing here is hard-coded. Point the app at a library annotated with some
    other naming and it will report that naming instead.
    """
    index: dict[str, dict[str, list[str]]] = {}

    def slot(overhang: str) -> dict[str, list[str]]:
        return index.setdefault(overhang, {"left": [], "right": [], "names": []})

    for entry in library.sorted_entries():
        if entry.connector_overhang:
            side = "left" if entry.call.part_type == "1" else "right"
            slot(entry.connector_overhang)[side].append(entry.name)

    votes: dict[str, dict[str, int]] = {}
    for entry in library.sorted_entries():
        span, overhangs = entry.cassette_span, entry.cassette_overhangs
        if not span or not overhangs or not entry.length:
            continue
        start, end = span
        for feature in entry.features:
            label = feature.label or ""
            match = _CONNECTOR_NAME.search(label)
            if not match or "scar" in label.lower():
                continue
            to_left = min((feature.end - start) % entry.length,
                          (start - feature.start) % entry.length)
            to_right = min((feature.end - end) % entry.length,
                           (end - feature.start) % entry.length)
            overhang = overhangs[0] if to_left <= to_right else overhangs[1]
            tally = votes.setdefault(overhang, {})
            tally[match.group(0)] = tally.get(match.group(0), 0) + 1

    for overhang, tally in votes.items():
        ranked = sorted(tally.items(), key=lambda pair: -pair[1])
        best, count = ranked[0]
        runner_up = ranked[1][1] if len(ranked) > 1 else 0
        if count >= 2 and count >= 2 * runner_up:
            slot(overhang)["names"].insert(0, best)

    return index


def describe_connector(overhang: str, index: dict[str, dict[str, list[str]]]) -> str:
    """A human label for a connector overhang, from whatever the files say."""
    slot = index.get(overhang)
    if not slot:
        return overhang
    if slot["names"]:
        return f"{overhang} ({'/'.join(slot['names'][:2])})"
    carriers = (slot["left"] + slot["right"])[:2]
    return f"{overhang} ({'/'.join(carriers)})" if carriers else overhang


def unit_features(entry: PlasmidEntry) -> list:
    """The features inside the released fragment - the cargo, not the vector."""
    if not entry.cassette_span:
        return list(entry.features)
    start, end = entry.cassette_span
    inside = []
    for feature in entry.features:
        offset = (feature.start - start) % entry.length
        length = (feature.end - feature.start) % entry.length or (feature.end - feature.start)
        span = (end - start) % entry.length or entry.length
        if offset + length <= span:
            inside.append(feature)
    return inside


def contents(entry: PlasmidEntry) -> str:
    """A one-line description of what a cassette carries, for the dropdown.

    Reads the features the plasmid already has, inside the fragment the reaction
    will actually move, and reports them the way a transcription unit is read:
    promoter, then whatever it drives, then terminator. Labels are the user's
    own - nothing here renames their parts.
    """
    inside = [
        f for f in unit_features(entry)
        if f.label and not _BACKBONE_NOISE.search(f.label)
    ]
    if not inside:
        return ""

    def first(kind: str) -> object | None:
        typed = [f for f in inside if f.type == kind]
        if typed:
            return typed[0]
        return next((f for f in inside if kind in f.label.lower()), None)

    promoter = first("promoter")
    terminator = first("terminator")
    skip = {id(f) for f in (promoter, terminator) if f}
    cargo = max(
        (f for f in inside if id(f) not in skip),
        key=lambda f: f.end - f.start,
        default=None,
    )

    parts = [f.label for f in (promoter, cargo, terminator) if f]
    seen: list[str] = []
    for label in parts:
        if label not in seen:
            seen.append(label)
    return " \u00b7 ".join(seen[:3])


def _within(position: int, span: tuple[int, int] | None, length: int) -> bool:
    """Circular containment: a span may run past the origin and wrap back."""
    if not span or not length:
        return False
    start, end = span
    return ((position - start) % length) < ((end - start) % length or length)


def feature_map(entry: PlasmidEntry) -> list[dict[str, object]]:
    """Every annotated feature, ready to draw as a map of the whole plasmid.

    ``kept`` marks the features on the fragment this plasmid contributes. For a
    destination vector that is the arm that survives the reaction; everything
    else is in the dropout, which the chosen cassettes replace.
    """
    return [
        {
            "label": feature.label or feature.type,
            "type": feature.type,
            "start": feature.start,
            "end": feature.end,
            "strand": feature.strand,
            "length": max(0, feature.end - feature.start),
            "color": FEATURE_COLORS.get(feature.type, DEFAULT_FEATURE_COLOR),
            "kept": _within(feature.start, entry.cassette_span, entry.length),
        }
        for feature in entry.features
        if feature.end > feature.start
    ]


def options(library: Library) -> dict[str, list[dict[str, object]]]:
    """The cassettes and backbones the screen can offer."""
    index = connector_names(library)

    def row(entry: PlasmidEntry, with_map: bool = False) -> dict[str, object]:
        left, right = entry.cassette_overhangs or ("", "")
        unit = entry.cassette_span
        payload: dict[str, object] = {
            "name": entry.name,
            "length": entry.length,
            "left_overhang": left,
            "right_overhang": right,
            "left_label": describe_connector(left, index),
            "right_label": describe_connector(right, index),
            "marker": entry.ecoli_marker,
            "linearizer_sites": entry.sites.linearizer_total,
            "contents": entry.component or contents(entry),
            "display": entry.display,
            "aliases": entry.aliases,
            "unit_length": (unit[1] - unit[0]) % entry.length if unit else 0,
            "unit_start": unit[0] if unit else None,
            "unit_end": unit[1] if unit else None,
            # the complement: what the reaction throws away and replaces
            "dropout_start": unit[1] if unit else None,
            "dropout_end": unit[0] if unit else None,
            "dropout_length": (unit[0] - unit[1]) % entry.length if unit else 0,
            "integration": entry.sites.linearizer_total >= 2,
        }
        if with_map:
            payload["features"] = feature_map(entry)
        return payload

    return {
        "cassettes": [row(e) for e in library.cassettes()],
        "backbones": [row(e, with_map=True) for e in library.multigene_vectors()],
    }


# --------------------------------------------------------------------------- #
# building
# --------------------------------------------------------------------------- #


def build(library: Library, design: MultigeneDesign) -> MultigeneResult:
    """Assemble the chosen cassettes into the chosen backbone."""
    scheme = library.scheme
    enzyme = scheme.multigene_enzyme
    issues: list[Issue] = []
    pieces: list[Piece] = []

    backbone_entry = library.get(design.backbone) if design.backbone else None
    if backbone_entry is None:
        issues.append(Issue(ERROR, "no_backbone", "choose a destination backbone"))
    elif not backbone_entry.is_multigene_vector:
        issues.append(
            Issue(ERROR, "not_a_backbone",
                  f"{backbone_entry.name} is not a destination vector: a backbone needs "
                  f"connector ends with a dropout between them")
        )
        backbone_entry = None

    if not design.transcription_units:
        issues.append(Issue(ERROR, "no_cassettes", "add at least one transcription unit"))

    chosen: list[PlasmidEntry] = []
    for index, name in enumerate(design.transcription_units, start=1):
        if not name:
            issues.append(Issue(ERROR, "empty_slot", f"TU{index} is empty"))
            continue
        entry = library.get(name)
        if entry is None:
            issues.append(Issue(ERROR, "unknown_cassette", f"no plasmid named {name}"))
            continue
        if not entry.is_cassette:
            issues.append(
                Issue(ERROR, "not_a_cassette",
                      f"{entry.name} is not a cassette: {enzyme.name} does not release a "
                      f"transcription unit with connector ends from it")
            )
            continue
        chosen.append(entry)

    for entry in ([backbone_entry] if backbone_entry else []) + chosen:
        try:
            pieces.append(library.piece(entry, enzyme))
        except (ValueError, OSError) as exc:
            issues.append(Issue(ERROR, "unreadable", f"{entry.name}: {exc}"))

    if any(i.level == ERROR for i in issues):
        failed = AssemblyResult(False, None, issues=issues, enzyme=enzyme.name)
        return MultigeneResult(failed, design, chain=_chain(backbone_entry, chosen, library))

    result = assemble(pieces, enzyme, name=design.name)
    result.issues = issues + result.issues
    result.ok = result.ok and not result.errors
    if not result.ok:
        result.product = None

    out = MultigeneResult(
        assembly=result,
        design=design,
        chain=_chain(backbone_entry, chosen, library),
    )
    if not result.ok:
        out.suggestions = _suggest(library, result, chosen)
        return out

    _annotate_scars(out, scheme)
    _integration_notes(out, backbone_entry, scheme)
    return out


def _chain(
    backbone: PlasmidEntry | None, chosen: list[PlasmidEntry], library: Library
) -> list[dict[str, str]]:
    """The connector chain as the screen draws it, whether or not it closes."""
    index = connector_names(library)
    chain = []
    if backbone:
        left, right = backbone.cassette_overhangs or ("", "")
        chain.append(
            {
                "role": "backbone",
                "name": backbone.name,
                "left": left,
                "right": right,
                "left_label": describe_connector(left, index),
                "right_label": describe_connector(right, index),
            }
        )
    for position, entry in enumerate(chosen, start=1):
        left, right = entry.cassette_overhangs or ("", "")
        chain.append(
            {
                "role": f"TU{position}",
                "name": entry.name,
                "left": left,
                "right": right,
                "left_label": describe_connector(left, index),
                "right_label": describe_connector(right, index),
            }
        )
    return chain


def _suggest(
    library: Library, result: AssemblyResult, chosen: list[PlasmidEntry]
) -> list[dict[str, str]]:
    """When the chain breaks, offer the cassettes that would bridge the gap.

    A ConL-spacer-ConR filler is just a cassette with the right ends, so this
    needs no special case: it looks for anything in the library that starts
    where the chain stopped.
    """
    dangling = set()
    for issue in result.issues:
        if issue.code in ("missing_overhang", "open_end"):
            dangling.update(re.findall(r"\b([ACGT]{4})\b", issue.message))
    if not dangling:
        return []

    used = {e.name for e in chosen}
    out = []
    for entry in library.cassettes():
        if entry.name in used or not entry.cassette_overhangs:
            continue
        left, right = entry.cassette_overhangs
        if left in dangling:
            out.append(
                {
                    "name": entry.name,
                    "left_overhang": left,
                    "right_overhang": right,
                    "why": f"starts at {left}, where the chain stops",
                }
            )
    return out[:8]


def _annotate_scars(out: MultigeneResult, scheme: Scheme) -> None:
    """Mark the connector barcode left at every junction."""
    product = out.assembly.product
    if product is None:
        return
    total = len(product.seq)
    for junction in out.assembly.junctions:
        start = max(0, junction.position - SCAR // 2)
        end = min(total, junction.position + SCAR // 2)
        product.features.append(
            SeqFeature(
                SimpleLocation(start, end, 1),
                type="misc_feature",
                qualifiers={
                    "label": [f"connector barcode {junction.overhang}"],
                    "note": [
                        f"{SCAR} bp connector scar joining {junction.upstream} to "
                        f"{junction.downstream}"
                    ],
                },
            )
        )
        out.scars.append(
            {
                "overhang": junction.overhang,
                "start": start,
                "end": end,
                "upstream": junction.upstream,
                "downstream": junction.downstream,
            }
        )


def _integration_notes(
    out: MultigeneResult, backbone: PlasmidEntry | None, scheme: Scheme
) -> None:
    """NotI linearization, and primers that verify each junction by colony PCR."""
    product = out.assembly.product
    if product is None or backbone is None:
        return

    from ..core.enzymes import find_sites

    sequence = str(product.seq).upper()
    sites = find_sites(sequence, scheme.linearizer)
    out.integration = len(sites) >= 2
    if out.integration:
        where = ", ".join(f"{s.start:,}" for s in sites)
        out.assembly.issues.append(
            Issue(INFO, "linearize",
                  f"{backbone.name} is an integration vector: cut the finished plasmid "
                  f"with {scheme.linearizer.name} ({len(sites)} sites, at {where}) before "
                  f"transforming, so the construct integrates rather than replicating")
        )
        out.assembly.issues.append(
            Issue(INFO, "locus_primers",
                  "the paper's locus-specific colony-PCR primers are in SI Table S4, which "
                  "is not bundled with this app; the junction primers below are designed "
                  "from the construct itself and verify the assembly, not the locus")
        )

    for junction in out.assembly.junctions:
        start = (junction.position - CHECK_PRIMER) % len(sequence)
        forward = sequence[start:start + CHECK_PRIMER]
        if len(forward) < CHECK_PRIMER:
            continue
        out.check_primers.append(
            CheckPrimer(
                name=f"check_{junction.upstream}_{junction.downstream}",
                sequence=forward,
                tm=round(float(MeltingTemp.Tm_NN(Seq(forward))), 1),
                position=start,
            )
        )


def save(library: Library, result: MultigeneResult, name: str) -> str | None:
    """Write the assembled multigene plasmid into the library folder."""
    from ..core import seqio

    if result.assembly.product is None:
        return None
    path = library.folder / f"{name}.gb"
    seqio.write_genbank(result.assembly.product, path)
    library.scan()
    return str(path)


# --------------------------------------------------------------------------- #
# a worked assembly to open on
# --------------------------------------------------------------------------- #


def _unit_preference(entry: PlasmidEntry) -> tuple:
    """How good a default a cassette makes: clean, described, simple."""
    return (
        bool(entry.internal_multigene_positions),
        not (entry.component or contents(entry)),
        entry.length,
        entry.name,
    )


def default_design(library: Library, name: str = "pMultigene") -> MultigeneDesign:
    """A closed, working multigene assembly to open the screen with.

    An empty Level 3 screen reports `no_backbone` and `empty_slot` the moment
    it loads, which reads as a broken app rather than as an invitation. Level 2
    already opens on a worked assembly; this is the same idea, and the same
    rule: nothing is hard-coded to a catalogue.

    The chain has to *close* - run from the backbone's right connector back
    round to its left - or the screen still opens on errors, so this walks the
    connectors greedily and only keeps a backbone it could close. Where no
    closing chain exists the backbone is still offered on its own, which is a
    starting point rather than a failure.
    """
    design = MultigeneDesign(name=name)
    cassettes = library.cassettes()
    if not cassettes:
        return design

    # index by the connector each cassette starts at, best default first
    starting: dict[str, list[PlasmidEntry]] = {}
    for entry in sorted(cassettes, key=_unit_preference):
        left, _ = entry.cassette_overhangs or ("", "")
        if left:
            starting.setdefault(left, []).append(entry)

    # Prefer a unit that expresses something. A spacer cassette assembles just
    # as well and shows nothing, so a screen seeded with two of them teaches
    # nothing about what the screen is for.
    #
    # Reading a unit's parts costs a tiling, and tiling every cassette in a
    # real library takes seconds - so it is done on demand and remembered,
    # which touches only the handful of candidates the walk actually weighs.
    expressive: dict[str, int] = {}

    def shows(entry: PlasmidEntry) -> int:
        """How many parts this unit can be drawn as: more is a better demo.

        A promoter, a gene and a terminator make the point the diagram exists
        to make; a single 2-3-4 composite is one block and a spacer is none.
        """
        if entry.name not in expressive:
            expressive[entry.name] = len(unit_contents(library, entry))
        return expressive[entry.name]

    for backbone in sorted(library.multigene_vectors(), key=_unit_preference):
        left, right = backbone.cassette_overhangs or ("", "")
        if not left or not right:
            continue
        chain = _walk(starting, start=right, stop=left, prefer=shows)
        if chain:
            design.backbone = backbone.name
            design.transcription_units = [e.name for e in chain]
            return design

    # nothing closed; open on the simplest backbone rather than on nothing
    vectors = sorted(library.multigene_vectors(), key=_unit_preference)
    if vectors:
        design.backbone = vectors[0].name
        design.transcription_units = [""]
    return design


def _walk(
    starting: dict[str, list[PlasmidEntry]],
    start: str,
    stop: str,
    limit: int = 4,
    prefer: Callable[[PlasmidEntry], int] | None = None,
) -> list[PlasmidEntry]:
    """A chain of cassettes from `start` to `stop`, or [] if none is short enough.

    Depth-first over the connector graph, refusing to use a cassette twice, so
    it cannot loop through a repeated connector forever.
    """
    def step(at: str, used: list[PlasmidEntry]) -> list[PlasmidEntry] | None:
        if at == stop and used:
            return used
        if len(used) >= limit:
            return None
        options = starting.get(at, [])
        if prefer is not None:
            # a stable partition: the ones worth showing first, each group
            # still in its existing preference order
            options = sorted(options, key=lambda e: -prefer(e))
        for entry in options:
            if entry in used:
                continue
            _, right = entry.cassette_overhangs or ("", "")
            if not right:
                continue
            found = step(right, [*used, entry])
            if found:
                return found
        return None

    return step(start, []) or []


# --------------------------------------------------------------------------- #
# what is inside each transcription unit
# --------------------------------------------------------------------------- #

#: The positions that make up a transcription unit: promoter, coding sequence,
#: terminator. Position 1 and 5 are the connectors that join units together and
#: 6-8 are the plasmid's own machinery - neither is what the unit *expresses*.
UNIT_POSITIONS = ("2", "3", "3a", "3b", "4", "4a", "4b", "234")


def unit_contents(library: Library, entry: PlasmidEntry) -> list[dict[str, object]]:
    """The promoter, CDS and terminator inside one cassette, in reading order.

    Level 3 assembles cassettes, so a chain of cassettes is all it normally
    shows - and that hides the thing you are actually checking, which is that
    *this* promoter is driving *that* gene. The parts are still in there as
    sequence, so they can be read back out and drawn above the unit they build.

    A part that is not in the library comes back as an entry with no name,
    since a gap in a diagram of a construct is a fact about the construct.
    """
    tiling = core_decompose.tile(library, entry)
    out: list[dict[str, object]] = []

    for match in tiling.matches:
        if match.part_type not in UNIT_POSITIONS:
            continue
        out.append({
            "name": match.name,
            "display": match.display,
            "component": match.component,
            "short": short_label(match.component or match.name, match.part_type),
            "part_type": match.part_type,
            "length": match.length,
            "start": match.start,
            "end": match.end,
            "color": color_for(match.part_type, library.scheme),
            "known": True,
        })

    for gap in tiling.unmatched:
        if gap.part_type not in UNIT_POSITIONS:
            continue
        out.append({
            "name": "",
            "display": f"{gap.length:,} bp not in the library",
            "component": "",
            "part_type": gap.part_type,
            "length": gap.length,
            "start": gap.start,
            "end": gap.end,
            "color": color_for(gap.part_type, library.scheme),
            "known": False,
        })

    out.sort(key=lambda p: p["start"])
    return out


#: Species prefixes the kit's own part names carry: ScTDH3, SpHIS5, KlLEU2.
_SPECIES = re.compile(r"^(Sc|Sp|Kl|Ec|At|Hs|Mm|Ca)(?=[A-Z0-9])")
_ROLE_WORD = re.compile(r"\b(promoter|terminator|prom|term)\b\.?", re.IGNORECASE)


def short_label(component: str, part_type: str) -> str:
    """A part's name in the short form a plasmid map has room for.

    `ScCCW12 Promoter` is sixteen characters that mean `prCCW12`, and at the
    width one fragment of a multigene construct gets, the difference decides
    whether the name is on one line or three. The species prefix goes because
    everything in a yeast kit is Sc unless it says otherwise, and the role goes
    into the affix that already means it.
    """
    text = _ROLE_WORD.sub("", component or "").strip(" ·-_")
    if not text:
        return component or ""
    # only the first thing named: `His3 Promoter · His3` is one part
    text = text.split("·")[0].strip()
    base = _SPECIES.sub("", text).strip()
    if not base:
        return text

    if part_type == "2":
        # `pTDH3` already says promoter; do not end up with `prpTDH3`
        gene = re.sub(r"^p(?=[A-Z])", "", base)
        return gene if gene.lower().startswith("pr") else f"pr{gene}"
    if part_type in ("4", "4b"):
        # likewise `tENo1` is the ENO1 terminator, not a gene called tENo1
        gene = re.sub(r"^t(?=[A-Z])", "", base)
        return gene if gene.lower().endswith("ter") else f"{gene}ter"
    return base


# --------------------------------------------------------------------------- #
# designing the cassettes a multigene construct needs
# --------------------------------------------------------------------------- #

#: The eight positions of a cassette, in circle order.
CIRCLE = ("1", "2", "3", "4", "5", "6", "7", "8")


@dataclass
class UnitSpec:
    """One transcription unit, as the person designing it thinks of it."""

    promoter: str = ""
    cds: str = ""
    terminator: str = ""
    cds_b: str = ""
    """A second half, when the coding sequence is split 3a + 3b."""
    terminator_b: str = ""
    """A second half, when position 4 is split 4a + 4b."""
    name: str = ""

    def slots(self) -> dict[str, str]:
        """The positions this unit fills, and with what."""
        out: dict[str, str] = {}
        if self.promoter:
            out["2"] = self.promoter
        if self.cds_b:
            out["3a"], out["3b"] = self.cds, self.cds_b
        elif self.cds:
            out["3"] = self.cds
        if self.terminator_b:
            out["4a"], out["4b"] = self.terminator, self.terminator_b
        elif self.terminator:
            out["4"] = self.terminator
        return out


@dataclass
class CassettePlan:
    """One eight-part plasmid that has to be built before the multigene step."""

    name: str
    parts: dict[str, str] = field(default_factory=dict)
    left: str = ""
    right: str = ""
    left_label: str = ""
    """`ConLS`, `ConL1` - which connector this unit begins with, by name."""
    right_label: str = ""
    """`ConR1`, `ConRE` - and which it ends with."""
    length: int = 0
    issues: list[Issue] = field(default_factory=list)
    record: object | None = None

    @property
    def ok(self) -> bool:
        return self.record is not None and not [i for i in self.issues if i.level == ERROR]


@dataclass
class DesignReport:
    """What to build, in the order it has to be built."""

    name: str
    backbone: str | None = None
    cassettes: list[CassettePlan] = field(default_factory=list)
    multigene_length: int = 0
    multigene: object | None = None
    issues: list[Issue] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return bool(self.cassettes) and all(c.ok for c in self.cassettes) and not self.errors

    @property
    def errors(self) -> list[Issue]:
        return [i for i in self.issues if i.level == ERROR]


def _connectors(library: Library) -> tuple[dict[str, PlasmidEntry], dict[str, PlasmidEntry]]:
    """Connector parts by the overhang they carry: the type 1s, then the type 5s.

    A part whose multigene site was domesticated away is no use here at all -
    it would make a cassette that could never be cut out again - so
    `level3_ready` filters them before anything is designed around them.
    """
    starts: dict[str, PlasmidEntry] = {}
    ends: dict[str, PlasmidEntry] = {}
    for entry in sorted(library.unique_entries(), key=lambda e: (len(e.name), e.name)):
        if not entry.connector_overhang or entry.level3_ready is False:
            continue
        table = starts if entry.call.part_type == "1" else (
            ends if entry.call.part_type == "5" else None
        )
        if table is not None:
            table.setdefault(entry.connector_overhang, entry)
    return starts, ends


#: `Con5`, `ConR1`, `ConLs` - the number a connector is known by, if it has one.
_CONNECTOR_NUMBER = re.compile(r"Con[LR]?(\d+)", re.IGNORECASE)


def connector_order(overhang: str, index: dict[str, dict[str, list[str]]]) -> tuple:
    """Where a connector falls in the kit's own numbering.

    Sorting the junctions by overhang instead - which is alphabetical DNA -
    is how a three-unit design came out running Con2 -> Con5 -> Con4 -> Con1.
    That assembles perfectly and reads as a mistake, and at the bench it is
    four chances to pick up the wrong tube.
    """
    names = index.get(overhang, {}).get("names", [])
    numbers = [
        int(match.group(1))
        for name in names
        if (match := _CONNECTOR_NUMBER.search(name))
    ]
    return (min(numbers) if numbers else 99, overhang)


def terminal_connectors(library: Library) -> tuple[set[str], set[str]]:
    """The two ends of the connector series: where a chain starts and closes.

    Structural, not a lookup. A connector that exists only as a type 1 part can
    begin a chain and can never continue one, so it is the series start -
    `ConLS` in the kit's naming. One that exists only as a type 5 part can only
    close - `ConRE`. Everything between exists in both directions, which is
    what lets it be a junction.

    This is what makes a full-range destination vector recognisable without
    hard-coding an overhang: it is the one whose ends are these two.
    """
    starts, ends = _connectors(library)
    return set(starts) - set(ends), set(ends) - set(starts)


def connector_label(
    overhang: str,
    side: str,
    index: dict[str, dict[str, list[str]]],
    first: set[str],
    last: set[str],
) -> str:
    """What to write on the tube: `ConL1`, `ConR2`, `ConLS`, `ConRE`.

    An overhang is one junction, and the plasmids either side of it are named
    for which side they sit on - the same `Con1` junction is `ConR1` as the
    right end of one unit and `ConL1` as the left end of the next. A part's own
    annotation only says `Con1`, because the plasmid does not know which of the
    two it is being used as; the position it was put in does.
    """
    if overhang in first:
        suffix = "S"
    elif overhang in last:
        suffix = "E"
    else:
        rank = connector_order(overhang, index)[0]
        suffix = str(rank) if rank != 99 else overhang
    return f"Con{side}{suffix}"


def connector_plan(
    library: Library,
    backbone: PlasmidEntry,
    count: int,
    index: dict[str, dict[str, list[str]]] | None = None,
) -> tuple[list[tuple[PlasmidEntry, PlasmidEntry]], list[Issue]]:
    """Which connector pair each unit needs, so the chain closes on the backbone.

    The chain runs from the backbone's right overhang round to its left, so
    unit *i* starts where unit *i-1* ended. Every junction between them needs a
    connector that exists in both directions - a type 5 to end one unit and a
    type 1 to begin the next - which is a much smaller set than the list of
    connectors, and the reason this can fail with a perfectly good backbone.
    """
    starts, ends = _connectors(library)
    left, right = backbone.cassette_overhangs or ("", "")
    issues: list[Issue] = []

    if not left or not right:
        return [], [Issue(ERROR, "no_backbone_ends",
                          f"{backbone.name} has no connector ends to design against")]
    if right not in starts:
        issues.append(Issue(ERROR, "no_first_connector",
                            f"nothing in the library begins a unit at {right}, where "
                            f"{backbone.name} expects the chain to start"))
    if left not in ends:
        issues.append(Issue(ERROR, "no_last_connector",
                            f"nothing in the library ends a unit at {left}, where "
                            f"{backbone.name} expects the chain to close"))
    if issues:
        return [], issues

    # junctions usable in the middle: a unit has to end there and the next
    # begin there, taken in the order the kit numbers them
    if index is None:
        index = connector_names(library)
    middles = sorted(
        (oh for oh in set(ends) & set(starts) if oh not in (left, right)),
        key=lambda oh: connector_order(oh, index),
    )
    if len(middles) < count - 1:
        return [], [Issue(
            ERROR, "not_enough_connectors",
            f"{count} units need {count - 1} junctions between them and only "
            f"{len(middles)} connector(s) in the library work in both directions",
        )]

    chain = [right, *middles[: count - 1], left]
    return [(starts[chain[i]], ends[chain[i + 1]]) for i in range(count)], []


#: Circle order including the split positions, so a design can mix them.
SLOT_ORDER = ("1", "2", "3", "3a", "3b", "4", "4a", "4b", "5", "6", "7", "8", "8a", "8b")


def _released(record, enzyme, name: str, component: str = "") -> Piece | None:
    """The site-free fragment `enzyme` cuts out of `record`, as a Piece.

    The library has this for plasmids it has indexed; a cassette designed a
    moment ago is not one of those, so it is done here from the record itself.
    """
    sequence = seqio.sequence(record)
    fragments = [f for f in digest(sequence, enzyme) if not f.has_sites]
    if len(fragments) != 1:
        return None
    fragment = fragments[0]
    return Piece(
        record=seqio.slice_record(record, fragment.start, fragment.end, name=name),
        left_overhang=fragment.left_overhang or "",
        right_overhang=fragment.right_overhang or "",
        source_name=name,
        part_type=None,
        component=component,
    )


def _default_part(library: Library, part_type: str) -> PlasmidEntry | None:
    """The simplest usable part at a position, for a slot nobody chose."""
    options = [
        e for e in library.parts_of_type(part_type)
        if not e.sites.part_enzyme_internal and not e.internal_multigene_positions
    ]
    return min(options, key=lambda e: (not e.component, len(e.name), e.name)) if options else None


def design(
    library: Library,
    units: Sequence[UnitSpec],
    backbone: str | None = None,
    shared: dict[str, str] | None = None,
    name: str = "pMultigene",
) -> DesignReport:
    """Work out the eight-part plasmids needed to build a multigene construct.

    The design direction. Level 2 asks "what does this cassette come out as";
    this asks the question a project actually starts from - *I want these genes
    expressed, in this order* - and answers it with the list of plasmids to
    build first, each one a complete eight-part assembly with its connectors
    already worked out so the cassettes chain in the order given.

    Everything not named is chosen: connectors from the backbone's ends, and
    the marker, origin and E. coli backbone from the simplest part at each
    position. Nothing is silently substituted at a position the caller filled.
    """
    report = DesignReport(name=name)
    if not units:
        report.issues.append(Issue(ERROR, "no_units", "a design needs at least one unit"))
        return report

    naming = connector_names(library)
    first_end, last_end = terminal_connectors(library)
    vectors = [v for v in library.multigene_vectors() if v.cassette_overhangs]
    if backbone:
        chosen = library.get(backbone)
        if chosen is None:
            report.issues.append(Issue(ERROR, "no_backbone", f"no plasmid named {backbone}"))
            return report
    else:
        # the first backbone the library can actually close a chain of this
        # length against, rather than the first one alphabetically
        # Prefer a full-range destination: one running from the series start
        # to the series end, so the units come out ConLS -> ConR1 / ConL1 ->
        # ConR2 ... -> ConRE, which is the order the kit is numbered in and
        # the order a bench protocol is written in. A short-range vector
        # assembles just as well and hands you Con2 -> Con3 -> ConR1.
        def canonical(vector: PlasmidEntry) -> tuple:
            left, right = vector.cassette_overhangs or ("", "")
            return (
                not (right in first_end and left in last_end),
                len(vector.name),
                vector.name,
            )

        chosen = next(
            (v for v in sorted(vectors, key=canonical)
             if not connector_plan(library, v, len(units), naming)[1]),
            None,
        )
        if chosen is None:
            report.issues.append(Issue(
                ERROR, "no_usable_backbone",
                f"no destination vector in the library can close a chain of "
                f"{len(units)} unit(s) with the connectors it has",
            ))
            return report
    report.backbone = chosen.name

    pairs, issues = connector_plan(library, chosen, len(units), naming)
    report.issues.extend(issues)
    if issues:
        # naming a backbone that cannot work is only half an answer; the other
        # half is which ones can, and the caller has no way to work that out
        workable = [
            v.name for v in vectors
            if v.name != chosen.name
            and not connector_plan(library, v, len(units), naming)[1]
        ]
        if workable:
            report.issues.append(Issue(
                INFO, "try_instead",
                f"{len(workable)} other destination vector(s) can close a chain of "
                f"{len(units)}: {', '.join(workable[:6])}"
                f"{'\u2026' if len(workable) > 6 else ''}",
            ))
        else:
            report.issues.append(Issue(
                INFO, "none_workable",
                f"no destination vector in the library can close a chain of "
                f"{len(units)} with the connector parts it has",
            ))
        return report

    fixed = dict(shared or {})
    for position in ("6", "7", "8"):
        if not fixed.get(position):
            found = _default_part(library, position)
            if found is None:
                report.issues.append(Issue(
                    ERROR, "no_part", f"the library has no usable type {position} part"
                ))
                return report
            fixed[position] = found.name

    unit_pieces: list[Piece] = []
    for index, (spec, (con_left, con_right)) in enumerate(zip(units, pairs), start=1):
        left_overhang = con_left.connector_overhang or ""
        right_overhang = con_right.connector_overhang or ""
        plan = CassettePlan(
            name=spec.name or f"{name}_TU{index}",
            left=left_overhang,
            right=right_overhang,
            left_label=connector_label(left_overhang, "L", naming, first_end, last_end),
            right_label=connector_label(right_overhang, "R", naming, first_end, last_end),
        )
        plan.parts = {"1": con_left.name, **spec.slots(), "5": con_right.name, **fixed}

        pieces: list[Piece] = []
        for slot in SLOT_ORDER:
            plasmid = plan.parts.get(slot)
            if not plasmid:
                continue
            entry = library.get(plasmid)
            if entry is None:
                plan.issues.append(Issue(ERROR, "no_part", f"no plasmid named {plasmid}"))
                continue
            try:
                pieces.append(library.piece(entry))
            except (ValueError, OSError) as error:
                plan.issues.append(Issue(ERROR, "unreadable_part", f"{plasmid}: {error}"))

        if not [i for i in plan.issues if i.level == ERROR]:
            built = assemble(pieces, library.scheme.part_enzyme, name=plan.name)
            plan.issues.extend(built.issues)
            plan.length = built.length
            plan.record = built.product
            if built.product is not None:
                unit = _released(
                    built.product, library.scheme.multigene_enzyme, plan.name,
                    component=spec.name,
                )
                if unit is None:
                    plan.issues.append(Issue(
                        ERROR, "no_unit",
                        f"{plan.name} assembles but "
                        f"{library.scheme.multigene_enzyme.name} does not release a "
                        f"single site-free fragment from it",
                    ))
                else:
                    unit_pieces.append(unit)
        report.cassettes.append(plan)

    if len(unit_pieces) == len(units):
        backbone_piece = library.piece(chosen, library.scheme.multigene_enzyme)
        multigene = assemble(
            [backbone_piece, *unit_pieces], library.scheme.multigene_enzyme, name=name
        )
        report.issues.extend(multigene.issues)
        report.multigene_length = multigene.length
        report.multigene = multigene.product
    return report


def design_report(report: DesignReport, library: Library) -> str:
    """The design as a build order: what to make first, and what it becomes.

    A design that did not work gets the reasons instead. Numbering the steps of
    something that cannot be built - "0 cassettes to build, then one multigene
    assembly" - reads as a result and is not one.
    """
    if not report.cassettes:
        lines = [
            f"{report.name} - cannot be designed",
            "=" * (len(report.name) + 22),
            "",
        ]
        if report.backbone:
            lines.append(f"Against destination backbone: {report.backbone}")
            lines.append("")
        for issue in report.issues:
            lines.append(f"  {'!' if issue.level == ERROR else '-'} {issue.message}")
        return "\n".join(lines) + "\n"

    lines = [
        f"{report.name} - build order",
        "=" * (len(report.name) + 14),
        "",
        f"{len(report.cassettes)} cassette(s) to build, then one multigene assembly.",
        f"Destination backbone: {report.backbone}",
        "",
    ]
    for index, plan in enumerate(report.cassettes, start=1):
        lines += [
            f"{index}. {plan.name}   {plan.length:,} bp   "
            f"connectors {plan.left_label} ({plan.left}) → "
            f"{plan.right_label} ({plan.right})",
            f"   {library.scheme.part_enzyme.name} assembly of:",
        ]
        for slot in SLOT_ORDER:
            if slot not in plan.parts:
                continue
            entry = library.get(plan.parts[slot])
            # a connector's own label says `Con1`; which side of the junction
            # it is being used as comes from the position it sits in
            if slot == "1":
                what = f" [{plan.left_label}]"
            elif slot == "5":
                what = f" [{plan.right_label}]"
            else:
                what = f" [{entry.component}]" if entry and entry.component else ""
            lines.append(f"      type {slot:<3} {plan.parts[slot]}{what}")
        for issue in plan.issues:
            if issue.level != INFO:
                lines.append(f"      ! {issue.level}: {issue.message}")
        lines.append("")

    lines += [
        f"{len(report.cassettes) + 1}. {report.name}   "
        f"{report.multigene_length:,} bp",
        f"   {library.scheme.multigene_enzyme.name} assembly of {report.backbone} "
        f"and the {len(report.cassettes)} cassette(s) above, in that order.",
        "",
    ]
    problems = [i for i in report.issues if i.level != INFO]
    if problems:
        lines += ["Before you start", "-" * 16]
        lines += [f"  ! {i.level}: {i.message}" for i in problems]
        lines.append("")
    return "\n".join(lines)
