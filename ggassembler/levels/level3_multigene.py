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

from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqUtils import MeltingTemp

from ..core.assembly import ERROR, INFO, WARNING, AssemblyResult, Issue, Piece, assemble
from ..core.library import Library, PlasmidEntry
from ..core.parttypes import YTK, Scheme

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

    for backbone in sorted(library.multigene_vectors(), key=_unit_preference):
        left, right = backbone.cassette_overhangs or ("", "")
        if not left or not right:
            continue
        chain = _walk(starting, start=right, stop=left)
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
        for entry in starting.get(at, []):
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
