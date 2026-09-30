"""Plasmid folder scan, index, mtime cache and manual type overrides.

What a plasmid *is* gets decided by simulating the digest, not by trusting its
name or its annotations, because in a real library both are often stale. The
enzymes settle it:

* the scheme's **part enzyme** (BsaI in YTK) releasing one fragment whose
  overhangs sit on the part circle means a type 1-8 part plasmid, ready for an
  eight-part assembly;
* the scheme's **multigene enzyme** (BsmBI) releasing a fragment with connector
  ends means a finished cassette, ready for a multi-TU assembly - or, with a
  dropout between those ends, the destination vector for one.

When the digest and an annotation disagree, the digest wins and the
disagreement is recorded as a conflict, never silently dropped. Every call
carries the evidence that produced it.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
from dataclasses import asdict, dataclass, field, replace
from pathlib import Path
from typing import Any, Sequence

from Bio.SeqRecord import SeqRecord

from . import parttypes, seqio
from .enzymes import Enzyme, Fragment, Site, digest, find_sites
from .assembly import Piece
from .parttypes import YTK, Scheme, color_for, span

CACHE_DIR = ".ggasm"
INDEX_FILE = "index.json"
OVERRIDES_FILE = "overrides.json"
DESIGNS_FILE = "designs.json"
BUILDS_FILE = "builds.jsonl"
CONFIG_FILE = "config.json"
INDEX_VERSION = 11

#: Confidence levels, best first.
DIGEST, ANNOTATION, FILENAME, NONE = "digest", "annotation", "filename", "none"

#: Roles a plasmid can hold, all decided by which enzyme cuts where.
PART, ENTRY_VECTOR, CONNECTOR, CASSETTE, MULTIGENE_VECTOR = (
    "part", "entry_vector", "connector", "cassette", "multigene_vector",
)

#: A folder name that means "kept, but not in play". Redundant copies of a
#: sequence go here rather than being deleted: the file is still on disk and
#: still in git history, it just stops appearing in every dropdown.
ARCHIVE_DIR = "_archive"

SKIP_DIRS = {
    CACHE_DIR, ARCHIVE_DIR, ".git", ".venv", "venv", "__pycache__",
    "node_modules", ".Rproj.user",
}

#: Names that say a sequence was checked by sequencing. Between two identical
#: files that is the one worth keeping, whatever it is called - the shortest
#: name is a tie-break, not a quality signal.
_SEQ_VERIFIED = re.compile(r"seq[ _-]*verified|sequence[ _-]*verified|verified", re.IGNORECASE)

#: Labels that describe where a part came from or how it was made, rather than
#: what it is. Everything else is fair game: for a type 8 part the E. coli
#: backbone *is* the component, so nothing is filtered on content grounds - the
#: fragment's own boundaries already exclude whatever is not part of it.
#: ``scar`` is matched with explicit boundaries rather than ``\b`` so that it
#: still catches ``scar_1`` and ``BsaI scar`` without eating ``mScarlet-I``.
_LABEL_NOISE = re.compile(
    r"^pytk\d+$|(?<![A-Za-z])scars?(?![A-Za-z])|\b(fp|rp)$|^con[lrs]?[0-9]*\s*(fp|rp)$|"
    r"^(promoter|terminator|cds|gene|insert|fragment)$|"
    r"benchling|amplicon|_bsai_|translation|barcode|"
    r"^(bsai|bsmbi|esp3i|noti|bbsi|sapi)$",
    re.IGNORECASE,
)

#: What a cloning pipeline appends to a fragment's name: the operation that
#: made it, not what it is.
_BOILERPLATE_TOKEN = re.compile(
    r"[ _\-]+(pcr|amplicon|bsa[i1]|bsmb[i1]|bbs[i1]|esp3i|not[i1]|sap[i1]|"
    r"largest|smallest|fragment|frag|insert|product|digest(?:ed)?|"
    r"seq|sequenced|verified|copy|clone|final|v\d+)$",
    re.IGNORECASE,
)

#: A library index code at the front of a file name - ``L114_``, ``C8_``,
#: ``GG212_``. Shaped exactly like a gene name (``TDH3_``, ``Cas9_``, ``ESA1_``),
#: so `_is_index_code` decides between them on the balance of letters to digits.
_NAME_PREFIX = re.compile(r"^([A-Za-z]{0,3})(\d{1,4})[_\-]")

#: What a file name says about the format rather than the contents.
_FORMAT_WORD = (
    r"gene|part|plasmid|vector|backbone|construct|yeast|"
    r"type[ _\-]?[1-8][ab]?|[1-8][ab]|t[1-8]|"
    r"moclo|ytk|seq|sequenced|verified|final|v\d+"
)
_NAME_SUFFIX = re.compile(rf"[ _\-]+(?:{_FORMAT_WORD})$", re.IGNORECASE)
_ONLY_FORMAT = re.compile(rf"^(?:{_FORMAT_WORD})$", re.IGNORECASE)


def _is_index_code(letters: str, digits: str) -> bool:
    """Whether ``L114``/``C8`` is a catalogue number rather than a gene name.

    A catalogue number is mostly digits - one or two letters then a serial.
    A gene name is mostly letters with one on the end: ``TDH3``, ``PHO5``,
    ``Cas9``. Getting this backwards renames ``TDH3_gene`` to ``gene``.
    """
    return len(letters) <= 2 or len(digits) >= len(letters)


def clean_label(label: str) -> str:
    """Strip the boilerplate a cloning pipeline appends, keeping the name.

    Benchling and tools like it name a fragment after the operation that made
    it: ``XYL2_amplicon_BsaI_largest`` is the largest BsaI fragment of the XYL2
    amplicon. ``XYL2`` is the part worth keeping, and discarding the whole
    label because it mentions BsaI throws away the only description the file
    has - which is why so many perfectly well annotated plasmids read as blank.
    """
    text = label.strip()
    while True:
        stripped = _BOILERPLATE_TOKEN.sub("", text).strip()
        if stripped == text or not stripped:
            break
        text = stripped
    return text or label.strip()


def name_component(name: str) -> str:
    """A description read off the file name, for a fragment with no labels.

    Some files carry nothing but Benchling's automatic translation, so the only
    statement of what is in them is what they are called: ``L114_ERG10_yeast``
    is the ERG10 gene. This is weaker evidence than an annotation and is
    recorded as such (`PlasmidEntry.component_source`), never used to decide a
    type - it only fills in a description that would otherwise be blank.

    Nothing is returned when the result would just repeat the file name, since
    the name is already on screen next to it.
    """
    def trim(text: str) -> str:
        while True:
            stripped = _NAME_SUFFIX.sub("", text).strip()
            if stripped == text or not stripped:
                break
            text = stripped
        return "" if _ONLY_FORMAT.match(text) else text

    stem = name.strip()
    prefix = _NAME_PREFIX.match(stem)
    without = stem[prefix.end():] if prefix and _is_index_code(*prefix.groups()) else stem
    body = trim(without) or trim(stem)
    body = body.replace("_", " ").strip()
    squash = lambda s: re.sub(r"[^a-z0-9]+", "", s.lower())  # noqa: E731
    if len(body) < 2 or squash(body) == squash(name):
        return ""
    return body

#: Feature types that annotate the work rather than the plasmid, never used.
_NOT_A_COMPONENT = {"primer", "primer_bind", "source"}

#: Binding sites: part of the molecule, but usually incidental to what a part
#: *is*, so they describe a fragment only when nothing else does. For a dropout
#: whose whole content is an I-SceI site, that site is the honest description.
_WEAK_COMPONENT = {"misc_binding", "protein_bind"}

_TYPE_IN_TEXT = re.compile(r"(?<![A-Za-z])type[\s_-]*([1-8][ab]?)(?![A-Za-z0-9])", re.IGNORECASE)
_YTK_NUMBER = re.compile(r"\bpYTK(\d+)", re.IGNORECASE)
_ECOLI_MARKERS = (
    ("AmpR", r"ampr|ampicillin|\bbla\b|beta-?lactamase"),
    ("KanR", r"kanr|kanamycin|neomycin|\bnptii\b|\baph\b"),
    ("CmR", r"\bcmr\b|chloramphenicol|\bcat\b"),
    ("SpecR", r"specr|spectinomycin|\baada\b"),
    ("TetR", r"\btetr\b|tetracycline"),
)


# --------------------------------------------------------------------------- #
# data model
# --------------------------------------------------------------------------- #


@dataclass
class FeatureInfo:
    type: str
    label: str
    start: int
    end: int
    strand: int | None


@dataclass
class TypeCall:
    """What the plasmid is, how sure we are, and why."""

    part_type: str | None = None
    confidence: str = NONE
    reason: str = "no evidence"
    source: str = "detected"
    """``detected`` or ``manual``."""
    five_prime: str | None = None
    three_prime: str | None = None
    reversed_sites: bool = False
    conflict: str | None = None
    assumed_circular: bool = False
    """The file said linear, but the digest says otherwise and we believed it."""

    @property
    def display_type(self) -> str:
        return self.part_type or "?"


@dataclass
class Destination:
    """A plasmid that can receive a new part, and what it expects of one.

    Two shapes, both common:

    * a **dropout vector** - ccdB or GFP between two outward-facing sites. The
      cloning enzyme's sites leave on the dropout, and the backbone that stays
      behind already carries the part enzyme's sites, so the insert needs only
      the part type's own overhangs.
    * an **existing part plasmid** - cut it with the part enzyme and its current
      part falls out, leaving the same backbone ready for a replacement of the
      same type.

    The one case that needs more of the insert is the universal entry vector
    (``pYTK001``), whose backbone has no part-enzyme sites: there the insert has
    to bring its own, which is what `supplies_part_sites` records.
    """

    name: str
    path: str
    cloning_enzyme: str
    accepts: tuple[str, str]
    """The overhangs an insert must present, 5' then 3'."""
    part_type: str | None
    """The part type those overhangs make, when they make a known one."""
    supplies_part_sites: bool
    """True when the backbone already carries the part enzyme's sites."""
    kind: str
    """``dropout`` or ``part_plasmid``."""
    backbone_length: int
    marker: str | None = None
    component: str = ""
    aliases: list[str] = field(default_factory=list)

    @property
    def display(self) -> str:
        what = f"type {self.part_type}" if self.part_type else "/".join(self.accepts)
        return f"{self.name} [{self.component}]" if self.component else f"{self.name} ({what})"


@dataclass
class SiteCounts:
    part_enzyme_total: int = 0
    part_enzyme_internal: int = 0
    multigene_enzyme_total: int = 0
    multigene_enzyme_internal: int = 0
    linearizer_total: int = 0


@dataclass
class PlasmidEntry:
    path: str
    name: str
    """Display name: the file stem, which is unique and is what the user sees."""
    record_id: str
    """The LOCUS/id inside the file, often blank, duplicated or '.'."""
    length: int
    topology: str
    call: TypeCall
    sites: SiteCounts
    scheme: str = YTK.name
    roles: list[str] = field(default_factory=list)
    connector_overhang: str | None = None
    """For a type 1 or 5 part: the connector-specific multigene-enzyme overhang,
    read from the sequence. ConLX and ConRX with the same X share it."""
    cassette_overhangs: tuple[str, str] | None = None
    """For a cassette or multigene vector: (ConL overhang, ConR overhang)."""
    cassette_span: tuple[int, int] | None = None
    """Where that released fragment sits in the plasmid, so callers can tell the
    transcription unit's own features from the backbone's without re-digesting."""
    ecoli_marker: str | None = None
    checksum: str = ""
    part_span: tuple[int, int] | None = None
    """Where the fragment this plasmid contributes sits, for the part enzyme."""
    component: str = ""
    """What the contributed fragment holds, read from its own annotations."""
    component_source: str = ""
    """``annotation`` when the description came from the fragment's own
    features, ``filename`` when the file had none and the name was read
    instead, ``""`` when there is no description at all. Never affects the type
    call - that is always the digest."""
    internal_multigene_positions: list[int] = field(default_factory=list)
    """Where the multigene enzyme cuts *inside* the fragment this plasmid has to
    contribute intact. Any of these would cut a cassette apart at Level 3."""
    multigene_overhangs: list[str] = field(default_factory=list)
    """Every overhang the multigene enzyme would leave, in order round the
    plasmid. Worth showing even when no connector pair could be made out."""
    aliases: list[str] = field(default_factory=list)
    """Other filenames holding this exact sequence, merged into this entry."""
    features: list[FeatureInfo] = field(default_factory=list)
    level3_ready: bool | None = None
    """For a type 1 or 5 part: whether it still carries the multigene enzyme's
    site, and so whether a cassette built with it could ever be released for a
    multigene assembly. `None` where the question does not arise - every other
    position contributes nothing to the cassette's ends."""
    conc_ng_ul: float | None = None
    """What a prep of this plasmid measured at. A property of the tube, not of
    the sequence, so it is kept in the overrides file rather than read from the
    record - and it is what makes an equimolar reaction setup possible."""

    @property
    def display(self) -> str:
        """``pYTK009 [ScTDH3 Promoter]`` - the name, then what is in it."""
        names = " / ".join([self.name, *self.aliases])
        return f"{names} [{self.component}]" if self.component else names

    @property
    def is_part(self) -> bool:
        return bool(self.call.part_type) and not self.call.reversed_sites

    @property
    def is_entry_vector(self) -> bool:
        return ENTRY_VECTOR in self.roles

    @property
    def is_cassette(self) -> bool:
        return CASSETTE in self.roles

    @property
    def is_multigene_vector(self) -> bool:
        return MULTIGENE_VECTOR in self.roles

    @property
    def is_assembled(self) -> bool:
        """True for anything the multigene enzyme reads as a finished construct."""
        return bool({CASSETTE, MULTIGENE_VECTOR, ENTRY_VECTOR} & set(self.roles))

    @property
    def overhangs(self) -> tuple[str | None, str | None]:
        return self.call.five_prime, self.call.three_prime

    def to_dict(self) -> dict[str, Any]:
        data = asdict(self)
        if self.cassette_overhangs:
            data["cassette_overhangs"] = list(self.cassette_overhangs)
        if self.cassette_span:
            data["cassette_span"] = list(self.cassette_span)
        if self.part_span:
            data["part_span"] = list(self.part_span)
        return data

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> PlasmidEntry:
        data = dict(data)
        call = TypeCall(**data.pop("call"))
        sites = SiteCounts(**data.pop("sites"))
        features = [FeatureInfo(**f) for f in data.pop("features", [])]
        overhangs = data.pop("cassette_overhangs", None)
        span_ = data.pop("cassette_span", None)
        part_span_ = data.pop("part_span", None)
        return cls(
            call=call,
            sites=sites,
            features=features,
            cassette_overhangs=tuple(overhangs) if overhangs else None,
            cassette_span=tuple(span_) if span_ else None,
            part_span=tuple(part_span_) if part_span_ else None,
            **data,
        )


# --------------------------------------------------------------------------- #
# detection
# --------------------------------------------------------------------------- #


@dataclass
class _Candidate:
    five: str
    three: str
    call: parttypes.SpanCall
    reversed_sites: bool
    extra_sites: int
    spans_origin: bool = False
    """True when one of the two flanking sites only exists because the ends were
    joined - the evidence a genuinely linear molecule must not be judged on."""


def _enclosed(site: Site, start: int, end: int, n: int) -> bool:
    offset = (site.start - start) % n
    return offset + (site.end - site.start) <= seqio.span_length(start, end, n)


def _candidates(seq: str, enzyme: Enzyme, scheme: Scheme) -> list[_Candidate]:
    """Every pair of cuts whose overhangs delimit a legal arc of the part circle."""
    n = len(seq)
    sites = find_sites(seq, enzyme)
    out: list[_Candidate] = []
    for left in sites:
        for right in sites:
            if left is right or left.top_cut == right.top_cut:
                continue
            call = span(left.overhang, right.overhang, scheme=scheme)
            if call is None or not call.legal or call.is_full_circle:
                continue
            inside = [s for s in sites if _enclosed(s, left.top_cut, right.top_cut, n)]
            flanks_inside = left in inside and right in inside
            out.append(
                _Candidate(
                    five=left.overhang,
                    three=right.overhang,
                    call=call,
                    reversed_sites=flanks_inside,
                    extra_sites=len(inside) - (2 if flanks_inside else 0),
                    spans_origin=left.end > n or right.end > n,
                )
            )
    # a part released site-free beats a dropout that carries its own sites;
    # among equals, the arc with the fewest stray sites in it wins.
    out.sort(key=lambda c: (c.reversed_sites, c.extra_sites, len(c.call.atoms)))
    return out


def detect_type(record: SeqRecord, path: Path | None = None, scheme: Scheme = YTK) -> TypeCall:
    """Classify one plasmid by digest, recording the evidence behind the call."""
    seq = seqio.sequence(record)
    circular = seqio.is_circular(record)
    enzyme = scheme.part_enzyme

    # A GenBank file can say "linear" and still be a plasmid - plenty of tools
    # write that by default. Rather than take its word and give up, try the
    # circular reading too, and believe it only when the evidence does not
    # depend on the ends being joined.
    assumed_circular = False
    candidates = _candidates(seq, enzyme, scheme) if seq else []
    if not circular:
        candidates = [c for c in candidates if not c.spans_origin]
        assumed_circular = bool(candidates)

    digest_call: TypeCall | None = None
    if True:
        if candidates:
            best = candidates[0]
            name = best.call.name + ("r" if best.reversed_sites else "")
            how = (
                f"{enzyme.name} sites reversed: fragment carries its own sites (dropout)"
                if best.reversed_sites
                else f"{enzyme.name} releases a site-free fragment"
            )
            extra = (
                f"; {best.extra_sites} further {enzyme.name} site(s) inside"
                if best.extra_sites
                else ""
            )
            digest_call = TypeCall(
                part_type=name,
                confidence=DIGEST,
                reason=f"{how}: {best.five} -> {best.three}{extra}",
                five_prime=best.five,
                three_prime=best.three,
                reversed_sites=best.reversed_sites,
                assumed_circular=assumed_circular,
                conflict=(
                    f"the file says this record is linear, but {enzyme.name} cuts it like a "
                    f"circular {name} part - read as circular"
                    if assumed_circular
                    else None
                ),
            )

    annotated = _type_from_annotation(record)
    named = _type_from_filename(path, record)

    if digest_call:
        for other, where in ((annotated, "annotation"), (named, "filename")):
            if digest_call.conflict:
                break
            if other and other[0] and _base_type(other[0]) != _base_type(digest_call.part_type):
                digest_call.conflict = (
                    f"{where} says type {other[0]} ({other[1]}), digest says "
                    f"{digest_call.part_type} - keeping the digest"
                )
                break
        return digest_call

    # Releasing a part takes two cuts, so a plasmid neither enzyme cuts cannot
    # contribute to any reaction, whatever it is called. That matters most for
    # a finished multigene construct: every site was consumed building it, but
    # it still carries the labels of the parts that went in - `L6_CYC1_type4`
    # on a 10.8 kb assembly - and believing one of those turns the whole
    # construct into a terminator.
    #
    # This is not the digest being silent, which is what lets an annotation
    # speak at all below. Finding no sites is positive evidence that nothing
    # can be cut out, so it is a disagreement, and the digest wins - with the
    # claim recorded as a conflict rather than dropped.
    #
    # A plasmid the *multigene* enzyme still cuts is left alone here: it is
    # probably a cassette, and `_reconcile` says so far more precisely once the
    # roles are known.
    claim = annotated or named
    uncuttable = (
        len(find_sites(seq, enzyme, circular)) < 2
        and len(find_sites(seq, scheme.multigene_enzyme, circular)) < 2
    )
    if claim and claim[0] and uncuttable:
        where = "annotation" if annotated else "file name"
        return TypeCall(
            part_type=None,
            confidence=NONE,
            reason=(
                f"neither {enzyme.name} nor {scheme.multigene_enzyme.name} cuts this "
                f"plasmid: nothing can be cut out of it"
            ),
            conflict=(
                f"{where} says type {claim[0]} ({claim[1]}), but no {enzyme.name} site is "
                f"left to cut it out - this looks like a finished construct, not a part"
            ),
        )

    if annotated:
        return TypeCall(part_type=annotated[0], confidence=ANNOTATION, reason=annotated[1])
    if named:
        return TypeCall(part_type=named[0], confidence=FILENAME, reason=named[1])

    reason = (
        f"no {enzyme.name} pair delimits a {scheme.name} part"
        if circular
        else "record is linear"
    )
    return TypeCall(part_type=None, confidence=NONE, reason=reason)


def _base_type(name: str | None) -> str | None:
    return name.rstrip("r") if name else None


def _type_from_annotation(record: SeqRecord) -> tuple[str | None, str] | None:
    """Look for a type named in a /note, /label or the source feature."""
    for feature in record.features:
        for key in ("label", "note", "standard_name", "product"):
            for value in feature.qualifiers.get(key, []):
                match = _TYPE_IN_TEXT.search(str(value))
                if match:
                    return match.group(1), f'{feature.type} /{key}="{value}"'
    comment = str(record.annotations.get("comment", ""))
    match = _TYPE_IN_TEXT.search(comment)
    if match:
        return match.group(1), f'comment "{comment.strip()[:60]}"'
    return None


def _type_from_filename(path: Path | None, record: SeqRecord) -> tuple[str | None, str] | None:
    names = [n for n in (path.name if path else None, record.id, record.name) if n]
    for name in names:
        match = _TYPE_IN_TEXT.search(name)
        if match:
            return match.group(1), f'name "{name}"'
    for name in names:
        match = _YTK_NUMBER.search(name)
        if match:
            return None, f'name "{name}" looks like a YTK plasmid, but names no type'
    return None


@dataclass
class _Roles:
    roles: list[str] = field(default_factory=list)
    connector_overhang: str | None = None
    cassette_overhangs: tuple[str, str] | None = None
    cassette_span: tuple[int, int] | None = None


def connector_pair(
    record: SeqRecord, vocabulary: set[str], scheme: Scheme
) -> tuple[tuple[str, str], tuple[int, int]] | None:
    """Find a transcription unit flanked by two known connector overhangs.

    An assembled cassette that still carries an internal site for the multigene
    enzyme does not cut into two clean pieces, so looking for exactly two
    fragments misses it. Instead, look for a *pair of cuts whose overhangs are
    connectors this library knows*, and take the arc between them.

    Which of the two arcs is the unit follows from the geometry, as everywhere
    else here: the flanking sites point inward, so neither of them lies within
    the unit they release. The E. coli marker breaks any remaining tie, since
    the cassette's backbone carries it and the unit does not.
    """
    seq = seqio.sequence(record)
    length = len(seq)
    sites = find_sites(seq, scheme.multigene_enzyme)
    cuts = sorted({s.top_cut: s for s in sites}.items())
    known = [(cut, site) for cut, site in cuts if site.overhang in vocabulary]
    if len(known) < 2:
        return None

    marker_at = _marker_position(record)
    best = None
    for left_cut, left in known:
        for right_cut, right in known:
            if left is right or left.overhang == right.overhang:
                continue
            width = (right_cut - left_cut) % length or length
            if width < 100 or width > length - 100:
                continue
            if _enclosed(left, left_cut, right_cut, length) or _enclosed(
                right, left_cut, right_cut, length
            ):
                continue  # a flanking site inside the arc means we have it backwards
            carries_marker = (
                marker_at is not None and ((marker_at - left_cut) % length) < width
            )
            rank = (carries_marker, width)
            if best is None or rank < best[0]:
                best = (rank, (left.overhang, right.overhang), (left_cut, right_cut))

    if best is None:
        return None
    return best[1], best[2]


def _marker_position(record: SeqRecord) -> int | None:
    """Where the E. coli selection marker sits, if it is annotated."""
    for feature in record.features:
        for values in feature.qualifiers.values():
            for value in values:
                for _, pattern in _ECOLI_MARKERS:
                    if re.search(pattern, str(value), re.IGNORECASE):
                        return int(feature.location.start)
    return None


def analyse_roles(record: SeqRecord, call: TypeCall, scheme: Scheme = YTK) -> _Roles:
    """Entry vector, connector and cassette detection, all via the multigene enzyme."""
    seq = seqio.sequence(record)
    info = _Roles()
    if not (seqio.is_circular(record) or call.assumed_circular):
        return info

    enzyme = scheme.multigene_enzyme
    fragments = digest(seq, enzyme)
    site_free = [f for f in fragments if not f.has_sites and f.left_overhang]

    for fragment in site_free:
        if fragment.overhangs == tuple(reversed(scheme.entry_vector_overhangs)):
            info.roles.append(ENTRY_VECTOR)
            break

    if ENTRY_VECTOR not in info.roles and len(fragments) == 2 and site_free:
        released = site_free[0]
        five, three = released.overhangs
        untouched_by_part_enzyme = span(five, three, scheme=scheme) is None
        if untouched_by_part_enzyme and five not in scheme.entry_vector_overhangs:
            # A plasmid with connector ends *and* a part-enzyme dropout between
            # them is a destination vector for the next level up, not a cassette
            # to put into one: it has nothing to contribute until its dropout is
            # replaced. Level 3 offers these as backbones and nothing else.
            info.roles.append(MULTIGENE_VECTOR if call.reversed_sites else CASSETTE)
            info.cassette_overhangs = (five, three)
            info.cassette_span = (released.start, released.end)

    if _base_type(call.part_type) in {"1", "5"} and call.five_prime:
        info.connector_overhang = _connector_overhang(seq, call, scheme)
        if info.connector_overhang:
            info.roles.append(CONNECTOR)

    if call.part_type and not call.reversed_sites:
        info.roles.insert(0, PART)

    return info


def _connector_overhang(seq: str, call: TypeCall, scheme: Scheme) -> str | None:
    """The multigene-enzyme overhang carried inside a type 1 or 5 part."""
    part = part_fragment(seq, call, scheme)
    if part is None:
        return None
    sites = [
        s
        for s in find_sites(seq, scheme.multigene_enzyme)
        if _enclosed(s, part.start, part.end, len(seq))
    ]
    return sites[0].overhang if len(sites) == 1 else None


def part_fragment(seq: str, call: TypeCall, scheme: Scheme = YTK) -> Fragment | None:
    """The fragment a part plasmid contributes to an assembly."""
    if not (call.five_prime and call.three_prime):
        return None
    for fragment in digest(seq, scheme.part_enzyme):
        if fragment.overhangs == (call.five_prime, call.three_prime):
            return fragment
    return None


def count_sites_for(
    record: SeqRecord, call: TypeCall, roles: _Roles, scheme: Scheme = YTK
) -> SiteCounts:
    """Site totals, and how many are *unexpected* - the ones that break a reaction."""
    seq = seqio.sequence(record)
    circular = seqio.is_circular(record) or call.assumed_circular
    part_total = len(find_sites(seq, scheme.part_enzyme, circular))
    multi_total = len(find_sites(seq, scheme.multigene_enzyme, circular))
    noti_total = len(find_sites(seq, scheme.linearizer, circular))

    expected_part = 2 if call.confidence == DIGEST else 0
    expected_multi = 0
    if {ENTRY_VECTOR, CASSETTE, MULTIGENE_VECTOR} & set(roles.roles):
        expected_multi += 2
    if roles.connector_overhang:
        expected_multi += 1

    return SiteCounts(
        part_enzyme_total=part_total,
        part_enzyme_internal=max(0, part_total - expected_part),
        multigene_enzyme_total=multi_total,
        multigene_enzyme_internal=max(0, multi_total - expected_multi),
        linearizer_total=noti_total,
    )


def _reconcile(call: TypeCall, roles: _Roles) -> TypeCall:
    """Refuse a weak part call for something the enzymes show is already assembled.

    An assembled cassette is full of part-derived features, so a ``/label`` like
    "type 4 terminator" inside it would otherwise type the whole 6 kb plasmid as
    a terminator part. A digest-backed call is never touched.
    """
    assembled = {CASSETTE, MULTIGENE_VECTOR, ENTRY_VECTOR} & set(roles.roles)
    if call.confidence in (DIGEST, NONE) or not assembled:
        return call
    role = next(r for r in (CASSETTE, MULTIGENE_VECTOR, ENTRY_VECTOR) if r in roles.roles)
    return TypeCall(
        part_type=None,
        confidence=NONE,
        reason=(
            f"{call.reason} names type {call.part_type}, but the enzymes show this "
            f"is an assembled {role.replace('_', ' ')}, not a part"
        ),
    )


def _ecoli_marker(record: SeqRecord) -> str | None:
    haystack = " ".join(
        str(value)
        for feature in record.features
        for values in feature.qualifiers.values()
        for value in values
    )
    for marker, pattern in _ECOLI_MARKERS:
        if re.search(pattern, haystack, re.IGNORECASE):
            return marker
    return None


def _features(record: SeqRecord, limit: int = 200) -> list[FeatureInfo]:
    out = []
    for feature in record.features[:limit]:
        label = ""
        for key in ("label", "gene", "product", "note", "standard_name"):
            if feature.qualifiers.get(key):
                label = str(feature.qualifiers[key][0])
                break
        out.append(
            FeatureInfo(
                type=feature.type,
                label=label,
                start=int(feature.location.start),
                end=int(feature.location.end),
                strand=feature.location.strand,
            )
        )
    return out


def _circular_overlap(offset: int, size: int, width: int, length: int) -> int:
    """How much of a feature at `offset` lies in the span ``[0, width)``.

    Both the feature and the span may run past the origin, so the feature's
    tail is intersected with the span a second time after it wraps.
    """
    head = max(0, min(offset + size, width) - max(offset, 0))
    tail = offset + size - length
    return head + (max(0, min(tail, width)) if tail > 0 else 0)


def features_within(
    features: list[FeatureInfo],
    span: tuple[int, int] | None,
    length: int,
    coverage: float = 0.8,
) -> list[FeatureInfo]:
    """The features lying inside a span, which may run past the origin.

    A feature counts as inside when most of it is - `coverage` of its length by
    default, not all of it. Annotations are drawn by hand and by tools that
    disagree about where a fragment ends: a fragment carries its left overhang
    and not its right, while an annotator commonly draws the CDS between the two
    overhangs or across both. Demanding total containment discards a label over
    a handful of bases, which is how a plasmid annotated ``mTurquoise`` across
    its whole insert ends up described as nothing at all.

    A whole-plasmid annotation is still excluded: only a fifth of it is in the
    part, so it fails the same test from the other side.
    """
    if not span or not length:
        return list(features)
    start, end = span
    width = (end - start) % length or length
    inside = []
    for feature in features:
        size = max(feature.end - feature.start, 1)
        offset = (feature.start - start) % length
        if _circular_overlap(offset, size, width, length) >= coverage * size:
            inside.append(feature)
    return inside


def remap_features(
    features: list[FeatureInfo],
    span: tuple[int, int] | None,
    length: int,
    offset: int,
    part_length: int,
    limit: int = 8,
    min_length: int = 60,
    max_coverage: float = 0.9,
) -> list[FeatureInfo]:
    """The features inside `span`, moved into the coordinates of a product.

    A part sits at one place in its own plasmid and somewhere else entirely in
    the construct built from it, so a map of the product cannot reuse either
    set of numbers. `span` says where the fragment sits in the source; `offset`
    says where that fragment starts in the product; everything in between is
    the same arithmetic `features_within` already does, kept here so there is
    one answer to "is this feature in this part" rather than two.

    What comes back is the part's *substructure*, which is why so much is
    dropped. A feature covering essentially the whole part says nothing the
    part does not already say, and drawn on a ring just inside the part band it
    reads as a second, paler copy of that band - two rings of the same thing.
    So `max_coverage` filters those out, leaving the promoter inside a cassette
    and the origin inside a backbone: the things the part band cannot show.

    Small features go too, and only the largest survive: at ring scale a 20 bp
    binding site is a sliver that cannot be read, and eight per part is already
    more than the annulus can hold.
    """
    if not span or not length or not part_length:
        return []
    start, _ = span
    out: list[FeatureInfo] = []
    seen: set[tuple[int, int]] = set()
    for feature in features_within(features, span, length):
        size = feature.end - feature.start
        if size < min_length or size > part_length * max_coverage:
            continue
        if feature.type in _NOT_A_COMPONENT:
            continue
        # a label that names the tool rather than the thing is not a feature
        label = clean_label(feature.label) if feature.label else ""
        if not label or _LABEL_NOISE.search(label):
            continue
        # the same span annotated twice draws twice and looks like a heavier arc
        if (feature.start, feature.end) in seen:
            continue
        seen.add((feature.start, feature.end))
        at = offset + ((feature.start - start) % length)
        # a feature drawn across the fragment's edge is clipped, not dropped
        out.append(
            replace(
                feature,
                label=label,
                start=at,
                end=min(at + size, offset + part_length),
            )
        )
    out.sort(key=lambda f: f.start - f.end)  # largest first
    out = out[:limit]

    # A lone annotation covering most of its part is a restatement of that part,
    # not structure within it - `ConS` across a 194 bp connector draws an arc
    # the same size as the arc above it. Two or more annotations *are*
    # structure, however big any one of them is, which is why `His3` inside a
    # marker part and `ARS4` inside an origin part still draw.
    if len(out) == 1 and (out[0].end - out[0].start) > part_length * 0.5:
        return []
    return sorted(out, key=lambda f: f.start)


def summarise(features: list[FeatureInfo], deprioritise: str | None = None) -> str:
    """What a fragment holds, in the order a transcription unit is read.

    Promoter, then whatever it drives, then terminator - and for a part that is
    only one of those, just that one. Labels are the user's own; nothing here
    renames anything.
    """
    named, weak = [], []
    for feature in features:
        if not feature.label or feature.type in _NOT_A_COMPONENT:
            continue
        # judge the label on what it names, not on how it was made
        label = clean_label(feature.label)
        if _LABEL_NOISE.search(label):
            continue
        cleaned = replace(feature, label=label) if label != feature.label else feature
        (weak if feature.type in _WEAK_COMPONENT else named).append(cleaned)
    named = named or weak
    if not named:
        return ""

    # The E. coli selection marker is how you keep the plasmid, not what it is -
    # unless it is all that is in the fragment, which is exactly a type 8 part.
    if deprioritise:
        without = [f for f in named if deprioritise.lower() not in f.label.lower()]
        if without:
            named = without

    def first(kind: str) -> FeatureInfo | None:
        typed = [f for f in named if f.type == kind]
        if typed:
            return typed[0]
        return next((f for f in named if kind in f.label.lower()), None)

    promoter, terminator = first("promoter"), first("terminator")
    skip = {id(f) for f in (promoter, terminator) if f}
    rest = sorted((f for f in named if id(f) not in skip),
                  key=lambda f: f.start - f.end)
    cargo = rest[:2] if not promoter and not terminator else rest[:1]

    out: list[str] = []
    for feature in [promoter, *cargo, terminator]:
        if feature and feature.label not in out:
            out.append(feature.label)
    return " \u00b7 ".join(out[:3])


def contributed_span(
    record: SeqRecord, call: TypeCall, roles: _Roles, scheme: Scheme
) -> tuple[int, int] | None:
    """Where the fragment this plasmid contributes to a reaction sits."""
    if roles.cassette_span:
        return roles.cassette_span
    fragment = part_fragment(seqio.sequence(record), call, scheme)
    return (fragment.start, fragment.end) if fragment else None


def internal_multigene_sites(
    record: SeqRecord, span: tuple[int, int] | None, roles: _Roles, scheme: Scheme
) -> list[int]:
    """Multigene-enzyme sites inside the fragment that has to survive intact.

    A type 1 or 5 part is supposed to carry exactly one - that is its connector -
    so that one is not counted against it.
    """
    if not span:
        return []
    seq = seqio.sequence(record)
    length = len(seq)
    inside = [
        site.start
        for site in find_sites(seq, scheme.multigene_enzyme)
        if _enclosed(site, span[0], span[1], length)
    ]
    if roles.connector_overhang and len(inside) == 1:
        return []
    return sorted(inside)


#: Why a plasmid could not be typed. Each is a different job to fix, which is
#: the point of separating them: "no BsaI pair" usually means a cassette or a
#: construct, "linear" is a file-format problem, and "uncuttable" means the
#: plasmid is finished and cannot contribute to anything.
TRIAGE_NONE = ""
TRIAGE_NO_PAIR = "no_part_pair"
TRIAGE_UNCUTTABLE = "uncuttable"
TRIAGE_LINEAR = "linear"


def level3_readiness(call: TypeCall, connector_overhang: str | None) -> bool | None:
    """Whether a connector part can still take part in a multigene assembly.

    A cassette is released from its plasmid by the multigene enzyme cutting
    inside the type 1 and type 5 parts at its ends. Some parts in a real
    library have had that site domesticated away - `L13_ConLS_BSMB1del` says so
    in its own name, and its annotation reads "Former Bsmb1". A cassette built
    from one assembles perfectly at Level 2 and then cannot be cut out at Level
    3, and nothing says so until the multigene step fails to find its ends.

    Only positions 1 and 5 are judged: no other part contributes to the ends,
    so for them the question does not arise and the answer is `None`.
    """
    if _base_type(call.part_type) not in {"1", "5"}:
        return None
    return bool(connector_overhang)


def triage(entry: PlasmidEntry) -> str:
    """Why this plasmid has no part type, as a key a filter can group on.

    Derived from what the detector recorded, not from re-reading the file, and
    not by matching on the reason text - the reason is written for a person.
    """
    if entry.call.part_type or entry.is_assembled:
        return TRIAGE_NONE
    if entry.topology != "circular" and not entry.call.assumed_circular:
        return TRIAGE_LINEAR
    if entry.sites.part_enzyme_total < 2 and entry.sites.multigene_enzyme_total < 2:
        return TRIAGE_UNCUTTABLE
    return TRIAGE_NO_PAIR


def describe(
    record: SeqRecord,
    path: Path | None = None,
    relpath: str | None = None,
    scheme: Scheme = YTK,
) -> PlasmidEntry:
    """Everything the library knows about one plasmid."""
    call = detect_type(record, path, scheme)
    roles = analyse_roles(record, call, scheme)
    call = _reconcile(call, roles)
    record_id = "" if record.id in (None, ".", "<unknown id>") else str(record.id)
    features = _features(record)
    span = contributed_span(record, call, roles, scheme)
    part_span = span if not roles.cassette_span else None
    name = (path.stem if path else None) or record_id or "unnamed"

    component = summarise(
        features_within(features, span, len(record.seq)),
        # position 8 *is* the E. coli backbone, so there the marker is the
        # component rather than the packaging
        deprioritise=(
            None
            if "8" in parttypes.positions(call.part_type or "", scheme)
            else _ecoli_marker(record)
        ),
    )
    component_source = "annotation" if component else ""
    if not component:
        component = name_component(name)
        component_source = "filename" if component else ""

    return PlasmidEntry(
        path=relpath or (str(path) if path else record_id),
        name=name,
        record_id=record_id,
        length=len(record.seq),
        topology="circular" if call.assumed_circular else seqio.topology(record),
        call=call,
        sites=count_sites_for(record, call, roles, scheme),
        scheme=scheme.name,
        roles=roles.roles,
        connector_overhang=roles.connector_overhang,
        level3_ready=level3_readiness(call, roles.connector_overhang),
        cassette_overhangs=roles.cassette_overhangs,
        cassette_span=roles.cassette_span,
        ecoli_marker=_ecoli_marker(record),
        checksum=hashlib.sha1(seqio.sequence(record).encode()).hexdigest()[:12],
        multigene_overhangs=[
            f.left_overhang
            for f in digest(seqio.sequence(record), scheme.multigene_enzyme)
            if f.left_overhang
        ] if (seqio.is_circular(record) or call.assumed_circular) else [],
        part_span=part_span,
        internal_multigene_positions=internal_multigene_sites(record, span, roles, scheme),
        component=component,
        component_source=component_source,
        features=features,
    )


# --------------------------------------------------------------------------- #
# the indexed folders
# --------------------------------------------------------------------------- #


class Library:
    """One or more folders of GenBank files, indexed and cached.

    Several roots can be scanned as one library - a project keeps its parts in
    more than one place. Entry paths are relative to the roots' common parent,
    so they stay readable and unique across roots.

    The cache lives in ``<base>/.ggasm/index.json`` and is keyed on each file's
    path, mtime and size, so a rescan only re-reads what changed. Manual
    assignments in ``overrides.json`` always win over detection.
    """

    def __init__(
        self,
        roots: str | Path | Sequence[str | Path],
        recursive: bool = True,
        scheme: Scheme = YTK,
        cache_dir: str | Path | None = None,
    ):
        if isinstance(roots, (str, Path)):
            roots = [roots]
        self.roots = [Path(r).expanduser().resolve() for r in roots]
        if not self.roots:
            raise ValueError("a library needs at least one root folder")
        self.recursive = recursive
        self.scheme = scheme
        self.base = self._common_base()
        self._cache_dir = Path(cache_dir).expanduser().resolve() if cache_dir else None
        self.entries: dict[str, PlasmidEntry] = {}
        self.errors: dict[str, str] = {}
        self._cache: dict[str, dict[str, Any]] = {}
        self._primary: dict[str, PlasmidEntry] = {}
        self._fragments: list[tuple[PlasmidEntry, str]] | None = None
        self.overrides: dict[str, dict[str, Any]] = {}

    def _common_base(self) -> Path:
        if len(self.roots) == 1:
            return self.roots[0]
        return Path(os.path.commonpath([str(r) for r in self.roots]))

    # -- paths ------------------------------------------------------------- #

    @property
    def folder(self) -> Path:
        """The first root - the folder a single-root library was pointed at."""
        return self.roots[0]

    @property
    def cache_dir(self) -> Path:
        return self._cache_dir or (self.base / CACHE_DIR)

    @property
    def index_path(self) -> Path:
        return self.cache_dir / INDEX_FILE

    @property
    def overrides_path(self) -> Path:
        return self.cache_dir / OVERRIDES_FILE

    def files(self) -> list[Path]:
        pattern = "**/*" if self.recursive else "*"
        seen: dict[Path, None] = {}
        for root in self.roots:
            for path in sorted(root.glob(pattern)):
                if not path.is_file() or path.suffix.lower() not in seqio.GENBANK_SUFFIXES:
                    continue
                if SKIP_DIRS & set(path.parts):
                    continue
                seen.setdefault(path.resolve(), None)
        return sorted(seen)

    def relpath(self, path: Path) -> str:
        try:
            return path.relative_to(self.base).as_posix()
        except ValueError:
            return path.as_posix()

    # -- scanning ---------------------------------------------------------- #

    def scan(self, force: bool = False) -> list[PlasmidEntry]:
        self._load_cache()
        self._load_overrides()
        fresh: dict[str, dict[str, Any]] = {}
        self.entries, self.errors = {}, {}

        for path in self.files():
            rel = self.relpath(path)
            stat = path.stat()
            key = {"mtime_ns": stat.st_mtime_ns, "size": stat.st_size}
            cached = self._cache.get(rel)
            if not force and cached and _same(cached, key) and cached.get("scheme") == self.scheme.name:
                entry = PlasmidEntry.from_dict(cached["entry"])
            else:
                try:
                    entry = self._describe_file(path, rel)
                except Exception as exc:  # a bad file must not stop the scan
                    self.errors[rel] = f"{type(exc).__name__}: {exc}"
                    continue
            fresh[rel] = {**key, "scheme": self.scheme.name, "entry": entry.to_dict()}
            self.entries[rel] = self._apply_override(entry)

        self._cache = fresh
        self._recover_cassettes()
        self._save_cache()
        self._merge_duplicates()
        return self.sorted_entries()

    @property
    def connector_vocabulary(self) -> set[str]:
        """Connector overhangs this library knows, from its type 1 and 5 parts."""
        return {e.connector_overhang for e in self.entries.values() if e.connector_overhang}

    def _recover_cassettes(self) -> None:
        """Second pass: find cassettes that a stray internal site hid.

        The first pass can only see one plasmid at a time. Once the whole folder
        is indexed the connector overhangs are known, and a construct whose ends
        are two of them is a cassette even if the enzyme also cuts it somewhere
        in the middle - which is worth knowing, and worth warning about.
        """
        vocabulary = self.connector_vocabulary
        if not vocabulary:
            return
        for rel, entry in self.entries.items():
            if entry.call.part_type or entry.is_assembled or entry.length == 0:
                continue
            if len(entry.multigene_overhangs) < 2:
                continue
            try:
                record = seqio.read_records(self.base / rel)[0]
            except (OSError, ValueError):
                continue
            found = connector_pair(record, vocabulary, self.scheme)
            if not found:
                continue
            (five, three), span = found
            entry.roles.append(CASSETTE)
            entry.cassette_overhangs = (five, three)
            entry.cassette_span = span
            # now that the cassette's own span is known, its annotations beat
            # anything read off the file name earlier
            recovered = summarise(features_within(entry.features, span, entry.length))
            if recovered and entry.component_source != "annotation":
                entry.component, entry.component_source = recovered, "annotation"
            elif not entry.component:
                entry.component = name_component(entry.name)
                entry.component_source = "filename" if entry.component else ""
            entry.internal_multigene_positions = internal_multigene_sites(
                record, span, _Roles(roles=[CASSETTE]), self.scheme
            )
            inside = entry.internal_multigene_positions
            extra = len(inside)
            where = ", ".join(f"{p:,}" for p in inside)
            entry.call = TypeCall(
                part_type=None,
                confidence=NONE,
                reason=(
                    f"no {self.scheme.part_enzyme.name} pair, but {five} and {three} are "
                    f"connector overhangs: an assembled cassette"
                    + (f" with {extra} internal {self.scheme.multigene_enzyme.name} "
                       f"site(s) at {where}" if extra > 0 else "")
                ),
                conflict=(
                    f"internal {self.scheme.multigene_enzyme.name} site(s) at {where} "
                    f"would cut this cassette apart in a multigene reaction - "
                    f"domesticate before using it at Level 3"
                    if extra > 0 else None
                ),
            )
            cached = self._cache.get(rel)
            if cached:
                cached["entry"] = entry.to_dict()

    def _describe_file(self, path: Path, rel: str) -> PlasmidEntry:
        records = seqio.read_records(path)
        if not records:
            raise ValueError("no records in file")
        if len(records) > 1:
            self.errors[rel] = f"{len(records)} records; only the first was indexed"
        return describe(records[0], path=path, relpath=rel, scheme=self.scheme)

    def sorted_entries(self) -> list[PlasmidEntry]:
        """Every indexed file, duplicates included."""
        return sorted(self.entries.values(), key=lambda e: e.path)

    def unique_entries(self) -> list[PlasmidEntry]:
        """One entry per distinct sequence, carrying the names of all its copies.

        The same plasmid filed twice is one plasmid. Rather than make the user
        pick a copy in every dropdown, the copies are merged: one entry keeps
        the others' names as aliases, and takes the best description any of them
        had - so a file named after its plate well and a file named after the
        catalogue both contribute what they know.
        """
        return sorted(self._primary.values(), key=lambda e: e.path)

    def _merge_duplicates(self) -> None:
        groups: dict[str, list[PlasmidEntry]] = {}
        for entry in self.sorted_entries():
            groups.setdefault(entry.checksum, []).append(entry)

        self._primary = {}
        self._fragments = None
        for checksum, group in groups.items():
            group.sort(key=self._primary_rank)
            keeper, rest = group[0], group[1:]
            keeper.aliases = sorted({e.name for e in rest if e.name != keeper.name})
            if rest and not keeper.component:
                # the canonical copy said nothing; let a twin's annotations speak
                keeper.component = next((e.component for e in group if e.component), "")
            self._primary[checksum] = keeper

    def _primary_rank(self, entry: PlasmidEntry) -> tuple:
        """Which copy speaks for the group: chosen by hand, then folder, .gb, name.

        The automatic order is a guess - shortest name in the first root - and
        it is often wrong about which name a lab actually uses. A copy marked
        canonical beats every heuristic, so a dropdown can stop reading
        `pYTK003 / A3_ConL1` the moment someone says which is which.
        """
        chosen = not (self.overrides.get(entry.path) or {}).get("canonical")
        in_first_root = not str(self.path_of(entry)).startswith(str(self.roots[0]))
        return (
            chosen,
            not _SEQ_VERIFIED.search(entry.name),
            in_first_root,
            not entry.path.endswith(".gb"),
            len(entry.name),
            entry.name,
        )

    def set_canonical(self, relpath: str, canonical: bool = True) -> None:
        """Make this copy the one that speaks for its group of identical files."""
        self._load_overrides()
        entry = self.entries.get(relpath)
        if entry is not None and canonical:
            # only one copy of a sequence can be the canonical one
            for other in self.sorted_entries():
                if other.checksum == entry.checksum and other.path != relpath:
                    record = dict(self.overrides.get(other.path) or {})
                    if record.pop("canonical", None) is not None:
                        self._store_override(other.path, record)
                        self._load_overrides()

        record = dict(self.overrides.get(relpath) or {})
        if canonical:
            record["canonical"] = True
        else:
            record.pop("canonical", None)
        self._store_override(relpath, record)

    # -- overrides --------------------------------------------------------- #

    def _apply_override(self, entry: PlasmidEntry) -> PlasmidEntry:
        override = self.overrides.get(entry.path)
        if not override:
            return entry

        # A concentration and a hand-assigned type live in the same record but
        # are independent: measuring a prep must not blank the type, and typing
        # a plasmid must not forget what it was measured at.
        if override.get("conc_ng_ul") is not None:
            entry.conc_ng_ul = float(override["conc_ng_ul"])
        if "part_type" not in override:
            return entry

        detected = entry.call.part_type
        wanted = override.get("part_type")
        entry.call = TypeCall(
            part_type=wanted,
            confidence=override.get("confidence", DIGEST),
            reason=override.get("reason", "set by hand"),
            source="manual",
            five_prime=entry.call.five_prime,
            three_prime=entry.call.three_prime,
            reversed_sites=entry.call.reversed_sites,
            conflict=f"detection said {detected}" if detected and detected != wanted else None,
        )
        return entry

    def set_override(self, relpath: str, part_type: str | None, reason: str = "set by hand") -> None:
        self._load_overrides()
        record = dict(self.overrides.get(relpath) or {})
        if part_type is None:
            record.pop("part_type", None)
            record.pop("reason", None)
        else:
            record.update({"part_type": part_type, "reason": reason})
        self._store_override(relpath, record)

    def set_concentration(self, relpath: str, conc_ng_ul: float | None) -> None:
        """Record what a prep of this plasmid measured at, in ng/µL.

        Equimolar Golden Gate needs a concentration per plasmid, and that is a
        property of the tube on your bench, not of the sequence - so it cannot
        be read out of the GenBank file and has to be kept beside it.
        """
        self._load_overrides()
        record = dict(self.overrides.get(relpath) or {})
        if conc_ng_ul is None:
            record.pop("conc_ng_ul", None)
        else:
            if conc_ng_ul <= 0:
                raise ValueError("a concentration must be greater than zero")
            record["conc_ng_ul"] = float(conc_ng_ul)
        self._store_override(relpath, record)

    def _store_override(self, relpath: str, record: dict[str, Any]) -> None:
        """Write one override record back, dropping it when nothing is left."""
        if record:
            self.overrides[relpath] = record
        else:
            self.overrides.pop(relpath, None)
        self._write_json(self.overrides_path, self.overrides)
        cached = self._cache.get(relpath)
        if relpath in self.entries and cached:
            self.entries[relpath] = self._apply_override(PlasmidEntry.from_dict(cached["entry"]))
            # `_primary` holds references into `self.entries`, and the line above
            # replaced one with a fresh object. Without re-merging, every screen
            # that reads the deduplicated view - which is all of them - keeps
            # showing the entry as it was before the edit.
            self._merge_duplicates()

    # -- queries ----------------------------------------------------------- #

    def get(self, name: str) -> PlasmidEntry | None:
        """Find a plasmid by name, path, record id, or the name of any copy."""
        for entry in self.unique_entries():
            if name in (entry.name, entry.path, entry.record_id, *entry.aliases):
                return entry
        for entry in self.sorted_entries():
            if name in (entry.name, entry.path, entry.record_id):
                return entry
        return None

    def path_of(self, entry: PlasmidEntry) -> Path:
        return self.base / entry.path

    def record(self, entry: PlasmidEntry) -> SeqRecord:
        """Re-read the plasmid behind an index entry."""
        return seqio.read_records(self.path_of(entry))[0]

    def piece(self, entry: PlasmidEntry, enzyme: Enzyme | None = None) -> Piece:
        """The fragment this plasmid contributes to a reaction, ready to assemble.

        Every level needs this and none of them should re-derive it: read the
        file, digest it, and take the fragment the reaction will actually carry -
        features and all.

        With the scheme's part enzyme (the default) that is the typed part. With
        the multigene enzyme it is whatever the plasmid releases site-free: a
        cassette's transcription unit, or a destination vector's backbone.
        """
        record = self.record(entry)
        if enzyme is not None and enzyme is not self.scheme.part_enzyme:
            return self._released_piece(entry, record, enzyme)
        fragment = part_fragment(seqio.sequence(record), entry.call, self.scheme)
        if fragment is None:
            raise ValueError(f"{entry.name} has no {self.scheme.part_enzyme.name} part to release")
        sub = seqio.slice_record(record, fragment.start, fragment.end, name=entry.name)
        return Piece(
            record=sub,
            left_overhang=fragment.left_overhang or "",
            right_overhang=fragment.right_overhang or "",
            source_name=entry.name,
            part_type=entry.call.part_type,
            color=color_for(entry.call.part_type, self.scheme),
            component=entry.component,
        )

    def _released_piece(self, entry: PlasmidEntry, record: SeqRecord, enzyme: Enzyme) -> Piece:
        """The one site-free fragment `enzyme` releases from this plasmid."""
        released = [f for f in digest(seqio.sequence(record), enzyme)
                    if not f.has_sites and f.left_overhang]
        if len(released) != 1:
            raise ValueError(
                f"{entry.name} releases {len(released)} site-free fragments with "
                f"{enzyme.name}; exactly one is needed"
            )
        fragment = released[0]
        sub = seqio.slice_record(record, fragment.start, fragment.end, name=entry.name)
        return Piece(
            record=sub,
            left_overhang=fragment.left_overhang or "",
            right_overhang=fragment.right_overhang or "",
            source_name=entry.name,
            part_type=entry.call.part_type,
            color=color_for(entry.call.part_type, self.scheme),
            component=entry.component,
        )

    def part_fragments(self) -> list[tuple[PlasmidEntry, str]]:
        """Every part in the library with the exact sequence it contributes.

        Cached, because the caller that needs it - reading an assembled
        construct back into parts - would otherwise re-read and re-digest every
        part file on the shelf for each construct it looks at. Cleared by
        `scan`, which is the only thing that can change the answer.
        """
        if self._fragments is not None:
            return self._fragments

        out: list[tuple[PlasmidEntry, str]] = []
        for entry in self.unique_entries():
            if not entry.is_part or not entry.part_span:
                continue
            try:
                sequence = str(self.piece(entry).record.seq).upper()
            except (ValueError, OSError):
                continue  # a file that moved, or a part that will not release
            if sequence:
                out.append((entry, sequence))
        self._fragments = out
        return out

    def labels_in_part(self, entry: PlasmidEntry, limit: int = 6) -> list[FeatureInfo]:
        """The annotations inside the fragment this plasmid contributes.

        What a picker should match on: the backbone's AmpR is in every plasmid
        in the kit and tells you nothing about which one you want.
        """
        span = entry.part_span or entry.cassette_span
        inside = features_within(entry.features, span, entry.length)
        out, seen = [], set()
        for feature in inside:
            if not feature.label or feature.type in _NOT_A_COMPONENT:
                continue
            label = clean_label(feature.label)
            if _LABEL_NOISE.search(label) or label in seen:
                continue
            seen.add(label)
            out.append(replace(feature, label=label))
            if len(out) >= limit:
                break
        return out

    def parts_of_type(self, part_type: str) -> list[PlasmidEntry]:
        return [e for e in self.unique_entries() if e.call.part_type == part_type]

    def parts_with_overhangs(self, five: str, three: str) -> list[PlasmidEntry]:
        """Every part whose released fragment has exactly this overhang pair."""
        return [
            e
            for e in self.unique_entries()
            if e.overhangs == (five.upper(), three.upper()) and not e.call.reversed_sites
        ]

    def destinations(self) -> list[Destination]:
        """Every plasmid a new part could be cloned into, and how.

        A part plasmid is itself a destination for another part of the same
        type - cut out what is in it and put the replacement in the same place -
        so the list is normally long, not one entry.
        """
        out: list[Destination] = []
        scheme = self.scheme
        for entry in self.unique_entries():
            try:
                record = self.record(entry)
            except (OSError, ValueError):
                continue
            out.extend(self._destinations_for(entry, record))
        out.sort(key=lambda d: (d.part_type or "~", d.name))
        return out

    #: A destination has to be a plasmid in its own right: it must keep an
    #: origin and a marker after the insert goes in, so a backbone much smaller
    #: than this is a fragment that happens to have two cuts, not a vector.
    MIN_BACKBONE = 1500
    #: And it has to give something up in exchange, not just nick itself.
    MIN_DROPOUT = 50

    def _destinations_for(self, entry: PlasmidEntry, record: SeqRecord) -> list[Destination]:
        scheme = self.scheme
        seq = seqio.sequence(record)
        if not (seqio.is_circular(record) or entry.call.assumed_circular):
            return []
        found: list[Destination] = []

        # a dropout vector, by whichever enzyme releases the dropout
        for enzyme in (scheme.multigene_enzyme, *scheme.alternate_cloning_enzymes):
            fragments = digest(seq, enzyme)
            if len(fragments) != 2:
                continue
            keep = [f for f in fragments if not f.has_sites and f.left_overhang]
            drop = [f for f in fragments if f.has_sites]
            if len(keep) != 1 or len(drop) != 1:
                continue
            if len(keep[0].seq) < self.MIN_BACKBONE or len(drop[0].seq) < self.MIN_DROPOUT:
                continue
            accepts = (drop[0].left_overhang or "", drop[0].right_overhang or "")
            call = span(*accepts, scheme=scheme)
            found.append(
                Destination(
                    name=entry.name,
                    path=entry.path,
                    cloning_enzyme=enzyme.name,
                    accepts=accepts,
                    part_type=call.name if call and call.legal else None,
                    supplies_part_sites=bool(
                        find_sites(keep[0].seq, scheme.part_enzyme, circular=False)
                    ),
                    kind="dropout",
                    backbone_length=len(keep[0].seq),
                    marker=entry.ecoli_marker,
                    component=entry.component,
                    aliases=entry.aliases,
                )
            )

        # an existing part plasmid, reopened with the part enzyme
        if entry.is_part and entry.call.five_prime and entry.call.three_prime:
            fragments = digest(seq, scheme.part_enzyme)
            backbone = [f for f in fragments if f.has_sites]
            if (
                len(fragments) == 2
                and len(backbone) == 1
                and len(backbone[0].seq) >= self.MIN_BACKBONE
            ):
                found.append(
                    Destination(
                        name=entry.name,
                        path=entry.path,
                        cloning_enzyme=scheme.part_enzyme.name,
                        accepts=(entry.call.five_prime, entry.call.three_prime),
                        part_type=entry.call.part_type,
                        supplies_part_sites=True,
                        kind="part_plasmid",
                        backbone_length=len(backbone[0].seq),
                        marker=entry.ecoli_marker,
                        component=entry.component,
                        aliases=entry.aliases,
                    )
                )
        return found

    def cassettes(self) -> list[PlasmidEntry]:
        """Multigene inputs: a ConL end, a ConR end, and a TU between them."""
        return [e for e in self.unique_entries() if e.is_cassette]

    def multigene_vectors(self) -> list[PlasmidEntry]:
        """Multigene backbones: connector ends with a dropout between them."""
        return [e for e in self.unique_entries() if e.is_multigene_vector]

    def entry_vectors(self) -> list[PlasmidEntry]:
        return [e for e in self.unique_entries() if e.is_entry_vector]

    def connector_overhangs(self) -> dict[str, list[str]]:
        """Connector overhang -> the parts carrying it, learned from the files."""
        out: dict[str, list[str]] = {}
        for entry in self.sorted_entries():
            if entry.connector_overhang:
                out.setdefault(entry.connector_overhang, []).append(entry.name)
        return out

    def duplicates(self) -> dict[str, list[str]]:
        """Identical sequences filed under more than one name."""
        by_sum: dict[str, list[str]] = {}
        for entry in self.sorted_entries():
            by_sum.setdefault(entry.checksum, []).append(entry.path)
        return {k: v for k, v in by_sum.items() if len(v) > 1}

    def duplicate_groups(self) -> list[dict[str, Any]]:
        """Each group of identical files, and which copy currently speaks for it."""
        out = []
        for checksum, paths in sorted(self.duplicates().items()):
            keeper = self._primary.get(checksum)
            out.append({
                "checksum": checksum,
                "primary": keeper.path if keeper else paths[0],
                "chosen_by_hand": bool(
                    keeper and (self.overrides.get(keeper.path) or {}).get("canonical")
                ),
                "copies": [
                    {
                        "path": path,
                        "name": self.entries[path].name,
                        "length": self.entries[path].length,
                        "component": self.entries[path].component,
                    }
                    for path in sorted(paths)
                    if path in self.entries
                ],
            })
        return out

    # -- persistence ------------------------------------------------------- #

    def _load_cache(self) -> None:
        data = self._read_json(self.index_path)
        self._cache = data.get("entries", {}) if data.get("version") == INDEX_VERSION else {}

    def _save_cache(self) -> None:
        self._write_json(
            self.index_path,
            {
                "version": INDEX_VERSION,
                "scheme": self.scheme.name,
                "roots": [str(r) for r in self.roots],
                "entries": self._cache,
            },
        )

    # -- the build log ----------------------------------------------------- #

    @property
    def builds_path(self) -> Path:
        return self.cache_dir / BUILDS_FILE

    def log_build(
        self,
        name: str,
        level: str,
        parts: Sequence[PlasmidEntry],
        length: int,
        issues: Sequence[str] = (),
    ) -> dict[str, Any]:
        """Record that a construct was exported, and exactly what went into it.

        Append-only, one JSON object per line, because the question it answers
        comes months later: *which* pYTK009 was in the thing I built in March,
        and has that file changed since. So each part is logged with the
        checksum of its sequence rather than only its name - a name is not
        evidence, and a file on a shared drive can be edited under you.
        """
        import datetime

        record = {
            "at": datetime.datetime.now().astimezone().isoformat(timespec="seconds"),
            "name": name,
            "level": level,
            "length": length,
            "parts": [
                {
                    "name": e.name,
                    "path": e.path,
                    "part_type": e.call.part_type,
                    "length": e.length,
                    "sha1": e.checksum,
                }
                for e in parts
            ],
            "issues": list(issues),
        }
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        with self.builds_path.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(record) + "\n")
        return record

    def builds(self, limit: int = 100) -> list[dict[str, Any]]:
        """The most recent builds, newest first. A damaged line is skipped."""
        if not self.builds_path.exists():
            return []
        out: list[dict[str, Any]] = []
        for line in self.builds_path.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if not line:
                continue
            try:
                out.append(json.loads(line))
            except json.JSONDecodeError:
                continue  # a torn write should not hide the rest of the log
        return list(reversed(out))[:limit]

    # -- saved designs ----------------------------------------------------- #

    @property
    def designs_path(self) -> Path:
        return self.cache_dir / DESIGNS_FILE

    def designs(self, level: str | None = None) -> list[dict[str, Any]]:
        """Every saved design, newest first.

        A design is a handful of plasmid names and view flags - it costs
        nothing to keep and everything to lose, since retyping eight choices
        after a restart is what makes a tool feel disposable. Stored beside the
        index rather than in it, because it describes what you are building
        rather than what is on the shelf.
        """
        saved = self._read_json(self.designs_path).get("designs", [])
        if level:
            saved = [d for d in saved if d.get("level") == level]
        return sorted(saved, key=lambda d: d.get("saved_at", ""), reverse=True)

    def save_design(self, level: str, name: str, design: dict[str, Any]) -> dict[str, Any]:
        """Keep a design under `name`, replacing any earlier one of that name."""
        name = name.strip()
        if not name:
            raise ValueError("a design needs a name")

        import datetime

        record = {
            "level": level,
            "name": name,
            "saved_at": datetime.datetime.now().astimezone().isoformat(timespec="seconds"),
            "design": design,
        }
        data = self._read_json(self.designs_path)
        kept = [
            d for d in data.get("designs", [])
            if not (d.get("level") == level and d.get("name") == name)
        ]
        data["designs"] = [record, *kept]
        self._write_json(self.designs_path, data)
        return record

    def delete_design(self, level: str, name: str) -> bool:
        data = self._read_json(self.designs_path)
        before = data.get("designs", [])
        after = [
            d for d in before
            if not (d.get("level") == level and d.get("name") == name)
        ]
        data["designs"] = after
        self._write_json(self.designs_path, data)
        return len(after) < len(before)

    def _load_overrides(self) -> None:
        self.overrides = self._read_json(self.overrides_path)

    @staticmethod
    def _read_json(path: Path) -> dict[str, Any]:
        try:
            return json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            return {}

    @staticmethod
    def _write_json(path: Path, data: dict[str, Any]) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(data, indent=2, sort_keys=True), encoding="utf-8")


def _same(cached: dict[str, Any], key: dict[str, Any]) -> bool:
    return cached.get("mtime_ns") == key["mtime_ns"] and cached.get("size") == key["size"]
