"""Level 1 - part-plasmid construction.

Take a stretch of template DNA and turn it into a part plasmid of a chosen
type: design the primers that add the flanking sites, scan for internal sites
that would ruin a later assembly, and predict the plasmid that comes out.

The two enzymes do different jobs here, and the primer tails encode both at
once. Reading the finished amplicon's top strand from the left:

    GCAT  CGTCTC n TCGG   GGTCTC a  [5' overhang]  ...part...
    pad   BsmBI    entry  BsaI       type flank
          (cuts out the amplicon)    (cuts the part out of the finished plasmid)

and the right end mirrors it, with both sites reversed so that each cuts
*inward*. So one PCR product carries the entry-vector ends that clone it, and,
one layer inside, the part-type ends that will release it again at Level 2.

`TCGG` ends in `GG` and `GGTCTC` begins with `GG`, so the two overlap by two
bases; every tail is built overlapped and then verified by simulating the
digest, falling back to the unoverlapped form if the shorter one does not cut
exactly right.

This module imports no other level.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace

from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import MeltingTemp

from ..core import parttypes, seqio
from ..core.assembly import ERROR, INFO, WARNING, AssemblyResult, Issue, Piece, assemble
from ..core.enzymes import Enzyme, find_sites
from ..core.enzymes import get_enzyme
from ..core.library import Destination, Library, PlasmidEntry, describe
from ..core.parttypes import YTK, Scheme

#: Pads, spacers and linkers taken from the L131 entry construct. The pads give
#: the outer enzyme somewhere to bind at the very end of a fragment, and the
#: linkers hold the outer cut away from the inner site so both enzymes can work.
PAD = "GCAT"
PAD_5 = "CGTAGTC"
PAD_3 = "GCTGCCA"
LINKER_5 = "CGGGAGGAAGTCTTTAGACC"
LINKER_3 = "GCCGAGGTTTCTCATTCAGCA"
SPACER = "A"
#: The spacer bases between a site and its overhang are arbitrary - the enzyme
#: only counts them - but these are the ones the L131 construct uses, so output
#: from this app reads the same as the constructs already on the bench.
SPACERS: dict[str, str] = {"BbsI": "TA"}
MIN_ANNEAL, MAX_ANNEAL = 18, 25
TARGET_TM = 60.0
TM_SPREAD = 3.0

#: How far either side of a site to look for a usable junction overhang. A
#: primer can only carry a mutation that falls inside its annealing region, so
#: the junction has to land within one primer's reach of the site it fixes.
JUNCTION_SEARCH = MAX_ANNEAL

#: Part types whose body is, by definition, coding sequence. Their 5' overhang
#: carries the frame: TATG supplies the ATG, and the 3a/3b junction GGT|TCT
#: lands on a codon boundary, so in every case the body starts in frame at 0.
CODING_TYPES = ("3", "3a", "3b")

STANDARD_CODONS = CodonTable.unambiguous_dna_by_id[1]
STOP_CODONS = ("TAA", "TAG", "TGA")


@dataclass
class Convention:
    """The YTK sequence conventions for particular part types."""

    gly_ser_linker: bool = True
    """Type 3/3b: append ``GG`` before ``ATCC`` so the junction reads ``GGATCC``."""
    strip_stop: bool = True
    """Type 3: drop a trailing stop codon, since the terminator follows."""
    stop_and_xhoi: bool = True
    """Type 4/4a: begin with ``TAA`` then ``CTCGAG``."""
    strip_start: bool = True
    """Type 3/3a: drop a leading ``ATG``, because the ``TATG`` overhang is one."""
    infer_from_sequence: bool = True
    """Let a terminal stop codon settle whether this CDS is meant to be fused.

    A CDS handed over *with* its stop codon, in frame, is a CDS that ends there.
    Stripping it and adding the Gly-Ser ``GG`` turns it into a fusable part -
    a different molecule, and one whose primers will not match the sequence that
    was pasted. So when a terminal in-frame stop is present, the fusion
    conventions stay out of the way unless they are asked for explicitly."""
    """Type 3/3a: drop a leading ``ATG``, because the ``TATG`` overhang is one.

    Paste a CDS straight out of a genome browser and it begins with its start
    codon; the part's own 5' overhang already spells it, so keeping both would
    put two methionines at the front."""


@dataclass
class PartRequest:
    """What the user asked for on the Level 1 screen."""

    part_type: str
    template: str | None = None
    """A plasmid name in the library."""
    sequence: str | None = None
    """A pasted sequence, used when `template` is None."""
    start: int = 1
    """1-based, inclusive."""
    end: int = 0
    """1-based, inclusive; 0 means to the end of the template."""
    destination: str | None = None
    """The plasmid to clone into. Any part plasmid works, as does a ccdB or GFP
    dropout vector; leaving it unset picks the best one in the library."""
    entry_vector: str | None = None
    """Deprecated name for `destination`."""
    mode: str = "pcr"
    """``pcr``, ``gblock`` or ``oligo``."""
    trim_to_insert: bool = True
    """Take the insert out of a sequence that is already a finished construct.

    Paste a whole entry construct and the sequence you actually want is the bit
    between its part-enzyme sites; wrapping new tails around the lot would build
    a construct around a construct.
    """
    domesticate: bool = False
    """Remove internal sites by silent mutation.

    Off by default, and deliberately so. In a coding sequence a synonymous codon
    is an obvious fix; in a promoter or a terminator there is no such thing as a
    silent change, and only the person who knows what the element does can say
    which base is safe to move. So the scan always runs and always reports what
    it would do - applying it is a decision, not a default.
    """
    name: str = "new_part"
    conventions: Convention = field(default_factory=Convention)


@dataclass
class Segment:
    """One labelled stretch of a primer, for the coloured display."""

    role: str
    seq: str


@dataclass
class Primer:
    name: str
    sequence: str
    annealing: str
    tm: float
    segments: list[Segment] = field(default_factory=list)

    @property
    def length(self) -> int:
        return len(self.sequence)


@dataclass
class Mutation:
    """A base change that destroys an internal site."""

    position: int
    """0-based, in the part body."""
    from_base: str
    to_base: str
    enzyme: str
    silent: bool
    gene: str = ""
    codon_from: str = ""
    codon_to: str = ""

    @property
    def description(self) -> str:
        where = f" in {self.gene}" if self.gene else ""
        how = (
            f"silent: {self.codon_from}->{self.codon_to}"
            if self.silent
            else "CHANGES THE PROTEIN" if self.codon_from else "not in a CDS"
        )
        return f"{self.enzyme} site at {self.position}{where}: {self.from_base}->{self.to_base} ({how})"


@dataclass
class FragmentPlan:
    """One PCR product; more than one means the part had to be domesticated."""

    index: int
    start: int
    end: int
    left_overhang: str
    right_overhang: str
    forward: Primer | None = None
    reverse: Primer | None = None
    mutations: list[Mutation] = field(default_factory=list)

    @property
    def length(self) -> int:
        return self.end - self.start


@dataclass
class Plan:
    """How a particular destination wants its insert built."""

    destination: Destination
    outer: Enzyme
    """The enzyme that cuts the amplicon into the destination."""
    inner: Enzyme | None
    """The part enzyme, when the insert has to bring its own sites; None when
    the destination's backbone already carries them."""
    pad_5: str = PAD_5
    pad_3: str = PAD_3
    linker_5: str = LINKER_5
    linker_3: str = LINKER_3

    @property
    def accepts(self) -> tuple[str, str]:
        return self.destination.accepts

    @property
    def carries_part_sites(self) -> bool:
        return self.inner is not None

    def forward_tail(self, five: str, added: str = "") -> str:
        """`added` is sequence the conventions introduced, not read off the
        template, so it belongs here rather than in the annealing region."""
        return _tail(self.accepts[0], five, self.outer, self.inner,
                     pad=self.pad_5, linker=self.linker_5 if self.inner else "",
                     added=added)

    def reverse_tail(self, three: str, added: str = "") -> str:
        return _tail(seqio.revcomp(self.accepts[1]), seqio.revcomp(three), self.outer,
                     self.inner, pad=seqio.revcomp(self.pad_3),
                     linker=seqio.revcomp(self.linker_3) if self.inner else "",
                     added=seqio.revcomp(added) if added else "")


@dataclass
class Level1Design:
    """Everything the Level 1 screen shows, and everything a bench user needs."""

    part_type: str
    five_prime: str
    three_prime: str
    body: str
    fragments: list[FragmentPlan] = field(default_factory=list)
    issues: list[Issue] = field(default_factory=list)
    product: SeqRecord | None = None
    entry_vector: str | None = None
    destination: str | None = None
    cloning_enzyme: str = ""
    accepts: tuple[str, str] | None = None
    destination_marker: str | None = None
    destination_universal: bool = False
    ligation: dict = field(default_factory=dict)
    binding: list = field(default_factory=list)
    validated_as: str | None = None
    oligos: list[Primer] = field(default_factory=list)
    gblock: str = ""
    mode: str = "pcr"

    @property
    def ok(self) -> bool:
        return self.product is not None and not self.errors

    @property
    def errors(self) -> list[Issue]:
        return [i for i in self.issues if i.level == ERROR]

    @property
    def warnings(self) -> list[Issue]:
        return [i for i in self.issues if i.level == WARNING]

    @property
    def primers(self) -> list[Primer]:
        out = []
        for fragment in self.fragments:
            out.extend(p for p in (fragment.forward, fragment.reverse) if p)
        return out

    def codes(self) -> set[str]:
        return {i.code for i in self.issues}


# --------------------------------------------------------------------------- #
# conventions
# --------------------------------------------------------------------------- #


def trim_to_insert(
    body: str, part_type: str, scheme: Scheme
) -> tuple[str, list[Issue]]:
    """If this is already a construct, return the part inside it.

    A finished entry construct carries two part-enzyme sites pointing inward at
    a known part type's overhangs. That is a fingerprint nothing else has, so
    when it is present the sequence is not an insert - it is something already
    built, and the insert is what sits between those sites.
    """
    sites = find_sites(body, scheme.part_enzyme, circular=False)
    forward = [s for s in sites if s.strand == "+"]
    reverse = [s for s in sites if s.strand == "-"]
    if not forward or not reverse:
        return body, []

    left, right = forward[0], reverse[-1]
    if right.top_cut <= left.top_cut:
        return body, []

    five = body[left.top_cut:left.top_cut + scheme.overhang_length]
    three = body[right.top_cut:right.top_cut + scheme.overhang_length]
    found = parttypes.span(five, three, scheme=scheme)
    if found is None or not found.legal:
        return body, []

    inner = body[left.top_cut + scheme.overhang_length:right.top_cut]
    trimmed = len(body) - len(inner)
    note = (
        f"this looks like a finished construct, not an insert: it already carries "
        f"{scheme.part_enzyme.name} sites releasing a type {found.name} part. Using "
        f"the {len(inner):,} bp between them and dropping {trimmed:,} bp of "
        f"flanking sequence"
    )
    if found.name != part_type:
        note += f". Note it was built as a type {found.name}, and you asked for {part_type}"
    return inner, [Issue(INFO, "trimmed_to_insert", note)]


def apply_conventions(
    body: str, part_type: str, convention: Convention
) -> tuple[str, str, str, list[Issue]]:
    """Fold in the YTK sequence rules for this part type.

    Returns the template-derived core and the bases the conventions added either
    side of it, kept separate on purpose. Anything added here does not exist on
    the template, so it has to be spelled by a primer *tail* - put it in the
    annealing region instead and the primer's 3' end cannot bind.
    """
    issues: list[Issue] = []
    body = body.upper()
    prefix = suffix = ""

    terminal_stop = (
        part_type in ("3", "3b")
        and len(body) % 3 == 0
        and body[-3:] in STOP_CODONS
    )
    if terminal_stop and convention.infer_from_sequence:
        convention = replace(convention, strip_stop=False, gly_ser_linker=False)
        issues.append(
            Issue(INFO, "terminal_cds",
                  f"this sequence ends in an in-frame {body[-3:]}, so it is being kept as "
                  f"a CDS that ends here: the stop stays and no Gly-Ser GG is added. "
                  f"Tick \u201cGly-Ser linker\u201d to make a fusable part instead")
        )

    if part_type in ("3", "3a") and convention.strip_start and body.startswith("ATG"):
        body = body[3:]
        issues.append(
            Issue(INFO, "start_stripped",
                  "dropped the leading ATG: the TATG overhang already supplies the "
                  "start codon")
        )

    if part_type in ("3", "3b") and convention.strip_stop:
        # Only strip a stop that really is one. The last three bases spelling
        # TAA does not make them a codon: unless the body is a whole number of
        # codons, those bases fall across a boundary and deleting them would
        # shift everything upstream of the junction.
        in_frame = len(body) % 3 == 0
        last = body[-3:]
        if last in STOP_CODONS and in_frame:
            body = body[:-3]
            issues.append(
                Issue(INFO, "stop_stripped",
                      f"dropped the trailing {last} stop codon; the terminator part "
                      f"supplies one. Turn this off for a CDS that should end here")
            )
        elif last in STOP_CODONS and not in_frame:
            issues.append(
                Issue(WARNING, "stop_out_of_frame",
                      f"the last three bases read {last}, but the sequence is "
                      f"{len(body):,} bp - not a whole number of codons - so they are "
                      f"not a stop codon and have been left alone")
            )

    if part_type in ("3", "3b") and convention.gly_ser_linker:
        suffix = "GG"
        issues.append(
            Issue(INFO, "gly_ser",
                  "appended GG before ATCC, giving the GGATCC (Gly-Ser / BamHI) junction")
        )
    elif part_type == "3a" and convention.gly_ser_linker:
        suffix = "GG"
        issues.append(Issue(INFO, "gly_ser", "appended GG before TTCT for the 3a/3b read-through"))

    if part_type in ("4", "4a") and convention.stop_and_xhoi:
        prefix = "TAACTCGAG"
        issues.append(Issue(INFO, "stop_xhoi", "prefixed TAA + CTCGAG (stop, then XhoI)"))

    return body, prefix, suffix, issues


# --------------------------------------------------------------------------- #
# primers
# --------------------------------------------------------------------------- #


def tm(seq: str) -> float:
    if len(seq) < 8:
        return 0.0
    return round(float(MeltingTemp.Tm_NN(Seq(seq))), 1)


def _pick_annealing(forward_source: str, reverse_source: str) -> tuple[str, str]:
    """Choose annealing lengths that are close to target and close to each other."""
    best, best_score = None, None
    for lf in range(MIN_ANNEAL, MAX_ANNEAL + 1):
        if lf > len(forward_source):
            break
        f = forward_source[:lf]
        tf = tm(f)
        for lr in range(MIN_ANNEAL, MAX_ANNEAL + 1):
            if lr > len(reverse_source):
                break
            r = reverse_source[:lr]
            tr = tm(r)
            score = abs(tf - tr) + 0.5 * (abs(tf - TARGET_TM) + abs(tr - TARGET_TM))
            if best_score is None or score < best_score:
                best, best_score = (f, r), score
    if best is None:
        return forward_source, reverse_source
    return best


def _spacer(enzyme: Enzyme) -> str:
    """The filler between a site and the overhang it cuts out."""
    width = enzyme.top_offset - len(enzyme.site)
    preset = SPACERS.get(enzyme.name, "")
    return preset[:width] if len(preset) >= width else SPACER * width


def _tail(
    outer_overhang: str,
    part_overhang: str,
    outer: Enzyme,
    inner: Enzyme | None,
    pad: str = PAD,
    linker: str = "",
    added: str = "",
) -> str:
    """The non-annealing head of a primer.

    Reading outwards from the annealing sequence: the part type's own overhang,
    then - only when the destination does not already carry them - the part
    enzyme's site that will release the part later, a linker, the overhang the
    destination accepts, and the outer enzyme that cuts the amplicon into it.

    With no linker the outer overhang and the inner site are adjacent and can
    share bases (``TCGG`` ends in ``GG``, ``GGTCTC`` begins with it), which is
    how the shorter pYTK001-style tail is built. The caller proves any shortened
    form still cuts correctly.
    """
    spacer = _spacer(outer)
    head = pad + outer.site + spacer + outer_overhang
    if inner is None:
        return head + added

    inner_spacer = _spacer(inner)
    inner_part = inner.site + inner_spacer + part_overhang
    if linker:
        return head + linker + inner_part + added
    for overlap in range(min(len(outer_overhang), len(inner.site)), 0, -1):
        if outer_overhang.endswith(inner.site[:overlap]):
            return head + inner_part[overlap:] + added
    return head + inner_part + added


def _segments(plan: Plan, outer_overhang: str, part_overhang: str, scheme: Scheme,
              tail: str, anneal: str) -> list[Segment]:
    """Label the pieces of a tail for the coloured display.

    Walks the tail left to right with a cursor rather than searching it, because
    the same four bases can legitimately appear twice - with the L131 vector the
    overhang it accepts and the type 3 flank are both ``TATG`` - and searching
    would label whichever came first as both.
    """
    order = [
        ("pad", tail[:len(plan.pad_5)] if tail.startswith(plan.pad_5) else PAD),
        ("multigene_enzyme", plan.outer.site),
        ("entry_overhang", outer_overhang),
        ("part_enzyme", plan.inner.site if plan.inner else ""),
        ("type_flank", part_overhang if plan.inner else ""),
    ]
    segments: list[Segment] = []
    cursor = 0
    for role, seq in order:
        if not seq:
            continue
        index = tail.find(seq, cursor)
        if index < 0:
            continue
        segments.append(Segment(role, seq))
        cursor = index + len(seq)
    segments.append(Segment("annealing", anneal))
    return segments


def design_primers(
    body: str,
    five: str,
    three: str,
    plan: Plan,
    scheme: Scheme,
    name: str,
    prefix: str = "",
    suffix: str = "",
) -> tuple[Primer, Primer]:
    """The pair that turns `body` into an amplicon ready for the destination.

    `body` here is the stretch of *template*; `prefix` and `suffix` are bases the
    conventions added, which no template carries and which therefore go in the
    tails. Annealing is chosen from the template bases alone, so both primers
    end on sequence that is really there to bind.
    """
    forward_tail = plan.forward_tail(five, prefix)
    # the right end is the mirror image, so the reverse primer is built on the
    # reverse complement of the three-prime overhang
    reverse_tail = plan.reverse_tail(three, suffix)

    forward_anneal, reverse_anneal = _pick_annealing(body, seqio.revcomp(body))

    forward = Primer(
        name=f"{name}_F",
        sequence=forward_tail + forward_anneal,
        annealing=forward_anneal,
        tm=tm(forward_anneal),
        segments=_segments(plan, plan.accepts[0], five, scheme, forward_tail, forward_anneal),
    )
    reverse = Primer(
        name=f"{name}_R",
        sequence=reverse_tail + reverse_anneal,
        annealing=reverse_anneal,
        tm=tm(reverse_anneal),
        segments=_segments(
            plan, seqio.revcomp(plan.accepts[1]), seqio.revcomp(three), scheme,
            reverse_tail, reverse_anneal,
        ),
    )
    return forward, reverse


def amplicon(body: str, five: str, three: str, plan: Plan) -> str:
    """The full PCR product: both primer tails with the body between them."""
    return plan.forward_tail(five) + body + seqio.revcomp(plan.reverse_tail(three))


def insert_fragment(body: str, five: str, three: str, plan: Plan) -> str:
    """The top strand the amplicon contributes once the entry enzyme has cut it.

    Rather than counting bases to find the cut, this digests the amplicon and
    takes the piece with the entry-vector ends - so the fragment is right by
    construction, and a mistake in a primer tail shows up here as a missing
    fragment rather than as an off-by-two nobody notices.
    """
    from ..core.enzymes import digest

    product = amplicon(body, five, three, plan)
    for fragment in digest(product, plan.outer, circular=False):
        if fragment.overhangs == plan.accepts:
            return fragment.seq
    raise ValueError(
        f"the primer tails do not present {plan.accepts[0]} / {plan.accepts[1]} ends "
        f"after {plan.outer.name} digestion"
    )


# --------------------------------------------------------------------------- #
# domestication
# --------------------------------------------------------------------------- #


def scanned_enzymes(scheme: Scheme, plan: "Plan | None" = None) -> list[Enzyme]:
    """Every enzyme whose site inside a part would cause trouble somewhere.

    The part enzyme and the cloning enzyme break *this* construction; the
    multigene enzyme breaks the assembly after it; the linearizer breaks
    integration. All four are worth reporting, so all four are scanned.
    """
    enzymes = [scheme.part_enzyme, scheme.multigene_enzyme, *scheme.alternate_cloning_enzymes,
               scheme.linearizer]
    if plan is not None and plan.outer not in enzymes:
        enzymes.append(plan.outer)
    seen, out = set(), []
    for enzyme in enzymes:
        if enzyme.name not in seen:
            seen.add(enzyme.name)
            out.append(enzyme)
    return out


def blocking_enzymes(scheme: Scheme, plan: "Plan | None" = None) -> set[str]:
    """The ones that would stop this construction working at all."""
    names = {scheme.part_enzyme.name}
    if plan is not None:
        names.add(plan.outer.name)
    return names


def internal_sites(body: str, scheme: Scheme, plan: "Plan | None" = None) -> dict[str, list[int]]:
    """Where each scanned enzyme cuts inside the part body."""
    out = {}
    for enzyme in scanned_enzymes(scheme, plan):
        positions = [s.start for s in find_sites(body, enzyme, circular=False)]
        if positions:
            out[enzyme.name] = positions
    return out


def _overhang_is_usable(candidate: str, taken: list[str], scheme: Scheme) -> bool:
    """A junction we choose ourselves must not collide with one we must honour."""
    if len(candidate) != scheme.overhang_length or set(candidate) - set("ACGT"):
        return False
    if seqio.revcomp(candidate) == candidate:
        return False
    if parttypes.conflicts_with_reserved(candidate, scheme):
        return False
    return all(not parttypes.near_match(candidate, other) for other in taken)


def choose_junction(body: str, near: int, taken: list[str], scheme: Scheme) -> int | None:
    """A split point near `near` whose 4 bases make a safe junction overhang."""
    width = scheme.overhang_length
    for offset in range(JUNCTION_SEARCH):
        for position in {near + offset, near - offset}:
            if position < width or position + width > len(body) - width:
                continue
            if _overhang_is_usable(body[position:position + width], taken, scheme):
                return position
    return None


def _cds_features(record: SeqRecord | None) -> list[SeqFeature]:
    return [f for f in (record.features if record else []) if f.type == "CDS"]


def propose_mutation(
    body: str, site_start: int, enzyme: Enzyme, frame_start: int | None, gene: str
) -> Mutation | None:
    """A single base change that kills a site, silent if the frame allows it."""
    site = body[site_start:site_start + len(enzyme.site)]
    for offset in range(len(site)):
        position = site_start + offset
        original = body[position]
        for base in "ACGT":
            if base == original:
                continue
            mutated = body[:position] + base + body[position + 1:]
            if find_sites(mutated[site_start:site_start + len(enzyme.site)], enzyme, circular=False):
                continue

            if frame_start is None:
                return Mutation(position, original, base, enzyme.name, silent=False, gene=gene)

            codon_index = (position - frame_start) // 3
            codon_start = frame_start + codon_index * 3
            if codon_start < 0 or codon_start + 3 > len(body):
                continue
            old_codon = body[codon_start:codon_start + 3]
            new_codon = mutated[codon_start:codon_start + 3]
            if len(old_codon) < 3 or set(old_codon) - set("ACGT"):
                continue
            if _translate(old_codon) == _translate(new_codon):
                return Mutation(position, original, base, enzyme.name, True, gene,
                                old_codon, new_codon)
    return None


def _translate(codon: str) -> str:
    if codon in STANDARD_CODONS.stop_codons:
        return "*"
    return STANDARD_CODONS.forward_table.get(codon, "?")


def domesticate(
    body: str,
    record: SeqRecord | None,
    scheme: Scheme,
    part_type: str,
    plan: "Plan | None" = None,
    apply_changes: bool = False,
) -> tuple[str, list[Mutation], list[Issue]]:
    """Scan for internal sites, propose a fix for each, and optionally apply it.

    The scan always runs: every internal BsaI, BsmBI, BbsI and NotI site is
    reported with the single base change that would remove it, and whether that
    change is silent. Applying them is the caller's decision - see
    `PartRequest.domesticate` for why.

    When changes are applied they are applied *here*, so everything downstream -
    the primers, the predicted plasmid, the validation - works on the sequence
    that will actually be built rather than the one that was pasted in.
    """
    issues: list[Issue] = []
    sites = internal_sites(body, scheme, plan)
    blocking = blocking_enzymes(scheme, plan)
    by_name = {e.name: e for e in scanned_enzymes(scheme, plan)}
    cds = _cds_features(record)

    for enzyme_name, positions in sites.items():
        where = ", ".join(f"{p:,}" for p in positions)
        blocks = enzyme_name in blocking
        issues.append(
            Issue(ERROR if (blocks and not apply_changes) else WARNING,
                  "blocking_site" if blocks else "internal_site",
                  f"{len(positions)} internal {enzyme_name} site(s) at {where}"
                  + (f" - {enzyme_name} releases this part, so the construct cannot "
                     f"be built until they are removed" if blocks and not apply_changes
                     else f" - would break {_what_breaks(enzyme_name, scheme)}"
                     if not blocks else ""))
        )

    mutations: list[Mutation] = []
    for enzyme_name, positions in sites.items():
        enzyme = by_name[enzyme_name]
        for site_start in sorted(positions):
            # positions shift as changes are applied, so re-find the site
            current = [s.start for s in find_sites(body, enzyme, circular=False)]
            here = min(current, key=lambda p: abs(p - site_start), default=None)
            if here is None:
                continue
            frame, gene = _frame_at(here, cds, part_type, len(body))
            mutation = propose_mutation(body, here, enzyme, frame, gene)
            if mutation is None:
                issues.append(
                    Issue(WARNING, "no_fix",
                          f"no single base change removes the {enzyme_name} site at "
                          f"{here:,}; this part needs resynthesis, not PCR")
                )
                continue
            mutations.append(mutation)
            if apply_changes:
                body = body[:mutation.position] + mutation.to_base + body[mutation.position + 1:]
                issues.append(
                    Issue(INFO if mutation.silent else WARNING,
                          "silent_mutation" if mutation.silent else "coding_mutation",
                          mutation.description)
                )
            else:
                issues.append(
                    Issue(INFO, "proposed_mutation",
                          f"could be removed by {mutation.description} - turn on "
                          f"\u201cremove internal sites\u201d to apply it")
                )

    if mutations and not apply_changes:
        silent = sum(1 for m in mutations if m.silent)
        issues.append(
            Issue(INFO, "domestication_available",
                  f"{len(mutations)} site(s) can be removed by a single base change each "
                  f"({silent} of them silent). Left off, because in a promoter or "
                  f"terminator no change is truly silent and the call is yours.")
        )

    return body, (mutations if apply_changes else []), issues


def _what_breaks(enzyme_name: str, scheme: Scheme) -> str:
    if enzyme_name == scheme.multigene_enzyme.name:
        return "a multigene assembly later"
    if enzyme_name == scheme.linearizer.name:
        return f"{scheme.linearizer.name} linearization before integration"
    return "a later step"


def plan_fragments(
    body: str,
    five: str,
    three: str,
    scheme: Scheme,
    mutations: list[Mutation],
) -> tuple[list[FragmentPlan], list[Issue]]:
    """Split the part so that every mutation lands inside a primer.

    A mutation near either end is carried by that end's own primer and needs no
    split. One in the middle does: the amplicon is cut there, and the two
    primers meeting at the junction spell the corrected bases between them.
    """
    issues: list[Issue] = []
    whole = [FragmentPlan(0, 0, len(body), five, three)]
    if not mutations:
        return whole, issues

    needs_split = [
        m for m in mutations
        if MAX_ANNEAL <= m.position < len(body) - MAX_ANNEAL
    ]
    carried = [m for m in mutations if m not in needs_split]
    if carried:
        issues.append(
            Issue(INFO, "carried_by_end_primer",
                  f"{len(carried)} change(s) sit within {MAX_ANNEAL} nt of an end and are "
                  f"spelled by that end's own primer")
        )
    if not needs_split:
        for mutation in mutations:
            whole[0].mutations.append(mutation)
        return whole, issues

    taken = [five, three]
    boundaries: list[int] = []
    for mutation in needs_split:
        position = choose_junction(body, mutation.position, taken, scheme)
        if position is None:
            issues.append(
                Issue(ERROR, "no_junction",
                      f"no safe junction overhang within {JUNCTION_SEARCH} nt of the site at "
                      f"{mutation.position:,}; every nearby 4-mer collides with "
                      f"{', '.join(sorted(set(taken)))}")
            )
            continue
        boundaries.append(position)
        taken.append(body[position:position + scheme.overhang_length])

    boundaries = sorted(set(boundaries))
    if not boundaries:
        return whole, issues

    edges = [0, *boundaries, len(body)]
    fragments = []
    for index, (start, end) in enumerate(zip(edges, edges[1:], strict=False)):
        left = five if index == 0 else body[start:start + scheme.overhang_length]
        right = three if end == len(body) else body[end:end + scheme.overhang_length]
        fragments.append(FragmentPlan(index, start, end, left, right))

    for mutation in mutations:
        owner = next((f for f in fragments if f.start <= mutation.position < f.end), fragments[0])
        owner.mutations.append(mutation)
        if not _covered_by_a_primer(mutation.position, owner):
            issues.append(
                Issue(WARNING, "mutation_unreachable",
                      f"the change at {mutation.position:,} is more than {MAX_ANNEAL} nt from "
                      f"either end of its fragment; order that fragment as a gBlock instead")
            )

    issues.append(
        Issue(INFO, "domesticated",
              f"split into {len(fragments)} fragments at junctions "
              f"{', '.join(f.right_overhang for f in fragments[:-1])}, chosen to avoid "
              f"{'/'.join(scheme.entry_vector_overhangs)} and any 3-of-4 match")
    )
    return fragments, issues


def _covered_by_a_primer(position: int, fragment: FragmentPlan) -> bool:
    return (
        position - fragment.start < MAX_ANNEAL
        or fragment.end - position <= MAX_ANNEAL
    )


def _frame_at(
    position: int, cds: list[SeqFeature], part_type: str = "", body_length: int = 0
) -> tuple[int | None, str]:
    """Where the reading frame containing `position` starts, and whose it is."""
    for feature in cds:
        start, end = int(feature.location.start), int(feature.location.end)
        if start <= position < end:
            label = feature.qualifiers.get("label") or feature.qualifiers.get("gene") or [""]
            return start, str(label[0])
    if part_type in CODING_TYPES and 0 <= position < body_length:
        return 0, f"the type {part_type} body"
    return None, ""


# --------------------------------------------------------------------------- #
# the whole design
# --------------------------------------------------------------------------- #


def choose_destination(
    library: Library, name: str | None, part_type: str, scheme: Scheme
) -> tuple[Destination | None, list[Issue]]:
    """Pick the plasmid to clone into, and say why if it will not do.

    Preference, when the user has not chosen: a dropout vector that takes this
    exact part type, then any universal entry vector, then an existing part
    plasmid of the same type - because dropping into a ccdB dropout is a cleaner
    selection than re-opening a working part plasmid.
    """
    issues: list[Issue] = []
    wanted = parttypes.type_overhangs(part_type, scheme)
    candidates = part_destinations(library, scheme)

    if name:
        chosen = [d for d in candidates if name in (d.name, d.path, *d.aliases)]
        if not chosen:
            return None, [Issue(ERROR, "unknown_destination",
                                f"no plasmid named {name} can take a part")]
        fits = [d for d in chosen if not d.supplies_part_sites or d.accepts == wanted]
        if not fits:
            return None, [
                Issue(ERROR, "destination_mismatch",
                      f"{name} takes {chosen[0].accepts[0]} \u2192 {chosen[0].accepts[1]} "
                      f"inserts, which is not a type {part_type} part "
                      f"({wanted[0]} \u2192 {wanted[1]})" if wanted else
                      f"{name} cannot take a type {part_type} part")
            ]
        return _best(fits), issues

    # A universal vector takes any part type, because the insert carries its own
    # part-enzyme sites and the vector's own ends never change - that is the
    # L131 system, and it is the default. A type-specific dropout needs a much
    # shorter primer, so it is offered next, then reopening a part plasmid.
    universal = [d for d in candidates if not d.supplies_part_sites]
    matching = [d for d in candidates if d.supplies_part_sites and d.accepts == wanted]
    dropout_first = [d for d in matching if d.kind == "dropout"]

    for group in (universal, dropout_first, matching):
        if group:
            return _best(group), issues
    return None, [
        Issue(ERROR, "no_destination",
              f"nothing in the library can take a type {part_type} part: it needs a "
              f"dropout vector or an existing type {part_type} part plasmid")
    ]


def part_destinations(library: Library, scheme: Scheme) -> list[Destination]:
    """Destinations that can hold a *part*, as opposed to a cassette.

    A plasmid with two cuts is not automatically somewhere a part can go: the
    ends it offers have to be a part type's own overhangs, or the fixed ends of
    a universal entry vector. Connector ends belong to Level 3.
    """
    out = []
    for destination in library.destinations():
        if destination.part_type is not None:
            out.append(destination)
        elif destination.accepts == scheme.entry_vector_overhangs:
            out.append(destination)
    return out


def _best(candidates: list[Destination]) -> Destination:
    """Prefer a ccdB dropout cut with the kit's alternate enzyme - the L131
    system - then any ccdB dropout, then the smallest backbone."""
    def rank(d: Destination) -> tuple:
        return (
            d.kind != "dropout",
            "ccdb" not in d.name.lower(),
            d.cloning_enzyme == "BsmBI",
            d.backbone_length,
            d.name,
        )
    return min(candidates, key=rank)


def make_plan(destination: Destination, scheme: Scheme) -> Plan:
    """How this destination wants its insert built."""
    outer = get_enzyme(destination.cloning_enzyme)
    inner = None if destination.supplies_part_sites else scheme.part_enzyme
    return Plan(destination=destination, outer=outer, inner=inner)


def design(library: Library, request: PartRequest) -> Level1Design:
    """Design the primers, scan for internal sites, and predict the plasmid."""
    scheme = library.scheme
    issues: list[Issue] = []

    pair = parttypes.type_overhangs(request.part_type, scheme)
    if pair is None:
        return Level1Design(
            part_type=request.part_type, five_prime="", three_prime="", body="",
            issues=[Issue(ERROR, "unknown_type", f"{request.part_type} is not a part type")],
        )
    five, three = pair

    record, source, error = _template(library, request)
    if error:
        return Level1Design(request.part_type, five, three, "", issues=[error])

    body = _slice(source, request)
    if len(body) < MIN_ANNEAL * 2:
        issues.append(
            Issue(ERROR, "too_short",
                  f"the chosen range is {len(body)} bp; a PCR product needs at least "
                  f"{MIN_ANNEAL * 2} bp for two primers")
        )
        return Level1Design(request.part_type, five, three, body, issues=issues)

    if request.trim_to_insert:
        body, trim_issues = trim_to_insert(body, request.part_type, scheme)
        issues += trim_issues

    core, prefix, suffix, convention_issues = apply_conventions(
        body, request.part_type, request.conventions
    )
    issues += convention_issues
    body = prefix + core + suffix

    destination, destination_issues = choose_destination(
        library, request.destination or request.entry_vector, request.part_type, scheme
    )
    issues += destination_issues
    if destination is None:
        return Level1Design(request.part_type, five, three, body, issues=issues)
    plan = make_plan(destination, scheme)

    body, mutations, domestication_issues = domesticate(
        body, record, scheme, request.part_type, plan, apply_changes=request.domesticate
    )
    issues += domestication_issues

    fragments, split_issues = plan_fragments(body, five, three, scheme, mutations)
    issues += split_issues

    for fragment in fragments:
        piece = body[fragment.start:fragment.end]
        # the convention bases sit at the very ends of the body, so they belong
        # to the first and last fragments' tails and to nothing else
        own_prefix = prefix if fragment.start == 0 else ""
        own_suffix = suffix if fragment.end == len(body) else ""
        piece = piece[len(own_prefix):len(piece) - len(own_suffix) or None]
        if len(piece) >= MIN_ANNEAL:
            forward, reverse = design_primers(
                piece, fragment.left_overhang, fragment.right_overhang, plan, scheme,
                f"{request.name}_{fragment.index + 1}" if len(fragments) > 1 else request.name,
                prefix=own_prefix, suffix=own_suffix,
            )
            fragment.forward, fragment.reverse = forward, reverse

    spread = [p.tm for p in _all_primers(fragments)]
    if spread and max(spread) - min(spread) > TM_SPREAD:
        issues.append(
            Issue(WARNING, "tm_spread",
                  f"primer Tm spans {min(spread):.1f}-{max(spread):.1f} \u00b0C; "
                  f"consider one annealing temperature per fragment")
        )

    design_obj = Level1Design(
        part_type=request.part_type,
        five_prime=five,
        three_prime=three,
        body=body,
        fragments=fragments,
        issues=issues,
        entry_vector=destination.name,
        destination=destination.name,
        cloning_enzyme=destination.cloning_enzyme,
        accepts=destination.accepts,
        destination_marker=destination.marker,
        destination_universal=not destination.supplies_part_sites,
        mode=request.mode,
    )

    _build_product(library, plan, design_obj, request, scheme)
    _alternate_formats(design_obj, plan, scheme)
    design_obj.ligation = ligation_view(library, plan, design_obj, scheme)
    design_obj.binding = binding_view(design_obj)
    return design_obj


def _all_primers(fragments: list[FragmentPlan]) -> list[Primer]:
    return [p for f in fragments for p in (f.forward, f.reverse) if p]


def _template(
    library: Library, request: PartRequest
) -> tuple[SeqRecord | None, str, Issue | None]:
    if request.template:
        entry = library.get(request.template)
        if entry is None:
            return None, "", Issue(ERROR, "unknown_template",
                                   f"no plasmid named {request.template}")
        record = library.record(entry)
        return record, seqio.sequence(record), None
    if request.sequence:
        cleaned = "".join(c for c in request.sequence.upper() if c in "ACGT")
        if not cleaned:
            return None, "", Issue(ERROR, "empty_sequence", "the pasted sequence has no DNA in it")
        return None, cleaned, None
    return None, "", Issue(ERROR, "no_template", "choose a template or paste a sequence")


def _slice(source: str, request: PartRequest) -> str:
    start = max(0, request.start - 1)
    end = request.end if request.end else len(source)
    end = min(end, len(source))
    if start < end:
        return source[start:end]
    return seqio.subseq(source, start, end, circular=True)  # a range across the origin


def _build_product(
    library: Library,
    plan: Plan,
    design_obj: Level1Design,
    request: PartRequest,
    scheme: Scheme,
) -> None:
    """Ligate the insert into the destination and check what comes out."""
    try:
        insert = insert_fragment(
            design_obj.body, design_obj.five_prime, design_obj.three_prime, plan
        )
        backbone = _destination_backbone(library, plan, scheme)
    except ValueError as exc:
        design_obj.issues.append(Issue(ERROR, "bad_destination", str(exc)))
        return

    insert_piece = Piece(
        record=SeqRecord(Seq(insert), id=request.name, annotations={"molecule_type": "DNA"}),
        left_overhang=plan.accepts[0],
        right_overhang=plan.accepts[1],
        source_name=request.name,
        part_type=design_obj.part_type,
        color=parttypes.color_for(design_obj.part_type, scheme),
        component=request.name,
    )

    result: AssemblyResult = assemble(
        [backbone, insert_piece], plan.outer, name=request.name
    )
    design_obj.issues += [i for i in result.issues if i.code != "internal_site"]

    if not result.ok or result.product is None:
        return

    product = result.product
    product.features.append(
        SeqFeature(
            SimpleLocation(0, len(product.seq), 1),
            type="misc_feature",
            qualifiers={"label": [f"type {design_obj.part_type} part"],
                        "note": [f"designed by GG Assembler from {request.template or 'a pasted sequence'}"]},
        )
    )
    design_obj.product = product

    # the promise of the whole screen: the thing we predict must read back as
    # the type that was asked for
    call = describe(product, scheme=scheme).call
    design_obj.validated_as = call.part_type
    if call.part_type != design_obj.part_type:
        design_obj.issues.append(
            Issue(ERROR, "validation_failed",
                  f"the predicted plasmid reads back as {call.part_type or 'nothing'} "
                  f"({call.reason}), not {design_obj.part_type}")
        )
    else:
        design_obj.issues.append(
            Issue(INFO, "validated",
                  f"the predicted plasmid re-digests as a type {call.part_type} part: "
                  f"{call.five_prime} → {call.three_prime}")
        )


def _destination_backbone(library: Library, plan: Plan, scheme: Scheme) -> Piece:
    """The destination minus whatever currently sits in it, ready for the insert.

    For a dropout vector that is the piece without the cloning enzyme's sites;
    for a part plasmid being reopened it is the piece *with* them, since there
    the sites live in the backbone and the old part is what leaves.
    """
    from ..core.enzymes import digest

    destination = plan.destination
    entry = library.get(destination.name)
    if entry is None:
        raise ValueError(f"{destination.name} is no longer in the library")
    record = library.record(entry)
    wanted = (destination.accepts[1], destination.accepts[0])

    for fragment in digest(seqio.sequence(record), plan.outer):
        keep = fragment.has_sites if destination.kind == "part_plasmid" else not fragment.has_sites
        if keep and fragment.overhangs == wanted:
            sub = seqio.slice_record(record, fragment.start, fragment.end, name=destination.name)
            return Piece(
                record=sub,
                left_overhang=fragment.left_overhang or "",
                right_overhang=fragment.right_overhang or "",
                source_name=destination.name,
                part_type="backbone",
                color="#8C8A86",
                component=destination.component or "destination backbone",
            )
    raise ValueError(
        f"{destination.name} does not leave a backbone with {wanted[0]} / {wanted[1]} "
        f"ends when cut with {plan.outer.name}"
    )


def _alternate_formats(design_obj: Level1Design, plan: Plan, scheme: Scheme) -> None:
    """The two non-PCR ways to get the same insert."""
    try:
        insert = insert_fragment(
            design_obj.body, design_obj.five_prime, design_obj.three_prime, plan
        )
    except ValueError:
        return
    design_obj.gblock = amplicon(
        design_obj.body, design_obj.five_prime, design_obj.three_prime, plan
    )

    entry = plan.accepts
    top = insert
    bottom = seqio.revcomp(insert[len(entry[0]):] + entry[1])
    design_obj.oligos = [
        Primer(name=f"{design_obj.part_type}_top", sequence=top, annealing=top, tm=tm(top[:25])),
        Primer(name=f"{design_obj.part_type}_bottom", sequence=bottom, annealing=bottom,
               tm=tm(bottom[:25])),
    ]


FLANK = 12
"""How many bases either side of a junction the ligation diagram shows."""


def ligation_view(
    library: Library, plan: Plan, design_obj: Level1Design, scheme: Scheme
) -> dict[str, object]:
    """The three ends that have to meet, with the bases around each one.

    Enough for a drawing of the reaction: where the primers sit on the template,
    what the finished amplicon is made of, and how its sticky ends pair with the
    ones the destination is left with once its dropout is cut out.
    """
    from ..core.enzymes import digest

    entry = library.get(plan.destination.name)
    if entry is None or not design_obj.body:
        return {}
    record = library.record(entry)
    seq = seqio.sequence(record)
    length = len(seq)

    wanted = (plan.accepts[1], plan.accepts[0])
    backbone = dropout = None
    for fragment in digest(seq, plan.outer):
        keep = (
            fragment.has_sites
            if plan.destination.kind == "part_plasmid"
            else not fragment.has_sites
        )
        if keep and fragment.overhangs == wanted:
            backbone = fragment
        else:
            dropout = fragment
    if backbone is None:
        return {}

    try:
        insert = insert_fragment(
            design_obj.body, design_obj.five_prime, design_obj.three_prime, plan
        )
    except ValueError:
        return {}

    width = scheme.part_enzyme.overhang_length
    return {
        "vector": {
            "name": plan.destination.name,
            "enzyme": plan.outer.name,
            "kind": plan.destination.kind,
            "length": entry.length,
            "backbone_length": len(backbone.seq),
            "dropout_length": len(dropout.seq) if dropout else 0,
            "dropout_component": plan.destination.component,
            # the backbone's own ends: it starts with the overhang the insert's
            # 3' end will pair with, and finishes just before the 5' one
            "left_overhang": backbone.left_overhang,
            "right_overhang": backbone.right_overhang,
            "left_flank": backbone.seq[width:width + FLANK],
            "right_flank": backbone.seq[-FLANK:] if len(backbone.seq) > FLANK else backbone.seq,
        },
        "insert": {
            "name": design_obj.part_type,
            "length": len(insert),
            "left_overhang": plan.accepts[0],
            "right_overhang": plan.accepts[1],
            "left_flank": insert[width:width + FLANK],
            "right_flank": insert[-FLANK:] if len(insert) > FLANK else insert,
        },
    }


def binding_view(design_obj: Level1Design) -> list[dict[str, object]]:
    """Where each primer sits on the template, and how much of it hangs off.

    A tail does not anneal on the first cycle - it is added by it - so the two
    are drawn differently, and this reports them separately rather than as one
    primer sequence.
    """
    out = []
    for fragment in design_obj.fragments:
        for primer, direction in ((fragment.forward, "forward"), (fragment.reverse, "reverse")):
            if primer is None:
                continue
            tail = primer.sequence[: len(primer.sequence) - len(primer.annealing)]
            out.append(
                {
                    "name": primer.name,
                    "direction": direction,
                    "tail": tail,
                    "annealing": primer.annealing,
                    "tm": primer.tm,
                    "length": primer.length,
                    "fragment": fragment.index,
                    "segments": [
                        {"role": segment.role, "seq": segment.seq}
                        for segment in primer.segments
                    ],
                }
            )
    return out


def save(library: Library, design_obj: Level1Design, name: str) -> str | None:
    """Write the predicted part plasmid into the library folder."""
    if design_obj.product is None:
        return None
    path = library.folder / f"{name}.gb"
    seqio.write_genbank(design_obj.product, path)
    library.scan()
    return str(path)
