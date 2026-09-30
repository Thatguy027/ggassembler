"""Type IIS site finding and overhang extraction (BsaI, BsmBI, NotI).

Cut geometry
------------
A Type IIS enzyme cuts at a fixed offset from its recognition site, leaving a
4-nt 5' overhang. For a site found on the top strand starting at index ``i``,
the top strand breaks before ``i + top_offset`` and the bottom strand before
``i + bottom_offset``. The four bases between those two cuts are the overhang;
by convention throughout this package an overhang is always written as the
**top-strand** sequence, 5'->3', which is how the YTK table in `parttypes` is
written.

Both fragments meeting at a cut carry that overhang - one on the top strand,
one on the bottom - which is exactly why they re-anneal. So `Fragment` records
its left overhang as the first four bases of its own top strand, and its right
overhang as the first four bases of the *next* fragment. Concatenating the
top strands of all fragments in order reconstitutes the molecule.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

from .seqio import revcomp, span_length, subseq, wrap

IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "[AG]", "Y": "[CT]", "S": "[GC]", "W": "[AT]",
    "K": "[GT]", "M": "[AC]", "B": "[CGT]", "D": "[AGT]",
    "H": "[ACT]", "V": "[ACG]", "N": ".",
}


@dataclass(frozen=True)
class Enzyme:
    """A restriction enzyme, described by where it cuts relative to its site."""

    name: str
    site: str
    top_offset: int
    bottom_offset: int

    @property
    def overhang_length(self) -> int:
        return self.bottom_offset - self.top_offset

    @property
    def is_palindromic(self) -> bool:
        return revcomp(self.site) == self.site

    def pattern(self, strand: str = "+") -> re.Pattern[str]:
        site = self.site if strand == "+" else revcomp(self.site)
        return re.compile("".join(IUPAC[b] for b in site))


#: BsaI, GGTCTC(1/5) - the Level 1 / Level 2 enzyme.
BSAI = Enzyme("BsaI", "GGTCTC", top_offset=7, bottom_offset=11)
#: BsmBI (Esp3I), CGTCTC(1/5) - the entry-vector and Level 3 enzyme.
BSMBI = Enzyme("BsmBI", "CGTCTC", top_offset=7, bottom_offset=11)
#: BbsI (BpiI), GAAGAC(2/6) - a third Type IIS cutter some destination vectors
#: use in place of BsmBI, so a part with an internal BsmBI site can still be
#: cloned without domesticating it first.
BBSI = Enzyme("BbsI", "GAAGAC", top_offset=8, bottom_offset=12)
#: NotI, GC^GGCCGC - linearizes integration vectors.
NOTI = Enzyme("NotI", "GCGGCCGC", top_offset=2, bottom_offset=6)

#: XhoI, C^TCGAG. Not a cloning enzyme here - the kit *writes* this site, at the
#: start of every Type 4 part, for BglBrick compatibility. Worth being able to
#: keep out of a synthetic CDS for exactly that reason: an extra copy inside the
#: gene makes the intended one useless.
XHOI = Enzyme("XhoI", "CTCGAG", top_offset=1, bottom_offset=5)
#: BamHI, G^GATCC. The kit writes this one too, as the GGATCC that the Gly-Ser
#: linker spells at the Type 3 / Type 4 junction.
BAMHI = Enzyme("BamHI", "GGATCC", top_offset=1, bottom_offset=5)

ENZYMES: dict[str, Enzyme] = {
    e.name.lower(): e for e in (BSAI, BSMBI, BBSI, NOTI, XHOI, BAMHI)
}


def get_enzyme(name: str) -> Enzyme:
    try:
        return ENZYMES[name.lower()]
    except KeyError:
        raise ValueError(f"unknown enzyme: {name}") from None


@dataclass(frozen=True)
class Site:
    """One recognition site and the cut it makes."""

    enzyme: str
    start: int
    """0-based index of the recognition sequence on the top strand."""
    end: int
    """Exclusive end of the recognition sequence (may wrap past the origin)."""
    strand: str
    """``+`` when the site cuts downstream of itself, ``-`` when upstream."""
    top_cut: int
    bottom_cut: int
    overhang: str
    """The 4-nt top-strand overhang, 5'->3'."""

    @property
    def cuts_downstream(self) -> bool:
        return self.strand == "+"


@dataclass
class Fragment:
    """A piece of DNA released by a digest."""

    start: int
    """Top-strand cut position where this fragment begins."""
    end: int
    """Top-strand cut position where the next fragment begins."""
    seq: str
    """Top strand from `start` (inclusive) to `end` (exclusive)."""
    left_overhang: str | None
    right_overhang: str | None
    sites: list[Site] = field(default_factory=list)
    """Recognition sites lying wholly inside this fragment."""
    source: str = ""

    def __len__(self) -> int:
        return len(self.seq)

    @property
    def has_sites(self) -> bool:
        return bool(self.sites)

    @property
    def overhangs(self) -> tuple[str | None, str | None]:
        return self.left_overhang, self.right_overhang


def find_sites(seq: str, enzyme: Enzyme, circular: bool = True) -> list[Site]:
    """Every site for `enzyme`, on both strands, sorted by position.

    On a circular molecule sites spanning the origin are found by searching the
    sequence concatenated with itself and mapping positions back.
    """
    seq = seq.upper()
    n = len(seq)
    if n == 0:
        return []
    haystack = seq + seq[: len(enzyme.site) - 1] if circular else seq

    found: dict[tuple[int, str], Site] = {}
    for strand in ("+", "-"):
        if strand == "-" and enzyme.is_palindromic:
            continue
        for match in enzyme.pattern(strand).finditer(haystack):
            i = match.start()
            if i >= n:
                continue
            if strand == "+":
                top = i + enzyme.top_offset
                bottom = i + enzyme.bottom_offset
            else:
                length = len(enzyme.site)
                bottom = i + length - enzyme.top_offset
                top = i + length - enzyme.bottom_offset
            if not circular and (top < 0 or bottom > n):
                continue  # cut falls off the end of a linear molecule
            top, bottom = (wrap(top, n), wrap(bottom, n)) if circular else (top, bottom)
            site = Site(
                enzyme=enzyme.name,
                start=i,
                end=i + len(enzyme.site),
                strand=strand,
                top_cut=top,
                bottom_cut=bottom,
                overhang=subseq(seq, top, top + enzyme.overhang_length, circular),
            )
            found[(i, strand)] = site

    return sorted(found.values(), key=lambda s: (s.start, s.strand))


def count_sites(seq: str, enzyme: Enzyme, circular: bool = True) -> int:
    return len(find_sites(seq, enzyme, circular))


def _site_inside(site: Site, start: int, end: int, n: int) -> bool:
    """True when a site's whole recognition sequence lies in the span."""
    offset = (site.start - start) % n
    return offset + (site.end - site.start) <= span_length(start, end, n)


def digest(
    seq: str,
    enzyme: Enzyme,
    circular: bool = True,
    source: str = "",
) -> list[Fragment]:
    """Cut `seq` with `enzyme` and return the fragments in top-strand order.

    A circular molecule with no sites comes back as one uncut fragment with no
    overhangs. Two sites give two fragments, and so on.
    """
    seq = seq.upper()
    n = len(seq)
    sites = find_sites(seq, enzyme, circular)

    if not circular:
        return _digest_linear(seq, sites, source)

    if not sites:
        return [Fragment(0, 0, seq, None, None, [], source)]

    cuts = sorted({s.top_cut: s for s in sites}.items())
    fragments = []
    for idx, (cut, site) in enumerate(cuts):
        next_cut, next_site = cuts[(idx + 1) % len(cuts)]
        frag = Fragment(
            start=cut,
            end=next_cut,
            seq=subseq(seq, cut, next_cut, circular=True),
            left_overhang=site.overhang,
            right_overhang=next_site.overhang,
            sites=[s for s in sites if _site_inside(s, cut, next_cut, n)],
            source=source,
        )
        fragments.append(frag)
    return fragments


def _digest_linear(seq: str, sites: list[Site], source: str) -> list[Fragment]:
    n = len(seq)
    cuts = sorted({s.top_cut: s for s in sites}.items())
    bounds = [(0, None), *[(c, s) for c, s in cuts], (n, None)]
    fragments = []
    for (start, left), (end, right) in zip(bounds, bounds[1:], strict=False):
        if start >= end:
            continue
        fragments.append(
            Fragment(
                start=start,
                end=end,
                seq=seq[start:end],
                left_overhang=left.overhang if left else None,
                right_overhang=right.overhang if right else None,
                sites=[s for s in sites if start <= s.start and s.end <= end],
                source=source,
            )
        )
    return fragments


def overhang_pair(fragment: Fragment) -> tuple[str, str] | None:
    """The fragment's (5', 3') overhangs, or None if either end is blunt."""
    left, right = fragment.overhangs
    if left is None or right is None:
        return None
    return left, right


def is_palindromic(overhang: str) -> bool:
    return revcomp(overhang) == overhang


def mismatch_positions(a: str, b: str) -> list[int]:
    """Where two equal-length overhangs differ; empty if they are identical."""
    if len(a) != len(b):
        return []
    return [i for i, (x, y) in enumerate(zip(a, b, strict=True)) if x != y]


def shares_three_of_four(a: str, b: str) -> bool:
    """True when two overhangs differ at exactly one position.

    The YTK SI names this as a misligation mode the overhang set was chosen to
    avoid, so it is a warning wherever two overhangs meet in one reaction.
    """
    return a != b and len(mismatch_positions(a, b)) == 1


def shares_noncontiguous_three(a: str, b: str) -> bool:
    """True when the three shared bases are not contiguous.

    For two 4-nt overhangs, sharing three bases *is* differing at one position,
    so this is not a separate condition but a sharper one: it asks whether the
    odd base out sits at an end or in the middle. ``ATCG``/``ATAG`` differ at
    position 2, leaving the matches split either side of it - the worse case,
    because a run broken in the middle still anneals across the whole overhang.
    A mismatch at either end leaves a contiguous run of three instead.
    """
    positions = mismatch_positions(a, b)
    if a == b or len(positions) != 1:
        return False
    return 0 < positions[0] < len(a) - 1
