"""Every row of the YTK table round-trips through `enzymes.py`.

Build a synthetic part plasmid for the row, digest it, and check that the
released fragment's overhangs are the pair the table promises.
"""

from __future__ import annotations

import pytest

from ggassembler.core import parttypes
from ggassembler.core.enzymes import BSAI, BSMBI, NOTI, digest, find_sites
from ggassembler.core.parttypes import span, type_overhangs
from ggassembler.core.seqio import revcomp

from . import synth

# The table exactly as printed in the paper's SI, transcribed independently of
# the BOUNDARIES tuple the code derives everything else from.
TABLE = [
    ("1", "CCCT", "AACG"),
    ("2", "AACG", "TATG"),
    ("3", "TATG", "ATCC"),
    ("3a", "TATG", "TTCT"),
    ("3b", "TTCT", "ATCC"),
    ("4", "ATCC", "GCTG"),
    ("4a", "ATCC", "TGGC"),
    ("4b", "TGGC", "GCTG"),
    ("5", "GCTG", "TACA"),
    ("6", "TACA", "GAGT"),
    ("7", "GAGT", "CCGA"),
    ("8", "CCGA", "CCCT"),
    ("8a", "CCGA", "CAAT"),
    ("8b", "CAAT", "CCCT"),
    ("234", "AACG", "GCTG"),
    ("678", "TACA", "CCCT"),
]


@pytest.mark.parametrize(("name", "five", "three"), TABLE)
def test_table_overhangs_match_the_paper(name, five, three):
    assert type_overhangs(name) == (five, three)


@pytest.mark.parametrize(("name", "five", "three"), TABLE)
def test_part_plasmid_releases_the_table_pair(name, five, three):
    record = synth.part_plasmid(name)
    fragments = digest(str(record.seq), BSAI)
    assert len(fragments) == 2

    released = [f for f in fragments if not f.has_sites]
    assert len(released) == 1, "the part fragment must be free of BsaI sites"
    assert released[0].overhangs == (five, three)


@pytest.mark.parametrize(("name", "five", "three"), TABLE)
def test_overhang_pair_resolves_back_to_the_type(name, five, three):
    call = span(five, three)
    assert call is not None and call.legal
    assert call.name == name


@pytest.mark.parametrize(("name", "five", "three"), TABLE)
def test_fragments_reassemble_into_the_original(name, five, three):
    record = synth.part_plasmid(name)
    seq = str(record.seq)
    fragments = digest(seq, BSAI)
    joined = "".join(f.seq for f in fragments)
    assert len(joined) == len(seq)
    # the digest starts at the first cut, so the join is a rotation of the input
    assert joined in seq + seq


def test_backbone_overhangs_are_an_illegal_span():
    """The leftover backbone is an arc that crosses the 8/1 boundary."""
    record = synth.part_plasmid("234")
    backbone = [f for f in digest(str(record.seq), BSAI) if f.has_sites][0]
    call = span(*backbone.overhangs)
    assert call is not None and not call.legal


def test_dropout_keeps_the_sites_and_the_span():
    """A 234r dropout has the 234 pair on the fragment that carries the sites."""
    record = synth.dropout_plasmid("234")
    fragments = digest(str(record.seq), BSAI)
    with_sites = [f for f in fragments if f.has_sites]
    assert len(with_sites) == 1
    assert with_sites[0].overhangs == type_overhangs("234")
    assert len(with_sites[0].sites) == 2

    site_free = [f for f in fragments if not f.has_sites][0]
    assert not span(*site_free.overhangs).legal


def test_entry_vector_overhangs():
    record = synth.entry_vector()
    fragments = digest(str(record.seq), BSMBI)
    site_free = [f for f in fragments if not f.has_sites][0]
    assert site_free.overhangs == tuple(reversed(parttypes.ENTRY_VECTOR_OVERHANGS))


def test_site_finding_wraps_the_origin():
    """A site split across the origin of a circular plasmid is still found."""
    record = synth.part_plasmid("2")
    seq = str(record.seq)
    rotated = seq[3:] + seq[:3]  # cut the leading BsaI site in half
    assert len(find_sites(rotated, BSAI)) == 2
    fragments = digest(rotated, BSAI)
    released = [f for f in fragments if not f.has_sites][0]
    assert released.overhangs == type_overhangs("2")


def test_reverse_site_cut_geometry():
    """A reverse site cuts upstream of itself, leaving the same top-strand overhang."""
    seq = "TTTT" + "ACGT" + "A" + revcomp(BSAI.site) + "TTTT"
    site = find_sites(seq, BSAI, circular=False)[0]
    assert site.strand == "-"
    assert site.overhang == "ACGT"
    assert site.bottom_cut - site.top_cut == 4


def test_noti_is_palindromic_and_cuts_inside_its_site():
    seq = "AAAA" + NOTI.site + "TTTT"
    sites = find_sites(seq, NOTI, circular=False)
    assert len(sites) == 1, "a palindromic site must not be reported twice"
    assert sites[0].overhang == "GGCC"


def test_uncut_plasmid_is_one_blunt_fragment():
    record = synth.plasmid_without_bsai()
    fragments = digest(str(record.seq), BSAI)
    assert len(fragments) == 1
    assert fragments[0].overhangs == (None, None)


def test_misligation_modes():
    from ggassembler.core.enzymes import shares_noncontiguous_three, shares_three_of_four

    assert shares_three_of_four("ATCG", "ATCA")
    assert not shares_three_of_four("ATCG", "ATCG")
    assert not shares_three_of_four("ATCG", "ATGA")

    assert shares_noncontiguous_three("ATCG", "ATAG")
    assert not shares_noncontiguous_three("ATCG", "ATGA")
    assert not shares_noncontiguous_three("ATCG", "ATCA")


def test_no_ytk_overhang_is_palindromic():
    from ggassembler.core.enzymes import is_palindromic

    for overhang in set(parttypes.BOUNDARIES) | set(parttypes.ENTRY_VECTOR_OVERHANGS):
        assert not is_palindromic(overhang), overhang
