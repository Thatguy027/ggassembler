"""Checking a sequenced clone against the plasmid it was meant to be.

Every case here is built by taking a real plasmid and doing to it exactly what
a vendor's assembler does - opening the circle somewhere else, handing it back
on the other strand - plus, separately, the mutations a clone actually carries.
The two must never be confused, and that is what most of this file is about.
"""

from __future__ import annotations

import random

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from ggassembler.core import align, seqio


def rotate(sequence: str, by: int) -> str:
    return sequence[by:] + sequence[:by]


@pytest.fixture(scope="module")
def reference() -> str:
    """A sequence with no long internal repeat, so anchors are unambiguous."""
    rng = random.Random(20250930)
    return "".join(rng.choice("ACGT") for _ in range(4000))


# --------------------------------------------------------------------------- #
# placing the clone
# --------------------------------------------------------------------------- #


def test_a_rotated_clone_is_not_a_difference(reference):
    """The single most important case: where the vendor opened the circle."""
    result = align.compare(reference, rotate(reference, 1500), name="c1")
    assert result.placement.placed
    assert result.placement.strand == "+"
    assert result.differences == []
    assert result.identity(reference) == 100.0


def test_a_reverse_complemented_clone_is_not_a_difference(reference):
    result = align.compare(reference, align.revcomp(rotate(reference, 900)), name="c1")
    assert result.placement.strand == "-"
    assert result.differences == []
    assert result.clean


@pytest.mark.parametrize("by", [0, 1, 7, 137, 1999, 3999])
def test_every_rotation_places_cleanly(reference, by):
    """A seam anywhere must not invent an indel where the circle was cut.

    Rotating by the voted diagonal alone puts the join at an arbitrary base; if
    the estimate is a few bases out, the alignment opens with a gap no clone
    carries. Cutting inside a shared exact match is what prevents that, and
    this is the test that would catch its removal.
    """
    assert align.compare(reference, rotate(reference, by)).differences == []


def test_an_unrelated_sequence_is_reported_as_unplaced(reference):
    """Better to say "this is not that plasmid" than to align it into noise."""
    rng = random.Random(11)
    other = "".join(rng.choice("ACGT") for _ in range(3000))
    result = align.compare(reference, other, name="wrong tube")
    assert not result.placement.placed
    assert not result.clean
    assert "does not match" in result.verdict()
    assert result.differences == []


# --------------------------------------------------------------------------- #
# the differences themselves
# --------------------------------------------------------------------------- #


def test_a_substitution_is_found_at_the_right_base(reference):
    mutant = list(reference)
    mutant[2000] = "A" if reference[2000] != "A" else "C"
    result = align.compare(reference, "".join(mutant))
    assert len(result.differences) == 1
    difference = result.differences[0]
    assert difference.kind == "substitution"
    assert difference.start == 2000
    assert difference.expected == reference[2000]
    assert difference.found == mutant[2000]
    assert difference.where == "2,001"


def test_a_substitution_survives_rotation_and_flipping(reference):
    """The coordinate reported is the reference's, whatever the clone's was."""
    mutant = list(reference)
    mutant[2000] = "A" if reference[2000] != "A" else "C"
    clone = align.revcomp(rotate("".join(mutant), 3333))
    result = align.compare(reference, clone)
    assert [(d.kind, d.start) for d in result.differences] == [("substitution", 2000)]


def test_adjacent_substitutions_are_one_event(reference):
    """Six changed bases in a row is one thing to look at, not six."""
    mutant = list(reference)
    for i in range(1000, 1006):
        mutant[i] = "A" if reference[i] != "A" else "C"
    result = align.compare(reference, "".join(mutant))
    assert len(result.differences) == 1
    assert result.differences[0].length == 6
    assert result.differences[0].where == "1,001-1,006"


def test_a_deletion_is_found_and_measured(reference):
    clone = reference[:1200] + reference[1212:]
    result = align.compare(reference, clone)
    deletions = [d for d in result.differences if d.kind == "deletion"]
    assert len(deletions) == 1
    assert deletions[0].length == 12
    assert deletions[0].start == 1200


def test_an_insertion_is_found_and_measured(reference):
    clone = reference[:2500] + "TTTTGGGGAAAA" + reference[2500:]
    result = align.compare(reference, clone)
    insertions = [d for d in result.differences if d.kind == "insertion"]
    assert len(insertions) == 1
    assert insertions[0].found == "TTTTGGGGAAAA"
    assert insertions[0].start == 2500


def test_a_large_dropout_comes_back_as_one_deletion(reference):
    """The case that matters: a clone that lost a whole part."""
    clone = reference[:800] + reference[1900:]
    result = align.compare(reference, rotate(clone, 400))
    deletions = [d for d in result.differences if d.kind == "deletion"]
    assert len(deletions) == 1
    assert deletions[0].length == 1100


def test_several_differences_in_one_clone(reference):
    mutant = list(reference)
    mutant[300] = "A" if reference[300] != "A" else "G"
    clone = "".join(mutant)
    clone = clone[:1500] + clone[1510:]
    clone = clone[:2500] + "CCCCC" + clone[2500:]
    result = align.compare(reference, clone)
    kinds = sorted(d.kind for d in result.differences)
    assert kinds == ["deletion", "insertion", "substitution"]
    assert result.verdict().startswith("3 differences")
    assert not result.clean


# --------------------------------------------------------------------------- #
# naming the place
# --------------------------------------------------------------------------- #


def test_a_difference_is_named_by_the_feature_it_lands_in(reference):
    regions = [
        align.Region(0, 500, "pTDH3", "promoter"),
        align.Region(500, 2000, "XI", "CDS"),
        align.Region(2000, 2300, "tCYC1", "terminator"),
    ]
    mutant = list(reference)
    mutant[1000] = "A" if reference[1000] != "A" else "G"
    result = align.compare(reference, "".join(mutant), regions=regions)
    assert result.differences[0].regions == ("XI",)
    assert "in XI" in result.differences[0].describe()


def test_a_difference_outside_every_feature_is_named_by_nothing(reference):
    regions = [align.Region(0, 500, "pTDH3", "promoter")]
    mutant = list(reference)
    mutant[3000] = "A" if reference[3000] != "A" else "G"
    result = align.compare(reference, "".join(mutant), regions=regions)
    assert result.differences[0].regions == ()


def test_a_feature_crossing_the_origin_still_covers_its_own_bases():
    """A circular plasmid's features wrap, and `end` past the length says so."""
    region = align.Region(3900, 4100, "wrapped")
    assert region.covers(50, 51, 4000)
    assert region.covers(3950, 3951, 4000)
    assert not region.covers(2000, 2001, 4000)


def test_regions_come_from_a_record_without_deciding_anything():
    record = SeqRecord(Seq("ACGT" * 100), id="x")
    record.features = [
        SeqFeature(SimpleLocation(0, 400), type="source", qualifiers={"label": ["whole"]}),
        SeqFeature(SimpleLocation(10, 90), type="promoter", qualifiers={"label": ["pTDH3"]}),
        SeqFeature(SimpleLocation(100, 200), type="CDS", qualifiers={"gene": ["XI"]}),
        SeqFeature(SimpleLocation(210, 220), type="CDS", qualifiers={}),
    ]
    regions = align.regions_from_record(record)
    # `source` spans everything and would label every difference with it
    assert [r.label for r in regions] == ["pTDH3", "XI"]


# --------------------------------------------------------------------------- #
# the rows the viewer draws
# --------------------------------------------------------------------------- #


def test_every_row_is_the_same_width(reference):
    """Columns only mean anything if they line up across every clone."""
    clones = [
        align.compare(reference, rotate(reference, 100), name="a"),
        align.compare(reference, reference[:1000] + "GGG" + reference[1000:], name="b"),
        align.compare(reference, reference[:2000] + reference[2005:], name="c"),
    ]
    merged = align.rows(reference, clones)
    assert len(merged.reference) == merged.columns
    assert all(len(row) == merged.columns for _, row in merged.rows)
    assert merged.columns == len(reference) + 3


def test_two_clones_with_different_insertions_share_the_widest_column(reference):
    """Otherwise the same screen column means a different base on each row."""
    short = align.compare(reference, reference[:1500] + "AA" + reference[1500:], name="short")
    long = align.compare(reference, reference[:1500] + "AAAAAA" + reference[1500:], name="long")
    merged = align.rows(reference, [short, long])
    # One column group, six wide - the wider of the two. Its exact coordinate
    # is deliberately not asserted: inserting A's next to an A is genuinely
    # ambiguous and the aligner left-shifts it, which is not a bug and not
    # something this module should pretend to resolve. What must hold is that
    # both clones are given the *same* group, or the columns lie.
    assert len(merged.insertions) == 1
    at, width = merged.insertions[0]
    assert width == 6
    assert all(len(row) == merged.columns for _, row in merged.rows)
    # the clone with the shorter insertion is padded, not truncated
    rows = dict(merged.rows)
    assert rows["short"][at : at + 6] == "AA----"
    assert rows["long"][at : at + 6] == "AAAAAA"


def test_a_clean_clone_reproduces_the_reference_exactly(reference):
    merged = align.rows(reference, [align.compare(reference, rotate(reference, 77), name="a")])
    assert merged.reference == reference
    assert merged.rows[0][1] == reference


# --------------------------------------------------------------------------- #
# reading the vendor's file
# --------------------------------------------------------------------------- #


def test_the_format_is_read_from_the_content_not_the_name():
    assert seqio.sniff_format(">x\nACGT\n") == "fasta"
    assert seqio.sniff_format("LOCUS       x  100 bp\n") == "genbank"
    assert seqio.sniff_format("ACGTACGT\n") == "plain"


def test_bare_sequence_is_accepted():
    records = seqio.read_text_records("acgt\nACGT\n", "pasted")
    assert str(records[0].seq) == "ACGTACGT"


def test_something_that_is_not_sequence_is_refused_rather_than_stripped():
    """A vendor's error page must not come back as a 40 bp plasmid."""
    with pytest.raises(ValueError, match="not DNA|not GenBank"):
        seqio.read_text_records("Sorry, your order could not be found.", "reply.txt")


def test_a_fasta_with_several_records_reads_as_several():
    text = ">clone1\nACGTACGTAC\n>clone2\nTTTTGGGGAA\n"
    assert len(seqio.read_text_records(text, "plate.fa")) == 2


def test_a_numbered_wrapped_paste_reads_as_sequence():
    """What copying out of a plasmid editor actually puts on the clipboard."""
    text = (
        "        1 atgtcgaaac tgtttaagga aaatcttagt\n"
        "       31 gaggttgcct ggaatatcgc caaaagctat\n"
    )
    records = seqio.read_text_records(text, "pasted")
    assert str(records[0].seq) == (
        "ATGTCGAAACTGTTTAAGGAAAATCTTAGTGAGGTTGCCTGGAATATCGCCAAAAGCTAT"
    )


def test_a_pasted_origin_block_reads_as_sequence():
    """An ORIGIN section alone is not a GenBank record, but it is a paste."""
    text = "ORIGIN\n        1 acgtacgtac gtacgtacgt\n//\n"
    records = seqio.read_text_records(text, "pasted")
    assert str(records[0].seq) == "ACGTACGTACGTACGTACGT"


def test_stripping_digits_does_not_start_swallowing_prose():
    """The guard that keeps the tolerance above from going too far.

    Digits and whitespace are never bases, so dropping them hides nothing.
    Letters are left alone on purpose - filtering a vendor's apology down to
    whichever of its letters spell DNA is how you get a 40 bp plasmid.
    """
    with pytest.raises(ValueError, match="not DNA"):
        seqio.read_text_records("Order 12345 could not be found.", "reply.txt")


def test_lower_case_sequence_comes_back_upper():
    assert str(seqio.read_text_records("acgt\n", "x")[0].seq) == "ACGT"


# --------------------------------------------------------------------------- #
# the alignment budget
# --------------------------------------------------------------------------- #
#
# A 20 kb multigene reference against a full-plasmid read is what this screen
# exists for, and also the pair most likely to ask the pairwise aligner for
# something quadratic in twenty thousand. The guard has to hold without the
# caller ever seeing an exception.


def big(n: int, seed: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def test_a_divergent_stretch_past_the_budget_is_reported_not_raised():
    """The failure has to arrive as data. An exception here reaches the browser
    as a 500 and the panel simply stays empty, which looks like the app hanging
    rather than like an answer about the clone."""
    reference = big(20_000, 9)
    clone = reference[:6_000] + big(8_000, 10) + reference[14_000:]

    result = align.compare(reference, clone, name="diverged")

    assert result.placement.placed
    unalignable = [d for d in result.differences if d.kind == "unaligned"]
    assert len(unalignable) == 1
    assert unalignable[0].length > align.MAX_WINDOW


def test_the_unalignable_message_names_the_limit():
    """Otherwise it reads as the app failing rather than as a fact about the
    clone, and the number is the only thing that says which."""
    reference = big(20_000, 9)
    clone = reference[:6_000] + big(8_000, 10) + reference[14_000:]
    note = next(
        d for d in align.compare(reference, clone).differences if d.kind == "unaligned"
    )
    assert f"{align.MAX_WINDOW:,}" in note.describe()
    assert "diverges" in note.describe()


def test_a_clone_that_is_half_another_plasmid_still_answers():
    reference = big(20_000, 9)
    clone = reference[:10_000] + big(10_000, 11)
    result = align.compare(reference, clone, name="chimera")
    assert [d.kind for d in result.differences] == ["unaligned"]
    assert "unalignable" in result.verdict()


def test_an_unalignable_stretch_is_not_counted_as_matching():
    reference = big(20_000, 9)
    clone = reference[:6_000] + big(8_000, 10) + reference[14_000:]
    result = align.compare(reference, clone)
    assert result.identity(reference) < 70.0


def test_a_twenty_kb_read_of_the_real_thing_is_still_clean():
    """The guard must not fire on the case the screen is actually for."""
    reference = big(20_000, 9)
    result = align.compare(reference, reference[7_000:] + reference[:7_000])
    assert result.differences == []
    assert result.identity(reference) == 100.0


def test_a_window_just_under_the_budget_is_aligned_properly():
    """The limit has to sit where real work still gets done base by base."""
    reference = big(20_000, 9)
    clone = reference[:6_000] + big(align.MAX_WINDOW - 200, 12) + reference[6_000 + align.MAX_WINDOW - 200:]
    kinds = {d.kind for d in align.compare(reference, clone).differences}
    assert "unaligned" not in kinds
    assert "substitution" in kinds
