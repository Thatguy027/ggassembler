"""Will a vendor actually make this gene?

A sequence can be clean of every restriction site and still be refused, quoted
at a premium, or delivered wrong. Synthesis houses screen on sequence features:
long single-base runs, exact internal repeats, windows of extreme GC. None of
that shows up in a translation check.

Worth testing rather than eyeballing, because the strategy that writes these
genes is exactly the one that produces them - always taking the most preferred
codon is a deterministic map from residue to bases, so a tract of one amino
acid becomes a tract of one codon.
"""

from __future__ import annotations

import random

from ggassembler.core.codons import (
    ECOLI,
    HOMOPOLYMER_LIMIT,
    SCEREVISIAE,
    back_translate,
    feasibility,
    gc_windows,
    longest_homopolymer,
    longest_repeat,
    translate,
)
from ggassembler.core.enzymes import BBSI, BSAI, BSMBI, NOTI

AVOID = (BSAI, BSMBI, BBSI, NOTI)
RESIDUES = "ACDEFGHIKLMNPQRSTVWY"


def filler(n: int, seed: int = 0) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice(RESIDUES) for _ in range(n))


def noise(n: int, seed: int = 0) -> str:
    """Random DNA, from one generator.

    Deliberately not `random.Random(seed).choice(...)` inside the loop: that
    builds a fresh generator per base and hands back the same one every time,
    so the "random" sequence is a homopolymer. It cost a confusing failure in
    this very file - the repeat finder was right and the fixture was not.
    """
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


# --------------------------------------------------------------------------- #
# homopolymer runs
# --------------------------------------------------------------------------- #


def test_the_longest_run_is_found_with_its_base_and_position():
    assert longest_homopolymer("ACGTTTTTTTTTACGT") == ("T", 9, 3)


def test_a_run_at_the_very_end_still_counts():
    """The scan has to close the open run when it falls off the end."""
    assert longest_homopolymer("ACGAAAA") == ("A", 4, 3)


def test_a_run_at_the_very_start_still_counts():
    assert longest_homopolymer("GGGGGACT") == ("G", 5, 0)


def test_the_first_of_two_equal_runs_wins():
    assert longest_homopolymer("AAATCCC") == ("A", 3, 0)


def test_no_sequence_is_not_an_error():
    """A report on nothing is a report, not a failure."""
    assert longest_homopolymer("") == ("", 0, -1)


def test_one_base_is_a_run_of_one():
    assert longest_homopolymer("A") == ("A", 1, 0)


# --------------------------------------------------------------------------- #
# exact repeats
# --------------------------------------------------------------------------- #


def test_an_exact_repeat_is_found_with_both_positions():
    motif = "ACGTACGTACGTACGTACGT"          # 20 bp
    dna = "TTTT" + motif + "GGGG" + motif + "AAAA"
    length, first, second = longest_repeat(dna)
    assert length >= 20
    assert dna[first : first + 20] == dna[second : second + 20]
    assert first < second


def test_the_longest_repeat_is_found_not_merely_a_long_one():
    """Binary search over the length has to land on the maximum, not a hit."""
    motif = noise(64, seed=1)
    length, first, second = longest_repeat(motif + "TT" + motif)
    assert length == 64
    assert (first, second) == (0, 66)


def test_a_sequence_with_no_long_repeat_says_so():
    assert longest_repeat(noise(600, seed=2)) == (0, -1, -1)


def test_a_sequence_too_short_to_hold_two_copies_says_so():
    assert longest_repeat("ACGTACGTACGTACGTACGT") == (0, -1, -1)


def test_overlapping_copies_still_count_as_a_repeat():
    """A tandem run is a repeat of itself, offset by one period."""
    length, first, second = longest_repeat("AC" * 40)
    assert length >= 20
    assert second > first


# --------------------------------------------------------------------------- #
# GC across windows
# --------------------------------------------------------------------------- #


def test_gc_is_measured_across_windows_not_over_the_whole_gene():
    """The case the whole-gene number hides: a comfortable average, a dead window."""
    dna = "AT" * 100 + "GC" * 100                 # 50% overall, 0% and 100% locally
    low, high, outside = gc_windows(dna)
    assert low == 0.0 and high == 100.0
    assert outside > 0


def test_an_even_sequence_has_no_window_outside_the_band():
    # a 50 bp window does not divide into ACGT, so it reads 48-52 rather than
    # a flat 50 - the point is that nothing leaves the band
    low, high, outside = gc_windows("ACGT" * 100)
    assert 45.0 <= low <= high <= 55.0
    assert outside == 0


def test_a_sequence_shorter_than_one_window_is_measured_whole():
    """Otherwise a short part reports nothing at all, which reads as passing."""
    low, high, outside = gc_windows("GGGGCCCC")
    assert (low, high) == (100.0, 100.0)
    assert outside == 1


def test_the_window_count_is_every_window_not_every_base():
    dna = "AT" * 50                                # 100 bp, all windows at 0%
    _, _, outside = gc_windows(dna, window=50)
    assert outside == len(dna) - 50 + 1


def test_no_sequence_is_not_an_error_here_either():
    assert gc_windows("") == (0.0, 0.0, 0)


# --------------------------------------------------------------------------- #
# the three together
# --------------------------------------------------------------------------- #


def test_feasibility_reports_every_check_with_its_own_verdict():
    # not "ACGT" * 100: that is a 396 bp tandem repeat, and the repeat check is
    # right to say so
    clean = feasibility(noise(600, seed=2))
    assert clean["homopolymer"]["ok"] and clean["gc"]["ok"] and clean["repeat"]["ok"]

    bad = feasibility("A" * 40 + "ACGT" * 40)
    assert not bad["homopolymer"]["ok"]
    assert bad["homopolymer"]["length"] == 41       # the run runs into the ACGT
    assert bad["homopolymer"]["limit"] == HOMOPOLYMER_LIMIT


def test_feasibility_carries_the_proof_not_just_the_verdict():
    """Every row has to be checkable: a position, or it is only an assertion."""
    motif = "ACGTACGTACGTACGTACGTAA"
    report = feasibility(motif + "TTT" + motif)
    assert report["repeat"]["first"] >= 0 and report["repeat"]["second"] >= 0
    assert report["homopolymer"]["at"] >= 0


# --------------------------------------------------------------------------- #
# the tie-break during codon choice
# --------------------------------------------------------------------------- #


def test_a_poly_lysine_tract_does_not_become_ninety_adenines():
    """The case that justifies smoothing, and it is not hypothetical.

    *E. coli*'s preferred lysine codon is AAA, so thirty lysines written by
    plain highest-preference choice are ninety consecutive adenines - a gene no
    vendor will make, produced by a strategy that never made a wrong decision.
    """
    protein = "M" + "K" * 30 + filler(150, seed=3)

    plain, _ = back_translate(protein, avoid=AVOID, table=ECOLI, smooth=False)
    assert longest_homopolymer(plain)[1] > 60, "the case being fixed has gone away"

    smoothed, _ = back_translate(protein, avoid=AVOID, table=ECOLI)
    assert longest_homopolymer(smoothed)[1] < longest_homopolymer(plain)[1] / 5


def test_smoothing_never_changes_the_protein():
    """It is a choice between synonyms, so this must hold by construction."""
    protein = "M" + "K" * 30 + "E" * 20 + filler(200, seed=5)
    dna, _ = back_translate(protein, avoid=AVOID, table=ECOLI)
    assert translate(dna).rstrip("*") == protein


def test_smoothing_never_introduces_a_site():
    """It only ever picks among codons that were already legal."""
    from ggassembler.core.enzymes import find_sites

    protein = "M" + "K" * 40 + filler(300, seed=7)
    dna, _ = back_translate(protein, avoid=AVOID, table=ECOLI)
    for enzyme in AVOID:
        assert not find_sites(dna, enzyme, circular=False)


def test_smoothing_leaves_a_sequence_that_needs_none_alone():
    """No churn where there is nothing to fix - the yeast table rarely runs."""
    protein = filler(300, seed=11)
    with_it, _ = back_translate(protein, avoid=AVOID, table=SCEREVISIAE)
    without, _ = back_translate(protein, avoid=AVOID, table=SCEREVISIAE, smooth=False)
    assert with_it == without


def test_an_unsmoothable_protein_still_gets_a_sequence():
    """Smoothing is a preference, never a requirement.

    Tryptophan and methionine have one codon each, so a run of them cannot be
    broken up at all. The answer is still a gene.
    """
    protein = "MWWWWWWWWWWMMMMMMMMMM" + filler(60, seed=13)
    dna, _ = back_translate(protein, avoid=AVOID, table=SCEREVISIAE)
    assert translate(dna).rstrip("*") == protein
