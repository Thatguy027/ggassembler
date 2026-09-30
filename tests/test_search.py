"""Annotation search, and the label cleaning it depends on.

The question this feature answers is "which plasmid has the TDH3 promoter in
it", so the tests are written the same way: build a library where the answer is
known, ask, and check the right plasmid comes back *first*. Ranking is the
feature - a search that returns the right plasmid ninth is a search nobody uses.
"""

from __future__ import annotations

import pytest
from Bio.SeqFeature import FeatureLocation, SeqFeature

from ggassembler.core import search
from ggassembler.core.library import (
    FeatureInfo,
    Library,
    clean_label,
    features_within,
    name_component,
    summarise,
)

from ggassembler.core.enzymes import BSAI
from ggassembler.core.seqio import write_genbank

from . import synth

#: Where `synth.part_plasmid` puts things: a forward site, then the 5' overhang,
#: then the insert. The released fragment starts at the overhang.
INSERT_START = len(synth.forward_site(BSAI)) + 4
INSERT_LENGTH = 120


def annotate(record, label, kind, start, end):
    record.features.append(
        SeqFeature(FeatureLocation(start, end), type=kind, qualifiers={"label": [label]})
    )
    return record


@pytest.fixture
def library(tmp_path):
    """Three parts labelled the way a real YTK library labels them."""
    made = {
        "pYTK009": ("2", "ScTDH3 Promoter", "promoter"),
        "pYTK056": ("4", "ScTDH1 Terminator", "terminator"),
        "pYTK033": ("3", "Venus", "CDS"),
    }
    for name, (part_type, label, kind) in made.items():
        record = synth.part_plasmid(part_type, name=name)
        record.features = []
        annotate(record, label, kind, INSERT_START, INSERT_START + INSERT_LENGTH)
        # every plasmid in the kit carries this, out in the backbone
        annotate(record, "AmpR", "CDS", INSERT_START + INSERT_LENGTH + 30, len(record.seq) - 10)
        write_genbank(record, tmp_path / f"{name}.gb")

    lib = Library(tmp_path)
    lib.scan()
    # the fixture is only worth trusting if the layout is what it assumes
    assert {e.name for e in lib.unique_entries()} == set(made)
    assert lib.get("pYTK009").component == "ScTDH3 Promoter"
    return lib


# --------------------------------------------------------------------------- #
# label cleaning
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(
    "raw, cleaned",
    [
        ("XYL2_amplicon_BsaI_largest", "XYL2"),
        ("gal2mut_BsaI_largest", "gal2mut"),
        ("XYL1_K270R_amplicon_BsaI_largest", "XYL1_K270R"),
        ("gre3_3primeHA_BsaI_largest", "gre3_3primeHA"),
        ("malibu_pcr_BsaI_largest", "malibu"),
        ("ScTDH3 Promoter", "ScTDH3 Promoter"),  # nothing to strip
        ("mScarlet-I", "mScarlet-I"),
    ],
)
def test_pipeline_boilerplate_is_stripped_not_the_gene(raw, cleaned):
    assert clean_label(raw) == cleaned


def test_a_label_that_is_all_boilerplate_survives_rather_than_vanishing():
    """Stripping must never return nothing: an odd label beats no label."""
    assert clean_label("BsaI_largest") == "BsaI"
    assert clean_label("amplicon") == "amplicon"


def test_mscarlet_is_not_mistaken_for_a_golden_gate_scar():
    """The old filter matched `scar` anywhere, which ate every mScarlet part."""
    assert summarise([FeatureInfo("CDS", "mScarlet-I", 0, 700, 1)]) == "mScarlet-I"
    assert summarise([FeatureInfo("misc_feature", "BsaI scar", 0, 4, 1)]) == ""
    assert summarise([FeatureInfo("misc_feature", "scar_1", 0, 4, 1)]) == ""


def test_a_benchling_translation_says_nothing_and_is_dropped():
    assert summarise([FeatureInfo("CDS", "Benchling translation", 0, 700, 1)]) == ""


def test_a_binding_site_describes_a_fragment_only_when_nothing_else_does():
    site = FeatureInfo("protein_bind", "I-SceI Recognition Site", 0, 18, 1)
    gene = FeatureInfo("CDS", "Venus", 0, 700, 1)
    assert summarise([site]) == "I-SceI Recognition Site"
    assert summarise([site, gene]) == "Venus"


# --------------------------------------------------------------------------- #
# the name fallback
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(
    "stem, expected",
    [
        ("L114_ERG10_yeast_gene", "ERG10"),
        ("L105_YBR110W_gene", "YBR110W"),
        ("L83_electra_4a", "electra"),
        ("L50_SpiRY_Type_3", "SpiRY"),
        ("L95_Cas9_Nickase", "Cas9 Nickase"),
        ("C9_Venus", "Venus"),
        ("GG212_Spobooster", "Spobooster"),
        # a gene name is shaped like an index code; the letters/digits balance
        # is what tells them apart, and getting it wrong renames this to "gene"
        ("TDH3_gene", "TDH3"),
        ("PHO5_gene", "PHO5"),
        ("HSP12_gene", "HSP12"),
        ("KAT5_Mutagenesis_Vector", "KAT5 Mutagenesis"),
    ],
)
def test_a_name_with_no_annotation_still_says_something(stem, expected):
    assert name_component(stem) == expected


@pytest.mark.parametrize("stem", ["pGNB043", "RAD21", "ACAT1", "L25_2", "pYTK009"])
def test_a_name_that_only_repeats_itself_adds_nothing(stem):
    """The name is already on screen; echoing it as a description is noise."""
    assert name_component(stem) == ""


def test_the_name_fallback_is_recorded_as_weaker_evidence(tmp_path):
    record = synth.part_plasmid("3")
    record.features = []
    write_genbank(record, tmp_path / "L114_ERG10_yeast_gene.gb")
    lib = Library(tmp_path)
    lib.scan()
    entry = lib.get("L114_ERG10_yeast_gene")
    assert entry.component == "ERG10"
    assert entry.component_source == "filename"
    # and it never touches the type call, which stays the digest's business
    assert entry.call.confidence == "digest"


# --------------------------------------------------------------------------- #
# span containment
# --------------------------------------------------------------------------- #


def test_an_annotation_drawn_across_the_overhang_still_counts_as_inside():
    """A fragment carries its left overhang and not its right; annotators often
    draw the insert between the two, so the feature runs a few bases past the
    end. That is the same feature, not a backbone one."""
    feature = FeatureInfo("misc_feature", "mTurquoise", 700, 1420, 1)
    assert features_within([feature], (699, 1419), 2383) == [feature]


def test_a_whole_plasmid_annotation_is_not_inside_the_part():
    """The same tolerance read from the other side: only a fifth of a
    plasmid-length feature is in the part, so it is not the part's label."""
    feature = FeatureInfo("rep_origin", "ColE1 origin", 0, 2383, 1)
    assert features_within([feature], (699, 1419), 2383) == []


def test_a_backbone_feature_stays_out():
    feature = FeatureInfo("CDS", "CmR", 1541, 2203, 1)
    assert features_within([feature], (699, 1419), 2383) == []


# --------------------------------------------------------------------------- #
# scoring
# --------------------------------------------------------------------------- #


def test_an_exact_label_beats_a_word_beats_a_fragment_of_one():
    assert search.score_term("venus", "Venus") == 3.0
    assert search.score_term("venus", "Venus CDS") == 2.0
    assert search.score_term("enus", "Venus") == 1.0
    assert search.score_term("nothing", "Venus") == 0.0


def test_a_species_prefix_does_not_hide_the_gene():
    """Half this library is named ScTDH3, KlLEU2, mScarlet. A camelCase break
    has to read as a word boundary or `tdh` ranks the promoter below every
    cassette that happens to spell the label without the prefix."""
    assert search.score_term("tdh", "ScTDH3 Promoter") == 2.0
    assert search.score_term("scarlet", "mScarlet-I") == 2.0
    # but a genuine mid-word match is still only a mid-word match
    assert search.score_term("dh3", "ScTDH3 Promoter") == 1.0


def test_punctuation_is_the_last_resort_not_the_first():
    assert search.score_term("tdh3", "TDH-3 promoter") == 0.75


# --------------------------------------------------------------------------- #
# searching
# --------------------------------------------------------------------------- #


def test_one_term_finds_both_the_promoter_and_the_terminator(library):
    hits = search.search(library, "tdh")[0]
    names = [h.name for h in hits]
    assert "pYTK009" in names and "pYTK056" in names
    assert "pYTK033" not in names


def test_a_hit_reports_the_part_type_it_would_fill(library):
    by_name = {h.name: h for h in search.search(library, "tdh")[0]}
    assert by_name["pYTK009"].part_type == "2"
    assert by_name["pYTK056"].part_type == "4"


def test_every_term_has_to_match(library):
    assert [h.name for h in search.search(library, "tdh terminator")[0]] == ["pYTK056"]
    assert search.search(library, "tdh venus")[0] == []


def test_a_part_match_outranks_a_backbone_one(library):
    """Every plasmid here carries AmpR, so searching for it must not look like
    a real hit - but it is still findable."""
    hits = search.search(library, "ampr")[0]
    assert len(hits) == 3
    assert all(h.where == "backbone" and not h.in_part for h in hits)

    venus = search.search(library, "venus")[0]
    assert venus[0].name == "pYTK033"
    assert venus[0].where in ("part", "component") and venus[0].in_part
    assert venus[0].score > hits[0].score


def test_an_empty_query_returns_nothing_rather_than_everything(library):
    assert search.search(library, "")[0] == []
    assert search.search(library, "   ")[0] == []


def test_a_search_can_be_narrowed_to_one_position(library):
    assert [h.name for h in search.search(library, "tdh", part_type="4")[0]] == ["pYTK056"]
    assert search.search(library, "tdh", part_type="1")[0] == []


def test_results_are_capped(library):
    assert len(search.search(library, "ampr", limit=2)[0]) == 2


def test_the_order_is_stable(library):
    first = [h.name for h in search.search(library, "ampr")[0]]
    for _ in range(5):
        assert [h.name for h in search.search(library, "ampr")[0]] == first


# --------------------------------------------------------------------------- #
# what a screen is allowed to offer
# --------------------------------------------------------------------------- #


@pytest.fixture
def library_with_a_finished_construct(library, tmp_path):
    """A multigene construct carrying the labels of the parts that built it.

    This is the real failure it guards against: a 10.8 kb assembly annotated
    `L6_CYC1_type4` was offered as a terminator part, because a label naming a
    component of a construct is not a statement about the construct.
    """
    record = synth.plasmid_without_bsai(name="xylb_gg50")
    annotate(record, "L6_CYC1_type4", "misc_feature", 10, 40)
    annotate(record, "xylb", "CDS", 60, 400)
    write_genbank(record, tmp_path / "xylb_gg50.gb")
    lib = Library(tmp_path)
    lib.scan()
    return lib


def test_a_finished_construct_is_found_but_never_typed(library_with_a_finished_construct):
    hits = search.search(library_with_a_finished_construct, "xylb")[0]
    assert [h.name for h in hits] == ["xylb_gg50"]
    assert hits[0].part_type is None
    assert not hits[0].usable_as_part


def test_the_cassette_screen_is_offered_nothing_it_cannot_use(
    library_with_a_finished_construct,
):
    hits, suppressed = search.search(library_with_a_finished_construct, "xylb", usable="part")
    assert hits == []
    assert suppressed == 1, "the match is withheld, and the screen is told so"


def test_a_real_part_survives_the_same_filter(library):
    hits, suppressed = search.search(library, "venus", usable="part")
    assert [h.name for h in hits] == ["pYTK033"]
    assert all(h.usable_as_part for h in hits)
    assert suppressed == 0


def test_a_part_is_not_offered_to_the_multigene_screen(library):
    hits, suppressed = search.search(library, "venus", usable="multigene")
    assert hits == []
    assert suppressed == 1


def test_an_unfiltered_search_still_returns_everything(library_with_a_finished_construct):
    hits, suppressed = search.search(library_with_a_finished_construct, "xylb")
    assert len(hits) == 1 and suppressed == 0
