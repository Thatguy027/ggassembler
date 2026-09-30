"""The generic simulator, and Level 2 on top of it."""

from __future__ import annotations

import pytest

from ggassembler.core.assembly import Piece, assemble
from ggassembler.core.enzymes import BSAI, BSMBI, digest, find_sites
from ggassembler.core.library import Library
from ggassembler.core.parttypes import YTK, type_overhangs
from ggassembler.core.seqio import write_genbank
from ggassembler.levels import level2_cassette as level2

from . import synth

#: The canonical cassette of the paper: ConLS, pTDH3, Venus, tADH1, ConR1,
#: URA3, CEN6/ARS4, AmpR-ColE1 - one part per position around the circle.
CANONICAL = {
    "1": "ConLS",
    "2": "pTDH3",
    "3": "Venus",
    "4": "tADH1",
    "5": "ConR1",
    "6": "URA3",
    "7": "CEN6_ARS4",
    "8": "AmpR_ColE1",
}


@pytest.fixture
def cassette_library(tmp_path):
    """A folder holding one part plasmid per position of the circle."""
    for part_type, name in CANONICAL.items():
        record = synth.part_plasmid(part_type, name=name, seed=synth.seed_for(name))
        write_genbank(record, tmp_path / f"{name}.gb")
    library = Library(tmp_path)
    library.scan()
    return library


def pieces_for(library, selections=CANONICAL):
    return [library.piece(library.get(name)) for name in selections.values()]


# --------------------------------------------------------------------------- #
# the simulator
# --------------------------------------------------------------------------- #


def test_eight_parts_close_into_one_circle(cassette_library):
    result = assemble(pieces_for(cassette_library), BSAI, name="pCassette")

    assert result.ok
    assert not result.errors
    assert result.product.annotations["topology"] == "circular"
    assert len(result.parts) == 8
    assert len(result.junctions) == 8


def test_product_length_is_the_sum_of_its_parts(cassette_library):
    pieces = pieces_for(cassette_library)
    result = assemble(pieces, BSAI)
    assert result.length == sum(len(p) for p in pieces)


def test_the_product_gives_back_the_same_eight_fragments(cassette_library):
    """Slicing the product at its junctions returns the inputs, base for base."""
    pieces = {p.source_name: p for p in pieces_for(cassette_library)}
    result = assemble(list(pieces.values()), BSAI)

    product = str(result.product.seq)
    assert len(result.parts) == len(pieces)
    for part in result.parts:
        assert product[part.start:part.end] == pieces[part.source_name].seq

    # and the junction overhangs are the first four bases of each next part
    for junction in result.junctions:
        assert (product + product)[junction.position:junction.position + 4] == junction.overhang


def test_the_assembly_enzyme_is_consumed(cassette_library):
    """A finished cassette has no part-enzyme sites left: they were all cut out."""
    result = assemble(pieces_for(cassette_library), BSAI)
    assert find_sites(str(result.product.seq), BSAI) == []


def test_order_of_input_does_not_matter(cassette_library):
    pieces = pieces_for(cassette_library)
    forward = assemble(list(pieces), BSAI)
    shuffled = assemble(list(reversed(pieces)), BSAI)

    # the same circle, read from a different starting point
    assert shuffled.length == forward.length
    assert str(shuffled.product.seq) in str(forward.product.seq) * 2


def test_features_survive_into_the_product(cassette_library):
    result = assemble(pieces_for(cassette_library), BSAI)
    labels = [
        q
        for feature in result.product.features
        for q in feature.qualifiers.get("label", [])
    ]
    for name in CANONICAL.values():
        assert name in labels, f"{name} lost its span feature"
    assert sum(1 for label in labels if label.startswith("junction ")) == 8


def test_every_junction_is_annotated_with_what_it_joins(cassette_library):
    result = assemble(pieces_for(cassette_library), BSAI)
    notes = [
        note
        for feature in result.product.features
        for note in feature.qualifiers.get("note", [])
        if "->" in note
    ]
    assert len(notes) == 8
    assert any("via AACG" in note for note in notes)


def test_the_product_records_its_inputs_and_enzyme(cassette_library):
    result = assemble(pieces_for(cassette_library), BSAI, name="pCassette")
    comment = result.product.annotations["comment"]
    assert "BsaI" in comment
    for name in CANONICAL.values():
        assert name in comment


def test_a_linear_product_can_be_requested(cassette_library):
    result = assemble(pieces_for(cassette_library), BSAI, circular=False)
    assert result.product.annotations["topology"] == "linear"


# --------------------------------------------------------------------------- #
# Level 2
# --------------------------------------------------------------------------- #


def test_default_design_has_eight_slots():
    design = level2.CassetteDesign()
    assert [s.key for s in level2.slots(design)] == ["1", "2", "3", "4", "5", "6", "7", "8"]


def test_view_modes_change_the_slots_not_the_model():
    split = level2.CassetteDesign(split_3=True, split_4=True, split_8=True)
    assert [s.key for s in level2.slots(split)] == [
        "1", "2", "3a", "3b", "4a", "4b", "5", "6", "7", "8a", "8b",
    ]

    composite = level2.CassetteDesign(composite_left=True, composite_right=True)
    assert [s.key for s in level2.slots(composite)] == ["1", "234", "5", "678"]


def test_slots_tile_the_circle_without_gaps():
    """Whatever the view mode, each slot's 3' overhang is the next one's 5'."""
    for design in (
        level2.CassetteDesign(),
        level2.CassetteDesign(split_3=True, split_8=True),
        level2.CassetteDesign(composite_left=True, composite_right=True),
    ):
        active = level2.slots(design)
        for here, nxt in zip(active, active[1:], strict=False):
            assert here.three_prime == nxt.five_prime
        assert active[-1].three_prime == active[0].five_prime


def test_slots_offer_only_parts_that_fit(cassette_library):
    design = level2.CassetteDesign()
    counts = level2.option_counts(cassette_library, design)
    assert counts == {str(i): 1 for i in range(1, 9)}

    promoter_slot = next(s for s in level2.slots(design) if s.key == "2")
    assert [e.name for e in level2.options(cassette_library, promoter_slot)] == ["pTDH3"]


def test_level2_build_assembles_the_canonical_cassette(cassette_library):
    design = level2.CassetteDesign(selections=dict(CANONICAL), name="pCassette")
    result = level2.build(cassette_library, design)

    assert result.ok, [str(i) for i in result.issues]
    assert len(result.parts) == 8
    assert result.product.id == "pCassette"
    assert [p.part_type for p in result.parts] == list(CANONICAL)


def test_level2_parts_carry_their_colour(cassette_library):
    design = level2.CassetteDesign(selections=dict(CANONICAL))
    result = level2.build(cassette_library, design)
    assert all(part.color for part in result.parts)
    assert len({part.color for part in result.parts}) > 1


def test_integration_configuration_is_flagged(tmp_path):
    for part_type, name in (*CANONICAL.items(), ("8a", "AmpR_ColE1_a"), ("8b", "HomL")):
        if part_type == "8":
            continue
        write_genbank(synth.part_plasmid(part_type, name=name, seed=synth.seed_for(name)),
                      tmp_path / f"{name}.gb")
    library = Library(tmp_path)
    library.scan()

    selections = {k: v for k, v in CANONICAL.items() if k != "8"}
    selections.update({"8a": "AmpR_ColE1_a", "8b": "HomL"})
    design = level2.CassetteDesign(selections=selections, split_8=True)
    assert design.is_integration

    result = level2.build(library, design)
    assert result.ok, [str(i) for i in result.issues]
    assert "integration" in result.codes()
    assert any("NotI" in i.message for i in result.issues)


def test_composite_slots_assemble_the_same_circle(tmp_path):
    """A 234 part and a 678 part build the same product as six separate parts."""
    for part_type, name in (("1", "ConLS"), ("234", "TU"), ("5", "ConR1"), ("678", "Backbone")):
        write_genbank(synth.part_plasmid(part_type, name=name, seed=synth.seed_for(name)),
                      tmp_path / f"{name}.gb")
    library = Library(tmp_path)
    library.scan()

    design = level2.CassetteDesign(
        selections={"1": "ConLS", "234": "TU", "5": "ConR1", "678": "Backbone"},
        composite_left=True,
        composite_right=True,
    )
    result = level2.build(library, design)
    assert result.ok, [str(i) for i in result.issues]
    assert len(result.parts) == 4
    assert find_sites(str(result.product.seq), BSAI) == []


def test_a_cassette_built_here_is_re_detected_as_a_part_free_construct(cassette_library):
    """The product must not look like a part plasmid when it lands in the library."""
    from ggassembler.core.library import describe

    result = level2.build(cassette_library, level2.CassetteDesign(selections=dict(CANONICAL)))
    entry = describe(result.product)
    assert entry.call.part_type is None
    assert entry.sites.part_enzyme_total == 0


def test_scheme_carries_the_enzymes_rather_than_hard_coding_them():
    assert YTK.part_enzyme is BSAI
    assert YTK.multigene_enzyme is BSMBI
    assert YTK.overhang_length == 4
    assert type_overhangs("234", YTK) == ("AACG", "GCTG")


def _connector_library(tmp_path, left_overhang: str, right_overhang: str) -> Library:
    write_genbank(synth.connector_plasmid("1", left_overhang, name="ConL"), tmp_path / "ConL.gb")
    write_genbank(synth.connector_plasmid("5", right_overhang, name="ConR"), tmp_path / "ConR.gb")
    for part_type, name in (("234", "TU"), ("678", "Backbone")):
        write_genbank(synth.part_plasmid(part_type, name=name, seed=synth.seed_for(name)),
                      tmp_path / f"{name}.gb")
    library = Library(tmp_path)
    library.scan()
    return library


def _composite_design():
    return level2.CassetteDesign(
        selections={"1": "ConL", "234": "TU", "5": "ConR", "678": "Backbone"},
        composite_left=True,
        composite_right=True,
    )


def test_a_cassette_with_one_connector_at_both_ends_is_flagged(tmp_path):
    """ConL1 with ConR1 leaves a cassette that cannot be ordered at Level 3."""
    library = _connector_library(tmp_path, "CCAA", "CCAA")
    result = level2.build(library, _composite_design())

    assert result.ok, [str(i) for i in result.issues]
    assert "same_connector" in result.codes()
    assert any("CCAA" in i.message for i in result.warnings)


def test_different_connectors_raise_no_such_warning(tmp_path):
    library = _connector_library(tmp_path, "CCAA", "GATG")
    result = level2.build(library, _composite_design())

    assert result.ok
    assert "same_connector" not in result.codes()


def test_a_finished_cassette_reports_its_connector_chain(tmp_path):
    """The product must come back as a Level 3 input, with its ends readable."""
    from ggassembler.core.library import describe

    library = _connector_library(tmp_path, "CCAA", "GATG")
    result = level2.build(library, _composite_design())

    entry = describe(result.product)
    assert entry.is_cassette
    assert entry.cassette_overhangs == ("CCAA", "GATG")
    assert entry.sites.part_enzyme_total == 0


# ------------------------------------------------- the screen's worked example --


def test_the_default_design_fills_every_slot(cassette_library):
    design = level2.default_design(cassette_library)
    keys = [s.key for s in level2.slots(design)]
    assert sorted(design.selections) == sorted(keys)


def test_the_default_design_actually_assembles(cassette_library):
    design = level2.default_design(cassette_library)
    result = level2.build(cassette_library, design)
    assert result.ok, [str(i) for i in result.issues]
    assert len(result.parts) == 8


def test_the_default_avoids_a_part_with_internal_sites(tmp_path):
    """Given a choice, the example should not open on a broken part."""
    write_genbank(synth.plasmid_with_internal_site("2", name="dirty"), tmp_path / "dirty.gb")
    for part_type, name in CANONICAL.items():
        if part_type == "2":
            continue
        write_genbank(synth.part_plasmid(part_type, name=name, seed=synth.seed_for(name)),
                      tmp_path / f"{name}.gb")
    write_genbank(synth.part_plasmid("2", name="clean", seed=synth.seed_for("clean")),
                  tmp_path / "clean.gb")
    library = Library(tmp_path)
    library.scan()

    design = level2.default_design(library)
    assert design.selections["2"] == "clean"


def test_the_default_picks_two_different_connectors(tmp_path):
    """A cassette with one connector at both ends is useless at Level 3."""
    write_genbank(synth.connector_plasmid("1", "CCAA", name="ConL1"), tmp_path / "ConL1.gb")
    write_genbank(synth.connector_plasmid("5", "CCAA", name="ConR1"), tmp_path / "ConR1.gb")
    write_genbank(synth.connector_plasmid("5", "GATG", name="ConR2"), tmp_path / "ConR2.gb")
    for part_type, name in CANONICAL.items():
        if part_type in ("1", "5"):
            continue
        write_genbank(synth.part_plasmid(part_type, name=name, seed=synth.seed_for(name)),
                      tmp_path / f"{name}.gb")
    library = Library(tmp_path)
    library.scan()

    design = level2.default_design(library)
    result = level2.build(library, design)
    assert design.selections["5"] == "ConR2"
    assert "same_connector" not in result.codes()


def test_a_part_plasmids_own_marker_does_not_warn(cassette_library):
    """Only the position 8 part contributes a marker to the product."""
    design = level2.CassetteDesign(selections=dict(CANONICAL))
    result = level2.build(cassette_library, design)
    assert "mixed_markers" not in result.codes()


def test_the_default_endpoint_returns_a_finished_assembly(cassette_library):
    from fastapi.testclient import TestClient

    from ggassembler.api.main import create_app

    with TestClient(create_app(cassette_library.folder)) as client:
        body = client.get("/api/level2/default").json()
        assert body["ok"] is True
        assert len(body["selections"]) == 8
        assert len(body["parts"]) == 8
        assert body["length"] > 0
        assert body["slots"][0]["selected"]


# ------------------------------------------------ what the product is labelled --


def test_parts_are_labelled_by_what_they_are_not_where_they_came_from(tmp_path):
    from Bio.SeqFeature import FeatureLocation, SeqFeature

    for part_type, name in CANONICAL.items():
        record = synth.part_plasmid(part_type, name=name, seed=synth.seed_for(name))
        record.features.append(
            SeqFeature(FeatureLocation(20, 90), type="promoter" if part_type == "2" else "CDS",
                       qualifiers={"label": [f"{name}_annotation"]})
        )
        write_genbank(record, tmp_path / f"{name}.gb")
    library = Library(tmp_path)
    library.scan()

    result = level2.build(library, level2.CassetteDesign(selections=dict(CANONICAL)))
    labels = [
        label
        for feature in result.product.features
        for label in feature.qualifiers.get("label", [])
    ]
    assert "pTDH3_annotation" in labels, "the span should carry the component name"
    assert "pTDH3" not in labels, "not the source plasmid name"

    # provenance survives, in the note
    notes = " ".join(
        note for f in result.product.features for note in f.qualifiers.get("note", [])
    )
    assert "from pTDH3" in notes


def test_the_product_does_not_describe_itself_as_a_part(cassette_library):
    """A note reading 'type 2 part from ...' would be read back as a type call."""
    from ggassembler.core.library import describe

    result = level2.build(cassette_library, level2.CassetteDesign(selections=dict(CANONICAL)))
    entry = describe(result.product)

    assert entry.call.part_type is None
    assert entry.call.confidence == "none"
    notes = " ".join(
        note for f in result.product.features for note in f.qualifiers.get("note", [])
    )
    assert "position 2 part" in notes
    assert "type 2 part" not in notes


def test_placed_parts_carry_their_component_to_the_api(cassette_library):
    result = level2.build(cassette_library, level2.CassetteDesign(selections=dict(CANONICAL)))
    for part in result.parts:
        assert part.label  # never blank: falls back to the source name


# --------------------------------------------------------------------------- #
# reading a finished cassette back into parts
# --------------------------------------------------------------------------- #


def test_a_cassette_decomposes_into_the_parts_that_built_it(tmp_path):
    """Build eight parts into a cassette, then read the cassette back.

    This is the property the whole feature rests on: what comes out has to be
    what went in. The assembly consumed every BsaI site, so the decomposition
    cannot digest - it matches sequence - and a round trip is the only honest
    check that the matching lines up with the boundaries.
    """
    for part_type, name in CANONICAL.items():
        write_genbank(
            synth.part_plasmid(part_type, name=name, seed=synth.seed_for(name)),
            tmp_path / f"{name}.gb",
        )
    library = Library(tmp_path)
    library.scan()

    design = level2.CassetteDesign(selections=dict(CANONICAL), name="built")
    built = level2.build(library, design)
    assert built.ok, [i.message for i in built.errors]

    # file the product back into the library and read it as a fresh plasmid
    write_genbank(built.product, tmp_path / "built.gb")
    library = Library(tmp_path)
    library.scan()

    found = level2.decompose(library, library.get("built"))
    assert found.complete, f"{found.covered} of {found.length} bp explained"
    assert found.selections == CANONICAL
    assert found.gaps == []


def test_the_decomposition_covers_every_base_exactly_once(tmp_path):
    for part_type, name in CANONICAL.items():
        write_genbank(
            synth.part_plasmid(part_type, name=name, seed=synth.seed_for(name)),
            tmp_path / f"{name}.gb",
        )
    library = Library(tmp_path)
    library.scan()
    built = level2.build(library, level2.CassetteDesign(selections=dict(CANONICAL)))
    write_genbank(built.product, tmp_path / "built.gb")
    library = Library(tmp_path)
    library.scan()

    found = level2.decompose(library, library.get("built"))
    assert found.covered == found.length
    # the matches abut: each one ends where the next begins, all the way round
    for before, after in zip(found.matches, found.matches[1:], strict=False):
        assert before.end == after.start
    assert found.matches[-1].end == found.matches[0].start


def test_a_cassette_with_a_part_from_outside_the_library_reports_the_gap(tmp_path):
    """Silence would be worse: the panels would look like the whole construct."""
    for part_type, name in CANONICAL.items():
        write_genbank(
            synth.part_plasmid(part_type, name=name, seed=synth.seed_for(name)),
            tmp_path / f"{name}.gb",
        )
    library = Library(tmp_path)
    library.scan()
    built = level2.build(library, level2.CassetteDesign(selections=dict(CANONICAL)))
    write_genbank(built.product, tmp_path / "built.gb")

    # take the type 3 part off the shelf; its stretch can no longer be explained
    (tmp_path / f"{CANONICAL['3']}.gb").unlink()
    library = Library(tmp_path)
    library.scan()

    found = level2.decompose(library, library.get("built"))
    assert not found.complete
    assert found.gaps, "an unexplained stretch must be reported"
    assert "3" not in found.selections
    assert found.covered < found.length


def test_a_plasmid_that_is_no_cassette_at_all_decomposes_to_nothing(tmp_path):
    write_genbank(synth.plasmid_without_bsai(name="mystery"), tmp_path / "mystery.gb")
    library = Library(tmp_path)
    library.scan()
    found = level2.decompose(library, library.get("mystery"))
    assert found.matches == []
    assert not found.complete


def test_a_split_part_turns_the_panel_split_on(tmp_path):
    """The panels a decomposition needs may not be the ones on screen."""
    selections = {**CANONICAL}
    del selections["3"]
    selections["3a"] = "part3a"
    selections["3b"] = "part3b"
    for part_type, name in selections.items():
        write_genbank(
            synth.part_plasmid(part_type, name=name, seed=synth.seed_for(name)),
            tmp_path / f"{name}.gb",
        )
    library = Library(tmp_path)
    library.scan()
    built = level2.build(
        library, level2.CassetteDesign(selections=dict(selections), split_3=True)
    )
    assert built.ok, [i.message for i in built.errors]
    write_genbank(built.product, tmp_path / "built.gb")
    library = Library(tmp_path)
    library.scan()

    found = level2.decompose(library, library.get("built"))
    assert found.split_3, "a 3a/3b pair has to open the split view"
    assert found.selections.get("3a") == "part3a"
    assert "3" not in found.selections


def test_an_unexplained_stretch_reports_the_position_it_occupies(tmp_path):
    """"1,750 bp of type 3 that isn't a part I know" beats "1,750 bp unexplained".

    The four bases at each end of the gap are the overhangs the missing part
    was ligated by, and they name its position even when the part itself is
    nowhere in the library.
    """
    for part_type, name in CANONICAL.items():
        write_genbank(
            synth.part_plasmid(part_type, name=name, seed=synth.seed_for(name)),
            tmp_path / f"{name}.gb",
        )
    library = Library(tmp_path)
    library.scan()
    built = level2.build(library, level2.CassetteDesign(selections=dict(CANONICAL)))
    write_genbank(built.product, tmp_path / "built.gb")

    # take the type 3 part off the shelf: its stretch is now unexplained
    (tmp_path / f"{CANONICAL['3']}.gb").unlink()
    library = Library(tmp_path)
    library.scan()

    found = level2.decompose(library, library.get("built"))
    assert len(found.unmatched) == 1
    gap = found.unmatched[0]
    assert gap.part_type == "3", f"the gap should name position 3, got {gap.part_type}"
    assert gap.left_overhang == "TATG" and gap.right_overhang == "ATCC"
    assert gap.length > 0


def test_the_old_gap_shape_still_works_for_callers_that_only_count(tmp_path):
    for part_type, name in CANONICAL.items():
        write_genbank(
            synth.part_plasmid(part_type, name=name, seed=synth.seed_for(name)),
            tmp_path / f"{name}.gb",
        )
    library = Library(tmp_path)
    library.scan()
    built = level2.build(library, level2.CassetteDesign(selections=dict(CANONICAL)))
    write_genbank(built.product, tmp_path / "built.gb")
    (tmp_path / f"{CANONICAL['3']}.gb").unlink()
    library = Library(tmp_path)
    library.scan()

    found = level2.decompose(library, library.get("built"))
    assert found.gaps == [(u.start, u.length) for u in found.unmatched]
    assert not found.complete
