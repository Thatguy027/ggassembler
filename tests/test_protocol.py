"""Reaction setup: the fmol maths, and refusing to invent what is missing.

The arithmetic here goes onto a bench, so the tests check it against numbers
you can verify by hand rather than against whatever the code happened to
produce. The other half of the file is about the gap: a concentration is a
property of a tube, and a protocol that quietly assumed one would be worse
than one that says it does not know.
"""

from __future__ import annotations

import pytest

from ggassembler.core import protocol
from ggassembler.core.enzymes import BBSI, BSAI, BSMBI
from ggassembler.core.library import Library
from ggassembler.core.seqio import write_genbank

from . import synth


# --------------------------------------------------------------------------- #
# the arithmetic
# --------------------------------------------------------------------------- #


def test_molecular_weight_of_a_kilobase_is_about_six_hundred_kilodaltons():
    assert protocol.molecular_weight(1000) == pytest.approx(618_000, rel=1e-3)


@pytest.mark.parametrize(
    "length_bp, ng",
    [
        # the numbers a cloning handbook quotes for 20 fmol
        (1000, 12.4),
        (2000, 24.7),
        (5000, 61.8),
        (10000, 123.6),
    ],
)
def test_twenty_femtomoles_weighs_what_the_handbook_says(length_bp, ng):
    assert protocol.ng_for(20, length_bp) == pytest.approx(ng, abs=0.1)


def test_mass_and_moles_are_inverses():
    for length in (700, 2400, 11000):
        ng = protocol.ng_for(20, length)
        assert protocol.fmol_for(ng, length) == pytest.approx(20)


def test_volume_is_mass_over_concentration():
    assert protocol.volume_for(50.0, 100.0) == pytest.approx(0.5)


def test_a_concentration_of_zero_is_refused_rather_than_dividing_by_it():
    with pytest.raises(ValueError):
        protocol.volume_for(50.0, 0.0)


def test_equimolar_means_the_longer_piece_gets_more_mass():
    """The whole point: equal ng would swamp the reaction with short pieces."""
    short = protocol.ng_for(20, 700)
    long = protocol.ng_for(20, 7000)
    assert long == pytest.approx(short * 10, rel=0.01)


# --------------------------------------------------------------------------- #
# the reaction
# --------------------------------------------------------------------------- #


def rx(**overrides):
    body = dict(
        name="pCassette",
        enzyme=BSAI,
        parts=[("promoter", 700, 100.0), ("cds", 2400, 50.0)],
        destination=("backbone", 4000, 200.0),
    )
    body.update(overrides)
    return protocol.reaction(**body)


def test_the_volumes_add_up_to_the_reaction():
    reaction = rx()
    total = sum(c.volume_ul for c in reaction.components)
    assert total == pytest.approx(reaction.total_ul, abs=0.05)


def test_every_piece_of_dna_gets_the_same_number_of_moles():
    reaction = rx()
    dna = [c for c in reaction.components if c.kind in ("part", "destination")]
    assert len(dna) == 3
    assert {c.fmol for c in dna} == {protocol.DEFAULT_FMOL}


def test_a_piece_with_no_concentration_has_no_volume_invented_for_it():
    reaction = rx(parts=[("promoter", 700, 100.0), ("cds", 2400, None)])
    missing = [c for c in reaction.components if c.name == "cds"][0]
    assert missing.volume_ul is None
    assert missing.ng is not None, "the mass is still known - only the volume is not"
    assert not reaction.ok
    assert reaction.missing == ["cds"]
    assert "no recorded concentration" in reaction.issues[0]


def test_the_reagents_are_there_and_scale_with_the_reaction():
    small = rx()
    big = rx(total_ul=20.0)

    def ligase(reaction):
        return [c for c in reaction.components if "ligase (" in c.name][0].volume_ul

    assert ligase(big) == pytest.approx(ligase(small) * 2)


def test_water_makes_up_the_difference():
    reaction = rx()
    water = [c for c in reaction.components if c.kind == "water"][0]
    other = sum(c.volume_ul for c in reaction.components if c.kind != "water")
    assert water.volume_ul == pytest.approx(reaction.total_ul - other, abs=0.02)


def test_dna_that_will_not_fit_is_reported_rather_than_a_negative_water_volume():
    """Dilute preps push the DNA past the reaction volume; that has to be said."""
    reaction = rx(
        parts=[(f"part{i}", 3000, 1.0) for i in range(6)],
        destination=("backbone", 8000, 1.0),
    )
    assert not reaction.ok
    assert any("more than the" in issue for issue in reaction.issues)
    water = [c for c in reaction.components if c.kind == "water"][0]
    assert water.volume_ul >= 0


# --------------------------------------------------------------------------- #
# cycling
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(
    "enzyme, temperature",
    [(BSAI, "37"), (BBSI, "37"), (BSMBI, "42")],
)
def test_the_cutting_step_uses_the_enzyme_s_own_temperature(enzyme, temperature):
    """BsmBI is cut at 42 °C, not 37 - the old UI hard-coded one number."""
    first = protocol.cycling(enzyme)[0]
    assert temperature in first.detail
    assert enzyme.name in first.detail


def test_the_programme_ends_held_cold():
    assert protocol.cycling(BSAI)[-1].detail.startswith("4")


# --------------------------------------------------------------------------- #
# the written protocol
# --------------------------------------------------------------------------- #


def test_the_text_version_carries_every_component_and_step():
    reaction = rx()
    text = protocol.as_text(reaction)
    for component in reaction.components:
        assert component.name in text
    for step in reaction.steps:
        assert step.label in text


def test_the_text_version_says_what_is_missing():
    text = protocol.as_text(rx(parts=[("cds", 2400, None)]))
    assert "?" in text
    assert "Before you set this up" in text


# --------------------------------------------------------------------------- #
# where the concentration is kept
# --------------------------------------------------------------------------- #


@pytest.fixture
def library(tmp_path):
    write_genbank(synth.part_plasmid("2", name="promoter"), tmp_path / "promoter.gb")
    lib = Library(tmp_path)
    lib.scan()
    return lib


def test_a_concentration_survives_a_rescan(library):
    library.set_concentration("promoter.gb", 132.5)
    assert library.get("promoter").conc_ng_ul == 132.5
    library.scan(force=True)
    assert library.get("promoter").conc_ng_ul == 132.5


def test_a_concentration_and_a_hand_assigned_type_do_not_erase_each_other(library):
    """They share one record in overrides.json, so this is worth pinning down."""
    library.set_concentration("promoter.gb", 88.0)
    library.set_override("promoter.gb", "3", reason="checked by sequencing")
    entry = library.get("promoter")
    assert entry.conc_ng_ul == 88.0
    assert entry.call.part_type == "3"

    library.set_concentration("promoter.gb", 99.0)
    entry = library.get("promoter")
    assert entry.call.part_type == "3", "measuring a prep must not blank the type"
    assert entry.conc_ng_ul == 99.0


def test_clearing_a_concentration_leaves_the_type_alone(library):
    library.set_concentration("promoter.gb", 88.0)
    library.set_override("promoter.gb", "3")
    library.set_concentration("promoter.gb", None)
    entry = library.get("promoter")
    assert entry.conc_ng_ul is None
    assert entry.call.part_type == "3"


def test_a_negative_concentration_is_refused(library):
    with pytest.raises(ValueError):
        library.set_concentration("promoter.gb", -5)


def test_an_unmeasured_plasmid_simply_has_none(library):
    assert library.get("promoter").conc_ng_ul is None
