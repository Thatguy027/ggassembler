"""A biological problem is a result, never an exception.

Every case here would be a reasonable thing for a user to try. Each must come
back as a specific `Issue` code with a message naming the overhang or part at
fault, so the UI can point at the panel that caused it.
"""

from __future__ import annotations

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from ggassembler.core.assembly import ERROR, Piece, assemble
from ggassembler.core.enzymes import BSAI
from ggassembler.core.library import Library
from ggassembler.core.seqio import write_genbank
from ggassembler.levels import level2_cassette as level2

from . import synth
from .test_assembly import CANONICAL, cassette_library, pieces_for  # noqa: F401


def piece(name: str, left: str, right: str, body: str = "") -> Piece:
    body = body or synth.filler(60, synth.seed_for(name))
    return Piece(
        record=SeqRecord(Seq(left + body), id=name),
        left_overhang=left,
        right_overhang=right,
        source_name=name,
    )


def test_a_missing_part_reports_the_overhang_that_has_no_partner():
    result = assemble([piece("a", "AACG", "TATG"), piece("b", "TATG", "ATCC")], BSAI)

    assert not result.ok
    assert result.product is None
    assert "missing_overhang" in result.codes()
    assert any("ATCC" in i.message for i in result.errors)
    assert any("AACG" in i.message for i in result.errors)


def test_a_duplicated_overhang_names_both_offenders():
    result = assemble(
        [
            piece("promoter_1", "AACG", "TATG"),
            piece("promoter_2", "AACG", "ATCC"),
            piece("cds", "TATG", "AACG"),
        ],
        BSAI,
    )

    assert not result.ok
    assert "duplicate_overhang" in result.codes()
    message = next(i.message for i in result.errors if i.code == "duplicate_overhang")
    assert "promoter_1" in message and "promoter_2" in message
    assert "AACG" in message


def test_a_three_of_four_near_match_warns_but_still_assembles():
    """AACG and AACT differ at one position: close enough to misligate."""
    result = assemble(
        [
            piece("a", "AACG", "AACT"),
            piece("b", "AACT", "TATG"),
            piece("c", "TATG", "AACG"),
        ],
        BSAI,
    )

    assert result.ok, [str(i) for i in result.issues]
    assert "near_match_3of4" in result.codes()
    warning = next(i for i in result.warnings if i.code == "near_match_3of4")
    assert "AACG" in warning.message and "AACT" in warning.message


def test_a_non_contiguous_three_match_warns():
    """ATCG and ATAG match at positions 1, 2 and 4 - the other misligation mode."""
    result = assemble(
        [
            piece("a", "ATCG", "ATAG"),
            piece("b", "ATAG", "TTCT"),
            piece("c", "TTCT", "ATCG"),
        ],
        BSAI,
    )
    assert result.ok
    assert "near_match_noncontiguous" in result.codes()


def test_a_palindromic_overhang_is_refused():
    result = assemble(
        [piece("a", "AATT", "TATG"), piece("b", "TATG", "AATT")],
        BSAI,
    )
    assert not result.ok
    assert "palindromic_overhang" in result.codes()


def test_an_overhang_of_the_wrong_width_is_refused():
    result = assemble([piece("a", "AAC", "TATG"), piece("b", "TATG", "AAC")], BSAI)
    assert not result.ok
    assert "bad_overhang" in result.codes()
    assert any("4 nt" in i.message for i in result.errors)


def test_two_separate_circles_do_not_count_as_one_product():
    result = assemble(
        [
            piece("a", "AACG", "TATG"),
            piece("b", "TATG", "AACG"),
            piece("c", "GCTG", "TACA"),
            piece("d", "TACA", "GCTG"),
        ],
        BSAI,
    )
    assert not result.ok
    assert "multiple_circles" in result.codes()
    assert any("c" in i.message and "d" in i.message for i in result.errors)


def test_an_internal_site_warns_without_blocking():
    body = synth.filler(40, 3) + BSAI.site + synth.filler(40, 4)
    result = assemble(
        [piece("undomesticated", "AACG", "TATG", body=body), piece("b", "TATG", "AACG")],
        BSAI,
    )
    assert result.ok
    assert "internal_site" in result.codes()
    assert any("undomesticated" in i.message for i in result.warnings)


def test_no_fragments_at_all():
    result = assemble([], BSAI)
    assert not result.ok
    assert "no_parts" in result.codes()


# --------------------------------------------------------------------------- #
# Level 2 refuses bad selections before it ever reaches the simulator
# --------------------------------------------------------------------------- #


def test_an_empty_slot_names_the_slot(cassette_library):  # noqa: F811
    selections = {k: v for k, v in CANONICAL.items() if k != "6"}
    result = level2.build(cassette_library, level2.CassetteDesign(selections=selections))

    assert not result.ok
    assert "empty_slot" in result.codes()
    assert any("TACA -> GAGT" in i.message for i in result.errors)


def test_an_unknown_part_name_is_reported(cassette_library):  # noqa: F811
    selections = dict(CANONICAL, **{"3": "does_not_exist"})
    result = level2.build(cassette_library, level2.CassetteDesign(selections=selections))

    assert not result.ok
    assert "unknown_part" in result.codes()


def test_a_part_in_the_wrong_slot_is_refused(cassette_library):  # noqa: F811
    selections = dict(CANONICAL, **{"3": "pTDH3"})  # a promoter in the CDS slot
    result = level2.build(cassette_library, level2.CassetteDesign(selections=selections))

    assert not result.ok
    assert "wrong_overhangs" in result.codes()
    message = next(i.message for i in result.errors if i.code == "wrong_overhangs")
    assert "pTDH3" in message and "TATG -> ATCC" in message


def test_a_dropout_cannot_be_used_as_a_part(tmp_path):
    write_genbank(synth.dropout_plasmid("234", name="GFP_dropout"), tmp_path / "drop.gb")
    for part_type, name in (("1", "ConLS"), ("5", "ConR1"), ("678", "Backbone")):
        write_genbank(synth.part_plasmid(part_type, name=name, seed=synth.seed_for(name)),
                      tmp_path / f"{name}.gb")
    library = Library(tmp_path)
    library.scan()

    design = level2.CassetteDesign(
        selections={"1": "ConLS", "234": "GFP_dropout", "5": "ConR1", "678": "Backbone"},
        composite_left=True,
        composite_right=True,
    )
    result = level2.build(library, design)

    assert not result.ok
    assert "dropout_as_part" in result.codes()
    assert any("consumed" in i.message for i in result.errors)


def test_nothing_raises_for_any_of_these(cassette_library):  # noqa: F811
    """The point of the whole module: none of this is an exception."""
    broken = [
        {},
        dict(CANONICAL, **{"2": "does_not_exist"}),
        dict(CANONICAL, **{"2": "Venus"}),
        {k: v for k, v in CANONICAL.items() if k != "1"},
    ]
    for selections in broken:
        result = level2.build(cassette_library, level2.CassetteDesign(selections=selections))
        assert not result.ok
        assert result.errors
        assert all(i.level == ERROR for i in result.errors)


# --------------------------------------------------------------------------- #
# ggasm check
# --------------------------------------------------------------------------- #


def test_check_needs_a_baseline_before_it_can_compare(tmp_path, capsys):
    from ggassembler.cli import main

    write_genbank(synth.part_plasmid("2", name="promoter"), tmp_path / "promoter.gb")
    code = main(["check", str(tmp_path), "--baseline", str(tmp_path / "base.json")])
    assert code == 2, "no baseline is neither a pass nor a failure"
    assert "--update" in capsys.readouterr().err


def test_a_clean_library_passes(tmp_path, capsys):
    from ggassembler.cli import main

    write_genbank(synth.part_plasmid("2", name="promoter"), tmp_path / "promoter.gb")
    baseline = tmp_path / "base.json"
    assert main(["check", str(tmp_path), "--baseline", str(baseline), "--update"]) == 0
    assert baseline.exists()
    assert main(["check", str(tmp_path), "--baseline", str(baseline)]) == 0
    assert "nothing new" in capsys.readouterr().out


def test_a_newly_added_problem_fails_the_check(tmp_path, capsys):
    """The whole point of the hook: this commit made things worse."""
    from ggassembler.cli import main

    write_genbank(synth.part_plasmid("2", name="promoter"), tmp_path / "promoter.gb")
    baseline = tmp_path / "base.json"
    main(["check", str(tmp_path), "--baseline", str(baseline), "--update"])

    write_genbank(synth.plasmid_without_bsai(name="mystery"), tmp_path / "mystery.gb")
    assert main(["check", str(tmp_path), "--baseline", str(baseline)]) == 1
    assert "new unrecognised: mystery" in capsys.readouterr().err


def test_a_problem_already_in_the_baseline_does_not_fail(tmp_path, capsys):
    """A library of several hundred files always has some. A hook that fails on
    what was already there is a hook people disable."""
    from ggassembler.cli import main

    write_genbank(synth.part_plasmid("2", name="promoter"), tmp_path / "promoter.gb")
    write_genbank(synth.plasmid_without_bsai(name="mystery"), tmp_path / "mystery.gb")
    baseline = tmp_path / "base.json"
    main(["check", str(tmp_path), "--baseline", str(baseline), "--update"])

    assert main(["check", str(tmp_path), "--baseline", str(baseline)]) == 0


def test_fixing_a_problem_is_reported_but_still_passes(tmp_path, capsys):
    from ggassembler.cli import main

    write_genbank(synth.part_plasmid("2", name="promoter"), tmp_path / "promoter.gb")
    write_genbank(synth.plasmid_without_bsai(name="mystery"), tmp_path / "mystery.gb")
    baseline = tmp_path / "base.json"
    main(["check", str(tmp_path), "--baseline", str(baseline), "--update"])

    (tmp_path / "mystery.gb").unlink()
    assert main(["check", str(tmp_path), "--baseline", str(baseline)]) == 0
    assert "1 fewer" in capsys.readouterr().out
