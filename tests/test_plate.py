"""The standing source plate.

It is a record of a physical object, and everything awkward about it follows
from that: it outlives the library index that was current when it was written,
it is edited by hand over weeks, and the volumes in it are the only warning you
get before a sweep runs a well dry halfway down the plate.
"""

from __future__ import annotations

import json
import re

import pytest
from fastapi.testclient import TestClient

from ggassembler.api.main import create_app
from ggassembler.core.library import Library
from ggassembler.core.seqio import write_genbank
from ggassembler.levels import plate_batch

from . import synth
from .test_assembly import CANONICAL


@pytest.fixture
def library(tmp_path) -> Library:
    for part_type, name in CANONICAL.items():
        write_genbank(
            synth.part_plasmid(part_type, name=name, seed=synth.seed_for(name)),
            tmp_path / f"{name}.gb",
        )
    lib = Library(tmp_path)
    lib.scan()
    return lib


@pytest.fixture
def client(tmp_path) -> TestClient:
    for part_type, name in CANONICAL.items():
        write_genbank(
            synth.part_plasmid(part_type, name=name, seed=synth.seed_for(name)),
            tmp_path / f"{name}.gb",
        )
    app = create_app(tmp_path)
    with TestClient(app) as client:
        client.library_dir = tmp_path
        yield client


def some_part(library: Library) -> str:
    return sorted(CANONICAL.values())[0]


# --------------------------------------------------------------------------- #
# well names
# --------------------------------------------------------------------------- #


def test_a_plate_has_ninety_six_wells_and_no_more():
    assert len(plate_batch.wells()) == 96
    assert plate_batch.wells()[0] == "A1"
    assert plate_batch.wells()[-1] == "H12"


def test_column_major_is_the_order_an_eight_channel_works_in():
    """One aspiration serves a whole column, so filling down the columns is
    what makes a sweep dispensable column-wise rather than well by well."""
    order = plate_batch.wells_column_major()
    assert order[:8] == ["A1", "B1", "C1", "D1", "E1", "F1", "G1", "H1"]
    assert order[8] == "A2"
    assert sorted(order) == sorted(plate_batch.wells())


@pytest.mark.parametrize("given,expected", [
    ("a1", "A1"), ("A01", "A1"), (" h12 ", "H12"), ("B09", "B9"),
])
def test_a_well_is_named_one_way_however_it_was_typed(given, expected):
    assert plate_batch.normalise_well(given) == expected


@pytest.mark.parametrize("bad", ["I1", "A13", "A0", "", "AA1", "1A", "A1.5"])
def test_something_that_is_not_a_well_is_refused(bad):
    with pytest.raises(ValueError):
        plate_batch.normalise_well(bad)


# --------------------------------------------------------------------------- #
# storing it
# --------------------------------------------------------------------------- #


def test_an_unwritten_plate_reads_as_an_empty_one(library):
    """First run: no file, and that is not an error."""
    plate = plate_batch.load(library)
    assert plate.wells == {}
    assert plate.labware == plate_batch.DEFAULT_LABWARE


def test_a_well_survives_a_round_trip(library):
    name = some_part(library)
    plate_batch.set_well(library, "a1", name, conc_ng_ul=52.0, volume_ul=180.0)

    again = plate_batch.load(library)
    assert set(again.wells) == {"A1"}
    assert again.wells["A1"].plasmid == name
    assert again.wells["A1"].conc_ng_ul == 52.0
    assert again.wells["A1"].volume_ul == 180.0


def test_the_file_is_where_the_rest_of_the_state_lives(library):
    plate_batch.set_well(library, "A1", some_part(library), volume_ul=100.0)
    written = json.loads(plate_batch.path_for(library).read_text())
    assert written["wells"]["A1"]["plasmid"] == some_part(library)
    assert plate_batch.path_for(library).parent == library.cache_dir


def test_a_well_can_be_emptied(library):
    plate_batch.set_well(library, "A1", some_part(library), volume_ul=100.0)
    plate_batch.set_well(library, "A1", "")
    assert plate_batch.load(library).wells == {}


def test_a_plasmid_that_is_not_in_the_library_is_refused_on_the_way_in(library):
    with pytest.raises(ValueError, match="pNope"):
        plate_batch.set_well(library, "A1", "pNope")


def test_editing_a_well_keeps_the_numbers_you_did_not_retype(library):
    """An empty box means "leave it", not "set it to zero"."""
    name = some_part(library)
    plate_batch.set_well(library, "A1", name, conc_ng_ul=52.0, volume_ul=180.0)
    plate_batch.set_well(library, "A1", name, volume_ul=120.0)

    well = plate_batch.load(library).wells["A1"]
    assert well.conc_ng_ul == 52.0, "the concentration was not retyped, not deleted"
    assert well.volume_ul == 120.0


def test_putting_a_different_plasmid_in_starts_the_well_over(library):
    """Keeping the old concentration here would be worse than losing it: it
    would be a measurement of something that is no longer in the well."""
    first, second = sorted(CANONICAL.values())[:2]
    plate_batch.set_well(library, "A1", first, conc_ng_ul=52.0, volume_ul=180.0)
    plate_batch.set_well(library, "A1", second)

    well = plate_batch.load(library).wells["A1"]
    assert well.plasmid == second
    assert well.volume_ul == 0.0


def test_a_well_falls_back_to_the_librarys_measured_concentration(library):
    """The number the Library screen records, when the well does not override."""
    name = some_part(library)
    # keyed by the file it came from, not by the display name
    library.set_concentration(library.get(name).path, 33.5)
    plate_batch.set_well(library, "B2", name, volume_ul=50.0)
    assert plate_batch.load(library).wells["B2"].conc_ng_ul == 33.5


def test_a_damaged_well_name_does_not_take_the_whole_plate_down(library):
    """A plate is hand-edited, and one bad key should cost one well."""
    plate_batch.set_well(library, "A1", some_part(library), volume_ul=10.0)
    path = plate_batch.path_for(library)
    raw = json.loads(path.read_text())
    raw["wells"]["ZZ9"] = {"plasmid": "pNope", "volume_ul": 1}
    path.write_text(json.dumps(raw))

    plate = plate_batch.load(library)
    assert set(plate.wells) == {"A1"}


# --------------------------------------------------------------------------- #
# what is wrong with it
# --------------------------------------------------------------------------- #


def codes(issues):
    return {i.code for i in issues}


def test_a_well_naming_a_plasmid_the_library_no_longer_has_is_an_error(library):
    """The plate outlives the index: a plasmid gets renamed or archived and the
    well goes on naming it. Caught on load, not when the tips are already on."""
    plate_batch.set_well(library, "A1", some_part(library), conc_ng_ul=10, volume_ul=100)
    raw = json.loads(plate_batch.path_for(library).read_text())
    raw["wells"]["A1"]["plasmid"] = "pGoneAway"
    plate_batch.path_for(library).write_text(json.dumps(raw))

    issues = plate_batch.check(plate_batch.load(library), library)
    assert "unknown_plasmid" in codes(issues)
    assert any("pGoneAway" in i.message for i in issues)


def test_a_well_with_no_concentration_is_flagged(library):
    plate_batch.set_well(library, "A1", some_part(library), conc_ng_ul=0, volume_ul=100)
    assert "no_concentration" in codes(plate_batch.check(plate_batch.load(library), library))


def test_a_well_at_dead_volume_is_flagged(library):
    plate_batch.set_well(library, "A1", some_part(library), conc_ng_ul=10, volume_ul=4.0)
    assert "at_dead_volume" in codes(plate_batch.check(plate_batch.load(library), library))


def test_the_same_plasmid_in_two_wells_is_worth_saying(library):
    name = some_part(library)
    plate_batch.set_well(library, "A1", name, conc_ng_ul=10, volume_ul=100)
    plate_batch.set_well(library, "B1", name, conc_ng_ul=10, volume_ul=100)
    issues = plate_batch.check(plate_batch.load(library), library)
    assert "duplicated" in codes(issues)


def test_a_healthy_plate_says_nothing(library):
    plate_batch.set_well(library, "A1", some_part(library), conc_ng_ul=52, volume_ul=180)
    assert plate_batch.check(plate_batch.load(library), library) == []


# --------------------------------------------------------------------------- #
# volume
# --------------------------------------------------------------------------- #


def test_a_run_that_fits_is_planned_against_the_wells_it_draws_from(library):
    name = some_part(library)
    plate_batch.set_well(library, "A1", name, conc_ng_ul=52, volume_ul=180)
    draw = plate_batch.plan_draw(plate_batch.load(library), {name: 58.0})
    assert draw.ok
    assert draw.per_well == {"A1": 58.0}


def test_a_run_that_would_take_a_well_below_dead_volume_is_refused(library):
    """This is the failure the whole feature exists to prevent: running dry at
    column 9 leaves a half-pipetted plate and a reaction with no ligase."""
    name = some_part(library)
    plate_batch.set_well(library, "A1", name, conc_ng_ul=52, volume_ul=60.0)

    plate = plate_batch.load(library)          # 60 uL, 5 of which is unreachable
    draw = plate_batch.plan_draw(plate, {name: 58.0})
    assert not draw.ok
    assert draw.short and draw.short[0][0] == "A1"
    assert draw.short[0][3] == pytest.approx(55.0)


def test_a_plasmid_that_is_not_on_the_plate_is_named(library):
    """"Add pYTK009 to a well" is the only useful thing to say here."""
    draw = plate_batch.plan_draw(plate_batch.load(library), {"pYTK009": 10.0})
    assert draw.missing == ["pYTK009"]
    assert not draw.ok


def test_a_draw_is_taken_from_one_well_rather_than_split(library):
    """Splitting across two wells is a different protocol. Quietly generating
    one nobody asked for is worse than saying the well is short."""
    name = some_part(library)
    plate_batch.set_well(library, "A1", name, conc_ng_ul=52, volume_ul=40.0)
    plate_batch.set_well(library, "B1", name, conc_ng_ul=52, volume_ul=200.0)

    draw = plate_batch.plan_draw(plate_batch.load(library), {name: 100.0})
    assert draw.ok
    assert draw.per_well == {"B1": 100.0}, "the first well that can cover it, not the first"


def test_applying_a_draw_decrements_the_plate(library):
    name = some_part(library)
    plate_batch.set_well(library, "A1", name, conc_ng_ul=52, volume_ul=180.0)

    draw = plate_batch.plan_draw(plate_batch.load(library), {name: 58.0})
    plate_batch.apply_draw(library, draw)

    assert plate_batch.load(library).wells["A1"].volume_ul == pytest.approx(122.0)


def test_a_draw_that_does_not_fit_is_never_applied(library):
    """The plate records a physical object; decrementing it for a run that was
    never pipetted is how it stops being one."""
    name = some_part(library)
    plate_batch.set_well(library, "A1", name, conc_ng_ul=52, volume_ul=20.0)
    draw = plate_batch.plan_draw(plate_batch.load(library), {name: 100.0})
    with pytest.raises(ValueError):
        plate_batch.apply_draw(library, draw)
    assert plate_batch.load(library).wells["A1"].volume_ul == 20.0


def test_needs_are_folded_per_plasmid():
    assert plate_batch.needs_from([("a", 1.0), ("b", 2.0), ("a", 0.5)]) == {"a": 1.5, "b": 2.0}


# --------------------------------------------------------------------------- #
# through the API
# --------------------------------------------------------------------------- #


def test_the_screen_gets_an_empty_plate_on_a_fresh_library(client):
    body = client.get("/api/plate/source").json()
    assert body["wells"] == [] and body["ok"] is True
    assert body["dead_volume_ul"] == plate_batch.DEFAULT_DEAD_VOLUME


def test_a_well_can_be_set_and_comes_back_in_the_summary(client):
    name = sorted(CANONICAL.values())[0]
    body = client.post("/api/plate/source/well", json={
        "well": "A1", "plasmid": name, "conc_ng_ul": 52.0, "volume_ul": 180.0,
    }).json()
    assert [w["well"] for w in body["wells"]] == ["A1"]
    assert body["wells"][0]["usable_ul"] == pytest.approx(175.0)
    assert body["wells"][0]["known"] is True
    assert body["wells"][0]["part_type"]


def test_an_unknown_plasmid_is_a_422_naming_it(client):
    answer = client.post("/api/plate/source/well", json={"well": "A1", "plasmid": "pNope"})
    assert answer.status_code == 422
    assert "pNope" in answer.json()["detail"]


def test_a_bad_well_name_is_a_422(client):
    answer = client.post("/api/plate/source/well", json={"well": "Z9", "plasmid": "x"})
    assert answer.status_code == 422


def test_the_plate_settings_can_be_changed(client):
    body = client.post("/api/plate/source/layout", json={
        "id": "SRC-07", "dead_volume_ul": 8.0,
    }).json()
    assert body["id"] == "SRC-07" and body["dead_volume_ul"] == 8.0


def test_a_negative_dead_volume_is_refused(client):
    answer = client.post("/api/plate/source/layout", json={"dead_volume_ul": -1})
    assert answer.status_code == 422


def test_candidates_can_be_narrowed_to_a_position(client):
    everything = client.get("/api/plate/candidates").json()
    promoters = client.get("/api/plate/candidates", params={"part_type": "2"}).json()
    assert everything and promoters
    assert len(promoters) < len(everything)
    assert all(c["part_type"] == "2" for c in promoters)


def test_the_geometry_endpoint_feeds_the_grid(client):
    body = client.get("/api/plate/wells").json()
    assert body["rows"] == list("ABCDEFGH")
    assert body["columns"] == list(range(1, 13))
    assert len(body["row_major"]) == 96 and body["row_major"][1] == "A2"
    assert body["column_major"][1] == "B1"


def test_the_screen_is_served(client):
    assert client.get("/plate").status_code == 200
    assert client.get("/static/plate/plate.js").status_code == 200
    assert client.get("/static/plate/plate.css").status_code == 200


# --------------------------------------------------------------------------- #
# the sweep
# --------------------------------------------------------------------------- #


def base_design(library: Library) -> dict:
    """The canonical eight-part cassette, as the base every well shares."""
    return dict(CANONICAL)


def factor(position: str, *names: str) -> plate_batch.Factor:
    return plate_batch.Factor(position=position, candidates=list(names))


@pytest.fixture
def wide(tmp_path) -> Library:
    """A library with several parts per position, so a sweep has room.

    The plain fixture holds one part of each type, which is enough to assemble
    a cassette and not enough to vary one - a 1x1 plate hides exactly the
    behaviour a sweep is about.
    """
    for part_type, name in CANONICAL.items():
        write_genbank(
            synth.part_plasmid(part_type, name=name, seed=synth.seed_for(name)),
            tmp_path / f"{name}.gb",
        )
    for part_type, count in (("2", 3), ("3", 3)):
        for n in range(count):
            name = f"alt{part_type}_{n}"
            write_genbank(
                synth.part_plasmid(part_type, name=name, seed=synth.seed_for(name)),
                tmp_path / f"{name}.gb",
            )
    lib = Library(tmp_path)
    lib.scan()
    return lib


def parts_for(library: Library, position: str, count: int = 2) -> list[str]:
    """Real parts that fit a position, for use as variants."""
    slot = next(s for s in plate_batch.slots_for(plate_batch.slot_keys(), library.scheme)
                if s.key == position)
    return [e.name for e in library.parts_with_overhangs(*slot.overhangs)][:count]


def two_of(library: Library, position: str) -> list[str]:
    return parts_for(library, position, 2)


# -- the cross product ---------------------------------------------------- #


def test_two_factors_make_the_factorial(library):
    """Row A is promoter 1, column 3 is CDS 3. Position carries identity, which
    is what frees the well fill to carry whether it assembles."""
    design = plate_batch.SweepDesign(
        base=base_design(library),
        row=factor("2", "p1", "p2", "p3"),
        column=factor("3", "c1", "c2"),
    )
    placed = plate_batch.layout(design)
    assert len(placed) == 6
    assert [w for w, *_ in placed] == ["A1", "B1", "C1", "A2", "B2", "C2"]
    assert dict((w, (r, c)) for w, _, _, r, c in placed)["B2"] == ("p2", "c2")


def test_a_full_plate_is_eight_by_twelve(library):
    design = plate_batch.SweepDesign(
        base=base_design(library),
        row=factor("2", *[f"p{n}" for n in range(8)]),
        column=factor("3", *[f"c{n}" for n in range(12)]),
    )
    placed = plate_batch.layout(design)
    assert len(placed) == 96
    assert placed[0][0] == "A1" and placed[-1][0] == "H12"


def test_one_factor_fills_column_major(library):
    """The order an 8-channel works in: one aspiration serves a column."""
    design = plate_batch.SweepDesign(
        base=base_design(library),
        row=factor("2", *[f"p{n}" for n in range(10)]),
    )
    assert [w for w, *_ in plate_batch.layout(design)][:10] == [
        "A1", "B1", "C1", "D1", "E1", "F1", "G1", "H1", "A2", "B2",
    ]


def test_a_partial_plate_stops_where_the_variants_do(library):
    design = plate_batch.SweepDesign(
        base=base_design(library),
        row=factor("2", "p1", "p2", "p3"),
        column=factor("3", "c1", "c2"),
    )
    assert len(plate_batch.layout(design)) == 6, "no phantom wells"


def test_a_one_factor_sweep_can_fill_the_whole_plate(library):
    design = plate_batch.SweepDesign(
        base=base_design(library),
        column=factor("3", *[f"c{n}" for n in range(96)]),
    )
    assert len(plate_batch.layout(design)) == 96


# -- selection-time rejection --------------------------------------------- #


def test_a_candidate_that_does_not_fit_is_rejected_before_assembly(library):
    """Ninety-six identical "wrong overhangs" errors is one fact said ninety-six
    times, and it is a fact about the candidate, not about the plate."""
    cds = two_of(library, "3")[0]
    design = plate_batch.SweepDesign(
        base=base_design(library),
        row=factor("2", cds),          # a type 3 part on the type 2 axis
    )
    issues = plate_batch.check_candidates(library, design)
    assert any(i.code == "bad_candidate" for i in issues)
    assert any("needs" in i.message for i in issues)


def test_a_dropout_is_refused_as_a_candidate(library):
    """It is consumed by the reaction and never ends up in the product."""
    drop = next((e.name for e in library.unique_entries() if e.call.reversed_sites), None)
    if drop is None:
        pytest.skip("no dropout in this fixture library")
    design = plate_batch.SweepDesign(base=base_design(library), row=factor("3", drop))
    assert any("dropout" in i.message for i in plate_batch.check_candidates(library, design))


def test_more_variants_than_the_axis_holds_is_refused(library):
    design = plate_batch.SweepDesign(
        base=base_design(library),
        row=factor("2", *[f"p{n}" for n in range(9)]),
    )
    assert any(i.code == "too_many" for i in plate_batch.check_candidates(library, design))
    wide = plate_batch.SweepDesign(
        base=base_design(library),
        column=factor("3", *[f"c{n}" for n in range(13)]),
    )
    assert any(i.code == "too_many" for i in plate_batch.check_candidates(library, wide))


def test_a_sweep_with_nothing_varying_is_not_a_sweep(library):
    design = plate_batch.SweepDesign(base=base_design(library))
    assert any(i.code == "no_factor" for i in plate_batch.check_candidates(library, design))


# -- building all of them -------------------------------------------------- #


def test_every_well_is_assembled_in_silico(wide):
    """The thing the app can do that nobody does by hand."""
    promoters = two_of(wide, "2")
    cdss = two_of(wide, "3")
    design = plate_batch.SweepDesign(
        base=base_design(wide),
        row=factor("2", *promoters),
        column=factor("3", *cdss),
        name="SW",
    )
    built = plate_batch.build_plate(wide, design)
    assert len(built) == len(promoters) * len(cdss)
    assert all(w.status == "ok" for w in built), [w.issues for w in built if w.issues]
    assert all(w.length > 0 for w in built)
    assert len({w.name for w in built}) == len(built), "every well needs its own name"


def test_a_bad_base_part_fails_each_well_rather_than_raising(wide):
    """One well's problem must never take the other ninety-five with it."""
    base = base_design(wide)
    base.pop("6")
    design = plate_batch.SweepDesign(
        base=base, row=factor("2", *two_of(wide, "2")),
    )
    built = plate_batch.build_plate(wide, design)
    assert built, "it still produced a plate"
    assert all(w.status == "err" for w in built)
    assert all(any("empty" in i.message for i in w.issues) for w in built)


def test_a_failure_shared_along_an_axis_is_named_on_the_axis(wide):
    """A factor is shared along an axis, so a bad part lights up its whole
    column. Eight identical warnings is the wrong report: it is one part."""
    promoters = two_of(wide, "2")
    cdss = two_of(wide, "3")
    base = base_design(wide)
    base.pop("7")                                  # breaks every well equally
    design = plate_batch.SweepDesign(
        base=base, row=factor("2", *promoters), column=factor("3", *cdss),
    )
    built = plate_batch.build_plate(wide, design)
    axes = plate_batch.axis_summary(built, design)
    assert all(a["err"] == a["wells"] for a in axes["columns"])
    assert all("empty_slot" in a["shared"] for a in axes["columns"])


def test_the_axis_summary_names_its_variant(wide):
    promoters = two_of(wide, "2")
    design = plate_batch.SweepDesign(
        base=base_design(wide), row=factor("2", *promoters),
        column=factor("3", *two_of(wide, "3")),
    )
    axes = plate_batch.axis_summary(plate_batch.build_plate(wide, design), design)
    assert [a["variant"] for a in axes["rows"]] == promoters


# -- what goes in the master mix ------------------------------------------- #


def test_the_shared_positions_are_everything_that_does_not_vary(wide):
    """768 transfers is not a protocol. Six of eight positions are identical in
    every well, and those are mixed once rather than pipetted ninety-six times."""
    design = plate_batch.SweepDesign(
        base=base_design(wide), row=factor("2", "a"), column=factor("3", "b"),
    )
    shared = plate_batch.shared_positions(design)
    assert "2" not in shared and "3" not in shared
    assert set(shared) == {"1", "4", "5", "6", "7", "8"}


def test_a_one_factor_sweep_shares_seven_positions(wide):
    design = plate_batch.SweepDesign(base=base_design(wide), row=factor("2", "a"))
    assert len(plate_batch.shared_positions(design)) == 7


def test_part_usage_counts_the_wells_each_plasmid_goes_into(wide):
    promoters = two_of(wide, "2")
    design = plate_batch.SweepDesign(
        base=base_design(wide), row=factor("2", *promoters),
        column=factor("3", *two_of(wide, "3")),
    )
    usage = plate_batch.part_usage(plate_batch.build_plate(wide, design))
    assert usage[CANONICAL["1"]] == 4, "a master-mix part is in every well"
    assert usage[promoters[0]] == 2, "a row variant is in one row"


# -- naming ---------------------------------------------------------------- #


def test_the_naming_pattern_is_followed(library):
    design = plate_batch.SweepDesign(pattern="{base}-{row}-{col}", name="SW")
    assert plate_batch.well_name(design, "SW", "pTDH3", "mRuby2", "A1") == "SW-pTDH3-mRuby2"


def test_an_unused_token_does_not_leave_a_dangling_separator(library):
    """A one-factor sweep has no column variant, and `SW_pTDH3_` is not a name."""
    design = plate_batch.SweepDesign(pattern="{base}_{row}_{col}", name="SW")
    assert plate_batch.well_name(design, "SW", "pTDH3", "", "A1") == "SW_pTDH3"


def test_the_well_token_is_available(library):
    design = plate_batch.SweepDesign(pattern="{base}_{well}", name="SW")
    assert plate_batch.well_name(design, "SW", "", "", "H12") == "SW_H12"


# -- through the API -------------------------------------------------------- #


def test_the_slots_endpoint_offers_positions_with_their_options(client):
    slots = client.get("/api/plate/slots").json()
    assert [s["key"] for s in slots] == ["1", "2", "3", "4", "5", "6", "7", "8"]
    assert all(s["options"] for s in slots)


def test_a_sweep_comes_back_well_by_well(client):
    slots = client.get("/api/plate/slots").json()
    base = {s["key"]: s["options"][0]["name"] for s in slots}
    body = client.post("/api/plate/sweep", json={
        "base": base, "name": "SW",
        "row": {"position": "2", "candidates": [base["2"]]},
    }).json()
    assert body["ok"] is True and body["blocked"] is False
    assert body["counts"]["total"] == 1
    assert body["wells"][0]["well"] == "A1"
    assert body["wells"][0]["length"] > 0
    assert body["shared_positions"] == ["1", "3", "4", "5", "6", "7", "8"]


def test_a_blocked_sweep_builds_nothing_and_says_why(client):
    """Refuse before assembling, not after ninety-six identical failures."""
    slots = client.get("/api/plate/slots").json()
    base = {s["key"]: s["options"][0]["name"] for s in slots}
    body = client.post("/api/plate/sweep", json={
        "base": base, "row": {"position": "2", "candidates": [base["3"]]},
    }).json()
    assert body["blocked"] is True
    assert body["wells"] == []
    assert any(i["code"] == "bad_candidate" for i in body["issues"])


# --------------------------------------------------------------------------- #
# the reaction
# --------------------------------------------------------------------------- #


def stocked(library: Library, names, conc=60.0, volume=150.0) -> None:
    """Put everything a run needs on the source plate."""
    for n, name in enumerate(sorted(set(names))):
        well = f"{plate_batch.ROWS[n % 8]}{n // 8 + 1}"
        plate_batch.set_well(library, well, name, conc_ng_ul=conc, volume_ul=volume)


def a_sweep(wide: Library, rows=2, columns=2) -> tuple:
    # asked for by count, not sliced from a fixed pair - slicing two items to
    # three gives two, and two "different" sweeps that are the same sweep
    promoters = parts_for(wide, "2", rows) or [CANONICAL["2"]]
    cdss = parts_for(wide, "3", columns) or [CANONICAL["3"]]
    assert len(promoters) == rows and len(cdss) == columns, "the library is too small"
    design = plate_batch.SweepDesign(
        base=dict(CANONICAL), name="SW",
        row=factor("2", *promoters), column=factor("3", *cdss),
    )
    built = plate_batch.build_plate(wide, design)
    return design, built, promoters, cdss


def test_the_mix_holds_only_what_every_well_shares(wide):
    """768 transfers is not a protocol. Six of eight positions are identical."""
    design, built, promoters, cdss = a_sweep(wide)
    stocked(wide, list(CANONICAL.values()) + promoters + cdss)

    plan = plate_batch.plan_reaction(wide, design, built, plate_batch.load(wide))
    parts = [c for c in plan.mix if c.kind.startswith("part")]
    assert {c.kind for c in parts} == {f"part {k}" for k in ("1", "4", "5", "6", "7", "8")}
    assert not any(c.kind in ("part 2", "part 3") for c in plan.mix)


def test_the_mix_carries_buffer_ligase_and_the_scheme_s_enzyme(wide):
    design, built, promoters, cdss = a_sweep(wide)
    stocked(wide, list(CANONICAL.values()) + promoters + cdss)
    plan = plate_batch.plan_reaction(wide, design, built, plate_batch.load(wide))
    names = {c.name for c in plan.mix}
    assert "T4 DNA ligase" in names and "T4 DNA ligase buffer" in names
    assert wide.scheme.part_enzyme.name in names


def test_the_mix_total_carries_the_overage(wide):
    """What the tips and the reservoir keep is not optional."""
    design, built, promoters, cdss = a_sweep(wide)
    stocked(wide, list(CANONICAL.values()) + promoters + cdss)
    setup = plate_batch.ReactionSetup(overage=1.2)
    plan = plate_batch.plan_reaction(wide, design, built, plate_batch.load(wide), setup)

    ligase = next(c for c in plan.mix if c.name == "T4 DNA ligase")
    assert ligase.total_ul == pytest.approx(ligase.per_well_ul * plan.wells * 1.2, rel=1e-3)


def test_the_volumes_are_the_shared_fmol_maths_not_a_second_copy(wide):
    """`core.protocol` owns this; a second implementation would drift."""
    from ggassembler.core import protocol

    design, built, promoters, cdss = a_sweep(wide)
    stocked(wide, list(CANONICAL.values()) + promoters + cdss, conc=60.0)
    setup = plate_batch.ReactionSetup(target_fmol=25.0)
    plan = plate_batch.plan_reaction(wide, design, built, plate_batch.load(wide), setup)

    part = next(c for c in plan.mix if c.kind == "part 1")
    expected = protocol.volume_for(protocol.ng_for(25.0, part.length), 60.0)
    # abs, not rel: the stored value is rounded to 0.001 uL, which is finer
    # than anything anyone pipettes and coarser than a 1e-3 relative tolerance
    assert part.per_well_ul == pytest.approx(expected, abs=0.001)


def test_a_part_with_no_concentration_is_costed_but_flagged(wide):
    """The app never presents an assumption as a measurement."""
    design, built, promoters, cdss = a_sweep(wide)
    stocked(wide, list(CANONICAL.values()) + promoters + cdss, conc=0.0)
    plan = plate_batch.plan_reaction(wide, design, built, plate_batch.load(wide))

    assert any(i.code == "assumed_concentration" for i in plan.issues)
    assert any(not c.measured and c.note for c in plan.mix)


def test_a_reaction_that_will_not_fit_the_tube_is_refused(wide):
    design, built, promoters, cdss = a_sweep(wide)
    stocked(wide, list(CANONICAL.values()) + promoters + cdss, conc=1.0)
    setup = plate_batch.ReactionSetup(total_ul=5.0, target_fmol=100.0)
    plan = plate_batch.plan_reaction(wide, design, built, plate_batch.load(wide), setup)
    assert any(i.code == "over_volume" for i in plan.issues)
    assert not plan.ok


def test_water_makes_up_the_difference(wide):
    design, built, promoters, cdss = a_sweep(wide)
    stocked(wide, list(CANONICAL.values()) + promoters + cdss)
    setup = plate_batch.ReactionSetup(total_ul=10.0, part_ul=1.0)
    plan = plate_batch.plan_reaction(wide, design, built, plate_batch.load(wide), setup)

    picked = len(design.factors) * setup.part_ul
    assert plan.mix_per_well_ul + picked == pytest.approx(10.0, abs=0.01)


def test_a_transfer_volume_far_off_the_fmol_target_is_flagged(wide):
    """A fixed per-part volume is what a protocol does; whether it delivers the
    target is a separate question, and the answer is worth saying."""
    design, built, promoters, cdss = a_sweep(wide)
    stocked(wide, list(CANONICAL.values()) + promoters + cdss, conc=500.0)
    plan = plate_batch.plan_reaction(wide, design, built, plate_batch.load(wide))
    assert any(i.code == "off_target_fmol" for i in plan.issues)


def test_every_varying_part_gets_one_transfer_per_well(wide):
    design, built, promoters, cdss = a_sweep(wide)
    stocked(wide, list(CANONICAL.values()) + promoters + cdss)
    plan = plate_batch.plan_reaction(wide, design, built, plate_batch.load(wide))
    assert len(plan.transfers) == len(built) * len(design.factors)
    assert {t.destination for t in plan.transfers} == {w.well for w in built}


def test_a_part_missing_from_the_source_plate_blocks_the_run(wide):
    """"Add pTDH3 to a well" before anything is pipetted, not after."""
    design, built, promoters, cdss = a_sweep(wide)
    stocked(wide, list(CANONICAL.values()))          # the variants are not stocked
    plan = plate_batch.plan_reaction(wide, design, built, plate_batch.load(wide))
    assert not plan.ok
    assert any(i.code == "not_on_plate" for i in plan.issues)


def test_a_well_too_low_to_cover_the_run_blocks_it(wide):
    design, built, promoters, cdss = a_sweep(wide)
    stocked(wide, list(CANONICAL.values()) + promoters + cdss, volume=6.0)
    plan = plate_batch.plan_reaction(wide, design, built, plate_batch.load(wide))
    assert not plan.ok
    assert any(i.code == "well_short" for i in plan.issues)


# --------------------------------------------------------------------------- #
# outputs
# --------------------------------------------------------------------------- #


def test_the_plate_map_has_a_column_per_position(wide):
    design, built, _p, _c = a_sweep(wide)
    csv = plate_batch.plate_map_csv(built, design, "SRC-01")
    header, first, *_ = csv.splitlines()
    assert header.split(",")[:3] == ["plate", "well", "construct"]
    assert "part_1" in header and "part_8" in header
    assert header.endswith("length_bp,status")
    assert first.startswith("SRC-01,A1,")
    assert len(csv.splitlines()) == len(built) + 1


def test_every_clean_well_yields_a_record(wide):
    design, built, _p, _c = a_sweep(wide)
    made = plate_batch.products(wide, design, built)
    assert len(made) == len([w for w in built if w.status != "err"])
    for well, record in made:
        assert len(record.seq) == well.length
        assert record.id


def test_a_failed_well_yields_no_record(wide):
    base = dict(CANONICAL)
    base.pop("6")
    design = plate_batch.SweepDesign(base=base, row=factor("2", *two_of(wide, "2")))
    built = plate_batch.build_plate(wide, design)
    assert plate_batch.products(wide, design, built) == []


def test_each_well_is_logged_with_the_plate_and_well_it_came_from(wide):
    """A colony picked off H7 three weeks later has to lead back to H7."""
    design, built, _p, _c = a_sweep(wide)
    written = plate_batch.log_builds(wide, design, built, "SRC-01")
    assert written == len(built)

    rows = wide.builds()
    assert len(rows) == len(built)
    assert all(row["level"] == "plate" for row in rows)
    assert any("[SRC-01 A1]" in row["name"] for row in rows)
    assert all(row["parts"] for row in rows)


# --------------------------------------------------------------------------- #
# through the API
# --------------------------------------------------------------------------- #


def api_sweep(client) -> dict:
    slots = client.get("/api/plate/slots").json()
    base = {s["key"]: s["options"][0]["name"] for s in slots if s["options"]}
    for n, name in enumerate(sorted(set(base.values()))):
        client.post("/api/plate/source/well", json={
            "well": f"{plate_batch.ROWS[n % 8]}{n // 8 + 1}",
            "plasmid": name, "conc_ng_ul": 60.0, "volume_ul": 150.0,
        })
    return {"base": base, "name": "SW",
            "row": {"position": "2", "candidates": [base["2"]]}}


def test_the_reaction_endpoint_returns_a_recipe_and_transfers(client):
    body = client.post("/api/plate/reaction", json=api_sweep(client)).json()
    assert body["wells"] == 1
    assert body["mix"] and body["transfers"]
    assert body["mix_total_ul"] > body["mix_per_well_ul"], "overage is applied"
    assert body["aspirations"] < body["wells"] * 8 + 1


def test_the_plate_map_downloads_as_csv(client):
    answer = client.post("/api/plate/map.csv", json=api_sweep(client))
    assert answer.status_code == 200
    assert "attachment" in answer.headers["content-disposition"]
    assert answer.text.startswith("plate,well,construct,")


def test_the_products_download_as_one_genbank_file(client):
    answer = client.post("/api/plate/products.gb", json=api_sweep(client))
    assert answer.status_code == 200
    assert answer.text.count("LOCUS") == 1


def test_the_products_also_download_as_a_zip_of_files(client):
    import io
    import zipfile

    answer = client.post("/api/plate/products.zip", json=api_sweep(client))
    assert answer.status_code == 200
    with zipfile.ZipFile(io.BytesIO(answer.content)) as archive:
        names = archive.namelist()
    assert any(n.endswith(".gb") for n in names)
    assert any(n.endswith(".csv") for n in names), "the map travels with the files"


def test_committing_logs_the_wells_and_draws_the_volumes(client):
    design = api_sweep(client)
    before = client.get("/api/plate/source").json()
    body = client.post("/api/plate/commit", json=design).json()

    assert body["ok"] and body["logged"] == 1
    assert len(client.get("/api/library/builds").json()) == 1

    after = {w["well"]: w["volume_ul"] for w in body["plate"]["wells"]}
    drawn = [w["well"] for w in before["wells"] if after[w["well"]] < w["volume_ul"]]
    assert drawn, "nothing was taken out of the plate"


def test_a_run_that_cannot_be_drawn_is_never_committed(client):
    """The plate records a physical object. Decrementing it for a run that was
    not pipetted is how it stops being one."""
    design = api_sweep(client)
    answer = client.post("/api/plate/commit", json={**design, "setup": {"part_ul": 400.0}})
    assert answer.status_code == 422
    assert client.get("/api/library/builds").json() == []


# --------------------------------------------------------------------------- #
# the generated protocol
# --------------------------------------------------------------------------- #
#
# The brief for this file is "test the template, not the output": every run
# ships code that has never been executed, so what has to hold is that the
# runner is the same every time and the data block is well formed.


import ast


def generated(wide: Library, rows=2, columns=3):
    design, built, promoters, cdss = a_sweep(wide, rows=rows, columns=columns)
    stocked(wide, list(CANONICAL.values()) + promoters + cdss)
    plate = plate_batch.load(wide)
    plan = plate_batch.plan_reaction(wide, design, built, plate)
    source = plate_batch.generate_protocol(wide, design, built, plate, plan, version="0.1.0")
    return design, built, plan, source


def literals(source: str) -> dict:
    """Every top-level literal assignment in the generated file."""
    tree = ast.parse(source)
    out = {}
    for node in tree.body:
        if isinstance(node, ast.Assign) and isinstance(node.targets[0], ast.Name):
            try:
                out[node.targets[0].id] = ast.literal_eval(node.value)
            except ValueError:
                out[node.targets[0].id] = None       # not a literal - see below
    return out


def test_the_generated_file_is_valid_python(wide):
    _design, _built, _plan, source = generated(wide)
    ast.parse(source)


def test_metadata_and_requirements_are_literals(wide):
    """The regression test for a bug only the simulator found.

    Opentrons reads `metadata` and `requirements` *statically*, walking the AST
    before the file is ever executed, and refuses anything that is not a
    literal. The first version built them from `PLATE["name"]` and
    `DECK["robot"]`, which is valid Python, reads perfectly well, and is
    rejected with "Could not read the contents of the metadata dict". No
    documentation page says so.
    """
    _design, _built, _plan, source = generated(wide)
    found = literals(source)

    assert found.get("metadata") is not None, "metadata is not a literal dict"
    assert found.get("requirements") is not None, "requirements is not a literal dict"
    assert found["requirements"]["robotType"] == "Flex"
    assert found["requirements"]["apiLevel"]
    assert "apiLevel" not in found["metadata"], "apiLevel belongs in requirements only"
    assert found["metadata"]["protocolName"]


def test_the_data_block_carries_the_whole_run(wide):
    _design, built, plan, source = generated(wide)
    found = literals(source)

    assert len(found["TRANSFERS"]) == len(plan.transfers)
    assert found["PLATE"]["wells"] == [w.well for w in built]
    assert found["PLATE"]["mix_per_well_ul"] == plan.mix_per_well_ul
    assert found["CYCLING"]["cycles"] == 30
    assert found["CYCLING"]["digest_c"] == 37, "BsaI cuts at 37"


def test_the_transfer_volumes_add_up(wide):
    """The arithmetic the protocol is a restatement of."""
    _design, built, plan, source = generated(wide)
    found = literals(source)

    picked = sum(move["volume_ul"] for move in found["TRANSFERS"])
    assert picked == pytest.approx(len(plan.transfers) * plan.setup.part_ul)

    per_well = found["PLATE"]["mix_per_well_ul"]
    assert per_well + plan.setup.part_ul * 2 == pytest.approx(plan.setup.total_ul, abs=0.01)
    assert found["PLATE"]["mix_total_ul"] == pytest.approx(
        per_well * len(built) * plan.setup.overage, rel=1e-3)


def test_every_well_gets_the_master_mix_by_column(wide):
    """Column-wise is what makes it an 8-channel job rather than 96 moves."""
    _design, built, _plan, source = generated(wide)
    found = literals(source)
    assert found["PLATE"]["columns"] == sorted({int(w.well[1:]) for w in built})


def test_the_runner_is_identical_between_two_different_runs(wide, tmp_path):
    """The point of splitting the file: only the data changes, so a diff
    between two plates shows the plates, and the code path is the tested one."""
    _d1, _b1, _p1, first = generated(wide, rows=2, columns=2)
    _d2, _b2, _p2, second = generated(wide, rows=2, columns=3)

    tail = plate_batch.RUNNER.strip()
    assert first.rstrip().endswith(tail)      # the file ends with a newline
    assert second.rstrip().endswith(tail)
    assert first[: first.index(tail)] != second[: second.index(tail)], "the data did change"


def test_the_header_says_where_this_run_came_from(wide):
    """A protocol found on a USB stick in three weeks has to be checkable."""
    design, built, _plan, source = generated(wide)
    head = source[: source.index("metadata")]

    assert "GG Assembler 0.1.0" in head
    assert "SRC-01" in head
    assert plate_batch.source_hash(plate_batch.load(wide)) in head
    assert design.row.candidates[0] in head, "the row factor is named"
    assert design.column.candidates[0] in head, "the column factor is named"
    assert f"Wells          {len(built)}" in head


def test_the_plate_hash_changes_when_the_plate_does(wide):
    """That is the whole point of it: a re-racked plate makes every source well
    in a saved protocol a guess."""
    stocked(wide, list(CANONICAL.values()))
    before = plate_batch.source_hash(plate_batch.load(wide))
    plate_batch.set_well(wide, "H12", CANONICAL["1"], conc_ng_ul=10, volume_ul=50)
    assert plate_batch.source_hash(plate_batch.load(wide)) != before


def test_a_tip_is_never_carried_between_two_plasmids(wide):
    """Cross-contamination here does not look like contamination three days
    later. It looks like a cloning failure."""
    assert 'out.setdefault((move["source"], move["plasmid"]), []).append(move)' \
        in plate_batch.RUNNER
    assert plate_batch.RUNNER.count("single.pick_up_tip()") >= 1
    assert plate_batch.RUNNER.count("single.drop_tip()") == \
        plate_batch.RUNNER.count("single.pick_up_tip()")


def test_the_deck_names_are_the_current_flex_ones(wide):
    """Checked against the Opentrons documentation rather than remembered."""
    deck = plate_batch.DECK
    assert deck["robot"] == "Flex"
    assert deck["thermocycler"] == "thermocyclerModuleV2"
    assert deck["single_pipette"].startswith("flex_1channel")
    assert deck["multi_pipette"].startswith("flex_8channel")
    assert deck["tiprack_labware"].startswith("opentrons_flex_96_tiprack")
    # Flex slots are coordinates, not numbers
    for slot in [deck["source_slot"], deck["reservoir_slot"], deck["trash_slot"],
                 *deck["tiprack_slots"]]:
        assert re.fullmatch(r"[A-D][1-4]", slot), slot
    # the Thermocycler occupies fixed slots and takes no slot argument
    assert "thermocycler_slot" not in deck
    assert "load_module(DECK[\"thermocycler\"])" in plate_batch.RUNNER


# --------------------------------------------------------------------------- #
# the simulate guard
# --------------------------------------------------------------------------- #


def test_a_missing_simulator_is_reported_as_missing_not_as_a_pass(wide, monkeypatch):
    """Treating "could not check" as "checked" is how a guard becomes worse
    than no guard at all."""
    import shutil as shutil_module

    monkeypatch.setattr(shutil_module, "which", lambda name: None)
    result = plate_batch.simulate("print('hello')")

    assert result.ran is False
    assert result.ok is False
    assert result.verdict == "not simulated"
    assert "not installed" in result.detail


def test_a_failing_simulation_carries_the_line_worth_reading(wide, monkeypatch):
    import shutil as shutil_module
    import subprocess

    monkeypatch.setattr(shutil_module, "which", lambda name: "/usr/bin/false")
    monkeypatch.setattr(subprocess, "run", lambda *a, **k: subprocess.CompletedProcess(
        a[0], 1, stdout="", stderr="Traceback\n  line\nMalformedPythonProtocolError: bad dict\n"))

    result = plate_batch.simulate("x = 1")
    assert result.ran and not result.ok
    assert result.verdict == "failed"
    assert "MalformedPythonProtocolError" in result.detail


def test_a_passing_simulation_says_so(wide, monkeypatch):
    import shutil as shutil_module
    import subprocess

    monkeypatch.setattr(shutil_module, "which", lambda name: "/usr/bin/true")
    monkeypatch.setattr(subprocess, "run", lambda *a, **k: subprocess.CompletedProcess(
        a[0], 0, stdout="Aspirating 8.0 uL\n", stderr=""))

    result = plate_batch.simulate("x = 1")
    assert result.ran and result.ok and result.verdict == "passed"
    assert "Aspirating" in result.output


# -- through the API -------------------------------------------------------- #


def test_the_protocol_endpoint_returns_the_source_and_the_verdict(client):
    body = client.post("/api/plate/protocol", json=api_sweep(client)).json()
    assert body["source"].startswith("# Generated by GG Assembler")
    assert body["lines"] > 100
    assert body["simulation"]["verdict"] in {"passed", "failed", "not simulated"}
    # the guard: only a protocol that actually passed may be downloaded
    assert body["downloadable"] is (body["simulation"]["ran"] and body["simulation"]["ok"])


def test_the_download_is_refused_unless_the_simulation_passed(client, monkeypatch):
    import shutil as shutil_module

    monkeypatch.setattr(shutil_module, "which", lambda name: None)
    answer = client.post("/api/plate/protocol.py", json=api_sweep(client))
    assert answer.status_code == 422
    assert "not simulated" in answer.json()["detail"]


def test_the_download_is_offered_when_it_passed(client, monkeypatch):
    monkeypatch.setattr(plate_batch, "simulate",
                        lambda source, timeout=180.0: plate_batch.Simulation(ran=True, ok=True))
    answer = client.post("/api/plate/protocol.py", json=api_sweep(client))
    assert answer.status_code == 200
    assert answer.text.startswith("# Generated by GG Assembler")
    assert "attachment" in answer.headers["content-disposition"]
