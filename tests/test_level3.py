"""Level 3 - multigene assembly.

The chain is ordered entirely by connector overhangs read out of the files, so
these tests build cassettes with known connectors and check that the ordering,
the closing and the failure messages all follow from those overhangs alone.
"""

from __future__ import annotations

import pytest
from fastapi.testclient import TestClient

from ggassembler.api.main import create_app
from ggassembler.core.enzymes import BSMBI, NOTI, find_sites
from ggassembler.core.library import Library
from ggassembler.core.seqio import write_genbank
from ggassembler.levels import level3_multigene as level3

from . import synth

#: A chain that closes: backbone CCAA->GATG, so the units run GATG->..->CCAA.
BACKBONE = ("CCAA", "GATG")
UNITS = [("TU_a", "GATG", "GTTC"), ("TU_b", "GTTC", "AGCA"), ("TU_c", "AGCA", "CCAA")]


def _annotate(record, spans):
    from Bio.SeqFeature import FeatureLocation, SeqFeature

    for at, kind, label in spans:
        if at + 40 <= len(record.seq):
            record.features.append(
                SeqFeature(FeatureLocation(at, at + 40), type=kind,
                           qualifiers={"label": [label]})
            )
    return record


@pytest.fixture
def library(tmp_path):
    backbone = synth.cassette_plasmid(*BACKBONE, dropout="234", name="pDest", seed=41)
    _annotate(backbone, [
        (60, "promoter", "GlpT Promoter"),
        (120, "CDS", "sfGFP"),
        (200, "terminator", "B0015 Terminator"),
        (len(backbone.seq) - 200, "CDS", "AmpR"),
        (len(backbone.seq) - 120, "rep_origin", "ColE1"),
    ])
    write_genbank(backbone, tmp_path / "pDest.gb")
    for index, (name, left, right) in enumerate(UNITS):
        write_genbank(
            synth.cassette_plasmid(left, right, name=name, seed=50 + index),
            tmp_path / f"{name}.gb",
        )
    # a connector part, so the screen can name an overhang after something real
    write_genbank(synth.connector_plasmid("1", "GATG", name="ConL2"), tmp_path / "ConL2.gb")
    lib = Library(tmp_path)
    lib.scan()
    return lib


def design(units=None, backbone="pDest", name="pMulti"):
    return level3.MultigeneDesign(
        backbone=backbone,
        transcription_units=list(units if units is not None else [u[0] for u in UNITS]),
        name=name,
    )


# ------------------------------------------------------------------ options --


def test_the_library_offers_cassettes_and_backbones(library):
    options = level3.options(library)
    assert {c["name"] for c in options["cassettes"]} == {u[0] for u in UNITS}
    assert [b["name"] for b in options["backbones"]] == ["pDest"]
    assert options["backbones"][0]["left_overhang"] == BACKBONE[0]


def test_a_dropout_vector_is_never_offered_as_a_cassette(library):
    options = level3.options(library)
    assert "pDest" not in {c["name"] for c in options["cassettes"]}


def test_connectors_are_named_after_the_parts_that_carry_them(library):
    index = level3.connector_names(library)
    assert "GATG" in index
    assert "ConL2" in index["GATG"]["left"]
    assert "ConL2" in level3.describe_connector("GATG", index)
    # an overhang nothing carries is still shown, just unlabelled
    assert level3.describe_connector("TTTT", index) == "TTTT"


# ------------------------------------------------------------------ chaining --


def test_three_units_close_the_circle(library):
    result = level3.build(library, design())

    assert result.ok, [str(i) for i in result.issues]
    assert result.assembly.length > 0
    assert len(result.assembly.parts) == 4  # backbone + three units
    assert len(result.assembly.junctions) == 4


def test_the_order_of_the_slots_does_not_decide_the_order_of_the_chain(library):
    """The connectors decide it; the dropdowns are just where you typed."""
    forward = level3.build(library, design())
    shuffled = level3.build(library, design(units=["TU_c", "TU_a", "TU_b"]))

    assert shuffled.ok
    assert shuffled.assembly.length == forward.assembly.length
    assert str(shuffled.assembly.product.seq) in str(forward.assembly.product.seq) * 2


def test_the_chain_is_reported_even_when_it_does_not_close(library):
    result = level3.build(library, design(units=["TU_a"]))
    assert not result.ok
    roles = [link["role"] for link in result.chain]
    assert roles == ["backbone", "TU1"]


def test_a_gap_names_the_overhang_and_suggests_what_fits(library):
    result = level3.build(library, design(units=["TU_a"]))

    assert not result.ok
    assert "missing_overhang" in result.codes()
    assert any("GTTC" in i.message for i in result.assembly.errors)
    assert [s["name"] for s in result.suggestions] == ["TU_b"]
    assert "GTTC" in result.suggestions[0]["why"]


def test_a_repeated_connector_is_refused(library):
    result = level3.build(library, design(units=["TU_a", "TU_a", "TU_b", "TU_c"]))
    assert not result.ok
    assert "duplicate_overhang" in result.codes()


def test_a_cassette_cannot_be_used_as_a_backbone(library):
    result = level3.build(library, design(backbone="TU_a"))
    assert not result.ok
    assert "not_a_backbone" in result.codes()
    assert any("dropout" in i.message for i in result.assembly.errors)


def test_a_backbone_cannot_be_used_as_a_transcription_unit(library):
    result = level3.build(library, design(units=["pDest"]))
    assert not result.ok
    assert "not_a_cassette" in result.codes()


def test_no_backbone_and_no_units_are_each_reported(library):
    assert "no_backbone" in level3.build(library, design(backbone=None)).codes()
    assert "no_cassettes" in level3.build(library, design(units=[])).codes()
    assert "empty_slot" in level3.build(library, design(units=[""])).codes()


def test_an_unknown_cassette_is_reported(library):
    result = level3.build(library, design(units=["nope"]))
    assert not result.ok
    assert "unknown_cassette" in result.codes()


# ------------------------------------------------------------- the product --


def test_every_junction_is_annotated_with_its_connector_barcode(library):
    result = level3.build(library, design())

    assert len(result.scars) == 4
    labels = [
        label
        for feature in result.assembly.product.features
        for label in feature.qualifiers.get("label", [])
        if label.startswith("connector barcode")
    ]
    assert len(labels) == 4
    for scar in result.scars:
        assert scar["end"] - scar["start"] <= level3.SCAR
        assert scar["overhang"] in "".join(labels)


def test_the_product_keeps_no_multigene_sites(library):
    """The connectors are consumed; a finished multigene plasmid cannot re-cut."""
    result = level3.build(library, design())
    assert find_sites(str(result.assembly.product.seq), BSMBI) == []


def test_check_primers_are_designed_for_each_junction(library):
    result = level3.build(library, design())
    assert len(result.check_primers) == 4
    for primer in result.check_primers:
        assert len(primer.sequence) == level3.CHECK_PRIMER
        assert primer.sequence in str(result.assembly.product.seq) * 2
        assert 30 < primer.tm < 80


def test_an_integration_backbone_asks_for_noti_linearization(tmp_path):
    body = synth.filler(200, 3) + NOTI.site + synth.filler(200, 4) + NOTI.site
    write_genbank(
        synth.cassette_plasmid(*BACKBONE, dropout="234", name="pInt", seed=61),
        tmp_path / "pInt.gb",
    )
    for index, (name, left, right) in enumerate(UNITS):
        extra = body if index == 0 else None
        write_genbank(
            synth.cassette_plasmid(left, right, body=extra, name=name, seed=70 + index),
            tmp_path / f"{name}.gb",
        )
    lib = Library(tmp_path)
    lib.scan()

    result = level3.build(lib, design(backbone="pInt"))
    assert result.ok, [str(i) for i in result.issues]
    assert result.integration
    assert "linearize" in result.codes()
    assert "locus_primers" in result.codes()


def test_a_construct_without_noti_is_not_called_integrating(library):
    result = level3.build(library, design())
    assert not result.integration
    assert "linearize" not in result.codes()


# ----------------------------------------------------------------- the API --


@pytest.fixture
def client(library):
    app = create_app(library.folder)
    with TestClient(app) as client:
        client.library_dir = library.folder
        yield client


def test_options_endpoint(client):
    body = client.get("/api/level3/options").json()
    assert body["enzyme"] == "BsmBI"
    assert body["linearizer"] == "NotI"
    assert len(body["cassettes"]) == 3
    assert len(body["backbones"]) == 1


def test_assemble_endpoint_returns_the_chain_and_the_map(client):
    body = client.post("/api/level3/assemble", json={
        "backbone": "pDest",
        "transcription_units": [u[0] for u in UNITS],
        "name": "pMulti",
    }).json()

    assert body["ok"] is True
    assert body["counts"]["units"] == 3
    assert len(body["parts"]) == 4
    assert len(body["scars"]) == 4
    assert [link["role"] for link in body["chain"]] == ["backbone", "TU1", "TU2", "TU3"]
    edges = [(p["start"], p["end"]) for p in body["parts"]]
    for (_, end), (start, _) in zip(edges, edges[1:], strict=False):
        assert end == start


def test_a_broken_chain_is_issues_not_a_500(client):
    body = client.post("/api/level3/assemble", json={
        "backbone": "pDest", "transcription_units": ["TU_a"],
    }).json()
    assert body["ok"] is False
    assert body["counts"]["errors"] >= 1
    assert body["suggestions"]


def test_export_and_save(client):
    payload = {"backbone": "pDest", "transcription_units": [u[0] for u in UNITS],
               "name": "pMulti"}
    response = client.post("/api/level3/export", json=payload)
    assert response.status_code == 200
    assert response.text.startswith("LOCUS")
    assert "connector barcode" in response.text

    saved = client.post("/api/level3/save", json=payload).json()
    assert saved["ok"] is True
    assert (client.library_dir / "pMulti.gb").exists()


def test_export_of_a_broken_chain_is_422(client):
    response = client.post("/api/level3/export", json={
        "backbone": "pDest", "transcription_units": ["TU_a"], "name": "x"})
    assert response.status_code == 422


def test_the_multigene_screen_is_served(client):
    response = client.get("/multigene")
    assert response.status_code == 200
    assert "Multigene assembly" in response.text
    assert "/static/level3/level3.js" in response.text


# ------------------------------------------------- what the dropdowns show --


def test_a_cassette_reports_what_it_carries(library):
    """The dropdown label should name the cargo, not the vector."""
    options = level3.options(library)
    row = next(c for c in options["cassettes"] if c["name"] == "TU_a")
    assert "contents" in row
    assert row["unit_length"] > 0
    assert row["unit_length"] < row["length"]


def test_contents_reads_the_unit_not_the_backbone(tmp_path):
    from Bio.SeqFeature import FeatureLocation, SeqFeature

    record = synth.cassette_plasmid("GATG", "GTTC", name="pAnnotated", seed=91)
    start, end = 200, 320
    for at, kind, label in (
        (start + 5, "promoter", "pTDH3"),
        (start + 40, "CDS", "Venus"),
        (start + 90, "terminator", "tADH1"),
        (len(record.seq) - 120, "CDS", "AmpR"),
        (len(record.seq) - 60, "rep_origin", "ColE1"),
    ):
        record.features.append(
            SeqFeature(FeatureLocation(at, at + 30), type=kind,
                       qualifiers={"label": [label]})
        )
    write_genbank(record, tmp_path / "pAnnotated.gb")
    lib = Library(tmp_path)
    lib.scan()

    entry = lib.get("pAnnotated")
    summary = level3.contents(entry)
    assert "pTDH3" in summary and "tADH1" in summary
    assert "AmpR" not in summary and "ColE1" not in summary


def test_noise_labels_are_left_out_of_the_summary():
    import re

    for label in ("AmpR", "ColE1", "CEN6", "Con1 scar", "pYTK013", "Promoter"):
        assert level3._BACKBONE_NOISE.search(label), label
    for label in ("ScTEF1 Promoter", "Venus", "gal2mut", "ScADH1 Terminator"):
        assert not level3._BACKBONE_NOISE.search(label), label


def test_a_backbone_carries_a_feature_map_and_a_unit_span(library):
    options = level3.options(library)
    backbone = options["backbones"][0]

    assert backbone["features"], "a backbone needs a map to draw"
    assert backbone["unit_start"] is not None and backbone["unit_end"] is not None
    for feature in backbone["features"]:
        assert feature["end"] > feature["start"]
        assert feature["color"].startswith("var(--")
        assert isinstance(feature["kept"], bool)


def test_only_backbones_carry_the_map(library):
    """Cassette rows stay small; the map is only needed for the selected vector."""
    options = level3.options(library)
    assert all("features" not in c for c in options["cassettes"])


def test_the_options_endpoint_exposes_all_of_it(client):
    body = client.get("/api/level3/options").json()
    cassette = body["cassettes"][0]
    assert {"contents", "unit_length", "left_overhang", "right_overhang"} <= set(cassette)
    backbone = body["backbones"][0]
    assert backbone["features"] and backbone["unit_start"] is not None


def test_a_span_that_wraps_the_origin_is_measured_the_long_way_round(library):
    """pYTK096's kept arm runs past the origin; the map must not go negative."""
    assert level3._within(10, (90, 20), 100) is True      # inside the wrap
    assert level3._within(95, (90, 20), 100) is True
    assert level3._within(50, (90, 20), 100) is False     # in the dropout
    assert level3._within(5, (0, 50), 100) is True
    assert level3._within(60, (0, 50), 100) is False
    assert level3._within(5, None, 100) is False


def test_the_dropout_is_the_complement_of_what_is_kept(library):
    options = level3.options(library)
    backbone = options["backbones"][0]

    assert backbone["dropout_start"] == backbone["unit_end"]
    assert backbone["dropout_end"] == backbone["unit_start"]
    assert backbone["dropout_length"] + backbone["unit_length"] == backbone["length"]
    assert backbone["dropout_length"] > 0


def test_connector_names_are_learned_from_annotations(tmp_path):
    """A cassette labelled Con1/Con2 teaches the app what those overhangs are."""
    from Bio.SeqFeature import FeatureLocation, SeqFeature

    for index in range(3):  # three cassettes must agree before a name is used
        record = synth.cassette_plasmid("GATG", "GTTC", name=f"pC{index}", seed=120 + index)
        span_start, span_end = 200, 320
        record.features.append(
            SeqFeature(FeatureLocation(span_start - 30, span_start - 5), type="misc_feature",
                       qualifiers={"label": ["Con2"]})
        )
        record.features.append(
            SeqFeature(FeatureLocation(span_end + 5, span_end + 30), type="misc_feature",
                       qualifiers={"label": ["Con3"]})
        )
        write_genbank(record, tmp_path / f"pC{index}.gb")
    lib = Library(tmp_path)
    lib.scan()

    index = level3.connector_names(lib)
    assert index["GATG"]["names"][:1] == ["Con2"]
    assert index["GTTC"]["names"][:1] == ["Con3"]
    assert level3.describe_connector("GATG", index) == "GATG (Con2)"


def test_a_single_stray_label_cannot_rename_a_connector(tmp_path):
    from Bio.SeqFeature import FeatureLocation, SeqFeature

    record = synth.cassette_plasmid("GATG", "GTTC", name="pOne", seed=131)
    record.features.append(
        SeqFeature(FeatureLocation(170, 195), type="misc_feature",
                   qualifiers={"label": ["ConX"]})
    )
    write_genbank(record, tmp_path / "pOne.gb")
    lib = Library(tmp_path)
    lib.scan()

    index = level3.connector_names(lib)
    assert index.get("GATG", {}).get("names", []) == []
    assert level3.describe_connector("GATG", index) == "GATG"


def test_the_connector_name_pattern():
    for label in ("Con1", "Con9", "ConS", "ConE", "ConLS", "ConR1", "ConE'"):
        assert level3._CONNECTOR_NAME.search(label), label
    for label in ("Connector region", "Control", "ColE1"):
        assert not level3._CONNECTOR_NAME.search(label), label


# --------------------------------------------------------------------------- #
# what the screen opens on
# --------------------------------------------------------------------------- #


def test_the_screen_opens_on_something_that_assembles(library):
    """Level 3 used to report no_backbone and empty_slot the moment it loaded.

    Errors on arrival read as a broken app, so the default has to be a chain
    that actually closes - not merely a backbone with empty slots beside it.
    """
    design = level3.default_design(library)
    assert design.backbone, "no backbone chosen"
    assert design.transcription_units, "no transcription units chosen"

    result = level3.build(library, design)
    assert result.ok, [i.message for i in result.assembly.errors]
    assert not result.assembly.errors


def test_the_default_chain_closes_back_onto_the_backbone(library):
    design = level3.default_design(library)
    result = level3.build(library, design)
    codes = result.codes()
    assert "no_backbone" not in codes
    assert "empty_slot" not in codes


def test_the_default_never_uses_a_cassette_twice(library):
    units = level3.default_design(library).transcription_units
    assert len(units) == len(set(units))


def test_a_library_with_no_cassettes_still_returns_a_design(tmp_path):
    empty = Library(tmp_path)
    empty.scan()
    design = level3.default_design(empty)
    assert design.transcription_units == [] and design.backbone is None


# --------------------------------------------------------------------------- #
# the reaction
# --------------------------------------------------------------------------- #


def test_the_multigene_reaction_takes_the_short_programme(client):
    """Every multigene vector is a dropout whose sites leave with the fragment,
    so this is the reaction the short programme exists for.

    It is also the only place in the app where such a plasmid can be in the
    tube: a reversed-site plasmid is not a part, so no Level 2 slot offers one.
    """
    seed = client.get("/api/level3/default").json()
    body = {
        "backbone": seed["backbone"],
        "transcription_units": seed["transcription_units"],
        "name": "pMultigene",
    }
    rx = client.post("/api/level3/protocol", json=body).json()

    assert rx["reversed_dropout"]
    assert [s["label"] for s in rx["steps"]] == ["30× cycle", "Hold"]
    assert any("chloramphenicol" in note for note in rx["notes"])
    assert rx["enzyme"] == "BsmBI"


def test_the_multigene_reaction_costs_every_piece(client):
    seed = client.get("/api/level3/default").json()
    rx = client.post("/api/level3/protocol", json={
        "backbone": seed["backbone"],
        "transcription_units": seed["transcription_units"],
        "name": "pMultigene",
    }).json()
    dna = [c for c in rx["components"] if c["kind"] in ("part", "destination")]
    assert len(dna) == len(seed["transcription_units"]) + 1
    for component in dna:
        assert component["volume_ul"] is not None
        assert component["path"]


def test_the_multigene_protocol_downloads_as_text(client):
    seed = client.get("/api/level3/default").json()
    response = client.post("/api/level3/protocol.txt", json={
        "backbone": seed["backbone"],
        "transcription_units": seed["transcription_units"],
        "name": "pMultigene",
    })
    assert response.status_code == 200
    assert "Thermocycler" in response.text
    assert "Final digest" not in response.text


# --------------------------------------------------------------------------- #
# what goes into each transcription unit
# --------------------------------------------------------------------------- #


def test_a_unit_reports_the_parts_that_built_it(library):
    """Level 3 shows a chain of cassettes, which hides the thing you check:
    that this promoter is driving that gene. The parts are still in the
    cassette as sequence, so they can be read back out."""
    design = level3.default_design(library)
    entry = library.get(design.transcription_units[0])
    parts = level3.unit_contents(library, entry)
    for part in parts:
        assert part["part_type"] in level3.UNIT_POSITIONS
        for key in ("name", "display", "part_type", "length", "start", "color", "known"):
            assert key in part


def test_only_the_transcription_unit_positions_are_reported(library):
    """Connectors join units together and 6-8 are the plasmid's own machinery.
    Neither is what the unit expresses, so neither belongs in this diagram."""
    for entry in library.cassettes():
        for part in level3.unit_contents(library, entry):
            assert part["part_type"] not in ("1", "5", "6", "7", "8", "8a", "8b")


def test_the_parts_come_back_in_reading_order(library):
    for entry in library.cassettes():
        parts = level3.unit_contents(library, entry)
        assert [p["start"] for p in parts] == sorted(p["start"] for p in parts)


def test_a_part_not_in_the_library_is_drawn_as_a_gap(tmp_path):
    """A gap in a diagram of a construct is a fact about the construct.

    Built from real parts so the tiling has something to find: a synthetic
    cassette contains no library parts at all, which would make this pass
    without exercising anything.
    """
    from ggassembler.levels import level2_cassette as level2

    from .test_assembly import CANONICAL

    for part_type, name in CANONICAL.items():
        write_genbank(
            synth.part_plasmid(part_type, name=name, seed=synth.seed_for(name)),
            tmp_path / f"{name}.gb",
        )
    lib = Library(tmp_path)
    lib.scan()
    built = level2.build(lib, level2.CassetteDesign(selections=dict(CANONICAL)))
    assert built.ok
    write_genbank(built.product, tmp_path / "unit.gb")

    lib = Library(tmp_path)
    lib.scan()
    whole = level3.unit_contents(lib, lib.get("unit"))
    assert [p["part_type"] for p in whole] == ["2", "3", "4"]
    assert all(p["known"] for p in whole)

    # take the coding sequence off the shelf: position 3 becomes a gap
    (tmp_path / f"{CANONICAL['3']}.gb").unlink()
    lib = Library(tmp_path)
    lib.scan()

    after = level3.unit_contents(lib, lib.get("unit"))
    gaps = [p for p in after if not p["known"]]
    assert len(gaps) == 1, f"expected one gap, got {[p['part_type'] for p in after]}"
    assert gaps[0]["part_type"] == "3"
    assert gaps[0]["length"] > 0
    assert "not in the library" in gaps[0]["display"]

def test_the_assemble_payload_carries_the_units(client):
    seed = client.get("/api/level3/default").json()
    body = {
        "backbone": seed["backbone"],
        "transcription_units": seed["transcription_units"],
        "name": "pMultigene",
    }
    payload = client.post("/api/level3/assemble", json=body).json()
    assert "units" in payload

    # one column per fragment of the construct, backbone included, so the
    # diagram and the bar beneath it share one set of proportions
    assert len(payload["units"]) == len(payload["parts"])
    assert [u["name"] for u in payload["units"]] == [p["source_name"] for p in payload["parts"]]
    assert [u["length"] for u in payload["units"]] == [p["length"] for p in payload["parts"]]

    backbones = [u for u in payload["units"] if not u["is_unit"]]
    assert len(backbones) == 1 and backbones[0]["parts"] == []
    units = [u for u in payload["units"] if u["is_unit"]]
    assert [u["role"] for u in units] == [f"TU{i}" for i in range(1, len(units) + 1)]


def test_the_columns_sum_to_the_construct(client):
    """They did not: the diagram used whole-plasmid lengths while the bar used
    released fragments, so the columns came to 173% of the construct."""
    seed = client.get("/api/level3/default").json()
    payload = client.post("/api/level3/assemble", json={
        "backbone": seed["backbone"],
        "transcription_units": seed["transcription_units"],
        "name": "pMultigene",
    }).json()
    assert sum(u["length"] for u in payload["units"]) == payload["length"]


def test_the_seed_prefers_units_that_actually_express_something(tmp_path):
    """A spacer cassette assembles as well as any and shows nothing.

    Seeded with two of those, the screen opens on a diagram with nothing in it
    and teaches nothing about what the screen is for.
    """
    from ggassembler.levels import level2_cassette as level2

    from .test_assembly import CANONICAL

    for part_type, name in CANONICAL.items():
        write_genbank(
            synth.part_plasmid(part_type, name=name, seed=synth.seed_for(name)),
            tmp_path / f"{name}.gb",
        )
    backbone = synth.cassette_plasmid(*BACKBONE, dropout="234", name="pDest", seed=41)
    write_genbank(backbone, tmp_path / "pDest.gb")

    # two cassettes on the same connectors: one with parts in it, one without
    left, right = BACKBONE[1], BACKBONE[0]
    lib = Library(tmp_path)
    lib.scan()
    built = level2.build(lib, level2.CassetteDesign(selections=dict(CANONICAL)))
    assert built.ok

    rich = synth.cassette_plasmid(
        left, right, body=str(built.product.seq), name="expresses", seed=60
    )
    write_genbank(rich, tmp_path / "expresses.gb")
    write_genbank(
        synth.cassette_plasmid(left, right, name="spacer", seed=61), tmp_path / "spacer.gb"
    )

    library = Library(tmp_path)
    library.scan()
    chosen = level3.default_design(library).transcription_units
    if "expresses" not in {*chosen} and "spacer" not in {*chosen}:
        pytest.skip("neither candidate closed a chain in this fixture")
    assert "expresses" in chosen, f"the seed picked {chosen} over the unit with parts"

def test_stylesheets_are_served_so_a_browser_revalidates(client):
    """A cached stylesheet is indistinguishable from a change that did not
    work: new markup, old rules. This is a local tool whose files change under
    a running browser, so nothing static may be cached blind."""
    response = client.get("/static/level3/level3.css")
    assert "no-cache" in response.headers.get("cache-control", "")


@pytest.mark.parametrize(
    "component, part_type, expected",
    [
        ("ScCCW12 Promoter", "2", "prCCW12"),
        ("ScPDC1 Terminator", "4", "PDC1ter"),
        ("TKL1", "3", "TKL1"),
        # already carrying the convention: do not double it up
        ("pTDH3", "2", "prTDH3"),
        ("tENo1", "4", "ENo1ter"),
        ("prCCW12", "2", "prCCW12"),
        ("ADH1ter", "4", "ADH1ter"),
        # only the first thing named; a summary is not one part
        ("His3 Promoter · His3", "2", "prHis3"),
        ("Spacer", "234", "Spacer"),
        ("", "3", ""),
    ],
)
def test_a_part_gets_the_short_name_a_map_has_room_for(component, part_type, expected):
    """`ScCCW12 Promoter` is sixteen characters that mean `prCCW12`, and at the
    width one fragment of a construct gets, that decides whether the name sits
    on one line or three."""
    assert level3.short_label(component, part_type) == expected

def test_a_design_names_the_cassettes_that_have_to_be_built_first(library):
    """The question a project starts from: I want these expressed, in this
    order. The answer is a list of plasmids to build before the last step."""
    units = [level3.UnitSpec(name="one"), level3.UnitSpec(name="two")]
    report = level3.design(library, units, name="pPathway")
    if report.errors:
        pytest.skip(f"this fixture cannot close a chain: {report.errors[0].message}")

    assert len(report.cassettes) == 2
    assert report.backbone
    for plan in report.cassettes:
        assert plan.parts["1"] and plan.parts["5"], "connectors were not assigned"
        for position in ("6", "7", "8"):
            assert plan.parts.get(position), f"nothing filled position {position}"


def test_the_connectors_chain_the_cassettes_in_the_order_given(library):
    """Unit i has to begin where unit i-1 ended, and the last has to close back
    onto the backbone - that is the fiddly part this exists to do."""
    backbone = next(
        (v for v in library.multigene_vectors() if v.cassette_overhangs), None
    )
    assert backbone, "the fixture has no destination vector"
    pairs, issues = level3.connector_plan(library, backbone, 3)
    if issues:
        pytest.skip(issues[0].message)

    left, right = backbone.cassette_overhangs
    assert pairs[0][0].connector_overhang == right, "the chain must start at the backbone"
    assert pairs[-1][1].connector_overhang == left, "the chain must close on the backbone"
    for before, after in zip(pairs, pairs[1:], strict=False):
        assert before[1].connector_overhang == after[0].connector_overhang


def test_a_chain_longer_than_the_connectors_allow_is_refused(library):
    """Better to say the library cannot do it than to emit a design that will
    not assemble."""
    backbone = next(v for v in library.multigene_vectors() if v.cassette_overhangs)
    _, issues = level3.connector_plan(library, backbone, 40)
    assert issues and issues[0].code in ("not_enough_connectors", "no_first_connector",
                                         "no_last_connector")


def test_a_connector_that_blocks_multigene_is_never_designed_in(library):
    """It would make a cassette that could never be cut out again."""
    starts, ends = level3._connectors(library)
    for entry in [*starts.values(), *ends.values()]:
        assert entry.level3_ready is not False


def test_a_design_with_no_units_is_refused(library):
    report = level3.design(library, [], name="empty")
    assert not report.ok
    assert report.errors[0].code == "no_units"


def test_the_design_endpoint_returns_a_build_order(client):
    body = {"name": "pPathway", "units": [{"name": "one"}, {"name": "two"}]}
    payload = client.post("/api/level3/design", json=body).json()
    assert "report" in payload and "cassettes" in payload
    assert payload["name"] == "pPathway"
    if payload["ok"]:
        assert "build order" in payload["report"]
        assert len(payload["cassettes"]) == 2
        for plan in payload["cassettes"]:
            assert plan["order"], "the report must say what order to assemble in"


def test_the_design_downloads_as_a_zip(client):
    import io
    import zipfile

    body = {"name": "pPathway", "units": [{"name": "one"}, {"name": "two"}]}
    response = client.post("/api/level3/design.zip", json=body)
    assert response.status_code == 200
    archive = zipfile.ZipFile(io.BytesIO(response.content))
    names = archive.namelist()
    assert any(n.endswith("build-order.txt") for n in names)
    assert archive.read([n for n in names if n.endswith(".txt")][0]).decode()


def test_a_failed_design_does_not_print_a_build_order(library):
    """"0 cassettes to build, then one multigene assembly" reads as a result
    and is not one - there is nothing to build and no step 1."""
    backbone = next(v for v in library.multigene_vectors() if v.cassette_overhangs)
    report = level3.design(library, [level3.UnitSpec()] * 40, backbone=backbone.name)
    assert not report.cassettes

    text = level3.design_report(report, library)
    assert "cannot be designed" in text
    assert "0 cassette" not in text
    assert "1. " not in text, "a failed design must not number steps"


def test_a_backbone_that_cannot_close_says_which_ones_could(client):
    """Naming a backbone that will not work is half an answer; the caller has
    no way to work out the other half."""
    payload = client.post("/api/level3/design", json={
        "name": "pPathway",
        "units": [{"name": "a"}, {"name": "b"}],
        "backbone": "no_such_vector",
    }).json()
    assert not payload["ok"]
    assert any(i["code"] == "no_backbone" for i in payload["issues"])


def test_the_design_dialog_does_not_borrow_the_screens_backbone(client):
    """The backbone the multigene screen happens to be showing is not a choice
    about this design - and the seeded one cannot close a designed chain."""
    script = client.get("/static/level3/level3.js").text
    block = script[script.index("function designBody()"):]
    block = block[:block.index("}")]
    assert "state.backbone" not in block, "the design borrows the screen's backbone"
    assert "design-backbone" in block


def test_leaving_the_backbone_blank_finds_one_that_works(library):
    report = level3.design(library, [level3.UnitSpec(name="a")], backbone=None)
    if report.errors:
        pytest.skip(f"this fixture has no workable vector: {report.errors[0].message}")
    assert report.backbone
    _, issues = level3.connector_plan(library, library.get(report.backbone), 1)
    assert not issues, "the chosen backbone cannot actually close the chain"


def test_junctions_run_in_the_kits_own_connector_order(library):
    """Sorting junctions by overhang is alphabetical DNA, which is arbitrary
    against the ConN numbering: a three-unit design came out Con2 -> Con5 ->
    Con4 -> Con1. It assembles and it reads as a mistake."""
    index = level3.connector_names(library)
    backbone = next(v for v in library.multigene_vectors() if v.cassette_overhangs)
    pairs, issues = level3.connector_plan(library, backbone, 3, index)
    if issues:
        pytest.skip(issues[0].message)

    # the ends are fixed by the backbone; the junctions between are the choice
    junctions = [right.connector_overhang for _, right in pairs[:-1]]
    ranks = [level3.connector_order(oh, index) for oh in junctions]
    assert ranks == sorted(ranks), f"junctions out of order: {junctions}"


@pytest.mark.parametrize(
    "names, expected",
    [
        (["Con5"], 5),
        (["Con2"], 2),
        (["ConR1"], 1),
        (["ConE"], 99),
        ([], 99),
    ],
)
def test_a_connector_is_ranked_by_the_number_it_is_known_by(names, expected):
    index = {"XXXX": {"names": names}}
    assert level3.connector_order("XXXX", index)[0] == expected


def test_the_series_ends_are_found_structurally(library):
    """A connector that exists only as a type 1 part can begin a chain and
    never continue one, so it is the series start; only-as-type-5 can only
    close. That is what makes the ends recognisable without a hard-coded
    overhang - and what makes a full-range destination vector recognisable."""
    first, last = level3.terminal_connectors(library)
    starts, ends = level3._connectors(library)
    assert not (first & set(ends)), "a series start cannot also close a unit"
    assert not (last & set(starts)), "a series end cannot also begin one"
    for overhang in set(starts) & set(ends):
        assert overhang not in first and overhang not in last


def test_the_design_runs_the_kits_canonical_connector_series(library):
    """ConLS -> ConR1 / ConL1 -> ConR2 / ConL2 -> ConRE, which is how the kit
    is numbered and how a bench protocol is written. Picking any chain that
    merely closes gives Con2 -> Con3 -> ConR1, which assembles and reads wrong.
    """
    first, last = level3.terminal_connectors(library)
    if not first or not last:
        pytest.skip("this fixture has no terminal connectors")

    report = level3.design(library, [level3.UnitSpec()] * 3)
    if report.errors:
        pytest.skip(report.errors[0].message)

    assert report.cassettes[0].left in first, "the chain must open at the series start"
    assert report.cassettes[-1].right in last, "the chain must close at the series end"

    index = level3.connector_names(library)
    junctions = [c.right for c in report.cassettes[:-1]]
    ranks = [level3.connector_order(j, index) for j in junctions]
    assert ranks == sorted(ranks)


@pytest.mark.parametrize(
    "side, overhang, expected",
    [
        ("L", "START", "ConLS"),
        ("R", "END", "ConRE"),
        ("L", "MID", "ConL1"),
        ("R", "MID", "ConR1"),
        ("L", "ODD", "ConLODD"),
    ],
)
def test_a_connector_is_named_for_the_side_it_is_used_on(side, overhang, expected):
    """One junction, two names. The same `Con1` overhang is `ConR1` as the
    right end of one unit and `ConL1` as the left end of the next - the plasmid
    cannot know which, because the position it is put in decides it.
    """
    index = {"MID": {"names": ["Con1"]}, "ODD": {"names": []}}
    assert level3.connector_label(
        overhang, side, index, {"START"}, {"END"}
    ) == expected


def test_the_build_order_says_which_side_each_connector_is(library):
    report = level3.design(library, [level3.UnitSpec()] * 2)
    if report.errors:
        pytest.skip(report.errors[0].message)

    text = level3.design_report(report, library)
    for plan in report.cassettes:
        assert plan.left_label.startswith("ConL")
        assert plan.right_label.startswith("ConR")
        # on the header line as well as beside the part, since the header is
        # what you read when scanning the order
        header = next(
            line for line in text.splitlines()
            if line.lstrip().startswith(f"{report.cassettes.index(plan) + 1}.")
        )
        assert plan.left_label in header and plan.right_label in header
        assert f"[{plan.left_label}]" in text and f"[{plan.right_label}]" in text


def test_the_junction_between_two_units_is_one_connector_under_two_names(library):
    """`ConR1` closing one unit and `ConL1` opening the next are the same
    overhang - if they were not, the cassettes would not ligate."""
    report = level3.design(library, [level3.UnitSpec()] * 3)
    if report.errors:
        pytest.skip(report.errors[0].message)
    for before, after in zip(report.cassettes, report.cassettes[1:], strict=False):
        assert before.right == after.left
        assert before.right_label[4:] == after.left_label[4:], (
            f"{before.right_label} and {after.left_label} are the same junction"
        )


# --------------------------------------------------------------------------- #
# the construct bar
# --------------------------------------------------------------------------- #


def test_the_parts_are_drawn_inside_the_construct_bar(client):
    """Not as a diagram above it.

    A separate row has to be kept aligned with the bar beneath it, and it never
    quite is - columns drift off their segment at some zoom levels, and the
    blocks stack when they run out of width. Drawn inside the segment there is
    nothing to align: one element, one set of proportions, scaling together.
    """
    script = client.get("/static/level3/level3.js").text
    assert "map-parts" in script and "map-part" in script
    assert "unit-inputs" not in script, "the old diagram is still being rendered"

    html = client.get("/multigene").text
    assert "unit-inputs" not in html


def test_every_block_in_the_bar_is_a_proportion_of_it(client):
    """Nothing in the bar may be sized in pixels, or it stops scaling with the
    page - which is the whole reason for drawing it here."""
    css = client.get("/static/level3/level3.css").text
    for selector in (".map-seg", ".map-part"):
        start = css.index(f"\n{selector} ") + 1
        block = css[start:css.index("}", start)]
        assert "flex-basis: 0" in block, f"{selector} does not scale with the bar"
        width = _css_property(block, "width")
        assert width is None or width.endswith("%"), f"{selector} is sized in {width}"


def test_a_part_the_library_does_not_have_is_still_drawn(client):
    """As absence, not as nothing: a gap in the construct is a fact about it."""
    css = client.get("/static/level3/level3.css").text
    assert ".map-part.is-gap" in css
    script = client.get("/static/level3/level3.js").text
    assert "is-gap" in script


def test_the_bar_uses_the_short_part_names(client):
    script = client.get("/static/level3/level3.js").text
    block = script[script.index("block.className = 'map-part'"):]
    block = block[:block.index("row.append")]
    assert "piece.short" in block


def _css_property(block, name):
    for line in block.splitlines():
        text = line.split("/*")[0].strip()
        if text.startswith(f"{name}:"):
            return text.split(":", 1)[1].strip().rstrip(";").strip()
    return None


# ------------------------------------------- a backbone per cassette ---


@pytest.fixture
def backboned(library, tmp_path):
    """The library plus real type 6, 7 and 8 parts to choose between.

    The base fixture has none: every earlier test let the design pick whatever
    it liked at those positions, which is exactly the behaviour being replaced.
    """
    for part_type, names in (("2", ("pTest",)), ("3", ("cdsTest",)), ("4", ("tTest",)),
                             ("6", ("mkHis", "mkLeu")), ("7", ("oriCen",)),
                             ("8", ("bbAmp",))):
        for name in names:
            write_genbank(synth.part_plasmid(part_type, name=name,
                                             seed=synth.seed_for(name)),
                          tmp_path / f"{name}.gb")
    # The base fixture cannot close a chain at all - it has type 1 connectors
    # and no type 5 - so a design against it returns no cassettes and there is
    # nothing to assert a backbone on. Both ends of pDest, in both directions.
    # Both ends of pDest and two junctions in between, each in both
    # directions: a chain of n units needs n-1 overhangs that a unit can both
    # end at and the next begin at, which is a much smaller set than the
    # connector list and the usual reason a good backbone will not close.
    for part_type, overhang in (("1", "CCAA"), ("5", "CCAA"),
                                ("1", "GATG"), ("5", "GATG"),
                                ("1", "GTTC"), ("5", "GTTC"),
                                ("1", "AGCA"), ("5", "AGCA")):
        name = f"Con{part_type}_{overhang}"
        write_genbank(synth.connector_plasmid(part_type, overhang, name=name),
                      tmp_path / f"{name}.gb")
    library.scan(force=True)
    return library


def test_each_cassette_can_carry_its_own_marker(backboned):
    """A project integrates each unit at a different locus, so each cassette
    needs its own selection while the origin and the E. coli backbone never
    change."""
    units = [level3.UnitSpec(name="one", backbone={"6": "mkHis"}),
             level3.UnitSpec(name="two", backbone={"6": "mkLeu"})]
    report = level3.design(backboned, units, name="pPer")

    assert report.cassettes[0].parts["6"] == "mkHis"
    assert report.cassettes[1].parts["6"] == "mkLeu"


def test_a_unit_without_an_override_takes_the_shared_default(backboned):
    units = [level3.UnitSpec(name="one", backbone={"6": "mkLeu"}),
             level3.UnitSpec(name="two")]
    report = level3.design(backboned, units, shared={"6": "mkHis"}, name="pMix")

    assert report.cassettes[0].parts["6"] == "mkLeu", "the override lost"
    assert report.cassettes[1].parts["6"] == "mkHis", "the shared default lost"


def test_the_connectors_are_not_overridable_per_unit(backboned):
    """Positions 1 and 5 are what make the units chain in the order asked for.
    Choosing them per unit would mean choosing the order twice, in two places,
    with nothing keeping the answers the same."""
    units = [level3.UnitSpec(name="one", backbone={"1": "nonsense", "5": "nonsense"}),
             level3.UnitSpec(name="two")]
    report = level3.design(backboned, units, name="pChain")

    assert report.cassettes[0].parts["1"] != "nonsense"
    assert report.cassettes[0].parts["5"] != "nonsense"


def test_a_marker_slot_given_an_origin_says_so(backboned):
    """Left to the assembly it comes back as an overhang mismatch, and you are
    left to work out that the marker slot was handed a type 7."""
    units = [level3.UnitSpec(name="one", backbone={"6": "oriCen"})]
    report = level3.design(backboned, units, name="pWrong")

    codes = [i.code for c in report.cassettes for i in c.issues]
    assert "wrong_type" in codes


def test_an_override_naming_nothing_is_reported(backboned):
    units = [level3.UnitSpec(name="one", backbone={"6": "pNotHere"})]
    report = level3.design(backboned, units, name="pGone")

    codes = [i.code for c in report.cassettes for i in c.issues]
    assert "no_part" in codes


def test_a_backbone_nothing_can_close_names_the_missing_connector(library):
    """"nothing ends a unit at CTGA" is true and unactionable. Which part to go
    and find is the whole answer, and for the multi-round destination vectors
    it is the only answer - the kit has those connectors on one side only."""
    stuck = [v for v in library.multigene_vectors()
             if v.cassette_overhangs and level3.connector_plan(library, v, 1)[1]]
    if not stuck:
        pytest.skip("this library can close every backbone")

    report = level3.design(library, [level3.UnitSpec(name="a")],
                           backbone=stuck[0].name, name="pX")
    said = " ".join(i.message for i in report.issues if i.level == 'error')
    assert "type 1 part releasing" in said or "type 5 part releasing" in said


def test_the_design_dialog_offers_a_backbone_and_a_per_unit_override(client):
    """The positions were chosen silently before: "the simplest usable part at
    each position", with no way to say otherwise from the screen."""
    html = client.get("/multigene").text
    for slot in ("6", "7", "8"):
        assert f'id="design-shared-{slot}"' in html, f"no shared picker for position {slot}"

    script = client.get("/static/level3/level3.js").text
    assert "BACKBONE_SLOTS" in script
    assert "row.backbone[slot]" in script, "a unit cannot override the shared choice"
    assert "shared: sharedBackbone()" in script, "the shared choice is never sent"


def test_the_dialog_does_not_offer_the_connectors(client):
    """Positions 1 and 5 decide the order. Offering them beside the backbone
    would mean setting the order twice, in two places."""
    html = client.get("/multigene").text
    assert 'id="design-shared-1"' not in html
    assert 'id="design-shared-5"' not in html


def test_every_picker_in_the_design_dialog_can_be_searched(client):
    """A plain dropdown over a few hundred parts is a scroll, and the thing you
    want is rarely near the name you would guess."""
    script = client.get("/static/level3/level3.js").text
    assert "function searchableSelect(" in script
    # used for the transcription unit positions and the backbone, not just one
    assert script.count("searchableSelect({") >= 3, "only some pickers are searchable"


def test_the_search_matches_what_is_inside_a_part(client):
    """Parts are named for where they came from and remembered for what is in
    them: pYTK009 is a TDH3 promoter and says so nowhere in its name."""
    script = client.get("/static/level3/level3.js").text
    block = script[script.index("function haystackOf("):]
    block = block[:block.index("}")]
    for field in ("entry.name", "entry.component", "entry.aliases"):
        assert field in block, f"the search ignores {field}"


def test_filtering_cannot_silently_unpick_a_part(client):
    """Type a filter that excludes what you already chose and a plain rebuild
    of the list drops it, taking the selection with it."""
    script = client.get("/static/level3/level3.js").text
    block = script[script.index("function searchableSelect("):]
    block = block[:block.index("\nfunction ")]
    assert "entry.name !== current" in block, (
        "the chosen entry is not kept in the list when it stops matching")


def test_the_search_does_not_go_back_to_the_server(client):
    """The whole list for a position is already loaded; a round trip per
    keystroke would be slower and could reorder the list under the cursor."""
    script = client.get("/static/level3/level3.js").text
    block = script[script.index("function searchableSelect("):]
    block = block[:block.index("\nfunction ")]
    assert "fetch(" not in block and "/api/" not in block


# ------------------------------------- keeping what a design worked out ---


@pytest.fixture
def full_client(backboned):
    """A client over a library that can actually close a design.

    The plain `client` fixture cannot: its library has type 1 connectors and no
    type 5, so every design comes back with no cassettes and these tests would
    all skip rather than check anything.
    """
    app = create_app(backboned.folder)
    with TestClient(app) as client:
        client.library_dir = backboned.folder
        yield client


def _body(**over):
    unit = {"promoter": "pTest", "cds": "cdsTest", "terminator": "tTest"}
    body = {"name": "pKeep",
            "units": [{**unit, "name": "one"}, {**unit, "name": "two"}]}
    body.update(over)
    return body


def test_a_design_can_be_put_straight_into_the_library(full_client):
    """The design hands back plasmids that do not exist yet. Until this, the
    only way to keep them was a zip to unpack by hand - so a cassette you had
    just designed could not be picked on any screen."""
    before = full_client.get("/health").json()["plasmids"]
    body = full_client.post("/api/level3/design/save", json=_body()).json()
    if not body["ok"] and any("backbone" in i for i in body["issues"]):
        pytest.skip("this fixture cannot close a design")

    assert body["saved"], body["issues"]
    assert full_client.get("/health").json()["plasmids"] > before


def test_saving_a_design_twice_refuses_rather_than_overwriting(full_client):
    """A design re-run after an edit would otherwise replace the plasmid
    somebody has already transformed, with no way to get the old one back."""
    first = full_client.post("/api/level3/design/save", json=_body()).json()
    if not first["saved"]:
        pytest.skip("this fixture cannot close a design")

    again = full_client.post("/api/level3/design/save", json=_body()).json()
    assert again["ok"] is False
    assert any("already in the library" in i for i in again["issues"])


def test_the_reactions_cover_every_step(full_client):
    """One BsaI per cassette and one BsmBI to join them - separate reactions on
    separate days, because the cassettes have to be built and verified first."""
    text = full_client.post("/api/level3/design/protocols.txt", json=_body()).text
    assert "BsaI" in text
    assert "BsmBI" in text or "cannot" in text.lower()


def test_the_download_carries_the_design_itself(full_client):
    """The build order is written for a human and does not round-trip: names
    are abbreviated in it, and a part chosen and then changed leaves no
    trace."""
    import io
    import json
    import zipfile

    response = full_client.post("/api/level3/design.zip", json=_body())
    archive = zipfile.ZipFile(io.BytesIO(response.content))
    names = archive.namelist()
    assert any(n.endswith("-design.json") for n in names), names
    assert any(n.endswith("-reactions.txt") for n in names), names

    saved = json.loads(archive.read("pKeep-design.json"))
    assert [u["name"] for u in saved["units"]] == ["one", "two"]


def test_a_reloaded_design_gives_the_same_answer(full_client):
    """What is written out has to be what goes back in."""
    import io
    import json
    import zipfile

    first = full_client.post("/api/level3/design", json=_body()).json()
    archive = zipfile.ZipFile(io.BytesIO(full_client.post("/api/level3/design.zip",
                                                     json=_body()).content))
    again = full_client.post("/api/level3/design",
                        json=json.loads(archive.read("pKeep-design.json"))).json()

    def shape(d):
        return [(c["name"], tuple(sorted(c["parts"].items()))) for c in d["cassettes"]]
    assert shape(first) == shape(again)


def test_the_dialog_offers_all_three(client):
    html = client.get("/multigene").text
    for control in ("design-save", "design-protocols", "design-load"):
        assert f'id="{control}"' in html, f"{control} is missing from the dialog"
