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
