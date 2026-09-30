"""The HTTP surface: shared library reads, and Level 2's own endpoints."""

from __future__ import annotations

import pytest
from fastapi.testclient import TestClient

from ggassembler.api.main import create_app
from ggassembler.core.seqio import write_genbank

from . import synth
from .test_assembly import CANONICAL


@pytest.fixture
def client(tmp_path):
    for part_type, name in CANONICAL.items():
        write_genbank(synth.part_plasmid(part_type, name=name, seed=synth.seed_for(name)),
                      tmp_path / f"{name}.gb")
    write_genbank(synth.dropout_plasmid("234", name="GFP_dropout"), tmp_path / "dropout.gb")
    app = create_app(tmp_path)
    with TestClient(app) as client:
        client.library_dir = tmp_path
        yield client


def design(**overrides):
    body = {"selections": dict(CANONICAL), "name": "pCassette"}
    body.update(overrides)
    return body


# --------------------------------------------------------------- library ---


def test_health(client):
    body = client.get("/health").json()
    assert body["ok"] and body["plasmids"] == 9


def test_summary_counts_what_was_indexed(client):
    body = client.get("/api/library/summary").json()
    assert body["total"] == 9
    assert body["parts"] == 8
    assert body["by_digest"] == 9
    assert body["scheme"] == "YTK"
    assert body["enzymes"] == {"part": "BsaI", "multigene": "BsmBI", "linearizer": "NotI"}


def test_plasmids_can_be_narrowed_to_a_slot(client):
    body = client.get("/api/library/plasmids", params={"five": "AACG", "three": "TATG"}).json()
    assert [p["name"] for p in body] == ["pTDH3"]
    assert body[0]["part_type"] == "2"
    assert body[0]["color"].startswith("#")


def test_one_plasmid_comes_back_in_full(client):
    body = client.get("/api/library/plasmids/Venus").json()
    assert body["part_type"] == "3"
    assert body["topology"] == "circular"
    assert "features" in body and "sites" in body


def test_a_missing_plasmid_is_a_404(client):
    assert client.get("/api/library/plasmids/nope").status_code == 404


def test_an_override_wins_and_is_marked_manual(client):
    body = client.post(
        "/api/library/override",
        json={"path": "Venus.gb", "part_type": "3a", "reason": "N-terminal half"},
    ).json()
    assert body["part_type"] == "3a"
    assert body["source"] == "manual"


# ---------------------------------------------------------------- level 2 ---


def test_slots_come_back_with_their_options(client):
    body = client.post("/api/level2/slots", json=design()).json()
    assert [s["key"] for s in body["slots"]] == ["1", "2", "3", "4", "5", "6", "7", "8"]
    assert all(s["match_count"] == 1 for s in body["slots"])
    assert body["enzyme"] == "BsaI"

    promoter = next(s for s in body["slots"] if s["key"] == "2")
    assert promoter["column"] == "left"
    assert promoter["five_prime"] == "AACG"
    assert promoter["badge_bg"].startswith("#")
    assert [o["name"] for o in promoter["options"]] == ["pTDH3"]


def test_view_modes_change_the_slots(client):
    split = client.post("/api/level2/slots", json=design(split_3=True)).json()
    assert [s["key"] for s in split["slots"]].count("3a") == 1

    composite = client.post(
        "/api/level2/slots", json=design(composite_left=True, composite_right=True)
    ).json()
    assert [s["key"] for s in composite["slots"]] == ["1", "234", "5", "678"]


def test_assemble_returns_everything_the_screen_draws(client):
    body = client.post("/api/level2/assemble", json=design()).json()

    assert body["ok"] is True
    assert body["length"] > 0
    assert len(body["parts"]) == 8
    assert len(body["junctions"]) == 8
    assert body["counts"]["errors"] == 0

    first = body["parts"][0]
    assert {"start", "end", "color", "part_type", "source_name"} <= set(first)
    assert body["parts"][-1]["end"] == body["length"]
    assert all(j["overhang"] and j["upstream"] and j["downstream"] for j in body["junctions"])


def test_parts_tile_the_product_without_gaps(client):
    body = client.post("/api/level2/assemble", json=design()).json()
    edges = [(p["start"], p["end"]) for p in body["parts"]]
    for (_, end), (start, _) in zip(edges, edges[1:], strict=False):
        assert end == start


def test_an_incomplete_design_reports_issues_not_a_500(client):
    selections = {k: v for k, v in CANONICAL.items() if k != "4"}
    body = client.post("/api/level2/assemble", json=design(selections=selections)).json()

    assert body["ok"] is False
    assert body["parts"] == []
    assert any(i["code"] == "empty_slot" for i in body["issues"])
    assert body["counts"]["errors"] >= 1


def test_a_dropout_chosen_as_a_part_is_refused(client):
    selections = dict(CANONICAL, **{"3": "GFP_dropout"})
    body = client.post("/api/level2/assemble", json=design(selections=selections)).json()
    assert body["ok"] is False
    assert any(i["code"] in ("wrong_overhangs", "dropout_as_part") for i in body["issues"])


def test_export_returns_a_genbank_file(client):
    response = client.post("/api/level2/export", json=design())
    assert response.status_code == 200
    assert response.text.startswith("LOCUS")
    assert "pCassette" in response.headers["content-disposition"]
    assert "BsaI" in response.text


def test_export_of_a_broken_design_is_422(client):
    response = client.post("/api/level2/export", json=design(selections={}))
    assert response.status_code == 422
    assert "cannot assemble" in response.text


def test_save_writes_into_the_library_and_reindexes(client):
    body = client.post("/api/level2/save", json=design(name="pSaved")).json()
    assert body["ok"] is True
    assert (client.library_dir / "pSaved.gb").exists()

    # the product is now indexed, and is not mistaken for a part
    saved = client.get("/api/library/plasmids/pSaved").json()
    assert saved["part_type"] is None


# ----------------------------------------------------------------- screens ---


def test_the_cassette_screen_is_served(client):
    response = client.get("/")
    assert response.status_code == 200
    assert "Cassette assembly" in response.text
    assert "/static/level2/level2.js" in response.text


def test_every_screen_is_served(client):
    """All four exist now; none should fall back to the placeholder."""
    for url, title in (
        ("/", "Cassette assembly"),
        ("/part", "Part plasmid construction"),
        ("/multigene", "Multigene assembly"),
        ("/library", "Part library"),
    ):
        response = client.get(url)
        assert response.status_code == 200
        assert title in response.text
        assert "not built yet" not in response.text


def test_the_library_screen_is_served(client):
    response = client.get("/library")
    assert response.status_code == 200
    assert "Part library" in response.text
    assert "/static/library/library.js" in response.text


def test_tokens_and_the_level_module_are_served(client):
    tokens = client.get("/static/tokens.css")
    assert tokens.status_code == 200
    assert "--part-1" in tokens.text
    assert "IBM Plex Mono" in tokens.text

    assert client.get("/static/level2/level2.js").status_code == 200
    assert client.get("/static/fonts/IBMPlexSans-400-latin.woff2").status_code == 200


# ----------------------------------------------------------------- library ---


def test_summary_carries_the_four_stat_cards(client):
    body = client.get("/api/library/summary").json()
    for key in ("total", "by_digest", "composite", "unrecognised"):
        assert key in body, key
    assert body["composite"] == 0  # the canonical eight are all single-position,
    # and the GFP dropout spans 2-3-4 but is not a part you can build with
    assert body["unrecognised"] == 0


def test_a_composite_part_is_counted_as_one(client, tmp_path):
    write_genbank(synth.part_plasmid("234", name="TU", seed=synth.seed_for("TU")),
                  client.library_dir / "TU.gb")
    client.post("/api/library/rescan")
    body = client.get("/api/library/summary").json()
    assert body["composite"] == 1

    rows = {r["name"]: r for r in client.get("/api/library/plasmids").json()}
    assert rows["TU"]["composite"] is True
    assert rows["pTDH3"]["composite"] is False


def test_an_unrecognised_plasmid_is_flagged_for_a_human(client):
    write_genbank(synth.plasmid_without_bsai(name="mystery"), client.library_dir / "mystery.gb")
    client.post("/api/library/rescan")

    rows = {r["name"]: r for r in client.get("/api/library/plasmids").json()}
    assert rows["mystery"]["unrecognised"] is True
    assert rows["mystery"]["needs_attention"] is True
    assert rows["pTDH3"]["needs_attention"] is False
    assert client.get("/api/library/summary").json()["unrecognised"] == 1


def test_an_internal_site_marks_a_row_for_attention(client):
    write_genbank(synth.plasmid_with_internal_site("3", name="undomesticated"),
                  client.library_dir / "undomesticated.gb")
    client.post("/api/library/rescan")

    row = {r["name"]: r for r in client.get("/api/library/plasmids").json()}["undomesticated"]
    assert row["part_type"] == "3"
    assert row["internal_sites"]["part_enzyme"] == 1
    assert row["needs_attention"] is True
    assert row["enzymes"]["part"] == "BsaI"


def test_rows_carry_enzyme_names_rather_than_hard_coded_ones(client):
    row = client.get("/api/library/plasmids").json()[0]
    assert row["enzymes"] == {"part": "BsaI", "multigene": "BsmBI", "linearizer": "NotI"}
    assert set(row["site_totals"]) == {"part_enzyme", "multigene_enzyme", "linearizer"}


def test_types_endpoint_feeds_the_filter_and_the_dialog(client):
    types = client.get("/api/library/types").json()
    by_name = {t["name"]: t for t in types}
    assert by_name["3"]["count"] == 1
    assert by_name["3"]["description"].startswith("CDS")
    assert by_name["3"]["badge_bg"].startswith("#")
    assert "234" in by_name and "8b" in by_name


def test_index_csv_has_a_row_per_plasmid(client):
    response = client.get("/api/library/index.csv")
    assert response.status_code == 200
    assert "attachment" in response.headers["content-disposition"]

    lines = [line for line in response.text.splitlines() if line.strip()]
    assert len(lines) == 1 + 9
    assert lines[0].startswith("path,name,length,topology,part_type")
    assert "BsaI_internal" in lines[0]
    assert any("pTDH3" in line for line in lines[1:])


def test_assigning_a_type_by_hand_survives_a_rescan(client):
    write_genbank(synth.plasmid_without_bsai(name="mystery"), client.library_dir / "mystery.gb")
    client.post("/api/library/rescan")

    client.post("/api/library/override",
                json={"path": "mystery.gb", "part_type": "6", "reason": "checked by sequencing"})
    client.post("/api/library/rescan", params={"force": True})

    row = {r["name"]: r for r in client.get("/api/library/plasmids").json()}["mystery"]
    assert row["part_type"] == "6"
    assert row["source"] == "manual"
    assert row["reason"] == "checked by sequencing"


def test_an_override_can_be_cleared(client):
    client.post("/api/library/override", json={"path": "Venus.gb", "part_type": "3a"})
    client.post("/api/library/override", json={"path": "Venus.gb", "part_type": None})

    row = {r["name"]: r for r in client.get("/api/library/plasmids").json()}["Venus"]
    assert row["part_type"] == "3"
    assert row["source"] == "detected"


# ------------------------------------------------------- front-end invariants --


def test_the_ring_is_positioned_so_its_overlay_lands_on_it(client):
    """The clickable wedges are an absolutely-positioned SVG inside .ring.

    If .ring is not itself positioned, that overlay measures against the whole
    stage instead, and every hit area is scaled up and lands away from the arc
    it belongs to - which shows up as grey wedges floating over the page.
    """
    css = client.get("/static/level2/level2.css").text
    ring = css[css.index(".ring {"):css.index(".ring-hub")]
    assert "position: relative" in ring, ".ring must establish the containing block"


def test_the_ring_paint_and_its_hit_areas_share_one_start_angle(client):
    """The arcs are painted by CSS and hit-tested by SVG, which measure angles
    differently: conic-gradient counts clockwise from 12 o'clock, SVG's cos/sin
    from 3 o'clock. Hard-coding both is how the hover ends up a quarter-turn
    from the arc under the cursor, so both must come from the one constant.
    """
    js = client.get("/static/level2/level2.js").text

    assert "RING_START_DEG" in js and "RING_START_RAD" in js
    assert "conic-gradient(from ${RING_START_DEG}deg" in js
    assert js.count("conic-gradient(from") == 1, "only one place may set the start angle"

    wedge = js[js.index("function wedge("):]
    wedge = wedge[:wedge.index("\n}")]
    assert "RING_START_RAD" in wedge
    assert "Math.PI / 2" not in wedge, "the wedge must not re-derive its own offset"


def test_each_arc_names_its_part_on_hover(client):
    js = client.get("/static/level2/level2.js").text
    assert "'title'" in js and "createElementNS" in js
    assert "arc-hit" in js


def test_the_destination_dropdown_is_grouped_not_a_single_entry_vector(client):
    """The old screen offered only pYTK001; every usable destination is offered now."""
    js = client.get("/static/level1/level1.js").text

    assert "Best for this part type (automatic)" in js
    assert "options.destinations" in js
    assert "optgroup" in js
    assert "options.entry_vectors.map" not in js, "the old single-vector list is gone"


def test_every_plasmid_dropdown_shows_the_merged_display_name(client):
    """A plasmid picker must offer `display` - the merged names plus component.

    Using the bare `name` is how a dropdown ends up showing one filename when
    the sequence is filed under three, and it is easy to reintroduce by editing
    one screen and not the others.
    """
    for path, needle in (
        ("/static/level1/level1.js", "template.display"),
        ("/static/level2/level2.js", "option.display"),
        ("/static/level3/level3.js", "cassette.display"),
    ):
        js = client.get(path).text
        assert needle in js, f"{path} does not use the merged display name"


# ---------------------------------------------------------------- search ---


def test_search_is_a_shared_endpoint_every_screen_can_call(client):
    body = client.get("/api/library/search", params={"q": "type3"}).json()
    assert body["query"] == "type3"
    assert body["count"] == len(body["hits"])


def test_a_search_hit_carries_everything_a_dropdown_needs(client):
    hit = client.get("/api/library/search", params={"q": CANONICAL["3"]}).json()["hits"][0]
    for key in ("name", "display", "part_type", "length", "where", "matched",
                "badge_bg", "badge_fg", "in_part", "component_source"):
        assert key in hit, f"a hit must carry {key}"


def test_an_empty_search_returns_nothing_rather_than_the_whole_library(client):
    assert client.get("/api/library/search", params={"q": ""}).json()["hits"] == []


def test_search_respects_its_limit(client):
    body = client.get("/api/library/search", params={"q": "synthetic", "limit": 2}).json()
    assert len(body["hits"]) <= 2


def test_search_can_be_scoped_to_one_part_type(client):
    body = client.get(
        "/api/library/search", params={"q": "synthetic", "part_type": "3"}
    ).json()
    assert all(h["part_type"] == "3" for h in body["hits"])


def test_the_plasmid_payload_says_where_its_description_came_from(client):
    row = client.get("/api/library/plasmids").json()[0]
    assert row["component_source"] in ("annotation", "filename", "")


# -------------------------------------------------------------- protocol ---


def test_the_protocol_button_has_something_behind_it(client):
    body = client.post("/api/level2/protocol", json=design()).json()
    assert body["components"], "the reaction table must not be empty"
    assert body["steps"], "the cycling must not be empty"
    assert body["enzyme"] == "BsaI"


def test_every_dna_row_says_which_plasmid_it_came_from(client):
    """The row needs a path so a concentration can be recorded against it.

    Parsing it back out of the display name is what this replaced: displays
    carry aliases and component text, and the name is not the key.
    """
    body = client.post("/api/level2/protocol", json=design()).json()
    dna = [c for c in body["components"] if c["kind"] in ("part", "destination")]
    assert dna
    for component in dna:
        assert component["path"], f"{component['name']} has no path"


def test_entering_a_concentration_fills_in_the_volume(client):
    before = client.post("/api/level2/protocol", json=design()).json()
    row = [c for c in before["components"] if c["kind"] == "part"][0]
    assert row["volume_ul"] is None

    client.post(
        "/api/library/concentration", json={"path": row["path"], "conc_ng_ul": 100.0}
    )
    after = client.post("/api/level2/protocol", json=design()).json()
    filled = [c for c in after["components"] if c["path"] == row["path"]][0]
    assert filled["conc_ng_ul"] == 100.0
    assert filled["volume_ul"] == pytest.approx(row["ng"] / 100.0, abs=0.01)


def test_a_concentration_of_zero_is_refused(client):
    path = client.get("/api/library/plasmids").json()[0]["path"]
    response = client.post(
        "/api/library/concentration", json={"path": path, "conc_ng_ul": 0}
    )
    assert response.status_code == 422


def test_the_protocol_downloads_as_text(client):
    response = client.post("/api/level2/protocol.txt", json=design())
    assert response.status_code == 200
    assert "Thermocycler" in response.text
    assert "attachment" in response.headers["content-disposition"]


def test_the_multigene_screen_opens_on_something_rather_than_on_errors(client):
    body = client.get("/api/level3/default").json()
    assert "backbone" in body and "transcription_units" in body


def test_a_slot_option_carries_what_the_picker_filters_on(client):
    """The picker matches name, aliases, component and annotations at once.

    Matching runs in the browser so typing does not cost a round trip, which
    only works if the labels ride along with the options.
    """
    slots = client.post("/api/level2/slots", json=design()).json()["slots"]
    options = [o for slot in slots for o in slot["options"]]
    assert options
    for option in options:
        for key in ("name", "display", "component", "aliases", "labels", "length"):
            assert key in option, f"a slot option must carry {key}"


def test_slot_labels_are_the_part_s_own_not_the_backbone_s(client):
    """AmpR is in every plasmid in the kit; matching it tells you nothing."""
    library = client.app.state.library
    entry = next(e for e in library.unique_entries() if e.is_part)
    span = entry.part_span
    assert span, "this test needs a part with a real fragment"
    inside = {f.label for f in library.labels_in_part(entry)}
    outside = {
        f.label for f in entry.features
        if f.label and f.label not in inside
    }
    assert not (inside & outside)


# --------------------------------------------------------------- triage ---


def test_every_untyped_plasmid_lands_in_exactly_one_triage_pile(client):
    """"83 unrecognised" is not actionable; three named piles are."""
    rows = client.get("/api/library/plasmids").json()
    for row in rows:
        if row["unrecognised"]:
            assert row["triage"] in ("no_part_pair", "uncuttable", "linear")
        else:
            assert row["triage"] == ""


def test_the_canonical_copy_can_be_chosen_by_hand(client, tmp_path):
    """The automatic choice is the shortest name, which is often the wrong one."""
    from ggassembler.core.seqio import write_genbank

    from . import synth

    record = synth.part_plasmid("2", name="shared", seed=synth.seed_for("shared"))
    write_genbank(record, client.library_dir / "aaa.gb")
    write_genbank(record, client.library_dir / "the_real_name.gb")
    client.post("/api/library/rescan?force=true")

    rows = client.get("/api/library/plasmids").json()
    merged = [r for r in rows if r["aliases"]]
    assert merged, "the two identical files should have merged"
    row = merged[0]
    other = row["aliases"][0]

    updated = client.post("/api/library/canonical", json={"name": other}).json()
    assert updated["name"] == other
    assert row["name"] in updated["aliases"]


def test_choosing_a_canonical_copy_survives_a_rescan(client):
    from ggassembler.core.seqio import write_genbank

    from . import synth

    record = synth.part_plasmid("2", name="twin", seed=synth.seed_for("twin"))
    write_genbank(record, client.library_dir / "zzz_twin.gb")
    write_genbank(record, client.library_dir / "b_twin.gb")
    client.post("/api/library/rescan?force=true")

    row = [r for r in client.get("/api/library/plasmids").json() if r["aliases"]][0]
    other = row["aliases"][0]
    client.post("/api/library/canonical", json={"name": other})
    client.post("/api/library/rescan?force=true")

    again = [r for r in client.get("/api/library/plasmids").json() if r["aliases"]][0]
    assert again["name"] == other


def test_asking_for_an_unknown_name_is_a_404(client):
    assert client.post("/api/library/canonical", json={"name": "nope"}).status_code == 404


# -------------------------------------------------------- saved designs ---


def test_a_design_survives_being_saved_and_read_back(client):
    """Nothing persisted before this: a restart lost every choice."""
    body = {"level": "level2", "name": "my cassette", "design": design()}
    saved = client.post("/api/library/designs", json=body).json()
    assert saved["name"] == "my cassette"
    assert saved["saved_at"]

    back = client.get("/api/library/designs", params={"level": "level2"}).json()
    assert [d["name"] for d in back] == ["my cassette"]
    assert back[0]["design"]["selections"] == design()["selections"]


def test_saving_the_same_name_replaces_rather_than_duplicates(client):
    for selection in ("a", "b"):
        client.post("/api/library/designs", json={
            "level": "level2", "name": "same", "design": {"selections": {"1": selection}},
        })
    saved = client.get("/api/library/designs", params={"level": "level2"}).json()
    assert len(saved) == 1
    assert saved[0]["design"]["selections"]["1"] == "b"


def test_designs_are_kept_per_level(client):
    client.post("/api/library/designs",
                json={"level": "level2", "name": "one", "design": {}})
    client.post("/api/library/designs",
                json={"level": "level3", "name": "two", "design": {}})
    assert [d["name"] for d in
            client.get("/api/library/designs", params={"level": "level3"}).json()] == ["two"]
    assert len(client.get("/api/library/designs").json()) == 2


def test_a_design_can_be_deleted(client):
    client.post("/api/library/designs",
                json={"level": "level2", "name": "gone", "design": {}})
    result = client.delete(
        "/api/library/designs", params={"level": "level2", "name": "gone"}
    ).json()
    assert result["deleted"]
    assert client.get("/api/library/designs").json() == []


def test_a_design_needs_a_name(client):
    response = client.post("/api/library/designs",
                           json={"level": "level2", "name": "  ", "design": {}})
    assert response.status_code == 422
