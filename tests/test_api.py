"""The HTTP surface: shared library reads, and Level 2's own endpoints."""

from __future__ import annotations

import json
import re

import pytest
from fastapi.testclient import TestClient

from ggassembler.api.main import ENV_FOLDERS, ENV_RECURSIVE, create_app, from_env
from ggassembler.cli import build_parser
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
    assert "badge_bg" not in promoter, "colour is the front end's job now"
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
    assert {"start", "end", "part_type", "source_name"} <= set(first)
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
    assert "badge_bg" not in by_name["3"], "colour is the front end's job now"
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


def test_the_map_is_positioned_so_its_overlay_lands_on_it(client):
    """Every layer of the map is absolutely positioned inside `.map`.

    If `.map` is not itself positioned, they measure against whatever ancestor
    is, and the clickable wedges end up scaled up and floating away from the
    arcs they belong to. This bit once already, when the overlay lived inside
    `.ring` instead.
    """
    css = client.get("/static/level2/level2.css").text
    block = css[css.index(".map {"):css.index(".map {") + 400]
    assert "position: relative" in block, ".map must establish the containing block"
    assert "aspect-ratio" in block, ".map must keep the ratio the geometry assumes"


def test_the_bands_are_masked_to_annuli_that_actually_show(client):
    """A radial-gradient with no size keyword sizes to farthest-corner.

    On a square box that radius is side/2 x sqrt(2), so a stop written as a
    fraction of the intended radius lands outside the circle and masks the
    whole band away - the ring simply disappears. `closest-side` is the radius
    the percentages are written against.
    """
    css = client.get("/static/level2/level2.css").text
    masks = [line for line in css.splitlines() if "radial-gradient" in line]
    assert masks, "the bands are not masked into annuli at all"
    for mask in masks:
        assert "closest-side" in mask, f"unsized radial-gradient: {mask.strip()}"


def test_no_layer_of_the_map_is_resized_in_pixels(client):
    """Each layer is a percentage of `.map` so the whole thing scales together.

    Overriding one in pixels - as the narrow-screen rules used to - moves it off
    the centre the callout geometry is computed around.
    """
    css = client.get("/static/level2/level2.css").text
    for selector in (".ring", ".ring-features", ".ring-hub"):
        for start in _occurrences(css, f"{selector} {{"):
            block = css[start:css.index("}", start)]
            for prop in ("width", "height"):
                value = _property(block, prop)
                if value is not None:
                    assert value.endswith("%"), (
                        f"{selector} sets {prop}: {value}, which will not scale"
                    )


def test_a_right_aligned_callout_still_caps_its_own_width(client):
    """`align-items: flex-end` switches off the width cap, and nothing says so.

    In a column flex container the cross axis is horizontal, so any
    `align-items` other than `stretch` sizes children to their content across
    it. The Type 8 callout's meta line rendered 306px wide inside a 143px box
    and escaped 78px past the left edge of the map - with `overflow: hidden`
    and `white-space: nowrap` both set on it, because neither can do anything
    without a width to overflow.

    `max-width: 100%` on the children restores the cap. It reads like a no-op
    and is not; the right-hand callouts, which are `stretch`, are genuinely
    unaffected by it, which is exactly why it looks removable.
    """
    css = client.get("/static/level2/level2.css").text

    alignment = [
        line for line in css.splitlines()
        if ".callout-text" in line and "align-items" in line and "stretch" not in line
    ]
    if not alignment:
        return  # no longer aligned that way, so the cap is not needed

    for start in _occurrences(css, ".callout-text > * {"):
        if "max-width: 100%" in css[start:css.index("}", start)]:
            break
    else:
        raise AssertionError(
            "`.callout-text` children are aligned to an edge "
            f"({alignment[0].strip()}) but never capped: add "
            "`.callout-text > * { max-width: 100%; }`"
        )


def test_the_callout_lines_still_ask_to_be_truncated(client):
    """The cap only helps alongside the overflow rules it exists to enable.

    Every block for the selector is searched, not the first one found: these
    names are also used inside a narrow-screen media query that sets only a
    font size, and taking that block as *the* rule is how a CSS assertion ends
    up passing or failing on text it was never about.
    """
    css = client.get("/static/level2/level2.css").text
    for selector in (".callout-name", ".callout-meta"):
        blocks = [
            css[start:css.index("}", start)]
            for start in _occurrences(css, f"{selector} {{")
        ]
        assert blocks, f"{selector} is not styled at all"
        assert any("nowrap" in b for b in blocks), f"{selector} would wrap, not truncate"
        assert any("overflow: hidden" in b for b in blocks), f"{selector} would spill"
        assert any("ellipsis" in b for b in blocks), f"{selector} truncates with no sign of it"


def test_no_endpoint_sends_a_literal_colour(client):
    """The other half of the theme rule, checked where it actually matters.

    A hex code in a payload is a colour chosen on the server, and the server
    has no idea which theme the page is wearing. Tokens are fine - `var(--x)`
    re-resolves - so the test is about literals, not about the word colour.

    Walked over whole responses rather than named fields: the failure mode is
    someone adding a new field, and a test that only knew the old names would
    pass through exactly the mistake it exists to catch.
    """
    hex_colour = re.compile(r"^#(?:[0-9a-fA-F]{3}|[0-9a-fA-F]{6})$")

    def literals(value, path="") -> list[str]:
        if isinstance(value, dict):
            return [b for k, v in value.items() for b in literals(v, f"{path}.{k}")]
        if isinstance(value, list):
            return [b for n, v in enumerate(value) for b in literals(v, f"{path}[{n}]")]
        if isinstance(value, str) and hex_colour.match(value.strip()):
            return [f"{path} = {value}"]
        return []

    endpoints = [
        client.get("/api/library/plasmids"),
        client.get("/api/library/types"),
        client.get("/api/library/search", params={"q": "TDH"}),
        client.get("/api/level2/default"),
        client.post("/api/level2/assemble", json=design()),
        client.post("/api/level2/slots", json=design()),
    ]
    offences = []
    for response in endpoints:
        assert response.status_code == 200, response.text[:200]
        offences += [f"{response.url.path}{bad}" for bad in literals(response.json())]

    assert not offences, "literal colours in API payloads:\n" + "\n".join(offences)


def test_a_part_still_says_what_type_it_is(client):
    """Removing the colour only works if what replaces it is there."""
    body = client.post("/api/level2/assemble", json=design()).json()
    assert body["parts"], "nothing was assembled, so this proves nothing"
    for part in body["parts"]:
        assert part.get("part_type"), f"{part['name']} has no part_type to colour by"


def _occurrences(text, needle):
    at, out = text.find(needle), []
    while at != -1:
        out.append(at)
        at = text.find(needle, at + 1)
    return out


def _property(block, name):
    """The value of one declaration, or None.

    Scanned line by line rather than with a lookbehind for `;`: these rules
    carry trailing `/* 300 / 680 */` comments, so the previous line does not
    end in a semicolon and a cleverer regex quietly matches nothing at all.
    """
    for line in block.splitlines():
        text = line.split("/*")[0].strip()
        if text.startswith(f"{name}:"):
            return text.split(":", 1)[1].strip().rstrip(";").strip()
    return None

def test_the_ring_paint_and_its_hit_areas_share_one_start_angle(client):
    """The arcs are painted by CSS and hit-tested by SVG, which measure angles
    differently: conic-gradient counts clockwise from 12 o'clock, SVG's cos/sin
    from 3 o'clock. Hard-coding both is how the hover ends up a quarter-turn
    from the arc under the cursor, so both must come from the one constant.
    """
    js = client.get("/static/level2/level2.js").text

    assert "RING_START_DEG" in js and "RING_START_RAD" in js
    assert "conic-gradient(from ${RING_START_DEG}deg" in js
    # two bands are painted now, the parts and their features. Both may exist;
    # what must not is a second source of truth for where the circle starts.
    starts = re.findall(r"conic-gradient\(from ([^,]+),", js)
    assert starts, "nothing paints the ring"
    assert set(starts) == {"${RING_START_DEG}deg"}, (
        f"a band sets its own start angle: {sorted(set(starts))}"
    )

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
                "in_part", "component_source"):
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


def test_entering_a_concentration_replaces_the_assumed_one(client):
    before = client.post("/api/level2/protocol", json=design()).json()
    row = [c for c in before["components"] if c["kind"] == "part"][0]
    assert not row["measured"], "an unmeasured prep should say so"
    assert row["volume_ul"] is not None, "it should still be pipettable"

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


# ------------------------------------------------------------ part features ---


def test_each_placed_part_carries_its_own_annotations(client):
    """The map's inner band draws these, so they must be per part, not per plasmid."""
    body = client.post("/api/level2/assemble", json=design()).json()
    assert body["ok"], [i["message"] for i in body["issues"]]
    for part in body["parts"]:
        assert "features" in part
        for feature in part["features"]:
            for key in ("label", "kind", "start", "end", "strand"):
                assert key in feature


def test_a_feature_is_reported_in_product_coordinates(client):
    """A part sits at one place in its own plasmid and elsewhere in the product.

    Reporting the source coordinates would put the feature band anywhere but
    over the part it belongs to.
    """
    from Bio.SeqFeature import FeatureLocation, SeqFeature

    from ggassembler.core.seqio import write_genbank

    from . import synth

    name = CANONICAL["3"]
    record = synth.part_plasmid(
        "3", insert=synth.filler(600, 3), name=name, seed=synth.seed_for(name)
    )
    # inside the released fragment: forward site + overhang, then the insert
    start = len(synth.forward_site(synth.BSAI)) + 4
    record.features.append(
        SeqFeature(FeatureLocation(start, start + 200), type="CDS",
                   qualifiers={"label": ["a marker gene"]})
    )
    record.features.append(
        SeqFeature(FeatureLocation(start + 220, start + 320), type="terminator",
                   qualifiers={"label": ["and its terminator"]})
    )
    write_genbank(record, client.library_dir / f"{name}.gb")
    client.post("/api/library/rescan?force=true")

    body = client.post("/api/level2/assemble", json=design()).json()
    part = next(p for p in body["parts"] if p["part_type"] == "3")
    cds = [f for f in part["features"] if f["kind"] == "CDS"]
    assert cds, "the CDS inside the type 3 part was not reported"
    for feature in cds:
        assert part["start"] <= feature["start"] < part["end"], "feature starts outside its part"
        assert feature["end"] <= part["end"], "feature runs past the end of its part"


def test_tiny_features_are_dropped_and_the_list_is_capped(client):
    body = client.post("/api/level2/assemble", json=design()).json()
    for part in body["parts"]:
        assert len(part["features"]) <= 8
        for feature in part["features"]:
            assert feature["end"] - feature["start"] >= 60


def test_a_feature_covering_the_whole_part_is_not_drawn(client):
    """It would restate the part band as a second, paler ring inside it.

    Most annotations in a real library do cover their whole part - a type 3
    plasmid's CDS is the type 3 part - so without this the feature band is a
    washed-out copy of the band outside it rather than a view of what is in it.
    """
    from Bio.SeqFeature import FeatureLocation, SeqFeature

    from ggassembler.core.seqio import write_genbank

    from . import synth

    name = CANONICAL["3"]
    record = synth.part_plasmid("3", name=name, seed=synth.seed_for(name))
    start = len(synth.forward_site(synth.BSAI)) + 4
    record.features.append(
        SeqFeature(FeatureLocation(start, start + 120), type="CDS",
                   qualifiers={"label": ["the whole insert"]})
    )
    write_genbank(record, client.library_dir / f"{name}.gb")
    client.post("/api/library/rescan?force=true")

    body = client.post("/api/level2/assemble", json=design()).json()
    part = next(p for p in body["parts"] if p["part_type"] == "3")
    assert not [f for f in part["features"] if f["label"] == "the whole insert"]


def test_a_tool_generated_label_is_not_drawn_as_a_feature(client):
    """`Benchling translation` names the tool, not the thing, and was being
    drawn on the ring - twice, since it duplicates the real annotation's span."""
    from Bio.SeqFeature import FeatureLocation, SeqFeature

    from ggassembler.core.seqio import write_genbank

    from . import synth

    name = CANONICAL["3"]
    record = synth.part_plasmid("3", name=name, seed=synth.seed_for(name))
    start = len(synth.forward_site(synth.BSAI)) + 4
    record.features.append(
        SeqFeature(FeatureLocation(start, start + 70), type="CDS",
                   qualifiers={"note": ["Benchling translation"]})
    )
    write_genbank(record, client.library_dir / f"{name}.gb")
    client.post("/api/library/rescan?force=true")

    body = client.post("/api/level2/assemble", json=design()).json()
    for part in body["parts"]:
        for feature in part["features"]:
            assert "benchling" not in feature["label"].lower()


def test_a_lone_annotation_spanning_its_part_is_not_drawn_either(client):
    """`ConS` across a 194 bp connector draws an arc the size of the arc above
    it. One annotation covering most of a part restates the part; two or more
    are structure, however large either one is."""
    from Bio.SeqFeature import FeatureLocation, SeqFeature

    from ggassembler.core.seqio import write_genbank

    from . import synth

    name = CANONICAL["3"]
    start = len(synth.forward_site(synth.BSAI)) + 4

    def rebuild(spans):
        record = synth.part_plasmid("3", name=name, seed=synth.seed_for(name))
        for offset, size, label in spans:
            record.features.append(
                SeqFeature(FeatureLocation(start + offset, start + offset + size),
                           type="misc_feature", qualifiers={"label": [label]})
            )
        write_genbank(record, client.library_dir / f"{name}.gb")
        client.post("/api/library/rescan?force=true")
        body = client.post("/api/level2/assemble", json=design()).json()
        return next(p for p in body["parts"] if p["part_type"] == "3")

    alone = rebuild([(0, 100, "most of the part")])
    assert alone["features"] == [], "a lone near-full annotation was drawn"

    # the same big feature, with a second one beside it, is structure
    paired = rebuild([(0, 100, "most of the part"), (0, 62, "and a bit of it")])
    assert len(paired["features"]) == 2


# ------------------------------------------------------- level 3 readiness ---


def test_a_slot_option_says_whether_it_blocks_the_level_above(client):
    slots = client.post("/api/level2/slots", json=design()).json()["slots"]
    for slot in slots:
        for option in slot["options"]:
            assert "level3_ready" in option
            if slot["key"] in ("1", "5"):
                assert option["level3_ready"] in (True, False)
            else:
                assert option["level3_ready"] is None


def test_the_library_row_carries_readiness_too(client):
    rows = client.get("/api/library/plasmids").json()
    assert any(r["level3_ready"] is not None for r in rows), "no connector parts indexed"
    for row in rows:
        if row["part_type"] in ("1", "5"):
            assert row["level3_ready"] in (True, False)


def test_a_cassette_built_on_a_blocked_connector_is_warned_about(client, tmp_path):
    """It assembles. That is the problem: nothing said so until Level 3 failed."""
    from ggassembler.core.enzymes import BSMBI, find_sites
    from ggassembler.core.seqio import write_genbank

    from . import synth

    name = CANONICAL["1"]
    record = synth.part_plasmid("1", insert=synth.filler(120, 42), name=name)
    assert not find_sites(str(record.seq), BSMBI)
    write_genbank(record, client.library_dir / f"{name}.gb")
    client.post("/api/library/rescan?force=true")

    body = client.post("/api/level2/assemble", json=design()).json()
    assert body["ok"], "the cassette still assembles - that is the whole point"
    codes = [i["code"] for i in body["issues"]]
    assert "blocks_multigene" in codes
    warning = next(i for i in body["issues"] if i["code"] == "blocks_multigene")
    assert warning["level"] == "warning"
    assert name in warning["message"]


# ------------------------------------------------------------------ sweep ---


def sweep_body(slot="2", candidates=None):
    return {
        "base_design": design(),
        "slot": slot,
        "candidates": candidates if candidates is not None else [CANONICAL["2"]],
    }


def test_a_sweep_builds_one_construct_per_candidate(client):
    rows = client.post("/api/level2/sweep", json=sweep_body()).json()
    assert rows["count"] == 1
    assert rows["rows"][0]["name"] == CANONICAL["2"]
    assert rows["rows"][0]["ok"]
    assert rows["rows"][0]["length"] > 0


def test_a_sweep_holds_every_other_position_still(client):
    """The point of a titration: one thing varies, everything else does not."""
    base = client.post("/api/level2/assemble", json=design()).json()
    rows = client.post("/api/level2/sweep", json=sweep_body()).json()["rows"]
    assert rows[0]["length"] == base["length"], "swapping in the same part changed the result"


def test_a_candidate_that_cannot_build_is_reported_not_dropped(client):
    """A part that fails is a result. Omitting it looks like it was never tried."""
    body = sweep_body(candidates=[CANONICAL["2"], "no_such_plasmid"])
    payload = client.post("/api/level2/sweep", json=body).json()
    assert payload["count"] == 2
    assert payload["built"] == 1
    failed = next(r for r in payload["rows"] if r["name"] == "no_such_plasmid")
    assert not failed["ok"]
    assert failed["errors"]


def test_a_sweep_reports_the_issues_of_each_candidate(client):
    rows = client.post("/api/level2/sweep", json=sweep_body()).json()["rows"]
    for row in rows:
        assert "errors" in row and "warnings" in row and "codes" in row


def test_writing_a_sweep_gives_a_genbank_each_plus_one_picklist(client):
    import io
    import zipfile

    body = sweep_body()
    response = client.post("/api/level2/sweep.zip", json=body)
    assert response.status_code == 200
    assert response.headers["content-type"] == "application/zip"

    archive = zipfile.ZipFile(io.BytesIO(response.content))
    names = archive.namelist()
    sheets = [n for n in names if n.endswith(".tsv")]
    records = [n for n in names if n.endswith(".gb")]
    assert len(sheets) == 1
    assert len(records) == 1

    sheet = archive.read(sheets[0]).decode()
    assert CANONICAL["2"] in sheet
    assert "assembles" in sheet
    assert archive.read(records[0]).decode().startswith("LOCUS")


def test_the_picklist_records_what_was_held_fixed(client):
    import io
    import zipfile

    archive = zipfile.ZipFile(
        io.BytesIO(client.post("/api/level2/sweep.zip", json=sweep_body()).content)
    )
    sheet = archive.read([n for n in archive.namelist() if n.endswith(".tsv")][0]).decode()
    assert "held fixed" in sheet
    assert CANONICAL["3"] in sheet, "the parts that did not vary should be named"


def test_an_empty_sweep_is_not_an_error(client):
    payload = client.post("/api/level2/sweep", json=sweep_body(candidates=[])).json()
    assert payload["count"] == 0 and payload["rows"] == []


# -------------------------------------------------------------- build log ---


def test_an_export_is_recorded_with_every_part_that_went_into_it(client):
    """The question this answers comes months later: which pYTK009 was in the
    thing I built in March, and has that file changed since. So a name is not
    enough - each part is logged with the checksum of its sequence."""
    assert client.get("/api/library/builds").json() == []

    response = client.post("/api/level2/export", json=design())
    assert response.status_code == 200

    log = client.get("/api/library/builds").json()
    assert len(log) == 1
    entry = log[0]
    assert entry["name"] == "pCassette"
    assert entry["level"] == "level2"
    assert entry["length"] > 0
    assert entry["at"]
    assert len(entry["parts"]) == len(design()["selections"])
    for part in entry["parts"]:
        assert part["sha1"], f"{part['name']} logged without a checksum"
        assert part["path"] and part["part_type"]


def test_the_log_is_append_only_and_newest_first(client):
    for name in ("first", "second", "third"):
        client.post("/api/level2/export", json=design(name=name))
    log = client.get("/api/library/builds").json()
    assert [e["name"] for e in log] == ["third", "second", "first"]


def test_a_failed_export_is_not_logged(client):
    """Nothing was built, so nothing was built."""
    client.post("/api/level2/export", json=design(selections={"1": "nope"}))
    assert client.get("/api/library/builds").json() == []


def test_the_log_records_the_issue_codes_the_build_carried(client):
    """A construct that built with a caveat should say so in the log.

    The fixtures' type 1 part carries no BsmBI site, so every cassette built
    from them warns that it could never be released for a Level 3 assembly -
    which is exactly the kind of thing worth finding in a build log later.
    """
    assembled = client.post("/api/level2/assemble", json=design()).json()
    expected = sorted({i["code"] for i in assembled["issues"] if i["level"] != "info"})
    assert expected, "this fixture no longer produces any warning to log"

    client.post("/api/level2/export", json=design())
    log = client.get("/api/library/builds").json()
    assert log[0]["issues"] == expected

def test_a_torn_line_does_not_hide_the_rest_of_the_log(client):
    client.post("/api/level2/export", json=design())
    library = client.app.state.library
    with library.builds_path.open("a", encoding="utf-8") as handle:
        handle.write('{"at": "truncated\n')
    client.post("/api/level2/export", json=design(name="after"))

    log = client.get("/api/library/builds").json()
    assert [e["name"] for e in log] == ["after", "pCassette"]


# ------------------------------------------------------------- map layout ---


def test_a_hidden_element_is_actually_hidden(client):
    """`[hidden]` in the UA sheet is a plain author-level rule, so any class
    that sets `display` beats it. `.clone-report { display: flex }` did, and
    left an empty teal bar on the page on every load."""
    tokens = client.get("/static/tokens.css").text
    rule = re.search(r"^\[hidden\]\s*\{([^}]*)\}", tokens, re.M)
    assert rule, "nothing anywhere makes the hidden attribute work"
    assert "display: none" in rule.group(1)
    assert "!important" in rule.group(1), "a class that sets display would beat it"


def test_every_panel_that_sets_display_can_still_be_hidden(client):
    """Anything toggled with .hidden in script must not out-specify [hidden]."""
    script = client.get("/static/level2/level2.js").text
    tokens = client.get("/static/tokens.css").text

    toggled = set(re.findall(r"el\('([a-z-]+)'\)\.hidden", script))
    toggled |= set(re.findall(r"(\w+)\.hidden = ", script))
    assert toggled, "nothing toggles hidden any more - is this test still needed?"
    # the global rule carries !important, so no display rule can out-specify it
    assert "display: none !important" in tokens


def test_left_hand_callouts_are_anchored_by_their_right_edge(client):
    """Positioned by `left` and pulled back with a transform, they escaped the
    container: a transform moves the box after layout, so shrink-to-fit never
    learns it has less room."""
    script = client.get("/static/level2/level2.js").text
    assert "box.style.right" in script, "left-hand callouts are not right-anchored"
    css = client.get("/static/level2/level2.css").text
    left_rule = css[css.index(".callout.is-left {"):]
    left_rule = left_rule[:left_rule.index("}")]
    assert "translateX(-100%)" not in left_rule


def test_the_callout_width_scales_with_the_map(client):
    """A width in px against a position in % runs off the end as the map narrows."""
    css = client.get("/static/level2/level2.css").text
    start = css.index("\n.callout {") + 1
    block = css[start:css.index("}", start)]
    width = _property(block, "max-width")
    assert width and width.endswith("%"), f"max-width is {width}, which will not scale"


def test_a_callout_does_not_print_its_component_twice(client):
    """`label` is the component when there is one, so a callout named by
    `label` with the component repeated below said the same thing twice."""
    script = client.get("/static/level2/level2.js").text
    block = script[script.index("name.className = 'callout-name'"):]
    block = block[:block.index("callout-meta")]
    assert "part.label" not in block, "the callout name is still the component"
    assert "part.source_name" in block


# --------------------------------------------------------------- reload ----


def test_serve_does_not_reload_unless_asked():
    """The default has to stay off: `--reload` costs a subprocess and throws
    the index away on every edit, which is wrong for the ordinary run."""
    args = build_parser().parse_args(["serve", "."])
    assert args.reload is False
    assert build_parser().parse_args(["serve", ".", "--reload"]).reload is True


def test_from_env_builds_the_library_the_cli_handed_over(client, monkeypatch):
    """Reload mode respawns the app in a subprocess, so the folders travel
    through the environment. If they did not arrive, the reloaded app would
    come back empty and look like the library had vanished."""
    monkeypatch.setenv(ENV_FOLDERS, json.dumps([str(client.library_dir)]))
    monkeypatch.setenv(ENV_RECURSIVE, "1")
    with TestClient(from_env()) as reloaded:
        assert reloaded.get("/health").json()["plasmids"] == 9


def test_from_env_honours_no_recursive(client, monkeypatch, tmp_path):
    """`--no-recursive` is part of what the subprocess needs to know: dropping
    it would index folders the user deliberately excluded."""
    nested = client.library_dir / "nested"
    nested.mkdir()
    write_genbank(synth.part_plasmid("3", name="pBuried"), nested / "pBuried.gb")

    monkeypatch.setenv(ENV_FOLDERS, json.dumps([str(client.library_dir)]))
    monkeypatch.setenv(ENV_RECURSIVE, "0")
    with TestClient(from_env()) as shallow:
        assert shallow.get("/health").json()["plasmids"] == 9

    monkeypatch.setenv(ENV_RECURSIVE, "1")
    with TestClient(from_env()) as deep:
        assert deep.get("/health").json()["plasmids"] == 10


def test_from_env_says_what_it_is_for_when_nothing_set_it(monkeypatch):
    """Imported by hand it would otherwise fail somewhere inside Library with
    a message about an empty path list."""
    monkeypatch.delenv(ENV_FOLDERS, raising=False)
    with pytest.raises(RuntimeError, match="ggasm serve --reload"):
        from_env()
