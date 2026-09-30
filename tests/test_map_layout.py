"""The Level 2 map's callout geometry, checked headlessly.

The drawing is not the interesting part; the collision pass is. Two 250 bp
parts next to each other put their labels 22 px apart, and the rule that
separates them without letting the leader lines cross is easy to get subtly
wrong in a way no screenshot makes obvious.

So `web/level2/layout.js` is pure arithmetic with no DOM, and this runs it
under node and checks the answers. Skipped where node is unavailable: this
guards a calculation, and a missing toolchain is not a broken one.
"""

from __future__ import annotations

import json
import shutil
import subprocess
import tempfile
from pathlib import Path

import pytest

LAYOUT = (
    Path(__file__).resolve().parent.parent
    / "ggassembler" / "web" / "level2" / "layout.js"
)

#: The cassette the app opens with, from the brief: eight parts round a
#: 6,304 bp circle. Types 1 and 2 land 27 px apart and types 4 and 5 land
#: 22 px apart, so both pairs have to be pushed to the 36 px minimum.
WORKED = [250, 670, 1450, 250, 250, 1300, 750, 1384]
TYPES = ["1", "2", "3", "4", "5", "6", "7", "8"]


def _parts(lengths, features=None):
    parts, at = [], 0
    for length, part_type in zip(lengths, TYPES, strict=False):
        parts.append({
            "part_type": part_type,
            "source_name": f"part{part_type}",
            "start": at,
            "end": at + length,
            "length": length,
            "color": "#6E9FC0",
            "features": (features or {}).get(part_type, []),
        })
        at += length
    return parts


def run_layout(parts, total, start_deg=-90):
    """Call layoutCallouts under node and bring the answer back as JSON."""
    node = shutil.which("node")
    if not node:
        pytest.skip("node is not installed")

    harness = f"""
import {{ layoutCallouts, RING, LABEL }} from {json.dumps(LAYOUT.as_posix())};
const parts = {json.dumps(parts)};
const out = layoutCallouts(parts, {total}, {start_deg});
console.log(JSON.stringify({{
  labelHeight: out.labelHeight,
  gap: out.gap,
  items: out.items.map((i) => ({{
    partType: i.part.part_type,
    theta: i.theta,
    side: i.side,
    naturalY: i.naturalY,
    y: i.y,
    stubStart: i.stubStart,
    stubEnd: i.stubEnd,
    labelX: i.labelX,
    anchor: i.anchor,
    displacement: i.displacement,
    needsElbow: i.needsElbow,
  }})),
  ring: RING,
  label: LABEL,
}}));
"""
    with tempfile.TemporaryDirectory() as tmp:
        script = Path(tmp) / "harness.mjs"
        script.write_text(harness, encoding="utf-8")
        result = subprocess.run(
            [node, str(script)], capture_output=True, text=True, timeout=30
        )
    if result.returncode:
        pytest.fail(f"layout.js failed:\n{result.stderr}")
    return json.loads(result.stdout)


@pytest.fixture(scope="module")
def worked():
    return run_layout(_parts(WORKED), sum(WORKED))


# --------------------------------------------------------------------------- #
# the rules that stop leader lines crossing
# --------------------------------------------------------------------------- #


def test_labels_stay_in_angular_order_on_each_side(worked):
    """Reordering is what makes leaders cross, so the order must be untouched."""
    for side in ("right", "left"):
        column = [i for i in worked["items"] if i["side"] == side]
        column.sort(key=lambda i: i["theta"] if side == "right" else -i["theta"])
        ys = [i["y"] for i in column]
        assert ys == sorted(ys), f"{side} column is out of angular order: {ys}"


def test_no_two_labels_on_a_side_are_closer_than_the_minimum(worked):
    gap = worked["gap"]
    for side in ("right", "left"):
        column = sorted(
            (i for i in worked["items"] if i["side"] == side), key=lambda i: i["y"]
        )
        for before, after in zip(column, column[1:], strict=False):
            room = after["y"] - before["y"]
            assert room >= gap - 0.01, (
                f"{before['partType']} and {after['partType']} are {room:.1f} px apart"
            )


def test_every_label_stays_inside_the_container(worked):
    ring, height = worked["ring"], worked["labelHeight"]
    top = ring["pad"] + height / 2
    bottom = ring["height"] - ring["pad"] - height / 2
    for item in worked["items"]:
        assert top - 0.5 <= item["y"] <= bottom + 0.5, (
            f"{item['partType']} sits at y={item['y']:.1f}, outside [{top}, {bottom}]"
        )


def test_the_worked_cassette_needs_no_correction_at_this_radius(worked):
    """The brief predicts types 1/2 land 27 px apart and 4/5 land 22 px apart.

    Neither follows from the geometry the same brief specifies. With the dot at
    r = outer(150) + stub(22) = 172, the real separations are 73.3 px and
    36.9 px, so nothing in this cassette is crowded and no label moves. The two
    predicted figures are not even consistent with one another: 22 px for 4/5
    implies r = 103, 27 px for 1/2 implies r = 63.

    Pinned here rather than worked around, so that if the artboard turns out to
    use different radii this test says exactly what changed.
    """
    by_type = {i["partType"]: i for i in worked["items"]}
    assert abs(by_type["2"]["naturalY"] - by_type["1"]["naturalY"]) == pytest.approx(73.3, abs=0.2)
    assert abs(by_type["5"]["naturalY"] - by_type["4"]["naturalY"]) == pytest.approx(36.9, abs=0.2)
    assert all(abs(i["displacement"]) < 0.01 for i in worked["items"])


def test_a_genuinely_crowded_pair_is_pushed_to_the_minimum():
    """Four parts where two are tiny and adjacent: the case callouts exist for."""
    lengths = [3000, 60, 60, 3000]
    out = run_layout(_parts(lengths), sum(lengths))
    by_type = {i["partType"]: i for i in out["items"]}
    natural = abs(by_type["3"]["naturalY"] - by_type["2"]["naturalY"])
    final = abs(by_type["3"]["y"] - by_type["2"]["y"])
    assert natural < out["gap"], f"the fixture is not crowded: {natural:.1f} px"
    assert final >= out["gap"] - 0.01, f"still overlapping at {final:.1f} px"


def test_crowding_spreads_both_ways_rather_than_piling_downward():
    """The backward pass is the point: without it a crowded group slides down
    the page and the last label leaves the container."""
    lengths = [2000, 50, 50, 50, 50, 2000]
    out = run_layout(_parts(lengths), sum(lengths))
    moved = [i for i in out["items"] if abs(i["displacement"]) > 0.01]
    assert moved, "nothing was crowded"
    assert any(i["displacement"] < 0 for i in moved), "every label was pushed down"
    assert any(i["displacement"] > 0 for i in moved), "every label was pushed up"


def test_a_label_that_barely_moved_needs_no_elbow(worked):
    for item in worked["items"]:
        if abs(item["displacement"]) <= worked["label"]["elbowThreshold"]:
            assert not item["needsElbow"]
        else:
            assert item["needsElbow"]


# --------------------------------------------------------------------------- #
# sides and geometry
# --------------------------------------------------------------------------- #


def test_the_right_half_anchors_left_and_the_left_half_anchors_right(worked):
    for item in worked["items"]:
        if item["theta"] < 180:
            assert item["side"] == "right" and item["anchor"] == "start"
        else:
            assert item["side"] == "left" and item["anchor"] == "end"


def test_a_stub_runs_outward_along_its_own_radius(worked):
    """The stub has to be radial, or the dot does not sit where the arc points."""
    ring = worked["ring"]
    for item in worked["items"]:
        for point, radius in (
            (item["stubStart"], ring["outer"]),
            (item["stubEnd"], ring["outer"] + ring["stub"]),
        ):
            dx = point["x"] - ring["cx"]
            dy = point["y"] - ring["cy"]
            assert abs((dx * dx + dy * dy) ** 0.5 - radius) < 0.01


def test_twelve_o_clock_is_the_top_of_the_circle():
    """The whole angle convention rests on this, and it is easy to invert."""
    parts = _parts([100, 100, 100, 100], )
    out = run_layout(parts, 400, start_deg=0)
    first = out["items"][0]
    ring = out["ring"]
    # the first part starts at theta 0 and spans 90 degrees, so its middle is
    # at 45: up and to the right of centre
    assert first["theta"] == pytest.approx(45)
    assert first["stubStart"]["x"] > ring["cx"]
    assert first["stubStart"]["y"] < ring["cy"]


def test_the_ring_starts_where_the_css_gradient_starts():
    """level2.js draws the band with conic-gradient(from -90deg); the callouts
    have to agree or every label points at the wrong arc."""
    out = run_layout(_parts([100, 100, 100, 100]), 400, start_deg=-90)
    ring = out["ring"]
    first = out["items"][0]
    # first part spans -90 to 0 degrees, middle at -45, i.e. up and to the left
    assert first["stubStart"]["x"] < ring["cx"]
    assert first["stubStart"]["y"] < ring["cy"]


# --------------------------------------------------------------------------- #
# crowding harder than the real case
# --------------------------------------------------------------------------- #


def test_many_tiny_parts_still_fit_inside_the_container():
    """Sixteen equal parts is worse than anything the kit produces; the layout
    should tighten rather than run off the page."""
    lengths = [100] * 16
    global TYPES
    kept = TYPES
    TYPES = [str(i) for i in range(16)]
    try:
        out = run_layout(_parts(lengths), sum(lengths))
    finally:
        TYPES = kept

    ring, height = out["ring"], out["labelHeight"]
    top = ring["pad"] + height / 2
    bottom = ring["height"] - ring["pad"] - height / 2
    for item in out["items"]:
        assert top - 0.5 <= item["y"] <= bottom + 0.5
    assert height <= out["label"]["height"]


def test_a_single_part_is_not_moved_at_all():
    out = run_layout(_parts([1000]), 1000)
    assert out["items"][0]["displacement"] == pytest.approx(0)
    assert not out["items"][0]["needsElbow"]
