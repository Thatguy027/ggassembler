"""The alignment viewer's column arithmetic, checked headlessly.

A screen column is not a reference base: an insertion in any clone opens
columns that belong to no reference position, and everything after them shifts.
Getting it wrong breaks nothing visibly - the bases render, the ruler counts -
it just points every label a few bases off, which is the one error a viewer
whose whole job is *where* must not make.

So `web/sequencing/columns.js` is arithmetic with no DOM, and this runs it
under node. Skipped where node is unavailable, like the map's layout test.
"""

from __future__ import annotations

import json
import shutil
import subprocess
import tempfile
from pathlib import Path

import pytest

COLUMNS = (
    Path(__file__).resolve().parent.parent
    / "ggassembler" / "web" / "sequencing" / "columns.js"
)


def run_js(body: str):
    """Run `body` against columns.js under node and return what it prints."""
    node = shutil.which("node")
    if not node:
        pytest.skip("node is not installed")

    with tempfile.TemporaryDirectory() as tmp:
        module = Path(tmp) / "columns.mjs"
        module.write_bytes(COLUMNS.read_bytes())
        script = Path(tmp) / "run.mjs"
        script.write_text(
            "import { index, classOf, runs, span, placed, ticks } from './columns.mjs';\n" + body,
            encoding="utf-8",
        )
        result = subprocess.run(
            [node, str(script)], capture_output=True, text=True, timeout=30
        )
    assert result.returncode == 0, result.stderr
    return json.loads(result.stdout)


# --------------------------------------------------------------------------- #
# columns and positions
# --------------------------------------------------------------------------- #


def test_with_no_insertions_a_column_is_its_own_base():
    out = run_js("""
      const { colOf, refOf, columns } = index([], 100);
      console.log(JSON.stringify({
        columns,
        colOf: [colOf[0], colOf[50], colOf[99]],
        refOf: [refOf[0], refOf[50], refOf[99]],
      }));
    """)
    assert out["columns"] == 100
    assert out["colOf"] == [0, 50, 99]
    assert out["refOf"] == [0, 50, 99]


def test_every_base_after_an_insertion_shifts_by_its_width():
    """The bug this exists for: labels that are right until the first indel."""
    out = run_js("""
      const { colOf, refOf, columns } = index([[40, 6]], 100);
      console.log(JSON.stringify({
        columns,
        before: colOf[39],
        insertionColumns: [refOf[40], refOf[45]],
        after: colOf[40],
        last: colOf[99],
        roundTrip: refOf[colOf[77]],
      }));
    """)
    assert out["columns"] == 106
    assert out["before"] == 39
    # the six opened columns belong to no reference base at all
    assert out["insertionColumns"] == [-1, -1]
    assert out["after"] == 46
    assert out["last"] == 105
    assert out["roundTrip"] == 77


def test_several_insertions_accumulate():
    out = run_js("""
      const { colOf, columns } = index([[10, 2], [20, 3], [30, 1]], 50);
      console.log(JSON.stringify({
        columns, at9: colOf[9], at10: colOf[10], at20: colOf[20], at30: colOf[30],
      }));
    """)
    assert out["columns"] == 56
    assert out["at9"] == 9
    assert out["at10"] == 12     # two columns opened before it
    assert out["at20"] == 25     # two more, plus three
    assert out["at30"] == 36     # and one more


def test_an_insertion_at_the_very_start_shifts_everything():
    out = run_js("""
      const { colOf, refOf, columns } = index([[0, 4]], 10);
      console.log(JSON.stringify({ columns, first: colOf[0], refAt0: refOf[0] }));
    """)
    assert out["columns"] == 14
    assert out["first"] == 4
    assert out["refAt0"] == -1


# --------------------------------------------------------------------------- #
# what a cell is
# --------------------------------------------------------------------------- #


def test_each_kind_of_cell_is_told_apart():
    out = run_js("""
      console.log(JSON.stringify([
        classOf('A', 'A'),   // match
        classOf('A', 'G'),   // mismatch
        classOf('A', '-'),   // the clone deleted this base
        classOf('-', 'T'),   // the clone inserted a base here
        classOf('-', '-'),   // padding: another clone's insertion, not this one's
        classOf('A', '.'),   // never covered by a read
        classOf('A', '?'),   // too divergent to align - nobody compared these
      ]));
    """)
    assert out == ["m", "x", "d", "i", "p", "u", "n"]


def test_an_unalignable_stretch_is_not_painted_as_a_mismatch():
    """Reading `?` as a mismatch paints eight thousand red cells across a row,
    which says the clone has eight thousand substitutions. What happened is
    that nobody compared those bases at all."""
    out = run_js("""
      const ref   = 'ACGTACGTAC';
      const clone = 'ACG??????C';
      console.log(JSON.stringify(runs(ref, clone, 0, 10).map(([kind, at]) => [kind, at])));
    """)
    assert out == [["m", 0], ["n", 3], ["m", 9]]


def test_cells_come_out_as_runs_not_one_element_each():
    out = run_js("""
      const ref   = 'ACGTACGTACGT';
      const clone = 'ACGTAGGTACGT';
      console.log(JSON.stringify(runs(ref, clone, 0, 12)));
    """)
    assert out == [["m", 0, "ACGTA"], ["x", 5, "G"], ["m", 6, "GTACGT"]]


def test_a_run_is_taken_only_over_the_window_asked_for():
    out = run_js("""
      const ref = 'AAAAAAAAAA';
      console.log(JSON.stringify(runs(ref, 'AAAAAAAAAA', 3, 6)));
    """)
    assert out == [["m", 3, "AAA"]]


# --------------------------------------------------------------------------- #
# the mark on the column you jumped to
# --------------------------------------------------------------------------- #


def test_a_substitution_is_marked_over_its_own_bases():
    out = run_js("""
      const { colOf } = index([], 1000);
      const events = [{ start: 500, end: 506, kind: 'substitution' }];
      console.log(JSON.stringify(span(500, events, colOf, [], 1000)));
    """)
    assert out == [500, 506]


def test_a_deletion_is_marked_over_the_bases_it_removed():
    out = run_js("""
      const { colOf } = index([], 1000);
      const events = [{ start: 200, end: 212, kind: 'deletion' }];
      console.log(JSON.stringify(span(200, events, colOf, [], 1000)));
    """)
    assert out == [200, 212]


def test_an_insertion_is_marked_backwards_over_the_columns_it_opened():
    """The one that goes the other way, and the reason this is a function.

    An insertion's bases live in the columns opened *before* the reference base
    it precedes. A mark running forward from `colOf[start]` would sit on the
    reference sequence after it and point at entirely the wrong bases.
    """
    out = run_js("""
      const { colOf } = index([[300, 8]], 1000);
      const events = [{ start: 300, end: 300, kind: 'insertion' }];
      console.log(JSON.stringify(span(300, events, colOf, [[300, 8]], 1000)));
    """)
    # base 300 sits at column 308; its eight inserted columns are 300-307
    assert out == [300, 309]


def test_a_substitution_past_an_insertion_is_marked_in_shifted_columns():
    out = run_js("""
      const { colOf } = index([[100, 5]], 1000);
      const events = [{ start: 400, end: 403, kind: 'substitution' }];
      console.log(JSON.stringify(span(400, events, colOf, [[100, 5]], 1000)));
    """)
    assert out == [405, 408]


def test_a_base_with_no_difference_gets_a_one_column_mark():
    out = run_js("""
      const { colOf } = index([], 1000);
      console.log(JSON.stringify(span(700, [], colOf, [], 1000)));
    """)
    assert out == [700, 701]


# --------------------------------------------------------------------------- #
# where a painted run is drawn
# --------------------------------------------------------------------------- #


def test_a_run_is_placed_at_its_own_column_not_its_offset_in_the_window():
    """The bug this exists for, and it reads exactly like truncated data.

    Only the columns on screen are in the DOM, but they sit inside a row as
    wide as the whole alignment. Positioning a run at its offset *within the
    painted window* is right only while the window starts at zero; every later
    repaint slides left by `from` columns and lands off screen, so the sequence
    appears to stop partway across and never comes back however far you scroll.
    """
    out = run_js("""
      const ref = 'ACGT'.repeat(500);
      const first = placed(ref, ref, 0, 100, 8);
      const later = placed(ref, ref, 900, 1000, 8);
      console.log(JSON.stringify({
        firstLeft: first[0].left,
        firstAt: first[0].at,
        laterLeft: later[0].left,
        laterAt: later[0].at,
      }));
    """)
    assert out["firstAt"] == 0 and out["firstLeft"] == 0
    # column 900 is drawn at 900 * 8, not at 0 because the window starts there
    assert out["laterAt"] == 900
    assert out["laterLeft"] == 7200


def test_the_same_column_lands_in_the_same_place_whatever_window_painted_it():
    """Scrolling must not move the sequence under the ruler."""
    out = run_js("""
      const ref   = 'ACGTACGTAC'.repeat(100);
      const clone = ref.slice(0, 504) + 'T' + ref.slice(505);
      const wide   = placed(ref, clone, 400, 700, 8).find((r) => r.at === 504);
      const narrow = placed(ref, clone, 500, 520, 8).find((r) => r.at === 504);
      console.log(JSON.stringify([wide.left, narrow.left, wide.kind, narrow.kind]));
    """)
    wide, narrow, wide_kind, narrow_kind = out
    assert wide == narrow == 504 * 8
    assert wide_kind == narrow_kind == "x"


def test_ruler_labels_are_placed_absolutely_too():
    out = run_js("""
      const { refOf } = index([], 1000);
      const early = ticks(refOf, 0, 60, 10, 8);
      const late  = ticks(refOf, 500, 560, 10, 8);
      console.log(JSON.stringify({
        early: early.map((t) => [t.label, t.left]),
        late: late.slice(0, 2).map((t) => [t.label, t.left]),
      }));
    """)
    # labels are 1-based; base 10 sits at column 9
    assert out["early"][0] == [10, 72]
    assert out["late"][0] == [510, 509 * 8]


def test_a_ruler_skips_the_columns_that_are_nobody_base():
    """An insertion column has no reference number to print."""
    out = run_js("""
      const { refOf } = index([[20, 6]], 1000);
      console.log(JSON.stringify(ticks(refOf, 0, 60, 10, 8).map((t) => [t.label, t.left / 8])));
    """)
    # The insertion opens six columns *before* reference index 20, so the base
    # labelled 20 - index 19 - is still ahead of it and does not move. Base 30
    # is index 29, past it, and shifts by six. Labels are 1-based; columns and
    # indices are not, which is exactly the confusion this module exists to
    # keep out of the rest of the screen.
    assert out == [[10, 9], [20, 19], [30, 35], [40, 45], [50, 55]]
