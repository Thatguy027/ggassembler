"""The Level 1 option explainers: one per option, and none left orphaned.

The failure this guards against is quiet. Adding a seventh convention checkbox
takes one line of markup and looks finished; the option then sits there with no
ⓘ beside it, or with an ⓘ that opens an empty panel, and nothing anywhere says
so. Both directions are checked, because an entry for an option that no longer
exists is dead copy that will be read as current.
"""

from __future__ import annotations

import json
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

import pytest

WEB = Path(__file__).resolve().parent.parent / "ggassembler" / "web" / "level1"
HTML = (WEB / "index.html").read_text(encoding="utf-8")
EXPLAINERS = (WEB / "explainers.js").read_text(encoding="utf-8")

#: The options that get an explanation: the sequence conventions and
#: domestication. Read out of the markup rather than listed here, so the test
#: cannot drift from the page.
OPTION_IDS = sorted(set(re.findall(r'<input type="checkbox" id="(c-[a-z]+)"', HTML)))


def run_js(body: str):
    node = shutil.which("node")
    if not node:
        pytest.skip("node is not installed")
    with tempfile.TemporaryDirectory() as tmp:
        (Path(tmp) / "explainers.mjs").write_bytes((WEB / "explainers.js").read_bytes())
        script = Path(tmp) / "run.mjs"
        script.write_text(
            "import { EXPLAINERS, MARKS, segments } from './explainers.mjs';\n" + body,
            encoding="utf-8",
        )
        result = subprocess.run([node, str(script)], capture_output=True, text=True, timeout=30)
    assert result.returncode == 0, result.stderr
    return json.loads(result.stdout)


def test_the_options_were_found_at_all():
    """Guard the guard: a regex that matches nothing would pass everything."""
    assert len(OPTION_IDS) >= 6
    assert "c-glyser" in OPTION_IDS and "c-domesticate" in OPTION_IDS


def test_every_option_has_an_info_button_pointing_at_the_shared_panel():
    missing = [
        option for option in OPTION_IDS
        if f'data-explains="{option}"' not in HTML
    ]
    assert not missing, f"options with no ⓘ: {missing}"


def test_every_info_button_is_a_real_button_with_a_label():
    """A div with a click handler is not a control, and is not reachable."""
    for option in OPTION_IDS:
        button = re.search(
            rf'<button[^>]*data-explains="{option}"[^>]*>', HTML
        ) or re.search(rf'<button[^>]*data-explains="{option}"(?:.|\n)*?>', HTML)
        assert button, f"{option}: no <button> carries data-explains"
        markup = button.group(0)
        assert 'type="button"' in markup, f"{option}: would submit, not explain"
        assert 'aria-controls="explainer"' in markup, f"{option}: points at no panel"
        assert 'aria-expanded=' in markup, f"{option}: never says whether it is open"
        assert 'aria-label="About:' in markup, f"{option}: unlabelled for a screen reader"


def test_the_info_buttons_sit_outside_the_labels_they_explain():
    """Inside a <label>, clicking the button also toggles the checkbox."""
    for row in re.findall(r'<label[^>]*>(?:.|\n)*?</label>', HTML):
        assert "data-explains" not in row, "an ⓘ is nested inside a <label>"


def test_the_panel_is_one_region_that_announces_itself():
    panel = re.search(r'<section[^>]*id="explainer"(?:.|\n)*?>', HTML)
    assert panel, "no shared explainer region"
    markup = panel.group(0)
    assert 'aria-live="polite"' in markup
    assert 'tabindex="-1"' in markup, "focus cannot be moved to it"
    assert HTML.count('id="explainer"') == 1, "more than one panel"


def test_the_panel_comes_after_the_predicted_part_plasmid_card():
    """It explains what you are looking at, so it goes below it, not above."""
    assert HTML.index("Predicted part plasmid") < HTML.index('id="explainer"')


def test_every_option_has_an_entry_and_no_entry_is_orphaned():
    keys = run_js("console.log(JSON.stringify(Object.keys(EXPLAINERS)));")
    assert sorted(keys) == OPTION_IDS


def test_every_entry_says_something():
    """An entry with a title and no prose is a panel that opens onto nothing.

    Counted across every field that holds prose, not just `body`: some options
    carry most of their substance in the worked example and the practice notes,
    and a threshold that only looked at `body` would be measuring layout rather
    than whether anything was explained.
    """
    entries = run_js("""
      const prose = (v) => [
        ...(v.body || []),
        v.example ? (v.example.caption || '') + ' ' + (v.example.after || '') : '',
        v.cost || '',
        ...(v.when || []),
      ].join(' ');
      console.log(JSON.stringify(Object.fromEntries(
        Object.entries(EXPLAINERS).map(([k, v]) => [k, {
          title: Boolean(v.title),
          words: prose(v).split(/\\s+/).filter(Boolean).length,
        }])
      )));
    """)
    for key, entry in entries.items():
        assert entry["title"], f"{key}: no title"
        assert entry["words"] >= 40, f"{key}: {entry['words']} words is not an explanation"


def test_the_marked_up_examples_survive_being_parsed():
    """Every mark used must be one the key can name, or it renders uncoloured."""
    out = run_js("""
      const used = new Set();
      let blank = [];
      for (const [key, entry] of Object.entries(EXPLAINERS)) {
        if (!entry.example) continue;
        const runs = segments(entry.example.text);
        if (runs.reduce((n, r) => n + r.text.length, 0) !== entry.example.text.replace(
              /\\[\\[([a-z]+)\\|([^\\]]*)\\]\\]/g, '$2').length) blank.push(key);
        for (const run of runs) if (run.mark) used.add(run.mark);
      }
      console.log(JSON.stringify({ used: [...used], marks: Object.keys(MARKS), lossy: blank }));
    """)
    assert not out["lossy"], f"an example lost text when parsed: {out['lossy']}"
    assert out["used"], "no example marks any bases"
    assert set(out["used"]) <= set(out["marks"]), "an example uses a mark with no key entry"


def test_the_marker_syntax_keeps_the_text_around_it():
    out = run_js("""
      console.log(JSON.stringify(segments('...NNN   [[added|GG]] [[overhang|ATCC]]   x')));
    """)
    assert out == [
        {"text": "...NNN   "},
        {"text": "GG", "mark": "added"},
        {"text": " "},
        {"text": "ATCC", "mark": "overhang"},
        {"text": "   x"},
    ]


def test_text_with_no_markers_comes_back_whole():
    assert run_js("console.log(JSON.stringify(segments('ATG CTA')));") == [{"text": "ATG CTA"}]


def test_the_gly_ser_entry_carries_the_interlock_and_the_cost():
    """Two things a reader has to be told, and the reason this copy exists."""
    entry = run_js("console.log(JSON.stringify(EXPLAINERS['c-glyser']));")
    assert "relatively innocuous" in entry["cost"]
    assert any("Strip a trailing stop codon" in line for line in entry["when"])
    assert "BamHI" in entry["example"]["after"]
