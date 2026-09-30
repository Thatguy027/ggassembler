# Claude Code prompt — GG Assembler, map treatment + next features

> Paste below the line, from the repo root.

---

Two pieces of work on GG Assembler: redesign the Level 2 plasmid map, then a set of features. Milestones are ordered; stop after each, run `uv run pytest`, and report.

**Before you start**, read the current state rather than trusting this document on it — `ggassembler/web/level2/level2.js`, `ggassembler/api/routes_level2.py`, `ggassembler/core/library.py`, `ggassembler/core/assembly.py`. Some of what follows may already exist. Keep the layering rules in `tests/test_layering.py` intact: `levels/*` never import each other, `core/*` imports neither `levels/` nor `api/`, each front-end level directory stays self-contained.

The design reference is the new **Map.dc.html** artboard on the design canvas, exported to `docs/design/`. Match it.

## M1 — API: features per placed part

`_result_payload` in `routes_level2.py` returns each `PlacedPart` with `start`, `end`, `color`, overhangs. Add a `features` list to each part: the annotated features of the source record that fall inside that part's span, remapped to product coordinates.

`core/library.py` already has `FeatureInfo` and `features_within(features, ...)` — use them, do not write a second implementation. Each feature entry needs `label`, `kind`, `start`, `end`, `strand`. Cap at the ~8 largest per part so the payload stays small; drop features under 60 bp.

Add a test asserting that a cassette built from the fixtures returns a CDS feature inside its Type 3 part with coordinates that fall within that part's `start`/`end`.

## M2 — The map

Replace the ring + legend-below with the callout treatment. Keep the existing clickable `.arc-hit` SVG wedge layer and its aria labels — that accessibility work stays.

### Geometry

Container 680 × 500, ring centre at (360, 270), outer radius 150.

Four concentric layers, outside in:

1. **Part band**, r 150 → 116. The existing conic gradient, unchanged, with its 1.6° hairline gaps.
2. **Gap**, r 116 → 110.
3. **Feature band**, r 110 → 96. A second conic gradient built from the M1 feature data, `transparent` wherever no feature sits. Fill each feature with its part's colour lightened ~25% toward white.
4. **Centre disc**, r 96. Length, `circular · N parts`, junction pill. Unchanged.

### Callouts

Per part, at its mid-arc angle θ (measured clockwise from 12 o'clock, consistent with `RING_START_DEG`):

- **Stub**: a 22 px line from r=150 outward along the radius. As a rotated element, `transform-origin: 0 50%; rotate(θ − 90deg)`.
- **Dot**: 4 px at the stub end (r=172).
- **Label**: 30 px tall, anchored at the stub end — left-anchored for θ < 180°, right-anchored (`flex-direction: row-reverse`, right-aligned text) for θ ≥ 180°. Contents: the type badge, the part name at 12 px/600, and a mono 10 px line `<length> bp · <component>` with a `→` when the part carries a stranded coding feature.

**Collision handling** — this is the part that will bite you. Compute each label's natural y from its stub end, then within each side:

1. Sort by angle, never by y, and never reorder. Reordering is what makes leaders cross.
2. Walk the list and push any label less than 36 px from its predecessor down to exactly 36 px.
3. Walk backwards from the bottom doing the same, so the correction spreads both ways instead of piling everything downward.
4. If the column still overflows the container, shrink the label height to 26 px and the gap to 30 px before you consider anything else.

Displacements under ~10 px need no elbow — the stub end still lands inside the label box. Above that, draw a vertical connector at the stub-end x from the stub end to the label's centre y.

Worked example, for a check: the cassette the app opens with is 6,304 bp with parts 250 / 670 / 1,450 / 250 / 250 / 1,300 / 750 / 1,384 bp. Types 1 and 2 land 27 px apart, and types 4 and 5 land 22 px apart — both pairs must be pushed to 36. Types 4 and 5 are only 14° of arc each, which is the whole reason for callouts: no label fits inside those arcs at any ring size.

### Chevrons

For each part whose largest feature is stranded, a small triangle at r=133, mid-arc, `rotate(θ)` for the forward strand and `rotate(θ + 180deg)` for the reverse.

### Linear view

Below the ring, a 40 px proportional track opened at the first part, `flex-grow` set to each part's bp. Parts whose largest feature is a stranded CDS get an arrow end via `clip-path: polygon(0 0, calc(100% - 14px) 0, 100% 50%, calc(100% - 14px) 100%, 0 100%)`. Name inside the block when it is wide enough, type number underneath.

### Toggles

The header's `Features` button toggles the feature band and the chevrons together; a second toggle shows/hides the linear view. Persist both in `localStorage` beside the existing recents key, wrapped in try/catch.

### Deliberately not doing

No bp ruler — it competes for the same annulus as the callouts and part lengths are already in the labels. No overhang labels inside the ring — eight 4-nt mono strings at that radius is clutter; they stay in the strip below, and hovering a ring boundary highlights the matching chip.

### Tests

The geometry is worth testing headlessly. Put the angle → stub/label maths in a pure module (`web/level2/layout.js` or similar, no DOM) and test it in Node via a small harness, or port the collision pass to Python under `core/` if that is cleaner given the layering rules. Assert: labels stay in angular order; no two labels on a side are closer than the minimum; every label stays inside the container.

## M3 — Decompose an existing cassette

The highest-value feature. 220 plasmids in the current library index as cassettes.

`POST /api/level2/decompose {name}` — take an assembled cassette from the library, digest it, match each released fragment back to a library part by overhang pair *and* sequence identity, and return a slot payload the Level 2 screen can load directly. Fragments with no library match come back as `unmatched` with their span and overhangs, so the screen shows "this position holds 1,450 bp that isn't a part I know" rather than silently dropping it.

On the screen: a "Load from cassette…" control in the header. This turns Level 2 from a build-from-scratch form into "clone this and swap the promoter", which is the actual workflow.

## M4 — Level 3 readiness

Some Type 1/5 parts in the library have had their internal BsmBI site removed — `L13_ConLS_BSMB1del` is the explicit case, and its own annotation says "Former Bsmb1". A cassette built from one assembles fine at Level 2 and then **cannot be used at Level 3**, and nothing surfaces that until the multigene step fails.

At index time, compute `level3_ready: bool` per part: does a Type 1 or Type 5 part still carry an internal BsmBI site yielding a resolvable connector overhang? Surface it in the Level 2 picker rows ("usable here · blocks multigene"), on the assembled-cassette summary, and as a Library filter.

## M5 — Part-swap matrix

`POST /api/level2/sweep {base_design, slot, candidates[]}` returns one assembly result per candidate. UI: pick a slot, multi-select candidates, get a results table (name, length, ok/issues) and a "write all" that exports every `.gbk` plus one combined picklist. The library has 47 Type 2 promoters indexed — a promoter titration against a fixed CDS is a standard experiment and this is two functions from working.

## M6 — Protocol, for real

`core/protocol.py` is currently a docstring, and `level2.js` only toggles `disabled` on `protocol-btn` — there is no click handler. Either implement or delete the button; a button that looks live and does nothing is worse than none.

To implement, you need a concentration per plasmid, which does not exist in the data model. Add `conc_ng_ul` to the library entry, persisted in `.ggasm/overrides.json`, editable from the Library row, defaulting to 50. Then emit, as JSON and copyable text: enzyme and buffer, equimolar 20 fmol per part with volumes from each plasmid's length and concentration, cycling (`30 × (37 °C 5 min / 16 °C 5 min)`, 60 °C 10 min, 80 °C 10 min), and the antibiotic plus colony-colour selection implied by the chosen dropout. When a reversed-site dropout (`234r`, `3r`, `8ar` — 34 such plasmids are indexed) is in the reaction, drop the final digest and heat-inactivation and add the note about elevated misassembly and chloramphenicol counter-screening.

## M7 — Smaller

- **Level 3 seeds empty** and shows `no_backbone` / `empty_slot` in red on arrival. Level 2 opens on a worked assembly; do the same here from one of the 31 indexed multigene vectors.
- **Level 1 has no library entry point.** Its template dropdown leads with "Paste a sequence instead…" and amplicon coordinates must be typed by hand. Let the user pick a library plasmid, then one of its annotated features, and fill start/end from it.
- **Type 1 panel is mislabelled** "5′ assembly connector" — the type is overhang-defined and the list includes `LexA_4 [LexA site]`. Lead with `Type 1 · CCCT→AACG`, demote the content description to a subtitle.
- **Library triage.** 83 plasmids are unrecognised and 26 conflict. Group the unrecognised by the `reason` string the detector already writes, and make those groups filter chips. For the 121 merged duplicate files, let the user pick the canonical name so pickers stop reading `pYTK003 / A3_ConL1`.
- **Build log.** Every export appends to `.ggasm/builds.jsonl`: timestamp, design name, each part with its file hash, predicted length, issue codes.
- **`ggasm check <folder>`** — run detection, exit nonzero on new conflicts or unrecognised files relative to a committed baseline. Usable as a pre-commit hook on the plasmid repo.
- **Verify the Library table is virtualized** before adding columns — 622 rows rendered eagerly behind a filter box will get slow.

## Constraints throughout

Keep it a no-build-step front end: vanilla ES modules, plain CSS, tokens from `web/tokens.css`. No new Python dependencies without saying why. Every new endpoint gets a test in the matching `tests/test_*.py`. Do not touch another level's files while working on one — and show me the changed-file list per milestone so that stays honest.
