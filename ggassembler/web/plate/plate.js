/* Plate - the standing source plate, and the sweeps built from it.
 *
 * The grid is the object, not a table of the same data. Ninety-six wells at a
 * glance is the only view that answers "which one will run out first", and a
 * list of ninety-six rows never does, however well sorted.
 *
 * So the grid carries state and the list carries detail, and they are two
 * views of one thing rather than two representations to keep in step: both are
 * rendered from the same payload on every change.
 */

const el = (id) => document.getElementById(id);

async function api(path, options) {
  const response = await fetch(path, options);
  const body = await response.json().catch(() => ({}));
  if (!response.ok) throw new Error(body.detail || `${response.status} ${response.statusText}`);
  return body;
}

const post = (path, payload) => api(path, {
  method: 'POST',
  headers: { 'Content-Type': 'application/json' },
  body: JSON.stringify(payload),
});

const state = {
  geometry: null,
  plate: null,
  selected: null,
  candidates: [],
  slots: [],
  base: {},
  row: { position: '', candidates: [] },
  column: { position: '', candidates: [] },
  sweep: null,
  timer: null,
};

const LIMIT = { row: 8, column: 12 };

/* ------------------------------------------------------------------ load */

async function load() {
  try {
    const [geometry, plate, candidates, slots] = await Promise.all([
      api('/api/plate/wells'),
      api('/api/plate/source'),
      api('/api/plate/candidates'),
      api('/api/plate/slots'),
    ]);
    state.geometry = geometry;
    state.plate = plate;
    state.candidates = candidates;
    state.slots = slots;
    // the base cassette starts from the first part that fits each position, so
    // the screen opens on something that assembles rather than on eight empty
    // dropdowns - the same choice the Cassette screen makes
    for (const slot of slots) {
      if (slot.options.length) state.base[slot.key] = slot.options[0].name;
    }
    state.row.position = pick(slots, '2');
    state.column.position = pick(slots, '3');
  } catch (error) {
    setStatus(`could not read the plate: ${error.message}`, true);
    return;
  }
  el('plasmid-list').replaceChildren(...state.candidates.map((c) => {
    const option = document.createElement('option');
    option.value = c.name;
    option.label = c.display || c.name;
    return option;
  }));
  renderFactors();
  renderAll();
  runSweep();
}

/** The first position with options, preferring the one asked for. */
function pick(slots, wanted) {
  const found = slots.find((s) => s.key === wanted && s.options.length);
  return found ? found.key : (slots.find((s) => s.options.length) || { key: '' }).key;
}

function renderAll() {
  renderAxes();
  renderGrid();
  renderList();
  renderIssues();
  renderSettings();
  const filled = state.plate.wells.length;
  el('source-note').textContent = filled
    ? `${state.plate.id} · ${filled} of 96 wells filled`
    : `${state.plate.id} · empty — click a well to put something in it`;
}

/* ------------------------------------------------------------------ grid */

function renderAxes() {
  el('source-cols').replaceChildren(...state.geometry.columns.map((n) => {
    const cell = document.createElement('span');
    cell.textContent = String(n);
    return cell;
  }));
  el('source-rows').replaceChildren(...state.geometry.rows.map((row) => {
    const cell = document.createElement('span');
    cell.textContent = row;
    return cell;
  }));
}

function byWell() {
  return new Map(state.plate.wells.map((w) => [w.well, w]));
}

function renderGrid() {
  const filled = byWell();
  const dead = state.plate.dead_volume_ul;
  const grid = el('source-grid');

  // the fill height is relative to the fullest well, not to the labware's
  // capacity: what matters here is which well runs out first, and against a
  // 200 uL maximum every well in a real plate looks equally empty
  const fullest = Math.max(10, ...state.plate.wells.map((w) => w.volume_ul));

  grid.replaceChildren(...state.geometry.row_major.map((name) => {
    const entry = filled.get(name);
    const cell = document.createElement('button');
    cell.type = 'button';
    cell.className = 'well';
    cell.dataset.well = name;
    cell.setAttribute('role', 'gridcell');

    if (!entry) {
      cell.classList.add('is-empty');
      cell.setAttribute('aria-label', `${name}, empty`);
    } else {
      const low = entry.volume_ul <= dead;
      if (!entry.known) cell.classList.add('is-bad');
      else if (low) cell.classList.add('is-low');
      cell.dataset.part = entry.part_type || '';

      const fill = document.createElement('i');
      fill.className = 'well-fill';
      fill.style.height = `${Math.min(100, (entry.volume_ul / fullest) * 100)}%`;
      cell.append(fill);

      cell.setAttribute('aria-label',
        `${name}, ${entry.plasmid}, ${entry.volume_ul} microlitres`
        + (entry.conc_ng_ul ? ` at ${entry.conc_ng_ul} nanograms per microlitre` : ', no concentration')
        + (entry.known ? '' : ', not in the library')
        + (low ? ', at or below dead volume' : ''));
      cell.title = `${name} · ${entry.plasmid}\n`
        + `${entry.volume_ul} µL${entry.conc_ng_ul ? ` at ${entry.conc_ng_ul} ng/µL` : ' · no concentration'}\n`
        + `${entry.usable_ul.toFixed(1)} µL reachable above the ${dead} µL dead volume`;
    }
    if (state.selected === name) cell.classList.add('is-selected');
    cell.addEventListener('click', () => select(name));
    cell.addEventListener('keydown', (event) => move(event, name));
    return cell;
  }));
}

/* Arrow keys walk the plate. A 96-cell grid where tab is the only way across
 * is 96 tab stops to reach H12, which is not navigation. */
function move(event, from) {
  const keys = { ArrowLeft: [0, -1], ArrowRight: [0, 1], ArrowUp: [-1, 0], ArrowDown: [1, 0] };
  const step = keys[event.key];
  if (!step) return;
  event.preventDefault();
  const rows = state.geometry.rows;
  const row = rows.indexOf(from[0]);
  const column = Number(from.slice(1));
  const next = `${rows[Math.min(rows.length - 1, Math.max(0, row + step[0]))]}`
    + `${Math.min(12, Math.max(1, column + step[1]))}`;
  const target = el('source-grid').querySelector(`[data-well="${next}"]`);
  if (target) target.focus();
}

/* ---------------------------------------------------------------- editor */

function select(name) {
  state.selected = name;
  const entry = byWell().get(name);

  el('editor').hidden = false;
  el('editor-well').textContent = name;
  el('editor-what').textContent = entry
    ? (entry.known ? `currently ${entry.plasmid}` : `${entry.plasmid} — not in the library`)
    : 'empty';
  el('editor-clear').hidden = !entry;
  el('e-plasmid').value = entry ? entry.plasmid : '';
  el('e-conc').value = entry && entry.conc_ng_ul ? entry.conc_ng_ul : '';
  el('e-volume').value = entry ? entry.volume_ul : '';
  el('editor-note').textContent = ' ';

  renderGrid();
  el('e-plasmid').focus();
}

async function saveWell() {
  if (!state.selected) return;
  const plasmid = el('e-plasmid').value.trim();
  if (!plasmid) {
    el('editor-note').textContent = 'name a plasmid, or use "Empty this well"';
    return;
  }
  const concentration = el('e-conc').value.trim();
  const volume = el('e-volume').value.trim();
  try {
    state.plate = await post('/api/plate/source/well', {
      well: state.selected,
      plasmid,
      // an empty box is "leave it as it was", not "set it to zero" - the two
      // are different and only one of them is ever what was meant
      conc_ng_ul: concentration === '' ? null : Number(concentration),
      volume_ul: volume === '' ? null : Number(volume),
    });
  } catch (error) {
    el('editor-note').textContent = error.message;
    return;
  }
  renderAll();
  select(state.selected);
  el('editor-note').textContent = `${state.selected} saved`;
}

async function clearWell() {
  if (!state.selected) return;
  try {
    state.plate = await post('/api/plate/source/well', { well: state.selected, plasmid: '' });
  } catch (error) {
    el('editor-note').textContent = error.message;
    return;
  }
  renderAll();
  select(state.selected);
}

/* ------------------------------------------------------------- the list */

function renderList() {
  const dead = state.plate.dead_volume_ul;
  el('well-list').replaceChildren(...state.plate.wells.map((entry) => {
    const row = document.createElement('li');
    row.className = 'well-row';
    if (!entry.known) row.classList.add('is-bad');
    else if (entry.volume_ul <= dead) row.classList.add('is-low');

    const where = document.createElement('b');
    where.textContent = entry.well;
    const what = document.createElement('span');
    what.className = 'w-name';
    what.textContent = entry.plasmid;
    const numbers = document.createElement('span');
    numbers.className = 'w-num';
    numbers.textContent = `${entry.volume_ul} µL`
      + (entry.conc_ng_ul ? ` · ${entry.conc_ng_ul} ng/µL` : ' · no conc.');
    const why = document.createElement('span');
    why.className = 'w-why';
    why.textContent = !entry.known
      ? 'not in the library'
      : entry.volume_ul <= dead ? `at or below the ${dead} µL dead volume` : '';

    row.append(where, what, numbers, why);
    row.addEventListener('click', () => {
      select(entry.well);
      el('editor').scrollIntoView({ behavior: 'smooth', block: 'nearest' });
    });
    return row;
  }));
}

function renderIssues() {
  el('plate-issues').replaceChildren(...state.plate.issues.map((issue) => {
    const row = document.createElement('li');
    row.className = `issue is-${issue.level}`;
    row.textContent = issue.message;
    return row;
  }));
}

/* ---------------------------------------------------------- the plate itself */

function renderSettings() {
  el('p-id').value = state.plate.id;
  el('p-labware').value = state.plate.labware;
  el('p-dead').value = state.plate.dead_volume_ul;
}

async function saveSettings() {
  try {
    state.plate = await post('/api/plate/source/layout', {
      id: el('p-id').value.trim(),
      labware: el('p-labware').value.trim(),
      dead_volume_ul: Number(el('p-dead').value) || 0,
    });
  } catch (error) {
    setStatus(error.message, true);
    return;
  }
  renderAll();
  setStatus('plate settings saved');
}

function setStatus(text, bad = false) {
  const status = el('status');
  status.textContent = text;
  status.classList.toggle('is-bad', Boolean(bad) && Boolean(text));
}

async function loadFolder() {
  try {
    const health = await fetch('/health').then((r) => r.json());
    el('folder').textContent = `${health.plasmids} plasmids`;
  } catch {
    el('folder').textContent = '';
  }
}

/* -------------------------------------------------------------------- wire */

el('editor-save').addEventListener('click', saveWell);
el('editor-clear').addEventListener('click', clearWell);
el('settings-save').addEventListener('click', saveSettings);
el('e-plasmid').addEventListener('keydown', (event) => {
  if (event.key === 'Enter') saveWell();
});
el('e-volume').addEventListener('keydown', (event) => {
  if (event.key === 'Enter') saveWell();
});

loadFolder();
load();

/* ------------------------------------------------------------- the sweep */

function slotFor(key) {
  return state.slots.find((s) => s.key === key) || null;
}

function renderFactors() {
  for (const axis of ['row', 'column']) {
    const short = axis === 'row' ? 'row' : 'col';
    const select = el(`${short}-position`);
    select.replaceChildren(...[
      Object.assign(document.createElement('option'), { value: '', textContent: 'nothing varies' }),
      ...state.slots.filter((s) => s.options.length).map((slot) => Object.assign(
        document.createElement('option'),
        { value: slot.key, textContent: `${slot.label} · ${slot.description}` },
      )),
    ]);
    select.value = state[axis].position;
    // the axis header wears the position's own part-type colour rather than a
    // new categorical ramp: a second 12-colour scale would fight the eight
    // part colours that already carry meaning everywhere else
    const badge = el(`${short}-badge`);
    badge.dataset.part = state[axis].position || '';
    badge.textContent = state[axis].position || '—';
    renderPicks(axis);
  }
  preview();
}

function renderPicks(axis) {
  const short = axis === 'row' ? 'row' : 'col';
  const host = el(`${short}-picks`);
  const slot = slotFor(state[axis].position);
  host.replaceChildren();
  if (!slot) {
    host.append(note('choose a position to vary'));
    return;
  }
  const chosen = new Set(state[axis].candidates);
  const limit = LIMIT[axis];

  const count = note(`${chosen.size} of ${limit} chosen · ${slot.options.length} fit this position`);
  count.className = 'picks-count';
  host.append(count);

  for (const option of slot.options) {
    const pick = document.createElement('button');
    pick.type = 'button';
    pick.className = 'pick';
    pick.textContent = option.component || option.name;
    pick.title = `${option.name}\n${option.length.toLocaleString()} bp`;
    if (chosen.has(option.name)) pick.classList.add('is-on');
    pick.setAttribute('aria-pressed', String(chosen.has(option.name)));
    pick.addEventListener('click', () => toggle(axis, option.name));
    host.append(pick);
  }
}

function note(text) {
  const span = document.createElement('span');
  span.className = 'card-note';
  span.textContent = text;
  return span;
}

function toggle(axis, name) {
  const list = state[axis].candidates;
  const at = list.indexOf(name);
  if (at >= 0) list.splice(at, 1);
  else if (list.length < LIMIT[axis]) list.push(name);
  else {
    setStatus(`the ${axis} axis holds ${LIMIT[axis]}; untick one first`, true);
    return;
  }
  setStatus('');
  renderPicks(axis);
  preview();
  schedule();
}

function preview() {
  const rowName = state.row.candidates[0] || 'row1';
  const colName = state.column.candidates[0] || 'col1';
  el('pattern-preview').textContent = 'A1 → ' + (el('sweep-pattern').value || '{base}_{row}_{col}')
    .replace('{base}', el('sweep-name').value || 'SWEEP')
    .replace('{row}', shortName(rowName))
    .replace('{col}', shortName(colName))
    .replace('{well}', 'A1');
}

function shortName(name) {
  const found = state.candidates.find((c) => c.name === name);
  const label = (found && found.component) || name;
  return label.split('·')[0].replace(/[^A-Za-z0-9]+/g, '').slice(0, 16) || name;
}

/* Debounced: ticking through eight promoters is eight clicks, and running 96
   in-silico assemblies after each one is 768 assemblies to show one answer. */
function schedule() {
  clearTimeout(state.timer);
  state.timer = setTimeout(runSweep, 220);
}

async function runSweep() {
  const design = {
    base: state.base,
    name: el('sweep-name').value || 'SWEEP',
    pattern: el('sweep-pattern').value || '{base}_{row}_{col}',
    row: { position: state.row.position, candidates: state.row.candidates },
    column: { position: state.column.position, candidates: state.column.candidates },
  };
  if (!state.row.candidates.length && !state.column.candidates.length) {
    state.sweep = null;
    renderDest();
    return;
  }
  try {
    state.sweep = await post('/api/plate/sweep', design);
  } catch (error) {
    setStatus(error.message, true);
    return;
  }
  renderDest();
}

/* ------------------------------------------------------ destination grid */

function renderDest() {
  const sweep = state.sweep;
  const blocker = el('blocker');
  const grid = el('dest-grid');

  if (!sweep) {
    blocker.hidden = false;
    blocker.textContent = 'Pick what varies, and the plate fills in.';
    grid.replaceChildren();
    el('dest-cols').replaceChildren();
    el('dest-rows').replaceChildren();
    el('strip').replaceChildren();
    el('generate').disabled = true;
    return;
  }

  const errors = sweep.issues.filter((i) => i.level === 'error');
  blocker.hidden = !errors.length;
  // the reason sits directly above the grid, not in a toast: it is the answer
  // to "why is Generate greyed out", and that question is asked here
  blocker.textContent = errors.map((i) => i.message).join(' · ');

  const byWell = new Map(sweep.wells.map((w) => [w.well, w]));
  const rows = state.row.candidates.length || 1;
  const columns = state.column.candidates.length || (state.row.candidates.length ? 1 : 0);

  el('dest-cols').style.gridTemplateColumns = `repeat(${Math.max(1, columns)}, minmax(0,1fr))`;
  grid.style.gridTemplateColumns = `repeat(${Math.max(1, columns)}, minmax(0,1fr))`;
  el('dest-rows').style.gridTemplateRows = `repeat(${Math.max(1, rows)}, minmax(0,1fr))`;

  el('dest-cols').replaceChildren(...(sweep.axes.columns || []).map((axis, n) => {
    const cell = document.createElement('span');
    cell.className = 'axis-label';
    cell.dataset.part = state.column.position || '';
    cell.textContent = axis.variant ? shortName(axis.variant) : `col ${n + 1}`;
    if (axis.shared.length) {
      cell.classList.add('is-flagged');
      cell.title = `column ${n + 1} · ${axis.variant} · ${axis.err} of ${axis.wells} failed · ${axis.shared.join(', ')}`;
    }
    return cell;
  }));

  el('dest-rows').replaceChildren(...(sweep.axes.rows || []).map((axis, n) => {
    const cell = document.createElement('span');
    cell.className = 'axis-row';
    cell.dataset.part = state.row.position || '';
    const letter = document.createElement('b');
    letter.textContent = state.geometry.rows[n];
    const name = document.createElement('span');
    name.textContent = axis.variant ? shortName(axis.variant) : '';
    cell.append(name, letter);
    if (axis.shared.length) {
      cell.classList.add('is-flagged');
      cell.title = `row ${state.geometry.rows[n]} · ${axis.variant} · ${axis.err} of ${axis.wells} failed · ${axis.shared.join(', ')}`;
    }
    return cell;
  }));

  const cells = [];
  for (let r = 0; r < rows; r += 1) {
    for (let c = 0; c < Math.max(1, columns); c += 1) {
      const name = `${state.geometry.rows[r]}${c + 1}`;
      const well = byWell.get(name);
      const cell = document.createElement('button');
      cell.type = 'button';
      cell.className = 'dest-well';
      cell.setAttribute('role', 'gridcell');
      cell.dataset.well = name;

      if (!well) {
        cell.classList.add('is-empty');
        cell.setAttribute('aria-label', `${name}, not used by this sweep`);
      } else {
        cell.classList.add(`is-${well.status}`);
        const trouble = well.issues.filter((i) => i.level !== 'info');
        cell.setAttribute('aria-label',
          `${name} — ${well.name}, ${well.length ? `${well.length} bp` : 'does not assemble'}`
          + (trouble.length ? `, ${trouble.length} ${trouble.length === 1 ? 'problem' : 'problems'}` : ''));
        cell.title = `${name} · ${well.name}\n${well.length ? `${well.length.toLocaleString()} bp` : 'does not assemble'}`
          + (trouble.length ? `\n${trouble.map((i) => i.message).join('\n')}` : '');
        cell.append(glyph(well.status));
        cell.addEventListener('click', () => showWell(well));
      }
      cell.addEventListener('keydown', (event) => moveDest(event, r, c, rows, columns));
      cells.push(cell);
    }
  }
  grid.replaceChildren(...cells);

  renderStrip(sweep);
  for (const id of ['generate', 'map-csv', 'products']) el(id).disabled = !sweep.ok;
  // a protocol shown for one design must not sit on screen looking current
  // after the design has changed under it
  el('protocol').hidden = true;
}

/* Never colour alone. A warning and an error are told apart by shape as well,
   for the same reason this app avoids a red/green pair anywhere else. */
function glyph(status) {
  const mark = document.createElement('i');
  mark.className = `glyph glyph-${status}`;
  mark.textContent = status === 'err' ? '×' : status === 'warn' ? '!' : '';
  return mark;
}

function moveDest(event, r, c, rows, columns) {
  const keys = { ArrowLeft: [0, -1], ArrowRight: [0, 1], ArrowUp: [-1, 0], ArrowDown: [1, 0] };
  const step = keys[event.key];
  if (!step) return;
  event.preventDefault();
  const nr = Math.min(rows - 1, Math.max(0, r + step[0]));
  const nc = Math.min(Math.max(1, columns) - 1, Math.max(0, c + step[1]));
  const target = el('dest-grid').querySelector(`[data-well="${state.geometry.rows[nr]}${nc + 1}"]`);
  if (target) target.focus();
}

function renderStrip(sweep) {
  const strip = el('strip');
  strip.replaceChildren();
  for (const [kind, label] of [['ok', 'clean'], ['warn', 'warning'], ['err', 'error']]) {
    const count = sweep.counts[kind];
    const chip = document.createElement('span');
    chip.className = `chip chip-${kind}`;
    chip.textContent = `${count} ${label}${count === 1 ? '' : 's'}`;
    strip.append(chip);
  }
  const mix = document.createElement('span');
  mix.className = 'card-note';
  mix.textContent = sweep.shared_positions && sweep.shared_positions.length
    ? `${sweep.shared_positions.length} of ${sweep.shared_positions.length + (state.row.candidates.length ? 1 : 0) + (state.column.candidates.length ? 1 : 0)} positions are the same in every well — those go in the master mix`
    : '';
  strip.append(mix);
}

function showWell(well) {
  const box = el('well-detail');
  box.hidden = false;
  box.replaceChildren();

  const head = document.createElement('div');
  head.className = 'detail-head';
  const where = document.createElement('b');
  where.textContent = well.well;
  const name = document.createElement('span');
  name.textContent = well.name;
  const size = document.createElement('span');
  size.className = 'card-note';
  size.textContent = well.length ? `${well.length.toLocaleString()} bp` : 'does not assemble';
  head.append(where, name, size);
  box.append(head);

  const parts = document.createElement('ul');
  parts.className = 'detail-parts';
  for (const slot of state.slots) {
    const chosen = well.selections[slot.key];
    if (!chosen) continue;
    const item = document.createElement('li');
    const badge = document.createElement('span');
    badge.className = 'badge';
    badge.dataset.part = slot.key;
    badge.textContent = slot.label;
    const what = document.createElement('span');
    what.textContent = chosen;
    const varies = document.createElement('span');
    varies.className = 'card-note';
    varies.textContent = slot.key === state.row.position ? 'row factor'
      : slot.key === state.column.position ? 'column factor' : 'master mix';
    item.append(badge, what, varies);
    parts.append(item);
  }
  box.append(parts);

  for (const issue of well.issues.filter((i) => i.level !== 'info')) {
    const row = document.createElement('div');
    row.className = `issue is-${issue.level}`;
    row.textContent = issue.message;
    box.append(row);
  }
}

el('row-position').addEventListener('change', () => {
  state.row.position = el('row-position').value;
  state.row.candidates = [];
  renderFactors();
  schedule();
});
el('col-position').addEventListener('change', () => {
  state.column.position = el('col-position').value;
  state.column.candidates = [];
  renderFactors();
  schedule();
});
el('sweep-name').addEventListener('input', () => { preview(); schedule(); });
el('sweep-pattern').addEventListener('input', () => { preview(); schedule(); });
el('revalidate').addEventListener('click', runSweep);
/* ------------------------------------------------------------- outputs */

function design() {
  return {
    base: state.base,
    name: el('sweep-name').value || 'SWEEP',
    pattern: el('sweep-pattern').value || '{base}_{row}_{col}',
    row: { position: state.row.position, candidates: state.row.candidates },
    column: { position: state.column.position, candidates: state.column.candidates },
  };
}

async function download(path, filename) {
  const response = await fetch(path, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(design()),
  });
  if (!response.ok) {
    const body = await response.json().catch(() => ({}));
    setStatus(body.detail || `${response.status} ${response.statusText}`, true);
    return;
  }
  const url = URL.createObjectURL(await response.blob());
  const link = document.createElement('a');
  link.href = url;
  link.download = filename;
  link.click();
  URL.revokeObjectURL(url);
  setStatus(`${filename} downloaded`);
}

/* Generating is where the protocol is simulated, because the point of
 * simulating is to decide whether a download should be offered at all. A
 * protocol that has never been executed is exactly the one that fails at the
 * bench with the tips already on. */
async function generate() {
  const button = el('generate');
  button.disabled = true;
  setStatus('generating and simulating\u2026');
  let result;
  try {
    result = await post('/api/plate/protocol', design());
  } catch (error) {
    setStatus(error.message, true);
    button.disabled = false;
    return;
  }
  setStatus('');
  button.disabled = false;

  el('protocol').hidden = false;
  el('protocol-source').textContent = result.source;

  const sim = result.simulation;
  const verdict = el('sim-verdict');
  verdict.textContent = `simulation ${sim.verdict}`;
  verdict.className = `sim sim-${sim.verdict.replace(/\s+/g, '-')}`;
  el('sim-detail').textContent = sim.detail
    || `${result.lines.toLocaleString()} lines \u00b7 ${result.transfers} transfers`;

  el('download').disabled = !result.downloadable;
  el('commit').disabled = !result.downloadable;
}

el('generate').addEventListener('click', generate);
el('map-csv').addEventListener('click',
  () => download('/api/plate/map.csv', `${el('sweep-name').value || 'sweep'}-plate.csv`));
el('products').addEventListener('click',
  () => download('/api/plate/products.zip', `${el('sweep-name').value || 'sweep'}-products.zip`));
el('download').addEventListener('click',
  () => download('/api/plate/protocol.py', `${el('sweep-name').value || 'sweep'}-protocol.py`));

el('commit').addEventListener('click', async () => {
  try {
    const result = await post('/api/plate/commit', design());
    state.plate = result.plate;
    renderAll();
    setStatus(`${result.logged} wells logged; the source plate has been drawn down`);
    el('commit').disabled = true;
  } catch (error) {
    setStatus(error.message, true);
  }
});
