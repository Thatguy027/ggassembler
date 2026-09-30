/* Part library.
 *
 * Reads only the shared /api/library/* endpoints - it is not a level and has no
 * endpoints of its own. Filtering happens in the browser over the indexed rows,
 * because the whole index is a few hundred rows and a round trip per keystroke
 * would be slower and worse.
 *
 * The one thing this screen writes is a manual type assignment, which the
 * server stores in .ggasm/overrides.json and which always beats detection.
 */

const state = { rows: [], types: [], filter: '', type: '', confidence: '' };

/* Why a plasmid has no part type. The detector already works this out; these
 * are the keys it groups under, so "83 unrecognised" becomes three separate,
 * actionable piles instead of one undifferentiated warning. */
const TRIAGE_FILTERS = ['no_part_pair', 'uncuttable', 'linear'];

let filterTimer = null;
const el = (id) => document.getElementById(id);

async function getJSON(path, options) {
  const response = await fetch(path, options);
  if (!response.ok) throw new Error(`${path}: ${response.status}`);
  return response.json();
}

// ------------------------------------------------------------------ stats ---

function renderStats(summary) {
  const cards = [
    [summary.unique ?? summary.total,
     summary.merged ? `unique plasmids (${summary.merged} duplicate files merged)` : '.gb / .gbk files scanned',
     ''],
    [summary.by_digest, `typed by ${summary.enzymes.part}/${summary.enzymes.multigene} digest`, 'stat-ok'],
    [summary.composite, 'composite (234, 678, 2·3·4…)', ''],
    [summary.unrecognised, 'unrecognised — need a manual type', summary.unrecognised ? 'stat-attention' : ''],
  ];
  el('stats').replaceChildren(
    ...cards.map(([value, label, cls]) => {
      const card = document.createElement('div');
      card.className = 'stat';
      const v = document.createElement('span');
      v.className = `stat-value ${cls}`.trim();
      v.textContent = value.toLocaleString();
      const l = document.createElement('span');
      l.className = 'stat-label';
      l.textContent = label;
      card.append(v, l);
      return card;
    }),
  );

  const roots = summary.roots.length === 1 ? summary.roots[0] : `${summary.roots.length} folders`;
  el('folder').textContent = roots;
  el('folder').title = summary.roots.join('\n');
  el('overhang-head').textContent = `Overhangs (${summary.enzymes.part})`;
}

// ------------------------------------------------------------------- rows ---

/** What the Internal sites column says, expected sites and all. */
function sitesText(row) {
  const parts = [];
  if (row.internal_multigene_positions && row.internal_multigene_positions.length) {
    const where = row.internal_multigene_positions.map((p) => p.toLocaleString()).join(', ');
    parts.push([`${row.enzymes.multigene} internal @ ${where}`, true]);
  }
  if (row.reversed_sites) parts.push([`${row.enzymes.part} reversed`, true]);
  if (row.internal_sites.part_enzyme) {
    parts.push([`${row.enzymes.part} ×${row.internal_sites.part_enzyme} internal`, true]);
  }
  if (row.connector_overhang) parts.push([`${row.enzymes.multigene} ×1 (expected)`, false]);
  else if (row.internal_sites.multigene_enzyme) {
    parts.push([`${row.enzymes.multigene} ×${row.internal_sites.multigene_enzyme}`, false]);
  }
  if (row.site_totals.linearizer) {
    parts.push([`${row.enzymes.linearizer} ×${row.site_totals.linearizer}`, false]);
  }
  if (!parts.length) return [['none', false]];
  return parts;
}

/** The Notes column: the most useful sentence we have about this row. */
function noteText(row) {
  if (row.conflict) return [row.conflict, true];
  if (row.internal_sites.part_enzyme) {
    return [`Internal ${row.enzymes.part} — domesticate before use`, true];
  }
  if (row.reversed_sites) return ['Reversed sites — omit the final digest step', false];
  if (row.unrecognised) return ['Assign a type manually', true];
  if (row.roles.includes('cassette')) {
    return [`Cassette ${row.cassette_overhangs.join(' → ')}`, false];
  }
  if (row.roles.includes('multigene_vector')) return ['Multigene destination vector', false];
  if (row.roles.includes('entry_vector')) return ['Part entry vector', false];
  if (row.connector_overhang) return [`Connector overhang ${row.connector_overhang}`, false];
  return [row.reason, false];
}

function evidenceCell(row) {
  const span = document.createElement('span');
  span.className = 'evidence';
  if (row.source === 'manual') {
    span.classList.add('evidence-manual');
    span.textContent = 'set by hand';
  } else if (row.confidence === 'digest') {
    span.classList.add('evidence-digest');
    span.textContent = row.internal_sites.part_enzyme ? 'digest, low conf.' : 'digest';
  } else if (row.confidence === 'none') {
    span.classList.add('evidence-weak');
    span.textContent = 'unclassified';
  } else {
    span.classList.add('evidence-weak');
    span.textContent = row.confidence;
  }
  span.title = row.reason;
  return span;
}

function typePill(row) {
  const pill = document.createElement('button');
  pill.type = 'button';
  pill.className = 'type-pill';
  if (row.source === 'manual') pill.classList.add('type-pill-manual');
  pill.textContent = row.part_type || '—';
  pill.style.background = row.part_type ? row.badge_bg : 'var(--unknown-bg)';
  pill.style.color = row.part_type ? row.badge_fg : 'var(--unknown-fg)';
  pill.title = `${row.reason}\n\nClick to assign a type by hand`;
  pill.setAttribute('aria-label', `Type ${row.part_type || 'unassigned'} for ${row.name}. Assign by hand.`);
  pill.addEventListener('click', () => openAssign(row));
  return pill;
}

function cell(text, className) {
  const td = document.createElement('td');
  if (className) td.className = className;
  td.textContent = text;
  return td;
}

/* What a prep measured at, in ng/uL.
 *
 * This is the one number a reaction setup needs that no GenBank file can
 * supply - it is a property of the tube, not of the sequence - so it is
 * entered here and kept in .ggasm/overrides.json beside the index. Without it
 * the Protocol button can show masses but not volumes.
 */
function concentrationCell(row) {
  const td = document.createElement('td');
  const input = document.createElement('input');
  input.className = 'conc-input';
  input.type = 'number';
  input.min = '0.1';
  input.step = '0.1';
  input.placeholder = '\u2014';
  input.value = row.conc_ng_ul === null || row.conc_ng_ul === undefined ? '' : row.conc_ng_ul;
  input.title = `what a prep of ${row.name} measured at, in ng/\u00b5L`;
  input.addEventListener('change', async () => {
    const raw = input.value.trim();
    const value = raw === '' ? null : Number(raw);
    if (value !== null && !(value > 0)) {
      input.value = row.conc_ng_ul ?? '';
      return;
    }
    input.disabled = true;
    try {
      const updated = await getJSON('/api/library/concentration', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ path: row.path, conc_ng_ul: value }),
      });
      row.conc_ng_ul = updated.conc_ng_ul ?? null;
    } finally {
      input.disabled = false;
    }
  });
  td.append(input);
  return td;
}

function renderRow(row) {
  const tr = document.createElement('tr');
  if (row.needs_attention) tr.classList.add('attention');

  tr.append(cell(row.path, 'cell-mono cell-name'));

  const name = document.createElement('td');
  name.className = 'cell-name';
  name.append(document.createTextNode(row.name));
  if (row.aliases && row.aliases.length) {
    /* The same sequence under several names. Which one speaks for the group
     * was a guess - shortest name in the first folder - and a guess is often
     * not the name the lab uses, so it can be chosen here instead and every
     * dropdown stops reading `pYTK003 / A3_ConL1`. */
    const pick = document.createElement('select');
    pick.className = 'cell-canonical';
    pick.title = 'which name this sequence goes by everywhere in the app';
    for (const option of [row.name, ...row.aliases]) {
      const item = document.createElement('option');
      item.value = option;
      item.textContent = option;
      if (option === row.name) item.selected = true;
      pick.append(item);
    }
    pick.addEventListener('change', async () => {
      pick.disabled = true;
      try {
        await getJSON('/api/library/canonical', {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({ name: pick.value }),
        });
        await load();
      } finally {
        pick.disabled = false;
      }
    });
    name.append(pick);
  }
  tr.append(name);

  const type = document.createElement('td');
  type.append(typePill(row));
  tr.append(type);

  const overhangs = document.createElement('td');
  overhangs.className = 'cell-mono';
  if (row.five_prime && row.three_prime) {
    overhangs.textContent = `${row.five_prime} → ${row.three_prime}`;
  } else if (row.cassette_overhangs) {
    overhangs.textContent = `${row.cassette_overhangs[0]} \u2192 ${row.cassette_overhangs[1]}`;
    overhangs.title = `${row.enzymes.multigene} connector ends`;
    const tag = document.createElement('span');
    tag.className = 'cell-dim';
    tag.style.display = 'block';
    tag.style.fontSize = '11px';
    tag.textContent = `${row.enzymes.multigene} connectors`;
    overhangs.append(tag);
  } else if (row.multigene_overhangs && row.multigene_overhangs.length) {
    overhangs.classList.add('cell-dim');
    overhangs.textContent = row.multigene_overhangs.join(' \u00b7 ');
    overhangs.title = `no ${row.enzymes.part} pair; these are the ${row.enzymes.multigene} overhangs`;
    const tag = document.createElement('span');
    tag.className = 'cell-attention';
    tag.style.display = 'block';
    tag.style.fontSize = '11px';
    tag.textContent = `no ${row.enzymes.part} pair`;
    overhangs.append(tag);
  } else {
    overhangs.classList.add('cell-attention');
    overhangs.textContent = `no ${row.enzymes.part} pair found`;
  }
  tr.append(overhangs);

  const sites = document.createElement('td');
  sites.className = 'cell-mono';
  for (const [text, warn] of sitesText(row)) {
    const span = document.createElement('span');
    span.className = warn ? 'cell-attention' : 'cell-dim';
    span.textContent = text;
    span.style.display = 'block';
    sites.append(span);
  }
  tr.append(sites);

  const evidence = document.createElement('td');
  evidence.append(evidenceCell(row));
  tr.append(evidence);

  // A description read off the file name is shown in italics: it is a guess
  // from the name, not something the file says, and these are exactly the
  // plasmids worth going back and annotating properly.
  const described = cell(row.component || '', 'cell-dim cell-component');
  if (row.component && row.component_source === 'filename') {
    described.classList.add('is-guessed');
    described.title = 'read from the file name — the GenBank file has no annotation here';
  }
  tr.append(described);

  tr.append(concentrationCell(row));

  const [note, warn] = noteText(row);
  tr.append(cell(note, warn ? 'cell-attention' : 'cell-dim'));
  tr.lastChild.style.fontSize = '12px';
  return tr;
}

function visible() {
  const needle = state.filter.trim().toLowerCase();
  return state.rows.filter((row) => {
    if (state.type && row.part_type !== state.type) return false;
    if (state.confidence === 'manual' && row.source !== 'manual') return false;
    if (state.confidence === 'attention' && !row.needs_attention) return false;
    if (state.confidence === 'internal-multigene'
        && !(row.internal_multigene_positions || []).length) return false;
    // the two description filters answer "what still needs annotating?"
    if (state.confidence === 'described-by-name'
        && row.component_source !== 'filename') return false;
    if (state.confidence === 'undescribed' && row.component) return false;
    // type 1 and 5 parts that lost their multigene site: they assemble at
    // Level 2 and make a cassette with no ends
    if (state.confidence === 'blocks-multigene' && row.level3_ready !== false) return false;
    // triage: why a plasmid has no type. Each is a different job to fix.
    if (TRIAGE_FILTERS.includes(state.confidence)
        && row.triage !== state.confidence) return false;
    if (state.confidence === 'assumed-circular' && !row.assumed_circular) return false;
    if (state.confidence === 'merged' && !(row.aliases || []).length) return false;
    if (
      state.confidence &&
      !['manual', 'attention', 'internal-multigene', 'described-by-name',
        'undescribed', 'assumed-circular', 'merged', 'blocks-multigene',
        ...TRIAGE_FILTERS]
        .includes(state.confidence) &&
      row.confidence !== state.confidence
    ) return false;
    if (!needle) return true;
    return [
      row.name, row.path, row.part_type, row.five_prime, row.three_prime,
      row.connector_overhang, row.ecoli_marker, row.component,
      ...(row.aliases || []), ...row.roles,
    ].some((field) => field && String(field).toLowerCase().includes(needle));
  });
}

/* How many rows go into the DOM at once.
 *
 * Every row is nine cells plus a concentration input, and merged rows add a
 * select on top of that, so the whole table is several thousand nodes. It used
 * to be rebuilt on every keystroke in the filter above it, which is exactly
 * the interaction you do most. So the table renders a chunk at a time and
 * grows as you scroll, and the filter is debounced: the work now scales with
 * what you can see rather than with what you own.
 */
const CHUNK = 120;

let shownCount = CHUNK;
let watcher = null;

function renderRows({ keepCount = false } = {}) {
  const rows = visible();
  if (!keepCount) shownCount = CHUNK;

  el('count').textContent =
    rows.length === state.rows.length
      ? `${state.rows.length} plasmids`
      : `${rows.length} of ${state.rows.length} plasmids`;

  if (watcher) { watcher.disconnect(); watcher = null; }

  if (!rows.length) {
    const tr = document.createElement('tr');
    const td = document.createElement('td');
    td.colSpan = 9;
    td.className = 'empty';
    td.textContent = 'Nothing matches those filters.';
    tr.append(td);
    el('rows').replaceChildren(tr);
    return;
  }

  const slice = rows.slice(0, shownCount);
  const body = el('rows');
  body.replaceChildren(...slice.map(renderRow));

  if (rows.length <= slice.length) return;

  // a sentinel row: when it scrolls into view, take the next chunk
  const tr = document.createElement('tr');
  const td = document.createElement('td');
  td.colSpan = 9;
  td.className = 'more';
  td.textContent = `${rows.length - slice.length} more — scroll to load`;
  tr.append(td);
  body.append(tr);

  if (typeof IntersectionObserver === 'undefined') {
    // no observer: make it a button rather than leaving rows unreachable
    td.classList.add('is-button');
    td.addEventListener('click', () => {
      shownCount += CHUNK;
      renderRows({ keepCount: true });
    });
    return;
  }
  watcher = new IntersectionObserver((entries) => {
    if (!entries.some((e) => e.isIntersecting)) return;
    shownCount += CHUNK;
    renderRows({ keepCount: true });
  }, { rootMargin: '200px' });
  watcher.observe(tr);
}

// ----------------------------------------------------------------- assign ---

let assigning = null;

function openAssign(row) {
  assigning = row;
  el('assign-title').textContent = row.name;
  el('assign-detected').textContent = row.part_type
    ? `Detected as ${row.part_type} — ${row.reason}`
    : row.reason;

  const select = el('assign-type');
  select.replaceChildren();
  const blank = document.createElement('option');
  blank.value = '';
  blank.textContent = 'No type';
  select.append(blank);
  for (const type of state.types) {
    const option = document.createElement('option');
    option.value = type.name;
    option.textContent = `${type.name} — ${type.description}`;
    if (type.name === row.part_type) option.selected = true;
    select.append(option);
  }
  el('assign-reason').value = row.source === 'manual' ? row.reason : '';
  el('assign-clear').hidden = row.source !== 'manual';
  el('assign-dialog').showModal();
}

async function saveAssign(partType) {
  if (!assigning) return;
  await getJSON('/api/library/override', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({
      path: assigning.path,
      part_type: partType,
      reason: el('assign-reason').value.trim() || 'set by hand',
    }),
  });
  el('assign-dialog').close();
  assigning = null;
  await load();
}

// ------------------------------------------------------------------ wiring ---

async function load() {
  const [summary, rows, types] = await Promise.all([
    getJSON('/api/library/summary'),
    getJSON('/api/library/plasmids'),
    getJSON('/api/library/types'),
  ]);
  state.rows = rows;
  state.types = types;
  renderStats(summary);

  const filter = el('type-filter');
  const chosen = filter.value;
  filter.replaceChildren();
  const all = document.createElement('option');
  all.value = '';
  all.textContent = 'All part types';
  filter.append(all);
  for (const type of types.filter((t) => t.count)) {
    const option = document.createElement('option');
    option.value = type.name;
    option.textContent = `Type ${type.name} (${type.count})`;
    filter.append(option);
  }
  filter.value = chosen;
  renderRows();
}

el('q').addEventListener('input', (event) => {
  state.filter = event.target.value;
  clearTimeout(filterTimer);
  filterTimer = setTimeout(renderRows, 120);
});
el('type-filter').addEventListener('change', (event) => {
  state.type = event.target.value;
  renderRows();
});
el('conf-filter').addEventListener('change', (event) => {
  state.confidence = event.target.value;
  renderRows();
});
el('rescan-btn').addEventListener('click', async () => {
  const button = el('rescan-btn');
  button.disabled = true;
  button.textContent = 'Scanning…';
  await fetch('/api/library/rescan?force=true', { method: 'POST' });
  await load();
  button.disabled = false;
  button.textContent = 'Re-scan';
});
el('assign-save').addEventListener('click', () => saveAssign(el('assign-type').value || null));
el('assign-clear').addEventListener('click', () => saveAssign(null));

load();
