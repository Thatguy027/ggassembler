/* Level 3 - multigene assembly.
 *
 * Talks only to /api/level3/* and the shared /api/library/*. Imports nothing
 * from a sibling level.
 *
 * The construct is drawn as a bar rather than a ring: at this level what you
 * are reading is an order, and a chain of connectors reads better straightened
 * out. The connector labels come from the server, which learned them from the
 * user's own part plasmids.
 */

const state = { backbone: null, transcription_units: [''], name: 'pMultigene' };
let catalogue = { cassettes: [], backbones: [], enzyme: 'BsmBI', linearizer: 'NotI' };

const el = (id) => document.getElementById(id);

async function post(path, body) {
  const response = await fetch(path, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(body),
  });
  if (!response.ok) throw new Error(`${path}: ${response.status}`);
  return response.json();
}

function note(text) {
  const span = document.createElement('span');
  span.className = 'card-note';
  span.textContent = text;
  return span;
}

// ------------------------------------------------------------------ slots ---

/** The overhang a slot has to start with, given what precedes it. */
function expectedLeft(index) {
  if (index === 0) {
    const backbone = catalogue.backbones.find((b) => b.name === state.backbone);
    return backbone ? backbone.right_overhang : null;
  }
  const previous = catalogue.cassettes.find(
    (c) => c.name === state.transcription_units[index - 1],
  );
  return previous ? previous.right_overhang : null;
}

function optionLabel(cassette) {
  if (cassette.display) return cassette.display;
  const what = cassette.contents ? ` [${cassette.contents}]` : '';
  return `${cassette.name}${what}`;
}

/** Group the cassettes by the connector they start with, fitting group first. */
function fillCassetteOptions(select, chosen, expected) {
  const blank = document.createElement('option');
  blank.value = '';
  blank.textContent = catalogue.cassettes.length
    ? 'Choose a cassette\u2026'
    : 'no cassettes in the library';
  select.append(blank);

  const groups = new Map();
  for (const cassette of catalogue.cassettes) {
    if (!groups.has(cassette.left_overhang)) groups.set(cassette.left_overhang, []);
    groups.get(cassette.left_overhang).push(cassette);
  }

  const keys = [...groups.keys()].sort((a, b) => {
    if (a === expected) return -1;
    if (b === expected) return 1;
    return a.localeCompare(b);
  });

  for (const overhang of keys) {
    const cassettes = groups.get(overhang);
    const group = document.createElement('optgroup');
    const fits = overhang === expected;
    const label = cassettes[0].left_label !== overhang ? ` ${cassettes[0].left_label.replace(overhang, '').trim()}` : '';
    group.label = `${fits ? '\u2713 fits here \u2014 ' : ''}starts at ${overhang}${label}`;
    for (const cassette of cassettes) {
      const option = document.createElement('option');
      option.value = cassette.name;
      option.textContent =
        `${optionLabel(cassette)} \u00b7 \u2192${cassette.right_overhang} \u00b7 ${cassette.unit_length || cassette.length} bp`;
      if (cassette.name === chosen) option.selected = true;
      option.title = [
        cassette.name,
        cassette.aliases && cassette.aliases.length
          ? `also filed as ${cassette.aliases.join(', ')}` : '',
        cassette.contents || 'no annotated contents',
        `${cassette.left_overhang} \u2192 ${cassette.right_overhang}`,
      ].filter(Boolean).join('\n');
      group.append(option);
    }
    select.append(group);
  }
}

function connectorLine(index, chosen, expected) {
  const line = document.createElement('div');
  line.className = 'tu-connectors';
  const cassette = catalogue.cassettes.find((c) => c.name === chosen);
  if (!cassette) {
    if (expected) {
      const hint = document.createElement('span');
      hint.textContent = `needs a cassette starting at `;
      const oh = document.createElement('span');
      oh.className = 'oh';
      oh.textContent = expected;
      line.append(hint, oh);
    }
    return line;
  }

  const left = document.createElement('span');
  left.className = 'oh';
  left.textContent = cassette.left_overhang;
  left.title = cassette.left_label;
  const arrow = document.createElement('span');
  arrow.textContent = '\u2192';
  const right = document.createElement('span');
  right.className = 'oh';
  right.textContent = cassette.right_overhang;
  right.title = cassette.right_label;
  line.append(left, arrow, right);

  if (expected) {
    const verdict = document.createElement('span');
    const fits = cassette.left_overhang === expected;
    verdict.className = fits ? 'tu-fit' : 'tu-misfit';
    verdict.textContent = fits
      ? index === 0 ? 'matches the backbone' : `matches TU${index}`
      : `does not follow ${expected}`;
    line.append(verdict);
  }
  return line;
}

function renderTus() {
  el('tu-list').replaceChildren(
    ...state.transcription_units.map((chosen, index) => {
      const row = document.createElement('div');
      row.className = 'tu';

      const label = document.createElement('span');
      label.className = 'tu-index';
      label.textContent = `TU${index + 1}`;

      const select = document.createElement('select');
      select.id = `tu-${index}`;
      select.setAttribute('aria-label', `Transcription unit ${index + 1}`);
      const expected = expectedLeft(index);
      fillCassetteOptions(select, chosen, expected);
      select.addEventListener('change', () => {
        state.transcription_units[index] = select.value;
        refresh();
      });

      const remove = document.createElement('button');
      remove.type = 'button';
      remove.className = 'tu-remove';
      remove.textContent = '\u00d7';
      remove.title = `Remove TU${index + 1}`;
      remove.setAttribute('aria-label', `Remove transcription unit ${index + 1}`);
      remove.addEventListener('click', () => {
        state.transcription_units.splice(index, 1);
        if (!state.transcription_units.length) state.transcription_units.push('');
        refresh();
      });

      row.append(label, select, remove, connectorLine(index, chosen, expected));
      return row;
    }),
  );
}

/** An annotated map of the selected backbone, with the dropout marked. */
function renderBackboneMap() {
  const box = el('backbone-map');
  const backbone = catalogue.backbones.find((b) => b.name === state.backbone);
  if (!backbone || !backbone.features) {
    box.hidden = true;
    return;
  }
  box.hidden = false;

  const track = document.createElement('div');
  track.className = 'pm-track';
  const pct = (value) => `${(value / backbone.length) * 100}%`;

  // the dropout - what the reaction replaces. It can run past the origin, in
  // which case it is drawn as the two arcs it really is.
  if (backbone.dropout_start !== null && backbone.dropout_end !== null) {
    const { dropout_start: from, dropout_end: to } = backbone;
    const spans = from <= to ? [[from, to]] : [[from, backbone.length], [0, to]];
    for (const [a, b] of spans) {
      if (b <= a) continue;
      const region = document.createElement('div');
      region.className = 'pm-unit';
      region.style.left = pct(a);
      region.style.width = pct(b - a);
      region.title = `dropout, ${backbone.dropout_length.toLocaleString()} bp \u2014 replaced by your cassettes`;
      track.append(region);
    }
  }

  for (const feature of backbone.features) {
    const block = document.createElement('div');
    block.className = 'pm-feature';
    block.style.left = pct(feature.start);
    block.style.width = pct(Math.max(feature.length, backbone.length / 300));
    block.style.background = feature.color;
    block.title = `${feature.label} \u00b7 ${feature.type} \u00b7 ${feature.start.toLocaleString()}\u2013${feature.end.toLocaleString()}`;
    track.append(block);
  }

  const scale = document.createElement('div');
  scale.className = 'pm-scale';
  const left = document.createElement('span');
  left.textContent = `${backbone.name} \u00b7 ${backbone.length.toLocaleString()} bp`;
  const right = document.createElement('span');
  right.textContent = backbone.integration
    ? `${backbone.marker || ''} \u00b7 integration vector`.trim()
    : backbone.marker || '';
  scale.append(left, right);

  const list = document.createElement('div');
  list.className = 'pm-list';
  const shown = [...backbone.features]
    .sort((a, b) => b.length - a.length)
    .slice(0, 10)
    .sort((a, b) => a.start - b.start);
  for (const feature of shown) {
    const item = document.createElement('span');
    item.className = `pm-item${feature.kept ? '' : ' is-unit'}`;
    const swatch = document.createElement('span');
    swatch.className = 'pm-swatch';
    swatch.style.background = feature.color;
    const text = document.createElement('span');
    text.textContent = feature.label;
    item.append(swatch, text);
    item.title = feature.kept
      ? 'kept: on the vector arm that survives the reaction'
      : 'in the dropout - replaced by your cassettes';
    list.append(item);
  }

  box.replaceChildren(track, scale, list);
}

function renderSuggestions(result) {
  const box = el('suggestions');
  if (!result.suggestions.length) {
    box.hidden = true;
    return;
  }
  box.hidden = false;
  box.replaceChildren(
    note('The chain stops before it closes. These cassettes start where it stops:'),
    ...result.suggestions.map((suggestion) => {
      const row = document.createElement('div');
      row.className = 'suggestion';
      const text = document.createElement('span');
      text.textContent = `${suggestion.name} · ${suggestion.left_overhang}→${suggestion.right_overhang}`;
      const button = document.createElement('button');
      button.type = 'button';
      button.textContent = 'Add';
      button.addEventListener('click', () => {
        const empty = state.transcription_units.indexOf('');
        if (empty >= 0) state.transcription_units[empty] = suggestion.name;
        else state.transcription_units.push(suggestion.name);
        refresh();
      });
      row.append(text, button);
      return row;
    }),
  );
}

// -------------------------------------------------------------------- map ---

function renderMap(result) {
  const map = el('map');
  if (!result.ok || !result.parts.length) {
    map.replaceChildren(note('choose a backbone and a cassette'));
    return;
  }

  /* The construct, with each transcription unit opened up into the parts that
   * built it. Drawn inside the bar rather than as a diagram above it: one
   * element means there is no alignment to keep, and the whole thing scales
   * with the page instead of a column of blocks drifting off its segment. */
  const contents = new Map((result.units || []).map((u) => [u.name, u]));

  const bar = document.createElement('div');
  bar.className = 'map-bar';
  for (const part of result.parts) {
    const seg = document.createElement('div');
    seg.className = 'map-seg';
    seg.style.flexGrow = String(part.length);
    seg.title = `${part.source_name} · ${part.length.toLocaleString()} bp · `
      + `${part.left_overhang}→${part.right_overhang}`;

    const unit = contents.get(part.source_name);
    const inside = (unit && unit.parts) || [];

    if (inside.length) {
      const row = document.createElement('div');
      row.className = 'map-parts';
      for (const piece of inside) {
        const block = document.createElement('div');
        block.className = 'map-part';
        if (!piece.known) block.classList.add('is-gap');
        // proportional within the segment, so a gene reads as a gene and a
        // terminator as a sliver - the same scale as the construct itself
        block.style.flexGrow = String(Math.max(piece.length, 1));
        block.dataset.part = piece.part_type || '';
        block.style.background = 'var(--part-color)';
        block.textContent = piece.short || piece.component || piece.name || '';
        block.title = `type ${piece.part_type} · `
          + `${piece.component || piece.name || 'not in the library'} · `
          + `${piece.length.toLocaleString()} bp`;
        row.append(block);
      }
      seg.append(row);
    } else {
      seg.classList.add('is-plain');
      seg.dataset.part = part.part_type || '';
      seg.style.background = 'var(--part-color)';
    }

    const label = document.createElement('span');
    label.className = 'map-seg-name';
    label.textContent = part.source_name;
    seg.append(label);
    bar.append(seg);
  }

  const scale = document.createElement('div');
  scale.className = 'map-scale';
  const left = document.createElement('span');
  left.textContent = '0';
  const right = document.createElement('span');
  right.textContent = `${result.length.toLocaleString()} bp · ${result.parts.length} fragments`;
  scale.append(left, right);

  map.replaceChildren(bar, scale);
}

function renderChain(result) {
  const strip = el('chain-strip');
  if (!result.chain.length) {
    strip.replaceChildren();
    return;
  }
  const nodes = [];
  result.chain.forEach((link, index) => {
    if (index) {
      const arrow = document.createElement('span');
      arrow.className = 'chain-arrow';
      arrow.textContent = '→';
      nodes.push(arrow);
    }
    const chip = document.createElement('span');
    chip.className = 'chain-link';
    const oh = document.createElement('span');
    oh.className = 'oh';
    oh.textContent = link.left;
    oh.title = link.left_label;
    const name = document.createElement('span');
    name.textContent = link.role;
    name.title = link.name;
    chip.append(oh, name);
    nodes.push(chip);
  });
  const last = result.chain[result.chain.length - 1];
  const close = document.createElement('span');
  close.className = 'oh';
  close.textContent = last.right;
  close.title = `${last.right_label} — closes back to the backbone`;
  nodes.push(document.createElement('span'), close);
  strip.replaceChildren(...nodes);
}

function renderIssues(result) {
  el('issues').replaceChildren(
    ...(result.issues.length
      ? result.issues.map((issue) => {
          const row = document.createElement('div');
          row.className = `issue issue-${issue.level}`;
          const code = document.createElement('span');
          code.className = 'issue-code';
          code.textContent = issue.code;
          const message = document.createElement('span');
          message.textContent = issue.message;
          row.append(code, message);
          return row;
        })
      : [note('nothing to report')]),
  );

  const box = el('primers-box');
  box.hidden = !result.check_primers.length;
  el('primers').replaceChildren(
    ...result.check_primers.map((primer) => {
      const row = document.createElement('tr');
      for (const [text, mono] of [
        [primer.name, false],
        [primer.sequence, true],
        [`${primer.tm} °C`, true],
      ]) {
        const cell = document.createElement('td');
        if (mono) cell.className = 'mono';
        cell.textContent = text;
        row.append(cell);
      }
      return row;
    }),
  );
}

// ---------------------------------------------------------------- wiring ---

let pending = null;

async function refresh() {
  state.name = el('name').value.trim() || 'multigene';
  state.backbone = el('backbone').value || null;
  renderTus();

  const mine = {};
  pending = mine;
  const result = await post('/api/level3/assemble', state);
  if (pending !== mine) return;

  renderMap(result);
  renderChain(result);
  renderIssues(result);
  renderSuggestions(result);

  el('subtitle').textContent = result.ok
    ? `${result.length.toLocaleString()} bp · ${result.counts.units} TUs · ${result.enzyme}${result.integration ? ' · integrating' : ''}`
    : ' ';
  el('export-btn').disabled = !result.ok;
  el('save-btn').disabled = !result.ok;

  const backbone = catalogue.backbones.find((b) => b.name === state.backbone);
  el('backbone-note').textContent = backbone
    ? `chain must run ${backbone.right_overhang} → … → ${backbone.left_overhang}`
    : '';
}

// ------------------------------------------------------------------ finder ---

/* The same shared search (/api/library/search), read the Level 3 way.
 *
 * Here a hit is not a part but a whole transcription unit, and where it goes is
 * decided by its connectors rather than by a position: it belongs in the first
 * slot whose expected left overhang it matches, and failing that at the end of
 * the chain. A multigene vector is offered as the backbone instead. Nothing is
 * changed until the question is answered.
 */

const finder = { hits: [], active: -1, chosen: null, seq: 0, suppressed: 0 };

function highlight(text, terms) {
  const fragment = document.createDocumentFragment();
  const pattern = terms
    .filter(Boolean)
    .map((t) => t.replace(/[.*+?^${}()|[\]\\]/g, '\\$&'))
    .join('|');
  if (!pattern) {
    fragment.append(text);
    return fragment;
  }
  let last = 0;
  for (const match of text.matchAll(new RegExp(`(${pattern})`, 'ig'))) {
    if (match.index > last) fragment.append(text.slice(last, match.index));
    const mark = document.createElement('mark');
    mark.textContent = match[0];
    fragment.append(mark);
    last = match.index + match[0].length;
  }
  if (last < text.length) fragment.append(text.slice(last));
  return fragment;
}

/** How a hit could be used here: as the backbone, as a TU, or not at all. */
function roleOf(hit) {
  if (catalogue.backbones.some((b) => b.name === hit.name)) return 'backbone';
  if (catalogue.cassettes.some((c) => c.name === hit.name)) return 'unit';
  return null;
}

/** The first empty-or-any slot whose expected left overhang this cassette fits. */
function slotFor(name) {
  const cassette = catalogue.cassettes.find((c) => c.name === name);
  if (!cassette) return state.transcription_units.length;
  for (let i = 0; i < state.transcription_units.length; i += 1) {
    if (expectedLeft(i) === cassette.left_overhang) return i;
  }
  const empty = state.transcription_units.indexOf('');
  return empty === -1 ? state.transcription_units.length : empty;
}

function hitRow(hit, index, terms) {
  const row = document.createElement('button');
  row.type = 'button';
  row.className = 'finder-hit';
  row.setAttribute('role', 'option');
  row.setAttribute('aria-selected', String(index === finder.active));

  const use = roleOf(hit);
  const mark = document.createElement('span');
  mark.className = 'badge';
  mark.textContent = use === 'backbone' ? 'BB' : use === 'unit' ? 'TU' : '—';
  mark.dataset.part = hit.part_type || '';
  row.append(mark);

  const what = document.createElement('span');
  what.className = 'finder-what';
  const name = document.createElement('span');
  name.className = 'finder-name';
  name.append(highlight(hit.display || hit.name, terms));
  const why = document.createElement('span');
  why.className = 'finder-why';
  why.append(highlight(hit.matched.slice(0, 3).join(' · ') || hit.component || hit.name, terms));
  const ends = hit.cassette_overhangs;
  const tail = document.createElement('span');
  tail.textContent = ends ? ` — ${ends[0]} → ${ends[1]}` : ' — not a multigene piece';
  why.append(tail);
  what.append(name, why);
  row.append(what);

  const len = document.createElement('span');
  len.className = 'finder-len';
  len.textContent = `${hit.length} bp`;
  row.append(len);

  if (!use) {
    row.classList.add('is-unusable');
    row.disabled = true;
    row.title = hit.part_type
      ? `a type ${hit.part_type} part — use it on the Cassette page`
      : 'no connector pair could be read out of this plasmid';
  } else {
    row.addEventListener('click', () => chooseHit(index));
  }
  return row;
}

function confirmBar(hit) {
  const bar = document.createElement('div');
  bar.className = 'finder-confirm';
  const use = roleOf(hit);

  const ask = document.createElement('span');
  ask.className = 'finder-ask';
  const name = document.createElement('b');
  name.textContent = hit.display || hit.name;
  if (use === 'backbone') {
    ask.append('Use ', name, ' as the destination backbone?');
  } else {
    const slot = slotFor(hit.name);
    const taken = state.transcription_units[slot];
    ask.append('Use ', name, ` in TU slot ${slot + 1}?`);
    if (taken) ask.append(` This replaces ${taken}.`);
    else if (slot >= state.transcription_units.length) ask.append(' A new slot is added for it.');
  }
  bar.append(ask);

  const cancel = document.createElement('button');
  cancel.type = 'button';
  cancel.className = 'btn';
  cancel.textContent = 'Cancel';
  cancel.addEventListener('click', () => { finder.chosen = null; renderFinder(); });

  const go = document.createElement('button');
  go.type = 'button';
  go.className = 'btn btn-primary';
  go.textContent = 'Use it';
  go.addEventListener('click', () => applyHit(hit));

  bar.append(cancel, go);
  return bar;
}

function chooseHit(index) {
  finder.active = index;
  finder.chosen = finder.hits[index];
  renderFinder();
  const go = el('finder-results').querySelector('.finder-confirm .btn-primary');
  if (go) go.focus();
}

function applyHit(hit) {
  if (roleOf(hit) === 'backbone') {
    state.backbone = hit.name;
    el('backbone').value = hit.name;
  } else {
    const slot = slotFor(hit.name);
    while (state.transcription_units.length <= slot) state.transcription_units.push('');
    state.transcription_units[slot] = hit.name;
  }
  closeFinder();
  refresh();
}

function renderFinder() {
  const box = el('finder-results');
  box.replaceChildren();
  const query = el('finder-input').value.trim();
  if (!query) { box.hidden = true; return; }
  box.hidden = false;

  /* Only finished cassettes and destination vectors are offered here; a part
   * plasmid belongs on the Cassette page. Withheld matches are counted rather
   * than silently dropped. */
  const withheld = finder.suppressed
    ? `${finder.suppressed} other match${finder.suppressed === 1 ? '' : 'es'} `
      + 'are part plasmids, not finished cassettes — build one on the Cassette page first.'
    : '';

  if (!finder.hits.length) {
    const hint = document.createElement('div');
    hint.className = 'finder-hint';
    hint.textContent = finder.suppressed
      ? `Nothing annotated “${query}” can be used here. ${withheld}`
      : `Nothing in the library is annotated “${query}”.`;
    box.append(hint);
    return;
  }
  const terms = query.split(/\s+/).filter(Boolean);
  finder.hits.forEach((hit, i) => box.append(hitRow(hit, i, terms)));
  if (withheld) {
    const note = document.createElement('div');
    note.className = 'finder-hint';
    note.textContent = withheld;
    box.append(note);
  }
  if (finder.chosen) box.append(confirmBar(finder.chosen));
}

function closeFinder() {
  finder.hits = [];
  finder.active = -1;
  finder.chosen = null;
  el('finder-input').value = '';
  el('finder-results').hidden = true;
  el('finder-results').replaceChildren();
}

async function runSearch() {
  const query = el('finder-input').value.trim();
  finder.chosen = null;
  finder.active = -1;
  if (query.length < 2) {
    finder.hits = [];
    renderFinder();
    return;
  }
  const seq = ++finder.seq;
  const response = await fetch(
    `/api/library/search?q=${encodeURIComponent(query)}&limit=20&usable=multigene`,
  );
  if (seq !== finder.seq) return;
  const result = await response.json();
  finder.hits = result.hits || [];
  finder.suppressed = result.suppressed || 0;
  renderFinder();
}

function wireFinder() {
  const input = el('finder-input');
  let timer = null;
  input.addEventListener('input', () => {
    clearTimeout(timer);
    timer = setTimeout(runSearch, 160);
  });
  input.addEventListener('keydown', (event) => {
    if (event.key === 'Escape') { closeFinder(); return; }
    if (!finder.hits.length) return;
    if (event.key === 'ArrowDown' || event.key === 'ArrowUp') {
      event.preventDefault();
      const step = event.key === 'ArrowDown' ? 1 : -1;
      let next = finder.active;
      for (let i = 0; i < finder.hits.length; i += 1) {
        next = (next + step + finder.hits.length) % finder.hits.length;
        if (roleOf(finder.hits[next])) break;
      }
      finder.active = finder.hits.some(roleOf) ? next : -1;
      renderFinder();
    } else if (event.key === 'Enter' && finder.active >= 0) {
      event.preventDefault();
      chooseHit(finder.active);
    }
  });
  document.addEventListener('click', (event) => {
    if (!event.target.closest('.finder')) el('finder-results').hidden = true;
  });
  input.addEventListener('focus', () => {
    if (finder.hits.length) el('finder-results').hidden = false;
  });
}

async function load() {
  const [summary, options] = await Promise.all([
    fetch('/api/library/summary').then((r) => r.json()),
    fetch('/api/level3/options').then((r) => r.json()),
  ]);
  catalogue = options;
  el('folder').textContent = summary.roots.length === 1
    ? summary.roots[0]
    : `${summary.roots.length} folders`;

  const backbones = el('backbone');
  backbones.replaceChildren();
  const blank = document.createElement('option');
  blank.value = '';
  blank.textContent = options.backbones.length
    ? 'Choose a backbone…'
    : 'no destination vectors in the library';
  backbones.append(blank);
  for (const vector of options.backbones) {
    const option = document.createElement('option');
    option.value = vector.name;
    option.textContent = `${vector.name} · ${vector.left_overhang}→${vector.right_overhang} · ${vector.length} bp`;
    backbones.append(option);
  }

  // Open on a real assembly rather than on two red errors. The server picks a
  // backbone whose connector chain actually closes, so the first thing the
  // screen shows is a finished construct.
  try {
    const seed = await (await fetch('/api/level3/default')).json();
    if (seed.backbone) {
      state.backbone = seed.backbone;
      backbones.value = seed.backbone;
      state.transcription_units = seed.transcription_units.length
        ? seed.transcription_units
        : [''];
      if (seed.name) {
        state.name = seed.name;
        el('name').value = seed.name;
      }
    }
  } catch {
    // an empty screen is a worse start than a seeded one, not a broken one
  }

  refresh();
}

el('backbone').addEventListener('change', refresh);
el('name').addEventListener('change', refresh);
wireFinder();
el('add-tu').addEventListener('click', () => {
  state.transcription_units.push('');
  refresh();
});

for (const [id, path] of [['export-btn', '/api/level3/export']]) {
  el(id).addEventListener('click', async () => {
    const response = await fetch(path, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(state),
    });
    if (!response.ok) return;
    const url = URL.createObjectURL(await response.blob());
    const link = document.createElement('a');
    link.href = url;
    link.download = `${state.name}.gb`;
    link.click();
    URL.revokeObjectURL(url);
  });
}

el('save-btn').addEventListener('click', async () => {
  // The name field keeps its default until someone changes it, so the save
  // asks for one: this is the moment the name is being decided, and there is
  // no way in the app to remove a plasmid saved under the wrong one.
  const chosen = prompt('Save to library as:', state.name || 'multigene');
  if (!chosen || !chosen.trim()) return;
  state.name = chosen.trim();
  el('name').value = state.name;

  const button = el('save-btn');
  button.disabled = true;
  const result = await post('/api/level3/save', state);
  button.textContent = result.ok ? 'Saved' : 'Could not save';
  setTimeout(() => { button.textContent = 'Save to library'; button.disabled = false; }, 2000);
});

load();

// ---------------------------------------------------------------- protocol ---

/* What to pipette for the multigene reaction.
 *
 * Every multigene vector is a dropout whose enzyme sites leave with the
 * fragment being replaced, so this reaction always takes the short programme:
 * cycling, then hold. No final digest, no heat inactivation, and the colony
 * colour cannot be trusted on its own - which the server says in its notes.
 */

let protocolText = '';

function rxRow(component) {
  const tr = document.createElement('tr');
  if (component.kind === 'reagent' || component.kind === 'water') {
    tr.classList.add('is-reagent');
  }
  if (!component.measured && component.kind !== 'reagent' && component.kind !== 'water') {
    tr.classList.add('is-missing');
  }
  const cells = [
    [component.name, ''],
    [component.volume_ul === null ? '' : component.volume_ul.toFixed(2), 'num'],
    [component.ng === null ? '' : component.ng.toFixed(1), 'num'],
    [component.length ? component.length.toLocaleString() : '', 'num'],
    [component.note, 'rx-note'],
  ];
  for (const [text, cls] of cells) {
    const td = document.createElement('td');
    if (cls) td.className = cls;
    td.textContent = text;
    tr.append(td);
  }
  return tr;
}

function renderProtocol(rx) {
  const body = el('protocol-body');
  body.replaceChildren();

  for (const [kind, messages] of [['rx-warn', rx.issues || []],
                                  ['rx-lead', rx.notes || []]]) {
    for (const message of messages) {
      const box = document.createElement('div');
      box.className = kind;
      box.textContent = message;
      body.append(box);
    }
  }

  const table = document.createElement('table');
  table.className = 'rx-table';
  const head = document.createElement('tr');
  for (const [label, cls] of [['Component', ''], ['µL', 'num'], ['ng', 'num'],
                              ['bp', 'num'], ['', '']]) {
    const th = document.createElement('th');
    th.textContent = label;
    if (cls) th.className = cls;
    head.append(th);
  }
  table.append(head, ...rx.components.map(rxRow));
  body.append(table);

  const steps = document.createElement('ul');
  steps.className = 'rx-steps';
  for (const step of rx.steps) {
    const li = document.createElement('li');
    const b = document.createElement('b');
    b.textContent = step.label;
    const detail = document.createElement('span');
    detail.textContent = step.detail;
    li.append(b, detail);
    steps.append(li);
  }
  body.append(steps);
}

function wireProtocol() {
  el('protocol-btn').addEventListener('click', async () => {
    const rx = await post('/api/level3/protocol', state);
    protocolText = rx.text || '';
    renderProtocol(rx);
    el('protocol-title').textContent = `${rx.name} · ${rx.enzyme} reaction`;
    el('protocol-dialog').showModal();
  });
  el('protocol-copy').addEventListener('click', async () => {
    if (!protocolText) return;
    await navigator.clipboard.writeText(protocolText);
    const button = el('protocol-copy');
    button.textContent = 'Copied';
    setTimeout(() => { button.textContent = 'Copy'; }, 1500);
  });
}

wireProtocol();

// ----------------------------------------------------------------- design ---

/* The design direction.
 *
 * Every other screen starts from parts you already have. A project starts from
 * the other end: *I want these genes expressed, in this order*. So this takes
 * promoter-CDS-terminator per unit and works backwards to the eight-part
 * plasmids that have to exist first, assigning the connectors so the cassettes
 * chain - which is the fiddly part, and the part worth automating.
 */

let designRows = [];
let designed = null;

/* A dropdown you can search, for the design dialog.
 *
 * A plain <select> over a few hundred parts is a scroll, and the thing you
 * want is rarely near the name you would guess: parts are named for where they
 * came from, and what you remember is what is in them. So the box matches the
 * name, the component the digest found inside, and any alias - the same three
 * the Cassette screen's picker matches on.
 *
 * Filtering happens here rather than over /api/library/search because the
 * whole list for a position is already loaded, and a round trip per keystroke
 * would be slower and could reorder the list under the cursor.
 */
function haystackOf(entry) {
  return [entry.name, entry.display, entry.component, ...(entry.aliases || [])]
    .filter(Boolean).join(' ').toLowerCase();
}

function searchableSelect({ entries, chosen, blank, onPick }) {
  const wrap = document.createElement('div');
  wrap.className = 'design-pick';

  const find = document.createElement('input');
  find.type = 'search';
  find.className = 'design-find';
  find.placeholder = 'search\u2026';
  find.setAttribute('aria-label', 'Search parts');

  const select = document.createElement('select');
  let current = chosen || '';

  function fill() {
    const terms = find.value.toLowerCase().split(/\s+/).filter(Boolean);
    select.replaceChildren();
    const none = document.createElement('option');
    none.value = '';
    none.textContent = blank;
    select.append(none);

    let shown = 0;
    for (const entry of entries) {
      const matches = !terms.length
        || terms.every((t) => haystackOf(entry).includes(t));
      // Whatever is currently chosen stays in the list even when it does not
      // match, so typing a filter never silently unpicks it.
      if (!matches && entry.name !== current) continue;
      const option = document.createElement('option');
      option.value = entry.name;
      option.textContent = entry.display || entry.name;
      select.append(option);
      if (matches) shown += 1;
    }
    select.value = current;
    find.dataset.found = terms.length ? String(shown) : '';
    wrap.dataset.empty = terms.length && !shown ? 'true' : 'false';
  }

  find.addEventListener('input', fill);
  select.addEventListener('change', () => {
    current = select.value;
    onPick(current);
  });

  fill();
  wrap.append(find, select);
  return wrap;
}

function partOptions(select, entries, chosen) {
  const blank = document.createElement('option');
  blank.value = '';
  blank.textContent = '—';
  select.append(blank);
  for (const entry of entries) {
    const option = document.createElement('option');
    option.value = entry.name;
    option.textContent = entry.display || entry.name;
    if (entry.name === chosen) option.selected = true;
    select.append(option);
  }
}

function designRow(row, index) {
  const wrap = document.createElement('div');
  wrap.className = 'design-unit';

  const head = document.createElement('div');
  head.className = 'design-unit-head';
  const label = document.createElement('span');
  label.className = 'lbl';
  label.textContent = `Unit ${index + 1}`;
  const nameInput = document.createElement('input');
  nameInput.className = 'design-unit-name';
  nameInput.placeholder = 'name (optional)';
  nameInput.value = row.name || '';
  nameInput.addEventListener('input', () => { row.name = nameInput.value; });
  const remove = document.createElement('button');
  remove.type = 'button';
  remove.className = 'btn';
  remove.textContent = 'Remove';
  remove.addEventListener('click', () => {
    designRows.splice(index, 1);
    renderDesignUnits();
  });
  head.append(label, nameInput, remove);
  wrap.append(head);

  const grid = document.createElement('div');
  grid.className = 'design-grid';
  for (const [slot, field, caption] of [
    ['2', 'promoter', 'Promoter'],
    ['3', 'cds', 'Coding sequence'],
    ['4', 'terminator', 'Terminator'],
  ]) {
    const cell = document.createElement('label');
    cell.className = 'design-cell';
    const text = document.createElement('span');
    text.className = 'lbl';
    text.textContent = caption;
    cell.append(text, searchableSelect({
      entries: designParts[slot] || [],
      chosen: row[field],
      blank: '\u2014',
      onPick: (name) => { row[field] = name; },
    }));
    grid.append(cell);
  }

  /* The backbone, per unit. Blank means "whatever the design is using", which
   * is the common case: a project usually varies the marker alone, because
   * each unit integrates at a different locus, and keeps one origin and one
   * E. coli backbone throughout. */
  row.backbone = row.backbone || {};
  for (const [slot, caption] of BACKBONE_SLOTS) {
    const cell = document.createElement('label');
    cell.className = 'design-cell';
    const text = document.createElement('span');
    text.className = 'lbl';
    text.textContent = caption;
    cell.append(text, searchableSelect({
      entries: designParts[slot] || [],
      chosen: row.backbone[slot] || '',
      blank: 'same as the rest',
      onPick: (name) => {
        if (name) row.backbone[slot] = name;
        else delete row.backbone[slot];
      },
    }));
    grid.append(cell);
  }
  wrap.append(grid);
  return wrap;
}

let designParts = {};

function renderDesignUnits() {
  const host = el('design-units');
  host.replaceChildren(...designRows.map(designRow));
}

function renderDesignReport(result) {
  designed = result;
  const host = el('design-report');
  host.replaceChildren();

  for (const issue of result.issues || []) {
    const box = document.createElement('div');
    box.className = issue.level === 'error' ? 'rx-warn' : 'rx-lead';
    box.textContent = issue.message;
    host.append(box);
  }

  const pre = document.createElement('pre');
  pre.className = 'design-report';
  pre.textContent = result.report;
  host.append(pre);
  el('design-write').disabled = !result.ok;
  el('design-protocols').disabled = !result.ok;
  el('design-save').disabled = !result.ok;
}

/* The cassette backbone: positions 6, 7 and 8, chosen once for the whole
 * design and overridable on any one unit.
 *
 * Not positions 1 and 5. Those are the connectors, and they are what makes the
 * units chain in the order the rows are in - offering them here would mean
 * setting the order twice, in two places, with nothing keeping the two
 * answers the same. */
const BACKBONE_SLOTS = [
  ['6', 'Marker'],
  ['7', 'Yeast origin / homology'],
  ['8', 'E. coli backbone'],
];

/* What the shared pickers are set to. Kept here rather than read back off the
 * element, because each one is replaced by a searchable control - a wrapper
 * round an input and a select - and a wrapper has no `value`. */
const sharedChoice = {};

/* The shared pickers get the same searchable control as the units. They are
 * replaced rather than filled, because the markup ships a plain <select> so
 * the dialog still reads sensibly before any script runs - and replaced again
 * on a reload, since a custom control cannot be set by assigning `.value`. */
function wireSharedPickers() {
  for (const [slot] of BACKBONE_SLOTS) {
    const picker = el(`design-shared-${slot}`);
    if (!picker) continue;
    const replacement = searchableSelect({
      entries: designParts[slot] || [],
      chosen: sharedChoice[slot] || '',
      blank: 'Choose for me',
      onPick: (name) => { sharedChoice[slot] = name; },
    });
    replacement.id = picker.id;
    picker.replaceWith(replacement);
  }
}

function sharedBackbone() {
  const out = {};
  for (const [slot] of BACKBONE_SLOTS) {
    if (sharedChoice[slot]) out[slot] = sharedChoice[slot];
  }
  return out;
}

function designBody() {
  return {
    name: el('design-title').value.trim() || 'pPathway',
    units: designRows,
    backbone: el('design-backbone').value || null,
    shared: sharedBackbone(),
  };
}

function wireDesign() {
  el('design-btn').addEventListener('click', async () => {
    if (!Object.keys(designParts).length) {
      const rows = await (await fetch('/api/library/plasmids')).json();
      for (const slot of ['2', '3', '4', '6', '7', '8']) {
        designParts[slot] = rows
          .filter((r) => r.part_type === slot)
          .sort((a, b) => (a.display || a.name).localeCompare(b.display || b.name));
      }
      wireSharedPickers();
    }
    const picker = el('design-backbone');
    if (picker.options.length <= 1) {
      for (const vector of catalogue.backbones) {
        const option = document.createElement('option');
        option.value = vector.name;
        option.textContent =
          `${vector.name} · ${vector.left_overhang}→${vector.right_overhang}`;
        picker.append(option);
      }
    }
    if (!designRows.length) designRows = [{}, {}];
    renderDesignUnits();
    el('design-report').replaceChildren();
    el('design-write').disabled = true;
    el('design-protocols').disabled = true;
    el('design-save').disabled = true;
    el('design-dialog').showModal();
  });

  el('design-add').addEventListener('click', () => {
    designRows.push({});
    renderDesignUnits();
  });

  el('design-run').addEventListener('click', async () => {
    const button = el('design-run');
    button.disabled = true;
    try {
      renderDesignReport(await post('/api/level3/design', designBody()));
    } finally {
      button.disabled = false;
    }
  });

  /* Download the reactions: one BsaI per cassette, then the BsmBI that joins
   * them. Separate from the plasmids because they are what you take to the
   * bench, and you want them on separate days. */
  el('design-protocols').addEventListener('click', async () => {
    const response = await fetch('/api/level3/design/protocols.txt', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(designBody()),
    });
    if (!response.ok) return;
    const url = URL.createObjectURL(await response.blob());
    const link = document.createElement('a');
    link.href = url;
    link.download = `${designBody().name}-reactions.txt`;
    link.click();
    URL.revokeObjectURL(url);
  });

  /* Put the lot in the library. Until this existed the only way to keep a
   * design was a zip to unpack by hand, which meant a cassette you had just
   * designed could not be picked on any screen. */
  el('design-save').addEventListener('click', async () => {
    const button = el('design-save');
    const was = button.textContent;
    button.disabled = true;
    button.textContent = 'Saving\u2026';
    try {
      const body = await post('/api/level3/design/save', designBody());
      if (body.ok) {
        button.textContent = `Saved ${body.saved.length}`;
        await refresh();
        setTimeout(() => { button.textContent = was; button.disabled = false; }, 2500);
        return;
      }
      window.alert('Not saved:\n\n' + (body.issues || []).join('\n'));
    } catch (error) {
      window.alert(`Not saved: ${error.message || error}`);
    }
    button.textContent = was;
    button.disabled = false;
  });

  /* Reload a design from the -design.json the download carries. The build
   * order is written for a human and does not round-trip: names are
   * abbreviated in it, and a part chosen and then changed leaves no trace. */
  el('design-load').addEventListener('change', async (event) => {
    const file = event.target.files && event.target.files[0];
    if (!file) return;
    try {
      const saved = JSON.parse(await file.text());
      if (!Array.isArray(saved.units) || !saved.units.length) {
        throw new Error('no transcription units in that file');
      }
      el('design-title').value = saved.name || 'pPathway';
      designRows = saved.units.map((u) => ({ ...u, backbone: { ...(u.backbone || {}) } }));
      for (const [slot] of BACKBONE_SLOTS) {
        sharedChoice[slot] = (saved.shared || {})[slot] || '';
      }
      el('design-backbone').value = saved.backbone || '';
      wireSharedPickers();
      renderDesignUnits();
      el('design-report').replaceChildren();
      el('design-write').disabled = true;
      el('design-protocols').disabled = true;
      el('design-save').disabled = true;
    } catch (error) {
      window.alert(
        `That is not a design file.\n\n${error.message || error}\n\n`
        + 'Use the -design.json from a Download all, not the build order.');
    } finally {
      event.target.value = '';
    }
  });

  el('design-write').addEventListener('click', async () => {
    const response = await fetch('/api/level3/design.zip', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(designBody()),
    });
    if (!response.ok) return;
    const url = URL.createObjectURL(await response.blob());
    const link = document.createElement('a');
    link.href = url;
    link.download = `${designBody().name}-design.zip`;
    link.click();
    URL.revokeObjectURL(url);
  });
}

wireDesign();
