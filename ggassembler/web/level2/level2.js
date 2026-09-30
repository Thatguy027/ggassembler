/* Level 2 - cassette assembly.
 *
 * This module owns the Level 2 screen's state and talks only to
 * /api/level2/* and the shared read-only /api/library/*. It imports nothing
 * from a sibling level, and nothing imports it.
 *
 * The state is exactly the design the API expects - view-mode flags plus one
 * chosen plasmid per slot - so every render is: post the design, draw what
 * comes back. The server decides what fits where; the browser never reasons
 * about overhangs on its own.
 */

import {
  RING, layoutCallouts, pointAt, strandOf,
} from './layout.js';

const state = {
  selections: {},
  split_3: false,
  split_4: false,
  split_8: false,
  composite_left: false,
  composite_right: false,
  name: 'pCassette',
  showDirection: true,
  showLinear: true,
};

const VIEW_KEY = 'ggasm.map-view';

/* The two map toggles, remembered between visits. Browser storage can throw
 * outright in a private window, so a failure here leaves the defaults standing
 * rather than taking the screen down with it. */
function loadView() {
  try {
    const saved = JSON.parse(localStorage.getItem(VIEW_KEY) || '{}');
    if (typeof saved.showDirection === 'boolean') state.showDirection = saved.showDirection;
    if (typeof saved.showLinear === 'boolean') state.showLinear = saved.showLinear;
  } catch {
    // defaults stand
  }
}

function saveView() {
  try {
    localStorage.setItem(VIEW_KEY, JSON.stringify({
      showDirection: state.showDirection, showLinear: state.showLinear,
    }));
  } catch {
    // remembering is a convenience, not a requirement
  }
}

function wireMapTools() {
  for (const [id, flag] of [['toggle-features', 'showDirection'],
                            ['toggle-linear', 'showLinear']]) {
    const button = el(id);
    button.setAttribute('aria-pressed', String(state[flag]));
    button.addEventListener('click', () => {
      state[flag] = !state[flag];
      button.setAttribute('aria-pressed', String(state[flag]));
      saveView();
      refresh();
    });
  }
}

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

// --------------------------------------------------------------- panels ---

/* The conic-gradient is the one place a part colour has to become a literal
 * string: a gradient stop cannot be a var() that changes underneath it, because
 * the whole background is one computed value. So it is read from the document
 * at paint time rather than baked in, and the ring is repainted on a theme
 * change - see `repaint` below. */
function partColor(partType) {
  const position = String(partType || '').charAt(0);
  const token = /[1-8]/.test(position) ? `--part-${position}` : '--unknown-bg';
  const value = getComputedStyle(document.documentElement).getPropertyValue(token);
  return value.trim() || 'var(--unknown-bg)';
}

function badge(slot) {
  const span = document.createElement('span');
  span.className = 'badge';
  span.textContent = slot.label;
  // the type, not a colour: tokens.css turns it into one, and a theme switch
  // re-resolves it without anything here running again
  span.dataset.part = slot.key || '';
  return span;
}

function segGroup(buttons) {
  const group = document.createElement('div');
  group.className = 'seg-group';
  group.setAttribute('role', 'group');
  for (const [label, pressed, onClick] of buttons) {
    const button = document.createElement('button');
    button.type = 'button';
    button.className = 'seg';
    button.textContent = label;
    button.setAttribute('aria-pressed', String(pressed));
    button.addEventListener('click', onClick);
    group.append(button);
  }
  return group;
}

/** The segmented control a panel carries, if its position can be split. */
function splitControl(slot) {
  if (slot.key === '3' || slot.key === '3a' || slot.key === '3b') {
    return segGroup([
      ['3', !state.split_3, () => setSplit('split_3', false, ['3a', '3b'])],
      ['3a+3b', state.split_3, () => setSplit('split_3', true, ['3'])],
    ]);
  }
  if (slot.key === '4' || slot.key === '4a' || slot.key === '4b') {
    return segGroup([
      ['4', !state.split_4, () => setSplit('split_4', false, ['4a', '4b'])],
      ['4a+4b', state.split_4, () => setSplit('split_4', true, ['4'])],
    ]);
  }
  if (slot.key === '8' || slot.key === '8a' || slot.key === '8b') {
    return segGroup([
      ['8', !state.split_8, () => setSplit('split_8', false, ['8a', '8b'])],
      ['8a+8b', state.split_8, () => setSplit('split_8', true, ['8'])],
    ]);
  }
  return null;
}

function setSplit(flag, value, clearKeys) {
  if (state[flag] === value) return;
  state[flag] = value;
  for (const key of clearKeys) delete state.selections[key];
  refresh();
}

// ------------------------------------------------------------------ picker ---

/* The part picker.
 *
 * A native <select> stops working somewhere around forty options, and this
 * library gives the type 3 slot ninety-one of them in ASCII order, so
 * `ACAT1` and `C10_mRuby2` come before every pYTK part. Worse, the detail
 * that tells them apart - the component, the aliases - lived in `title`,
 * which a mouse reveals only on a hover and a keyboard never reveals at all.
 *
 * So: a text box that filters, over four fields at once, with the detail on
 * the row where it can be read. The <select> stays in the DOM, hidden and in
 * step, so the panel still submits and still works without this script.
 *
 * Matching runs in the browser against data the slot already carries - no
 * round trip per keystroke - which is why `labels` is in the slot payload.
 */

const RECENT_KEY = 'ggasm.recent-parts';
const RECENT_MAX = 8;

function recentParts() {
  try {
    return JSON.parse(localStorage.getItem(RECENT_KEY) || '[]');
  } catch {
    return []; // a private window, or a cleared store: no history is fine
  }
}

function rememberPart(name) {
  try {
    const kept = [name, ...recentParts().filter((n) => n !== name)].slice(0, RECENT_MAX);
    localStorage.setItem(RECENT_KEY, JSON.stringify(kept));
  } catch {
    // remembering is a convenience; failing to is not an error
  }
}

/** Everything about an option worth typing at. */
function haystack(option) {
  return [
    option.name,
    option.display,
    option.component,
    ...(option.aliases || []),
    ...(option.labels || []),
  ].filter(Boolean).join(' \u0001 ').toLowerCase();
}

function matches(option, terms) {
  if (!terms.length) return true;
  const hay = haystack(option);
  return terms.every((term) => hay.includes(term));
}

/** Recently used first, then your own parts, then the YTK reference set. */
function groupOptions(options) {
  const recent = recentParts();
  const groups = [
    { label: 'Recently used', items: [] },
    { label: 'Your parts', items: [] },
    { label: 'YTK reference', items: [] },
  ];
  for (const option of options) {
    if (recent.includes(option.name)) groups[0].items.push(option);
    else if (/^pytk\d+/i.test(option.name)) groups[2].items.push(option);
    else groups[1].items.push(option);
  }
  groups[0].items.sort((a, b) => recent.indexOf(a.name) - recent.indexOf(b.name));
  return groups.filter((g) => g.items.length);
}

function optionRow(option, terms, chosen) {
  const item = document.createElement('button');
  item.type = 'button';
  item.className = 'pick-row';
  item.setAttribute('role', 'option');
  item.dataset.name = option.name;
  if (option.name === chosen) item.setAttribute('aria-selected', 'true');

  const what = document.createElement('span');
  what.className = 'pick-what';

  const name = document.createElement('span');
  name.className = 'pick-name';
  name.append(highlight(option.display || option.name, terms));
  what.append(name);

  // the detail that used to hide in `title`, where nobody could read it
  const detail = [];
  if (option.aliases && option.aliases.length) {
    detail.push(`also ${option.aliases.join(', ')}`);
  }
  if (option.labels && option.labels.length) detail.push(option.labels.join(' · '));
  if (detail.length) {
    const sub = document.createElement('span');
    sub.className = 'pick-sub';
    sub.append(highlight(detail.join(' — '), terms));
    if (option.component_source === 'filename') sub.classList.add('is-guessed');
    what.append(sub);
  }
  item.append(what);

  if (option.internal_sites) {
    const warn = document.createElement('span');
    warn.className = 'pick-warn';
    warn.textContent = `${option.internal_sites} internal site${
      option.internal_sites === 1 ? '' : 's'}`;
    warn.title = 'this part still carries a site for the assembly enzyme';
    item.append(warn);
  }

  // A connector whose multigene site was domesticated away works perfectly
  // here and makes a cassette that can never be cut out again. Said on the row
  // rather than after the build, which is the only point it can still change
  // the choice.
  if (option.level3_ready === false) {
    const blocks = document.createElement('span');
    blocks.className = 'pick-warn is-blocking';
    blocks.textContent = 'blocks multigene';
    blocks.title =
      'usable here, but this part carries no site for the multigene enzyme, '
      + 'so the finished cassette could not be released for a Level 3 assembly';
    item.append(blocks);
  }

  const len = document.createElement('span');
  len.className = 'pick-len';
  len.textContent = `${option.length} bp`;
  item.append(len);
  return item;
}

function picker(slot) {
  const wrap = document.createElement('div');
  wrap.className = 'pick';

  // kept in step and hidden: the panel still works without this script
  const select = document.createElement('select');
  select.id = `slot-${slot.key}`;
  select.className = 'pick-fallback';
  select.tabIndex = -1;
  const blank = document.createElement('option');
  blank.value = '';
  blank.textContent = slot.match_count ? 'Choose a part…' : 'No part in the library fits';
  select.append(blank);
  for (const option of slot.options) {
    const item = document.createElement('option');
    item.value = option.name;
    item.textContent = `${option.display || option.name} · ${option.length} bp`;
    if (option.name === slot.selected) item.selected = true;
    select.append(item);
  }
  select.addEventListener('change', () => choose(select.value));

  const chosen = slot.options.find((o) => o.name === slot.selected);
  const input = document.createElement('input');
  input.className = 'pick-input';
  input.type = 'text';
  input.autocomplete = 'off';
  input.spellcheck = false;
  input.setAttribute('role', 'combobox');
  input.setAttribute('aria-expanded', 'false');
  input.setAttribute('aria-controls', `pick-list-${slot.key}`);
  input.placeholder = slot.match_count
    ? `Search ${slot.match_count} part${slot.match_count === 1 ? '' : 's'}…`
    : 'No part in the library fits';
  input.value = chosen ? chosen.display || chosen.name : '';
  input.disabled = slot.match_count === 0;
  input.title = chosen ? `${chosen.name} · ${chosen.length} bp` : '';

  const list = document.createElement('div');
  list.className = 'pick-list';
  list.id = `pick-list-${slot.key}`;
  list.setAttribute('role', 'listbox');
  list.hidden = true;

  let active = -1;
  let shown = [];

  function close(revert = true) {
    list.hidden = true;
    input.setAttribute('aria-expanded', 'false');
    active = -1;
    if (revert) {
      const still = slot.options.find((o) => o.name === state.selections[slot.key]);
      input.value = still ? still.display || still.name : '';
    }
  }

  function choose(name) {
    if (name) {
      state.selections[slot.key] = name;
      rememberPart(name);
    } else {
      delete state.selections[slot.key];
    }
    close(false);
    refresh();
  }

  function render() {
    const terms = input.value.trim().toLowerCase().split(/\s+/).filter(Boolean);
    // a box still showing the current selection is not a filter
    const typing = !(chosen && input.value === (chosen.display || chosen.name));
    const filtered = typing ? slot.options.filter((o) => matches(o, terms)) : slot.options;
    const marks = typing ? terms : [];

    list.replaceChildren();
    shown = [];
    if (!filtered.length) {
      const hint = document.createElement('div');
      hint.className = 'pick-hint';
      hint.textContent = `No part here matches “${input.value.trim()}”.`;
      list.append(hint);
      return;
    }

    for (const group of groupOptions(filtered)) {
      const head = document.createElement('div');
      head.className = 'pick-group';
      head.textContent = group.label;
      list.append(head);
      for (const option of group.items) {
        const item = optionRow(option, marks, state.selections[slot.key]);
        item.addEventListener('mousedown', (event) => {
          event.preventDefault(); // keep focus so blur does not revert first
          choose(option.name);
        });
        list.append(item);
        shown.push(item);
      }
    }
    if (state.selections[slot.key]) {
      const clear = document.createElement('button');
      clear.type = 'button';
      clear.className = 'pick-clear';
      clear.textContent = 'Clear this slot';
      clear.addEventListener('mousedown', (event) => {
        event.preventDefault();
        choose('');
      });
      list.append(clear);
    }
  }

  function open() {
    if (input.disabled) return;
    render();
    list.hidden = false;
    input.setAttribute('aria-expanded', 'true');
  }

  function move(step) {
    if (list.hidden) { open(); return; }
    if (!shown.length) return;
    active = (active + step + shown.length) % shown.length;
    for (const [i, item] of shown.entries()) {
      item.classList.toggle('is-active', i === active);
    }
    shown[active].scrollIntoView({ block: 'nearest' });
  }

  input.addEventListener('focus', () => { input.select(); open(); });
  input.addEventListener('input', () => { active = -1; open(); });
  input.addEventListener('blur', () => close());
  input.addEventListener('keydown', (event) => {
    if (event.key === 'ArrowDown') { event.preventDefault(); move(1); }
    else if (event.key === 'ArrowUp') { event.preventDefault(); move(-1); }
    else if (event.key === 'Escape') { event.preventDefault(); close(); }
    else if (event.key === 'Enter') {
      event.preventDefault();
      if (active >= 0 && shown[active]) choose(shown[active].dataset.name);
      else if (shown.length === 1) choose(shown[0].dataset.name);
    }
  });

  wrap.append(input, list, select);
  return wrap;
}

function panel(slot) {
  const section = document.createElement('section');
  section.className = 'card panel';
  section.dataset.slot = slot.key;
  if (!slot.selected) section.classList.add('is-empty');

  const head = document.createElement('div');
  head.className = 'card-head';
  head.append(badge(slot));

  /* What a slot *is* is its pair of overhangs; what usually goes there is a
   * convention. Leading with "5' assembly connector" while the list offers
   * LexA_4 [LexA site] invites picking a binding site in the belief it is a
   * connector, so the position leads and the description follows it. */
  const title = document.createElement('span');
  title.className = 'card-title';
  title.textContent = `Type ${slot.key}`;
  head.append(title);

  const ends = document.createElement('span');
  ends.className = 'card-ends oh';
  ends.textContent = `${slot.five_prime}→${slot.three_prime}`;
  head.append(ends);

  const spacer = document.createElement('span');
  spacer.style.flexGrow = '1';
  head.append(spacer);

  const control = splitControl(slot);
  if (control) head.append(control);
  section.append(head);

  const label = document.createElement('label');
  label.className = 'vh';
  label.setAttribute('for', `slot-${slot.key}`);
  label.textContent = `Type ${slot.key} part`;
  section.append(label);

  const row = document.createElement('div');
  row.className = 'panel-overhangs';

  row.append(picker(slot));
  section.append(row);

  const note = document.createElement('span');
  note.className = 'card-note';
  note.textContent =
    `${slot.description} · ${slot.match_count} part${
      slot.match_count === 1 ? '' : 's'} fit`;
  note.title = 'what usually goes here; the overhangs above are what actually fits';
  section.append(note);

  return section;
}

// ------------------------------------------------------------------ ring ---

/* Where the first part starts on the ring, as CSS measures it: clockwise from
 * 12 o'clock. The design canvas begins at 9 o'clock, so this is -90.
 *
 * The painted ring and the invisible wedges that catch the clicks have to agree
 * on this, and they use different conventions to get there: CSS conic-gradient
 * counts clockwise from the top, while SVG's cos/sin count from 3 o'clock. Both
 * are derived from this one constant so they cannot drift apart again.
 */
const RING_START_DEG = -90;
const RING_START_RAD = ((RING_START_DEG - 90) * Math.PI) / 180;

/* The map.
 *
 * The part band, the centre disc, and one callout per part round the outside.
 * A second band drawn from each part's annotations was tried and taken out:
 * most annotations in a real library span their whole part, so it read as a
 * paler copy of the band outside it rather than as a view into one. What the
 * feature data is still worth is direction - the chevrons below.
 *
 * Callouts rather than labels on the arcs, because of parts like a 250 bp
 * connector in a 6 kb cassette: fourteen degrees of arc holds no text at any
 * ring size this page can afford. The geometry and the collision pass live in
 * layout.js, which has no DOM in it and is tested directly.
 */

const SVG_NS = 'http://www.w3.org/2000/svg';

function svgEl(tag, attrs = {}) {
  const node = document.createElementNS(SVG_NS, tag);
  for (const [k, v] of Object.entries(attrs)) node.setAttribute(k, v);
  return node;
}

/** Paint both bands, the clickable wedges, the chevrons and the callouts. */
/* The ring is the one thing on this screen with a colour baked into a computed
 * value rather than a var(), so it is the one thing a theme switch cannot fix
 * on its own. Everything else re-resolves through `data-part`. */
let lastDrawn = null;

document.addEventListener('themechange', () => {
  if (lastDrawn) drawRing(lastDrawn.parts, lastDrawn.length);
});

function drawRing(parts, total) {
  const band = el('ring');
  const layer = el('map-layer');

  if (!parts.length || !total) {
    band.style.background = 'var(--chip)';
    layer.replaceChildren();
    el('callouts').replaceChildren();
    return;
  }

  // --- part band: unchanged, hairline gaps and all
  const stops = [];
  const hairline = 1.6;
  for (const part of parts) {
    const from = (part.start / total) * 360;
    const to = (part.end / total) * 360;
    stops.push(`${partColor(part.part_type)} ${from}deg ${Math.max(from, to - hairline)}deg`);
    stops.push(`var(--surface) ${Math.max(from, to - hairline)}deg ${to}deg`);
  }
  band.style.background = `conic-gradient(from ${RING_START_DEG}deg, ${stops.join(', ')})`;

  layer.replaceChildren();

  // the invisible wedges that catch clicks and keyboard focus
  for (const part of parts) {
    const path = svgEl('path', { d: wedge(part.start / total, part.end / total) });
    path.classList.add('arc-hit');
    path.setAttribute('tabindex', '0');
    path.setAttribute('role', 'button');
    path.setAttribute('aria-label', `${part.part_type}: ${part.label || part.source_name}`);
    const title = svgEl('title');
    title.textContent = [
      `Type ${part.part_type} · ${part.label || part.source_name}`,
      `${part.length.toLocaleString()} bp`,
      `joins ${part.left_overhang} → ${part.right_overhang}`,
    ].join('\n');
    path.append(title);
    const focus = () => focusPanel(part.part_type);
    path.addEventListener('click', focus);
    path.addEventListener('keydown', (event) => {
      if (event.key === 'Enter' || event.key === ' ') { event.preventDefault(); focus(); }
    });
    layer.append(path);
  }

  const placed = layoutCallouts(parts, total, RING_START_DEG);

  if (state.showDirection) {
    for (const item of placed.items) {
      const strand = strandOf(item.part);
      if (!strand) continue;
      const at = pointAt(RING.chevron, item.theta);
      layer.append(svgEl('path', {
        d: 'M -5 4 L 0 -4 L 5 4 Z',
        class: 'ring-chevron',
        transform: `translate(${at.x} ${at.y}) rotate(${strand > 0 ? item.theta : item.theta + 180})`,
      }));
    }
  }

  for (const item of placed.items) {
    layer.append(svgEl('line', {
      class: 'callout-stub',
      x1: item.stubStart.x, y1: item.stubStart.y,
      x2: item.stubEnd.x, y2: item.stubEnd.y,
    }));
    layer.append(svgEl('circle', {
      class: 'callout-dot', cx: item.dot.x, cy: item.dot.y, r: 2,
    }));
    // only where the label moved far enough that the stub end would otherwise
    // fall outside its own box
    if (item.needsElbow) {
      layer.append(svgEl('line', {
        class: 'callout-elbow',
        x1: item.stubEnd.x, y1: item.stubEnd.y, x2: item.stubEnd.x, y2: item.y,
      }));
    }
  }

  renderCallouts(placed);
}

/** One label per part, placed where the collision pass put it. */
function renderCallouts(placed) {
  const host = el('callouts');
  host.replaceChildren();

  for (const item of placed.items) {
    const part = item.part;
    const box = document.createElement('button');
    box.type = 'button';
    box.className = `callout is-${item.side}`;
    // right-hand labels grow rightwards from the stub; left-hand ones grow
    // leftwards *to* it, which is `right`, not `left` plus a transform. A
    // transform moves the box after layout, so shrink-to-fit never knows it
    // has less room and the label walks off the left of the container.
    const across = (item.labelX / RING.width) * 100;
    if (item.side === 'left') box.style.right = `${100 - across}%`;
    else box.style.left = `${across}%`;
    box.style.top = `${(item.y / RING.height) * 100}%`;
    box.style.height = `${placed.labelHeight}px`;
    box.addEventListener('click', () => focusPanel(part.part_type));

    const badge = document.createElement('span');
    badge.className = 'callout-badge';
    badge.textContent = part.part_type || '?';
    badge.dataset.part = part.part_type || '';
    badge.style.background = 'var(--part-color)';

    // `label` is the component wherever there is one, so naming the callout by
    // it and then repeating the component below printed the same string twice -
    // once whole, once truncated. The plasmid is what you picked, what you
    // order by, and what the panel on the left shows, so that is the name.
    const name = document.createElement('span');
    name.className = 'callout-name';
    name.textContent = part.source_name;

    const meta = document.createElement('span');
    meta.className = 'callout-meta';
    const strand = strandOf(part);
    meta.textContent = `${part.length.toLocaleString()} bp`
      + (part.component ? ` · ${part.component}` : '')
      + (strand ? (strand > 0 ? ' →' : ' ←') : '');

    const text = document.createElement('span');
    text.className = 'callout-text';
    text.append(name, meta);
    box.append(badge, text);
    host.append(box);
  }
}

/** The clickable wedge for one arc, in the map's own coordinates. */
function wedge(fromFraction, toFraction) {
  const a0 = fromFraction * 2 * Math.PI + RING_START_RAD;
  const a1 = toFraction * 2 * Math.PI + RING_START_RAD;
  const large = a1 - a0 > Math.PI ? 1 : 0;
  const p = (radius, angle) =>
    `${RING.cx + radius * Math.cos(angle)} ${RING.cy + radius * Math.sin(angle)}`;
  return [
    `M ${p(RING.outer, a0)}`,
    `A ${RING.outer} ${RING.outer} 0 ${large} 1 ${p(RING.outer, a1)}`,
    `L ${p(RING.bandInner, a1)}`,
    `A ${RING.bandInner} ${RING.bandInner} 0 ${large} 0 ${p(RING.bandInner, a0)}`,
    'Z',
  ].join(' ');
}

function focusPanel(partType) {
  // the combobox input, not the fallback <select> beside it - that one is
  // display:none, and focusing it silently does nothing
  const panelEl = document.querySelector(`.panel[data-slot="${partType}"] .pick-input`);
  if (panelEl) {
    panelEl.focus();
    panelEl.closest('.panel').scrollIntoView({ block: 'nearest', behavior: 'smooth' });
  }
}

// --------------------------------------------------------------- rendering ---

function renderPanels(slots) {
  const left = slots.filter((s) => s.column === 'left');
  const right = slots.filter((s) => s.column === 'right');
  el('panels-left').replaceChildren(...left.map(panel));
  el('panels-right').replaceChildren(...right.map(panel));

  for (const button of document.querySelectorAll('[data-view]')) {
    const view = button.dataset.view;
    const on =
      (view === 'composite-left' && state.composite_left) ||
      (view === 'parts-left' && !state.composite_left) ||
      (view === 'composite-right' && state.composite_right) ||
      (view === 'parts-right' && !state.composite_right);
    button.setAttribute('aria-pressed', String(on));
  }
}


function renderJunctions(junctions, issues) {
  const strip = el('junction-strip');
  if (!junctions.length) {
    strip.replaceChildren(note('pick a part in every slot to see the junctions'));
    return;
  }
  const risky = new Set();
  for (const issue of issues) {
    if (issue.code.startsWith('near_match') || issue.code === 'duplicate_overhang') {
      for (const match of issue.message.matchAll(/\b([ACGT]{4})\b/g)) risky.add(match[1]);
    }
  }
  const chips = junctions.map((junction) => {
    const chip = document.createElement('span');
    chip.className = risky.has(junction.overhang) ? 'oh oh-bad' : 'oh oh-ok';
    chip.textContent = junction.overhang;
    chip.title = `${junction.upstream} → ${junction.downstream}`;
    return chip;
  });
  chips.push(
    note(risky.size ? 'some overhangs can cross-match' : 'no repeats, no 3/4-nt cross-matches'),
  );
  strip.replaceChildren(...chips);
}

function note(text) {
  const span = document.createElement('span');
  span.className = 'card-note';
  span.style.alignSelf = 'center';
  span.textContent = text;
  return span;
}

function renderIssues(issues) {
  const box = el('issues');
  const shown = issues.filter((i) => i.code !== 'empty_slot');
  box.hidden = shown.length === 0;
  box.replaceChildren(
    ...shown.map((issue) => {
      const row = document.createElement('div');
      row.className = `issue issue-${issue.level}`;
      const code = document.createElement('span');
      code.className = 'issue-code';
      code.textContent = issue.code;
      const message = document.createElement('span');
      message.textContent = issue.message;
      row.append(code, message);
      return row;
    }),
  );
}

function renderStage(result) {
  const empty = result.issues.filter((i) => i.code === 'empty_slot').length;
  const pill = el('ring-pill');

  el('ring-bp').textContent = result.ok ? `${result.length.toLocaleString()} bp` : '—';
  el('ring-sub').textContent = result.ok
    ? `circular · ${result.parts.length} parts`
    : empty
      ? `${empty} slot${empty === 1 ? '' : 's'} still empty`
      : 'cannot assemble yet';

  if (result.ok) {
    pill.hidden = false;
    pill.className = result.counts.warnings ? 'pill pill-warn' : 'pill';
    pill.textContent = result.counts.warnings
      ? `${result.counts.junctions} junctions · ${result.counts.warnings} warning${result.counts.warnings === 1 ? '' : 's'}`
      : `${result.counts.junctions} / ${result.counts.junctions} junctions compatible`;
  } else if (result.counts.errors && !empty) {
    pill.hidden = false;
    pill.className = 'pill pill-error';
    pill.textContent = `${result.counts.errors} problem${result.counts.errors === 1 ? '' : 's'}`;
  } else {
    pill.hidden = true;
  }

  const screening = result.issues.find((i) => i.code === 'screening');
  el('stage-subtitle').textContent = result.ok
    ? `Level 2 cassette · ${result.enzyme}${screening ? ' · ' + screening.message.split(': ')[1] : ''}`
    : ' ';

  el('rx-enzyme').textContent = result.enzyme || '—';
  el('rx-selection').textContent = screening
    ? screening.message.split('screen ')[1] || '—'
    : result.ok ? 'white colonies' : '—';

  el('export-btn').disabled = !result.ok;
  el('protocol-btn').disabled = !result.ok;
}

// ------------------------------------------------------------------ finder ---

/* Search the library by annotation and drop the result into the slot it fits.
 *
 * The search itself is shared (/api/library/search); what a hit *means* is not,
 * and this is the Level 2 answer: a hit is a part type, a part type is a panel,
 * and some panels only exist in one view. A type 3a part has nowhere to go
 * while the 3 panel is whole, so choosing one flips that panel to 3a+3b rather
 * than refusing. Nothing is filled without being asked first.
 */

const finder = { hits: [], active: -1, chosen: null, seq: 0, suppressed: 0 };

/** The view flag a part type needs turned on or off before its panel exists. */
const VIEW_FOR_TYPE = {
  3: ['split_3', false], '3a': ['split_3', true], '3b': ['split_3', true],
  4: ['split_4', false], '4a': ['split_4', true], '4b': ['split_4', true],
  8: ['split_8', false], '8a': ['split_8', true], '8b': ['split_8', true],
  234: ['composite_left', true], 678: ['composite_right', true],
};

/** Which column a type sits in, so the composite toggles can be cleared. */
const LEFT_TYPES = new Set(['2', '3', '3a', '3b', '4', '4a', '4b', '234']);

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
  const re = new RegExp(`(${pattern})`, 'ig');
  let last = 0;
  for (const match of text.matchAll(re)) {
    if (match.index > last) fragment.append(text.slice(last, match.index));
    const mark = document.createElement('mark');
    mark.textContent = match[0];
    fragment.append(mark);
    last = match.index + match[0].length;
  }
  if (last < text.length) fragment.append(text.slice(last));
  return fragment;
}

/** One line of explanation: what matched, and where it sits in the plasmid. */
function whyLine(hit, terms) {
  const line = document.createElement('span');
  line.className = 'finder-why';
  const where = {
    part: 'in the part', cassette: 'in the cassette', component: 'in the part',
    name: 'in the name', backbone: 'in the backbone only',
  }[hit.where] || '';
  const text = hit.matched.slice(0, 3).join(' · ') || hit.component || hit.name;
  line.append(highlight(text, terms));
  const tail = document.createElement('span');
  if (hit.component_source === 'filename' && hit.where === 'name') {
    tail.className = 'is-guessed';
    tail.textContent = ` — ${where}, no annotation in the file`;
  } else {
    tail.textContent = ` — ${where}`;
  }
  line.append(tail);
  return line;
}

function hitRow(hit, index, terms) {
  const row = document.createElement('button');
  row.type = 'button';
  row.className = 'finder-hit';
  row.setAttribute('role', 'option');
  row.setAttribute('aria-selected', String(index === finder.active));

  const mark = document.createElement('span');
  mark.className = 'badge';
  mark.textContent = hit.part_type ? hit.part_type : '—';
  mark.dataset.part = hit.part_type || '';
  row.append(mark);

  const what = document.createElement('span');
  what.className = 'finder-what';
  const name = document.createElement('span');
  name.className = 'finder-name';
  name.append(highlight(hit.display || hit.name, terms));
  if (hit.internal_sites) name.append(' ⚠');
  what.append(name, whyLine(hit, terms));
  row.append(what);

  const len = document.createElement('span');
  len.className = 'finder-len';
  len.textContent = `${hit.length} bp`;
  row.append(len);

  if (!hit.part_type) {
    row.classList.add('is-unusable');
    row.disabled = true;
    row.title = hit.roles.includes('cassette')
      ? 'a finished cassette — use it on the Multigene page'
      : 'no BsaI part could be read out of this plasmid';
  } else {
    row.addEventListener('click', () => chooseHit(index));
  }
  return row;
}

/** Step two: ask before changing the design. */
function confirmBar(hit) {
  const bar = document.createElement('div');
  bar.className = 'finder-confirm';

  const ask = document.createElement('span');
  ask.className = 'finder-ask';
  const name = document.createElement('b');
  name.textContent = hit.display || hit.name;
  ask.append('Use ', name, ` as the type ${hit.part_type} part?`);
  const view = VIEW_FOR_TYPE[hit.part_type];
  if (view && state[view[0]] !== view[1]) {
    const note = document.createElement('span');
    note.textContent = ` This switches that column to ${
      view[0] === 'composite_left' ? '2·3·4'
      : view[0] === 'composite_right' ? '6·7·8'
      : view[1] ? hit.part_type.replace(/[ab]$/, '') + 'a+' + hit.part_type.replace(/[ab]$/, '') + 'b'
      : 'a single ' + hit.part_type + ' panel'
    }.`;
    ask.append(note);
  }
  bar.append(ask);

  const cancel = document.createElement('button');
  cancel.type = 'button';
  cancel.className = 'btn';
  cancel.textContent = 'Cancel';
  cancel.addEventListener('click', () => { finder.chosen = null; renderFinder(); });

  const use = document.createElement('button');
  use.type = 'button';
  use.className = 'btn btn-primary';
  use.textContent = 'Use it';
  use.addEventListener('click', () => applyHit(hit));

  bar.append(cancel, use);
  return bar;
}

function chooseHit(index) {
  finder.active = index;
  finder.chosen = finder.hits[index];
  renderFinder();
  const use = el('finder-results').querySelector('.finder-confirm .btn-primary');
  if (use) use.focus();
}

/** Put the part in its slot, opening the view that slot lives in. */
function applyHit(hit) {
  const view = VIEW_FOR_TYPE[hit.part_type];
  if (view) {
    const [flag, wanted] = view;
    if (state[flag] !== wanted) {
      state[flag] = wanted;
      // the panels this replaces are gone, so their choices go with them
      const dropped = flag === 'composite_left'
        ? ['2', '3', '3a', '3b', '4', '4a', '4b']
        : flag === 'composite_right'
          ? ['6', '7', '8', '8a', '8b']
          : wanted
            ? [flag.slice(-1)]
            : [`${flag.slice(-1)}a`, `${flag.slice(-1)}b`];
      for (const key of dropped) delete state.selections[key];
    }
  }
  // a composite part needs its column collapsed; a plain one needs it open
  if (hit.part_type === '234') state.composite_left = true;
  else if (LEFT_TYPES.has(hit.part_type)) state.composite_left = false;
  if (hit.part_type === '678') state.composite_right = true;
  else if (['6', '7', '8', '8a', '8b'].includes(hit.part_type)) state.composite_right = false;

  state.selections[hit.part_type] = hit.name;
  closeFinder();
  refresh();
}

function renderFinder() {
  const box = el('finder-results');
  box.replaceChildren();
  const query = el('finder-input').value.trim();

  if (!query) { box.hidden = true; return; }
  box.hidden = false;

  /* Only part plasmids are offered here. A match that cannot be cut into a
   * part - a finished construct still carrying the labels of the parts that
   * built it, or a destination vector - is withheld rather than shown greyed,
   * but it is counted out loud: silently returning nothing for a name that is
   * plainly in the library reads as a broken search. */
  const withheld = finder.suppressed
    ? `${finder.suppressed} other match${finder.suppressed === 1 ? '' : 'es'} `
      + 'are not part plasmids — finished constructs or destination vectors. '
      + 'The Library tab lists them.'
    : '';

  if (!finder.hits.length) {
    const hint = document.createElement('div');
    hint.className = 'finder-hint';
    hint.textContent = finder.suppressed
      ? `Nothing annotated “${query}” can be used as a part. ${withheld}`
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
    `/api/library/search?q=${encodeURIComponent(query)}&limit=20&usable=part`,
  );
  if (seq !== finder.seq) return; // a later keystroke already answered
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
    const usable = finder.hits.filter((h) => h.part_type).length;
    if (event.key === 'Escape') { closeFinder(); return; }
    if (!finder.hits.length) return;
    if (event.key === 'ArrowDown' || event.key === 'ArrowUp') {
      event.preventDefault();
      const step = event.key === 'ArrowDown' ? 1 : -1;
      let next = finder.active;
      for (let i = 0; i < finder.hits.length; i += 1) {
        next = (next + step + finder.hits.length) % finder.hits.length;
        if (finder.hits[next].part_type) break;
      }
      finder.active = usable ? next : -1;
      renderFinder();
    } else if (event.key === 'Enter' && finder.active >= 0) {
      event.preventDefault();
      chooseHit(finder.active);
    }
  });
  document.addEventListener('click', (event) => {
    if (!event.target.closest('.finder')) {
      el('finder-results').hidden = true;
    }
  });
  input.addEventListener('focus', () => {
    if (finder.hits.length) el('finder-results').hidden = false;
  });
}

// ------------------------------------------------------------------- clone ---

/* Reverse mode: read a finished cassette back into its parts.
 *
 * The screen was only ever "build eight from scratch", but the thing you
 * actually do is take a construct that works and change one part of it. The
 * assembly consumed every BsaI site, so the server cannot digest it back - it
 * matches each library part's fragment against the sequence instead, and the
 * stretches nothing explains come back as gaps rather than being guessed at.
 */

let cassettes = [];

function cloneRow(entry) {
  const row = document.createElement('button');
  row.type = 'button';
  row.className = 'pick-row';
  const what = document.createElement('span');
  what.className = 'pick-what';
  const name = document.createElement('span');
  name.className = 'pick-name';
  name.textContent = entry.display || entry.name;
  what.append(name);
  if (entry.component) {
    const sub = document.createElement('span');
    sub.className = 'pick-sub';
    sub.textContent = entry.component;
    what.append(sub);
  }
  const len = document.createElement('span');
  len.className = 'pick-len';
  len.textContent = `${entry.length} bp`;
  row.append(what, len);
  row.addEventListener('click', () => cloneFrom(entry.name));
  return row;
}

function renderCloneList() {
  const needle = el('clone-filter').value.trim().toLowerCase();
  const terms = needle.split(/\s+/).filter(Boolean);
  const shown = cassettes.filter((entry) => {
    if (!terms.length) return true;
    const hay = [entry.name, entry.display, entry.component, ...(entry.aliases || [])]
      .filter(Boolean).join(' ').toLowerCase();
    return terms.every((t) => hay.includes(t));
  });
  const list = el('clone-list');
  list.replaceChildren();
  if (!shown.length) {
    const hint = document.createElement('div');
    hint.className = 'pick-hint';
    hint.textContent = 'No cassette matches that.';
    list.append(hint);
    return;
  }
  for (const entry of shown.slice(0, 200)) list.append(cloneRow(entry));
  if (shown.length > 200) {
    const more = document.createElement('div');
    more.className = 'pick-hint';
    more.textContent = `${shown.length - 200} more — keep typing to narrow it down.`;
    list.append(more);
  }
}

async function cloneFrom(name) {
  const result = await post('/api/level2/decompose', { cassette: name });
  el('clone-dialog').close();

  state.selections = result.selections || {};
  state.split_3 = result.split_3;
  state.split_4 = result.split_4;
  state.split_8 = result.split_8;
  state.composite_left = result.composite_left;
  state.composite_right = result.composite_right;
  state.name = result.name;
  el('construct-name').value = result.name;

  renderPanels(result.slots);
  renderStage(result);
  lastDrawn = result;
  drawRing(result.parts, result.length);
  renderLinear(result.parts);
  renderJunctions(result.junctions, result.issues);
  renderIssues(result.issues);
  reportClone(result.source);
}

/* What the match found, said plainly. A cassette built with a part that is not
 * in this library cannot be filled in completely, and pretending otherwise
 * would hand you a design that quietly differs from the thing on your bench. */
function reportClone(source) {
  if (!source) return;
  const box = el('clone-report');
  box.hidden = false;
  box.replaceChildren();
  const text = document.createElement('span');
  if (source.complete) {
    text.textContent =
      `Every base of ${source.display} is accounted for by ${source.matches.length} `
      + 'library parts. Change any panel and rebuild.';
  } else {
    box.classList.add('is-partial');
    // the overhangs at each end of a gap say which position the unknown part
    // occupies, which is a far more useful thing to report than a bare length
    const where = (source.unmatched || [])
      .map((u) => `${u.length.toLocaleString()} bp`
        + (u.part_type ? ` at position ${u.part_type}` : '')
        + (u.left_overhang ? ` (${u.left_overhang}→${u.right_overhang})` : ''))
      .join(', ');
    text.textContent =
      `${source.matches.length} of ${source.display}'s parts were found in the library, `
      + `covering ${source.covered.toLocaleString()} of ${source.length.toLocaleString()} bp. `
      + `Unaccounted for: ${where} — no part on the shelf has that sequence, so the `
      + 'panels below are what could be identified, not the whole construct.';
  }
  box.append(text);
  const dismiss = document.createElement('button');
  dismiss.type = 'button';
  dismiss.className = 'btn';
  dismiss.textContent = 'Dismiss';
  dismiss.addEventListener('click', () => { box.hidden = true; });
  box.append(dismiss);
}

async function openClone() {
  if (!cassettes.length) {
    const response = await fetch('/api/library/plasmids?role=cassette');
    cassettes = (await response.json()).sort((a, b) => a.name.localeCompare(b.name));
  }
  renderCloneList();
  el('clone-dialog').showModal();
  el('clone-filter').focus();
}

function wireClone() {
  el('clone-btn').addEventListener('click', openClone);
  el('clone-filter').addEventListener('input', renderCloneList);
}

// ---------------------------------------------------------------- protocol ---

/* What to pipette, and the cycling.
 *
 * The volumes are equimolar, which needs a concentration per plasmid - the one
 * number no GenBank file can supply. Rather than assume one, a piece with no
 * recorded concentration is shown with an input in place of its volume, so the
 * gap is fixable at the moment you notice it; entering a number writes it to
 * the library and the table refills.
 */

let protocolText = '';

function num(value, digits = 2) {
  const td = document.createElement('td');
  td.className = 'num';
  td.textContent = value === null || value === undefined ? '' : value.toFixed(digits);
  return td;
}

/** The volume cell, or an input when the concentration is not known yet. */
function volumeCell(component) {
  if (component.volume_ul !== null && component.volume_ul !== undefined) {
    return num(component.volume_ul);
  }
  const td = document.createElement('td');
  td.className = 'num';
  const input = document.createElement('input');
  input.className = 'rx-conc';
  input.type = 'number';
  input.min = '0.1';
  input.step = '0.1';
  input.placeholder = 'ng/µL';
  input.title = `measure ${component.name} and enter it here`;
  input.addEventListener('change', async () => {
    const value = Number(input.value);
    if (!(value > 0)) return;
    input.disabled = true;
    await post('/api/library/concentration', {
      path: component.path,
      conc_ng_ul: value,
    });
    await showProtocol();
  });
  td.append(input);
  return td;
}

function renderProtocol(rx) {
  const body = el('protocol-body');
  body.replaceChildren();

  for (const issue of rx.issues || []) {
    const warn = document.createElement('div');
    warn.className = 'rx-warn';
    warn.textContent = issue;
    body.append(warn);
  }

  const lead = document.createElement('p');
  lead.className = 'rx-lead';
  lead.textContent =
    `${rx.fmol_each} fmol of each piece in ${rx.total_ul} µL. `
    + 'Equimolar, not equal by mass — a short part and a long backbone at the '
    + 'same ng/µL differ several-fold in copy number.';
  body.append(lead);

  const table = document.createElement('table');
  table.className = 'rx-table';
  const head = document.createElement('tr');
  for (const [label, cls] of [['Component', ''], ['µL', 'num'], ['ng', 'num'],
                              ['bp', 'num'], ['ng/µL', 'num'], ['', '']]) {
    const th = document.createElement('th');
    th.textContent = label;
    if (cls) th.className = cls;
    head.append(th);
  }
  table.append(head);

  for (const c of rx.components) {
    const tr = document.createElement('tr');
    if (c.kind === 'reagent' || c.kind === 'water') tr.classList.add('is-reagent');
    if (c.kind !== 'reagent' && c.kind !== 'water' && c.volume_ul === null) {
      tr.classList.add('is-missing');
    }
    const name = document.createElement('td');
    name.textContent = c.name;
    tr.append(name, volumeCell(c), num(c.ng, 1));
    const bp = document.createElement('td');
    bp.className = 'num';
    bp.textContent = c.length ? c.length.toLocaleString() : '';
    tr.append(bp, num(c.conc_ng_ul, 1));
    const note = document.createElement('td');
    note.className = 'rx-note';
    note.textContent = c.note;
    tr.append(note);
    table.append(tr);
  }

  const total = document.createElement('tr');
  total.className = 'rx-total';
  const label = document.createElement('td');
  label.textContent = 'Total';
  const sum = rx.components.reduce((a, c) => a + (c.volume_ul || 0), 0);
  total.append(label, num(sum));
  for (let i = 0; i < 4; i += 1) total.append(document.createElement('td'));
  table.append(total);
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

  if (rx.selection) {
    const sel = document.createElement('p');
    sel.className = 'rx-lead';
    sel.style.marginTop = '12px';
    sel.textContent = `Selection: ${rx.selection}`;
    body.append(sel);
  }
}

async function showProtocol() {
  const rx = await post('/api/level2/protocol', state);
  if (rx.ok === false && !rx.components) {
    el('protocol-body').textContent =
      `This design does not assemble yet: ${(rx.issues || []).join('; ')}`;
    protocolText = '';
  } else {
    protocolText = rx.text || '';
    renderProtocol(rx);
  }
  el('protocol-title').textContent = `${rx.name || state.name} · ${rx.enzyme || ''} reaction`;
  const dialog = el('protocol-dialog');
  if (!dialog.open) dialog.showModal();
}

function wireProtocol() {
  el('protocol-btn').addEventListener('click', showProtocol);
  el('protocol-copy').addEventListener('click', async () => {
    if (!protocolText) return;
    await navigator.clipboard.writeText(protocolText);
    const button = el('protocol-copy');
    button.textContent = 'Copied';
    setTimeout(() => { button.textContent = 'Copy'; }, 1500);
  });
  el('protocol-download').addEventListener('click', async () => {
    const response = await fetch('/api/level2/protocol.txt', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(state),
    });
    if (!response.ok) return;
    const url = URL.createObjectURL(await response.blob());
    const link = document.createElement('a');
    link.href = url;
    link.download = `${state.name}-protocol.txt`;
    link.click();
    URL.revokeObjectURL(url);
  });
}

// ------------------------------------------------------------------ wiring ---

let pending = null;

async function refresh() {
  state.name = el('construct-name').value.trim() || 'cassette';
  const mine = {};
  pending = mine;
  const result = await post('/api/level2/assemble', state);
  if (pending !== mine) return; // a newer request already answered

  renderPanels(result.slots);
  renderStage(result);
  drawRing(result.parts, result.length);
  renderLinear(result.parts);
  renderJunctions(result.junctions, result.issues);
  renderIssues(result.issues);
}

function wireViewToggles() {
  for (const button of document.querySelectorAll('[data-view]')) {
    button.addEventListener('click', () => {
      const view = button.dataset.view;
      if (view.endsWith('-left')) {
        const on = view.startsWith('composite');
        if (state.composite_left === on) return;
        state.composite_left = on;
        for (const key of ['2', '3', '3a', '3b', '4', '4a', '4b', '234']) {
          delete state.selections[key];
        }
      } else {
        const on = view.startsWith('composite');
        if (state.composite_right === on) return;
        state.composite_right = on;
        for (const key of ['6', '7', '8', '8a', '8b', '678']) delete state.selections[key];
      }
      refresh();
    });
  }
}

/** Open on a real assembly rather than eight empty dropdowns. */
async function loadDefault() {
  try {
    const result = await (await fetch('/api/level2/default')).json();
    if (!result.selections || !Object.keys(result.selections).length) return false;
    Object.assign(state.selections, result.selections);
    if (result.name) {
      state.name = result.name;
      el('construct-name').value = result.name;
    }
    renderPanels(result.slots);
    renderStage(result);
    drawRing(result.parts, result.length);
    renderLegend(result.parts);
  renderLinear(result.parts);
    renderJunctions(result.junctions, result.issues);
    renderIssues(result.issues);
    return true;
  } catch {
    return false;
  }
}

async function loadSummary() {
  const response = await fetch('/api/library/summary');
  const summary = await response.json();
  const roots = summary.roots.length === 1 ? summary.roots[0] : `${summary.roots.length} folders`;
  el('folder').textContent = `${roots} · ${summary.parts} parts indexed`;
  el('rx-enzyme').textContent = summary.enzymes.part;
}

function wireExport() {
  el('export-btn').addEventListener('click', async () => {
    const response = await fetch('/api/level2/export', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(state),
    });
    if (!response.ok) return;
    const blob = await response.blob();
    const url = URL.createObjectURL(blob);
    const link = document.createElement('a');
    link.href = url;
    link.download = `${state.name}.gb`;
    link.click();
    URL.revokeObjectURL(url);
  });
}

el('construct-name').addEventListener('change', refresh);
wireViewToggles();
wireFinder();
loadView();
wireMapTools();
wireClone();
wireSaved();
wireProtocol();
wireExport();
loadSummary();
loadSaved();
loadDefault().then((loaded) => {
  if (!loaded) refresh();
});

// ------------------------------------------------------------------- saved ---

/* Saved designs.
 *
 * A design is eight plasmid names and five flags. Nothing persisted before
 * this, so closing the tab threw away every choice - which is what makes a
 * tool feel disposable. Kept server-side in .ggasm/designs.json rather than in
 * localStorage, because the design belongs with the library it refers to, and
 * a name that no longer exists there should be visible as such.
 */

async function loadSaved() {
  const select = el('saved-select');
  const keep = select.value;
  let saved = [];
  try {
    saved = await (await fetch('/api/library/designs?level=level2')).json();
  } catch {
    return;
  }
  select.replaceChildren(select.firstElementChild);
  for (const record of saved) {
    const option = document.createElement('option');
    option.value = record.name;
    option.textContent = record.name;
    option.title = `saved ${record.saved_at}`;
    select.append(option);
  }
  select.value = keep;
  savedDesigns = saved;
}

let savedDesigns = [];

function wireSaved() {
  el('save-btn').addEventListener('click', async () => {
    const name = prompt('Save this design as:', state.name || 'cassette');
    if (!name || !name.trim()) return;
    await post('/api/library/designs', {
      level: 'level2', name: name.trim(), design: state,
    });
    await loadSaved();
    el('saved-select').value = name.trim();
  });

  el('saved-select').addEventListener('change', async () => {
    const record = savedDesigns.find((d) => d.name === el('saved-select').value);
    if (!record) return;
    Object.assign(state, record.design);
    el('construct-name').value = state.name || 'pCassette';
    await refresh();
  });
}

// ------------------------------------------------------------------ linear ---

/* The construct straightened out.
 *
 * The ring is proportional and correct, but a circle is a poor place to read
 * an order off. What you check before ordering reads as a line: does the promoter drive the
 * CDS, does the terminator follow it. Blocks are sized by bp, and a part whose
 * largest feature is stranded gets an arrow end, so direction is shown rather
 * than implied.
 */

function renderLinear(parts) {
  const box = el('linear');
  box.replaceChildren();
  box.hidden = !state.showLinear || !parts.length;
  if (box.hidden) return;

  const track = document.createElement('div');
  track.className = 'linear-track';

  for (const part of parts) {
    const cell = document.createElement('div');
    cell.className = 'linear-cell';
    cell.style.flexGrow = String(part.length);

    const block = document.createElement('button');
    block.type = 'button';
    block.className = 'linear-block';
    block.dataset.part = part.part_type || '';
    block.style.background = 'var(--part-color)';
    if (strandOf(part)) block.classList.add('is-directional');
    block.title = `Type ${part.part_type} · ${part.label || part.source_name} · `
      + `${part.length.toLocaleString()} bp`;
    block.addEventListener('click', () => focusPanel(part.part_type));

    const name = document.createElement('span');
    name.className = 'linear-name';
    name.textContent = part.label || part.source_name;
    block.append(name);

    const tick = document.createElement('span');
    tick.className = 'linear-tick';
    tick.textContent = part.part_type || '';

    cell.append(block, tick);
    track.append(cell);
  }

  box.append(track);

  const caption = document.createElement('span');
  caption.className = 'card-note';
  caption.textContent = '5′ → 3′ · click a block to jump to its panel';
  box.append(caption);
}

// ------------------------------------------------------------------ sweep ---

/* Vary one position across many parts, holding the rest of the design still.
 *
 * A promoter titration against a fixed coding sequence is the standard shape
 * of this: forty-seven type 2 parts, one construct each. Doing it by hand is
 * forty-seven passes through this screen to learn the same two things every
 * time - does it assemble, and how long does it come out.
 */

let sweepSlots = [];
let sweepRows = [];

function sweepBody() {
  const chosen = [...document.querySelectorAll('#sweep-candidates input:checked')]
    .map((box) => box.value);
  return { base_design: state, slot: el('sweep-slot').value, candidates: chosen };
}

function renderCandidates() {
  const slot = sweepSlots.find((s) => s.key === el('sweep-slot').value);
  const host = el('sweep-candidates');
  host.replaceChildren();
  if (!slot) return;

  for (const option of slot.options) {
    const row = document.createElement('label');
    row.className = 'sweep-candidate';
    const box = document.createElement('input');
    box.type = 'checkbox';
    box.value = option.name;
    box.checked = true;
    box.addEventListener('change', updateSweepCount);
    const text = document.createElement('span');
    text.textContent = `${option.display || option.name} · ${option.length} bp`;
    row.append(box, text);
    host.append(row);
  }
  updateSweepCount();
}

function updateSweepCount() {
  const n = document.querySelectorAll('#sweep-candidates input:checked').length;
  el('sweep-count').textContent = `${n} selected`;
  el('sweep-write').disabled = !sweepRows.length;
}

function renderSweepResults(payload) {
  sweepRows = payload.rows;
  const host = el('sweep-results');
  host.replaceChildren();

  const summary = document.createElement('p');
  summary.className = 'rx-lead';
  summary.textContent =
    `${payload.built} of ${payload.count} assemble at position ${payload.slot}.`;
  host.append(summary);

  const table = document.createElement('table');
  table.className = 'rx-table';
  const head = document.createElement('tr');
  for (const [label, cls] of [['Part', ''], ['Part bp', 'num'], ['Construct bp', 'num'],
                              ['Builds', ''], ['Notes', '']]) {
    const th = document.createElement('th');
    th.textContent = label;
    if (cls) th.className = cls;
    head.append(th);
  }
  table.append(head);

  for (const row of payload.rows) {
    const tr = document.createElement('tr');
    if (!row.ok) tr.classList.add('is-missing');
    const cells = [
      [row.display || row.name, ''],
      [row.part_length.toLocaleString(), 'num'],
      [row.ok ? row.length.toLocaleString() : '', 'num'],
      [row.ok ? 'yes' : 'no', ''],
      [[...row.errors, ...row.warnings].join('; '), 'rx-note'],
    ];
    for (const [text, cls] of cells) {
      const td = document.createElement('td');
      if (cls) td.className = cls;
      td.textContent = text;
      tr.append(td);
    }
    table.append(tr);
  }
  host.append(table);
  updateSweepCount();
}

function wireSweep() {
  el('sweep-btn').addEventListener('click', async () => {
    const payload = await post('/api/level2/slots', state);
    sweepSlots = payload.slots.filter((s) => s.options.length);
    const picker = el('sweep-slot');
    picker.replaceChildren();
    for (const slot of sweepSlots) {
      const option = document.createElement('option');
      option.value = slot.key;
      option.textContent = `Type ${slot.key} · ${slot.options.length} parts`;
      picker.append(option);
    }
    sweepRows = [];
    el('sweep-results').replaceChildren();
    renderCandidates();
    el('sweep-dialog').showModal();
  });

  el('sweep-slot').addEventListener('change', () => {
    sweepRows = [];
    el('sweep-results').replaceChildren();
    renderCandidates();
  });
  el('sweep-all').addEventListener('click', () => {
    for (const box of document.querySelectorAll('#sweep-candidates input')) box.checked = true;
    updateSweepCount();
  });
  el('sweep-none').addEventListener('click', () => {
    for (const box of document.querySelectorAll('#sweep-candidates input')) box.checked = false;
    updateSweepCount();
  });

  el('sweep-run').addEventListener('click', async () => {
    const button = el('sweep-run');
    button.disabled = true;
    button.textContent = 'Building…';
    try {
      renderSweepResults(await post('/api/level2/sweep', sweepBody()));
    } finally {
      button.disabled = false;
      button.textContent = 'Build all';
    }
  });

  el('sweep-write').addEventListener('click', async () => {
    const response = await fetch('/api/level2/sweep.zip', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(sweepBody()),
    });
    if (!response.ok) return;
    const url = URL.createObjectURL(await response.blob());
    const link = document.createElement('a');
    link.href = url;
    link.download = `${state.name}-sweep.zip`;
    link.click();
    URL.revokeObjectURL(url);
  });
}

wireSweep();
