/* Sequencing - check clones against the plasmid they were meant to be.
 *
 * The files never leave this machine: they are read with FileReader and posted
 * as text to the local server, which is the same thing every other screen here
 * does with the folder on disk.
 *
 * Two views of the same answer, because they answer different questions.
 * "Around each difference" is what you read when you want to know whether a
 * clone is usable: a short window per event, in place, named by the feature it
 * falls in. "The whole alignment" is what you read when you want to see the
 * construct end to end - ten thousand columns of it, so only the ones on
 * screen are ever in the DOM.
 */

import { index as columnIndex, placed, span as columnSpan, ticks } from './columns.js';

const el = (id) => document.getElementById(id);
const api = (path, options) => fetch(path, options).then(async (response) => {
  const body = await response.json().catch(() => ({}));
  if (!response.ok) throw new Error(body.detail || `${response.status} ${response.statusText}`);
  return body;
});

/* ---------------------------------------------------------------- state */

const state = {
  source: 'library',
  references: { library: [], expected: [] },
  chosen: '',
  pasted: 0,
  repaint: null,
  referenceFile: null,
  clones: [],
  result: null,
  jumpAt: -1,
};

/* Column geometry. Measured rather than assumed: the monospace face is
 * vendored, and a fallback font would make every tick label point at the wrong
 * base if this were hard-coded. */
let CH = 8;

function measure() {
  const probe = document.createElement('span');
  probe.className = 'aln-probe';
  probe.textContent = 'A'.repeat(100);
  document.body.appendChild(probe);
  const width = probe.getBoundingClientRect().width / 100;
  probe.remove();
  if (width > 1) CH = width;
  document.documentElement.style.setProperty('--ch', `${CH}px`);
}

/* ------------------------------------------------------------ references */

async function loadReferences() {
  try {
    state.references = await api('/api/sequencing/references');
  } catch (error) {
    setStatus(`could not read the library: ${error.message}`, true);
    return;
  }
  renderReferences();
}

function currentList() {
  return state.source === 'expected' ? state.references.expected : state.references.library;
}

const LIST_LIMIT = 500;

function renderReferences() {
  const list = currentList();
  const select = el('reference');
  const kept = select.value;

  select.replaceChildren(...list.slice(0, LIST_LIMIT).map((r) => option(r)));
  if (list.some((r) => r.name === kept)) select.value = kept;
  // a plasmid picked out of the search may sit past the browse list's cut-off;
  // it still has to be selectable, so it is added rather than the list grown
  if (state.chosen && !list.slice(0, LIST_LIMIT).some((r) => r.name === state.chosen)) {
    const found = list.find((r) => r.name === state.chosen);
    if (found) select.prepend(option(found));
  }
  if (state.chosen && list.some((r) => r.name === state.chosen)) select.value = state.chosen;

  const shown = Math.min(list.length, LIST_LIMIT);
  el('ref-note').textContent = state.source === 'expected'
    ? (list.length
        ? `${list.length} kept from the Cassette and Multigene screens`
        : 'nothing kept yet — export a construct from Cassette or Multigene and it appears here')
    : (list.length > shown
        ? `${shown.toLocaleString()} of ${list.length.toLocaleString()} plasmids — search to reach the rest`
        : `${list.length.toLocaleString()} plasmids`);
  updateRunButton();
}

function option(r) {
  const node = document.createElement('option');
  node.value = r.name;
  const what = r.part_type ? `type ${r.part_type}` : ((r.roles || [])[0] || '');
  node.textContent = `${r.name} · ${r.length.toLocaleString()} bp${what ? ` · ${what}` : ''}`;
  return node;
}

function setSource(source) {
  state.source = source;
  for (const [id, value] of [['src-library', 'library'], ['src-expected', 'expected'], ['src-upload', 'upload']]) {
    el(id).setAttribute('aria-pressed', String(value === source));
  }
  el('pick-named').hidden = source === 'upload';
  el('pick-file').hidden = source !== 'upload';
  closeFinder();
  if (source === 'upload') {
    el('ref-note').textContent = state.referenceFile ? state.referenceFile.name : 'no file chosen';
    updateRunButton();
  } else {
    renderReferences();
  }
}

/* ---------------------------------------------------------------- finder */

/* The same shared search (/api/library/search) the Cassette and Multigene
 * screens wear, read the way this screen needs it.
 *
 * Two deliberate differences from those two, both because a reference is
 * compared against rather than assembled with:
 *
 * No `usable=` scope. There, a hit has to fit a slot or close a chain, so a
 * finished construct offered as a part is a real mistake and gets withheld.
 * Here anything with a sequence is a legitimate thing to have sequenced - a
 * part plasmid, a cassette, a 12 kb multigene construct, even a plasmid the
 * digest could not place. Withholding any of them would only hide the one you
 * actually sent off.
 *
 * And no confirm step. There, picking a hit rewrites a design you have been
 * building, so it asks first. Here it selects a reference, which is one click
 * to change and destroys nothing, and a confirm bar would be friction
 * guarding against no risk.
 */

const finder = { hits: [], active: -1, seq: 0 };

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

function hitRow(hit, position, terms) {
  const row = document.createElement('button');
  row.type = 'button';
  row.className = 'finder-hit';
  row.setAttribute('role', 'option');
  row.setAttribute('aria-selected', String(position === finder.active));

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
  const why = document.createElement('span');
  why.className = 'finder-why';
  why.append(highlight(
    (hit.matched || []).slice(0, 3).join(' · ') || hit.component || hit.name, terms,
  ));
  what.append(name, why);
  row.append(what);

  const len = document.createElement('span');
  len.className = 'finder-len';
  len.textContent = `${hit.length.toLocaleString()} bp`;
  row.append(len);

  row.addEventListener('click', () => chooseHit(position));
  return row;
}

function chooseHit(position) {
  const hit = finder.hits[position];
  if (!hit) return;
  state.chosen = hit.name;
  closeFinder();
  renderReferences();
  setStatus(`reference: ${hit.display || hit.name}`);
}

function renderFinder() {
  const box = el('finder-results');
  box.replaceChildren();
  const query = el('finder-input').value.trim();
  if (!query) { box.hidden = true; return; }
  box.hidden = false;

  if (!finder.hits.length) {
    const hint = document.createElement('div');
    hint.className = 'finder-hint';
    hint.textContent = state.source === 'expected'
      ? `No construct you have exported is called “${query}”.`
      : `Nothing in the library is named or annotated “${query}”.`;
    box.append(hint);
    return;
  }
  const terms = query.split(/\s+/).filter(Boolean);
  finder.hits.forEach((hit, i) => box.append(hitRow(hit, i, terms)));
}

function closeFinder() {
  finder.hits = [];
  finder.active = -1;
  el('finder-input').value = '';
  el('finder-results').hidden = true;
  el('finder-results').replaceChildren();
}

async function runSearch() {
  const query = el('finder-input').value.trim();
  finder.active = -1;
  if (query.length < 2) {
    finder.hits = [];
    renderFinder();
    return;
  }

  // The kept constructs are a handful and are already here, so searching them
  // over the network would be a request to be told what this page knows.
  if (state.source === 'expected') {
    const term = query.toLowerCase();
    finder.hits = state.references.expected.filter((r) => r.name.toLowerCase().includes(term));
    renderFinder();
    return;
  }

  const seq = ++finder.seq;
  const response = await fetch(`/api/library/search?q=${encodeURIComponent(query)}&limit=20`);
  if (seq !== finder.seq) return;   // an older keystroke, answered late
  const result = await response.json();
  finder.hits = result.hits || [];
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
      finder.active = (finder.active + step + finder.hits.length) % finder.hits.length;
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

/* ----------------------------------------------------------------- files */

function readFile(file) {
  return new Promise((resolve, reject) => {
    const reader = new FileReader();
    reader.onerror = () => reject(new Error(`could not read ${file.name}`));
    reader.onload = () => resolve({ name: file.name, text: String(reader.result || '') });
    reader.readAsText(file);
  });
}

async function addClones(files) {
  const incoming = await Promise.all([...files].map(readFile));
  const names = new Set(state.clones.map((c) => c.name));
  for (const file of incoming) {
    if (!names.has(file.name)) state.clones.push(file);
  }
  renderFiles();
}

/* What to call a sequence that arrived without a file name.
 *
 * A file name is the best identifier there is - it is what is written on the
 * tube - and pasted text has none. What it often does have is a header the
 * person already typed, so a FASTA `>` line or a GenBank LOCUS is read before
 * falling back to counting. The name shown in the list is the name sent, so
 * what you see here is what the report calls it.
 */
function nameForPaste(text, typed) {
  if (typed.trim()) return unique(typed.trim());
  const first = text.trim().split('\n', 1)[0].trim();
  if (first.startsWith('>')) {
    const header = first.slice(1).trim().split(/\s+/)[0];
    if (header) return unique(header);
  }
  if (first.startsWith('LOCUS')) {
    const locus = first.split(/\s+/)[1];
    if (locus) return unique(locus);
  }
  state.pasted += 1;
  return unique(`pasted ${state.pasted}`);
}

/** The same name twice would let one clone quietly replace another. */
function unique(name) {
  const taken = new Set(state.clones.map((c) => c.name));
  if (!taken.has(name)) return name;
  for (let n = 2; ; n += 1) {
    if (!taken.has(`${name} (${n})`)) return `${name} (${n})`;
  }
}

function addPasted() {
  const text = el('paste-text').value;
  const note = el('paste-note');
  if (!text.trim()) {
    note.textContent = 'nothing to add';
    return;
  }
  const name = nameForPaste(text, el('paste-name').value);
  state.clones.push({ name, text });
  el('paste-text').value = '';
  el('paste-name').value = '';
  // whether it is really sequence is the server's call, and it names the
  // problem precisely; guessing here would only be a second, worse answer
  note.textContent = `added as ${name}`;
  renderFiles();
}

function renderFiles() {
  el('file-list').replaceChildren(...state.clones.map((file, index) => {
    const li = document.createElement('li');
    const name = document.createElement('span');
    name.className = 'file-name';
    name.textContent = file.name;
    const size = document.createElement('span');
    size.className = 'file-size';
    size.textContent = `${file.text.length.toLocaleString()} chars`;
    const drop = document.createElement('button');
    drop.className = 'file-x';
    drop.type = 'button';
    drop.title = `Remove ${file.name}`;
    drop.textContent = '×';
    drop.addEventListener('click', () => {
      state.clones.splice(index, 1);
      renderFiles();
    });
    li.append(name, size, drop);
    return li;
  }));
  el('clone-note').textContent = state.clones.length
    ? `${state.clones.length} file${state.clones.length > 1 ? 's' : ''}`
    : 'one file per clone';
  el('clear-btn').disabled = !state.clones.length;
  updateRunButton();
}

function updateRunButton() {
  const haveReference = state.source === 'upload'
    ? Boolean(state.referenceFile)
    : Boolean(el('reference').value);
  el('run-btn').disabled = !(haveReference && state.clones.length);
}

/* ------------------------------------------------------------------- run */

async function run() {
  el('run-btn').disabled = true;
  setStatus('aligning…');
  try {
    state.result = await api('/api/sequencing/run', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        source: state.source,
        reference: state.source === 'upload' ? '' : el('reference').value,
        reference_file: state.referenceFile,
        clones: state.clones,
      }),
    });
  } catch (error) {
    setStatus(error.message, true);
    updateRunButton();
    return;
  }
  setStatus('');
  index(state.result);
  render();
  updateRunButton();
  el('report').scrollIntoView({ behavior: 'smooth', block: 'start' });
}

/* The column arithmetic lives in columns.js, with no DOM in it, so it can be
 * run under node - see tests/test_columns.py. */
function index(result) {
  const { colOf, refOf } = columnIndex(result.alignment.insertions, result.reference.length);
  result.colOf = colOf;
  result.refOf = refOf;

  result.events = [];
  result.clones.forEach((clone, row) => {
    for (const d of clone.differences) result.events.push({ ...d, clone: clone.name, row });
  });
  result.events.sort((a, b) => a.start - b.start);
  state.jumpAt = -1;
}

/* -------------------------------------------------------------- verdicts */

function render() {
  const result = state.result;
  el('report').hidden = false;
  el('headline').textContent = result.summary.headline;
  el('headline-note').textContent =
    `${result.reference.name} · ${result.reference.length.toLocaleString()} bp · `
    + `${result.summary.checked} clone${result.summary.checked === 1 ? '' : 's'} checked`;

  el('verdicts').replaceChildren(
    ...result.clones.map(clone),
    ...result.rejected.map((bad) => row('is-rejected', bad.name, `not read — ${bad.why}`, '')),
  );
  el('open-all-btn').hidden = !result.clones.some((c) => c.differences.length);
  el('open-all-btn').textContent = 'Open all';

  renderOverview();
  renderViewer();
}

/** One clone's line: name, what it is, and the numbers, in that order. */
function row(kind, name, what, meta) {
  const line = document.createElement('div');
  line.className = `verdict ${kind}`;
  line.append(
    span('v-name', name),
    // the name is already the chip beside this; repeating it here is how a
    // line ends up reading "clone3 clone3 has no mismatches"
    span('v-what', what),
    span('v-meta', meta),
  );
  return line;
}

function span(className, text) {
  const node = document.createElement('span');
  node.className = className;
  node.textContent = text;
  return node;
}

function clone(c) {
  const kind = c.clean ? 'is-clean' : c.placed ? 'is-changed' : 'is-unplaced';
  const what = c.placed ? `has ${c.verdict}` : c.verdict;

  let meta = c.note;
  if (c.placed) {
    const bits = [`${c.identity}% identity`, `${c.length.toLocaleString()} bp`];
    // the read's own start and strand are reported once and then never again:
    // they are how the vendor happened to open the circle, not a finding
    if (c.strand === '-') bits.push('read on the opposite strand');
    if (c.offset) bits.push(`opens at ${(c.offset + 1).toLocaleString()}`);
    if (c.length_difference) {
      bits.push(`${c.length_difference > 0 ? '+' : ''}${c.length_difference.toLocaleString()} bp`);
    }
    meta = bits.join(' · ');
  }

  // A clone with nothing to show stays a plain line. Wrapping it in a
  // disclosure that opens onto nothing is a control that looks live and is
  // not, and "no mismatches" is the whole answer for that clone anyway.
  if (!c.differences.length) return row(kind, c.name, what, meta);

  const box = document.createElement('details');
  box.className = `verdict ${kind}`;
  const head = document.createElement('summary');
  head.append(span('v-name', c.name), span('v-what', what), span('v-meta', meta));
  box.append(head);

  const list = document.createElement('ul');
  list.className = 'v-diffs';
  for (const d of c.differences) {
    const item = document.createElement('li');
    item.className = `d-${d.kind}`;
    const button = document.createElement('button');
    button.type = 'button';
    button.className = 'd-jump';
    button.textContent = d.text;
    button.addEventListener('click', () => jumpTo(d.start));
    item.append(button);
    list.append(item);
  }
  box.append(list);
  return box;
}

/* -------------------------------------------------------------- overview */

function renderOverview() {
  const result = state.result;
  const length = result.reference.length;
  const box = el('overview');
  box.replaceChildren();

  const lane = document.createElement('div');
  lane.className = 'ov-features';
  for (const region of result.reference.regions) {
    const block = document.createElement('span');
    block.className = 'ov-feature';
    const start = Math.max(0, region.start);
    const end = Math.min(length, region.end);
    block.style.left = `${(100 * start) / length}%`;
    block.style.width = `${Math.max(0.15, (100 * (end - start)) / length)}%`;
    block.title = `${region.label} · ${(start + 1).toLocaleString()}–${end.toLocaleString()}`;
    block.textContent = region.label;
    lane.append(block);
  }
  box.append(lane);

  for (const c of result.clones) {
    const row = document.createElement('div');
    row.className = 'ov-row' + (c.clean ? ' is-clean' : c.placed ? '' : ' is-unplaced');

    const name = document.createElement('span');
    name.className = 'ov-name';
    name.textContent = c.name;

    const track = document.createElement('span');
    track.className = 'ov-track';
    for (const d of c.differences) {
      const tick = document.createElement('button');
      tick.type = 'button';
      tick.className = `ov-tick d-${d.kind}`;
      tick.style.left = `${(100 * d.start) / length}%`;
      tick.title = `${c.name}: ${d.text}`;
      tick.addEventListener('click', () => jumpTo(d.start));
      track.append(tick);
    }
    if (!c.placed) {
      const none = document.createElement('span');
      none.className = 'ov-none';
      none.textContent = 'not this plasmid';
      track.append(none);
    }
    row.append(name, track);
    box.append(row);
  }
}

/* ---------------------------------------------------------------- viewer */

const RULER = 10;    // a tick label every this many columns

function escape(text) {
  return text.replace(/&/g, '&amp;').replace(/</g, '&lt;');
}

function seqHtml(refRow, row, from, to) {
  return placed(refRow, row, from, to, CH)
    .map((run) => `<span class="c-${run.kind}" style="left:${run.left}px">${escape(run.text)}</span>`)
    .join('');
}

function rowElement(name, inner, width, extra = '') {
  const row = document.createElement('div');
  row.className = `aln-row ${extra}`;
  row.innerHTML =
    `<span class="aln-name"></span><span class="aln-seq" style="width:${width}px">${inner}</span>`;
  row.querySelector('.aln-name').textContent = name;
  return row;
}

function rulerHtml(from, to) {
  return ticks(state.result.refOf, from, to, RULER, CH)
    .map((tick) => `<span class="aln-tick" style="left:${tick.left}px">${tick.label.toLocaleString()}</span>`)
    .join('');
}

/* One alignment, end to end, scrolled through.
 *
 * There is deliberately no second view. An earlier draft also offered short
 * windows cut around each difference, and it made the screen harder to read
 * rather than easier: two layouts of the same data, and the one you were
 * looking at was never obviously the one you wanted. The differences are
 * already findable three other ways - the ticks on the overview, the buttons
 * in each clone's row, and Next/Previous here - all of which scroll this one
 * alignment rather than replacing it with a different shape.
 */
function renderViewer() {
  const viewer = el('viewer');
  viewer.replaceChildren();
  viewer.onscroll = null;

  const result = state.result;
  const alignment = result.alignment;

  if (!alignment.rows.length) {
    const empty = document.createElement('p');
    empty.className = 'aln-empty';
    empty.textContent = 'No clone could be placed on this reference, so there is nothing to align.';
    viewer.append(empty);
    el('viewer-note').textContent = '';
    return;
  }

  const width = alignment.columns * CH;
  const differences = result.events.length;
  el('viewer-note').textContent =
    `${alignment.columns.toLocaleString()} columns · `
    + (differences
        ? `${differences} difference${differences === 1 ? '' : 's'}, marked in the rows`
        : 'no differences anywhere');

  const scroll = document.createElement('div');
  scroll.className = 'aln-scroll';

  const cursor = document.createElement('div');
  cursor.className = 'aln-cursor';
  cursor.hidden = true;
  scroll.append(cursor);

  const rows = [
    rowElement('', '', width, 'aln-ruler'),
    rowElement(result.reference.name, '', width, 'is-ref'),
    ...alignment.rows.map(({ name }) => rowElement(name, '', width)),
  ];
  for (const row of rows) scroll.append(row);
  viewer.append(scroll);

  // Only the columns on screen are ever in the DOM. A 12 kb construct across
  // eight clones is over a hundred thousand cells; rendering them all is what
  // turns a viewer into a frozen tab.
  let drawn = [-1, -1];
  const paint = () => {
    const first = Math.max(0, Math.floor(viewer.scrollLeft / CH) - 200);
    const last = Math.min(alignment.columns, first + Math.ceil(viewer.clientWidth / CH) + 400);
    if (first >= drawn[0] && last <= drawn[1]) return;
    drawn = [first, last];
    const sequences = rows.map((row) => row.querySelector('.aln-seq'));
    sequences[0].innerHTML = rulerHtml(first, last);
    sequences[1].innerHTML = seqHtml(alignment.reference, alignment.reference, first, last);
    alignment.rows.forEach(({ row }, n) => {
      sequences[n + 2].innerHTML = seqHtml(alignment.reference, row, first, last);
    });
  };
  paint();
  viewer.onscroll = paint;
  // and on resize: a wider window needs more columns than the last paint drew,
  // and nothing else would ask for them until the next scroll
  state.repaint = paint;
  if (state.jumpAt >= 0) jumpTo(state.jumpAt);
}

/* ----------------------------------------------------------------- jumps */

/** Scroll the alignment to a reference base and mark the column it landed on.
 *
 * The mark matters: scrolling to the middle of ten thousand identical-looking
 * columns and leaving the reader to spot which one was meant is how a jump
 * button ends up feeling broken.
 */
function jumpTo(at) {
  state.jumpAt = at;
  const { colOf, events, alignment, reference } = state.result;
  const [from, to] = columnSpan(at, events, colOf, alignment.insertions, reference.length);
  const viewer = el('viewer');
  viewer.scrollLeft = Math.max(0, from * CH - viewer.clientWidth / 2);

  const cursor = viewer.querySelector('.aln-cursor');
  if (!cursor) return;
  cursor.hidden = false;
  cursor.style.left = `calc(var(--aln-name) + ${from * CH}px)`;
  cursor.style.width = `${(to - from) * CH}px`;
}

function step(direction) {
  const events = state.result && state.result.events;
  if (!events || !events.length) return;
  let next;
  if (state.jumpAt < 0) {
    next = direction > 0 ? events[0] : events[events.length - 1];
  } else if (direction > 0) {
    next = events.find((d) => d.start > state.jumpAt) || events[0];
  } else {
    next = [...events].reverse().find((d) => d.start < state.jumpAt) || events[events.length - 1];
  }
  jumpTo(next.start);
}

/* ------------------------------------------------------------------ chrome */

function setStatus(text, bad = false) {
  const status = el('status');
  status.textContent = text;
  status.classList.toggle('is-bad', Boolean(bad) && Boolean(text));
}

function dropTarget(node, onFiles) {
  for (const type of ['dragenter', 'dragover']) {
    node.addEventListener(type, (event) => {
      event.preventDefault();
      node.classList.add('is-over');
    });
  }
  for (const type of ['dragleave', 'drop']) {
    node.addEventListener(type, (event) => {
      event.preventDefault();
      node.classList.remove('is-over');
      if (type === 'drop' && event.dataTransfer) onFiles(event.dataTransfer.files);
    });
  }
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

el('src-library').addEventListener('click', () => setSource('library'));
el('src-expected').addEventListener('click', () => setSource('expected'));
el('src-upload').addEventListener('click', () => setSource('upload'));

el('reference').addEventListener('change', () => {
  state.chosen = el('reference').value;
  updateRunButton();
});

el('ref-file').addEventListener('change', async (event) => {
  const [file] = event.target.files || [];
  if (!file) return;
  state.referenceFile = await readFile(file);
  el('ref-file-name').textContent = file.name;
  el('ref-note').textContent = file.name;
  updateRunButton();
});

el('clone-files').addEventListener('change', (event) => addClones(event.target.files || []));
el('paste-add-btn').addEventListener('click', addPasted);

el('ref-paste-btn').addEventListener('click', () => {
  const text = el('ref-paste-text').value;
  if (!text.trim()) {
    el('ref-paste-note').textContent = 'nothing to use';
    return;
  }
  const name = nameForPaste(text, '');
  state.referenceFile = { name, text };
  el('ref-paste-note').textContent = `using ${name}`;
  el('ref-file-name').textContent = name;
  el('ref-note').textContent = name;
  updateRunButton();
});

el('clear-btn').addEventListener('click', () => {
  state.clones = [];
  state.pasted = 0;
  el('clone-files').value = '';
  el('paste-text').value = '';
  el('paste-name').value = '';
  el('paste-note').textContent = '';
  renderFiles();
});
el('run-btn').addEventListener('click', run);

el('open-all-btn').addEventListener('click', () => {
  const boxes = [...el('verdicts').querySelectorAll('details')];
  const opening = boxes.some((box) => !box.open);
  for (const box of boxes) box.open = opening;
  el('open-all-btn').textContent = opening ? 'Close all' : 'Open all';
});
el('prev-btn').addEventListener('click', () => step(-1));
el('next-btn').addEventListener('click', () => step(1));

dropTarget(el('clone-drop'), addClones);
dropTarget(el('ref-drop'), async (files) => {
  const [file] = files;
  if (!file) return;
  state.referenceFile = await readFile(file);
  el('ref-file-name').textContent = file.name;
  el('ref-note').textContent = file.name;
  updateRunButton();
});

window.addEventListener('resize', () => { if (state.repaint) state.repaint(); });

wireFinder();
measure();
loadFolder();
loadReferences();
renderFiles();
