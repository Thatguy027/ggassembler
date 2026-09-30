/* Level 1 - part plasmid construction.
 *
 * Talks only to /api/level1/* and the shared /api/library/*. Imports nothing
 * from a sibling level.
 *
 * The screen holds the request; the server does every piece of reasoning,
 * including the one that matters most - feeding the predicted plasmid back
 * through the detector to prove it reads as the type that was asked for.
 */

const state = {
  part_type: '3',
  template: '',
  sequence: '',
  start: 1,
  end: 0,
  entry_vector: null,
  mode: 'pcr',
  name: 'new_part',
  conventions: { gly_ser_linker: true, strip_stop: true, stop_and_xhoi: true },
  domesticate: false,
};

let enzymes = { part: 'BsaI', multigene: 'BsmBI', linearizer: 'NotI' };
let entryOverhangs = ['TCGG', 'GACC'];

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

// ---------------------------------------------------------------- primers ---

const ROLE_LABELS = {
  pad: 'pad',
  multigene_enzyme: () => enzymes.multigene,
  entry_overhang: () => entryOverhangs.join('/') + ' entry-vector OH',
  part_enzyme: () => `${enzymes.part}, part-side`,
  type_flank: 'type flanks',
};

function renderPrimers(design) {
  const box = el('primers');
  const primers = design.fragments.flatMap((f) => [f.forward, f.reverse]).filter(Boolean);
  if (!primers.length) {
    box.replaceChildren(note('no primers for this insert source'));
    return;
  }

  box.replaceChildren(
    ...primers.map((primer) => {
      const wrap = document.createElement('div');
      wrap.className = 'primer';

      const head = document.createElement('div');
      head.className = 'primer-head';
      const name = document.createElement('span');
      name.className = 'primer-name';
      name.textContent = primer.name;
      const meta = document.createElement('span');
      meta.className = 'primer-meta';
      meta.textContent = `${primer.length} nt · Tm ${primer.tm} °C`;
      head.append(name, meta);

      const seq = document.createElement('div');
      seq.className = 'primer-seq';
      let consumed = 0;
      for (const segment of primer.segments) {
        const index = primer.sequence.indexOf(segment.seq, consumed);
        if (index < 0) continue;
        if (index > consumed) {
          seq.append(document.createTextNode(primer.sequence.slice(consumed, index)));
        }
        const span = document.createElement('span');
        span.className = `seg-${segment.role}`;
        span.textContent = segment.seq;
        span.title = typeof ROLE_LABELS[segment.role] === 'function'
          ? ROLE_LABELS[segment.role]()
          : ROLE_LABELS[segment.role] || segment.role;
        seq.append(span);
        consumed = index + segment.seq.length;
      }
      if (consumed < primer.sequence.length) {
        seq.append(document.createTextNode(primer.sequence.slice(consumed)));
      }

      wrap.append(head, seq);
      return wrap;
    }),
  );

  el('legend-keys').replaceChildren(
    ...[
      ['multigene_enzyme', enzymes.multigene, enzymes.multigene],
      ['entry_overhang', entryOverhangs.join('/'), 'entry-vector OH'],
      ['part_enzyme', `${enzymes.part} site`, `${enzymes.part}, part-side`],
      ['type_flank', `${design.five_prime}/${design.three_prime}`, `Type ${design.part_type} flanks`],
    ].map(([role, code, label]) => {
      const key = document.createElement('span');
      key.className = 'legend-key';
      const mono = document.createElement('span');
      mono.className = `mono seg-${role}`;
      mono.textContent = code;
      key.append(mono, document.createTextNode(' ' + label));
      return key;
    }),
  );
}

// ------------------------------------------------- primers on the template ---

const COMPLEMENT = { A: 'T', C: 'G', G: 'C', T: 'A', N: 'N' };
const complement = (seq) => [...seq].map((b) => COMPLEMENT[b] || b).join('');
const reverse = (seq) => [...seq].reverse().join('');

function span(text, className) {
  const node = document.createElement('span');
  if (className) node.className = className;
  node.textContent = text;
  return node;
}

/** Each primer over the stretch of template it anneals to, tail hanging off. */
function renderBinding(design) {
  const box = el('binding');
  const binding = design.binding || [];
  if (!binding.length) {
    box.replaceChildren(note('nothing designed yet'));
    return;
  }

  const rows = [];
  for (const primer of binding) {
    const row = document.createElement('div');
    row.className = 'binding-row';

    const label = document.createElement('div');
    label.className = 'binding-label';
    label.append(
      span(primer.name),
      span(`${primer.length} nt · Tm ${primer.tm} °C · ${primer.annealing.length} nt anneals`),
    );

    const strand = document.createElement('div');
    strand.className = 'binding-strand';
    if (primer.direction === 'forward') {
      // tail flares off the 5' end, annealing sits on the template
      strand.append(
        span("5'-", 'binding-tick'),
        span(primer.tail, 'binding-tail'),
        span(primer.annealing, 'binding-anneal'),
        span("-3'", 'binding-tick'),
      );
      const template = document.createElement('div');
      template.className = 'binding-strand binding-template';
      template.append(
        span(' '.repeat(3 + primer.tail.length) + "3'-", 'binding-tick'),
        span(complement(primer.annealing)),
        span('···  template', 'binding-tick'),
      );
      row.append(label, strand, template);
    } else {
      // the reverse primer reads along the bottom strand, so it is drawn under
      // the template and right to left
      const template = document.createElement('div');
      template.className = 'binding-strand binding-template';
      template.append(
        span("  template  ···5'-", 'binding-tick'),
        span(reverse(complement(primer.annealing))),
        span("-3'", 'binding-tick'),
      );
      strand.append(
        span(' '.repeat(16) + "3'-", 'binding-tick'),
        span(reverse(primer.annealing), 'binding-anneal'),
        span(reverse(primer.tail), 'binding-tail'),
        span("-5'", 'binding-tick'),
      );
      row.append(label, template, strand);
    }
    rows.push(row);
  }

  const key = document.createElement('div');
  key.className = 'binding-label';
  key.append(span('tail added by the reaction', 'binding-tail'),
             span('anneals to the template', 'binding-anneal'));
  rows.push(key);
  box.replaceChildren(...rows);
}

// ------------------------------------------------- the construct going in ---

const SEG_COLORS = {
  pad: 'var(--badge-5-bg)',
  multigene_enzyme: 'var(--badge-2-bg)',
  entry_overhang: 'var(--badge-6-bg)',
  part_enzyme: 'var(--badge-1-bg)',
  type_flank: 'var(--badge-4-bg)',
  annealing: 'var(--chip)',
  body: 'var(--part-3)',
};

function renderConstruct(design) {
  const box = el('construct');
  const first = (design.fragments || [])[0];
  if (!first || !first.forward) {
    box.replaceChildren(note('nothing designed yet'));
    return;
  }

  // the amplicon, read outwards from the middle: forward tail, body, mirrored
  const parts = [];
  for (const segment of first.forward.segments) {
    if (segment.role === 'annealing') continue;
    parts.push({ role: segment.role, seq: segment.seq, side: "5'" });
  }
  parts.push({ role: 'body', seq: `${design.body_length} bp insert`, side: '' });
  const tailSegments = first.reverse ? [...first.reverse.segments] : [];
  for (const segment of tailSegments.reverse()) {
    if (segment.role === 'annealing') continue;
    parts.push({ role: segment.role, seq: segment.seq, side: "3'" });
  }

  const bar = document.createElement('div');
  bar.className = 'construct-bar';
  for (const part of parts) {
    const seg = document.createElement('div');
    seg.className = 'construct-seg';
    seg.style.background = SEG_COLORS[part.role] || 'var(--chip)';
    seg.style.flexGrow = String(part.role === 'body' ? 8 : Math.max(part.seq.length, 4));
    seg.textContent = part.role === 'body' ? part.seq : part.seq;
    seg.title = `${part.role.replace(/_/g, ' ')}${part.side ? ' (' + part.side + ' end)' : ''}: ${part.seq}`;
    bar.append(seg);
  }

  const key = document.createElement('div');
  key.className = 'construct-key';
  const labels = {
    pad: 'pad',
    multigene_enzyme: `${design.cloning_enzyme} (cuts it into the vector)`,
    entry_overhang: 'end the vector accepts',
    part_enzyme: `${enzymes.part} (releases the part later)`,
    type_flank: `type ${design.part_type} flank`,
    body: 'your sequence',
  };
  for (const [role, label] of Object.entries(labels)) {
    if (!parts.some((p) => p.role === role)) continue;
    const item = document.createElement('span');
    const swatch = document.createElement('span');
    swatch.className = 'construct-swatch';
    swatch.style.background = SEG_COLORS[role];
    item.append(swatch, document.createTextNode(label));
    key.append(item);
  }
  box.replaceChildren(bar, key);
}

// ------------------------------------------------------------ the ligation ---

/** One junction as an annealed duplex, with the sticky end in the middle. */
function junction(title, leftSeq, leftClass, overhang, rightSeq, rightClass) {
  const wrap = document.createElement('div');
  wrap.className = 'junction';

  const heading = document.createElement('div');
  heading.className = 'junction-title';
  heading.textContent = title;

  const top = document.createElement('div');
  top.className = 'junction-duplex';
  top.append(
    span("5'-", 'junction-nick'),
    span(leftSeq, leftClass),
    span(overhang, 'junction-oh'),
    span(rightSeq, rightClass),
    span("-3'", 'junction-nick'),
  );

  const bottom = document.createElement('div');
  bottom.className = 'junction-duplex';
  bottom.append(
    span("3'-", 'junction-nick'),
    span(complement(leftSeq), leftClass),
    span(complement(overhang), 'junction-oh'),
    span(complement(rightSeq), rightClass),
    span("-5'", 'junction-nick'),
  );

  wrap.append(heading, top, bottom);
  return wrap;
}

function renderLigation(design) {
  const box = el('ligation');
  const view = design.ligation || {};
  if (!view.vector || !view.insert) {
    box.replaceChildren(note('nothing designed yet'));
    return;
  }
  const { vector, insert } = view;

  const summary = document.createElement('div');
  summary.className = 'junction-title';
  summary.textContent =
    `${vector.name} cut with ${vector.enzyme}: ${vector.dropout_length.toLocaleString()} bp ` +
    `${vector.dropout_component ? 'dropout (' + vector.dropout_component + ')' : 'dropout'} out, ` +
    `${vector.backbone_length.toLocaleString()} bp backbone kept, ` +
    `${insert.length.toLocaleString()} bp insert in.`;

  box.replaceChildren(
    summary,
    junction(
      `vector 3' end → insert 5' end · ${insert.left_overhang}`,
      vector.right_flank, 'junction-vector',
      insert.left_overhang,
      insert.left_flank, 'junction-insert',
    ),
    junction(
      `insert 3' end → vector 5' end · ${insert.right_overhang}`,
      insert.right_flank, 'junction-insert',
      insert.right_overhang,
      vector.left_flank, 'junction-vector',
    ),
  );
}

// ---------------------------------------------------------- domestication ---

function scanRow(ok, text) {
  const row = document.createElement('div');
  row.className = `scan-row ${ok ? 'scan-ok' : 'scan-warn'}`;
  const mark = document.createElement('span');
  mark.className = 'scan-mark';
  mark.textContent = ok ? '✓' : '!';
  const label = document.createElement('span');
  label.textContent = text;
  row.append(mark, label);
  return row;
}

function renderDomestication(design) {
  const rows = [];
  const internal = design.issues.filter((i) => i.code === 'internal_site');
  const mutations = design.fragments.flatMap((f) => f.mutations);

  if (internal.length) {
    for (const issue of internal) {
      const extra = design.fragments.length > 1
        ? ` — split into ${design.fragments.length} fragments${mutations.length ? ', mutations proposed' : ''}`
        : '';
      rows.push(scanRow(false, issue.message + extra));
    }
  } else {
    rows.push(scanRow(true, `No internal ${enzymes.part}, ${enzymes.multigene} or ${enzymes.linearizer} sites`));
  }

  for (const mutation of mutations) {
    rows.push(scanRow(mutation.silent, mutation.description));
  }
  el('domestication').replaceChildren(...rows);

  const junctions = design.fragments.slice(0, -1).map((f) => f.right_overhang);
  el('junction-note').textContent = junctions.length
    ? `Internal junction overhangs ${junctions.join(', ')} chosen to avoid ${entryOverhangs.join(' / ')} and any 3-of-4 match.`
    : '';
}

// --------------------------------------------------------------- results ---

function note(text) {
  const span = document.createElement('span');
  span.className = 'card-note';
  span.textContent = text;
  return span;
}

function renderResults(design) {
  const verdict = el('verdict');
  if (design.product_length) {
    const line = document.createElement('div');
    line.className = 'verdict-line';
    const bp = document.createElement('span');
    bp.className = 'verdict-bp';
    bp.textContent = `${design.product_length.toLocaleString()} bp`;
    const pill = document.createElement('span');
    const good = design.validated_as === design.part_type;
    pill.className = good ? 'pill' : 'pill pill-error';
    pill.textContent = good
      ? `re-digests as type ${design.validated_as}`
      : `reads back as ${design.validated_as || 'nothing'}`;
    line.append(bp, pill);
    verdict.replaceChildren(line, note(`in ${design.entry_vector} · ${design.five_prime} → ${design.three_prime}`));
  } else {
    verdict.replaceChildren(note('nothing designed yet'));
  }

  el('issues').replaceChildren(
    ...design.issues.map((issue) => {
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

  el('f-enzyme').textContent = design.cloning_enzyme || enzymes.multigene;
  el('f-selection').textContent = design.destination_marker
    ? `${design.destination_marker} · ccdB counter-selection`
    : 'ccdB counter-selection';
  el('f-fragments').textContent = design.fragments.length
    ? `${design.fragments.length} + vector`
    : '—';
  el('export-btn').disabled = !design.product_length;
  el('save-btn').disabled = !design.product_length;

  const alt = el('alt-field');
  if (design.mode === 'gblock') {
    alt.hidden = false;
    el('alt').value = design.gblock;
  } else if (design.mode === 'oligo') {
    alt.hidden = false;
    el('alt').value = design.oligos.map((o) => `>${o.name}\n${o.sequence}`).join('\n');
  } else {
    alt.hidden = true;
  }
}

// ----------------------------------------------------------------- wiring ---

let pending = null;

async function refresh() {
  read();
  const mine = {};
  pending = mine;
  const design = await post('/api/level1/design', state);
  if (pending !== mine) return;
  renderPrimers(design);
  renderBinding(design);
  renderConstruct(design);
  renderLigation(design);
  renderDomestication(design);
  renderResults(design);
}

// ------------------------------------------------------------------ region ---

/* Where the amplicon starts and ends, read off the template's own annotations.
 *
 * The screen used to ask for the coordinates as two numbers, which meant
 * opening the .gbk in another tool to look them up - for a screen whose whole
 * job is to save you that trip. Every file here already carries its CDS
 * features, so the range can come from the file: pick the gene, get the
 * numbers. The boxes stay editable, because a range is sometimes not a
 * feature.
 */

/** Feature types worth offering as an amplicon, in the order they read. */
const REGION_TYPES = ['CDS', 'gene', 'promoter', 'terminator', 'misc_feature', 'mRNA'];

async function loadRegions(name) {
  const field = el('region-field');
  const select = el('region');
  select.replaceChildren(select.firstElementChild);
  if (!name) {
    field.hidden = true;
    return;
  }

  let plasmid;
  try {
    plasmid = await (await fetch(`/api/library/plasmids/${encodeURIComponent(name)}`)).json();
  } catch {
    field.hidden = true;
    return;
  }

  const usable = (plasmid.features || []).filter(
    (f) => f.label && f.end > f.start && REGION_TYPES.includes(f.type),
  );
  // the part this plasmid already contributes is the most likely answer of all
  if (plasmid.part_type && plasmid.five_prime) {
    const option = document.createElement('option');
    option.value = 'part';
    option.textContent = `The type ${plasmid.part_type} part itself`;
    select.append(option);
  }
  for (const feature of usable) {
    const option = document.createElement('option');
    option.value = `${feature.start + 1}:${feature.end}`;
    option.textContent =
      `${feature.label} · ${feature.type} · ${feature.end - feature.start} bp`;
    option.title = `${feature.start + 1}–${feature.end} on ${plasmid.name}`;
    select.append(option);
  }
  field.hidden = !usable.length && !plasmid.part_type;
  select.disabled = !usable.length && !plasmid.part_type;
}

function applyRegion() {
  const value = el('region').value;
  if (!value) {
    el('start').value = '1';
    el('end').value = '0';
  } else if (value === 'part') {
    // the server already knows where the part sits; 0/0 asks it to use that
    el('start').value = '1';
    el('end').value = '0';
  } else {
    const [start, end] = value.split(':');
    el('start').value = start;
    el('end').value = end;
  }
  refresh();
}

function read() {
  state.part_type = el('part-type').value || '3';
  state.template = el('template').value || null;
  state.sequence = el('sequence').value || null;
  state.start = Number(el('start').value) || 1;
  state.end = Number(el('end').value) || 0;
  state.entry_vector = el('entry-vector').value || null;
  state.destination = state.entry_vector;
  state.domesticate = el('c-domesticate').checked;
  state.name = el('name').value.trim() || 'new_part';
  state.conventions = {
    gly_ser_linker: el('c-glyser').checked,
    strip_stop: el('c-stop').checked,
    strip_start: el('c-start').checked,
    infer_from_sequence: el('c-infer').checked,
    stop_and_xhoi: el('c-xhoi').checked,
  };
  el('paste-field').hidden = Boolean(state.template);
}

async function load() {
  const [summary, options] = await Promise.all([
    fetch('/api/library/summary').then((r) => r.json()),
    fetch('/api/level1/options').then((r) => r.json()),
  ]);
  enzymes = options.enzymes;
  entryOverhangs = options.entry_overhangs;

  el('folder').textContent = summary.roots.length === 1
    ? summary.roots[0]
    : `${summary.roots.length} folders`;

  const types = el('part-type');
  types.replaceChildren(
    ...options.types.map((type) => {
      const option = document.createElement('option');
      option.value = type.name;
      option.textContent = `Type ${type.name} — ${type.description.split(' (')[0]}`;
      if (type.name === state.part_type) option.selected = true;
      return option;
    }),
  );

  // Destinations, grouped so the universal vectors - the ones that take any
  // part type - come first, then the vectors already cut for this type, then
  // the part plasmids that can simply be reopened.
  const vectors = el('entry-vector');
  vectors.replaceChildren();
  const auto = document.createElement('option');
  auto.value = '';
  auto.textContent = 'Best for this part type (automatic)';
  vectors.append(auto);

  const groups = [
    ['Universal — takes any part type', (d) => d.universal],
    ['Dropout vectors cut for one type', (d) => !d.universal && d.kind === 'dropout'],
    ['Part plasmids to reopen', (d) => d.kind === 'part_plasmid'],
  ];
  for (const [label, test] of groups) {
    const matching = (options.destinations || []).filter(test);
    if (!matching.length) continue;
    const group = document.createElement('optgroup');
    group.label = `${label} (${matching.length})`;
    for (const d of matching.slice(0, 150)) {
      const option = document.createElement('option');
      option.value = d.name;
      const what = d.component ? ` [${d.component}]` : '';
      option.textContent = `${d.name}${what} · ${d.cloning_enzyme} · ${d.accepts.join('→')}`;
      option.title = [
        d.kind === 'dropout' ? 'dropout vector' : 'existing part plasmid',
        `cut with ${d.cloning_enzyme}, accepts ${d.accepts.join(' → ')}`,
        d.universal
          ? 'the insert brings its own BsaI sites, so it works for any part type'
          : 'this vector supplies the BsaI sites, so the primers are much shorter',
      ].join('\n');
      group.append(option);
    }
    vectors.append(group);
  }
  if (!(options.destinations || []).length) {
    const option = document.createElement('option');
    option.value = '';
    option.textContent = 'nothing in this library can take a part';
    vectors.append(option);
  }

  const templates = el('template');
  templates.replaceChildren(templates.firstElementChild);
  for (const template of options.templates) {
    const option = document.createElement('option');
    option.value = template.name;
    option.textContent =
      `${template.display || template.name} · ${template.length} bp${template.part_type ? ', type ' + template.part_type : ''}`;
    option.title = template.component || '';
    templates.append(option);
  }

  refresh();
}

for (const id of ['part-type', 'sequence', 'start', 'end', 'name',
                  'entry-vector', 'c-glyser', 'c-stop', 'c-start', 'c-xhoi', 'c-infer', 'c-domesticate']) {
  const node = el(id);
  node.addEventListener(node.tagName === 'SELECT' || node.type === 'checkbox' ? 'change' : 'input',
    () => refresh());
}

// choosing a template reloads the regions it offers, then redesigns
el('template').addEventListener('change', async () => {
  el('region').value = '';
  el('start').value = '1';
  el('end').value = '0';
  await loadRegions(el('template').value);
  refresh();
});
el('region').addEventListener('change', applyRegion);

for (const button of document.querySelectorAll('[data-mode]')) {
  button.addEventListener('click', () => {
    state.mode = button.dataset.mode;
    for (const other of document.querySelectorAll('[data-mode]')) {
      other.setAttribute('aria-pressed', String(other === button));
    }
    refresh();
  });
}

el('export-btn').addEventListener('click', async () => {
  const response = await fetch('/api/level1/export', {
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

el('save-btn').addEventListener('click', async () => {
  const button = el('save-btn');
  button.disabled = true;
  const result = await post('/api/level1/save', state);
  button.textContent = result.ok ? 'Saved' : 'Could not save';
  setTimeout(() => { button.textContent = 'Save to library'; button.disabled = false; }, 2000);
});

load();
