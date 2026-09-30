/* The Proteins screen.
 *
 * The one screen that uses the network. Everything else in this app reads a
 * local folder and would work on a plane; this asks UniProt for a protein and
 * writes it back as DNA that can be ordered as a part.
 *
 * It talks only to /api/uniprot/* and the shared /api/library/summary, and
 * nothing imports it.
 */

const el = (id) => document.getElementById(id);
let made = null;

function say(text, bad = false) {
  const node = el('status');
  node.textContent = text;
  node.classList.toggle('is-bad', bad);
}

async function getJSON(path, options) {
  const response = await fetch(path, options);
  if (!response.ok) {
    let detail = `${response.status}`;
    try {
      detail = (await response.json()).detail || detail;
    } catch {
      // a failure with no JSON body is still a failure worth naming
    }
    throw new Error(detail);
  }
  return response.json();
}

function resultRow(entry) {
  const row = document.createElement('button');
  row.type = 'button';
  row.className = 'result';

  const badge = document.createElement('span');
  badge.className = entry.reviewed ? 'pill is-reviewed' : 'pill';
  badge.textContent = entry.reviewed ? 'reviewed' : 'unreviewed';
  badge.title = entry.reviewed
    ? 'curated by UniProt'
    : 'automatically annotated, not curated';

  const what = document.createElement('span');
  what.className = 'result-what';
  const name = document.createElement('span');
  name.className = 'result-name';
  name.textContent = entry.label;
  const meta = document.createElement('span');
  meta.className = 'result-meta';
  meta.textContent = `${entry.accession} · ${entry.length} aa · ${entry.organism}`;
  what.append(name, meta);

  row.append(badge, what);
  row.addEventListener('click', () => build(entry));
  return row;
}

async function search() {
  const query = el('q').value.trim();
  if (!query) return;
  say('Searching UniProt…');
  el('results').replaceChildren();
  try {
    const found = await getJSON(
      `/api/uniprot/search?q=${encodeURIComponent(query)}`
      + `&organism=${encodeURIComponent(el('organism').value)}`,
    );
    say(found.count
      ? `${found.count} result${found.count === 1 ? '' : 's'} — pick one to write it as DNA`
      : `Nothing in UniProt matches “${query}”.`);
    el('results').replaceChildren(...found.results.map(resultRow));
  } catch (error) {
    // being offline is ordinary for a tool that otherwise never asks
    say(`Could not search UniProt: ${error.message}`, true);
  }
}

function summaryRow(label, value, note) {
  const row = document.createElement('div');
  row.className = 'sum-row';
  const key = document.createElement('span');
  key.textContent = label;
  const val = document.createElement('span');
  val.className = 'mono';
  val.textContent = value;
  row.append(key, val);
  if (note) {
    const extra = document.createElement('span');
    extra.className = 'card-note';
    extra.textContent = note;
    row.append(extra);
  }
  return row;
}

async function build(entry, protein = null) {
  say(`Writing ${entry.accession || entry.label} as DNA…`);
  try {
    made = await getJSON('/api/uniprot/part', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        accession: protein ? '' : entry.accession,
        protein: protein || '',
        part_type: el('part-type').value,
        codon_table: el('codon-table').value,
        avoid: chosenEnzymes(),
      }),
    });
  } catch (error) {
    // the server names a bad enzyme or site precisely; showing it beside the
    // field that caused it beats a status line at the top of the page
    if (extraSites().length) el('extra-note').textContent = error.message;
    say(`Could not write ${entry.accession || entry.label}: ${error.message}`, true);
    return;
  }

  say(`${entry.label} — ${made.length.toLocaleString()} bp`);
  el('part-title').textContent = `${made.name} · coding sequence`;
  el('part-card').hidden = false;

  const box = el('part-summary');
  box.replaceChildren(
    summaryRow('Protein', `${made.protein_length} aa`, entry.organism),
    summaryRow('Coding sequence', `${made.length.toLocaleString()} bp`, made.part_type_label),
    summaryRow('Stop codon', made.stop ? 'yes' : 'no', made.stop_reason),
    summaryRow('Codon usage', made.codon_table_name, `GC ${made.gc}%`),
    summaryRow('Sites avoided', made.avoided.join(', ') || 'none',
      made.avoided.length
        ? 'none of these appears in the sequence'
        : 'nothing was kept out - check before ordering'),
    summaryRow('Codons passed over', String(made.compromised.length),
      made.compromised.length
        ? 'the preferred codon would have spelled a site here'
        : 'every residue got its preferred codon'),
    summaryRow('Translates back', made.verified ? 'yes' : 'NO',
      made.verified ? 'the DNA spells the protein exactly' : 'do not use this'),
    ...feasibilityRows(made),
  );

  el('dna').textContent = (made.dna.match(/.{1,60}/g) || []).join('\n');
}

/* The two choices that change the sequence: which organism's codon usage to
 * write it for, and which sites to keep out of it. Both are offered rather
 * than fixed - a part bound for E. coli wants E. coli codons, and an enzyme
 * the part will never meet is not worth a run of rare codons to avoid. */
async function loadOptions() {
  let options;
  try {
    options = await getJSON('/api/uniprot/options');
  } catch {
    return; // the page still works on the defaults
  }

  const targets = el('part-type');
  targets.replaceChildren();
  for (const target of options.part_types || []) {
    const option = document.createElement('option');
    option.value = target.key;
    option.textContent = target.label;
    option.dataset.why = target.why;
    if (target.key === '3') option.selected = true;   // the usual reason to be here
    targets.append(option);
  }
  showWhy();

  const tables = el('codon-table');
  tables.replaceChildren();
  for (const table of options.codon_tables) {
    const option = document.createElement('option');
    option.value = table.key;
    option.textContent = table.name;
    if (table.key === 'scerevisiae') option.selected = true;
    tables.append(option);
  }

  const host = el('enzymes');
  host.replaceChildren();
  for (const enzyme of [...options.enzymes, ...options.other_enzymes]) {
    const label = document.createElement('label');
    label.className = 'enzyme';
    const box = document.createElement('input');
    box.type = 'checkbox';
    box.value = enzyme.name;
    box.checked = enzyme.default;
    box.addEventListener('change', rebuild);
    const text = document.createElement('span');
    text.textContent = enzyme.name;
    const site = document.createElement('span');
    site.className = 'enzyme-site';
    site.textContent = enzyme.site;
    label.append(box, text, site);
    host.append(label);
  }
}

/* The part type decides the stop codon, so the reason is shown beside it
 * rather than discovered afterwards in the summary. */
function showWhy() {
  const chosen = el('part-type').selectedOptions[0];
  el('part-type-why').textContent = chosen ? chosen.dataset.why || '' : '';
}

function chosenEnzymes() {
  const ticked = [...document.querySelectorAll('#enzymes input:checked')].map((b) => b.value);
  return [...ticked, ...extraSites()];
}

/* Anything the person typed: a name the app knows, or the site itself. The
 * server decides which - it owns the enzyme list and the IUPAC alphabet, and
 * two places deciding what counts as a site is one place too many. */
function extraSites() {
  return (el('extra-sites').value || '')
    .split(/[,;\s]+/)
    .map((s) => s.trim())
    .filter(Boolean);
}

/* Will a vendor make it? None of this shows up in a translation check, and
 * the strategy that writes these genes is the one that produces it: always
 * taking the most preferred codon maps a repeated protein motif straight onto
 * repeated bases, and a tract of one residue onto a tract of one codon.
 *
 * Every row carries its own proof - a position, a base, a count - because a
 * verdict you cannot check is an assertion. */
function feasibilityRows(made) {
  const f = made.feasibility;
  if (!f) return [];
  const rows = [];

  const h = f.homopolymer;
  rows.push(verdictRow('Longest single-base run', h.ok,
    h.length ? `${h.length} \u00d7 ${h.base}` : 'none',
    h.ok
      ? `under the ${h.limit} that vendors query`
      : `${h.length} ${h.base}s from base ${(h.at + 1).toLocaleString()} \u2014 at or over ${h.limit}`));

  const r = f.repeat;
  rows.push(verdictRow('Longest exact repeat', r.ok,
    r.length ? `${r.length} bp` : 'none',
    r.ok
      ? `nothing repeated at ${r.limit} bp or more`
      : `identical at ${(r.first + 1).toLocaleString()} and ${(r.second + 1).toLocaleString()}`));

  const g = f.gc;
  rows.push(verdictRow('GC across ' + g.window + ' bp windows', g.ok,
    `${g.min}\u2013${g.max}%`,
    g.ok
      ? `every window inside ${g.low}\u2013${g.high}%`
      : `${g.outside.toLocaleString()} window${g.outside === 1 ? '' : 's'} outside ${g.low}\u2013${g.high}%`));

  return rows;
}

function verdictRow(label, ok, value, note) {
  const row = summaryRow(label, value, note);
  row.classList.add(ok ? 'is-ok' : 'is-flagged');
  return row;
}

/* Applied on a button rather than as you type: half a typed enzyme name is a
 * 422 from the server, and rebuilding on every keystroke would spend the whole
 * time showing errors about sequences nobody asked for. */
async function applyExtraSites() {
  const sites = extraSites();
  el('extra-note').textContent = sites.length
    ? `also keeping out ${sites.join(', ')}`
    : '\u00a0';
  if (made) await rebuild();
}

/** Rewrite the sequence under the current options, keeping the same protein. */
async function rebuild() {
  if (!made) return;
  await build(made.source || { accession: made.accession, label: made.name },
              made.source ? null : made.protein);
}

function wire() {
  el('go').addEventListener('click', search);
  el('codon-table').addEventListener('change', rebuild);
  el('part-type').addEventListener('change', () => { showWhy(); rebuild(); });
  el('extra-apply').addEventListener('click', applyExtraSites);
  el('extra-sites').addEventListener('keydown', (event) => {
    if (event.key === 'Enter') applyExtraSites();
  });
  el('q').addEventListener('keydown', (event) => {
    if (event.key === 'Enter') search();
  });
  el('organism').addEventListener('change', () => {
    if (el('q').value.trim()) search();
  });

  el('copy').addEventListener('click', async () => {
    if (!made) return;
    await navigator.clipboard.writeText(made.dna);
    const button = el('copy');
    button.textContent = 'Copied';
    setTimeout(() => { button.textContent = 'Copy DNA'; }, 1500);
  });

  /* Hand the sequence to the Part screen. It reads its own state from the URL
   * on load, so the sequence travels in the fragment rather than a store -
   * nothing is persisted and a reload of that screen keeps working. */
  el('to-part').addEventListener('click', () => {
    if (!made) return;
    const parameters = new URLSearchParams({
      sequence: made.dna, name: made.name, part_type: '3',
    });
    window.location.href = `/part#${parameters.toString()}`;
  });
}

async function loadSummary() {
  try {
    const summary = await (await fetch('/api/library/summary')).json();
    el('folder').textContent = `${summary.parts} parts indexed`;
  } catch {
    el('folder').textContent = '';
  }
}

wire();
loadOptions();
loadSummary();
