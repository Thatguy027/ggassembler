/* What each Level 1 convention actually does, keyed by the checkbox it explains.
 *
 * Content, not markup. The panel that shows it is built from these entries by
 * `explainerNode` below, so there is one copy of every sentence and adding an
 * option means adding an entry rather than another block of HTML that has to
 * be kept in step with the rest.
 *
 * Sequence examples are written as plain text so they stay readable here, with
 * `[[mark|text]]` marking the bases worth telling apart:
 *
 *   added     this convention put these bases there
 *   overhang  the 4 nt junction itself
 *   supplied  another part brings these
 *
 * The blocks are rendered into a <pre>, so the alignment in the string is the
 * alignment on screen: keep the spacing when editing.
 */

export const MARKS = {
  added: 'added by this convention',
  overhang: 'the 4 nt overhang',
  supplied: 'supplied by another part',
};

export const EXPLAINERS = {
  'c-glyser': {
    title: 'Gly-Ser linker (GG before the 3′ flank)',
    body: [
      'An overhang is 4 nt, which is not a multiple of three, so a CDS that ended '
      + 'flush against ATCC would frameshift everything downstream. Adding two bases '
      + 'before the overhang makes the junction 6 nt — two whole codons — and the '
      + 'reading frame survives into whatever comes next. GG is the pair the kit chose.',
    ],
    example: {
      caption: 'With the stop stripped and a Type 4 terminator downstream:',
      text:
        '...NNN   [[added|GG]] [[overhang|ATCC]]   [[supplied|TAA CTCGAG...]]\n'
        + '  CDS   added  overhang   from the Type 4 part\n'
        + '\n'
        + '...NNN | GGA | TCC | TAA\n'
        + '  Leu    Gly   Ser   stop',
      after:
        'GGATCC is a Gly-Ser linker and also a BamHI site, which is what gives the '
        + 'assembled cassette BglBrick compatibility.',
    },
    cost:
      'Those two residues are on your protein whether or not you wanted a fusion. '
      + 'With a plain Type 4 terminator you get Gly-Ser-stop — two amino acids longer '
      + 'than the native sequence. The paper calls this “relatively innocuous”, not free.',
    when: [
      'Turn it on when a Type 4a C-terminal fusion follows (tag, fluorescent protein, '
      + 'localisation signal), or when you are splitting a CDS across 3a and 3b.',
      'Leave it off when the native C-terminus matters. Then keep your own stop codon '
      + 'in the part.',
      'Requires “Strip a trailing stop codon”. With the stop still in, the junction '
      + 'reads …TAA GGATCC and the linker is never translated.',
      'On the split types: 3a appends GG before TTCT → GGTTCT, also Gly-Ser. 3b appends '
      + 'GG before ATCC, same as Type 3.',
    ],
  },

  'c-start': {
    title: 'Type 3: drop a leading ATG (TATG supplies it)',
    body: [
      'The Type 2 → Type 3 junction overhang is TATG, and its ATG is the start codon '
      + '— the promoter part’s downstream overhang doubles as the start. Coding '
      + 'sequences in this standard therefore begin at the second codon.',
    ],
    example: {
      text:
        'overhang            part begins here\n'
        + '  [[overhang|TATG]]              TCTGAAGAATCT...\n'
        + '  [[supplied|ATG]] | TCT | GAA | GAA | TCT\n'
        + '  Met   Ser   Glu   Glu   Ser',
      after:
        'If your source sequence carries its own ATG and this is off, you get ATGATG… '
        + '— an extra methionine. Frame is fine, the protein has a residue you '
        + 'didn’t ask for.',
    },
  },

  'c-stop': {
    title: 'Strip a trailing stop codon',
    body: [
      'The Type 4 part supplies TAA immediately after the ATCC overhang, so the stop '
      + 'is the terminator part’s job, not yours. The Type 3 convention omits it so '
      + 'that read-through into a Type 4 or 4a part is possible at all.',
    ],
    when: [
      'Leave it on unless you deliberately want a part that terminates on its own and '
      + 'cannot be fused — in which case also leave the Gly-Ser linker off.',
    ],
  },

  'c-xhoi': {
    title: 'Type 4: lead with TAA + CTCGAG',
    body: [
      'A Type 4 part begins with the stop codon the Type 3 part omitted, then an XhoI '
      + 'site, then the terminator sequence:',
    ],
    example: {
      text:
        '[[overhang|ATCC]]  [[added|TAA]]  [[added|CTCGAG]]  NNN...terminator\n'
        + '      stop  XhoI',
      after: 'The XhoI site is there for BglBrick compatibility.',
    },
    when: [
      'Type 4a is the exception: it carries a C-terminal coding sequence and ends with '
      + 'its own TAA + XhoI, rather than enabling read-through of its TGGC flank. If you '
      + 'are building a 4a part, this convention applies at its 3′ end, not its 5′ end.',
    ],
  },

  'c-infer': {
    title: 'Let a terminal stop codon mean “this CDS ends here”',
    body: [
      'When the pasted or selected sequence ends in TAA, TAG or TGA in frame, treat '
      + 'that as the end of the coding sequence instead of asking you for explicit '
      + 'coordinates.',
    ],
    when: [
      'Where it misfires: a sequence that happens to end in those three bases out of '
      + 'frame, or one where you pasted the 3′ UTR along with the CDS. Turn it off and '
      + 'give coordinates when you are not sure what you pasted.',
    ],
  },

  'c-domesticate': {
    title: 'Remove internal BsaI / BsmBI / BbsI / NotI sites by silent mutation',
    body: [
      'Off by default, and it should stay that way.',
      '“Silent” only means anything inside a coding sequence. In a promoter, '
      + 'terminator, origin or homology arm there is no such thing as a silent base '
      + 'change — the app cannot know which bases carry function, so it reports what it '
      + 'would change and leaves the decision to you.',
      'When on: a site inside an annotated CDS gets a synonymous codon; a site outside '
      + 'one is reported and refused, not guessed at.',
      'The alternative is to split the amplicon at the offending site and order two '
      + 'primer pairs — the fragment count on the destination panel is showing you the '
      + 'cost of that route.',
      'Related: the Proteins screen avoids these sites at codon-choice time instead, '
      + 'which is strictly better when you are ordering a synthetic gene — nothing has '
      + 'to be mutated because nothing bad was ever written.',
    ],
  },
};

/** Split `[[mark|text]]` out of an example into plain and marked runs. */
export function segments(text) {
  const out = [];
  const pattern = /\[\[([a-z]+)\|([^\]]*)\]\]/g;
  let last = 0;
  for (const match of text.matchAll(pattern)) {
    if (match.index > last) out.push({ text: text.slice(last, match.index) });
    out.push({ text: match[2], mark: match[1] });
    last = match.index + match[0].length;
  }
  if (last < text.length) out.push({ text: text.slice(last) });
  return out;
}

function paragraph(text, className = '') {
  const node = document.createElement('p');
  if (className) node.className = className;
  node.textContent = text;
  return node;
}

function heading(text) {
  const node = document.createElement('h4');
  node.className = 'ex-head';
  node.textContent = text;
  return node;
}

function exampleNode(example) {
  const box = document.createElement('div');
  box.className = 'ex-example';
  if (example.caption) box.append(paragraph(example.caption, 'ex-caption'));

  const pre = document.createElement('pre');
  pre.className = 'ex-seq';
  const used = new Set();
  for (const run of segments(example.text)) {
    if (!run.mark) {
      pre.append(run.text);
      continue;
    }
    used.add(run.mark);
    const span = document.createElement('span');
    span.className = `ex-${run.mark}`;
    span.textContent = run.text;
    pre.append(span);
  }
  box.append(pre);

  if (used.size) {
    const key = document.createElement('ul');
    key.className = 'ex-key';
    for (const mark of Object.keys(MARKS)) {
      if (!used.has(mark)) continue;
      const item = document.createElement('li');
      const swatch = document.createElement('span');
      swatch.className = `ex-swatch ex-${mark}`;
      item.append(swatch, document.createTextNode(MARKS[mark]));
      key.append(item);
    }
    box.append(key);
  }

  if (example.after) box.append(paragraph(example.after));
  return box;
}

/** One explainer as DOM, ready to drop into the shared panel. */
export function explainerNode(entry) {
  const fragment = document.createDocumentFragment();

  const title = document.createElement('h3');
  title.className = 'ex-title';
  title.id = 'explainer-title';
  title.textContent = entry.title;
  fragment.append(title);

  for (const text of entry.body || []) fragment.append(paragraph(text));
  if (entry.example) fragment.append(exampleNode(entry.example));

  if (entry.cost) {
    fragment.append(heading('What it costs'), paragraph(entry.cost));
  }
  if (entry.when && entry.when.length) {
    fragment.append(heading('In practice'));
    const list = document.createElement('ul');
    list.className = 'ex-when';
    for (const text of entry.when) {
      const item = document.createElement('li');
      item.textContent = text;
      list.append(item);
    }
    fragment.append(list);
  }
  return fragment;
}
