/* Column arithmetic for the alignment viewer. No DOM, on purpose.
 *
 * A screen column is not a reference base. Wherever any clone carries an
 * insertion the merged alignment opens extra columns that belong to no
 * reference position at all, so the two coordinate systems drift apart by
 * however many bases every clone before this point inserted.
 *
 * Get that wrong and nothing looks broken: the bases still render, the ruler
 * still counts, and every label quietly points a few bases off - which is the
 * one kind of error a viewer must not make, because the whole reason to look
 * at it is to find out *where* something is. So it lives here as arithmetic
 * and is tested under node rather than eyeballed.
 */

/** Both directions between reference positions and screen columns.
 *
 * `colOf[r]` is the column showing reference base `r`; `refOf[c]` is the base
 * column `c` shows, or -1 where the column is an insertion belonging to no
 * reference base. They are not inverses of each other and both are needed.
 */
export function index(insertions, referenceLength) {
  const colOf = new Int32Array(referenceLength + 1);
  const widths = new Map(insertions.map(([at, width]) => [at, width]));

  let column = 0;
  for (let r = 0; r < referenceLength; r += 1) {
    column += widths.get(r) || 0;
    colOf[r] = column;
    column += 1;
  }
  colOf[referenceLength] = column;

  const refOf = new Int32Array(column).fill(-1);
  for (let r = 0; r < referenceLength; r += 1) refOf[colOf[r]] = r;

  return { colOf, refOf, columns: column };
}

/** What one cell is: match, mismatch, deletion, insertion, uncovered, padding. */
export function classOf(refChar, cloneChar) {
  if (refChar === '-') return cloneChar === '-' ? 'p' : 'i';
  if (cloneChar === '.') return 'u';
  if (cloneChar === '-') return 'd';
  return cloneChar === refChar ? 'm' : 'x';
}

/** Cells of one kind, merged into runs, over the half-open column range.
 *
 * A window of 8,000 columns across six clones is 48,000 cells and almost all
 * of them match; one element per cell is what turns a viewer like this into a
 * frozen tab. Nearly all of them come out as a handful of runs instead.
 */
export function runs(referenceRow, row, from, to) {
  if (to <= from) return [];
  const out = [];
  let start = from;
  let kind = classOf(referenceRow[from], row[from]);
  for (let i = from + 1; i < to; i += 1) {
    const next = classOf(referenceRow[i], row[i]);
    if (next !== kind) {
      out.push([kind, start, row.slice(start, i)]);
      kind = next;
      start = i;
    }
  }
  out.push([kind, start, row.slice(start, to)]);
  return out;
}

/** The columns one difference occupies, as `[from, to)`.
 *
 * Not simply `colOf[start]` plus its length, and the two kinds go opposite
 * ways. A substitution or a deletion runs forward from its base, so its far
 * end is `colOf[end]` - which already counts any insertion columns falling
 * inside it. An insertion goes backwards: its bases sit in the columns opened
 * immediately *before* the reference base it precedes, so a mark that ran
 * forward from `colOf[start]` would point at the wrong bases entirely.
 */
export function span(at, events, colOf, insertions, referenceLength) {
  let from = colOf[at];
  let to = colOf[at] + 1;
  for (const event of events) {
    if (event.start !== at) continue;
    if (event.kind === 'insertion') {
      const opened = insertions.find(([where]) => where === at);
      from = Math.min(from, colOf[at] - (opened ? opened[1] : 0));
    } else {
      to = Math.max(to, colOf[Math.min(event.end, referenceLength)]);
    }
  }
  return [from, Math.max(to, from + 1)];
}

/** The runs in `[from, to)`, each with the x it must be drawn at.
 *
 * The x is absolute - `column * charWidth` from the start of the whole
 * alignment, never from the start of the window being painted. That
 * distinction is the whole reason this is a function.
 *
 * Only the columns on screen are ever in the DOM, but they sit inside a row
 * as wide as the entire alignment. A run positioned at its offset *within the
 * window* therefore lands correctly only while the window starts at zero, and
 * silently slides left by `from` columns once it does not: the sequence looks
 * like it stops partway across and the rest of the row stays blank however far
 * you scroll. It shipped that way, and reads exactly like truncated data.
 */
export function placed(referenceRow, row, from, to, charWidth) {
  return runs(referenceRow, row, from, to).map(([kind, at, text]) => ({
    kind,
    at,
    left: at * charWidth,
    text,
  }));
}

/** Ruler labels for `[from, to)`: every `every`-th reference base, placed.
 *
 * Absolute for the same reason, and skipping insertion columns, which belong
 * to no reference base and so have no number to show.
 */
export function ticks(refOf, from, to, every, charWidth) {
  const out = [];
  for (let i = from; i < to; i += 1) {
    const at = refOf[i];
    if (at < 0 || (at + 1) % every) continue;
    out.push({ at, left: i * charWidth, label: at + 1 });
  }
  return out;
}
