/* Level 2 map geometry. Pure arithmetic - no DOM, no fetch, no state.
 *
 * It lives apart from level2.js because the interesting part is not the
 * drawing but the collision pass, and that is worth testing rather than
 * eyeballing: a cassette with two 250 bp parts side by side gives two labels
 * 22 px apart, and the rule that separates them without letting their leader
 * lines cross is easy to get subtly wrong.
 *
 * Angles are degrees clockwise from 12 o'clock, which is what RING_START_DEG
 * in level2.js measures too. Screen coordinates follow from
 *   x = cx + r sin(theta),  y = cy - r cos(theta)
 * so theta 0 is the top, 90 is 3 o'clock, 180 the bottom.
 */

export const RING = {
  width: 680,
  height: 500,
  cx: 360,
  cy: 270,
  outer: 150,
  bandInner: 116,
  featureOuter: 110,
  featureInner: 96,
  chevron: 133,
  stub: 22,
  pad: 12,
};

export const LABEL = {
  height: 30,
  gap: 36,
  tightHeight: 26,
  tightGap: 30,
  elbowThreshold: 10,
};

/** Where `theta` degrees clockwise from 12 lands at radius `r`. */
export function pointAt(r, theta, ring = RING) {
  const rad = (theta * Math.PI) / 180;
  return { x: ring.cx + r * Math.sin(rad), y: ring.cy - r * Math.cos(rad) };
}

/** Degrees into the circle, normalised to [0, 360). */
export function normalise(theta) {
  return ((theta % 360) + 360) % 360;
}

/** The angle at the middle of a part's arc. */
export function midAngle(part, total, startDeg) {
  const mid = (part.start + part.length / 2) / total;
  return normalise(startDeg + mid * 360);
}

/** The largest annotated feature in a part, or null. */
export function primaryFeature(part) {
  const features = part.features || [];
  if (!features.length) return null;
  return features.reduce((best, f) =>
    (f.end - f.start) > (best.end - best.start) ? f : best);
}

/** Whether a part's biggest feature says which way it is read. */
export function strandOf(part) {
  const feature = primaryFeature(part);
  return feature && feature.strand ? feature.strand : 0;
}

/**
 * Place one callout label per part, separated enough to read.
 *
 * Two rules matter and both come from leader lines crossing:
 *   - labels are ordered by angle and that order is never changed, so a label
 *     pushed past its neighbour would drag its own leader across the other's;
 *   - the correction is applied forwards and then backwards, so a crowded
 *     group spreads in both directions instead of sliding down the page.
 */
export function layoutCallouts(parts, total, startDeg, ring = RING, label = LABEL) {
  const placed = parts.map((part, index) => {
    const theta = midAngle(part, total, startDeg);
    const stubStart = pointAt(ring.outer, theta, ring);
    const stubEnd = pointAt(ring.outer + ring.stub, theta, ring);
    return {
      index,
      part,
      theta,
      side: theta < 180 ? 'right' : 'left',
      stubStart,
      stubEnd,
      naturalY: stubEnd.y,
      y: stubEnd.y,
    };
  });

  let height = label.height;
  let gap = label.gap;
  for (let attempt = 0; attempt < 2; attempt += 1) {
    for (const side of ['right', 'left']) {
      const column = placed.filter((item) => item.side === side);
      // top to bottom *by angle*: clockwise on the right, anticlockwise on the
      // left. Sorting by y instead is what lets leaders cross.
      column.sort((a, b) => (side === 'right' ? a.theta - b.theta : b.theta - a.theta));
      for (const item of column) item.y = item.naturalY;
      separate(column, gap);
      clamp(column, ring, height);
    }
    if (!overflows(placed, ring, height)) break;
    // one step tighter before giving up on the layout entirely
    height = label.tightHeight;
    gap = label.tightGap;
  }

  for (const item of placed) {
    item.displacement = item.y - item.naturalY;
    item.needsElbow = Math.abs(item.displacement) > label.elbowThreshold;
    item.labelHeight = height;
    item.anchor = item.side === 'right' ? 'start' : 'end';
    item.dot = item.stubEnd;
    item.labelX = item.stubEnd.x;
  }
  return { items: placed, labelHeight: height, gap };
}

/**
 * Push neighbours apart to `gap`, spreading the correction both ways.
 *
 * The obvious version - one pass down setting `y[i] = y[i-1] + gap`, then one
 * pass up - does not work: after the downward pass every gap is already at
 * least `gap`, so the upward pass has nothing left to do and the whole crowded
 * group has slid toward the bottom of the page.
 *
 * Splitting each overlap between the two labels that share it gives what that
 * was reaching for. A pair is always separated symmetrically, so their order
 * is never swapped; a label shoved into its own neighbour is settled on the
 * next pass, and the group ends up centred near where it started.
 */
function separate(column, gap, passes = 40) {
  for (let pass = 0; pass < passes; pass += 1) {
    let worst = 0;
    for (let i = 1; i < column.length; i += 1) {
      const overlap = gap - (column[i].y - column[i - 1].y);
      if (overlap > 0.01) {
        column[i - 1].y -= overlap / 2;
        column[i].y += overlap / 2;
        worst = Math.max(worst, overlap);
      }
    }
    if (worst <= 0.01) return;
  }
}

/** Keep the column inside the container without losing its spacing. */
function clamp(column, ring, height) {
  if (!column.length) return;
  const top = ring.pad + height / 2;
  const bottom = ring.height - ring.pad - height / 2;
  const first = column[0];
  const last = column[column.length - 1];
  const shift = Math.max(0, top - first.y) - Math.max(0, last.y - bottom);
  if (shift) for (const item of column) item.y += shift;
}

function overflows(placed, ring, height) {
  const top = ring.pad + height / 2;
  const bottom = ring.height - ring.pad - height / 2;
  return placed.some((item) => item.y < top - 0.5 || item.y > bottom + 0.5);
}
