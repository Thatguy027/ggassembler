/* Sharing the library with the rest of the lab.
 *
 * The control lives on this screen and nowhere else, because this is the screen
 * about the plasmid library - and because `web/` may hold no shared script, so
 * a badge on all eight screens would mean eight copies of this file drifting
 * apart. Saving a plasmid on any screen still publishes it: that happens on the
 * server, and this is where you come to see whether it worked.
 *
 * Every piece of state shown here is server state. Nothing is cached in the
 * page, because a badge that says "up to date" from memory while the push
 * failed is worse than no badge.
 */

const el = (id) => document.getElementById(id);

async function call(path, options) {
  const response = await fetch(path, options);
  const body = await response.json().catch(() => ({}));
  if (!response.ok) throw new Error(body.detail || `${path}: ${response.status}`);
  return body;
}

const post = (path, payload) => call(path, {
  method: 'POST',
  headers: { 'Content-Type': 'application/json' },
  body: JSON.stringify(payload || {}),
});

/* The dot, and what it means. Amber is "there is something to do", not "there
 * is something wrong" - an unshared plasmid is the normal state five seconds
 * after you save one. */
function tone(state) {
  if (!state.available) return 'off';
  if (state.last && state.last.ok === false) return 'bad';
  if (state.incoming || state.pending) return 'due';
  return 'ok';
}

function render(state, note) {
  const box = el('sync');
  if (!box) return;
  box.dataset.tone = tone(state);

  el('sync-summary').textContent = state.available
    ? state.summary
    : (state.detail || 'not set up for sharing');

  const where = el('sync-where');
  where.textContent = state.available ? `${state.remote}/${state.branch}` : '';

  el('sync-toggle').checked = !!state.enabled;
  el('sync-now').disabled = !state.available;

  /* The last outcome outlives the request that produced it, so a background
   * push that failed while you were on another screen is still visible when
   * you come back. */
  const last = el('sync-last');
  const message = note || (state.last ? state.last.message : '');
  last.textContent = message || '';
  last.hidden = !message;
  last.dataset.ok = state.last ? String(state.last.ok !== false) : 'true';

  const problems = el('sync-problems');
  const found = (state.verdict && state.verdict.introduced) || [];
  problems.hidden = found.length === 0;
  if (found.length) {
    el('sync-problem-count').textContent =
      `${found.length} new problem${found.length === 1 ? '' : 's'} stopped this from being shared`;
    el('sync-problem-list').innerHTML = '';
    for (const line of found) {
      const item = document.createElement('li');
      item.textContent = line;
      el('sync-problem-list').appendChild(item);
    }
  }
}

async function refresh(fetchRemote) {
  try {
    render(await call(`/api/sync?fetch=${fetchRemote ? 'true' : 'false'}`));
  } catch (error) {
    render({ available: false, detail: String(error.message || error) });
  }
}

export function startSync(onLibraryChanged) {
  if (!el('sync')) return;

  el('sync-toggle').addEventListener('change', async (event) => {
    const on = event.target.checked;
    try {
      render(await post('/api/sync/enabled', { on }),
        on ? 'Your plasmids will be shared with the lab.'
           : 'Your plasmids stay on this machine. You will still receive the lab’s.');
    } catch (error) {
      render({ available: false, detail: String(error.message || error) });
    }
  });

  el('sync-now').addEventListener('click', async () => {
    const button = el('sync-now');
    button.disabled = true;
    button.textContent = 'Syncing…';
    try {
      await post('/api/sync/pull');
      const after = await post('/api/sync/share', { message: 'update plasmids' });
      render(after, after.blocked ? null : 'Up to date with your lab.');
      if (onLibraryChanged) onLibraryChanged();
    } catch (error) {
      render({ available: false, detail: String(error.message || error) });
    } finally {
      button.disabled = false;
      button.textContent = 'Sync now';
    }
  });

  refresh(true);
  /* Slow on purpose. This is a lab library, not a chat window: a colleague's
   * plasmid arriving within the minute is soon enough, and a tighter poll
   * spends somebody's battery asking git the same question. */
  setInterval(() => refresh(true), 60000);
}
