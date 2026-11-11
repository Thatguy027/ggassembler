"""The app in a native window, with no browser and no terminal.

Same server as ``ggasm serve``: a thread runs uvicorn on a loopback port, and a
`pywebview` window is pointed at it. What the window adds is the two things a
browser was standing in for - a way to choose the library, and somewhere to
remember the choice - because a bundled app is launched by double-clicking it
and gets no arguments at all.

The window is created *before* the library is known, showing a splash, and the
real URL is loaded into that same window once the server answers. This is not
an aesthetic choice: the folder dialog is a method on a window, so there has to
be a window before anything can be picked, and the server cannot start before a
library is picked. One window that changes what it is showing avoids both a
second window flashing up and a dialog with nothing behind it.
"""

from __future__ import annotations

import json
import os
import shutil
import socket
import sys
import threading
from pathlib import Path
from typing import Sequence

APP_NAME = "GG Assembler"

#: The kits the bundle carries, laid down on first launch. Only a bundle has
#: them - `packaging/fetch_seed.py` clones them in at build time - so a pip
#: install finds nothing here and asks for a folder instead.
SEED = Path(__file__).resolve().parent / "seed"

#: Where they are laid down. Documents rather than Application Support: this
#: is the person's own library from the moment it exists, they will be asked
#: to find it in a file dialog, and they will want to put it under git one day.
SEEDED = Path.home() / "Documents" / f"{APP_NAME} Library"

#: How long to wait for uvicorn to claim its port before loading the window
#: anyway. Generous, because the first launch against an unindexed library
#: scans every file in it; a window that gives up early is worse than one that
#: sits on the splash a moment longer.
STARTUP_TIMEOUT = 120.0

#: Where the chosen library is remembered between launches.
if sys.platform == "darwin":
    SETTINGS = Path.home() / "Library" / "Application Support" / APP_NAME / "settings.json"
elif os.name == "nt":
    SETTINGS = Path(os.environ.get("APPDATA", Path.home())) / APP_NAME / "settings.json"
else:
    SETTINGS = Path(
        os.environ.get("XDG_CONFIG_HOME", Path.home() / ".config")
    ) / "gg-assembler" / "settings.json"

SPLASH = """
<!doctype html><meta charset="utf-8">
<style>
  :root { color-scheme: light dark; }
  html, body { height: 100%; margin: 0; }
  body {
    display: flex; align-items: center; justify-content: center;
    font: 15px/1.5 -apple-system, BlinkMacSystemFont, "Segoe UI", system-ui, sans-serif;
    background: Canvas; color: CanvasText;
  }
  p { opacity: .65; }
</style>
<p id="m">Starting GG Assembler &hellip;</p>
"""


class Chosen:
    """What the bootstrap thread works out, in the order it works it out."""

    def __init__(self) -> None:
        self.folders: list[Path] = []
        self.server = None


def remembered() -> list[Path]:
    """The library from last time, if it is still a folder.

    A library that has been moved or deleted is treated as no library rather
    than as an error: the next thing that happens is the picker, which is what
    someone whose folder moved needs anyway.
    """
    try:
        saved = json.loads(SETTINGS.read_text())
    except (OSError, ValueError):
        return []
    folders = [Path(p) for p in saved.get("folders", [])]
    return folders if folders and all(f.is_dir() for f in folders) else []


def remember(folders: Sequence[Path]) -> None:
    try:
        SETTINGS.parent.mkdir(parents=True, exist_ok=True)
        SETTINGS.write_text(json.dumps(
            {"folders": [str(Path(f).resolve()) for f in folders]}, indent=2))
    except OSError:
        # Not being able to remember is worth none of the app: the window is
        # already up and the library already works for this session.
        pass


def seed() -> list[Path]:
    """Lay the bundled kits down as a library of the person's own, once.

    The bundle is for someone who will not open a terminal, and that person
    cannot run ``ggasm init --starter``: without this they would meet a file
    dialog asking for plasmids they have no way to have got. So the app opens
    on 192 published plasmids and the picker becomes something they reach for
    when they have their own library, rather than a wall on first launch.

    Copied out rather than served from inside the bundle. Everything about a
    library is written to - the index, saved constructs, the build log, the
    designs - and the bundle is both read-only in principle and replaced
    wholesale by the next version, which would take every saved construct with
    it.

    Returns nothing when there is no seed (a pip install) or when the folder
    already exists, which is the ordinary case from the second launch on: a
    remembered library is found before this is ever called.
    """
    if not SEED.is_dir() or SEEDED.exists():
        return []
    try:
        # copytree to a temporary name and then rename, so a copy interrupted
        # half way does not leave a folder that looks seeded and is missing
        # half the kit - the next launch would index it and quietly offer
        # fewer parts.
        partial = SEEDED.with_name(SEEDED.name + " (copying)")
        if partial.exists():
            shutil.rmtree(partial)
        partial.parent.mkdir(parents=True, exist_ok=True)
        shutil.copytree(SEED, partial)
        partial.rename(SEEDED)
    except OSError:
        # A read-only or missing Documents folder is not worth refusing to
        # start over: the picker is still there, and still works.
        return []
    return [SEEDED]


def free_port(host: str = "127.0.0.1") -> int:
    """A port the OS picked, so two copies of the app never collide.

    ``ggasm serve`` defaults to 8737 and refuses to start when it is taken,
    which is right for a server someone is pointing a browser at by hand. A
    window carries its own address, nobody types it, and a lab member opening
    the app twice should get two windows rather than an error.
    """
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as probe:
        probe.bind((host, 0))
        return probe.getsockname()[1]


def serve_in_thread(folders: Sequence[Path], host: str, port: int, **kwargs: object):
    """Start uvicorn on a daemon thread and return the server once it is up."""
    import uvicorn

    from .api.main import create_app

    app = create_app(folders, **kwargs)  # type: ignore[arg-type]
    config = uvicorn.Config(app, host=host, port=port, log_level="warning")
    server = uvicorn.Server(config)
    threading.Thread(target=server.run, daemon=True).start()

    # Wait for the socket, not for a fixed delay. Indexing a few hundred
    # plasmids is quick with a warm cache and slow without one, and a window
    # loaded too early shows a connection error that never retries.
    waited = 0.0
    while not server.started and waited < STARTUP_TIMEOUT:
        threading.Event().wait(0.05)
        waited += 0.05
    return server


def run(
    folders: Sequence[str | Path] | None = None,
    recursive: bool = True,
    cache_dir: str | Path | None = None,
    data_dir: str | Path | None = None,
    host: str = "127.0.0.1",
) -> int:
    """Open the app in a native window. Returns once the window is closed."""
    try:
        import webview
        import webview.menu
    except ImportError:
        print(
            "the desktop window needs pywebview, which is not installed.\n\n"
            "  pip install 'ggassembler[desktop]'\n\n"
            "Or use the browser:  ggasm serve <library>",
            file=sys.stderr,
        )
        return 2

    chosen = Chosen()
    chosen.folders = [Path(f) for f in (folders or [])]
    window = webview.create_window(APP_NAME, html=SPLASH, width=1280, height=860)

    def serve(where: Sequence[Path]) -> None:
        """Point the window at these folders, replacing whatever it was showing.

        A new server on a new port rather than a reconfigured one: the library
        is built at ``create_app`` time and reaches the routers through the
        app, so there is nothing to swap underneath a running server. The old
        one is told to exit and its port is simply abandoned.
        """
        previous = chosen.server
        port = free_port(host)
        chosen.server = serve_in_thread(
            where, host, port,
            recursive=recursive, cache_dir=cache_dir, data_dir=data_dir,
        )
        chosen.folders = list(where)
        remember(where)
        window.load_url(f"http://{host}:{port}/")
        if previous is not None:
            previous.should_exit = True

    def ask() -> list[Path]:
        picked = window.create_file_dialog(webview.FOLDER_DIALOG)
        return [Path(p) for p in picked] if picked else []

    def open_library() -> None:
        """File ▸ Open library. The only way to switch without a terminal."""
        picked = ask()
        if picked:
            serve(picked)

    def bootstrap(_window: object = None) -> None:
        where = chosen.folders or remembered() or seed() or ask()
        if not where:
            # Nothing to show and nothing to ask: closing the window is the
            # honest response to a cancelled picker on a first launch.
            window.destroy()
            return
        serve(where)

    menu = [webview.menu.Menu("File", [
        webview.menu.MenuAction("Open library…", open_library),
    ])]

    # Eight screens keep something in localStorage - the theme, which columns
    # are shown, what was last typed into a picker. pywebview defaults to a
    # private session that is thrown away on quit, which would make the window
    # forget between launches everything the same app remembers in a browser.
    webview.start(
        bootstrap, window, menu=menu,
        private_mode=False, storage_path=str(SETTINGS.parent / "webview"),
    )

    if chosen.server is not None:
        chosen.server.should_exit = True
    return 0
