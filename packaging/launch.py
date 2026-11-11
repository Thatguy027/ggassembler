"""Entry point for the bundled app.

A double-clicked app gets no arguments, so there is almost nothing here to
parse: the window asks for a library the first time and remembers it after
that. The CLI keeps its own entry point - this one exists so the bundle has
something to start that is not ``argparse``.

The one thing it does read is ``GGASM_SELF_TEST``, which CI sets. A bundle that
cannot start is the failure that does not show up in a build log: PyInstaller
reports success for an app whose first import is missing, and the only way to
find out is to run it. A window cannot be run on a build machine, so the
self-test takes the same path up to the window and stops there.

An environment variable rather than a flag because a windowed bundle never
sees one. PyInstaller's GUI bootloader hands Python an empty ``sys.argv``, so
``--self-test`` was silently ignored and the build machine opened the app
instead of testing it - which looks exactly like a hang.
"""

from __future__ import annotations

import multiprocessing
import os
import sys


def self_test(folder: str) -> int:
    """Start the server the way the window would, then exit.

    Exercises the parts a bundle actually breaks: the package imports, the
    frontend files being where ``api.main`` looks for them, and uvicorn
    binding a port.
    """
    import urllib.request

    from ggassembler import desktop

    port = desktop.free_port()
    server = desktop.serve_in_thread([folder], "127.0.0.1", port)
    if not server.started:
        print("the server did not start", file=sys.stderr)
        return 1

    with urllib.request.urlopen(f"http://127.0.0.1:{port}/health", timeout=30) as response:
        if response.status != 200:
            print(f"/health answered {response.status}", file=sys.stderr)
            return 1

    # The frontend is data, not code, so it is the thing most likely to be
    # missing from a bundle while every import still succeeds.
    with urllib.request.urlopen(f"http://127.0.0.1:{port}/", timeout=30) as response:
        body = response.read()
    if b"<title>" not in body:
        print("the app served no page", file=sys.stderr)
        return 1

    server.should_exit = True
    print(f"bundle starts: served {len(body)} bytes from {desktop.APP_NAME}")
    return 0


#: Set to a library folder to start the server, check it answers, and exit.
ENV_SELF_TEST = "GGASM_SELF_TEST"


def main() -> int:
    # uvicorn's reloader and anything else that spawns would otherwise restart
    # the bundle itself, which on Windows opens a second app instead of a
    # worker. Harmless here, required the moment anything spawns.
    multiprocessing.freeze_support()

    folder = os.environ.get(ENV_SELF_TEST)
    if folder:
        return self_test(folder)

    from ggassembler.desktop import run

    return run()


if __name__ == "__main__":
    sys.exit(main())
