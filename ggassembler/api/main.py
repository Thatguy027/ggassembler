"""FastAPI app factory: mounts the routers and serves ``web/`` as static files.

The app holds exactly one `Library`, on ``app.state.library``. Routers reach it
through the request rather than through a shared module, so no router has to
import another one.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Sequence

from fastapi import FastAPI
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles

from ..core.library import Library
from ..core.parttypes import YTK, Scheme
from . import (
    routes_level1,
    routes_level2,
    routes_level3,
    routes_library,
    routes_plate,
    routes_sequencing,
    routes_sync,
    routes_uniprot,
)

WEB = Path(__file__).resolve().parent.parent / "web"

#: The package itself - what ``--reload`` watches.
PACKAGE = Path(__file__).resolve().parent.parent

#: How ``ggasm serve --reload`` hands the library to the subprocess.
ENV_FOLDERS = "GGASM_FOLDERS"
ENV_RECURSIVE = "GGASM_RECURSIVE"
ENV_CACHE_DIR = "GGASM_CACHE_DIR"
ENV_DATA_DIR = "GGASM_DATA_DIR"

#: Each screen: its URL, the folder under web/, and the tab it lights up.
SCREENS = {
    "/": "level2",
    "/cassette": "level2",
    "/part": "level1",
    "/multigene": "level3",
    "/library": "library",
    "/proteins": "uniprot",
    "/sequencing": "sequencing",
    "/plate": "plate",
}


def create_app(
    roots: str | Path | Sequence[str | Path],
    recursive: bool = True,
    scheme: Scheme = YTK,
    scan: bool = True,
    cache_dir: str | Path | None = None,
    data_dir: str | Path | None = None,
) -> FastAPI:
    app = FastAPI(title="GG Assembler", version="0.1.0")

    library = Library(roots, recursive=recursive, scheme=scheme,
                      cache_dir=cache_dir, data_dir=data_dir)
    if scan:
        library.scan()
    app.state.library = library

    app.include_router(routes_library.router)
    app.include_router(routes_uniprot.router)
    app.include_router(routes_sequencing.router)
    app.include_router(routes_sync.router)
    app.include_router(routes_plate.router)
    app.include_router(routes_level2.router)

    # Level 1 and Level 3 come online when their routers exist; until then the
    # modules are empty and the app runs without them. Adding one touches no
    # other level's file.
    for module in (routes_level1, routes_level3):
        router = getattr(module, "router", None)
        if router is not None:
            app.include_router(router)

    # Always revalidate. This is a local tool whose stylesheets and scripts
    # change under a running browser, and a cached one is indistinguishable
    # from a change that did not work: the page renders with new markup and
    # old rules, which is how a flex diagram comes out stacked vertically.
    # ETags make revalidation nearly free, so nothing is actually re-sent
    # unless it changed.
    class FreshStatic(StaticFiles):
        def file_response(self, *args: object, **kwargs: object):  # type: ignore[override]
            response = super().file_response(*args, **kwargs)
            response.headers["Cache-Control"] = "no-cache, must-revalidate"
            return response

    app.mount("/static", FreshStatic(directory=WEB), name="static")

    for url, folder in SCREENS.items():
        app.get(url, include_in_schema=False)(_screen(folder))

    @app.get("/health", include_in_schema=False)
    def health() -> dict[str, object]:
        return {"ok": True, "plasmids": len(library.entries), "scheme": library.scheme.name}

    return app


def from_env() -> FastAPI:
    """Build the app from the environment, for ``ggasm serve --reload``.

    Reload mode runs the app in a subprocess that uvicorn imports by name, so
    there is no way to pass it the folders as arguments: nothing built in the
    parent survives the respawn. The CLI puts them in the environment and this
    reads them back.

    JSON rather than a separator-joined string because a path may legally
    contain whatever separator one picks, and a library folder with a colon in
    its name should not quietly index half of itself.
    """
    raw = os.environ.get(ENV_FOLDERS)
    if not raw:
        raise RuntimeError(
            f"{ENV_FOLDERS} is not set. This factory exists so that "
            f"`ggasm serve --reload` can respawn the app; to build one "
            f"yourself, call create_app() with the folders."
        )
    return create_app(
        [Path(folder) for folder in json.loads(raw)],
        recursive=os.environ.get(ENV_RECURSIVE, "1") != "0",
        cache_dir=os.environ.get(ENV_CACHE_DIR) or None,
        data_dir=os.environ.get(ENV_DATA_DIR) or None,
    )


def _screen(folder: str):
    """Serve one level's index.html, or a placeholder until it is built."""

    def handler() -> FileResponse:
        index = WEB / folder / "index.html"
        if not index.exists():
            index = WEB / "not_built.html"
        return FileResponse(index)

    handler.__name__ = f"screen_{folder}"
    return handler
