"""FastAPI app factory: mounts the routers and serves ``web/`` as static files.

The app holds exactly one `Library`, on ``app.state.library``. Routers reach it
through the request rather than through a shared module, so no router has to
import another one.
"""

from __future__ import annotations

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
    routes_uniprot,
)

WEB = Path(__file__).resolve().parent.parent / "web"

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
) -> FastAPI:
    app = FastAPI(title="GG Assembler", version="0.1.0")

    library = Library(roots, recursive=recursive, scheme=scheme)
    if scan:
        library.scan()
    app.state.library = library

    app.include_router(routes_library.router)
    app.include_router(routes_uniprot.router)
    app.include_router(routes_sequencing.router)
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


def _screen(folder: str):
    """Serve one level's index.html, or a placeholder until it is built."""

    def handler() -> FileResponse:
        index = WEB / folder / "index.html"
        if not index.exists():
            index = WEB / "not_built.html"
        return FileResponse(index)

    handler.__name__ = f"screen_{folder}"
    return handler
