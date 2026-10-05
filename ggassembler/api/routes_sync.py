"""Sharing the library with the rest of an organisation (``/api/sync/*``).

Thin on purpose: every decision about what may be published lives in
`core.sync` and `core.baseline`, so the rule is the same whether it was reached
through this screen or through the command line.

One thing is decided here rather than there, because it is about a request and
not about a library: sharing after a save runs as a background task. A push
over a VPN can take the better part of a minute, and a save that blocks on it
is a save that feels broken. The outcome is remembered so the badge can report
a failure rather than letting it disappear - the worst version of this feature
is one where a push quietly fails and you believe your work is in the library.
"""

from __future__ import annotations

from typing import Any

from fastapi import APIRouter, BackgroundTasks, HTTPException, Request
from pydantic import BaseModel

from ..core import baseline, sync
from ..core.library import Library

router = APIRouter(prefix="/api/sync", tags=["sync"])


def get_library(request: Request) -> Library:
    return request.app.state.library


class ToggleRequest(BaseModel):
    on: bool


class ShareRequest(BaseModel):
    message: str = "update plasmids"
    force: bool = False


def _state(library: Library, found: sync.Status, last: dict[str, Any] | None) -> dict[str, Any]:
    return {
        "enabled": sync.enabled(library),
        "available": found.available,
        "detail": found.detail,
        "summary": found.summary(),
        "remote": found.remote,
        "branch": found.branch,
        "incoming": found.incoming,
        "outgoing": found.outgoing,
        "unshared": found.unshared,
        "pending": found.pending,
        "last": last,
    }


def _last(request: Request) -> dict[str, Any] | None:
    return getattr(request.app.state, "last_share", None)


def _remember(request: Request, ok: bool, message: str) -> dict[str, Any]:
    record = {"ok": ok, "message": message}
    request.app.state.last_share = record
    return record


@router.get("")
def read(request: Request, fetch: bool = False) -> dict[str, Any]:
    """What sharing would do. `fetch` asks the remote; without it, this is local."""
    library = get_library(request)
    return _state(library, sync.status(library, fetch=fetch), _last(request))


@router.post("/enabled")
def toggle(request: Request, body: ToggleRequest) -> dict[str, Any]:
    """Turn publishing off for this person, on this machine.

    Receiving is deliberately not affected. Somebody who has stopped publishing
    still needs their colleagues' plasmids: working against a library you know
    to be stale is how an assembly gets repeated.
    """
    library = get_library(request)
    sync.set_enabled(library, body.on)
    return _state(library, sync.status(library), _last(request))


@router.post("/pull")
def pull(request: Request) -> dict[str, Any]:
    library = get_library(request)
    try:
        found = sync.pull(library)
    except sync.SyncError as exc:
        raise HTTPException(status_code=409, detail=str(exc)) from exc
    return _state(library, found, _last(request))


@router.post("/share")
def share(request: Request, body: ShareRequest) -> dict[str, Any]:
    library = get_library(request)
    try:
        found, verdict = sync.share(library, body.message, force=body.force)
    except sync.SyncError as exc:
        _remember(request, False, str(exc))
        raise HTTPException(status_code=409, detail=str(exc)) from exc

    blocked = verdict is not None and not verdict.ok and not body.force
    out = _state(library, found, _remember(
        request, not blocked,
        verdict.summary() if blocked else "shared with your lab",
    ))
    out["blocked"] = blocked
    out["verdict"] = {
        "ok": bool(verdict.ok),
        "introduced": list(verdict.introduced),
        "summary": verdict.summary(),
        "baselined": verdict.baselined,
    }
    return out


@router.post("/baseline")
def rebaseline(request: Request) -> dict[str, Any]:
    """Agree that the library's current problems are the known ones.

    Deliberately an action someone takes on purpose. It is the one way to make
    the gate stop objecting, so it should be a commit with a name on it that a
    colleague can ask about, not something the app does to get out of the way.
    """
    library = get_library(request)
    state = baseline.write(library)
    return {
        "total": state["total"],
        "conflicts": len(state["conflicts"]),
        "unrecognised": len(state["unrecognised"]),
    }


def after_save(request: Request, tasks: BackgroundTasks, what: str) -> None:
    """Publish a just-saved plasmid, if this person shares and the gate allows.

    Called by the screens that save. Failures are remembered rather than
    raised: the save itself succeeded, and losing the plasmid because the
    network was down would be the wrong trade.
    """
    library = get_library(request)
    if not sync.enabled(library):
        return
    tasks.add_task(_share_quietly, request.app, library, what)


def _share_quietly(app: Any, library: Library, what: str) -> None:
    try:
        _, verdict = sync.share(library, what)
    except sync.SyncError as exc:
        app.state.last_share = {"ok": False, "message": str(exc)}
        return
    ok = verdict is None or verdict.ok
    app.state.last_share = {
        "ok": ok,
        "message": f"shared {what}" if ok else verdict.summary(),
    }
