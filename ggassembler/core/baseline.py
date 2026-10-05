"""What the library is known to be broken about, and what is newly broken.

Lifted out of `cli.py` when sharing needed the same gate the pre-commit hook
uses. The rule it encodes is the one worth having in exactly one place: judge a
change by what it *adds*, never by the absolute state.

A library of several hundred files always carries some conflicts and some
plasmids nothing can classify. A gate that refuses to let anyone publish until
all of them are fixed is a gate that gets switched off in a week. A gate that
only objects to problems this change introduced is one people keep, because it
only ever fires about something they just did.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .library import Library

#: Lives in the shared directory, so an organisation agrees on one notion of
#: "known problems" and re-baselining is a commit somebody can question.
BASELINE_FILE = "baseline.json"


def state(library: Library) -> dict[str, Any]:
    """The facts a baseline records: what is broken, and what it is broken about."""
    entries = library.unique_entries()
    return {
        "total": len(entries),
        "conflicts": {e.name: e.call.conflict for e in entries if e.call.conflict},
        "unrecognised": {
            e.name: e.call.reason
            for e in entries
            if e.call.part_type is None and not e.is_assembled
        },
    }


@dataclass
class Verdict:
    """Whether a change is safe to publish, and what to say if it is not."""

    ok: bool = True
    introduced: list[str] = field(default_factory=list)
    fixed: int = 0
    total: int = 0
    baselined: bool = True
    """False when no baseline has been recorded yet, so nothing was compared."""

    def summary(self) -> str:
        if not self.baselined:
            return "no baseline recorded yet, so nothing was checked"
        if self.ok:
            fixed = f", {self.fixed} fewer than before" if self.fixed else ""
            return f"{self.total} plasmids, nothing new{fixed}"
        n = len(self.introduced)
        return f"{n} new problem{'' if n == 1 else 's'} this change would publish"


def path_for(library: Library) -> Path:
    return library.data_dir / BASELINE_FILE


def read(library: Library) -> dict[str, Any] | None:
    try:
        return json.loads(path_for(library).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None


def write(library: Library) -> dict[str, Any]:
    """Record the current state as the agreed starting point."""
    current = state(library)
    path = path_for(library)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(current, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return current


def compare(library: Library) -> Verdict:
    """Judge the library against its baseline, counting only what is new."""
    current = state(library)
    before = read(library)
    if before is None:
        return Verdict(ok=True, total=current["total"], baselined=False)

    introduced = []
    for kind, label in (("conflicts", "conflict"), ("unrecognised", "unrecognised")):
        known = set(before.get(kind, {}))
        for name in sorted(set(current[kind]) - known):
            introduced.append(f"{name}: {label} - {current[kind][name]}")

    fixed = sum(
        len(set(before.get(kind, {})) - set(current[kind]))
        for kind in ("conflicts", "unrecognised")
    )
    return Verdict(
        ok=not introduced,
        introduced=introduced,
        fixed=fixed,
        total=current["total"],
    )
