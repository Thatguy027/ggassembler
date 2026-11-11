"""Fetch the published kits so the bundle can ship a library inside it.

Run before PyInstaller:

    python packaging/fetch_seed.py

The bundle exists for someone who will not open a terminal, and that person
cannot run ``ggasm init --starter`` to get plasmids. So the app carries a copy
of the starter library and lays it down on first launch.

Cloned at build time rather than committed. The kits already live in their own
repository, and a second copy in this one would be a fork nobody remembers to
update - the bundle would ship whatever the plasmids looked like on the day
somebody pasted them in.

``.git`` is left behind on purpose. What is laid down is the person's own
working library, not a clone of a public repository: they will save their own
constructs into it, and a repository whose remote is a kit everyone reads is
the one place those must not be pushed.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from ggassembler.core.sync import STARTER  # noqa: E402

SEED = Path(__file__).resolve().parent / "seed"

#: Everything else in the kit repository is the kit itself.
SKIP = {".git", ".github", ".ggasm"}


def main() -> int:
    if SEED.exists():
        shutil.rmtree(SEED)

    scratch = SEED.with_name("seed-clone")
    if scratch.exists():
        shutil.rmtree(scratch)

    print(f"cloning {STARTER} ...")
    done = subprocess.run(
        ["git", "clone", "--depth", "1", STARTER, str(scratch)],
        capture_output=True, text=True,
    )
    if done.returncode != 0:
        print(done.stderr.strip(), file=sys.stderr)
        return 1

    shutil.copytree(
        scratch, SEED,
        ignore=shutil.ignore_patterns(*SKIP),
    )
    shutil.rmtree(scratch)

    plasmids = [p for p in SEED.rglob("*") if p.suffix in {".gb", ".gbk", ".genbank"}]
    if not plasmids:
        print(f"::error::{STARTER} produced no plasmids", file=sys.stderr)
        return 1

    print(f"seed ready: {len(plasmids)} plasmids in {SEED}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
