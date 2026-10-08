"""Sharing a plasmid library with the rest of an organisation, over git.

The library is a folder of GenBank files that is already a git repository, and
git already answers every question this feature asks. It replicates a whole
copy to every machine, so the app keeps working on a plane. It merges two
people's additions without being told how. It records who changed a sequence
and when, which for a plasmid library is worth as much as the sequence. And a
private repository per organisation gives privacy between organisations with no
server to run, no accounts to store and - the part that matters for adoption -
nobody's unpublished constructs sitting on infrastructure somebody else owns.

So this module adds no sync protocol. It drives `git` through a subprocess and
spends its effort on the two things git does not provide: deciding what is safe
to publish, and saying what happened in words a bench scientist can act on.

Shelling out rather than taking a library dependency is deliberate. Every
machine that cloned the library already has git, and anything this module does
can be reproduced by typing the same commands by hand - which is what you want
at 11pm when someone's push is rejected and the app is the only thing standing
between them and their data.
"""

from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

from .library import CACHE_DIR, DATA_DIR, Library

#: Long enough for a slow push over a VPN, short enough that a hung network
#: does not freeze the screen the user is looking at.
TIMEOUT = 120

#: What the shared directory needs git to know. `merge=union` is the whole
#: reason the build log is one JSON object per line: two people logging builds
#: on the same afternoon append to the same file, and union keeps both sets of
#: lines instead of raising a conflict over work that does not overlap.
ATTRIBUTES = f"{DATA_DIR}/builds.jsonl merge=union\n"

#: Derived state must never be committed. It is keyed on mtimes, so it differs
#: on every machine and would conflict on every pull while meaning nothing.
IGNORES = (f"{CACHE_DIR}/", ".DS_Store")


class SyncError(RuntimeError):
    """Something went wrong that the person at the keyboard has to decide about."""


@dataclass
class Status:
    """What sharing would do right now, in terms the screen can show."""

    repo: Path | None = None
    """The library's git repository, or None if the library is not in one."""
    remote: str = ""
    branch: str = ""
    incoming: int = 0
    """Commits on the remote this copy has not got."""
    outgoing: int = 0
    """Commits here that nobody else has."""
    unshared: list[str] = field(default_factory=list)
    """Files changed or added locally and never committed."""
    detail: str = ""
    """Why sharing is unavailable, when it is."""

    @property
    def available(self) -> bool:
        return self.repo is not None and bool(self.remote)

    @property
    def pending(self) -> int:
        return self.outgoing + len(self.unshared)

    def summary(self) -> str:
        """One line for the badge, saying the thing worth acting on."""
        if not self.available:
            return self.detail or "not set up for sharing"
        parts = []
        if self.incoming:
            parts.append(f"{self.incoming} update{_s(self.incoming)} from your lab")
        if self.pending:
            parts.append(f"{self.pending} of yours not shared")
        return ", ".join(parts) or "up to date"


def _s(n: int) -> str:
    return "" if n == 1 else "s"


def git_available() -> bool:
    return shutil.which("git") is not None


def run(repo: Path, *args: str, check: bool = True) -> str:
    """One git command, with its stderr kept for the message if it fails."""
    if not git_available():
        raise SyncError("git is not installed, so the library cannot be shared")
    try:
        done = subprocess.run(
            ["git", *args], cwd=repo, capture_output=True, text=True, timeout=TIMEOUT,
        )
    except subprocess.TimeoutExpired:
        raise SyncError(
            f"git {args[0]} did not finish within {TIMEOUT}s - check the network "
            f"or your access to the remote"
        ) from None
    if check and done.returncode != 0:
        raise SyncError((done.stderr or done.stdout).strip() or f"git {args[0]} failed")
    return done.stdout.strip()


def repo_of(library: Library) -> Path | None:
    """The git repository holding the library, if there is exactly one.

    A library can be several folders. Shared, they may be several repositories
    - a public base library plus a private one - and publishing into "the"
    repository would then be a guess. Refusing to guess is the point: the
    wrong guess publishes an organisation's private plasmid to a public repo.
    """
    if not git_available():
        return None
    found: set[Path] = set()
    for root in library.roots:
        try:
            top = run(root, "rev-parse", "--show-toplevel")
        except (SyncError, OSError):
            continue
        if top:
            found.add(Path(top).resolve())
    return found.pop() if len(found) == 1 else None


#: The setting, kept per person per machine. Default on: a lab library that
#: only fills up when people remember to share is a lab library that does not.
SETTING = "share_plasmids"


#: Marks a library you read from rather than contribute to. Set when `ggasm
#: init --starter` clones the published kits, which everyone shares and nobody
#: should be publishing their own constructs into.
REFERENCE = "reference_library"


def enabled(library: Library) -> bool:
    return bool(library.settings().get(SETTING, True))


def is_reference(library: Library) -> bool:
    return bool(library.settings().get(REFERENCE, False))


def mark_reference(library: Library) -> None:
    """Make a library read-only as far as publishing is concerned.

    A starter library's remote is a public repository that everyone clones.
    Default-on sharing plus a save made while pointed at it would push one
    lab's construct into the reference every other lab reads - and the save
    that did it would look exactly like every other save.
    """
    library.set_setting(REFERENCE, True)
    library.set_setting(SETTING, False)


def set_enabled(library: Library, on: bool) -> bool:
    library.set_setting(SETTING, bool(on))
    return bool(on)


def paths_to_share(library: Library, repo: Path) -> list[str]:
    """The library's own paths inside the repository, and nothing else.

    A plasmid repository is rarely only plasmids - this one sits inside an R
    package. Staging the whole worktree because somebody saved a plasmid would
    sweep up a half-finished script and publish it under a commit message about
    a plasmid. So sharing names its own paths and leaves the rest of the
    repository to whoever is working on it.
    """
    out = []
    for folder in [*library.roots, library.data_dir]:
        try:
            out.append(folder.resolve().relative_to(repo).as_posix() or ".")
        except ValueError:
            continue  # outside this repository: not ours to publish
    return sorted(set(out))


def prepare(library: Library, repo: Path) -> list[str]:
    """Teach the repository the two things it has to know. Idempotent.

    Returns what it changed, so setting up can say so rather than happening
    invisibly.
    """
    changed = []
    attrs = repo / ".gitattributes"
    existing = attrs.read_text(encoding="utf-8") if attrs.exists() else ""
    if "builds.jsonl" not in existing:
        prefix = "" if (not existing or existing.endswith("\n")) else "\n"
        attrs.write_text(existing + prefix + ATTRIBUTES, encoding="utf-8")
        changed.append(".gitattributes")

    ignore = repo / ".gitignore"
    existing = ignore.read_text(encoding="utf-8") if ignore.exists() else ""
    lines = {line.strip() for line in existing.splitlines()}
    missing = [rule for rule in IGNORES if rule not in lines]
    if missing:
        prefix = "" if (not existing or existing.endswith("\n")) else "\n"
        ignore.write_text(existing + prefix + "\n".join(missing) + "\n", encoding="utf-8")
        changed.append(".gitignore")
    return changed


def status(library: Library, fetch: bool = False) -> Status:
    """What sharing would do right now.

    `fetch` touches the network, so the screen can ask for a cheap answer from
    what is already known and a true one only when the person asks for it.
    """
    if not git_available():
        return Status(detail="git is not installed")
    repo = repo_of(library)
    if repo is None:
        return Status(detail="the library is not in a single git repository")

    found = Status(repo=repo)
    try:
        # `symbolic-ref` rather than `rev-parse --abbrev-ref HEAD`, which cannot
        # name a branch before the first commit and fails with "ambiguous
        # argument 'HEAD'" - an unhelpful thing to show someone who has just
        # run `git init` on a folder of plasmids.
        found.branch = run(repo, "symbolic-ref", "--short", "HEAD")
        found.remote = run(repo, "remote", check=False).splitlines()[0] if run(
            repo, "remote", check=False) else ""
    except SyncError as exc:
        return Status(repo=repo, detail=str(exc))

    if not found.remote:
        found.detail = "the repository has no remote, so there is nobody to share with"
        return found

    if fetch:
        try:
            run(repo, "fetch", "--quiet", found.remote)
        except SyncError as exc:
            found.detail = str(exc)

    upstream = f"{found.remote}/{found.branch}"
    counts = run(repo, "rev-list", "--left-right", "--count", f"{upstream}...HEAD",
                 check=False)
    if counts and len(counts.split()) == 2:
        found.incoming, found.outgoing = (int(n) for n in counts.split())

    mine = paths_to_share(library, repo)
    porcelain = run(repo, "status", "--porcelain", "--", *mine, check=False)
    found.unshared = [
        name for name in (line[3:].strip() for line in porcelain.splitlines())
        # Filtered here as well as ignored in .gitignore, because the ignore
        # rule is only written when sharing is set up: until then a fresh clone
        # would report its own index as work waiting to be published.
        if name and f"{CACHE_DIR}/" not in f"{name}/"
    ]
    return found


def pull(library: Library) -> Status:
    """Take what the rest of the organisation has. Always allowed.

    Receiving is never gated and never switched off by the toggle. Someone who
    has stopped publishing still needs to see their colleagues' plasmids -
    working against a library you know to be stale is a way to repeat an
    assembly that somebody already did.
    """
    found = status(library, fetch=True)
    if not found.available:
        raise SyncError(found.detail or "sharing is not set up")
    run(found.repo, "pull", "--rebase", "--autostash", found.remote, found.branch)
    library.scan()
    return status(library)


def share(library: Library, message: str, force: bool = False) -> tuple[Status, "object"]:
    """Publish what is here, if the library is no worse than its baseline.

    The gate is the whole reason this is not just `git push`. A plasmid that
    will not digest, or one nothing can classify, is a plasmid that breaks a
    colleague's dropdown rather than yours - and by then it is in the history
    of everybody's clone.
    """
    from . import baseline  # local: baseline reads a Library, as this does

    found = status(library)
    if is_reference(library):
        raise SyncError(
            f"{library.folder} is a reference library - the published kits that "
            f"everyone clones - so your own constructs do not belong in it. "
            f"Point the app at your lab's library and save there instead."
        )
    if not found.available:
        raise SyncError(found.detail or "sharing is not set up")

    verdict = baseline.compare(library)
    if not verdict.ok and not force:
        return found, verdict

    prepare(library, found.repo)
    run(found.repo, "add", "--", *paths_to_share(library, found.repo))
    staged = run(found.repo, "diff", "--cached", "--name-only", check=False)
    if staged:
        run(found.repo, "commit", "--quiet", "-m", message)
    run(found.repo, "push", found.remote, f"HEAD:{found.branch}")
    return status(library), verdict


#: What git says when the remote is there but you are not allowed in. Matched
#: so the app can answer the actual question - "set up your GitHub access" -
#: rather than passing on a sentence about public keys and shell access.
DENIED = (
    "permission denied",
    "could not read from remote repository",
    "authentication failed",
    "repository not found",
    "access rights",
)


#: The two published kits, as a repository anyone can read without being in a
#: lab. `ggasm init` falls back to it so that installing the app and having
#: nothing to open is not the first thing that happens to a new user. HTTPS
#: rather than SSH on purpose: this one needs no key, and asking someone to
#: set up GitHub access before they have seen the app work is the wrong order.
STARTER = "https://github.com/Thatguy027/ytk-starter-library.git"


def clone(url: str, into: Path) -> Path:
    """Clone an organisation's library, or say why not in words.

    This is where a new person meets git for the first time, and it is the
    step most likely to stop them: a lab member who has never pushed anything
    has no key on their account, and git's own message for that talks about
    public keys and shell access rather than about what to do next.
    """
    if not git_available():
        raise SyncError(
            "git is not installed. Install it (on macOS, `xcode-select --install`) "
            "and run this again."
        )
    into = Path(into).expanduser().resolve()
    if into.exists() and any(into.iterdir()):
        raise SyncError(f"{into} already exists and is not empty")

    into.parent.mkdir(parents=True, exist_ok=True)
    try:
        done = subprocess.run(
            ["git", "clone", url, str(into)],
            capture_output=True, text=True, timeout=TIMEOUT * 5,
        )
    except subprocess.TimeoutExpired:
        raise SyncError(f"cloning {url} took too long - check the network") from None

    if done.returncode != 0:
        raise SyncError(clone_failure(url, (done.stderr or done.stdout).strip()))
    return into


def clone_failure(url: str, message: str) -> str:
    """Turn git's refusal into something a bench scientist can act on.

    Separated from `clone` so the wording is reachable without a network and a
    repository nobody can read: the message is the feature here, and a feature
    only testable by failing to reach GitHub is a feature nobody tests.
    """
    if not any(phrase in message.lower() for phrase in DENIED):
        return message or f"could not clone {url}"
    return (
        f"GitHub would not let you read {url}.\n\n"
        f"  Either you have no access to it - ask whoever runs the library "
        f"to add you - or this machine has no key on your account.\n"
        f"  To check: `ssh -T git@github.com`.\n"
        f"  To fix:   `gh auth login`, or add an SSH key at "
        f"https://github.com/settings/keys\n\n"
        f"git said: {message}"
    )
