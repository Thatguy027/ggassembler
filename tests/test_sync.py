"""Sharing a library over git: what gets published, and what is stopped first.

These tests drive real repositories in a temporary directory rather than a
mocked git. The reason is the module's own premise - that anything it does can
be reproduced by typing the same commands by hand - and a mock asserts only
that this file and that one agree about what git does.
"""

from __future__ import annotations

import json
import subprocess

import pytest

from ggassembler.core import baseline, sync
from ggassembler.core.library import CACHE_DIR, Library
from ggassembler.core.seqio import write_genbank

from . import synth
from .test_assembly import CANONICAL

pytestmark = pytest.mark.skipif(not sync.git_available(), reason="git is not installed")


def git(repo, *args):
    done = subprocess.run(["git", *args], cwd=repo, capture_output=True, text=True)
    assert done.returncode == 0, done.stderr
    return done.stdout.strip()


def plasmids_into(folder):
    folder.mkdir(parents=True, exist_ok=True)
    for part_type, name in CANONICAL.items():
        write_genbank(synth.part_plasmid(part_type, name=name, seed=synth.seed_for(name)),
                      folder / f"{name}.gb")


@pytest.fixture
def lab(tmp_path):
    """A shared remote, and one person's clone of it with plasmids in place."""
    remote = tmp_path / "remote.git"
    git(tmp_path, "init", "--quiet", "--bare", "--initial-branch=main", str(remote))

    work = tmp_path / "mine"
    git(tmp_path, "clone", "--quiet", str(remote), str(work))
    git(work, "config", "user.email", "test@example.invalid")
    git(work, "config", "user.name", "Test")

    plasmids_into(work / "plasmids")
    git(work, "add", "-A")
    git(work, "commit", "--quiet", "-m", "the plasmids we started with")
    git(work, "push", "--quiet", "origin", "HEAD:main")

    library = Library(work / "plasmids")
    library.scan()
    baseline.write(library)

    # The baseline is shared data, so leaving it uncommitted would correctly
    # show up as work waiting to be published and make every count in these
    # tests one higher than the thing being tested.
    sync.prepare(library, work)
    git(work, "add", "-A")
    git(work, "commit", "--quiet", "-m", "agree a baseline")
    git(work, "push", "--quiet", "origin", "HEAD:main")
    return library


def colleague(tmp_path, remote_name="remote.git", as_name="theirs"):
    """A second clone, standing in for someone else in the lab."""
    work = tmp_path / as_name
    git(tmp_path, "clone", "--quiet", str(tmp_path / remote_name), str(work))
    git(work, "config", "user.email", "them@example.invalid")
    git(work, "config", "user.name", "Them")
    return work


# -------------------------------------------------------------- the toggle ---


def test_sharing_is_on_unless_someone_turns_it_off(lab):
    """A library that only fills up when people remember to share does not."""
    assert sync.enabled(lab) is True


def test_the_toggle_persists_and_stays_on_this_machine(lab):
    sync.set_enabled(lab, False)
    assert sync.enabled(lab) is False
    assert sync.enabled(Library(lab.roots[0])) is False
    assert lab.settings_path.parent.name == CACHE_DIR


# -------------------------------------------------------------- the status ---


def test_a_library_outside_git_says_so_rather_than_failing(tmp_path):
    plasmids_into(tmp_path / "loose")
    found = sync.status(Library(tmp_path / "loose"))
    assert not found.available
    assert "git repository" in found.detail


def test_a_repository_with_no_remote_says_there_is_nobody_to_share_with(tmp_path):
    work = tmp_path / "solo"
    plasmids_into(work / "plasmids")
    git(tmp_path, "init", "--quiet", "--initial-branch=main", str(work))
    found = sync.status(Library(work / "plasmids"))
    assert not found.available and "no remote" in found.detail


def test_a_fresh_clone_has_nothing_to_do(lab):
    found = sync.status(lab, fetch=True)
    assert found.available and found.pending == 0 and found.incoming == 0
    assert found.summary() == "up to date"


def test_a_new_plasmid_shows_up_as_unshared(lab):
    write_genbank(synth.part_plasmid("3", name="pNew"), lab.roots[0] / "pNew.gb")
    found = sync.status(lab)
    assert found.pending == 1
    assert "1 of yours not shared" in found.summary()


# --------------------------------------------------------------- the gate ---


def test_a_plasmid_nothing_can_classify_is_not_published(lab):
    """The gate's reason for existing: a plasmid that will not digest breaks a
    colleague's dropdown, not yours, and by then it is in everyone's history."""
    write_genbank(synth.plasmid_without_bsai(name="pJunk"), lab.roots[0] / "pJunk.gb")
    lab.scan()

    found, verdict = sync.share(lab, "add pJunk")
    assert verdict.ok is False
    assert any("pJunk" in line for line in verdict.introduced)
    assert found.pending >= 1, "it was published despite the gate"


def test_problems_the_library_already_had_do_not_block_anyone(lab):
    """A gate that demands a clean library is a gate switched off in a week."""
    write_genbank(synth.plasmid_without_bsai(name="pOld"), lab.roots[0] / "pOld.gb")
    lab.scan()
    baseline.write(lab)        # agreed as known

    write_genbank(synth.part_plasmid("4", name="pFine"), lab.roots[0] / "pFine.gb")
    lab.scan()
    _, verdict = sync.share(lab, "add pFine")
    assert verdict.ok is True


def test_force_publishes_anyway(lab):
    """The gate advises; the person at the keyboard decides."""
    write_genbank(synth.plasmid_without_bsai(name="pJunk"), lab.roots[0] / "pJunk.gb")
    lab.scan()
    found, verdict = sync.share(lab, "add pJunk anyway", force=True)
    assert verdict.ok is False and found.pending == 0


# ------------------------------------------------------- publish and receive ---


def test_sharing_reaches_the_rest_of_the_lab(lab, tmp_path):
    write_genbank(synth.part_plasmid("3", name="pMine"), lab.roots[0] / "pMine.gb")
    lab.scan()
    found, verdict = sync.share(lab, "add pMine")
    assert verdict.ok and found.pending == 0

    theirs = colleague(tmp_path)
    assert (theirs / "plasmids" / "pMine.gb").exists()


def test_a_colleagues_plasmid_arrives_on_the_next_sync(lab, tmp_path):
    theirs = colleague(tmp_path)
    plasmid = synth.part_plasmid("4", name="pTheirs")
    write_genbank(plasmid, theirs / "plasmids" / "pTheirs.gb")
    git(theirs, "add", "-A")
    git(theirs, "commit", "--quiet", "-m", "add pTheirs")
    git(theirs, "push", "--quiet", "origin", "HEAD:main")

    before = len(lab.entries)
    assert sync.status(lab, fetch=True).incoming == 1
    sync.pull(lab)
    assert len(lab.entries) == before + 1
    assert any(e.name == "pTheirs" for e in lab.entries.values())


def test_curation_travels_with_the_plasmid(lab, tmp_path):
    """The whole reason overrides moved into the shared directory. Before this,
    a colleague's plasmid arrived and their reading of it did not."""
    name = sorted(lab.entries)[0]
    lab.set_concentration(name, 84.0)
    sync.share(lab, "measured a prep")

    theirs = colleague(tmp_path)
    elsewhere = Library(theirs / "plasmids")
    elsewhere.scan()
    assert elsewhere.entries[name].conc_ng_ul == 84.0


# ------------------------------------------------- what must not be published ---


def test_saving_a_plasmid_does_not_publish_unrelated_work(lab):
    """The repository is rarely only plasmids - this one sits in an R package.
    Staging the worktree would publish a half-written script under a commit
    message about a plasmid."""
    repo = sync.repo_of(lab)
    (repo / "analysis.R").write_text("# half-finished\n", encoding="utf-8")

    write_genbank(synth.part_plasmid("3", name="pMine"), lab.roots[0] / "pMine.gb")
    lab.scan()
    sync.share(lab, "add pMine")

    assert "analysis.R" in git(repo, "status", "--porcelain"), "someone's script was published"


def test_the_index_is_never_committed(lab):
    """Keyed on mtime, so it differs per machine and would conflict on every
    pull while meaning nothing."""
    write_genbank(synth.part_plasmid("3", name="pMine"), lab.roots[0] / "pMine.gb")
    lab.scan()
    sync.share(lab, "add pMine")
    tracked = git(sync.repo_of(lab), "ls-files")
    assert CACHE_DIR not in tracked


# ----------------------------------------------------------------- setup ---


def test_preparing_teaches_the_repository_both_rules(lab):
    repo = sync.repo_of(lab)
    (repo / ".gitattributes").unlink()      # as a repository that has never
    (repo / ".gitignore").unlink()          # been shared from would look
    changed = sync.prepare(lab, repo)
    assert set(changed) == {".gitattributes", ".gitignore"}
    assert "merge=union" in (repo / ".gitattributes").read_text()
    assert f"{CACHE_DIR}/" in (repo / ".gitignore").read_text()


def test_preparing_twice_changes_nothing_the_second_time(lab):
    repo = sync.repo_of(lab)
    sync.prepare(lab, repo)
    assert sync.prepare(lab, repo) == []


def test_preparing_keeps_what_the_repository_already_said(lab):
    repo = sync.repo_of(lab)
    (repo / ".gitignore").write_text("*.Rproj.user/\n", encoding="utf-8")
    sync.prepare(lab, repo)
    assert "*.Rproj.user/" in (repo / ".gitignore").read_text()


def test_the_build_log_is_set_to_merge_rather_than_conflict(lab):
    """Two people logging builds on the same afternoon append to one file."""
    repo = sync.repo_of(lab)
    sync.prepare(lab, repo)
    assert "builds.jsonl merge=union" in (repo / ".gitattributes").read_text()


# --------------------------------------------------------------- onboarding ---


def test_cloning_a_library_gets_someone_from_nothing_to_working(lab, tmp_path):
    """The step that stops people: installed the app, no plasmids to open."""
    into = sync.clone(str(sync.repo_of(lab)), tmp_path / "newcomer")
    arrived = Library(into / "plasmids")
    arrived.scan()
    assert len(arrived.entries) == len(lab.entries)


def test_cloning_refuses_to_write_into_somebody_elses_folder(lab, tmp_path):
    busy = tmp_path / "busy"
    busy.mkdir()
    (busy / "notes.txt").write_text("mine", encoding="utf-8")
    with pytest.raises(sync.SyncError, match="not empty"):
        sync.clone(str(sync.repo_of(lab)), busy)
    assert (busy / "notes.txt").read_text() == "mine"


@pytest.mark.parametrize("refusal", [
    "ERROR: Repository not found.\nfatal: Could not read from remote repository.",
    "git@github.com: Permission denied (publickey).",
    "fatal: Authentication failed for 'https://github.com/org/library.git/'",
    "Please make sure you have the correct access rights",
])
def test_every_way_github_says_no_leads_to_the_same_advice(refusal):
    """git's own wording for this talks about public keys and shell access. A
    lab member who has never pushed anything needs to be told to set up their
    GitHub access, which is the one thing that message never says."""
    said = sync.clone_failure("git@github.com:org/library.git", refusal)
    assert "no key on your account" in said
    assert "gh auth login" in said
    assert refusal.splitlines()[0] in said, "git's own words were thrown away"


def test_an_ordinary_failure_is_passed_through_unchanged(lab, tmp_path):
    """Not every failure is an access problem, and dressing a disk error up as
    one would send somebody to go and fix their SSH keys for nothing."""
    said = sync.clone_failure("u", "fatal: destination path exists and is not an empty directory")
    assert "gh auth login" not in said
    assert said.startswith("fatal:")


def run_init(*argv) -> int:
    from ggassembler.cli import build_parser, cmd_init
    return cmd_init(build_parser().parse_args(["init", *argv]))


def test_init_leaves_a_plain_folder_ready_to_share(lab, capsys):
    """Run against a library that is already there, it is the three commands
    nobody should have to be told: ignore the index, merge the build log,
    agree a baseline."""
    repo = sync.repo_of(lab)
    (repo / ".gitattributes").unlink()
    baseline.path_for(lab).unlink()

    assert run_init(str(lab.roots[0])) == 0
    assert "merge=union" in (repo / ".gitattributes").read_text()
    assert baseline.read(lab) is not None
    assert "Start the app with" in capsys.readouterr().out


def test_init_is_safe_to_run_twice(lab):
    assert run_init(str(lab.roots[0])) == 0
    assert run_init(str(lab.roots[0])) == 0


def test_init_does_not_overwrite_an_agreed_baseline(lab):
    """Re-baselining is how the gate stops objecting, so it has to stay a
    thing somebody did on purpose - never a side effect of running setup."""
    import json
    before = json.loads(baseline.path_for(lab).read_text())
    write_genbank(synth.plasmid_without_bsai(name="pNew"), lab.roots[0] / "pNew.gb")
    run_init(str(lab.roots[0]))
    assert json.loads(baseline.path_for(lab).read_text()) == before


def test_init_on_a_folder_that_is_not_there_says_how_to_clone_one(tmp_path, capsys):
    assert run_init(str(tmp_path / "nothing-here")) == 2
    assert "--from" in capsys.readouterr().err


def test_the_starter_library_is_a_public_https_url():
    """A new user has no SSH key on their account yet, and asking them to set
    GitHub access up before they have seen the app work is the wrong order."""
    assert sync.STARTER.startswith("https://")
    assert sync.STARTER.endswith(".git")


def test_init_refuses_two_sources_at_once(capsys):
    from ggassembler.cli import build_parser, cmd_init
    args = build_parser().parse_args(["init", "x", "--starter", "--from", "u"])
    assert cmd_init(args) == 2
    assert "not both" in capsys.readouterr().err


def test_init_offers_the_starter_when_there_is_nothing_to_open(tmp_path, capsys):
    """The first thing that happens to someone who installed the app and has
    no plasmids must not be a bare error."""
    assert run_init(str(tmp_path / "nowhere")) == 2
    said = capsys.readouterr().err
    assert "--starter" in said and "--from" in said
