"""The split between state that syncs and state that must not.

A library folder holds two unrelated things. One is derived - an index keyed on
mtimes, rebuildable from the files in seconds, and different on every machine.
The other is what people decided: that this plasmid is really a Type 3, that
this prep measured at 84 ng/uL, that this design is worth keeping. Nothing can
recompute the second kind.

Before these tests they lived in the same hidden directory, which meant the
gitignore rule covering the index also covered the curation. A colleague's
plasmid file would reach you and their reading of it would not.
"""

from __future__ import annotations

import json

import pytest

from ggassembler.core.library import CACHE_DIR, DATA_DIR, Library
from ggassembler.core.seqio import write_genbank

from . import synth
from .test_assembly import CANONICAL


@pytest.fixture
def folder(tmp_path):
    for part_type, name in CANONICAL.items():
        write_genbank(synth.part_plasmid(part_type, name=name, seed=synth.seed_for(name)),
                      tmp_path / f"{name}.gb")
    return tmp_path


def a_library(folder):
    library = Library(folder)
    library.scan()
    return library


def first_part(library):
    return next(iter(sorted(library.entries)))


# ------------------------------------------------------- which side it is on ---


def test_the_index_is_derived_and_stays_out_of_the_shared_directory(folder):
    """Keyed on mtime, so it differs on every machine. Shared, it would
    conflict on every single pull."""
    library = a_library(folder)
    assert library.index_path.parent.name == CACHE_DIR
    assert library.index_path.exists()
    assert not (folder / DATA_DIR / "index.json").exists()


def test_curation_lands_in_the_shared_directory(folder):
    """The whole point: a colleague's reading of a plasmid travels with it."""
    library = a_library(folder)
    library.set_concentration(first_part(library), 84.0)
    written = list((folder / DATA_DIR).rglob("*.json"))
    assert written, "the override was not written anywhere shareable"


def test_a_persons_own_settings_are_never_shared(folder):
    """Whether this person syncs is theirs. In the shared directory, one person
    switching sync off would switch it off for everyone."""
    library = a_library(folder)
    assert library.settings_path.parent.name == CACHE_DIR
    assert DATA_DIR not in library.settings_path.parts


def test_predicted_constructs_are_not_indexed_as_plasmids(folder):
    """`expected/` moved out of the hidden directory, so its .gb files now sit
    under a scanned root. Without DATA_DIR in SKIP_DIRS the library indexes its
    own predictions as though someone had cloned them."""
    library = a_library(folder)
    before = len(library.entries)
    library.remember_expected("pPredicted", synth.part_plasmid("3", name="pPredicted"))
    assert library.expected_dir.is_dir(), "nothing was written, so this proves nothing"

    again = a_library(folder)
    assert len(again.entries) == before, "a predicted construct was indexed as a plasmid"


# ------------------------------------------------------------- merge shape ---


def test_two_people_curating_two_plasmids_touch_two_files(folder):
    """The merge property, stated as a property of the filesystem. One shared
    dict would put both edits on the same line of the same file."""
    library = a_library(folder)
    names = sorted(library.entries)[:2]
    library.set_concentration(names[0], 10.0)
    library.set_concentration(names[1], 20.0)

    files = {p.read_text(encoding="utf-8") for p in library.overrides_dir.rglob("*.json")}
    assert len(files) == 2, "both overrides went into one file; git will conflict on it"


def test_an_override_file_is_named_after_its_plasmid(folder):
    """So `git log` on one file is the history of what this plasmid was thought
    to be. A hash would be unique and tell a reviewer nothing."""
    library = a_library(folder)
    name = first_part(library)
    library.set_concentration(name, 42.0)
    assert library._override_path(name).exists()
    assert name.split("/")[-1] in library._override_path(name).name


def test_clearing_an_override_removes_its_file(folder):
    """Left behind, an empty record is a file in everyone's next pull saying
    nothing."""
    library = a_library(folder)
    name = first_part(library)
    library.set_concentration(name, 42.0)
    path = library._override_path(name)
    library.set_concentration(name, None)
    assert not path.exists()


def test_an_override_knows_its_own_plasmid(folder):
    """The path is stored in the file, not only encoded in the filename, so
    reading back never depends on reversing the slug."""
    library = a_library(folder)
    name = first_part(library)
    library.set_concentration(name, 42.0)
    assert json.loads(library._override_path(name).read_text())["path"] == name


def test_overrides_survive_a_reload(folder):
    library = a_library(folder)
    name = first_part(library)
    library.set_concentration(name, 42.0)
    assert a_library(folder).entries[name].conc_ng_ul == 42.0


# ---------------------------------------------------------------- designs ---


def test_each_design_is_its_own_file(folder):
    library = a_library(folder)
    library.save_design("level2", "first", {"selections": {}})
    library.save_design("level2", "second", {"selections": {}})
    assert len(list(library.designs_dir.rglob("*.json"))) == 2


def test_saving_a_design_twice_replaces_it(folder):
    library = a_library(folder)
    library.save_design("level2", "same", {"selections": {"1": "a"}})
    library.save_design("level2", "same", {"selections": {"1": "b"}})
    saved = library.designs("level2")
    assert len(saved) == 1 and saved[0]["design"]["selections"]["1"] == "b"


def test_deleting_a_design_removes_the_file_and_says_so(folder):
    library = a_library(folder)
    library.save_design("level2", "doomed", {})
    assert library.delete_design("level2", "doomed") is True
    assert library.delete_design("level2", "doomed") is False
    assert library.designs("level2") == []


def test_two_levels_may_use_the_same_design_name(folder):
    library = a_library(folder)
    library.save_design("level1", "shared name", {})
    library.save_design("level2", "shared name", {})
    assert len(library.designs()) == 2


# -------------------------------------------------------------- migration ---


def test_an_old_overrides_file_is_split_into_one_file_each(folder):
    """Nobody should lose their curation to a change in storage layout."""
    library = Library(folder)
    library.cache_dir.mkdir(parents=True, exist_ok=True)
    legacy = library.legacy_overrides_path
    legacy.write_text(json.dumps({
        "pOne.gb": {"conc_ng_ul": 11.0},
        "pTwo.gb": {"part_type": "3"},
    }), encoding="utf-8")

    library.scan()
    assert library.overrides["pOne.gb"]["conc_ng_ul"] == 11.0
    assert library.overrides["pTwo.gb"]["part_type"] == "3"
    assert len(list(library.overrides_dir.rglob("*.json"))) == 2


def test_the_old_file_is_kept_rather_than_deleted(folder):
    """It holds concentrations measured at a bench, which nothing can
    recompute. A migration that loses those is worse than a stray file."""
    library = Library(folder)
    library.cache_dir.mkdir(parents=True, exist_ok=True)
    legacy = library.legacy_overrides_path
    legacy.write_text(json.dumps({"pOne.gb": {"conc_ng_ul": 11.0}}), encoding="utf-8")

    library.scan()
    assert not legacy.exists()
    assert legacy.with_suffix(".json.migrated").exists()


def test_migrating_twice_does_not_undo_an_edit_made_in_between(folder):
    """The second scan must not resurrect the old value over the new one."""
    library = Library(folder)
    library.cache_dir.mkdir(parents=True, exist_ok=True)
    library.legacy_overrides_path.write_text(
        json.dumps({"pOne.gb": {"conc_ng_ul": 11.0}}), encoding="utf-8")
    library.scan()

    library.set_concentration("pOne.gb", 99.0)
    again = a_library(folder)
    assert again.overrides["pOne.gb"]["conc_ng_ul"] == 99.0


def test_old_designs_migrate_too(folder):
    library = Library(folder)
    library.cache_dir.mkdir(parents=True, exist_ok=True)
    library.legacy_designs_path.write_text(json.dumps({"designs": [
        {"level": "level2", "name": "old one", "saved_at": "2026-01-01T00:00:00", "design": {}},
    ]}), encoding="utf-8")

    assert [d["name"] for d in library.designs("level2")] == ["old one"]
    assert library.legacy_designs_path.with_suffix(".json.migrated").exists()
