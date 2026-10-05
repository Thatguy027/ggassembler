"""The sharing endpoints, and the control on the Library screen.

The test that matters most here is the dullest one: that saving a plasmid still
works when the library is not in a git repository at all. Sharing is a feature
some people will never switch on, and a save that fails because a push could
not be attempted would make this change strictly worse than not making it.
"""

from __future__ import annotations

import pytest
from fastapi.testclient import TestClient

from ggassembler.api.main import create_app
from ggassembler.core.library import CACHE_DIR, Library
from ggassembler.core.seqio import write_genbank

from . import synth
from .test_assembly import CANONICAL


@pytest.fixture
def client(tmp_path):
    for part_type, name in CANONICAL.items():
        write_genbank(synth.part_plasmid(part_type, name=name, seed=synth.seed_for(name)),
                      tmp_path / f"{name}.gb")
    with TestClient(create_app(tmp_path)) as client:
        client.library_dir = tmp_path
        yield client


def design(**overrides):
    body = {"selections": dict(CANONICAL), "name": "pCassette"}
    body.update(overrides)
    return body


# ------------------------------------------------------------- the default ---


def test_sharing_is_on_out_of_the_box(client):
    """A lab library that fills up only when people remember to share does not."""
    assert client.get("/api/sync").json()["enabled"] is True


def test_a_library_outside_git_is_honest_about_it(client):
    body = client.get("/api/sync").json()
    assert body["available"] is False
    assert body["detail"] and body["summary"] == body["detail"]


def test_the_toggle_persists(client):
    assert client.post("/api/sync/enabled", json={"on": False}).json()["enabled"] is False
    assert client.get("/api/sync").json()["enabled"] is False
    assert client.post("/api/sync/enabled", json={"on": True}).json()["enabled"] is True


def test_the_toggle_is_this_persons_business_and_not_the_labs(client):
    """In the shared directory it would switch sharing off for everyone, which
    is invisible until somebody's month of work turns out never to have left
    their laptop."""
    client.post("/api/sync/enabled", json={"on": False})
    library = Library(client.library_dir)
    assert library.settings_path.parent.name == CACHE_DIR
    assert library.settings_path.exists()


# --------------------------------------------- it must never break a save ---


def test_saving_works_when_there_is_nothing_to_share_with(client):
    """The dull test that matters: sharing is optional, saving is not."""
    response = client.post("/api/level2/save", json=design())
    assert response.status_code == 200, response.text
    assert response.json()["ok"] is True
    assert (client.library_dir / "pCassette.gb").exists()


def test_saving_works_with_sharing_switched_off(client):
    client.post("/api/sync/enabled", json={"on": False})
    assert client.post("/api/level2/save", json=design(name="pQuiet")).json()["ok"] is True


def test_a_failed_share_is_remembered_rather_than_lost(client):
    """A push that quietly fails while you believe your work is in the library
    is the worst version of this feature."""
    client.post("/api/level2/save", json=design(name="pNoted"))
    last = client.get("/api/sync").json()["last"]
    assert last is None or last["ok"] is False


# ---------------------------------------------------------- not set up yet ---


def test_pulling_without_a_remote_explains_itself(client):
    response = client.post("/api/sync/pull")
    assert response.status_code == 409
    assert "git" in response.json()["detail"] or "remote" in response.json()["detail"]


def test_sharing_without_a_remote_explains_itself(client):
    response = client.post("/api/sync/share", json={"message": "x"})
    assert response.status_code == 409
    assert response.json()["detail"]


def test_a_baseline_can_be_agreed_from_the_app(client):
    body = client.post("/api/sync/baseline").json()
    assert body["total"] == len(CANONICAL)
    assert (client.library_dir / "ggasm" / "baseline.json").exists()


# ------------------------------------------------------------- the control ---


def test_the_toggle_is_checked_in_the_markup(client):
    """It defaults on, and the markup has to agree with the server or the
    switch flips itself the first time the page loads."""
    html = client.get("/library").text
    toggle = html[html.index('id="sync-toggle"') - 60:html.index('id="sync-toggle"') + 40]
    assert "checked" in toggle


def test_the_control_says_what_it_shares_with(client):
    html = client.get("/library").text
    assert "Share my plasmids with the lab" in html


def test_colour_is_never_the_only_signal(client):
    """Someone who cannot tell amber from green still has to be able to read
    what state sharing is in."""
    script = client.get("/static/library/sync.js").text
    assert "sync-summary" in script, "the dot changes but nothing says so in words"


def test_the_screen_rereads_the_library_after_a_sync(client):
    """A sync pulls in what colleagues added; without a reload the table shows
    a library one sync out of date."""
    script = client.get("/static/library/library.js").text
    assert "startSync(load)" in script


def test_the_guide_no_longer_points_at_the_old_storage(client):
    """It named .ggasm/overrides.json, which is now neither the path nor the
    point - overrides are shared, and that is worth saying."""
    html = client.get("/library").text
    assert ".ggasm/overrides.json" not in html
    assert "ggasm/overrides/" in html
