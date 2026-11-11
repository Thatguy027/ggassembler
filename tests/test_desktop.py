"""The native window: everything about it that does not need a window.

The window itself is not testable without a display, and the parts worth
testing are not the window anyway - they are what a double-clicked app has to
work out for itself, because it is handed no arguments. Choosing a port,
remembering the library, and noticing that the remembered one is gone.
"""

from __future__ import annotations

import json
import socket

import pytest

from ggassembler import desktop


@pytest.fixture
def settings(tmp_path, monkeypatch):
    """Point the remembered-library file somewhere disposable."""
    path = tmp_path / "state" / "settings.json"
    monkeypatch.setattr(desktop, "SETTINGS", path)
    return path


def test_free_port_is_actually_free():
    port = desktop.free_port()
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as probe:
        probe.bind(("127.0.0.1", port))


def test_free_port_differs_between_windows():
    """Two copies of the app open two windows rather than one error.

    ``ggasm serve`` refuses a taken port, which is right for something a
    browser is pointed at by hand. A window carries its own address.
    """
    assert desktop.free_port() != desktop.free_port()


def test_nothing_remembered_at_first(settings):
    assert desktop.remembered() == []


def test_a_library_is_remembered(settings, tmp_path):
    folder = tmp_path / "plasmids"
    folder.mkdir()

    desktop.remember([folder])

    assert desktop.remembered() == [folder.resolve()]
    assert json.loads(settings.read_text())["folders"] == [str(folder.resolve())]


def test_several_roots_are_remembered_together(settings, tmp_path):
    first, second = tmp_path / "ytk", tmp_path / "mylab"
    first.mkdir()
    second.mkdir()

    desktop.remember([first, second])

    assert desktop.remembered() == [first.resolve(), second.resolve()]


def test_a_library_that_moved_is_not_remembered(settings, tmp_path):
    """A folder that is gone has to reach the picker, not an error dialog.

    Someone whose library moved needs to be asked again, and that is exactly
    what an empty answer here produces.
    """
    folder = tmp_path / "plasmids"
    folder.mkdir()
    desktop.remember([folder])
    folder.rmdir()

    assert desktop.remembered() == []


def test_one_missing_root_discards_the_set(settings, tmp_path):
    """Half a library is the wrong library.

    Serving only the roots that survived would index a subset and say nothing
    about it - and the screens would quietly offer fewer parts than the design
    they are reloading was built from.
    """
    first, second = tmp_path / "ytk", tmp_path / "mylab"
    first.mkdir()
    second.mkdir()
    desktop.remember([first, second])
    second.rmdir()

    assert desktop.remembered() == []


def test_unreadable_settings_are_not_fatal(settings):
    settings.parent.mkdir(parents=True)
    settings.write_text("{ this is not json")

    assert desktop.remembered() == []


def test_remembering_survives_an_unwritable_home(monkeypatch, tmp_path):
    """Not being able to remember costs the next launch, not this one."""
    monkeypatch.setattr(desktop, "SETTINGS", tmp_path / "nope" / "settings.json")
    monkeypatch.setattr(
        desktop.Path, "mkdir",
        lambda *a, **k: (_ for _ in ()).throw(OSError("read-only")),
    )

    desktop.remember([tmp_path])  # must not raise
