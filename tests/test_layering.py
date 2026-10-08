"""The architecture test.

The three assembly levels must stay independently modifiable: changing Level 2
must never require editing a Level 1 or Level 3 file. That is a property of the
import graph, so it is checked here rather than trusted.

Python rules
    core/*       may not import levels/*, api/* or web/*
    levels/*     may import core/*, never a sibling level, never api/*
    api/*        may not import web/*; only routes_levelN.py may import levels,
                 and only levels/levelN_*
    cli.py, __init__.py   may import anything (they are the wiring)

Front-end rules
    web/<level>/ imports nothing from a sibling level directory and fetches only
                 its own /api/<level>/* plus the shared endpoints in SHARED_API
    web/         holds no shared script or stylesheet other than tokens.css
"""

from __future__ import annotations

import ast
import re
from pathlib import Path

PKG = Path(__file__).resolve().parent.parent / "ggassembler"
WEB = PKG / "web"
LEVEL_DIRS = ("level1", "level2", "level3")

#: Endpoints any screen may call. They qualify by the same test the Library
#: endpoints always passed: the router behind them imports no level, so calling
#: one cannot couple two screens to each other. Sharing a plasmid library is a
#: property of the library, not of whatever screen you happened to be on.
SHARED_API = ("/api/library", "/api/sync")


# --------------------------------------------------------------------------- #
# import extraction
# --------------------------------------------------------------------------- #


def python_modules() -> list[tuple[str, Path]]:
    """Every module in the package, as (dotted name relative to the package, path)."""
    out = []
    for path in sorted(PKG.rglob("*.py")):
        rel = path.relative_to(PKG)
        parts = list(rel.parts)
        if parts[-1] == "__init__.py":
            parts = parts[:-1]
        else:
            parts[-1] = parts[-1][: -len(".py")]
        out.append((".".join(parts), path))
    return out


def _package_of(dotted: str, path: Path) -> list[str]:
    parts = dotted.split(".") if dotted else []
    return parts if path.name == "__init__.py" else parts[:-1]


def internal_imports(dotted: str, path: Path) -> list[tuple[str, int]]:
    """Package-internal import targets in `path`, as (dotted target, line number).

    Targets are relative to the package root, e.g. ``core.enzymes``. Absolute
    (``ggassembler.core.enzymes``) and relative (``from ..core import enzymes``)
    forms both resolve to the same string.
    """
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    pkg = _package_of(dotted, path)
    targets: list[tuple[str, int]] = []

    def record(name: str, lineno: int) -> None:
        if name:
            targets.append((name, lineno))

    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                if alias.name == "ggassembler" or alias.name.startswith("ggassembler."):
                    record(alias.name[len("ggassembler.") :], node.lineno)
        elif isinstance(node, ast.ImportFrom):
            if node.level:
                # `from . import x` -> current package; each extra dot goes up one.
                base = pkg[: len(pkg) - (node.level - 1)] if node.level > 1 else pkg
                head = ".".join([*base, *(node.module.split(".") if node.module else [])])
            elif node.module == "ggassembler":
                head = ""
            elif node.module and node.module.startswith("ggassembler."):
                head = node.module[len("ggassembler.") :]
            else:
                continue
            record(head, node.lineno)
            for alias in node.names:
                record(f"{head}.{alias.name}" if head else alias.name, node.lineno)

    return targets


def layer(dotted: str) -> str:
    """``core``, ``levels``, ``api``, ``web``, or ``""`` for a top-level module."""
    head = dotted.split(".")[0] if dotted else ""
    return head if head in {"core", "levels", "api", "web"} else ""


def area_of(dotted: str) -> str | None:
    """The screen a module belongs to, from its file name.

    ``levels/level2_cassette.py`` and ``api/routes_level2.py`` are both
    ``level2``; ``levels/plate_batch.py`` and ``api/routes_plate.py`` are both
    ``plate``. Named by convention rather than by a list, so a new screen that
    follows the convention is governed by the same rule as the three assembly
    levels without this file having to learn about it.
    """
    tail = dotted.split(".")[1] if "." in dotted else ""
    if not tail:
        return None
    if tail.startswith("routes_"):
        return tail[len("routes_"):] or None
    return tail.split("_")[0] or None


#: Kept under the old name: most of this file reads better saying "level".
level_of = area_of


def violation(src: str, dst: str) -> str | None:
    """The rule broken by `src` importing `dst`, or None if the import is legal."""
    src_layer, dst_layer = layer(src), layer(dst)
    if dst_layer == "":
        return None  # the root package itself

    if src_layer == "core":
        if dst_layer != "core":
            return f"core/ must not import {dst_layer}/"

    elif src_layer == "levels":
        if dst_layer in {"api", "web"}:
            return f"levels/ must not import {dst_layer}/"
        if dst_layer == "levels" and area_of(dst) != area_of(src):
            return "levels must not import one another"

    elif src_layer == "api":
        if dst_layer == "web":
            return "api/ must not import web/"
        if dst_layer == "levels":
            own = area_of(src)
            if own is None:
                return f"only a routes_<screen>.py may import levels/ (not {src})"
            if dst == "levels":
                # the bare package, from `from ..levels import levelN_x`; the
                # submodule that import also records is what gets judged
                return None
            if area_of(dst) != own:
                return f"api/routes_{own}.py may import only levels/{own}_*"

    return None


# --------------------------------------------------------------------------- #
# tests
# --------------------------------------------------------------------------- #


def test_package_layout_present():
    assert (PKG / "core").is_dir()
    assert (PKG / "levels").is_dir()
    assert (PKG / "api").is_dir()
    assert WEB.is_dir()
    assert python_modules(), "no modules found - is the package path right?"


def test_python_import_boundaries():
    offences = []
    for dotted, path in python_modules():
        for target, lineno in internal_imports(dotted, path):
            why = violation(dotted, target)
            if why:
                rel = path.relative_to(PKG.parent)
                offences.append(f"{rel}:{lineno}: imports {target} - {why}")
    assert not offences, "import-boundary violations:\n" + "\n".join(offences)


def test_web_levels_are_self_contained():
    """No level's front end may reach into a sibling level's directory."""
    offences = []
    import_re = re.compile(r"""(?:from|import)\s+['"]([^'"]+)['"]""")
    for own in LEVEL_DIRS:
        for path in sorted((WEB / own).rglob("*.js")):
            text = path.read_text(encoding="utf-8")
            for spec in import_re.findall(text):
                for other in LEVEL_DIRS:
                    if other != own and re.search(rf"\b{other}\b", spec):
                        rel = path.relative_to(WEB)
                        offences.append(f"web/{rel}: imports {spec} from sibling {other}/")
    assert not offences, "front-end level coupling:\n" + "\n".join(offences)


def test_the_library_screen_uses_only_shared_endpoints():
    """The Library screen is not a level: it owns no /api/levelN endpoints."""
    offences = []
    api_re = re.compile(r"""['"`](/api/[A-Za-z0-9_\-/{}$.]*)""")
    for path in sorted((WEB / "library").rglob("*.js")):
        for url in api_re.findall(path.read_text(encoding="utf-8")):
            if not url.startswith(SHARED_API):
                offences.append(f"web/library/{path.name}: calls {url}, not a shared endpoint")
    assert not offences, "library screen reached into a level:\n" + "\n".join(offences)


def test_the_sequencing_screen_uses_only_its_own_and_shared_endpoints():
    """Not a level either: it checks clones, it does not design anything.

    It may reach `/api/library*`, which is the shared surface every screen is
    allowed - that is how it wears the same search box as the others - but
    never a level's own endpoints.
    """
    offences = []
    api_re = re.compile(r"""['"`](/api/[A-Za-z0-9_\-/{}$.]*)""")
    for path in sorted((WEB / "sequencing").rglob("*.js")):
        for url in api_re.findall(path.read_text(encoding="utf-8")):
            if url.startswith(("/api/sequencing", *SHARED_API)):
                continue
            offences.append(f"web/sequencing/{path.name}: calls {url}")
    assert not offences, "sequencing screen reached outside its own API:\n" + "\n".join(offences)


def test_the_sequencing_finder_withholds_nothing():
    """Every other screen scopes the search; this one must not.

    On the Cassette and Multigene screens a hit has to fit a slot or close a
    chain, so `usable=` withholds the ones that cannot. Here a reference is
    compared against rather than assembled with, and anything with a sequence
    is a legitimate thing to have sequenced - including the finished construct
    those screens deliberately hide. Adding a scope here would hide the very
    plasmid the person sent off, so its absence is asserted rather than left
    to be tidied up by someone matching the other screens.
    """
    script = (WEB / "sequencing" / "sequencing.js").read_text(encoding="utf-8")
    # the URL, not the word: the module explains this choice in a comment, and
    # a test that searched the whole file would match its own explanation -
    # which has caught me out on this repo before
    urls = re.findall(r"/api/library/search\?[^`'\"\s]*", script)
    assert urls, "the sequencing screen lost its search box"
    assert not any("usable" in url for url in urls), \
        f"the sequencing search must not withhold anything: {urls}"


def test_the_alignment_viewer_does_not_render_every_column_at_once():
    """A 12 kb construct across eight clones is over a hundred thousand cells.

    Guarded rather than measured, like the library table: the obvious edit -
    drawing the whole width once and letting the browser scroll it - looks
    simpler and is what turns the screen into a frozen tab. Two things have to
    survive: a window on the columns, and runs instead of one span per base.
    """
    viewer = (WEB / "sequencing" / "sequencing.js").read_text(encoding="utf-8")
    columns = (WEB / "sequencing" / "columns.js").read_text(encoding="utf-8")
    assert "scrollLeft" in viewer, "the whole alignment is not windowed on scroll"
    assert re.search(r"viewer\.onscroll\s*=\s*paint", viewer), "nothing repaints on scroll"
    assert "clientWidth" in viewer, "the window is not sized to what is on screen"
    assert "export function runs" in columns, "cells are no longer merged into runs"


def test_web_levels_fetch_only_their_own_endpoints():
    offences = []
    api_re = re.compile(r"""['"`](/api/[A-Za-z0-9_\-/{}$.]*)""")
    for own in LEVEL_DIRS:
        for path in sorted((WEB / own).rglob("*.js")):
            for url in api_re.findall(path.read_text(encoding="utf-8")):
                if url.startswith((*SHARED_API, f"/api/{own}")):
                    continue
                rel = path.relative_to(WEB)
                offences.append(f"web/{rel}: calls {url}, outside /api/{own} and /api/library")
    assert not offences, "front-end endpoint violations:\n" + "\n".join(offences)


def test_rule_checker_catches_known_violations():
    """Guard against a checker that passes because it checks nothing."""
    assert violation("core.assembly", "levels.level2_cassette")
    assert violation("core.library", "api.routes_library")
    assert violation("levels.level2_cassette", "levels.level1_part")
    assert violation("levels.level3_multigene", "api.main")
    assert violation("api.routes_level2", "levels.level3_multigene")
    assert violation("core.library", "levels")
    assert violation("api.main", "levels.level2_cassette")

    assert violation("levels.level2_cassette", "core.assembly") is None
    assert violation("api.routes_level2", "levels.level2_cassette") is None
    # a screen that is not one of the three assembly levels is governed the
    # same way, by name, with no list of screens anywhere in this file
    assert violation("api.routes_plate", "levels.plate_batch") is None
    assert violation("api.routes_plate", "levels.level2_cassette")
    assert violation("levels.plate_batch", "levels.level3_multigene")
    assert violation("api.routes_library", "levels.plate_batch")
    assert violation("api.routes_level2", "levels") is None  # `from ..levels import ...`
    assert violation("api.routes_level2", "core.parttypes") is None
    assert violation("core.assembly", "core.enzymes") is None
    assert violation("cli", "api.main") is None


def test_relative_and_absolute_imports_resolve_alike(tmp_path):
    src = "from ..core import assembly\nfrom ggassembler.core import assembly as a2\n"
    path = tmp_path / "routes_level2.py"
    path.write_text(src, encoding="utf-8")
    targets = {t for t, _ in internal_imports("api.routes_level2", path)}
    assert "core" in targets and "core.assembly" in targets


def test_no_screen_takes_a_badge_colour_from_the_api():
    """Colour is resolved from the part type in CSS, never sent as hex.

    It used to be sent: every payload carried `color`, `badge_bg` and
    `badge_fg`. That is invisible until the page can be themed, and then it
    fails in the worst way - a theme switch restyles everything except the ring
    arcs and the type badges, which are the only places colour carries meaning.

    The matching half of this - that no endpoint sends a literal colour at all -
    is checked against live responses in `test_api.py`, because whether a value
    is a token or a hex code is not something a grep over the front end can
    see. What it *can* see is the two field names that used to carry hex.
    """
    offences = []
    for path in sorted(WEB.rglob("*.js")):
        source = path.read_text(encoding="utf-8")
        for field in ("badge_bg", "badge_fg"):
            if field in source:
                offences.append(f"web/{path.relative_to(WEB)}: reads {field} from a payload")
    assert not offences, "colour taken from the API:\n" + "\n".join(offences)


def test_the_part_palette_is_mapped_in_the_one_shared_file():
    """Guard the other half: the mapping has to exist for the above to work."""
    tokens = (WEB / "tokens.css").read_text(encoding="utf-8")
    for position in range(1, 9):
        assert f'[data-part^="{position}"]' in tokens, f"no rule for part {position}"
        assert f"--part-{position})" in tokens
    assert "--part-color" in tokens and "--badge-bg" in tokens and "--badge-fg" in tokens


def test_tokens_css_is_the_only_shared_front_end_file():
    shared = [
        p.name
        for p in WEB.iterdir()
        if p.is_file() and p.suffix in {".js", ".css", ".mjs"} and p.name != "tokens.css"
    ]
    assert not shared, f"web/ may share only tokens.css, found: {sorted(shared)}"


def test_the_finder_is_wired_on_the_screens_that_offer_it():
    """A search box in the markup with no script behind it is a dead control.

    This has bitten before: an edit that silently missed the file leaves a
    screen looking finished and doing nothing.
    """
    offences = []
    for own in ("level2", "level3", "sequencing"):
        html = (WEB / own / "index.html").read_text(encoding="utf-8")
        if 'id="finder-input"' not in html:
            offences.append(f"web/{own}/index.html: no search box")
            continue
        script = (WEB / own / f"{own}.js").read_text(encoding="utf-8")
        if "/api/library/search" not in script:
            offences.append(f"web/{own}/{own}.js: never calls the search endpoint")
        # the call, not the declaration - `function wireFinder()` proves nothing
        if not re.search(r"^\s*wireFinder\(\);", script, re.MULTILINE):
            offences.append(f"web/{own}/{own}.js: defines the finder but never calls wireFinder()")
    assert not offences, "search box not wired:\n" + "\n".join(offences)


def test_the_finder_asks_before_it_changes_the_design():
    """Picking a result must offer a confirm step, not silently fill a slot.

    The two design screens only. The Sequencing screen is deliberately not in
    this list: there a hit selects a reference to compare against, which is one
    click to change and overwrites no work, so asking first would be friction
    guarding against nothing.
    """
    offences = []
    for own in ("level2", "level3"):
        script = (WEB / own / f"{own}.js").read_text(encoding="utf-8")
        if "finder-confirm" not in script:
            offences.append(f"web/{own}/{own}.js: no confirm step before applying a hit")
    assert not offences, "\n".join(offences)


def test_each_screen_only_searches_for_what_it_can_use():
    """The cassette screen must not offer a finished construct as a part.

    A 10.8 kb multigene assembly annotated `L6_CYC1_type4` was offered as a
    terminator; the fix is the `usable` scope on the shared search, and a
    screen that forgets to pass it silently regains the bug.
    """
    wanted = {"level2": "usable=part", "level3": "usable=multigene"}
    offences = []
    for own, scope in wanted.items():
        script = (WEB / own / f"{own}.js").read_text(encoding="utf-8")
        if scope not in script:
            offences.append(f"web/{own}/{own}.js: searches without {scope}")
    assert not offences, "unscoped search:\n" + "\n".join(offences)


def test_no_button_in_the_markup_is_left_without_a_handler():
    """A disabled-toggle with no click handler looks finished and does nothing.

    The Protocol button shipped that way: `disabled` was set from the assembly
    result, so it lit up on a valid design and then did nothing when pressed.
    """
    offences = []
    for directory in sorted(p for p in WEB.iterdir() if p.is_dir()):
        html = directory / "index.html"
        scripts = sorted(directory.glob("*.js"))
        if not html.exists() or not scripts:
            continue
        # Every script in the screen's folder, not just the one named after it.
        # A screen may split its code across files - the handler for a button
        # is as likely to be in the module that owns that feature.
        source = "\n".join(p.read_text(encoding="utf-8") for p in scripts)
        for button_id in re.findall(r'<button[^>]*\bid="([^"]+)"', html.read_text("utf-8")):
            # The id has to appear *near* an addEventListener. Merely being
            # mentioned is what the dead Protocol button already managed:
            # `el('protocol-btn').disabled = !result.ok` names it and wires
            # nothing, so a test that only looked for the name passed.
            # A wider window than it looks: ids are often wired through a
            # table and a loop, so the listener is several statements below
            # the name. Verified to still catch a button with no handler at
            # all, which is the failure this exists for.
            near = rf"""['"]{re.escape(button_id)}['"](?:.|\n){{0,420}}?addEventListener"""
            if re.search(near, source, re.DOTALL):
                continue
            offences.append(f"web/{directory.name}/index.html: #{button_id} has no handler")
    assert not offences, "dead controls:\n" + "\n".join(offences)


def test_every_front_end_script_parses():
    """A JS syntax error takes a whole screen down and no Python test sees it.

    Skipped where node is unavailable rather than failing: this guards against
    a typo, and a missing toolchain is not one.
    """
    import shutil
    import subprocess

    node = shutil.which("node")
    if not node:
        import pytest

        pytest.skip("node is not installed")

    import tempfile

    broken = []
    with tempfile.TemporaryDirectory() as tmp:
        for path in sorted(WEB.rglob("*.js")):
            # the suffix has to be .mjs, or node parses it as CommonJS and
            # every `import` reads as a syntax error
            copy = Path(tmp) / f"{path.stem}.mjs"
            copy.write_bytes(path.read_bytes())
            result = subprocess.run(
                [node, "--check", str(copy)], capture_output=True, text=True, timeout=30
            )
            if result.returncode:
                detail = next(
                    (line for line in result.stderr.splitlines() if "Error" in line),
                    result.stderr.strip()[:200],
                )
                broken.append(f"{path.relative_to(WEB)}: {detail.strip()}")
    assert not broken, "front-end syntax errors:\n" + "\n".join(broken)


def test_the_library_table_does_not_render_every_row_at_once():
    """622 rows of nine cells rebuilt on each keystroke is the filter's own cost.

    Guarded rather than measured: the check is that the chunking is still
    there, since the obvious edit - going back to `rows.map(renderRow)` - looks
    harmless and is not.
    """
    source = (WEB / "library" / "library.js").read_text(encoding="utf-8")
    assert "CHUNK" in source, "the table renders every matching row at once"
    assert "slice(0, shownCount)" in source, "no window on the rendered rows"
    assert re.search(r"setTimeout\(renderRows", source), "the filter is not debounced"


def _without_comments(source: str) -> str:
    """JavaScript with /* block */ and // line comments removed."""
    source = re.sub(r"/\*.*?\*/", "", source, flags=re.DOTALL)
    return re.sub(r"^\s*//.*$", "", source, flags=re.MULTILINE)


def test_every_save_endpoint_is_reachable_from_its_screen():
    """An endpoint nothing calls is a feature that does not exist.

    /api/level2/save wrote the cassette into the library, was tested, and was
    wired to no button at all: the Cassette screen's "Save" saved a *design*
    through /api/library/designs, so pressing it looked like it had added a
    plasmid and had not. Level 1 and Level 3 both said "Save to library" and
    did it. This checks each screen reaches its own save.
    """
    offences = []
    for level in LEVEL_DIRS:
        source = "\n".join(_without_comments(p.read_text(encoding="utf-8"))
                           for p in sorted((WEB / level).glob("*.js")))
        # Comments stripped first: the comment explaining this bug names the
        # endpoint, and a test that reads its own explanation as evidence
        # passes while the button is still dead.
        if f"/api/{level}/save" not in source:
            offences.append(f"web/{level}/ never calls /api/{level}/save")
    assert not offences, "unreachable save endpoints:\n" + "\n".join(offences)


def test_saving_a_design_and_saving_a_plasmid_are_different_controls():
    """They were one button on Level 2, and the one it did was the lesser."""
    html = (WEB / "level2" / "index.html").read_text(encoding="utf-8")
    assert 'id="save-design-btn"' in html and 'id="save-btn"' in html
    assert "Save to library" in html, "nothing says where a construct goes"


def test_saving_to_the_library_asks_for_a_name():
    """The name field holds a default - pCassette, new_part - until somebody
    changes it, and nothing in the app removes a plasmid saved under the wrong
    one. The first construct saved on this project was called pCassette, and
    the fix was three git commands in a terminal.
    """
    offences = []
    for level in LEVEL_DIRS:
        source = _without_comments(
            (WEB / level / f"{level}.js").read_text(encoding="utf-8"))
        block = source[source.index("el('save-btn').addEventListener"):]
        block = block[:1200]
        if "prompt(" not in block:
            offences.append(f"web/{level}/ saves without asking for a name")
    assert not offences, "\n".join(offences)
