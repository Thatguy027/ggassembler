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
                 its own /api/<level>/* plus the shared /api/library* endpoints
    web/         holds no shared script or stylesheet other than tokens.css
"""

from __future__ import annotations

import ast
import re
from pathlib import Path

PKG = Path(__file__).resolve().parent.parent / "ggassembler"
WEB = PKG / "web"
LEVEL_DIRS = ("level1", "level2", "level3")


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


def level_of(dotted: str) -> str | None:
    """``level1`` for ``levels.level1_part`` or ``api.routes_level1``, else None."""
    m = re.search(r"level([123])", dotted.split(".")[1] if "." in dotted else "")
    return f"level{m.group(1)}" if m else None


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
        if dst_layer == "levels" and level_of(dst) != level_of(src):
            return "levels must not import one another"

    elif src_layer == "api":
        if dst_layer == "web":
            return "api/ must not import web/"
        if dst_layer == "levels":
            own = level_of(src)
            if own is None:
                return f"only api/routes_levelN.py may import levels/ (not {src})"
            if dst == "levels":
                # the bare package, from `from ..levels import levelN_x`; the
                # submodule that import also records is what gets judged
                return None
            if level_of(dst) != own:
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
            if not url.startswith("/api/library"):
                offences.append(f"web/library/{path.name}: calls {url}, not a shared endpoint")
    assert not offences, "library screen reached into a level:\n" + "\n".join(offences)


def test_web_levels_fetch_only_their_own_endpoints():
    offences = []
    api_re = re.compile(r"""['"`](/api/[A-Za-z0-9_\-/{}$.]*)""")
    for own in LEVEL_DIRS:
        for path in sorted((WEB / own).rglob("*.js")):
            for url in api_re.findall(path.read_text(encoding="utf-8")):
                if url.startswith("/api/library") or url.startswith(f"/api/{own}"):
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
    for own in ("level2", "level3"):
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
    """Picking a result must offer a confirm step, not silently fill a slot."""
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
        script = directory / f"{directory.name}.js"
        if not html.exists() or not script.exists():
            continue
        source = script.read_text(encoding="utf-8")
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
