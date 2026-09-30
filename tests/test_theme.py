"""Two themes, and the contrast that has to hold in both.

Colour is the one part of a design that fails silently. A token that is only
defined in Paper renders as *nothing* in Slate - not wrong, absent - and a
foreground that is legible on white can be invisible on #1A1D21 while the page
still looks broadly fine in a screenshot. Neither shows up in any other test.

So: every token defined in one theme must be defined in the other, and every
text pair that carries meaning is measured against WCAG rather than eyeballed.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

WEB = Path(__file__).resolve().parent.parent / "ggassembler" / "web"
TOKENS = (WEB / "tokens.css").read_text(encoding="utf-8")

#: WCAG AA: 4.5:1 for body text, 3:1 at 24px or bold 18.66px and above.
BODY, LARGE = 4.5, 3.0


def block(selector: str) -> dict[str, str]:
    """The custom properties declared in one rule."""
    start = TOKENS.index(selector + " {")
    body = TOKENS[start : TOKENS.index("\n}", start)]
    return dict(re.findall(r"(--[a-z0-9-]+)\s*:\s*([^;]+);", body))


PAPER = block(":root")
SLATE = block(':root[data-theme="slate"]')


def resolve(name: str, theme: dict[str, str]) -> str:
    """A token's value, following one level of `var(--other)` indirection."""
    value = (theme.get(name) or PAPER.get(name, "")).strip()
    reference = re.fullmatch(r"var\((--[a-z0-9-]+)\)", value)
    return resolve(reference.group(1), theme) if reference else value


def rgb(value: str) -> tuple[float, float, float]:
    value = value.strip().lstrip("#")
    if len(value) == 3:
        value = "".join(c * 2 for c in value)
    return tuple(int(value[i : i + 2], 16) / 255 for i in (0, 2, 4))


def luminance(value: str) -> float:
    def channel(c: float) -> float:
        return c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4
    r, g, b = (channel(c) for c in rgb(value))
    return 0.2126 * r + 0.7152 * g + 0.0722 * b


def contrast(fg: str, bg: str) -> float:
    a, b = luminance(fg), luminance(bg)
    high, low = max(a, b), min(a, b)
    return (high + 0.05) / (low + 0.05)


def themes():
    return [("Paper", PAPER), ("Slate", SLATE)]


# --------------------------------------------------------------------------- #
# both themes define the same things
# --------------------------------------------------------------------------- #


def test_the_themes_were_parsed_at_all():
    """Guard the guard: an empty block would make every check below vacuous."""
    assert len(PAPER) > 40, f"only found {len(PAPER)} tokens in :root"
    assert len(SLATE) > 30, f"only found {len(SLATE)} tokens in Slate"


#: Type, geometry and the fixed header live in Paper only on purpose - they do
#: not change between themes, and repeating them in Slate would be two places
#: to edit a radius.
SHARED_ON_PURPOSE = re.compile(
    r"^--(font|r|h|col$|on-dark|header$|header-active|header-line|bg$|part-[1-8]$"
    r"|badge-[1-8]-(bg|fg)$)"
)


def test_every_colour_token_exists_in_both_themes():
    """A token missing from Slate renders as nothing, not as a fallback."""
    missing = [
        name for name in PAPER
        if name not in SLATE
        and not SHARED_ON_PURPOSE.match(name)
        and "shadow" not in name
    ]
    assert not missing, f"defined in Paper but not Slate: {missing}"


def test_slate_redefines_the_part_palette_and_the_badges():
    """These carry meaning, so they are the ones a theme must not skip."""
    for position in range(1, 9):
        assert f"--part-{position}" in SLATE
        assert f"--badge-{position}-bg" in SLATE
        assert f"--badge-{position}-fg" in SLATE


def test_both_themes_define_all_three_elevations():
    for name, theme in themes():
        for level in (1, 2, 3):
            assert f"--shadow-{level}" in theme, f"{name} has no --shadow-{level}"


def test_a_well_recedes_in_both_themes():
    """The one that is easy to get backwards.

    `--surface-sunk` is lighter than the card on Paper and *darker* on Slate. A
    well recedes on dark by going down; lightening it makes it read as raised
    and the form turns inside out.
    """
    assert luminance(resolve("--surface-sunk", PAPER)) < luminance(resolve("--surface", PAPER))
    assert luminance(resolve("--surface-sunk", SLATE)) < luminance(resolve("--surface", SLATE))


def test_a_card_sits_above_the_page_in_both_themes():
    """The recalibration: a white card on #F7F5F0 barely reads as a card."""
    for name, theme in themes():
        # a ratio, not a luminance gap: at the dark end every absolute
        # difference is tiny and a fixed gap would demand a Slate page far
        # lighter than a dark theme should be
        ratio = contrast(resolve("--surface", theme), resolve("--ground", theme))
        assert ratio > 1.05, f"{name}: cards do not separate from the page ({ratio:.3f})"


def test_slate_carries_elevation_by_lightness_not_by_shadow():
    """Shadows barely register on #121417, so raised has to be visibly lighter."""
    raised = luminance(resolve("--surface-raised", SLATE))
    surface = luminance(resolve("--surface", SLATE))
    assert raised > surface * 1.3, "a raised surface must read as raised on dark"


# --------------------------------------------------------------------------- #
# contrast
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("theme_name", ["Paper", "Slate"])
@pytest.mark.parametrize("ink", ["--ink", "--ink-2"])
@pytest.mark.parametrize("ground", ["--surface", "--surface-sunk", "--ground"])
def test_body_text_is_legible(theme_name, ink, ground):
    theme = dict(themes())[theme_name]
    ratio = contrast(resolve(ink, theme), resolve(ground, theme))
    assert ratio >= BODY, f"{theme_name}: {ink} on {ground} is {ratio:.2f}:1"


@pytest.mark.parametrize("theme_name", ["Paper", "Slate"])
@pytest.mark.parametrize("ground", ["--surface", "--surface-sunk"])
def test_the_quietest_ink_is_legible(theme_name, ground):
    """`--ink-3` is captions and meta lines - small, so it is held to 4.5:1."""
    theme = dict(themes())[theme_name]
    ratio = contrast(resolve("--ink-3", theme), resolve(ground, theme))
    assert ratio >= BODY, f"{theme_name}: --ink-3 on {ground} is {ratio:.2f}:1"


@pytest.mark.parametrize("theme_name", ["Paper", "Slate"])
@pytest.mark.parametrize("position", range(1, 9))
def test_every_badge_is_legible_on_its_own_ground(theme_name, position):
    """Eight pairs per theme, and one bad pair is one part type nobody can read."""
    theme = dict(themes())[theme_name]
    fg = resolve(f"--badge-{position}-fg", theme)
    bg = resolve(f"--badge-{position}-bg", theme)
    ratio = contrast(fg, bg)
    assert ratio >= BODY, f"{theme_name}: badge {position} is {ratio:.2f}:1"


@pytest.mark.parametrize("theme_name", ["Paper", "Slate"])
@pytest.mark.parametrize("pair", [("--ok", "--ok-bg"), ("--warn", "--warn-bg"),
                                  ("--error", "--error-bg"),
                                  ("--attention", "--attention-bg"),
                                  ("--unknown-fg", "--unknown-bg")])
def test_status_text_is_legible_on_its_own_wash(theme_name, pair):
    theme = dict(themes())[theme_name]
    fg, bg = (resolve(name, theme) for name in pair)
    ratio = contrast(fg, bg)
    assert ratio >= BODY, f"{theme_name}: {pair[0]} on {pair[1]} is {ratio:.2f}:1"


@pytest.mark.parametrize("theme_name", ["Paper", "Slate"])
def test_the_primary_button_is_legible(theme_name):
    """`--on-accent` on `--accent`, which is what `.btn-primary` renders.

    White works on Paper's rust and is 2.6:1 on Slate's peach, so this is one
    of the pairs a theme cannot simply inherit."""
    theme = dict(themes())[theme_name]
    ratio = contrast(resolve("--on-accent", theme), resolve("--accent", theme))
    assert ratio >= LARGE, f"{theme_name}: button text is {ratio:.2f}:1"


@pytest.mark.parametrize("theme_name", ["Paper", "Slate"])
def test_what_sits_on_a_part_colour_is_legible_on_every_one(theme_name):
    """`--on-part` is white on Paper's mid-tones and near-black on Slate's
    pastels. One value has to work across all eight arcs."""
    theme = dict(themes())[theme_name]
    on_part = resolve("--on-part", theme)
    worst = min(
        (contrast(on_part, resolve(f"--part-{n}", theme)), n) for n in range(1, 9)
    )
    # 3:1 rather than 4.5, and knowingly: one token serves eight arcs, and the
    # worst pair is Type 7 brown at about 3.9. What this replaced was white,
    # which was 1.75:1 on Type 3 yellow - a character nobody could read.
    assert worst[0] >= LARGE, f"{theme_name}: part {worst[1]} is {worst[0]:.2f}:1"


@pytest.mark.parametrize("theme_name", ["Paper", "Slate"])
def test_the_eight_part_colours_stay_separable(theme_name):
    """Type 7 brown sank into the ground at its Paper value and stopped being
    tellable from Type 8 grey, which is why Slate lifts it furthest."""
    theme = dict(themes())[theme_name]
    values = [luminance(resolve(f"--part-{n}", theme)) for n in range(1, 9)]
    ground = luminance(resolve("--ground", theme))
    for n, value in enumerate(values, start=1):
        assert abs(value - ground) > 0.02, f"{theme_name}: part {n} merges into the page"


# --------------------------------------------------------------------------- #
# the inline resolver
# --------------------------------------------------------------------------- #


PAGES = sorted(WEB.glob("*/index.html"))


def theme_script(html: str) -> str:
    start = html.index("<script>")
    return html[start : html.index("</script>", start)]


def test_every_screen_carries_the_resolver():
    assert len(PAGES) >= 6, "screens have moved"
    for page in PAGES:
        assert "ggasmTheme" in page.read_text(encoding="utf-8"), f"{page.parent.name} has none"


def test_the_copies_have_not_drifted():
    """It is duplicated because it must be inline - so this is what keeps the
    copies honest. Editing one and not the others is the obvious failure and
    would leave one screen on a stale theme."""
    scripts = {page: theme_script(page.read_text(encoding="utf-8")) for page in PAGES}
    first = scripts[PAGES[0]]
    drifted = [p.parent.name for p, s in scripts.items() if s != first]
    assert not drifted, f"the theme script differs on: {drifted}"


def test_the_resolver_runs_before_the_stylesheet():
    """Otherwise every load of a dark page flashes white first."""
    for page in PAGES:
        html = page.read_text(encoding="utf-8")
        # the <link>, not the word: the script's own comment says "tokens.css"
        # and matching that would compare the block against itself
        link = html.index('href="/static/tokens.css"')
        assert html.index("ggasmTheme") < link, page.parent.name


def test_every_storage_call_is_guarded():
    """localStorage throws in a private window rather than returning nothing,
    and a colour scheme is not worth a broken page."""
    script = theme_script(PAGES[0].read_text(encoding="utf-8"))
    for call in ("localStorage.getItem", "localStorage.setItem"):
        at = script.index(call)
        window = script[max(0, at - 220) : at + 120]
        assert "try {" in window and "catch" in window, f"{call} is unguarded"


def test_the_toggle_offers_three_states_and_defaults_to_auto():
    script = theme_script(PAGES[0].read_text(encoding="utf-8"))
    assert "'auto', 'light', 'dark'" in script
    assert "return 'auto'" in script, "a failed read must fall back to auto"


def test_the_toggle_is_a_labelled_button_in_the_header():
    script = theme_script(PAGES[0].read_text(encoding="utf-8"))
    assert "aria-label" in script
    assert "topbar-right" in script, "the control does not live in the header"
    assert ".theme-toggle" in TOKENS, "the control has no styling"
