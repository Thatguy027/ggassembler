"""Draw the app icon: an eight-part cassette, the thing the app is for.

Run when the colours or the shape change:

    python packaging/make_icon.py

Writes ``icon.icns`` (macOS) and ``icon.ico`` (Windows) next to this file, and
both are committed, so an ordinary build needs neither Pillow nor this script.

The ring is the same object the Cassette screen draws - eight arcs, one per YTK
position, in the same eight colours from ``tokens.css`` and in the same order.
Someone who has used the app for an afternoon already knows what it means,
which is more than a letter in a circle would say.

Drawn rather than rendered from SVG because macOS ships no SVG rasteriser and
the one honest alternative is a headless browser in CI. Pillow is a build-time
dependency that does the whole job in thirty lines.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

from PIL import Image, ImageDraw

HERE = Path(__file__).resolve().parent

#: The light theme's part arcs, in position order. Light rather than dark
#: because an icon sits on the user's wallpaper and a dock of its own, not on
#: the app's background.
PARTS = [
    "#6E9FC0", "#8FB861", "#D9C45A", "#D97F92",
    "#8E73C4", "#D9AB5E", "#9E6A3C", "#8C8A86",
]

#: The app's dark surface, as the tile the ring sits on.
TILE = "#1A1D21"

#: The master everything is resampled from.
CANVAS = 1024

#: Drawn this many times larger and reduced. Pillow has no antialiased arc, so
#: every edge in the drawing is a staircase and the reduction is the only
#: thing that smooths it - at 1:1 the segment boundaries come out visibly
#: ragged. It is also why nothing here is a hairline: a 1px gap at 1024
#: disappears entirely at 16.
OVERSAMPLE = 4

#: Every size macOS asks for, and the ones Windows uses. 16 is the reason the
#: ring is thick and the gaps are wide.
ICNS = [16, 32, 64, 128, 256, 512, 1024]
ICO = [16, 24, 32, 48, 64, 128, 256]


def draw(size: int = CANVAS) -> Image.Image:
    """The icon at `size`, drawn large and reduced so the edges are smooth."""
    big = _draw(size * OVERSAMPLE)
    return big.resize((size, size), Image.LANCZOS)


def _draw(size: int) -> Image.Image:
    image = Image.new("RGBA", (size, size), (0, 0, 0, 0))
    pen = ImageDraw.Draw(image)

    # A rounded square, which is what macOS and Windows both expect of an app;
    # a bare ring on the wallpaper reads as a loading spinner.
    margin = size * 0.06
    pen.rounded_rectangle(
        (margin, margin, size - margin, size - margin),
        radius=size * 0.20, fill=TILE,
    )

    outer = size * 0.33
    width = size * 0.115
    middle = size / 2
    box = (middle - outer, middle - outer, middle + outer, middle + outer)

    # Eight arcs with a gap between them, starting at twelve o'clock and going
    # clockwise, the way the Cassette screen numbers its slots.
    gap = 4.0
    span = 360 / len(PARTS)
    for position, colour in enumerate(PARTS):
        start = -90 + position * span + gap / 2
        pen.arc(box, start, start + span - gap, fill=colour, width=round(width))

    return image


def main() -> int:
    master = draw()

    iconset = HERE / "icon.iconset"
    if iconset.exists():
        shutil.rmtree(iconset)
    iconset.mkdir()
    for size in ICNS:
        frame = master.resize((size, size), Image.LANCZOS)
        frame.save(iconset / f"icon_{size}x{size}.png")
        # Retina: macOS wants the 2x of each size under the smaller name.
        if size * 2 in ICNS or size <= 512:
            master.resize((size * 2, size * 2), Image.LANCZOS).save(
                iconset / f"icon_{size}x{size}@2x.png")

    master.resize((256, 256), Image.LANCZOS).save(
        HERE / "icon.ico", sizes=[(n, n) for n in ICO])
    print(f"wrote {HERE / 'icon.ico'}")

    if sys.platform != "darwin":
        print("not macOS: leaving icon.icns alone, icon.ico is written")
        return 0

    done = subprocess.run(
        ["iconutil", "-c", "icns", str(iconset), "-o", str(HERE / "icon.icns")],
        capture_output=True, text=True,
    )
    if done.returncode != 0:
        print(done.stderr.strip(), file=sys.stderr)
        return 1
    shutil.rmtree(iconset)
    print(f"wrote {HERE / 'icon.icns'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
