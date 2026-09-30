# GG Assembler

A local desktop-in-browser app for designing Golden Gate assemblies in the Yeast
Toolkit (YTK) standard of Lee, DeLoache, Cervantes & Dueber, *ACS Synth. Biol.*
2015. Point it at a folder of GenBank files; it infers what each plasmid is by
simulating the restriction digest, then drives three independent assembly
designers (part, cassette, multigene).

Built from `../goldengateR/gui/gg-assembler-claude-code-prompt.md`.

## Install

```sh
uv sync                 # dev
uv run ggasm --help
```

`pip install -e .` works too.

## Test

```sh
uv run pytest
```

## Layout

```
ggassembler/
  core/     pure logic - no web framework, no level-specific knowledge
  levels/   one module per assembly level; levels never import each other
  api/      FastAPI routes; routes_levelN touches only levels/levelN_* and core
  web/      static front end; one self-contained directory per level
tests/
docs/design/   design canvas exports + tokens.css
```

`tests/test_layering.py` walks the AST of every module and fails on any import
that crosses those boundaries, so the three levels stay independently
modifiable.

## Run it

```sh
uv run ggasm serve ../goldengateR/inst/extdata/plasmids
```

Indexes the folder, starts the app on <http://127.0.0.1:8737> and opens a tab.
Several folders scan as one library:

```sh
uv run ggasm serve ..          # every plasmid in the project
```

All four screens exist: **Part plasmid** (Level 1), **Cassette** (Level 2),
**Multigene** (Level 3) and **Library**.

The Library screen lists what the digest found, flags the rows that need a
human - a conflict, an internal site, or no call at all - and lets you assign a
type by hand. Manual assignments live in `.ggasm/overrides.json` and always beat
detection.

## Command line

```sh
# index every plasmid in the project (several folders scan as one library)
uv run ggasm scan ..
uv run ggasm scan ../goldengateR/inst/extdata/plasmids --no-recursive

# assemble a cassette from named library parts
uv run ggasm build ../goldengateR/inst/extdata/plasmids \
  -p pYTK003 -p pMYT001 -p pMYT016 -p pMYT025 \
  -p pYTK068 -p pYTK074 -p pYTK081 -p pYTK084 \
  --name pMyCassette -o pMyCassette.gb
```

The index is cached in `<common parent>/.ggasm/index.json`; manual type
assignments go in `.ggasm/overrides.json` and always beat detection.

## Design

`docs/design/` holds the canvas this is built to: one HTML page per screen,
extracted from `goldengateR/gui/GG Assembler.html`. Every colour, type size and
radius in the app comes from `ggassembler/web/tokens.css`, which is transcribed
from those pages - including the eight part-arc hues, which `core/parttypes.py`
also serves to the API so the map, the legend and the panel badges cannot drift
apart. The three typefaces (IBM Plex Sans, IBM Plex Mono, Space Grotesk) are
vendored in `ggassembler/web/fonts/`, so the app needs no network.

## Assembly schemes

`core/parttypes.Scheme` binds a part enzyme, a multigene enzyme, a linearizer
and an overhang table into one object. `YTK` is the implemented kit - BsaI for
the eight-part cassette, BsmBI for multigene and entry vectors, NotI to
linearize. Another kit, or another enzyme pair chosen to dodge an internal
site, is a new `Scheme`, not a change to the code that uses one: every function
that cares takes `scheme=` and defaults to `YTK`.
