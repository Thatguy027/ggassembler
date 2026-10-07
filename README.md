# GG Assembler

A local desktop-in-browser app for designing Golden Gate assemblies in the Yeast
Toolkit (YTK) standard of Lee, DeLoache, Cervantes & Dueber, *ACS Synth. Biol.*
2015.

Point it at a folder of GenBank files. It works out what each plasmid **is** by
simulating the restriction digest, then drives seven screens: part construction,
cassette assembly, multigene assembly, a library index, a protein search that
writes a gene you can order, a plate that builds ninety-six cassettes at once,
and a check of what came back from sequencing.

## Get started

Three commands, from nothing to a running app with 192 published plasmids in it:

```sh
uv tool install ggassembler
ggasm init ytk --starter
ggasm serve ytk
```

The last command opens the app in your browser. That is the whole setup —
everything below is detail you can read when you need it.

| | |
|---|---|
| **Already have plasmids?** | `ggasm serve /path/to/your/genbank/folder` |
| **Joining a lab?** | `ggasm init plasmids --from git@github.com:<org>/<library>.git` |
| **Working on the app itself?** | `uv sync && uv run ggasm serve <folder> --reload` |
| **Need Python?** | 3.11 or newer. [uv](https://docs.astral.sh/uv/) installs it for you. |
| **No uv?** | `pipx install ggassembler`, or `pip install ggassembler` into a virtualenv. |

Nothing is uploaded anywhere. The app runs on your machine, reads your files,
and only talks to the network if you ask it to search UniProt or share a
library with your lab.

---

<details>
<summary>The longer version of the above</summary>

```sh
uv sync
uv run ggasm serve ../goldengateR/inst/extdata/plasmids
```

With nothing on your machine yet, start from the two published kits:

```sh
ggasm init ytk --starter
ggasm serve ytk
```

That clones [ytk-starter-library](https://github.com/Thatguy027/ytk-starter-library)
— pYTK001–096 and pMYT001–096, 121 usable parts covering every position of the
part circle, no detection conflicts — which is enough to design and simulate a
complete cassette without adding anything of your own.

Joining a lab that already has its own library:

```sh
ggasm init plasmids --from git@github.com:<org>/<library>.git
ggasm serve plasmids
```

A public base library and a private one scan as a single library, so you can
have both:

```sh
ggasm serve ~/ggasm/ytk ~/ggasm/mylab --cache-dir ~/.ggasm
```

`init` clones it, tells the repository to ignore the index and merge the build
log, and agrees a baseline if there is not one. It is the three commands nobody
should have to be told twice, and the one place where GitHub access is likely
to stop someone - so when it does, it says to run `gh auth login` rather than
passing on git's sentence about public keys.

Serves on `http://127.0.0.1:8737`. Several folders scan as one library:

```sh
uv run ggasm serve ..          # every plasmid in the project
```

Files under `web/` are served with `no-cache`, so an edited stylesheet or
script needs nothing but a browser refresh. Python is not like that - the
routers register when the app is built, so a running server keeps serving the
routes it started with. `--reload` restarts it when a `.py` file changes:

```sh
uv run ggasm serve ../goldengateR/inst/extdata/plasmids --reload
```

</details>

## Sharing a library with a lab

A library folder holds two unrelated things, and they belong on opposite sides
of a `.gitignore`:

| | What it is | Shared |
|---|---|---|
| `.ggasm/` | the index, keyed on mtimes, and this person's settings | no - it differs on every machine and would conflict on every pull |
| `ggasm/`  | overrides, saved designs, the build log, the baseline | **yes** - nothing can recompute a concentration measured at a bench |

Both sit beside the library unless `--cache-dir` / `--data-dir` say otherwise,
which matters when several roots are scanned as one library: their default
location is whatever ancestor the roots share.

Sharing itself is git. The library folder is already a repository - for this
lab, [ggassembler-library](https://github.com/Thatguy027/ggassembler-library) -
so one private repository per organisation gives privacy between organisations with no
server to run and nobody's unpublished constructs on someone else's
infrastructure. The Library screen carries the control - **Share my plasmids
with the lab**, on by default - and switching it off stops publishing without
stopping you receiving, because working against a library you know to be stale
is how an assembly gets repeated.

Before anything is published it is checked against `ggasm/baseline.json`, the
same gate `ggasm check` applies: a change is judged by the problems it *adds*,
never by the absolute state. A library of several hundred files always carries
some conflicts, and a gate that demands they all be fixed first is a gate that
gets switched off in a week.

Overrides and designs are one file each so that two people curating two
plasmids touch two files, and the build log is one JSON object per line with
`merge=union` set, so concurrent appends merge instead of conflicting.

## The idea it rests on

**What a plasmid is gets decided by the enzymes, never by its name or its
annotations.** In a real library both go stale: a file called `pYTK009.gb` may
have been edited, and a label reading `type 4 terminator` may describe a part
that went *into* the construct rather than the construct itself.

So the index simulates the digest. A plasmid whose part enzyme releases one
site-free fragment with overhangs on the part circle is a type 1–8 part. One
whose multigene enzyme releases a fragment with connector ends is a finished
cassette — or, with a dropout between those ends, the destination vector for
one. Where the digest and an annotation disagree, the digest wins and the
disagreement is **recorded as a conflict, never silently dropped**.

This is load-bearing. On one real library of 622 plasmids it recovers 20
cassettes hidden by stray internal sites, finds 6 files labelled `linear` that
are plainly circular, and rejects 23 assembled constructs that carried a part's
label and would otherwise have been offered as that part.

## The screens

**Part plasmid** — paste a sequence, pick a position, get the PCR primers that
add the right BsaI ends. Handles the lab's ccdB/L131/BbsI entry system. Removing
internal sites is offered, not assumed: on a promoter or a terminator it is
often unclear which change is safe, and that is the user's call. Shows the
primers on the template, the predicted construct, and the digested destination
so the ligation can be read.

**Cassette** — the eight positions of a transcription unit, each offering only
the parts whose overhangs fit. Opens on a worked assembly from your own library
rather than eight empty dropdowns. A ring map with callout labels, a linear
view, a junction strip, and a part picker that filters on name, aliases,
component and the annotations inside the part. Two directions beyond that:

- **Start from a cassette** — pick a construct you already have and get its
  panels filled in, ready to swap one part. The assembly consumed every BsaI
  site, so this is not a digest: each library part's fragment is looked for in
  the sequence and the ones that abut are chained into a tiling.
- **Sweep a position** — build the same design once per candidate in one
  position. A promoter titration against a fixed CDS is 48 constructs; you get
  a table and a zip of every `.gb` plus a picklist.

**Multigene** — cassettes chained into one construct by their connectors. The
bar shows each unit opened into the promoter, CDS and terminator that built it.
**Design** works the other way round: say what you want expressed, and it
returns the eight-part plasmids you have to build first, with connectors
assigned so they chain — `ConLS → ConR1 / ConL1 → ConR2 … → ConRE`, the kit's
own series — plus a build order and a `.gb` for every plasmid in it.

**Library** — what the digest found, with the rows that need a human flagged:
a conflict, an internal site, or no call at all. Unrecognised plasmids are
grouped by *why*. Identical files are merged, and you choose which name the
sequence goes by. Manual type assignments and measured concentrations live in
`.ggasm/` and always beat detection.

**Proteins** — the one screen that uses the network. Search UniProt, and the
protein comes back as DNA. Codons are chosen one at a time, never using one
that would spell a site — which is not the same as domesticating afterwards,
where a synonymous change made to remove one site can spell another two codons
downstream. Two choices change the sequence: codon usage (*S. cerevisiae* or
*E. coli*, the classical preferred-codon orderings for highly expressed genes —
an ordering, not a frequency table, so no CAI is claimed) and which enzymes to
keep out, defaulting to the scheme's own. Unticking one is a real choice: if
the part will never meet that enzyme, avoiding it is not worth a run of rare
codons. The result is a gene to order, not the organism's own sequence.

**Plate** — the same cassette built ninety-six ways. Pick a factor for the rows
and one for the columns, and every well is a real assembly checked before it is
costed, not a label. Carries a source-plate map of what is physically in which
well, the reaction maths, a plate map .csv, and an Opentrons Flex protocol that
was validated by running `opentrons_simulate` against it rather than by reading
the docs.

**Sequencing** — the clones a vendor sent back, checked against the plasmid
they were meant to be. A whole-plasmid read starts at an arbitrary base and is
as often as not on the opposite strand; neither is a difference, so each clone
is *placed* first — shared k-mers vote on a diagonal, and both circles are then
cut inside a shared exact match so the join cannot invent an indel that no
clone carries. What is left is read off base by base and named by the feature
it lands in. A clone that is not this plasmid is said to be unplaced rather
than aligned into a wall of false mismatches. Constructs exported from the
Cassette and Multigene screens are kept as references, so a design made today
can be checked against reads that arrive next week.

## What it will not do

- **Guess.** A part that cannot be read is reported as unreadable; a stretch of
  a construct that no library part explains is reported as a gap, with the
  overhangs that name the position it occupies.
- **Assume a concentration.** A reaction costed at the 50 ng/µL default says so
  on the row rather than presenting the volume as a measurement.
- **Silently domesticate.** Removing an internal site changes the sequence you
  will order, so it is always a choice.
- **Align something that is not the plasmid.** A sequencing result that shares
  no stretch with the reference is reported as unplaced. Forcing it through a
  global alignment would produce thousands of mismatches and bury the one
  clone that is genuinely wrong by a single base.

## Command line

```sh
uv run ggasm scan ..                       # index and print what it found
uv run ggasm serve ../goldengateR          # the app
uv run ggasm build <folder> -p pYTK003 -p pMYT001 ... -o out.gb

# fail only on problems this commit introduced, for a pre-commit hook
uv run ggasm check <folder> --update       # record the baseline once, commit it
uv run ggasm check <folder>                # exit 1 on anything new
```

`check` is deliberately quiet about problems that were already there. A library
of several hundred files always has some, and a hook that fails on all of them
is a hook people disable in a week.

## Layout

```
ggassembler/
  core/     enzymes, the part circle, the library index, assembly, search,
            decomposition, protocol, codons, UniProt, alignment
  levels/   one module per assembly level; levels never import each other
  api/      FastAPI; routes_levelN touches only levels/levelN_* and core
  web/      static front end, one self-contained directory per screen
tests/
docs/design/
```

`tests/test_layering.py` walks the AST of every module and fails on any import
that crosses those boundaries, so the three levels stay independently
modifiable: changing Level 2 can never require editing a Level 1 or Level 3
file. Work that both levels need — reading a construct back into parts, for
instance — goes to `core/` rather than being imported sideways or duplicated.

The front end is vanilla ES modules and plain CSS with no build step. Only
`tokens.css` is shared between screens. Two pieces of arithmetic are kept in
files with no DOM in them and run directly under node, because both are easy to
get subtly wrong and impossible to eyeball: `web/level2/layout.js`, the rule
that separates two crowded map labels without letting their leader lines cross,
and `web/sequencing/columns.js`, the mapping between a screen column and a
reference base — which drifts by however many bases the clones before it
inserted, and when it is wrong nothing looks broken; every label just points a
few bases off.

## Test

```sh
uv run pytest        # 993 tests, on 3.11, 3.12 and 3.13 in CI
```

Several exist because of a specific bug and say so. A few examples: an
annealing region must exist on its own template (a convention base placed in
the annealing region instead of the tail matched 1,581 bp upstream and would
have amplified a wrong product that looked fine); a button in the markup must
have a real handler; no stylesheet may be sized in pixels where the layout
scales; `ConR1` closing one unit must be the same overhang as `ConL1` opening
the next; and a sequencing result rotated to start at every one of six
different bases must still come back with no differences at all.

## Schemes

`core/parttypes.Scheme` binds a part enzyme, a multigene enzyme, a linearizer,
alternate cloning enzymes and an overhang table into one object. `YTK` is the
implemented kit — BsaI for the eight-part cassette, BsmBI for multigene and
entry vectors, BbsI for the ccdB entry system, NotI to linearize. Another kit,
or another enzyme chosen to dodge an internal site, is a new `Scheme` rather
than a change to the code that uses one: every function that cares takes
`scheme=` and defaults to `YTK`.

## Design

`docs/design/` holds the canvas this is built to. Every colour, type size and
radius comes from `web/tokens.css`, including the eight part-arc hues, which
`core/parttypes.py` also serves to the API so the map, the callouts and the
panel badges cannot drift apart. The three typefaces are vendored in
`web/fonts/`, so every screen but Proteins works with no network at all.

## Licence

MIT. Copyright (c) 2026 Stefan Zdraljevic. See [LICENSE](LICENSE).

The bundled fonts are not covered by it. IBM Plex and Space Grotesk are
redistributed under the SIL Open Font License 1.1, whose text travels with them
in [`ggassembler/web/fonts/LICENSE-fonts.txt`](ggassembler/web/fonts/LICENSE-fonts.txt).

No plasmid sequences are in this repository, now or in its history. The test
fixtures are built from the part-type definitions rather than checked in, so
the tests describe the standard rather than any particular lab's collection.
A library is something you point the app at.
