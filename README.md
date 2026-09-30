# GG Assembler

A local desktop-in-browser app for designing Golden Gate assemblies in the Yeast
Toolkit (YTK) standard of Lee, DeLoache, Cervantes & Dueber, *ACS Synth. Biol.*
2015.

Point it at a folder of GenBank files. It works out what each plasmid **is** by
simulating the restriction digest, then drives five screens: part construction,
cassette assembly, multigene assembly, a library index, and a protein search
that writes a gene you can order.

```sh
uv sync
uv run ggasm serve ../goldengateR/inst/extdata/plasmids
```

Opens on <http://127.0.0.1:8737>. Several folders scan as one library:

```sh
uv run ggasm serve ..          # every plasmid in the project
```

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

## What it will not do

- **Guess.** A part that cannot be read is reported as unreadable; a stretch of
  a construct that no library part explains is reported as a gap, with the
  overhangs that name the position it occupies.
- **Assume a concentration.** A reaction costed at the 50 ng/µL default says so
  on the row rather than presenting the volume as a measurement.
- **Silently domesticate.** Removing an internal site changes the sequence you
  will order, so it is always a choice.

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
            decomposition, protocol, codons, UniProt
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
`tokens.css` is shared between screens. The map's geometry and its label
collision pass live in `web/level2/layout.js` with no DOM in them, so they are
tested directly under node — the rule that separates two crowded labels without
letting their leader lines cross is easy to get subtly wrong and impossible to
eyeball.

## Test

```sh
uv run pytest        # 599 tests
```

Several exist because of a specific bug and say so. A few examples: an
annealing region must exist on its own template (a convention base placed in
the annealing region instead of the tail matched 1,581 bp upstream and would
have amplified a wrong product that looked fine); a button in the markup must
have a real handler; no stylesheet may be sized in pixels where the layout
scales; and `ConR1` closing one unit must be the same overhang as `ConL1`
opening the next.

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
