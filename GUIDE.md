# Using GG Assembler

Nine things, in the order you will need them. Each one is the shortest version
that works; open the fold under it when you want to know why, or when the short
version is not enough.

1. [Install the app](#1-install-the-app)
2. [Get a plasmid library](#2-get-a-plasmid-library)
3. [Start the app](#3-start-the-app)
4. [Share your lab's library](#4-share-your-labs-library)
5. [Add new plasmids](#5-add-new-plasmids)
6. [Build a part plasmid](#6-build-a-part-plasmid)
7. [Build an eight-part cassette](#7-build-an-eight-part-cassette)
8. [Build a multi-TU plasmid](#8-build-a-multi-tu-plasmid)
9. [Check a sequencing result](#9-check-a-sequencing-result)

---

## 1. Install the app

```sh
uv tool install ggassembler
```

<details>
<summary>No uv, or it is not on your PATH</summary>

`uv` is the fastest route because it installs Python for you. Without it:

```sh
pipx install ggassembler          # or
pip install ggassembler           # inside a virtualenv
```

Needs Python 3.11 or newer. If `ggasm` is not found after installing, the tool
directory is not on your PATH — `uv tool update-shell` fixes that for uv, then
open a new terminal.

To upgrade later: `uv tool upgrade ggassembler`.
</details>

---

## 2. Get a plasmid library

The app does not ship with plasmids. It reads a folder of GenBank files, and
you need one before anything else is useful.

```sh
ggasm init ytk --starter
```

That clones the published YTK and MYT kits — 192 records, enough to build a
complete cassette.

<details>
<summary>Joining a lab that has its own library</summary>

```sh
ggasm init plasmids --from git@github.com:<org>/<library>.git
```

This is the usual case once you are in a lab. It clones the repository, tells
git to ignore the index and merge the build log, and records a baseline if
there is not one.

If it stops with a message about access, you have no key on your GitHub account
or no access to that repository. `gh auth login` sets the first up; the second
is a conversation with whoever runs the library.
</details>

<details>
<summary>Already have a folder of GenBank files</summary>

Nothing to clone. Point the app at it directly in step 3. Files may be `.gb` or
`.gbk`, in one folder or nested — names and annotations are ignored, because
what each plasmid *is* gets decided by simulating the digest.

To make an existing folder shareable later: `ggasm init /path/to/folder`.
</details>

---

## 3. Start the app

```sh
ggasm serve ytk
```

It opens in your browser. Leave the terminal running — closing it stops the app.
Press `Ctrl-C` there when you are finished.

<details>
<summary>Several libraries at once, a different port, and other options</summary>

A public base library and your lab's own scan as a single library:

```sh
ggasm serve ~/ggasm/ytk ~/ggasm/mylab --cache-dir ~/.ggasm
```

`--cache-dir` matters when the folders are far apart: without it the index
lands in whatever parent directory they happen to share.

| | |
|---|---|
| `--port 9000` | if 8737 is taken |
| `--no-browser` | do not open a tab |
| `--reload` | restart on `.py` changes — only useful if you are editing the app |

It serves on `http://127.0.0.1:8737`, which means *this machine*. Nobody else
can reach it, and nothing is uploaded anywhere.

The first start indexes every file and prints how many it found. Later starts
re-read only what changed.
</details>

---

## 4. Share your lab's library

Open the **Library** tab. Sharing is on by default.

![The Library screen, with the sharing bar at the top](docs/images/04-library-sync.png)

Saving a plasmid on any screen publishes it to your lab. **Sync now** pulls in
what colleagues have added.

<details>
<summary>What the toggle does, and what it does not</summary>

**Share my plasmids with the lab** controls publishing only. Switch it off and
your plasmids stay on this machine — but you still receive everyone else's,
because working against a library you know is stale is how an assembly gets
repeated.

The setting is yours alone. It lives in `.ggasm/` beside the library, which is
never shared, so switching it off does not switch it off for the lab.

Before anything is published it is checked: a plasmid that will not digest, or
one nothing can classify, is refused rather than pushed. Problems the library
already had never block you — only ones your change introduces.

Sharing is git. Your library folder is a repository, and one private repository
per lab is what keeps one lab's plasmids out of another's.
</details>

---

## 5. Add new plasmids

Copy the `.gb` or `.gbk` files into your library folder, then press **Re-scan**
on the Library screen.

![The Re-scan button, top right of the Library screen](docs/images/05-rescan.png)

They are available immediately, and still there next time you start the app.

<details>
<summary>They are there but not usable, or classified wrongly</summary>

Everything is decided by the digest, so a file whose name says `promoter` and
whose sequence says otherwise is listed as what it is. Rows needing a human are
flagged amber: a conflict between digest and annotations, an internal site, or
nothing the digest could call.

Click the type pill on any row to assign a type by hand. Manual assignments beat
detection, are stored one file per plasmid under `ggasm/overrides/`, and are
shared with your lab alongside the plasmid itself.

Re-scan reads only files whose timestamp or size changed, so it stays quick on a
library of several hundred.

Plasmids saved from inside the app land in the library folder on their own — no
copying needed.
</details>

---

## 6. Build a part plasmid

**Part plasmid** tab. Paste a sequence, pick what it should become, save.

![The Part screen with the sequence box, part type and Save circled](docs/images/06-part.png)

1. Paste your sequence into the sequence box
2. Choose a **part type** — 3 for a CDS, 2 for a promoter, 4 for a terminator
3. Name it
4. **Save to library**

It designs the primers that add the right BsaI ends, and shows them on the
template with the predicted construct below.

<details>
<summary>Sequence conventions, and when to change them</summary>

The checkboxes under **Sequence conventions** decide what gets added or removed
at the ends. Each has an ⓘ explaining what it does and when it matters. The ones
that catch people out:

- **Type 3: drop a leading ATG** — the TATG overhang already supplies it, so a
  CDS pasted with its own ATG would otherwise get two.
- **Strip a trailing stop codon** — needed if anything is fused downstream.
- **Type 4: lead with TAA + CTCGAG** — supplies the stop a Type 3 part had
  removed.
- **Gly-Ser linker** — off by default. Only for a fusion where the two halves
  need spacing.

**Removing internal sites** is offered, never assumed. On a promoter or
terminator it is often unclear which change is safe, and that is your call, not
the app's.

Pick the **entry vector** to match the ends you want. **Export .gbk** gives you
the file without adding it to the library.
</details>

---

## 7. Build an eight-part cassette

**Cassette** tab. It opens on a worked assembly from your own library, so you
can see a complete design before changing anything.

![The Cassette screen, with a part panel and the Save button circled](docs/images/07-cassette.png)

1. Click a panel and pick a part — each offers only parts whose overhangs fit
   that position
2. Name the construct
3. **Save**

The ring map redraws as you go. If the assembly cannot close, it says which
junction is the problem rather than just failing.

<details>
<summary>Splitting positions, composite parts, and the other two directions</summary>

Panels carry a segmented toggle where the standard allows a split — `3 | 3a+3b`
on Type 3, the same on 4 and 8. Switch it and the single panel becomes two.
A `2·3·4` or `6·7·8` composite collapses three positions into one.

**Start from a cassette…** takes a construct you already have and fills the
panels in, ready to swap one part. The assembly consumed every BsaI site, so
this is not a digest — each library part's fragment is looked for in the
sequence and the ones that abut are chained together.

**Sweep a position…** builds the same design once per candidate — every
promoter in your library against a fixed CDS, say — and writes them all out.

**Protocol** gives the pipetting steps and volumes. **Direction** and
**Linear view** change how the map is drawn, not what is built.
</details>

---

## 8. Build a multi-TU plasmid

**Multigene** tab. Chain cassettes into one construct.

![The Multigene screen with Add TU and the vector selector circled](docs/images/08-multigene.png)

1. **Add TU** for each transcription unit, and pick a cassette for each
2. Choose a **multigene vector**
3. Name it, then **Save to library**

Order comes from the connector pair each cassette carries, not from the order
you added them.

<details>
<summary>Designing from the other end, when you do not have the cassettes yet</summary>

**Design…** works backwards: say what you want expressed and it works out which
cassettes you would need, which of them you already have, and what is missing.
**Work it out** runs it; **Download all** writes the lot.

This is the more useful direction when starting a new construct — you rarely
begin with the cassettes already in hand.

**Protocol** gives the assembly steps for the BsmBI reaction.
</details>

---

## 9. Check a sequencing result

**Sequencing** tab. Compare clones a vendor sent back against what they should
have been.

![The Sequencing screen with the reference picker, Add clone and Align circled](docs/images/09-sequencing.png)

1. Pick the reference — **A plasmid in the library**, **A construct you
   designed**, or **A file or paste**
2. **Add clone** for each sequenced file
3. **Align**

Each clone gets a verdict — matches, or differs here — with a scrollable
alignment below it.

<details>
<summary>Whole-plasmid reads, and reading the alignment</summary>

Built for whole-plasmid nanopore results, where the clone starts at an arbitrary
position and may be the reverse complement. It finds where the two circles
correspond and opens them at the same place, so a rotated read does not look
like a rearrangement.

**Next difference →** and **← Previous difference** jump between mismatches.
**Open all** expands every clone's detail at once.

Differences are listed below the alignment, folded, so a clone with forty
mismatches does not bury one with a single mismatch.

A long stretch the app will not align base by base is shown grey and labelled,
rather than painted as thousands of mismatches.

A construct you saved from the Cassette or Multigene screen is available as
**A construct you designed** — the predicted sequence is kept precisely so that
sequencing it later has something to check against.
</details>

---

## When something is wrong

| | |
|---|---|
| `ggasm: command not found` | first check it is actually installed: `uv tool list` should name `ggassembler`. If it is there, zsh has cached the miss — run `rehash` or open a new terminal. Only if it is still missing is PATH the problem: `uv tool update-shell` |
| `not a folder: ytk` | the library has not been created yet — `ggasm init ytk --starter` |
| The page will not load | the terminal running `ggasm serve` has stopped; start it again |
| A plasmid is missing | **Re-scan** on the Library screen; if still missing, it is not a GenBank file the app could read |
| A part is not offered | its overhangs do not fit that position — the Library screen says what it was detected as |
| Sharing says it cannot reach GitHub | `ssh -T git@github.com` to check access, `gh auth login` to set it up |
| Changes to the app's code do nothing | restart `ggasm serve`, or run it with `--reload` |
