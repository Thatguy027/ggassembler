# Screenshots for GUIDE.md

Six images, referenced by exact filename from [`../../GUIDE.md`](../../GUIDE.md).
Until they exist the guide still reads correctly — each one renders as its alt
text, which says what it would have shown.

Take them with the app running (`ggasm serve ytk`) at a window width of roughly
1400px, in whichever theme you prefer. Circle the named control in red.

| File | Screen | Capture | Circle |
|---|---|---|---|
| `04-library-sync.png` | `/library` | the sharing bar across the top | the **Share my plasmids with the lab** switch and **Sync now** |
| `05-rescan.png` | `/library` | the top bar, right-hand side | **Re-scan** |
| `06-part.png` | `/part` | a filled-in design, not an empty form | the sequence box, the **part type** selector, **Save to library** |
| `07-cassette.png` | `/` | the worked assembly it opens on, ring map visible | one part panel and **Save** |
| `08-multigene.png` | `/multigene` | two or three TUs added | **Add TU** and the **Multigene vector** selector |
| `09-sequencing.png` | `/sequencing` | after an alignment has run, a verdict visible | the reference picker, **Add clone**, **Align** |

Two things worth getting right, because they are what makes a screenshot useful
rather than decorative:

**Fill the screen in first.** An empty form photographs badly and teaches
nothing. The Cassette screen opens on a real assembly from your own library
already, so it is ready as-is.

**Nothing unpublished in frame.** These go in a public repository. The starter
library is a safe thing to screenshot; your lab's collection is not, so run
`ggasm serve ytk` against the starter library rather than the full one.
