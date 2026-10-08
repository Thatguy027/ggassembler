# Screenshots for GUIDE.md

| File | Screen | Shows | Marked |
|---|---|---|---|
| `part_demo.png` | `/part` | the demo CDS designed as a Type 3 part | sequence box, part type, Save to library |
| `cassette.png` | `/` | a complete eight-part assembly | one part panel, Save |
| `multi.png` | `/multigene` | four transcription units | vector selector, Add TU |
| `align_demo.png` | `/sequencing` | three clones aligned - clean, deletion, substitution | reference sources, drop zone, Align |
| `protein.png` | `/proteins` | not referenced by the guide yet | - |

Still missing, for sections 4 and 5, which read fine without them:

| | Screen | Would show |
|---|---|---|
| sharing | `/library` | the **Share my plasmids with the lab** switch and **Sync now** |
| re-scan | `/library` | **Re-scan**, top right |

## If you retake any of these

Run the app against the starter library, not your own:

```sh
ggasm serve ~/ytk
```

Check the header afterwards. If it names any library but `ytk`, something else
was already holding port 8737 and the browser was showing *that* server - which
is how the first set came to be shot against a private library despite the
right command being typed. `ggasm serve` now refuses to start in that case
rather than printing an address it is not serving.

Demo material for filling the screens is in `~/ggasm-demo/` - a synthetic CDS
for the Part screen, and three clones of pYTK047 for Sequencing.

The paths in the current set were painted out after the fact. Retaking against
`~/ytk` avoids needing that.
