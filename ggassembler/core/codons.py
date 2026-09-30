"""Turn a protein back into DNA that can actually be used as a part.

UniProt, and any other protein database, hands you amino acids. A type 3 part
is DNA, so something has to bridge that - and the useful way to bridge it is
not simply to pick a codon per residue but to pick codons that *avoid the
kit's restriction sites while doing so*.

Domesticating a sequence afterwards means finding a site, changing a base, and
hoping the change is silent and does not create a new site two codons along.
Choosing the codon in the first place cannot create a site at all: at each
residue the most preferred codon that keeps the growing sequence clean is
taken, and only the last few bases can be affected, so the check is local and
exact.

The ordering is the classical set of *S. cerevisiae* preferred codons - the
ones over-represented in highly expressed genes. It is a preference order, not
a frequency table: this module does not compute a CAI and does not claim one.
Where the first choice would spell a site, the next synonymous codon is used,
so a run of rare codons only ever appears where the sequence forces it.
"""

from __future__ import annotations

from .enzymes import Enzyme

#: Synonymous codons per amino acid, most preferred first.
#:
#: These are the classical preferred-codon orderings for highly expressed genes
#: - the sets that have been in use since Bennetzen & Hall and Sharp & Cowe.
#: They are an *ordering*, not a frequency table: this module does not compute
#: a CAI and does not claim one. What the ordering has to be right about is
#: which codon to reach for first and what to fall back to, which is all the
#: site-avoidance below needs.
SCEREVISIAE: dict[str, tuple[str, ...]] = {
    "A": ("GCT", "GCC", "GCA", "GCG"),
    "R": ("AGA", "AGG", "CGT", "CGC", "CGA", "CGG"),
    "N": ("AAC", "AAT"),
    "D": ("GAC", "GAT"),
    "C": ("TGT", "TGC"),
    "Q": ("CAA", "CAG"),
    "E": ("GAA", "GAG"),
    "G": ("GGT", "GGC", "GGA", "GGG"),
    "H": ("CAC", "CAT"),
    "I": ("ATT", "ATC", "ATA"),
    "L": ("TTG", "CTA", "CTT", "CTG", "TTA", "CTC"),
    "K": ("AAG", "AAA"),
    "M": ("ATG",),
    "F": ("TTC", "TTT"),
    "P": ("CCA", "CCT", "CCC", "CCG"),
    "S": ("TCT", "TCC", "TCA", "AGT", "TCG", "AGC"),
    "T": ("ACT", "ACC", "ACA", "ACG"),
    "W": ("TGG",),
    "Y": ("TAC", "TAT"),
    "V": ("GTT", "GTC", "GTA", "GTG"),
    "*": ("TAA", "TGA", "TAG"),
}

#: The same for *E. coli*, for a part that will be expressed there rather than
#: in yeast - a marker, or a protein being made before it is moved across.
ECOLI: dict[str, tuple[str, ...]] = {
    "A": ("GCG", "GCT", "GCC", "GCA"),
    "R": ("CGT", "CGC", "CGG", "CGA", "AGA", "AGG"),
    "N": ("AAC", "AAT"),
    "D": ("GAT", "GAC"),
    "C": ("TGC", "TGT"),
    "Q": ("CAG", "CAA"),
    "E": ("GAA", "GAG"),
    "G": ("GGC", "GGT", "GGG", "GGA"),
    "H": ("CAT", "CAC"),
    "I": ("ATC", "ATT", "ATA"),
    "L": ("CTG", "TTA", "TTG", "CTC", "CTT", "CTA"),
    "K": ("AAA", "AAG"),
    "M": ("ATG",),
    "F": ("TTT", "TTC"),
    "P": ("CCG", "CCA", "CCT", "CCC"),
    "S": ("AGC", "TCT", "TCC", "AGT", "TCG", "TCA"),
    "T": ("ACC", "ACG", "ACT", "ACA"),
    "W": ("TGG",),
    "Y": ("TAT", "TAC"),
    "V": ("GTG", "GTT", "GTC", "GTA"),
    "*": ("TAA", "TGA", "TAG"),
}

#: The tables a caller can choose between, by the name the screen shows.
TABLES: dict[str, dict[str, tuple[str, ...]]] = {
    "scerevisiae": SCEREVISIAE,
    "ecoli": ECOLI,
}

TABLE_NAMES: dict[str, str] = {
    "scerevisiae": "S. cerevisiae",
    "ecoli": "E. coli",
}

#: The default, this being a yeast toolkit.
CODONS = SCEREVISIAE

#: Residues that carry no information about which codon to use.
AMBIGUOUS = {"X", "B", "Z", "J", "U", "O"}


def get_table(name: str | None) -> dict[str, tuple[str, ...]]:
    """One codon table by name, defaulting to yeast."""
    if not name:
        return SCEREVISIAE
    try:
        return TABLES[name.lower().replace(".", "").replace(" ", "")]
    except KeyError:
        raise ValueError(
            f"unknown codon table: {name}. Known: {', '.join(sorted(TABLES))}"
        ) from None


class BackTranslationError(ValueError):
    """The protein cannot be written as DNA that avoids every site."""


def clean_protein(sequence: str) -> str:
    """Strip whitespace, FASTA headers and numbering from a pasted protein."""
    out = []
    for line in sequence.splitlines():
        if line.startswith(">"):
            continue
        out.append("".join(c for c in line if c.isalpha() or c == "*"))
    return "".join(out).upper()


def _makes_site(sequence: str, codon: str, avoid: tuple[Enzyme, ...]) -> bool:
    """Whether adding `codon` spells a site for anything in `avoid`.

    Only the join needs checking: everything earlier was already clean, and a
    site cannot be longer than its own recognition sequence, so it is enough to
    look at the last few bases either side of the new codon - on both strands,
    since these enzymes cut whichever way round their site appears.
    """
    for enzyme in avoid:
        # The *recognition sequence*, not the cut. `find_sites` drops a site
        # whose cut falls outside the sequence it was given, and BbsI cuts
        # eight to twelve bases away - so in a window this short it reports
        # nothing at all, and twelve BbsI sites went straight through.
        window = len(enzyme.site) + len(codon)
        tail = (sequence + codon)[-window:]
        if enzyme.pattern("+").search(tail):
            return True
        if not enzyme.is_palindromic and enzyme.pattern("-").search(tail):
            return True
    return False


def _trailing_run(text: str) -> int:
    """How many of the same base the sequence currently ends on."""
    if not text:
        return 0
    last = text[-1]
    n = 1
    while n < len(text) and text[-1 - n] == last:
        n += 1
    return n


def back_translate(
    protein: str,
    avoid: tuple[Enzyme, ...] = (),
    stop: str | None = None,
    table: dict[str, tuple[str, ...]] | None = None,
    smooth: bool = True,
) -> tuple[str, list[int]]:
    """Write `protein` as DNA, never spelling a site for anything in `avoid`.

    Returns the sequence and the positions (in residues) where the preferred
    codon had to be passed over - which is worth reporting, because those are
    the only places the result differs from a plain codon-optimised gene. Both
    reasons for passing one over land in that list; what a reader needs from it
    is where the sequence is not the textbook answer, not why.

    `stop` appends a stop codon when the protein does not carry one.

    `smooth` lets a long single-base run break the tie between codons that are
    equally legal. Always taking the most preferred codon is a deterministic
    map from residue to bases, so a tract of one amino acid becomes a tract of
    one codon: a poly-lysine stretch under the *E. coli* table is AAA repeated,
    and thirty lysines are ninety adenines that no vendor will synthesise. It
    is a preference and never a requirement - see the retry below.
    """
    try:
        return _write(protein, avoid, stop, table, smooth)
    except BackTranslationError:
        if not smooth:
            raise
        # Preferring a synonym can walk the search into a corner that the plain
        # order would have got through. A gene with an awkward run in it beats
        # no gene at all, so the unsmoothed answer is still the right one.
        return _write(protein, avoid, stop, table, False)


def _write(
    protein: str,
    avoid: tuple[Enzyme, ...],
    stop: str | None,
    table: dict[str, tuple[str, ...]] | None,
    smooth: bool,
) -> tuple[str, list[int]]:
    codons = table or SCEREVISIAE
    residues = clean_protein(protein)
    if not residues:
        raise BackTranslationError("no protein sequence given")

    unknown = sorted({r for r in residues if r not in codons})
    if unknown:
        raise BackTranslationError(
            f"cannot write {', '.join(unknown)} as DNA"
            + (" - an ambiguous residue has no codon" if set(unknown) & AMBIGUOUS else "")
        )

    # A forward pass alone is not enough. The codons already written decide
    # which are legal next, so a residue can arrive with every one of its
    # codons spelling a site - and the fix is not there but one position back.
    # So this walks forward greedily and reverses a step when it is stuck,
    # which is a search rather than a scan and needs a budget to stay one.
    out: list[str] = []
    picked: list[int] = []
    budget = 40 * len(residues) + 1000
    index = 0

    while index < len(residues):
        if budget <= 0:
            raise BackTranslationError(
                "gave up finding codons that avoid every site; the protein may "
                "need a site removed from the avoid list"
            )
        budget -= 1

        options = codons[residues[index]]
        start = picked[index] if index < len(picked) else 0
        sequence = "".join(out)
        legal = [
            r for r in range(start, len(options))
            if not _makes_site(sequence, options[r], avoid)
        ]
        choice = legal[0] if legal else None

        if smooth and len(legal) > 1:
            # Chosen only from codons that were already legal, so this cannot
            # introduce a site; and only when a kinder one exists, so it cannot
            # remove the answer.
            kinder = next(
                (r for r in legal
                 if _trailing_run(sequence + options[r]) < HOMOPOLYMER_LIMIT),
                None,
            )
            if kinder is not None:
                choice = kinder

        if choice is None:
            if index == 0:
                raise BackTranslationError(
                    f"every codon for {residues[0]} spells a site for one of "
                    f"{', '.join(e.name for e in avoid)}"
                )
            # undo this position and make the one before it choose again
            del picked[index:]
            index -= 1
            picked[index] += 1
            out.pop()
            continue

        if index < len(picked):
            picked[index] = choice
        else:
            picked.append(choice)
        out.append(options[choice])
        index += 1

    compromised = [i for i, rank in enumerate(picked) if rank]

    if stop and residues[-1] != "*":
        out.append(stop)
    return "".join(out), compromised


def translate(dna: str) -> str:
    """The protein a coding sequence spells, for checking a round trip."""
    table = {
        codon: residue
        for source in (SCEREVISIAE, ECOLI)
        for residue, group in source.items()
        for codon in group
    }
    # the synonymous table above is not the whole genetic code; fill the rest
    for codon, residue in _REST.items():
        table.setdefault(codon, residue)
    dna = dna.upper().replace("U", "T")
    return "".join(
        table.get(dna[i:i + 3], "X") for i in range(0, len(dna) - len(dna) % 3, 3)
    )


#: The codons the preference table above happens not to list, so `translate`
#: can read a sequence this module did not write.
_REST: dict[str, str] = {
    "CTC": "L", "AGC": "S", "ATA": "I", "GTG": "V", "CCG": "P", "ACG": "T",
    "GCG": "A", "CGG": "R", "GGG": "G", "TAG": "*", "TGA": "*",
}


# --------------------------------------------------------------------------- #
# will a vendor actually make this?
# --------------------------------------------------------------------------- #
#
# A gene that is clean of restriction sites can still be refused, or quoted at
# a premium, or silently delivered wrong. Synthesis houses screen on sequence
# features rather than on biology: long single-base runs, exact internal
# repeats, and windows of extreme GC. None of that is visible in a translation
# check, and all of it is made *more likely* by the strategy above - always
# taking the most preferred codon is a deterministic map from residue to bases,
# so a run of one amino acid becomes a run of one codon, and a repeated motif
# in the protein becomes an exact repeat in the DNA.
#
# These are reported rather than enforced. Thresholds differ between vendors
# and change; what does not change is that you want to see the numbers before
# you pay for the gene.

#: A single-base run at or above this is worth flagging. Below it, nothing any
#: mainstream vendor objects to; at 8 the quotes start carrying notes.
HOMOPOLYMER_LIMIT = 8

#: An exact internal repeat at or above this length is what makes an assembly
#: mis-prime and a synthesis house ask questions.
REPEAT_LIMIT = 20

#: GC is measured across a sliding window, not over the whole gene: a sequence
#: can sit at a perfectly comfortable 45% overall and still carry a 50 bp
#: stretch at 15% that will not amplify.
GC_WINDOW = 50
GC_LOW, GC_HIGH = 25.0, 75.0


def longest_homopolymer(dna: str) -> tuple[str, int, int]:
    """The longest single-base run, as ``(base, length, 0-based start)``.

    An empty sequence has no run, reported as ``("", 0, -1)`` rather than
    raising: this is a report, and a report on nothing is not an error.
    """
    dna = dna.upper()
    if not dna:
        return "", 0, -1
    best = (dna[0], 1, 0)
    start = 0
    for i in range(1, len(dna) + 1):
        if i < len(dna) and dna[i] == dna[start]:
            continue
        if i - start > best[1]:
            best = (dna[start], i - start, start)
        start = i
    return best


def longest_repeat(dna: str, minimum: int = REPEAT_LIMIT) -> tuple[int, int, int]:
    """The longest exact repeat of at least `minimum`, as ``(length, first, second)``.

    Both positions are 0-based starts of the two copies; ``(0, -1, -1)`` means
    there is no repeat that long.

    Found by binary search over the length with a hash set per length, which is
    O(n log n) on the sequence rather than the O(n^2) a pairwise scan would
    cost. A 4 kb CDS is checked in a few milliseconds, which matters because
    this runs on every keystroke-triggered rebuild.
    """
    dna = dna.upper()
    if len(dna) < minimum * 2:
        return 0, -1, -1

    def found(length: int) -> tuple[int, int] | None:
        seen: dict[str, int] = {}
        for i in range(len(dna) - length + 1):
            chunk = dna[i : i + length]
            if chunk in seen:
                return seen[chunk], i
            seen[chunk] = i
        return None

    if found(minimum) is None:
        return 0, -1, -1

    low, high, best = minimum, len(dna) // 2 + 1, (minimum, *found(minimum))
    while low <= high:
        middle = (low + high) // 2
        hit = found(middle)
        if hit is None:
            high = middle - 1
        else:
            best = (middle, *hit)
            low = middle + 1
    return best


def gc_windows(dna: str, window: int = GC_WINDOW) -> tuple[float, float, int]:
    """GC across sliding windows, as ``(lowest %, highest %, windows outside)``.

    A sequence shorter than one window is measured whole, so the answer is
    still about the sequence rather than about there being no windows.
    """
    dna = dna.upper()
    if not dna:
        return 0.0, 0.0, 0
    size = min(window, len(dna))

    running = sum(1 for base in dna[:size] if base in "GC")
    lowest = highest = 100.0 * running / size
    outside = 1 if not (GC_LOW <= lowest <= GC_HIGH) else 0

    for i in range(size, len(dna)):
        running += (dna[i] in "GC") - (dna[i - size] in "GC")
        fraction = 100.0 * running / size
        lowest = min(lowest, fraction)
        highest = max(highest, fraction)
        if not (GC_LOW <= fraction <= GC_HIGH):
            outside += 1

    return round(lowest, 1), round(highest, 1), outside


def feasibility(dna: str) -> dict[str, object]:
    """Every synthesis check in one shape, for a panel to render row by row."""
    base, run, at = longest_homopolymer(dna)
    length, first, second = longest_repeat(dna)
    low, high, outside = gc_windows(dna)
    return {
        "homopolymer": {
            "base": base, "length": run, "at": at,
            "limit": HOMOPOLYMER_LIMIT, "ok": run < HOMOPOLYMER_LIMIT,
        },
        "repeat": {
            "length": length, "first": first, "second": second,
            "limit": REPEAT_LIMIT, "ok": length < REPEAT_LIMIT,
        },
        "gc": {
            "min": low, "max": high, "outside": outside, "window": GC_WINDOW,
            "low": GC_LOW, "high": GC_HIGH, "ok": outside == 0,
        },
    }
