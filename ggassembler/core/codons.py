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


def back_translate(
    protein: str,
    avoid: tuple[Enzyme, ...] = (),
    stop: str | None = None,
    table: dict[str, tuple[str, ...]] | None = None,
) -> tuple[str, list[int]]:
    """Write `protein` as DNA, never spelling a site for anything in `avoid`.

    Returns the sequence and the positions (in residues) where the preferred
    codon had to be passed over - which is worth reporting, because those are
    the only places the result differs from a plain codon-optimised gene.

    `stop` appends a stop codon when the protein does not carry one.
    """
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
        choice = next(
            (r for r in range(start, len(options))
             if not _makes_site(sequence, options[r], avoid)),
            None,
        )

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
