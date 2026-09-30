"""GenBank read/write via Biopython ``SeqRecord``; circular topology handling.

Everything downstream assumes plasmids are circular, so the one job that
matters here is making a circular sequence searchable as if it were linear:
`doubled()` concatenates the sequence with itself, and `wrap()` maps a position
in that doubled coordinate space back onto the real molecule.
"""

from __future__ import annotations

import warnings
from pathlib import Path

from Bio import BiopythonParserWarning, BiopythonWarning, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

GENBANK_SUFFIXES = (".gb", ".gbk", ".genbank")
FASTA_SUFFIXES = (".fa", ".fasta", ".fna", ".seq")


def read_records(path: str | Path) -> list[SeqRecord]:
    """Read every record in a GenBank or FASTA file."""
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix in GENBANK_SUFFIXES:
        fmt = "genbank"
    elif suffix in FASTA_SUFFIXES:
        fmt = "fasta"
    else:
        raise ValueError(f"unsupported sequence file: {path}")
    with warnings.catch_warnings():
        # plenty of real-world GenBank files have a non-standard LOCUS line;
        # Biopython still parses them correctly and the noise buries real errors
        warnings.simplefilter("ignore", BiopythonParserWarning)
        return list(SeqIO.parse(str(path), fmt))


def read_record(path: str | Path) -> SeqRecord:
    """Read a file expected to hold exactly one record."""
    records = read_records(path)
    if len(records) != 1:
        raise ValueError(f"{path}: expected 1 record, found {len(records)}")
    return records[0]


def write_genbank(record: SeqRecord | list[SeqRecord], path: str | Path) -> Path:
    """Write one or more records as GenBank, creating parent directories."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    records = record if isinstance(record, list) else [record]
    for rec in records:
        rec.annotations.setdefault("molecule_type", "DNA")
    with warnings.catch_warnings():
        # long ApE qualifier keys are common in real files and harmless here
        warnings.simplefilter("ignore", BiopythonWarning)
        SeqIO.write(records, str(path), "genbank")
    return path


def is_circular(record: SeqRecord) -> bool:
    """True unless the record explicitly says it is linear.

    GenBank files written by many tools omit the topology field; for a plasmid
    library, circular is the safer default and the caller can override.
    """
    return str(record.annotations.get("topology", "circular")).lower() != "linear"


def topology(record: SeqRecord) -> str:
    return "circular" if is_circular(record) else "linear"


def sequence(record: SeqRecord) -> str:
    return str(record.seq).upper()


def doubled(seq: str) -> str:
    """The sequence concatenated with itself, for searching across the origin."""
    return seq + seq


def wrap(position: int, length: int) -> int:
    """Map a position in doubled coordinates back onto a circle of `length`."""
    return position % length if length else 0


def revcomp(seq: str) -> str:
    return str(Seq(seq).reverse_complement()).upper()


def subseq(seq: str, start: int, end: int, circular: bool = True) -> str:
    """Sequence from `start` up to `end`, wrapping over the origin when needed.

    `start == end` on a circular molecule means the whole molecule, not the
    empty string: it is the fragment produced by a single cut.
    """
    n = len(seq)
    if not circular:
        return seq[start:end]
    start, end = wrap(start, n), wrap(end, n)
    if start < end:
        return seq[start:end]
    return seq[start:] + seq[:end]


def span_length(start: int, end: int, length: int) -> int:
    """Length of the circular span from `start` to `end`."""
    return (end - start) % length or length


def slice_record(
    record: SeqRecord,
    start: int,
    end: int,
    circular: bool = True,
    name: str | None = None,
) -> SeqRecord:
    """Cut a sub-record out of a record, carrying its features across.

    Coordinates are the same top-strand cut positions `enzymes.digest` uses, so
    a fragment's record can be built straight from a `Fragment`. Features that
    fall wholly inside come across intact; ones that straddle an end are
    clipped and marked, because a clipped CDS is still worth seeing on a map.
    """
    seq = sequence(record)
    n = len(seq)
    length = span_length(start, end, n) if circular else end - start
    sub = SeqRecord(
        Seq(subseq(seq, start, end, circular)),
        id=name or record.id,
        name=(name or record.name or "fragment")[:16],
        description=record.description,
        annotations={"molecule_type": "DNA", "topology": "linear"},
    )

    for feature in record.features:
        mapped = _map_feature(feature, start, n, length, circular)
        if mapped is not None:
            sub.features.append(mapped)
    return sub


def _map_feature(feature, start: int, n: int, length: int, circular: bool):
    """Move a feature into fragment coordinates, or None if it falls outside."""
    from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation

    pieces = []
    clipped = False
    for part in feature.location.parts:
        begin = (int(part.start) - start) % n if circular else int(part.start) - start
        finish = begin + (int(part.end) - int(part.start))
        if begin >= length or finish <= 0:
            continue
        if begin < 0 or finish > length:
            clipped = True
            begin, finish = max(0, begin), min(length, finish)
        if finish > begin:
            pieces.append(SimpleLocation(begin, finish, part.strand))

    if not pieces:
        return None

    location = pieces[0] if len(pieces) == 1 else CompoundLocation(pieces)
    qualifiers = {k: list(v) for k, v in feature.qualifiers.items()}
    if clipped:
        qualifiers.setdefault("note", []).append("clipped by digest")
    return SeqFeature(location, type=feature.type, qualifiers=qualifiers)
