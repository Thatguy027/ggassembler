"""Programmatic fixtures.

Golden fixtures are built from the part-type definitions rather than checked in
as sequence files, so the tests describe the standard rather than a particular
vendor's plasmid. Anything in ``tests/fixtures/`` is a real file supplied by the
user and is read as-is.
"""

from __future__ import annotations

import random
import zlib

from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from ggassembler.core import parttypes
from ggassembler.core.enzymes import BSAI, BSMBI, Enzyme, find_sites
from ggassembler.core.seqio import revcomp

SPACER = "A"


def seed_for(name: str) -> int:
    """A stable seed from a name.

    Python's `hash()` is salted per process, so using it here made fixtures - and
    therefore the tests - differ between runs.
    """
    return zlib.crc32(name.encode()) % 100_000


def _expected_sites(record: SeqRecord, enzyme: Enzyme, count: int) -> bool:
    return len(find_sites(str(record.seq), enzyme, circular=True)) == count


def filler(length: int, seed: int = 0) -> str:
    """Random-but-deterministic sequence carrying no BsaI, BsmBI or NotI site."""
    rng = random.Random(seed)
    for attempt in range(200):
        seq = "".join(rng.choice("ACGT") for _ in range(length))
        if not _has_any_site(seq):
            return seq
        rng = random.Random(seed + 1000 * (attempt + 1))
    raise RuntimeError("could not generate a site-free filler")


def _has_any_site(seq: str) -> bool:
    from ggassembler.core.enzymes import NOTI

    return any(find_sites(seq, e, circular=False) for e in (BSAI, BSMBI, NOTI))


def _spacer(enzyme: Enzyme) -> str:
    """The filler between a site and its overhang - one base for BsaI and
    BsmBI, two for BbsI. Taken from the enzyme rather than assumed."""
    return SPACER * (enzyme.top_offset - len(enzyme.site))


def forward_site(enzyme: Enzyme) -> str:
    """A site that cuts to its right; the overhang follows immediately after."""
    return enzyme.site + _spacer(enzyme)


def reverse_site(enzyme: Enzyme) -> str:
    """A site that cuts to its left; the overhang sits immediately before."""
    return _spacer(enzyme) + revcomp(enzyme.site)


def part_plasmid(
    part_type: str,
    insert: str | None = None,
    enzyme: Enzyme = BSAI,
    backbone: str | None = None,
    name: str | None = None,
    seed: int = 1,
) -> SeqRecord:
    """A canonical part plasmid: two sites pointing inward, releasing the part.

    Layout, left to right:
        [site ->][5' overhang][insert][3' overhang][<- site][backbone]
    so that digesting releases a site-free fragment whose overhangs are the
    part type's pair.
    """
    five, three = _overhangs(part_type)
    insert = insert if insert is not None else filler(120, seed)
    backbone = backbone if backbone is not None else filler(300, seed + 7)
    label = name or f"synthetic_type{part_type}"
    for attempt in range(50):
        pad = filler(300, seed + 7 + 31 * attempt) if backbone is None else backbone
        body = insert if insert is not None else filler(120, seed + 101 * attempt)
        record = _record(
            forward_site(enzyme) + five + body + three + reverse_site(enzyme) + pad,
            label,
        )
        # a segment boundary can spell an extra site by accident; reject and retry
        if _expected_sites(record, enzyme, 2) or insert is not None:
            return record
    raise RuntimeError(f"could not build a clean part plasmid for {part_type}")


def dropout_plasmid(
    part_type: str,
    enzyme: Enzyme = BSAI,
    dropout: str | None = None,
    backbone: str | None = None,
    name: str | None = None,
    seed: int = 2,
) -> SeqRecord:
    """A destination vector: sites point outward, so the dropout carries them.

    The released fragment is the one holding the recognition sites; it is
    consumed by the reaction and never appears in the product. This is the
    ``234r`` geometry, and the same geometry the part entry vector uses.
    """
    five, three = _overhangs(part_type)
    label = name or f"synthetic_{part_type}r"
    for attempt in range(50):
        middle = dropout if dropout is not None else filler(150, seed + 101 * attempt)
        pad = backbone if backbone is not None else filler(400, seed + 7 + 31 * attempt)
        record = _record(
            five + reverse_site(enzyme) + middle + forward_site(enzyme) + three + pad,
            label,
        )
        if _expected_sites(record, enzyme, 2) or dropout is not None:
            return record
    raise RuntimeError(f"could not build a clean dropout plasmid for {part_type}")


def entry_vector(name: str = "synthetic_entry_vector", seed: int = 3) -> SeqRecord:
    """A pYTK001-style part entry vector: BsmBI dropout with TCGG/GACC ends.

    The backbone is kept plasmid-sized on purpose: a destination has to carry an
    origin and a marker, and the library will not offer a fragment as a vector.
    """
    five, three = parttypes.ENTRY_VECTOR_OVERHANGS
    dropout = filler(150, seed)
    backbone = filler(1800, seed + 7)
    seq = (
        five
        + reverse_site(BSMBI)
        + dropout
        + forward_site(BSMBI)
        + three
        + backbone
    )
    return _record(seq, name)


def dropout_vector(
    five: str,
    three: str,
    enzyme: Enzyme = BSMBI,
    name: str = "synthetic_dropout_vector",
    backbone: int = 1800,
    seed: int = 12,
) -> SeqRecord:
    """A ccdB/GFP-style destination: sites point outward, the dropout carries them.

    Cutting it releases the dropout and leaves a site-free backbone whose ends
    accept an insert presenting `five` -> `three`. The backbone is plasmid-sized,
    because the library will not offer a fragment as a vector.
    """
    seq = (
        five
        + reverse_site(enzyme)
        + filler(400, seed)
        + forward_site(enzyme)
        + three
        + filler(backbone, seed + 5)
    )
    return _record(seq, name)


def connector_plasmid(
    part_type: str,
    connector_overhang: str,
    name: str,
    seed: int = 4,
) -> SeqRecord:
    """A type 1 or 5 part whose insert carries an internal BsmBI site.

    The BsmBI overhang is what orders a Level 3 assembly; ConLX and ConRX with
    the same X share it. It is written into the sequence here, never assumed.
    """
    lead = filler(20, seed)
    tail = filler(20, seed + 1)
    if part_type.startswith("1"):
        # points rightward, out of the connector and into the cassette
        insert = lead + forward_site(BSMBI) + connector_overhang + tail
    else:
        insert = lead + connector_overhang + reverse_site(BSMBI) + tail
    record = part_plasmid(part_type, insert=insert, name=name, seed=seed)
    record.features.append(
        SeqFeature(FeatureLocation(0, len(record.seq)), type="misc_feature",
                   qualifiers={"label": [name]})
    )
    return record


def cassette_plasmid(
    connector_five: str,
    connector_three: str,
    body: str | None = None,
    dropout: str | None = None,
    name: str = "synthetic_cassette",
    seed: int = 8,
) -> SeqRecord:
    """An assembled cassette: BsmBI sites in the connectors point inward.

    Digesting with BsmBI releases the transcription unit flanked by the two
    connector overhangs, site-free - which is what makes it a Level 3 input.
    Pass `dropout` to get a destination vector instead: the same connector ends
    with a BsaI dropout of that type between them.
    """
    if dropout:
        middle = str(dropout_plasmid(dropout, seed=seed).seq)
    else:
        middle = body if body is not None else filler(400, seed)
    seq = (
        forward_site(BSMBI)
        + connector_five
        + middle
        + connector_three
        + reverse_site(BSMBI)
        + filler(300, seed + 1)
    )
    return _record(seq, name)


def plasmid_without_bsai(name: str = "synthetic_no_sites", seed: int = 5) -> SeqRecord:
    return _record(filler(800, seed), name)


def plasmid_with_internal_site(
    part_type: str = "3",
    name: str = "synthetic_internal_bsai",
    seed: int = 6,
) -> SeqRecord:
    """A part plasmid whose insert still carries an undomesticated BsaI site."""
    insert = filler(60, seed) + BSAI.site + filler(60, seed + 1)
    return part_plasmid(part_type, insert=insert, name=name, seed=seed)


def _overhangs(part_type: str) -> tuple[str, str]:
    pair = parttypes.type_overhangs(part_type)
    if pair is None:
        raise ValueError(f"not a YTK part type: {part_type}")
    return pair


def _record(seq: str, name: str) -> SeqRecord:
    record = SeqRecord(
        Seq(seq),
        id=name,
        name=name[:16],
        description=f"{name} (synthetic fixture)",
        annotations={"molecule_type": "DNA", "topology": "circular"},
    )
    return record


def plasmid_with_unreadable_sites(
    name: str = "synthetic_odd_overhangs",
    enzyme: Enzyme = BSAI,
    seed: int = 9,
) -> SeqRecord:
    """A plasmid the part enzyme cuts, but not into anything the scheme knows.

    Two inward-facing sites, so a fragment *is* released - it just does not sit
    on the part circle, so no type can be read off it. This is the case the
    annotation fallback exists for: the digest genuinely has no answer, rather
    than answering "nothing can be cut out of this at all".
    """
    return _record(
        forward_site(enzyme) + "AAAA" + filler(120, seed)
        + "TTTT" + reverse_site(enzyme) + filler(300, seed + 7),
        name,
    )
