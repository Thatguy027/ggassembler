"""Detection: what the digest says a plasmid is, and how sure it is."""

from __future__ import annotations

import json

import pytest
from Bio.SeqFeature import FeatureLocation, SeqFeature

from ggassembler.core.library import INDEX_VERSION, Library, describe, detect_type
from ggassembler.core.seqio import write_genbank

from . import synth

PART_TYPES = ["1", "2", "3", "3a", "3b", "4", "4a", "4b", "5", "6", "7", "8", "8a", "8b"]


@pytest.mark.parametrize("part_type", PART_TYPES)
def test_part_plasmids_are_typed_by_digest(part_type):
    record = synth.part_plasmid(part_type)
    call = detect_type(record)
    assert call.part_type == part_type
    assert call.confidence == "digest"
    assert not call.reversed_sites
    assert call.conflict is None
    assert "BsaI" in call.reason


def test_dropout_is_typed_with_the_r_suffix():
    call = detect_type(synth.dropout_plasmid("234"))
    assert call.part_type == "234r"
    assert call.confidence == "digest"
    assert call.reversed_sites
    assert (call.five_prime, call.three_prime) == ("AACG", "GCTG")


def test_composite_parts_are_named_by_their_span():
    for name in ("234", "678"):
        call = detect_type(synth.part_plasmid(name))
        assert call.part_type == name
        assert call.confidence == "digest"


def test_entry_vector_is_detected_via_bsmbi():
    entry = describe(synth.entry_vector())
    assert entry.is_entry_vector
    assert not entry.is_cassette
    # it has no BsaI part to release
    assert entry.call.part_type is None


def test_connector_overhang_is_read_from_the_file():
    """ConLX and ConRX with the same X share a connector overhang."""
    conl = describe(synth.connector_plasmid("1", "CCAA", name="ConL1"))
    conr = describe(synth.connector_plasmid("5", "CCAA", name="ConR1"))

    assert conl.call.part_type == "1"
    assert conr.call.part_type == "5"
    assert conl.connector_overhang == "CCAA"
    assert conr.connector_overhang == "CCAA"
    assert "connector" in conl.roles


def test_internal_bsai_site_is_counted_but_does_not_break_the_call():
    entry = describe(synth.plasmid_with_internal_site("3"))
    assert entry.call.part_type == "3"
    assert entry.call.confidence == "digest"
    assert entry.sites.part_enzyme_total == 3
    assert entry.sites.part_enzyme_internal == 1
    assert "further BsaI site" in entry.call.reason


def test_plasmid_with_no_bsai_pair_is_unclassified():
    call = detect_type(synth.plasmid_without_bsai())
    assert call.part_type is None
    assert call.confidence == "none"
    assert call.reason


def test_annotation_is_used_only_when_the_digest_fails():
    """The digest having no answer is what lets an annotation speak."""
    record = synth.plasmid_with_unreadable_sites(name="mystery")
    record.features.append(
        SeqFeature(FeatureLocation(0, 10), type="misc_feature",
                   qualifiers={"note": ["YTK type 6 marker"]})
    )
    call = detect_type(record)
    assert call.part_type == "6"
    assert call.confidence == "annotation"
    assert "note" in call.reason


def test_filename_is_the_last_resort():
    record = synth.plasmid_with_unreadable_sites(name="pMYT_type_7_origin")
    call = detect_type(record)
    assert call.part_type == "7"
    assert call.confidence == "filename"


def test_a_plasmid_no_enzyme_cuts_is_never_a_part_however_it_is_labelled():
    """A finished construct keeps the labels of the parts that built it.

    `L6_CYC1_type4` sitting inside a 10.8 kb assembly is a record of what went
    in, not a statement about the whole plasmid - and since every site was
    consumed building it, nothing can be cut out of it at all. Believing the
    label offers a multigene construct as a terminator.
    """
    record = synth.plasmid_without_bsai(name="assembled")
    record.features.append(
        SeqFeature(FeatureLocation(10, 40), type="misc_feature",
                   qualifiers={"label": ["L6_CYC1_type4"]})
    )
    call = detect_type(record)
    assert call.part_type is None
    assert call.confidence == "none"
    # the claim is recorded, not dropped
    assert "type 4" in call.conflict
    assert "not a part" in call.conflict


def test_the_same_holds_when_only_the_file_name_makes_the_claim():
    call = detect_type(synth.plasmid_without_bsai(name="pMYT_type_7_origin"))
    assert call.part_type is None
    assert "file name" in call.conflict


def test_digest_beats_a_disagreeing_annotation_and_records_the_conflict():
    record = synth.part_plasmid("2", name="mislabelled")
    record.features.append(
        SeqFeature(FeatureLocation(0, 10), type="misc_feature",
                   qualifiers={"label": ["type 6 marker"]})
    )
    call = detect_type(record)
    assert call.part_type == "2"
    assert call.confidence == "digest"
    assert call.conflict and "type 6" in call.conflict
    assert "keeping the digest" in call.conflict


def test_agreeing_annotation_raises_no_conflict():
    record = synth.part_plasmid("6", name="marker")
    record.features.append(
        SeqFeature(FeatureLocation(0, 10), type="misc_feature",
                   qualifiers={"label": ["type 6"]})
    )
    assert detect_type(record).conflict is None


def test_ecoli_marker_and_features_are_recorded():
    record = synth.part_plasmid("8", name="backbone")
    record.features.append(
        SeqFeature(FeatureLocation(0, 30), type="CDS", qualifiers={"label": ["AmpR"]})
    )
    entry = describe(record)
    assert entry.ecoli_marker == "AmpR"
    assert any(f.label == "AmpR" for f in entry.features)
    assert entry.topology == "circular"


# --------------------------------------------------------------------------- #
# the indexed folder
# --------------------------------------------------------------------------- #


@pytest.fixture
def folder(tmp_path):
    write_genbank(synth.part_plasmid("2", name="promoter"), tmp_path / "promoter.gb")
    write_genbank(synth.part_plasmid("3", name="cds"), tmp_path / "cds.gb")
    write_genbank(synth.dropout_plasmid("234", name="dropout"), tmp_path / "sub" / "dropout.gb")
    return tmp_path


def test_scan_indexes_recursively(folder):
    entries = Library(folder).scan()
    assert {e.call.part_type for e in entries} == {"2", "3", "234r"}


def test_no_recursive_stops_at_the_top_level(folder):
    entries = Library(folder, recursive=False).scan()
    assert {e.call.part_type for e in entries} == {"2", "3"}


def test_index_is_cached_and_reused(folder):
    library = Library(folder)
    library.scan()
    assert library.index_path.exists()

    cached = json.loads(library.index_path.read_text())
    assert cached["version"] == INDEX_VERSION
    assert "promoter.gb" in cached["entries"]

    # a second scan reads the cache rather than the file
    (folder / "promoter.gb").rename(folder / "hidden.gb")
    (folder / "hidden.gb").rename(folder / "promoter.gb")
    again = Library(folder).scan()
    assert len(again) == 3


def test_changed_file_is_rescanned(folder):
    library = Library(folder)
    library.scan()
    write_genbank(synth.part_plasmid("6", name="promoter"), folder / "promoter.gb")

    entries = {e.path: e for e in Library(folder).scan()}
    assert entries["promoter.gb"].call.part_type == "6"


def test_manual_override_wins_and_is_marked(folder):
    library = Library(folder)
    library.scan()
    library.set_override("cds.gb", "3a", reason="it is really the N-terminal half")

    entry = {e.path: e for e in Library(folder).scan()}["cds.gb"]
    assert entry.call.part_type == "3a"
    assert entry.call.source == "manual"
    assert entry.call.conflict and "detection said 3" in entry.call.conflict


def test_override_can_be_cleared(folder):
    library = Library(folder)
    library.scan()
    library.set_override("cds.gb", "3a")
    library.set_override("cds.gb", None)
    entry = {e.path: e for e in Library(folder).scan()}["cds.gb"]
    assert entry.call.part_type == "3"
    assert entry.call.source == "detected"


def test_unreadable_file_is_reported_not_fatal(folder):
    (folder / "broken.gb").write_text("this is not a GenBank file")
    library = Library(folder)
    entries = library.scan()
    assert len(entries) == 3
    assert "broken.gb" in library.errors


def test_queries(folder):
    library = Library(folder)
    library.scan()
    assert [e.name for e in library.parts_of_type("2")] == ["promoter"]
    assert [e.name for e in library.parts_with_overhangs("AACG", "TATG")] == ["promoter"]
    # a dropout is not offered as a part, even though its span matches
    assert library.parts_with_overhangs("AACG", "GCTG") == []


# --------------------------------------------------------------------------- #
# assembled constructs are not parts
# --------------------------------------------------------------------------- #


def test_cassette_is_detected_by_its_connector_ends():
    entry = describe(synth.cassette_plasmid("CCAA", "GATG"))
    assert entry.is_cassette
    assert entry.cassette_overhangs == ("CCAA", "GATG")
    assert entry.call.part_type is None


def test_a_label_inside_a_cassette_does_not_type_the_whole_plasmid():
    record = synth.cassette_plasmid("CCAA", "GATG")
    record.features.append(
        SeqFeature(FeatureLocation(10, 40), type="terminator",
                   qualifiers={"label": ["type 4 terminator"]})
    )
    entry = describe(record)
    assert entry.call.part_type is None
    assert entry.call.confidence == "none"
    assert "assembled cassette, not a part" in entry.call.reason


def test_a_dropout_between_connectors_is_a_multigene_vector():
    entry = describe(synth.cassette_plasmid("AGCA", "CTGA", dropout="234"))
    assert entry.call.part_type == "234r"
    assert entry.is_multigene_vector
    assert not entry.is_cassette
    assert entry.cassette_overhangs == ("AGCA", "CTGA")


# --------------------------------------------------------------------------- #
# a file can be wrong about its own topology
# --------------------------------------------------------------------------- #


def test_a_plasmid_labelled_linear_is_still_typed():
    """Plenty of tools write 'linear' by default; the digest knows better."""
    record = synth.part_plasmid("3", name="says_linear")
    record.annotations["topology"] = "linear"

    call = detect_type(record)
    assert call.part_type == "3"
    assert call.confidence == "digest"
    assert call.assumed_circular
    assert call.conflict and "says this record is linear" in call.conflict


def test_a_reinterpreted_record_is_reported_as_circular():
    record = synth.part_plasmid("6", name="says_linear")
    record.annotations["topology"] = "linear"
    entry = describe(record)

    assert entry.topology == "circular"
    assert entry.sites.part_enzyme_total == 2
    assert entry.call.conflict


def test_a_genuinely_linear_fragment_is_not_forced_circular():
    """A PCR product must not be typed on a site that only the join creates."""
    from ggassembler.core.enzymes import BSAI
    from ggassembler.core.seqio import revcomp

    # ends that only form a pair if you glue them together
    seq = ("A" + revcomp(BSAI.site) + synth.filler(200, 9)
           + "TATG" + synth.filler(150, 10) + BSAI.site + "A")
    record = synth._record(seq, "pcr_product")
    record.annotations["topology"] = "linear"

    call = detect_type(record)
    assert not call.assumed_circular
    assert call.part_type is None


def test_a_linear_record_with_no_sites_stays_unclassified():
    record = synth.plasmid_without_bsai(name="linear_blob")
    record.annotations["topology"] = "linear"
    call = detect_type(record)
    assert call.part_type is None
    assert not call.assumed_circular


def test_an_explicitly_circular_record_carries_no_topology_conflict():
    assert detect_type(synth.part_plasmid("3")).conflict is None


# --------------------------------------------------------------------------- #
# the same plasmid filed twice is one plasmid
# --------------------------------------------------------------------------- #


def test_duplicate_files_merge_into_one_entry_keeping_both_names(tmp_path):
    record = synth.part_plasmid("2", name="pTDH3")
    write_genbank(record, tmp_path / "pYTK009.gb")
    write_genbank(record, tmp_path / "A9_TDH3_promoter.gbk")
    write_genbank(synth.part_plasmid("3", name="cds"), tmp_path / "cds.gb")

    library = Library(tmp_path)
    library.scan()

    assert len(library.entries) == 3, "both files are still indexed"
    assert len(library.unique_entries()) == 2, "but they are one plasmid"

    merged = library.get("pYTK009")
    assert merged.aliases == ["A9_TDH3_promoter"]
    assert "pYTK009" in merged.display and "A9_TDH3_promoter" in merged.display


def test_either_name_finds_the_merged_plasmid(tmp_path):
    record = synth.part_plasmid("2", name="pTDH3")
    write_genbank(record, tmp_path / "pYTK009.gb")
    write_genbank(record, tmp_path / "A9_TDH3_promoter.gbk")
    library = Library(tmp_path)
    library.scan()

    assert library.get("pYTK009") is library.get("A9_TDH3_promoter")


def test_a_merged_plasmid_appears_once_in_a_dropdown(tmp_path):
    record = synth.part_plasmid("2", name="pTDH3")
    for name in ("pYTK009.gb", "A9_TDH3_promoter.gbk", "spare_copy.gb"):
        write_genbank(record, tmp_path / name)
    library = Library(tmp_path)
    library.scan()

    matches = library.parts_with_overhangs("AACG", "TATG")
    assert len(matches) == 1
    assert sorted(matches[0].aliases) == ["A9_TDH3_promoter", "spare_copy"]


def test_the_canonical_copy_speaks_for_the_group(tmp_path):
    """A .gb in the first root outranks a .gbk copy elsewhere."""
    record = synth.part_plasmid("2", name="pTDH3")
    write_genbank(record, tmp_path / "pYTK009.gb")
    (tmp_path / "other").mkdir()
    write_genbank(record, tmp_path / "other" / "A9_TDH3_promoter.gbk")

    library = Library([tmp_path, tmp_path / "other"])
    library.scan()
    assert library.unique_entries()[0].name == "pYTK009"


def test_a_twins_annotations_fill_in_a_blank_description(tmp_path):
    from Bio.SeqFeature import FeatureLocation, SeqFeature

    bare = synth.part_plasmid("2", name="bare")
    described = synth.part_plasmid("2", name="described")
    described.features.append(
        SeqFeature(FeatureLocation(10, 120), type="promoter",
                   qualifiers={"label": ["ScTDH3 Promoter"]})
    )
    write_genbank(bare, tmp_path / "pYTK009.gb")
    write_genbank(described, tmp_path / "A9_TDH3_promoter.gb")

    library = Library(tmp_path)
    library.scan()
    merged = library.get("pYTK009")
    assert merged.component == "ScTDH3 Promoter"


def test_a_component_is_read_for_every_part_type(tmp_path):
    from Bio.SeqFeature import FeatureLocation, SeqFeature

    record = synth.part_plasmid("4", name="tADH1")
    record.features.append(
        SeqFeature(FeatureLocation(20, 90), type="terminator",
                   qualifiers={"label": ["ScADH1 Terminator"]})
    )
    write_genbank(record, tmp_path / "tADH1.gb")
    library = Library(tmp_path)
    library.scan()

    entry = library.get("tADH1")
    assert entry.component == "ScADH1 Terminator"
    assert entry.display == "tADH1 [ScADH1 Terminator]"


def test_provenance_and_primer_labels_are_not_components():
    from ggassembler.core.library import FeatureInfo, summarise

    noise = [
        FeatureInfo("misc_feature", "pYTK013", 0, 50, 1),
        FeatureInfo("misc_feature", "Con1 scar", 0, 20, 1),
        FeatureInfo("primer", "GPDpro-F", 0, 20, 1),
        FeatureInfo("protein_bind", "BsmBI", 0, 6, 1),
    ]
    assert summarise(noise) == ""
    assert summarise([*noise, FeatureInfo("promoter", "ScTDH3 Promoter", 0, 700, 1)]) == (
        "ScTDH3 Promoter"
    )


# --------------------------------------------------------------------------- #
# cassettes that a stray internal site would otherwise hide
# --------------------------------------------------------------------------- #


def _with_internal_bsmbi(left, right, name, seed=200):
    """A cassette whose transcription unit also carries a stray BsmBI site."""
    from ggassembler.core.enzymes import BSMBI

    from .synth import filler, forward_site
    body = filler(300, seed) + forward_site(BSMBI) + "TTAG" + filler(300, seed + 1)
    return synth.cassette_plasmid(left, right, body=body, name=name, seed=seed)


def test_a_cassette_with_a_stray_internal_site_is_still_recognised(tmp_path):
    """One extra site used to hide the whole cassette; now it is a warning."""
    write_genbank(synth.connector_plasmid("1", "CCAA", name="ConL1"), tmp_path / "ConL1.gb")
    write_genbank(synth.connector_plasmid("5", "GATG", name="ConR2"), tmp_path / "ConR2.gb")
    write_genbank(_with_internal_bsmbi("CCAA", "GATG", "pStray"), tmp_path / "pStray.gb")

    library = Library(tmp_path)
    library.scan()

    entry = library.get("pStray")
    assert entry.is_cassette
    assert entry.cassette_overhangs == ("CCAA", "GATG")
    assert entry.internal_multigene_positions, "the stray site must be located"
    assert entry.call.conflict and "cut this cassette apart" in entry.call.conflict
    assert "domesticate" in entry.call.conflict


def test_the_stray_site_position_is_reported(tmp_path):
    write_genbank(synth.connector_plasmid("1", "CCAA", name="ConL1"), tmp_path / "ConL1.gb")
    write_genbank(synth.connector_plasmid("5", "GATG", name="ConR2"), tmp_path / "ConR2.gb")
    record = _with_internal_bsmbi("CCAA", "GATG", "pStray")
    write_genbank(record, tmp_path / "pStray.gb")

    library = Library(tmp_path)
    library.scan()
    entry = library.get("pStray")

    from ggassembler.core.enzymes import BSMBI, find_sites
    real = {s.start for s in find_sites(str(record.seq), BSMBI)}
    for position in entry.internal_multigene_positions:
        assert position in real
    assert str(entry.internal_multigene_positions[0]) in entry.call.reason.replace(",", "")


def test_recovery_needs_a_library_that_knows_its_connectors(tmp_path):
    """With no connector parts indexed, there is no vocabulary to match against."""
    write_genbank(_with_internal_bsmbi("CCAA", "GATG", "pStray"), tmp_path / "pStray.gb")
    library = Library(tmp_path)
    library.scan()
    assert not library.get("pStray").is_cassette


def test_a_clean_cassette_reports_no_internal_sites(tmp_path):
    write_genbank(synth.connector_plasmid("1", "CCAA", name="ConL1"), tmp_path / "ConL1.gb")
    write_genbank(synth.cassette_plasmid("CCAA", "GATG", name="pClean", seed=7),
                  tmp_path / "pClean.gb")
    library = Library(tmp_path)
    library.scan()

    entry = library.get("pClean")
    assert entry.is_cassette
    assert entry.internal_multigene_positions == []
    assert entry.call.conflict is None


def test_a_connector_part_is_not_faulted_for_its_own_site(tmp_path):
    """A type 1 part is supposed to carry exactly one BsmBI site."""
    write_genbank(synth.connector_plasmid("1", "CCAA", name="ConL1"), tmp_path / "ConL1.gb")
    library = Library(tmp_path)
    library.scan()

    entry = library.get("ConL1")
    assert entry.connector_overhang == "CCAA"
    assert entry.internal_multigene_positions == []


def test_multigene_overhangs_are_recorded_even_without_a_call(tmp_path):
    write_genbank(_with_internal_bsmbi("CCAA", "GATG", "pStray"), tmp_path / "pStray.gb")
    library = Library(tmp_path)
    library.scan()

    entry = library.get("pStray")
    assert len(entry.multigene_overhangs) >= 3
    assert "CCAA" in entry.multigene_overhangs


def test_a_hand_assigned_type_reaches_the_view_every_screen_reads(tmp_path):
    """Overrides used to update `entries` but not the deduplicated view.

    `_primary` holds references into `entries`, so replacing an entry left the
    merged view pointing at the old object - and every dropdown, the library
    table and the search all read the merged view. The override was written to
    disk, reported back by the API, and invisible everywhere it mattered.
    """
    write_genbank(synth.part_plasmid("2", name="promoter"), tmp_path / "promoter.gb")
    library = Library(tmp_path)
    library.scan()
    assert library.get("promoter").call.part_type == "2"

    library.set_override("promoter.gb", "3", reason="checked by sequencing")
    assert library.get("promoter").call.part_type == "3"
    assert [e.call.part_type for e in library.unique_entries()] == ["3"]
    assert library.parts_of_type("3"), "the override must reach the type index too"


def test_clearing_an_override_puts_the_detected_type_back(tmp_path):
    write_genbank(synth.part_plasmid("2", name="promoter"), tmp_path / "promoter.gb")
    library = Library(tmp_path)
    library.scan()
    library.set_override("promoter.gb", "3")
    library.set_override("promoter.gb", None)
    assert library.get("promoter").call.part_type == "2"


def test_an_archive_folder_is_not_indexed(tmp_path):
    """Redundant copies are moved aside rather than deleted, so the indexer has
    to agree not to look there - otherwise archiving changes nothing."""
    from ggassembler.core.library import ARCHIVE_DIR

    record = synth.part_plasmid("2", name="promoter")
    write_genbank(record, tmp_path / "promoter.gb")
    archived = tmp_path / ARCHIVE_DIR / "plasmids"
    archived.mkdir(parents=True)
    write_genbank(record, archived / "promoter_old_name.gb")

    library = Library(tmp_path)
    library.scan()
    assert len(library.entries) == 1
    assert library.get("promoter_old_name") is None
    assert library.duplicates() == {}


def test_moving_a_file_back_out_of_the_archive_restores_it(tmp_path):
    from ggassembler.core.library import ARCHIVE_DIR

    record = synth.part_plasmid("2", name="promoter")
    archived = tmp_path / ARCHIVE_DIR
    archived.mkdir()
    write_genbank(record, archived / "promoter.gb")

    library = Library(tmp_path)
    library.scan()
    assert library.entries == {}

    (archived / "promoter.gb").rename(tmp_path / "promoter.gb")
    library = Library(tmp_path)
    library.scan()
    assert library.get("promoter") is not None


def test_a_sequence_verified_copy_speaks_for_its_group(tmp_path):
    """Between two identical files the verified one is the one worth keeping.

    The old rule preferred the shortest name, which archived
    `gal2_cassette_xylb_Seq_verified.gbk` in favour of `gal2_cassette.gb` -
    keeping the unverified file and discarding the provenance.
    """
    record = synth.part_plasmid("2", name="thing")
    write_genbank(record, tmp_path / "thing.gb")
    write_genbank(record, tmp_path / "thing_seq_verified.gb")

    library = Library(tmp_path)
    library.scan()
    entries = library.unique_entries()
    assert len(entries) == 1
    assert entries[0].name == "thing_seq_verified"
    assert entries[0].aliases == ["thing"]


def test_a_hand_chosen_canonical_still_beats_a_verified_name(tmp_path):
    record = synth.part_plasmid("2", name="thing")
    write_genbank(record, tmp_path / "thing.gb")
    write_genbank(record, tmp_path / "thing_seq_verified.gb")
    library = Library(tmp_path)
    library.scan()
    library.set_canonical("thing.gb")
    assert library.unique_entries()[0].name == "thing"


# --------------------------------------------------------------------------- #
# readiness for the level above
# --------------------------------------------------------------------------- #


def test_a_connector_that_kept_its_multigene_site_is_level3_ready(tmp_path):
    from ggassembler.core.enzymes import BSMBI

    insert = (
        synth.filler(40, 11) + synth.forward_site(BSMBI) + "CCAA" + synth.filler(40, 12)
    )
    record = synth.part_plasmid("1", insert=insert, name="ConLS")
    write_genbank(record, tmp_path / "ConLS.gb")
    library = Library(tmp_path)
    library.scan()
    entry = library.get("ConLS")
    assert entry.connector_overhang, "this fixture should carry a connector"
    assert entry.level3_ready is True


def test_a_connector_with_no_multigene_site_blocks_the_level_above(tmp_path):
    """`L13_ConLS_BSMB1del` is the real case: its own label reads "Former Bsmb1".

    It assembles at Level 2 like any other type 1 part and produces a cassette
    that cannot be cut out of its plasmid, which nothing surfaced until the
    multigene step failed to find the ends.
    """
    from ggassembler.core.enzymes import BSMBI, find_sites

    record = synth.part_plasmid("1", insert=synth.filler(120, 42), name="ConLS_del")
    assert not find_sites(str(record.seq), BSMBI), "fixture must have no BsmBI site"
    write_genbank(record, tmp_path / "ConLS_del.gb")

    library = Library(tmp_path)
    library.scan()
    entry = library.get("ConLS_del")
    assert entry.call.part_type == "1"
    assert entry.connector_overhang is None
    assert entry.level3_ready is False


def test_the_question_does_not_arise_for_other_positions(tmp_path):
    """Only positions 1 and 5 contribute to a cassette's ends, so only they can
    block the level above. A promoter is neither ready nor blocking."""
    for part_type, name in (("2", "promoter"), ("3", "cds"), ("6", "marker")):
        write_genbank(synth.part_plasmid(part_type, name=name), tmp_path / f"{name}.gb")
    library = Library(tmp_path)
    library.scan()
    for name in ("promoter", "cds", "marker"):
        assert library.get(name).level3_ready is None


def test_readiness_survives_the_cache(tmp_path):
    write_genbank(synth.part_plasmid("5", name="ConR1"), tmp_path / "ConR1.gb")
    first = Library(tmp_path)
    first.scan()
    expected = first.get("ConR1").level3_ready
    assert expected is not None

    second = Library(tmp_path)
    second.scan()  # reads the index written above rather than re-digesting
    assert second.get("ConR1").level3_ready == expected


# ------------------------------------------------- acceptors with no BsaI ---


def test_a_backbone_with_no_part_sites_is_still_a_backbone(tmp_path):
    """The rule asked whether the *BsaI* sites were reversed. A destination
    vector with no BsaI sites at all cannot answer that, so every one of them
    was filed as a cassette and never offered as somewhere to assemble into -
    which is what made eleven acceptors in a real library invisible.

    What separates the two is which piece comes out. A cassette releases its
    transcription unit and keeps the E. coli backbone; an acceptor releases the
    backbone and keeps a stuffer.
    """
    acceptor = synth.cassette_plasmid(
        "AGCA", "CTGA", body=synth.filler(3200, seed=7), name="pAcceptor", seed=12)
    write_genbank(acceptor, tmp_path / "pAcceptor.gb")
    lib = Library(tmp_path)
    lib.scan()

    entry = lib.get("pAcceptor")
    assert entry.sites.part_enzyme_total == 0, "this test is about plasmids with no BsaI"
    assert entry.is_multigene_vector, "a backbone that releases its backbone is a backbone"


def test_a_small_plasmid_is_never_read_as_a_backbone(tmp_path):
    """A 1.3 kb synthetic cassette carrying 318 bp of backbone releases 76% of
    itself, which is acceptor-shaped by fraction alone. A real destination
    vector carries a marker, an origin and an E. coli backbone and cannot be
    that small, so size is the check that tells them apart."""
    from ggassembler.core.library import ACCEPTOR_MINIMUM

    # acceptor-shaped by fraction, and far too small to be one
    small = synth.cassette_plasmid(
        "AGCA", "CTGA", body=synth.filler(900, seed=5), name="tiny", seed=61)
    write_genbank(small, tmp_path / "tiny.gb")
    lib = Library(tmp_path)
    lib.scan()

    entry = lib.get("tiny")
    released = (entry.cassette_span[1] - entry.cassette_span[0]) % entry.length
    assert released / entry.length > 0.70, "fixture does not reach the fraction rule"
    assert entry.length < ACCEPTOR_MINIMUM
    assert not entry.is_multigene_vector, "size is what stops this being a backbone"


def test_a_cassette_that_keeps_its_backbone_stays_a_cassette(tmp_path):
    """The other side of the same rule: a large insert does not make a cassette
    into a backbone. Forty-eight of them crossed a 50% threshold on insert size
    alone before it was raised."""
    cassette = synth.cassette_plasmid(
        "CTGA", "CCAA", body=synth.filler(900, seed=3), name="pBig", seed=13)
    write_genbank(cassette, tmp_path / "pBig.gb")
    # and a backbone's worth of sequence on the retained side
    lib = Library(tmp_path)
    lib.scan()
    entry = lib.get("pBig")
    released = (entry.cassette_span[1] - entry.cassette_span[0]) % entry.length
    assert released / entry.length < 0.95, "fixture does not exercise the rule"


def test_a_cached_index_does_not_outlive_a_change_to_detection(tmp_path):
    """The cache is keyed on each file's mtime and size, so a rule change
    applied to unchanged files never runs unless the version says so. The
    acceptor rule shipped once without this, which would have meant installing
    the new code and still reading the old classification."""
    import json

    from ggassembler.core.library import INDEX_VERSION

    write_genbank(synth.part_plasmid("3", name="pOne"), tmp_path / "pOne.gb")
    lib = Library(tmp_path)
    lib.scan()
    assert lib.index_path.exists()

    stale = json.loads(lib.index_path.read_text())
    stale["version"] = INDEX_VERSION - 1
    lib.index_path.write_text(json.dumps(stale))

    again = Library(tmp_path)
    again.scan()
    assert json.loads(again.index_path.read_text())["version"] == INDEX_VERSION, (
        "an index written by an older detection rule was kept")
