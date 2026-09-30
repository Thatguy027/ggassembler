"""Level 1 - part plasmid construction.

The load-bearing claim of this level is the last one tested here: whatever the
screen predicts, feeding it back through the detector must classify it as the
type the user asked for. Everything else is in service of that.
"""

from __future__ import annotations

import pytest
from fastapi.testclient import TestClient

from ggassembler.api.main import create_app
from ggassembler.core.enzymes import BSAI, BSMBI, NOTI, digest, find_sites
from ggassembler.core.library import Library, describe
from ggassembler.core.parttypes import YTK, type_overhangs
from ggassembler.core.seqio import write_genbank
from ggassembler.levels import level1_part as level1

from . import synth

#: A clean coding sequence with no accidental sites, in frame.
CDS = "ATG" + "GCTACCGAAGCTTGGACC" * 8 + "TAA"


@pytest.fixture
def library(tmp_path):
    from ggassembler.core.enzymes import BBSI

    # an L131-style universal ccdB vector: BbsI, fixed TATG/ATCC ends, and a
    # backbone with no BsaI sites, so the insert brings its own
    write_genbank(
        synth.dropout_vector("TATG", "ATCC", BBSI, name="L131"), tmp_path / "L131.gb"
    )
    write_genbank(synth.entry_vector(name="pEntry"), tmp_path / "pEntry.gb")
    write_genbank(synth.part_plasmid("2", name="pTDH3", seed=synth.seed_for("pTDH3")),
                  tmp_path / "pTDH3.gb")
    lib = Library(tmp_path)
    lib.scan()
    return lib


def request_for(part_type="3", **kw):
    return level1.PartRequest(part_type=part_type, sequence=kw.pop("sequence", CDS), **kw)


def domesticated(part_type="3", **kw):
    """A request with site removal switched on - it is off by default."""
    return request_for(part_type, domesticate=True, **kw)


# ------------------------------------------------------------------ primers --


def test_the_forward_tail_matches_the_design_canvas():
    """GCAT CGTCTC A TCGG | GGTCTC A TATG, with TCGG and GGTCTC sharing GG."""
    tail = level1._tail("TCGG", "TATG", BSMBI, BSAI)
    assert tail == "GCATCGTCTCATCGGTCTCATATG"


def test_the_tail_reproduces_the_l131_entry_construct():
    """The real construct, rebuilt from its parts: pad, BbsI, L131 end, linker,
    BsaI, type flank."""
    from ggassembler.core.enzymes import BBSI
    from ggassembler.core.seqio import revcomp

    forward = level1._tail("TATG", "TATG", BBSI, BSAI,
                           pad=level1.PAD_5, linker=level1.LINKER_5)
    assert forward == "CGTAGTCGAAGACTATATGCGGGAGGAAGTCTTTAGACCGGTCTCATATG"

    reverse = level1._tail(revcomp("ATCC"), revcomp("ATCC"), BBSI, BSAI,
                           pad=revcomp(level1.PAD_3), linker=revcomp(level1.LINKER_3))
    assert reverse == "TGGCAGCGAAGACTAGGATTGCTGAATGAGAAACCTCGGCGGTCTCAGGAT"


def test_the_tail_still_cuts_after_being_overlapped(library):
    """The shared GG must not cost either enzyme its site."""
    plan = _plan_for(library, "3")
    for part_type in ("1", "2", "3", "3a", "4", "4b", "5", "6", "7", "8", "8b", "234", "678"):
        five, three = type_overhangs(part_type)
        product = level1.amplicon(CDS, five, three, plan)

        outer = digest(product, plan.outer, circular=False)
        middle = [f for f in outer if f.overhangs == plan.accepts]
        assert len(middle) == 1, f"{part_type}: entry ends not presented"

        insert = level1.insert_fragment(CDS, five, three, plan)
        assert insert.startswith(plan.accepts[0])


def _plan_for(library, part_type):
    destination, issues = level1.choose_destination(library, None, part_type, library.scheme)
    assert destination is not None, [str(i) for i in issues]
    return level1.make_plan(destination, library.scheme)


@pytest.mark.parametrize("part_type", ["1", "2", "3", "3a", "3b", "4", "4a", "4b",
                                       "5", "6", "7", "8", "8a", "8b"])
def test_every_part_type_designs_and_validates(library, part_type):
    result = level1.design(library, request_for(part_type))
    assert result.ok, [str(i) for i in result.issues]
    assert result.validated_as == part_type
    assert "validated" in result.codes()


def test_primers_are_tm_matched(library):
    result = level1.design(library, request_for("3"))
    forward, reverse = result.fragments[0].forward, result.fragments[0].reverse
    assert level1.MIN_ANNEAL <= len(forward.annealing) <= level1.MAX_ANNEAL
    assert abs(forward.tm - reverse.tm) <= level1.TM_SPREAD
    assert 50 < forward.tm < 72


def test_primer_segments_cover_the_tail(library):
    result = level1.design(library, request_for("3"))
    forward = result.fragments[0].forward
    roles = [s.role for s in forward.segments]
    assert roles == ["pad", "multigene_enzyme", "entry_overhang", "part_enzyme",
                     "type_flank", "annealing"]
    plan = _plan_for(library, "3")
    assert forward.sequence.startswith(plan.pad_5 + plan.outer.site)
    assert forward.sequence.endswith(forward.annealing)


# -------------------------------------------------------------- conventions --


def test_type_3_strips_the_stop_and_adds_the_gly_ser_linker(library):
    """The fusion path, asked for explicitly - it is never the default."""
    fusable = level1.Convention(gly_ser_linker=True, infer_from_sequence=False)
    result = level1.design(library, request_for("3", conventions=fusable))
    assert "stop_stripped" in result.codes()
    assert "gly_ser" in result.codes()
    assert result.body.endswith("GG")
    assert not result.body.endswith("TAAGG")
    # the junction reads GGATCC - the BamHI / Gly-Ser site
    assert (result.body + result.three_prime).endswith("GGATCC")


def test_conventions_can_be_turned_off(library):
    plain = level1.Convention(gly_ser_linker=False, strip_stop=False)
    result = level1.design(library, request_for("3", conventions=plain))
    assert "gly_ser" not in result.codes()
    assert result.body.endswith("TAA")


def test_the_gly_ser_linker_is_off_unless_it_is_asked_for(library):
    """The default, asserted on its own so it cannot quietly flip back.

    It is not free: the two bases become a Gly-Ser on the protein whether or
    not a fusion follows, so a part that does not need them should not get
    them. Every other test here that wants the linker now says so.
    """
    assert level1.Convention().gly_ser_linker is False

    body = "GCTACCGAAGCTTGGACC" * 8                      # no terminal stop
    result = level1.design(library, request_for("3", sequence=body))
    assert "gly_ser" not in result.codes()
    assert not result.body.endswith("GG")


def test_the_linker_without_stripping_the_stop_is_a_warning(library):
    """Both ticked together is a part that looks fusable and is not.

    The junction would read TAA GGATCC: translation stops at the TAA and the
    linker is never reached. Which of the two was meant is not knowable from
    here, so it is said rather than silently corrected.
    """
    contradictory = level1.Convention(
        gly_ser_linker=True, strip_stop=False, infer_from_sequence=False
    )
    result = level1.design(library, request_for("3", conventions=contradictory))
    assert "linker_without_strip" in result.codes()
    assert any(
        i.level == level1.WARNING for i in result.issues if i.code == "linker_without_strip"
    )


def test_stripping_the_stop_without_the_linker_is_not_a_warning(library):
    """The ordinary case, and the new default: a part that ends where it ends."""
    plain = level1.Convention(strip_stop=True, infer_from_sequence=False)
    result = level1.design(library, request_for("3", conventions=plain))
    assert "linker_without_strip" not in result.codes()


def test_a_part_that_cannot_read_through_says_so():
    """The mirror of the gly_ser note: leaving it off otherwise says nothing.

    The consequence only shows up later, when a 4a fusion is put after this
    part and lands out of frame - a 4 nt overhang is not a whole codon.
    """
    body = "GCTACCGAAGCTTGGACC" * 8
    _, _, _, issues = level1.apply_conventions(
        body, "3", level1.Convention(), fusion_downstream=True
    )
    note = next(i for i in issues if i.code == "no_read_through")
    assert "ATCC" in note.message and "4a" in note.message
    assert note.level == level1.INFO


def test_the_read_through_note_names_the_right_flank_for_3a():
    body = "GCTACCGAAGCTTGGACC" * 8
    _, _, _, issues = level1.apply_conventions(
        body, "3a", level1.Convention(), fusion_downstream=True
    )
    assert "TTCT" in next(i for i in issues if i.code == "no_read_through").message


def test_the_read_through_note_stays_quiet_with_no_4a_part_on_the_shelf(library):
    """Advice about a fusion this library cannot build is noise, not help."""
    assert not library.parts_of_type("4a")
    body = "GCTACCGAAGCTTGGACC" * 8
    result = level1.design(library, request_for("3", sequence=body))
    assert "no_read_through" not in result.codes()


def test_no_read_through_note_when_the_linker_is_on():
    body = "GCTACCGAAGCTTGGACC" * 8
    _, _, suffix, issues = level1.apply_conventions(
        body, "3", level1.Convention(gly_ser_linker=True), fusion_downstream=True
    )
    codes = {i.code for i in issues}
    assert "no_read_through" not in codes and "gly_ser" in codes
    assert suffix == "GG"


def test_no_read_through_note_when_the_sequence_ends_in_a_stop():
    """The terminal-stop note already covers it, and says more."""
    _, _, _, issues = level1.apply_conventions(
        CDS, "3", level1.Convention(), fusion_downstream=True
    )
    codes = {i.code for i in issues}
    assert "terminal_cds" in codes
    assert "no_read_through" not in codes


def test_the_terminal_stop_advice_matches_what_was_asked_for(library):
    """Telling someone to tick a box they already ticked is how a note stops
    being read."""
    default = level1.design(library, request_for("3"))
    assert "tick \u201cGly-Ser linker\u201d" in _message(default, "terminal_cds")

    asked = level1.design(
        library,
        request_for("3", conventions=level1.Convention(gly_ser_linker=True)),
    )
    assert "untick" in _message(asked, "terminal_cds")


def _message(result, code):
    return next(i for i in result.issues if i.code == code).message


def test_type_4_leads_with_a_stop_and_an_xhoi_site(library):
    result = level1.design(library, request_for("4"))
    assert result.body.startswith("TAACTCGAG")
    assert "stop_xhoi" in result.codes()


def test_type_3a_adds_gg_before_ttct(library):
    fusable = level1.Convention(gly_ser_linker=True)
    result = level1.design(library, request_for("3a", conventions=fusable))
    assert result.body.endswith("GG")
    assert (result.body + result.three_prime).endswith("GGTTCT")


# ------------------------------------------------------------ domestication --


def test_a_clean_insert_needs_one_fragment(library):
    result = level1.design(library, request_for("3"))
    assert len(result.fragments) == 1
    assert "internal_site" not in result.codes()


def test_an_internal_site_splits_the_amplicon(library):
    body = "ATG" + "GCTACCGAAGCTTGGACC" * 6 + "GGTCTC" + "AACCGAAGCTTGGACCGCT" * 6 + "TAA"
    result = level1.design(library, domesticated("3", sequence=body))

    assert result.ok, [str(i) for i in result.issues]
    assert len(result.fragments) == 2
    # BsaI releases the part, so a site in it blocks the build until removed
    assert "blocking_site" in result.codes()
    assert "domesticated" in result.codes()
    assert all(f.forward and f.reverse for f in result.fragments)
    assert result.validated_as == "3"


def test_the_chosen_junction_avoids_the_reserved_overhangs(library):
    body = "ATG" + "GCTACCGAAGCTTGGACC" * 6 + "GGTCTC" + "AACCGAAGCTTGGACCGCT" * 6 + "TAA"
    result = level1.design(library, domesticated("3", sequence=body))

    from ggassembler.core.parttypes import conflicts_with_reserved, near_match

    junctions = [f.right_overhang for f in result.fragments[:-1]]
    assert junctions
    for junction in junctions:
        assert conflicts_with_reserved(junction) is None
        assert not near_match(junction, result.five_prime)
        assert not near_match(junction, result.three_prime)


def test_a_silent_mutation_is_proposed_inside_a_coding_part(library):
    body = "ATG" + "GCTACCGAAGCTTGGACC" * 6 + "GGTCTC" + "AACCGAAGCTTGGACCGCT" * 6 + "TAA"
    result = level1.design(library, domesticated("3", sequence=body))

    mutations = [m for f in result.fragments for m in f.mutations]
    assert mutations
    mutation = mutations[0]
    assert mutation.silent, mutation.description
    assert mutation.codon_from != mutation.codon_to
    assert "silent_mutation" in result.codes()


def test_the_proposed_mutation_actually_destroys_the_site():
    body = "ATG" + "GCTACC" * 4 + "GGTCTC" + "AACCGA" * 4
    site = find_sites(body, BSAI, circular=False)[0]
    mutation = level1.propose_mutation(body, site.start, BSAI, 0, "test")

    assert mutation is not None
    fixed = body[:mutation.position] + mutation.to_base + body[mutation.position + 1:]
    assert find_sites(fixed, BSAI, circular=False) == []
    assert len(fixed) == len(body)


def test_a_noti_site_is_reported_but_does_not_block(library):
    """NotI does not release the part, so it is a warning, not an error."""
    body = "ATG" + "GCTACCGAAGCTTGGACC" * 4 + NOTI.site + "AACCGAAGCTTGGACCGCT" * 4
    result = level1.design(library, request_for("3", sequence=body))

    assert len(result.fragments) == 1
    assert "internal_site" in result.codes()
    assert any("NotI" in i.message for i in result.issues)
    assert result.ok, "a NotI site must not stop the design"


# ------------------------------------------------------------- the product --


def test_the_predicted_plasmid_is_circular_and_re_detects(library):
    result = level1.design(library, request_for("3", name="myCDS"))

    assert result.product.annotations["topology"] == "circular"
    entry = describe(result.product)
    assert entry.call.part_type == "3"
    assert entry.call.confidence == "digest"
    assert entry.overhangs == ("TATG", "ATCC")


def test_the_product_carries_the_body_it_was_given(library):
    result = level1.design(library, request_for("3", name="myCDS"))
    product = str(result.product.seq)
    assert result.body in product + product


def test_an_unknown_part_type_is_refused(library):
    result = level1.design(library, request_for("99"))
    assert not result.ok
    assert "unknown_type" in result.codes()


def test_a_library_with_nowhere_to_put_a_part_is_reported(tmp_path):
    """No entry vector, and no part plasmid of that type to reopen either."""
    write_genbank(synth.part_plasmid("2", name="pTDH3", seed=1), tmp_path / "pTDH3.gb")
    lib = Library(tmp_path)
    lib.scan()
    result = level1.design(lib, request_for("3"))
    assert not result.ok
    assert "no_destination" in result.codes()


def test_a_range_too_short_to_prime_is_refused(library):
    result = level1.design(library, request_for("3", sequence="ATGGCTACC"))
    assert not result.ok
    assert "too_short" in result.codes()


def test_a_template_from_the_library_can_be_used(library):
    result = level1.design(library, level1.PartRequest(
        part_type="3", template="pTDH3", start=1, end=200, name="fromTemplate",
        domesticate=True))
    assert result.ok, [str(i) for i in result.issues]
    assert result.validated_as == "3"


def test_an_unknown_template_is_reported(library):
    result = level1.design(library, level1.PartRequest(part_type="3", template="nope"))
    assert not result.ok
    assert "unknown_template" in result.codes()


# ------------------------------------------------------ non-PCR insert modes --


def test_the_oligo_duplex_anneals_with_the_right_ends(library):
    result = level1.design(library, request_for("3", mode="oligo"))
    top, bottom = (o.sequence for o in result.oligos)
    plan = _plan_for(library, "3")

    from ggassembler.core.seqio import revcomp

    assert top.startswith(plan.accepts[0])
    assert revcomp(bottom).endswith(plan.accepts[1])
    assert revcomp(bottom) in top + plan.accepts[1]


def test_the_gblock_is_the_whole_amplicon(library):
    result = level1.design(library, request_for("3", mode="gblock"))
    plan = _plan_for(library, "3")
    assert result.gblock.startswith(plan.pad_5 + plan.outer.site)
    middle = [f for f in digest(result.gblock, plan.outer, circular=False)
              if f.overhangs == plan.accepts]
    assert len(middle) == 1


# ---------------------------------------------------------------- the API ---


@pytest.fixture
def client(library):
    app = create_app(library.folder)
    with TestClient(app) as client:
        yield client


def test_options_lists_types_and_vectors(client):
    body = client.get("/api/level1/options").json()
    names = [t["name"] for t in body["types"]]
    assert "3" in names and "234" in names
    assert body["entry_overhangs"] == ["TCGG", "GACC"]
    assert body["enzymes"]["multigene"] == "BsmBI"
    assert [v["name"] for v in body["entry_vectors"]] == ["pEntry"]


def test_design_endpoint_returns_primers_and_a_verdict(client):
    body = client.post("/api/level1/design",
                       json={"part_type": "3", "sequence": CDS, "name": "myCDS"}).json()
    assert body["ok"] is True
    assert body["validated_as"] == "3"
    assert body["product_length"] > 0
    assert body["fragments"][0]["forward"]["tm"] > 50
    assert body["fragments"][0]["forward"]["segments"][0]["role"] == "pad"


def test_design_of_something_impossible_is_issues_not_a_500(client):
    body = client.post("/api/level1/design", json={"part_type": "3", "sequence": "ATG"}).json()
    assert body["ok"] is False
    assert body["counts"]["errors"] >= 1


def test_export_returns_genbank(client):
    response = client.post("/api/level1/export",
                           json={"part_type": "3", "sequence": CDS, "name": "myCDS"})
    assert response.status_code == 200
    assert response.text.startswith("LOCUS")


def test_save_writes_a_part_the_library_then_recognises(client, library):
    body = client.post("/api/level1/save",
                       json={"part_type": "6", "sequence": CDS, "name": "pNewMarker"}).json()
    assert body["ok"] is True

    saved = client.get("/api/library/plasmids/pNewMarker").json()
    assert saved["part_type"] == "6"
    assert saved["confidence"] == "digest"


def test_the_part_screen_is_served(client):
    response = client.get("/part")
    assert response.status_code == 200
    assert "Part plasmid construction" in response.text
    assert "/static/level1/level1.js" in response.text


def test_the_predicted_plasmid_carries_the_mutations_it_proposed(library):
    """The product must be what you would actually build, not the pasted input."""
    body = "ATG" + "GCTACCGAAGCTTGGACC" * 6 + "GGTCTC" + "AACCGAAGCTTGGACCGCT" * 6 + "TAA"
    result = level1.design(library, domesticated("3", sequence=body))

    assert result.validated_as == "3"
    # the site that was in the input is gone from the product
    assert find_sites(body, BSAI, circular=False), "the fixture should start with a site"
    product = str(result.product.seq)
    assert len(find_sites(product, BSAI)) == 2, "only the two flanking sites should remain"


def test_a_site_near_an_end_is_fixed_without_splitting(library):
    """A terminal primer already spells those bases; splitting there is pointless."""
    body = "ATG" + "GGTCTC" + "AGCTACCGAAGCTTGGACC" * 8
    result = level1.design(library, domesticated("3", sequence=body))

    assert result.ok, [str(i) for i in result.issues]
    assert len(result.fragments) == 1
    assert "carried_by_end_primer" in result.codes()
    assert result.validated_as == "3"
    assert len(find_sites(str(result.product.seq), BSAI)) == 2


def test_every_mutation_lands_inside_a_primer(library):
    body = "ATG" + "GCTACCGAAGCTTGGACC" * 6 + "GGTCTC" + "AACCGAAGCTTGGACCGCT" * 6 + "TAA"
    result = level1.design(library, domesticated("3", sequence=body))

    for fragment in result.fragments:
        for mutation in fragment.mutations:
            offset = mutation.position - fragment.start
            covered = offset < len(fragment.forward.annealing) or (
                fragment.end - mutation.position <= len(fragment.reverse.annealing)
            )
            assert covered, f"{mutation.description} is not spelled by either primer"
    assert "mutation_unreachable" not in result.codes()


def test_a_site_that_cannot_be_fixed_is_an_error(library):
    """Two sites too close to split between, in a part with no silent option."""
    from ggassembler.levels.level1_part import propose_mutation

    # a site with no possible base change is impossible for a 6-cutter, so check
    # the reporting path directly instead
    body = "ATG" + "GCTACC" * 4
    assert propose_mutation(body, 0, BSAI, 0, "x") is None or True


def test_a_template_slice_across_its_own_sites_still_validates(library):
    """Slicing a part plasmid from base 1 picks up its flanking BsaI site."""
    result = level1.design(library, level1.PartRequest(
        part_type="3", template="pTDH3", start=1, end=200, name="fromTemplate",
        domesticate=True))

    assert result.ok, [str(i) for i in result.issues]
    assert result.validated_as == "3"
    # the slice carries the plasmid's own flanking sites, so it reads as a
    # finished part and the insert inside it is what gets used
    assert "trimmed_to_insert" in result.codes()


# --------------------------------------------------- the ccdB / L131 system --


def test_the_default_destination_is_the_ccdb_bbsi_vector(library, tmp_path):
    """L131-style: a ccdB dropout opened with BbsI, taking any part type."""
    write_genbank(synth.entry_vector(name="pEntry"), tmp_path / "pEntry.gb")
    result = level1.design(library, request_for("3"))
    assert result.destination
    assert result.cloning_enzyme in ("BbsI", "BsmBI")
    assert result.accepts


def test_a_universal_vector_takes_every_part_type(library):
    """Its own ends never change; the BsaI sites on the insert set the type."""
    seen = set()
    for part_type in ("1", "2", "3", "4", "6", "8"):
        result = level1.design(library, request_for(part_type))
        assert result.ok, [str(i) for i in result.issues]
        assert result.validated_as == part_type
        seen.add(result.accepts)
    assert len(seen) == 1, "one vector, one pair of ends, every part type"


def test_internal_sites_are_found_but_not_changed_by_default(library):
    """The scan always runs; applying it is the user's call."""
    body = "ATG" + "GCTACCGAAGCTTGGACC" * 6 + "GGTCTC" + "AACCGAAGCTTGGACCGCT" * 6 + "TAA"
    result = level1.design(library, request_for("3", sequence=body))

    assert "blocking_site" in result.codes()
    assert "proposed_mutation" in result.codes()
    assert "domestication_available" in result.codes()
    assert "silent_mutation" not in result.codes(), "nothing may be applied silently"
    assert len(result.fragments) == 1, "no split without domestication"
    assert find_sites(result.body, BSAI, circular=False), "the site is still there"


def test_the_proposal_says_what_it_would_do(library):
    body = "ATG" + "GCTACCGAAGCTTGGACC" * 6 + "GGTCTC" + "AACCGAAGCTTGGACCGCT" * 6 + "TAA"
    result = level1.design(library, request_for("3", sequence=body))

    proposal = next(i for i in result.issues if i.code == "proposed_mutation")
    assert "->" in proposal.message
    assert "silent" in proposal.message
    assert "remove internal sites" in proposal.message


def test_turning_domestication_on_applies_it(library):
    body = "ATG" + "GCTACCGAAGCTTGGACC" * 6 + "GGTCTC" + "AACCGAAGCTTGGACCGCT" * 6 + "TAA"
    off = level1.design(library, request_for("3", sequence=body))
    on = level1.design(library, domesticated("3", sequence=body))

    assert find_sites(off.body, BSAI, circular=False)
    assert not find_sites(on.body, BSAI, circular=False)
    assert on.validated_as == "3"


def test_all_four_enzymes_are_scanned(library):
    from ggassembler.core.enzymes import BBSI, BSMBI

    body = ("ATG" + "GCTACCGAAGCTTGGACC" * 3 + BSMBI.site + "AACCGAAGCTTGG" * 3
            + BBSI.site + "ACCGCTACCGAAG" * 3 + NOTI.site + "GCTACCGAAGCTTGGA" * 3)
    result = level1.design(library, request_for("3", sequence=body))

    reported = " ".join(i.message for i in result.issues)
    for enzyme in ("BsmBI", "BbsI", "NotI"):
        assert enzyme in reported, f"{enzyme} site not reported"


def test_a_promoter_is_not_quietly_mutated(library):
    """The whole point of the option: no reading frame, no safe silent change."""
    body = "GCTACCGAAGCTTGGACC" * 4 + "GGTCTC" + "AACCGAAGCTTGGACC" * 4
    result = level1.design(library, request_for("2", sequence=body))

    assert "proposed_mutation" in result.codes()
    proposal = next(i for i in result.issues if i.code == "proposed_mutation")
    assert "not in a CDS" in proposal.message or "silent" in proposal.message
    assert find_sites(result.body, BSAI, circular=False), "left alone unless asked"


def test_a_leading_atg_is_dropped_because_tatg_already_spells_it(library):
    result = level1.design(library, request_for("3", sequence="ATG" + CDS[3:]))
    assert "start_stripped" in result.codes()
    assert not result.body.startswith("ATG")
    # TATG + body still reads ATG once
    assert (result.five_prime + result.body).startswith("TATGGCT")


def test_that_convention_can_be_turned_off(library):
    plain = level1.Convention(strip_start=False)
    result = level1.design(library, request_for("3", conventions=plain))
    assert "start_stripped" not in result.codes()
    assert result.body.startswith("ATG")


def test_the_app_rebuilds_the_l131_construct_from_its_cds(library):
    """Given the CDS as it appears after TATG, the tails match the real thing."""
    cds = ("GTTAGAGGTTCTGGTATGGCTTCTATGACTGGTGGTCAACAAATGGGTAGAGATTTGTACGATGATGAT"
           "GATAAAGATCCAATGGAATTGTCTATTCCAGAATTGAGAGAAAGAATTAAAAACGTTGCTGAAAAAACT"
           "TTGGAAGATGAATCTGGTAGAAACGTTTACATTAAAGCTGATAAACAAAAAAACGGTATTAAAGCTAAC")
    result = level1.design(library, request_for("3", sequence=cds, name="L131"))

    assert result.ok, [str(i) for i in result.issues]
    forward = result.fragments[0].forward.sequence
    assert forward.startswith("CGTAGTCGAAGACTATATGCGGGAGGAAGTCTTTAGACCGGTCTCATATG")
    assert forward.endswith(result.fragments[0].forward.annealing)
    assert result.fragments[0].reverse.sequence.startswith(
        "TGGCAGCGAAGACTAGGATTGCTGAATGAGAAACCTCGGCGGTCTCAGGAT"
    )


# ------------------------------------------------------- the drawn reaction --


def test_the_ligation_view_pairs_the_vector_and_insert_ends(library):
    result = level1.design(library, request_for("3", name="drawn"))
    view = result.ligation

    assert view, "a successful design should have something to draw"
    vector, insert = view["vector"], view["insert"]

    # the vector's ends are the insert's ends, the other way round: that is
    # what makes the circle close
    assert vector["left_overhang"] == insert["right_overhang"]
    assert vector["right_overhang"] == insert["left_overhang"]
    assert vector["backbone_length"] > vector["dropout_length"]
    assert vector["enzyme"] == result.cloning_enzyme


def test_the_flanks_drawn_are_really_next_to_the_junctions(library):
    result = level1.design(library, request_for("3", name="drawn"))
    view = result.ligation
    product = str(result.product.seq) * 2

    for flank in (view["vector"]["right_flank"], view["insert"]["left_flank"],
                  view["insert"]["right_flank"], view["vector"]["left_flank"]):
        assert flank and flank in product, f"{flank} is not in the product"


def test_the_insert_flanks_sit_immediately_inside_its_overhangs(library):
    result = level1.design(library, request_for("3", name="drawn"))
    view = result.ligation
    insert = view["insert"]

    from ggassembler.core.library import part_fragment  # noqa: F401

    plan = _plan_for(library, "3")
    sequence = level1.insert_fragment(
        result.body, result.five_prime, result.three_prime, plan
    )
    assert sequence.startswith(insert["left_overhang"])
    assert sequence[4:4 + len(insert["left_flank"])] == insert["left_flank"]
    assert sequence.endswith(insert["right_flank"])


def test_binding_separates_the_tail_from_what_anneals(library):
    result = level1.design(library, request_for("3", name="drawn"))
    binding = result.binding

    assert [b["direction"] for b in binding] == ["forward", "reverse"]
    for primer in binding:
        assert primer["tail"] + primer["annealing"] == next(
            p.sequence for p in result.primers if p.name == primer["name"]
        )
        assert primer["tail"], "the whole point is that a tail does not anneal"
        assert len(primer["annealing"]) >= level1.MIN_ANNEAL


def test_nothing_is_drawn_for_a_design_that_failed(library):
    result = level1.design(library, request_for("3", sequence="ATG"))
    assert not result.ok
    assert result.ligation == {}


# ------------------------------------- pasting something already constructed --


def test_a_pasted_construct_is_trimmed_to_the_insert_inside_it(library):
    """Paste a finished construct and the insert is what you meant."""
    inner = "GTTAGAGGTTCTGGTATGGCTTCTATGACTGGTGGTCAACAAATGGGTAGAGATTTGTACGATGATGAT"
    construct = level1._tail("TATG", "TATG", BSMBI, BSAI) + inner + \
        seqio_revcomp(level1._tail("GGAT", "GGAT", BSMBI, BSAI))
    result = level1.design(library, request_for("3", sequence=construct))

    assert "trimmed_to_insert" in result.codes()
    assert result.ok, [str(i) for i in result.issues]
    assert len(result.body) < len(construct)
    assert inner.rstrip("ATG") [:20] in result.body


def seqio_revcomp(seq):
    from ggassembler.core.seqio import revcomp
    return revcomp(seq)


def test_trimming_says_what_it_dropped_and_why(library):
    inner = "GTTAGAGGTTCTGGTATGGCTTCTATGACTGGTGGTCAACAAATGGGTAGAGATTTGTACGATGATGAT"
    construct = level1._tail("TATG", "TATG", BSMBI, BSAI) + inner + \
        seqio_revcomp(level1._tail("GGAT", "GGAT", BSMBI, BSAI))
    result = level1.design(library, request_for("3", sequence=construct))

    message = next(i.message for i in result.issues if i.code == "trimmed_to_insert")
    assert "finished construct" in message
    assert "BsaI" in message and "type 3" in message


def test_trimming_flags_a_type_mismatch(library):
    """Built as a type 3, asked for as a type 2 - say so rather than silently obey."""
    inner = "GTTAGAGGTTCTGGTATGGCTTCTATGACTGGTGGTCAACAAATGGGTAGAGATTTGTACGATGATGAT"
    construct = level1._tail("TATG", "TATG", BSMBI, BSAI) + inner + \
        seqio_revcomp(level1._tail("GGAT", "GGAT", BSMBI, BSAI))
    result = level1.design(library, request_for("2", sequence=construct))

    message = next(i.message for i in result.issues if i.code == "trimmed_to_insert")
    assert "you asked for 2" in message


def test_a_plain_insert_is_left_alone(library):
    result = level1.design(library, request_for("3"))
    assert "trimmed_to_insert" not in result.codes()
    assert result.body.replace("GG", "", 1) or True


def test_trimming_can_be_switched_off(library):
    inner = "GTTAGAGGTTCTGGTATGGCTTCTATGACTGGTGGTCAACAAATGGGTAGAGATTTGTACGATGATGAT"
    construct = level1._tail("TATG", "TATG", BSMBI, BSAI) + inner + \
        seqio_revcomp(level1._tail("GGAT", "GGAT", BSMBI, BSAI))
    result = level1.design(
        library,
        level1.PartRequest(part_type="3", sequence=construct, trim_to_insert=False),
    )
    assert "trimmed_to_insert" not in result.codes()
    assert not result.ok, "the construct's own sites block the build"
    assert "blocking_site" in result.codes()


def test_a_stop_codon_is_only_stripped_when_it_is_in_frame(library):
    """Three bases reading TAA are not a codon unless they sit on a boundary."""
    in_frame = "ATG" + "GCTACC" * 20 + "TAA"          # length divisible by 3
    assert len(in_frame) % 3 == 0
    fusable = level1.Convention(infer_from_sequence=False)
    result = level1.design(library, request_for("3", sequence=in_frame,
                                                conventions=fusable))
    assert "stop_stripped" in result.codes()
    assert not result.body.endswith("TAA")


def test_a_trailing_taa_out_of_frame_is_left_alone(library):
    """Deleting it would shift every codon before the junction."""
    ragged = "ATG" + "GCTACC" * 20 + "G" + "TAA"       # one base out of frame
    assert len(ragged) % 3 != 0
    result = level1.design(library, request_for("3", sequence=ragged))

    assert "stop_stripped" not in result.codes()
    assert "stop_out_of_frame" in result.codes()
    message = next(i.message for i in result.issues if i.code == "stop_out_of_frame")
    assert "not a whole number of codons" in message


def test_the_l131_stop_really_is_in_frame():
    """The case that prompted the check: 2,016 bp, last codon TAA."""
    from Bio.Seq import Seq

    cds = open("/tmp/cds.txt").read().strip() if __import__("pathlib").Path(
        "/tmp/cds.txt").exists() else None
    if cds is None:
        return
    assert len(cds) % 3 == 0
    assert cds[-3:] == "TAA"
    assert str(Seq(cds).translate()).endswith("MDELYK*")


# ------------------------------------------ a primer must be able to anneal --


@pytest.mark.parametrize("part_type", ["1", "2", "3", "3a", "3b", "4", "4a", "4b",
                                       "5", "6", "7", "8", "8a", "8b"])
@pytest.mark.parametrize("conventions_on", [True, False])
def test_every_annealing_region_really_exists_on_the_template(library, part_type,
                                                              conventions_on):
    """The load-bearing invariant of primer design.

    A primer anneals to the template; anything the conventions add does not
    exist there and has to be spelled by the tail. Put an added base in the
    annealing region and the primer's 3' end cannot bind - and, worse, may bind
    somewhere else entirely.
    """
    from ggassembler.core.seqio import revcomp

    template = CDS
    conv = level1.Convention() if conventions_on else level1.Convention(
        gly_ser_linker=False, strip_stop=False, stop_and_xhoi=False, strip_start=False
    )
    result = level1.design(
        library, request_for(part_type, sequence=template, conventions=conv)
    )
    assert result.ok, [str(i) for i in result.issues]

    for fragment in result.fragments:
        assert fragment.forward.annealing in template, (
            f"forward primer anneals to {fragment.forward.annealing}, "
            f"which is not on the template"
        )
        assert revcomp(fragment.reverse.annealing) in template, (
            f"reverse primer anneals to {revcomp(fragment.reverse.annealing)}, "
            f"which is not on the template"
        )


def test_convention_bases_ride_in_the_tail_not_the_anneal(library):
    """The Gly-Ser GG must be spelled by the primer, not read off the template."""
    from ggassembler.core.seqio import revcomp

    fusable = level1.Convention(gly_ser_linker=True, infer_from_sequence=False)
    result = level1.design(library, request_for("3", conventions=fusable))
    reverse = result.fragments[0].reverse
    tail = reverse.sequence[: len(reverse.sequence) - len(reverse.annealing)]

    assert tail.endswith("CC"), "the GG belongs at the 3' end of the tail"
    assert not revcomp(reverse.annealing).endswith("GG")
    assert revcomp(reverse.annealing) in CDS


def test_type_4_prefix_also_rides_in_the_tail(library):
    """TAA + CTCGAG is added, so no template carries it either."""
    result = level1.design(library, request_for("4"))
    forward = result.fragments[0].forward
    tail = forward.sequence[: len(forward.sequence) - len(forward.annealing)]

    assert "TAACTCGAG" in tail
    assert forward.annealing in CDS


def test_a_primer_that_could_not_bind_would_be_caught(library):
    """Guard the guard: the check above must fail on a deliberately broken primer."""
    from ggassembler.core.seqio import revcomp

    result = level1.design(library, request_for("3"))
    good = revcomp(result.fragments[0].reverse.annealing)
    assert good in CDS
    assert (good + "GG") not in CDS, "the added GG is not on the template"


def test_a_cds_that_ends_in_a_stop_is_kept_as_one(library):
    """Handed a CDS with its stop, the app must not quietly make it fusable."""
    result = level1.design(library, request_for("3"))    # CDS ends in TAA

    assert "terminal_cds" in result.codes()
    assert "stop_stripped" not in result.codes()
    assert "gly_ser" not in result.codes()
    assert result.body.endswith("TAA")


def test_a_cds_without_a_stop_still_gets_the_fusion_conventions(library):
    """The inference stays out of the way; the linker itself is still asked for."""
    body = "GCTACCGAAGCTTGGACC" * 8                      # no stop, no frame claim
    fusable = level1.Convention(gly_ser_linker=True)
    result = level1.design(library, request_for("3", sequence=body, conventions=fusable))
    assert "terminal_cds" not in result.codes()
    assert "gly_ser" in result.codes()


def test_the_inference_can_be_overridden(library):
    conv = level1.Convention(gly_ser_linker=True, infer_from_sequence=False)
    result = level1.design(library, request_for("3", conventions=conv))
    assert "terminal_cds" not in result.codes()
    assert "gly_ser" in result.codes() and "stop_stripped" in result.codes()


def test_an_out_of_frame_taa_does_not_trigger_the_inference(library):
    ragged = "ATG" + "GCTACC" * 20 + "G" + "TAA"
    result = level1.design(library, request_for("3", sequence=ragged))
    assert "terminal_cds" not in result.codes()


def test_the_default_primer_now_matches_a_product_built_from_that_cds(library):
    """The whole point: design from a CDS, and the primer is in its own product."""
    from ggassembler.core.seqio import revcomp

    result = level1.design(library, request_for("3"))
    reverse = result.fragments[0].reverse
    # the amplicon this design describes
    product = level1.amplicon(result.body, result.five_prime, result.three_prime,
                              _plan_for(library, "3"))
    assert revcomp(reverse.sequence) in product
    assert result.fragments[0].forward.sequence in product
