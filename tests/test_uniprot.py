"""Finding a protein and writing it back as DNA.

Nothing here touches the network. The UniProt client takes its fetcher as an
argument precisely so these can run on a plane and can't go red because a
curator renamed an entry - what they check is the parsing and the codon
choice, which is the part this app is responsible for.

The one thing they cannot check is whether the live service still answers in
this shape. That is what `UniProtError` and the 503 are for.
"""

from __future__ import annotations

import re

import pytest

from ggassembler.core import uniprot
from ggassembler.core.codons import (
    SCEREVISIAE,
    get_table,
    CODONS,
    BackTranslationError,
    back_translate,
    clean_protein,
    translate,
)
from ggassembler.core.enzymes import BBSI, BSAI, BSMBI, NOTI

AVOID = (BSAI, BSMBI, BBSI, NOTI)

#: One search result, shaped the way UniProt shapes them.
RAW = {
    "primaryAccession": "P42826",
    "uniProtkbId": "XKS1_YEAST",
    "entryType": "UniProtKB reviewed (Swiss-Prot)",
    "proteinDescription": {"recommendedName": {"fullName": {"value": "Xylulose kinase"}}},
    "genes": [{"geneName": {"value": "XKS1"}}],
    "organism": {"scientificName": "Saccharomyces cerevisiae"},
    "sequence": {"length": 600, "value": "MLCSVIQRQTREVSNTMSLDSYYLGFDLSTQQLKCLAINQDLKIVHSETVEF"},
}


def fake(results):
    """A fetcher that answers with `results` and never opens a socket."""
    def fetch(url: str):
        return {"results": results} if "/search?" in url else results[0]
    return fetch


# --------------------------------------------------------------------------- #
# reading what UniProt sends
# --------------------------------------------------------------------------- #


def test_a_result_is_reduced_to_what_a_part_designer_needs():
    found = uniprot.parse_entry(RAW)
    assert found.accession == "P42826"
    assert found.protein == "Xylulose kinase"
    assert found.genes == ["XKS1"]
    assert found.reviewed
    assert found.label == "XKS1 - Xylulose kinase"


def test_a_sparse_entry_does_not_break_the_parse():
    """UniProt omits fields rather than nulling them, and an unreviewed entry
    routinely has no recommended name."""
    found = uniprot.parse_entry({"primaryAccession": "A0A000"})
    assert found.accession == "A0A000"
    assert found.protein == "A0A000", "an entry with no name still needs one"
    assert found.genes == [] and not found.reviewed and found.length == 0


def test_a_submitted_name_is_used_when_there_is_no_recommended_one():
    found = uniprot.parse_entry({
        "primaryAccession": "X1",
        "proteinDescription": {"submissionNames": [{"fullName": {"value": "Putative kinase"}}]},
    })
    assert found.protein == "Putative kinase"


def test_reviewed_entries_come_first():
    draft = {**RAW, "primaryAccession": "Q1", "entryType": "UniProtKB unreviewed (TrEMBL)"}
    found = uniprot.search("kinase", fetch=fake([draft, RAW]))
    assert [e.accession for e in found] == ["P42826", "Q1"]


def test_an_organism_narrows_the_query_rather_than_the_results():
    seen = {}

    def fetch(url: str):
        seen["url"] = url
        return {"results": []}

    uniprot.search("kinase", organism="559292", fetch=fetch)
    assert "organism_id%3A559292" in seen["url"] or "organism_id:559292" in seen["url"]


def test_an_empty_query_asks_nothing():
    def fetch(url: str):
        raise AssertionError("an empty query should not reach the network")

    assert uniprot.search("   ", fetch=fetch) == []


def test_an_entry_with_no_sequence_is_an_error():
    with pytest.raises(uniprot.UniProtError):
        uniprot.entry("P42826", fetch=lambda url: {"primaryAccession": "P42826"})


# --------------------------------------------------------------------------- #
# writing it back as DNA
# --------------------------------------------------------------------------- #


def sites_in(dna: str) -> int:
    total = 0
    for enzyme in AVOID:
        total += len(enzyme.pattern("+").findall(dna))
        if not enzyme.is_palindromic:
            total += len(enzyme.pattern("-").findall(dna))
    return total


def test_the_dna_spells_the_protein_exactly():
    protein = "".join(k for k in CODONS if k != "*")
    dna, _ = back_translate(protein, avoid=AVOID)
    assert translate(dna) == protein


def test_no_site_is_ever_spelled():
    """Chosen codon by codon, so a site cannot appear - unlike domesticating
    afterwards, where a synonymous change can spell another one downstream."""
    import random

    random.seed(11)
    residues = "".join(k for k in CODONS if k != "*")
    protein = "".join(random.choice(residues) for _ in range(4000))
    dna, _ = back_translate(protein, avoid=AVOID, stop="TAA")
    assert sites_in(dna) == 0
    assert translate(dna).rstrip("*") == protein


def test_a_site_spanning_two_codons_is_caught():
    """The only sites that matter are the ones a codon join creates, and a
    cut-position finder will not see them: BbsI cuts twelve bases away, so in
    a window this short it reports nothing at all."""
    protein = "EDRRGGSLEDLK" * 40
    dna, _ = back_translate(protein, avoid=(BBSI,))
    assert not BBSI.pattern("+").search(dna)
    assert not BBSI.pattern("-").search(dna)


def test_the_preferred_codon_is_used_wherever_it_can_be():
    dna, compromised = back_translate("MGGGGGGK", avoid=AVOID)
    assert compromised == [], "nothing here forces a second-choice codon"
    assert dna.startswith("ATG")
    assert "GGT" in dna, "the preferred glycine codon should be used"


def test_the_positions_that_had_to_compromise_are_reported():
    """They are the only places the result differs from a plain optimised gene."""
    import random

    random.seed(3)
    residues = "".join(k for k in CODONS if k != "*")
    protein = "".join(random.choice(residues) for _ in range(2000))
    _, compromised = back_translate(protein, avoid=AVOID)
    assert compromised == sorted(set(compromised))
    assert all(0 <= i < len(protein) for i in compromised)


def test_backtracking_recovers_when_a_residue_has_no_legal_codon():
    """The codons already written decide which are legal next, so a residue can
    arrive with every one of its codons spelling a site. The fix is not there
    but one position back, which a forward-only pass cannot reach."""
    import random

    random.seed(7)
    residues = "".join(k for k in CODONS if k != "*")
    # long enough that a greedy pass provably gets stuck: it did, at residue 183
    protein = "".join(random.choice(residues) for _ in range(3000))
    dna, _ = back_translate(protein, avoid=AVOID)
    assert translate(dna) == protein
    assert sites_in(dna) == 0


def test_a_stop_is_added_only_when_the_protein_lacks_one():
    with_stop, _ = back_translate("MGK*", avoid=AVOID, stop="TAA")
    without, _ = back_translate("MGK", avoid=AVOID, stop="TAA")
    assert len(with_stop) == len(without) == 12


def test_an_ambiguous_residue_is_refused_rather_than_guessed():
    with pytest.raises(BackTranslationError, match="X"):
        back_translate("MGXK", avoid=AVOID)


def test_a_fasta_header_and_its_whitespace_are_stripped():
    assert clean_protein(">sp|P42826|XKS1_YEAST\nMLCS VIQR\nQTRE\n") == "MLCSVIQRQTRE"


def test_an_empty_protein_is_an_error():
    with pytest.raises(BackTranslationError):
        back_translate("", avoid=AVOID)


# --------------------------------------------------------------------------- #
# the endpoints
# --------------------------------------------------------------------------- #


@pytest.fixture
def client(tmp_path):
    from fastapi.testclient import TestClient

    from ggassembler.api.main import create_app
    from ggassembler.core.seqio import write_genbank

    from . import synth

    write_genbank(synth.part_plasmid("2", name="promoter"), tmp_path / "promoter.gb")
    with TestClient(create_app(tmp_path)) as client:
        yield client


def test_a_pasted_protein_needs_no_network(client):
    """The screen is useful offline too: paste a sequence, get a part."""
    body = {"protein": "MGKALEDLRQ", "name": "test"}
    made = client.post("/api/uniprot/part", json=body).json()
    assert made["verified"], "the DNA must spell the protein back"
    assert made["length"] == 33, "ten residues plus a stop"
    assert made["stop_added"]
    assert sites_in(made["dna"]) == 0
    assert set(made["avoided"]) >= {"BsaI", "BsmBI", "NotI"}


def test_the_round_trip_is_checked_on_the_server(client):
    """A part that does not spell its protein is not a part, however clean its
    sites are - so the check is made rather than assumed."""
    made = client.post("/api/uniprot/part", json={"protein": "MGK"}).json()
    assert "verified" in made and made["verified"] is True


def test_asking_for_nothing_is_refused(client):
    assert client.post("/api/uniprot/part", json={}).status_code == 422


def test_an_unwritable_protein_is_a_422_not_a_crash(client):
    response = client.post("/api/uniprot/part", json={"protein": "MGXK"})
    assert response.status_code == 422
    assert "X" in response.json()["detail"]


def test_being_offline_is_reported_as_such(client, monkeypatch):
    """Ordinary for a tool that otherwise never asks anything of the network."""
    def refuse(*args, **kwargs):
        raise uniprot.UniProtError("could not reach UniProt: offline")

    monkeypatch.setattr(uniprot, "search", refuse)
    response = client.get("/api/uniprot/search", params={"q": "kinase"})
    assert response.status_code == 503
    assert "UniProt" in response.json()["detail"]


def test_the_screen_is_served_and_linked(client):
    assert client.get("/proteins").status_code == 200
    for page in ("/", "/part", "/multigene", "/library"):
        assert 'href="/proteins"' in client.get(page).text, f"{page} has no link to it"


# --------------------------------------------------------------------------- #
# the two choices that change the sequence
# --------------------------------------------------------------------------- #


def test_a_codon_table_can_be_chosen(client):
    """A part bound for E. coli wants E. coli codons."""
    from ggassembler.core.codons import ECOLI, SCEREVISIAE, get_table

    assert get_table(None) is SCEREVISIAE
    assert get_table("ecoli") is ECOLI
    assert get_table("S. cerevisiae") is SCEREVISIAE, "the display name should resolve"
    with pytest.raises(ValueError, match="unknown codon table"):
        get_table("nonesuch")


def test_the_two_tables_give_genuinely_different_sequences():
    protein = "MGKALEDLRQAGGSRPWLEEK" * 6
    yeast, _ = back_translate(protein, avoid=AVOID, table=get_table("scerevisiae"))
    coli, _ = back_translate(protein, avoid=AVOID, table=get_table("ecoli"))
    assert yeast != coli
    assert translate(yeast) == translate(coli) == protein

    def gc(dna):
        return sum(dna.count(b) for b in "GC") / len(dna)

    # E. coli's preferred codons are markedly GC-richer; if this ever stops
    # holding, one of the tables has been mangled
    assert gc(coli) > gc(yeast) + 0.10


def test_every_table_covers_every_residue():
    from ggassembler.core.codons import TABLES

    residues = set(SCEREVISIAE)
    for name, table in TABLES.items():
        assert set(table) == residues, f"{name} is missing residues"
        for residue, codons in table.items():
            assert codons, f"{name} has no codon for {residue}"
            for codon in codons:
                assert len(codon) == 3 and set(codon) <= set("ACGT")
                assert translate(codon) == residue, f"{name}: {codon} is not {residue}"


def test_the_enzymes_to_avoid_can_be_chosen(client):
    """Unticking one is a choice - the part may be going somewhere that never
    sees that enzyme, and the alternative is a run of rare codons."""
    protein = "MGKALEDLRQ" * 8
    default = client.post("/api/uniprot/part", json={"protein": protein}).json()
    assert set(default["avoided"]) >= {"BsaI", "BsmBI", "NotI"}

    one = client.post(
        "/api/uniprot/part", json={"protein": protein, "avoid": ["BsaI"]}
    ).json()
    assert one["avoided"] == ["BsaI"]


def test_avoiding_nothing_is_a_choice_not_an_omission(client):
    """`None` means the scheme's own enzymes; `[]` means avoid nothing. They
    have to be different or one of them is unreachable."""
    protein = "MGKALEDLRQ" * 8
    none_given = client.post("/api/uniprot/part", json={"protein": protein}).json()
    explicit = client.post(
        "/api/uniprot/part", json={"protein": protein, "avoid": []}
    ).json()
    assert none_given["avoided"] and explicit["avoided"] == []


def test_an_unknown_enzyme_is_refused(client):
    response = client.post(
        "/api/uniprot/part", json={"protein": "MGK", "avoid": ["EcoRI"]}
    )
    assert response.status_code == 422


def test_an_unknown_codon_table_is_refused(client):
    response = client.post(
        "/api/uniprot/part", json={"protein": "MGK", "codon_table": "nonesuch"}
    )
    assert response.status_code == 422


def test_the_options_endpoint_lists_what_can_be_picked(client):
    options = client.get("/api/uniprot/options").json()
    assert {t["key"] for t in options["codon_tables"]} == {"scerevisiae", "ecoli"}
    assert [e["name"] for e in options["enzymes"]], "no default enzymes offered"
    for enzyme in options["enzymes"]:
        assert enzyme["default"] and enzyme["site"]


def test_the_page_offers_both_controls(client):
    html = client.get("/proteins").text
    assert 'id="codon-table"' in html and 'id="enzymes"' in html
    script = client.get("/static/uniprot/uniprot.js").text
    assert "codon_table" in script and "avoid: chosenEnzymes()" in script


# ----------------------------------------------------------------- ranking ---


def _raw(accession, gene, length, reviewed=True, name="", protein="a protein"):
    return {
        "primaryAccession": accession,
        "uniProtkbId": name or f"{gene}_YEAST",
        "genes": [{"geneName": {"value": gene}}] if gene else [],
        "proteinDescription": {"recommendedName": {"fullName": {"value": protein}}},
        "organism": {"scientificName": "Saccharomyces cerevisiae"},
        "sequence": {"length": length, "value": "M" * length},
        "entryType": "UniProtKB reviewed (Swiss-Prot)" if reviewed else "UniProtKB unreviewed (TrEMBL)",
    }


def _search(query, results, **kw):
    """Search against a canned answer, so ranking is tested and not the network."""
    return uniprot.search(query, fetch=lambda url: {"results": results}, **kw)


def test_the_gene_you_named_comes_first():
    """UniProt put GRE3 above XKS1 for the query `XKS1`.

    Defensible as text search and useless as an answer: someone typing a gene
    name has already decided which gene they want.
    """
    found = _search("XKS1", [
        _raw("P38715", "GRE3", 327, protein="NADPH-dependent aldose reductase"),
        _raw("P42826", "XKS1", 600, protein="Xylulose kinase"),
    ])
    assert [e.genes[0] for e in found] == ["XKS1", "GRE3"]


def test_an_accession_matches_exactly_too():
    found = _search("P38715", [
        _raw("P42826", "XKS1", 600),
        _raw("P38715", "GRE3", 327),
    ])
    assert found[0].accession == "P38715"


def test_the_match_is_not_case_sensitive():
    found = _search("xks1", [_raw("P38715", "GRE3", 327), _raw("P42826", "XKS1", 600)])
    assert found[0].genes == ["XKS1"]


def test_a_prefix_match_beats_an_unrelated_one_but_loses_to_an_exact_one():
    found = _search("XKS", [
        _raw("Q00000", "OTHER", 100),
        _raw("Q00001", "XKS1B", 900),
        _raw("Q00002", "XKS", 400),
    ])
    assert [e.genes[0] for e in found] == ["XKS", "XKS1B", "OTHER"]


def test_reviewed_still_wins_among_equally_good_matches():
    found = _search("XKS1", [
        _raw("A00000", "XKS1", 600, reviewed=False),
        _raw("P42826", "XKS1", 600, reviewed=True),
    ])
    assert found[0].accession == "P42826"


def test_an_exact_match_outranks_a_reviewed_entry_that_is_not_one():
    """The order that was wrong: curation is a tie-break, not the first key."""
    found = _search("XKS1", [
        _raw("P38715", "GRE3", 327, reviewed=True),
        _raw("A00000", "XKS1", 600, reviewed=False),
    ])
    assert found[0].genes == ["XKS1"]


def test_length_breaks_the_last_tie():
    """Among near-duplicates the shorter record is usually the canonical one."""
    found = _search("XKS1", [
        _raw("A00001", "XKS1", 900),
        _raw("A00002", "XKS1", 600),
    ])
    assert [e.length for e in found] == [600, 900]


def test_an_entry_with_no_gene_name_still_ranks():
    """Unreviewed records often have none; it must not raise."""
    found = _search("XKS1", [_raw("A00003", "", 200), _raw("P42826", "XKS1", 600)])
    assert found[0].genes == ["XKS1"]


# ------------------------------------------------------------- part types ---


PROTEIN = "MGKALEDLRQFATVSNW"


def make(client, **body):
    answer = client.post("/api/uniprot/part", json={"protein": PROTEIN, **body})
    assert answer.status_code == 200, answer.json()
    return answer.json()


def test_a_type_3_target_gets_no_stop_codon(client):
    """The contradiction this fixes.

    A Type 3 part omits its stop - the Type 4 terminator supplies TAA right
    after the ATCC overhang. Handing Level 1 a sequence with one on it means
    Level 1's own conventions strip it straight back off, and the two screens
    disagree about the same part.
    """
    made = make(client, part_type="3")
    assert made["stop"] is False
    assert not made["stop_added"]
    assert made["length"] == len(PROTEIN) * 3
    assert made["verified"]
    assert "Type 4 terminator" in made["stop_reason"]


@pytest.mark.parametrize("part_type", ["3", "3a", "3b"])
def test_every_coding_type_omits_the_stop(client, part_type):
    assert make(client, part_type=part_type)["stop"] is False


def test_a_4a_fusion_keeps_its_stop(client):
    """4a ends the protein itself, rather than reading through its TGGC flank."""
    made = make(client, part_type="4a")
    assert made["stop"] is True
    assert made["dna"].endswith("TAA")
    assert made["length"] == (len(PROTEIN) + 1) * 3


def test_a_plain_gene_to_order_keeps_its_stop(client):
    made = make(client, part_type="")
    assert made["stop"] is True
    assert made["dna"].endswith("TAA")
    assert "open reading frame" in made["stop_reason"]


def test_the_convention_can_be_overridden_but_says_so(client):
    """Following the type is the default, not a rule with no way out."""
    made = make(client, part_type="3", stop=True)
    assert made["stop"] is True
    assert made["stop_followed_convention"] is False
    assert make(client, part_type="3")["stop_followed_convention"] is True


def test_a_protein_that_arrives_with_a_stop_loses_it_for_a_type_3(client):
    """Otherwise the answer depends on what the source record happened to hold."""
    made = client.post(
        "/api/uniprot/part", json={"protein": PROTEIN + "*", "part_type": "3"}
    ).json()
    assert made["stop"] is False
    assert made["stop_removed"] is True
    assert made["length"] == len(PROTEIN) * 3
    assert made["verified"]


def test_an_unknown_part_type_is_refused_by_name(client):
    answer = client.post("/api/uniprot/part", json={"protein": PROTEIN, "part_type": "9"})
    assert answer.status_code == 422
    assert "9" in answer.json()["detail"]


def test_the_options_endpoint_offers_the_part_types(client):
    types = client.get("/api/uniprot/options").json()["part_types"]
    keys = {t["key"]: t for t in types}
    assert keys["3"]["stop"] is False and keys["4a"]["stop"] is True
    assert "" in keys, "there must be a way to ask for a plain gene"
    assert all(t["why"] for t in types), "every type has to say why"


# ------------------------------------------------------ sites to keep out ---


def test_a_site_can_be_given_by_name_including_the_bglbrick_pair(client):
    """The kit writes XhoI and BamHI itself, so a CDS carrying one is a problem."""
    made = make(client, avoid=["XhoI", "BamHI"])
    assert made["avoided"] == ["XhoI", "BamHI"]
    assert "CTCGAG" not in made["dna"]
    assert "GGATCC" not in made["dna"]


def test_a_site_can_be_spelled_out_instead_of_named(client):
    """The only way to ask for something outside the list without shipping REBASE."""
    made = make(client, avoid=["GAATTC"])
    assert made["avoided"] == ["GAATTC"]
    assert "GAATTC" not in made["dna"]


def test_an_iupac_site_is_accepted_and_honoured(client):
    made = make(client, avoid=["GGNCC"])
    assert made["avoided"] == ["GGNCC"]
    assert not re.search("GG[ACGT]CC", made["dna"])


def test_names_and_sites_can_be_mixed(client):
    made = make(client, avoid=["BsaI", "CTCGAG"])
    assert made["avoided"] == ["BsaI", "CTCGAG"]


def test_something_that_is_neither_is_refused_by_name(client):
    answer = client.post(
        "/api/uniprot/part", json={"protein": PROTEIN, "avoid": ["NotAnEnzyme"]}
    )
    assert answer.status_code == 422
    assert "NotAnEnzyme" in answer.json()["detail"]


def test_a_site_too_short_to_mean_anything_is_refused(client):
    """`AT` would forbid most of the genetic code and look like a hang."""
    answer = client.post("/api/uniprot/part", json={"protein": PROTEIN, "avoid": ["AT"]})
    assert answer.status_code == 422
    assert "short" in answer.json()["detail"]


# ------------------------------------------------- synthesis feasibility ---


def test_the_part_payload_carries_the_feasibility_checks(client):
    report = make(client, part_type="3")["feasibility"]
    assert set(report) == {"homopolymer", "repeat", "gc"}
    for check in report.values():
        assert "ok" in check and "limit" in check or "window" in check


def test_the_checks_carry_their_own_proof(client):
    """A row that cannot be checked is an assertion, not a report."""
    report = make(client, part_type="3")["feasibility"]
    assert report["homopolymer"]["at"] >= 0
    assert report["gc"]["window"] == 50
    assert report["gc"]["low"] == 25.0 and report["gc"]["high"] == 75.0


def test_a_repeated_protein_motif_shows_up_as_a_dna_repeat(client):
    """Highest-preference choice maps a repeated motif straight onto the bases."""
    made = client.post("/api/uniprot/part", json={
        "protein": "M" + "ACDEFGHIKLMNPQRST" * 4, "part_type": "3",
    }).json()
    assert made["feasibility"]["repeat"]["length"] >= 20
    assert made["feasibility"]["repeat"]["ok"] is False


def test_a_poly_lysine_gene_in_the_ecoli_table_is_not_ninety_adenines(client):
    """The tie-break, end to end through the API."""
    made = client.post("/api/uniprot/part", json={
        "protein": "M" + "K" * 30 + "ACDEFGHIKLMNPQRSTVWY" * 4,
        "codon_table": "ecoli",
    }).json()
    assert made["verified"]
    assert made["feasibility"]["homopolymer"]["length"] < 20
