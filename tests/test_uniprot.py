"""Finding a protein and writing it back as DNA.

Nothing here touches the network. The UniProt client takes its fetcher as an
argument precisely so these can run on a plane and can't go red because a
curator renamed an entry - what they check is the parsing and the codon
choice, which is the part this app is responsible for.

The one thing they cannot check is whether the live service still answers in
this shape. That is what `UniProtError` and the 503 are for.
"""

from __future__ import annotations

import pytest

from ggassembler.core import uniprot
from ggassembler.core.codons import (
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
