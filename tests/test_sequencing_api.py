"""``/api/sequencing/*`` - the screen that checks clones against a reference."""

from __future__ import annotations

import pytest
from fastapi.testclient import TestClient

from ggassembler.api.main import create_app
from ggassembler.core import align
from ggassembler.core.seqio import write_genbank

from . import synth
from .test_assembly import CANONICAL


@pytest.fixture
def client(tmp_path):
    for part_type, name in CANONICAL.items():
        write_genbank(
            synth.part_plasmid(part_type, name=name, seed=synth.seed_for(name)),
            tmp_path / f"{name}.gb",
        )
    app = create_app(tmp_path)
    with TestClient(app) as client:
        client.library_dir = tmp_path
        yield client


def reference_of(client, name: str = "pTDH3") -> str:
    entry = client.app.state.library.get(name)
    return str(client.app.state.library.record(entry).seq).upper()


def fasta(name: str, sequence: str) -> dict:
    return {"name": f"{name}.fa", "text": f">{name}\n{sequence}\n"}


def run(client, clones, reference="pTDH3", source="library", **extra):
    body = {"source": source, "reference": reference, "clones": clones}
    body.update(extra)
    return client.post("/api/sequencing/run", json=body)


# ------------------------------------------------------------- references ---


def test_references_lists_the_library(client):
    body = client.get("/api/sequencing/references").json()
    assert {r["name"] for r in body["library"]} >= set(CANONICAL.values())
    assert body["expected"] == []


def test_an_exported_cassette_becomes_a_reference_you_can_check_against(client):
    """The point of keeping it: the design and the reads are a week apart.

    Exporting a construct is the moment its predicted sequence exists; if it is
    only ever a download, the file is missing on the day the clones come back.
    """
    export = client.post(
        "/api/level2/export", json={"selections": dict(CANONICAL), "name": "pCassette"}
    )
    assert export.status_code == 200

    body = client.get("/api/sequencing/references").json()
    assert [r["name"] for r in body["expected"]] == ["pCassette"]

    # and it works as a reference, which is the part that could silently rot
    kept = client.app.state.library.expected_record("pCassette")
    sequence = str(kept.seq).upper()
    answer = run(client, [fasta("c1", sequence)], reference="pCassette", source="expected")
    assert answer.status_code == 200
    assert answer.json()["clones"][0]["clean"]


# -------------------------------------------------------------------- run ---


def test_a_rotated_clone_comes_back_clean(client):
    sequence = reference_of(client)
    rotated = sequence[600:] + sequence[:600]
    body = run(client, [fasta("clone1", rotated)]).json()
    clone = body["clones"][0]
    assert clone["name"] == "clone1"
    assert clone["clean"] and clone["placed"]
    assert clone["verdict"] == "no mismatches"
    assert body["summary"]["clean"] == ["clone1"]
    assert "clone1 has no mismatches" == f"{clone['name']} has {clone['verdict']}"


def test_a_flipped_clone_says_so_without_calling_it_a_difference(client):
    sequence = reference_of(client)
    body = run(client, [fasta("clone2", align.revcomp(sequence))]).json()
    clone = body["clones"][0]
    assert clone["strand"] == "-"
    assert clone["clean"]


def test_a_mismatch_is_reported_with_its_position(client):
    sequence = reference_of(client)
    mutant = list(sequence)
    mutant[120] = "A" if sequence[120] != "A" else "G"
    body = run(client, [fasta("clone3", "".join(mutant))]).json()
    clone = body["clones"][0]
    assert not clone["clean"]
    assert clone["differences"][0]["kind"] == "substitution"
    assert clone["differences"][0]["start"] == 120
    assert clone["differences"][0]["where"] == "121"
    assert body["summary"]["changed"] == ["clone3"]


def test_a_clone_that_is_not_this_plasmid_is_unplaced_not_wrong(client):
    """Anything else would bury a real result under thousands of mismatches."""
    other = reference_of(client, "pTEF1") if client.app.state.library.get("pTEF1") else "ACGT" * 300
    body = run(client, [{"name": "clone4.fa", "text": f">x\n{'ACGT' * 400}\n"}]).json()
    clone = body["clones"][0]
    assert not clone["placed"]
    assert clone["differences"] == []
    assert "does not match" in clone["verdict"]
    assert body["summary"]["unplaced"] == ["clone4"]
    assert other  # the fixture library varies; this only guards the setup above


def test_several_clones_share_one_set_of_columns(client):
    """A column on screen must be the same reference base on every row."""
    sequence = reference_of(client)
    clones = [
        fasta("a", sequence),
        fasta("b", sequence[:200] + "TTTT" + sequence[200:]),
        fasta("c", sequence[:300] + sequence[306:]),
    ]
    body = run(client, clones).json()
    alignment = body["alignment"]
    assert len(alignment["reference"]) == alignment["columns"]
    assert all(len(row["row"]) == alignment["columns"] for row in alignment["rows"])
    assert [row["name"] for row in alignment["rows"]] == ["a", "b", "c"]


def test_the_headline_names_the_clones_worth_growing_up(client):
    sequence = reference_of(client)
    mutant = list(sequence)
    mutant[50] = "A" if sequence[50] != "A" else "C"
    body = run(client, [fasta("good", sequence), fasta("bad", "".join(mutant))]).json()
    assert body["summary"]["headline"] == "good has no mismatches"


def test_the_file_name_becomes_the_clone_name(client):
    sequence = reference_of(client)
    body = run(client, [{"name": "P12345_A03_consensus.gbk", "text": f">x\n{sequence}\n"}]).json()
    assert body["clones"][0]["name"] == "P12345_A03"


def test_a_file_that_is_not_sequence_is_rejected_by_name(client):
    sequence = reference_of(client)
    body = run(
        client,
        [fasta("good", sequence), {"name": "oops.txt", "text": "no such order"}],
    ).json()
    assert body["rejected"] == [{"name": "oops", "why": body["rejected"][0]["why"]}]
    assert "not DNA" in body["rejected"][0]["why"]
    assert [c["name"] for c in body["clones"]] == ["good"]


def test_an_uploaded_reference_needs_no_library_entry(client):
    sequence = reference_of(client)
    answer = client.post(
        "/api/sequencing/run",
        json={
            "source": "upload",
            "reference_file": {"name": "expected.fa", "text": f">expected\n{sequence}\n"},
            "clones": [fasta("clone1", sequence[900:] + sequence[:900])],
        },
    )
    assert answer.status_code == 200
    body = answer.json()
    assert body["reference"]["name"] == "expected"
    assert body["clones"][0]["clean"]


def test_a_reference_that_does_not_exist_is_a_404(client):
    answer = run(client, [fasta("c", "ACGT" * 100)], reference="pNope")
    assert answer.status_code == 404


def test_no_files_is_refused_rather_than_answered_emptily(client):
    assert run(client, []).status_code == 422


# ------------------------------------------------------------------ pasted ---


def test_bare_pasted_sequence_needs_no_header_or_file(client):
    """Someone copying a consensus out of a vendor's web page has neither."""
    sequence = reference_of(client)
    body = run(client, [{"name": "pasted 1", "text": sequence[300:] + sequence[:300]}]).json()
    assert body["clones"][0]["name"] == "pasted 1"
    assert body["clones"][0]["clean"]


def test_pasted_sequence_survives_the_whitespace_a_copy_brings_with_it(client):
    """Numbered, wrapped, upper and lower case - what a web page hands you."""
    sequence = reference_of(client)
    messy = "\n".join(
        f"  {i + 1:>6} " + " ".join(sequence[i:i + 60][j:j + 10].lower() for j in range(0, 60, 10))
        for i in range(0, len(sequence), 60)
    )
    body = run(client, [{"name": "A01", "text": messy}]).json()
    assert body["rejected"] == []
    assert body["clones"][0]["name"] == "A01"
    assert body["clones"][0]["clean"]


def test_a_pasted_multi_record_fasta_becomes_several_clones(client):
    """One paste off a plate can hold every colony on it."""
    sequence = reference_of(client)
    mutant = list(sequence)
    mutant[80] = "A" if sequence[80] != "A" else "C"
    text = f">A01\n{sequence}\n>A02\n{''.join(mutant)}\n"
    body = run(client, [{"name": "pasted 1", "text": text}]).json()
    assert [c["name"] for c in body["clones"]] == ["pasted 1 \u00b7 A01", "pasted 1 \u00b7 A02"]
    assert body["clones"][0]["clean"]
    assert not body["clones"][1]["clean"]


def test_a_pasted_reference_needs_no_file_either(client):
    sequence = reference_of(client)
    answer = client.post(
        "/api/sequencing/run",
        json={
            "source": "upload",
            "reference_file": {"name": "pasted 1", "text": sequence},
            "clones": [{"name": "A01", "text": sequence[120:] + sequence[:120]}],
        },
    )
    assert answer.status_code == 200
    assert answer.json()["clones"][0]["clean"]
