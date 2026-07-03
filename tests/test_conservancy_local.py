"""Tests for the local Epitope Conservancy Analysis reimplementation.

Includes a round-trip integration test proving the generated CSVs are consumed correctly by the
existing POA2 conservancy parser (compatibility with the IEDB output format).
"""
from poa.core.parsers import conservancy as conservancy_parser
from poa.services import conservancy_client as cc


def test_max_identity_exact_and_mismatch():
    assert cc.max_identity("AYI", "GGAYIGG") == 100.0          # exact window present
    assert round(cc.max_identity("AYI", "GGAYLGG"), 2) == round(2 / 3 * 100, 2)  # one mismatch
    assert cc.max_identity("AYI", "GG") == 0.0                 # protein shorter than epitope


def test_conservancy_rows_structure_and_values():
    epitopes = [("SARS_SPIKE_Bepipred3.0_4_6", "AYI")]
    proteins = [("P1", "GGAYIGG"), ("P2", "GGAYLGG")]  # identities 100 and 66.67
    rows = cc.conservancy_rows(epitopes, proteins, threshold=70)
    assert len(rows) == 1
    row = rows[0]
    assert list(row.keys()) == cc.CSV_COLUMNS
    assert row["Epitope name"] == "SARS_SPIKE_Bepipred3.0_4_6"
    assert row["Epitope length"] == 3
    assert row["Percent of protein sequence matches at identity <= 100%"] == "50.00% (1/2)"  # only P1 >= 70
    assert row["Maximum identity"] == "100.00%"
    assert row["Minimum identity"].startswith("66.6")


def test_roundtrip_with_poa2_parser(tmp_path):
    # simulate a POA1 per-species epitope FASTA
    poa1_dir = tmp_path / "poa1" / "Conservancy Analysis"
    poa1_dir.mkdir(parents=True)
    (poa1_dir / "SARS_epitopes.fasta").write_text(">SARS_SPIKE_Bepipred3.0_4_6\nAYI\n")

    proteins = [("P1", "GGAYIGG"), ("P2", "GGAYLGG")]
    out_dir = tmp_path / "csvs"

    written = cc.run_conservancy_for_dir(str(poa1_dir), proteins, threshold=70, out_dir=str(out_dir))
    assert len(written) == 1

    # POA2's parser must accept the generated CSV. Percent match is 50% (1/2).
    kept = conservancy_parser.EpitConservAnalysis(70, ">=", 40, 100, 60, str(out_dir))
    assert len(kept) == 1
    assert kept.iloc[0]["Epitope name"] == "SARS_SPIKE_Bepipred3.0_4_6"
    assert "Protein(s) sequence match(es) at Sequence identity threshold >=70" in kept.columns

    # With a stricter seq-match requirement (60%), the 50% row is filtered out.
    empty = conservancy_parser.EpitConservAnalysis(70, ">=", 60, 100, 60, str(out_dir))
    assert len(empty) == 0
