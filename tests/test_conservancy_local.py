"""Tests for the local Epitope Conservancy Analysis reimplementation.

Includes a round-trip integration test proving the generated CSVs are consumed correctly by the
existing POA2 conservancy parser (compatibility with the IEDB output format), and the checks that
keep the sequence-identity threshold honest end to end: it must change the result, travel with the
data, and stop a POA2 run that was told a different value.
"""
import logging

import pytest

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
    # the header spells out the threshold the row was computed at, exactly as IEDB writes it
    assert list(row.keys()) == cc.csv_columns(">=", 70)
    assert row["Epitope name"] == "SARS_SPIKE_Bepipred3.0_4_6"
    assert row["Epitope length"] == 3
    assert row["Percent of protein sequence matches at identity >= 70%"] == "50.00% (1/2)"  # only P1
    assert row["Maximum identity"] == "100.00%"
    assert row["Minimum identity"].startswith("66.6")


def test_threshold_changes_the_result():
    """The regression this suite previously could not catch: varying -t must move the numbers."""
    epitopes = [("SP_E_M_1_3", "AYI")]
    proteins = [("P1", "GGAYIGG"), ("P2", "GGAYLGG")]  # identities 100 and 66.67

    at_50 = cc.conservancy_rows(epitopes, proteins, threshold=50)[0]
    at_70 = cc.conservancy_rows(epitopes, proteins, threshold=70)[0]
    assert at_50[cc.percent_column(">=", 50)] == "100.00% (2/2)"
    assert at_70[cc.percent_column(">=", 70)] == "50.00% (1/2)"


def test_unique_objective_uses_the_less_than_operator():
    """'<' used to be ignored here: the count was hardcoded to '>=', so unique == conserved."""
    epitopes = [("SP_E_M_1_3", "AYI")]
    proteins = [("P1", "GGAYIGG"), ("P2", "GGAYLGG")]  # identities 100 and 66.67

    conserved = cc.conservancy_rows(epitopes, proteins, threshold=70, operator=">=")[0]
    unique = cc.conservancy_rows(epitopes, proteins, threshold=70, operator="<")[0]
    assert conserved[cc.percent_column(">=", 70)] == "50.00% (1/2)"   # P1 only
    assert unique[cc.percent_column("<", 70)] == "50.00% (1/2)"       # P2 only — the complement


def test_matches_threshold_rejects_unknown_operator():
    with pytest.raises(ValueError):
        cc.matches_threshold(80.0, 70.0, operator="~=")


# --------------------------------------------------------------------------- header plumbing
def test_percent_column_round_trip():
    assert conservancy_parser.percent_column(">=", 70) == (
        "Percent of protein sequence matches at identity >= 70%")
    assert conservancy_parser.parse_percent_column(
        "Percent of protein sequence matches at identity >= 70%") == (">=", 70.0)
    assert conservancy_parser.parse_percent_column(
        "Percent of protein sequence matches at identity <= 100%") == ("<=", 100.0)
    assert conservancy_parser.parse_percent_column("Minimum identity") is None


def test_find_percent_column_rejects_mixed_thresholds():
    columns = [conservancy_parser.percent_column(">=", 70),
               conservancy_parser.percent_column(">=", 90)]
    with pytest.raises(ValueError, match="different sequence identity thresholds"):
        conservancy_parser.find_percent_column(columns)


def test_operator_family_groups_iedb_and_poa_spellings():
    assert conservancy_parser.operator_family(">=") == conservancy_parser.operator_family(">") == "ge"
    assert conservancy_parser.operator_family("<") == conservancy_parser.operator_family("<=") == "lt"
    assert conservancy_parser.operator_family("?") == ""


# --------------------------------------------------------------------------- sidecar metadata
def test_metadata_round_trip(tmp_path):
    conservancy_parser.write_metadata(str(tmp_path), 70, ">=", source="local",
                                      species={"DENV1": {"n_proteins": 12, "discriminating": True}})
    meta = conservancy_parser.read_metadata(str(tmp_path))
    assert meta.threshold == 70.0
    assert meta.operator == ">="
    assert meta.family == "ge"
    assert meta.species["DENV1"]["n_proteins"] == 12
    assert meta.generated_at


def test_read_metadata_ignores_corrupt_file(tmp_path):
    (tmp_path / conservancy_parser.METADATA_FILENAME).write_text("{not json")
    assert conservancy_parser.read_metadata(str(tmp_path)) is None


# --------------------------------------------------------------------------- integration
def _write_poa1_epitopes(tmp_path, name="SARS_epitopes.fasta", body=">SARS_SPIKE_Bepipred3.0_4_6\nAYI\n"):
    poa1_dir = tmp_path / "poa1" / "Conservancy Analysis"
    poa1_dir.mkdir(parents=True, exist_ok=True)
    (poa1_dir / name).write_text(body)
    return poa1_dir


def test_roundtrip_with_poa2_parser(tmp_path):
    poa1_dir = _write_poa1_epitopes(tmp_path)
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


def test_run_for_dir_records_metadata(tmp_path):
    poa1_dir = _write_poa1_epitopes(tmp_path)
    out_dir = tmp_path / "csvs"
    cc.run_conservancy_for_dir(str(poa1_dir), [("P1", "GGAYIGG"), ("P2", "GGAYLGG")],
                               threshold=70, out_dir=str(out_dir), operator=">=")

    meta = conservancy_parser.read_metadata(str(out_dir))
    assert (meta.threshold, meta.operator, meta.source) == (70.0, ">=", "local")
    assert meta.species["SARS"]["n_proteins"] == 2
    assert meta.species["SARS"]["n_epitopes"] == 1
    assert meta.species["SARS"]["discriminating"] is True


def test_single_protein_comparison_set_is_flagged(tmp_path, caplog):
    """One protein cannot produce a conservancy — the epitope came from it, so it always matches."""
    poa1_dir = _write_poa1_epitopes(tmp_path)
    out_dir = tmp_path / "csvs"
    with caplog.at_level(logging.WARNING):
        cc.run_conservancy_for_dir(str(poa1_dir), [("P1", "GGAYIGG")],
                                   threshold=70, out_dir=str(out_dir))

    assert "cannot discriminate" in caplog.text
    meta = conservancy_parser.read_metadata(str(out_dir))
    assert meta.species["SARS"]["discriminating"] is False


def test_poa2_refuses_a_threshold_the_csvs_were_not_built_with(tmp_path):
    poa1_dir = _write_poa1_epitopes(tmp_path)
    out_dir = tmp_path / "csvs"
    cc.run_conservancy_for_dir(str(poa1_dir), [("P1", "GGAYIGG"), ("P2", "GGAYLGG")],
                               threshold=70, out_dir=str(out_dir), operator=">=")

    with pytest.raises(ValueError, match="threshold"):
        conservancy_parser.EpitConservAnalysis(90, ">=", 40, 100, 60, str(out_dir))

    with pytest.raises(ValueError, match="operator"):
        conservancy_parser.EpitConservAnalysis(70, "<", 40, 100, 60, str(out_dir))

    # the values the CSVs were actually produced with still go through
    assert len(conservancy_parser.EpitConservAnalysis(70, ">=", 40, 100, 60, str(out_dir))) == 1


def test_non_default_header_is_readable(tmp_path):
    """A real IEDB CSV at any threshold but the default used to crash the parser."""
    csv = (
        "Epitope #,Epitope name,Epitope sequence,Epitope length,"
        "Percent of protein sequence matches at identity >= 70%,Minimum identity,Maximum identity,"
        "View details\n"
        "1,SP_E_M_1_3,AYI,3,100.00% (2/2),80.00%,100.00%,details\n"
    )
    d = tmp_path / "csvs"
    d.mkdir()
    (d / "sp.csv").write_text(csv)

    df = conservancy_parser.EpitConservAnalysis(70, ">=", 60, 100, 60, str(d))
    assert len(df) == 1
    assert df.iloc[0]["Minimum identity(%)"] == 80.0


def test_header_mismatch_warns_without_metadata(tmp_path, caplog):
    """Legacy/manual CSVs must stay usable — the header disagreement is a warning, not a stop."""
    csv = (
        "Epitope #,Epitope name,Epitope sequence,Epitope length,"
        "Percent of protein sequence matches at identity <= 100%,Minimum identity,Maximum identity,"
        "View details\n"
        "1,SP_E_M_1_3,AYI,3,100.00% (2/2),80.00%,100.00%,details\n"
    )
    d = tmp_path / "csvs"
    d.mkdir()
    (d / "sp.csv").write_text(csv)

    with caplog.at_level(logging.WARNING):
        df = conservancy_parser.EpitConservAnalysis(70, ">=", 60, 100, 60, str(d))
    assert len(df) == 1
    assert "does not match" in caplog.text

    # ...and -strict turns the same disagreement into a hard stop
    with pytest.raises(ValueError, match="does not match"):
        conservancy_parser.EpitConservAnalysis(70, ">=", 60, 100, 60, str(d), strict=True)


def test_malformed_identity_cells_do_not_crash(tmp_path):
    """'-' / blank identity cells used to raise AttributeError inside the cleaning step."""
    csv = (
        "Epitope #,Epitope name,Epitope sequence,Epitope length,"
        "Percent of protein sequence matches at identity >= 70%,Minimum identity,Maximum identity,"
        "View details\n"
        "1,SP_E_M_1_3,AYI,3,100.00% (2/2),80.00%,100.00%,d\n"
        "2,SP_E_M_4_6,LMN,3,100.00% (2/2),-,,d\n"
    )
    d = tmp_path / "csvs"
    d.mkdir()
    (d / "sp.csv").write_text(csv)

    df = conservancy_parser.EpitConservAnalysis(70, ">=", 60, 100, 60, str(d))
    assert set(df["Epitope name"]) == {"SP_E_M_1_3"}   # the unusable row drops out, no exception


# --------------------------------------------------------------------------- identity filter
def _two_epitope_dir(tmp_path):
    csv = (
        "Epitope #,Epitope name,Epitope sequence,Epitope length,"
        "Percent of protein sequence matches at identity >= 70%,Minimum identity,Maximum identity,"
        "View details\n"
        "1,SP_E_M_1_3,AYI,3,100.00% (2/2),90.00%,100.00%,d\n"    # conserved everywhere
        "2,SP_E_M_4_6,LMN,3,100.00% (2/2),60.00%,65.00%,d\n"     # below the threshold everywhere
    )
    d = tmp_path / "csvs"
    d.mkdir()
    (d / "sp.csv").write_text(csv)
    return d


def test_identity_filter_conserved_keeps_only_epitopes_above_the_threshold(tmp_path):
    d = _two_epitope_dir(tmp_path)
    both = conservancy_parser.EpitConservAnalysis(70, ">=", 60, 100, 0, str(d))
    assert set(both["Epitope name"]) == {"SP_E_M_1_3", "SP_E_M_4_6"}

    filtered = conservancy_parser.EpitConservAnalysis(70, ">=", 60, 100, 0, str(d),
                                                      identity_filter=True)
    assert set(filtered["Epitope name"]) == {"SP_E_M_1_3"}   # minimum identity 90 >= 70


def test_identity_filter_unique_keeps_only_epitopes_below_the_threshold(tmp_path, caplog):
    csv = (
        "Epitope #,Epitope name,Epitope sequence,Epitope length,"
        "Percent of protein sequence matches at identity < 70%,Minimum identity,Maximum identity,"
        "View details\n"
        "1,SP_E_M_1_3,AYI,3,100.00% (2/2),90.00%,100.00%,d\n"
        "2,SP_E_M_4_6,LMN,3,100.00% (2/2),60.00%,65.00%,d\n"
    )
    d = tmp_path / "csvs"
    d.mkdir()
    (d / "sp.csv").write_text(csv)

    filtered = conservancy_parser.EpitConservAnalysis(70, "<", 60, 100, 0, str(d),
                                                      identity_filter=True)
    assert set(filtered["Epitope name"]) == {"SP_E_M_4_6"}   # maximum identity 65 < 70
