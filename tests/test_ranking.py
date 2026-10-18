"""Tests for the ranking/filtering logic (length limits, IC50/HLA, conservancy thresholds)."""
import pandas as pd

from poa.core.parsers import bepipred, conservancy, mhcii, others


def _bepipred_df():
    return pd.DataFrame({
        "Specie": ["SARS", "SARS", "SARS"],
        "Protein": ["SPIKE", "SPIKE", "SPIKE"],
        "ID_Sequence": ["NP1", "NP1", "NP1"],
        "Initial Position": [1, 5, 20],
        "Final Position": [3, 12, 24],
        "Peptide Sequence": ["ABC", "ABCDEFGH", "ABCDE"],  # lengths 3, 8, 5
    })


def test_bepipred_length_min_max():
    df = bepipred.finalresultsBepipred(_bepipred_df(), 4, 6, 0)
    assert set(df["Peptide Sequence"]) == {"ABCDE"}  # only len 5 in [4,6]


def test_bepipred_length_min_only():
    df = bepipred.finalresultsBepipred(_bepipred_df(), 5, 0, 0)
    assert set(df["Peptide Sequence"]) == {"ABCDEFGH", "ABCDE"}  # len >= 5


def test_bepipred_no_length_filter_keeps_all():
    df = bepipred.finalresultsBepipred(_bepipred_df(), 0, 0, 0)
    assert len(df) == 3


def test_others_length_filter(tmp_path):
    f = tmp_path / "x.fasta"
    f.write_text(
        ">P_SARS_M_NP1_1_3\nABC\n"
        ">P_SARS_M_NP1_5_12\nABCDEFGH\n"
        ">P_SARS_M_NP1_20_24\nABCDE\n"
    )
    df = others.fasta_epitopes(str(f), 4, 6)
    assert set(df["Peptide Sequence"]) == {"ABCDE"}


def test_mhcii_ic50_threshold_variation(tmp_path):
    header = ["allele", "seq_num", "start", "end", "length", "core_peptide", "peptide",
              "method", "smm_align_ic50", "nn_align_ic50", "nn_align_rank", "nn_align_adjusted_rank"]
    rows = [
        ["HLA-DRB1*01:01", "1", "3", "17", "15", "c1", "PEPTIDE_STRONG", "nn_align", "10", "20", "1", "1"],
        ["HLA-DRB1*01:01", "1", "3", "17", "15", "c2", "PEPTIDE_MEDIUM", "nn_align", "10", "300", "5", "5"],
    ]
    df_raw = pd.DataFrame(rows, columns=header)

    # threshold 50 -> only the 20 nM peptide
    strong = mhcii.MHCIIAntigenEpitopes(df_raw, "SARS", "SPIKE", "DR", 50)
    assert set(strong["Peptide Sequence"]) == {"PEPTIDE_STRONG"}

    # threshold 500 -> both survive
    both = mhcii.MHCIIAntigenEpitopes(df_raw, "SARS", "SPIKE", "DR", 500)
    assert set(both["Peptide Sequence"]) == {"PEPTIDE_STRONG", "PEPTIDE_MEDIUM"}


def test_conservancy_max_identity_and_seqmatch(tmp_path):
    csv = (
        "Epitope #,Epitope name,Epitope sequence,Epitope length,"
        "Percent of protein sequence matches at identity <= 100%,Minimum identity,Maximum identity,View details\n"
        "1,A_B_M_1_3,AYI,3,100.00% (2/2),80.00%,100.00%,d\n"   # passes at imax 100
        "2,A_B_M_4_6,LMN,3,100.00% (2/2),70.00%,95.00%,d\n"    # excluded when imax=90
        "3,A_B_M_7_9,QRS,3,20.00% (1/5),80.00%,100.00%,d\n"    # excluded by seq_match 60
    )
    d = tmp_path / "csvs"
    d.mkdir()
    (d / "sp.csv").write_text(csv)

    # imax=90 removes row1 (max 100) and row3 (max 100); row2 max 95 also > 90 -> all removed
    df_low_max = conservancy.EpitConservAnalysis(70, "<", 60, 90, 60, str(d))
    assert df_low_max.empty

    # imax=100, seq_match=60 -> row3 removed (20% < 60); rows 1 and 2 remain
    df = conservancy.EpitConservAnalysis(70, ">=", 60, 100, 60, str(d))
    assert set(df["Epitope name"]) == {"A_B_M_1_3", "A_B_M_4_6"}
