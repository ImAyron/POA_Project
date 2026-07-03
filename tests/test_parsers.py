"""Tests for every POA1 parser: each must produce the standardized epitope columns."""
import json

import pytest

from poa.core.models import STANDARD_COLUMNS
from poa.core.parsers import bepipred, conservancy, mhcii, netctl, others, papimed


# --------------------------------------------------------------------------- BepiPred 2.0 (JSON)
def test_bepipred2_json(tmp_path):
    data = {
        "antigens": {
            "SPIKE_SARS_NP1": {
                "AA": ["M", "K", "T", "A", "Y", "I", "A"],
                "PRED": [0.1, 0.2, 0.9, 0.8, 0.7, 0.2, 0.1],
            }
        }
    }
    f = tmp_path / "bp2.json"
    f.write_text(json.dumps(data))

    raw = bepipred.bp2_JsonAnalysis(str(f))
    df = bepipred.bp2_AntigenEpitopes(raw)
    df = bepipred.finalresultsBepipred(df, 0, 0, 0)

    assert list(df.columns) == STANDARD_COLUMNS
    assert len(df) == 1
    row = df.iloc[0]
    assert row["Method"] == "Bepipred2.0"
    assert row["Specie"] == "SARS"
    assert row["Protein"] == "SPIKE"
    assert row["ID_Sequence"] == "NP1"
    assert row["Peptide Sequence"] == "TAY"
    assert int(row["Initial Position"]) == 3
    assert int(row["Final Position"]) == 5


# --------------------------------------------------------------------------- BepiPred 3.0 (FASTA)
def test_bepipred3_fasta(tmp_path):
    f = tmp_path / "bp3.fasta"
    f.write_text(">SARS_SPIKE_NP1\nmktAYIamkgvLMNkqrst\n")

    df = bepipred.bp3_FastaAnalysis(str(f))
    df = bepipred.finalresultsBepipred(df, 0, 0, 1)

    assert list(df.columns) == STANDARD_COLUMNS
    assert set(df["Peptide Sequence"]) == {"AYI", "LMN"}
    first = df[df["Peptide Sequence"] == "AYI"].iloc[0]
    assert first["Method"] == "Bepipred3.0"
    assert first["Specie"] == "SARS"
    assert first["Protein"] == "SPIKE"
    assert int(first["Initial Position"]) == 4
    assert int(first["Final Position"]) == 6


# --------------------------------------------------------------------------- PAP / IMED (TXT)
def test_papimed_txt(tmp_path):
    content = (
        ">SPIKE_SARS_NP1\n"
        "n\tStart Position\tSequence\tEnd Position\n"
        "1\t3\tTAYIA\t7\n"
        "2\t10\tGVLMNK\t15\n"
    )
    f = tmp_path / "pap.txt"
    f.write_text(content)

    df = papimed.PAPepitopes(str(f), 0, 0)
    assert list(df.columns) == STANDARD_COLUMNS
    assert len(df) == 2
    assert set(df["Peptide Sequence"]) == {"TAYIA", "GVLMNK"}
    r = df.iloc[0]
    assert r["Method"] == "PAP/IMED"
    assert r["Specie"] == "SARS"
    assert r["Protein"] == "SPIKE"
    assert r["ID_Sequence"] == "NP1"


# --------------------------------------------------------------------------- NetCTL (HTML)
@pytest.mark.xfail(
    reason="Known bug in the original netctl parser (preserved verbatim): rows are appended "
    "with the file line number as a non-sequential .loc index (rejected by pandas>=2), and the "
    "marker/append logic ('if last != E: append -') is inconsistent with selecting '<-E', "
    "overflowing the 16 columns. Awaiting user approval to fix the selection logic.",
    strict=True,
    raises=ValueError,
)
def test_netctl_html(tmp_path):
    lines = [
        "NetCTL 1.2 predictions",
        "1 id SPIKE_SARS pep AADEFGHIK aff 0.1 0.2 cle 0.3 tap 0.4 t2 COMB 0.9 <-E",
        "2 id SPIKE_SARS pep BBCDEFGHI aff 0.1 0.2 cle 0.3 tap 0.4 t2 COMB 0.1",
    ]
    f = tmp_path / "netctl.html"
    f.write_text("\n".join(lines) + "\n")

    df = netctl.netctlAntigenEpitopes(str(f))
    assert list(df.columns) == STANDARD_COLUMNS
    assert len(df) == 1  # only the '<-E' line is an epitope
    row = df.iloc[0]
    assert row["Method"] == "NetCTL"
    assert row["Specie"] == "SARS"
    assert row["Protein"] == "SPIKE"
    assert row["Peptide Sequence"] == "AADEFGHIK"
    assert int(row["Initial Position"]) == 1
    assert int(row["Final Position"]) == 9  # 1 + len(9) - 1


# --------------------------------------------------------------------------- MHC-II (dir of TSV)
def test_mhcii_dir(tmp_path):
    directory = tmp_path / "mhcii"
    directory.mkdir()
    # write a clean TSV file directly
    header = ["allele", "seq_num", "start", "end", "length", "core_peptide", "peptide",
              "method", "smm_align_ic50", "nn_align_ic50", "nn_align_rank", "nn_align_adjusted_rank"]
    rows = [
        ["HLA-DRB1*01:01", "1", "3", "17", "15", "AAADEFGHI", "AAADEFGHIKLMNOP", "nn_align", "100", "20", "1.0", "1.0"],
        ["HLA-DRB1*01:01", "1", "20", "34", "15", "QRSTUVWXY", "QRSTUVWXYZABCDE", "nn_align", "600", "300", "9.0", "9.0"],
        ["HLA-DQB1*02:01", "1", "5", "19", "15", "GGGHHHIII", "GGGHHHIIIJJJKKK", "nn_align", "30", "10", "0.5", "0.5"],
    ]
    (directory / "SPIKE_SARS.html").write_text(
        "\t".join(header) + "\n" + "\n".join("\t".join(r) for r in rows) + "\n"
    )

    specie, protein, files = mhcii.files_map_MHCIIBD(str(directory))
    assert specie == ["SARS"] and protein == ["SPIKE"]

    df_raw = mhcii.MHCIIHTMLconverter(files[0])
    out = mhcii.MHCIIAntigenEpitopes(df_raw, specie[0], protein[0], "DR", 50)

    # Only the DR allele with nn_align_ic50 <= 50 survives (row 1)
    assert len(out) == 1
    row = out.iloc[0]
    assert row["Method"] == "MHCII-Binding"
    assert row["Specie"] == "SARS"
    assert row["Protein"] == "SPIKE"
    assert row["Allele"] == "HLA-DRB1*01:01"
    assert row["Peptide Sequence"] == "AAADEFGHIKLMNOP"


# --------------------------------------------------------------------------- Others (FASTA)
def test_others_fasta(tmp_path):
    f = tmp_path / "x.fasta"
    f.write_text(">SPIKE_SARS_CTLpred_NP1_3_9\nTAYIAMK\n")

    df = others.fasta_epitopes(str(f), 0, 0)
    assert list(df.columns) == STANDARD_COLUMNS
    assert len(df) == 1
    row = df.iloc[0]
    assert row["Method"] == "CTLPRED"
    assert row["Specie"] == "SARS"
    assert row["Protein"] == "SPIKE"
    assert row["ID_Sequence"] == "NP1"
    assert row["Peptide Sequence"] == "TAYIAMK"
    assert row["Initial Position"] == "3"
    assert row["Final Position"] == "9"


# --------------------------------------------------------------------------- Conservancy (CSV)
def test_conservancy_filter(tmp_path):
    csv = (
        "Epitope #,Epitope name,Epitope sequence,Epitope length,"
        "Percent of protein sequence matches at identity <= 100%,Minimum identity,Maximum identity,View details\n"
        "1,SARS_SPIKE_Bepipred3.0_4_6,AYI,3,100.00% (2/2),80.00%,100.00%,details\n"
        "2,SARS_SPIKE_Bepipred3.0_12_14,LMN,3,50.00% (1/2),40.00%,60.00%,details\n"
    )
    d = tmp_path / "csvs"
    d.mkdir()
    (d / "SARS.csv").write_text(csv)

    df = conservancy.EpitConservAnalysis(70, ">=", 60, 100, 60, str(d))
    # Row 2 excluded (Minimum identity 40 < 60)
    assert len(df) == 1
    row = df.iloc[0]
    assert row["Epitope name"] == "SARS_SPIKE_Bepipred3.0_4_6"
    assert "Protein(s) sequence match(es) at Sequence identity threshold >=70" in df.columns
    assert row["Minimum identity(%)"] == 80.0
    assert row["Maximum identity(%)"] == 100.0
