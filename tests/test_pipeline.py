"""End-to-end POA1 test: consolidation + report + per-species FASTA output."""
import os
from types import SimpleNamespace

from poa.core import pipeline
from poa.core.models import STANDARD_COLUMNS


def _poa1_args(tmp_path, b3, x, f, out, export="n"):
    return SimpleNamespace(
        b2="", b3=str(b3), bmin=0, bmax=0,
        p="", pmin=0, pmax=0,
        n="", m="", mhla="DR", mic=50,
        x=str(x), xmin=0, xmax=0,
        d=str(out), f=str(f), e=export,
    )


def test_run_poa1_end_to_end(tmp_path):
    # -b3 input (uppercase = epitope): two epitopes AYI (4-6) and LMN (12-14)
    b3 = tmp_path / "bp3.fasta"
    b3.write_text(">SARS_SPIKE_NP1\nmktAYIamkgvLMNkqrst\n")

    # -x input (other predictor)
    x = tmp_path / "x.fasta"
    x.write_text(">SPIKE_SARS_CTLpred_NP1_3_9\nTAYIAMK\n")

    # -f proteins (long enough that no epitope exceeds the shortest protein)
    f = tmp_path / "proteins.fasta"
    f.write_text(">SPIKE_SARS_NP1\nMKTAYIAMKGVLMNKQRSTAYIAMK\n")

    out = tmp_path / "results"
    out.mkdir()

    result = pipeline.run_poa1(_poa1_args(tmp_path, b3, x, f, out))

    # consolidated predictions: 2 from BepiPred3 + 1 from others
    assert list(result.predictions.columns) == STANDARD_COLUMNS
    assert len(result.predictions) == 3
    assert set(result.predictions["Method"]) == {"Bepipred3.0", "CTLPRED"}

    # report written (UTF-8, so Portuguese accents are preserved)
    assert os.path.exists(result.report_path)
    report_text = open(result.report_path, encoding="utf-8").read()
    assert "POA - Pipeline de Otimização de Antígenos" in report_text
    assert "Total:\t3" in report_text

    # per-species FASTA created (all epitopes are from species SARS)
    conservancy_fasta = os.path.join(result.conservancy_dir, "SARS_epitopes.fasta")
    assert os.path.exists(conservancy_fasta)
    fasta = open(conservancy_fasta).read()
    # headers must be Specie_Protein_Method_Init_Final (exactly 5 fields for POA2)
    assert ">SARS_SPIKE_Bepipred3.0_4_6" in fasta
    assert "AYI" in fasta


def test_run_poa1_requires_a_method(tmp_path):
    f = tmp_path / "proteins.fasta"
    f.write_text(">SPIKE_SARS_NP1\nMKTAYIAMK\n")
    out = tmp_path / "results"
    out.mkdir()
    args = _poa1_args(tmp_path, b3="", x="", f=f, out=out)
    try:
        pipeline.run_poa1(args)
        assert False, "expected an exception when no prediction method is provided"
    except Exception as exc:
        assert "required" in str(exc).lower()


def test_run_poa1_creates_missing_output_dir(tmp_path):
    # -d points to a directory that does not exist yet; the pipeline must create it (even with -e y).
    b3 = tmp_path / "bp3.fasta"
    b3.write_text(">SARS_SPIKE_NP1\nmktAYIamk\n")
    f = tmp_path / "proteins.fasta"
    f.write_text(">SPIKE_SARS_NP1\nMKTAYIAMKQRST\n")
    out = tmp_path / "not_created_yet" / "out"  # nested, does not exist
    args = _poa1_args(tmp_path, b3=b3, x="", f=f, out=out, export="y")
    result = pipeline.run_poa1(args)
    assert out.exists()
    assert os.path.exists(result.report_path)
    assert os.path.exists(os.path.join(str(out), "Bepipred_Epitopes.xlsx"))


def test_run_poa1_xlsx_export(tmp_path):
    b3 = tmp_path / "bp3.fasta"
    b3.write_text(">SARS_SPIKE_NP1\nmktAYIamk\n")
    f = tmp_path / "proteins.fasta"
    f.write_text(">SPIKE_SARS_NP1\nMKTAYIAMKQRST\n")
    out = tmp_path / "results"
    out.mkdir()
    args = _poa1_args(tmp_path, b3=b3, x="", f=f, out=out, export="y")
    pipeline.run_poa1(args)
    assert os.path.exists(os.path.join(str(out), "Bepipred_Epitopes.xlsx"))
