"""Tests for the external-tool clients (pure adapters + graceful unavailability)."""
import pytest

from poa.core.parsers import mhcii as mhcii_parser
from poa.core.parsers import papimed
from poa.services import antigenic_client, bepipred_client, mhcii_client, tmhmm_client
from poa.services.base import ServiceUnavailable


# --------------------------------------------------------------------------- MHC-II TSV adapter
def test_mhcii_normalize_maps_aliases_and_fills_missing(tmp_path):
    tsv = (
        "allele\tseq_num\tstart\tend\tlength\tcore_peptide\tpeptide\tic50\trank\n"
        "HLA-DRB1*01:01\t1\t3\t17\t15\tCOREPEPTID\tAAADEFGHIKLMNOP\t20\t1.0\n"
    )
    norm = mhcii_client.normalize_tsv(tsv, requested_method="nn_align")
    header = norm.splitlines()[0].split("\t")
    assert header == mhcii_client.REQUIRED_COLUMNS

    # write normalized file and confirm the existing parser accepts it end-to-end
    path = mhcii_client.write_result_file(tsv, "SPIKE", "SARS", str(tmp_path))
    df_raw = mhcii_parser.MHCIIHTMLconverter(path)
    out = mhcii_parser.MHCIIAntigenEpitopes(df_raw, "SARS", "SPIKE", "DR", 50)
    assert len(out) == 1
    row = out.iloc[0]
    assert row["Peptide Sequence"] == "AAADEFGHIKLMNOP"
    assert row["Allele"] == "HLA-DRB1*01:01"


def test_mhcii_submit_bad_response_detected(monkeypatch):
    class _Resp:
        status_code = 200
        text = "<html>error page</html>"

    import poa.services.mhcii_client as mc

    def fake_post(url, data=None, timeout=None):
        return _Resp()

    import requests
    monkeypatch.setattr(requests, "post", fake_post)
    with pytest.raises(Exception):
        mc.submit("SEQ", ["HLA-DRB1*01:01"], 15)


# --------------------------------------------------------------------------- EMBOSS antigenic adapter
def test_antigenic_gff_to_pap_roundtrip(tmp_path):
    fasta = tmp_path / "prot.fasta"
    fasta.write_text(">SPIKE_SARS_NP1\nMKTAYIAMKGVLMNK\n")  # positions 3..7 => TAYIA
    gff = "SPIKE_SARS_NP1\tantigenic\tSO:0001067\t3\t7\t1.20\t.\t.\tID=1\n"

    feats = antigenic_client.parse_gff(gff)
    assert feats["SPIKE_SARS_NP1"] == [(3, 7, 1.2)]

    pap_txt = antigenic_client.to_pap_txt(gff, str(fasta))
    pap_file = tmp_path / "pap.txt"
    pap_file.write_text(pap_txt)

    df = papimed.PAPepitopes(str(pap_file), 0, 0)
    assert len(df) == 1
    row = df.iloc[0]
    assert row["Peptide Sequence"] == "TAYIA"
    assert row["Protein"] == "SPIKE"
    assert row["Specie"] == "SARS"
    assert row["Method"] == "PAP/IMED"


# --------------------------------------------------------------------------- graceful unavailability
def test_tmhmm_unavailable_raises_clearly():
    if not tmhmm_client.is_available():
        with pytest.raises(ServiceUnavailable):
            tmhmm_client.predict_topology("whatever.fasta")


def test_bepipred_cli_missing_raises():
    import importlib.util
    import shutil

    installed = (shutil.which("bp3") or shutil.which("bepipred3_CLI.py")
                 or importlib.util.find_spec("bp3") is not None)
    if not installed:
        with pytest.raises(ServiceUnavailable):
            bepipred_client._resolve_cli()
