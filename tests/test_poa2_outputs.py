"""POA2 output naming: two different questions must not land on one file.

The POA2 result file used to be named after the sequence-identity threshold alone. But the
threshold is only half of the criterion — ``-g`` selects *conserved* epitopes and ``-l`` selects
*unique* ones, which are opposite questions with opposite answers. Running both at the same
threshold into the same ``-r`` silently replaced the first result with the second, and ``-rf``
(which region of the membrane) collapsed the same way.
"""
from types import SimpleNamespace

import pytest

from poa.core import pipeline
from poa.services import conservancy_client as cc

EPITOPES = ">SARS_SPIKE_Bepipred3.0_4_6\nAYI\n>SARS_SPIKE_Bepipred3.0_1_3\nGGA\n"
PROTEINS = [("SPIKE_SARS_ref", "GGAYIGG"), ("SPIKE_SARS_v2", "GGAYLGG")]


@pytest.fixture
def prepared(tmp_path, monkeypatch):
    """Conservancy CSVs plus a reference FASTA, with the topology step stubbed out."""
    poa1_dir = tmp_path / "poa1" / "Conservancy Analysis"
    poa1_dir.mkdir(parents=True)
    (poa1_dir / "SARS_epitopes.fasta").write_text(EPITOPES, encoding="utf-8")

    csv_dir = tmp_path / "csvs"
    cc.run_conservancy_for_dir(str(poa1_dir), PROTEINS, threshold=70,
                               out_dir=str(csv_dir), operator=">=")

    fasta = tmp_path / "ref.fasta"
    fasta.write_text(">SPIKE_SARS_ref\nGGAYIGG\n", encoding="utf-8")

    # pyTMHMM is an optional native dependency; the naming has nothing to do with it.
    monkeypatch.setattr(pipeline.topology, "tmhmmAnalysis",
                        lambda args, df: df.assign(Portion_Outside=1, Portion_TM=0,
                                                   Portion_Inside=0))
    return csv_dir, fasta, tmp_path


def _args(csv_dir, fasta, out, *, g="", l="", rf=None, t=70):      # noqa: E741 - mirrors the CLI
    return SimpleNamespace(d=str(csv_dir), f=str(fasta), t=t, g=g, l=l, r=str(out),
                           imin=0, imax=100, m=0, rf=rf, idf=False, strict=False)


def test_conserved_and_unique_results_do_not_overwrite_each_other(prepared):
    csv_dir, fasta, tmp_path = prepared
    out = tmp_path / "poa2"

    conserved = pipeline.run_poa2(_args(csv_dir, fasta, out, g="y"))
    # The CSVs were built at '>=', so a '<' run has to be given its own set; reusing the same
    # directory is what the old naming allowed. Build the matching CSVs first.
    poa1_dir = tmp_path / "poa1" / "Conservancy Analysis"
    lt_dir = tmp_path / "csvs_lt"
    cc.run_conservancy_for_dir(str(poa1_dir), PROTEINS, threshold=70,
                               out_dir=str(lt_dir), operator="<")
    unique = pipeline.run_poa2(_args(lt_dir, fasta, out, l="y"))

    assert conserved.xlsx_path != unique.xlsx_path
    assert "ge70" in conserved.xlsx_path and "lt70" in unique.xlsx_path
    from pathlib import Path
    assert Path(conserved.xlsx_path).is_file() and Path(unique.xlsx_path).is_file()


def test_region_filters_each_get_their_own_fasta(prepared):
    """-rf 1 (outside) and -rf 3 (inside) select disjoint sets; one name cannot hold both."""
    csv_dir, fasta, tmp_path = prepared
    out = tmp_path / "poa2"

    outside = pipeline.run_poa2(_args(csv_dir, fasta, out, g="y", rf=1))
    inside = pipeline.run_poa2(_args(csv_dir, fasta, out, g="y", rf=3))

    assert outside.fasta_path != inside.fasta_path
    assert "rf1" in outside.fasta_path and "rf3" in inside.fasta_path


def test_threshold_still_appears_in_the_name(prepared):
    """The threshold stays in the name — the operator is added to it, not substituted for it."""
    csv_dir, fasta, tmp_path = prepared

    result = pipeline.run_poa2(_args(csv_dir, fasta, tmp_path / "poa2", g="y"))

    assert "70" in result.xlsx_path
