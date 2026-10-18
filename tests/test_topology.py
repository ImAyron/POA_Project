"""Tests for the TMHMM per-epitope region math (no pyTMHMM install required)."""
from types import SimpleNamespace

import pandas as pd
import pytest

from poa.core import topology


def test_epit_region_fractions_mixed():
    # annotation o=outer, m=transmembrane, i=inner
    out, tm, ins = topology.epitTMHMMcaract("oommii", 0, 6)
    assert (out, tm, ins) == (round(2 / 6, 4), round(2 / 6, 4), round(2 / 6, 4))


def test_epit_region_fractions_offset():
    # only look at residues [2:6] = "mmii"
    out, tm, ins = topology.epitTMHMMcaract("oommii", 2, 4)
    assert out == 0.0
    assert tm == 0.5
    assert ins == 0.5


def test_range_verify():
    assert topology.rangeVerify(5, 10) is True
    assert topology.rangeVerify(0, 10) is False
    assert topology.rangeVerify(-1, 10) is False
    assert topology.rangeVerify(11, 10) is False


# --------------------------------------------------------------------------- virus matching
def test_same_virus_uses_the_species_field():
    # 'DENV1' is a substring of both, but only one has DENV1 as its species field
    assert topology.sameVirus("DENV1", "E_DENV1_REF") is True
    assert topology.sameVirus("DENV1", "DENV1_NS1") is False   # this is Specie_Protein, inverted
    assert topology.sameVirus("DENV1", "E_DENV2_REF") is False
    # a header with no separator keeps the original substring behaviour
    assert topology.sameVirus("DENV1", "DENV1") is True


# --------------------------------------------------------------------------- row alignment
def _fake_predict(monkeypatch, annotations):
    """Stand in for pyTMHMM: map each -f record id to a canned o/m/i annotation."""
    monkeypatch.setattr(topology, "pyTMHMMpredict",
                        lambda f: (list(annotations), [annotations[k] for k in annotations]))


def _frame(*names_and_seqs):
    return pd.DataFrame({"Epitope name": [n for n, _ in names_and_seqs],
                         "Epitope sequence": [s for _, s in names_and_seqs]})


def test_tmhmm_analysis_one_row_per_epitope_with_duplicated_records(tmp_path, monkeypatch):
    """
    Regression: an epitope matching several -f records used to append several values, so the
    column no longer lined up with the frame ('Length of values (32) does not match (16)').
    """
    f = tmp_path / "proteins_all.fasta"
    # the same protein twice: once with the pipeline convention, once with inverted headers
    f.write_text(">DENV1_NS1\nMKTAYIAKQR\n>E_DENV1_ref\nMKTAYIAKQR\n")
    _fake_predict(monkeypatch, {"DENV1_NS1": "oooooooooo", "E_DENV1_REF": "iiiiiiiiii"})

    df = _frame(("DENV1_E_Bepipred2.0_3_5", "TAY"), ("DENV1_E_Bepipred2.0_6_8", "IAK"))
    out = topology.tmhmmAnalysis(SimpleNamespace(f=str(f)), df)

    assert len(out) == 2
    # only E_DENV1_ref has DENV1 as its species field, so its 'inside' annotation is the one used
    assert list(out["Portion_Inside"]) == [1.0, 1.0]
    assert list(out["Portion_Outside"]) == [0.0, 0.0]


def test_tmhmm_analysis_warns_and_uses_first_on_genuine_duplicates(tmp_path, monkeypatch):
    """Two records of the SAME species still yield one value per row, with a warning."""
    f = tmp_path / "dup.fasta"
    f.write_text(">E_DENV1_ref\nMKTAYIAKQR\n>E_DENV1_copy\nMKTAYIAKQR\n")
    _fake_predict(monkeypatch, {"E_DENV1_REF": "oooooooooo", "E_DENV1_COPY": "iiiiiiiiii"})

    df = _frame(("DENV1_E_Bepipred2.0_3_5", "TAY"))
    with pytest.warns(UserWarning, match="matched 2 protein records"):
        out = topology.tmhmmAnalysis(SimpleNamespace(f=str(f)), df)

    assert len(out) == 1
    assert list(out["Portion_Outside"]) == [1.0]   # the first record wins


def test_tmhmm_analysis_keeps_placeholder_when_nothing_matches(tmp_path, monkeypatch):
    """An unmatched epitope must consume its own row instead of shifting every later one."""
    f = tmp_path / "proteins.fasta"
    f.write_text(">E_DENV1_ref\nMKTAYIAKQR\n")
    _fake_predict(monkeypatch, {"E_DENV1_REF": "ooooommmmm"})

    df = _frame(("DENV1_E_Bepipred2.0_1_3", "MKT"),
                ("DENV2_E_Bepipred2.0_1_3", "WWW"),   # other species, absent from the FASTA
                ("DENV1_E_Bepipred2.0_9_10", "QR"))
    with pytest.warns(UserWarning):
        out = topology.tmhmmAnalysis(SimpleNamespace(f=str(f)), df)

    assert len(out) == 3
    assert list(out["Portion_Outside"]) == [1.0, "-", 0.0]   # row 2 keeps the placeholder
    assert list(out["Portion_TM"]) == [0.0, "-", 1.0]        # row 3 still gets ITS own values
