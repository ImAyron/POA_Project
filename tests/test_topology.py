"""Tests for the TMHMM per-epitope region math (no pyTMHMM install required)."""
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
