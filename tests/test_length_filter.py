"""The length bound that reaches every method (POA1 ``-lmin``/``-lmax``).

``-bmin``/``-pmin``/``-xmin`` each filter one method, and NetCTL (``-n``) and MHC-II (``-m``) have
no bound of their own — so a minimum length asked for through those fields left every short
epitope coming from those two in the consolidated table. Measured on a real PAP/IMED + MHC-II run
(``results/15112026-TESTE_15_11``): ``-pmin 16`` kept all 677 MHC-II 15-mers.
"""
import pandas as pd
import pytest
from types import SimpleNamespace

from poa.core import pipeline
from poa.gui import backend


# --------------------------------------------------------------------------- unit
def _frame(*sequences):
    return pd.DataFrame({
        "Method": ["M"] * len(sequences),
        "Specie": ["SARS"] * len(sequences),
        "Protein": ["SPIKE"] * len(sequences),
        "ID_Sequence": ["NP1"] * len(sequences),
        "Initial Position": range(1, len(sequences) + 1),
        "Final Position": range(1, len(sequences) + 1),
        "Peptide Sequence": list(sequences),
    })


def test_filter_by_length_bounds():
    frame = _frame("AA", "AAAAA", "AAAAAAAAAA")

    assert list(pipeline.filter_by_length(frame, 5, 0)["Peptide Sequence"]) == [
        "AAAAA", "AAAAAAAAAA"]
    assert list(pipeline.filter_by_length(frame, 0, 5)["Peptide Sequence"]) == ["AA", "AAAAA"]
    assert list(pipeline.filter_by_length(frame, 5, 5)["Peptide Sequence"]) == ["AAAAA"]


def test_filter_by_length_zero_disables_each_bound():
    frame = _frame("AA", "AAAAA")
    assert len(pipeline.filter_by_length(frame, 0, 0)) == 2


def test_filter_by_length_reindexes():
    # Downstream code addresses rows positionally (bepipred's assembly loop does), so a filtered
    # frame must not carry holes in its index.
    out = pipeline.filter_by_length(_frame("AA", "AAAAA", "AAA"), 3, 0)
    assert list(out.index) == [0, 1]


def test_filter_by_length_handles_empty_frame():
    empty = _frame()
    assert pipeline.filter_by_length(empty, 5, 0).empty


# --------------------------------------------------------------------------- the regression
def _netctl_html(path, *peptides):
    """A NetCTL result page with one '<-E' epitope line per peptide."""
    lines = [
        f"{i} id SPIKE_SARS pep {pep} aff 0.1 0.2 cle 0.3 tap 0.4 t2 COMB 0.9 <-E"
        for i, pep in enumerate(peptides, start=1)
    ]
    path.write_text("NetCTL 1.2 predictions\n" + "\n".join(lines) + "\n")
    return str(path)


def _args(tmp_path, **overrides):
    out = tmp_path / "results"
    out.mkdir(exist_ok=True)
    f = tmp_path / "proteins.fasta"
    f.write_text(">SPIKE_SARS_NP1\n" + "MKTAYIAMKGVLMNKQRSTAYIAMKDEFGHIKLMNPQRSTVWY" * 3 + "\n")
    base = dict(b2="", b3="", bmin=0, bmax=0, p="", pmin=0, pmax=0,
                n="", m="", mhla="DR", mic=50, x="", xmin=0, xmax=0,
                lmin=0, lmax=0, d=str(out), f=str(f), e="n")
    base.update(overrides)
    return SimpleNamespace(**base)


def test_per_method_bound_does_not_reach_netctl(tmp_path):
    """The bug: -xmin is a bound on 'others', and NetCTL epitopes walk past it."""
    netctl = _netctl_html(tmp_path / "netctl.html", "AAA", "DEFGHIKLM")

    df = pipeline.consolidate_predictions(_args(tmp_path, n=netctl, xmin=5))

    assert sorted(df["Peptide Sequence"]) == ["AAA", "DEFGHIKLM"]


def test_global_bound_reaches_netctl(tmp_path):
    """The fix: -lmin applies to every method, NetCTL included."""
    netctl = _netctl_html(tmp_path / "netctl.html", "AAA", "DEFGHIKLM")

    df = pipeline.consolidate_predictions(_args(tmp_path, n=netctl, lmin=5))

    assert list(df["Peptide Sequence"]) == ["DEFGHIKLM"]


def test_global_bound_reaches_mhcii(tmp_path):
    """The other method with no bound of its own — the one in the user's run."""
    directory = tmp_path / "mhcii"
    directory.mkdir()
    header = ["allele", "seq_num", "start", "end", "length", "core_peptide", "peptide",
              "method", "smm_align_ic50", "nn_align_ic50", "nn_align_rank",
              "nn_align_adjusted_rank"]
    rows = [
        ["HLA-DRB1*01:01", "1", "3", "17", "15", "AAADEFGHI", "AAADEFGHIKLMNOP",
         "nn_align", "100", "20", "1.0", "1.0"],
        ["HLA-DRB1*01:01", "1", "20", "22", "3", "QRS", "QRS",
         "nn_align", "30", "10", "0.5", "0.5"],
    ]
    (directory / "SPIKE_SARS.html").write_text(
        "\t".join(header) + "\n" + "\n".join("\t".join(r) for r in rows) + "\n")

    unfiltered = pipeline.consolidate_predictions(_args(tmp_path, m=str(directory)))
    assert len(unfiltered) == 2

    filtered = pipeline.consolidate_predictions(_args(tmp_path, m=str(directory), lmin=5))
    assert list(filtered["Peptide Sequence"]) == ["AAADEFGHIKLMNOP"]


def test_global_bound_adds_to_the_per_method_one(tmp_path):
    """Both bounds apply; the surviving set is the stricter of the two."""
    x = tmp_path / "x.fasta"
    x.write_text(">SPIKE_SARS_CTLPRED_NP1_1_3\nAAA\n"
                 ">SPIKE_SARS_CTLPRED_NP1_1_6\nAAADEF\n"
                 ">SPIKE_SARS_CTLPRED_NP1_1_9\nAAADEFGHI\n")

    # -xmin 4 would keep the 6- and 9-mers; -lmax 8 removes the 9-mer as well.
    df = pipeline.consolidate_predictions(_args(tmp_path, x=str(x), xmin=4, lmax=8))

    assert list(df["Peptide Sequence"]) == ["AAADEF"]


def test_global_bound_filters_the_per_method_xlsx_export(tmp_path):
    """'-e' writes from the per-method frames, which the filter must reach too."""
    netctl = _netctl_html(tmp_path / "netctl.html", "AAA", "DEFGHIKLM")
    args = _args(tmp_path, n=netctl, lmin=5, e="y")

    pipeline.consolidate_predictions(args)

    exported = pd.read_excel(f"{args.d}/NetCTL_Epitopes.xlsx")
    assert list(exported["Peptide Sequence"]) == ["DEFGHIKLM"]


def test_run_poa1_without_the_attributes_still_works(tmp_path):
    """An args namespace predating -lmin/-lmax (the CLI of an older script, a stored manifest)."""
    netctl = _netctl_html(tmp_path / "netctl.html", "AAA")
    args = _args(tmp_path, n=netctl)
    del args.lmin, args.lmax

    assert len(pipeline.consolidate_predictions(args)) == 1


# --------------------------------------------------------------------------- GUI plumbing
def test_build_poa1_args_carries_the_global_bound():
    args = backend.build_poa1_args({"n": "netctl.html"},
                                   {"lmin": 5, "lmax": 25}, "proteins.fasta", "out")
    assert (args.lmin, args.lmax) == (5, 25)


def test_build_poa1_args_defaults_the_global_bound_to_off():
    args = backend.build_poa1_args({"n": "netctl.html"}, {}, "proteins.fasta", "out")
    assert (args.lmin, args.lmax) == (0, 0)


def test_cli_accepts_the_global_bound():
    from poa.cli.poa1 import build_parser

    parsed = build_parser().parse_args(
        ["-n", "netctl.html", "-d", "out", "-f", "proteins.fasta", "-lmin", "5", "-lmax", "25"])
    assert (parsed.lmin, parsed.lmax) == (5, 25)
    assert build_parser().parse_args(
        ["-n", "n.html", "-d", "out", "-f", "p.fasta"]).lmin == 0


@pytest.mark.parametrize("key", ["len_min_all", "len_max_all"])
def test_help_text_exists(key):
    from poa.gui.help_texts import HELP

    assert key in HELP and "método" in HELP[key]


# --------------------------------------------------------------------------- step 1 as rendered
streamlit_testing = pytest.importorskip("streamlit.testing.v1",
                                        reason="streamlit not installed")


def _step1(tmp_path, monkeypatch, sources=None, **state):
    """Step 1 open, on a throwaway analysis folder."""
    from pathlib import Path

    monkeypatch.setenv("POA_RESULTS_DIR", str(tmp_path))
    app = str(Path(__file__).resolve().parents[1] / "poa" / "gui" / "app.py")
    at = streamlit_testing.AppTest.from_file(app, default_timeout=90)
    at.session_state["run_path"] = str(tmp_path)
    if sources is not None:
        at.session_state["poa1_sources"] = sources
    for key, value in state.items():
        at.session_state[key] = value
    return at.run()


def _number(at, label):
    return next(n for n in at.number_input if label in n.label)


def test_step1_offers_the_global_bound_and_records_it(tmp_path, monkeypatch):
    at = _step1(tmp_path, monkeypatch)

    for wanted in (5, 9, 7):
        at = _number(at, "Mínimo geral").set_value(wanted).run()
        assert at.session_state["params_poa1"]["lmin"] == wanted
    at = _number(at, "Máximo geral").set_value(25).run()
    assert at.session_state["params_poa1"]["lmax"] == 25
    assert not at.exception


def test_step1_warns_when_a_per_method_bound_cannot_reach_the_active_methods(
        tmp_path, monkeypatch):
    at = _step1(tmp_path, monkeypatch, sources={"m": ["mhcii"]})

    at = _number(at, "Bepipred min").set_value(5).run()

    assert any("não é alcançado" in str(w.value) for w in at.warning)
    assert any("MHC-II" in str(w.value) for w in at.warning)


def test_step1_drops_the_warning_once_the_global_bound_is_set(tmp_path, monkeypatch):
    at = _step1(tmp_path, monkeypatch, sources={"m": ["mhcii"]})
    at = _number(at, "Bepipred min").set_value(5).run()

    at = _number(at, "Mínimo geral").set_value(5).run()

    assert not any("não é alcançado" in str(w.value) for w in at.warning)


def test_step1_does_not_warn_when_every_active_method_has_its_own_bound(tmp_path, monkeypatch):
    at = _step1(tmp_path, monkeypatch, sources={"b3": ["bp3.fasta"]})

    at = _number(at, "Bepipred min").set_value(5).run()

    assert not any("não é alcançado" in str(w.value) for w in at.warning)


def test_step1_refuses_an_impossible_range(tmp_path, monkeypatch):
    at = _step1(tmp_path, monkeypatch)
    at = _number(at, "Mínimo geral").set_value(20).run()

    at = _number(at, "Máximo geral").set_value(5).run()

    assert any("maior que o máximo geral" in str(e.value) for e in at.error)
