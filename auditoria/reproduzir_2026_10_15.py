r"""Reproducoes isoladas da auditoria; nao altera dados reais nem consulta a rede.

Executar da raiz: .venv\Scripts\python.exe auditoria/reproduzir_2026_10_15.py
As assercoes documentam o comportamento defeituoso observado nesta revisao.
Nao sao testes de aceitacao: devem mudar quando os defeitos forem corrigidos.
"""
from __future__ import annotations

import contextlib
import io
import sys
import tempfile
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd
from Bio import SeqIO

from poa.cli.poa0 import build_parser
from poa.core import dbsearch, pipeline, structure_map
from poa.core.parsers import conservancy as cp, mhcii
from poa.gui import backend, viz
from poa.services import conservancy_client as cc, mhcii_client


def put(root, name, content):
    path = root / name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")
    return str(path)


def pdb(sequence, chain):
    names = {"A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU"}
    return "\n".join(
        f"ATOM  {i:>5}  CA  {names[aa]} {chain}{i:>4}    "
        f"{float(i):8.3f}{0.:8.3f}{0.:8.3f}  1.00  0.00           C"
        for i, aa in enumerate(sequence, 1)
    ) + "\n"


def run(root):
    # 01: Uma reexecucao POA1 remove B da tabela, mas nao do disco/etapa seguinte.
    refs = put(root, "refs.fasta", ">E_A_ref\nACDE\n>E_B_ref\nACDE\n")
    pred = put(root, "pred.fasta", ">E_A_M_1_4\nACDE\n>E_B_M_1_4\nACDE\n")
    out = root / "poa1"
    args = backend.build_poa1_args({"x": pred}, {}, refs, str(out))
    pipeline.run_poa1(args)
    put(root, "pred.fasta", ">E_A_M_1_4\nACDE\n")
    result = pipeline.run_poa1(args)
    csv_dir = root / "csvs"
    cc.run_conservancy_for_dir(result.conservancy_dir, refs, 70, str(csv_dir))
    assert result.predictions.Specie.tolist() == ["A"]
    assert (csv_dir / "B_conservancy.csv").exists()
    assert cp.stale_csvs(str(csv_dir)) == []
    print("01 CONFIRMADO: B removida do POA1 reaparece na conservancia; detector de sobras vazio.")

    # 02: Metadados antigos anulam a validacao dos cabecalhos ate com strict=True.
    meta_dir = root / "metadata"
    rows = cc.conservancy_rows([("A_E_M_1_4", "ACDE")], [("E_A_ref", "ACDE")], 90)
    put(meta_dir, "A.csv", cc.rows_to_csv(rows))
    cp.write_metadata(str(meta_dir), 70, ">=", species={"A": {"csv": "A.csv"}})
    selected = cp.EpitConservAnalysis(70, ">=", 0, 100, 0, str(meta_dir), strict=True)
    assert any("threshold >=70" in col for col in selected.columns)
    assert len(selected) == 1
    print("02 CONFIRMADO: CSV >=90 foi aceito e rotulado >=70 mesmo com strict=True.")

    # 03: O nome do arquivo MHC-II nao distingue IDs da mesma proteina/especie.
    a = mhcii_client.write_result_file("resultado ref1", "E", "A", str(root / "mhcii"), normalize=False)
    b = mhcii_client.write_result_file("resultado ref2", "E", "A", str(root / "mhcii"), normalize=False)
    assert a == b and Path(a).read_text() == "resultado ref2"
    print("03 CONFIRMADO: segunda referencia MHC-II substitui integralmente a primeira.")

    # 04: O merge de referencias elimina conflitos sem avisar.
    a = put(root, "a.fasta", ">E_A_ref\nAAAA\n")
    b = put(root, "b.fasta", ">E_A_ref\nCCCC\n")
    merged = backend.merge_fastas([a, b], str(root / "merged.fasta"))
    records = list(SeqIO.parse(merged, "fasta"))
    assert len(records) == 1 and str(records[0].seq) == "AAAA"
    print("04 CONFIRMADO: duas sequencias distintas com mesmo ID viram uma sem erro.")

    # 05: A referencia 3D e escolhida por ordem alfabetica entre TODOS os FASTAs.
    ctx = backend.WorkContext(root / "viz")
    put(ctx.inputs_dir, "comparison_A.fasta", ">E_A_variant\nCCCC\n")
    put(ctx.inputs_dir, "ref_0_original.fasta", ">E_A_ref\nAAAA\n")
    chosen = backend.reference_sequence(ctx, "E", "A")
    assert chosen["sequence"] == "CCCC"
    print("05 CONFIRMADO: referencia 3D escolhida foi comparison_A.fasta, nao a referencia original.")

    # 06: POA0 aceita coordenadas fora da sequencia e grava FASTA vazio/truncado.
    entry = {"sequence": {"value": "ACDE", "length": 4}, "organism": {"scientificName": "Species A"}}
    lengths = []
    with patch.object(dbsearch.db, "get_entry", return_value=entry):
        for region in ("0-3", "4-2", "2-999"):
            args = build_parser().parse_args(["-a", "TEST", "-r", region, "--nome", "E", "-d", str(root / region)])
            report = dbsearch.run_poa0(args)
            lengths.append(report["sequencia"]["tamanho"])
    assert lengths == [0, 0, 3]
    print("06 CONFIRMADO: intervalos 0-3, 4-2, 2-999 aceitos; comprimentos gravados =", lengths)

    # 07: os.walk localiza arquivos aninhados mas monta caminhos na raiz.
    nested = put(root / "nested_mhcii", "sub/E_A.html", "allele\n")
    _, _, paths = mhcii.files_map_MHCIIBD(str(root / "nested_mhcii"))
    assert paths and paths[0] != nested and not Path(paths[0]).exists()
    print("07 CONFIRMADO: MHC-II em subpasta devolve caminho inexistente.")

    # 08: Upload com ../ sai da pasta designada (apenas dentro do temporario de auditoria).
    dest = root / "uploads"
    saved = Path(backend.save_upload(dest, "../fora.txt", b"audit"))
    assert saved.resolve() == (root / "fora.txt").resolve()
    print("08 CONFIRMADO: save_upload aceita ../ e grava fora da pasta designada.")

    # 09: Nova analise com nome equivalente reutiliza a pasta com resultados anteriores.
    first = backend.create_run(root / "runs", "Teste A")
    put(first.poa2_out_dir, "anterior.txt", "resultado anterior")
    second = backend.create_run(root / "runs", "Teste-A")
    assert first.root == second.root and (second.poa2_out_dir / "anterior.txt").exists()
    print("09 CONFIRMADO: criar 'Teste-A' reutiliza 'Teste A', inclusive suas saidas.")

    # 10: Mapeia cadeia A, mas o seletor do visualizador destaca todas as cadeias.
    text = pdb("ACDE", "A") + pdb("CCCC", "B")
    report = structure_map.map_epitopes(text, "ACDE", [(1, 4)], chain=None)
    assert report.ok and report.chain == "A"
    with patch("py3Dmol.view") as view:
        viz.build_3dmol_view_html(text, report.ranges, chain=None)
        selector = view.return_value.addStyle.call_args.args[0]
        assert "chain" not in selector
    print("10 CONFIRMADO: apenas A foi validada, mas seletor 3D nao restringe cadeia.")

    # 11: Uma proteina curta NAO relacionada bloqueia epitopo valido de outra especie.
    refs = put(root, "viability.fasta", ">E_A_ref\nACDE\n>E_B_ref\nAC\n")
    pred = put(root, "viability_pred.fasta", ">E_A_M_1_4\nACDE\n")
    args = backend.build_poa1_args({"x": pred}, {}, refs, str(root / "viability"))
    try:
        pipeline.run_poa1(args)
    except Exception as exc:
        assert "length greater" in str(exc)
    else:
        raise AssertionError("Esperava bloqueio indevido da referencia curta")
    print("11 CONFIRMADO: referencia curta de B bloqueia epitopo valido de A.")

    # 12: Metadados manuais sem inventario classificam TODO CSV como sobra.
    cp.write_metadata(str(meta_dir), 90, ">=", source="iedb")
    assert cp.stale_csvs(str(meta_dir)) == ["A.csv"]
    try:
        cp.EpitConservAnalysis(90, ">=", 0, 100, 0, str(meta_dir), strict=True)
    except ValueError as exc:
        assert "not produced" in str(exc)
    else:
        raise AssertionError("Esperava falso positivo de CSV residual")
    print("12 CONFIRMADO: CSV manual com metadados corretos e recusado como sobra em strict.")

    # 13: NaN faz a comparacao de diferencas retornar False e o validador aprovar.
    from tools import validate_conservancy

    ep = put(root, "validate_ep.fasta", ">A_E_M_1_4\nACDE\n")
    prot = put(root, "validate_prot.fasta", ">E_A_ref\nACDE\n")
    rows = cc.conservancy_rows([("A_E_M_1_4", "ACDE")], [("E_A_ref", "ACDE")], 70)
    for column in (cp.percent_column(">=", 70), "Minimum identity", "Maximum identity"):
        rows[0][column] = "NaN"
    invalid_csv = put(root, "invalid.csv", cc.rows_to_csv(rows))
    with contextlib.redirect_stdout(io.StringIO()):
        code = validate_conservancy.main(["-i", invalid_csv, "-e", ep, "-p", prot, "-t", "70"])
    assert code == 0
    print("13 CONFIRMADO: validador declara equivalencia mesmo com tres numeros NaN no CSV.")


def run_gui(root):
    from streamlit.testing.v1 import AppTest
    from poa import logging_conf
    from poa.core.models import STANDARD_COLUMNS

    app = str(Path(__file__).resolve().parents[1] / "poa/gui/app.py")
    ctx = backend.WorkContext(root / "gui")
    refs = put(ctx.inputs_dir, "ref.fasta", ">E_A_ref\nACDE\n")
    pred = put(ctx.inputs_dir, "pred.fasta", ">E_A_M_1_4\nACDE\n")
    first = backend.run_poa1(ctx, {"x": pred}, {}, refs)
    prior = pipeline.Poa2Result(pd.DataFrame(), "anterior.xlsx", None, ">=")

    # Impede handlers persistentes do Streamlit de prender o diretorio temporario no Windows.
    with patch.object(logging_conf, "use_log_file", return_value=None):
        at = AppTest.from_file(app, default_timeout=90)
        for key, value in {"run_path": str(ctx.root), "proteins_path": refs,
                           "poa1_sources": {"x": [pred]}, "poa1_result": first,
                           "poa2_result": prior, "conservancy_ready": True}.items():
            at.session_state[key] = value
        at.run()
        at.sidebar.radio[0].set_value("3 · POA1 & Relatório").run()
        assert not at.exception
        next(b for b in at.button if "Executar POA1" in b.label).click().run()
        assert not at.exception
        assert at.session_state["poa2_result"].xlsx_path == "anterior.xlsx"
        assert at.session_state["conservancy_ready"] is False
        print("14 CONFIRMADO: reexecutar POA1 mantem o resultado POA2 da execucao anterior.")

        structure = put(ctx.inputs_dir, "structure.pdb", pdb("ACDE", "A"))
        at.session_state["poa1_result"] = pipeline.Poa1Result(
            pd.DataFrame(columns=STANDARD_COLUMNS), first.report_path, first.conservancy_dir)
        at.session_state["pdb_path"] = structure
        at.sidebar.radio[0].set_value("7 · Visualização (2D/3D)").run()
        assert any("split" in str(e.message) for e in at.exception)
        print("15 CONFIRMADO: POA1 vazio com PDB salvo quebra etapa 7 com None.split.")


if __name__ == "__main__":
    with tempfile.TemporaryDirectory(prefix="poa_auditoria_") as temp:
        run(Path(temp))
        if "--gui" in sys.argv:
            run_gui(Path(temp))
    print("Reproducoes concluidas. Dados sinteticos temporarios removidos.")
