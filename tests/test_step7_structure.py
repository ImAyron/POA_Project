"""A etapa 7 buscando a estrutura e resolvendo a numeração pelo SIFTS.

Duas coisas que o repositório tinha pela metade:

* ``map_epitopes`` aceitava ``sifts_offset`` e **nenhum chamador passava um**. O POA0 calculava e
  gravava o mapeamento (``estrutura.offset_regiao``), a etapa 1 dizia ao usuário que a etapa 7 o
  usava, e a etapa 7 caía no alinhamento — que funciona, mas é estimativa onde havia exatidão.
* Baixar a estrutura do RCSB só existia na etapa 1, enquanto gerava a referência. Quem enviou o
  próprio FASTA, ou quis trocar de estrutura depois, tinha que ir ao site buscar o arquivo.
"""
import json

import pytest

from poa.gui import backend


@pytest.fixture
def ctx(tmp_path, monkeypatch):
    monkeypatch.setenv("POA_RESULTS_DIR", str(tmp_path))
    return backend.create_run(tmp_path, "teste etapa 7")


def _poa0_report(ctx, *, pdb_name="structure_6H3X.pdb", offset=481, covers=True, length=222):
    """O relatório que o POA0 grava, reduzido ao que a etapa 7 lê."""
    (ctx.inputs_dir / pdb_name).write_text("ATOM\nEND\n", encoding="utf-8")
    report = {
        "proteina": {"accession": "A0A0D4BSW3"},
        "regiao": {"start": 482, "end": 703},
        "sequencia": {"inicio": 482, "fim": 703, "tamanho": length},
        "arquivos": {"estrutura": str(ctx.inputs_dir / pdb_name)},
        "estrutura": {
            "pdb_id": pdb_name[len("structure_"):-len(".pdb")],
            "offset_regiao": offset,
            "cobre_regiao": covers,
            "sifts": {"chain": "A", "unp_start": 482, "unp_end": 703, "offset": 0},
        },
    }
    (ctx.inputs_dir / backend.POA0_REPORT_FILENAME).write_text(
        json.dumps(report, ensure_ascii=False), encoding="utf-8")
    return str(ctx.inputs_dir / pdb_name)


# --------------------------------------------------------------- o mapeamento SIFTS
def test_the_recorded_sifts_mapping_is_read_back(ctx):
    """O POA0 grava; até agora ninguém lia."""
    pdb = _poa0_report(ctx)

    mapping = backend.structure_mapping(ctx, pdb)

    assert mapping is not None
    assert mapping["offset"] == 481
    assert mapping["span"] == (1, 222)
    assert mapping["pdb_id"] == "6H3X"
    assert mapping["covers_region"] is True


def test_a_different_structure_file_does_not_borrow_the_mapping(ctx):
    """
    O relatório descreve **uma** estrutura.

    Aplicar o offset dela a outro arquivo colocaria os epítopos com confiança sobre os resíduos
    errados — o mesmo tipo de falha silenciosa do 4C23, só que pior, porque viria com o selo de
    "exato".
    """
    _poa0_report(ctx)
    outra = ctx.inputs_dir / "structure_4C23.pdb"
    outra.write_text("ATOM\nEND\n", encoding="utf-8")

    assert backend.structure_mapping(ctx, str(outra)) is None


def test_no_poa0_report_means_no_mapping(ctx):
    """Sem POA0 não há região registrada, e o alinhamento continua sendo a resposta certa."""
    pdb = ctx.inputs_dir / "structure_6H3X.pdb"
    pdb.write_text("ATOM\nEND\n", encoding="utf-8")

    assert backend.structure_mapping(ctx, str(pdb)) is None
    assert backend.structure_mapping(ctx, None) is None


def test_a_report_without_sifts_falls_back_to_alignment(ctx):
    """PDBe nem sempre tem mapeamento; o POA0 avisa e segue, e a etapa 7 também deve."""
    pdb = _poa0_report(ctx)
    report = json.loads((ctx.inputs_dir / backend.POA0_REPORT_FILENAME).read_text(encoding="utf-8"))
    del report["estrutura"]["offset_regiao"]
    (ctx.inputs_dir / backend.POA0_REPORT_FILENAME).write_text(json.dumps(report), encoding="utf-8")

    assert backend.structure_mapping(ctx, pdb) is None


def test_a_partial_structure_is_flagged(ctx):
    pdb = _poa0_report(ctx, covers=False)

    assert backend.structure_mapping(ctx, pdb)["covers_region"] is False


def test_an_unreadable_report_is_ignored_not_fatal(ctx):
    pdb = ctx.inputs_dir / "structure_6H3X.pdb"
    pdb.write_text("ATOM\nEND\n", encoding="utf-8")
    (ctx.inputs_dir / backend.POA0_REPORT_FILENAME).write_text("{nao é json", encoding="utf-8")

    assert backend.poa0_report(ctx) is None
    assert backend.structure_mapping(ctx, str(pdb)) is None


# --------------------------------------------------------------- o accession do cabeçalho
@pytest.mark.parametrize("header, expected", [
    ("NS1_DENV1_P29990", "P29990"),            # o caso do POA0
    ("GCHEAD_OROV_A0A0D4BSW3", "A0A0D4BSW3"),  # accession de 10 caracteres
    ("NS1_DENV1_P29990_extra", "P29990"),      # campo composto, primeiro pedaço
    ("NS1_DENV1_ref", ""),                     # 'ref' não é accession
    ("NS1_DENV1", ""),                         # sem terceiro campo
    ("NS1_DENV1_12345", ""),                   # não começa com letra
])
def test_accession_is_only_taken_when_it_looks_like_one(tmp_path, header, expected):
    """
    Um ID que não é accession é melhor tratado como ausente do que chutado.

    O PDBe responde por accession; um valor errado listaria estruturas de outra proteína, e a
    pessoa não teria como saber — a lista viria preenchida e plausível.
    """

    fasta = tmp_path / "ref.fasta"
    fasta.write_text(f">{header}\nMKTAYIAKQRQ\n", encoding="utf-8")

    protein, specie = header.split("_")[0], header.split("_")[1]
    assert backend.accession_from_header(f"{protein}|{specie}", str(fasta)) == expected


def test_no_reference_means_no_accession():
    """Sem referência não há cabeçalho, e portanto nada a perguntar ao PDBe."""
    assert backend.accession_from_header("NS1|DENV1", None) == ""
