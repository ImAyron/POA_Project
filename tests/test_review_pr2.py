"""Regressões encontradas revisando o PR #2 (POA0, guard de estrutura, etapas 4 e 5).

Cada caso aqui foi **medido** na branch antes da correção, não deduzido: o comentário de cada
teste diz o que acontecia. O guard da estrutura é o mais sério — ele aprovava justamente a
estrutura errada que o PR existe para barrar.
"""
import random

import pytest

from poa.core import dbsearch, structure_map
from poa.core.structure_map import DEFAULT_MIN_COVERAGE, DEFAULT_MIN_IDENTITY
from poa.services import conservancy_client as cc

AMINOACIDS = "ACDEFGHIKLMNPQRSTVWY"


def _random_protein(length: int, seed: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice(AMINOACIDS) for _ in range(length))


def _pdb(sequence: str, chain: str = "A", first: int = 1) -> str:
    """Um PDB mínimo: um CA por resíduo, que é o que ``chain_residues`` lê."""
    three = {one: code for code, one in structure_map.THREE_TO_ONE.items() if len(one) == 1}
    lines = []
    for offset, letter in enumerate(sequence):
        lines.append(
            f"ATOM  {offset + 1:5d}  CA  {three.get(letter, 'ALA'):>3s} {chain}"
            f"{first + offset:4d}       0.000   0.000   0.000  1.00  0.00           C")
    return "\n".join(lines) + "\nEND\n"


# --------------------------------------------------------------- guard da estrutura
def test_structure_sharing_a_few_residues_is_refused():
    """
    Identidade alta sobre pouquíssimas colunas não é evidência de nada.

    Os dois gaps terminais são livres, o que é certo para uma estrutura que cobre parte de uma
    proteína maior, mas deixa o alinhador descartar quase toda a estrutura como sobra e guardar um
    trecho curto bem pareado. Como a aprovação só exigia ``columns > 0``, isso passava.

    Medido nesta branch antes da correção: de 216 pares de sequências aleatórias sem relação, um
    foi aceito com 100% de identidade sobre **uma única coluna** — 0,7% de uma estrutura de 150
    resíduos contra uma referência de 500 aa.
    """
    # As sementes são as do caso que a varredura encontrou, não escolhidas à toa: com elas o
    # alinhador devolve exatamente 1 coluna a 100%. Outras sementes dão identidade baixa e o
    # guard já as recusava pelo critério de identidade — este teste existe para o caso em que a
    # identidade **passa**.
    reference = _random_protein(500, seed=15)
    structure = _random_protein(150, seed=16)

    index_map, identity, columns = structure_map.align_sequences(reference, structure)
    assert identity >= DEFAULT_MIN_IDENTITY and columns == 1, (
        "a premissa do teste mudou: o alinhador nao devolve mais 1 coluna a 100% para este par")

    report = structure_map.map_epitopes(_pdb(structure), reference, [(251, 255)])

    assert not report.ok, "uma estrutura sem relacao foi aceita"
    assert report.aligned_columns < DEFAULT_MIN_COVERAGE * len(structure)
    assert any("alinharam" in w for w in report.warnings), report.warnings


def test_a_structure_of_the_right_protein_still_passes():
    """
    A correção não pode recusar o caso legítimo: um cristal que cobre parte da proteína.

    É exatamente o que o projeto tem — poliproteínas virais com estrutura resolvida de um domínio.
    """
    reference = _random_protein(400, seed=3)
    domain = reference[100:300]            # 200 de 400 resíduos, idênticos

    report = structure_map.map_epitopes(_pdb(domain, first=101), reference, [(120, 130)])

    assert report.ok, report.warnings
    assert report.identity == pytest.approx(1.0)
    assert report.aligned_columns == len(domain)
    assert report.ranges


def test_min_identity_is_honoured_on_the_sifts_path():
    """
    ``min_identity`` era ignorado quando havia mapeamento SIFTS.

    A comparação usava a constante do módulo, então ``check_structure.py -i 90`` não tinha efeito
    nenhum nesse caminho — a única rota que o POA0 pretende usar para estruturas que ele baixa.
    """
    reference = "MKTAYIAKQR"
    structure = "MKTAYIAKQR"[:5] + "WWWWW"     # 50% de identidade sobre o trecho coberto

    lenient = structure_map.map_epitopes(_pdb(structure), reference, [(1, 5)],
                                         sifts_offset=0, sifts_span=(1, 10), min_identity=0.4)
    strict = structure_map.map_epitopes(_pdb(structure), reference, [(1, 5)],
                                        sifts_offset=0, sifts_span=(1, 10), min_identity=0.9)

    assert lenient.ok, lenient.warnings
    assert not strict.ok
    assert strict.identity == pytest.approx(lenient.identity)


# --------------------------------------------------------------- POA0
def test_region_outside_the_sequence_is_refused(monkeypatch, tmp_path):
    """
    ``-r`` fora dos limites gerava um FASTA vazio e saía com código 0.

    Fatiamento em Python não reclama: ``sequence[4999:6000]`` de uma proteína de 1420 aa é ``""``.
    O POA0 escrevia o cabeçalho sem sequência embaixo, registrava ``tamanho: 0``, e o POA1 rodava
    sobre uma referência vazia. Alcançável tanto por engano no ``-r`` quanto por uma região do
    InterPro descrita sobre outra versão da sequência do UniProt.
    """
    sequence = "M" * 1420
    entry = {"primaryAccession": "P00001", "sequence": {"value": sequence},
             "organism": {"scientificName": "Teste virus"},
             "proteinDescription": {"recommendedName": {"fullName": {"value": "Prot"}}}}
    monkeypatch.setattr(dbsearch.db, "get_entry", lambda acc: entry)
    monkeypatch.setattr(dbsearch.db, "get_sequence", lambda acc, entry=None: sequence)
    monkeypatch.setattr(dbsearch, "resolve_region", lambda *a, **k: dbsearch.db.Region(
        name="X", start=5000, end=6000, source="manual", description=""))

    args = type("Args", (), {"a": "P00001", "s": None, "r": "5000-6000", "d": str(tmp_path),
                             "nome": None, "sigla": "TST", "estrutura": None, "limite": 25})()

    with pytest.raises(SystemExit, match="nao cabe"):
        dbsearch.run_poa0(args)

    assert not list(tmp_path.glob("*.fasta")), "nenhum FASTA deve ser escrito"


# --------------------------------------------------------------- conservância
def test_world_fasta_with_an_empty_record_is_skipped(tmp_path, caplog):
    """
    ``max_identity`` passou a devolver ``None`` e um chamador não foi atualizado.

    ``best_world_fasta`` fazia ``max(ident for ...)``; com um ``None`` na lista isso levantava
    ``TypeError: '>' not supported between instances of 'float' and 'NoneType'``. Dispara quando o
    primeiro registro de um world FASTA está vazio — caso que a etapa 1 já avisa que acontece.
    """
    from poa.services.realdata_import import best_world_match

    (tmp_path / "vazio.fasta").write_text(">SO_CABECALHO\n", encoding="utf-8")
    (tmp_path / "bom.fasta").write_text(">P_SP_x\nMKTAYIAKQRQ\n", encoding="utf-8")

    path, identity = best_world_match("MKTAYIAKQRQ", str(tmp_path))

    assert path is not None and path.endswith("bom.fasta")
    assert identity == pytest.approx(100.0)


def test_metadata_records_the_denominator_the_csv_used(tmp_path):
    """
    O sidecar contradizia o CSV ao lado dele.

    Proteínas mais curtas que o epítopo passaram a sair do denominador, por epítopo, mas o
    sidecar só gravava ``n_proteins``. Um conjunto de 2 com uma proteína curta escrevia
    ``100.00% (1/1)`` no CSV e ``n_proteins: 2`` no sidecar, e é do sidecar que a grade de 100/n
    do aviso do ``-m`` e da tabela da etapa 5 é derivada.
    """
    epitopes = [("EP_SP_M_1_9", "AYIAMKTAY")]
    proteins = [("P_SP_a", "GGAYIAMKTAYGG"), ("P_SP_curta", "MKT")]
    rows = cc.conservancy_rows(epitopes, proteins, 70, ">=")

    assert cc._comparable_count(rows[0]) == 1, "o CSV comparou 1 proteína, não 2"

    epitope_dir = tmp_path / "epitopes"
    epitope_dir.mkdir()
    (epitope_dir / "SP_epitopes.fasta").write_text(">EP_SP_M_1_9\nAYIAMKTAY\n", encoding="utf-8")
    cc.run_conservancy_for_dir(str(epitope_dir), proteins, threshold=70,
                               out_dir=str(tmp_path / "csv"), operator=">=")

    from poa.core.parsers.conservancy import read_metadata
    meta = read_metadata(str(tmp_path / "csv"))
    assert meta.species["SP"]["n_proteins"] == 2          # o conjunto oferecido
    assert meta.species["SP"]["n_comparable"] == 1        # o que o CSV de fato usou
