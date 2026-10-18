"""Regressoes do code-review das correcoes da auditoria — nucleo: backend, report, fasta_out."""
import json
import os
import time
from pathlib import Path

import pandas as pd
import pytest

from poa.core import fasta_out, report
from poa.gui import backend


def put(root, name, text):
    path = Path(root) / name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return str(path)


def settle(path):
    """Backdate mtime past the settle window so the digest cache is allowed to remember it."""
    old = time.time() - 3600
    os.utime(path, (old, old))


def epitope_row(protein, specie, sequence, id_sequence=""):
    return pd.DataFrame([{
        "Method": "M", "Specie": specie, "Protein": protein, "ID_Sequence": id_sequence,
        "Initial Position": "1", "Final Position": str(len(sequence)),
        "Peptide Sequence": sequence,
    }])


# ----------------------------------------------------- achado 1: checkViability silently skipping
def test_01_epitope_with_no_matching_reference_is_not_waved_through():
    """Trava o achado 1: sem referencia para a Proteina/Especie, o epitopo passava sem checagem."""
    # The reference calls the protein ENV; the predictor calls it E — so nothing matches and the
    # scoped limit is empty. A 15 aa epitope against a 10 aa reference set must still be refused.
    with pytest.raises(ValueError, match="no matching reference"):
        report.checkViability({"ENV_DENV1_REF": 10},
                              epitope_row("E", "DENV1", "ACDEFGHIKLMNPQR"),
                              "Peptide Sequence")


def test_01_full_poa1_run_stops_instead_of_reporting_0_of_0(tmp_path):
    """Trava o achado 1 de ponta a ponta: a execucao seguia e a conservancia saia '0.00% (0/0)'."""
    ctx = backend.WorkContext(tmp_path)
    refs = put(ctx.inputs_dir, "ref.fasta", ">ENV_DENV1_ref\nACDEFGHIKL\n")
    pred = put(ctx.inputs_dir, "pred.fasta", ">E_DENV1_M_1_15\nACDEFGHIKLMNPQR\n")
    with pytest.raises(ValueError, match="no matching reference"):
        backend.run_poa1(ctx, {"x": pred}, {}, refs)


def test_01_unmatched_epitope_that_fits_every_reference_still_passes():
    """O achado 1 nao autoriza recusar o que a regra antiga aceitava: cabe na menor referencia."""
    # Pre-audit the rule was 'longer than ANY reference'; a 2 aa epitope breaks no reference here,
    # so restoring that fallback must not invent a new refusal for species without a reference.
    report.checkViability({"ENV_DENV1_REF": 10, "NS1_DENV1_REF": 4},
                          epitope_row("E", "DENV1", "AC"), "Peptide Sequence")


def test_01_matching_reference_still_scopes_the_limit():
    """A correcao da auditoria continua valendo: so a referencia da propria especie limita."""
    # 4 aa fits E_A (4 aa) and must not be refused by the unrelated 2 aa E_B reference.
    report.checkViability({"E_A_REF": 4, "E_B_REF": 2},
                          epitope_row("E", "A", "ACDE"), "Peptide Sequence")
    with pytest.raises(ValueError, match="matching reference"):
        report.checkViability({"E_A_REF": 4, "E_B_REF": 2},
                              epitope_row("E", "B", "ACDE"), "Peptide Sequence")


def test_01_empty_reference_set_has_nothing_to_check():
    """Sem referencias nao ha limite a aplicar — nao pode estourar em min() de conjunto vazio."""
    report.checkViability({}, epitope_row("E", "DENV1", "ACDEFGHIKLMNPQR"), "Peptide Sequence")


# ----------------------------------------------------- achado 2: fingerprint re-hashing everything
def count_hashed_files(monkeypatch):
    """Record every path whose bytes ``input_fingerprint`` actually reads."""
    hashed = []
    real = backend.file_hash
    monkeypatch.setattr(backend, "file_hash",
                        lambda p: (hashed.append(str(p)), real(p))[1])
    return hashed


def test_02_fingerprint_does_not_rehash_an_unchanged_file(tmp_path, monkeypatch):
    """Trava o achado 2: cada chamada relia e fazia SHA-256 de todo byte do arquivo de novo."""
    big = Path(put(tmp_path, "bepipred2.json", "{}" * 4096))
    settle(big)
    hashed = count_hashed_files(monkeypatch)

    first = backend.input_fingerprint([str(big)])
    # app.py fingerprints the same inputs up to four times per Streamlit rerun.
    for _ in range(3):
        assert backend.input_fingerprint([str(big)]) == first
    assert len(hashed) == 1


def test_02_fingerprint_does_not_rehash_an_unchanged_directory(tmp_path, monkeypatch):
    """Trava o achado 2 para diretorios: o rglob re-hasheava todo CSV a cada interacao."""
    csv_dir = tmp_path / "conservancy_csv"
    for i in range(5):
        settle(put(csv_dir, f"S{i}_conservancy.csv", "Epitope #\n1\n"))
    hashed = count_hashed_files(monkeypatch)

    first = backend.input_fingerprint([str(csv_dir)])
    assert backend.input_fingerprint([str(csv_dir)]) == first
    assert len(hashed) == 5


def test_02_cached_digest_still_notices_an_edited_file(tmp_path):
    """Guarda do achado 2: o cache nao pode mascarar uma edicao de mesmo tamanho."""
    target = Path(put(tmp_path, "ref.fasta", ">E_A_ref\nACDE\n"))
    settle(target)
    before = backend.input_fingerprint([str(target)])
    target.write_text(">E_A_ref\nWWWW\n", encoding="utf-8")  # byte count unchanged, mtime is not
    assert backend.input_fingerprint([str(target)]) != before


def test_02_a_file_rewritten_in_the_same_tick_is_not_served_from_cache(tmp_path):
    """Guarda do achado 2: mtime grosseiro esconderia a reescrita no mesmo instante do hash."""
    target = Path(put(tmp_path, "fresh.fasta", ">E_A_ref\nACDE\n"))
    stamp = target.stat().st_mtime_ns
    before = backend.input_fingerprint([str(target)])
    # Same byte count and the very same mtime: all a (size, mtime) key can see. Only a file that
    # has been still for a while may be trusted from cache, and this one has not.
    target.write_text(">E_A_ref\nWWWW\n", encoding="utf-8")
    os.utime(target, ns=(stamp, stamp))
    assert backend.input_fingerprint([str(target)]) != before


# ----------------------------------------------------- achado 3: dead import
def test_03_fasta_out_no_longer_carries_a_dead_os_import():
    """Trava o achado 3: `import os` sobreviveu a troca de os.makedirs por staged_outputs."""
    assert not hasattr(fasta_out, "os")
    assert "os.makedirs" not in (fasta_out.__doc__ or "")


# ----------------------------------------------------- achado 4: invariante de create_run
def test_04_create_run_always_writes_its_own_metadata(tmp_path):
    """Fixa a invariante que o ramo morto do achado 4 escondia: run.json nunca e herdado."""
    # A folder bearing the name create_run would pick is stepped over, not reopened, so the
    # metadata it writes always describes the run the caller just asked for.
    squatter = tmp_path / backend.run_dir_name("Teste A")
    squatter.mkdir(parents=True)
    (squatter / backend.RUN_METADATA_FILENAME).write_text(
        json.dumps({"label": "OUTRA ANALISE", "folder": squatter.name}), encoding="utf-8")

    ctx = backend.create_run(tmp_path, "Teste A")

    assert ctx.root != squatter
    meta = json.loads((ctx.root / backend.RUN_METADATA_FILENAME).read_text(encoding="utf-8"))
    assert meta["label"] == "Teste A"
    assert meta["folder"] == ctx.root.name
    assert meta["created_at"]
    assert json.loads((squatter / backend.RUN_METADATA_FILENAME).read_text(
        encoding="utf-8"))["label"] == "OUTRA ANALISE"
