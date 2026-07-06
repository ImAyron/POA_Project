# Guia de testes na WSL (Ubuntu)

Passo a passo para validar o POA de ponta a ponta no WSL2 — incluindo o que **não** pôde ser
verificado no Windows: **POA2 + pyTMHMM**, **EMBOSS antigenic** (substituto do PAP/IMED),
**API MHC-II do IEDB** e (opcional) **BepiPred-3.0**.

> Convenções: comandos em blocos são para colar no terminal da WSL. Onde aparecer
> `>>> esperado:` é o que você deve ver se deu certo. Os blocos `python - <<'PY' … PY`
> resolvem o seu `$HOME` sozinhos (via `expanduser`), então não precisa editar caminhos.

---

## 0. Pré-requisitos

- Windows com **WSL2 + Ubuntu** (`wsl --install -d Ubuntu` no PowerShell, se ainda não tiver).
- **conda/miniconda dentro da WSL** (não é o conda do Windows). Se não tiver:
  ```bash
  cd ~
  wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
  bash Miniconda3-latest-Linux-x86_64.sh -b -p ~/miniconda3
  ~/miniconda3/bin/conda init bash
  exec bash        # recarrega o shell
  ```

> ⚠️ **Não reutilize o `.venv` criado no Windows** (tem binários `.exe`, não roda no Linux).
> Na WSL usamos um **ambiente conda** próprio.

---

## 1. Trazer o código para a WSL

Recomendado: **clonar o repositório local (com a branch da automação) para o home da WSL**.
Evita a lentidão do `/mnt/c`, problemas de fim de linha (CRLF) e o OneDrive.

```bash
git clone -b feature/automation-gui \
  "/mnt/c/Users/pichau/OneDrive/Desktop/POA/POA_Project" ~/POA_Project
cd ~/POA_Project
git log --oneline -3
```
`>>> esperado:` os commits `fix(netctl)…`, `docs(readme)…`, etc.

*(Alternativa, sem copiar: `cd /mnt/c/Users/pichau/OneDrive/Desktop/POA/POA_Project` — funciona,
mas é mais lento e o `.venv` do Windows fica visível; apenas não o use.)*

---

## 2. Criar o ambiente e instalar dependências

```bash
# ambiente (Python 3.11 = ampla compatibilidade com pyTMHMM/biopython)
# usamos --override-channels -c conda-forge para NÃO tocar nos canais 'default' da Anaconda
# (que hoje exigem aceitar Termos de Serviço). Assim evitamos o CondaToSNonInteractiveError.
conda create -n poa --override-channels -c conda-forge python=3.11 pip -y
conda activate poa

# dependências Python do pipeline + GUI
pip install pandas numpy biopython openpyxl pytest requests streamlit plotly

# EMBOSS (substituto do PAP/IMED) — só existe no Linux, via bioconda
conda install -n poa --override-channels -c conda-forge -c bioconda emboss -y

# pyTMHMM (topologia de membrana no POA2) — precisa de compilador C E numpy<2.
# O pyTMHMM 1.3.6 não compila contra numpy 2.x (usa np.int_t) e seu setup.py importa numpy
# em tempo de build (por isso --no-build-isolation, para usar o numpy já instalado no env).
conda install -n poa --override-channels -c conda-forge cython c-compiler -y
pip install "numpy<2"
pip install --no-build-isolation pyTMHMM
python -c "import pyTMHMM; print('pyTMHMM OK')"
```

> ⚠️ Mantenha `numpy<2` neste env. Se um `pip install` futuro subir o numpy para 2.x, o pyTMHMM
> (compilado contra numpy 1.x) quebra com `numpy.core.multiarray failed to import` — basta rodar
> `pip install "numpy<2"` de novo (não precisa recompilar o pyTMHMM). E **não** use
> `--force-reinstall pyTMHMM`, pois isso faz o pip re-subir o numpy para 2.x.

Verifique:
```bash
python -c "import pandas, numpy, Bio, streamlit, plotly, requests; print('py deps OK')"
which antigenic && antigenic -help 2>&1 | head -n 3
python -c "import pyTMHMM; print('pyTMHMM OK')"
```
`>>> esperado:` `py deps OK`, o caminho do `antigenic` + uso, e `pyTMHMM OK`.

*(Se o `pip install pyTMHMM` falhar por causa do numpy, veja **Solução de problemas** no fim.)*

---

## 3. Rodar a suíte de testes automatizados

```bash
cd ~/POA_Project
pytest -q
```
`>>> esperado:` `43 passed`. Valida parsers, ranqueamento/filtragem, conservancy local,
clientes de serviço e o backend da GUI **na sua máquina**.

---

## 4. Teste end-to-end OFFLINE (POA1 → Conservancy local → POA2)

Exercita o pipeline inteiro **sem internet** e valida o **POA2 + pyTMHMM**. Crie o mini-dataset:

```bash
mkdir -p ~/poa_test

# Proteínas (-f): cabeçalho Proteína_Espécie_ID. SARS e MERS idênticas => alta conservância.
cat > ~/poa_test/proteins.fasta <<'FASTA'
>SPIKE_SARS_NP1
MKTAYIAMKGVLMNKQRSTAAILLVVGGAAWLFFIICCMMNNQQRRSSTTVVYY
>SPIKE_MERS_NP2
MKTAYIAMKGVLMNKQRSTAAILLVVGGAAWLFFIICCMMNNQQRRSSTTVVYY
FASTA

# Entrada estilo BepiPred-3.0 (-b3): MAIÚSCULAS = epítopo; devem ser substring da proteína.
cat > ~/poa_test/bp3.fasta <<'FASTA'
>SARS_SPIKE_NP1
mktAYIAMkgvLMNKqrstaailLVVGGaawlffiiccmmnnqqrrssttvvyy
FASTA
```

**4.1 — POA1** (consolida, gera relatório, FASTAs por espécie e .xlsx):
```bash
cd ~/POA_Project
python POA1_v1.0.py -b3 ~/poa_test/bp3.fasta -f ~/poa_test/proteins.fasta -d ~/poa_test/out -e y
ls -R ~/poa_test/out
cat ~/poa_test/out/Analysis_report.txt
cat "$HOME/poa_test/out/Conservancy Analysis/SARS_epitopes.fasta"
```
`>>> esperado:` `Analysis_report.txt`, `Bepipred_Epitopes.xlsx`, e
`Conservancy Analysis/SARS_epitopes.fasta` com cabeçalhos `>SARS_SPIKE_Bepipred3.0_<ini>_<fim>`.

**4.2 — Conservancy local** (substitui a submissão manual ao IEDB):
```bash
cd ~/POA_Project
python - <<'PY'
import os
from poa.services import conservancy_client
base = os.path.expanduser("~/poa_test")
written = conservancy_client.run_conservancy_for_dir(
    os.path.join(base, "out", "Conservancy Analysis"),
    os.path.join(base, "proteins.fasta"),
    threshold=70,
    out_dir=os.path.join(base, "csvs"),
)
print("CSVs gerados:", written)
PY
cat ~/poa_test/csvs/*.csv
```
`>>> esperado:` um `SARS_conservancy.csv` com colunas do IEDB e ~`100.00% (2/2)`.

**4.3 — POA2** (filtra conservância + topologia TMHMM; gera .xlsx e .fasta):
```bash
cd ~/POA_Project
python POA2_v1.0.py -g True -t 70 -d ~/poa_test/csvs -f ~/poa_test/proteins.fasta -r ~/poa_test/poa2 -rf 0
ls ~/poa_test/poa2
```
`>>> esperado:` `POA2_analysis_70.xlsx` (com colunas Portion_Outside/TM/Inside) e
`POA2_analysis.fasta`. Se aparecer `pyTMHMM não instalado`, volte ao passo 2.

---

## 5. Testar as ferramentas externas REAIS

### 5.1 EMBOSS antigenic (substituto do PAP/IMED)
```bash
cd ~/POA_Project
python - <<'PY'
import os
from poa.services import antigenic_client
base = os.path.expanduser("~/poa_test")
fasta = os.path.join(base, "proteins.fasta")
res = antigenic_client.predict(fasta)
with open(os.path.join(base, "pap.txt"), "w") as fh:
    fh.write(antigenic_client.to_pap_txt(res.content, fasta))
print("fonte:", res.source, "-> pap.txt gerado")
PY
cat ~/poa_test/pap.txt
# Alimenta o POA1 pelo caminho PAP/IMED:
python POA1_v1.0.py -p ~/poa_test/pap.txt -f ~/poa_test/proteins.fasta -d ~/poa_test/out_pap
```
`>>> esperado:` `pap.txt` com blocos `>SPIKE_SARS_NP1` + linhas `n<TAB>ini<TAB>peptídeo<TAB>fim`,
e o POA1 rodando sem erro.

### 5.2 API MHC-II do IEDB (precisa de internet) — **confirmar as colunas**
```bash
cd ~/POA_Project
python - <<'PY'
from poa.services import mhcii_client
tsv = mhcii_client.submit("MKTAYIAMKGVLMNKQRSTAAILLVVGG", "HLA-DRB1*01:01", length=15, method="nn_align")
linhas = tsv.splitlines()
print("=== CABEÇALHO REAL DA API ===")
print(linhas[0])
print("=== 1a linha de dados ===")
print(linhas[1] if len(linhas) > 1 else "(sem dados)")
PY
```
`>>> esperado:` uma linha de cabeçalho separada por TAB. **Compare** os nomes de colunas com
`REQUIRED_COLUMNS` em `poa/services/mhcii_client.py` (`allele, seq_num, start, end, method,
peptide, smm_align_ic50, nn_align_ic50, nn_align_rank, nn_align_adjusted_rank`). Se algum nome
diferir (ex.: `ic50` em vez de `nn_align_ic50`), me mande o cabeçalho real que eu ajusto o
mapeamento em `_ALIASES` — o normalizador já cobre variações comuns, mas convém confirmar.

### 5.3 BepiPred-3.0 local (opcional, pesado — baixa pesos ESM-2 ~GB)
```bash
pip install bp3
cd ~/POA_Project
python - <<'PY'
import os
from poa.services import bepipred_client
base = os.path.expanduser("~/poa_test")
res = bepipred_client.predict(os.path.join(base, "proteins.fasta"))
with open(os.path.join(base, "bp3_out.fasta"), "w") as fh:
    fh.write(res.content)
print("fonte:", res.source)
PY
python POA1_v1.0.py -b3 ~/poa_test/bp3_out.fasta -f ~/poa_test/proteins.fasta -d ~/poa_test/out_bp3
```

### 5.4 Validar a Conservancy local contra o IEDB (recomendado antes de produção)
1. Pegue um `*_epitopes.fasta` de `~/poa_test/out/Conservancy Analysis/`.
2. Submeta manualmente em https://tools.iedb.org/conservancy/ com o mesmo conjunto de proteínas
   e o mesmo limiar (70), baixe o CSV.
3. Compare os valores (Minimum/Maximum identity e "% of protein sequence matches") com o CSV em
   `~/poa_test/csvs/`. Devem bater. Se divergirem, me mande os dois CSVs que eu ajusto o cálculo.

---

## 6. Rodar a interface gráfica (Streamlit)

```bash
cd ~/POA_Project
streamlit run poa/gui/app.py
```
Abra **http://localhost:8501** no navegador do **Windows** (o WSL2 encaminha o localhost).
Percorra as 6 etapas: envie `proteins.fasta`, automatize/faça upload das predições, rode o POA1
(veja tabela + gráficos), calcule a Conservancy local, rode o POA2 e baixe os resultados.
Para parar: `Ctrl+C` no terminal.

---

## 7. Solução de problemas

- **`CondaToSNonInteractiveError` (Terms of Service dos canais da Anaconda)** → use
  `--override-channels -c conda-forge` (como nos comandos acima) para não usar os canais
  `pkgs/main`/`pkgs/r`. Alternativa: aceitar uma vez com
  `conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main` (e idem
  para `.../pkgs/r`).
- **`antigenic: command not found`** → EMBOSS não instalou ou o env não está ativo. Rode
  `conda activate poa` e reinstale: `conda install -n poa --override-channels -c conda-forge -c bioconda emboss -y`.
- **`pip install pyTMHMM` falha** — três sintomas comuns e a receita que funciona:
  1. `No module named 'numpy'` no build → o setup.py importa numpy; use `--no-build-isolation`.
  2. `hmm.pyx: np.int_t ... Invalid type` → numpy 2.x removeu `np.int_t`; use `numpy<2`.
  3. `numpy.core.multiarray failed to import` ao importar → foi compilado com numpy 1.x mas o
     runtime está com numpy 2.x → volte o numpy para 1.x.
  Sequência que resolve tudo (com `cython`+`c-compiler` já instalados no env):
  ```bash
  pip install "numpy<2"
  pip install --no-build-isolation pyTMHMM   # NÃO use --force-reinstall (re-sobe o numpy)
  python -c "import pyTMHMM; print('pyTMHMM OK')"
  ```
  O pandas funciona com `numpy>=1.26,<2`.
- **`ModuleNotFoundError: poa`** → rode a partir da raiz do repo (`cd ~/POA_Project`);
  o `pytest.ini` já define `pythonpath = .`.
- **Fim de linha / permissão estranhos** → você está em `/mnt/c` (OneDrive). Prefira o clone em
  `~/POA_Project` (passo 1).
- **API MHC-II sem resposta** → pode ser instabilidade do IEDB; o cliente levanta
  `ServiceUnavailable` e a GUI oferece o upload manual como plano B.

---

## Resumo do que cada teste cobre

| Passo | Valida |
|---|---|
| 3 | Parsers, ranqueamento/filtragem, conservancy local, clientes, backend GUI (43 testes) |
| 4 | Pipeline completo offline: POA1 → Conservancy local → **POA2 + pyTMHMM** |
| 5.1 | **EMBOSS antigenic** → PAP/IMED → POA1 |
| 5.2 | **API MHC-II do IEDB** + conferência das colunas do TSV |
| 5.3 | **BepiPred-3.0** local (bp3) → POA1 |
| 5.4 | Equivalência da **Conservancy local vs IEDB** |
| 6 | **GUI** ponta a ponta |
