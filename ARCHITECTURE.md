# ARQUITETURA & PLANO DE AUTOMAÇÃO — POA (Pipeline de Otimização de Antígenos)

> Documento da **Fase 0 (Exploração)**. Descreve o pipeline atual *exatamente como está*,
> marca cada ponto de intervenção manual, apresenta a investigação da disponibilidade
> **atual (2026)** de cada ferramenta externa e propõe um plano de automação + GUI.
> **Nada de código de produção foi alterado ainda.** A lógica científica de
> ranqueamento/seleção está documentada aqui e será **preservada** na refatoração.

Data da análise: 2026-07-03 · Base: commit `62cf8f4` (branch `main`).

---

## 1. Visão geral

O POA é um pipeline **semiautomático** de vacinologia reversa em duas etapas:

- **POA1** (`POA1_v1.0.py`): lê resultados de predição de epítopos de várias ferramentas web
  (formatos heterogêneos: `.json`, `.txt`, `.html`, `.fasta`), normaliza tudo numa matriz
  única, aplica filtros por método, gera um relatório estatístico e **exporta FASTAs de
  epítopos agrupados por espécie** para a etapa seguinte (+ `.xlsx` opcional).
- **POA2** (`POA2_v1.0.py`): recebe os CSVs da **Epitope Conservancy Analysis (IEDB)**, filtra
  por limiares de conservação e usa **pyTMHMM** para classificar a topologia de membrana de
  cada epítopo (externo / transmembranar / interno); gera `.xlsx` e `.fasta`.

Entre POA1 e POA2 existe uma etapa **inteiramente manual**: submeter cada FASTA de POA1 à
ferramenta web de Conservancy Analysis e baixar os CSVs.

---

## 2. Fluxo de dados exato

```
                          ┌─────────────────── ETAPA MANUAL (web) ───────────────────┐
                          │  Bepipred · PAP/IMED · NetCTL · MHC-II BP                 │
  proteins.fasta ───────► │  (usuário submete, baixa, renomeia, formata resultados)  │
                          └───────────────┬──────────────────────────────────────────┘
                                          │  b2.json / b3.fasta / pap.txt / netctl.html / dir de *.html / x.fasta
                                          ▼
   ┌──────────────────────────────── POA1_v1.0.py (main) ───────────────────────────────┐
   │  checkIntegrity(-f)  → rejeita 'x' nas sequências                                    │
   │  SeqIO → dict {ID: len}  (para checar viabilidade)                                   │
   │  PrEpiAn.runningPrEpiAn(args):                                                       │
   │     bepipred.py  ─┐                                                                  │
   │     papImed.py    ├─► cada parser → DataFrame padronizado ─► pd.concat(join=inner)   │
   │     netctl.py     │      colunas: Method, Specie, Protein, ID_Sequence,             │
   │     mhcii.py      │               Initial Position, Final Position, Peptide Sequence │
   │     othersPred.py ─┘                                                                 │
   │  checkViability()  → nenhum epítopo maior que a menor proteína                       │
   │  writereport()     → {-d}/Analysis_report.txt                                        │
   │  prepianResultsforConservancyAnalysis.filesforEptConsAnalysis():                     │
   │      cria {-d}/Conservancy Analysis/ e escreve {Specie}_epitopes.fasta por espécie   │
   │      header: >{Specie}_{Protein}_{Method}_{Init}_{Final}                             │
   │  (opcional -e y) → *.xlsx por método                                                 │
   └───────────────────────────────────────────┬─────────────────────────────────────────┘
                                                │  {Specie}_epitopes.fasta (1 por espécie)
                                                ▼
                          ┌──────────── ETAPA MANUAL (web IEDB) ────────────┐
                          │  Epitope Conservancy Analysis                   │
                          │  usuário submete cada FASTA + conjunto de       │
                          │  proteínas comparadoras, baixa CSVs → 1 pasta   │
                          └───────────────┬─────────────────────────────────┘
                                          │  *.csv
                                          ▼
   ┌──────────────────────────────── POA2_v1.0.py (main) ───────────────────────────────┐
   │  conservancyAnalysis.EpitConservAnalysis(t, símbolo, m, imax, imin, -d):            │
   │      lê/concatena CSVs, limpa '%', filtra por Min/Max identity e % de matches       │
   │  PepiTMHMM.tmhmmAnalysis(args, df):                                                  │
   │      TMHMM.pyTMHMMpredict(-f) → topologia por proteína                              │
   │      casa cada epítopo à proteína (por substring) → % de resíduos O/M/I             │
   │  df.to_excel(POA2_analysis_{t}.xlsx)                                                 │
   │  (opcional -rf) getfastafile() → FASTA por região de membrana                       │
   └─────────────────────────────────────────────────────────────────────────────────────┘
```

### Módulos e responsabilidades

| Arquivo | Papel | Camada-alvo (§9) |
|---|---|---|
| `POA1_v1.0.py` | CLI da etapa 1 + relatório + orquestração | interface/CLI |
| `POA2_v1.0.py` | CLI da etapa 2 + escrita xlsx/fasta | interface/CLI |
| `PrEpiAn.py` | orquestra os parsers, consolida DataFrame, exporta xlsx | core (lógica) |
| `bepipred.py` | parser Bepipred 2.0 (JSON) e 3.0 (FASTA) + filtro de tamanho | core (parser) |
| `papImed.py` | parser PAP/IMED (TXT copiado da web) | core (parser) |
| `netctl.py` | parser NetCTL 1.2 (HTML) | core (parser) |
| `mhcii.py` | parser MHC-II BP (pasta de HTMLs) + filtro IC50/HLA | core (parser) |
| `othersPred.py` | parser genérico (FASTA padronizado) | core (parser) |
| `prepianResultsforConservancyAnalysis.py` | escreve FASTAs por espécie p/ Conservancy | core (saída) |
| `conservancyAnalysis.py` | parser + filtro dos CSVs da Conservancy Analysis | core (parser/filtro) |
| `PepiTMHMM.py` | casa epítopo↔proteína e calcula topologia | core (lógica) |
| `TMHMM.py` | wrapper do `pyTMHMM` | services (externo local) |

---

## 3. Argumentos de linha de comando (contrato a preservar)

### POA1 (`python POA1_v1.0.py ...`)
Obrigatórios: **`-f`** (FASTA de proteínas), **`-d`** (pasta de saída) e **pelo menos um** método.

| Flag | Significado | Formato | Default |
|---|---|---|---|
| `-b2` | Bepipred 2.0 | `.json` (JSON Summary) | `''` |
| `-b3` | Bepipred 3.0 | `.fasta` (maiúsculas = epítopo) | `''` |
| `-bmin`/`-bmax` | tamanho min/max Bepipred | int | 0/0 |
| `-p` | PAP/IMED | `.txt` (tabela copiada) | `''` |
| `-pmin`/`-pmax` | tamanho min/max PAP | int | 0/0 |
| `-n` | NetCTL 1.2 | `.html` | `''` |
| `-m` | MHC-II BP | **diretório** de `.html` | `''` |
| `-mhla` | classe HLA (DR/DP/DQ) | str | `DR` |
| `-mic` | limiar IC50 (nM) | int | 50 |
| `-x` | outros preditores | `.fasta` padronizado | `''` |
| `-xmin`/`-xmax` | tamanho min/max outros | int | 0/0 |
| `-e` | exportar `.xlsx` (y/n) | str | `n` |

> **Nota:** `-b2` e `-b3` são mutuamente exclusivos (exceção lançada se ambos).

### POA2 (`python POA2_v1.0.py ...`)
Obrigatórios: **`-d`** (pasta de CSVs), **`-t`** (limiar de identidade), **`-f`** (FASTA de
proteínas), e **`-g` OU `-l`** (grupo mutuamente exclusivo, obrigatório).

| Flag | Significado | Default |
|---|---|---|
| `-g` / `-l` | conservados (`>=`) / únicos (`<`) | — (obrigatório escolher 1) |
| `-t` | limiar de identidade de sequência | (obrigatório) |
| `-d` | pasta com CSVs da Conservancy | (obrigatório) |
| `-f` | FASTA de proteínas (p/ TMHMM) | (obrigatório) |
| `-r` | pasta de saída | `''` (cwd) |
| `-imin`/`-imax` | identidade mín/máx (%) | 60/100 |
| `-m` | % de sequências com match ≥ limiar | 60 |
| `-rf` | FASTA por região: 0=todos 1=externo 2=TM 3=interno | `None` |

---

## 4. Formatos de entrada/saída por parser (contrato interno)

Todos os parsers de POA1 devolvem um DataFrame com as **mesmas 7 colunas**:
`Method, Specie, Protein, ID_Sequence, Initial Position, Final Position, Peptide Sequence`
(MHC-II também produz `Allele`, mas ela é **descartada** na consolidação — ver §6).

- **Bepipred 2.0** (`bp2_*`): lê JSON `antigens[id] = {AA[], PRED[], ...}`; resíduo é
  "Epitope" se `PRED > 0.5`; agrupa resíduos contíguos em regiões. Header do antígeno é
  parseado por regex `(\w+?)_(\w+?)_(\w+)` → Protein, Specie, ID_Sequence.
- **Bepipred 3.0** (`bp3_*`): lê FASTA onde **letras maiúsculas = epítopo**; regex `[A-Z]+`
  extrai segmentos; header `split('_', 2)` → Specie, Protein, ID_Sequence.
- **PAP/IMED** (`PAPepitopes`): lê `.txt` com blocos `>ID` seguidos de linhas
  `n<TAB>start<TAB>sequence<TAB>end`; header `split('_')` → Protein, Specie, ID_Sequence.
- **NetCTL** (`netctlAntigenEpitopes`): lê `.html`, pega linhas que começam com dígito;
  seleciona onde a última coluna é `<-E` (ligante identificado); header `split('_')`
  → [0]=Protein, [1]=Specie.
- **MHC-II** (`mhcii.py`): varre **diretório** de `.html`; nome do arquivo
  `Protein_Specie.html` dá metadados; converte a tabela TSV embutida (linha `allele...`
  = cabeçalho, linhas `H...` = dados); usa **NN_align IC50**; filtra por classe HLA e
  `nn_align_ic50 <= -mic`.
- **Outros** (`fasta_epitopes`): FASTA com header
  `>Protein_Specie_Method_[NP_]ID_Init_Final`.
- **Saída para Conservancy** (`prepianResultsforConservancyAnalysis`): um FASTA por espécie,
  header `>{Specie}_{Protein}_{Method}_{Init}_{Final}` (**exatamente 5 campos** — POA2 depende disso).
- **Conservancy CSV** (`conservancyAnalysis`): espera colunas `Epitope #, Epitope name,
  Epitope sequence, Epitope length, "Percent of protein sequence matches at identity <= 100%",
  Minimum identity, Maximum identity, View details`.
- **TMHMM** (`TMHMM.py` + `PepiTMHMM.py`): `pyTMHMM.predict()` devolve string `o/m/i` por
  resíduo; conta a fração do epítopo em cada região.

---

## 5. Pontos de intervenção MANUAL (o que vamos automatizar)

| # | Ponto manual atual | Esforço | Estratégia de automação (§8) |
|---|---|---|---|
| **M1** | Submeter proteínas ao **Bepipred**, esperar a fila, baixar JSON/FASTA | alto | pacote local `bp3` (3.0) / standalone (2.0); fallback upload |
| **M2** | Submeter ao **PAP/IMED** *uma proteína por vez* e **copiar/colar a tabela** num `.txt` | muito alto (frágil) | **EMBOSS `antigenic` local** (mesmo algoritmo K&T); fallback upload |
| **M3** | Submeter ao **NetCTL 1.2**, esperar carregar, "salvar como" `.html` | alto | automação de navegador (fila assíncrona) ou standalone WSL; fallback upload |
| **M4** | Submeter ao **MHC-II BP** *proteína por proteína*, salvar `.html` com nome exato | muito alto | **API REST do IEDB** (síncrona); fallback upload |
| **M5** | Pegar FASTAs de POA1 → **Conservancy Analysis** (web) → baixar CSVs → 1 pasta | muito alto | **reimplementação local do algoritmo** (recomendado) ou navegador; fallback upload |
| **M6** | Nomear/organizar arquivos e passar POA1→POA2 na mão | médio | orquestração pela GUI (fluxo guiado) |

O **plano B semiautomático** (upload manual do arquivo de resultado) é **preservado para todas
as ferramentas** — se um serviço externo falhar ou estiver indisponível, o usuário faz o upload
e o pipeline segue exatamente como hoje.

---

## 6. Fragilidades e bugs observados (documentar; **não** alterar sem aviso)

Estes pontos afetam a robustez, **não** a lógica científica. Serão tratados só com sua aprovação
e cobertos por testes antes de qualquer mudança:

1. **Coluna `Allele` é perdida** — `PrEpiAn` consolida com `pd.concat(join="inner")`, e como só
   o DataFrame do MHC-II tem `Allele`, ela some do resultado final. Se o alelo importa no
   relatório, isso é uma perda silenciosa.
2. **`SB`/`WB` fixos em 50/500** — em `mhcii.MHCIIAntigenEpitopes` a coluna `Prediction_NN` usa
   limiares fixos, independente de `-mic` (o filtro real por `-mic` funciona; a *rotulagem* não).
3. **`PepiTMHMM` pode desalinhar tamanhos** — se um epítopo não casar com nenhuma proteína, as
   listas `Out/TM/Ins` ficam menores que o DataFrame e a atribuição de colunas quebra.
4. **`files_map_MHCIIBD`** usa `os.path.join(directory, file)` em vez de `folders` do `os.walk`
   (falha se houver subpastas). Também assume nome `Protein_Specie.html` estritamente.
5. **Acoplamento por nome** — o casamento epítopo↔proteína em POA2 depende de `Specie` ser
   substring do header da proteína; nomes inconsistentes silenciam resultados (só `warnings.warn`).
6. **Paths com `/` hardcoded** e `os.mkdir` (não `makedirs`) — frágil no Windows e falha se a
   pasta pai não existir.
7. **Sem `requirements.txt`, sem testes, sem logging** — adicionaremos.

---

## 7. Disponibilidade ATUAL das ferramentas externas (investigação 2026)

> Verificado por busca web + fetch das páginas/repos reais em 2026-07. **Não** assumimos que
> endpoints antigos funcionam. Legenda de estratégia: 🟢 API/pacote programável · 🟡 automação
> de navegador · 🔴 sem automação viável (só upload manual).

### 7.1 Bepipred 2.0 / 3.0 — DTU Health Tech
- **Web**: vivos em `https://services.healthtech.dtu.dk/services/BepiPred-2.0/` e `.../BepiPred-3.0/`
  (os URLs antigos `service.php?...` redirecionam). **Sem API REST**; fila assíncrona com Job ID.
- **BepiPred-3.0 🟢**: distribuído como **código**. `pip install bp3` (v0.0.12.7, 2024) **ou** clone
  `github.com/UberClifford/BepiPred-3.0` → CLI `python bepipred3_CLI.py -i in.fasta -o out/ -pred vt_pred`.
  Requer **PyTorch + ESM-2** (baixa pesos no 1º uso; pesado, ~GB). Saída inclui **FASTA com
  maiúsculas = epítopo** → *diretamente compatível* com o parser `bp3_FastaAnalysis`. Também há
  `pybiolib` (`biolib run DTU/BepiPred-3`).
- **BepiPred-2.0 🟡**: só via **standalone DTU** (v2.0c) atrás de cadastro acadêmico (exige
  e-mail institucional). Sem PyPI. Saída JSON só pela web.
- **Licença**: acadêmica/não-comercial pela DTU (o PyPI `bp3` declara MIT — **conflito a
  esclarecer** antes de uso comercial). Sem CAPTCHA; `robots.txt` = 404.

### 7.2 NetCTL 1.2 — DTU Health Tech
- **Web 🟡**: vivo em `https://services.healthtech.dtu.dk/services/NetCTL-1.2/`. **Sem API,
  sem PyPI, sem binário Windows**. Standalone só **Linux/IRIX/Darwin** (tcsh, chama netMHC/netChop),
  cadastro acadêmico. Fluxo web = fila assíncrona (POST → Job ID → poll). Sem CAPTCHA.
- **Consequência**: no Windows, a automação local exige **WSL/Docker**; caso contrário, a via
  realista é **automação de navegador** (Playwright) — ou manter **upload manual**.

### 7.3 MHC-II Binding Predictions — IEDB
- **Web**: `https://tools.iedb.org/mhcii/` vivo (com banner de migração p/ o "next-gen").
- **API REST legada 🟢 (confirmada viva)**: `https://tools-cluster-interface.iedb.org/tools_api/mhcii/`
  — `POST` **síncrono**, corpo `method=...&sequence_text=...&allele=...&length=...`, resposta
  **TSV** com colunas `allele, seq_num, start, end, length, core_peptide, peptide` + colunas de
  score por método (inclui `nn_align_ic50`, `nn_align_rank`). Sem login/chave. (GET devolve 405 →
  endpoint existe, é POST-only.)
- **API next-gen 🟢**: `POST https://api-nextgen-tools.iedb.org/api/v1/pipeline` (assíncrona,
  poll `results_uri`) — mais durável/futuro.
- **Standalone 🟢**: `IEDB_MHC_II-3.1.12.tar.gz` (~1.5 GB, Linux). 
- **Impacto no POA**: a API devolve **TSV** com as mesmas colunas que o parser HTML já lê. Basta
  um adaptador TSV (ou salvar a resposta e reusar a lógica de `MHCIIHTMLconverter`).

### 7.4 Epitope Conservancy Analysis — IEDB  ⚠️ ponto crítico
- **Web 🟡**: `https://tools.iedb.org/conservancy/` vivo, **sem banner de migração**.
- **SEM API REST** (não está na lista de `tools_api/`) e **SEM standalone** (a própria página de
  download declara que não há versão baixável). **Confirmado.**
- **Algoritmo é simples** (Bui et al., *BMC Bioinformatics* 2007): para cada epítopo, desliza-o
  contra cada proteína e calcula a **identidade máxima (%)**; a conservância é a fração de
  proteínas com match ≥ limiar. → **Reimplementação local viável (~algumas dezenas de linhas)**,
  o que **remove o maior gargalo manual (M5)**. *Isto toca a "lógica científica" — precisa da sua
  aprovação e será validado contra saídas reais do IEDB.*

### 7.5 PAP/IMED (antigenicidade)  ⚠️ tool original fora do ar
- **Original 🔴**: `http://imed.med.ucm.es/Tools/antigenic.pl` responde, mas devolve **página de
  acesso restrito da UCM** (só via VPN institucional). Não é utilizável publicamente.
- **Substituto local 🟢**: o PAP é um front-end do programa **EMBOSS `antigenic`** (mesmo método
  Kolaskar & Tongaonkar, 1990). Instalável via **Bioconda** (`conda install -c bioconda emboss`,
  EMBOSS 6.6.0). No Windows: **WSL2 + conda** (ou Docker). Uso:
  `antigenic -sequence prot.fasta -minlen 6 -rformat gff -outfile out.gff -auto` (saída estruturada:
  start/end/score/sequência). Sem wrapper Biopython → chamar por `subprocess`.
- **EBI REST**: o Job Dispatcher da EBI existe, mas **não** oferece `antigenic` → via morta.

### 7.6 pyTMHMM (topologia de membrana) — já local
- Usado por POA2 via `TMHMM.py`. Continua **local** (`pip install pyTMHMM`), sem dependência web.
  Manter; apenas isolar atrás da camada de serviços e tratar erros de instalação.

### Resumo executivo por ferramenta

| Ferramenta | API pública? | Pacote local? | Estratégia primária | Fallback |
|---|---|---|---|---|
| Bepipred 3.0 | não | **sim** (`bp3`/GitHub) | 🟢 pacote local | upload |
| Bepipred 2.0 | não | standalone (cadastro) | 🟡 navegador / manual | upload |
| NetCTL 1.2 | não | Linux-only | 🟡 navegador (ou WSL) | upload |
| MHC-II BP | **sim (REST)** | Linux 1.5 GB | 🟢 **API REST IEDB** | upload |
| Conservancy | **não** | **não** | 🟢 **reimpl. local** (a aprovar) | navegador / upload |
| PAP/IMED | fora do ar | **EMBOSS antigenic** | 🟢 EMBOSS local (WSL) | upload |
| pyTMHMM | — | **sim** | 🟢 já local | — |

---

## 8. Plano de automação por ferramenta

Padrão comum a **todas** as integrações (camada `services/`, §9):

```
predict(input_fasta, params) -> ResultArtifact
  1. calcula chave de cache = hash(ferramenta + versão + input + params)
  2. se houver no cache → retorna sem resubmeter
  3. tenta a estratégia primária (pacote local / API / navegador)
  4. em falha/indisponibilidade → loga erro claro e sinaliza FALLBACK
  5. FALLBACK = usuário fornece o arquivo de resultado (comportamento atual)
  6. em qualquer caminho, o resultado passa pelo MESMO parser de core/
```

- **Bepipred 3.0** → `bp3` local; saída FASTA já compatível com o parser existente.
- **Bepipred 2.0** → manter parser JSON; automação por navegador é opcional (baixo ROI); default = upload.
- **NetCTL 1.2** → cliente de navegador (Playwright) para a fila assíncrona; default seguro = upload.
- **MHC-II** → cliente da **API REST IEDB** (POST síncrono); adaptador TSV→DataFrame reusando a
  lógica de colunas atual. Fallback upload.
- **Conservancy** → **reimplementação local** do cálculo de identidade (a aprovar); produz um CSV
  no **mesmo esquema de colunas** que POA2 já consome, garantindo compatibilidade. Alternativa:
  navegador. Fallback upload dos CSVs.
- **PAP** → EMBOSS `antigenic` via `subprocess`; adaptador da saída para o esquema de `papImed`.
  Fallback upload do `.txt`.
- **pyTMHMM** → manter; isolar e tratar `ImportError`/falhas com mensagem clara.

---

## 9. Arquitetura-alvo (3 camadas, CLI preservada)

```
poa/
├── core/            # (a) parsing + LÓGICA do pipeline — puro, sem rede. Reutilizável por CLI e GUI.
│   ├── models.py            # dataclasses/constantes das 7 colunas padrão
│   ├── parsers/             # bepipred, papimed, netctl, mhcii, others, conservancy  (extraídos dos .py atuais)
│   ├── ranking.py           # filtros por método (mesma lógica científica de hoje)
│   ├── report.py            # writereport()
│   ├── conservancy.py       # filtro dos CSVs (atual) + (novo) cálculo local opcional
│   ├── topology.py          # PepiTMHMM/TMHMM (lógica de casamento e % O/M/I)
│   └── pipeline.py          # run_poa1(config)->artefatos ; run_poa2(config)->artefatos
├── services/        # (b) INTEGRAÇÃO com serviços externos (rede / subprocess / pacotes)
│   ├── base.py              # interface Predictor + cache + fallback + logging
│   ├── bepipred_client.py   # bp3 local
│   ├── netctl_client.py     # navegador (Playwright)
│   ├── mhcii_client.py      # API REST IEDB
│   ├── conservancy_client.py# local/navegador
│   ├── antigenic_client.py  # EMBOSS antigenic (subprocess)
│   └── tmhmm_client.py      # pyTMHMM
├── cache/           # resultados em disco, indexados por hash (não re-submeter)
├── cli/             # POA1/POA2 como hoje (argumentos idênticos)  → chamam core.pipeline
├── gui/             # app Streamlit (§11) — só orquestra core + services
└── tests/           # pytest (§12)

# Compat: POA1_v1.0.py e POA2_v1.0.py permanecem como wrappers finos → cli/, sem quebrar uso atual.
```

**Princípios**: (1) `core/` não importa nada de rede; (2) a CLI continua funcionando com os
mesmos argumentos; (3) a GUI e a CLI compartilham exatamente o mesmo `core`; (4) toda ferramenta
externa tem fallback manual.

---

## 10. Erros, logging e cache

- **Logging** (`logging` stdlib): um logger por camada; níveis INFO (progresso) / WARNING
  (fallback acionado) / ERROR (falha). A GUI mostra o log por etapa.
- **Erros claros**: cada cliente externo captura timeout/HTTP/indisponibilidade e devolve uma
  mensagem acionável + oferta de fallback manual. Nunca "trava" o pipeline.
- **Cache**: `cache/<ferramenta>/<hash>.json|tsv|fasta`. Chave = SHA-256 de
  (nome+versão da ferramenta, parâmetros, conteúdo do input). Botão "ignorar cache" na GUI.

---

## 11. Interface gráfica — recomendação: **Streamlit** ✅

**Recomendo Streamlit** (web local), não PyQt/Tkinter. Justificativa:

| Critério | Streamlit | PyQt/Tkinter |
|---|---|---|
| Esforço p/ o fluxo (upload, tabelas, download, gráficos) | baixo (widgets prontos) | alto (layout + event loop manual) |
| Manutenção (seu requisito explícito) | **simples** (Python puro, 1 arquivo por página) | maior (sinais/slots, threading p/ não travar UI) |
| Tabelas e gráficos científicos | `st.dataframe`, `st.line_chart`, Plotly nativo | precisa integrar matplotlib/Qt manualmente |
| Progresso/status por etapa | `st.status`, `st.progress` nativos | implementar na mão |
| Multiplataforma / zero instalação de GUI | roda no navegador local | depende de Qt/Tk no SO |
| Empacotar como desktop (se quiser depois) | possível (stlite/electron/`streamlit` + atalho) | nativo |

Contras do Streamlit: reexecuta o script a cada interação (contornável com `st.session_state` e
cache) e não é um "app nativo". Para um **fluxo linear guiado** com uploads, tabelas e downloads —
exatamente o caso do POA — os prós superam. Tarefas longas (predições/navegador) rodam em
*background* com atualização de status para não bloquear a UI.

**Fluxo da GUI (6 passos, conforme pedido):**
1. Upload do FASTA de proteínas + parâmetros de POA1 (tamanhos, HLA, IC50, etc.).
2. Disparo/monitoramento das predições por ferramenta (barra de status) **ou upload manual** (fallback).
3. Executar POA1 + visualizar o relatório estatístico e a tabela de epítopos.
4. Submeter à Conservancy Analysis (local/navegador) **ou upload dos CSVs**.
5. Executar POA2 com parâmetros de conservação/topologia.
6. Visualizar e baixar resultados finais (`.xlsx`/`.fasta`) com tabelas e gráficos dos epítopos.

---

## 12. Testes (pytest)

- **Parsers**: fixtures de exemplos pequenos de cada formato (JSON b2, FASTA b3, TXT PAP,
  HTML NetCTL, HTML MHC-II, FASTA "outros", CSV Conservancy) → asserts nas 7 colunas.
- **Ranqueamento/filtragem**: limites de tamanho (bmin/bmax…), filtro IC50/HLA do MHC-II,
  filtros de identidade/`-m` da Conservancy, símbolo `>=`/`<`.
- **Topologia**: `epitTMHMMcaract` com strings `o/m/i` conhecidas → frações corretas.
- **Conservancy local (se aprovado)**: comparar contra um CSV real baixado do IEDB (teste de
  equivalência) antes de confiar na reimplementação.
- **Regressão de CLI**: rodar POA1/POA2 num dataset mínimo e comparar saídas.

---

## 13. `requirements.txt` planejado (com motivos)

```
pandas            # matrizes de epítopos (já usado)
numpy             # suporte numérico ao parser Bepipred (já usado)
biopython         # SeqIO / parsing FASTA (já usado)
openpyxl          # escrita .xlsx via pandas.to_excel (dependência implícita hoje)
pyTMHMM           # topologia de membrana em POA2 (já usado)
requests          # cliente HTTP p/ API REST do IEDB (MHC-II)
streamlit         # GUI (novo)
plotly            # gráficos dos epítopos na GUI (novo)
playwright        # automação de navegador p/ NetCTL / Conservancy (novo, opcional)
pytest            # testes (dev)
# bp3            — Bepipred 3.0 local (pesado: PyTorch+ESM-2). Instalação opcional/documentada.
# EMBOSS         — não é pip: instalar via conda/WSL/Docker (documentar no README).
```

---

## 14. Roadmap incremental (1 commit por etapa)

1. **Fase 0** — este documento + `requirements.txt` + esqueleto de pastas (sem quebrar nada).
2. **Refatoração core/** — extrair parsers e lógica para `core/`, mantendo `POA1/POA2` funcionando;
   adicionar testes dos parsers (rede nenhuma).
3. **Camada services + cache + logging** — clientes com fallback manual; começar pelo mais seguro
   (MHC-II API) e por EMBOSS/Bepipred local.
4. **Conservancy** — (após sua aprovação) reimplementação local validada contra o IEDB.
5. **NetCTL/Bepipred-2.0** — automação de navegador (opcional) com fallback.
6. **GUI Streamlit** — fluxo de 6 passos, progresso, gráficos, downloads.
7. **Polimento** — testes de regressão de CLI, docs no README, tratamento de bugs da §6 (com aviso).

---

## 15. Pontos de decisão que preciso confirmar antes de implementar

1. **GUI**: confirmo **Streamlit**? (recomendação acima).
2. **Conservancy Analysis**: autoriza a **reimplementação local** do cálculo (validada contra o
   IEDB) — já que **não há API nem standalone**? Sem isso, a automação de M5 fica só por navegador.
   *(Toca a lógica científica → só faço com seu OK.)*
3. **PAP/IMED**: aceita substituir o tool (fora do ar) por **EMBOSS `antigenic` local**? Você tem
   **WSL2/conda ou Docker** disponível no Windows para EMBOSS e (opcionalmente) NetCTL?
4. **Bepipred**: prioriza **3.0 via `bp3` local** (dependências pesadas: PyTorch+ESM-2) como via
   automatizada, mantendo 2.0 como upload manual? Ou quer 2.0 automatizado por navegador também?
5. **NetCTL 1.2**: prefere **automação de navegador** (Playwright) ou apenas manter **upload
   manual** (dado que não há API nem binário Windows)?
6. **Licença/uso**: o uso é **acadêmico/não-comercial**? (as licenças DTU e do standalone MHC-II
   restringem uso comercial — relevante para a estratégia de distribuição).
```
