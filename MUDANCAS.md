# POA v1.1 — Relatório de Mudanças Implementadas

**Projeto:** POA — Pipeline de Otimização de Antígenos
**Branch:** `feature/viz-realdata` (derivada de `feature/automation-gui`, base `main`)
**Data:** 17 de agosto de 2026
**Base de comparação:** commit `62cf8f4` (`main`)

---

## 1. Resumo executivo

O POA era um pipeline semiautomático composto por scripts Python soltos na raiz do repositório,
com todas as submissões às ferramentas web feitas manualmente pelo usuário. Esta entrega
reorganiza o código em um pacote Python em camadas, adiciona uma camada de automação para as
ferramentas externas, uma interface gráfica web local, visualização 2D/3D dos epítopos e uma
suíte de testes — **sem alterar a lógica científica de ranqueamento e seleção de epítopos** e
**sem quebrar a linha de comando original**.

| Indicador | Antes (`main`) | Depois (`feature/viz-realdata`) |
|---|---|---|
| Organização do código | 13 scripts na raiz | pacote `poa/` em 4 camadas (44 arquivos) |
| Testes automatizados | nenhum | 82 testes (`pytest`), todos passando |
| Interface | apenas linha de comando | CLI (inalterada) + GUI Streamlit de 7 etapas |
| Submissão às ferramentas web | 100% manual | automatizada ou reimplementada, com fallback manual |
| Análise de conservação | manual no site do IEDB | reimplementação local + upload manual |
| Visualização dos resultados | planilhas e FASTAs | mapa 2D de epítopos + estrutura 3D interativa |
| Documentação | README | README + `ARCHITECTURE.md` + `TESTING_WSL.md` + este relatório |

Total consolidado do diff `main → feature/viz-realdata`: **63 arquivos alterados, 4.934 linhas
adicionadas, 1.326 removidas** (as remoções são os scripts antigos migrados para o pacote).

---

## 2. Reorganização do código em pacote (`refactor`)

Os scripts da raiz foram convertidos em um pacote com separação clara de responsabilidades. Os
arquivos originais `POA1_v1.0.py` e `POA2_v1.0.py` permaneceram como pontos de entrada finos
(*shims*) que apenas chamam o pacote, de modo que todos os comandos documentados no README
continuam funcionando exatamente como antes.

| Camada | Pacote | Responsabilidade |
|---|---|---|
| Núcleo | `poa/core` | Parsing dos formatos de predição, ranqueamento, relatório, topologia. Sem rede, sem I/O externo. |
| Serviços | `poa/services` | Comunicação com ferramentas externas (API, subprocesso, navegador), cache e fallback manual. |
| Linha de comando | `poa/cli` | `poa1`, `poa2` e o novo `import_realdata`. Mesmos argumentos de antes. |
| Interface | `poa/gui` | Aplicação Streamlit e helpers de visualização. |

**Mapeamento dos scripts antigos:**

| Script original (removido) | Novo local |
|---|---|
| `PrEpiAn.py` | `poa/core/pipeline.py` (`consolidate_predictions`) |
| `bepipred.py` | `poa/core/parsers/bepipred.py` |
| `papImed.py` | `poa/core/parsers/papimed.py` |
| `netctl.py` | `poa/core/parsers/netctl.py` |
| `mhcii.py` | `poa/core/parsers/mhcii.py` |
| `othersPred.py` | `poa/core/parsers/others.py` |
| `conservancyAnalysis.py` | `poa/core/parsers/conservancy.py` |
| `TMHMM.py`, `PepiTMHMM.py` | `poa/core/topology.py` |
| `prepianResultsforConservancyAnalysis.py` | `poa/core/fasta_out.py` |

Além disso, as chamadas `print` de diagnóstico foram substituídas por *logging* configurável
(`poa/logging_conf.py`), e as colunas padronizadas do pipeline foram centralizadas em
`poa/core/models.py`.

---

## 3. Camada de automação das ferramentas externas

Cada ferramenta externa ganhou um cliente próprio em `poa/services`. Todos seguem o mesmo
contrato: quando o serviço não está disponível, lançam `ServiceUnavailable`, o pipeline registra
uma mensagem clara e oferece o **upload manual do arquivo de resultado** — preservando
integralmente o comportamento semiautomático original.

| Ferramenta | Estratégia implementada | Fallback |
|---|---|---|
| MHC-II Binding (IEDB) | API REST legada síncrona (`tools-cluster-interface.iedb.org/tools_api/mhcii/`) | upload manual |
| BepiPred-3.0 | pacote local `bp3` (PyTorch + ESM-2), saída FASTA já compatível com o parser | upload manual |
| PAP/IMED | EMBOSS `antigenic` local via subprocesso (mesmo método Kolaskar–Tongaonkar do site, que saiu do ar) | upload manual |
| Epitope Conservancy (IEDB) | **reimplementação local** do algoritmo | upload manual do CSV |
| NetCTL 1.2 / BepiPred-2.0 | automação de navegador (Playwright), *best effort* | upload manual (padrão na GUI) |
| pyTMHMM | execução local, inalterada | — |

Pontos de projeto relevantes:

* **Saída normalizada.** Cada cliente grava o resultado no mesmo formato que o parser
  correspondente já consumia. O cliente MHC-II, por exemplo, converte o TSV da API para o layout
  `<Proteína>_<Espécie>.html`, e o cliente EMBOSS converte o GFF para o layout `.txt` do PAP/IMED.
  Isso significa que **o caminho automático e o caminho manual usam exatamente o mesmo parser** —
  não há divergência de comportamento entre eles.
* **Cache em disco** (`poa/services/cache.py`): chave SHA-256 derivada de ferramenta, versão,
  parâmetros e conteúdo de entrada. Evita ressubmeter as mesmas sequências. Os arquivos ficam em
  `.poa_cache/<ferramenta>/<chave>.<ext>`, legíveis e fáceis de limpar.
* **Automação de navegador marcada como experimental.** NetCTL 1.2 e BepiPred-2.0 (DTU) só
  expõem formulários web assíncronos (submeter → id do job → *polling* → resultado), cujo HTML
  pode mudar sem aviso. Por isso a GUI mantém essas duas ferramentas em modo manual por padrão.

### 3.1 Reimplementação local da Epitope Conservancy Analysis

Verificamos que o IEDB **não oferece API nem versão standalone** para essa ferramenta, o que a
tornava o principal gargalo manual entre POA1 e POA2. Com aprovação explícita do usuário — por
tocar em lógica científica —, o algoritmo foi reimplementado localmente em
`poa/services/conservancy_client.py`, seguindo Bui et al. (BMC Bioinformatics, 2007):

> Para um epítopo linear de comprimento L, desliza-se uma janela sem *gaps* de comprimento L ao
> longo de cada proteína, calculando a identidade percentual (correspondências / L × 100) em cada
> deslocamento; a identidade do epítopo com aquela proteína é o **máximo** entre todos os
> deslocamentos. Dado um limiar de identidade T, reporta-se por epítopo o percentual (e a
> contagem) de proteínas com identidade ≥ T, além da identidade mínima e máxima observadas.

O CSV gerado usa exatamente os mesmos cabeçalhos de coluna que o parser do POA2 já consumia,
portanto **o POA2 funciona de forma idêntica** com CSVs vindos do IEDB ou da reimplementação.

> **Recomendação de validação:** antes de usar em produção, compare a saída local com uma
> exportação real do IEDB para o seu conjunto de dados. As verificações de equivalência estão em
> `tests/test_conservancy_local.py`.

---

## 4. Interface gráfica (Streamlit)

`poa/gui/app.py` implementa um fluxo guiado que orquestra o pipeline inteiro em **7 etapas**:

1. **Proteínas** — upload do FASTA de referência (`-f`).
2. **Predições** — dispara cada ferramenta automatizada ou aceita upload manual do resultado.
3. **POA1** — consolidação, ranqueamento e geração dos FASTAs por espécie.
4. **Conservancy Analysis** — execução local ou upload dos CSVs do IEDB.
5. **POA2** — filtro por conservação e classificação de topologia de membrana.
6. **Resultados** — tabelas, estatísticas e download dos `.xlsx` / `.fasta`.
7. **Visualização (2D/3D)** — mapa de epítopos e estrutura tridimensional (ver seção 5).

A orquestração ficou em `poa/gui/backend.py`, propositalmente **sem dependência do Streamlit**,
o que permite testá-la sem subir a interface (`tests/test_gui_backend.py`).

Execução:

```bash
streamlit run poa/gui/app.py      # ou:  python run_gui.py
```

---

## 5. Visualização 2D e 3D dos epítopos *(novo nesta branch)*

Nova etapa da GUI, apoiada em `poa/gui/viz.py` (também livre de imports do Streamlit, portanto
testável — `tests/test_viz.py`).

**Mapa 2D de epítopos (Plotly).** Cada epítopo é desenhado como um segmento horizontal sobre o
eixo de resíduos, com uma trilha por `Proteína_Espécie` e cor por método de predição. Permite ver
de relance onde os métodos concordam, quais regiões concentram epítopos e qual a extensão de cada
um. Não exige estrutura tridimensional — funciona assim que o POA1 roda.

**Visualizador 3D (py3Dmol / 3Dmol.js).** O usuário envia um arquivo `.pdb` e escolhe de qual
`Proteína_Espécie` destacar os epítopos; eles são realçados em vermelho sobre a estrutura, em
representação *cartoon*, *stick* ou *sphere*, com superfície de van der Waals opcional — no
estilo de visualização do Discovery Studio.

Duas ressalvas foram tratadas explicitamente na interface:

* **Numeração de resíduos.** A numeração do PDB precisa corresponder à posição na sequência usada
  nas predições. A interface exibe esse alerta junto da renderização, já que *gaps* ou
  deslocamentos de numeração produziriam destaques incorretos sem nenhum erro aparente.
* **Dependência de rede.** O 3Dmol.js é carregado de um CDN, então a visualização 3D exige
  internet no momento da renderização. Falhas são exibidas com `st.error`, sem derrubar a
  aplicação.

Nova dependência: `py3Dmol>=2.0` (adicionada ao `requirements.txt`).

---

## 6. Adaptador para os arquivos originais da pesquisa *(novo nesta branch)*

Ao testar o pipeline com os arquivos reais da pesquisa (`2_fastas_mundo/`,
`5_predicao_de_epitopo/`), constatamos que eles **não seguem as convenções de cabeçalho** que os
parsers esperam. O problema central: as exportações JSON do BepiPred-2.0 trazem uma única chave
de antígeno genérica, `"Sequence"`. O parser extrai proteína/espécie/ID dessa chave via expressão
regular `(\w+?)_(\w+?)_(\w+)`, de modo que uma chave genérica produzia `NaN` — e os arquivos de
conservação saíam nomeados `nan_epitopes.fasta`.

`poa/services/realdata_import.py` resolve isso adaptando **apenas a nomenclatura, nunca os dados
preditos**:

* Reescreve a chave do antígeno para o formato `Proteína_Espécie_ID` (ex.: `E_DENV1_ref`),
  copiando os escores por resíduo literalmente.
* Reconstrói o FASTA de referência (`-f`) a partir do próprio array `AA` do JSON — ou seja, a
  sequência exata sobre a qual o BepiPred rodou, eliminando qualquer risco de descompasso entre a
  referência e a predição.
* `best_world_match` mapeia automaticamente cada vírus ao seu conjunto de diversidade para a
  etapa de conservação, por identidade com o primeiro registro. Empates são desfeitos pelo
  comprimento mais próximo, o que faz o CHIKV escolher o conjunto de E1 em vez da poliproteína
  E2-E1-6K que meramente o contém.
* Referências que contêm o resíduo ambíguo `X` são sinalizadas com `has_x`, pois a checagem de
  integridade do POA1 as recusa como `-f`.

**Validação cruzada:** o epítopo `IGNRDFVE` (posições 6–13) reconstruído para DENV1 coincide
exatamente com o `bepipred_denv1.fasta` produzido manualmente pelo pesquisador.

Duas formas de uso:

```bash
# linha de comando: prepara e, com --run, executa POA1 + Conservancy local
python -m poa.cli.import_realdata \
    --pred-dir  ".../5_predicao_de_epitopo" \
    --world-dir ".../2_fastas_mundo" \
    --out       "./realdata_run" \
    --protein E --chikv-protein E1 --threshold 70 --run
```

Na GUI, a etapa 2 oferece o mesmo adaptador como caminho de upload, podendo adotar a referência
reconstruída como `-f` caso ainda não tenha sido fornecida.

**Resultado do teste ponta a ponta:** POA1 + Conservancy local rodaram para CHIKV, DENV1, DENV2 e
DENV4, com 13 a 16 epítopos por vírus e percentuais de conservação coerentes. DENV3 foi pulado
por conter `X` na referência.

---

## 7. Correções

### 7.1 Parser do NetCTL 1.2 (correção científica, aprovada pelo usuário)

Bug **pré-existente** em `netctlAntigenEpitopes`, com duas causas somadas:

1. As linhas eram inseridas usando o número da linha do arquivo como índice não sequencial
   (`.loc[i]`), operação rejeitada pelo pandas ≥ 2.
2. A lógica `if new_df_line[-1] != 'E': append('-')` estourava as 16 colunas nas linhas de epítopo
   reais, de modo que o critério de seleção `Identified_MHC_ligands == '<-E'` **nunca podia
   corresponder** — isto é, nenhum epítopo do NetCTL era selecionado.

Correção aplicada (índice sequencial; preenchimento apenas de linhas curtas até a contagem de
colunas, preservando o marcador `<-E` nas linhas de epítopo; guarda para entrada vazia retornando
as 7 colunas padrão). Como isso muda quais linhas são selecionadas, a alteração foi submetida e
**aprovada explicitamente pelo usuário**. Coberta por testes de caminho feliz, ausência de dados e
exclusão de linhas não-epítopo.

### 7.2 Análise de mais de uma espécie na mesma execução (correção)

Ao testar **duas espécies ao mesmo tempo**, as etapas 3 (POA1) e 4 (Conservancy) não conseguiam
identificar as espécies. Três causas independentes:

1. **Parser do BepiPred-3.0 lia o cabeçalho invertido.** `bp3_FastaAnalysis` fazia
   `sp, prot, idSeqNumber = header.split("_", 2)`, mas a convenção documentada — e usada pelo `-f`,
   pelo BepiPred-2.0, pelo PAP/IMED, pelo NetCTL, pelo MHC-II e pelo parser de "outros" — é
   `Proteína_Espécie_ID`. O resultado era espécie e proteína trocadas em toda a análise (a espécie
   virava o nome da proteína). Corrigido para `prot, sp, idSeqNumber`, com erro explícito quando o
   cabeçalho não tem os três campos.
2. **A GUI guardava um único arquivo por método.** Cada `st.file_uploader` (e cada importação de
   JSON do BepiPred-2.0) sobrescrevia `poa1_files[chave]`, então o segundo organismo **substituía
   silenciosamente** o primeiro; além disso o `-f` continuava sendo a referência da primeira
   espécie. Agora a GUI acumula uma *lista* de arquivos por método e os funde antes do POA1
   (`backend.merge_bepipred2` / `merge_fastas` / `merge_text`, orquestrados por
   `consolidate_sources`), e todas as referências reconstruídas são unidas em um único `-f`.
3. **A conservância comparava todas as espécies com o mesmo conjunto de proteínas.** O POA1 grava
   um `<ESPÉCIE>_epitopes.fasta` por espécie, mas `run_conservancy_for_dir` usava o mesmo conjunto
   de comparação para todos — epítopos de DENV1 medidos contra proteínas de CHIKV dão identidades
   sem sentido. Agora cada espécie é comparada com **o conjunto dela**: um conjunto próprio
   informado por espécie (ex.: o *world set* de diversidade), ou o subconjunto do `-f` cujo
   cabeçalho tem aquela espécie; o conjunto completo continua sendo o *fallback*, com aviso.

A etapa 3 passou a mostrar um diagnóstico de espécies (`backend.species_diagnostics`): quais
espécies estão no `-f`, quais aparecem nos epítopos, quais têm epítopos **sem proteína
correspondente no `-f`** (a falha que produzia conservância sem sentido) e quais ficaram sem
nenhum epítopo. A etapa 1 também passou a listar as espécies detectadas nos cabeçalhos.

### 7.3 Cabeçalhos fora da convenção — espécie não identificada

Duas formas dos arquivos reais em que **nenhuma espécie era extraída** (todos os epítopos com
`Specie = NaN`, gerando um único `nan_epitopes.fasta`):

* **Um JSON do BepiPred-2.0 com várias espécies.** Submeter `denv1_ns1` + `denv2_ns1` numa só
  execução devolve as chaves de antígeno `denv1` e `denv2` — um único campo, sem proteína e sem ID.
  A regex `(\w+?)_(\w+?)_(\w+)` do parser não casa e as duas espécies caem fora. O adaptador
  passou a aceitar um **mapeamento por antígeno** (`prepare_bepipred2(..., mapping=...)`), e a
  etapa 2 mostra um formulário com uma linha por antígeno (espécie + proteína), reconstruindo o
  `-f` de todos eles a partir dos arrays `AA`.
* **`-f` com cabeçalhos `Espécie_Proteína`.** `denv1_ns1` é lido como Proteína=DENV1,
  Espécie=NS1 — a ordem documentada é `Proteína_Espécie_ID`. A etapa 1 passou a exibir uma tabela
  de como cada cabeçalho foi interpretado e oferece a inversão dos dois primeiros campos
  (`normalize_fasta_headers`), que reescreve apenas os cabeçalhos, nunca as sequências.

O upload do BepiPred-2.0 passou a ser **um só controle, validado no envio**: as chaves de
antígeno são conferidas contra `Protein_Specie_ID` (`bepipred.BP2_ANTIGEN_KEY`, agora exposto como
constante). Se conformarem, o arquivo é usado direto; se não, ele é recusado com a lista das chaves
problemáticas e o formulário de mapeamento é aberto ali mesmo — antes era possível enviar o JSON
pelo uploader simples e só descobrir o problema duas etapas adiante, com todos os epítopos sem
espécie. A etapa 3 também avisa quando nenhuma espécie do `-f` coincide com as das predições,
apontando a caixa de inversão de cabeçalhos.

### 7.4 Saída da interface em `results/` (antes era `%TEMP%`)

A GUI criava um diretório de trabalho novo em `%TEMP%\poa_gui_XXXX` a cada inicialização — difícil
de achar, diferente a cada execução e sujeito à limpeza de disco do Windows. Passou a gravar tudo
em **`POA_Project/results/`** (`backend.default_results_root()`, sobrescrevível pela variável de
ambiente `POA_RESULTS_DIR`), mantendo a mesma árvore: `inputs/`, `mhcii/`, `poa1_out/`
(com `Conservancy Analysis/`), `conservancy_csv/`, `poa2_out/`, `cache/`. A pasta já estava no
`.gitignore`. O caminho aparece na barra lateral e na etapa 4, que também ganhou botões de download
dos CSVs de conservância. Como a pasta é fixa, uma nova execução sobrescreve arquivos de mesmo nome
— o aviso está na barra lateral.

### 7.5 POA2 no Windows — orientação correta e retomada da análise

A etapa 5 dizia "Instale com `pip install pyTMHMM`", o que **não funciona no Windows**: não há
*wheel* publicada (só macOS arm64/cp311), o build exige o Microsoft C++ Build Tools e o pyTMHMM
1.3.6 usa `np.int_t`, removido no numpy 2.x — verificado nesta máquina, o build falha em
`hmm.pyx:40: Invalid type`. A mensagem passou a ser específica por plataforma: no Windows aponta a
WSL (onde o ambiente conda `poa` do `TESTING_WSL.md` já tem o pyTMHMM compilado); nos demais
sistemas mostra o comando que realmente funciona (`numpy<2` + `--no-build-isolation`).

Como a sessão do Streamlit vive no processo e os arquivos não, trocar de ambiente significaria
refazer tudo. A etapa 5 ganhou **"Retomar com estes arquivos"**: detecta os CSVs em
`results/conservancy_csv/` (via `backend.resumable_run`), deixa escolher o `-f` entre os FASTAs de
`results/inputs/` e segue direto para o POA2. Junto veio a correção de um bug introduzido com a
tabela de cabeçalhos: `_refresh_reference()` zerava o `-f` sempre que a etapa 1 era renderizada
sem uploads, apagando uma referência adotada em outro lugar; agora só limpa quando explicitamente
pedido (`allow_clear=True`, usado pelo botão "Limpar espécies importadas").

### 7.6 `tmhmmAnalysis` — colunas de topologia desalinhadas (correção científica)

Bug **pré-existente** em `poa/core/topology.py`. As listas `Portion_Outside/TM/Inside` eram
preenchidas com um `append` **por registro do `-f` que casava** com o epítopo, e depois atribuídas
como colunas do DataFrame — o que só funciona se cada epítopo casar com exatamente um registro.
Duas consequências:

* **≠ 1 casamento → erro.** Com `proteins_all.fasta` (4 registros: `denv1_ns1`, `denv2_ns1`,
  `E_DENV1_ref`, `E_DENV2_ref`), cada epítopo casava com 2 → `ValueError: Length of values (32)
  does not match length of index (16)`.
* **Pior: corrupção silenciosa.** Se num mesmo lote alguns epítopos casassem 2 vezes e outros
  nenhuma, os totais podiam coincidir e as colunas eram gravadas **deslocadas** — a topologia de
  um epítopo atribuída a outro, sem erro nenhum.

A causa do casamento múltiplo era o teste `if EpitopeVirus not in seq_polyprot_id`, um *substring*:
`DENV1` está contido tanto em `DENV1_NS1` quanto em `E_DENV1_REF`.

Correção: (a) `sameVirus()` compara o **campo de espécie** do cabeçalho (`Protein_Specie_ID`),
mantendo o comportamento de substring só para cabeçalhos sem `_`; (b) as listas passaram a ter
**exatamente um valor por linha** — o primeiro registro casado, com aviso quando há mais de um
(`-f` duplicado), e o marcador `-` preservado na própria linha quando não há nenhum, em vez de
deslocar todas as seguintes. Como isso altera resultados nos casos de 0 ou ≥2 casamentos (antes:
erro ou desalinhamento), a mudança precisa da sua confirmação. Coberta por 4 testes novos.

A etapa 5 também passou a rotular cada `-f` candidato na retomada com as espécies que ele contém
(✅/⚠️/❌ conforme cubra as espécies dos CSVs), para não escolher por engano um FASTA que juntou
cabeçalhos não convertidos.

### 7.7 Retomada reconstrói também a tabela do POA1

A retomada da etapa 5 restaurava só o `-f` e os CSVs, então as etapas 6 (resultados) e 7
(visualização 2D/3D) continuavam bloqueadas por `poa1_result is None`. Agora
`backend.poa1_result_from_disk()` reconstrói a tabela de epítopos a partir dos
`<ESPÉCIE>_epitopes.fasta` que o próprio POA1 gravou — os cabeçalhos
`Espécie_Proteína_Método_Início_Fim` carregam todas as colunas padrão exceto `ID_Sequence` — e a
sessão retomada abre os mapas 2D/3D normalmente.

### 7.8 Robustez

* Criação automática dos diretórios de saída do POA1 (`-d`) e do POA2, evitando falha quando a
  pasta ainda não existe (sem impacto científico).
* Inclusão da raiz do repositório no `sys.path` para que `streamlit run poa/gui/app.py` resolva o
  pacote `poa` corretamente.
* Adequações ao pandas 3.0 encontradas durante a refatoração (`df.loc[índice_não_sequencial]` e
  `apply(axis=1)` sobre frames vazios agora levantam exceção; *Copy-on-Write* passou a ser padrão).

---

## 8. Testes

Suíte com **82 testes**, todos passando (`pytest`, configuração em `pytest.ini`):

| Arquivo | Cobertura |
|---|---|
| `tests/test_parsers.py` | os seis parsers de formato de predição |
| `tests/test_ranking.py` | filtros e ranqueamento (lógica científica) |
| `tests/test_pipeline.py` | orquestração POA1/POA2 |
| `tests/test_topology.py` | classificação de topologia de membrana |
| `tests/test_conservancy_local.py` | equivalência da conservação local com o formato IEDB |
| `tests/test_clients.py`, `test_services.py` | clientes de serviço, cache e fallback |
| `tests/test_browser_client.py` | automação de navegador (com mocks) |
| `tests/test_gui_backend.py` | orquestração da GUI e adaptador de dados reais |
| `tests/test_viz.py` | mapa 2D e geração do HTML do visualizador 3D |
| `tests/test_realdata_import.py` | adaptador dos arquivos originais da pesquisa |

```bash
.venv/Scripts/python.exe -m pytest      # Windows
python -m pytest                        # Linux/WSL
```

---

## 9. Documentação adicionada

* **`ARCHITECTURE.md`** (442 linhas) — descreve o pipeline original exatamente como estava, marca
  cada ponto de intervenção manual, apresenta a investigação de disponibilidade de cada ferramenta
  externa em 2026 e detalha o plano de automação. É a referência para entender *por que* cada
  estratégia de automação foi escolhida.
* **`TESTING_WSL.md`** (282 linhas) — guia passo a passo para rodar as partes que dependem de
  Linux (EMBOSS `antigenic`, NetCTL standalone, pyTMHMM) em WSL2 + conda, incluindo a receita que
  efetivamente funciona para instalar o pyTMHMM (`numpy<2` + `--no-build-isolation`) e o uso de
  `conda-forge --override-channels` para evitar o bloqueio dos termos de uso da Anaconda.
* **`README.md`** — nova seção 0 documentando a camada de automação, a GUI, a instalação dos
  extras e o status por ferramenta.
* **`requirements.txt`** e **`.gitignore`** — criados.

---

## 10. Compatibilidade — o que **não** mudou

Esta é a garantia central da entrega:

* **A lógica científica de ranqueamento e seleção de epítopos é a mesma**, com a única exceção da
  correção do NetCTL descrita em 7.1, que foi aprovada explicitamente por corrigir um
  comportamento comprovadamente errado.
* **A linha de comando é idêntica.** Os mesmos argumentos, os mesmos nomes de arquivo de saída:

  ```bash
  python POA1_v1.0.py -f proteins.fasta -d results -b3 bepipred3.fasta -e y
  python POA2_v1.0.py -g True -t 70 -d conservancy_csvs -f proteins.fasta -rf 1
  ```

* **Toda automação é opcional.** Qualquer ferramenta pode continuar sendo usada manualmente, via
  upload do arquivo de resultado, exatamente como no fluxo original.
* **A licença MIT foi mantida.**

---

## 11. Limitações e pendências conhecidas

| Item | Situação |
|---|---|
| Conservação local | Precisa ser validada contra uma exportação real do IEDB para o conjunto de dados em uso antes do uso em produção. |
| Automação de NetCTL / BepiPred-2.0 | Experimental. Os seletores foram escritos contra as páginas da DTU e podem quebrar; o upload manual é o caminho confiável. |
| POA2 no Windows | Bloqueado — o pyTMHMM não instala nativamente. Use o WSL conforme `TESTING_WSL.md`. |
| EMBOSS, BepiPred-3.0, Playwright | Não instalados na máquina de desenvolvimento; esses caminhos não puderam ser verificados em execução real no Windows. |
| DENV3 | Não processável enquanto o resíduo ambíguo `X` na referência não for resolvido a montante. |
| Visualização 3D | Depende de internet (3Dmol.js via CDN) e de numeração de resíduos do PDB compatível com as posições das predições. |
| Branch | `feature/viz-realdata` ainda **não foi enviada** ao repositório remoto. |

---

## 12. Histórico de commits (`main` → `feature/viz-realdata`)

| Commit | Descrição |
|---|---|
| `5b8bd7e` | docs: plano de arquitetura/automação, requirements e gitignore |
| `366c0b8` | refactor: extração do pipeline para o pacote `poa` (core/cli) com suíte pytest |
| `f908aa8` | feat(services): estrutura de cache/fallback manual + Conservancy Analysis local |
| `28798e5` | feat(services): clientes MHC-II (API IEDB), EMBOSS antigenic, pyTMHMM e BepiPred-3.0 |
| `5618e9f` | feat(gui): aplicação Streamlit orquestrando o pipeline (fluxo de 6 etapas) |
| `25d6931` | feat(services): automação Playwright *best-effort* para NetCTL/BepiPred-2.0 |
| `8b142d8` | docs(readme): camada de automação, GUI, instalação (EMBOSS/WSL) e CLI |
| `bf1e055` | fix(netctl): correção da seleção de epítopos (correção científica aprovada) |
| `424d620` | docs: guia passo a passo de testes em WSL (`TESTING_WSL.md`) |
| `24c049f` | docs: `conda-forge --override-channels` no guia WSL |
| `0a7434c` | docs: receita funcional de instalação do pyTMHMM em WSL |
| `3ab5555` | fix(pipeline): criação automática dos diretórios de saída de POA1/POA2 |
| `1d441a0` | fix(gui): raiz do repositório no `sys.path` para `streamlit run` |
| `7729348` | feat(realdata): adaptador para os arquivos originais do fluxo manual |
| `9e3490d` | feat(gui): mapa 2D de epítopos, visualizador 3D e etapa de importação de dados reais |

---

## 13. Como executar

```bash
# Ambiente
python -m venv .venv
.venv\Scripts\activate          # Windows
source .venv/bin/activate       # Linux/macOS
pip install -r requirements.txt

# Interface gráfica (fluxo completo de 7 etapas)
streamlit run poa/gui/app.py

# Linha de comando (inalterada)
python POA1_v1.0.py -f proteins.fasta -d results -b3 bepipred3.fasta -e y
python POA2_v1.0.py -g True -t 70 -d conservancy_csvs -f proteins.fasta -rf 1

# Dados originais da pesquisa
python -m poa.cli.import_realdata --pred-dir <5_...> --world-dir <2_...> \
    --out ./realdata_run --protein E --chikv-protein E1 --threshold 70 --run

# Testes
python -m pytest
```

Extras opcionais (não são pacotes pip ou são pesados — instale apenas se quiser aquela automação):

* **EMBOSS `antigenic`** (substitui o site PAP/IMED, fora do ar): `conda install -c bioconda emboss`
  (no Windows: WSL2 + conda, ou Docker).
* **BepiPred-3.0** local (PyTorch + ESM-2): `pip install bp3`.
* **Automação de navegador** para NetCTL 1.2 / BepiPred-2.0: `pip install playwright && playwright install chromium`.

---

*Documento gerado em 17 de agosto de 2026 · POA — Pipeline de Otimização de Antígenos · Licença MIT*
