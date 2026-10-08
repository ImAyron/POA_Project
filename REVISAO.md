# Revisão do POA — 08/10/2026

Leitura completa do repositório (objetivo, camadas, documentação e suíte) seguida de duas
correções. A validação **não** foi executada: esta máquina Windows não tem ambiente Python do
projeto montado, e o fluxo de trabalho adotado é editar no Windows e testar no Linux.

## Correções aplicadas

| Problema | Comportamento após a correção |
|---|---|
| O parser de outros preditores exigia 6 campos no cabeçalho, mas o README documenta o ID do NCBI como opcional (`se_houver`) | Cabeçalho de 5 campos (`Proteína_Espécie_Método_Início_Fim`) é aceito, com `ID_Sequence` vazio; a mensagem de erro passa a descrever as duas formas |
| A etapa 5 aplicava o limiar 70 quando os CSVs não o declaravam: ela caía em `ss.threshold`, que é o valor da etapa 4 apenas se a etapa 4 rodou na mesma sessão — numa análise retomada do disco era o default inicial | O limiar e o critério passam a ser pedidos na própria etapa 5 quando não há como lê-los dos CSVs, e o valor informado é o aplicado |
| Retomar uma análise não restaurava o limiar nem o critério, deixando a sessão em 70/`>=` enquanto os CSVs no disco tinham outro valor | `_resume_from_results` adota o limiar e o critério gravados junto dos CSVs, de modo que voltar à etapa 4 mostra o valor real em vez de 70 |
| Sair de uma etapa e voltar apagava o que havia sido preenchido nela: os filtros do POA2 voltavam a 60/100/60, os parâmetros do MHC-II e os controles 3D aos valores iniciais, e o modo da etapa 4 à reimplementação local | Política de estado explícita no topo de `poa/gui/app.py`: todo campo é espelhado numa chave de sessão comum, lida como padrão do widget e reescrita a cada execução |
| O PDB da etapa 7 era lido direto do uploader, sem ser gravado — sair da etapa descartava a estrutura e a visualização 3D sumia | O arquivo é gravado em `inputs/` como todos os outros uploads e guardado por caminho; a etapa 7 renderiza a partir dele |
| Não havia como começar do zero: só existia "Trocar de análise", que preserva os parâmetros de filtragem | Botão **"Nova análise"** na barra lateral — limpa toda a sessão, inclusive formulários e seletor de etapa, e desvincula a pasta para que uma nova seja criada |
| Dos 43 campos da interface, apenas 1 explicava o que esperava; o significado de `-t`, `-m`, `imin`/`imax` e das convenções de cabeçalho só existia no README | 34 textos de ajuda em `poa/gui/help_texts.py`, exibidos como "?" ao lado do campo (`help=` do Streamlit), cobrindo os campos não triviais das etapas 1, 2, 4, 5 e 7 |
| O projeto tinha 23 chamadas de log, **nenhuma** de erro: uma falha não deixava registro, as etapas não tinham fronteira e nenhuma linha dizia quantos epítopos saíram de cada método | `log_step` em `poa/logging_conf.py` delimita cada etapa com início, resultado contado e tempo, e registra `FAIL` com o tipo da exceção antes de repropagá-la; POA1 e POA2 passam a reportar contagens por método e por filtro |
| As mensagens não tinham forma comum: cada módulo escrevia à sua maneira e nada identificava se a linha era um passo, um resultado ou um erro | Gramática única `<TAG> <etapa> \| <mensagem>` com seis tags (`STEP`, `DONE`, `FAIL`, `NOTE`, `WARN`, `DATA`) e etapas nomeadas como caminho do backend (`POA1/bepipred-2.0`, `POA2/topology`, `service/iedb-mhcii`), aplicada também aos clientes de serviço |
| Os símbolos `▶ ✔ ✖` da primeira versão podiam levantar `UnicodeEncodeError` num console Windows em página de código legada — uma falha no log derrubando a execução | Tags ASCII, com teste que falha se deixarem de ser |
| A dependência `py3Dmol` era importada sem tratamento: faltando o pacote, a etapa 7 exibia "No module named 'py3Dmol'" | `RuntimeError` com a instrução de instalação, e seção nova no `README.md` com os três requisitos da visualização 3D (pacote, internet para o CDN do 3Dmol.js, numeração de resíduo do PDB) |
| As orientações de instalação listavam só os extras de instalação separada; o `py3Dmol` não era citado em lugar nenhum, e a ausência deixava ambíguo se a visualização 3D exigia algum passo a mais | A seção de instalação passa a dizer o que o `requirements.txt` cobre, nomeando `streamlit`, `plotly` e `py3Dmol`, e a verificação do `COMO_TESTAR.md` importa os dois últimos |
| `TESTING_WSL.md` e `COMO_TESTAR.md` repetem a lista de dependências à mão, e nada impedia que divergissem do `requirements.txt` | Teste estrutural compara as duas listas e falha quando um pacote obrigatório não aparece no guia |
| Não havia registro em arquivo: o log existia só no terminal de quem executou | Toda análise grava `poa.log` na própria pasta, em DEBUG, independente do nível do console — o registro viaja com os resultados entre Windows e Linux |
| `poa/core/parsers/conservancy.py` e `poa/services/conservancy_client.py` usavam o mesmo nome de logger (`conservancy`), tornando as mensagens indistinguíveis | O cliente passa a usar `conservancy_client`; um teste estrutural impede nomes repetidos |
| A etapa 5 mantinha a tabela do POA2 de uma execução anterior quando a execução falhava ou quando os parâmetros mudavam sem reexecutar | A execução descarta o resultado anterior antes de começar e retorna em caso de falha; se os parâmetros na tela não forem os da tabela, um aviso diz que ela é de execução anterior |
| `MUDANCAS.md` afirmava que a branch não havia sido enviada ao remoto | A linha descreve o estado real: publicada em `origin` e sincronizada |

A correção do parser **passa a aceitar** epítopos que eram recusados: um arquivo de outro preditor
sem ID do NCBI que antes abortava o POA1 agora é processado — é mudança de resultado no caso
afetado, no sentido de deixar de descartar dados válidos.

A correção do limiar **pode alterar a seleção** em análises retomadas cujos CSVs não declaravam o
limiar: o POA2 aplicava 70 e passa a aplicar o valor informado. Com `-idf` ligado o limiar é
critério de seleção, então nesses casos o conjunto de epítopos muda. Quando os CSVs declaram o
limiar — o caso normal, tanto da reimplementação local quanto de um download do IEDB — nada muda:
o dado continua tendo precedência sobre o formulário, que é a garantia introduzida em `2bd54d6`.

## Encontrado na autorrevisão das próprias alterações

| Problema introduzido | Correção |
|---|---|
| `use_log_file` (então `add_log_file`) **somava** handlers em vez de trocá-los: ao abrir uma segunda análise na mesma sessão, as linhas dela também eram gravadas no `poa.log` da primeira — desfazendo exatamente o que o log por análise existe para garantir | A função passa a remover a destinação anterior antes de instalar a nova; coberto por teste que escreve em duas análises e confere o isolamento |
| A validação de métodos do POA1 passou a usar teste de verdade (`if value`) onde o resto do pipeline usa `!= ""` — divergência silenciosa se algum argumento chegasse como `None` | Voltou a `!= ""`, o mesmo teste que `consolidate_predictions` aplica |

Nenhum dos dois era detectado pelos testes existentes: o primeiro passou por 20 casos de logging
porque nenhum abria dois arquivos diferentes.

## Removido por não se pagar

- Teste de comprimento mínimo dos textos de ajuda (60 caracteres): limiar arbitrário que não
  verifica se o texto é útil.
- Teste de formato dos nomes de etapa: policiamento de convenção com risco de falso positivo,
  já que `test_pipeline_stages_are_bracketed_by_log_step` cobre o que importa.
- Parametrização do teste de gramática sobre as seis tags: seis casos para uma propriedade só,
  colapsados em um.
- Asserção de que as tags têm quatro caracteres: restrição cosmética sem efeito real.

## Limpeza do terminal (após o primeiro teste manual no Linux)

| Problema | Correção |
|---|---|
| APIs depreciadas do Streamlit imprimiam um bloco "Please replace ..." a **cada reexecução** — dezenas de linhas por execução, mais saída que a própria análise | 17 `use_container_width=True` → `width="stretch"`; `components.html` → `_embed_html()`, que resolve `st.iframe` por atributo e cai no antigo quando não existe |
| `gio: http://localhost:8501: Operation not supported` a cada inicialização na WSL, da tentativa de abrir o navegador | `.streamlit/config.toml` versionado com `headless = true`, mais `logger.level = "warning"` e `gatherUsageStats = false` |
| `poa.topology` imprimia um banner e uma linha por sequência de referência em INFO, sem a gramática e competindo com o resultado da etapa | Banner removido (a fronteira `POA2/topology` já anuncia a etapa) e o detalhe por sequência passou a `DATA` em DEBUG |
| As linhas do próprio `logging_conf` não seguiam a gramática que ele define | `NOTE logging \| log file: ...` e `WARN logging \| ...` |
| `requirements.txt` pedia `streamlit>=1.30`, que não tem o parâmetro `width` que o código agora usa | Piso em `>=1.50`; os docs do Streamlit descrevem o parâmetro mas não dizem em que versão entrou, então é o lado conservador da onda de depreciação, não a versão exata |

Dois testes novos impedem o retorno: um falha se qualquer API depreciada reaparecer no `app.py`,
outro confere as três chaves do `config.toml`.

### Lacuna conhecida, deliberadamente não mexida

`poa/core/topology.py` usa `warnings.warn()` em quatro pontos para avisar sobre epítopos que não
casaram com nenhuma proteína ou que casaram com várias. Isso sai no formato do `warnings` do Python,
em stderr, e **não chega ao `poa.log`** — o registro da execução fica incompleto justamente nos
casos que mais importam. Converter para `logger.warning` quebraria `tests/test_topology.py`, que
depende de `pytest.warns(UserWarning)`, então é mudança de contrato e não cabia nesta limpeza.

## Verificado e correto (sem alteração)

Três pontos examinados em profundidade que **não** eram defeitos:

- A fronteira entre antígenos no BepiPred-2.0 e a ordem dos campos na topologia estão corretas; os
  nomes de epítopo são `Espécie_Proteína_Método_Início_Fim`, gerados por `fasta_out.py`, e a
  comparação de proteína em `topology.py` usa o campo certo.
- O guard de divisão por zero em `report.py` é necessário: sem ele um relatório vazio imprimiria
  `inf` e `-inf`.
- `stage_counts` não chamar `check_threshold` não produz funil enganoso: a etapa 5 já bloqueia a
  execução quando o operador dos CSVs discorda do objetivo, o limiar é lido dos próprios CSVs em
  vez do formulário, e limiares misturados fazem o funil sair vazio com a mensagem correta.

## Validação

```powershell
# parcial nesta máquina (só os módulos sem pandas/biopython/streamlit):
python -m pytest tests/test_logging.py tests/test_gui_help.py -q   # 53 passaram
# completo — rodar no ambiente Linux:
python -m pytest -q
git diff --check
```

- **A suíte completa não foi executada.** Não há `.venv` no Windows e a WSL não tem conda; o
  Anaconda local tem pytest mas não biopython/pandas/streamlit, então 15 dos 19 arquivos de teste
  não chegam a ser importados.
- **Os 45 casos que não dependem dessas bibliotecas foram executados e passaram** aqui
  (`test_logging.py` e `test_gui_help.py`). Um deles apontou um problema real na primeira
  tentativa — o `caplog` não via o aviso porque `setup_logging` desliga a propagação — corrigido
  no teste.
- Foram adicionados 7 testes: um em `tests/test_review_regressions.py` (cabeçalho de 5 campos),
  um em `tests/test_step_outputs.py` (aviso de tabela remanescente), dois em
  `tests/test_gui_app_threshold.py` (limiar informado na etapa 5; retomada adotando o limiar dos
  CSVs) e sete no novo `tests/test_gui_state.py` (parâmetros sobrevivendo à navegação; forma
  completa espelhada; objetivo seguindo os CSVs e não o formulário salvo; "Nova análise" limpando
  tudo e desvinculando a pasta; "Trocar de análise" preservando os parâmetros; nenhum upload lido
  sem ser gravado; toda chave de sessão declarada antes do uso). Os dois últimos leem a AST do
  `app.py`, então valem também para campos e uploads adicionados depois.
- Mais 4 casos (um deles parametrizado sobre 30 campos) no novo `tests/test_gui_help.py`: cada
  campo não trivial tem texto e o texto está ligado a ele, nenhum texto fica órfão, nenhum é curto
  demais para ser útil, e nenhuma ajuda é escrita direto no widget em vez da tabela compartilhada.
- Mais 20 casos no novo `tests/test_logging.py` (um parametrizado sobre as seis tags): fronteira
  de etapa com contagem e tempo, erro registrado e repropagado, `detail()` só em DEBUG, toda linha
  com tag e etapa, tags ASCII de comprimento fixo, `setup_logging` idempotente, nível vindo do
  ambiente (inclusive valor inválido), arquivo de log em DEBUG com console em INFO, arquivo
  adicionado uma vez por caminho, caminho inválido avisando sem interromper, nomes de logger
  únicos, toda etapa do pipeline delimitada por `log_step` e nomes de etapa dentro do padrão.
- A contagem de 155 testes citada no `README.md` e na revisão de 03/10 **precisa ser remedida**
  após estas adições. O repositório declara 150 funções `def test_`; o total coletado é maior por
  causa dos `@pytest.mark.parametrize`.

---

# Revisão do POA — 03/10/2026

Revisão das camadas de parsing, pipeline, serviços, backend, interface e documentação,
incluindo as alterações locais de prévias, funil e exemplos que já estavam em andamento.
A validação foi local, no Windows, usando o ambiente `.venv` do projeto.

## Correções aplicadas

| Problema | Comportamento após a correção |
|---|---|
| BepiPred-2.0 unia resíduos de antígenos diferentes quando as posições eram consecutivas | A troca de antígeno encerra o epítopo, preservando a espécie e a sequência de cada registro |
| O parser de outros preditores tratava cada linha do FASTA como um epítopo | Linhas são unidas por registro; IDs completos são preservados e cabeçalhos incompletos recebem erro explicativo |
| POA1 sem epítopos dividia por zero no relatório | Relatório válido com total e comprimento médio zero, sem valores infinitos |
| A topologia podia usar outra proteína da mesma espécie ou a primeira ocorrência de um peptídeo repetido | A associação considera também a proteína; usa as posições declaradas quando correspondem ao peptídeo |
| POA2 calculava topologia mesmo para uma seleção vazia | Retorna uma tabela vazia com as colunas de topologia, sem invocar o preditor |
| CSVs vazios, incompletos ou com limiares misturados interrompiam as prévias | A interface informa a falha; o funil não mostra uma contagem parcial enganosa; a validação do POA2 continua rejeitando limiares misturados |
| Busca de CSV aceitava `arquivo.csv.bak` | Aceita apenas a extensão CSV, sem distinguir maiúsculas, em ordem determinística |
| Funil e execução usavam valores padrão diferentes | Ambos usam 60% para identidade mínima e percentual mínimo de proteínas com match |
| Instalação básica exigia compilar pyTMHMM | pyTMHMM passa a ser instalado separadamente, seguindo o guia WSL |
| Guia de testes sugeria substituir o código local com `git checkout` | Os comandos executam diretamente as regressões sem sobrescrever arquivos |

As correções de parsing e associação de topologia podem alterar resultados nos casos afetados.
Os limiares de seleção existentes foram preservados. Quando as coordenadas não correspondem à
sequência, a topologia mantém a busca histórica pela primeira ocorrência, agora com aviso.
Referências duplicadas ainda usam o primeiro registro correspondente, com aviso.

## Validação

```powershell
.venv\Scripts\python.exe -m pytest -q
git diff --check
```

- 143 testes existentes passaram antes das correções.
- 155 testes passaram após as correções, incluindo 12 novos casos em
  `tests/test_review_regressions.py` e os testes de interface Streamlit.
- A sandbox negou acesso aos diretórios temporários do pytest; a suíte passou após execução
  com a permissão necessária para usar os diretórios temporários normais.
- Não foram executados os serviços remotos, o modelo BepiPred local nem o pyTMHMM real.
  Os testes dessas integrações usam simulações; o roteiro WSL cobre a validação manual.
- A instalação em ambiente limpo e a regeneração de `MUDANCAS.pdf` não foram executadas.

## Melhorias propostas para as próximas entregas

1. **Invalidar resultados dependentes ao mudar entradas ou parâmetros.** Hoje as pastas são
   reutilizadas e os geradores gravam apenas os arquivos atuais. Espécies removidas podem deixar
   FASTAs/CSVs antigos disponíveis. Usar um manifesto de entradas e saídas por execução, com
   geração em pasta temporária e publicação somente após sucesso.
2. **Conferir metadados contra todos os cabeçalhos CSV.** `inspect_directory` dá prioridade ao
   JSON de metadados. Validar sua consistência com cada arquivo evita que uma substituição manual
   preserve um limiar antigo e rotule dados incorretamente.
3. **Explicitar a política para referências ambíguas.** `merge_fastas` mantém o primeiro cabeçalho
   duplicado. Rejeitar sequências conflitantes e carregar o ID da referência até o POA2 permitiria
   distinguir cópias da mesma proteína e reduzir a necessidade de buscas por ocorrência.
4. **Centralizar validações e leitura de arquivos.** Compartilhar a cadeia de filtros entre
   execução e funil, validar limites e cabeçalhos antes de gravar artefatos e restringir nomes de
   arquivos de uploads/adaptadores às pastas de destino.
5. **Automatizar a matriz de testes e instalações.** Configurar CI em Windows e Linux, testar
   instalação limpa e separar dependências de execução, desenvolvimento e automação opcional.
   Validar pyTMHMM real no Linux com um conjunto pequeno e resultados de referência.

Essas propostas não estão implementadas neste commit. Esta revisão e a suíte de testes não
substituem a validação dos resultados científicos contra conjuntos de referência externos.
