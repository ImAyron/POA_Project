# A inversão `Espécie_Proteína` alcançando as predições — 26/10/2026

Relatado a partir de um caso real: um consenso de E2 de alphavirus, cabeçalho
`alphavirus_e2_consenso` — **espécie alphavirus, proteína E2**. Com a caixa da etapa 1 marcada, a
interface ainda mostrava `Espécies: E2` nas etapas 2 e 3.

Suíte: de 526 para **533 passando**, no ambiente conda `poa` da WSL.

## A raiz

A inversão era aplicada **só ao FASTA de referência** (`normalize_fasta_headers` sobre
`ss.ref_raw`). Nada a aplicava às predições, e a mesma string era lida de duas formas:

| Origem | Cabeçalho | Resultado |
|---|---|---|
| Referência (`-f`) | `alphavirus_e2_consenso` | proteína `E2`, espécie `ALPHAVIRUS` ✓ |
| Predição (chave do JSON `-b2`) | `alphavirus_e2_consenso` | proteína `ALPHAVIRUS`, espécie `E2` ✗ |

Daí os três diagnósticos simultâneos da etapa 3, todos corretos: epítopos com espécie `E2`, `-f`
com espécie `ALPHAVIRUS`, nenhuma coincidência.

**O adaptador da etapa 2 não ajudava**, por um motivo sutil: ele só é oferecido quando as chaves de
antígeno falham no regex da convenção, `(\w+?)_(\w+?)_(\w+)`. A chave `alphavirus_e2_consenso` tem
três campos e **passa** — a verificação é sintática, o erro é semântico, e `alphavirus_e2_consenso`
e `e2_alphavirus_consenso` são indistinguíveis para um regex.

## A correção

Num ponto só: logo após `consolidate_predictions` e antes do relatório, dos arquivos por espécie e
de tudo que vem depois. Cada método carrega a convenção à sua maneira — chave de JSON no `-b2`,
cabeçalho FASTA no `-b3` e no `-x`, linha `>ID` no `-p`, identificador dentro da página no `-n`,
**nome de arquivo** no `-m` — e tratar os seis formatos separadamente seriam seis lugares para
manter em sincronia, cada parser tendo que conhecer uma configuração da interface.

Exposto também na linha de comando (`-swap`), para a interface e a CLI não divergirem, e desligado
por padrão.

Verificado com os arquivos reais do caso:

```
-f reescrito : >E2_ALPHAVIRUS_CONSENSO
swap=False -> Especie ['e2']         faltando no -f ['E2']  FASTA e2_epitopes.fasta
swap=True  -> Especie ['alphavirus'] faltando no -f []      FASTA alphavirus_epitopes.fasta
```

## Duas falhas de comunicação, corrigidas junto

- **O rótulo da caixa mentia por omissão:** "Meus cabeçalhos estão em `Espécie_Proteína`", no
  plural e sem ressalva, agindo em um lugar só. Passa a dizer "na referência e nas predições".
- **O diagnóstico da etapa 3 prescrevia o remédio errado:** mandava marcar a caixa quando ela já
  estava marcada, porque o desencontro que a inversão pela metade produzia era idêntico ao de não
  tê-la marcado. Agora distingue os dois casos.

## Consequência para análises já feitas

Qualquer análise rodada com a caixa marcada e cabeçalhos `Espécie_Proteína` agrupou a conservância
pelo nome da proteína. **Precisa ser reexecutada** a partir da etapa 3.

---

# Etapa 7: buscar a estrutura e usar o SIFTS — 19/10/2026

Duas coisas que o repositório tinha pela metade, apontadas na revisão do PR #2 e ainda abertas.
Suíte: de 513 para **526 passando**, no ambiente conda `poa` da WSL.

## O SIFTS era código morto

`map_epitopes` aceitava `sifts_offset` e **nenhum chamador passava um**. O POA0 calculava o
mapeamento e gravava em `estrutura.offset_regiao`; a etapa 1 dizia ao usuário que "a estrutura
baixada aqui já vem com o mapeamento SIFTS — que é o que a etapa 7 usa"; e a etapa 7 caía no
alinhamento. A afirmação era falsa, e a diferença não é cosmética: SIFTS declara a correspondência
resíduo a resíduo, o alinhamento a estima.

`backend.structure_mapping` lê o registro de volta e a etapa 7 o repassa. A interface agora diz
qual das duas rotas resolveu a numeração, em vez de deixar implícito.

Uma guarda que vale registrar: o relatório descreve **uma** estrutura. Se o arquivo no visualizador
for outro, o mapeamento é ignorado — aplicar o offset de uma estrutura a outra colocaria os
epítopos com confiança sobre os resíduos errados, que é a falha do 4C23 agravada pelo selo de
"exato".

## Buscar a estrutura na etapa 7

Baixar do RCSB só existia na etapa 1, enquanto o POA0 gerava a referência. Quem enviou o próprio
FASTA, ou quis trocar de estrutura depois, tinha que ir ao site. Agora há um expansor na etapa 7
que lê o accession do terceiro campo do cabeçalho, pergunta ao PDBe quais entradas cobrem a
proteína e baixa a escolhida — reaproveitando `list_structures` e `fetch_structure`, sem código de
rede novo. O arquivo entra pelo mesmo `save_upload` com `role="structure"`, então já é gravado na
análise e recarregado depois.

Sem accession no cabeçalho, a etapa diz isso e deixa o upload como caminho — em vez de buscar algo
que não consegue identificar. E uma estrutura baixada aqui ainda é conferida por **alinhamento**:
o offset exato depende de qual trecho da sequência a referência cobre, e só o POA0 registra isso.
A interface diz essa diferença em vez de insinuá-la.

## Um defeito que o teste pegou

O regex de accession do UniProt começava com `[A-NR-Z]`, que **exclui O, P e Q** — e `P29990`,
a NS1 de dengue, é exatamente a forma que ele recusava. Trocado pelos dois padrões oficiais do
registro, com os dois formatos (6 e 10 caracteres) cobertos por teste parametrizado.

## Um achado da revisão do PR #2 que se verificou improcedente

"O destaque é pintado em todas as cadeias" **não é mais um defeito vivo**: a etapa 7 atribui
`chain = report.chain` após o mapeamento, e nos caminhos em que não há verificação os intervalos
ficam vazios — nada é destacado. Deixo registrado para não ser reinvestigado.

## Ainda em aberto

- **`.docx` versionado** (`Relatorio_Alternativas_Ferramentas_POA.docx`, 40 KB, fora do
  `.gitignore`). Não removido: apagar arquivo versionado de outro autor é decisão do time.
- Cadeia única escolhida por ordem de arquivo; *outage* de API indistinguível de "sem dados";
  apagar análise falhando no Windows; saída do POA2 sobrescrevendo quando só `-m`/`-imin`/`-imax`/
  `-idf` mudam; e o denominador da conservância versus o IEDB, que continua sendo decisão
  científica.

---

# Recarregar uma análise — 16/10/2026

Reabrir uma análise passou a devolver **todos** os arquivos dela à etapa em que foram usados.
Suíte: de 318 para **343 passando**, executada no ambiente conda `poa` da WSL.

## O problema

`_resume_from_results` cobria só a etapa 5: CSVs de conservância, o `-f` e a tabela do POA1, atrás
de um botão. Tudo o mais ficava no disco com a interface mostrando campo vazio — os arquivos de
predição da etapa 2, os conjuntos de comparação por espécie da etapa 4, a estrutura da etapa 7.
Continuar uma análise significava reenviar o que ela já continha.

## A decisão de desenho

O nome do arquivo **não diz de que etapa ele veio**. Um FASTA de referência e um de "outros
preditores" são ambos FASTAs em `inputs/`; um arquivo que a pessoa nomeou não segue convenção
nenhuma; e o POA0 escreve `<PROT>_<ESPECIE>.fasta`, que é indistinguível de uma referência enviada
à mão. Inferir pelo nome seria adivinhação.

Por isso a procedência é **gravada na chegada**, em `input_history.json` na pasta da análise:
papel, método, espécie, nome original e quando entrou. Os oito pontos de upload da interface
passam a registrar, e o POA0 registra o que ele mesmo escreve.

Dois pontos que valem registro:

- **Caminhos relativos à pasta da análise.** A mesma pasta é aberta do Windows e da WSL, onde o
  caminho absoluto é outro; gravar o absoluto de um lado tornaria o registro inútil do outro.
- **O histórico nunca é condição para salvar.** O arquivo é gravado primeiro; falhar ao registrar
  a linha de procedência vira `WARN` e segue. Perder uma linha de histórico vale menos que perder
  o upload.

## Análises que já estavam no disco

Essas não têm registro, e são justamente as que o pedido quer recarregar. Para elas há
classificação por nome (`_classify_legacy_input`), que é **melhor esforço** — a interface avisa
quando foi esse o caminho e pede para conferir. O histórico sempre vence o nome quando ambos
existem; há teste para esse caso, com um arquivo cujo nome diria uma coisa e o registro diz outra.

## Um defeito encontrado pelo próprio teste

O MHC-II é o único método que o POA1 consome como **pasta** de `Proteína_Espécie.html`, não como
arquivo. Cada HTML é registrado individualmente — isso é a procedência correta — mas a primeira
versão devolvia os caminhos individuais, e o POA1 teria recebido um arquivo onde espera um
diretório. Corrigido: a etapa recebe a pasta de volta.

## Testes

- `tests/test_restore_analysis.py` (23 casos): o registro em si, caminhos relativos, sem
  duplicatas, histórico ilegível ignorado, salvar sem papel não grava nada; a restauração por
  etapa, arquivo apagado do disco sai da restauração, análise vazia; e 12 casos parametrizados da
  classificação por nome, mais a precedência do histórico sobre o nome.
- `tests/test_gui_app_runs.py` (2 casos novos, via `AppTest`): reabrir pela interface e encontrar
  cada arquivo na sua etapa; e a restauração acontecendo **uma vez**, não a cada interação — repetir
  desfaria o que a pessoa mexeu desde então.

---

# Revisão do PR #2 — 16/10/2026

Revisão do POA0, da verificação estrutural e das correções das etapas 4 e 5. Diferente das
revisões anteriores, esta foi **medida**: a suíte e os casos abaixo rodaram no ambiente conda
`poa` da WSL (`conda run -n poa python -m pytest`), não apenas lidos.

Suíte antes: 312 passando. Depois das correções: **318 passando**.

## O que o PR acerta

- A descoberta do **PDB 4C23** é real: duas análises de DENV usaram como estrutura a L-fuculose
  quinase de *S. pneumoniae*, com 0 dos 19 epítopos presentes. Há ferramenta para reproduzir.
- O desenho do guard está certo onde importa — identidade por alinhamento e não janela
  deslizante, numeração real do arquivo em vez de 1-based presumido, recorte ao observado em vez
  de extrapolação, falha fechada em cadeia vazia, e na recusa a estrutura é mostrada **sem**
  destaque em vez de escondida.
- `validate_conservancy.py` endereça a pendência mais antiga do projeto.
- A exclusão de análise está corretamente contida: caminhos resolvidos, `..`/symlink/absoluto
  barrados, `cache` e raiz recusados, confirmação por nome.
- Os `warnings.warn` do `topology.py` passaram a chegar ao `poa.log` via `captureWarnings` —
  fechando uma lacuna que a revisão de 08/10 tinha deixado registrada.

## Correções aplicadas

| Problema (medido) | Correção |
|---|---|
| **O guard aprovava a estrutura errada.** A identidade é uma razão sobre as colunas alinhadas, e a aprovação só exigia `columns > 0`. Com os dois gaps terminais livres, o alinhador descarta quase toda a estrutura como sobra e guarda um trecho curto bem pareado. Medição: de 216 pares de sequências sem relação, **1 foi aceito a 100% de identidade sobre uma única coluna** — 0,7% de uma estrutura de 150 resíduos contra referência de 500 aa | Piso de cobertura `DEFAULT_MIN_COVERAGE = 0.50`: a identidade só vale se a estrutura tiver alinhado em pelo menos metade do seu comprimento. O caso legítimo — cristal de um domínio de poliproteína — continua passando, e há teste para os dois lados |
| `min_identity` era ignorado na rota SIFTS: a comparação usava a constante do módulo, então `check_structure.py -i 90` não tinha efeito na única rota que o POA0 usa para estruturas que ele baixa | O argumento é propagado e honrado nas duas rotas |
| **`-r` fora dos limites gerava FASTA vazio e saía com código 0.** `sequence[4999:6000]` de uma proteína de 1420 aa é `""`; o POA0 escrevia cabeçalho sem sequência, registrava `tamanho: 0`, e o POA1 rodava sobre referência vazia | Validação `1 <= start <= end <= len(sequence)` com mensagem que distingue erro de digitação de região descrita sobre outra versão da sequência |
| **`max_identity` passou a devolver `None` e um chamador não foi atualizado:** `best_world_match` fazia `max(ident for ...)` e levantava `TypeError` quando o primeiro registro de um world FASTA estava vazio — caso que a etapa 1 já avisa que acontece | O arquivo não pontuável é registrado em `WARN` e pulado |
| Falha no download da estrutura destruía o relatório: `_add_structure` rodava antes do `json.dump`, então o FASTA ficava no disco sem registro da região e do offset de onde veio | O relatório é gravado sempre, com `estrutura_erro`, e a exceção é repropagada depois para o código de saída continuar dizendo que algo falhou |
| O sidecar contradizia o CSV ao lado dele: gravava `n_proteins: 2` enquanto o CSV dizia `100.00% (1/1)`, e é do sidecar que a grade de 100/n do aviso do `-m` e da tabela da etapa 5 é derivada | Passa a gravar também `n_comparable` — o denominador que o CSV de fato usou |

## Pendências deixadas para decisão, não corrigidas

**O denominador da conservância diverge do IEDB.** Excluir do denominador as proteínas mais curtas
que o epítopo é defensável, e a justificativa no docstring é boa — mas o IEDB usa
`matched / |conjunto|`. Isso torna o resultado local **não equivalente** ao do site, e a própria
`tools/validate_conservancy.py` vai reportar divergência em qualquer conjunto com fragmentos
curtos, imprimindo "não use os números locais". Não alterei: é decisão científica, não defeito.
Ou se alinha ao IEDB, ou se documenta a divergência e se ajusta o validador.

Residual relacionado: quando `comparable == 0` a linha sai `0.00% (0/0)` com
`Minimum identity 0.00%` — a mesma falha que o docstring diz corrigir.

## Fragilidades apontadas e não corrigidas neste commit

Confirmadas por leitura, fora do escopo destas correções:

- **O caminho SIFTS é código morto.** `dbsearch` calcula e grava `offset_regiao`, mas nem a
  interface nem `check_structure.py` passam `sifts_offset=`. Enquanto isso a etapa 1 afirma ao
  usuário que "a estrutura baixada aqui já vem com o mapeamento SIFTS… que é o que a etapa 7 usa".
  A afirmação é falsa como está ligado.
- **O destaque é pintado em todas as cadeias.** Verifica-se uma cadeia, mas o `chain` é omitido da
  seleção do 3Dmol quando o campo está vazio — os números de resíduo vão para todo parceiro de um
  complexo. É um resíduo do próprio problema que o PR corrige.
- **Só uma cadeia é examinada, escolhida por ordem no arquivo.** Num complexo antígeno-anticorpo
  com H/L primeiro, o antígeno nunca é tentado e uma estrutura boa é recusada.
- **Outage de API é indistinguível de "sem dados"**: 404, 429, 500, 503 e timeout viram a mesma
  exceção, e três chamadores a convertem em `[]` logando "nenhuma estrutura conhecida".
- **Apagar a análise recém-fechada falha no Windows** — o `FileHandler` do `poa.log` fica aberto e
  o `rmtree` bate em `WinError 32`. No Linux o unlink funciona, então o teste não pega.
- **Quatro avisos novos quebram a gramática de log** (`logger.warning` sem o prefixo
  `WARN <etapa> |`), então `grep 'WARN'` e `grep 'POA2/'` não os encontram.
- **Nome de saída do POA2 ainda sobrescreve** quando só `-m`, `-imin`, `-imax` ou `-idf` mudam.
- `pdb_id` é interpolado em caminho e URL sem validar formato; `--limite` é ignorado por três dos
  quatro modos de listagem; `database_client` não usa a camada de cache e a interface dispara duas
  requisições idênticas a `best_structures` por clique; nenhum teste exercita `_get_json`, então
  nada verifica que um 4xx/5xx vira `ServiceUnavailable`.
- **Higiene:** um `.docx` de 40 KB foi commitado sem estar no `.gitignore`; README e MUDANCAS ainda
  dizem "155 testes".

---

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
