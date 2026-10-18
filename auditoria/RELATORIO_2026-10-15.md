# Auditoria do POA — 15/10/2026

Revisão da base `de92f62`, no Windows, cobrindo os pontos de entrada, parsers,
POA0/POA1/POA2, serviços, cache, persistência, interface Streamlit, visualização,
ferramentas auxiliares e documentação de execução.

**Resultado: 312 testes existentes passaram; 15 problemas foram reproduzidos em
verificações adicionais.** A suíte atual não cobre esses cenários. Alguns riscos
de persistência e metadados já constavam como propostas em `REVISAO.md`; foram
reproduzidos novamente no código atual, não tratados como descobertas inéditas.

Esta entrega acrescenta o relatório e o script de reprodução. Não aplica correções
à aplicação nem modifica os dados das análises em `results/`.

## Verificação executada e limites

```powershell
.\.venv\Scripts\python.exe -m pytest -q -p no:cacheprovider
# 312 passed in 79.63s

.\.venv\Scripts\python.exe -m pip check
# No broken requirements found.

.\.venv\Scripts\python.exe auditoria/reproduzir_2026_10_15.py --gui
# 15 casos confirmados; inclui a interface real via Streamlit AppTest.

git diff --check
```

O Python do PATH não tinha pytest; foi usado o ambiente `.venv` já existente.
As primeiras tentativas na sandbox falharam por acesso aos temporários/cache;
a execução completa com acesso aos temporários passou. Esses erros de ambiente
não foram contabilizados como bugs do aplicativo.

Também foram verificados os três comandos `POA0_v1.0.py --help`,
`POA1_v1.0.py --help` e `POA2_v1.0.py --help`, e executado o POA1 pela CLI sobre
`exemplos/limiar/`, com relatório e exportação XLSX. A conservância local e os
filtros selecionaram os dois epítopos a 70% e apenas DENV2 a 90%, conforme o exemplo.

Ambiente: pandas 3.0.6, numpy 2.5.3, Biopython 1.88, Streamlit 1.64.0,
pytest 9.1.1, Plotly 7.1.0 e py3Dmol 2.5.5.

**Limites:** pyTMHMM, bp3 e Playwright não estão instalados neste ambiente.
Não foram executados os preditores externos reais, o EMBOSS, chamadas a APIs,
o fluxo Linux/WSL nem uma instalação limpa. Os testes dessas integrações usam
simulações. O teste 3D verifica o seletor entregue ao visualizador; não renderiza
WebGL em um navegador. Não houve nova validação científica contra exportação real
do IEDB. Os arquivos de resultados anteriores não foram usados como evidência de
correção científica.

## Problemas prioritários — P1

P1 significa risco de perda de dados ou interpretação incorreta de resultados em
fluxos suportados. Os números abaixo correspondem aos casos do script.

### 01 — Espécies removidas voltam pela pasta de saída do POA1

- **Local:** `poa/core/fasta_out.py:21–42`; consumo em
  `poa/services/conservancy_client.py:260` em diante.
- **Reprodução:** executar POA1 com A+B e repetir na mesma pasta somente com A.
  A tabela nova contém A, mas `B_epitopes.fasta` permanece. A conservância processa
  ambos e registra B nos metadados novos; `stale_csvs()` retorna lista vazia.
- **Impacto:** resultados antigos voltam a integrar a análise, inclusive sem o
  aviso de CSV residual. O mesmo problema ocorre quando a nova seleção fica vazia.
- **Correção:** publicar um conjunto completo de saídas por execução e consumir
  exclusivamente o inventário correspondente. Usar pasta de preparação e substituir
  as saídas geradas somente após sucesso, preservando arquivos enviados pelo usuário.

### 02 — Metadados antigos prevalecem sobre CSV incompatível até em modo estrito

- **Local:** `poa/core/parsers/conservancy.py:200–204`.
- **Reprodução:** CSV com cabeçalho `>=90%`, sidecar com `>=70%`, ambos apontando
  para o mesmo arquivo. `EpitConservAnalysis(..., strict=True)` aceita a entrada
  e devolve a coluna rotulada como `>=70`.
- **Impacto:** a validação afirma um critério diferente daquele usado no cálculo;
  substituir manualmente um CSV pode gerar essa divergência.
- **Correção:** comparar sidecar, cabeçalhos de todos os CSVs e parâmetros solicitados;
  rejeitar contradições. Vincular os metadados aos hashes dos arquivos.

### 03 — Automação MHC-II sobrescreve referências da mesma proteína/espécie

- **Local:** `poa/gui/app.py:627–635`; `poa/services/mhcii_client.py:113–121`.
- **Reprodução:** duas referências `E_A_ref1` e `E_A_ref2` são processadas
  separadamente, mas ambas gravam `E_A.html`; a segunda substitui a primeira.
- **Impacto:** a interface pode informar duas predições concluídas, enquanto só
  a última chega ao POA1. O parser também descarta a identificação da referência.
- **Correção:** incluir o ID da referência no artefato e preservá-lo na tabela;
  registrar todas as respostas em um manifesto, sem colisões de nome.

### 04 — Merge de FASTAs descarta sequências conflitantes silenciosamente

- **Local:** `poa/gui/backend.py:448–466`.
- **Reprodução:** dois FASTAs com `E_A_ref`, contendo `AAAA` e `CCCC`, produzem
  somente `AAAA`, sem erro ou aviso.
- **Impacto:** referências e resultados BepiPred-3 podem desaparecer. A validação
  de cabeçalhos duplicados na topologia não percebe o conflito, pois o merge já
  eliminou uma das sequências.
- **Correção:** deduplicar apenas pares ID+sequência idênticos; rejeitar IDs iguais
  com conteúdo diferente e informar os arquivos envolvidos.

### 05 — Visualização usa a primeira sequência encontrada, não a referência ativa

- **Local:** `poa/gui/backend.py:888–920`; chamada em `poa/gui/app.py:1479–1489`.
- **Reprodução:** `comparison_A.fasta` e `ref_0_original.fasta` têm a mesma
  proteína/espécie e sequências diferentes. O primeiro, por ordem alfabética,
  é escolhido para conferir/mapear a estrutura.
- **Impacto:** rejeição de uma estrutura correta ou mapeamento das coordenadas
  dos epítopos contra outra variante. A busca também ignora referências com
  extensões `.fa`, `.faa` ou `.txt`, embora a etapa 1 as aceite.
- **Correção:** usar a referência efetivamente associada à predição, com ID e hash
  persistidos; evitar inferir sua origem procurando qualquer FASTA em `inputs/`.

### 10 — Valida uma cadeia, mas destaca resíduos em todas as cadeias

- **Local:** `poa/gui/app.py:1488–1508`; `poa/gui/viz.py:129–133`.
- **Reprodução:** estrutura com A=`ACDE` e B=`CCCC`, deixando o campo cadeia vazio.
  O mapeamento valida A, mas a chamada de renderização recebe `chain=None`;
  o seletor de destaque contém apenas `resi`.
- **Impacto:** a numeração validada em A também colore B, apesar de B não ter
  sido validada e poder representar outra proteína.
- **Correção:** transmitir `report.chain` ao visualizador. Para destacar várias
  cadeias, mapear e validar cada uma separadamente.

### 13 — Validador científico aprova métricas NaN

- **Local:** `tools/validate_conservancy.py:62–65,116`.
- **Reprodução:** CSV de comparação com `NaN` nas métricas percentual, identidade
  mínima e máxima. O programa retorna código 0 e declara equivalência.
- **Causa:** `abs(NaN - valor) > tolerancia` é falso; não há verificação de finitude.
- **Correção:** exigir métricas presentes, finitas e dentro do domínio válido;
  validar também a tolerância. Distinguir dados inválidos de resultados equivalentes
  e rejeitar validação sem epítopos comparáveis.

### 14 — Reexecutar POA1 preserva o resultado anterior do POA2

- **Local:** `poa/gui/app.py:893–894`; apresentação em `step_results()`.
- **Reprodução:** AppTest com resultado POA2 anterior; clicar em executar POA1.
  `conservancy_ready` fica falso, mas `poa2_result` ainda aponta para a planilha
  anterior. A etapa 6 continua usando esse objeto.
- **Impacto:** resultados finais antigos aparecem ao lado de entradas/saídas novas.
  A checagem de parâmetros do POA2 não identifica troca de conteúdo de arquivos.
- **Correção:** invalidar resultados dependentes quando mudarem entradas ou uma
  etapa anterior for executada; associar cada resultado à versão/hash das entradas.

## Outros problemas confirmados — P2

| Caso | Local | Reprodução e impacto | Correção proposta |
|---|---|---|---|
| **06** | `poa/core/dbsearch.py:191` | POA0 aceita `0-3`, `4-2` e `2-999` para sequência de 4 resíduos. Grava comprimentos 0, 0 e 3 sem rejeitar coordenadas inválidas; o relatório mantém o fim 999. | Exigir `1 <= início <= fim <= comprimento` antes de criar saídas. |
| **07** | `poa/core/parsers/mhcii.py:51–54` | `os.walk` encontra `sub/E_A.html`, mas o caminho devolvido é `raiz/E_A.html`, que não existe. | Usar a pasta devolvida pelo `os.walk`; ordenar arquivos e validar o nome. |
| **08** | `poa/gui/backend.py:403–407` | `save_upload(dest, '../fora.txt', ...)` escreve fora de `dest`. A prova usa apenas dados sintéticos dentro do temporário. | Rejeitar caminhos absolutos, separadores e componentes de travessia; resolver e conferir contenção antes de gravar. Aplicar a mesma regra a nomes derivados de espécies/adaptadores. |
| **09** | `poa/gui/backend.py:215–230` | Criar `Teste A` e depois `Teste-A` no mesmo dia devolve a mesma pasta, com as saídas anteriores. A tela promete uma nova análise independente. | Criar ID único ou recusar colisão explicitamente, oferecendo reabertura separada. |
| **11** | `poa/core/report.py:13–27` | Um epítopo válido de 4 resíduos de A é rejeitado porque a referência de B tem 2 resíduos, mesmo sendo outra espécie. | Conferir viabilidade contra a referência correspondente. Tratar proteínas não comparáveis na conservância por epítopo, como o serviço local já faz. |
| **12** | `poa/core/parsers/conservancy.py:368–373`; `poa/gui/app.py:967,984` | Metadados de upload manual são gravados sem inventário `species`. Todos os CSVs passam a ser classificados como sobras; `strict=True` rejeita um conjunto válido. | Diferenciar inventário ausente de inventário vazio; gravar os arquivos do upload e validar sua proveniência. |
| **15** | `poa/gui/app.py:1472` | POA1 vazio com PDB já carregado deixa o seletor de proteína vazio. A etapa 7 chama `None.split()` e encerra com `AttributeError`. Confirmado via AppTest. | Tratar seleção vazia antes dos controles de destaque e exibir a estrutura sem epítopos. |

No caso 08 foi demonstrada a falha de contenção da função de escrita. Não foi
demonstrado ataque remoto pelo componente de upload do Streamlit; a exposição
depende dos chamadores e de como os nomes chegam à função.

## Melhorias propostas e ordem de implementação

1. **Integridade das execuções:** resolver 01, 02 e 14 juntos. Um manifesto deve
   registrar entradas, parâmetros, hashes, saídas e dependências entre etapas.
   Publicar resultados somente após sucesso e conservar versões anteriores para
   consulta. Isso também evita misturar arquivos após uma falha parcial.
2. **Identificação das referências:** resolver 03, 04 e 05 usando IDs explícitos
   desde a importação até a topologia e a visualização. Hoje os nomes de epítopos
   exportados não preservam `ID_Sequence`, o que dificulta distinguir variantes.
3. **Validação centralizada:** corrigir 06, 08, 11, 12 e 13; validar limites,
   coordenadas, cabeçalhos e dados numéricos antes de gravar saídas. Tornar o
   contrato de `<`/`<=` e `>`/`>=` explícito, sem tratar igualdade como equivalente
   apenas para compatibilizar rótulos. Preservar limiares decimais em vez de
   convertê-los para `int` na interface/backend.
4. **Visualização vinculada aos dados:** corrigir 10 e 15, utilizar a referência
   ativa e não destacar posições cruas quando faltar referência verificável.
   Representar resíduos por cadeia, número e código de inserção, em vez de somente
   intervalos numéricos. Integrar à GUI o mapeamento SIFTS persistido pelo POA0.
5. **Testes de transição:** converter os cenários desta auditoria em regressões
   que exijam o comportamento correto, especialmente reexecução, falha parcial,
   troca de arquivos, referências conflitantes e resultado vazio. Os testes atuais
   cobrem vários fluxos isolados, mas não todas as transições entre eles.
6. **Integração contínua e instalação reproduzível:** executar a suíte em Windows
   e Linux; verificar instalação limpa; separar dependências básicas, GUI,
   desenvolvimento e preditores opcionais. A lista atual inclui Playwright mesmo
   descrito como opcional. Atualizar a contagem de testes da documentação: o README
   ainda cita 155, enquanto a execução atual coletou e aprovou 312.
7. **Serviços e concorrência:** versionar cache pelo modelo/binário efetivo e gravar
   entradas de forma atômica; esperar um marcador de resultado final no navegador,
   pois `networkidle` sozinho não comprova conclusão do trabalho remoto. Isolar
   logs por execução/sessão: o handler global atual é substituído a cada análise e
   não representa sessões simultâneas de maneira independente.
8. **Desempenho:** medir a conservância com conjuntos representativos antes de
   otimizar. Deduplicar sequências para o cálculo e reutilizar identidades, mantendo
   a multiplicidade original no denominador. Evitar carregar modelos/predizer
   topologia repetidamente para referências que não mudaram.

As propostas 6–8 incluem melhorias de arquitetura e cobertura, não novos bugs
contabilizados entre os 15 casos reproduzidos.

## Como reproduzir

O arquivo [reproduzir_2026_10_15.py](reproduzir_2026_10_15.py) usa entradas sintéticas,
não consulta serviços e apaga sua pasta temporária ao terminar. Sem `--gui` executa
os casos 01–13; com `--gui`, também 14–15. O traceback de `None.split()` no último
caso é esperado: é precisamente a falha observada pela auditoria.

As asserções registram o comportamento defeituoso desta revisão. **Não devem ser
incorporadas à suíte como requisitos do produto:** depois das correções, substituí-las
por expectativas do comportamento correto. Os testes do visualizador usam uma
simulação do objeto py3Dmol; as verificações da interface usam AppTest sem navegador.
