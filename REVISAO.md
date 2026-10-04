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
