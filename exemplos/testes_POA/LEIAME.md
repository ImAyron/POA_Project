# Teste manual com os dados reais de `testes_POA/`

Roteiro para percorrer a interface com os arquivos que já estavam em `testes_POA/` — saída real
do BepiPred-2.0 sobre a NS1 de DENV1 e DENV2, mais a NS1 de seis flavivírus.

Os números abaixo foram medidos nesta máquina rodando o mesmo caminho sem navegador; se a tela
mostrar outra coisa, é diferença de verdade.

## Os arquivos originais e o que atrapalha

| Arquivo | Conteúdo | Convenção `Proteína_Espécie_ID` |
|---|---|---|
| `dv1_dv2.json` | BepiPred-2.0, antígenos `denv1` e `denv2` | ❌ sem `_` — precisa do adaptador da etapa 2 |
| `dv1_dv2.fasta` | as duas NS1, `denv1_ns1` / `denv2_ns1` | ❌ invertido — é `Espécie_Proteína` |
| `DENV1..4`, `YFV`, `ZIKA` `.fasta` | uma NS1 cada (352 aa), `>Dv1`, `>yfv`, … | ❌ sem `_` nenhum |

Nenhum serve como `-f` do jeito que está. Dois caminhos cobrem isso, e os dois já existem na
interface:

* **adaptador do JSON** (etapa 2) — reconstrói o `-f` a partir do próprio array `AA`;
* **checkbox "meus cabeçalhos estão em `Espécie_Proteína`"** (etapa 1) — para o `dv1_dv2.fasta`.

O `diversidade_NS1_flavivirus.fasta` desta pasta é derivado: são as seis NS1 reunidas num só
arquivo, com cabeçalho na convenção, porque o campo de comparação da etapa 4 aceita **um**
arquivo por espécie.

## Caminho A — pelo adaptador do JSON (recomendado)

1. **Etapa 1:** não envie nada. O `-f` vem do adaptador.
2. **Etapa 2:** expansor do BepiPred-2.0 → `testes_POA/dv1_dv2.json`. A interface recusa as
   chaves `denv1`/`denv2` e abre o mapeamento. Preencha **Espécie** `DENV1` / `DENV2` e
   **Proteína** `NS1` nos dois, e clique em **Adaptar e adicionar**.
   Esperado: `NS1_DENV1_ref` e `NS1_DENV2_ref`, 352 aa cada, sem resíduo `X`.
3. **Volte à etapa 1:** a tabela mostra as duas sequências, 704 resíduos no total, coluna
   *Ambíguos* em `—`.
4. **Etapa 3:** executar o POA1 → **19 epítopos** (DENV1: 11, DENV2: 8).
5. **Etapa 4, primeiro sem conjunto de comparação.** Esperado: o aviso de que DENV1 e DENV2 têm
   **1 proteína** de comparação e o limiar não discrimina nada. Rode assim mesmo e confira na
   tabela de resultado: *Identidade mín.* = `100.0` para as duas, com qualquer limiar.
6. **Etapa 4, agora com diversidade.** Envie `diversidade_NS1_flavivirus.fasta` nos **dois**
   campos (DENV1 e DENV2) e recalcule. A *Identidade mín.* cai para `23.26` (DENV1) e `20.00`
   (DENV2) — o limiar passou a ver variação de verdade.

## O limiar fazendo efeito

Mesmo epítopo de DENV1, variando só o limiar com que os CSVs são gerados:

| Limiar | 1º epítopo do CSV | 6º epítopo |
|---|---|---|
| 50 % | `83.33% (5/6)` | `100.00% (6/6)` |
| 70 % | `33.33% (2/6)` | `100.00% (6/6)` |
| 90 % | `16.67% (1/6)` | `66.67% (4/6)` |

E o cabeçalho da coluna acompanha: `Percent of protein sequence matches at identity >= 90%`.

## Etapa 5 — o funil

Dos 19 epítopos, quantos sobrevivem ao filtro `-m` (com `imin 0`, `imax 100`, sem filtro de
identidade):

| Limiar dos CSVs | `-m 0` | `-m 20` | `-m 40` | `-m 60` | `-m 80` | `-m 100` |
|---|---|---|---|---|---|---|
| 50 % | 19 | 19 | 18 | 18 | 16 | 10 |
| 70 % | 19 | 16 | 13 | 10 | 2 | 2 |
| 90 % | 19 | 7 | 5 | 3 | 0 | 0 |

A linha de 90 % com `-m 80` é o caso que vale olhar: **0 epítopos**, e a etapa 5 diz isso
*antes* de rodar o TMHMM, apontando qual filtro derrubou. No Windows o botão de executar fica
desabilitado (falta o `pyTMHMM`), mas o funil é calculado e mostrado assim mesmo.

## Caminho B — pelo FASTA, sem o adaptador

Na etapa 1, envie `testes_POA/dv1_dv2.fasta` e **marque** o checkbox
*"Meus cabeçalhos estão em `Espécie_Proteína`"*. Sem ele, `denv1_ns1` é lido como proteína
`DENV1` da espécie `NS1`, e as duas sequências caem na mesma "espécie" `NS1` — erro silencioso,
não há mensagem de erro, só um resultado errado. Vale marcar e desmarcar uma vez para ver a
coluna *Espécie* da tabela mudar.

Esse caminho dá o `-f`, mas não os epítopos: a etapa 2 ainda precisa do `dv1_dv2.json`.
