# Exemplo — o limiar de identidade fazendo efeito

Conjunto mínimo para conferir à mão, pela interface, que o limiar das etapas 4 e 5 passou a
mudar o resultado. Duas espécies de propósito: uma com conjunto de diversidade e outra sem.

| Arquivo | Onde entra |
|---|---|
| `proteinas.fasta` | etapa 1, FASTA de proteínas (`-f`) — `E_DENV1_ref` e `E_DENV2_ref` |
| `epitopos_outros.fasta` | etapa 2, expansor **"➕ Outros preditores"** — um epítopo por espécie |
| `diversidade_DENV1.fasta` | etapa 4, **"FASTA de comparação para DENV1"** — 6 variantes |

As 6 variantes de DENV1 têm identidade 100, 100, 87.5, 75, 62.5 e 62.5 % contra o epítopo
`AYIAMKGQ`. O resultado da etapa 4 muda com o limiar:

| Limiar | DENV1 (com diversidade) | DENV2 (sem diversidade) |
|---|---|---|
| 50 % | `100.00% (6/6)` | `100.00% (1/1)` |
| 70 % | `66.67% (4/6)` | `100.00% (1/1)` |
| 90 % | `33.33% (2/6)` | `100.00% (1/1)` |

**DENV2 é o caso que parecia "o limiar não está sendo respeitado".** Sem um conjunto de
diversidade, a única proteína de comparação é aquela de onde o epítopo foi predito: ele casa
com ela a 100 % e a conservância sai 100 % para qualquer limiar. A etapa 4 agora avisa isso em
vez de deixar passar calado.

No POA2, com `-m 60` e `-imin 60`: a `t=70` sobrevivem os dois epítopos; a `t=90` só o de
DENV2 — o de DENV1 cai, porque 33.33 % < 60 %.
