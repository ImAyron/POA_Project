# Como testar o POA

Guia para conferir, por conta própria, que o pipeline está funcionando — e em particular que o
**limiar de identidade** das etapas 4 e 5 faz efeito, e que cada análise fica na sua própria pasta.

Tudo aqui foi executado nesta máquina (Windows 11, Python 3.12); os números das tabelas são
medidos, não estimados.

---

## 1. Preparar o ambiente (uma vez só)

```powershell
python -m venv .venv
.venv\Scripts\python.exe -m pip install -r requirements.txt
```

> O `pyTMHMM` é instalado separadamente, conforme `TESTING_WSL.md`, para evitar que seu build
> nativo bloqueie a instalação básica. Sem ele, as etapas 1 a 4 e o painel de filtros da etapa 5
> funcionam; o botão de executar a topologia fica desabilitado.

Confira:

```powershell
.venv\Scripts\python.exe -c "import pandas, Bio, streamlit; print('ok')"
```

---

## 2. Testes automatizados

```powershell
.venv\Scripts\python.exe -m pytest -q
```

Esperado nesta revisão: **155 passando**. Se algum falhar, o nome do teste já diz o que quebrou.

Para rodar só o que cobre as mudanças recentes:

```powershell
.venv\Scripts\python.exe -m pytest -q -k "conservancy or threshold or run"
```

Os testes de interface dirigem o Streamlit sem navegador, incluindo os filtros da **etapa 5**
com e sem pyTMHMM disponível.

---

## 3. Provar que as correções são regressões reais

Os casos de regressão usam entradas que reproduzem falhas específicas e conferem o resultado
esperado. Execute-os sem substituir arquivos da sua cópia de trabalho:

```powershell
.venv\Scripts\python.exe -m pytest tests/test_review_regressions.py tests/test_step_outputs.py tests/test_conservancy_local.py -q
```

Para comparar com uma versão antiga, use uma cópia isolada do repositório. Os casos novos e
os limites da validação estão descritos em `REVISAO.md`.

---

## 4. Teste manual pela interface

```powershell
.venv\Scripts\streamlit.exe run poa\gui\app.py
```

Abre em `http://localhost:8501`. O conjunto de exemplo está em [`exemplos/limiar/`](exemplos/limiar/).

### 4.1 Pastas por análise

Ao abrir, a interface pergunta qual análise usar.

1. Crie uma com um nome qualquer → deve surgir `results/<DDMMAAAA>-<NOME>/`.
2. Clique em **"Trocar de análise"** e crie uma segunda → pasta separada.
3. Em **"Abrir existente"**, as duas aparecem com data e progresso, mais recente primeiro.
4. Reabra a primeira: o que você tinha carregado nela continua no disco, e o estado da
   segunda **não** vaza para ela.

### 4.2 Percorrer o pipeline

| Etapa | O que enviar |
|---|---|
| 1 | `exemplos/limiar/proteinas.fasta` → 2 sequências, espécies DENV1 e DENV2 |
| 2 | expansor **"➕ Outros preditores"** → `exemplos/limiar/epitopos_outros.fasta` |
| 3 | executar o POA1 → 2 epítopos, um por espécie |

### 4.3 O limiar (o teste principal)

Na etapa 4, **primeiro rode sem enviar nada** no campo de comparação. Esperado: o aviso de que
DENV1 e DENV2 têm só 1 proteína de comparação e o limiar não discrimina nada.

Agora envie `exemplos/limiar/diversidade_DENV1.fasta` em **"FASTA de comparação para DENV1"** e
recalcule variando o limiar. Abra `results/<sua análise>/conservancy_csv/DENV1_conservancy.csv`:

| Limiar | DENV1 (com diversidade) | DENV2 (sem diversidade) |
|---|---|---|
| 50 % | `100.00% (6/6)` | `100.00% (1/1)` |
| 70 % | `66.67% (4/6)` | `100.00% (1/1)` |
| 90 % | `33.33% (2/6)` | `100.00% (1/1)` |

Confira também:

* o **cabeçalho** do CSV traz o limiar usado — `Percent of protein sequence matches at identity >= 70%`;
* apareceu um **`conservancy_meta.json`** ao lado dos CSVs, com limiar, critério e o tamanho do
  conjunto de comparação de cada espécie.

A coluna DENV2 é o comportamento que parecia "o limiar não está sendo respeitado": sem conjunto
de diversidade, a única proteína de comparação é aquela de onde o epítopo foi predito, então ele
casa a 100 % e a conservância sai 100 % com qualquer limiar. Isso continua acontecendo — é
matemática, não bug —, mas agora a interface avisa em vez de deixar passar calado.

---

## 5. Executar a topologia com pyTMHMM

A etapa 5 mostra os controles de limiar e o funil dos filtros mesmo sem `pyTMHMM`. No Windows,
você pode ajustar os parâmetros e verificar os CSVs; a dependência bloqueia apenas a execução.

Esse comportamento está coberto por `tests/test_gui_app_threshold.py` e
`tests/test_step_outputs.py`, com a disponibilidade do preditor simulada.

Para executar a topologia, use a WSL — veja [`TESTING_WSL.md`](TESTING_WSL.md). Resumo do caminho:

```powershell
wsl -l -v              # se só aparecer 'docker-desktop', você não tem distro utilizável
wsl --install -d Ubuntu
```

Depois, dentro do Ubuntu, siga o `TESTING_WSL.md` §0 a §2 (miniconda, ambiente `poa`, pyTMHMM com
`numpy<2` e `--no-build-isolation`). Para continuar na WSL uma análise começada no Windows, rode a
interface direto sobre a pasta do Windows (§6.1) ou exporte `POA_RESULTS_DIR` apontando para ela
(§1) — senão a WSL cria a sua própria `results/` e as análises não se encontram. Na tela inicial,
**abra a mesma análise**; a etapa 5 então oferece "Retomar com estes arquivos".

---

## 6. Testar com dados reais seus

Os testes provam que o limiar passou a ter efeito; eles **não** dizem se os epítopos que sobram
nos seus dados fazem sentido biológico. Para isso:

1. Crie uma análise nova com um nome que identifique o teste.
2. Rode o fluxo com os seus arquivos.
3. Na etapa 4, confira o aviso sobre o tamanho do conjunto de comparação de cada espécie. Se
   aparecer, o resultado de conservância daquela espécie não significa nada ainda — falta o
   conjunto de diversidade (*world set*), não é defeito do código.
4. Compare com a análise anterior: como cada uma tem a sua pasta, as duas coexistem.
