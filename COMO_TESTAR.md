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

> O `pyTMHMM` é o **único** pacote instalado separadamente (conforme `TESTING_WSL.md`), porque seu
> build nativo bloquearia a instalação básica no Windows. Sem ele, as etapas 1 a 4 e o painel de
> filtros da etapa 5 funcionam; o botão de executar a topologia fica desabilitado. Todo o resto —
> inclusive o `py3Dmol` da visualização 3D — vem no `requirements.txt`.

Confira:

```powershell
.venv\Scripts\python.exe -c "import pandas, Bio, streamlit, plotly, py3Dmol; print('ok')"
```

Se esse comando falhar em `py3Dmol`, o ambiente foi montado antes de o visualizador 3D existir ou
pacote a pacote; rode o `pip install -r requirements.txt` de novo. As etapas 1 a 6 funcionam sem
ele — só a estrutura 3D da etapa 7 não renderiza.

---

## 2. Testes automatizados

```powershell
.venv\Scripts\python.exe -m pytest -q
```

Esperado nesta revisão: **254 passando**. Se algum falhar, o nome do teste já diz o que quebrou.

Para rodar só o que cobre as mudanças recentes:

```powershell
.venv\Scripts\python.exe -m pytest -q -k "structure_map or dbsearch"      # 36, POA0 + guarda 3D
.venv\Scripts\python.exe -m pytest -q -k "conservancy or threshold or run"
```

Nenhum teste da suíte acessa a rede: as respostas do UniProt, do InterPro e do PDBe são
simuladas, com cópias reduzidas das respostas reais. Por isso a suíte roda em ~11 s e funciona
offline. Os testes do §8 abaixo é que exercitam os bancos de verdade.

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

---

## 7. Ler o log de uma execução

Toda análise grava um `poa.log` na própria pasta, sempre em `DEBUG` — mesmo quando o terminal está
em `INFO`. É o primeiro lugar a olhar quando um resultado não bate com o esperado, porque registra
o que cada etapa produziu, e não apenas que ela rodou.

As linhas têm a forma `<TAG> <etapa> | <mensagem>`, então dá para filtrar por parte do pipeline:

```powershell
# o que cada método do POA1 selecionou
Select-String 'DONE POA1/' "results\08102026-MEU_TESTE\poa.log"

# tudo que não terminou
Select-String 'FAIL' "results\08102026-MEU_TESTE\poa.log"

# a etapa POA2 inteira, sub-etapas incluídas
Select-String 'POA2/' "results\08102026-MEU_TESTE\poa.log"
```

Três verificações que valem fazer no log antes de confiar num resultado:

1. **`DONE POA1 | epitopes=N`** — confira se `N` é o número que você esperava. Zero aqui significa
   que os filtros de tamanho da etapa 1 descartaram tudo.
2. **`species=`** nas linhas `DONE` — se aparecer uma espécie que você não reconhece, ou faltar
   uma que deveria estar, é cabeçalho fora da convenção `Proteína_Espécie_ID`.
3. **Linhas `WARN`** — cada uma marca um resultado que ficou suspeito, não um erro. A mais
   importante é a do conjunto de comparação com menos de duas proteínas: a conservância daquela
   espécie sai 100% com qualquer limiar.

Para ver o detalhe por item também no terminal, suba o nível:

```powershell
$env:POA_LOG_LEVEL = "DEBUG"
streamlit run poa/gui/app.py
```

Na linha de comando, `POA_LOG_FILE` aponta o arquivo, já que não há pasta de análise:

```powershell
$env:POA_LOG_FILE = "poa1.log"
python POA1_v1.0.py -f proteins.fasta -d results -b3 bepipred3.fasta
```

---

## 8. Testar o POA0 e a verificação da estrutura 3D

Esta parte precisa de **internet** (UniProt, InterPro, PDBe, RCSB). Nada aqui grava dentro de
`results/` — escolha uma pasta de rascunho.

### 8.1 Da espécie ao FASTA, para qualquer vírus

As quatro perguntas, cada uma mostrando as opções e encerrando:

```powershell
.venv\Scripts\python.exe POA0_v1.0.py -s "Oropouche virus"                 # qual táxon?
.venv\Scripts\python.exe POA0_v1.0.py -s "Oropouche virus" -l proteinas    # qual proteína?
.venv\Scripts\python.exe POA0_v1.0.py -a A0A0D4BSW3 -l regioes             # qual região?
.venv\Scripts\python.exe POA0_v1.0.py -a A0A0D4BSW3 -l estruturas          # qual estrutura?
```

O que esperar, medido nesta máquina:

| Comando | Resultado |
|---|---|
| `-s "Oropouche virus"` | táxon **118655**, e `Oropouche like virus` como segunda opção |
| `-l proteinas` | **A0A0D4BSW3** em primeiro — 1.420 aa, TrEMBL, a única com PDB (6H3X) |
| `-a A0A0D4BSW3 -l regioes` | 3 regiões: `Domain+Pfam` 27–306, `PDB` 482–703, `Domain+Pfam` 514–1359 |

A coluna de fonte mostra `Domain+Pfam` quando duas fontes descrevem o mesmo intervalo — num
registro TrEMBL o `Domain` do UniProt é derivado do Pfam, então as duas concordam.

Troque o nome por qualquer outro vírus. `-s "Zika virus"` dá o táxon 64320 e, nas regiões, as
**12 proteínas maduras** anotadas como `Chain` (é entrada SwissProt), entre elas
`Envelope protein E` em 291–790.

### 8.2 Gerar os arquivos

```powershell
.venv\Scripts\python.exe POA0_v1.0.py -a A0A0D4BSW3 -r 482-703 --nome GcHead `
                                      --sigla OROV --estrutura 6H3X -d results\teste_poa0
```

Confira três coisas em `results\teste_poa0\`:

1. `GcHead_OROV.fasta` — **222 aa**, cabeçalho `>GcHead_OROV_A0A0D4BSW3` (três campos:
   proteína, espécie, acessão).
2. `structure_6H3X.pdb` — baixado do RCSB.
3. `POA0_relatorio.json` — procedência. O campo `estrutura.sifts.pdb_start` deve ser **482**
   (numeração do autor) e `entity_start` deve ser **1** (índice da entidade). São números
   diferentes para o mesmo resíduo; confundi-los desloca todo epítopo em 481 posições.

**Regressão:** pedir a região pelo id PDB dá o mesmo arquivo que pedir por coordenadas.

```powershell
.venv\Scripts\python.exe POA0_v1.0.py -a A0A0D4BSW3 -r 6H3X --nome GcHead --sigla OROV -d results\teste_pdbid
Compare-Object (Get-Content results\teste_poa0\GcHead_OROV.fasta) (Get-Content results\teste_pdbid\GcHead_OROV.fasta)
```

O `Compare-Object` não deve imprimir nada — os arquivos são idênticos.

### 8.3 A verificação da estrutura, fora da interface

`tools\check_structure.py` faz a mesma conferência da etapa 7 em linha de comando. Sai com
código **0** quando a estrutura serve e **1** quando não serve.

**O caso que deve passar:**

```powershell
.venv\Scripts\python.exe tools\check_structure.py -p results\teste_poa0\structure_6H3X.pdb `
                                                  -f results\teste_poa0\GcHead_OROV.fasta -c A
```

Esperado: `POLYPROTEIN · OROPOUCHE VIRUS`, **100,0% de identidade**, 221 resíduos observados,
`RESULTADO: a estrutura corresponde`.

Vem junto um aviso de "1 epítopo parcialmente ausente no cristal", e ele está certo: sem `-e` a
ferramenta confere a referência inteira (1–222) como se fosse um epítopo só, e o resíduo **703**
(o último) não foi resolvido na estrutura — a cadeia A tem 221 dos 222. É exatamente o tipo de
detalhe que a figura sozinha não mostra.

**O caso que deve falhar** — e este você já tem no repositório, sem precisar de internet:

```powershell
.venv\Scripts\python.exe tools\check_structure.py `
    -p "results\08102026-TESTANDO\inputs\structure_4C23.pdb" `
    -f "results\08102026-TESTANDO\inputs\ref_0_DENV1.fasta" `
    -e "results\08102026-TESTANDO\poa1_out\Conservancy Analysis\DENV1_epitopes.fasta" -c A
```

Esperado: `L-FUCULOSE KINASE FUCK · STREPTOCOCCUS PNEUMONIAE`, **25,0% de identidade**,
**0 de 11 epítopos mapeados**, `RESULTADO: a estrutura NAO corresponde`.

Esse é o defeito que motivou a verificação. O 4C23 é uma quinase bacteriana de açúcar, sem
relação com dengue, e foi usado como estrutura em duas análises de DENV1/DENV2. Nenhuma das 19
sequências de epítopo existe nele, mas as 19 posições caem dentro da numeração da cadeia A (que
vai de −5 a 467) — então, antes desta revisão, as 19 eram pintadas de vermelho sem erro algum.

### 8.4 O POA0 pela interface (etapa 1)

```powershell
.venv\Scripts\streamlit.exe run poa\gui\app.py
```

Crie uma análise nova e abra a **etapa 1**. No alto, o expansor
**"🔎 Não tenho o FASTA — buscar pela espécie (POA0)"**. Dentro dele, as mesmas quatro perguntas
da linha de comando, uma destravando a seguinte:

1. Digite `Oropouche virus` → **1. Buscar organismo** → deve listar o táxon **118655**.
2. **Listar proteínas deste táxon** → **A0A0D4BSW3** em primeiro, com `PDB: 6H3X`.
3. **Listar regiões desta proteína** → três regiões, e um seletor de estrutura aparece com
   `6H3X  UniProt 482-703  15.6% da sequência  2.09 A  cadeia(s) A`.
4. Escolha a região, ajuste os campos do cabeçalho (ex.: `GcHead` e `OROV`), escolha a estrutura
   e clique **⬇️ Gerar o FASTA**.

O que deve acontecer: o arquivo aparece em `inputs/`, a tabela de cabeçalhos logo abaixo passa a
mostrá-lo com espécie `OROV`, e a estrutura baixada já fica carregada na etapa 7. Nenhuma
dessas telas pede upload — é o caminho de quem não tem o FASTA.

> Deixe a caixa **"Meus cabeçalhos estão em `Espécie_Proteína`"** desmarcada: o POA0 já grava na
> convenção certa, e inverter os campos corromperia o que ele escreveu.

O painel é coberto por `tests/test_gui_poa0.py` (5 testes, com os bancos simulados), inclusive o
caso de o UniProt estar fora do ar — a interface mostra o erro e continua utilizável.

### 8.5 A verificação da estrutura pela interface (etapa 7)

```powershell
.venv\Scripts\streamlit.exe run poa\gui\app.py
```

Abra a análise `08102026-TESTANDO` (ou rode o fluxo das etapas 1 a 3 com os seus dados) e vá à
etapa 7. Envie `results\08102026-TESTANDO\inputs\structure_4C23.pdb`. Deve acontecer, nesta ordem:

1. Uma faixa azul dizendo o que o arquivo contém: **L-FUCULOSE KINASE FUCK · STREPTOCOCCUS
   PNEUMONIAE · X-RAY DIFFRACTION · 2.00 ANGSTROMS**.
2. Uma faixa vermelha de **estrutura incompatível**, com a identidade medida.
3. A estrutura desenhada **em cinza, sem nenhum destaque** — ela continua visível, o que é
   recusado é a afirmação de que aqueles resíduos são os epítopos.

Depois envie `results\teste_poa0\structure_6H3X.pdb` numa análise cuja referência seja
`GcHead_OROV.fasta`: aí sim os epítopos aparecem em vermelho, com uma faixa verde dizendo
quantos foram mapeados, a identidade e por qual método.

> A etapa 7 precisa de internet para desenhar (o 3Dmol.js vem de um CDN) e do `py3Dmol`
> instalado. A verificação em si não precisa de nenhum dos dois — por isso o
> `tools\check_structure.py` funciona offline.

---

## 9. Validar a conservância local contra o IEDB

A reimplementação local reproduz o método de Bui et al. (2007) que o IEDB aplica, mas o IEDB
**não** expõe essa ferramenta por API (reconfirmado em 12/10/2026). Então a equivalência não pode
ser verificada automaticamente: ela depende de uma exportação que você faz uma vez, à mão, e
guarda como evidência.

**É o passo que separa "o pipeline roda" de "o número é defensável em publicação.**

1. Rode a etapa 4 localmente e guarde o FASTA de epítopos e o de comparação.
2. No [Epitope Conservancy Analysis](http://tools.iedb.org/conservancy/), submeta **os mesmos**
   epítopos, **o mesmo** conjunto de proteínas e **o mesmo** limiar/operador.
3. Baixe o CSV e compare:

```powershell
.venv\Scripts\python.exe tools\validate_conservancy.py `
    -i iedb_export.csv `
    -e "results\<analise>\poa1_out\Conservancy Analysis\<SP>_epitopes.fasta" `
    -p results\<analise>\inputs\comparison_<SP>.fasta `
    -t 70 --operador ">="
```

A ferramenta recalcula localmente e compara os três campos numéricos de cada epítopo (percentual,
`Minimum identity`, `Maximum identity`), com tolerância de 0,011 pp — o bastante para absorver o
arredondamento de 2 casas do IEDB e nada mais. Sai com código **0** se tudo bate e **1** se diverge.

Se divergir, o relatório diz em qual epítopo e em qual campo, com o delta. Duas causas prováveis
antes de suspeitar do algoritmo:

* **limiar ou operador diferentes** entre o site e o comando — o cabeçalho das duas colunas é
  impresso no alto justamente para essa conferência. O IEDB escreve `<=` onde o POA escreve `<`;
* **conjunto de comparação diferente** — epítopos "só no IEDB" ou "só no local" aparecem listados.

Guarde o CSV do IEDB junto da análise. Ele é a evidência de que, **naquele conjunto de dados**,
os números locais podem substituir os do site.

### 9.1 O algoritmo, sem precisar do site

Independente da validação acima, 16 testes fixam o método contra valores calculados à mão —
janela sem gaps, máximo sobre os offsets, a fronteira exata do limiar e a semântica das colunas
de identidade mínima e máxima:

```powershell
.venv\Scripts\python.exe -m pytest tests\test_conservancy_algorithm.py -q
```

Se algum deles quebrar, o cálculo mudou — e qualquer resultado publicado com a versão anterior
deixa de valer.

---

## 10. Apagar uma análise

Na tela inicial, aba **"Abrir existente"**, selecione a análise e abra o expansor
**🗑️ Apagar esta análise**. Ele mostra quantos arquivos e quantos MB serão removidos.

A confirmação é **digitar o nome da pasta**, não um segundo botão. Um par de botões vira reflexo;
digitar o nome é uma decisão — e torna impossível apagar uma análise diferente daquela que você
está lendo na tela. O botão só habilita quando o texto confere exatamente.

**Não tem desfazer e nada vai para a lixeira.** Some a pasta inteira: predições (que custaram
chamadas de API e minutos), CSVs de conservância, planilhas do POA2 e o log.

O backend recusa quatro alvos mesmo que chamado com os argumentos errados:

| Alvo | Por quê |
|---|---|
| A raiz `results/` | só uma análise por vez |
| Qualquer caminho fora da raiz | um argumento errado não pode sair da árvore de resultados |
| A pasta `cache/` | é compartilhada por todas as análises; apagá-la custaria horas de API |
| Pasta que não é análise (sem `run.json` nem as subpastas) | `results/` também guarda saídas do POA0 e rascunhos |

Para conferir as salvaguardas:

```powershell
.venv\Scripts\python.exe -m pytest tests\test_run_delete.py -q
```

São 15 testes, dos quais 9 verificam recusas — a parte que importa numa operação irreversível.
