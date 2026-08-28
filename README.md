# POA: Pipeline de Otimização de Antígenos (Antigen Optimization Pipeline)

Esta é a versão em **português (PT-BR)**, a principal deste repositório.
A versão em inglês está [ao final do documento](#english).

## 0. AUTOMAÇÃO E INTERFACE GRÁFICA (v1.1)

> Esta versão acrescenta uma camada opcional de automação e uma interface gráfica **sem alterar a
> lógica científica de ranqueamento e seleção de epítopos**. A linha de comando original continua
> funcionando exatamente como antes. O projeto completo e a estratégia de automação de cada
> ferramenta estão documentados em [`ARCHITECTURE.md`](ARCHITECTURE.md). O relatório detalhado de
> tudo que mudou está em [`MUDANCAS.md`](MUDANCAS.md) — também disponível como
> [`MUDANCAS.pdf`](MUDANCAS.pdf).

**Organização do pacote** (refatorado em camadas; os scripts da raiz tornaram-se pontos de entrada
finos): `poa/core` (parsing + lógica do pipeline, sem rede) · `poa/services` (ferramentas externas) ·
`poa/cli` (POA1/POA2, mesmos argumentos) · `poa/gui` (aplicação Streamlit).

**Instalação**
```
python -m venv .venv
# Windows: .venv\Scripts\activate   |   Linux/macOS: source .venv/bin/activate
pip install -r requirements.txt
```
Extras que não são pacotes pip ou são pesados (instale apenas se quiser aquela automação):
* **EMBOSS `antigenic`** (substitui o site PAP/IMED, que saiu do ar; mesmo método Kolaskar–Tongaonkar):
  `conda install -c bioconda emboss` (no Windows: WSL2 + conda, ou Docker).
* Modelo local do **BepiPred-3.0** (PyTorch + ESM-2): `pip install bp3`.
* **Automação de navegador** para NetCTL 1.2 / BepiPred-2.0: `pip install playwright && playwright install chromium`.

**Executar a interface gráfica** (fluxo guiado de 7 etapas: proteínas → predições → POA1 →
conservação → POA2 → resultados → visualização):
```
streamlit run poa/gui/app.py      # ou:  python run_gui.py
```

**Executar pela linha de comando** (inalterada):
```
python POA1_v1.0.py -f proteins.fasta -d results -b3 bepipred3.fasta -e y
python POA2_v1.0.py -g True -t 70 -d conservancy_csvs -f proteins.fasta -rf 1
```

**Arquivos originais da pesquisa** (adaptador de nomenclatura para os JSONs do BepiPred-2.0 do fluxo
manual, que trazem a chave de antígeno genérica `Sequence`):
```
python -m poa.cli.import_realdata --pred-dir <5_...> --world-dir <2_...> \
    --out ./realdata_run --protein E --chikv-protein E1 --threshold 70 --run
```

**Status da automação por ferramenta** (verificado em 2026 — detalhes em `ARCHITECTURE.md`):

| Ferramenta | Estratégia | Fallback |
|---|---|---|
| MHC-II Binding (IEDB) | API REST (síncrona) | upload manual |
| BepiPred-3.0 | pacote local `bp3` | upload manual |
| PAP/IMED | EMBOSS `antigenic` local | upload manual |
| Epitope Conservancy (IEDB) | **reimplementação local** (não existe API nem standalone) | upload manual do CSV |
| NetCTL 1.2 / BepiPred-2.0 | automação de navegador (*best effort*) | upload manual |
| pyTMHMM | local (inalterado) | — |

Em todos os casos, se um serviço estiver indisponível o pipeline registra uma mensagem clara e
permite o upload manual do arquivo de resultado, preservando o comportamento semiautomático original.

**Visualização** (etapa 7 da interface): mapa 2D dos epítopos ao longo da sequência — uma trilha por
proteína/espécie, colorida por método de predição — e visualizador 3D sobre uma estrutura `.pdb`
enviada pelo usuário, com os epítopos destacados (py3Dmol / 3Dmol.js; a renderização 3D requer
internet, e a numeração de resíduos do PDB precisa corresponder às posições das predições).

**Testes:** `pytest` — 82 testes cobrindo parsers, ranqueamento/filtros, conservação local, clientes
de serviço, backend da interface, visualização e o adaptador dos dados reais.

> **Corrigido:** o parser do NetCTL 1.2 tinha um bug **pré-existente** (inserção de linhas por índice
> não sequencial, incompatível com pandas ≥ 2; e lógica de marcador/preenchimento inconsistente com a
> seleção `<-E` dos epítopos, de modo que nenhum epítopo do NetCTL era selecionado). Foi corrigido
> — com aprovação explícita do usuário, por alterar quais linhas são selecionadas — para usar índice
> sequencial e preencher apenas as linhas curtas até a contagem de colunas, mantendo a seleção `<-E`
> documentada neste README. Coberto por testes unitários.

> **Análise de várias espécies na mesma execução (corrigido):** a interface guardava um único
> arquivo por método, então ao testar duas espécies ao mesmo tempo a segunda substituía a primeira
> silenciosamente e o `-f` continuava sendo a referência da primeira; além disso o parser do
> BepiPred-3.0 lia o cabeçalho como `Espécie_Proteína_ID` (invertido em relação à convenção
> `Proteína_Espécie_ID` usada por todo o resto do pipeline), e a análise de conservância comparava
> todas as espécies com o mesmo conjunto de proteínas. Agora: cada método aceita **um arquivo por
> espécie** (unidos automaticamente antes do POA1), as referências são unidas em um único `-f`,
> cada espécie é comparada com o conjunto de proteínas dela na etapa 4, e a etapa 3 exibe um
> diagnóstico das espécies identificadas — incluindo as que aparecem nos epítopos mas faltam no
> `-f`. Coberto por `tests/test_multispecies.py`.

> **Cabeçalhos fora da convenção (corrigido):** se o JSON do BepiPred-2.0 traz várias espécies em
> um só arquivo (chaves como `denv1`/`denv2`, sem proteína nem ID) ou se o `-f` usa
> `Espécie_Proteína` (`denv1_ns1`), nenhuma espécie era extraída e todos os epítopos saíam sem
> espécie. Agora o adaptador da etapa 2 mapeia **cada antígeno** para sua espécie e proteína, e a
> etapa 1 mostra como cada cabeçalho foi interpretado, com opção de inverter os dois primeiros
> campos.

> **Saída da interface:** a GUI grava tudo em `POA_Project/results/` (antes: um diretório temporário
> por sessão). Use a variável `POA_RESULTS_DIR` para apontar outra pasta. O caminho aparece na barra
> lateral; como é fixo, uma nova execução sobrescreve arquivos de mesmo nome.

> **POA2 no Windows:** `pip install pyTMHMM` não funciona (sem *wheel* para Windows, precisa do
> Microsoft C++ Build Tools e o pacote não compila contra numpy 2.x). Rode a interface pela WSL,
> apontando para a mesma pasta do projeto — veja `TESTING_WSL.md` §6.1. A etapa 5 oferece
> **"Retomar com estes arquivos"** a partir dos CSVs já em `results/`, sem repetir as etapas 1 a 4.

> **Corrigido (topologia do POA2):** `tmhmmAnalysis` preenchia as colunas `Portion_*` com um valor
> por *registro casado* do `-f`, não por epítopo. Um `-f` com registros redundantes fazia cada
> epítopo casar duas vezes (o teste de espécie era um *substring*: `DENV1` casa com `DENV1_NS1` e
> com `E_DENV1_REF`) e a atribuição falhava; pior, com casamentos a mais e a menos no mesmo lote os
> totais podiam bater e as colunas eram gravadas deslocadas, atribuindo a topologia de um epítopo a
> outro. Agora a comparação usa o campo de espécie do cabeçalho e há exatamente um valor por linha.

## 1. INTRODUÇÃO

O Pipeline de Otimização de Antígenos (POA) é um pipeline semiautomático em python para análise, organização, recuperação e triagem de epítopos resultantes da predição de peptídeos antigênicos em ferramentas web. 

A partir dos resultados obtidos em ferramentas de predição para células B, células T-citotóxicas e células T-auxiliares, o POA identifica os melhores alvos e cria um banco de dados com estes antígenos. Além disso, o POA auxilia na identificação dos peptídeos com alta similaridade de resíduos entre espécies de um determinado grupo de organismos, com o auxílio da ferramenta [Epitope Conservancy Analysis](http://tools.iedb.org/conservancy/), e na caracterização da localização dos epítopos em porções transmembranares, com o auxílio da ferramenta [pyTMHMM](https://github.com/bosborne/pyTMHMM).


### 1.1 Predição de Peptídeos Antigênicos

O screening de alto rendimento, o estudo de sequências genômicas de organismos patogênicos e a aplicação de técnicas de vacinologia reversa têm possibilitado a construção de vacinas e testes diagnósticos de forma rápida e eficiente. Neste contexto, o uso de abordagens in silico que integram diversas técnicas e conhecimentos, como a construção de multi-epitopos, tem se tornando cada vez mais presente na produção de produtos biotecnológicos e intervenções imuno-terapêuticas contra muitas doenças virais, bacterianas e parasitárias ([Dar et al., 2021](https://www.nature.com/articles/s41598-021-90868-2); [Enayatkhani et al., 2021](https://www.tandfonline.com/doi/full/10.1080/07391102.2020.1756411?casa_token=Ctu83G9J4aEAAAAA%3Aocl_27U5qwcFn-WXSFZIaaSmChfq-hC5ggb281Z1hGuocFuOXH65acPngTeO6HavrexguLxkdEbfAVQ); [Shey et al., 2019](https://www.nature.com/articles/s41598-019-40833-x)).


Com o POA você organiza os resultados das predições in silico de epítopos de proteínas imunogênicas, obtém uma melhor visualização destes resultados (podendo ser um enorme volume de dados) e otimiza o tempo de análise destas predições. Além disso, direciona a análise fazendo a triagem de epítopos únicos (ou compartilhados) num grupo de organismos e epítopos expostos na estrutura membranar.


Para epítopos de células B, o POA possui suporte para análise do resultado de predição da ferramenta [Bepipred 2.0](https://services.healthtech.dtu.dk/service.php?BepiPred-2.0) e [Predicted Antigenic Peptides/IMED](http://imed.med.ucm.es/Tools/antigenic.pl).
Para epítopos de células T-citotóxica, o POA possui suporte para análise do resultado de predição da ferramenta [NetCTL - 1.2](https://services.healthtech.dtu.dk/service.php?NetCTL-1.2).
Para epítopos de células T-auxiliares, o POA possui suporte para análise do resultado de predição da ferramenta [MHC-II Binding Predictions](http://tools.iedb.org/mhcii/).
Resultados de outros preditores podem também ser analisados pelo POA quando previamente organizados em arquivos fasta padronizados (ver seção 3.1.1 - Arquivos de entrada no POA1).



## 2. INSTALAÇÃO

No momento há apenas uma maneira de instalar o POA v1.0, instalando-o manualmente. Para que ele funcione é necessário que todas as dependências também sejam instaladas (instalação manual).
As dependências necessárias estão listadas abaixo:

Python version = 3.7.7

Pandas version 1.3.4

Numpy version 1.21.4

Biopython version 1.79

pyTMHMM


Depois que as dependências forem instaladas, instale o POA v1.0 (Linux) com: git clone https://github.com/UbiratanBatista/POA_Project.git



## 3. USO DO POA


Por se tratar de um pipeline semi-automático, algumas das etapas ainda precisarão ser feitas de forma manual pelo usuário.


A primeira parte do pipeline (POA1) fará a análise dos resultados de predição e criará um banco de dados de antígenos, com base nos valores de rankeamento estabelecidos por cada ferramenta de predição. Ao final, o POA1 organiza os resultados em múltiplos arquivos fasta para serem utilizados na etapa seguinte.


Os arquivos resultantes da primeira fase do POA (POA1) deverão ser submetidos na plataforma web da ferramenta [Epitope Conservancy Analysis](http://tools.iedb.org/conservancy/). Os resultados da análise de conservação deverão ser organizados em uma pasta e esta será submetida a segunda parte do pipeline (POA2).


O POA2 analisará os epítopos de maior ou menor conservação (a depender do objetivo da análise) e aplicará sobre estes epítopos a análise de predição de estruturas transmembranares, TMHMM. Por fim retornará os epítopos resultantes da análise de conservação e sua caracterização quanto à posição na membrana (externa, transmembranar e interna).


### 3.1 POA1


O POA1 recebe os resultados da predição das ferramentas web em diferentes formatos e utiliza um algoritmo que seleciona os peptídeos de interesse e os organiza em matrizes padronizadas. 


#### 3.1.1 Arquivos de entrada (Arquivos dos resultados de predição)


Como pipeline semi-automático, a primeira coisa a se fazer é realizar a predição nas ferramentas web de predição de epítopos (ver seção 1.1). Para isto deve-se submeter, via web, as sequências das proteínas de interesse nestas ferramentas de predição e seguir os protocolos e parâmetros definidos para cada método e objetivo. O POA1 recebe os diferentes formatos de arquivos de saída de cada método e executa a análise, indicando e organizando os peptídeos preditamente imunogênicos.


###### 3.1.1.1 Bepipred


O arquivo fasta contendo todas as proteínas para a análise deve ser submetido à ferramenta Bepipred. O resultado esperado é algo como isto:

![Captura de tela de 2022-01-26 17-45-11](https://user-images.githubusercontent.com/72517648/151360191-20879962-26e7-4df5-a79e-1f9d79082c8d.png)


O resultado da análise do Bepipred 2.0 é uma página web (html) contendo as proteínas e as regiões preditamente antigênicas (marcadas em E) em cada sequência. Na página de resultados tem a opção de fazer o Download deste resultado em arquivo .json (JSON Summary). O arquivo json contendo toda esta informação será utilizado pelo POA1.


###### 3.1.1.2 Predicting Antigenic Peptides/IMED


Diferentemente do Bepipred, a ferramenta do IMED não recebe mais que uma proteína por vez. Neste caso, as sequências das proteínas deverão ser submetidas uma por vez, definindo sempre o cabeçalho. O resultado da predição é algo como: 

![Captura de tela de 2022-01-26 17-56-13](https://user-images.githubusercontent.com/72517648/151360188-21f130dd-ae97-4067-b7c9-8f138bff95bf.png)


Esta tabela deve ser copiada (ctrl+c) e colada (ctrl+v) em um arquivo de texto (.txt), seguindo o seguinte formato:

![Captura de tela de 2022-01-26 17-58-37](https://user-images.githubusercontent.com/72517648/151360185-4a3b4863-5eab-4024-99d9-5105ac705e53.png)


O arquivo .txt formatado e contendo os epítopos preditos pelo Predicting Antigenic Peptides/IMED será utilizado pelo POA1.


###### 3.1.1.3 NetCTL


O arquivo fasta contendo todas as proteínas para a análise deve ser submetido à ferramenta NetCTL 1.2. O resultado esperado é algo como isto:

![Captura de tela de 2022-01-26 18-02-17](https://user-images.githubusercontent.com/72517648/151360184-6091362a-bcb6-455f-b363-64e575096806.png)


O resultado de saída da análise do NetCTL 1.2 é uma página web (html) contendo os peptídeos resultantes da análise e aqueles preditamente antigênicos (indicados com < - E). Esta página de resultados deverá ser salva em formato .html, O download pode ser feito com a opção “salvar como” através do click do botão direito do mouse sobre a página. ATENÇÃO: Aguardar a página carregar toda a informação para depois efetuar o Download. O arquivo .html contendo toda esta informação será utilizado pelo POA1.


###### 3.1.1.4 MHC-II Binding Predictions


A ferramenta MHC-II Binding Predictions não avalia separadamente as proteínas quando submetidas ao servidor em um único arquivo .fasta. As proteínas deverão ser submetidas separadamente para a análise. O resultado será algo parecido com isto:

![Captura de tela de 2022-01-26 18-07-50](https://user-images.githubusercontent.com/72517648/151360180-72b31616-f40f-4bce-a78c-25124cffc675.png)


Assim como para o NetCTL, o resultado da predição do MHC-II Binding Predictions é uma página .html (escolha a opção “Text file” na seleção “Output format”) contendo os peptídeos, os scores definidos pelos algoritmos usados no servidor e algumas outras informações. Esta página html deverá ser baixada através da opção “salvar como” (click do botão direito do mouse sobre a página). O usuário deverá nomear o arquivo .html resultante da seguinte forma:

> < proteína >_< organismo >.html
 
 ou
 
> polyp_< organismo >.html (caso esteja trabalhando com poliproteínas completas)


ATENÇÃO: Esperar a página de resultados carregar por completo para fazer o download.
Como a predição ocorre com uma proteína por vez, todos os arquivos .html deverão ser organizados em uma única pasta (ou diretório). O POA1 utilizará o caminho (path) desta pasta para encontrar os arquivos.


###### 3.1.1.5 Outros preditores


Esta seção é aberta aos resultados de predição de qualquer outro preditor de peptídeos antigênicos, além dos que foram mencionados acima. Os resultados destas outras predições (não se restringe a um único método) devem ser organizados em um arquivo formato fasta. Antes de ser submetido ao POA1, o arquivo fasta contendo os epítopos preditos deverá ser organizado no seguinte formato:

```
><proteína>_<organismo>_<método_preditor>_<ID_NCBI(se_houver)>_<posição_inicial>_<posição_final>

Epitopo1
 
><proteína>_<organismo>_<método_preditor>_<ID_NCBI(se_houver)>_<posição_inicial>_<posição_final>
 
Epitopo2
…
```

#### 3.1.2 Argumentos obrigatórios

 ![3d](https://user-images.githubusercontent.com/72517648/151359262-9553ca42-63ad-4c9e-976f-0ac6704f49df.png)


 __parâmetro -b:__ Arquivo .Json contendo o resultado da predição do Bepipred 2.0 (seção 3.1.1.1);
 
 __parâmetro -p:__ Arquivo .fasta contendo os resultados da predição do PAP/IMED (seção 3.1.1.2);
 
 __parâmetro -m:__ Diretório/pasta contendo todos os arquivos .html resultados da análise do MHC-II Binding Predictions (seção 3.1.1.4);
 
 __parâmetro -n:__ Arquivo .html contendo o resultado da análise do NetCTL 1.2 (seção 3.1.1.3);
 
 __parâmetro -x:__ Arquivo .fasta contendo o resultado da análise de algum outro preditor escolhido pelo usuário, seguindo formatação pré-definida (ver seção 3.1.1.5);

 
 ATENÇÃO: Deverão ser submetidos os resultados de pelo menos um método de predição de epítopos (ver seção 3.1.1).


 Além dos resultados da predição de epitopos (__argumentos n, m, b, p, e x__), para utilizar o POA1 é necessário submeter ao pipeline  (__parâmetro -f__) um arquivo fasta contendo todas as poliproteínas (ou proteínas, caso tenham sido utilizadas em sequências separadas) usadas para a predição de epítopos. Este arquivo deve seguir a seguinte formatação:


 Para poliproteínas completas:


 ```
 >polyp_<nome_do_organismo1>_<identificador1_NCBI(se houver)>...

 Sequencia de aminoácidos1 ...

 >polyp_<nome_do_organismo2>_<identificador2_NCBI(se houver)>...

 Sequencia de aminoácidos2 ...
 ...

 ```
 
 Para proteínas completas:
 
 ```
 ><nome_da_proteína1>_<nome_do_organismo1>_<identificador1_NCBI(se houver)>...

 Sequencia de aminoácidos1 ...

 ><nome_da_proteína2>_<nome_do_organismo2>_<identificador2_NCBI(se houver)>...

 Sequencia de aminoácidos2 ...
 ...

 ```
 

 Para as análises realizadas, a informação do cabeçalho será identificada a partir do símbolo separador (“_”).


 Além disso, deverá ser designada obrigatoriamente uma pasta para receber os resultados da análise (__parâmetro -d__).



#### 3.1.3 Argumentos opcionais

  ![3d](https://user-images.githubusercontent.com/72517648/151359262-9553ca42-63ad-4c9e-976f-0ac6704f49df.png)


 __parâmetro -mhla:__ Tipo de alelo de HLA que estará sendo investigado na triagem do resultado da análise de ligantes do MHC-II, alelos humanos HLA-DP, HLA-DQ e HLA-DR (default = DR).
 
 __parâmetro -mic:__ Threshold do IC-50 do algoritmo NN_align 2.3 usado como referência do POA1 para seleção dos epitopos. Valores de IC50 < 50 nM admitem alta afinidade de ligação do epitopo com o MHCII, IC50 < 500 nM é considerado uma ligação de afinidade intermediária e IC50 < 5000 nM caracterizam ligações de baixa afinidade (default = 50).
 
 __parâmetro -bmin:__ Valor mínimo para tamanho de epitopos preditos pela ferramenta Bepipred (default = 0, ou seja não há seleção de tamanho mínimo para os epitopos Bepipred).
 
 __parâmetro -bmax:__ Valor máximo para tamanho dos epitopos preditos pela ferramenta Bepipred (default = 0, ou seja não há seleção de tamanho máximo para os epitopos Bepipred).
 
 __parâmetro -pmin:__ Valor mínimo para tamanho de epitopos preditos pela ferramenta PAP/IMED (default = 0, ou seja não há seleção de tamanho mínimo para os epitopos PAP/IMED).
 
 __parâmetro -pmax:__ Valor máximo para tamanho dos epitopos preditos pela ferramenta PAP/IMED  (default = 0, ou seja não há seleção de tamanho máximo para os epitopos PAP/IMED).
 
 __parâmetro -xmin:__ Valor mínimo para tamanho de epitopos preditos por outros preditores que não aqueles mencionados na seção 3.1.1  (default = 0, ou seja não há seleção de tamanho para estes epítopos).
 
 __parâmetro -xmax:__ Valor máximo para tamanho dos epitopos preditos por outros preditores que não aqueles mencionados na seção 3.1.1 (default = 0, ou seja não há seleção de tamanho para estes epítopos).
 
 __parâmetro -e:__ Opção para que os resultados também sejam organizados em planilha Excel (.xslx), além dos arquivos fastas que são gerados obrigatoriamente (default = n).



#### 3.1.4 Rodando o POA1

Com os arquivos dos resultados das predições em mãos, pode-se executar o POA1.

![2](https://user-images.githubusercontent.com/72517648/151359259-0473b328-cb54-423b-ac3c-f49cefc94213.png)

#### 3.1.5 Arquivos de Saída 

Como resultado da análise do POA1, os epítopos serão organizados em arquivos fasta para serem submetidos a ferramenta web Epitope Analysis Conservancy. O fasta de epítopos serão separados pela espécie do organismo, ou seja, ao final da análise todos os epítopos da espécie X serão organizados em um único fasta. O fasta contendo estes epitopos será nomeado:


X_epitope.fasta (<sp_do_organismo>_epitope.fasta)

Na análise de conservação, estes epítopos deverão ser comparados às sequências de proteínas de outras espécies diferentes da qual eles foram retirados. Desta forma, para evitar redundâncias na comparação entre epitopos-proteínas, sugere-se criar um fasta contendo todas as proteínas presentes na análise, exceto aquelas do organismo em que os próprios epítopos foram retirados. Sugere-se ainda utilizar poliproteínas completas para fazer a comparação dos epitopos.


O POA1 constrói ainda um relatório de análises contendo algumas estatísticas dos epitopos analisados. Este relatório contém algumas informações como número de proteínas, número de epítopos e quantidade maior e menor de epitopos dentre as proteínas. Além disso, o pipeline permite a exportação dos epitopos em formato .xslx (Excel) (ver seção 3.1.3 - argumentos opcionais).


### 3.2 POA2


Após a análise de conservação feita pela ferramenta Epitope Conservancy Analysis (IEDB) com os arquivos resultantes do POA1, os dados de conservação para cada peptídeo já podem ser trabalhados. O POA2 recebe os resultados da predição de análise de conservação e seleciona os peptídeos com base em valores de threshold estabelecidos. Além disso, utiliza a ferramenta pyTMHMM  para classificar o resultado quanto à topologia de membrana, a partir da predição de hélices transmembranares. Os peptídeos caracterizados são armazenados em arquivos de planilha (.xlsx) e em arquivos fasta.


#### 3.2.1 Arquivos de entrada


Como discutido, os resultados da primeira metade do pipeline (POA1) devem ser submetidos a análise de conservação na ferramenta online [Epitope Conservancy Analysis (IEDB)](http://tools.iedb.org/conservancy/). A análise de conservação deverá ser feita para todos os fastas de epítopos (separados por espécie do organismo no POA1). Os resultados são tabelas (formato .csv) que deverão ser baixadas e armazenadas em um único diretório. O POA2 recebe e avalia cada um dos arquivos fasta presentes no diretório, seleciona os epitopos com base no valor de conservação definido pelo usuário e organiza em uma planilha os epítopos que são mais ou menos conservados (ver seção 3.2.2) e sua classificação quanto a topologia da membrana. 


#### 3.2.2 Argumentos Obrigatórios

 ![Design sem nome](https://user-images.githubusercontent.com/72517648/151359269-a50cbc97-443b-4dc4-9082-b1b52527a566.png)


Para utilizar o POA2 será necessário organizar todos os resultados da análise de conservação em um único diretório, o caminho (path) deste diretório será utilizado pelo POA2 (__parâmetro -d__). 

O usuário deverá informar ao POA2 qual foi o objetivo da análise, se esta foi feita visando encontrar epítopos conservados (>=) ou epítopos únicos (<), __parâmetro -g ou -l__, respectivamente, e qual o threshold de identidade de sequência (Sequence identity threshold) foi utilizado para fazer a predição (__parâmetro -t__).

Por último, deverá apresentar um arquivo fasta que contenha todas as polyproteínas/proteínas usadas para a predição de epítopos (__parâmetro -f__). Este arquivo poderá ser igual àquele submetido no POA1, seguindo a mesma formatação (ver seção 3.1.2).



#### 3.2.3 Argumentos Opcionais

![Design sem nome](https://user-images.githubusercontent.com/72517648/151359269-a50cbc97-443b-4dc4-9082-b1b52527a566.png)


 __parâmetro -r:__ Diretório onde os resultados da análise do POA2 (planilhas e ou arquivos fasta) serão organizados;
 
 __parâmetro -imin:__ Threshold mínimo de identidade aceito no resultado da análise de conservação entre os peptídeos e as proteínas avaliadas (default = 60).
 
 __parâmetro -imax:__ Threshold máximo de identidade aceito no resultado da análise de conservação entre os peptídeos e as proteínas avaliadas (default = 100).
 
 __parâmetro -m:__ porcentagem das sequências que tiveram match de identidade acima do threshold de identidade de sequência avaliado (parâmtero -t) (default = 60).
 
 __parâmetro -rf:__ Opções de arquivo fasta contendo os resultados da análise de conservação (default = None). 
 [0] Todos os epitopos resultantes da análise de conservação, sem considerar a classificação por topologia de membrana.
 [1] Os epitopos resultantes da análise de conservação que estiverem em porções expostas da membrana.
 [2] Os epitopos resultantes da análise de conservação que estiverem em porções transmembranares.
 [3] Os epitopos resultantes da análise de conservação que estiverem em porções internas à membrana.


#### 3.2.4 Rodando o POA2


Com os arquivos dos resultados das predições em mãos, pode-se executar o pipeline.
Ex: Análise de conservação para epítopos conservados (-g True)

![dfd](https://user-images.githubusercontent.com/72517648/151359274-2b88b315-0f56-4494-84a8-7025be6e21df.png)


#### 3.2.5 Arquivos de Saída


Como resultado da análise do POA2, todos os epítopos selecionados na análise de conservação são organizados em um único arquivo excel (.xlsx). Além de conter as informações já disponibilizadas na planilha de resultado da ferramenta do IEDB, a tabela conta com colunas que identificam a porções do epitopo quanto a topologia da membrana, em porções externas, transmembranares e internas (representação da porcentagem, valores de 0 a 1). Esta análise de topologia da membrana é feita usando um algoritmo derivado do TMHMM ([Sonnhammer et al., 1998](https://pubmed.ncbi.nlm.nih.gov/9783223/)), o [pyTMHMM](https://github.com/bosborne/pyTMHMM). 
Os epitopos triados podem também ser organizados em arquivos fasta, como demonstrado na seção 3.2.3, em __parâmetro -rf__.

## CITAÇÕES

.....

--------------------------------------------------------------

## English

This is the English version of the README. The Portuguese (PT-BR) version is at the top of this
document.

## 0. AUTOMATION & GRAPHICAL INTERFACE (v1.1)

> This release adds an optional automation layer and a graphical interface **without changing the
> scientific ranking/selection logic**. The original command line still works exactly as before.
> The full design and per-tool automation strategy are documented in [`ARCHITECTURE.md`](ARCHITECTURE.md).
> A complete report of everything that changed (PT-BR) is in [`MUDANCAS.md`](MUDANCAS.md) —
> also available as [`MUDANCAS.pdf`](MUDANCAS.pdf).

**Package layout** (refactored into layers; the flat scripts became thin shims):
`poa/core` (parsing + pipeline logic, no network) · `poa/services` (external tools) ·
`poa/cli` (POA1/POA2, same arguments) · `poa/gui` (Streamlit app).

**Install**
```
python -m venv .venv
# Windows: .venv\Scripts\activate   |   Linux/macOS: source .venv/bin/activate
pip install -r requirements.txt
```
Extras that are not pip packages / are heavy (install only if you want that automation):
* **EMBOSS `antigenic`** (replaces the offline PAP/IMED site, same Kolaskar–Tongaonkar method):
  `conda install -c bioconda emboss` (on Windows: WSL2 + conda, or Docker).
* **BepiPred-3.0** local model (PyTorch + ESM-2): `pip install bp3`.
* **Browser automation** for NetCTL 1.2 / BepiPred-2.0: `pip install playwright && playwright install chromium`.

**Run the GUI** (guided 6-step flow: proteins → predictions → POA1 → conservancy → POA2 → results):
```
streamlit run poa/gui/app.py      # or:  python run_gui.py
```

**Run the CLI** (unchanged):
```
python POA1_v1.0.py -f proteins.fasta -d results -b3 bepipred3.fasta -e y
python POA2_v1.0.py -g True -t 70 -d conservancy_csvs -f proteins.fasta -rf 1
```

**Automation status per tool** (verified 2026 — details in `ARCHITECTURE.md`):

| Tool | Strategy | Fallback |
|---|---|---|
| MHC-II Binding (IEDB) | REST API (synchronous) | manual upload |
| BepiPred-3.0 | local `bp3` package | manual upload |
| PAP/IMED | local EMBOSS `antigenic` | manual upload |
| Epitope Conservancy (IEDB) | **local reimplementation** (no API/standalone exists) | manual CSV upload |
| NetCTL 1.2 / BepiPred-2.0 | best-effort browser automation | manual upload |
| pyTMHMM | local (unchanged) | — |

In every case, if a service is unavailable the pipeline logs a clear message and lets you upload
the result file manually (preserving the original semi-automatic behavior).

**Tests:** `pytest` (parsers, ranking/filtering, local conservancy, service clients, GUI backend).

> **Fixed:** the NetCTL 1.2 HTML parser had a pre-existing bug (row append via a non-sequential
> index, incompatible with pandas ≥ 2; and marker/append logic inconsistent with the `<-E`
> epitope selection). It was corrected (user-approved) to use a sequential index and to pad only
> short rows to the column count, keeping the `<-E` selection documented in the README; covered
> by unit tests.

## 1. INTRODUCTION

High-throughput screening, genomic sequence analysis of pathogenic organisms, and reverse vaccinology techniques have revolutionized the rapid development of vaccines and diagnostic tests. In this context, in silico approaches, such as multi-epitope construction, have become indispensable for designing biotechnological products and immunotherapeutic interventions against viral, bacterial, and parasitic diseases ([Dar et al., 2021](https://www.nature.com/articles/s41598-021-90868-2); [Enayatkhani et al., 2021](https://www.tandfonline.com/doi/full/10.1080/07391102.2020.1756411?casa_token=Ctu83G9J4aEAAAAA%3Aocl_27U5qwcFn-WXSFZIaaSmChfq-hC5ggb281Z1hGuocFuOXH65acPngTeO6HavrexguLxkdEbfAVQ); [Shey et al., 2019](https://www.nature.com/articles/s41598-019-40833-x)).

The Antigen Optimization Pipeline (POA, in Portuguese version) is a semi-automated Python-based pipeline designed to analyze, organize, retrieve, and screen epitopes derived from antigenic peptide predictions using various web-based tools. POA integrates results from B cell, cytotoxic T cell (Tc), and T helper (Th) cell epitope prediction tools to identify optimal antigenic targets and construct a comprehensive antigen database. Additionally, POA enables the identification of conserved epitopes across species using the [Epitope Conservancy Analysis - IEDB](http://tools.iedb.org/conservancy) tool and characterizes epitope localization in transmembrane regions through the [pyTMHMM](https://github.com/bosborne/pyTMHMM) tool. By integrating these approaches, POA identifies high-potential candidates within viral proteins, driving advancements in integrative therapeutic interventions. This capability provides innovative strategies to tackle the global burden of zoonotic and vector-borne diseases, paving the way for the development of targeted vaccines and immunotherapies.

### 1.1 Antigenic peptide prediction

POA optimizes the organization and visualization of in silico epitope prediction results, enabling more efficient analysis of large datasets. It also optimizes the screening process by identifying unique or shared epitopes within a group of organisms and characterizing epitopes exposed on membrane structures (Batista et al., 2025).

Supported Prediction Tools:

* B cell epitopes: [Bepipred 2.0](https://services.healthtech.dtu.dk/service.php?BepiPred-2.0) and [Predicted Antigenic Peptides (PAP) - IMED](http://imed.med.ucm.es/Tools/antigenic.pl).

* Tc epitopes: [NetCTL - 1.2](https://services.healthtech.dtu.dk/service.php?NetCTL-1.2).

* Th epitopes: [MHC-II Binding Predictions](http://tools.iedb.org/mhcii/).

POA is also capable of analyzing results from additional prediction tools, provided the data is formatted into standardized FASTA files (see Section 3.1.1 POA1 - Input files).

## 2. INSTALLATION

POA v1.0 must be installed manually. Ensure all dependencies are installed for proper functionality. The required dependencies are:

* Python version > 3.6.0
* Pandas version >= 1.4.0
* Numpy version >= 1.21.4
* Argparse version 1.4.0
* Datetime version 4.4
* Biopython version 1.79
* Html version 1.16
* Json
* pyTMHMM (This package needs manual installation. Please follow the instructions in [pyTMHMM tool](https://github.com/bosborne/pyTMHMM)) 

To install POA v1.0 on Linux, use the following command:

```
git clone https://github.com/UbiratanBatista/POA_Project.git
```

## 3. USING THE POA PIPELINE

POA is a semi-automated pipeline, requiring some manual steps by the user. The pipeline consists of two main stages: POA1 and POA2.

The first stage of the pipeline, POA1, analyzes prediction results from epitope prediction tools (see Section 1.1 - Antigenic Peptide Prediction), creates an antigen database based on epitope ranking values, and organizes the results into multiple FASTA files for use in the next stage.

Following this, the output files from POA1 should be submitted to the [Epitope Conservancy Analysis - IEDB](http://tools.iedb.org/conservancy) web platform. The conservancy analysis results are then processed by POA2, the second stage of the pipeline. POA2 evaluates epitopes with higher or lower conservancy (depending on the analysis objective) and performs transmembrane structure prediction using TMHMM. Ultimately, POA2 returns the conservancy analysis results along with the characterization of epitope localization in the membrane (external, transmembrane, or internal).

### 3.1 POA1

POA1 analyses prediction results from web tools, which may be in various formats, ranks epitopes based on tool-specific criteria, and generates a curated antigen database. The results are stored in standardized matrices and subsequently organized into multiple FASTA files for downstream analysis.

#### 3.1.1 POA1 - Input files

The first step is to perform predictions in the epitope prediction web tools (see section 1.1), following the protocols and parameters defined for each method and objective. POA1 accepts output files from various epitope prediction tools. Below are the specific requirements for each tool:

##### 3.1.1.1 Bepipred 2.0

Submit a ```.fasta``` file containing all proteins for analysis.

Download the results as a ```.json``` file (JSON Summary).

Example output:
![Captura de tela de 2022-01-26 17-45-11](https://user-images.githubusercontent.com/72517648/151360191-20879962-26e7-4df5-a79e-1f9d79082c8d.png)

*The result of the Bepipred 2.0 analysis is a web page (html) containing the proteins and the predicted antigenic regions (marked with E) in each sequence.

##### 3.1.1.2 Predicting Antigenic Peptides (PAP) - IMED

For the PAP tool, ensure each protein sequence is submitted separately, with the header formatted as " ```><protein>_<organism>...``` ". The prediction result will appear as follows:
![Captura de tela de 2022-01-26 17-56-13](https://user-images.githubusercontent.com/72517648/151360188-21f130dd-ae97-4067-b7c9-8f138bff95bf.png)

Copy the prediction table and paste it into a ```.txt``` file, adhering to the following format:
![Captura de tela de 2022-01-26 17-58-37](https://user-images.githubusercontent.com/72517648/151360185-4a3b4863-5eab-4024-99d9-5105ac705e53.png)

##### 3.1.1.3 NetCTL 1.2

Submit a protein ```.fasta`` file for the NetCTL 1.2 analysis. The prediction result will appear as follows:
![Captura de tela de 2022-01-26 18-02-17](https://user-images.githubusercontent.com/72517648/151360184-6091362a-bcb6-455f-b363-64e575096806.png)

*The NetCTL 1.2 output is a web page (html) containing many peptides and those predicted to be antigenic (indicated with < - E). Save the results as an ```.html``` file (ATTENTION: ensure the page is fully loaded before saving).

##### 3.1.1.4 MHC-II Binding Predictions

Submit proteins individually to the MHCIIBP tool. The output will resemble the following:
![Captura de tela de 2022-01-26 18-07-50](https://user-images.githubusercontent.com/72517648/151360180-72b31616-f40f-4bce-a78c-25124cffc675.png)

The MHCIIBP result is an web page (select “Text file” under “Output format”) containing the predicted peptides, their scores as defined by the prediction algorithms, and additional relevant information. Save the results as ```.html``` files, ensuring the page is fully loaded before saving. Name the files as ```<protein>_<organism>.html``` or, for complete polyproteins, ```polyp_<organism>.html```. Organize all ```.html``` files into a **single folder (or directory)** for further processing.

##### 3.1.1.5 Other Predictors

Results from additional prediction tools can also be analyzed by POA, provided they are formatted into a standardized ```FASTA``` file. The file should be organized as follows:
```
><protein>_<organism>_<predictor_method>_<NCBI_ID_if_any>_<start_position>_<end_position>
Epitope 1

><protein>_<organism>_<predictor_method>_<NCBI_ID_if_any>_<start_position>_<end_position>
Epitope 2
```

#### 3.1.2 POA1 - Mandatory Arguments

![3d](https://user-images.githubusercontent.com/72517648/151359262-9553ca42-63ad-4c9e-976f-0ac6704f49df.png)

* **parameter -b2**: ```.json``` file containing the Bepipred 2.0 prediction results (see Section 3.1.1.1 Bepipred 2.0).
* **parameter -b3**: ```.fasta``` file containing the Bepipred 3.0 prediction results.
* **parameter -p**: ```.fasta``` file containing the PAP-IMED prediction results (see Section 3.1.1.2 Predicting Antigenic Peptides (PAP) - IMED).
* **parameter -n**: ```.html``` file containing the NetCTL 1.2 analysis results (see Section 3.1.1.3 NetCTL 1.2).
* **parameter -m**: Directory containing all ```.html``` files from the MHC-II Binding Predictions analysis (see Section 3.1.1.4 MHC-II Binding Predictions).
* **parameter -x**: ```.fasta``` file containing results from any other predictor, formatted as specified (see Section 3.1.1.5 Other Predictors).

**Note:** At least one epitope prediction result must be submitted (see Section 3.1.1 POA1 - Input files).

In addition to the prediction results (parameters ```-b```, ```-p```, ```-m```,```-n```, and ```-x```), a ```.fasta``` file containing all proteins (or polyproteins) used for epitope prediction must be submitted (**parameter -f**). The file should follow the format below:

For complete polyproteins:
```
>polyp_<organism_name1>_<NCBI_identifier1_if_any>...

Amino_acid_sequence 1

>polyp_<organism_name2>_<NCBI_identifier2_if_any>...

Amino_acid_sequence 2
```

For individual proteins:
```
><protein_name1>_<organism_name1>_<NCBI_identifier1_if_any>...

Amino_acid_sequence 1

><protein_name2>_<organism_name2>_<NCBI_identifier2_if_any>...

Amino_acid_sequence 2
```
The header information is parsed using the separator symbol ("_") for analysis.

Additionally, a folder MUST be designated to store the analysis results (**parameter -d**).

#### 3.1.3 POA1 - Optional Arguments

![3d](https://user-images.githubusercontent.com/72517648/151359262-9553ca42-63ad-4c9e-976f-0ac6704f49df.png)

* **parameter -bmin**: Minimum size for epitopes predicted by Bepipred (default = 0, no minimum size).
* **parameter -bmax**: Maximum size for epitopes predicted by Bepipred (default = 0, no maximum size).
* **parameter -pmin**: Minimum size for epitopes predicted by PAP - IMED (default = 0, no minimum size).
* **parameter -pmax**: Maximum size for epitopes predicted by PAP - IMED (default = 0, no maximum size).
* **parameter -mhla**: Type of HLA allele investigated in the MHC-II ligand analysis (options: HLA-DP, HLA-DQ, HLA-DR; default = DR).
* **parameter -mic**: IC-50 threshold for the NN_align 2.3 algorithm used for epitope selection: IC50 < 50 nM (High-affinity binding to MHCII), IC50 < 500 nM (Intermediate affinity binding to MHCII), and IC50 < 5000 nM (Low-affinity binding to MHCII); Default = 50.
* **parameter -xmin**: Minimum size for epitopes predicted by other predictors (default = 0, no minimum size).
* **parameter -xmax**: Maximum size for epitopes predicted by other predictors (default = 0, no maximum size).
* **parameter -e**: Option to export results as an Excel spreadsheet ```(.xlsx)``` in addition to the default ```.fasta``` files (default = n).

#### 3.1.4 Running POA1

Once the prediction result files are prepared, POA1 can be executed.

![2](https://user-images.githubusercontent.com/72517648/151359259-0473b328-cb54-423b-ac3c-f49cefc94213.png)

#### 3.1.5 POA1 - Output files

POA1 organizes the epitopes into ```.fasta``` files for submission to the Epitope Conservancy Analysis web tool. The epitopes are grouped by organism species, with all epitopes from a species stored in a single ```.fasta``` file named:

```<SpeciesID_epitope>.fasta```

For the conservancy analysis, these epitopes should be compared to protein sequences from other species. To avoid redundancy, it is recommended to create a ```.fasta``` file containing all proteins except those from the organism from which the epitopes were extracted. Using complete polyproteins for comparison is also suggested.

POA1 also generates an analysis report containing statistics on the analyzed epitopes, including: the total number of proteins, the total number of epitopes, and the proteins with the highest and lowest number of epitopes. Additionally, the pipeline provides the option to export the epitopes in ```.xlsx``` format (see Section 3.1.3 POA1 - Optional Arguments).

### 3.2 POA2

After the [Epitope Conservancy Analysis - IEDB](http://tools.iedb.org/conservancy) tool processes the files generated by POA1, POA2 analyzes the conservancy data for each peptide. POA2 selects peptides based on user-defined conservancy thresholds and uses the [pyTMHMM](https://github.com/bosborne/pyTMHMM) tool to classify their membrane topology based on transmembrane helix predictions. The results are stored in ```.xlsx``` (spreadsheet) and ```.fasta``` files for further use.

#### 3.2.1 POA2 - Input files

The results from POA1 must first undergo conservancy analysis using the [Epitope Conservancy Analysis - IEDB](http://tools.iedb.org/conservancy) tool. This analysis should be performed for all ```.fasta``` files (organized by organism species in POA1). The results, in ```.csv``` format, should be downloaded and stored in a single directory. POA2 processes these files, selects epitopes based on the user-defined conservancy threshold, and organizes them into a spreadsheet. The spreadsheet includes epitopes with higher or lower conservancy (see Section 3.2.2 POA2 - Mandatory Arguments) and their membrane topology classification.

#### 3.2.2 POA2 - Mandatory Arguments

![Design sem nome](https://user-images.githubusercontent.com/72517648/151359269-a50cbc97-443b-4dc4-9082-b1b52527a566.png)

* **parameter -d**: Directory containing the conservancy analysis results (.csv files).
* **parameter -g or -l:**: Specifies the analysis objective: Use ```-g``` to identify conserved epitopes (>= threshold), or use ```-l``` to identify unique epitopes (< threshold).
* **parameter -t**: Sequence identity threshold used for the conservancy analysis.
* **parameter -f**: ```.fasta``` file containing all polyproteins/proteins used for epitope prediction (same file as submitted in POA1; see Section 3.1.2 POA1 - Mandatory Arguments for formatting details).

#### 3.2.3 POA2 - Optional Arguments

![Design sem nome](https://user-images.githubusercontent.com/72517648/151359269-a50cbc97-443b-4dc4-9082-b1b52527a566.png)

* **parameter -r**: Directory to store POA2 analysis results (```.xlsx``` and/or ```.fasta``` files).
* **parameter -imin**: Minimum identity threshold for conservancy analysis (default = 60).
* **parameter -imax**: Maximum identity threshold for conservancy analysis (default = 100).
* **parameter -m**: Percentage of sequences with identity matches above the threshold (default = 60).
* **parameter -rf**: Options for organizing ```.fasta``` files based on membrane topology: [0]: All epitopes from the conservancy analysis (no membrane topology classification); [1]: Epitopes located in exposed membrane portions; [2]: Epitopes located in transmembrane portions; [3]: Epitopes located in internal membrane portions. Default = 0

#### 3.2.4 Running POA2

Once the conservancy analysis files are prepared, POA2 can be executed. Below is an example command for analyzing conserved epitopes (-g True):

![dfd](https://user-images.githubusercontent.com/72517648/151359274-2b88b315-0f56-4494-84a8-7025be6e21df.png)

#### 3.2.5 POA2 - Output files

POA2 generates an Excel file ```.(.xlsx)```. containing all selected epitopes from the conservancy analysis. The spreadsheet includes:

* Conservancy data from the Epitope Conservancy Analysis (IEDB) tool.

* Membrane topology classification (external, transmembrane, or internal) as a percentage (values from 0 to 1), derived from the [pyTMHMM](https://github.com/bosborne/pyTMHMM) algorithm ([Sonnhammer et al., 1998](https://pubmed.ncbi.nlm.nih.gov/9783223/)).

Additionally, the screened epitopes can be organized into ```.fasta``` files based on the options specified in the -rf parameter (see Section 3.2.3 POA2 - Optional Arguments).

## 4. CITATIONS

......................

*END*
