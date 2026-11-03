"""Tooltip text for the GUI fields, in one reviewable place.

Streamlit's ``help=`` renders a small "?" next to a field's label and shows the text on hover.
Keeping the strings here instead of inline does two things: the forms in :mod:`poa.gui.app` stay
readable, and the wording — which is the pipeline's user-facing documentation — can be reviewed
and corrected without touching the layout.

Each entry answers, in this order: **what the field is**, **what to put in it**, and, when the
answer is not obvious, **what goes wrong if it is wrong**. Only fields whose meaning is not
evident from the label are covered; a trivial field with a tooltip is noise that trains people to
ignore the "?" on the fields that need one.

The three conservancy numbers (``t``, ``m``, ``imin``/``imax``) get the most care: they are
different quantities with similar names, and confusing them is the most common way to read a POA2
result as meaning something it does not.
"""
from __future__ import annotations

HELP = {
    # ----------------------------------------------------------------- step 1: POA0 search
    "poa0_descendants": (
        "Desmarcado, a busca usa `organism_id`, que casa **só** com o táxon exato. Marcado, usa "
        "`taxonomy_id`, que inclui cepas e isolados descendentes. Para vírus a diferença costuma "
        "ser de uma ordem de grandeza — vale marcar quando o táxon exato devolve pouca coisa, e "
        "deixar desmarcado quando você quer a espécie de referência e não os isolados."
    ),
    "poa0_token_name": (
        "O campo de proteína do cabeçalho `Proteína_Espécie_ID` (ex.: `GcHead`). **Não pode "
        "conter `_`**: o pipeline separa os campos por underscore, então um underscore aqui "
        "desloca todos os campos seguintes e a espécie passa a ser lida do meio do nome da "
        "proteína. Pontuação e espaços são removidos automaticamente."
    ),
    "poa0_token_specie": (
        "O campo de espécie do cabeçalho (ex.: `OROV`, `DENV1`). É por ele que o POA2 liga cada "
        "epítopo à sua sequência de referência, então use a mesma sigla em todas as proteínas da "
        "mesma espécie. Também não pode conter `_`."
    ),
    "poa0_structure": (
        "Baixa a estrutura escolhida do RCSB e grava junto o mapeamento SIFTS — a correspondência "
        "exata entre os resíduos da sequência e a numeração que o arquivo usa. É isso que permite "
        "à etapa 7 destacar os epítopos na posição certa em vez de supor que a estrutura começa "
        "no resíduo 1. A lista vem ordenada por cobertura da sequência e, em seguida, resolução."
    ),

    # ----------------------------------------------------------------- step 1: proteins & params
    "proteins_fasta": (
        "O conjunto de proteínas sobre o qual os epítopos foram preditos — o `-f` do POA1. "
        "Cabeçalho na convenção `Proteína_Espécie_ID` (ex.: `NS1_DENV1_ref`); o separador é o "
        "`_`. É daqui que o pipeline tira a espécie de cada epítopo, então um cabeçalho fora da "
        "convenção faz os epítopos saírem sem espécie. Pode enviar um arquivo por espécie: eles "
        "são unidos num único `-f`."
    ),
    "swap_headers": (
        "Marque se os seus cabeçalhos estão na ordem inversa, `Espécie_Proteína` (ex.: "
        "`alphavirus_e2`). Sem isso, `alphavirus_e2` é lido como proteína `ALPHAVIRUS` da espécie "
        "`E2`, e todas as sequências caem numa 'espécie' que é na verdade o nome da proteína — sem "
        "mensagem de erro, só um resultado errado.\n\n"
        "Esta caixa define a convenção **da análise inteira**: é o `-f` que sai daqui que diz, "
        "daí em diante, qual campo é a proteína e qual é a espécie.\n\n"
        "Os arquivos de predição não precisam seguir a mesma ordem. O POA1 compara cada epítopo "
        "com o `-f` e corrige os que vierem invertidos — o que importa porque eles não são "
        "uniformes: um arquivo que **você** envia traz a ordem da ferramenta que o gerou, "
        "enquanto um arquivo que a interface **gera** (EMBOSS antigenic, BepiPred-3.0, API do "
        "MHC-II) já sai na ordem do `-f`, porque é dele que foi gerado."
    ),
    "len_min": (
        "Comprimento mínimo, em resíduos, para um epítopo deste método ser mantido. "
        "**0 desliga o filtro.** Útil para descartar fragmentos curtos demais para serem alvos."
    ),
    "len_max": (
        "Comprimento máximo, em resíduos, para um epítopo deste método ser mantido. "
        "**0 desliga o filtro.**"
    ),
    "mhla": (
        "Classe do alelo HLA considerada na triagem do resultado do MHC-II: DR, DP ou DQ. "
        "Só as linhas do alelo desta classe entram na análise. DR é o padrão e o mais estudado."
    ),
    "mic": (
        "Limiar de IC50 (nM) do NN_align para um peptídeo ser aceito como ligante do MHC-II: "
        "mantém as linhas com IC50 **menor ou igual** a este valor. Referência: < 50 nM é alta "
        "afinidade, < 500 nM intermediária, < 5000 nM baixa. Quanto maior o valor, mais "
        "permissiva a seleção."
    ),

    # ----------------------------------------------------------------- step 2: predictions
    "mhcii_alleles": (
        "Alelos a submeter à API do IEDB, separados por vírgula, na notação do próprio IEDB — "
        "ex.: `HLA-DRB1*01:01, HLA-DRB1*03:01`. Cada alelo é uma predição a mais, então a "
        "consulta fica mais lenta à medida que você adiciona."
    ),
    "mhcii_length": (
        "Comprimento, em resíduos, dos peptídeos que o MHC-II vai avaliar. 15 é o usual para "
        "células T-auxiliares, porque é a faixa que a fenda do MHC classe II acomoda."
    ),
    "mhcii_upload": (
        "Plano B para quando a API do IEDB não responder: os `.html` salvos do site, **um por "
        "proteína**, nomeados `Proteína_Espécie.html` (ou `polyp_Espécie.html` para "
        "poliproteínas). O nome do arquivo é de onde saem a proteína e a espécie — um nome fora "
        "do padrão não é lido."
    ),
    "b2_json": (
        "O *JSON Summary* baixado do BepiPred-2.0. Se as chaves de antígeno já seguirem "
        "`Proteína_Espécie_ID`, o arquivo é usado direto; se vierem genéricas (`Sequence`, "
        "`denv1`), a interface abre um mapeamento para você dizer a espécie e a proteína de cada "
        "antígeno, e reconstrói as proteínas de referência a partir do próprio JSON."
    ),
    # one per method, since each tool's file format is the thing people get wrong
    "upload_b3": (
        "O FASTA de saída do BepiPred-3.0, em que **letras maiúsculas marcam o epítopo** e as "
        "minúsculas o resto da sequência — é assim que o parser identifica as regiões. Cabeçalho "
        "em `Proteína_Espécie_ID`. Pode enviar um arquivo por espécie."
    ),
    "upload_p": (
        "A tabela de resultado do PAP/IMED colada num `.txt`: blocos `>ID` seguidos de linhas "
        "com número, posição inicial, sequência e posição final separados por tabulação. O site "
        "original saiu do ar — o caminho recomendado hoje é o EMBOSS `antigenic` acima, que usa "
        "o mesmo método (Kolaskar & Tongaonkar)."
    ),
    "upload_n": (
        "A página de resultado do NetCTL 1.2 salva como `.html`. Espere a página carregar por "
        "completo antes de salvar. São selecionados os peptídeos marcados com `<-E`, que é como "
        "o NetCTL indica o ligante identificado."
    ),
    "upload_x": (
        "Resultados de qualquer outro preditor, organizados num FASTA com cabeçalho "
        "`Proteína_Espécie_Método_ID_Início_Fim`. O ID do NCBI é opcional — sem ele o cabeçalho "
        "tem cinco campos, `Proteína_Espécie_Método_Início_Fim`. Uma sequência pode ocupar "
        "várias linhas."
    ),
    "map_specie": (
        "A espécie deste antígeno, como ela deve aparecer nos resultados — ex.: `DENV1`. É a "
        "chave que liga o epítopo à proteína dele no `-f` e que separa os arquivos da etapa 4, "
        "então use o mesmo nome nos dois lugares."
    ),
    "map_protein": (
        "A proteína deste antígeno — ex.: `NS1`, `E`, ou `polyp` para uma poliproteína inteira."
    ),

    # ----------------------------------------------------------------- step 4: conservancy
    "cons_threshold": (
        "A identidade de sequência (%) a partir da qual uma proteína conta como *match* de um "
        "epítopo. É com este valor que a conservância é **calculada** — ele fica gravado junto "
        "dos CSVs e a etapa 5 o lê de volta. Não é um filtro de resultado: é o critério que "
        "define o que a coluna de % de matches significa."
    ),
    "cons_operator": (
        "O que você está procurando. **Conservados (≥)**: conta como match a proteína cuja "
        "identidade *alcança* o limiar — epítopos presentes em todo o grupo, candidatos a vacina "
        "ampla. **Únicos (<)**: conta a proteína cuja identidade fica *abaixo* do limiar — "
        "epítopos exclusivos de uma espécie, candidatos a diagnóstico diferencial."
    ),
    "cons_mode": (
        "**Reimplementação local** calcula a conservância aqui mesmo (Bui et al., 2007), sem "
        "depender do site do IEDB — que não oferece API nem versão para baixar. **Upload manual** "
        "é para quando você já rodou no site e tem os CSVs."
    ),
    "comparison_set": (
        "Com quais proteínas os epítopos desta espécie serão comparados — o *world set* de "
        "diversidade. O padrão é usar as proteínas da própria espécie que estão no `-f`, o que "
        "costuma ser pouco: **com uma única proteína de comparação todo epítopo sai com 100% de "
        "conservância**, qualquer que seja o limiar, porque o epítopo foi predito nela mesma. "
        "Envie aqui um conjunto com variação real (várias cepas ou espécies) para o limiar "
        "discriminar algo."
    ),
    "cons_csv_upload": (
        "Os `.csv` baixados do Epitope Conservancy Analysis do IEDB, um por espécie. Todos "
        "precisam ter sido gerados com o **mesmo** limiar — misturar limiares na mesma pasta "
        "filtraria cada espécie por um critério diferente, e o POA2 recusa."
    ),
    "cons_declared_threshold": (
        "Com que limiar estes CSVs foram realmente gerados. O cabeçalho deles não traz o valor, "
        "então sem esta informação a etapa 5 não tem como validar o que está filtrando."
    ),

    # ----------------------------------------------------------------- step 5: POA2
    "poa2_objective": (
        "O mesmo critério da etapa 4, agora aplicado na seleção. Ele **precisa** coincidir com o "
        "critério com que os CSVs foram gerados: filtrar com um critério e rotular o resultado "
        "com outro é erro científico silencioso, e o POA2 bloqueia a execução se divergirem."
    ),
    "poa2_rf": (
        "Gera também um FASTA dos epítopos selecionados, separado por região de membrana — útil "
        "para levar só os expostos (**externos**) adiante, que são os alcançáveis por anticorpo. "
        "A topologia vem do pyTMHMM. `não gerar` produz apenas a planilha."
    ),
    "poa2_imin": (
        "Descarta o epítopo cuja **menor** identidade entre as proteínas comparadas fique abaixo "
        "deste valor. É um piso de semelhança: garante que não exista nenhuma proteína do "
        "conjunto em que o epítopo seja muito diferente. Não confunda com o limiar (t), que "
        "define o que conta como match; este filtra o resultado já calculado."
    ),
    "poa2_imax": (
        "Descarta o epítopo cuja **maior** identidade entre as proteínas comparadas passe deste "
        "valor. Serve para procurar epítopos que *não* sejam idênticos a nada do conjunto — "
        "baixe-o quando o objetivo for epítopos únicos. Em 100 não descarta ninguém."
    ),
    "poa2_m": (
        "Percentual mínimo das proteínas comparadas que precisam ter dado match para o epítopo "
        "ser mantido. É o filtro principal de conservância: com `-m 80`, o epítopo só passa se "
        "ao menos 80% das proteínas do conjunto satisfizeram o limiar. O funil abaixo mostra "
        "quantos epítopos sobrevivem a ele."
    ),
    "poa2_idf": (
        "Por padrão o limiar só define a coluna de % de matches, como no POA original. Marcado, o "
        "POA2 também descarta epítopos cuja identidade não satisfaça o critério: **≥** exige que "
        "a identidade *mínima* alcance o limiar (conservado em todo o conjunto); **<** exige que "
        "a *máxima* fique abaixo dele (único). Ligar isto muda quais epítopos são selecionados."
    ),
    "poa2_declared_t": (
        "O limiar com que estes CSVs foram calculados. Não foi possível descobri-lo "
        "automaticamente, então informe-o — é ele que o POA2 vai aplicar e registrar no "
        "resultado."
    ),
    "poa2_declared_op": (
        "O critério com que estes CSVs foram calculados: se as proteínas contadas como match "
        "foram as que alcançaram o limiar (≥) ou as que ficaram abaixo dele (<)."
    ),

    # ----------------------------------------------------------------- step 7: visualization
    "viz_pdb": (
        "Uma estrutura `.pdb` da proteína escolhida acima, para ver os epítopos sobre ela. "
        "**Atenção:** a numeração de resíduo do PDB precisa corresponder à posição na sequência "
        "usada nas predições — estruturas experimentais costumam ter gaps e começar em outro "
        "resíduo, e nesse caso o destaque sai deslocado."
    ),
    "viz_style": (
        "Como a cadeia principal é desenhada. `cartoon` mostra a dobra (hélices e folhas) e é a "
        "escolha usual; `stick` e `sphere` mostram os átomos, úteis em regiões pequenas."
    ),
    "viz_surface": (
        "Desenha a superfície acessível ao solvente por cima da estrutura. Ajuda a julgar se um "
        "epítopo está de fato exposto — o que importa para acesso de anticorpo — mas deixa a "
        "renderização mais lenta."
    ),
    "viz_chain": (
        "Restringe o destaque a uma cadeia do PDB (ex.: `A`). Deixe vazio para considerar todas. "
        "Use quando a estrutura tiver várias cópias da proteína e o destaque em todas poluir a "
        "visualização."
    ),
}
