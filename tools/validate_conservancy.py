#!/usr/bin/env python
"""Comparar a conservância local com uma exportação real do IEDB.

A reimplementação local em :mod:`poa.services.conservancy_client` reproduz o método de
Bui et al. (BMC Bioinformatics 2007) que o Epitope Conservancy Analysis do IEDB aplica. O IEDB
**não** expõe esse cálculo por API (verificado em 12/10/2026: a página da ferramenta responde,
mas a lista de APIs não menciona conservância), então a equivalência não pode ser verificada
automaticamente — ela depende de uma exportação que você faz à mão, uma vez, e guarda.

Fluxo:

1. No site, rode a análise com **os mesmos** epítopos, **o mesmo** conjunto de proteínas e **o
   mesmo** limiar/operador que o POA usará.
2. Baixe o CSV.
3. Rode esta ferramenta. Ela recalcula localmente e compara linha a linha.

    python tools/validate_conservancy.py -i iedb.csv -e epitopos.fasta -p proteinas.fasta \\
                                         -t 70 --operador ">="

Código de saída 0 quando tudo bate dentro da tolerância, 1 quando há divergência, 2 em erro de
uso. Guarde o CSV do IEDB junto do projeto: ele é a evidência de que os números locais podem
substituir os do site naquele conjunto de dados.
"""
import argparse
import csv
import io
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from poa.services import conservancy_client as cc

# Tolerância em pontos percentuais. O IEDB arredonda para 2 casas; divergências acima disso são
# diferença de método, não de formatação.
TOLERANCIA = 0.011


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        add_help=False,
        description="Compara a conservancia local com um CSV exportado do IEDB.")
    parser._action_groups.pop()
    required = parser.add_argument_group("required arguments")
    optional = parser.add_argument_group("optional arguments")

    required.add_argument("-h", "--help", action="help", default=argparse.SUPPRESS,
                          help="Print this help message.")
    required.add_argument("-i", help="CSV exportado do IEDB", required=True, type=str)
    required.add_argument("-e", help="FASTA de epitopos usado na analise", required=True, type=str)
    required.add_argument("-p", help="FASTA de proteinas (conjunto de comparacao)",
                          required=True, type=str)
    required.add_argument("-t", help="Limiar de identidade usado no site", required=True,
                          type=float)
    optional.add_argument("--operador", help="'>=' (conservados) ou '<' (unicos). Padrao: >=",
                          type=str, default=">=", choices=[">=", "<"])
    optional.add_argument("--tolerancia", help=f"Pontos percentuais (padrao = {TOLERANCIA})",
                          type=float, default=TOLERANCIA)
    return parser


def _num(valor) -> float:
    """Converte '67.05%' ou '50.00% (1/2)' no número da frente."""
    texto = str(valor).strip().replace("%", "")
    return float(texto.split("(")[0].split()[0]) if texto else float("nan")


def ler_csv(texto: str):
    """Indexa as linhas de um CSV de conservância pela sequência do epítopo."""
    linhas = {}
    leitor = csv.DictReader(io.StringIO(texto))
    coluna_pct = next((c for c in (leitor.fieldnames or [])
                       if c and c.startswith("Percent of protein")), None)
    if coluna_pct is None:
        raise SystemExit("O CSV nao tem a coluna 'Percent of protein sequence matches...'. "
                         "Confira se e mesmo uma exportacao do Epitope Conservancy Analysis.")
    for linha in leitor:
        chave = str(linha.get("Epitope sequence", "")).strip().upper()
        if chave:
            linhas[chave] = {
                "nome": linha.get("Epitope name", ""),
                "percent": _num(linha.get(coluna_pct)),
                "min": _num(linha.get("Minimum identity")),
                "max": _num(linha.get("Maximum identity")),
            }
    return linhas, coluna_pct


def main(argv=None):
    args = build_parser().parse_args(argv)

    iedb, coluna_iedb = ler_csv(open(args.i, encoding="utf-8", errors="replace").read())
    local, coluna_local = ler_csv(cc.run_conservancy(args.e, args.p, args.t, args.operador))

    print("COLUNAS")
    print(f"  IEDB  : {coluna_iedb}")
    print(f"  local : {coluna_local}")
    if coluna_iedb != coluna_local:
        print("  AVISO: os cabecalhos diferem. Confira se o limiar e o operador sao os mesmos;")
        print("         o IEDB escreve '<=' onde o POA escreve '<'.")

    so_iedb = sorted(set(iedb) - set(local))
    so_local = sorted(set(local) - set(iedb))
    comuns = sorted(set(iedb) & set(local))

    print(f"\nEPITOPOS\n  no IEDB: {len(iedb)} | local: {len(local)} | em ambos: {len(comuns)}")
    for seq in so_iedb[:5]:
        print(f"  so no IEDB : {seq[:40]}")
    for seq in so_local[:5]:
        print(f"  so no local: {seq[:40]}")

    divergentes = []
    for seq in comuns:
        a, b = iedb[seq], local[seq]
        for campo in ("percent", "min", "max"):
            if abs(a[campo] - b[campo]) > args.tolerancia:
                divergentes.append((seq, a["nome"], campo, a[campo], b[campo]))

    print(f"\nCOMPARACAO (tolerancia {args.tolerancia} pp)")
    if divergentes:
        print(f"  {len(divergentes)} divergencia(s) em {len(comuns)} epitopo(s):\n")
        print(f"  {'epitopo':<26} {'campo':<8} {'IEDB':>9} {'local':>9} {'delta':>9}")
        print("  " + "-" * 64)
        for seq, nome, campo, v_iedb, v_local in divergentes[:30]:
            print(f"  {(nome or seq)[:25]:<26} {campo:<8} {v_iedb:>9.2f} {v_local:>9.2f} "
                  f"{v_local - v_iedb:>+9.2f}")
        if len(divergentes) > 30:
            print(f"  ... e mais {len(divergentes) - 30}")
    else:
        print(f"  nenhuma: os {len(comuns)} epitopos batem nos tres campos.")

    print()
    if divergentes or so_iedb or so_local:
        print("RESULTADO: a reimplementacao local DIVERGE do IEDB neste conjunto. Nao use os "
              "numeros locais sem entender a causa.")
        return 1
    print("RESULTADO: equivalente ao IEDB neste conjunto de dados. Guarde o CSV do IEDB como "
          "evidencia.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
