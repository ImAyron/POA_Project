"""POA0 command-line interface — database search, in the argument style of POA1/POA2.

Four listing modes answer the four questions that used to be answered by hand, each printing the
options and stopping, so the choice is always explicit:

    python POA0_v1.0.py -s "Oropouche virus"                      # qual taxon?
    python POA0_v1.0.py -s "Oropouche virus" -l proteinas         # qual proteina?
    python POA0_v1.0.py -a A0A0D4BSW3 -l regioes                  # qual regiao?
    python POA0_v1.0.py -a A0A0D4BSW3 -l estruturas               # qual estrutura?

and then the run that writes the files:

    python POA0_v1.0.py -a A0A0D4BSW3 -r 482-703 --nome GcHead \\
                        --sigla OROV --estrutura 6H3X -d entradas
"""
from __future__ import annotations

import argparse
import sys

from ..core import dbsearch
from ..logging_conf import setup_logging

LIST_MODES = ("taxons", "proteinas", "regioes", "estruturas")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        add_help=False,
        description="POA0 - busca em bancos de dados: de uma especie ao FASTA que o POA1 consome.")
    parser._action_groups.pop()
    required = parser.add_argument_group("required arguments (-s ou -a)")
    optional = parser.add_argument_group("optional arguments")

    required.add_argument("-h", "--help", action="help", default=argparse.SUPPRESS,
                          help="Print this help message.")
    required.add_argument("-s", help="Nome do organismo, ex.: \"Oropouche virus\"",
                          type=str, default="")
    required.add_argument("-a", help="Acessao UniProt, ex.: A0A0D4BSW3 (dispensa -s)",
                          type=str, default="")
    required.add_argument("-d", help="Diretorio de saida", type=str, default="")

    optional.add_argument("-l", help=f"Modo listagem: {', '.join(LIST_MODES)}. "
                                     f"Mostra as opcoes e encerra.",
                          type=str, choices=LIST_MODES, default="")
    optional.add_argument("-p", help="Proteina escolhida (acessao ou parte do nome)",
                          type=str, default="")
    optional.add_argument("-r", help="Regiao escolhida: nome, id PDB ou intervalo 'inicio-fim'. "
                                     "Sem -r, usa a sequencia inteira.",
                          type=str, default="")
    optional.add_argument("--nome", help="Campo de proteina no cabecalho, ex.: GcHead "
                                         "(padrao: o nome da regiao, sem pontuacao)",
                          type=str, default=None)
    optional.add_argument("--sigla", help="Campo de especie no cabecalho, ex.: OROV "
                                          "(padrao: a primeira palavra do organismo)",
                          type=str, default="")
    optional.add_argument("--estrutura", help="Id PDB a baixar do RCSB, ex.: 6H3X",
                          type=str, default="")
    optional.add_argument("--descendentes", help="Incluir taxons descendentes (cepas e isolados)",
                          action="store_true")
    optional.add_argument("--limite", help="Maximo de itens nas listagens (padrao = 25)",
                          type=int, default=25)
    return parser


def _print_list(title: str, items) -> int:
    """Print a listing and report how many items it had."""
    print(f"\n{title}\n" + "-" * len(title))
    if not items:
        print("  (nada encontrado)")
        return 0
    for item in items:
        print(f"  {item}")
    print(f"\n{len(items)} item(ns).")
    return len(items)


def _run_listing(args) -> int:
    """Handle the ``-l`` modes. Returns the process exit code."""
    if args.l == "taxons":
        if not args.s:
            print("Para listar taxons e preciso -s com o nome do organismo.", file=sys.stderr)
            return 2
        _print_list(f"Taxons para '{args.s}'", dbsearch.list_taxa(args.s))
        return 0

    if args.l == "proteinas":
        if not args.s:
            print("Para listar proteinas e preciso -s com o nome do organismo.", file=sys.stderr)
            return 2
        taxon = dbsearch.db.resolve_taxon(args.s)
        print(f"Organismo: {taxon}")
        _print_list(f"Proteinas de {taxon.scientific_name}",
                    dbsearch.list_proteins(taxon.taxon_id,
                                           include_descendants=args.descendentes,
                                           limit=args.limite))
        return 0

    accession = args.a
    if not accession:
        print("Para listar regioes ou estruturas e preciso -a com a acessao UniProt "
              "(use -l proteinas para descobri-la).", file=sys.stderr)
        return 2

    if args.l == "regioes":
        _print_list(f"Regioes de {accession}", dbsearch.list_regions(accession))
    else:
        _print_list(f"Estruturas de {accession}", dbsearch.list_structures(accession))
    return 0


def main(argv=None):
    setup_logging()
    parser = build_parser()
    args = parser.parse_args(argv)

    if not args.s and not args.a:
        parser.error("informe -s (nome do organismo) ou -a (acessao UniProt).")

    if args.l:
        return _run_listing(args)

    if not args.d:
        parser.error("-d e obrigatorio para gravar os arquivos "
                     "(use -l para apenas listar as opcoes).")

    report = dbsearch.run_poa0(args)

    print("\nPOA0 - arquivos gerados:")
    for label, path in report["arquivos"].items():
        print(f"  {label:<10} {path}")
    for warning in report.get("avisos", []):
        print(f"  AVISO: {warning}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
