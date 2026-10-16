#!/usr/bin/env python
"""Conferir se uma estrutura .pdb corresponde à sequência usada nas predições.

A mesma verificação que a etapa 7 da interface faz, em linha de comando, para poder ser usada
sobre arquivos já gravados sem subir o Streamlit:

    python tools/check_structure.py -p estrutura.pdb -f referencia.fasta
    python tools/check_structure.py -p estrutura.pdb -f referencia.fasta -e epitopos.fasta -c A

Sai com código 0 quando a estrutura serve, 1 quando não serve e 2 em erro de uso — então também
funciona dentro de um script.
"""
import argparse
import os
import sys

# Com PYTHONSAFEPATH=1 o diretório do script não entra no sys.path, e a raiz do repositório
# (um nível acima de tools/) é que contém o pacote ``poa``.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from poa.core import structure_map
from poa.services import database_client


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        add_help=False,
        description="Confere se um .pdb e a sequencia das predicoes sao da mesma proteina.")
    parser._action_groups.pop()
    required = parser.add_argument_group("required arguments")
    optional = parser.add_argument_group("optional arguments")

    required.add_argument("-h", "--help", action="help", default=argparse.SUPPRESS,
                          help="Print this help message.")
    required.add_argument("-p", help="Arquivo .pdb da estrutura", required=True, type=str)
    required.add_argument("-f", help="FASTA da sequencia usada nas predicoes (o -f do POA)",
                          required=True, type=str)
    optional.add_argument("-e", help="FASTA de epitopos do POA1 (nomes Especie_Prot_Metodo_ini_fim)",
                          type=str, default="")
    optional.add_argument("-c", help="Cadeia a conferir (padrao: a primeira do arquivo)",
                          type=str, default="")
    optional.add_argument("-i", help="Identidade minima em %% (padrao = 70)", type=float,
                          default=structure_map.DEFAULT_MIN_IDENTITY * 100)
    return parser


def read_ranges(path):
    """Intervalos (inicio, fim) a partir dos nomes de um FASTA de epitopos do POA1."""
    from Bio import SeqIO

    ranges = []
    for record in SeqIO.parse(path, "fasta"):
        fields = str(record.id).split("_")
        try:
            ranges.append((int(fields[-2]), int(fields[-1])))
        except (IndexError, ValueError):
            print(f"  aviso: nome fora da convencao, ignorado: {record.id}")
    return ranges


def main(argv=None):
    from Bio import SeqIO

    args = build_parser().parse_args(argv)

    pdb_text = open(args.p, encoding="utf-8", errors="replace").read()
    records = list(SeqIO.parse(args.f, "fasta"))
    if not records:
        print(f"ERRO: {args.f} nao tem nenhum registro.", file=sys.stderr)
        return 2
    if len(records) > 1:
        print(f"Aviso: {args.f} tem {len(records)} registros; conferindo o primeiro "
              f"({records[0].id}). Use um FASTA por proteina.")
    reference = str(records[0].seq).upper()

    summary = database_client.structure_summary(pdb_text)
    print("ESTRUTURA")
    print(f"  arquivo  : {args.p}")
    for label, key in (("molecula", "molecule"), ("organismo", "organism"),
                       ("metodo", "method"), ("resolucao", "resolution"),
                       ("cadeias", "chains")):
        if summary.get(key):
            print(f"  {label:<9}: {summary[key]}")

    print("\nREFERENCIA")
    print(f"  {records[0].id}  ({len(reference)} aa)")

    ranges = read_ranges(args.e) if args.e else [(1, len(reference))]
    report = structure_map.map_epitopes(pdb_text, reference, ranges,
                                        chain=args.c or None,
                                        min_identity=args.i / 100.0)

    print("\nVERIFICACAO")
    print(f"  cadeia conferida     : {report.chain} ({report.structure_length} residuos observados)")
    print(f"  identidade           : {report.identity_pct:.1f}%  (minimo {args.i:.0f}%)")
    print(f"  metodo               : {report.method}")
    if args.e:
        print(f"  epitopos             : {len(ranges)}")
        print(f"  mapeados             : {len(report.ranges)}")
        if report.partial:
            print(f"  parciais no cristal  : {len(report.partial)}")
        if report.unmapped:
            print(f"  sem correspondencia  : {len(report.unmapped)}")

    for warning in report.warnings:
        print(f"\n  AVISO: {warning}")

    print()
    if report.ok:
        print("RESULTADO: a estrutura corresponde a sequencia. O destaque pode ser confiado.")
        return 0
    print("RESULTADO: a estrutura NAO corresponde. Destacar epitopos nela produziria uma "
          "figura sem significado.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
