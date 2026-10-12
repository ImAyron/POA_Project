"""POA0 — from an organism name to the files the rest of the pipeline consumes.

The stage that used to be done by hand: find the organism, pick the protein, pick the region
worth predicting on, and write the FASTA. Optionally download the structure that covers that
region, so the 3D view in the GUI shows the *right* protein instead of whatever file was at hand.

Header convention, which the rest of the pipeline depends on:

    >Protein_Specie_Accession

:func:`poa.core.topology.sameVirus` compares an epitope's species against field 1 of this header
and the protein name against field 0, and every field is obtained by splitting on ``_`` — so
neither token may contain an underscore. :func:`sanitize_token` enforces that.

One file is written for the predictions *and* for ``-f``. Predictions report positions relative
to the sequence they ran on, so feeding POA1/POA2 a different sequence (the whole polyprotein,
say) would offset every epitope. Using one file for both makes that mismatch impossible.
"""
from __future__ import annotations

import json
import os
import re
from dataclasses import asdict, is_dataclass
from typing import List, Optional

from ..logging_conf import get_logger, log_step
from ..services import database_client as db

logger = get_logger("poa0")

REPORT_NAME = "POA0_relatorio.json"


def sanitize_token(text: str, fallback: str = "SEQ") -> str:
    """
    Turn free text into a header field: no underscores, no spaces, no punctuation.

    The pipeline splits headers and epitope names on ``_``, so an underscore inside a protein or
    species name silently shifts every other field — the species would be read from the middle
    of the protein's name.

    Args:
        text: the raw name, e.g. "Envelope protein E".
        fallback: used when nothing usable remains.

    Returns:
        A token safe to place in a header, e.g. "EnvelopeproteinE".
    """
    cleaned = re.sub(r"[^0-9A-Za-z]+", "", str(text or ""))
    return cleaned or fallback


def species_token(organism: str) -> str:
    """Short, underscore-free species token from an organism name ("Oropouche virus" -> OROPOUCHE)."""
    first = str(organism or "").strip().split()[0] if str(organism or "").strip() else ""
    return sanitize_token(first, "ORG").upper()


def _json_default(value):
    if is_dataclass(value):
        return asdict(value)
    return str(value)


def write_fasta(path: str, header: str, sequence: str, width: int = 60) -> str:
    """Write a single-record FASTA and return its path."""
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write(f">{header}\n")
        for start in range(0, len(sequence), width):
            handle.write(sequence[start:start + width] + "\n")
    return path


# --------------------------------------------------------------------------- listing modes
def list_taxa(name: str) -> List[db.Taxon]:
    """Candidate organisms for a name."""
    return db.search_taxa(name)


def list_proteins(taxon_id: int, *, include_descendants: bool = False,
                  limit: int = 25) -> List[db.Protein]:
    """Protein entries of an organism."""
    return db.search_proteins(taxon_id, include_descendants=include_descendants, limit=limit)


def list_regions(accession: str) -> List[db.Region]:
    """Candidate sub-sequences of a protein."""
    return db.list_regions(accession)


def list_structures(accession: str) -> List[db.Structure]:
    """Structures covering part of a protein."""
    return db.list_structures(accession)


# --------------------------------------------------------------------------- generation
def resolve_region(accession: str, wanted: str,
                   entry: Optional[dict] = None) -> Optional[db.Region]:
    """
    Find the region a user named.

    Accepts a region name ("Envelope protein E"), a PDB identifier ("6H3X") or explicit
    coordinates ("482-703"), because which of those a user has depends on where they are
    starting from.

    Args:
        accession: UniProt accession.
        wanted: the region name, PDB id or ``start-end`` span.
        entry: a previously fetched UniProt entry.

    Returns:
        The matching region, or None.
    """
    span = re.fullmatch(r"\s*(\d+)\s*-\s*(\d+)\s*", wanted or "")
    if span:
        start, end = int(span.group(1)), int(span.group(2))
        return db.Region(name=f"{start}-{end}", start=start, end=end, source="manual",
                         description="intervalo informado na linha de comando")

    regions = db.list_regions(accession, entry=entry)
    target = str(wanted or "").strip().lower()
    for region in regions:
        if region.name.lower() == target:
            return region
    for region in regions:                      # partial match, longest region first
        if target and target in region.name.lower():
            return region
    return None


def run_poa0(args):
    """
    Run POA0: resolve the organism/protein/region and write the pipeline's input files.

    Expects the namespace :func:`poa.cli.poa0.build_parser` produces.

    Returns:
        dict: the report that is also written as ``POA0_relatorio.json``.
    """
    report = {"entrada": {k: v for k, v in vars(args).items() if v not in (None, "", False)},
              "arquivos": {}, "avisos": []}

    # --- organism and protein ------------------------------------------------------------
    accession = args.a
    organism = args.s or ""
    if not accession:
        with log_step(logger, "POA0/taxonomia"):
            taxon = db.resolve_taxon(args.s)
            organism = taxon.scientific_name
            report["taxon"] = asdict(taxon)
            logger.info("organismo: %s (taxon %s)", taxon.scientific_name, taxon.taxon_id)

        with log_step(logger, "POA0/proteinas"):
            proteins = db.search_proteins(taxon.taxon_id,
                                          include_descendants=args.descendentes, limit=25)
            if not proteins:
                raise SystemExit(f"Nenhuma proteina encontrada para '{args.s}'.")
            chosen = _pick_protein(proteins, args.p)
            if chosen is None:
                raise SystemExit(
                    f"Nenhuma proteina de '{args.s}' corresponde a '{args.p}'. "
                    f"Use -l proteinas para ver as opcoes.")
            accession = chosen.accession
            report["proteina"] = asdict(chosen)
            logger.info("proteina: %s (%s, %d aa)", chosen.accession, chosen.name, chosen.length)

    with log_step(logger, "POA0/sequencia"):
        entry = db.get_entry(accession)
        sequence = db.get_sequence(accession, entry=entry)
        if not sequence:
            raise SystemExit(f"A entrada {accession} nao tem sequencia.")
        organism = organism or _organism_of(entry)
        report.setdefault("proteina", {"accession": accession,
                                       "name": db._entry_name(entry),
                                       "length": len(sequence)})

    # --- region --------------------------------------------------------------------------
    region = None
    if args.r:
        region = resolve_region(accession, args.r, entry=entry)
        if region is None:
            raise SystemExit(
                f"Nenhuma regiao de {accession} corresponde a '{args.r}'. "
                f"Use -l regioes para ver as opcoes.")
    if region is None:
        region = db.Region(name=db._entry_name(entry), start=1, end=len(sequence),
                           source="completa", description="sequencia inteira")

    sub_sequence = sequence[region.start - 1:region.end]
    report["regiao"] = asdict(region)

    # --- tokens and output ---------------------------------------------------------------
    specie = sanitize_token(args.sigla or species_token(organism)).upper()
    protein_token = sanitize_token(args.nome or region.name, "PROT")
    header = f"{protein_token}_{specie}_{accession}"
    report["cabecalho"] = header

    if args.nome is None and protein_token != str(region.name):
        report["avisos"].append(
            f"O nome da regiao '{region.name}' virou o campo '{protein_token}' no cabecalho "
            f"(o pipeline divide os campos por '_'). Use --nome para escolher outro.")

    os.makedirs(args.d, exist_ok=True)
    fasta_path = os.path.join(args.d, f"{protein_token}_{specie}.fasta")
    write_fasta(fasta_path, header, sub_sequence)
    report["arquivos"]["fasta"] = fasta_path
    report["sequencia"] = {"inicio": region.start, "fim": region.end,
                           "tamanho": len(sub_sequence)}
    logger.info("FASTA: %s (%d aa, %s %d-%d)", fasta_path, len(sub_sequence),
                region.source, region.start, region.end)

    # --- structure -----------------------------------------------------------------------
    if args.estrutura:
        _add_structure(args, accession, region, report)

    report_path = os.path.join(args.d, REPORT_NAME)
    with open(report_path, "w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, ensure_ascii=False, default=_json_default)
    report["arquivos"]["relatorio"] = report_path
    return report


def _pick_protein(proteins: List[db.Protein], wanted: Optional[str]) -> Optional[db.Protein]:
    """Choose a protein by accession or name; the best-ranked one when nothing was asked for."""
    if not wanted:
        return proteins[0]
    target = wanted.strip().lower()
    for protein in proteins:
        if protein.accession.lower() == target:
            return protein
    for protein in proteins:
        if target in protein.name.lower():
            return protein
    return None


def _organism_of(entry: dict) -> str:
    return entry.get("organism", {}).get("scientificName", "")


def _add_structure(args, accession: str, region: db.Region, report: dict) -> None:
    """
    Download the chosen structure and record how its numbering relates to the region.

    The mapping is the point. A structure file carries its own residue numbering, which may or
    may not be the one the predictions use; recording the SIFTS correspondence here means the
    viewer never has to guess it.
    """
    with log_step(logger, "POA0/estrutura"):
        pdb_id = args.estrutura.upper()
        text = db.fetch_structure(pdb_id)
        path = os.path.join(args.d, f"structure_{pdb_id}.pdb")
        with open(path, "w", encoding="utf-8", newline="\n") as handle:
            handle.write(text)
        report["arquivos"]["estrutura"] = path

        summary = db.structure_summary(text)
        report["estrutura"] = {"pdb_id": pdb_id, "cabecalho": summary}
        logger.info("estrutura %s: %s | %s", pdb_id,
                    summary.get("molecule", "?"), summary.get("organism", "?"))

        mappings = [m for m in db.sifts_mappings(pdb_id, accession=accession)]
        if not mappings:
            report["avisos"].append(
                f"A estrutura {pdb_id} nao tem mapeamento SIFTS para {accession}: a numeracao "
                f"dos residuos sera resolvida por alinhamento na hora de visualizar.")
            return

        # Positions in the written FASTA are relative to the region, so the offset that matters
        # is from the region's own coordinates, not the full sequence's.
        chosen = max(mappings, key=lambda m: m.unp_end - m.unp_start)
        report["estrutura"]["sifts"] = asdict(chosen)
        report["estrutura"]["offset_regiao"] = chosen.offset + (region.start - 1)
        report["estrutura"]["cobre_regiao"] = (chosen.unp_start <= region.start
                                               and chosen.unp_end >= region.end)
        if not report["estrutura"]["cobre_regiao"]:
            report["avisos"].append(
                f"A estrutura {pdb_id} cobre {chosen.unp_start}-{chosen.unp_end} da sequencia, "
                f"mas a regiao escolhida e {region.start}-{region.end}: parte dos epitopos nao "
                f"tera correspondencia estrutural.")
