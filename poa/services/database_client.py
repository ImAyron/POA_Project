"""Database search: a species name in, the FASTA the pipeline consumes out.

This is the stage that precedes POA1. The rest of the pipeline starts from a protein FASTA that
somebody had to find by hand; this module finds it, for any organism, from four public services:

* **UniProt** (``rest.uniprot.org``) — taxonomy lookup, protein entries, sequences and the
  ``Chain``/``Domain`` features that delimit mature proteins and domains.
* **InterPro** (``ebi.ac.uk/interpro/api``) — Pfam domains with coordinates. UniProt cross-
  references Pfam but does not carry its coordinates, so the boundaries come from here.
* **PDBe SIFTS** (``ebi.ac.uk/pdbe/api/mappings``) — the residue-level correspondence between a
  UniProt sequence and a PDB entry. This is what makes a structure usable: it says *exactly*
  which UniProt residues a structure covers and how they are numbered in the file.
* **RCSB** (``files.rcsb.org``) — the structure file itself.

Nothing here writes files; :mod:`poa.cli.poa0` does that. Every call raises
:class:`ServiceUnavailable` on a network problem so the caller can fall back to a manual FASTA,
exactly like the prediction clients.

Scope note: UniProt's ``organism_id`` matches one taxon, ``taxonomy_id`` includes its descendants.
The difference is large for virus species with many isolates, so it is an explicit parameter
rather than a hidden default.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence

from ..logging_conf import get_logger
from .base import ServiceUnavailable

logger = get_logger("dbsearch")

UNIPROT_API = "https://rest.uniprot.org"
INTERPRO_API = "https://www.ebi.ac.uk/interpro/api"
PDBE_SIFTS_API = "https://www.ebi.ac.uk/pdbe/api/mappings/uniprot"
PDBE_BEST_STRUCTURES = "https://www.ebi.ac.uk/pdbe/api/mappings/best_structures"
RCSB_FILES = "https://files.rcsb.org/download"

USER_AGENT = "POA/1.0 (Pipeline de Otimizacao de Antigenos)"

# Region sources, most to least specific. A 'Chain' is a mature protein cleaved from a
# polyprotein and is what a vaccine target usually is; 'PDB' is the only source that guarantees
# an experimentally solved structure exists for the region.
REGION_SOURCES = ("Chain", "Domain", "Pfam", "PDB")


# --------------------------------------------------------------------------- data types
@dataclass
class Taxon:
    """A candidate organism from the taxonomy search."""

    taxon_id: int
    scientific_name: str
    common_name: str = ""
    rank: str = ""

    def __str__(self) -> str:
        extra = f" ({self.common_name})" if self.common_name else ""
        return f"{self.taxon_id}  {self.scientific_name}{extra}"


@dataclass
class Protein:
    """A UniProt entry, without its sequence (listings stay cheap)."""

    accession: str
    name: str
    length: int
    reviewed: bool
    pdb_ids: List[str] = field(default_factory=list)

    @property
    def database(self) -> str:
        return "SwissProt" if self.reviewed else "TrEMBL"

    def __str__(self) -> str:
        pdb = f"  PDB: {','.join(self.pdb_ids[:3])}" if self.pdb_ids else ""
        return f"{self.accession}  {self.length:>5} aa  {self.database:<10} {self.name}{pdb}"


@dataclass
class Region:
    """A sub-sequence of a protein worth predicting epitopes on."""

    name: str
    start: int
    end: int
    source: str                 # one of REGION_SOURCES
    description: str = ""

    @property
    def length(self) -> int:
        return self.end - self.start + 1

    def __str__(self) -> str:
        return f"[{self.source:<12}] {self.start:>5}-{self.end:<5} {self.length:>5} aa  {self.name}"


@dataclass
class ChainMapping:
    """
    SIFTS correspondence between one PDB chain and a stretch of a UniProt sequence.

    PDBe reports two different numbers for the same residue and they are easy to confuse:
    ``residue_number`` is the 1-based index within the entity, while
    ``author_residue_number`` is what the coordinate file actually says. For PDB 6H3X they are
    1 and 482 for the same residue — taking the former for the latter shifts every epitope by
    481 positions. ``pdb_start``/``pdb_end`` here are always the author numbering, the one a
    viewer needs; the entity indices are kept separately.
    """

    pdb_id: str
    chain_id: str
    unp_start: int
    unp_end: int
    pdb_start: int
    pdb_end: int
    identity: float
    accession: str = ""
    entity_start: int = 0
    entity_end: int = 0
    derived_end: bool = False       # pdb_end was computed, not reported

    @property
    def offset(self) -> int:
        """Add this to a UniProt position to get the author residue number in the file."""
        return self.pdb_start - self.unp_start


@dataclass
class Structure:
    """
    A PDB entry covering part of a UniProt sequence, as returned by the *listing*.

    Deliberately does not carry PDB residue numbers: obtaining those costs one request per
    entry, and an entry like Zika's polyprotein has hundreds. Call :func:`sifts_mappings` for
    the one structure actually chosen.
    """

    pdb_id: str
    accession: str
    chains: List[str] = field(default_factory=list)
    unp_start: int = 0
    unp_end: int = 0
    coverage: float = 0.0           # fraction of the UniProt sequence covered
    resolution: Optional[float] = None

    @property
    def span(self) -> str:
        return f"{self.unp_start}-{self.unp_end}"

    def __str__(self) -> str:
        res = f"{self.resolution:.2f} A" if self.resolution else "n/d"
        return (f"{self.pdb_id}  UniProt {self.span:<12} {self.coverage * 100:>5.1f}% da sequencia"
                f"  {res:>8}  cadeia(s) {','.join(self.chains)}")


# --------------------------------------------------------------------------- HTTP
def _get_json(url: str, params: Optional[dict] = None, timeout: int = 60):
    """
    GET a JSON document.

    Args:
        url: absolute URL.
        params: query string parameters.
        timeout: seconds.

    Returns:
        The decoded JSON document.

    Raises:
        ServiceUnavailable: the host could not be reached or did not return JSON.
    """
    try:
        import requests
    except ImportError as exc:  # pragma: no cover - requests is a declared dependency
        raise ServiceUnavailable("The 'requests' package is required for database search.") from exc

    headers = {"User-Agent": USER_AGENT, "Accept": "application/json"}
    try:
        resp = requests.get(url, params=params, headers=headers, timeout=timeout)
        resp.raise_for_status()
        return resp.json()
    except requests.RequestException as exc:
        raise ServiceUnavailable(f"Could not reach {url}: {exc}") from exc
    except ValueError as exc:
        raise ServiceUnavailable(f"{url} did not return JSON: {exc}") from exc


def _get_text(url: str, timeout: int = 120) -> str:
    """GET a plain-text document (a structure file). Raises ServiceUnavailable on failure."""
    try:
        import requests
    except ImportError as exc:  # pragma: no cover
        raise ServiceUnavailable("The 'requests' package is required for database search.") from exc

    try:
        resp = requests.get(url, headers={"User-Agent": USER_AGENT}, timeout=timeout)
        resp.raise_for_status()
        return resp.text
    except requests.RequestException as exc:
        raise ServiceUnavailable(f"Could not download {url}: {exc}") from exc


# --------------------------------------------------------------------------- taxonomy
def search_taxa(name: str, limit: int = 10) -> List[Taxon]:
    """
    Find candidate organisms by name.

    A name is almost never unique — "Zika virus" matches the species, isolates and a synthetic
    construct — so this returns candidates for the caller to choose from rather than guessing.

    Args:
        name: free-text organism name, e.g. "Oropouche virus".
        limit: maximum candidates.

    Returns:
        Candidates, exact scientific-name matches first.
    """
    data = _get_json(f"{UNIPROT_API}/taxonomy/search",
                     {"query": name, "size": min(limit, 25), "format": "json"})
    taxa = [
        Taxon(taxon_id=t["taxonId"],
              scientific_name=t.get("scientificName", ""),
              common_name=t.get("commonName", "") or "",
              rank=t.get("rank", "") or "")
        for t in data.get("results", [])
    ]
    # An exact name match is what the user typed; UniProt does not always rank it first.
    wanted = name.strip().lower()
    taxa.sort(key=lambda t: (t.scientific_name.lower() != wanted, len(t.scientific_name)))
    logger.debug("taxonomy '%s' -> %d candidate(s)", name, len(taxa))
    return taxa[:limit]


def resolve_taxon(name: str) -> Taxon:
    """
    Return the single best taxon for a name.

    Raises:
        ServiceUnavailable: nothing matched.
    """
    taxa = search_taxa(name, limit=10)
    if not taxa:
        raise ServiceUnavailable(f"No organism in UniProt matches '{name}'.")
    return taxa[0]


# --------------------------------------------------------------------------- proteins
def search_proteins(taxon_id: int, *, include_descendants: bool = False, limit: int = 25,
                    reviewed_first: bool = True) -> List[Protein]:
    """
    List the protein entries of an organism.

    Args:
        taxon_id: NCBI/UniProt taxon identifier.
        include_descendants: ``taxonomy_id`` (strains and isolates too) instead of the exact
            ``organism_id``. For viruses this can multiply the result count by an order of
            magnitude, so it is off by default.
        limit: maximum entries.
        reviewed_first: put SwissProt entries before TrEMBL ones. SwissProt entries carry the
            mature proteins as ``Chain`` features, which TrEMBL entries do not.

    Returns:
        The entries found.
    """
    field_name = "taxonomy_id" if include_descendants else "organism_id"
    # Ranking can only reorder what was fetched: asking for exactly `limit` entries and then
    # sorting them would hide the one entry with a structure behind an arbitrary page boundary.
    # That is how A0A0D4BSW3 — the only Oropouche polyprotein with a PDB — fell off a 4-row list.
    fetch_size = min(max(limit * 5, 50), 500) if reviewed_first else min(limit, 500)
    data = _get_json(f"{UNIPROT_API}/uniprotkb/search", {
        "query": f"{field_name}:{taxon_id}",
        "size": fetch_size,
        "format": "json",
        "fields": "accession,protein_name,length,reviewed,xref_pdb",
    })

    proteins = []
    for entry in data.get("results", []):
        proteins.append(Protein(
            accession=entry["primaryAccession"],
            name=_entry_name(entry),
            length=entry.get("sequence", {}).get("length", 0),
            reviewed="unreviewed" not in entry.get("entryType", "").lower(),
            pdb_ids=[x["id"] for x in entry.get("uniProtKBCrossReferences", [])
                     if x.get("database") == "PDB"],
        ))

    if reviewed_first:
        # A structure is what makes a target usable, so entries that have one come first among
        # equals; otherwise the choice would be the arbitrary order UniProt returned.
        proteins.sort(key=lambda p: (not p.reviewed, not p.pdb_ids, -p.length))
    logger.debug("taxon %s -> %d protein(s)", taxon_id, len(proteins))
    return proteins[:limit]


def _entry_name(entry: dict) -> str:
    """Best available protein name: recommended, else submitted, else the accession."""
    desc = entry.get("proteinDescription", {})
    rec = desc.get("recommendedName", {}).get("fullName", {}).get("value")
    if rec:
        return rec
    for sub in desc.get("submissionNames", []):
        value = sub.get("fullName", {}).get("value")
        if value:
            return value
    return entry.get("primaryAccession", "?")


def get_entry(accession: str) -> dict:
    """Return the full UniProt JSON entry for an accession."""
    return _get_json(f"{UNIPROT_API}/uniprotkb/{accession}.json")


def get_sequence(accession: str, entry: Optional[dict] = None) -> str:
    """Return the amino-acid sequence of an entry as a plain string."""
    entry = entry or get_entry(accession)
    return entry.get("sequence", {}).get("value", "")


# --------------------------------------------------------------------------- regions
def list_regions(accession: str, entry: Optional[dict] = None,
                 sources: Sequence[str] = REGION_SOURCES,
                 structure_limit: int = 10) -> List[Region]:
    """
    Collect the candidate sub-sequences of a protein from every source that has one.

    A polyprotein is not a vaccine target; one of its mature proteins or domains is. Which
    source carries that information depends on the entry: SwissProt entries annotate mature
    proteins as ``Chain``, while TrEMBL entries (most virus isolates) have only ``Domain``,
    Pfam and whatever a structure covers.

    Args:
        accession: UniProt accession.
        entry: a previously fetched entry, to avoid a second request.
        sources: which of :data:`REGION_SOURCES` to consult.
        structure_limit: how many structures the ``PDB`` source may contribute. Well-studied
            proteins have hundreds, which would bury the handful of biologically meaningful
            regions under near-duplicate crystal forms.

    Returns:
        Regions sorted by start position, deduplicated on (start, end, name).
    """
    entry = entry or get_entry(accession)
    full_length = entry.get("sequence", {}).get("length", 0)
    regions: List[Region] = []

    for feature in entry.get("features", []):
        ftype = feature.get("type")
        if ftype not in ("Chain", "Domain") or ftype not in sources:
            continue
        start, end = _feature_span(feature)
        if start is None:
            continue
        # A SwissProt polyprotein annotates itself as a Chain spanning the whole sequence; that
        # is the protein we are already looking at, not a sub-region of it.
        if ftype == "Chain" and start == 1 and end == full_length:
            continue
        regions.append(Region(name=feature.get("description") or ftype,
                              start=start, end=end, source=ftype,
                              description=feature.get("description", "")))

    if "Pfam" in sources:
        regions.extend(_pfam_regions(accession))
    if "PDB" in sources:
        for structure in list_structures(accession, limit=structure_limit):
            regions.append(Region(
                name=structure.pdb_id,
                start=structure.unp_start, end=structure.unp_end, source="PDB",
                description=f"estrutura {structure.pdb_id}, cadeia(s) "
                            f"{','.join(structure.chains)}, "
                            f"{structure.coverage * 100:.1f}% da sequencia"))

    # Two sources frequently describe the same span: a TrEMBL entry's Domain features are derived
    # from Pfam, so both report it. Collapsing them to one region but keeping both source names
    # avoids claiming a region has less support than it has.
    merged: Dict[tuple, Region] = {}
    for region in sorted(regions, key=lambda r: (r.start, r.end, r.source)):
        key = (region.start, region.end, region.name.lower())
        if key in merged:
            previous = merged[key]
            if region.source not in previous.source.split("+"):
                previous.source = f"{previous.source}+{region.source}"
        else:
            merged[key] = region
    unique = list(merged.values())
    logger.debug("%s -> %d region(s)", accession, len(unique))
    return unique


def _feature_span(feature: dict):
    """Return (start, end) of a UniProt feature, or (None, None) when either end is unknown."""
    location = feature.get("location", {})
    start = location.get("start", {}).get("value")
    end = location.get("end", {}).get("value")
    if start is None or end is None:
        return None, None
    return int(start), int(end)


def _pfam_regions(accession: str) -> List[Region]:
    """Pfam domains with coordinates, via InterPro. Returns [] when the entry has none."""
    url = f"{INTERPRO_API}/entry/pfam/protein/uniprot/{accession}/"
    try:
        data = _get_json(url)
    except ServiceUnavailable:
        # InterPro answers 404 for a protein with no Pfam match, which is an answer, not an outage.
        logger.debug("no Pfam entries for %s", accession)
        return []

    regions = []
    for result in data.get("results", []):
        meta = result.get("metadata", {})
        for protein in result.get("proteins", []):
            for location in protein.get("entry_protein_locations", []):
                for fragment in location.get("fragments", []):
                    regions.append(Region(
                        name=meta.get("name") or meta.get("accession", "Pfam"),
                        start=int(fragment["start"]), end=int(fragment["end"]),
                        source="Pfam", description=meta.get("accession", "")))
    return regions


# --------------------------------------------------------------------------- structures
def list_structures(accession: str, limit: int = 25) -> List[Structure]:
    """
    Return the PDB entries that cover part of a UniProt sequence.

    Uses PDBe's ``best_structures`` endpoint, which answers for the whole accession in one
    request. Walking the entry's PDB cross-references instead would cost one request per
    structure — for Zika's polyprotein (281 entries, 600 chain mappings) that is minutes of
    sequential HTTP for a list the user only wants to glance at.

    Args:
        accession: UniProt accession.
        limit: maximum entries returned.

    Returns:
        One :class:`Structure` per PDB entry, best coverage first and better resolution breaking
        ties. Empty when PDBe knows no structure for the accession.
    """
    try:
        data = _get_json(f"{PDBE_BEST_STRUCTURES}/{accession}")
    except ServiceUnavailable:
        logger.debug("no structures known for %s", accession)
        return []

    by_entry: Dict[str, Structure] = {}
    for hits in data.values():
        for hit in hits:
            pdb_id = str(hit.get("pdb_id", "")).upper()
            if not pdb_id:
                continue
            start, end = int(hit.get("unp_start", 0)), int(hit.get("unp_end", 0))
            structure = by_entry.get(pdb_id)
            if structure is None:
                structure = Structure(pdb_id=pdb_id, accession=accession,
                                      unp_start=start, unp_end=end,
                                      coverage=float(hit.get("coverage") or 0.0),
                                      resolution=hit.get("resolution"))
                by_entry[pdb_id] = structure
            else:
                # Several chains of one entry: the entry covers the union of their spans.
                structure.unp_start = min(structure.unp_start, start)
                structure.unp_end = max(structure.unp_end, end)
            chain = hit.get("chain_id")
            if chain and chain not in structure.chains:
                structure.chains.append(chain)

    structures = sorted(by_entry.values(),
                        key=lambda s: (-s.coverage, s.resolution or 99.0, s.pdb_id))
    for structure in structures:
        structure.chains.sort()
    logger.debug("%s -> %d structure(s)", accession, len(structures))
    return structures[:limit]


def sifts_mappings(pdb_id: str, accession: Optional[str] = None) -> List[ChainMapping]:
    """
    Residue-level correspondence between a PDB entry and the UniProt sequence(s) it covers.

    This is the authoritative answer to "which residue of my sequence is residue N of this
    structure". Crystal structures routinely start their numbering somewhere other than 1 —
    at a construct tag, at the first ordered residue, or at the position in the parent
    polyprotein — and guessing the offset silently mislabels every residue downstream.

    Args:
        pdb_id: four-character PDB identifier, case-insensitive.
        accession: keep only the mappings to this UniProt entry.

    Returns:
        One mapping per (chain, segment). Empty when PDBe has no UniProt mapping for the entry,
        which is normal for structures with no UniProt counterpart.
    """
    try:
        data = _get_json(f"{PDBE_SIFTS_API}/{pdb_id.lower()}")
    except ServiceUnavailable:
        logger.debug("no SIFTS mapping for %s", pdb_id)
        return []

    mappings = []
    for entry_id, body in data.items():
        for acc, info in body.get("UniProt", {}).items():
            if accession and acc != accession:
                continue
            for m in info.get("mappings", []):
                try:
                    start, end = m["start"], m["end"]
                    unp_start, unp_end = int(m["unp_start"]), int(m["unp_end"])
                    # Author numbering is what the coordinate file uses; the entity index is a
                    # fallback only when PDBe omits the author number.
                    author_start = start.get("author_residue_number")
                    pdb_start = int(author_start if author_start is not None
                                    else start["residue_number"])
                    author_end = end.get("author_residue_number")
                    derived = author_end is None
                    # PDBe leaves the author number of the last residue null for some entries
                    # (6H3X among them). Author numbering is contiguous across a SIFTS segment
                    # by construction, so the end follows from the start and the segment length.
                    pdb_end = int(author_end if not derived
                                  else pdb_start + (unp_end - unp_start))
                    mapping = ChainMapping(
                        pdb_id=entry_id.upper(),
                        chain_id=m.get("chain_id") or m.get("struct_asym_id", "A"),
                        unp_start=unp_start, unp_end=unp_end,
                        pdb_start=pdb_start, pdb_end=pdb_end,
                        identity=float(m.get("identity", 0.0)),
                        accession=acc,
                        entity_start=int(start.get("residue_number") or 0),
                        entity_end=int(end.get("residue_number") or 0),
                        derived_end=derived,
                    )
                except (KeyError, TypeError, ValueError):
                    continue
                mappings.append(mapping)
    return mappings


def fetch_structure(pdb_id: str) -> str:
    """
    Download a structure from RCSB in PDB format.

    Args:
        pdb_id: four-character PDB identifier.

    Returns:
        The file contents.

    Raises:
        ServiceUnavailable: the download failed. Very large structures exist only as mmCIF;
        those fail here and the caller should offer the mmCIF URL instead.
    """
    text = _get_text(f"{RCSB_FILES}/{pdb_id.upper()}.pdb")
    logger.debug("downloaded %s (%d bytes)", pdb_id, len(text))
    return text


def structure_summary(pdb_text: str) -> Dict[str, str]:
    """
    Read the identity of a structure out of its own header.

    Used to tell the user what a structure actually *is* before it is used as if it were their
    antigen. A PDB identifier carries no meaning on sight, so a file can be mistaken for an
    unrelated protein with no visible symptom.

    Args:
        pdb_text: the contents of a PDB file.

    Returns:
        Whatever of ``title``, ``molecule``, ``organism``, ``method``, ``resolution``,
        ``chains`` could be read. Missing fields are absent from the dict.
    """
    out: Dict[str, str] = {}
    chains = set()
    for line in pdb_text.splitlines():
        record = line[:6].strip()
        if record == "TITLE" and "title" not in out:
            out["title"] = line[10:].strip()
        elif record == "COMPND" and "MOLECULE:" in line:
            out.setdefault("molecule", line.split("MOLECULE:", 1)[1].strip().rstrip(";"))
        elif record == "SOURCE" and "ORGANISM_SCIENTIFIC:" in line:
            out.setdefault("organism",
                           line.split("ORGANISM_SCIENTIFIC:", 1)[1].strip().rstrip(";"))
        elif record == "EXPDTA" and "method" not in out:
            out["method"] = line[10:].strip()
        elif record == "REMARK" and "RESOLUTION." in line and "resolution" not in out:
            out["resolution"] = line.split("RESOLUTION.", 1)[1].strip()
        elif record == "SEQRES":
            chains.add(line[11])
        elif record == "ATOM":
            break          # the header is over; the rest is coordinates
    if chains:
        out["chains"] = ",".join(sorted(chains))
    return out
