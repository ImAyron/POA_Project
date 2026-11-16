"""Map epitope positions onto the residue numbering of a structure, or refuse to.

Highlighting an epitope on a structure means coloring residue numbers. Those numbers are only
meaningful if the structure is of the protein the predictions were run on, numbered the way that
protein is numbered — and neither is guaranteed:

* A crystal structure of a domain often restarts numbering at 1 (PDB 6H3X numbers the Oropouche
  Gc head 1-222, where the polyprotein has it at 482-703).
* Expression constructs add tags, so chain A of PDB 4C23 begins at residue -5.
* Residues disordered in the crystal are absent from the file entirely.
* And a structure can simply be of a different protein.

The last case is the dangerous one, because nothing fails: 3Dmol.js colors whatever residue
numbers it is given. Running the 19 DENV epitopes of ``results/08102026-TESTANDO`` against 4C23
(L-fuculose kinase, *Streptococcus pneumoniae*) finds **0 of 19** epitope sequences in the
structure, yet all 19 position ranges fall inside chain A's numbering — so all 19 would be
painted red on an unrelated bacterial enzyme, with no error and no visible symptom.

This module answers "which residue of this structure is position N of my sequence" by comparing
the sequences themselves, and reports an unusable answer as unusable instead of a plausible
picture.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

from ..logging_conf import get_logger

logger = get_logger("structmap")

# Below this fraction of identical residues the structure is treated as a different protein and
# the mapping is refused rather than offered with a caveat.
DEFAULT_MIN_IDENTITY = 0.70

# ...and below this fraction of the structure actually aligning, the identity is not evidence of
# anything. Identity alone is a ratio over the aligned columns, so one matching residue out of a
# 150-residue chain scores 100%. Both end gaps are free (see _aligner), which is right for a
# structure covering part of a longer protein, but it also lets the aligner discard almost the
# whole structure as overhang and keep a tiny high-scoring overlap. Measured on this branch: of
# 216 pairs of unrelated random sequences, one was accepted at 100% identity over a single
# aligned column (0.7% of a 150-residue structure against a 500-aa reference) — exactly the
# wrong-structure case this module exists to refuse.
DEFAULT_MIN_COVERAGE = 0.50

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
    "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
    "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    # Modified residues common in crystal structures, mapped to their parent amino acid.
    "MSE": "M", "SEC": "U", "PYL": "O", "HYP": "P", "PCA": "Q", "SEP": "S", "TPO": "T",
    "PTR": "Y", "CSO": "C", "CME": "C", "KCX": "K", "LLP": "K", "MLY": "K",
}


@dataclass
class Residue:
    """One observed residue of a chain: what it is and how the file numbers it."""

    number: int
    icode: str
    letter: str

    @property
    def label(self) -> str:
        return f"{self.number}{self.icode.strip()}"


@dataclass
class MappingReport:
    """The outcome of mapping epitope positions onto a structure."""

    chain: str
    method: str                                     # 'sifts' | 'alignment'
    identity: float                                 # 0.0-1.0 over aligned columns
    ranges: List[Tuple[int, int]] = field(default_factory=list)      # in PDB numbering
    # The reference-sequence range each entry of ``ranges`` came from, index for index. ``ranges``
    # alone cannot be traced back: unmapped epitopes are left out of it, so position N of the
    # output is not epitope N of the input, and a viewer that wants to say *which* epitope a
    # colored stretch is had nothing to look it up by.
    sources: List[Tuple[int, int]] = field(default_factory=list)     # in reference numbering
    unmapped: List[Tuple[int, int]] = field(default_factory=list)    # no counterpart in structure
    partial: List[Tuple[int, int]] = field(default_factory=list)     # mapped, but gapped in the file
    warnings: List[str] = field(default_factory=list)
    ok: bool = True
    structure_length: int = 0
    aligned_columns: int = 0

    @property
    def identity_pct(self) -> float:
        return self.identity * 100.0

    def summary(self) -> str:
        """One line in Portuguese, for the interface and the logs."""
        if not self.ok:
            return (f"Estrutura incompativel: identidade {self.identity_pct:.1f}% com a sequencia "
                    f"das predicoes (cadeia {self.chain}). Destaque recusado.")
        method = "mapeamento SIFTS" if self.method == "sifts" else "alinhamento"
        parts = [f"cadeia {self.chain}", f"identidade {self.identity_pct:.1f}%", f"via {method}",
                 f"{len(self.ranges)} epitopo(s) mapeado(s)"]
        if self.partial:
            parts.append(f"{len(self.partial)} com residuos ausentes no cristal")
        if self.unmapped:
            parts.append(f"{len(self.unmapped)} sem correspondencia na estrutura")
        return "; ".join(parts) + "."


# --------------------------------------------------------------------------- parsing
def structure_format(text: str) -> str:
    """
    ``'cif'`` for an mmCIF file, ``'pdb'`` for anything else.

    The two are told apart by what only mmCIF has — a ``data_`` block holding an ``_atom_site``
    table — rather than by the file name: a structure arrives here as text, and the name it was
    uploaded under is somebody's choice.
    """
    return "cif" if "_atom_site." in text and "data_" in text[:4096] else "pdb"


def _cif_atoms(text: str):
    """
    Yield ``(group, atom, altloc, resname, chain, number, icode)`` for each atom of an mmCIF file.

    mmCIF carries two numberings side by side. ``label_seq_id`` counts residues of the entity
    from 1; ``auth_seq_id`` is the numbering the authors deposited, the one a PDB-format file
    of the same entry uses. For 6H3X they are 1 and 482 for the same residue. The ``auth_``
    fields are read — chain and residue number both — because that is what 3Dmol.js reads, and a
    residue number that means one thing here and another in the viewer colors the wrong residues
    with nothing to show for it.

    Only the first model is read, which is what the viewer draws.
    """
    columns: List[str] = []
    rows = False
    index: Dict[str, int] = {}
    first_model = None
    for line in text.splitlines():
        if not rows:
            if line.startswith("_atom_site."):
                columns.append(line.split()[0][len("_atom_site."):])
                continue
            if not columns:
                continue
            rows = True
            index = {name: position for position, name in enumerate(columns)}
        if not line or line[0] in "#_" or line.startswith(("loop_", "data_")):
            break
        fields = line.split()            # no value of this table contains a space
        if len(fields) != len(columns):
            continue

        def value(*names, default=""):
            for name in names:
                if name in index:
                    return fields[index[name]].strip("'\"")
            return default

        model = value("pdbx_PDB_model_num")
        if first_model is None:
            first_model = model
        elif model != first_model:
            continue
        try:
            number = int(value("auth_seq_id", "label_seq_id"))
        except ValueError:
            continue                     # '.', as on the rows of a ligand with no auth number
        altloc, icode = value("label_alt_id"), value("pdbx_PDB_ins_code")
        yield (value("group_PDB"), value("label_atom_id", "auth_atom_id"),
               " " if altloc in (".", "?", "") else altloc,
               value("label_comp_id", "auth_comp_id").upper(),
               value("auth_asym_id", "label_asym_id"), number,
               " " if icode in (".", "?", "") else icode)


def _ca_atoms(text: str):
    """Yield ``(chain, number, icode, resname, altloc)`` for each alpha-carbon record, either format."""
    if structure_format(text) == "cif":
        for _group, atom, altloc, resname, chain, number, icode in _cif_atoms(text):
            if atom == "CA":
                yield chain, number, icode, resname, altloc
        return
    for line in text.splitlines():
        if not line.startswith(("ATOM  ", "HETATM")) or len(line) < 54:
            continue
        if line[12:16].strip() != "CA":
            continue
        try:
            number = int(line[22:26])
        except ValueError:
            continue
        yield line[21], number, line[26], line[17:20].strip().upper(), line[16]


def chain_ids(pdb_text: str) -> List[str]:
    """Return the chain identifiers that have observed residues, in order of appearance."""
    seen: List[str] = []
    if structure_format(pdb_text) == "cif":
        for record in _cif_atoms(pdb_text):
            if record[4] not in seen:
                seen.append(record[4])
        return seen
    for line in pdb_text.splitlines():
        if line.startswith(("ATOM  ", "HETATM")) and len(line) > 21:
            chain = line[21]
            if chain not in seen:
                seen.append(chain)
    return seen


def chain_residues(pdb_text: str, chain: Optional[str] = None) -> List[Residue]:
    """
    Return the observed residues of a chain, in file order.

    Reads one entry per residue from the alpha-carbon records, so the result contains exactly
    the residues a viewer can color — residues missing from the crystal are absent here too,
    and that is the point: they cannot be highlighted and must not be counted as mapped.

    Args:
        pdb_text: contents of a structure file, PDB or mmCIF.
        chain: chain identifier; the first chain found when omitted.

    Returns:
        The residues, each carrying its number, insertion code and one-letter code.
    """
    if chain is None:
        ids = chain_ids(pdb_text)
        chain = ids[0] if ids else "A"
    return residues_by_chain(pdb_text).get(chain, [])


def chain_sequence(pdb_text: str, chain: Optional[str] = None) -> str:
    """The observed sequence of a chain, as a plain string."""
    return "".join(r.letter for r in chain_residues(pdb_text, chain))


# --------------------------------------------------------------------------- alignment
def _aligner():
    """A semi-global aligner: internal gaps are penalised, overhanging ends are free."""
    try:
        from Bio import Align
        from Bio.Align import substitution_matrices
    except ImportError as exc:  # pragma: no cover - biopython is a declared dependency
        raise RuntimeError(
            "O alinhamento estrutura-sequencia precisa do biopython "
            "(ja listado em requirements.txt): pip install biopython"
        ) from exc

    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    try:
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    except Exception:  # pragma: no cover - fall back to plain identity scoring
        aligner.match_score, aligner.mismatch_score = 1.0, -1.0
    aligner.open_gap_score = -11.0
    aligner.extend_gap_score = -1.0
    # A structure usually covers part of a longer sequence; penalising the overhang would push
    # the alignment to spread the structure across the whole protein instead of placing it.
    # Biopython renamed these (target_end_gap_score -> end_insertion_score); the old names still
    # work but warn. Probing with hasattr is not an option: reading a gap score raises
    # ValueError("gap scores are different") whenever the two sides differ, which they do here.
    # Writing is safe, and an unknown attribute raises AttributeError, so try new then old.
    for new_name, old_name in (("end_insertion_score", "target_end_gap_score"),
                               ("end_deletion_score", "query_end_gap_score")):
        try:
            setattr(aligner, new_name, 0.0)
        except AttributeError:
            setattr(aligner, old_name, 0.0)
    return aligner


def align_sequences(ref_seq: str, struct_seq: str) -> Tuple[Dict[int, int], float, int]:
    """
    Align a reference sequence to a structure's observed sequence.

    Args:
        ref_seq: the sequence the predictions were run on.
        struct_seq: the observed sequence of the structure's chain.

    Returns:
        ``(index_map, identity, aligned_columns)`` where ``index_map`` sends a 0-based index in
        ``ref_seq`` to a 0-based index in ``struct_seq``, ``identity`` is the fraction of aligned
        columns holding the same residue, and ``aligned_columns`` is how many columns had a
        residue on both sides.
    """
    if not ref_seq or not struct_seq:
        return {}, 0.0, 0

    alignment = _aligner().align(ref_seq, struct_seq)[0]
    index_map: Dict[int, int] = {}
    matches = columns = 0

    # ``aligned`` gives the ungapped blocks as ((ref_start, ref_end), (struct_start, struct_end)).
    for (ref_start, ref_end), (st_start, st_end) in zip(*alignment.aligned):
        for offset in range(ref_end - ref_start):
            ref_i, st_i = ref_start + offset, st_start + offset
            index_map[ref_i] = st_i
            columns += 1
            if ref_seq[ref_i] == struct_seq[st_i]:
                matches += 1

    identity = matches / columns if columns else 0.0
    return index_map, identity, columns


# --------------------------------------------------------------------------- mapping
def map_epitopes(
    pdb_text: str,
    ref_seq: str,
    ranges: Sequence[Tuple[int, int]],
    *,
    chain: Optional[str] = None,
    min_identity: float = DEFAULT_MIN_IDENTITY,
    sifts_offset: Optional[int] = None,
    sifts_span: Optional[Tuple[int, int]] = None,
) -> MappingReport:
    """
    Translate epitope positions into the residue numbers of a structure.

    Args:
        pdb_text: contents of the structure file.
        ref_seq: the sequence the predictions were run on. Its first residue is position 1,
            matching the ``Initial Position``/``Final Position`` columns POA1 emits.
        ranges: epitope ``(start, end)`` pairs, 1-based and inclusive, in ``ref_seq`` coordinates.
        chain: which chain to map onto; the first observed chain when omitted.
        min_identity: below this fraction of identical residues the structure is judged to be a
            different protein and the mapping is refused.
        sifts_offset: add this to a reference position to get the PDB residue number, when a
            SIFTS mapping is available. Exact, and skips the alignment.
        sifts_span: the ``(start, end)`` of ``ref_seq`` that the SIFTS mapping covers.

    Returns:
        A :class:`MappingReport`. When ``ok`` is false the caller must not highlight anything:
        the ranges would land on residues of an unrelated protein.
    """
    residues = chain_residues(pdb_text, chain)
    available = chain_ids(pdb_text)
    resolved_chain = chain or (available or ["A"])[0]
    return _map_onto(residues, resolved_chain, available, ref_seq, ranges,
                     min_identity=min_identity, sifts_offset=sifts_offset, sifts_span=sifts_span)


def _map_onto(residues, chain, available, ref_seq, ranges, *, min_identity,
              sifts_offset=None, sifts_span=None) -> MappingReport:
    """Map onto residues that were already read: the part of :func:`map_epitopes` after parsing."""
    report = MappingReport(chain=chain, method="alignment", identity=0.0,
                           structure_length=len(residues))

    if not residues:
        report.ok = False
        report.warnings.append(
            f"A cadeia '{chain}' nao tem residuos de proteina no arquivo. "
            f"Cadeias disponiveis: {', '.join(available) or 'nenhuma'}.")
        return report

    struct_seq = "".join(r.letter for r in residues)

    if sifts_offset is not None:
        _map_by_sifts(report, residues, ref_seq, struct_seq, ranges, sifts_offset, sifts_span,
                      min_identity)
    else:
        _map_by_alignment(report, residues, ref_seq, struct_seq, ranges, min_identity)

    logger.debug("structure mapping: %s", report.summary())
    return report


def residues_by_chain(pdb_text: str) -> Dict[str, List[Residue]]:
    """
    The observed protein residues of every chain, read in one pass over the file.

    Asking chain by chain would read the file once per chain: a 2.4 MB entry with a dozen chains
    would be split into lines a dozen times on every rerun of the interface. Chains with no
    protein residue (waters, ligands, nucleic acids) are absent from the result.
    """
    found: Dict[str, List[Residue]] = {}
    seen: set = set()
    for chain, number, icode, resname, altloc in _ca_atoms(pdb_text):
        if altloc not in (" ", "A"):
            continue
        letter = THREE_TO_ONE.get(resname)
        if letter is None:
            continue                      # a ligand or a nucleotide, not part of the chain
        if (chain, number, icode) in seen:
            continue
        seen.add((chain, number, icode))
        found.setdefault(chain, []).append(Residue(number=number, icode=icode, letter=letter))
    return found


def map_epitopes_best_chain(
    pdb_text: str,
    ref_seq: str,
    ranges: Sequence[Tuple[int, int]],
    *,
    min_identity: float = DEFAULT_MIN_IDENTITY,
) -> MappingReport:
    """
    Map onto whichever chain of the file *is* the protein, instead of assuming it is the first.

    :func:`map_epitopes` takes the first chain when none is named. That holds for a structure of
    the antigen alone and fails for the entries a person is most likely to pick by hand: in a
    complex the antibody chains often come first, and a file can open with a nucleic acid or a
    peptide. The first chain is then compared, found to be another protein, and the structure is
    refused as incompatible while the right chain sits further down the same file.

    Every distinct chain is tried and the best accepted one is used — the most identical residues
    over the aligned stretch, so a full-length copy beats a fragment. Chains with the same
    observed sequence are compared once.

    Returns:
        The report of the chosen chain. When no chain is accepted, the report of the closest one,
        with the first warning saying how many were tried — so the refusal names the best
        candidate rather than whichever chain happened to be first.
    """
    by_chain = residues_by_chain(pdb_text)
    available = chain_ids(pdb_text)
    if not by_chain:
        return map_epitopes(pdb_text, ref_seq, ranges, min_identity=min_identity)

    def matches(report: MappingReport) -> float:
        return report.identity * report.aligned_columns

    tried: Dict[str, MappingReport] = {}          # observed sequence -> its report
    for chain, residues in by_chain.items():
        sequence = "".join(r.letter for r in residues)
        if sequence not in tried:
            tried[sequence] = _map_onto(residues, chain, available, ref_seq, ranges,
                                        min_identity=min_identity)

    accepted = [report for report in tried.values() if report.ok]
    if accepted:
        return max(accepted, key=matches)         # max() keeps the first of equals: file order
    closest = max(tried.values(), key=matches)
    if len(tried) > 1 and closest.warnings:
        closest.warnings[0] += (f" Nenhuma das {len(tried)} cadeias distintas do arquivo "
                                f"corresponde; a cadeia {closest.chain} e a mais proxima.")
    return closest


def chain_copies(pdb_text: str, chain: str, min_shared: float = 0.9) -> List[str]:
    """
    The other chains that are copies of ``chain`` numbered the same way.

    A crystal of a dimer has the protein twice, and coloring the epitopes on one copy only draws
    a picture in which half of the antigen appears to have none. The residue numbers found for
    ``chain`` can be reused on another chain exactly when that chain calls the same residues by
    the same numbers, which is what this checks: every residue number the two chains share must
    be the same amino acid, and they must share most of their residues — copies routinely differ
    by a disordered loop, so identity of the observed lists would be too strict.

    Args:
        pdb_text: contents of the structure file.
        chain: the chain the epitopes were mapped onto.
        min_shared: fraction of the shorter chain's residues that must be present in both.

    Returns:
        Chain identifiers, in file order, excluding ``chain`` itself.
    """
    by_chain = residues_by_chain(pdb_text)
    target = {(r.number, r.icode): r.letter for r in by_chain.get(chain, [])}
    if not target:
        return []
    copies = []
    for other, residues in by_chain.items():
        if other == chain:
            continue
        theirs = {(r.number, r.icode): r.letter for r in residues}
        shared = target.keys() & theirs.keys()
        if (len(shared) >= min_shared * min(len(target), len(theirs))
                and all(target[key] == theirs[key] for key in shared)):
            copies.append(other)
    return copies


def _map_by_sifts(report, residues, ref_seq, struct_seq, ranges, offset, span,
                  min_identity=DEFAULT_MIN_IDENTITY):
    """Exact path: SIFTS already states the correspondence, so only verify and shift."""
    report.method = "sifts"
    by_number = {r.number: r for r in residues}
    lo, hi = span if span else (1, len(ref_seq))

    matches = columns = 0
    for position in range(lo, hi + 1):
        residue = by_number.get(position + offset)
        if residue is None or position > len(ref_seq):
            continue
        columns += 1
        if ref_seq[position - 1] == residue.letter:
            matches += 1
    report.identity = matches / columns if columns else 0.0
    report.aligned_columns = columns

    if report.identity < min_identity:
        # SIFTS said these correspond, yet the residues disagree: the reference sequence is not
        # the one the structure was solved from (a different isolate, or the wrong region).
        report.ok = False
        report.warnings.append(
            f"O mapeamento SIFTS aponta deslocamento {offset:+d}, mas so {report.identity_pct:.1f}% "
            f"dos residuos coincidem. A sequencia usada nas predicoes nao e a da estrutura.")
        return

    for start, end in ranges:
        if start < lo or end > hi:
            report.unmapped.append((start, end))
            continue
        present = sum(1 for position in range(start, end + 1)
                      if position + offset in by_number)
        if not present:
            # SIFTS says where these residues *would* be numbered; the file has none of them. An
            # epitope outside the stretch the crystal resolves was being reported as mapped and
            # merely "partial" — a range of residue numbers with nothing to color, counted among
            # the highlighted ones.
            report.unmapped.append((start, end))
            continue
        report.ranges.append((start + offset, end + offset))
        report.sources.append((start, end))
        if present < end - start + 1:
            report.partial.append((start, end))

    if report.unmapped:
        report.warnings.append(
            f"{len(report.unmapped)} epitopo(s) caem fora do trecho que a estrutura cobre "
            f"({lo}-{hi} da sequencia de referencia).")
    if report.partial:
        report.warnings.append(
            f"{len(report.partial)} epitopo(s) estao parcialmente ausentes no cristal: "
            f"a estrutura mostra menos residuos do que o epitopo tem.")


def _map_by_alignment(report, residues, ref_seq, struct_seq, ranges, min_identity,
                      min_coverage=DEFAULT_MIN_COVERAGE):
    """Fallback path: no SIFTS mapping, so derive the correspondence from the sequences."""
    index_map, identity, columns = align_sequences(ref_seq, struct_seq)
    report.identity = identity
    report.aligned_columns = columns
    coverage = columns / len(struct_seq) if struct_seq else 0.0

    if identity < min_identity or not columns:
        report.ok = False
        report.warnings.append(
            f"A estrutura tem {identity * 100:.1f}% de identidade com a sequencia usada nas "
            f"predicoes (minimo {min_identity * 100:.0f}%). Provavelmente e outra proteina: "
            f"destacar os epitopos produziria uma figura sem significado.")
        return

    # Identity passed, but over how much? The ratio is taken over the aligned columns alone, so a
    # handful of them can read 100% while the rest of the chain was discarded as end gap. A
    # structure of the right protein aligns over most of its length; one that only overlaps in a
    # few residues is not evidence, whatever the ratio says.
    if coverage < min_coverage:
        report.ok = False
        report.warnings.append(
            f"Apenas {columns} dos {len(struct_seq)} residuos da estrutura "
            f"({coverage * 100:.1f}%) alinharam com a sequencia usada nas predicoes, abaixo do "
            f"minimo de {min_coverage * 100:.0f}%. A identidade de {identity * 100:.1f}% e "
            f"calculada so sobre esses residuos, entao nao diz que a estrutura e desta proteina: "
            f"destacar os epitopos produziria uma figura sem significado.")
        return

    for start, end in ranges:
        # Only positions that actually aligned are used. Where an end is missing the range is
        # clipped to the observed part rather than extrapolated: without SIFTS there is no way
        # to know the residue number of a residue the file does not contain, and numbering with
        # gaps or insertion codes makes "last observed + 1" a guess that can land anywhere.
        observed = [index_map[position - 1] for position in range(start, end + 1)
                    if position - 1 in index_map]
        if not observed:
            report.unmapped.append((start, end))
            continue
        report.ranges.append((residues[observed[0]].number, residues[observed[-1]].number))
        report.sources.append((start, end))
        # Residues disordered in the crystal are absent from the file, so an epitope can be
        # placed correctly and still be only partly visible. The viewer draws what exists; the
        # report says how much that is, which the picture alone cannot.
        if len(observed) < end - start + 1:
            report.partial.append((start, end))

    if report.unmapped:
        report.warnings.append(
            f"{len(report.unmapped)} epitopo(s) nao tem correspondencia na estrutura "
            f"(fora do trecho resolvido ou com as extremidades ausentes).")
    if report.partial:
        report.warnings.append(
            f"{len(report.partial)} epitopo(s) estao parcialmente ausentes no cristal: "
            f"a estrutura mostra menos residuos do que o epitopo tem.")
    if identity < 0.95:
        report.warnings.append(
            f"Identidade de {identity * 100:.1f}%: a estrutura e de uma variante da proteina, "
            f"nao da sequencia exata usada nas predicoes.")
