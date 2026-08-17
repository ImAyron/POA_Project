"""UI-agnostic orchestration for the GUI: working dirs, argument building, stage runners.

Kept free of Streamlit imports so it can be unit-tested and reused. The Streamlit app
(:mod:`poa.gui.app`) is a thin layer that calls into this module.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace
from typing import Dict, List, Optional, Sequence, Tuple

import pandas as pd

from ..core import pipeline
from ..services import realdata_import
from ..services.cache import Cache


@dataclass
class WorkContext:
    """Per-session working directory tree."""

    root: Path

    def __post_init__(self):
        self.root = Path(self.root)
        for d in (self.inputs_dir, self.mhcii_dir, self.poa1_out_dir,
                  self.conservancy_csv_dir, self.poa2_out_dir, self.cache_dir):
            d.mkdir(parents=True, exist_ok=True)

    @property
    def inputs_dir(self) -> Path:
        return self.root / "inputs"

    @property
    def mhcii_dir(self) -> Path:
        return self.root / "mhcii"

    @property
    def poa1_out_dir(self) -> Path:
        return self.root / "poa1_out"

    @property
    def conservancy_epitopes_dir(self) -> Path:
        # where POA1 writes the per-species epitope FASTAs
        return self.poa1_out_dir / "Conservancy Analysis"

    @property
    def conservancy_csv_dir(self) -> Path:
        return self.root / "conservancy_csv"

    @property
    def poa2_out_dir(self) -> Path:
        return self.root / "poa2_out"

    @property
    def cache_dir(self) -> Path:
        return self.root / "cache"

    def cache(self) -> Cache:
        return Cache(self.cache_dir)


# --------------------------------------------------------------------------- argument builders
def build_poa1_args(files: Dict[str, str], params: Dict, proteins_path: str, out_dir: str) -> SimpleNamespace:
    """
    Assemble a POA1 args namespace.

    Parameters:
        files: optional keys 'b2','b3','p','n','m','x' mapping to result file/dir paths.
        params: optional keys 'bmin','bmax','pmin','pmax','mhla','mic','xmin','xmax','e'.
        proteins_path: path to the protein FASTA (-f).
        out_dir: analysis output directory (-d).
    """
    return SimpleNamespace(
        b2=files.get("b2", ""),
        b3=files.get("b3", ""),
        bmin=int(params.get("bmin", 0)),
        bmax=int(params.get("bmax", 0)),
        p=files.get("p", ""),
        pmin=int(params.get("pmin", 0)),
        pmax=int(params.get("pmax", 0)),
        n=files.get("n", ""),
        m=files.get("m", ""),
        mhla=str(params.get("mhla", "DR")),
        mic=int(params.get("mic", 50)),
        x=files.get("x", ""),
        xmin=int(params.get("xmin", 0)),
        xmax=int(params.get("xmax", 0)),
        d=out_dir,
        f=proteins_path,
        e=str(params.get("e", "n")),
    )


def build_poa2_args(params: Dict, csv_dir: str, proteins_path: str, out_dir: str) -> SimpleNamespace:
    """
    Assemble a POA2 args namespace.

    Parameters:
        params: keys 'objective' ('conserved'|'unique'), 't', 'imin','imax','m','rf'.
    """
    conserved = params.get("objective", "conserved") == "conserved"
    return SimpleNamespace(
        g="True" if conserved else None,
        l=None if conserved else "True",
        t=int(params["t"]),
        d=csv_dir,
        r=out_dir,
        imin=int(params.get("imin", 60)),
        imax=int(params.get("imax", 100)),
        m=int(params.get("m", 60)),
        f=proteins_path,
        rf=(None if params.get("rf", None) in (None, "") else int(params["rf"])),
    )


# --------------------------------------------------------------------------- stage runners
def run_poa1(ctx: WorkContext, files: Dict[str, str], params: Dict, proteins_path: str) -> pipeline.Poa1Result:
    args = build_poa1_args(files, params, proteins_path, str(ctx.poa1_out_dir))
    return pipeline.run_poa1(args)


def run_poa2(ctx: WorkContext, params: Dict, proteins_path: str) -> pipeline.Poa2Result:
    args = build_poa2_args(params, str(ctx.conservancy_csv_dir), proteins_path, str(ctx.poa2_out_dir))
    return pipeline.run_poa2(args)


# --------------------------------------------------------------------------- summaries for charts
def predictions_summary(df: pd.DataFrame) -> Dict[str, pd.DataFrame]:
    """Return small aggregate tables (by method, by species) for charting the POA1 output."""
    out: Dict[str, pd.DataFrame] = {}
    if df is None or df.empty:
        return {"by_method": pd.DataFrame(columns=["Method", "count"]),
                "by_species": pd.DataFrame(columns=["Specie", "count"]),
                "length_hist": pd.DataFrame(columns=["length", "count"])}
    out["by_method"] = df.groupby("Method").size().reset_index(name="count")
    out["by_species"] = df.groupby("Specie").size().reset_index(name="count")
    lengths = df["Peptide Sequence"].astype(str).str.len()
    out["length_hist"] = lengths.value_counts().sort_index().reset_index()
    out["length_hist"].columns = ["length", "count"]
    return out


def topology_summary(df: pd.DataFrame) -> pd.DataFrame:
    """Aggregate POA2 topology fractions (mean outer/TM/inner) for charting."""
    cols = ["Portion_Outside", "Portion_TM", "Portion_Inside"]
    if df is None or df.empty or not all(c in df.columns for c in cols):
        return pd.DataFrame(columns=["region", "mean_fraction"])
    means = {
        "Outside": pd.to_numeric(df["Portion_Outside"], errors="coerce").mean(),
        "Transmembrane": pd.to_numeric(df["Portion_TM"], errors="coerce").mean(),
        "Inside": pd.to_numeric(df["Portion_Inside"], errors="coerce").mean(),
    }
    return pd.DataFrame({"region": list(means.keys()), "mean_fraction": list(means.values())})


def save_upload(dest_dir: Path, filename: str, data: bytes) -> str:
    """Persist an uploaded file's bytes to ``dest_dir/filename`` and return the path."""
    dest_dir.mkdir(parents=True, exist_ok=True)
    path = dest_dir / filename
    path.write_bytes(data)
    return str(path)


# --------------------------------------------------------------------------- real-data adapter
def import_bepipred2(ctx: WorkContext, filename: str, data: bytes,
                     specie: str, protein: str = "E") -> realdata_import.Prepared:
    """
    Adapt an original-workflow BepiPred-2.0 JSON upload to pipeline conventions.

    Persists the upload, then rewrites its antigen key to ``Protein_Specie_ID`` and rebuilds the
    reference protein FASTA from the JSON's own ``AA`` array (see :mod:`poa.services.realdata_import`).
    Returns the :class:`~poa.services.realdata_import.Prepared` record (paths + ``has_x`` flag);
    the caller wires ``b2_json`` into the POA1 files and may adopt ``reference_fasta`` as ``-f``.
    """
    raw_path = save_upload(ctx.inputs_dir, filename, data)
    return realdata_import.prepare_bepipred2(raw_path, specie, str(ctx.inputs_dir), protein=protein)
