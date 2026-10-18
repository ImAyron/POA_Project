"""CLI: run the pipeline directly on the original research folders.

Adapts every ``bepipred_*.json`` in a predictions folder to the pipeline's conventions
(:mod:`poa.services.realdata_import`) and, with ``--run``, executes POA1 + local Conservancy
for each virus, auto-mapping it to its world/diversity FASTA in ``--world-dir``.

Example:
    python -m poa.cli.import_realdata \
        --pred-dir  ".../5_predicao_de_epitopo" \
        --world-dir ".../2_fastas_mundo" \
        --out       "./realdata_run" \
        --protein E --chikv-protein E1 --threshold 70 --run
"""
from __future__ import annotations

import argparse
import os
from types import SimpleNamespace

from ..core.pipeline import run_poa1
from ..logging_conf import setup_logging
from ..services import conservancy_client, realdata_import


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Adapt and run the pipeline on the original research folders.")
    p.add_argument("--pred-dir", required=True, help="Folder with bepipred_*.json (e.g. 5_predicao_de_epitopo)")
    p.add_argument("--world-dir", default="", help="Folder with world/diversity FASTAs for conservancy (e.g. 2_fastas_mundo)")
    p.add_argument("--out", required=True, help="Output root directory")
    p.add_argument("--protein", default="E", help="Default protein label for headers (default: E)")
    p.add_argument("--chikv-protein", default="E1", help="Protein label override for CHIKV (default: E1)")
    p.add_argument("--threshold", type=int, default=70, help="Conservancy sequence-identity threshold (default: 70)")
    p.add_argument("--run", action="store_true", help="Run POA1 + Conservancy (otherwise only prepare inputs)")
    return p


def main(argv=None):
    setup_logging()
    args = build_parser().parse_args(argv)

    prep_dir = os.path.join(args.out, "_prepared")
    prepared = realdata_import.prepare_bepipred2_dir(
        args.pred_dir, prep_dir, protein=args.protein,
        protein_map={"CHIKV": args.chikv_protein},
    )
    print(f"Prepared {len(prepared)} virus input(s) in {prep_dir}")
    for pr in prepared:
        flag = "  [X in reference — POA1 will refuse]" if pr.has_x else ""
        print(f"  {pr.specie:6s} protein={pr.protein:3s} ref_len={len(pr.sequence):4d}{flag}")

    if not args.run:
        print("\n(dry run — pass --run to execute POA1 + Conservancy)")
        return prepared

    if not args.world_dir:
        raise SystemExit("--run requires --world-dir for the conservancy step")

    print("\n=== POA1 + Conservancy ===")
    for pr in prepared:
        vdir = os.path.join(args.out, pr.specie)
        if pr.has_x:
            print(f"[{pr.specie}] SKIPPED — reference contains 'X'")
            continue

        poa1_args = SimpleNamespace(
            b2=pr.b2_json, b3="", bmin=0, bmax=0, p="", pmin=0, pmax=0,
            n="", m="", mhla="DR", mic=50, x="", xmin=0, xmax=0,
            d=vdir, f=pr.reference_fasta, e="y",
        )
        r1 = run_poa1(poa1_args)

        world_path, ident = (None, 0.0)
        if args.world_dir:
            world_path, ident = realdata_import.best_world_match(pr.sequence, args.world_dir)
        if not world_path:
            print(f"[{pr.specie}] POA1 ok ({len(r1.predictions)} epitopes) — no world match found, conservancy skipped")
            continue

        csv_dir = os.path.join(vdir, "csvs")
        written = conservancy_client.run_conservancy_for_dir(
            r1.conservancy_dir, world_path, threshold=args.threshold, out_dir=csv_dir,
        )
        print(f"[{pr.specie}] epitopes={len(r1.predictions):3d}  world={os.path.basename(world_path)} "
              f"(id~{ident:.0f}%)  csvs={len(written)}")

    print("\nDONE. POA2 (TMHMM topology) must run where pyTMHMM is installed (see TESTING_WSL.md).")
    return prepared


if __name__ == "__main__":
    main()
