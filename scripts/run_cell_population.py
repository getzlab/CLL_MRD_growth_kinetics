#!/usr/bin/env python3
"""Utility to fetch Cell Population inputs from Terra and run the solver."""

import argparse
import subprocess
import sys
from pathlib import Path
from typing import Dict, Iterable, Optional

from dalmatian.wmanager import WorkspaceManager

FILE_FIELD_CANDIDATES: Dict[str, Iterable[str]] = {
    "cluster_ccfs": (
        "cluster_ccfs",
        "cluster_ccfs_Aug2",
        "cluster_ccfs_tsv",
    ),
    "mut_ccfs": (
        "mut_ccfs",
        "mut_ccfs_Aug2",
        "mut_ccfs_tsv",
    ),
    "sif": (
        "sif_out",
        "sif_file",
    ),
    "build_tree_posteriors": (
        "build_tree_posteriors",
        "tree_tsv",
        "tree_tsv_Aug2",
    ),
}

DEFAULT_FILENAMES: Dict[str, str] = {
    "cluster_ccfs": "{patient_id}.cluster_ccfs.txt",
    "mut_ccfs": "{patient_id}.mut_ccfs.txt",
    "sif": "{patient_id}.sif",
    "build_tree_posteriors": "{patient_id}_build_tree_posteriors.tsv",
}


def pick_participant_value(row, candidates: Iterable[str], kind: str) -> str:
    for key in candidates:
        if key in row and isinstance(row[key], str) and row[key]:
            return row[key]
    raise KeyError(
        f"Unable to locate a value for '{kind}'. Tried columns: {', '.join(candidates)}"
    )


def download_gcs(uri: str, dest: Path, *, dry_run: bool = False) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dry_run:
        print(f"[DRY-RUN] gsutil cp {uri} {dest}")
        return
    cmd = ["gsutil", "cp", uri, str(dest)]
    subprocess.run(cmd, check=True)


def run_cell_population(
    patient_id: str,
    sif: Path,
    mut_ccf: Path,
    cluster_ccf: Path,
    tree_tsv: Path,
    tree_number: int,
    *,
    dry_run: bool = False,
) -> None:
    root = Path(__file__).resolve().parents[1]
    cell_pop_dir = root / "Cell_Population"
    cli = [
        sys.executable,
        str(cell_pop_dir / "PhylogicNDT.py"),
        "CellPopulation",
        "-i",
        patient_id,
        "-sif",
        str(sif.resolve()),
        "-m",
        str(mut_ccf.resolve()),
        "-c",
        str(cluster_ccf.resolve()),
        "-t",
        str(tree_tsv.resolve()),
        "--tree_number",
        str(tree_number),
    ]
    if dry_run:
        print("[DRY-RUN]", " ".join(cli))
        return
    subprocess.run(cli, check=True, cwd=cell_pop_dir)


def main() -> None:
    parser = argparse.ArgumentParser(description="Fetch inputs and run Cell Population analysis")
    parser.add_argument("patient_id", help="Participant identifier (e.g. CLL14-1063)")
    parser.add_argument("workspace", help="Terra workspace in the form namespace/name")
    parser.add_argument(
        "--tree-number",
        type=int,
        required=True,
        help="Tree number to pass to the CellPopulation module",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("Cell_Population"),
        help="Directory where input files should be written (default: Cell_Population)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the actions that would be taken without downloading or running",
    )
    args = parser.parse_args()

    output_dir: Path = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    wm = WorkspaceManager(args.workspace)
    participants = wm.get_participants()
    if args.patient_id not in participants.index:
        raise SystemExit(f"Patient {args.patient_id} not found in workspace {args.workspace}")
    participant_row = participants.loc[args.patient_id]

    uris: Dict[str, str] = {}
    for kind, candidates in FILE_FIELD_CANDIDATES.items():
        try:
            uris[kind] = pick_participant_value(participant_row, candidates, kind)
        except KeyError as err:
            raise SystemExit(str(err)) from err

    local_paths: Dict[str, Path] = {}
    for kind, uri in uris.items():
        filename = DEFAULT_FILENAMES[kind].format(patient_id=args.patient_id)
        dest = output_dir / filename
        download_gcs(uri, dest, dry_run=args.dry_run)
        local_paths[kind] = dest

    run_cell_population(
        args.patient_id,
        sif=local_paths["sif"],
        mut_ccf=local_paths["mut_ccfs"],
        cluster_ccf=local_paths["cluster_ccfs"],
        tree_tsv=local_paths["build_tree_posteriors"],
        tree_number=args.tree_number,
        dry_run=args.dry_run,
    )

    if not args.dry_run:
        print(
            f"Completed CellPopulation for {args.patient_id}. "
            f"Inputs downloaded to {output_dir.resolve()}"
        )


if __name__ == "__main__":
    main()
