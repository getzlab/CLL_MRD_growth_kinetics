#!/usr/bin/env python3
import argparse
import json
import random
from datetime import datetime
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from clonokinetix import model_helper, plotly_helper


@dataclass
class PatientConfig:
    patient_id: str
    wbc_file: Path
    sample_ids: Sequence[str]
    output_html: Path

    umi_start: int
    umi_end: int
    treatment_start: int
    treatment_end: int
    tree_choice: int

    treatment_file: Optional[Path] = None
    branch_annotations: Dict[int, str] = field(default_factory=dict)
    fludarabine_windows: Sequence[Tuple[int, int]] = field(default_factory=tuple)

    # Data locations; fall back to Cell_Population layout if omitted
    cluster_ccf_path: Optional[Path] = None
    tree_posteriors_path: Optional[Path] = None
    abundance_path: Optional[Path] = None
    mcmc_path: Optional[Path] = None

    # Terra support (falls back to local files if not set)
    workspace: Optional[str] = None

    # noise injection options
    noise_low: float = 0.0
    noise_high: float = 0.01
    noise_seed: int = 42
    model_iteration: int = 2
    index_samples_model: Optional[Tuple[int, Optional[int]]] = None
    fi_modeling: bool = False
    

def _default_cell_population_path(patient_id: str, suffix: str) -> Path:
    # Files with dots: cluster_ccfs, mut_ccfs, sif
    # Files with underscores: build_tree_posteriors, cell_population_abundances, etc.
    if suffix.startswith("cluster_ccfs") or suffix.startswith("mut_ccfs") or suffix == "sif":
        return Path("Cell_Population") / f"{patient_id}.{suffix}"
    else:
        return Path("Cell_Population") / f"{patient_id}_{suffix}"


def add_uniform_noise(abundance: Dict[int, Dict[int, List[float]]],
                      low: float,
                      high: float,
                      seed: int) -> Dict[int, Dict[int, List[float]]]:
    rng = random.Random(seed)
    noisy = {}
    for cluster, iter_map in abundance.items():
        noisy[cluster] = {}
        for iteration, values in iter_map.items():
            perturbed = [x + rng.uniform(low, high) for x in values]
            total = sum(perturbed) or 1.0
            noisy[cluster][iteration] = [x / total for x in perturbed]
    return noisy


def transpose_iterations(noisy: Dict[int, Dict[int, List[float]]]) -> Dict[int, Dict[int, List[float]]]:
    per_iter: Dict[int, Dict[int, List[float]]] = {}
    for cluster, iter_map in noisy.items():
        for iteration, values in iter_map.items():
            per_iter.setdefault(iteration, {})[cluster] = values
    # transpose back so each cluster -> iteration -> list stays aligned
    recomposed: Dict[int, Dict[int, List[float]]] = {}
    for cluster in noisy.keys():
        recomposed[cluster] = {iteration: per_iter[iteration][cluster] for iteration in per_iter}
    return recomposed


def _validate_treatment_window(cfg: PatientConfig, wbc_df: pd.DataFrame) -> None:
    patient_rows = wbc_df[wbc_df["Patient"] == cfg.patient_id]
    if patient_rows.empty:
        raise ValueError(f"No rows found for patient {cfg.patient_id} in WBC file {cfg.wbc_file}")

    start_values = patient_rows["Start of treatment"].dropna().astype(str).str.strip()
    end_values = patient_rows["End of treatment"].dropna().astype(str).str.strip()

    if start_values.empty or end_values.empty:
        # If either date is missing we cannot validate; surface a clear error
        raise ValueError(
            f"Missing treatment start/end dates for patient {cfg.patient_id} in WBC file {cfg.wbc_file}"
        )

    start_date = pd.to_datetime(start_values.iloc[0], errors="coerce")
    end_date = pd.to_datetime(end_values.iloc[0], errors="coerce")
    if pd.isna(start_date) or pd.isna(end_date):
        raise ValueError(
            f"Unable to parse treatment dates for patient {cfg.patient_id}: "
            f"start='{start_values.iloc[0]}', end='{end_values.iloc[0]}'"
        )

    observed_duration = int((end_date - start_date).days)
    if observed_duration != cfg.treatment_end:
        raise ValueError(
            f"Configured treatment_end ({cfg.treatment_end}) for patient {cfg.patient_id} "
            f"does not match duration between Start of treatment ({start_date.date()}) and "
            f"End of treatment ({end_date.date()}) in '{cfg.wbc_file}' ({observed_duration} days)."
        )


def _validate_sample_counts(cfg: PatientConfig, patient_df: pd.DataFrame) -> None:
    """Ensure JSON-configured samples match the number of WBC samples."""
    if not cfg.sample_ids:
        return
    wbc_samples = (
        patient_df.loc[patient_df["Sample"].notna(), "Sample"]
        .astype(str)
        .tolist()
    )
    if len(cfg.sample_ids) != len(wbc_samples):
        raise ValueError(
            f"Sample count mismatch for patient {cfg.patient_id}: "
            f"{len(cfg.sample_ids)} identifiers in config vs {len(wbc_samples)} WBC samples."
        )


def _parse_tree_edges(edges_str: str) -> set:
    edges = set()
    for item in edges_str.split(","):
        item = item.strip()
        if not item:
            continue
        parent_str, child_str = item.split("-")
        parent = None if parent_str.strip() == "None" else int(parent_str.strip())
        child = None if child_str.strip() == "None" else int(child_str.strip())
        edges.add((parent, child))
    return edges


def _derive_edges_from_abundance(abundance_df: pd.DataFrame) -> set:
    edges = set()
    cell_pops = abundance_df["Cell_population"].dropna().unique()
    for label in cell_pops:
        parts = [part for part in label.split("_") if part.startswith("CL")]
        if not parts:
            continue
        cluster_ids = [int(part[2:]) for part in parts if part[2:].isdigit()]
        if not cluster_ids:
            continue
        edges.add((None, cluster_ids[0]))
        for parent, child in zip(cluster_ids, cluster_ids[1:]):
            edges.add((parent, child))
    return edges


def _validate_tree_choice(cfg: PatientConfig, tree_df: pd.DataFrame, abundance_df: pd.DataFrame) -> None:
    if tree_df.empty or abundance_df.empty:
        return
    tree_index = cfg.tree_choice - 1 if cfg.tree_choice > 0 else 0
    tree_index = max(0, min(tree_index, len(tree_df) - 1))
    selected_edges_str = tree_df.iloc[tree_index]["edges"]
    selected_edges = _parse_tree_edges(selected_edges_str)
    abundance_edges = _derive_edges_from_abundance(abundance_df)
    if selected_edges != abundance_edges:
        def _sort_key(edge):
            parent = -1 if edge[0] is None else edge[0]
            child = -1 if edge[1] is None else edge[1]
            return (parent, child)
        raise ValueError(
            f"Tree mismatch for patient {cfg.patient_id}: "
            f"config tree_choice={cfg.tree_choice} yields edges {sorted(selected_edges, key=_sort_key)} "
            f"but CellPopulation abundances imply {sorted(abundance_edges, key=_sort_key)}."
        )


def _resolve_index_slice(cfg: PatientConfig, times_sample: Sequence[int]) -> slice:
    if cfg.index_samples_model is not None:
        start, stop = cfg.index_samples_model
        return slice(start, stop)

    post_tx_indices = [i for i, t in enumerate(times_sample) if t > cfg.treatment_end]
    if post_tx_indices:
        return slice(post_tx_indices[0], post_tx_indices[-1] + 1)

    return slice(1, len(times_sample))


def load_patient_data(cfg: PatientConfig):
    wbc_df = pd.read_csv(cfg.wbc_file)
    if cfg.treatment_file:
        treatment_df = pd.read_csv(cfg.treatment_file, sep="\t")
    else:
        treatment_df = pd.DataFrame(
            [
                {
                    "tx": "Treatment",
                    "tx_start": cfg.treatment_start,
                    "tx_end": cfg.treatment_end,
                }
            ]
        )

    _validate_treatment_window(cfg, wbc_df)

    if cfg.workspace:
        from dalmatian.wmanager import WorkspaceManager  # optional
        wm = WorkspaceManager(cfg.workspace)
        participants = wm.get_participants()
        participant_row = participants.loc[cfg.patient_id]

        cluster_ccf_df = pd.read_csv(participant_row["cluster_ccfs"], sep="\t")
        # Optional tree from workspace; default to local fallback if not provided
        tree_df = pd.read_csv(
            cfg.tree_posteriors_path or _default_cell_population_path(cfg.patient_id, "build_tree_posteriors.tsv"),
            sep="\t",
        )
        abundance_df = pd.read_csv(
            cfg.abundance_path or _default_cell_population_path(cfg.patient_id, "cell_population_abundances.tsv"),
            sep="\t",
        )
        mcmc_df = pd.read_csv(
            cfg.mcmc_path or _default_cell_population_path(cfg.patient_id, "cell_population_mcmc_trace.tsv"),
            sep="\t",
        )
    else:
        cluster_ccf_df = pd.read_csv(
            cfg.cluster_ccf_path or _default_cell_population_path(cfg.patient_id, "cluster_ccfs.txt"),
            sep="\t",
        )
        tree_df = pd.read_csv(
            cfg.tree_posteriors_path or _default_cell_population_path(cfg.patient_id, "build_tree_posteriors.tsv"),
            sep="\t",
        )
        abundance_df = pd.read_csv(
            cfg.abundance_path or _default_cell_population_path(cfg.patient_id, "cell_population_abundances.tsv"),
            sep="\t",
        )
        mcmc_df = pd.read_csv(
            cfg.mcmc_path or _default_cell_population_path(cfg.patient_id, "cell_population_mcmc_trace.tsv"),
            sep="\t",
        )

    return wbc_df, treatment_df, cluster_ccf_df, abundance_df, mcmc_df, tree_df


def filter_patient_data(wbc_df: pd.DataFrame, patient_id: str):
    patient_df = wbc_df[wbc_df["Patient"] == patient_id].reset_index(drop=True)
    times_sample = [int(v) for v in patient_df.loc[patient_df.Sample.notna(), "Time_since_start_tx"]]
    wbc_all = [float(v) for v in patient_df["WBC"]]
    cll_count = [float(v) for v in patient_df.loc[patient_df.Sample.notna(), "CLL count estm"]]
    all_times = [int(v) for v in patient_df["Time_since_start_tx"]]
    return patient_df, times_sample, cll_count, wbc_all, all_times


def _extract_treatment_label(patient_df: pd.DataFrame) -> Optional[str]:
    if "Treatment" not in patient_df.columns:
        return None

    treatments = (
        patient_df["Treatment"]
        .dropna()
        .astype(str)
        .str.strip()
    )
    unique_non_empty = [value for value in treatments if value]
    if not unique_non_empty:
        return None

    # Preserve first occurrence order while collapsing duplicates
    return next(iter(dict.fromkeys(unique_non_empty)))


def build_patient_report(cfg: PatientConfig) -> None:
    wbc_df, treatment_df, cluster_ccf_df, abundance_df, mcmc_df, tree_df = load_patient_data(cfg)
    wbc_patient_df, times_sample, cll_counts, wbc_counts, all_times = filter_patient_data(wbc_df, cfg.patient_id)
    _validate_sample_counts(cfg, wbc_patient_df)
    _validate_tree_choice(cfg, tree_df, abundance_df)

    treatment_label = _extract_treatment_label(wbc_patient_df)
    if treatment_label:
        treatment_df = treatment_df.copy()
        treatment_df["tx"] = treatment_label

    metadata_table = plotly_helper.plot_metadata_table(wbc_patient_df, cfg.patient_id)

    # tree_choice is 1-based for user convenience; clamp to valid 0-based index
    tree_index = cfg.tree_choice - 1 if cfg.tree_choice > 0 else 0
    tree_index = max(0, min(tree_index, len(tree_df) - 1))

    ccf_plot = plotly_helper.plot_ccf(cluster_ccf_df, times_sample, treatment_df)
    combined_ccf_tree = plotly_helper.plot_ccf_tree_combined(
        tree_df=tree_df,
        tree_selected=tree_index,
        ccf_df=cluster_ccf_df,
        times_sample=times_sample,
        treatment_df=treatment_df,
        branch_annotations=cfg.branch_annotations,
    )
    cll_plot = plotly_helper.plot_CLL_count(
        cfg.patient_id,
        times_sample,
        cll_counts,
        cfg.umi_start,
        cfg.umi_end,
        cfg.treatment_start,
        cfg.treatment_end,
        treatment_label=treatment_label,
    )
    cluster_list, cluster_abundance = model_helper.get_abundance(abundance_df, mcmc_df, cfg.sample_ids)
    _, log_subclone = model_helper.calc_subclone(cll_counts, cluster_abundance, cluster_list)

    mcmc_abundance = model_helper.get_all_abundance(cluster_list, mcmc_df, cfg.sample_ids)
    noisy_mcmc_abundance = add_uniform_noise(mcmc_abundance, cfg.noise_low, cfg.noise_high, cfg.noise_seed)
    noisy_mcmc_abundance = transpose_iterations(noisy_mcmc_abundance)

    _, log_subclone_mcmc = model_helper.calc_subclone(
        cll_counts,
        noisy_mcmc_abundance,
        cluster_list,
        input_type="mcmc",
    )

    times_after_treatment = [t for t in all_times if t > cfg.treatment_end]
    times_after_treatment.insert(0, cfg.treatment_end)
    extrapolate_start_idx = 1

    subclone_plot = plotly_helper.plot_subclones(
        cluster_list,
        times_sample,
        cll_counts,
        log_subclone,
        extrapolate_start_idx,
        times_after_treatment,
        treatment_df,
        cfg.treatment_end,
    )
    linear_model_plot = plotly_helper.plot_linear_model_mcmc(
        cluster_list,
        times_sample,
        cll_counts,
        log_subclone_mcmc,
        extrapolate_start_idx,
        times_after_treatment,
        treatment_df,
        cfg.treatment_end,
    )

    wbc_patient_df["CLL count estm"] = pd.to_numeric(
        wbc_patient_df["CLL count estm"], errors="coerce"
    )
    post_treatment_mask = (
        (wbc_patient_df["CLL count estm"] > 0)
        & (wbc_patient_df["Time_since_start_tx"] > cfg.treatment_end)
    )
    times_sliced_after = (
        wbc_patient_df.loc[post_treatment_mask, "Time_since_start_tx"]
        .dropna()
        .astype(int)
        .tolist()
    )
    wbc_model = (
        wbc_patient_df.loc[post_treatment_mask, "CLL count estm"]
        .dropna()
        .astype(float)
        .tolist()
    )

    extra_plots: List[str] = []
    if wbc_model and times_sliced_after:
        index_slice = _resolve_index_slice(cfg, times_sample)
        start = max(0, index_slice.start or 0)
        stop = index_slice.stop if index_slice.stop is not None else len(times_sample)
        stop = max(start + 1, min(stop, len(times_sample)))
        index_slice = slice(start, stop)

        try:
            X, y = model_helper.create_inputs(
                times_sliced_after,
                log_subclone_mcmc,
                cfg.model_iteration,
                index_slice,
                times_sample,
            )
            logsumexp_points = np.log(wbc_model)
            multi_cluster_model = model_helper.MultiClusterLinearRegression(
                len(cluster_list),
                X,
                y,
            )
            multi_cluster_model.fit(logsumexp_points)
        except Exception as exc:  # pragma: no cover - defensive
            print(f"[WARN] skipping new-model plots for {cfg.patient_id}: {exc}")
        else:
            extra_plots.append(
                plotly_helper.plot_subclones_new_model(
                    cluster_list,
                    times_sample,
                    wbc_model,
                    log_subclone,
                    extrapolate_start_idx,
                    times_after_treatment,
                    times_sliced_after,
                    treatment_df,
                    cfg.treatment_end,
                    multi_cluster_model,
                    cfg.fi_modeling,
                )
            )
            extra_plots.append(
                plotly_helper.plot_mcmc_model(
                    cluster_list,
                    index_slice,
                    times_after_treatment,
                    times_sample,
                    times_sliced_after,
                    cfg.sample_ids,
                    wbc_model,
                    log_subclone_mcmc,
                    treatment_df,
                    cfg.treatment_end,
                    cfg.fi_modeling,
                )
            )

    plots = [
        metadata_table,
        ccf_plot,
        combined_ccf_tree,
        cll_plot,
        subclone_plot,
        linear_model_plot,
    ]
    plots.extend(extra_plots)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    requested_path = Path(cfg.output_html)
    if requested_path.suffix:
        output_dir = requested_path.parent
    else:
        output_dir = requested_path

    output_dir.mkdir(parents=True, exist_ok=True)
    timestamped_output = output_dir / f"{cfg.patient_id}_report_{timestamp}.html"

    plotly_helper.create_html_file(
        plots,
        output_file=str(timestamped_output),
    )


def load_config(path: Path) -> Iterable[PatientConfig]:
    data = json.loads(path.read_text())
    if isinstance(data, dict):
        data = [data]
    for entry in data:
        entry["wbc_file"] = Path(entry["wbc_file"])
        entry["output_html"] = Path(entry["output_html"])
        entry["sample_ids"] = entry.get("sample_ids") or []
        if entry.get("treatment_file"):
            entry["treatment_file"] = Path(entry["treatment_file"])
        else:
            entry["treatment_file"] = None
        for key in ("cluster_ccf_path", "tree_posteriors_path", "abundance_path", "mcmc_path"):
            if entry.get(key):
                entry[key] = Path(entry[key])
        entry["fludarabine_windows"] = [tuple(win) for win in entry.get("fludarabine_windows", [])]
        
        if entry.get("index_samples_model") is not None:
            values = entry["index_samples_model"]
            if isinstance(values, (list, tuple)):
                if len(values) == 1:
                    values = (values[0], None)
                entry["index_samples_model"] = tuple(values)
            else:
                raise ValueError("index_samples_model must be a list or tuple")
        yield PatientConfig(**entry)


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate clonokinetics patient report.")
    parser.add_argument("--config", type=Path, required=True, help="JSON configuration file (single patient or list).")
    args = parser.parse_args()

    for patient_cfg in load_config(args.config):
        build_patient_report(patient_cfg)


if __name__ == "__main__":
    main()
