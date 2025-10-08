#!/usr/bin/env python3
import argparse
import json
import random
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
    treatment_file: Path
    sample_ids: Sequence[str]
    output_html: Path

    umi_start: int
    umi_end: int
    treatment_start: int
    treatment_end: int
    tree_choice: int

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
    treatment_df = pd.read_csv(cfg.treatment_file, sep="\t")

    if cfg.workspace:
        import dalmatian  # optional
        wm = dalmatian.WorkspaceManager(cfg.workspace)
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


def build_patient_report(cfg: PatientConfig) -> None:
    wbc_df, treatment_df, cluster_ccf_df, abundance_df, mcmc_df, tree_df = load_patient_data(cfg)
    wbc_patient_df, times_sample, cll_counts, wbc_counts, all_times = filter_patient_data(wbc_df, cfg.patient_id)

    metadata_table = plotly_helper.plot_metadata_table(wbc_patient_df, cfg.patient_id)

    ccf_plot = plotly_helper.plot_ccf(cluster_ccf_df, times_sample, treatment_df)
    combined_ccf_tree = plotly_helper.plot_ccf_tree_combined(
        tree_df=tree_df,
        tree_selected=cfg.tree_choice,
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
    )
    cluster_list, cluster_abundance = model_helper.get_abundance(abundance_df, mcmc_df, cfg.sample_ids)
    _, log_subclone = model_helper.calc_subclone(cll_counts, cluster_abundance, cluster_list)

    mcmc_abundance = model_helper.get_all_abundance(cluster_list, mcmc_df, cfg.sample_ids, times_sample)
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
    times_sliced_after = [
        int(value)
        for value in wbc_patient_df.loc[
            wbc_patient_df["CLL count estm"] > 0, "Time_since_start_tx"
        ]
        if int(value) > 0
    ]
    wbc_model = [
        float(val)
        for val in wbc_patient_df.loc[
            wbc_patient_df["Time_since_start_tx"] > cfg.treatment_end,
            "CLL count estm",
        ]
        if not pd.isna(val) and val > 0
    ]

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

    plotly_helper.create_html_file(
        plots,
        output_file=str(cfg.output_html),
    )


def load_config(path: Path) -> Iterable[PatientConfig]:
    data = json.loads(path.read_text())
    if isinstance(data, dict):
        data = [data]
    for entry in data:
        entry["wbc_file"] = Path(entry["wbc_file"])
        entry["treatment_file"] = Path(entry["treatment_file"])
        entry["output_html"] = Path(entry["output_html"])
        entry["sample_ids"] = entry.get("sample_ids") or []
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
