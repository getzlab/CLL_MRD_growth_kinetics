#!/usr/bin/env python3
import argparse
import json
import random
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

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

    cll_plot = plotly_helper.plot_CLL_count(
        cfg.patient_id,
        times_sample,
        cll_counts,
        cfg.umi_start,
        cfg.umi_end,
        cfg.treatment_start,
        cfg.treatment_end,
    )
    
    metadata_table = plotly_helper.plot_metadata_table(wbc_patient_df, cfg.patient_id)
    tree_plot = plotly_helper.plot_tree_plotly(tree_df, cfg.tree_choice)
    ccf_plot = plotly_helper.plot_ccf(cluster_ccf_df, times_sample, treatment_df)
    combined_ccf_tree = plotly_helper.plot_ccf_tree_combined(
        tree_df=tree_df,
        tree_selected=cfg.tree_choice,
        ccf_df=cluster_ccf_df,
        times_sample=times_sample,
        treatment_df=treatment_df,
        branch_annotations=cfg.branch_annotations,
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

    plotly_helper.create_html_file(
        [
            cll_plot,
            metadata_table,
            tree_plot,
            ccf_plot,
            combined_ccf_tree,
            subclone_plot,
            linear_model_plot,
        ],
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
        yield PatientConfig(**entry)


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate clonokinetics patient report.")
    parser.add_argument("--config", type=Path, required=True, help="JSON configuration file (single patient or list).")
    args = parser.parse_args()

    for patient_cfg in load_config(args.config):
        build_patient_report(patient_cfg)


if __name__ == "__main__":
    main()