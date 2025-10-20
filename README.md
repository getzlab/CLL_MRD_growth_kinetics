# Fludarabine MRD Analysis

This repository collects the exploratory and production material for the CLL MRD growth kinetics projects.

## Repository layout
- `notebooks/` – all Jupyter notebooks grouped by project (`crc`, `fi`, `gcll`, `misc`, and `autogen`).
  - `gcll/active/` holds the working set of GCLL notebooks.
  - `gcll/archive/` retains legacy notebook snapshots for reference.
- `data/` – shared inputs arranged by study: `gcll/` (cell counts, treatment logs, spreadsheets), `crc/inputs/` & `crc/cnv/`, `fi/`, `rnaqc/`, and `wbc_variance/`.
- `outputs/` – generated growth kinetics reports exported from notebooks
- `external/` – third-party tools vendored into the repo (`comut`, `signature_analyzer`).
- `Cell_Population/` plus helper scripts (`helper.py`, `model_helper.py`, `plotly_helper.py`) remain at the top level for backward compatibility with existing notebooks.




## Scripts

- `scripts/run_cell_population.py` – fetches the Cell Population inputs (cluster CCFs, mutation CCFs, sif, build_tree_posteriors) for a patient from a Terra workspace via `dalmatian`, then runs `python run_cell_population.py`. Example:
  ```bash
  scripts/run_cell_population.py CLL14-1056 broad-firecloud-ibmwatson/TAG_CLL_Clonal_Kinetic_UMI_PrAN --tree-number 4
  ```
  

- `scripts/build_patient_report.py` – renders the Plotly report from a JSON configuration. Point it at a config (single patient or list) with your data paths and options, then run:
  ```bash
  python scripts/build_patient_report.py --config data/patient_config.json
  ```
  Reports land in `outputs/current/html_reports/`. 

