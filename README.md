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




