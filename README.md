# Fludarabine MRD Analysis

This repository collects the exploratory and production material for the CLL MRD growth kinetics projects.

## Repository layout
- `notebooks/` – all Jupyter notebooks grouped by project (`crc`, `fi`, `gcll`, `misc`, and `autogen`).
  - `gcll/active/` holds the working set of GCLL notebooks.
  - `gcll/archive/` retains legacy notebook snapshots for reference.
- `outputs/` – generated artefacts exported from notebooks, split into `current/html_reports` and the archived `archive/html_archive` snapshots.
- `external/` – third-party tools vendored into the repo (`comut`, `signature_analyzer`).
- `Cell_Population/`, `CRC_cnv/`, `CRC inputs/`, `wbc_files_variance/`, etc. – experimental or source data that still live at the root. See "Next steps" below for consolidation ideas.
- Helper scripts such as `helper.py`, `model_helper.py`, and `plotly_helper.py` remain at the top level for backward compatibility with existing notebooks.

## Working with notebooks
Launch Jupyter from the repository root so that relative imports (for the helper scripts) and data references continue to resolve:

```bash
cd /Users/lil/PycharmProjects/Fludarabine
jupyter lab
```

All new notebooks should live inside one of the sub-folders under `notebooks/`. When branching a new analysis, consider duplicating an existing notebook into the relevant project folder.

## Next steps
- Gather the various CSV/XLSX/txt inputs under a dedicated `data/` hierarchy and update notebook references accordingly.
- Package the helper scripts under a `src/` module (e.g., via `pip install -e .`) so that imports no longer rely on the repository root.
- Add lightweight instructions for reproducing plots (dependencies, environment setup, etc.).
- Capture data provenance for files housed in `CRC_cnv`, `CRC inputs`, `wbc_files_variance`, and similar directories.

Keeping these follow-ups in mind will make it easier to share or on-board new collaborators without relying on tribal knowledge.
