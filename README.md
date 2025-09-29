# Fludarabine MRD Analysis

This repository collects the exploratory and production material for the CLL MRD growth kinetics projects.

## Repository layout
- `notebooks/` – all Jupyter notebooks grouped by project (`crc`, `fi`, `gcll`, `misc`, and `autogen`).
  - `gcll/active/` holds the working set of GCLL notebooks.
  - `gcll/archive/` retains legacy notebook snapshots for reference.
- `data/` – shared inputs arranged by study: `gcll/` (cell counts, treatment logs, spreadsheets), `crc/inputs/` & `crc/cnv/`, `fi/`, `rnaqc/`, and `wbc_variance/`.
- `outputs/` – generated artefacts exported from notebooks, split into `current/html_reports` and the archived `archive/html_archive` snapshots.
- `src/fludarabine/` – Python package that exposes the helper modules (`helper`, `model_helper`, `plotly_helper`) for reuse across notebooks.
- `external/` – third-party tools vendored into the repo (`comut`, `signature_analyzer`).
- `Cell_Population/` stores the tree/posterior assets pulled from Terra; they remain referenced directly by notebooks for now.

## Working with notebooks
Set up and activate a virtual environment from the repository root, then install the helper package in editable mode so notebooks can import it without manipulating `sys.path`:

```bash
cd /Users/lil/PycharmProjects/Fludarabine
python -m venv venv
source venv/bin/activate
pip install --upgrade pip
pip install -e .
```

After the environment is prepared, launch Jupyter from the repository root so that relative data paths resolve correctly:

```bash
jupyter lab
```

All new notebooks should live inside one of the sub-folders under `notebooks/`. When branching a new analysis, consider duplicating an existing notebook into the relevant project folder.

## Next steps
- Add lightweight instructions for reproducing plots (dependencies, environment setup, etc.).
- Capture data provenance for files housed in `data/` (source, refresh cadence, ownership) so collaborators know what can be regenerated.
- Audit remaining large artefacts (e.g., HTML exports living beside notebooks) and decide which belong under `outputs/` versus per-notebook scratch space.

Keeping these follow-ups in mind will make it easier to share or on-board new collaborators without relying on tribal knowledge.
