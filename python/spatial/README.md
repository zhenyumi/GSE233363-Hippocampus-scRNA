# python/spatial/ — Optional Spatial Python Workflows

This directory is reserved for future Python spatial transcriptomics workflows.

Current status:

- No Python analysis is implemented here yet.
- Python is expected only for optional cross-modality reproduction work such as Tangram/scanpy, if Phase Spatial-07 decides it is needed.
- The current DG/Hippo RDS-based reproduction stages run in R and do not require Python.

## Rules

- Do not create a Python environment or install packages without an approved plan.
- Do not commit Python virtual environments, caches, notebook checkpoints, AnnData/H5AD files, or generated outputs.
- Keep Python source code here; generated outputs belong in `data/processed/spatial/`, `figures/spatial/`, `results/spatial/`, `reports/spatial/`, or `cache/spatial/`.
- If Tangram is implemented later, record the exact Python version, package versions, CPU/GPU mode, input files, and deviations from the author notebook.

## Planned Use

The only currently planned Python route is optional:

- Phase Spatial-08: Tangram / scanpy cross-modality mapping, based on `docs/original_code_from_paper/Scripts/Python/notebook_tangram.ipynb`.

If Phase Spatial-07 concludes Tangram is unnecessary for the paper-reproduction handoff or the CA1/CA3/DG mitochondrial question, this directory should remain source-only and unused.
