# Repository Structure

This repository has two active R analysis branches and one optional Python branch:

- scRNA-seq ferroDG stage-2 reproduction.
- Spatial transcriptomics reproduction and later CA1/CA3/DG extensions.
- Optional Python/Tangram cross-modality mapping, only if a future validation phase approves it.

The R branches share one R/renv project. Optional Python work must use its own approved environment plan. Source code, documentation, environments, and generated outputs must stay separated.

## Canonical Layout

| Path | Role | Git policy |
|------|------|------------|
| `R/` | R source root; contains branch-specific subdirectories only | Tracked |
| `R/scrna/` | Executable R source code for the scRNA-seq pipeline | Tracked |
| `R/spatial/` | Executable R source code for spatial transcriptomics stages, using `s##_` script names | Tracked |
| `python/spatial/` | Optional future Python source code for Tangram/scanpy spatial workflows | Tracked source only |
| `run_pipeline.R` | scRNA-seq pipeline runner only | Tracked |
| `docs/` | Planning documents, reference registries, result summaries, and copied author code | Tracked, except copyrighted/local PDFs |
| `docs/original_code_from_paper/` | Read-only author-code reference copied for reproduction | Tracked when license/copyright allows |
| `gene_lists/` | Small reviewed gene-list source files used by scripts | Tracked |
| `data/` | Raw input data and generated processed tables/RDS files | Ignored |
| `figures/` | Generated figures from scRNA-seq or spatial scripts | Ignored |
| `results/` | Generated analysis result exports | Ignored |
| `reports/` | Generated local reports and logs | Ignored |
| `cache/` | Short-lived scratch files and temporary outputs | Ignored |
| `renv.lock`, `.Rprofile` | Reproducible R environment entrypoints | Tracked |
| `renv/library/`, `renv/staging/`, `renv/sandbox/` | Local package libraries/caches | Ignored |
| `.venv/`, `venv/`, `__pycache__/`, `.ipynb_checkpoints/` | Local Python environments/caches | Ignored |

## Retired `analysis/` Directory

`analysis/` is not a current code location in this repository.

The historical ferroDG source used paths such as `analysis/stage-2/`, but this reproduction stores maintained scripts under `R/scrna/` and `R/spatial/`. If an empty local `analysis/` directory appears, delete it. Do not add new code, scratch files, or outputs under `analysis/`.

## Output Conventions

scRNA-seq outputs:

- Tables and RDS checkpoints: `data/processed/`
- Figures: `figures/stage-2/`

Spatial outputs:

- Tables and RDS checkpoints: `data/processed/spatial/`
- Figures: `figures/spatial/`
- Larger result exports: `results/spatial/`
- Generated reports/logs: `reports/spatial/`
- Scratch files: `cache/spatial/`

Generated files should not be committed unless they are intentionally promoted into a small reviewed document under `docs/`.

## Temporary File Rules

- Do not write temporary files into the repository root.
- Do not write generated CSV/RDS/PDF/PNG/log files into `R/`, `R/scrna/`, `R/spatial/`, `python/`, `docs/`, or `gene_lists/`.
- Use `cache/` or `cache/spatial/` for short-lived scratch files, and remove them when they are no longer needed.
- Keep `Rplots.pdf` out of the repository root; scripts that create PDF output should use explicit output paths and close graphics devices.

## R Script Rules

- Do not place executable scripts directly in `R/`.
- scRNA-seq scripts remain in `R/scrna/`.
- Spatial scripts remain in `R/spatial/`.
- Add only scRNA-seq scripts to `run_pipeline.R`; that runner is reserved for the scRNA-seq branch.
- A future `run_spatial_pipeline.R` may be added only after the spatial R workflow stabilizes.

## Optional Python Rules

- Python is expected only for optional Tangram/scanpy cross-modality mapping.
- Python source belongs in `python/spatial/`.
- Do not create or commit Python environments, package caches, notebook checkpoints, H5AD/AnnData files, or generated results.
- Do not install Python packages without an approved dependency plan.
- If Python is used later, document Python version, package versions, CPU/GPU mode, input paths, output paths, and deviations from the author notebook.

## Environment Rules

- Use the project renv environment.
- Do not edit `renv.lock` unless packages were intentionally installed, updated, or recorded.
- Local renv package libraries stay untracked.
- Python environment details are not managed by renv and must be documented separately before use.
