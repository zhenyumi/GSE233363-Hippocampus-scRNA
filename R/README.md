# R/ — R Analysis Source

R source code is organized by analysis branch:

| Path | Purpose |
|------|---------|
| `R/scrna/` | scRNA-seq ferroDG stage-2 reproduction scripts |
| `R/spatial/` | spatial transcriptomics reproduction and extension scripts |

Do not place executable analysis scripts directly in `R/`. Add them to the appropriate subdirectory and update the relevant README, plan, and runner.

Generated outputs do not belong under `R/`; use `data/`, `figures/`, `results/`, `reports/`, or `cache/` according to `docs/repository_structure.md`.
