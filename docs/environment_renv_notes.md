# Environment and renv Notes — Spatial Extension

## Current State

**No packages have been added to renv.lock for spatial analysis.**

The current lockfile supports the existing scRNA-seq pipeline.
Exact package presence must be verified by future audit.

## Future Spatial Dependencies

Future dependencies TBD after object inspection and ref-bio routing.
This section will be updated once actual analysis needs are determined.

## renv Rules

1. **Do NOT run `renv::snapshot()` during scaffolding** — no new packages have been added
2. Only snapshot after intentional package installation
3. Never conflict with Seurat v5.5.0 (the scRNA-seq pipeline's primary dependency)
4. When adding spatial packages, consult official documentation for installation; update renv.lock and this file accordingly
5. Document any new package additions in this file

## System Environment

| Component | Version | Notes |
|-----------|---------|-------|
| OS | macOS (darwin, aarch64) | Apple Silicon |
| R | 4.5.1 | |
| Seurat | 5.5.0 | Primary dependency |
| renv | Active | Manages project dependencies |
| Memory | 16GB | See `docs/memory_management.md` |

## When Ready to Add Spatial Packages

Consult official Seurat and Bioconductor documentation for installation instructions.
Document any new package additions in this file with actual installed versions.
