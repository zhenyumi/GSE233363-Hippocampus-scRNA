# Memory Management — Spatial Transcriptomics

## System Constraints

| Component | Value |
|-----------|-------|
| OS | macOS (darwin, aarch64) |
| RAM | 16GB |
| Chip | Apple Silicon |
| Seurat version | v5.5.0 |

## Rules for Spatial Analysis

1. **Load one spatial object at a time** — do not hold multiple Visium objects in memory simultaneously
2. **Use lightweight summaries** — do not create duplicate RDS copies in memory
3. **Call `gc()` after large operations** — especially after loading, subsetting, or computing
4. **Memory management specifics**: Deferred to official Seurat documentation and actual object structure at first load
5. **Save intermediate results to disk** — then release the in-memory object

## Spatial Object Memory Behavior

Spatial object file sizes are unverified and may be large; measure file size and memory behavior at first load.

When spatial RDS files become available:
1. Record file size on disk before loading
2. Record R session memory usage after loading (use available base R/system tools or documented packages after dependency review)
3. Record memory usage after subsetting or other operations
4. Record memory usage after releasing large objects
5. Document findings in the inspection script output

## Memory Management Pattern for Spatial Scripts

```r
# Load
obj <- readRDS("data/raw/spatial/sample.rds")
gc()

# Inspect (lightweight)
cat("Dimensions:", dim(obj), "\n")
cat("Assays:", Assays(obj), "\n")

# Record provenance/object structure
# ... inspection ...

# Save lightweight CSV/text inspection summary
write.csv(summary_df, "data/processed/spatial/inspection_summary.csv", row.names = FALSE)

# Release
rm(obj, summary_df)
gc()
```

## Known Memory Patterns (scRNA-seq Pipeline)

The scRNA-seq pipeline uses ~4-6GB peak memory. Spatial objects with images may use more.
See `docs/reproduction_plan.md` Section 7 for scRNA-seq memory management details.

## Troubleshooting

- If R runs out of memory: restart session and consult official Seurat documentation for memory optimization
- If loading fails: check file size first; very large files (>2GB) may need chunked processing
- If `gc()` doesn't free enough memory: restart R session between major operations
