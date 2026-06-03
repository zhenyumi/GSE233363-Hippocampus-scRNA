#!/usr/bin/env Rscript
# s02_dg_hippo_metadata.R
# Phase Spatial-02: DG + Hippo metadata/structure verification
# Dependencies: Seurat + base R only (no ggplot2, no dplyr)
# Output: data/processed/spatial/phase02_metadata/

library(Seurat)

# ---- 0. Setup ----
start_time <- Sys.time()
r_version <- paste0(R.version$major, ".", R.version$minor)
seurat_version <- as.character(packageVersion("Seurat"))

base_dir <- "data/raw/GSE233363_official"
out_dir <- "data/processed/spatial/phase02_metadata"
dg_dir <- file.path(out_dir, "DG")
hippo_dir <- file.path(out_dir, "Hippo")
dir.create(dg_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(hippo_dir, recursive = TRUE, showWarnings = FALSE)

cat("=== Phase Spatial-02: DG + Hippo Metadata Verification ===\n")
cat("Start time:", format(start_time), "\n")
cat("R version:", r_version, "\n")
cat("Seurat version:", seurat_version, "\n\n")

all_warnings <- list()
all_checks <- list()  # Collect validation checks across objects

# ---- Helper functions ----

# Record a validation check
add_check <- function(obj_name, check_name, status, details) {
  all_checks[[length(all_checks) + 1]] <<- data.frame(
    object = obj_name,
    check = check_name,
    status = status,
    details = details,
    stringsAsFactors = FALSE
  )
}

# Safely get image spot counts
get_image_spot_counts <- function(sobj) {
  imgs <- Images(sobj)
  results <- data.frame(
    image_name = character(0),
    n_spots = integer(0),
    coord_columns = character(0),
    x_min = numeric(0),
    x_max = numeric(0),
    y_min = numeric(0),
    y_max = numeric(0),
    stringsAsFactors = FALSE
  )
  for (img_name in imgs) {
    tryCatch({
      coords <- GetTissueCoordinates(sobj, image = img_name)
      if (!is.null(coords) && nrow(coords) > 0) {
        coord_cols <- paste(colnames(coords), collapse = ", ")
        # Coordinates are typically imagecol, imagerow
        x_vals <- coords[, 1]
        y_vals <- coords[, 2]
        results <- rbind(results, data.frame(
          image_name = img_name,
          n_spots = nrow(coords),
          coord_columns = coord_cols,
          x_min = min(x_vals, na.rm = TRUE),
          x_max = max(x_vals, na.rm = TRUE),
          y_min = min(y_vals, na.rm = TRUE),
          y_max = max(y_vals, na.rm = TRUE),
          stringsAsFactors = FALSE
        ))
      }
    }, error = function(e) {
      results <<- rbind(results, data.frame(
        image_name = img_name,
        n_spots = NA_integer_,
        coord_columns = paste0("ERROR: ", e$message),
        x_min = NA_real_, x_max = NA_real_,
        y_min = NA_real_, y_max = NA_real_,
        stringsAsFactors = FALSE
      ))
    })
  }
  return(results)
}

# Metadata summary: all columns
get_metadata_summary <- function(meta_df) {
  result <- data.frame(
    column = character(0),
    class = character(0),
    n_unique = integer(0),
    n_na = integer(0),
    levels_preview = character(0),
    stringsAsFactors = FALSE
  )
  for (col_name in colnames(meta_df)) {
    col_data <- meta_df[[col_name]]
    n_u <- length(unique(col_data))
    n_na <- sum(is.na(col_data))
    if (is.factor(col_data)) {
      lvl <- paste(levels(col_data)[1:min(5, length(levels(col_data)))], collapse = ", ")
      lvl <- paste0("factor(", length(levels(col_data)), "): ", lvl)
    } else if (is.character(col_data)) {
      uvals <- unique(col_data)
      lvl <- paste(uvals[1:min(5, length(uvals))], collapse = ", ")
    } else {
      lvl <- paste0("range: ", suppressWarnings(min(col_data, na.rm = TRUE)),
                     " - ", suppressWarnings(max(col_data, na.rm = TRUE)))
    }
    result <- rbind(result, data.frame(
      column = col_name,
      class = class(col_data)[1],
      n_unique = n_u,
      n_na = n_na,
      levels_preview = lvl,
      stringsAsFactors = FALSE
    ))
  }
  return(result)
}

# Distribution table
get_distribution <- function(x) {
  tb <- table(x, useNA = "ifany")
  data.frame(
    level = names(tb),
    count = as.integer(tb),
    stringsAsFactors = FALSE
  )
}

# QC metric audit
audit_qc_metric <- function(x, name) {
  data.frame(
    metric = name,
    class = class(x)[1],
    n_na = sum(is.na(x)),
    min = suppressWarnings(min(x, na.rm = TRUE)),
    max = suppressWarnings(max(x, na.rm = TRUE)),
    median = suppressWarnings(median(x, na.rm = TRUE)),
    mean = suppressWarnings(mean(x, na.rm = TRUE)),
    stringsAsFactors = FALSE
  )
}

# Object structure summary
get_object_structure <- function(sobj, expected) {
  assays <- names(sobj@assays)
  imgs <- Images(sobj)
  reductions <- names(sobj@reductions)
  data.frame(
    field = c("n_genes", "n_spots", "default_assay", "assay_names",
              "n_images", "reduction_names"),
    expected = c(
      as.character(expected$n_genes),
      as.character(expected$n_spots),
      expected$default_assay,
      expected$assays,
      as.character(expected$n_images),
      expected$reductions
    ),
    observed = c(
      as.character(nrow(sobj)),
      as.character(ncol(sobj)),
      DefaultAssay(sobj),
      paste(assays, collapse = ", "),
      as.character(length(imgs)),
      paste(sort(reductions), collapse = ", ")
    ),
    stringsAsFactors = FALSE
  )
}

# ---- Expected values from Stage-01 ----
dg_expected <- list(
  n_genes = 32285,
  n_spots = 2364,
  default_assay = "Spatial",
  assays = "Spatial, SCT",
  n_images = 16,
  reductions = "pca, tsne, umap"
)

hippo_expected <- list(
  n_genes = 17096,
  n_spots = 9921,
  default_assay = "SCT",
  assays = "Spatial, SCT",
  n_images = 16,
  reductions = "pca, tsne, umap"
)

# Sort expected reductions alphabetically for order-independent comparison
dg_expected$reductions <- paste(sort(strsplit(dg_expected$reductions, ", ")[[1]]), collapse = ", ")
hippo_expected$reductions <- paste(sort(strsplit(hippo_expected$reductions, ", ")[[1]]), collapse = ", ")

# Stage-01 Age distributions
dg_age_expected <- c(Young = 797L, Middle = 519L, Old = 1048L)
hippo_age_expected <- c(Young = 2483L, Middle = 2402L, Old = 5036L)

# ---- 1. Process DG ----
cat("=== Processing DG ===\n")
dg_file <- file.path(base_dir, "seurat_Visium_DG_All.rds")
dg_size <- file.info(dg_file)$size
cat("File:", dg_file, "\n")
cat("Size:", format(dg_size, big.mark = ","), "bytes\n")

gc_before_dg <- gc(verbose = FALSE)
dg_t0 <- Sys.time()

cat("Loading DG object...\n")
sobj <- readRDS(dg_file)
dg_t1 <- Sys.time()
gc_after_dg <- gc(verbose = FALSE)

cat("Loaded. Dimensions:", nrow(sobj), "genes x", ncol(sobj), "spots\n")

# 1.2 Object structure
dg_structure <- get_object_structure(sobj, dg_expected)
write.csv(dg_structure, file.path(dg_dir, "dg_object_structure_summary.csv"), row.names = FALSE)

# Check structure
for (i in seq_len(nrow(dg_structure))) {
  if (dg_structure$expected[i] == dg_structure$observed[i]) {
    add_check("DG", paste0("structure_", dg_structure$field[i]), "PASS",
              paste0("expected=", dg_structure$expected[i]))
  } else {
    add_check("DG", paste0("structure_", dg_structure$field[i]), "FAIL",
              paste0("expected=", dg_structure$expected[i], " observed=", dg_structure$observed[i]))
  }
}

# 1.3 Extract metadata
meta <- sobj@meta.data

# 1.4 Image spot counts
cat("Extracting image spot counts...\n")
dg_img_counts <- get_image_spot_counts(sobj)
dg_img_sum <- sum(dg_img_counts$n_spots, na.rm = TRUE)
dg_img_counts_summary <- rbind(
  dg_img_counts,
  data.frame(image_name = "TOTAL_SUM", n_spots = dg_img_sum,
             coord_columns = "", x_min = NA, x_max = NA, y_min = NA, y_max = NA,
             stringsAsFactors = FALSE),
  data.frame(image_name = "OBJECT_NCOL", n_spots = ncol(sobj),
             coord_columns = "", x_min = NA, x_max = NA, y_min = NA, y_max = NA,
             stringsAsFactors = FALSE)
)
write.csv(dg_img_counts_summary, file.path(dg_dir, "dg_image_spot_counts.csv"), row.names = FALSE)

if (dg_img_sum == ncol(sobj)) {
  add_check("DG", "image_spot_count_consistency", "PASS",
            paste0("sum=", dg_img_sum, " == ncol=", ncol(sobj)))
} else {
  add_check("DG", "image_spot_count_consistency", "WARNING",
            paste0("sum=", dg_img_sum, " != ncol=", ncol(sobj),
                   " diff=", abs(dg_img_sum - ncol(sobj))))
}

# 1.5 Metadata summary
dg_meta_summary <- get_metadata_summary(meta)
write.csv(dg_meta_summary, file.path(dg_dir, "dg_metadata_summary.csv"), row.names = FALSE)

# 1.6 Core column existence checks
for (col in c("Age", "Region", "GEMgroup", "Position",
              "percentMito", "percentRibo", "nCount_Spatial", "nFeature_Spatial")) {
  if (col %in% colnames(meta)) {
    add_check("DG", paste0("column_exists_", col), "PASS", paste0("class=", class(meta[[col]])[1]))
  } else {
    add_check("DG", paste0("column_exists_", col), ifelse(col %in% c("Age", "Region", "GEMgroup"), "FAIL", "WARNING"),
              "MISSING")
  }
}

# SCT_snn_res.1 / seurat_clusters (non-blocking)
for (col in c("SCT_snn_res.1", "seurat_clusters")) {
  if (col %in% colnames(meta)) {
    add_check("DG", paste0("column_exists_", col), "PASS",
              paste0("n_levels=", length(unique(meta[[col]])), " (non-blocking)"))
  } else {
    add_check("DG", paste0("column_exists_", col), "PASS", "absent (non-blocking)")
  }
}

# 1.7 Age distribution
cat("Computing distributions...\n")
dg_age_dist <- get_distribution(meta$Age)
dg_age_dist$expected <- dg_age_expected[dg_age_dist$level]
dg_age_dist$diff <- dg_age_dist$count - dg_age_dist$expected
write.csv(dg_age_dist, file.path(dg_dir, "dg_age_distribution.csv"), row.names = FALSE)

# Age check
for (i in seq_len(nrow(dg_age_dist))) {
  lvl <- dg_age_dist$level[i]
  if (!is.na(dg_age_expected[lvl]) && dg_age_dist$count[i] == dg_age_expected[lvl]) {
    add_check("DG", paste0("age_", lvl), "PASS",
              paste0("expected=", dg_age_expected[lvl]))
  } else if (is.na(dg_age_expected[lvl])) {
    add_check("DG", paste0("age_", lvl), "WARNING",
              paste0("unexpected level, count=", dg_age_dist$count[i]))
  } else {
    add_check("DG", paste0("age_", lvl), "FAIL",
              paste0("expected=", dg_age_expected[lvl],
                     " observed=", dg_age_dist$count[i],
                     " diff=", dg_age_dist$diff[i]))
  }
}

# 1.8 Region distribution
dg_region_dist <- get_distribution(meta$Region)
write.csv(dg_region_dist, file.path(dg_dir, "dg_region_distribution.csv"), row.names = FALSE)

# 1.9 Cross-tabs
dg_gemgroup_age <- as.data.frame.matrix(table(meta$GEMgroup, meta$Age))
dg_gemgroup_age <- cbind(GEMgroup = rownames(dg_gemgroup_age), dg_gemgroup_age)
write.csv(dg_gemgroup_age, file.path(dg_dir, "dg_gemgroup_age_crosstab.csv"), row.names = FALSE)

dg_gemgroup_pos <- as.data.frame.matrix(table(meta$GEMgroup, meta$Position))
dg_gemgroup_pos <- cbind(GEMgroup = rownames(dg_gemgroup_pos), dg_gemgroup_pos)
write.csv(dg_gemgroup_pos, file.path(dg_dir, "dg_gemgroup_position_crosstab.csv"), row.names = FALSE)

dg_age_region <- as.data.frame.matrix(table(meta$Age, meta$Region))
dg_age_region <- cbind(Age = rownames(dg_age_region), dg_age_region)
write.csv(dg_age_region, file.path(dg_dir, "dg_age_region_crosstab.csv"), row.names = FALSE)

# 1.10 GEMgroup mapping validation
# Expected: GEMgroup 1-4 → Young, 5-8 → Middle, 9-16 → Old
# Expected: odd GEMgroup → Posterior, even → Anterior
gemgroup_vals <- levels(meta$GEMgroup)
for (gg in gemgroup_vals) {
  gg_num <- suppressWarnings(as.integer(gg))
  if (is.na(gg_num)) next
  # Age check
  expected_age <- if (gg_num <= 4) "Young" else if (gg_num <= 8) "Middle" else "Old"
  actual_ages <- unique(as.character(meta$Age[meta$GEMgroup == gg]))
  if (length(actual_ages) == 1 && actual_ages == expected_age) {
    add_check("DG", paste0("gemgroup_", gg, "_age"), "PASS",
              paste0("expected=", expected_age))
  } else {
    add_check("DG", paste0("gemgroup_", gg, "_age"), "FAIL",
              paste0("expected=", expected_age, " observed=", paste(actual_ages, collapse = "/")))
  }
  # Position check
  expected_pos <- if (gg_num %% 2 == 1) "Posterior" else "Anterior"
  actual_pos <- unique(as.character(meta$Position[meta$GEMgroup == gg]))
  if (length(actual_pos) == 1 && actual_pos == expected_pos) {
    add_check("DG", paste0("gemgroup_", gg, "_position"), "PASS",
              paste0("expected=", expected_pos))
  } else {
    add_check("DG", paste0("gemgroup_", gg, "_position"), "FAIL",
              paste0("expected=", expected_pos, " observed=", paste(actual_pos, collapse = "/")))
  }
}

# 1.11 QC metric audit
cat("Auditing QC metrics...\n")
qc_cols <- intersect(c("percentMito", "percentRibo", "nCount_Spatial", "nFeature_Spatial"), colnames(meta))
dg_qc_audit <- do.call(rbind, lapply(qc_cols, function(cn) audit_qc_metric(meta[[cn]], cn)))
write.csv(dg_qc_audit, file.path(dg_dir, "dg_qc_metric_summary.csv"), row.names = FALSE)

# Flag impossible values
for (cn in intersect(c("percentMito", "percentRibo"), colnames(meta))) {
  vals <- meta[[cn]]
  if (any(vals < 0, na.rm = TRUE) || any(vals > 100, na.rm = TRUE)) {
    add_check("DG", paste0("qc_impossible_", cn), "WARNING",
              paste0("min=", min(vals, na.rm = TRUE), " max=", max(vals, na.rm = TRUE),
                     " (expected 0-100)"))
  } else {
    add_check("DG", paste0("qc_impossible_", cn), "PASS",
              paste0("range OK: ", min(vals, na.rm = TRUE), " - ", max(vals, na.rm = TRUE)))
  }
}

# 1.12 Write DG validation checks
dg_checks <- do.call(rbind, all_checks[which(sapply(all_checks, function(x) x$object == "DG"))])
write.csv(dg_checks, file.path(dg_dir, "dg_validation_checks.csv"), row.names = FALSE)

dg_n_fail <- sum(dg_checks$status == "FAIL")
dg_n_warn <- sum(dg_checks$status == "WARNING")
cat("DG complete:", nrow(dg_checks), "checks,", dg_n_fail, "FAIL,", dg_n_warn, "WARNING\n")

# Cleanup DG
rm(sobj)
gc()
cat("DG object released.\n\n")

# ---- 2. Process Hippo ----
cat("=== Processing Hippo ===\n")
hippo_file <- file.path(base_dir, "seurat_Visium_Hippo_All.rds")
hippo_size <- file.info(hippo_file)$size
cat("File:", hippo_file, "\n")
cat("Size:", format(hippo_size, big.mark = ","), "bytes\n")

gc_before_hippo <- gc(verbose = FALSE)
hippo_t0 <- Sys.time()

cat("Loading Hippo object...\n")
sobj <- readRDS(hippo_file)
hippo_t1 <- Sys.time()
gc_after_hippo <- gc(verbose = FALSE)

cat("Loaded. Dimensions:", nrow(sobj), "genes x", ncol(sobj), "spots\n")

# 2.2 Object structure
hippo_structure <- get_object_structure(sobj, hippo_expected)
write.csv(hippo_structure, file.path(hippo_dir, "hippo_object_structure_summary.csv"), row.names = FALSE)

for (i in seq_len(nrow(hippo_structure))) {
  if (hippo_structure$expected[i] == hippo_structure$observed[i]) {
    add_check("Hippo", paste0("structure_", hippo_structure$field[i]), "PASS",
              paste0("expected=", hippo_structure$expected[i]))
  } else {
    add_check("Hippo", paste0("structure_", hippo_structure$field[i]), "FAIL",
              paste0("expected=", hippo_structure$expected[i], " observed=", hippo_structure$observed[i]))
  }
}

# 2.3 Extract metadata
meta <- sobj@meta.data

# 2.4 Image spot counts
cat("Extracting image spot counts...\n")
hippo_img_counts <- get_image_spot_counts(sobj)
hippo_img_sum <- sum(hippo_img_counts$n_spots, na.rm = TRUE)
hippo_img_counts_summary <- rbind(
  hippo_img_counts,
  data.frame(image_name = "TOTAL_SUM", n_spots = hippo_img_sum,
             coord_columns = "", x_min = NA, x_max = NA, y_min = NA, y_max = NA,
             stringsAsFactors = FALSE),
  data.frame(image_name = "OBJECT_NCOL", n_spots = ncol(sobj),
             coord_columns = "", x_min = NA, x_max = NA, y_min = NA, y_max = NA,
             stringsAsFactors = FALSE)
)
write.csv(hippo_img_counts_summary, file.path(hippo_dir, "hippo_image_spot_counts.csv"), row.names = FALSE)

if (hippo_img_sum == ncol(sobj)) {
  add_check("Hippo", "image_spot_count_consistency", "PASS",
            paste0("sum=", hippo_img_sum, " == ncol=", ncol(sobj)))
} else {
  add_check("Hippo", "image_spot_count_consistency", "WARNING",
            paste0("sum=", hippo_img_sum, " != ncol=", ncol(sobj),
                   " diff=", abs(hippo_img_sum - ncol(sobj))))
}

# 2.5 Metadata summary
hippo_meta_summary <- get_metadata_summary(meta)
write.csv(hippo_meta_summary, file.path(hippo_dir, "hippo_metadata_summary.csv"), row.names = FALSE)

# 2.6 Core column existence checks
for (col in c("Age", "Region", "GEMgroup", "Position",
              "percentMito", "percentRibo", "nCount_Spatial", "nFeature_Spatial")) {
  if (col %in% colnames(meta)) {
    add_check("Hippo", paste0("column_exists_", col), "PASS", paste0("class=", class(meta[[col]])[1]))
  } else {
    add_check("Hippo", paste0("column_exists_", col), ifelse(col %in% c("Age", "Region", "GEMgroup"), "FAIL", "WARNING"),
              "MISSING")
  }
}

for (col in c("SCT_snn_res.1", "seurat_clusters")) {
  if (col %in% colnames(meta)) {
    add_check("Hippo", paste0("column_exists_", col), "PASS",
              paste0("n_levels=", length(unique(meta[[col]])), " (non-blocking)"))
  } else {
    add_check("Hippo", paste0("column_exists_", col), "PASS", "absent (non-blocking)")
  }
}

# 2.7 Age distribution
cat("Computing distributions...\n")
hippo_age_dist <- get_distribution(meta$Age)
hippo_age_dist$expected <- hippo_age_expected[hippo_age_dist$level]
hippo_age_dist$diff <- hippo_age_dist$count - hippo_age_dist$expected
write.csv(hippo_age_dist, file.path(hippo_dir, "hippo_age_distribution.csv"), row.names = FALSE)

for (i in seq_len(nrow(hippo_age_dist))) {
  lvl <- hippo_age_dist$level[i]
  if (!is.na(hippo_age_expected[lvl]) && hippo_age_dist$count[i] == hippo_age_expected[lvl]) {
    add_check("Hippo", paste0("age_", lvl), "PASS",
              paste0("expected=", hippo_age_expected[lvl]))
  } else if (is.na(hippo_age_expected[lvl])) {
    add_check("Hippo", paste0("age_", lvl), "WARNING",
              paste0("unexpected level, count=", hippo_age_dist$count[i]))
  } else {
    add_check("Hippo", paste0("age_", lvl), "FAIL",
              paste0("expected=", hippo_age_expected[lvl],
                     " observed=", hippo_age_dist$count[i],
                     " diff=", hippo_age_dist$diff[i]))
  }
}

# 2.8 Region distribution
hippo_region_dist <- get_distribution(meta$Region)
write.csv(hippo_region_dist, file.path(hippo_dir, "hippo_region_distribution.csv"), row.names = FALSE)

# 2.9 Cross-tabs
hippo_gemgroup_age <- as.data.frame.matrix(table(meta$GEMgroup, meta$Age))
hippo_gemgroup_age <- cbind(GEMgroup = rownames(hippo_gemgroup_age), hippo_gemgroup_age)
write.csv(hippo_gemgroup_age, file.path(hippo_dir, "hippo_gemgroup_age_crosstab.csv"), row.names = FALSE)

hippo_gemgroup_pos <- as.data.frame.matrix(table(meta$GEMgroup, meta$Position))
hippo_gemgroup_pos <- cbind(GEMgroup = rownames(hippo_gemgroup_pos), hippo_gemgroup_pos)
write.csv(hippo_gemgroup_pos, file.path(hippo_dir, "hippo_gemgroup_position_crosstab.csv"), row.names = FALSE)

hippo_age_region <- as.data.frame.matrix(table(meta$Age, meta$Region))
hippo_age_region <- cbind(Age = rownames(hippo_age_region), hippo_age_region)
write.csv(hippo_age_region, file.path(hippo_dir, "hippo_age_region_crosstab.csv"), row.names = FALSE)

# 2.10 GEMgroup mapping validation
gemgroup_vals <- levels(meta$GEMgroup)
for (gg in gemgroup_vals) {
  gg_num <- suppressWarnings(as.integer(gg))
  if (is.na(gg_num)) next
  expected_age <- if (gg_num <= 4) "Young" else if (gg_num <= 8) "Middle" else "Old"
  actual_ages <- unique(as.character(meta$Age[meta$GEMgroup == gg]))
  if (length(actual_ages) == 1 && actual_ages == expected_age) {
    add_check("Hippo", paste0("gemgroup_", gg, "_age"), "PASS",
              paste0("expected=", expected_age))
  } else {
    add_check("Hippo", paste0("gemgroup_", gg, "_age"), "FAIL",
              paste0("expected=", expected_age, " observed=", paste(actual_ages, collapse = "/")))
  }
  expected_pos <- if (gg_num %% 2 == 1) "Posterior" else "Anterior"
  actual_pos <- unique(as.character(meta$Position[meta$GEMgroup == gg]))
  if (length(actual_pos) == 1 && actual_pos == expected_pos) {
    add_check("Hippo", paste0("gemgroup_", gg, "_position"), "PASS",
              paste0("expected=", expected_pos))
  } else {
    add_check("Hippo", paste0("gemgroup_", gg, "_position"), "FAIL",
              paste0("expected=", expected_pos, " observed=", paste(actual_pos, collapse = "/")))
  }
}

# 2.11 QC metric audit
qc_cols <- intersect(c("percentMito", "percentRibo", "nCount_Spatial", "nFeature_Spatial"), colnames(meta))
hippo_qc_audit <- do.call(rbind, lapply(qc_cols, function(cn) audit_qc_metric(meta[[cn]], cn)))
write.csv(hippo_qc_audit, file.path(hippo_dir, "hippo_qc_metric_summary.csv"), row.names = FALSE)

for (cn in intersect(c("percentMito", "percentRibo"), colnames(meta))) {
  vals <- meta[[cn]]
  if (any(vals < 0, na.rm = TRUE) || any(vals > 100, na.rm = TRUE)) {
    add_check("Hippo", paste0("qc_impossible_", cn), "WARNING",
              paste0("min=", min(vals, na.rm = TRUE), " max=", max(vals, na.rm = TRUE),
                     " (expected 0-100)"))
  } else {
    add_check("Hippo", paste0("qc_impossible_", cn), "PASS",
              paste0("range OK: ", min(vals, na.rm = TRUE), " - ", max(vals, na.rm = TRUE)))
  }
}

# 2.12 Write Hippo validation checks
hippo_checks <- do.call(rbind, all_checks[which(sapply(all_checks, function(x) x$object == "Hippo"))])
write.csv(hippo_checks, file.path(hippo_dir, "hippo_validation_checks.csv"), row.names = FALSE)

hippo_n_fail <- sum(hippo_checks$status == "FAIL")
hippo_n_warn <- sum(hippo_checks$status == "WARNING")
cat("Hippo complete:", nrow(hippo_checks), "checks,", hippo_n_fail, "FAIL,", hippo_n_warn, "WARNING\n")

# Cleanup Hippo
rm(sobj)
gc()
cat("Hippo object released.\n\n")

# ---- 3. Cross-Object Comparison ----
cat("=== Cross-Object Comparison ===\n")

# Region count consistency: Hippo ML+GCL+Hilus vs DG total
hippo_ml <- as.integer(hippo_region_dist$count[hippo_region_dist$level == "ML"])
hippo_gcl <- as.integer(hippo_region_dist$count[hippo_region_dist$level == "GCL"])
hippo_hilus <- as.integer(hippo_region_dist$count[hippo_region_dist$level == "Hilus"])
hippo_dg_sum <- sum(c(hippo_ml, hippo_gcl, hippo_hilus), na.rm = TRUE)
dg_total <- 2364L  # from Stage-01 inspection

cross_compare <- data.frame(
  check = c("Hippo_ML_plus_GCL_plus_Hilus_vs_DG_total",
            "Hippo_ML_count",
            "Hippo_GCL_count",
            "Hippo_Hilus_count",
            "Hippo_dg_region_sum",
            "DG_total_spots"),
  value = c(
    ifelse(hippo_dg_sum == dg_total, "CONSISTENT", "DISCREPANCY"),
    as.character(hippo_ml),
    as.character(hippo_gcl),
    as.character(hippo_hilus),
    as.character(hippo_dg_sum),
    as.character(dg_total)
  ),
  note = c(
    "count-level consistency only; not spot-level subset proof",
    "", "", "", "", ""
  ),
  stringsAsFactors = FALSE
)
write.csv(cross_compare, file.path(out_dir, "cross_object_comparison.csv"), row.names = FALSE)

if (hippo_dg_sum == dg_total) {
  add_check("CrossObject", "region_count_consistency", "PASS",
            paste0("Hippo ML+GCL+Hilus=", hippo_dg_sum, " == DG total=", dg_total))
} else {
  add_check("CrossObject", "region_count_consistency", "WARNING",
            paste0("Hippo ML+GCL+Hilus=", hippo_dg_sum, " != DG total=", dg_total,
                   " diff=", abs(hippo_dg_sum - dg_total)))
}

# ---- 4. Validation Summary ----
cat("\n=== Writing Validation Summary ===\n")

all_checks_df <- do.call(rbind, all_checks)
n_total <- nrow(all_checks_df)
n_pass <- sum(all_checks_df$status == "PASS")
n_fail <- sum(all_checks_df$status == "FAIL")
n_warn <- sum(all_checks_df$status == "WARNING")
blocking_cols <- c("Age", "Region", "GEMgroup")
blocking_fails <- all_checks_df[all_checks_df$status == "FAIL" &
                                  sapply(all_checks_df$check, function(x) any(sapply(blocking_cols, function(bc) grepl(bc, x)))), ]
phase_status <- if (nrow(blocking_fails) > 0) "failed/blocking" else if (n_fail > 0) "completed_with_failures" else if (n_warn > 0) "completed_with_warnings" else "passed"

end_time <- Sys.time()
dg_load_time <- as.numeric(difftime(dg_t1, dg_t0, units = "secs"))
hippo_load_time <- as.numeric(difftime(hippo_t1, hippo_t0, units = "secs"))

sink(file.path(out_dir, "phase02_validation_summary.txt"))
cat("=== Phase Spatial-02 Validation Summary ===\n")
cat("Date:", format(end_time), "\n\n")
cat("Overall status:", phase_status, "\n")
cat("Total checks:", n_total, "\n")
cat("PASS:", n_pass, "\n")
cat("FAIL:", n_fail, "\n")
cat("WARNING:", n_warn, "\n\n")
cat("Blocking failures (core columns):", nrow(blocking_fails), "\n")
if (nrow(blocking_fails) > 0) {
  cat("  Blocking checks:\n")
  for (i in seq_len(nrow(blocking_fails))) {
    cat("    -", blocking_fails$check[i], ":", blocking_fails$details[i], "\n")
  }
}
cat("\n--- DG Summary ---\n")
cat("Checks:", nrow(dg_checks), "\n")
cat("FAIL:", sum(dg_checks$status == "FAIL"), "\n")
cat("WARNING:", sum(dg_checks$status == "WARNING"), "\n")
cat("Load time:", round(dg_load_time, 1), "s\n\n")
cat("--- Hippo Summary ---\n")
cat("Checks:", nrow(hippo_checks), "\n")
cat("FAIL:", sum(hippo_checks$status == "FAIL"), "\n")
cat("WARNING:", sum(hippo_checks$status == "WARNING"), "\n")
cat("Load time:", round(hippo_load_time, 1), "s\n\n")
cat("--- Cross-Object ---\n")
cat("Region count consistency:", cross_compare$value[1], "\n")
cat("Note:", cross_compare$note[1], "\n\n")
cat("--- All FAIL/WARNING checks ---\n")
fw <- all_checks_df[all_checks_df$status != "PASS", ]
if (nrow(fw) > 0) {
  for (i in seq_len(nrow(fw))) {
    cat(sprintf("[%s] %s: %s — %s\n", fw$status[i], fw$object[i], fw$check[i], fw$details[i]))
  }
} else {
  cat("None.\n")
}
sink()

# Provenance
sink(file.path(out_dir, "phase02_provenance.txt"))
cat("=== Phase Spatial-02 Provenance ===\n")
cat("Date:", format(end_time), "\n")
cat("R version:", r_version, "\n")
cat("Seurat version:", seurat_version, "\n")
cat("Script: R/spatial/s02_dg_hippo_metadata.R\n")
cat("DG file:", dg_file, "\n")
cat("DG size:", format(dg_size, big.mark = ","), "bytes\n")
cat("DG load time:", round(dg_load_time, 1), "s\n")
cat("Hippo file:", hippo_file, "\n")
cat("Hippo size:", format(hippo_size, big.mark = ","), "bytes\n")
cat("Hippo load time:", round(hippo_load_time, 1), "s\n")
cat("Total checks:", n_total, "\n")
cat("PASS:", n_pass, " FAIL:", n_fail, " WARNING:", n_warn, "\n")
cat("Overall status:", phase_status, "\n")
if (length(all_warnings) > 0) {
  cat("\nWarnings collected:\n")
  for (w in all_warnings) cat(" -", w, "\n")
}
sink()

cat("\n=== Phase Spatial-02 Complete ===\n")
cat("Status:", phase_status, "\n")
cat("Checks:", n_total, "(PASS:", n_pass, "FAIL:", n_fail, "WARNING:", n_warn, ")\n")
cat("Output:", out_dir, "\n")
cat("End time:", format(end_time), "\n")
