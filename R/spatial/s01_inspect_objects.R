## =========================================================
## Script: s01_inspect_objects.R
## Purpose: Inspect Visium spatial Seurat objects — record assays,
##          spatial slots, metadata, dimensions, and provenance
## Usage:
##   Rscript R/spatial/s01_inspect_objects.R [object_key] [--confirm-large]
##
##   object_key: DG (default), Hippo, WholeBrain
##   --confirm-large: required for WholeBrain (4.6 GB)
##
## Output: 7 lightweight TXT/CSV files in
##         data/processed/spatial/inspection/{object_key}/
##         + object_inventory.csv in data/processed/spatial/inspection/
## =========================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(jsonlite)
})

## =========================================================
## 0. Argument parsing
## =========================================================
args <- commandArgs(trailingOnly = TRUE)

# Extract flags
confirm_large <- "--confirm-large" %in% args
positional     <- args[!grepl("^--", args)]

object_key <- if (length(positional) >= 1) positional[1] else "DG"

# Reject Chromium explicitly
if (tolower(object_key) == "chromium") {
  stop("Chromium is a scRNA-seq object and must NOT be loaded by spatial scripts. ",
       "Use the scRNA-seq pipeline (R/01_load_and_subset.R) for Chromium.")
}

# Validate object_key
valid_keys <- c("DG", "Hippo", "WholeBrain")
if (!object_key %in% valid_keys) {
  stop("Invalid object_key '", object_key, "'. Must be one of: ",
       paste(valid_keys, collapse = ", "))
}

# WholeBrain requires --confirm-large
if (object_key == "WholeBrain" && !confirm_large) {
  stop("WholeBrain RDS is 4.6 GB. Loading it may exceed available memory (16 GB). ",
       "Re-run with: Rscript R/spatial/s01_inspect_objects.R WholeBrain --confirm-large")
}

## =========================================================
## 0b. Object configuration
## =========================================================
object_config <- list(
  DG = list(
    filename    = "seurat_Visium_DG_All.rds",
    meta_key    = "seurat_Visium_DG_All.rds",
    label       = "Visium DG"
  ),
  Hippo = list(
    filename    = "seurat_Visium_Hippo_All.rds",
    meta_key    = "seurat_Visium_Hippo_All.rds",
    label       = "Visium Hippo"
  ),
  WholeBrain = list(
    filename    = "seurat_Visium_WholeBrain.rds",
    meta_key    = "seurat_Visium_WholeBrain.rds",
    label       = "Visium WholeBrain"
  )
)

cfg <- object_config[[object_key]]
rds_filename <- cfg$filename
rds_label    <- cfg$label

message("=== Spatial s01: Inspect ", rds_label, " Object ===")
message("Object key: ", object_key)
message("Start time: ", Sys.time())

## =========================================================
## 0c. Paths
## =========================================================
project_dir  <- getwd()
data_dir     <- file.path(project_dir, "data", "raw", "GSE233363_official")
inspect_root <- file.path(project_dir, "data", "processed", "spatial", "inspection")
output_dir   <- file.path(inspect_root, object_key)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

rds_file <- file.path(data_dir, rds_filename)

if (!file.exists(rds_file)) {
  stop("Spatial RDS not found: ", rds_file,
       "\nPlace ", rds_filename, " in data/raw/GSE233363_official/")
}

## =========================================================
## 0d. Seurat version validation
## =========================================================
seurat_version <- tryCatch(packageVersion("Seurat"), error = function(e) NULL)
if (is.null(seurat_version)) {
  stop("Cannot determine Seurat version. ",
       "Ensure Seurat is installed (use renv::restore()).")
}
seurat_vstr <- as.character(seurat_version)
if (package_version(seurat_vstr) < "5.0.0") {
  stop("Seurat version ", seurat_vstr, " is too old. ",
       "This script requires Seurat >= 5.0.0. ",
       "Use renv::restore() to install the correct version.")
}
message("  Seurat version: ", seurat_vstr, " (>= 5.0.0 OK)")

## =========================================================
## 1. File info + METADATA.json
## =========================================================
message("[1/9] Checking file info...")

fi <- file.info(rds_file)

meta_file <- file.path(data_dir, "METADATA.json")
meta <- NULL
expected_md5 <- NA_character_
if (file.exists(meta_file)) {
  meta <- tryCatch(fromJSON(meta_file), error = function(e) NULL)
  if (!is.null(meta) && !is.null(meta$files) && is.data.frame(meta$files) &&
      "filename" %in% names(meta$files) && "md5" %in% names(meta$files)) {
    md5_row <- tryCatch(
      meta$files[meta$files$filename == cfg$meta_key, ],
      error = function(e) NULL
    )
    if (!is.null(md5_row) && nrow(md5_row) == 1) {
      expected_md5 <- md5_row$md5
    }
  }
}

## =========================================================
## 2. Load RDS + Seurat class check
## =========================================================
message("[2/9] Loading ", rds_filename, " (", round(fi$size / 1e6, 1), " MB)...")

load_result <- tryCatch({
  sobj <- readRDS(rds_file)
  list(ok = TRUE, obj = sobj, error = NULL)
}, error = function(e) {
  list(ok = FALSE, obj = NULL, error = conditionMessage(e))
})

if (!load_result$ok) {
  stop("readRDS() failed on ", object_key, ": ", load_result$error,
       "\nCannot proceed. Check RDS file integrity.")
}

sobj <- load_result$obj

# Class check: not Seurat → write minimal output, stop
obj_class <- tryCatch(class(sobj)[1], error = function(e) "ERROR")

if (!inherits(sobj, "Seurat")) {
  message("  ERROR: Object is not Seurat. Actual class: ", obj_class)
  writeLines(c(
    paste("=== Provenance:", rds_filename, "==="),
    paste("File path:", rds_file),
    paste("File size:", fi$size, "bytes"),
    paste("File mtime:", fi$mtime),
    paste("Actual class:", obj_class),
    "STATUS: NOT A SEURAT OBJECT — inspection stopped",
    paste("Inspection date:", Sys.time())
  ), file.path(output_dir, "provenance.txt"))
  writeLines(c(
    "=== Inspection Summary ===",
    paste("Object class:", obj_class),
    "STATUS: NOT A SEURAT OBJECT — inspection stopped"
  ), file.path(output_dir, "inspection_summary.txt"))
  stop("Loaded object is not Seurat (class: ", obj_class, "). ",
       "Cannot call Seurat API. Stopped.")
}

message("  Loaded successfully. Class: Seurat.")

## =========================================================
## 3. Basic structure
## =========================================================
message("[3/9] Recording basic structure...")

obj_dims      <- tryCatch(dim(sobj), error = function(e) c(NA_integer_, NA_integer_))
n_genes       <- obj_dims[1]
n_cells       <- obj_dims[2]
default_assay <- tryCatch(DefaultAssay(sobj), error = function(e) "ERROR")

# Assays() first, fallback to names(sobj@assays)
all_assays <- tryCatch(Assays(sobj), error = function(e) {
  tryCatch(names(sobj@assays), error = function(e2) character(0))
})

# capture.output(str(...))
str_lines <- tryCatch(
  capture.output(str(sobj, max.level = 1, list.len = 20)),
  error = function(e) paste("str() failed:", e$message)
)

slot_names <- tryCatch(slotNames(sobj), error = function(e) "ERROR")

## =========================================================
## 4. Assay inspection
## =========================================================
message("[4/9] Inspecting assays...")

assay_rows <- list()
for (aname in all_assays) {
  a <- tryCatch(sobj[[aname]], error = function(e) NULL)
  if (is.null(a)) {
    assay_rows[[length(assay_rows) + 1]] <- data.frame(
      assay = aname, class = "ERROR", n_features = NA_integer_,
      layers = "ERROR", counts_status = "ERROR",
      data_status = "ERROR", scale_data_status = "ERROR",
      stringsAsFactors = FALSE
    )
    next
  }

  a_class <- tryCatch(class(a)[1], error = function(e) "ERROR")
  n_feat  <- tryCatch(nrow(a), error = function(e) NA_integer_)

  lyr <- tryCatch(Layers(a), error = function(e) NULL)
  if (!is.null(lyr)) {
    lyr_str <- paste(lyr, collapse = ", ")
  } else {
    lyr_str <- "UNKNOWN"
  }

  counts_status <- tryCatch({
    GetAssayData(a, layer = "counts"); "present"
  }, error = function(e) "missing")

  data_status <- tryCatch({
    GetAssayData(a, layer = "data"); "present"
  }, error = function(e) "missing")

  scale_data_status <- tryCatch({
    GetAssayData(a, layer = "scale.data"); "present"
  }, error = function(e) "missing")

  assay_rows[[length(assay_rows) + 1]] <- data.frame(
    assay = aname, class = a_class, n_features = n_feat,
    layers = lyr_str, counts_status = counts_status,
    data_status = data_status, scale_data_status = scale_data_status,
    stringsAsFactors = FALSE
  )
}

assay_df <- do.call(rbind, assay_rows)
write.csv(assay_df, file.path(output_dir, "assay_info.csv"), row.names = FALSE)
message("  Saved: assay_info.csv (", nrow(assay_df), " assays)")

## =========================================================
## 5. Spatial info
## =========================================================
message("[5/9] Inspecting spatial info...")

known_coord_pairs <- list(
  c("x", "y"),
  c("imagecol", "imagerow"),
  c("pxl_col_in_fullres", "pxl_row_in_fullres"),
  c("X", "Y")
)

spatial_names <- tryCatch(Images(sobj), error = function(e) character(0))

if (length(spatial_names) == 0) {
  message("  WARNING: No spatial images found via Images().")
  spatial_df <- data.frame(
    image_name = character(0), image_class = character(0),
    spatial_key = character(0), coord_columns = character(0),
    x_min = numeric(0), x_max = numeric(0),
    y_min = numeric(0), y_max = numeric(0),
    n_spots = integer(0), stringsAsFactors = FALSE
  )
} else {
  spatial_rows <- list()
  for (img_name in spatial_names) {
    img_obj <- tryCatch(sobj@images[[img_name]], error = function(e) NULL)

    img_class <- if (!is.null(img_obj)) {
      tryCatch(class(img_obj)[1], error = function(e) "ERROR")
    } else {
      "ERROR"
    }

    # 3-strategy SpatialKey probe
    sk <- tryCatch(SpatialKey(sobj, image = img_name), error = function(e) NULL)
    if (is.null(sk) && !is.null(img_obj)) {
      sk <- tryCatch(Key(img_obj), error = function(e) NULL)
    }
    if (is.null(sk)) {
      sk <- tryCatch(Key(sobj[[default_assay]]), error = function(e) NULL)
    }
    spatial_key <- if (!is.null(sk)) sk else NA_character_

    # Coordinates + range detection
    coord_cols <- NA_character_
    n_spots <- NA_integer_
    x_min <- x_max <- y_min <- y_max <- NA_real_
    if (!is.null(img_obj)) {
      coords <- tryCatch(GetTissueCoordinates(img_obj), error = function(e) NULL)
      if (!is.null(coords)) {
        coord_cols <- paste(colnames(coords), collapse = ", ")
        n_spots <- nrow(coords)

        for (pair in known_coord_pairs) {
          x_name <- pair[1]
          y_name <- pair[2]
          if (x_name %in% colnames(coords) && y_name %in% colnames(coords)) {
            x_vals <- tryCatch(as.numeric(coords[[x_name]]), error = function(e) NULL)
            y_vals <- tryCatch(as.numeric(coords[[y_name]]), error = function(e) NULL)
            if (!is.null(x_vals) && !is.null(y_vals)) {
              x_min <- min(x_vals, na.rm = TRUE)
              x_max <- max(x_vals, na.rm = TRUE)
              y_min <- min(y_vals, na.rm = TRUE)
              y_max <- max(y_vals, na.rm = TRUE)
              break
            }
          }
        }
      }
    }

    spatial_rows[[length(spatial_rows) + 1]] <- data.frame(
      image_name = img_name, image_class = img_class,
      spatial_key = spatial_key, coord_columns = coord_cols,
      x_min = x_min, x_max = x_max, y_min = y_min, y_max = y_max,
      n_spots = n_spots, stringsAsFactors = FALSE
    )
  }
  spatial_df <- do.call(rbind, spatial_rows)
}

write.csv(spatial_df, file.path(output_dir, "spatial_info.csv"), row.names = FALSE)
message("  Saved: spatial_info.csv (", nrow(spatial_df), " images)")

## =========================================================
## 6. Metadata inspection
## =========================================================
message("[6/9] Inspecting metadata...")

md <- tryCatch(sobj@meta.data, error = function(e) NULL)

if (is.null(md)) {
  message("  WARNING: Could not access @meta.data.")
  writeLines("ERROR: meta.data not accessible", file.path(output_dir, "metadata_columns.csv"))
} else {
  col_classes <- sapply(md, function(x) class(x)[1])
  col_nunique <- sapply(md, function(x) length(unique(x)))

  skip_levels <- sapply(names(md), function(col) {
    cls <- col_classes[col]
    if (cls %in% c("factor", "logical")) return(FALSE)
    if (cls == "character" && col_nunique[col] <= 100) return(FALSE)
    return(TRUE)
  })

  md_summary <- data.frame(
    column = names(md), class = col_classes,
    n_unique = col_nunique,
    is_categorical = col_classes %in% c("factor", "character", "logical"),
    skipped_levels = skip_levels,
    stringsAsFactors = FALSE
  )
  write.csv(md_summary, file.path(output_dir, "metadata_columns.csv"), row.names = FALSE)
  message("  Saved: metadata_columns.csv (", nrow(md_summary), " columns)")

  expand_cols <- names(skip_levels)[!skip_levels]
  if (length(expand_cols) > 0) {
    levels_rows <- list()
    for (col in expand_cols) {
      tbl <- tryCatch(table(md[[col]], useNA = "ifany"), error = function(e) NULL)
      if (!is.null(tbl)) {
        for (i in seq_along(tbl)) {
          levels_rows[[length(levels_rows) + 1]] <- data.frame(
            column = col, level = names(tbl)[i], count = as.integer(tbl[i]),
            stringsAsFactors = FALSE
          )
        }
      }
    }
    if (length(levels_rows) > 0) {
      levels_df <- do.call(rbind, levels_rows)
      write.csv(levels_df, file.path(output_dir, "metadata_levels.csv"), row.names = FALSE)
      message("  Saved: metadata_levels.csv (", nrow(levels_df), " level entries)")
    } else {
      writeLines("No level entries generated", file.path(output_dir, "metadata_levels.csv"))
    }
  } else {
    writeLines("No columns eligible for level expansion", file.path(output_dir, "metadata_levels.csv"))
  }
}

## =========================================================
## 7. Reductions
## =========================================================
message("[7/9] Inspecting reductions...")

red_names <- tryCatch(names(sobj@reductions), error = function(e) character(0))

if (length(red_names) == 0) {
  message("  No reductions found.")
  writeLines("No reductions found", file.path(output_dir, "reductions.csv"))
} else {
  red_rows <- list()
  for (rname in red_names) {
    r <- tryCatch(sobj@reductions[[rname]], error = function(e) NULL)
    if (!is.null(r)) {
      r_class <- tryCatch(class(r)[1], error = function(e) "ERROR")
      r_dims  <- tryCatch(dim(Embeddings(r)), error = function(e) c(NA_integer_, NA_integer_))
      red_rows[[length(red_rows) + 1]] <- data.frame(
        reduction = rname, class = r_class,
        n_cells = r_dims[1], n_dimensions = r_dims[2],
        stringsAsFactors = FALSE
      )
    }
  }
  red_df <- do.call(rbind, red_rows)
  write.csv(red_df, file.path(output_dir, "reductions.csv"), row.names = FALSE)
  message("  Saved: reductions.csv (", nrow(red_df), " reductions)")
}

## =========================================================
## 8. Inspection summary + provenance
## =========================================================
message("[8/9] Writing summary and provenance...")

summary_lines <- c(
  paste("=== Seurat", rds_label, "Inspection Summary ==="),
  paste("Date:", Sys.time()),
  "",
  paste("Seurat version:", seurat_vstr),
  paste("Class:", obj_class),
  paste("Dimensions:", n_genes, "genes x", n_cells, "spots"),
  paste("Default assay:", default_assay),
  paste("All assays:", paste(all_assays, collapse = ", ")),
  paste("Spatial images:", if (length(spatial_names) > 0) paste(spatial_names, collapse = ", ") else "NONE"),
  paste("Reductions:", if (length(red_names) > 0) paste(red_names, collapse = ", ") else "NONE"),
  "",
  "--- str() output (max.level=1, list.len=20) ---",
  str_lines,
  "",
  "--- slotNames() ---",
  paste(slot_names, collapse = ", ")
)
writeLines(summary_lines, file.path(output_dir, "inspection_summary.txt"))

md5_status <- if (!is.na(expected_md5)) expected_md5 else "NOT_IN_METADATA"

get_meta <- function(field) {
  if (is.null(meta) || is.null(meta[[field]])) "NOT_FOUND" else as.character(meta[[field]])
}

provenance_lines <- c(
  paste("=== Provenance:", rds_filename, "==="),
  paste("File path:", rds_file),
  paste("File size:", fi$size, "bytes (", round(fi$size / 1e6, 1), "MB)"),
  paste("File mtime:", fi$mtime),
  paste("Expected MD5:", md5_status),
  paste("Actual MD5: NOT COMPUTED"),
  "",
  paste("Dataset:", get_meta("dataset")),
  paste("Source:", get_meta("source")),
  paste("Zenodo DOI:", get_meta("doi")),
  paste("Paper DOI:", get_meta("paper_doi")),
  paste("GitHub:", get_meta("github")),
  paste("Species:", get_meta("species")),
  paste("Technology:", get_meta("technology")),
  paste("Description:", get_meta("description")),
  paste("Downloaded date:", get_meta("downloaded_date")),
  "",
  paste("Seurat version loaded:", seurat_vstr),
  paste("R version:", R.version.string),
  paste("Loaded from:", rds_file),
  paste("Inspection date:", Sys.time())
)
writeLines(provenance_lines, file.path(output_dir, "provenance.txt"))

## =========================================================
## 9. Cleanup
## =========================================================
message("[9/9] Cleaning up...")
rm(sobj)
gc_out <- gc(verbose = FALSE)

message("\n=== Spatial s01 Inspection Complete (", object_key, ") ===")
message("Output directory: ", output_dir)
message("Files written:")
for (f in list.files(output_dir)) {
  message("  ", f)
}
message("\nEnd time: ", Sys.time())
