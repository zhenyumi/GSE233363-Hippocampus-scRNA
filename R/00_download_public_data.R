## =========================================================
## Script: 00_download_public_data.R
## Purpose: Download official author-provided RDS files from Zenodo
## Input: Zenodo (DOI: 10.5281/zenodo.13893852)
## Output: data/raw/GSE233363_official/*.rds
## =========================================================

suppressPackageStartupMessages({
  library(httr)
  library(tools)
})

ZENODO_BASE <- "https://zenodo.org/api/records/13893852/files"

FILES <- data.frame(
  filename = c(
    "seurat_Chromium_All.rds",
    "seurat_Visium_DG_All.rds",
    "seurat_Visium_Hippo_All.rds",
    "seurat_Visium_WholeBrain.rds"
  ),
  size = c(
    "1.6 GB",
    "424 MB",
    "1.2 GB",
    "4.6 GB"
  ),
  md5 = c(
    "6b1f5ece1be8d119df3300e6fe280633",
    "518e4d786cf816dda10ef09886ee5821",
    "b87a53ea521a546e22a91bb2fdaff041",
    "abfd050382cf8b61af3d37100af2e0e3"
  ),
  stringsAsFactors = FALSE
)

project_dir <- getwd()
output_dir  <- file.path(project_dir, "data", "raw", "GSE233363_official")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

verify_md5 <- function(filepath, expected_md5) {
  if (!file.exists(filepath)) return(FALSE)
  actual <- tolower(md5sum(filepath))
  expected <- tolower(expected_md5)
  if (is.na(actual)) return(FALSE)
  actual == expected
}

download_file <- function(url, destfile) {
  message("  Downloading to: ", destfile)
  response <- GET(url, write_disk(destfile, overwrite = TRUE), progress())
  return(response$status_code)
}

message("=== Download Official GSE233363 Author-Provided RDS Files ===")
message("Start time: ", Sys.time())
message("Source: Zenodo DOI 10.5281/zenodo.13893852")
message("Target: ", output_dir)
message("")

total_size <- sum(c(1728987813, 445053868, 1302821630, 4890330884))
message("Total download size: ~", round(total_size / 1e9, 1), " GB")
message("Estimated time: varies by network speed (10 min - 2 hours)")
message("")

for (i in seq_len(nrow(FILES))) {
  fname <- FILES$filename[i]
  fsize <- FILES$size[i]
  fmd5  <- FILES$md5[i]
  dest  <- file.path(output_dir, fname)

  cat(sprintf("[%d/%d] %s (%s)\n", i, nrow(FILES), fname, fsize))

  if (verify_md5(dest, fmd5)) {
    message("  Already downloaded and verified (MD5 matches). Skipping.")
    next
  }

  if (file.exists(dest)) {
    message("  File exists but MD5 mismatch. Re-downloading...")
    file.remove(dest)
  }

  file_url <- paste0(ZENODO_BASE, "/", fname, "/content")
  status <- download_file(file_url, dest)

  if (status == 200 && verify_md5(dest, fmd5)) {
    message("  Download complete and verified (MD5 match).")
  } else if (status == 200) {
    warning("  Download complete but MD5 mismatch! Expected: ", fmd5)
  } else {
    warning("  Download failed with status: ", status)
  }

  gc()
  message("")
}

message("=== Download Complete ===")
message("End time: ", Sys.time())

existing <- FILES$filename[file.exists(file.path(output_dir, FILES$filename))]
message("Files present: ", length(existing), "/", nrow(FILES))
if (length(existing) > 0) {
  for (f in existing) {
    fi <- file.info(file.path(output_dir, f))
    message(sprintf("  %s (%.1f MB)", f, fi$size / 1e6))
  }
}

missing_files <- setdiff(FILES$filename, existing)
if (length(missing_files) > 0) {
  warning("Missing files: ", paste(missing_files, collapse = ", "))
}
