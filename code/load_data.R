# =============================================================================
# load_data.R — single data entry for the campus-soundscape-zones release.
#
# Data source: the Zenodo workbook
#   data/Campus_soundwalk_summer_2026_perception_acoustics_visual.xlsx  (CC BY 4.0)
#
# Each sheet is read with its numeric cells as numbers, written once to a CSV
# cache and read back with readr, so every analysis receives the column types
# that readr infers from a CSV, and both "NA" and "" map to NA.
#
# code/helpers.R and three analysis modules read their tables through
# read_frozen_csv(), which names the table's file in the analysis pipeline; the
# mapping below resolves each name to its workbook sheet. Nothing else in the
# repository touches the workbook.
# =============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(readr)
})

XLSX_PATH <- file.path("data", "Campus_soundwalk_summer_2026_perception_acoustics_visual.xlsx")
CACHE_DIR <- file.path("output", "_workbook_cache")

# pipeline file name -> workbook sheet
SHEET_FOR_FILE <- c(
  "p35_records.csv"       = "records",
  "p35_people.csv"        = "people",
  "p35_env.csv"           = "environments",
  "p35_sessions.csv"      = "sessions",
  "p35_sites.csv"         = "sites",
  "indices_by_env.csv"    = "ecoacoustic_indices",
  "detections_by_env.csv" = "detected_sources"
)

read_workbook_sheet <- function(sheet, guess_max) {
  if (!file.exists(XLSX_PATH)) {
    stop("Workbook not found at ", XLSX_PATH, "\n",
         "Download it from Zenodo and place it in data/ — see data/README.md.",
         call. = FALSE)
  }
  dir.create(CACHE_DIR, recursive = TRUE, showWarnings = FALSE)
  cached <- file.path(CACHE_DIR, paste0(sheet, ".csv"))
  if (!file.exists(cached)) {
    raw <- readxl::read_excel(XLSX_PATH, sheet = sheet, guess_max = 1e5)
    readr::write_csv(raw, cached, na = "")
  }
  readr::read_csv(cached, show_col_types = FALSE, progress = FALSE, guess_max = guess_max)
}

#' Read a deposited table by its file name in the analysis pipeline.
read_frozen_csv <- function(filename, guess_max = 1000) {
  sheet <- SHEET_FOR_FILE[[filename]]
  if (is.null(sheet)) {
    stop("No workbook sheet is registered for ", filename, call. = FALSE)
  }
  read_workbook_sheet(sheet, guess_max)
}
