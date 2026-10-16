# =============================================================================
# run_all.R — reproduce every figure and table of the campus soundscape paper
# from the Zenodo workbook.
#
#   Rscript run_all.R
#
# Reads   data/Campus_soundwalk_summer_2026_perception_acoustics_visual.xlsx
#         data/osm/osm_campus_features.csv  (optional; Figure 1 only)
# Writes  output/<analysis>/   one folder per analysis module
#         output/figures/      Figures 1-6 (PNG, 600 dpi)
#         output/tables/       Tables 1-2 (LaTeX fragments)
#         output/tables/si/    Supplementary Tables S1-S34 (LaTeX fragments)
# =============================================================================

start_time <- Sys.time()

if (!file.exists(file.path("code", "helpers.R"))) {
  stop("Run this script from the repository root: Rscript run_all.R", call. = FALSE)
}

# Analysis modules. Each reads the deposited tables and writes its own result
# tables; a06b also reads the result table of a06, so it runs after it. The
# order is otherwise for readability only.
MODULES <- c(
  "a01_inventory",
  "a02_acoustic_windows",
  "a03_zone_profiles",
  "a04_acoustic",
  "a05_sources",
  "a06_visual",
  "a06b_without_w04",
  "a07_daytype_slot",
  "a08_person_level",
  "a09_variance",
  "a09b_random_slope_permutation_all",
  "a10b_zone_by_comfort",
  "a11_feature_families",
  "a12_robustness",
  "a13_compensation",
  "a14b_zone_slopes",
  "a15_richness_diversity",
  "a16_quadrants",
  "a16b_comfort_by_quadrant",
  "a17_models",
  "a18_models",
  "a19_source_structure",
  "a20_core_stability",
  "a20_storyline_checks",
  "a21_cowalker"
)

# Each script runs in a fresh R session, as in the working pipeline, so no
# object, option or random-number state passes between them.
rscript <- file.path(R.home("bin"), "Rscript")
run <- function(script) {
  t0 <- Sys.time()
  status <- system2(rscript, file.path("code", paste0(script, ".R")))
  if (status != 0) stop(script, " failed (exit status ", status, ")", call. = FALSE)
  message("    ", script, " done in ",
          round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), " s")
}

dir.create("output", showWarnings = FALSE)
for (module in MODULES) {
  message("\n=== ", module, " ===")
  run(module)
}

message("\n=== display items ===")
run("build_figures")
run("build_tables")

message("\nComplete in ",
        round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 1),
        " min. Figures are in output/figures and table fragments in output/tables.")
