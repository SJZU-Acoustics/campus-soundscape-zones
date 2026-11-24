# A24 — same-session listeners needed to estimate a place mean.
source(file.path("code", "helpers.R"))
A09 <- file.path("output", "analysis_09_variance_partitioning")
OUT <- file.path("output", "analysis_24_listeners_per_point")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
components <- read_csv(file.path(A09, "recording_mean_reliability.csv"), show_col_types = FALSE)
reliability <- function(place, session, independent, k) {
  place / (place + session + independent / k)
}
listeners_needed <- function(place, session, independent, target) {
  allowance <- place * (1 / target - 1) - session
  ifelse(allowance > 0, pmax(1, ceiling(independent / allowance)), NA_real_)
}
stopifnot(with(components,
  all(abs(reliability(icc1_place, session_share, independent_share, k_raters_per_recording) -
          icc_k_place) < 1e-9)))
proj <- components %>% select(outcome, icc1_place, session_share, independent_share)
for (k in c(3, 5, 8, 10)) {
  proj[[paste0("rel_", k)]] <- with(proj, reliability(icc1_place, session_share, independent_share, k))
}
for (target in c(0.70, 0.80, 0.90)) {
  proj[[sprintf("k_for_%.2f", target)]] <- with(proj, listeners_needed(icc1_place, session_share, independent_share, target))
}
write_csv(proj, file.path(OUT, "listeners_per_point.csv"))
