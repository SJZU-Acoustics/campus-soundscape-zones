# A24 — listeners needed per stop. Spearman-Brown projection from the four-way variance partition of A09
# (990 ratings): the reliability of a stop mean of k listeners, with site and stop variance as signal and
# evaluator, session and residual variance as error. Runs after A09.
source(file.path("code", "helpers.R"))
A09 <- file.path("output", "analysis_09_variance_partitioning")
OUT <- file.path("output", "analysis_24_listeners_per_point")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
vp <- read_csv(file.path(A09, "variance_partition_990.csv"), show_col_types = FALSE)
icc <- vp %>% group_by(outcome) %>% summarise(icc1_place = sum(share[grp %in% c("site", "env_id")]), .groups = "drop")
rel_k <- function(icc, k) k * icc / (1 + (k - 1) * icc)
k_for <- function(icc, target) ceiling(target * (1 - icc) / (icc * (1 - target)))
proj <- icc %>% mutate(rel_3 = rel_k(icc1_place, 3), rel_5 = rel_k(icc1_place, 5), rel_8 = rel_k(icc1_place, 8), rel_10 = rel_k(icc1_place, 10),
                       k_for_0.70 = k_for(icc1_place, 0.70), k_for_0.80 = k_for(icc1_place, 0.80), k_for_0.90 = k_for(icc1_place, 0.90))
# check against the reliability that A09 reports at the observed mean number of raters per stop
chk <- read_csv(file.path(A09, "recording_mean_reliability.csv"), show_col_types = FALSE)
stopifnot(all(abs(rel_k(chk$icc1_place, chk$k_raters_per_recording) - chk$icc_k_place) < 1e-9),
          all(abs(icc$icc1_place[match(chk$outcome, icc$outcome)] - chk$icc1_place) < 1e-9))
write_csv(proj, file.path(OUT, "listeners_per_point.csv"))
