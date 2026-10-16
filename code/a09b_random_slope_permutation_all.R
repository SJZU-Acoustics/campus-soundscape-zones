# A09b — within-evaluator permutation test of the level random slope for all four outcomes: the same model pair as
# A09 (ML; uncorrelated slope, `||`), seed 35, Comfort and noisiness first so their nulls reproduce A09 exactly, then
# Richness and appropriateness. Slope descriptives are taken from the same ML random-slope fit.
# One pre-specified test per outcome, reported with its own permutation p.
source(file.path("code", "helpers.R"))
OUT <- file.path("output", "analysis_09_variance_partitioning")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
d <- load_p35("wav")
set.seed(35)
NPERM <- 1000
res <- list()
for (y in c("comfort", "rating_noisiness", "richness", "rating_appropriateness")) {
  f0 <- as.formula(paste0(y, " ~ bk_laeq_z + (1 | evaluator) + (1 | env_id)"))
  f1 <- as.formula(paste0(y, " ~ bk_laeq_z + (1 + bk_laeq_z || evaluator) + (1 | env_id)"))
  lrt_of <- function(dd) { m0 <- checked_resample_lmer(f0, dd, REML = FALSE); m1 <- checked_resample_lmer(f1, dd, REML = FALSE); as.numeric(2 * (logLik(m1) - logLik(m0))) }
  obs <- lrt_of(d)
  null <- vapply(seq_len(NPERM), function(i) { dp <- d %>% group_by(evaluator) %>% mutate(bk_laeq_z = sample(bk_laeq_z)) %>% ungroup(); lrt_of(dp) }, numeric(1))
  p_perm <- (sum(null >= obs) + 1) / (NPERM + 1)
  m1 <- checked_resample_lmer(f1, d, REML = FALSE)
  vc <- as.data.frame(VarCorr(m1)) %>% filter(grepl("^evaluator", grp), var1 == "bk_laeq_z", is.na(var2))
  sl <- fixef(m1)[["bk_laeq_z"]] + ranef(m1)$evaluator[["bk_laeq_z"]]
  res[[y]] <- tibble(outcome = y, mean_slope = fixef(m1)[["bk_laeq_z"]], sd_slope = sqrt(vc$vcov[1]),
                     blup_slope_min = min(sl), blup_slope_max = max(sl), n_positive_slopes = sum(sl > 0), n_evaluators = length(sl),
                     lrt_observed = obs, null_q95 = quantile(null, 0.95), null_q99 = quantile(null, 0.99), p_permutation = p_perm, n_perm = NPERM)
  write_csv(tibble(null_lrt = null), file.path(OUT, paste0("perm_null_all_", y, ".csv")))
  cat(y, ": LRT", round(obs, 2), "| null 95%", round(quantile(null, 0.95), 2), "| p_perm", signif(p_perm, 3), "| slope", round(fixef(m1)[["bk_laeq_z"]], 3), "SD", round(sqrt(vc$vcov[1]), 3), "| positive", sum(sl > 0), "\n")
}
res <- bind_rows(res); write_csv(res, file.path(OUT, "random_slope_permutation_all.csv"))
cat("done\n")
