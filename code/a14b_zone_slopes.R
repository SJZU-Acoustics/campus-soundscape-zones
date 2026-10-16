# A14b — zone-specific slopes (the common-slope assumption is rejected: zone × level LRT p 0.001 / 0.03): per-zone level slope,
# and the neutral level only for zones whose slope is negative and significant; cluster bootstrap over recordings (1,000) on the B&K layer.
source(file.path("code", "helpers.R"))
suppressPackageStartupMessages(library(emmeans))
OUT <- file.path("output", "analysis_14_neutral_levels"); w <- function(x, f) write_csv(x, file.path(OUT, f))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
REF <- 55; d <- load_p35("wav") %>% mutate(L = bk_laeq - REF)
rng <- d %>% distinct(env_id, zone, L) %>% group_by(zone) %>% summarise(n_env = n(), lo = min(L) + REF, hi = max(L) + REF, median_level = median(L) + REF, .groups = "drop")
thr_zone <- function(m, neutral) { f <- lme4::fixef(m); int <- f["(Intercept)"] + c(0, f[paste0("zone", ZONE_LEVELS[-1])]); sl <- f["L"] + c(0, f[paste0("zone", ZONE_LEVELS[-1], ":L")]); list(thr = REF + (neutral - int) / sl, slope = sl) }
out <- list()
for (y in c("rating_appropriateness", "comfort")) {
  neutral <- if (y == "comfort") 0 else 3
  m <- fit_backbone("zone * L", y, d); tr <- as_tibble(emtrends(m, "zone", var = "L")); est <- thr_zone(m, neutral)
  set.seed(141); envs <- unique(d$env_id)
  boot <- replicate(1000, { s <- sample(envs, replace = TRUE); db <- resample_stops(d, s)
    thr_zone(checked_resample_lmer(as.formula(paste0(y, " ~ zone * L + (1 | evaluator) + (1 | env_id)")), data = db), neutral)$thr })
  out[[y]] <- tibble(outcome = y, zone = ZONE_LEVELS, slope_per_dB = tr$L.trend, slope_se = tr$SE, slope_lo = tr$lower.CL, slope_hi = tr$upper.CL,
                     neutral_level = est$thr, ci_lo = apply(boot, 1, quantile, 0.025), ci_hi = apply(boot, 1, quantile, 0.975)) %>% left_join(rng, by = "zone") %>%
    mutate(slope_significant_negative = slope_hi < 0, headroom_dB = neutral_level - median_level, within_observed_range = neutral_level >= lo & neutral_level <= hi,
           reportable = slope_significant_negative & within_observed_range)
}
out <- bind_rows(out); w(out, "neutral_levels_zone_slopes.csv")
print(out %>% mutate(across(c(neutral_level, ci_lo, ci_hi, lo, hi, median_level, headroom_dB), ~round(.x, 1)), across(c(slope_per_dB, slope_se, slope_lo, slope_hi), ~round(.x, 4))) %>%
        select(outcome, zone, slope_per_dB, slope_lo, slope_hi, neutral_level, ci_lo, ci_hi, lo, hi, median_level, headroom_dB, reportable), n = 16, width = 220)
