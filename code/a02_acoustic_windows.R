# A02 — the two acoustic windows (XL2 full measurement vs B&K 30-s region) and the
# structure of the psychoacoustic indicator set.
source(file.path("code", "helpers.R"))
suppressPackageStartupMessages(library(performance))
OUT <- file.path("output", "analysis_02_acoustic_windows"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
w <- function(x, f) write_csv(x, file.path(OUT, f))
d <- load_p35("wav"); e <- load_p35_env() %>% filter(wav_recorded == 1)

# ---- 1. agreement of the two windows (143 recordings)
e <- e %>% mutate(diff = bk_laeq - xl2_laeq, zone = factor(ZONE_LEVELS[site], levels = ZONE_LEVELS))
fit <- lm(bk_laeq ~ xl2_laeq, data = e)
cat("bias (30-s − full): mean", round(mean(e$diff), 2), "sd", round(sd(e$diff), 2), "median", round(median(e$diff), 2), "| r =", round(cor(e$bk_laeq, e$xl2_laeq), 3),
    "| slope", round(coef(fit)[2], 3), "| |diff|>5 dB:", sum(abs(e$diff) > 5), "of", nrow(e), "\n")
by_zone <- e %>% group_by(zone) %>% summarise(n = n(), bias = mean(diff), sd = sd(diff), r = cor(bk_laeq, xl2_laeq), .groups = "drop")
by_sess <- e %>% group_by(session) %>% summarise(n = n(), bias = mean(diff), sd = sd(diff), .groups = "drop")
w(by_zone, "agreement_by_zone.csv"); w(by_sess, "agreement_by_session.csv"); print(by_zone %>% mutate(across(where(is.numeric), ~round(.x, 2))))
# is the bias related to the level or to transients (L5-L95 in the full window)?
cat("r(diff, xl2_laeq) =", round(cor(e$diff, e$xl2_laeq), 3), " r(diff, xl2_dl5_95) =", round(cor(e$diff, e$xl2_dl5_95), 3), " r(diff, bk_dl5_95) =", round(cor(e$diff, e$bk_dl5_95), 3), "\n")

# ---- 2. which window predicts the ratings (666 records, crossed RE)
res <- list()
for (y in c("comfort", "richness", "rating_noisiness", "rating_appropriateness")) {
  m_bk <- fit_backbone("bk_laeq_z", y, d); m_xl <- fit_backbone("xl2_laeq_z", y, d); m_both <- fit_backbone("bk_laeq_z + xl2_laeq_z", y, d)
  cb <- summary(m_bk)$coefficients["bk_laeq_z", ]; cx <- summary(m_xl)$coefficients["xl2_laeq_z", ]; bb <- summary(m_both)$coefficients
  res[[y]] <- tibble(outcome = y,
                     beta_bk = cb[1], se_bk = cb[2], p_bk = cb[5], r2m_bk = r2_nakagawa(m_bk)$R2_marginal, aic_bk = AIC(m_bk),
                     beta_xl2 = cx[1], se_xl2 = cx[2], p_xl2 = cx[5], r2m_xl2 = r2_nakagawa(m_xl)$R2_marginal, aic_xl2 = AIC(m_xl),
                     both_beta_bk = bb["bk_laeq_z", 1], both_p_bk = bb["bk_laeq_z", 5], both_beta_xl2 = bb["xl2_laeq_z", 1], both_p_xl2 = bb["xl2_laeq_z", 5])
}
res <- bind_rows(res); w(res, "window_predictive_comparison.csv")
print(res %>% mutate(across(where(is.numeric), ~signif(.x, 3))), width = 220)

# ---- 3. indicator structure (143 recordings)
ind <- c("bk_laeq", "bk_l5", "bk_l95", "bk_dl5_95", "n", "n5", "n5_n95", "s", "r", "f", "tonality")
cm <- cor(e[, ind]); w(as.data.frame(round(cm, 3)) %>% mutate(indicator = ind) %>% relocate(indicator), "indicator_correlations.csv")
cat("\nr with bk_laeq:", paste(ind[-1], round(cm[1, -1], 2), sep = "=", collapse = " "), "\n")
pc <- prcomp(e[, ind], scale. = TRUE); ve <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)
cat("indicator PCA variance %:", ve[1:5], "\n")
load <- as.data.frame(round(pc$rotation[, 1:3], 2)) %>% mutate(indicator = ind) %>% relocate(indicator); w(load, "indicator_pca_loadings.csv"); print(load)
vif <- sapply(ind, function(v) { others <- setdiff(ind, v); 1 / (1 - summary(lm(reformulate(others, v), data = e))$r.squared) })
cat("VIF (all 11 together):", paste(ind, round(vif, 1), sep = "=", collapse = " "), "\n")
# env-level mean ratings vs indicators, and partial r given bk_laeq
em <- d %>% group_by(env_id) %>% summarise(comfort = mean(comfort), richness = mean(richness), noisiness = mean(rating_noisiness), appropriateness = mean(rating_appropriateness), .groups = "drop") %>%
  inner_join(e %>% select(env_id, all_of(ind)), by = "env_id")
pr <- function(y, x, z) { ry <- resid(lm(em[[y]] ~ em[[z]])); rx <- resid(lm(em[[x]] ~ em[[z]])); cor(ry, rx) }
tab <- expand_grid(outcome = c("comfort", "richness", "noisiness", "appropriateness"), indicator = ind) %>% rowwise() %>%
  mutate(r_env = cor(em[[outcome]], em[[indicator]]), partial_r_given_laeq = if (indicator == "bk_laeq") NA_real_ else pr(outcome, indicator, "bk_laeq")) %>% ungroup()
w(tab, "env_level_indicator_correlations.csv")
print(tab %>% mutate(across(where(is.numeric), ~round(.x, 2))) %>% pivot_wider(names_from = outcome, values_from = c(r_env, partial_r_given_laeq)), width = 220)
