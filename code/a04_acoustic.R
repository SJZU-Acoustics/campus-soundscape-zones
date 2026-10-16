# A04 — acoustic indicators → the four outcomes (666 WAV records; crossed evaluator + recording RE).
# Families (A02): one level indicator (bk_laeq; N5 alternate) plus dynamics/spectral/temporal indicators.
source(file.path("code", "helpers.R"))
suppressPackageStartupMessages(library(performance))
OUT <- file.path("output", "analysis_04_acoustic_to_perception"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
w <- function(x, f) write_csv(x, file.path(OUT, f))
d <- load_p35("wav"); dall <- load_p35("all")
OUTS <- c("comfort", "richness", "rating_noisiness", "rating_appropriateness")
IND <- c("bk_laeq", "n5", "bk_dl5_95", "n5_n95", "s", "r", "f", "tonality")

# ---- 1. single-indicator screen (family = 4 outcomes × 8 indicators = 32)
scan <- expand_grid(outcome = OUTS, indicator = IND) %>% rowwise() %>% mutate(
  fit = list(fit_backbone(paste0(indicator, "_z"), outcome, d)),
  beta = fixef(fit)[2], se = summary(fit)$coefficients[2, 2], p = summary(fit)$coefficients[2, 5],
  r2m = r2_nakagawa(fit)$R2_marginal) %>% ungroup() %>% select(-fit) %>% bh_family(p_col = "p")
w(scan, "single_indicator_screen.csv")
print(scan %>% mutate(across(c(beta, se, r2m), ~round(.x, 3)), p = signif(p, 2), q = signif(q, 2)) %>% arrange(outcome, q), n = 40, width = 150)

# ---- 2. beyond level: each non-level indicator added to LAeq (family = 4 × 6 = 24)
beyond <- expand_grid(outcome = OUTS, indicator = setdiff(IND, c("bk_laeq", "n5"))) %>% rowwise() %>% mutate(
  fit = list(fit_backbone(paste0("bk_laeq_z + ", indicator, "_z"), outcome, d)),
  beta_ind = fixef(fit)[3], p_ind = summary(fit)$coefficients[3, 5], beta_laeq = fixef(fit)[2],
  r2m = r2_nakagawa(fit)$R2_marginal, r2m_laeq_only = r2_nakagawa(fit_backbone("bk_laeq_z", outcome, d))$R2_marginal) %>% ungroup() %>% select(-fit) %>%
  mutate(delta_r2m = r2m - r2m_laeq_only) %>% bh_family(p_col = "p_ind")
w(beyond, "beyond_level_screen.csv")
print(beyond %>% mutate(across(c(beta_ind, beta_laeq, r2m, delta_r2m), ~round(.x, 3)), p_ind = signif(p_ind, 2), q = signif(q, 2)) %>% arrange(outcome, q), n = 30, width = 150)

# ---- 3. multi-indicator model: level + tonality + sharpness (+ dynamics); leave-one-site-out generality of the slopes
multi <- list()
for (y in OUTS) {
  m <- fit_backbone("bk_laeq_z + tonality_z + s_z + bk_dl5_95_z", y, d)
  co <- summary(m)$coefficients
  multi[[y]] <- as_tibble(co, rownames = "term") %>% mutate(outcome = y, r2m = r2_nakagawa(m)$R2_marginal, r2c = r2_nakagawa(m)$R2_conditional)
}
multi <- bind_rows(multi); w(multi, "multi_indicator_models.csv")
print(multi %>% filter(term != "(Intercept)") %>% mutate(across(c(Estimate, `Std. Error`, r2m, r2c), ~round(.x, 3)), p = signif(`Pr(>|t|)`, 2)) %>% select(outcome, term, Estimate, `Std. Error`, p, r2m, r2c), n = 20)

loso <- expand_grid(outcome = c("comfort", "rating_noisiness"), left_out = levels(d$site)) %>% rowwise() %>% mutate(
  fit = list(fit_backbone("bk_laeq_z + tonality_z + s_z", outcome, filter(d, site != left_out))),
  b_laeq = fixef(fit)["bk_laeq_z"], b_ton = fixef(fit)["tonality_z"], b_s = fixef(fit)["s_z"]) %>% ungroup() %>% select(-fit)
w(loso, "loso_slopes.csv")
print(loso %>% group_by(outcome) %>% summarise(across(c(b_laeq, b_ton, b_s), list(min = min, max = max)), .groups = "drop") %>% mutate(across(where(is.numeric), ~round(.x, 3))), width = 150)

# ---- 4. whole-sample check with the XL2 level (990 records) and the level slope on the 2 July-only rows
xl <- tibble(outcome = OUTS) %>% rowwise() %>% mutate(
  fit = list(fit_backbone("xl2_laeq_z", outcome, dall)), beta = fixef(fit)[2], p = summary(fit)$coefficients[2, 5], r2m = r2_nakagawa(fit)$R2_marginal,
  fit02 = list(fit_backbone("xl2_laeq_z", outcome, filter(dall, wav_recorded == 0))), beta_0702 = fixef(fit02)[2], p_0702 = summary(fit02)$coefficients[2, 5]) %>% ungroup() %>% select(-fit, -fit02)
w(xl, "xl2_whole_sample.csv"); print(xl %>% mutate(across(c(beta, r2m, beta_0702), ~round(.x, 3)), p = signif(p, 2), p_0702 = signif(p_0702, 2)))

# ---- 5. shape: is the level effect linear? quadratic term and a 3-band comparison
for (y in c("comfort", "rating_noisiness")) {
  mq <- fit_backbone("bk_laeq_z + I(bk_laeq_z^2)", y, d)
  cat(y, ": quadratic term beta", round(fixef(mq)[3], 3), "p", signif(summary(mq)$coefficients[3, 5], 2), "\n")
}
