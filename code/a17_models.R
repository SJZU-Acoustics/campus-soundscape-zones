# A17 — ecoacoustic indices → perceived richness (and comfort, noisiness, appropriateness).
# Indices per recording are the deposited ecoacoustic_indices table (scikit-maad; 30-s region window primary, full recording
# companion), joined to the record table by env_id. Predictions, stated before any audio was analysed: ACI, H, ADI and BI
# correlate positively with richness after LAeq adjustment; NDSI positively with comfort/appropriateness and negatively with noisiness;
# the indices lift richness LOSO recording-mean R² above 0.2; they add little to comfort beyond level + tonality + source groups.
source(file.path("code", "helpers.R"))
suppressPackageStartupMessages({library(performance)})
OUT <- file.path("output", "analysis_17_ecoacoustic_indices"); w <- function(x, f) write_csv(x, file.path(OUT, f))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
idx <- read_frozen_csv("indices_by_env.csv")
PRIMARY <- c("ACI", "ADI", "AEI", "BI", "NDSI", "Ht", "Hf", "H_total", "NDSI_kasten", "BI_2_8k", "ADI_10k", "AEI_10k")
front <- c("env_id", "window", "duration_s", "region_start_s")
ALL_IDX <- setdiff(names(idx), c(front, "LAeq_recomputed", "LZeq_recomputed"))
SECONDARY <- setdiff(ALL_IDX, PRIMARY)
reg <- idx %>% filter(window == "region") %>% select(-window, -duration_s, -region_start_s)
ful <- idx %>% filter(window == "full") %>% select(-window, -duration_s, -region_start_s)
# drop indices that are constant or NA on the region window
keep <- ALL_IDX[sapply(ALL_IDX, function(v) { x <- reg[[v]]; sum(!is.na(x)) >= 130 && sd(x, na.rm = TRUE) > 0 })]
dropped <- setdiff(ALL_IDX, keep); cat("indices dropped (NA/constant):", paste(dropped, collapse = ", "), "\n")
SECONDARY <- intersect(SECONDARY, keep); PRIMARY <- intersect(PRIMARY, keep)
d <- load_p35("wav") %>% mutate(env_id = as.character(env_id)) %>% inner_join(reg, by = "env_id") %>% mutate(env_id = factor(env_id))
stopifnot(nrow(d) == 666)
src_cols <- paste0("src_", SOURCES)
d <- d %>% mutate(total_sal = rowSums(across(all_of(src_cols))), n_sources = rowSums(across(all_of(src_cols)) > 0))
for (v in keep) d[[paste0(v, "_z")]] <- as.numeric(scale(d[[v]]))
OUTS <- c("richness", "comfort", "rating_noisiness", "rating_appropriateness")

# ---- 1. descriptives, and agreement between the two windows (recording level, 143)
desc <- reg %>% summarise(across(all_of(PRIMARY), list(mean = ~mean(.x, na.rm = TRUE), sd = ~sd(.x, na.rm = TRUE), min = ~min(.x, na.rm = TRUE), max = ~max(.x, na.rm = TRUE)))) %>%
  pivot_longer(everything(), names_to = c("index", "stat"), names_pattern = "(.*)_(mean|sd|min|max)$") %>% pivot_wider(names_from = stat, values_from = value)
agree <- tibble(index = keep) %>% rowwise() %>% mutate(r_region_full = cor(reg[[index]], ful[[index]], use = "pairwise", method = "pearson"),
  rho_region_full = cor(reg[[index]], ful[[index]], use = "pairwise", method = "spearman"), median_region_minus_full = median(reg[[index]] - ful[[index]], na.rm = TRUE)) %>% ungroup()
w(desc, "descriptives_primary.csv"); w(agree, "window_agreement.csv")
cat("\n== window agreement (primary) ==\n"); print(agree %>% filter(index %in% PRIMARY) %>% mutate(across(where(is.numeric), ~round(.x, 3))), n = 20)

# ---- 2. correlation structure: primary indices with each other and with the frozen acoustic layer (recording level)
envl <- d %>% group_by(env_id) %>% summarise(across(c(all_of(keep), bk_laeq, tonality, s, n, bk_dl5_95, src_natural, src_human, src_mechanical, total_sal, n_sources, richness, comfort), ~mean(.x, na.rm = TRUE)), .groups = "drop")
cm <- cor(envl[, c(PRIMARY, "bk_laeq", "tonality", "s", "bk_dl5_95")], use = "pairwise", method = "spearman")
w(as.data.frame(round(cm, 2)) %>% mutate(var = rownames(cm)) %>% relocate(var), "index_correlations_env.csv")
cat("\n== Spearman with LAeq / tonality / sharpness (recording level) ==\n")
print(as_tibble(round(cm[PRIMARY, c("bk_laeq", "tonality", "s", "bk_dl5_95")], 2)) %>% mutate(index = PRIMARY) %>% relocate(index), n = 20)
pc <- prcomp(scale(envl[, PRIMARY])); ve <- summary(pc)$importance[2, 1:4]
cat("PCA of the 12 primary indices (143 recordings): variance PC1-4", round(ve, 2), "\n")
load <- as_tibble(round(pc$rotation[, 1:3], 2)) %>% mutate(index = PRIMARY) %>% relocate(index); w(load, "pca_loadings_primary.csv"); print(load, n = 20)

# ---- 3. single-index screen, level adjusted, crossed RE: primary family (12 × 4 = 48) and secondary family (k × 4), BH within family
screen <- function(vars, fam) expand_grid(outcome = OUTS, index = vars) %>% rowwise() %>% mutate(
  fit = list(fit_backbone(paste0("bk_laeq_z + ", index, "_z"), outcome, d)), n = nrow(d), beta = fixef(fit)[3], se = summary(fit)$coefficients[3, 2],
  p = summary(fit)$coefficients[3, 5], r2m = r2_nakagawa(fit)$R2_marginal, r2m_level_only = r2_nakagawa(fit_backbone("bk_laeq_z", outcome, d))$R2_marginal) %>%
  ungroup() %>% select(-fit) %>% mutate(family = fam) %>% bh_family(p_col = "p")
sc1 <- screen(PRIMARY, "primary"); sc2 <- screen(SECONDARY, "secondary")
w(sc1, "single_index_screen_primary.csv"); w(sc2, "single_index_screen_secondary.csv")
cat("\n== primary screen (β per SD, level adjusted; family 48) ==\n")
print(sc1 %>% mutate(across(c(beta, r2m, r2m_level_only), ~round(.x, 3)), across(c(p, q), ~signif(.x, 2))) %>% select(outcome, index, beta, p, q, r2m, r2m_level_only) %>% arrange(outcome, q), n = 48, width = 150)
cat("\n== secondary screen: q < 0.05 ==\n")
print(sc2 %>% filter(q < 0.05) %>% mutate(across(c(beta, r2m), ~round(.x, 3)), across(c(p, q), ~signif(.x, 2))) %>% select(outcome, index, beta, p, q, r2m, family_size) %>% arrange(outcome, q), n = 100, width = 150)
# replication on the full-window indices (same family, region → full)
dF <- load_p35("wav") %>% mutate(env_id = as.character(env_id)) %>% inner_join(ful, by = "env_id") %>% mutate(env_id = factor(env_id))
for (v in PRIMARY) dF[[paste0(v, "_z")]] <- as.numeric(scale(dF[[v]]))
scF <- expand_grid(outcome = OUTS, index = PRIMARY) %>% rowwise() %>% mutate(fit = list(fit_backbone(paste0("bk_laeq_z + ", index, "_z"), outcome, dF)), beta_full = fixef(fit)[3], p_full = summary(fit)$coefficients[3, 5]) %>% ungroup() %>% select(-fit)
w(sc1 %>% left_join(scF, by = c("outcome", "index")), "single_index_screen_primary_with_fullwindow.csv")
cat("\nsign agreement region vs full window (primary, 48):", sum(sign(sc1$beta) == sign(scF$beta_full[match(paste(sc1$outcome, sc1$index), paste(scF$outcome, scF$index))])), "/ 48\n")

# ---- 4. increment beyond the five-feature model (A11): each primary index added; and the whole primary set (PCA scores) added
base <- "bk_laeq_z + tonality_z + src_natural + src_human + src_mechanical"
inc <- expand_grid(outcome = c("richness", "comfort", "rating_appropriateness", "rating_noisiness"), index = PRIMARY) %>% rowwise() %>% mutate(
  m0 = list(fit_backbone(base, outcome, d)), m1 = list(fit_backbone(paste0(base, " + ", index, "_z"), outcome, d)),
  beta = fixef(m1)[[paste0(index, "_z")]], p = summary(m1)$coefficients[paste0(index, "_z"), 5], r2m_base = r2_nakagawa(m0)$R2_marginal, r2m_with = r2_nakagawa(m1)$R2_marginal,
  d_aic = AIC(refitML(m1)) - AIC(refitML(m0))) %>% ungroup() %>% select(-m0, -m1) %>% group_by(outcome) %>% mutate(q = p.adjust(p, "BH"), family_size = n()) %>% ungroup()
w(inc, "increment_beyond_five_feature.csv")
cat("\n== increment beyond LAeq + tonality + 3 groups (per outcome family 12) ==\n")
print(inc %>% mutate(across(c(beta, r2m_base, r2m_with), ~round(.x, 3)), d_aic = round(d_aic, 1), across(c(p, q), ~signif(.x, 2))) %>% select(outcome, index, beta, p, q, r2m_base, r2m_with, d_aic) %>% arrange(outcome, q), n = 48, width = 150)

# ---- 5. cross-validated prediction (A11 machinery): LOSO / LOEO recording-mean R²
cv <- function(y, feats, fold_var, data) { folds <- unique(data[[fold_var]]); pred <- rep(NA_real_, nrow(data))
  for (f in folds) { tr <- data[[fold_var]] != f; m <- lm(reformulate(feats, y), data = data[tr, ]); pred[!tr] <- predict(m, newdata = data[!tr, ]) }
  env <- tibble(env_id = data$env_id, y = data[[y]], p = pred) %>% group_by(env_id) %>% summarise(y = mean(y), p = mean(p), .groups = "drop")
  c(r2_record = 1 - sum((data[[y]] - pred)^2) / sum((data[[y]] - mean(data[[y]]))^2), r2_env_mean = 1 - sum((env$y - env$p)^2) / sum((env$y - mean(env$y))^2)) }
IDX_CORE <- c("ACI", "ADI", "BI", "NDSI_kasten", "Ht", "Hf")  # the six named in the plan's prediction (one NDSI, one H pair)
FAM <- list(level_only = "bk_laeq", level_tonality = c("bk_laeq", "tonality"), indices_core = IDX_CORE, indices_primary = PRIMARY,
            level_indices_core = c("bk_laeq", "tonality", IDX_CORE), level_indices_primary = c("bk_laeq", "tonality", PRIMARY),
            five_feature = c("bk_laeq", "tonality", "src_natural", "src_human", "src_mechanical"),
            five_feature_indices_core = c("bk_laeq", "tonality", "src_natural", "src_human", "src_mechanical", IDX_CORE),
            groups_indices_core = c("src_natural", "src_human", "src_mechanical", IDX_CORE))
cvres <- expand_grid(outcome = OUTS, family = names(FAM), scheme = c("site", "evaluator")) %>% rowwise() %>%
  mutate(k = length(FAM[[family]]), r = list(cv(outcome, FAM[[family]], scheme, d)), r2_record = r[["r2_record"]], r2_env_mean = r[["r2_env_mean"]]) %>% ungroup() %>% select(-r)
w(cvres, "cv_r2_by_family.csv")
cat("\n== cross-validated recording-mean R² (site = new place, evaluator = new person) ==\n")
print(cvres %>% mutate(across(c(r2_record, r2_env_mean), ~round(.x, 3))) %>% select(-r2_record) %>% pivot_wider(names_from = scheme, values_from = r2_env_mean) %>% arrange(outcome, desc(site)), n = 40, width = 150)

# ---- 6. are the indices measuring what people report? recording-level Spearman with perceived group means and total salience (family 12 × 5)
perc <- expand_grid(index = PRIMARY, target = c("src_natural", "src_human", "src_mechanical", "total_sal", "n_sources")) %>% rowwise() %>% mutate(
  rho = cor(envl[[index]], envl[[target]], method = "spearman"), p = cor.test(envl[[index]], envl[[target]], method = "spearman", exact = FALSE)$p.value) %>% ungroup() %>% bh_family(p_col = "p")
w(perc, "indices_vs_perceived_env.csv")
cat("\n== indices vs perceived salience (143 recordings; family 60) ==\n")
print(perc %>% mutate(rho = round(rho, 2), across(c(p, q), ~signif(.x, 2))) %>% select(-family_size) %>% pivot_wider(names_from = target, values_from = c(rho, q)) %>% select(index, starts_with("rho"), starts_with("q_")), n = 20, width = 200)
# the same partialling LAeq (recording level): residual correlations
part <- expand_grid(index = PRIMARY, target = c("src_natural", "src_mechanical", "total_sal", "richness", "comfort")) %>% rowwise() %>% mutate(
  rho_partial_laeq = { rx <- resid(lm(envl[[index]] ~ envl$bk_laeq)); ry <- resid(lm(envl[[target]] ~ envl$bk_laeq)); cor(rx, ry, method = "spearman") }) %>% ungroup()
w(part, "indices_vs_perceived_env_partial_laeq.csv")
cat("\n== partial (given LAeq) Spearman, recording level ==\n"); print(part %>% mutate(rho_partial_laeq = round(rho_partial_laeq, 2)) %>% pivot_wider(names_from = target, values_from = rho_partial_laeq), n = 20, width = 150)

# ---- 7. flagged-row sensitivity for any primary result at q < 0.05 (richness / comfort)
hits <- sc1 %>% filter(q < 0.05, outcome %in% c("richness", "comfort"))
if (nrow(hits)) { sens <- hits %>% rowwise() %>% mutate(beta_noflag = fixef(fit_backbone(paste0("bk_laeq_z + ", index, "_z"), outcome, filter(d, !flagged)))[3]) %>% ungroup() %>% select(outcome, index, beta, beta_noflag)
  w(sens, "sensitivity_no_flagged.csv"); cat("\n== no-flagged-rows sensitivity ==\n"); print(sens %>% mutate(across(c(beta, beta_noflag), ~round(.x, 3)))) }
