# A18 — detected sources (PANNs CNN14; the deposited detected_sources table) vs perceived sources, and detected sources →
# the four outcomes, LOSO/LOEO, and the detected exchange rate. The AudioSet-to-source mapping (deposited audioset_source_map)
# and the score definitions were fixed before any detection was correlated with a rating. Predictions: detected bird/insect/water
# presence vs perceived natural r > 0.5 at recording level; vehicles/engines vs perceived mechanical r > 0.5; speech/footsteps vs human
# 0.3–0.5; detected groups predict comfort/appropriateness with A05's signs and comfort LOSO recording-mean R² ≥ 0.35; the 2.5 : 1
# exchange rate reproduced in sign.
source(file.path("code", "helpers.R"))
suppressPackageStartupMessages({library(performance)})
OUT <- file.path("output", "analysis_18_detected_sources"); w <- function(x, f) write_csv(x, file.path(OUT, f))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
det <- read_frozen_csv("detections_by_env.csv")
reg <- det %>% filter(window == "region") %>% select(-window, -duration_s, -n_seconds, -top5)
ful <- det %>% filter(window == "full") %>% select(-window, -duration_s, -n_seconds, -top5)
GROUPS <- c("natural", "human", "mechanical")
d <- load_p35("wav") %>% mutate(env_id = as.character(env_id)) %>% inner_join(reg, by = "env_id") %>% mutate(env_id = factor(env_id))
stopifnot(nrow(d) == 666)
src_cols <- paste0("src_", SOURCES)
d <- d %>% mutate(total_sal = rowSums(across(all_of(src_cols))))
OUTS <- c("comfort", "rating_appropriateness", "rating_noisiness", "richness")

# ---- 1. recording-level agreement: detected vs perceived, per source (24) and per group (3); three score types
envl <- d %>% group_by(env_id, site, zone) %>% summarise(across(c(all_of(src_cols), src_natural, src_human, src_mechanical, total_sal, comfort, richness, rating_appropriateness, rating_noisiness, bk_laeq, tonality, naturalness,
  starts_with("det_"), starts_with("detmax_"), starts_with("sedmean_"), starts_with("sedfrac_"), starts_with("ctx_")), ~mean(.x, na.rm = TRUE)), .groups = "drop")
cat("recordings:", nrow(envl), "\n")
agree <- expand_grid(source = c(SOURCES, GROUPS), score = c("det", "sedmean", "sedfrac")) %>% rowwise() %>% mutate(
  perceived = if (source %in% GROUPS) paste0("src_", source) else paste0("src_", source), detected = paste0(score, "_", source),
  perceived_mean = mean(envl[[perceived]]), perceived_prevalence = mean(envl[[perceived]] > 0), detected_mean = mean(envl[[detected]]),
  rho = suppressWarnings(cor(envl[[detected]], envl[[perceived]], method = "spearman")),
  p = suppressWarnings(cor.test(envl[[detected]], envl[[perceived]], method = "spearman", exact = FALSE)$p.value), n = nrow(envl)) %>% ungroup() %>%
  mutate(level = ifelse(source %in% GROUPS, "group", "source")) %>% group_by(score) %>% mutate(q = p.adjust(p, "BH"), family_size = n()) %>% ungroup()
w(agree, "agreement_detected_vs_perceived.csv")
cat("\n== agreement (Spearman, 143 recordings; BH within score type, family 27) — clip-probability score ==\n")
print(agree %>% filter(score == "det") %>% mutate(across(c(perceived_mean, detected_mean, rho), ~round(.x, 2)), perceived_prevalence = round(perceived_prevalence, 2), across(c(p, q), ~signif(.x, 2))) %>%
  select(source, level, perceived_mean, perceived_prevalence, detected_mean, rho, p, q) %>% arrange(level, desc(rho)), n = 30, width = 150)
cat("\n== the three score types, groups only ==\n")
print(agree %>% filter(level == "group") %>% mutate(rho = round(rho, 2), q = signif(q, 2)) %>% select(source, score, rho, q) %>% pivot_wider(names_from = score, values_from = c(rho, q)), width = 150)
# the group score built as max over sources (secondary)
gm <- tibble(group = GROUPS) %>% rowwise() %>% mutate(rho_detmax = cor(envl[[paste0("detmax_", group)]], envl[[paste0("src_", group)]], method = "spearman"),
  rho_det_mean = cor(envl[[paste0("det_", group)]], envl[[paste0("src_", group)]], method = "spearman")) %>% ungroup(); print(gm %>% mutate(across(where(is.numeric), ~round(.x, 2))))
# balance and context
bal <- tibble(pair = c("det_balance vs perceived natural−mechanical", "ctx_outside_natural vs perceived natural", "ctx_outside_natural vs visual naturalness", "ctx_outside_urban vs perceived mechanical", "ctx_outside_urban − ctx_outside_natural vs perceived mech−nat"),
  rho = c(cor(envl$det_balance, envl$src_natural - envl$src_mechanical, method = "spearman"), cor(envl$ctx_outside_natural, envl$src_natural, method = "spearman"),
          cor(envl$ctx_outside_natural, envl$naturalness, method = "spearman"), cor(envl$ctx_outside_urban, envl$src_mechanical, method = "spearman"),
          cor(envl$ctx_outside_urban - envl$ctx_outside_natural, envl$src_mechanical - envl$src_natural, method = "spearman")))
w(bal, "balance_and_context.csv"); print(bal %>% mutate(rho = round(rho, 2)))
# top disagreements for the three headline sources (perceived rank vs detected rank)
dis <- list()
for (s in c("birds", "vehicles", "conversation", "water", "ac_units")) {
  e <- envl %>% transmute(env_id, site, zone, perceived = .data[[paste0("src_", s)]], detected = .data[[paste0("det_", s)]]) %>%
    mutate(rank_gap = rank(perceived) - rank(detected)) %>% arrange(desc(abs(rank_gap))) %>% head(4) %>% mutate(source = s)
  dis[[s]] <- e }
dis <- bind_rows(dis) %>% left_join(det %>% filter(window == "region") %>% select(env_id, top5) %>% mutate(env_id = factor(env_id)), by = "env_id")
w(dis, "top_disagreements.csv"); cat("\n== largest rank disagreements (perceived vs detected) ==\n"); print(dis %>% mutate(across(c(perceived, detected), ~round(.x, 2))) %>% select(source, env_id, zone, perceived, detected, rank_gap, top5), n = 25, width = 200)
# by-zone agreement for the groups (does the tagger work everywhere?)
bz <- bind_rows(lapply(split(envl, envl$zone), function(e) tibble(zone = e$zone[1], n = nrow(e), rho_natural = suppressWarnings(cor(e$det_natural, e$src_natural, method = "spearman")),
  rho_human = suppressWarnings(cor(e$det_human, e$src_human, method = "spearman")), rho_mechanical = suppressWarnings(cor(e$det_mechanical, e$src_mechanical, method = "spearman")))))
w(bz, "agreement_by_zone.csv"); print(bz %>% mutate(across(starts_with("rho"), ~round(.x, 2))))
# region vs full window agreement of the detected scores
wa <- tibble(score = c(paste0("det_", GROUPS), paste0("det_", c("birds", "vehicles", "conversation", "water", "insects", "footsteps")))) %>% rowwise() %>%
  mutate(rho_region_full = cor(reg[[score]], ful[[score]], method = "spearman")) %>% ungroup(); w(wa, "window_agreement.csv"); print(wa %>% mutate(rho_region_full = round(rho_region_full, 2)))

# ---- 2. models: detected groups → four outcomes (LAeq adjusted, crossed RE); family 12; then head-to-head with perceived groups
dz <- d %>% mutate(across(c(det_natural, det_human, det_mechanical, src_natural, src_human, src_mechanical, sedfrac_natural, sedfrac_human, sedfrac_mechanical), ~as.numeric(scale(.x)), .names = "{.col}_z"))
grp <- list(); hh <- list()
for (y in OUTS) {
  m0 <- fit_backbone("bk_laeq_z", y, dz); m1 <- fit_backbone("bk_laeq_z + det_natural_z + det_human_z + det_mechanical_z", y, dz)
  m2 <- fit_backbone("bk_laeq_z + tonality_z + det_natural_z + det_human_z + det_mechanical_z", y, dz)
  mp <- fit_backbone("bk_laeq_z + tonality_z + src_natural_z + src_human_z + src_mechanical_z", y, dz)
  m3 <- fit_backbone("bk_laeq_z + tonality_z + det_natural_z + det_human_z + det_mechanical_z + src_natural_z + src_human_z + src_mechanical_z", y, dz)
  ms <- fit_backbone("bk_laeq_z + sedfrac_natural_z + sedfrac_human_z + sedfrac_mechanical_z", y, dz)
  grp[[y]] <- as_tibble(summary(m1)$coefficients, rownames = "term") %>% filter(grepl("^det_", term)) %>% mutate(outcome = y, r2m_level_only = r2_nakagawa(m0)$R2_marginal, r2m_detected = r2_nakagawa(m1)$R2_marginal,
    r2m_level_tonality_detected = r2_nakagawa(m2)$R2_marginal, r2m_level_tonality_perceived = r2_nakagawa(mp)$R2_marginal, r2m_both = r2_nakagawa(m3)$R2_marginal, r2m_sedfrac = r2_nakagawa(ms)$R2_marginal)
  hh[[y]] <- as_tibble(summary(m3)$coefficients, rownames = "term") %>% filter(grepl("_z$", term), term != "bk_laeq_z", term != "tonality_z") %>% mutate(outcome = y)
}
grp <- bind_rows(grp) %>% bh_family(p_col = "Pr(>|t|)"); hh <- bind_rows(hh)
w(grp, "detected_groups_models.csv"); w(hh, "head_to_head_detected_vs_perceived.csv")
cat("\n== detected groups → outcomes (β per SD, LAeq adjusted; family 12) ==\n")
print(grp %>% mutate(across(c(Estimate, starts_with("r2m")), ~round(.x, 3)), p = signif(`Pr(>|t|)`, 2), q = signif(q, 2)) %>% select(outcome, term, Estimate, p, q, r2m_level_only, r2m_detected, r2m_level_tonality_detected, r2m_level_tonality_perceived, r2m_both, r2m_sedfrac), n = 12, width = 200)
cat("\n== head-to-head: detected and perceived groups in one model ==\n")
print(hh %>% mutate(Estimate = round(Estimate, 3), p = signif(`Pr(>|t|)`, 2)) %>% select(outcome, term, Estimate, p) %>% pivot_wider(names_from = outcome, values_from = c(Estimate, p)), n = 12, width = 200)

# ---- 3. cross-validated prediction (A11 machinery): LOSO / LOEO recording-mean R²
cv <- function(y, feats, fold_var, data) { folds <- unique(data[[fold_var]]); pred <- rep(NA_real_, nrow(data))
  for (f in folds) { tr <- data[[fold_var]] != f; m <- lm(reformulate(feats, y), data = data[tr, ]); pred[!tr] <- predict(m, newdata = data[!tr, ]) }
  env <- tibble(env_id = data$env_id, y = data[[y]], p = pred) %>% group_by(env_id) %>% summarise(y = mean(y), p = mean(p), .groups = "drop")
  c(r2_record = 1 - sum((data[[y]] - pred)^2) / sum((data[[y]] - mean(data[[y]]))^2), r2_env_mean = 1 - sum((env$y - env$p)^2) / sum((env$y - mean(env$y))^2)) }
DET24 <- paste0("det_", SOURCES)
FAM <- list(level_only = "bk_laeq", level_tonality = c("bk_laeq", "tonality"), detected_groups = c("det_natural", "det_human", "det_mechanical"),
            level_tonality_detected = c("bk_laeq", "tonality", "det_natural", "det_human", "det_mechanical"),
            level_tonality_detected24 = c("bk_laeq", "tonality", DET24),
            level_tonality_sedfrac = c("bk_laeq", "tonality", "sedfrac_natural", "sedfrac_human", "sedfrac_mechanical"),
            five_feature_perceived = c("bk_laeq", "tonality", "src_natural", "src_human", "src_mechanical"),
            perceived_plus_detected = c("bk_laeq", "tonality", "src_natural", "src_human", "src_mechanical", "det_natural", "det_human", "det_mechanical"))
cvres <- expand_grid(outcome = OUTS, family = names(FAM), scheme = c("site", "evaluator")) %>% rowwise() %>%
  mutate(k = length(FAM[[family]]), r = list(cv(outcome, FAM[[family]], scheme, d)), r2_record = r[["r2_record"]], r2_env_mean = r[["r2_env_mean"]]) %>% ungroup() %>% select(-r)
w(cvres, "cv_r2_by_family.csv")
cat("\n== cross-validated recording-mean R² ==\n")
print(cvres %>% mutate(across(c(r2_record, r2_env_mean), ~round(.x, 3))) %>% select(-r2_record) %>% pivot_wider(names_from = scheme, values_from = r2_env_mean) %>% arrange(outcome, desc(site)), n = 40, width = 150)

# ---- 4. exchange rate from detected scores (A13 method): −b_mech / b_nat on the raw 0–1 scale and per SD; cluster bootstrap over recordings (1,000)
cen <- function(x) x %>% mutate(dnat_c = det_natural - mean(det_natural), dmech_c = det_mechanical - mean(det_mechanical), dhum_c = det_human - mean(det_human))
dc <- cen(d)
ex <- list()
for (y in c("comfort", "rating_appropriateness")) {
  ma <- fit_backbone("bk_laeq_z + dnat_c + dmech_c + dhum_c", y, dc); mi <- fit_backbone("bk_laeq_z + dnat_c * dmech_c + dhum_c", y, dc); ca <- summary(ma)$coefficients; ci <- summary(mi)$coefficients
  ex[[y]] <- tibble(outcome = y, b_nat = ca["dnat_c", 1], p_nat = ca["dnat_c", 5], b_mech = ca["dmech_c", 1], p_mech = ca["dmech_c", 5], ratio_raw = -ca["dmech_c", 1] / ca["dnat_c", 1],
                    ratio_per_sd = (-ca["dmech_c", 1] * sd(d$det_mechanical)) / (ca["dnat_c", 1] * sd(d$det_natural)), b_int = ci["dnat_c:dmech_c", 1], p_int = ci["dnat_c:dmech_c", 5]) }
ex <- bind_rows(ex)
set.seed(18); B <- 1000; envs <- unique(dc$env_id)
boot <- replicate(B, { s <- sample(envs, replace = TRUE); dd <- resample_stops(dc, s)
  sapply(c("comfort", "rating_appropriateness"), function(y) { f <- lme4::fixef(checked_resample_lmer(as.formula(paste0(y, " ~ bk_laeq_z + dnat_c + dmech_c + dhum_c + (1 | evaluator) + (1 | env_id)")), data = dd)); (-f["dmech_c"] * sd(d$det_mechanical)) / (f["dnat_c"] * sd(d$det_natural)) }) })
ex <- ex %>% mutate(ratio_per_sd_ci_lo = apply(boot, 1, quantile, 0.025), ratio_per_sd_ci_hi = apply(boot, 1, quantile, 0.975), n_boot = B,
                    boot_frac_positive = rowMeans(boot > 0))
w(ex, "detected_exchange_rate.csv"); cat("\n== detected exchange rate (per SD of the detected scores; compare A13 perceived-group ratios in their own raw units) ==\n"); print(ex %>% mutate(across(where(is.numeric), ~signif(.x, 3))), width = 200)
# the perceived ratio per SD for the same comparison
ps <- list(); for (y in c("comfort", "rating_appropriateness")) { f <- fixef(fit_backbone("bk_laeq_z + src_natural + src_mechanical + src_human", y, d)); ps[[y]] <- tibble(outcome = y, perceived_ratio_per_sd = (-f["src_mechanical"] * sd(d$src_mechanical)) / (f["src_natural"] * sd(d$src_natural))) }
print(bind_rows(ps) %>% mutate(perceived_ratio_per_sd = round(perceived_ratio_per_sd, 2)))

# ---- 5. replication on the full-window detections and without flagged rows (group model, comfort/appropriateness)
dF <- load_p35("wav") %>% mutate(env_id = as.character(env_id)) %>% inner_join(ful, by = "env_id") %>% mutate(env_id = factor(env_id)) %>%
  mutate(across(c(det_natural, det_human, det_mechanical), ~as.numeric(scale(.x)), .names = "{.col}_z"))
rep <- list()
for (y in OUTS) { mF <- fit_backbone("bk_laeq_z + det_natural_z + det_human_z + det_mechanical_z", y, dF); mN <- fit_backbone("bk_laeq_z + det_natural_z + det_human_z + det_mechanical_z", y, filter(dz, !flagged))
  rep[[y]] <- tibble(outcome = y, term = c("det_natural_z", "det_human_z", "det_mechanical_z"), beta_full_window = fixef(mF)[3:5], p_full_window = summary(mF)$coefficients[3:5, 5], beta_no_flagged = fixef(mN)[3:5]) }
rep <- bind_rows(rep); w(rep, "replication_fullwindow_noflagged.csv"); cat("\n== replication: full window / no flagged rows ==\n"); print(rep %>% mutate(across(where(is.numeric), ~signif(.x, 3))), n = 12)
