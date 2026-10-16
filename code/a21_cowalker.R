# A21 — co-walkers' source salience. Each rating's source predictors are replaced by the mean salience reported by the other
# three or four evaluators at the same stop (leave-self-out), so the predictor no longer comes from the rater's own questionnaire.
# Declared before execution and bounded to the established headline results:
#   1. five-feature models of the four outcomes with co-walkers' group salience (666 ratings; 990/XL2 companion);
#   2. the natural–mechanical exchange ratio for Comfort and appropriateness (A13 method; 1,000 cluster bootstraps over stops);
#   3. the natural × mechanical Richness interaction (A20 primary specification; LRT on ML refits);
#   4. the zone-dependent Comfort and Richness slopes on appropriateness (A10 specification) with co-walkers' mean Comfort and Richness;
#   5. the site-held-out hierarchy of stop-mean R²: level → level + tonality → + detected groups → + co-walkers' groups → + own groups,
#      for the four outcomes, with the same folds and fixed-effect predictor as A11/A18; Richness also with the interaction.
# Single pre-specified tests report their own p; the twelve co-walker group coefficients of part 1 form one declared BH family.
source(file.path("code", "helpers.R"))
suppressPackageStartupMessages(library(performance))
OUT <- file.path("output", "analysis_21_cowalker_salience")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
w <- function(x, f) write_csv(x, file.path(OUT, f))
GROUPS <- c("natural", "human", "mechanical")
OUTS <- c("comfort", "richness", "rating_appropriateness", "rating_noisiness")

# leave-self-out mean of a column within each stop
loo_mean <- function(x, g) { s <- ave(x, g, FUN = sum); n <- ave(x, g, FUN = length); (s - x) / (n - 1) }
add_cowalker <- function(d) {
  stopifnot(all(table(droplevels(d$env_id)) >= 4))
  for (g in GROUPS) d[[paste0("cw_", g)]] <- loo_mean(d[[paste0("src_", g)]], d$env_id)
  d$cw_comfort <- loo_mean(d$comfort, d$env_id); d$cw_richness <- loo_mean(d$richness, d$env_id)
  d
}
det <- read_frozen_csv("detections_by_env.csv") %>%
  filter(window == "region") %>% select(env_id, det_natural, det_human, det_mechanical)
d <- load_p35("wav") %>% add_cowalker() %>% mutate(env_id = as.character(env_id)) %>% left_join(det, by = "env_id") %>% mutate(env_id = factor(env_id))
da <- load_p35("all") %>% add_cowalker()
stopifnot(nrow(d) == 666, nrow(da) == 990, all(complete.cases(d[, c("det_natural", "det_human", "det_mechanical")])))
z <- function(x) as.numeric(scale(x))
d <- d %>% mutate(across(c(paste0("cw_", GROUPS), paste0("src_", GROUPS)), z, .names = "{.col}_z"))
da <- da %>% mutate(across(c(paste0("cw_", GROUPS), paste0("src_", GROUPS)), z, .names = "{.col}_z"))

# agreement between own and co-walkers' salience (record level)
agree <- bind_rows(lapply(GROUPS, function(g) tibble(group = g, r_own_cowalker = cor(d[[paste0("src_", g)]], d[[paste0("cw_", g)]]),
  r_own_cowalker_990 = cor(da[[paste0("src_", g)]], da[[paste0("cw_", g)]]))))
w(agree, "own_vs_cowalker_agreement.csv"); print(agree)

# ---- 1. five-feature models: own vs co-walkers' group salience
coefs <- list(); r2 <- list()
for (y in OUTS) for (src in c("own", "cowalker")) for (layer in c("wav", "xl2_990")) {
  pre <- if (src == "own") "src_" else "cw_"
  dat <- if (layer == "wav") d else da
  phys <- if (layer == "wav") "bk_laeq_z + tonality_z" else "xl2_laeq_z"
  m <- fit_backbone(paste(phys, "+", paste0(pre, GROUPS, "_z", collapse = " + ")), y, dat)
  coefs[[paste(y, src, layer)]] <- as_tibble(summary(m)$coefficients, rownames = "term") %>% filter(term != "(Intercept)") %>%
    mutate(outcome = y, salience = src, layer = layer, n = nobs(m),
           lo = Estimate - qt(0.975, df) * `Std. Error`, hi = Estimate + qt(0.975, df) * `Std. Error`)
  r2[[paste(y, src, layer)]] <- tibble(outcome = y, salience = src, layer = layer, r2m = r2_nakagawa(m)$R2_marginal)
}
coefs <- bind_rows(coefs); r2 <- bind_rows(r2)
fam <- coefs %>% filter(salience == "cowalker", layer == "wav", grepl("^cw_", term)) %>% bh_family("Pr(>|t|)") %>% select(outcome, term, salience, layer, q, family_size)
coefs <- coefs %>% left_join(fam, by = c("outcome", "term", "salience", "layer"))
w(coefs, "five_feature_own_vs_cowalker.csv"); w(r2, "five_feature_r2.csv")
cat("\n== five-feature group coefficients (per SD), own vs co-walkers ==\n")
print(coefs %>% filter(grepl("natural|human|mechanical", term)) %>% mutate(across(c(Estimate, lo, hi), ~round(.x, 3)), p = signif(`Pr(>|t|)`, 2), q = signif(q, 2)) %>%
  select(outcome, layer, salience, term, Estimate, lo, hi, p, q), n = 60, width = 200)

# ---- 2. exchange ratio with co-walkers' salience (A13 method)
cen <- function(x) x %>% mutate(nat_c = cw_natural - mean(cw_natural), mech_c = cw_mechanical - mean(cw_mechanical), hum_c = cw_human - mean(cw_human))
dc <- cen(d); dac <- cen(da)
ex <- bind_rows(lapply(c("comfort", "rating_appropriateness"), function(y) {
  ca <- summary(fit_backbone("bk_laeq_z + nat_c + mech_c + hum_c", y, dc))$coefficients
  ci <- summary(fit_backbone("bk_laeq_z + nat_c * mech_c + hum_c", y, dc))$coefficients
  cx <- summary(fit_backbone("xl2_laeq_z + nat_c + mech_c + hum_c", y, dac))$coefficients
  tibble(outcome = y, b_nat = ca["nat_c", 1], p_nat = ca["nat_c", 5], b_mech = ca["mech_c", 1], p_mech = ca["mech_c", 5],
         ratio = -ca["mech_c", 1] / ca["nat_c", 1], b_int = ci["nat_c:mech_c", 1], p_int = ci["nat_c:mech_c", 5],
         ratio_990_xl2 = -cx["mech_c", 1] / cx["nat_c", 1])
}))
set.seed(21); B <- 1000; envs <- unique(as.character(dc$env_id))
boot <- replicate(B, { s <- sample(envs, replace = TRUE); dd <- resample_stops(dc, s)
  sapply(c("comfort", "rating_appropriateness"), function(y) {
    f <- lme4::fixef(checked_resample_lmer(as.formula(paste0(y, " ~ bk_laeq_z + nat_c + mech_c + hum_c + (1 | evaluator) + (1 | env_id)")), data = dd))
    -f["mech_c"] / f["nat_c"] }) })
ex <- ex %>% mutate(ci_lo = apply(boot, 1, quantile, 0.025), ci_hi = apply(boot, 1, quantile, 0.975), n_boot = B)
w(ex, "exchange_ratio_cowalker.csv")
cat("\n== exchange ratio, co-walkers' salience ==\n"); print(ex %>% mutate(across(where(is.numeric), ~signif(.x, 3))), width = 200)

# ---- 3. Richness interaction (A20 primary specification) with co-walkers' salience
ri <- list()
for (src in c("own", "cowalker")) {
  pre <- if (src == "own") "src_" else "cw_"
  base <- paste0("bk_laeq + tonality + ", paste0(pre, GROUPS, collapse = " + "))
  m0 <- fit_backbone(base, "richness", d); m1 <- fit_backbone(paste0(base, " + ", pre, "natural:", pre, "mechanical"), "richness", d)
  lrt <- anova(m0, m1, refit = TRUE); cf <- summary(m1)$coefficients; term <- paste0(pre, "natural:", pre, "mechanical")
  lim <- quantile(d[[paste0(pre, "mechanical")]], c(0.05, 0.95))
  sl <- as_tibble(emmeans::emtrends(m1, as.formula(paste0("~", pre, "mechanical")), var = paste0(pre, "natural"),
                  at = setNames(list(unname(lim)), paste0(pre, "mechanical"))))
  ri[[src]] <- tibble(salience = src, b_int = cf[term, 1], lo = cf[term, 1] - qt(0.975, cf[term, 3]) * cf[term, 2],
                      hi = cf[term, 1] + qt(0.975, cf[term, 3]) * cf[term, 2], chisq = lrt$Chisq[2], p = lrt$`Pr(>Chisq)`[2],
                      slope_p05 = sl[[2]][1], slope_p95 = sl[[2]][2])
}
ri <- bind_rows(ri); w(ri, "richness_interaction_own_vs_cowalker.csv")
cat("\n== Richness natural × mechanical interaction ==\n"); print(ri %>% mutate(across(where(is.numeric), ~signif(.x, 3))), width = 200)

# ---- 4. zone-dependent appropriateness slopes with co-walkers' mean Comfort and Richness (A10 specification, 990/XL2)
za <- list()
for (src in c("own", "cowalker")) {
  cvar <- if (src == "own") "comfort" else "cw_comfort"; rvar <- if (src == "own") "richness" else "cw_richness"
  common <- paste0(cvar, " + ", rvar, " + xl2_laeq_z + zone")
  m0 <- fit_backbone(common, "rating_appropriateness", da)
  mr <- fit_backbone(paste0(common, " + zone:", rvar), "rating_appropriateness", da)
  mc <- fit_backbone(paste0(common, " + zone:", cvar), "rating_appropriateness", da)
  lr_r <- anova(m0, mr, refit = TRUE); lr_c <- anova(m0, mc, refit = TRUE); c0 <- summary(m0)$coefficients
  sr <- as_tibble(emmeans::emtrends(mr, ~zone, var = rvar)); sc <- as_tibble(emmeans::emtrends(mc, ~zone, var = cvar))
  za[[src]] <- tibble(salience = src, b_comfort = c0[cvar, 1], se_comfort = c0[cvar, 2], b_richness = c0[rvar, 1], se_richness = c0[rvar, 2],
                      chisq_zone_richness = lr_r$Chisq[2], p_zone_richness = lr_r$`Pr(>Chisq)`[2],
                      chisq_zone_comfort = lr_c$Chisq[2], p_zone_comfort = lr_c$`Pr(>Chisq)`[2],
                      min_comfort_slope = min(sc[[2]]), n_comfort_slopes_positive_ci = sum(sc$lower.CL > 0),
                      richness_slope_sport = sr[[2]][sr$zone == "sport"], richness_slope_residential = sr[[2]][sr$zone == "residential"])
  names(sr)[2] <- "trend"; names(sc)[2] <- "trend"
  w(bind_rows(sr %>% mutate(slope = "richness"), sc %>% mutate(slope = "comfort")) %>% mutate(salience = src), paste0("zone_slopes_", src, ".csv"))
}
za <- bind_rows(za); w(za, "appropriateness_zone_own_vs_cowalker.csv")
cat("\n== appropriateness: axes and zone interactions ==\n"); print(za %>% mutate(across(where(is.numeric), ~signif(.x, 3))), width = 250)

# ---- 5. site-held-out hierarchy (A11/A18 procedure: lm on seven sites, predict the eighth, stop-mean R² pooled)
cv <- function(y, feats) {
  pred <- rep(NA_real_, nrow(d))
  for (s in unique(d$site)) { tr <- d$site != s; m <- lm(reformulate(feats, y), data = d[tr, ]); pred[!tr] <- predict(m, newdata = d[!tr, ]) }
  e <- tibble(env_id = d$env_id, site = d$site, y = d[[y]], p = pred) %>% group_by(env_id, site) %>% summarise(y = mean(y), p = mean(p), .groups = "drop")
  list(r2 = 1 - sum((e$y - e$p)^2) / sum((e$y - mean(e$y))^2),
       folds = e %>% group_by(site) %>% summarise(n_stops = n(), mse = mean((y - p)^2), .groups = "drop"))
}
LT <- c("bk_laeq", "tonality")
H <- list(level = "bk_laeq", level_tonality = LT,
          detected = c(LT, paste0("det_", GROUPS)),
          cowalker = c(LT, paste0("cw_", GROUPS)),
          own = c(LT, paste0("src_", GROUPS)),
          cowalker_interaction = c(LT, paste0("cw_", GROUPS), "cw_natural:cw_mechanical"),
          own_interaction = c(LT, paste0("src_", GROUPS), "src_natural:src_mechanical"))
hier <- list(); folds <- list()
for (y in OUTS) for (h in names(H)) {
  if (grepl("interaction", h) && y != "richness") next
  r <- cv(y, H[[h]]); hier[[paste(y, h)]] <- tibble(outcome = y, features = h, k = length(H[[h]]), r2_stop_mean = r$r2)
  folds[[paste(y, h)]] <- r$folds %>% mutate(outcome = y, features = h)
}
hier <- bind_rows(hier); w(hier, "heldout_hierarchy.csv"); w(bind_rows(folds), "heldout_hierarchy_site_errors.csv")
cat("\n== site-held-out stop-mean R² hierarchy ==\n")
print(hier %>% mutate(r2_stop_mean = round(r2_stop_mean, 3)) %>% select(-k) %>% pivot_wider(names_from = outcome, values_from = r2_stop_mean), width = 200)
writeLines(c(capture.output(sessionInfo())), file.path(OUT, "session_info.txt"))

# ---- 6. ADDED AFTER parts 1–5 had run (16 October 2026), not in the declared set: the A16 quadrant contrast with quadrants defined
#         from the co-walkers' mean Comfort and Richness instead of the rater's own (990/XL2; level and zone adjusted; Holm over six pairs).
QL <- c("calm", "vibrant", "monotonous", "chaotic")
# A leave-self-out mean that is exactly zero can come out as +/-1e-17 in floating point, and its sign would then decide the
# quadrant. Values within 1e-9 of zero are set to zero, so neutral means are excluded as exact zeros are for the rater's own
# ratings (the smallest non-zero |mean| is about 0.003).
quad <- function(c, r) {
  c[abs(c) < 1e-9] <- 0; r[abs(r) < 1e-9] <- 0
  factor(case_when(c == 0 | r == 0 ~ NA_character_, c > 0 & r <= 0 ~ "calm", c > 0 & r > 0 ~ "vibrant",
                   c <= 0 & r <= 0 ~ "monotonous", TRUE ~ "chaotic"), levels = QL)
}
qd <- list()
for (src in c("own", "cowalker")) {
  dq <- da %>% mutate(quadrant = if (src == "own") quad(comfort, richness) else quad(cw_comfort, cw_richness)) %>% filter(!is.na(quadrant))
  m0 <- fit_backbone("xl2_laeq_z + zone + quadrant", "rating_appropriateness", dq)
  pw <- as_tibble(pairs(emmeans::emmeans(m0, "quadrant"), adjust = "holm")) %>% mutate(salience = src, n = nrow(dq))
  qd[[src]] <- pw
}
qd <- bind_rows(qd); w(qd, "quadrant_pairwise_own_vs_cowalker.csv")
cat("\n== appropriateness by quadrant (Holm, six pairs) ==\n"); print(qd %>% mutate(across(c(estimate, SE), ~round(.x, 2)), p.value = signif(p.value, 2)) %>% select(salience, n, contrast, estimate, SE, p.value), n = 20)
