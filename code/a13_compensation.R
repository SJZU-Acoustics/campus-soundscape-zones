# A13 — does natural salience offset the mechanical penalty? Additive offset ratio (natural points per mechanical point) with a
# cluster bootstrap over recordings (1,000), the natural × mechanical interaction, the birds × vehicles single-source case,
# a level three-way check, and the 990/XL2 replication. Salience groups on their raw 0–10 scale (group means).
source(file.path("code", "helpers.R"))
suppressPackageStartupMessages({library(performance); library(emmeans)})
OUT <- file.path("output", "analysis_13_compensation"); w <- function(x, f) write_csv(x, file.path(OUT, f))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
cen <- function(x) x %>% mutate(nat_c = src_natural - mean(src_natural), mech_c = src_mechanical - mean(src_mechanical), hum_c = src_human - mean(src_human))
d <- cen(load_p35("wav")); dall <- cen(load_p35("all"))
cat("group means (0-10): natural", round(mean(d$src_natural), 2), "human", round(mean(d$src_human), 2), "mechanical", round(mean(d$src_mechanical), 2), "\n")
OUTS <- c("comfort", "rating_appropriateness", "rating_noisiness", "richness")
# ---- 1. additive model and offset ratio; interaction model
res <- list()
for (y in OUTS) {
  ma <- fit_backbone("bk_laeq_z + nat_c + mech_c + hum_c", y, d); mi <- fit_backbone("bk_laeq_z + nat_c * mech_c + hum_c", y, d)
  ca <- summary(ma)$coefficients; ci <- summary(mi)$coefficients
  res[[y]] <- tibble(outcome = y, b_nat = ca["nat_c", 1], b_mech = ca["mech_c", 1], p_nat = ca["nat_c", 5], p_mech = ca["mech_c", 5],
                     offset_ratio = -ca["mech_c", 1] / ca["nat_c", 1], b_int = ci["nat_c:mech_c", 1], p_int = ci["nat_c:mech_c", 5],
                     b_nat_int = ci["nat_c", 1], b_mech_int = ci["mech_c", 1], r2m_add = r2_nakagawa(ma)$R2_marginal, r2m_int = r2_nakagawa(mi)$R2_marginal)
}
res <- bind_rows(res) %>% mutate(q_int = p.adjust(p_int, "BH"), family_size = n())
w(res, "additive_and_interaction.csv"); print(res %>% mutate(across(where(is.numeric), ~signif(.x, 3))), width = 200)
# ---- 2. cluster bootstrap (recordings) of the offset ratio for comfort and appropriateness
set.seed(13); B <- 1000; envs <- unique(d$env_id)
boot <- replicate(B, {
  s <- sample(envs, replace = TRUE); dd <- resample_stops(d, s)
  sapply(c("comfort", "rating_appropriateness"), function(y) { f <- lme4::fixef(checked_resample_lmer(as.formula(paste0(y, " ~ bk_laeq_z + nat_c + mech_c + hum_c + (1 | evaluator) + (1 | env_id)")), data = dd)); -f["mech_c"] / f["nat_c"] })
})
bt <- tibble(outcome = c("comfort", "rating_appropriateness"), ratio = res$offset_ratio[match(outcome, res$outcome)],
             ci_lo = apply(boot, 1, quantile, 0.025), ci_hi = apply(boot, 1, quantile, 0.975), n_boot = B)
w(bt, "offset_ratio_bootstrap.csv"); print(bt %>% mutate(across(where(is.numeric), ~round(.x, 2))))
# ---- 3. single-source case: birds × vehicles (level adjusted)
for (y in c("comfort", "rating_appropriateness")) {
  ma <- fit_backbone("bk_laeq_z + src_birds + src_vehicles", y, d); mi <- fit_backbone("bk_laeq_z + src_birds * src_vehicles", y, d)
  ca <- summary(ma)$coefficients; ci <- summary(mi)$coefficients
  cat(y, ": birds", round(ca["src_birds", 1], 4), "vehicles", round(ca["src_vehicles", 1], 4), "offset ratio", round(-ca["src_vehicles", 1] / ca["src_birds", 1], 2),
      "| birds:vehicles", round(ci["src_birds:src_vehicles", 1], 4), "p", signif(ci["src_birds:src_vehicles", 5], 2), "\n")
}
# ---- 4. does the compensation depend on level? (one three-way test per outcome)
for (y in c("comfort", "rating_appropriateness")) {
  m3 <- fit_backbone("bk_laeq_z * nat_c * mech_c + hum_c", y, d); c3 <- summary(m3)$coefficients
  cat(y, ": nat×mech×level", round(c3["bk_laeq_z:nat_c:mech_c", 1], 4), "p", signif(c3["bk_laeq_z:nat_c:mech_c", 5], 2), "| nat×level", round(c3["bk_laeq_z:nat_c", 1], 4), "p", signif(c3["bk_laeq_z:nat_c", 5], 2), "| mech×level", round(c3["bk_laeq_z:mech_c", 1], 4), "p", signif(c3["bk_laeq_z:mech_c", 5], 2), "\n")
}
# ---- 5. prediction grid for display (comfort, appropriateness at natural 1/3/5 × mechanical 0/2/4, mean level)
mi <- fit_backbone("bk_laeq_z + nat_c * mech_c + hum_c", "comfort", d); ma2 <- fit_backbone("bk_laeq_z + nat_c * mech_c + hum_c", "rating_appropriateness", d)
grid <- expand_grid(src_natural = c(1, 3, 5), src_mechanical = c(0, 2, 4)) %>% mutate(nat_c = src_natural - mean(d$src_natural), mech_c = src_mechanical - mean(d$src_mechanical), hum_c = 0, bk_laeq_z = 0)
grid$comfort <- predict(mi, newdata = grid, re.form = NA); grid$appropriateness <- predict(ma2, newdata = grid, re.form = NA)
w(grid, "prediction_grid.csv"); print(grid %>% select(src_natural, src_mechanical, comfort, appropriateness) %>% mutate(across(c(comfort, appropriateness), ~round(.x, 2))) %>% pivot_wider(names_from = src_mechanical, values_from = c(comfort, appropriateness)))
# ---- 6. replication on 990 with XL2
for (y in c("comfort", "rating_appropriateness")) {
  ma <- fit_backbone("xl2_laeq_z + nat_c + mech_c + hum_c", y, dall); mi <- fit_backbone("xl2_laeq_z + nat_c * mech_c + hum_c", y, dall); ca <- summary(ma)$coefficients; ci <- summary(mi)$coefficients
  cat("990/XL2", y, ": ratio", round(-ca["mech_c", 1] / ca["nat_c", 1], 2), "| interaction", round(ci["nat_c:mech_c", 1], 4), "p", signif(ci["nat_c:mech_c", 5], 2), "\n")
}
