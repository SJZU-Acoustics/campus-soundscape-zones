# A03 — how the eight functional zones differ in perception, level and source make-up.
source(file.path("code", "helpers.R"))
suppressPackageStartupMessages({ library(emmeans); library(performance) })
OUT <- file.path("output", "analysis_03_zone_profiles"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
w <- function(x, f) write_csv(x, file.path(OUT, f))
d <- load_p35("all"); e <- load_p35_env() %>% mutate(zone = factor(ZONE_LEVELS[site], levels = ZONE_LEVELS))
emm_options(lmer.df = "satterthwaite")

# ---- 1. outcomes by zone (all 990; evaluator + env random intercepts)
em_all <- list(); contrasts <- list(); omni <- list()
for (y in c("comfort", "richness", "rating_appropriateness", "rating_noisiness")) {
  m <- fit_backbone("zone", y, d); m0 <- fit_backbone("1", y, d)
  # Different fixed effects require ML likelihoods; retain REML for the EMMs.
  lrt <- anova(m0, m, refit = TRUE)
  omni[[y]] <- tibble(outcome = y, chisq = lrt$Chisq[2], df = lrt$Df[2], p_lrt = lrt$`Pr(>Chisq)`[2], r2_marginal = r2_nakagawa(m)$R2_marginal)
  em <- emmeans(m, ~ zone)
  em_all[[y]] <- as_tibble(em) %>% mutate(outcome = y)
  contrasts[[y]] <- as_tibble(pairs(em, adjust = "holm")) %>% mutate(outcome = y)
}
omni <- bind_rows(omni); w(omni, "zone_omnibus.csv"); print(omni %>% mutate(across(where(is.numeric), ~signif(.x, 3))))
em_all <- bind_rows(em_all); w(em_all, "zone_emmeans.csv")
print(em_all %>% select(outcome, zone, emmean, SE) %>% mutate(across(c(emmean, SE), ~round(.x, 3))) %>% pivot_wider(names_from = outcome, values_from = c(emmean, SE)), width = 200)
contrasts <- bind_rows(contrasts); w(contrasts, "zone_pairwise_holm.csv")
sig <- contrasts %>% filter(p.value < 0.05) %>% group_by(outcome) %>% summarise(n_sig = n(), pairs = paste(contrast, collapse = "; "))
cat("\nHolm-significant pairwise zone contrasts (28 per outcome):\n"); print(sig, width = 300)

# ---- 2. level and visual by zone (env level)
lev <- list()
for (v in c("xl2_laeq", "bk_laeq", "n5", "s", "naturalness", "artificiality", "vdi")) {
  dd <- e %>% filter(!is.na(.data[[v]]))
  a <- anova(lm(reformulate("zone", v), data = dd))
  lev[[v]] <- tibble(variable = v, n_env = nrow(dd), F = a$`F value`[1], p = a$`Pr(>F)`[1], eta2 = a$`Sum Sq`[1] / sum(a$`Sum Sq`))
}
lev <- bind_rows(lev); w(lev, "env_variables_by_zone_anova.csv"); print(lev %>% mutate(across(where(is.numeric), ~signif(.x, 3))))
zm <- e %>% group_by(zone) %>% summarise(n_env = n(), xl2_laeq = mean(xl2_laeq), bk_laeq = mean(bk_laeq, na.rm = TRUE), n5 = mean(n5, na.rm = TRUE), naturalness = mean(naturalness), .groups = "drop")
w(zm, "env_means_by_zone.csv")

# ---- 3. does comfort follow level across zones and within zones?
zc <- em_all %>% filter(outcome == "comfort") %>% select(zone, comfort = emmean) %>% inner_join(zm, by = "zone")
cat("\nacross the 8 zones: r(zone comfort, zone XL2 LAeq) =", round(cor(zc$comfort, zc$xl2_laeq), 3), "; r(zone comfort, zone naturalness) =", round(cor(zc$comfort, zc$naturalness), 3), "\n")
dw <- load_p35("wav")
m_int <- fit_backbone("bk_laeq_z * zone", "comfort", dw); m_add <- fit_backbone("bk_laeq_z + zone", "comfort", dw)
cat("within-zone level slopes differ? LRT interaction: p =", signif(anova(m_add, m_int, refit = TRUE)$`Pr(>Chisq)`[2], 3), "\n")
tr <- as_tibble(emtrends(m_int, ~ zone, var = "bk_laeq_z")); w(tr, "comfort_level_slope_by_zone.csv")
print(tr %>% mutate(across(where(is.numeric), ~round(.x, 3))))

# ---- 4. source fingerprints
src_cols <- paste0("src_", SOURCES)
zs <- d %>% group_by(zone) %>% summarise(across(all_of(src_cols), mean), natural = mean(src_natural), human = mean(src_human), mechanical = mean(src_mechanical), .groups = "drop")
fp <- zs %>% pivot_longer(all_of(src_cols), names_to = "source", values_to = "mean") %>% mutate(source = str_remove(source, "src_")) %>%
  group_by(source) %>% mutate(overall = mean(mean), excess = mean - overall) %>% ungroup()
w(fp, "source_fingerprint_long.csv")
top <- fp %>% group_by(zone) %>% slice_max(excess, n = 3) %>% summarise(distinctive = paste0(source, " (+", round(excess, 1), ")", collapse = "; "))
print(top, width = 200)
print(zs %>% select(zone, natural, human, mechanical) %>% mutate(across(where(is.numeric), ~round(.x, 2))))
# appropriateness by zone adjusted for noisiness and comfort
m_app <- fit_backbone("zone + rating_noisiness + comfort", "rating_appropriateness", d)
cat("\nappropriateness ~ zone + noisiness + comfort: noisiness beta", round(fixef(m_app)["rating_noisiness"], 3), " comfort beta", round(fixef(m_app)["comfort"], 3), "\n")
print(as_tibble(emmeans(m_app, ~ zone)) %>% mutate(across(where(is.numeric), ~round(.x, 3))))
