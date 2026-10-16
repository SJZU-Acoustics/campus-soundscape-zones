# A08 — person-level variables (noise sensitivity, WHO-5, rater type, sex, education) as main effects and as moderators of the level slope.
# 30 evaluators: person main effects are tests on 30 person means (evaluator RE); moderation is a cross-level interaction.
source(file.path("code", "helpers.R"))
suppressPackageStartupMessages(library(performance))
OUT <- file.path("output", "analysis_08_person_level")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
w <- function(x, f) write_csv(x, file.path(OUT, f))
prep <- function(x) x %>% mutate(education = factor(education), sex = factor(sex),
                                 ns_z = as.numeric(scale(ns_sum)), who5_z = as.numeric(scale(who5_sum)))
d <- prep(load_p35("wav")); dall <- prep(load_p35("all"))
OUTS <- c("comfort", "richness", "rating_noisiness", "rating_appropriateness")
PV <- c("ns_z", "who5_z", "rater_type", "sex", "education")
term_of <- function(pv) switch(pv, ns_z = "ns_z", who5_z = "who5_z", rater_type = "rater_typetrained", sex = paste0("sex", levels(d$sex)[2]), education = paste0("education", levels(d$education)[2]))

# ---- 1. main effects (family 20), level adjusted, 666 records; replicated on 990 with XL2
main <- expand_grid(outcome = OUTS, pv = PV) %>% rowwise() %>% mutate(
  fit = list(fit_backbone(paste0("bk_laeq_z + ", pv), outcome, d)),
  beta = fixef(fit)[term_of(pv)], p = summary(fit)$coefficients[term_of(pv), 5],
  fit_all = list(fit_backbone(paste0("xl2_laeq_z + ", pv), outcome, dall)),
  beta_990 = fixef(fit_all)[term_of(pv)], p_990 = summary(fit_all)$coefficients[term_of(pv), 5]) %>% ungroup() %>% select(-fit, -fit_all) %>% bh_family(p_col = "p")
w(main, "person_main_effects.csv"); print(main %>% mutate(across(c(beta, beta_990), ~round(.x, 3)), across(c(p, q, p_990), ~signif(.x, 2))) %>% arrange(q), n = 20)

# ---- 2. moderation of the level slope (family 20)
mod <- expand_grid(outcome = OUTS, pv = PV) %>% rowwise() %>% mutate(
  fit = list(fit_backbone(paste0("bk_laeq_z * ", pv), outcome, d)),
  beta_int = fixef(fit)[paste0("bk_laeq_z:", term_of(pv))], p_int = summary(fit)$coefficients[paste0("bk_laeq_z:", term_of(pv)), 5],
  fit_all = list(fit_backbone(paste0("xl2_laeq_z * ", pv), outcome, dall)),
  beta_int_990 = fixef(fit_all)[paste0("xl2_laeq_z:", term_of(pv))], p_int_990 = summary(fit_all)$coefficients[paste0("xl2_laeq_z:", term_of(pv)), 5]) %>% ungroup() %>% select(-fit, -fit_all) %>% bh_family(p_col = "p_int")
w(mod, "person_moderation_of_level_slope.csv"); print(mod %>% mutate(across(c(beta_int, beta_int_990), ~round(.x, 3)), across(c(p_int, q, p_int_990), ~signif(.x, 2))) %>% arrange(q), n = 20)

# ---- 3. noise sensitivity × source type (does NS sharpen the mechanical penalty?) — 4 outcomes × 2 = 8
nsx <- expand_grid(outcome = OUTS, src = c("src_mechanical", "src_natural")) %>% rowwise() %>% mutate(
  fit = list(fit_backbone(paste0("bk_laeq_z + ns_z * ", src), outcome, d)),
  beta_int = fixef(fit)[paste0("ns_z:", src)], p_int = summary(fit)$coefficients[paste0("ns_z:", src), 5]) %>% ungroup() %>% select(-fit) %>% bh_family(p_col = "p_int")
w(nsx, "ns_by_source_interaction.csv"); print(nsx %>% mutate(beta_int = round(beta_int, 3), across(c(p_int, q), ~signif(.x, 2))))

# ---- 4. evaluator random slopes for level: how heterogeneous are the slopes? (LRT reference here; permutation null in A09)
rs <- list()
for (y in OUTS) {
  m0 <- lme4::lmer(as.formula(paste0(y, " ~ bk_laeq_z + (1 | evaluator) + (1 | env_id)")), data = d, REML = FALSE)
  m1 <- lme4::lmer(as.formula(paste0(y, " ~ bk_laeq_z + (1 + bk_laeq_z | evaluator) + (1 | env_id)")), data = d, REML = FALSE)
  vc <- as.data.frame(VarCorr(m1)); sd_slope <- vc$sdcor[vc$grp == "evaluator" & vc$var1 == "bk_laeq_z" & is.na(vc$var2)]
  cor_is <- vc$sdcor[vc$grp == "evaluator" & !is.na(vc$var2)]
  lrt <- 2 * (logLik(m1) - logLik(m0)); p_mix <- 0.5 * pchisq(lrt, 1, lower.tail = FALSE) + 0.5 * pchisq(lrt, 2, lower.tail = FALSE)
  bl <- ranef(m1)$evaluator$bk_laeq_z + fixef(m1)["bk_laeq_z"]
  rs[[y]] <- tibble(outcome = y, mean_slope = fixef(m1)["bk_laeq_z"], sd_slope = sd_slope, cor_int_slope = cor_is, lrt = as.numeric(lrt), p_mixture = p_mix,
                    blup_slope_min = min(bl), blup_slope_max = max(bl), n_positive_slopes = sum(bl > 0))
}
rs <- bind_rows(rs); w(rs, "evaluator_random_slopes.csv"); print(rs %>% mutate(across(where(is.numeric), ~signif(.x, 3))), width = 200)

# ---- 5. person-level descriptives of the ratings: mean comfort per evaluator vs NS, WHO-5 (30 points)
pe <- d %>% group_by(evaluator, ns_sum, who5_sum, rater_type, sex) %>% summarise(across(all_of(OUTS), mean), n = n(), .groups = "drop")
w(pe, "evaluator_means.csv")
cat("person-level r: comfort~NS", round(cor(pe$comfort, pe$ns_sum), 2), "| noisiness~NS", round(cor(pe$rating_noisiness, pe$ns_sum), 2), "| comfort~WHO5", round(cor(pe$comfort, pe$who5_sum), 2), "| richness~WHO5", round(cor(pe$richness, pe$who5_sum), 2), "\n")
