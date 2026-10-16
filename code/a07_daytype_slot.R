# A07 — weekday vs weekend within the common slots (09-12 … 18-21) and time-of-day effects on level, sources and ratings.
# Design caveat: each session has its own evaluator group, so daytype/slot are confounded with evaluator cohorts —
# the evaluator RE absorbs person means; a session RE is added as a check (14 sessions).
source(file.path("code", "helpers.R"))
suppressPackageStartupMessages({library(performance); library(emmeans)})
OUT <- file.path("output", "analysis_07_daytype_slot")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
w <- function(x, f) write_csv(x, file.path(OUT, f))
d <- load_p35("all")
common <- d %>% distinct(slot, daytype) %>% count(slot) %>% filter(n == 2) %>% pull(slot) %>% as.character()
dc <- d %>% filter(slot %in% common) %>% droplevels()
cat("common slots:", common, "| records", nrow(dc), "| weekday", sum(dc$daytype == "weekday"), "weekend", sum(dc$daytype == "weekend"), "\n")
OUTS <- c("comfort", "richness", "rating_noisiness", "rating_appropriateness")

# ---- 1. environment level: level and source profile by daytype and slot (env-level, common slots)
e <- dc %>% group_by(env_id, site, zone, slot, daytype, session) %>%
  summarise(xl2_laeq = first(xl2_laeq), bk_laeq = first(bk_laeq), src_natural = mean(src_natural), src_human = mean(src_human), src_mechanical = mean(src_mechanical),
            n_sources = mean(rowSums(across(all_of(paste0("src_", SOURCES))) > 0)), .groups = "drop")
env_tests <- list()
for (v in c("xl2_laeq", "bk_laeq", "src_natural", "src_human", "src_mechanical", "n_sources")) {
  m <- lm(reformulate(c("zone", "slot", "daytype"), v), data = e)
  a <- anova(m); em <- emmeans(m, "daytype") %>% as_tibble()
  env_tests[[v]] <- tibble(variable = v, n = sum(!is.na(e[[v]])), p_daytype = a["daytype", "Pr(>F)"], p_slot = a["slot", "Pr(>F)"],
                           weekday = em$emmean[em$daytype == "weekday"], weekend = em$emmean[em$daytype == "weekend"],
                           slot_emm = paste(sprintf("%s %.2f", levels(e$slot), as_tibble(emmeans(m, "slot"))$emmean), collapse = "; "))
}
env_tests <- bind_rows(env_tests) %>% mutate(q_daytype = p.adjust(p_daytype, "BH"), q_slot = p.adjust(p_slot, "BH"), family_size = n())
w(env_tests, "env_level_daytype_slot.csv"); print(env_tests %>% mutate(across(c(weekday, weekend), ~round(.x, 2)), across(c(p_daytype, p_slot, q_daytype, q_slot), ~signif(.x, 2))), width = 200)

# ---- 2. ratings: daytype and slot with zone fixed, evaluator + env RE (common slots); session RE as check
rt <- list()
for (y in OUTS) {
  m <- fit_backbone("zone + slot + daytype", y, dc)
  a <- anova(m)
  emd <- as_tibble(emmeans(m, "daytype")); ems <- as_tibble(emmeans(m, "slot"))
  ms <- lmerTest::lmer(as.formula(paste0(y, " ~ zone + slot + daytype + (1 | evaluator) + (1 | env_id) + (1 | session)")), data = dc)
  as_ <- anova(ms)
  rt[[y]] <- tibble(outcome = y, p_daytype = a["daytype", "Pr(>F)"], p_slot = a["slot", "Pr(>F)"],
                    weekday = emd$emmean[emd$daytype == "weekday"], weekend = emd$emmean[emd$daytype == "weekend"],
                    slot_emm = paste(sprintf("%s %.2f", ems$slot, ems$emmean), collapse = "; "),
                    p_daytype_sessionRE = as_["daytype", "Pr(>F)"], p_slot_sessionRE = as_["slot", "Pr(>F)"],
                    var_session = as.data.frame(VarCorr(ms))$vcov[as.data.frame(VarCorr(ms))$grp == "session"])
}
rt <- bind_rows(rt) %>% mutate(q_daytype = p.adjust(p_daytype, "BH"), q_slot = p.adjust(p_slot, "BH"), family_size = n())
w(rt, "ratings_daytype_slot.csv"); print(rt %>% mutate(across(c(weekday, weekend, var_session), ~round(.x, 3)), across(starts_with("p_"), ~signif(.x, 2)), across(starts_with("q_"), ~signif(.x, 2))), width = 220)

# ---- 3. does the weekday–weekend difference in ratings survive adjusting for level and sources? (mechanism check)
for (y in OUTS) {
  m <- fit_backbone("zone + slot + daytype + xl2_laeq_z + src_natural + src_human + src_mechanical", y, dc)
  emd <- as_tibble(emmeans(m, "daytype"))
  cat(y, ": daytype adjusted for level+sources — weekday", round(emd$emmean[1], 3), "weekend", round(emd$emmean[2], 3), "p", signif(anova(m)["daytype", "Pr(>F)"], 2), "\n")
}
# ---- 4. slot effect across the full weekday set (5 slots incl. 06-09), level and comfort
dw <- d %>% filter(daytype == "weekday")
ew <- dw %>% distinct(env_id, zone, slot, xl2_laeq)
print(as_tibble(emmeans(lm(xl2_laeq ~ zone + slot, ew), "slot")) %>% mutate(across(where(is.numeric), ~round(.x, 2))))
mw <- fit_backbone("zone + slot", "comfort", dw); print(as_tibble(emmeans(mw, "slot")) %>% mutate(across(where(is.numeric), ~round(.x, 3)))); cat("weekday slot effect on comfort p", signif(anova(mw)["slot", "Pr(>F)"], 2), "\n")
# ---- 5. by-zone daytype contrasts (which zones change at the weekend?) — comfort and level, descriptive
zc <- dc %>% group_by(zone, daytype) %>% summarise(comfort = mean(comfort), richness = mean(richness), xl2_laeq = mean(xl2_laeq), n = n(), .groups = "drop") %>%
  pivot_wider(names_from = daytype, values_from = c(comfort, richness, xl2_laeq, n)) %>%
  mutate(d_comfort = comfort_weekend - comfort_weekday, d_richness = richness_weekend - richness_weekday, d_laeq = xl2_laeq_weekend - xl2_laeq_weekday)
w(zc, "zone_by_daytype_descriptives.csv"); print(zc %>% select(zone, starts_with("d_"), n_weekday, n_weekend) %>% mutate(across(where(is.numeric), ~round(.x, 2))))
