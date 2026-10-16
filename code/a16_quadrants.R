# A16 — ISO-style quadrants from the signs of comfort and richness (calm: C>0,R<0; vibrant: C>0,R>0; monotonous: C<0,R<0; chaotic: C<0,R>0).
# Exact zeros on either axis are excluded in the primary analysis (counted) and folded into the ≤0 side as a sensitivity check.
source(file.path("code", "helpers.R"))
suppressPackageStartupMessages({library(performance); library(emmeans)})
OUT <- file.path("output", "analysis_16_quadrants"); w <- function(x, f) write_csv(x, file.path(OUT, f))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
QL <- c("calm", "vibrant", "monotonous", "chaotic")
# Means that are exactly zero can come out as +/-1e-17 in floating point; values within 1e-9 of zero are set to zero first.
quad <- function(c, r, strict = TRUE) { c[abs(c) < 1e-9] <- 0; r[abs(r) < 1e-9] <- 0
  cp <- if (strict) c > 0 else c > 0; rp <- if (strict) r > 0 else r > 0
  factor(case_when(strict & (c == 0 | r == 0) ~ NA_character_, cp & !rp ~ "calm", cp & rp ~ "vibrant", !cp & !rp ~ "monotonous", TRUE ~ "chaotic"), levels = QL) }
d <- load_p35("all") %>% mutate(quadrant = quad(comfort, richness), quadrant_all = quad(comfort, richness, strict = FALSE))
cat("records:", nrow(d), "| exact zero on comfort", sum(d$comfort == 0), "on richness", sum(d$richness == 0), "| excluded (either)", sum(is.na(d$quadrant)), "\n")
print(table(d$quadrant, useNA = "ifany")); print(round(100 * prop.table(table(d$quadrant)), 1))
# ---- 1. distribution by zone (records, and recordings by co-walker majority)
tz <- d %>% filter(!is.na(quadrant)) %>% count(zone, quadrant) %>% group_by(zone) %>% mutate(pct = round(100 * n / sum(n), 1)) %>% ungroup()
w(tz, "quadrant_by_zone_records.csv"); print(tz %>% select(zone, quadrant, pct) %>% pivot_wider(names_from = quadrant, values_from = pct, values_fill = 0))
ct <- chisq.test(table(d$zone[!is.na(d$quadrant)], d$quadrant[!is.na(d$quadrant)])); cat("zone × quadrant chi-square", round(ct$statistic, 1), "df", ct$parameter, "p", signif(ct$p.value, 2), "\n")
env_q <- d %>% group_by(env_id, zone) %>% summarise(comfort = mean(comfort), richness = mean(richness), .groups = "drop") %>% mutate(quadrant = quad(comfort, richness, strict = FALSE))
w(env_q, "quadrant_by_recording.csv"); print(env_q %>% count(zone, quadrant) %>% group_by(zone) %>% mutate(pct = round(100 * n / sum(n))) %>% select(-n) %>% pivot_wider(names_from = quadrant, values_from = pct, values_fill = 0))
# ---- 2. appropriateness by quadrant within zone
dq <- d %>% filter(!is.na(quadrant))
m0 <- fit_backbone("xl2_laeq_z + zone + quadrant", "rating_appropriateness", dq); m1 <- fit_backbone("xl2_laeq_z + zone * quadrant", "rating_appropriateness", dq)
cat("quadrant main effect: "); print(anova(m0)["quadrant", ]); cat("zone × quadrant LRT p", signif(anova(m0, m1, refit = TRUE)$`Pr(>Chisq)`[2], 2), "\n")
em_all <- as_tibble(emmeans(m0, "quadrant")); print(em_all %>% mutate(across(where(is.numeric), ~round(.x, 2))))
em <- as_tibble(emmeans(m1, ~ quadrant | zone)); w(em, "appropriateness_emmeans_quadrant_by_zone.csv")
print(em %>% select(zone, quadrant, emmean) %>% mutate(emmean = round(emmean, 2)) %>% pivot_wider(names_from = quadrant, values_from = emmean))
best <- em %>% group_by(zone) %>% slice_max(emmean, n = 1) %>% select(zone, best_quadrant = quadrant, emmean); print(best %>% mutate(emmean = round(emmean, 2)))
vc <- as_tibble(contrast(emmeans(m1, ~ quadrant | zone), method = list(vibrant_minus_calm = c(-1, 1, 0, 0)))) %>% mutate(p_holm = p.adjust(p.value, "holm"))
w(vc, "vibrant_minus_calm_by_zone.csv"); print(vc %>% select(zone, estimate, SE, p.value, p_holm) %>% mutate(across(c(estimate, SE), ~round(.x, 2)), across(c(p.value, p_holm), ~signif(.x, 2))))
# sensitivity: zeros folded in
m1b <- fit_backbone("xl2_laeq_z + zone * quadrant_all", "rating_appropriateness", d)
vcb <- as_tibble(contrast(emmeans(m1b, ~ quadrant_all | zone), method = list(vibrant_minus_calm = c(-1, 1, 0, 0)))) %>% mutate(p_holm = p.adjust(p.value, "holm"))
print(vcb %>% select(zone, estimate, p_holm) %>% mutate(estimate = round(estimate, 2), p_holm = signif(p_holm, 2)))
# ---- 3. noisiness and the auditory-environment rating by quadrant (check the quadrant labels behave)
for (y in c("rating_noisiness", "rating_auditory_env")) { m <- fit_backbone("xl2_laeq_z + zone + quadrant", y, dq); e <- as_tibble(emmeans(m, "quadrant")); cat(y, ":", paste(e$quadrant, round(e$emmean, 2), collapse = " | "), "\n") }
# ---- 4. addendum: overall pairwise quadrant contrasts (Holm) for appropriateness — is monotonous worse than chaotic?
pw <- as_tibble(pairs(emmeans(m0, "quadrant"), adjust = "holm")); w(pw, "quadrant_pairwise_overall.csv"); print(pw %>% mutate(across(c(estimate, SE), ~round(.x, 2)), p.value = signif(p.value, 2)) %>% select(contrast, estimate, SE, p.value))
mc <- fit_backbone("xl2_laeq_z + zone + quadrant", "comfort", dq); ec <- as_tibble(emmeans(mc, "quadrant")); cat("comfort by quadrant:", paste(ec$quadrant, round(ec$emmean, 2), collapse = " | "), "\n")
# monotonous vs chaotic within zone
mcz <- as_tibble(contrast(emmeans(m1, ~ quadrant | zone), method = list(chaotic_minus_monotonous = c(0, 0, -1, 1)))) %>% mutate(p_holm = p.adjust(p.value, "holm"))
w(mcz, "chaotic_minus_monotonous_by_zone.csv"); print(mcz %>% select(zone, estimate, SE, p.value, p_holm) %>% mutate(across(c(estimate, SE), ~round(.x, 2)), across(c(p.value, p_holm), ~signif(.x, 2))))
