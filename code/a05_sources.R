# A05 — perceived source salience → the four outcomes, adjusting for level (666 WAV records, crossed RE).
source(file.path("code", "helpers.R"))
suppressPackageStartupMessages(library(performance))
OUT <- file.path("output", "analysis_05_sources_to_perception"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
w <- function(x, f) write_csv(x, file.path(OUT, f))
d <- load_p35("wav") %>% mutate(n_sources = rowSums(across(all_of(paste0("src_", SOURCES))) > 0),
                                 across(c(src_natural, src_human, src_mechanical, n_sources), ~as.numeric(scale(.x)), .names = "{.col}_z"))
dall <- load_p35("all") %>% mutate(n_sources = rowSums(across(all_of(paste0("src_", SOURCES))) > 0),
                                    across(c(src_natural, src_human, src_mechanical, n_sources), ~as.numeric(scale(.x)), .names = "{.col}_z"))
OUTS <- c("comfort", "richness", "rating_noisiness", "rating_appropriateness")

# ---- 1. the three source groups together, adjusting for level (12 coefficients = one family)
grp <- list()
for (y in OUTS) {
  m <- fit_backbone("bk_laeq_z + src_natural_z + src_human_z + src_mechanical_z", y, d)
  co <- as_tibble(summary(m)$coefficients, rownames = "term") %>% filter(term != "(Intercept)") %>% mutate(outcome = y, r2m = r2_nakagawa(m)$R2_marginal,
        r2m_level_only = r2_nakagawa(fit_backbone("bk_laeq_z", y, d))$R2_marginal)
  grp[[y]] <- co
}
grp <- bind_rows(grp); grp_fam <- grp %>% filter(term != "bk_laeq_z") %>% bh_family(p_col = "Pr(>|t|)")
w(bind_rows(grp_fam, grp %>% filter(term == "bk_laeq_z")), "source_groups_models.csv")
print(grp_fam %>% mutate(across(c(Estimate, `Std. Error`, r2m, r2m_level_only), ~round(.x, 3)), p = signif(`Pr(>|t|)`, 2), q = signif(q, 2)) %>% select(outcome, term, Estimate, p, q, r2m, r2m_level_only), n = 12)
# whole-sample replication with XL2 level
grp_all <- list()
for (y in OUTS) {
  m <- fit_backbone("xl2_laeq_z + src_natural_z + src_human_z + src_mechanical_z", y, dall)
  grp_all[[y]] <- as_tibble(summary(m)$coefficients, rownames = "term") %>% filter(term != "(Intercept)") %>% mutate(outcome = y)
}
grp_all <- bind_rows(grp_all); w(grp_all, "source_groups_models_all990_xl2.csv")
print(grp_all %>% filter(term != "xl2_laeq_z") %>% mutate(Estimate = round(Estimate, 3), p = signif(`Pr(>|t|)`, 2)) %>% select(outcome, term, Estimate, p), n = 12)

# ---- 2. single-source screen adjusting for level (family = 4 × 24 = 96)
single <- expand_grid(outcome = OUTS, source = SOURCES) %>% rowwise() %>% mutate(
  fit = list(fit_backbone(paste0("bk_laeq_z + src_", source), outcome, d)),
  beta_per_point = fixef(fit)[3], se = summary(fit)$coefficients[3, 2], p = summary(fit)$coefficients[3, 5]) %>% ungroup() %>% select(-fit) %>% bh_family(p_col = "p")
w(single, "single_source_screen.csv")
cat("\nsingle-source screen: q<0.05 in", sum(single$q < 0.05), "of", nrow(single), "\n")
print(single %>% filter(q < 0.05) %>% mutate(across(c(beta_per_point, se), ~round(.x, 3)), p = signif(p, 2), q = signif(q, 2)) %>% arrange(outcome, q), n = 60)

# ---- 3. number of sources heard → richness (the "eventful" hypothesis), adjusting for level
for (y in OUTS) {
  m <- fit_backbone("bk_laeq_z + n_sources_z", y, d)
  cat(y, ": n_sources beta", round(fixef(m)[3], 3), "p", signif(summary(m)$coefficients[3, 5], 2), "| r2m", round(r2_nakagawa(m)$R2_marginal, 3), "\n")
}
# ---- 4. birds: the campus's dominant source — dose response with level adjusted (comfort, appropriateness)
for (y in c("comfort", "rating_appropriateness", "richness")) {
  m <- fit_backbone("bk_laeq_z + src_birds + src_vehicles + src_conversation", y, d)
  co <- summary(m)$coefficients
  cat(y, ": birds", round(co["src_birds", 1], 3), "(p", signif(co["src_birds", 5], 2), ") vehicles", round(co["src_vehicles", 1], 3), "(p", signif(co["src_vehicles", 5], 2), ") conversation", round(co["src_conversation", 1], 3), "(p", signif(co["src_conversation", 5], 2), ")\n")
}
