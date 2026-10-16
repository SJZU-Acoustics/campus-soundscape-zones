# A12 — robustness of the headline model (outcome ~ LAeq + tonality + natural + human + mechanical, evaluator + env RE):
# full 666; without flagged rows; without trained raters; with session RE; XL2 level instead of B&K; all 990 with XL2; all 990 without 2 July;
# and cumulative-link mixed models for the two 5-point items. Zone omnibus repeated on the subsets.
source(file.path("code", "helpers.R"))
suppressPackageStartupMessages({library(performance); library(ordinal)})
OUT <- file.path("output", "analysis_12_robustness")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
w <- function(x, f) write_csv(x, file.path(OUT, f))
zs <- function(x) x %>% mutate(across(c(src_natural, src_human, src_mechanical), ~as.numeric(scale(.x)), .names = "{.col}_z"))
d <- zs(load_p35("wav")); dall <- zs(load_p35("all"))
OUTS <- c("comfort", "richness", "rating_noisiness", "rating_appropriateness")
RHS_BK <- "bk_laeq_z + tonality_z + src_natural_z + src_human_z + src_mechanical_z"; RHS_XL <- "xl2_laeq_z + src_natural_z + src_human_z + src_mechanical_z"
variants <- list(
  full_666 = list(data = d, rhs = RHS_BK), no_flagged = list(data = filter(d, !flagged), rhs = RHS_BK), no_trained = list(data = filter(d, rater_type == "naive"), rhs = RHS_BK),
  session_RE = list(data = d, rhs = RHS_BK, extra = " + (1 | session)"), xl2_on_666 = list(data = d, rhs = RHS_XL),
  all_990_xl2 = list(data = dall, rhs = RHS_XL), no_0702_xl2 = list(data = filter(dall, date != "2026-07-02"), rhs = RHS_XL), no_flagged_990 = list(data = filter(dall, !flagged), rhs = RHS_XL))
rows <- list()
for (vn in names(variants)) for (y in OUTS) {
  v <- variants[[vn]]; extra <- if (is.null(v$extra)) "" else v$extra
  m <- lmerTest::lmer(as.formula(paste0(y, " ~ ", v$rhs, " + (1 | evaluator) + (1 | env_id)", extra)), data = v$data)
  co <- summary(m)$coefficients; terms <- rownames(co)[-1]
  rows[[paste(vn, y)]] <- as_tibble(co[-1, , drop = FALSE], rownames = "term") %>% mutate(variant = vn, outcome = y, n = nrow(v$data), r2m = r2_nakagawa(m)$R2_marginal) %>%
    mutate(term = str_replace(term, "^(bk|xl2)_laeq_z$", "level_z"))
}
rows <- bind_rows(rows); w(rows, "headline_model_variants.csv")
print(rows %>% filter(term %in% c("level_z", "tonality_z", "src_natural_z", "src_mechanical_z")) %>% mutate(b = sprintf("%.3f%s", Estimate, ifelse(`Pr(>|t|)` < 0.001, "***", ifelse(`Pr(>|t|)` < 0.01, "**", ifelse(`Pr(>|t|)` < 0.05, "*", ""))))) %>%
        select(outcome, term, variant, b) %>% pivot_wider(names_from = variant, values_from = b), n = 20, width = 250)
print(rows %>% distinct(variant, outcome, n, r2m) %>% mutate(r2m = round(r2m, 3)) %>% pivot_wider(names_from = variant, values_from = r2m) %>% select(-n), width = 250)
# ---- zone omnibus on the subsets
zo <- list()
for (vn in c("full_666", "no_flagged", "no_trained", "all_990_xl2", "no_0702_xl2")) for (y in OUTS) {
  dd <- variants[[vn]]$data; m1 <- fit_backbone("zone", y, dd); m0 <- fit_backbone("1", y, dd)
  zo[[paste(vn, y)]] <- tibble(variant = vn, outcome = y, n = nrow(dd), lrt_p = anova(m0, m1, refit = TRUE)$`Pr(>Chisq)`[2], r2m = r2_nakagawa(m1)$R2_marginal)
}
zo <- bind_rows(zo); w(zo, "zone_omnibus_variants.csv"); print(zo %>% mutate(lrt_p = signif(lrt_p, 2), r2m = round(r2m, 3)) %>% pivot_wider(names_from = variant, values_from = c(lrt_p, r2m)) %>% select(-n), width = 250)
# ---- ordinal check for the two 5-point items (CLMM, 666)
for (y in c("rating_noisiness", "rating_appropriateness")) {
  dd <- d %>% mutate(yo = factor(.data[[y]], ordered = TRUE))
  cm <- clmm(as.formula(paste0("yo ~ ", RHS_BK, " + (1 | evaluator) + (1 | env_id)")), data = dd, link = "logit")
  co <- coef(summary(cm)); co <- co[grepl("_z$", rownames(co)), , drop = FALSE]
  cat("\nCLMM", y, "\n"); print(round(co[, c(1, 4)], 3))
  w(as_tibble(co, rownames = "term") %>% mutate(outcome = y), paste0("clmm_", y, ".csv"))
}
