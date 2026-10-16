# A11 — which feature family generalises: cross-validated prediction of the four outcomes from acoustic, source, visual and combined
# feature sets, leave-one-site-out (new place) and leave-one-evaluator-out (new person). Fixed-effect (lm) predictions on held-out records;
# out-of-sample R² at the record level and at the recording-mean level (evaluator noise averaged out).
source(file.path("code", "helpers.R"))
OUT <- file.path("output", "analysis_11_feature_families")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
w <- function(x, f) write_csv(x, file.path(OUT, f))
d <- load_p35("wav") %>% mutate(n_sources = rowSums(across(all_of(paste0("src_", SOURCES))) > 0))
OUTS <- c("comfort", "richness", "rating_noisiness", "rating_appropriateness")
FAM <- list(level_only = "bk_laeq", acoustic = c("bk_laeq", "tonality", "s", "bk_dl5_95", "n5_n95", "r", "f"), level_tonality = c("bk_laeq", "tonality"),
            source_groups = c("src_natural", "src_human", "src_mechanical"), sources_24 = paste0("src_", SOURCES), n_sources = "n_sources",
            visual = c("naturalness", "vdi", "vri", "vis_people", "vis_vehicles", "vis_sky"),
            acoustic_visual = c("bk_laeq", "tonality", "s", "naturalness", "vdi", "vri", "vis_people", "vis_vehicles"),
            acoustic_sourcegroups = c("bk_laeq", "tonality", "src_natural", "src_human", "src_mechanical"),
            all = c("bk_laeq", "tonality", "s", "bk_dl5_95", "src_natural", "src_human", "src_mechanical", "naturalness", "vdi", "vri", "vis_people", "vis_vehicles"))
cv <- function(y, feats, fold_var) {
  folds <- unique(d[[fold_var]]); pred <- rep(NA_real_, nrow(d))
  for (f in folds) { tr <- d[[fold_var]] != f; m <- lm(reformulate(feats, y), data = d[tr, ]); pred[!tr] <- predict(m, newdata = d[!tr, ]) }
  r2_rec <- 1 - sum((d[[y]] - pred)^2) / sum((d[[y]] - mean(d[[y]]))^2)
  env <- tibble(env_id = d$env_id, y = d[[y]], p = pred) %>% group_by(env_id) %>% summarise(y = mean(y), p = mean(p), .groups = "drop")
  r2_env <- 1 - sum((env$y - env$p)^2) / sum((env$y - mean(env$y))^2)
  c(r2_record = r2_rec, r2_env_mean = r2_env)
}
res <- expand_grid(outcome = OUTS, family = names(FAM), scheme = c("site", "evaluator")) %>% rowwise() %>%
  mutate(k = length(FAM[[family]]), r = list(cv(outcome, FAM[[family]], scheme)), r2_record = r[["r2_record"]], r2_env_mean = r[["r2_env_mean"]]) %>% ungroup() %>% select(-r)
w(res, "cv_r2_by_family.csv")
print(res %>% mutate(across(c(r2_record, r2_env_mean), ~round(.x, 3))) %>% pivot_wider(names_from = scheme, values_from = c(r2_record, r2_env_mean)) %>% arrange(outcome, desc(r2_env_mean_site)), n = 60, width = 200)
# in-sample reference (lm on all 666) for the same families
ins <- expand_grid(outcome = OUTS, family = names(FAM)) %>% rowwise() %>% mutate(r2_insample = summary(lm(reformulate(FAM[[family]], outcome), data = d))$r.squared) %>% ungroup()
w(ins, "insample_r2_by_family.csv")
