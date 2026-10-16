# A15 — what richness is: count, salient count, total salience, maximum, Shannon, Simpson, evenness, dominance of the 24-source salience
# vector → richness (level adjusted; crossed RE), combinations, LOSO/LOEO recording-mean R², comfort as contrast, response-style check.
source(file.path("code", "helpers.R"))
suppressPackageStartupMessages(library(performance))
OUT <- file.path("output", "analysis_15_richness_diversity"); w <- function(x, f) write_csv(x, file.path(OUT, f))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
src_cols <- paste0("src_", SOURCES)
idx <- function(x) {
  S <- as.matrix(x[, src_cols]); tot <- rowSums(S); P <- S / ifelse(tot > 0, tot, NA)
  H <- -rowSums(ifelse(P > 0, P * log(P), 0), na.rm = TRUE); nsrc <- rowSums(S > 0)  # NB: not `n` — the table has a loudness column `n`
  x %>% mutate(n_sources = nsrc, n_salient = rowSums(S >= 5), total_sal = tot, max_sal = apply(S, 1, max), shannon = ifelse(tot > 0, H, NA), simpson = ifelse(tot > 0, 1 - rowSums(P^2, na.rm = TRUE), NA),
               evenness = ifelse(nsrc > 1, H / log(nsrc), NA), dominance = ifelse(tot > 0, max_sal / tot, NA))
}
d <- idx(load_p35("wav")); dall <- idx(load_p35("all"))
IDX <- c("n_sources", "n_salient", "total_sal", "max_sal", "shannon", "simpson", "evenness", "dominance")
cat("records with zero total salience:", sum(d$total_sal == 0), "(wav) /", sum(dall$total_sal == 0), "(all); evenness NA:", sum(is.na(d$evenness)), "\n")
print(d %>% summarise(across(all_of(IDX), list(mean = ~mean(.x, na.rm = TRUE), sd = ~sd(.x, na.rm = TRUE)))) %>% pivot_longer(everything(), names_to = c("index", "stat"), names_pattern = "(.*)_(mean|sd)$") %>% pivot_wider(names_from = stat, values_from = value) %>% mutate(across(c(mean, sd), ~round(.x, 2))))
cm <- cor(d[, IDX], use = "pairwise"); w(as.data.frame(round(cm, 2)) %>% mutate(index = IDX) %>% relocate(index), "index_correlations.csv")
dz <- d %>% mutate(across(all_of(IDX), ~as.numeric(scale(.x)), .names = "{.col}_z")); dallz <- dall %>% mutate(across(all_of(IDX), ~as.numeric(scale(.x)), .names = "{.col}_z"))
# ---- 1. single-index screen, richness and comfort (family 8 each)
scan <- expand_grid(outcome = c("richness", "comfort"), index = IDX) %>% rowwise() %>% mutate(
  fit = list(fit_backbone(paste0("bk_laeq_z + ", index, "_z"), outcome, dz)), n = nobs(fit),
  beta = fixef(fit)[3], p = summary(fit)$coefficients[3, 5], r2m = r2_nakagawa(fit)$R2_marginal, aic = AIC(refitML(fit))) %>% ungroup() %>% select(-fit) %>%
  group_by(outcome) %>% mutate(q = p.adjust(p, "BH"), family_size = n()) %>% ungroup()
w(scan, "single_index_screen.csv"); print(scan %>% mutate(across(c(beta, r2m), ~round(.x, 3)), aic = round(aic, 1), across(c(p, q), ~signif(.x, 2))) %>% arrange(outcome, desc(r2m)), n = 16, width = 150)
# ---- 2. combinations for richness: is diversity informative beyond the count, and vice versa? total salience as the rating-intensity alternative
combos <- c(count = "n_sources_z", shannon = "shannon_z", count_shannon = "n_sources_z + shannon_z", count_total = "n_sources_z + total_sal_z", count_evenness = "n_sources_z + evenness_z",
            count_max = "n_sources_z + max_sal_z", groups3 = "src_natural + src_human + src_mechanical", groups3_count = "src_natural + src_human + src_mechanical + n_sources_z", groups3_shannon = "src_natural + src_human + src_mechanical + shannon_z")
dd <- dz %>% filter(!is.na(evenness_z))
cb <- tibble(model = names(combos), rhs = combos) %>% rowwise() %>% mutate(fit = list(fit_backbone(paste0("bk_laeq_z + ", rhs), "richness", dd)), r2m = r2_nakagawa(fit)$R2_marginal, aic = AIC(refitML(fit)),
  terms = paste(sprintf("%s %.3f (p %.2g)", names(fixef(fit))[-(1:2)], fixef(fit)[-(1:2)], summary(fit)$coefficients[-(1:2), 5]), collapse = "; ")) %>% ungroup()
coefficients <- bind_rows(lapply(seq_len(nrow(cb)), function(i) {
  co <- as_tibble(summary(cb$fit[[i]])$coefficients, rownames = "term")
  co %>% mutate(model = cb$model[i], n = nobs(cb$fit[[i]]),
                lo = Estimate - qt(0.975, df) * `Std. Error`,
                hi = Estimate + qt(0.975, df) * `Std. Error`)
}))
w(coefficients, "richness_combination_coefficients.csv")
cb <- cb %>% mutate(n = vapply(fit, nobs, numeric(1))) %>% select(-fit)
w(cb, "richness_combinations.csv"); print(cb %>% mutate(r2m = round(r2m, 3), aic = round(aic, 1)) %>% select(model, r2m, aic, terms), n = 12, width = 250)
# ---- 3. cross-validated (lm) recording-mean R² for richness: new site / new person
cv <- function(feats, fold_var, data) { folds <- unique(data[[fold_var]]); pred <- rep(NA_real_, nrow(data))
  for (f in folds) { tr <- data[[fold_var]] != f; m <- lm(reformulate(feats, "richness"), data = data[tr, ]); pred[!tr] <- predict(m, newdata = data[!tr, ]) }
  env <- tibble(env_id = data$env_id, y = data$richness, p = pred) %>% group_by(env_id) %>% summarise(y = mean(y), p = mean(p), .groups = "drop")
  c(record = 1 - sum((data$richness - pred)^2) / sum((data$richness - mean(data$richness))^2), env_mean = 1 - sum((env$y - env$p)^2) / sum((env$y - mean(env$y))^2)) }
FS <- list(count = c("bk_laeq", "n_sources"), shannon = c("bk_laeq", "shannon"), count_shannon = c("bk_laeq", "n_sources", "shannon"), groups3 = c("bk_laeq", "src_natural", "src_human", "src_mechanical"), groups3_shannon = c("bk_laeq", "src_natural", "src_human", "src_mechanical", "shannon"))
cvres <- expand_grid(family = names(FS), scheme = c("site", "evaluator")) %>% rowwise() %>% mutate(r = list(cv(FS[[family]], scheme, dd)), r2_record = r[["record"]], r2_env_mean = r[["env_mean"]]) %>% ungroup() %>% select(-r)
w(cvres, "cv_richness.csv"); print(cvres %>% mutate(across(c(r2_record, r2_env_mean), ~round(.x, 3))) %>% pivot_wider(names_from = scheme, values_from = c(r2_record, r2_env_mean)))
# ---- 4. response-style check: between-person association of mean total salience / n_sources with mean richness (30 points), and within-person slopes
pe <- d %>% group_by(evaluator) %>% summarise(richness = mean(richness), total_sal = mean(total_sal), n_sources = mean(n_sources), shannon = mean(shannon, na.rm = TRUE), .groups = "drop")
cat("between-person r: richness~total_sal", round(cor(pe$richness, pe$total_sal), 2), "| richness~n_sources", round(cor(pe$richness, pe$n_sources), 2), "| richness~shannon", round(cor(pe$richness, pe$shannon), 2), "\n")
dw <- dz %>% group_by(evaluator) %>% mutate(across(c(n_sources_z, shannon_z, total_sal_z), ~.x - mean(.x, na.rm = TRUE), .names = "{.col}_w")) %>% ungroup() %>% filter(!is.na(shannon_z_w))
mw <- fit_backbone("bk_laeq_z + n_sources_z_w + shannon_z_w + total_sal_z_w", "richness", dw); print(round(summary(mw)$coefficients[-1, c(1, 5)], 4))
# ---- 5. replication on 990 (XL2 level) of the single-index screen for richness
rep990 <- tibble(index = IDX) %>% rowwise() %>% mutate(fit = list(fit_backbone(paste0("xl2_laeq_z + ", index, "_z"), "richness", dallz)),
  beta = fixef(fit)[3], p = summary(fit)$coefficients[3, 5], r2m = r2_nakagawa(fit)$R2_marginal) %>% ungroup() %>% select(-fit)
w(rep990, "single_index_screen_990_xl2.csv"); print(rep990 %>% mutate(across(c(beta, r2m), ~round(.x, 3)), p = signif(p, 2)) %>% arrange(desc(r2m)))
