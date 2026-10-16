# A19 — is the a priori natural / human / mechanical grouping the structure the data show?
# PCA and hierarchical clustering of the 24 source saliences (record level, 990; recording level, 208), agreement with the a priori
# groups (adjusted Rand index), within- vs between-group coherence, and a head-to-head of the three a priori group means against
# three data-driven components as predictors (R²m, AIC, LOSO/LOEO recording-mean R²).
# Predicted before running: the first data-driven axis separates natural from mechanical sources; human sources split between axes
# rather than forming their own; three components do not beat the three a priori groups out of sample.
source(file.path("code", "helpers.R"))
suppressPackageStartupMessages(library(performance))
OUT <- file.path("output", "analysis_19_source_structure"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
w <- function(x, f) write_csv(x, file.path(OUT, f))
src_cols <- paste0("src_", SOURCES); grp <- rep(names(SOURCE_GROUPS), lengths(SOURCE_GROUPS)); names(grp) <- unlist(SOURCE_GROUPS)
d <- load_p35("all"); dw <- load_p35("wav")
S <- as.matrix(d[, src_cols]); colnames(S) <- SOURCES
# ---- 1. correlation structure and coherence
R <- cor(S); w(as.data.frame(round(R, 2)) %>% mutate(source = SOURCES) %>% relocate(source), "source_correlations.csv")
coh <- expand_grid(a = SOURCES, b = SOURCES) %>% filter(a < b) %>% mutate(r = R[cbind(a, b)], same = grp[a] == grp[b], pair = ifelse(same, grp[a], paste(sort(c(grp[a], grp[b]))[1], sort(c(grp[a], grp[b]))[2], sep = "–"))) %>% rowwise() %>% mutate(pair = ifelse(same, grp[a], paste(sort(c(grp[a], grp[b])), collapse = "–"))) %>% ungroup()
cs <- coh %>% group_by(pair, same) %>% summarise(n_pairs = n(), mean_r = mean(r), .groups = "drop") %>% arrange(desc(same), desc(mean_r)); w(cs, "group_coherence.csv")
cat("== mean inter-source correlation within / between a priori groups ==\n"); print(cs %>% mutate(mean_r = round(mean_r, 3)))
# ---- 2. PCA (record level) and hierarchical clustering of sources
pc <- prcomp(S, scale. = TRUE); ve <- summary(pc)$importance[2, 1:5]; cat("PCA variance PC1-5:", round(ve, 3), "\n")
load <- as_tibble(round(pc$rotation[, 1:4], 2)) %>% mutate(source = SOURCES, group = grp[SOURCES]) %>% relocate(source, group) %>% arrange(group, desc(PC1)); w(load, "pca_loadings.csv"); print(load, n = 24)
ari <- function(a, b) { t <- table(a, b); n <- sum(t); s <- function(x) sum(choose(x, 2)); si <- sum(choose(t, 2)); sa <- s(rowSums(t)); sb <- s(colSums(t)); e <- sa * sb / choose(n, 2); (si - e) / ((sa + sb) / 2 - e) }
hc <- hclust(as.dist(1 - R), method = "ward.D2"); cl <- list()
for (k in 2:5) { m <- cutree(hc, k); cl[[as.character(k)]] <- tibble(k = k, source = SOURCES, group = grp[SOURCES], cluster = m); cat(sprintf("k = %d: ARI vs a priori groups %.3f | ", k, ari(m, grp[SOURCES]))); for (j in 1:k) cat(sprintf("[%d: %s] ", j, paste(SOURCES[m == j], collapse = " "))); cat("\n") }
w(bind_rows(cl), "hclust_memberships.csv")
# a "data-driven 3" assignment: k = 3 Ward clusters; and the PCA sign-split of PC1
# ---- 3. head-to-head as predictors: a priori 3 groups vs PC1–3 scores vs Ward-3 cluster means (record-level; 666 WAV records for the level layer)
Sw <- as.matrix(dw[, src_cols]); colnames(Sw) <- SOURCES; sc <- predict(pc, newdata = Sw)[, 1:3]; m3 <- cutree(hc, 3)
dw <- dw %>% mutate(pc1 = sc[, 1], pc2 = sc[, 2], pc3 = sc[, 3], cl1 = rowMeans(Sw[, m3 == 1, drop = FALSE]), cl2 = rowMeans(Sw[, m3 == 2, drop = FALSE]), cl3 = rowMeans(Sw[, m3 == 3, drop = FALSE]))
OUTS <- c("comfort", "richness", "rating_noisiness", "rating_appropriateness")
SETS <- list(apriori3 = c("src_natural", "src_human", "src_mechanical"), pca3 = c("pc1", "pc2", "pc3"), ward3 = c("cl1", "cl2", "cl3"), sources24 = src_cols)
hh <- expand_grid(outcome = OUTS, set = names(SETS)) %>% rowwise() %>% mutate(fit = list(fit_backbone(paste(c("bk_laeq_z", SETS[[set]]), collapse = " + "), outcome, dw)), r2m = r2_nakagawa(fit)$R2_marginal, aic = AIC(refitML(fit)), k = length(SETS[[set]])) %>% ungroup() %>% select(-fit)
w(hh, "head_to_head_insample.csv"); cat("\n== in-sample: level + set (R²m, AIC) ==\n"); print(hh %>% mutate(r2m = round(r2m, 3), aic = round(aic, 1)) %>% pivot_wider(names_from = set, values_from = c(r2m, aic, k)), width = 200)
cv <- function(y, set, fold_var, data) { folds <- unique(data[[fold_var]]); pred <- rep(NA_real_, nrow(data))
  feats <- c("bk_laeq", "tonality", SETS[[set]])
  for (f in folds) {
    tr <- data[[fold_var]] != f
    fold_data <- data
    # Learned representations must not see the held-out site or evaluator.
    # The pooled PCA/clustering above remains a descriptive analysis only.
    if (set %in% c("pca3", "ward3")) {
      mat <- as.matrix(data[, src_cols])
      keep <- apply(mat[tr, , drop = FALSE], 2, sd) > 0
      mat <- mat[, keep, drop = FALSE]
      stopifnot(ncol(mat) >= 3)
      if (set == "pca3") {
        transform <- prcomp(mat[tr, , drop = FALSE], scale. = TRUE)
        scores <- predict(transform, mat)[, 1:3, drop = FALSE]
        fold_data[, c("pc1", "pc2", "pc3")] <- as_tibble(scores)
      } else {
        membership <- cutree(hclust(as.dist(1 - cor(mat[tr, , drop = FALSE])), method = "ward.D2"), 3)
        for (j in 1:3) fold_data[[paste0("cl", j)]] <- rowMeans(mat[, membership == j, drop = FALSE])
      }
    }
    m <- lm(reformulate(feats, y), data = fold_data[tr, ])
    pred[!tr] <- predict(m, newdata = fold_data[!tr, ])
  }
  env <- tibble(env_id = data$env_id, y = data[[y]], p = pred) %>% group_by(env_id) %>% summarise(y = mean(y), p = mean(p), .groups = "drop")
  1 - sum((env$y - env$p)^2) / sum((env$y - mean(env$y))^2) }
cvres <- expand_grid(outcome = OUTS, set = names(SETS), scheme = c("site", "evaluator")) %>% rowwise() %>% mutate(r2_env_mean = cv(outcome, set, scheme, dw)) %>% ungroup() %>% mutate(preprocessing = "training fold only")
w(cvres, "cv_r2_by_set.csv"); cat("\n== LOSO / LOEO recording-mean R² (level + tonality + set) ==\n"); print(cvres %>% mutate(r2_env_mean = round(r2_env_mean, 3)) %>% pivot_wider(names_from = c(set, scheme), values_from = r2_env_mean), width = 200)
# ---- 4. is "human" a group at all? drop it / merge it: comfort and richness models with natural + mechanical only, and with human split into speech-like vs activity-like
alt <- tibble(model = c("nat+mech", "nat+hum+mech", "nat+mech+speech+activity"), rhs = c("src_natural + src_mechanical", "src_natural + src_human + src_mechanical", "src_natural + src_mechanical + src_conversation + I((src_children + src_playing + src_footsteps + src_running)/4)"))
for (y in c("comfort", "richness", "rating_appropriateness")) { r <- alt %>% rowwise() %>% mutate(fit = list(fit_backbone(paste0("bk_laeq_z + ", rhs), y, dw)), r2m = r2_nakagawa(fit)$R2_marginal, aic = AIC(refitML(fit))) %>% ungroup() %>% select(-fit); cat("\n", y, ": "); print(r %>% mutate(r2m = round(r2m, 3), aic = round(aic, 1)) %>% select(model, r2m, aic)) }
# ---- 5. recording-level structure (208 env means) — does the three-group structure hold where the ratings are averaged?
E <- d %>% group_by(env_id) %>% summarise(across(all_of(src_cols), mean), .groups = "drop") %>% select(-env_id) %>% as.matrix(); colnames(E) <- SOURCES
Re <- cor(E); hce <- hclust(as.dist(1 - Re), method = "ward.D2"); me <- cutree(hce, 3); cat("\nrecording level (208): k = 3 ARI", round(ari(me, grp[SOURCES]), 3), "|"); for (j in 1:3) cat(sprintf(" [%d: %s]", j, paste(SOURCES[me == j], collapse = " "))); cat("\n")
pce <- prcomp(E, scale. = TRUE); cat("recording-level PCA variance PC1-3:", round(summary(pce)$importance[2, 1:3], 3), "\n"); w(as_tibble(round(pce$rotation[, 1:3], 2)) %>% mutate(source = SOURCES, group = grp[SOURCES]) %>% relocate(source, group), "pca_loadings_recording_level.csv")
