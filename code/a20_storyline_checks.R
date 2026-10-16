# A20 — bounded checks of Richness declared before execution: total salience, diversity, the source groups and the
# natural × mechanical interaction, the natural slope across mechanical salience, and site-held-out prediction with
# perceived and detected sources, under the same sample variants as the core checks.
source(file.path("code", "helpers.R"))
suppressPackageStartupMessages(library(performance))
OUT <- file.path("output", "analysis_20_storyline_checks")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
w <- function(x, f) write_csv(x, file.path(OUT, f))
prepare <- function(d) {
  S <- as.matrix(d[, paste0("src_", SOURCES)]); total <- rowSums(S)
  P <- S / total
  d$total_sal <- total
  d$shannon <- -rowSums(ifelse(P > 0, P * log(P), 0))
  stopifnot(all(total > 0), !anyNA(d$shannon))
  d
}
d <- prepare(load_p35("wav")); da <- prepare(load_p35("all"))
det <- read_frozen_csv("detections_by_env.csv")
join_detection <- function(x, window_name) {
  x %>% mutate(env_id = as.character(env_id)) %>% left_join(det %>% filter(window == window_name) %>% select(env_id, det_natural, det_human, det_mechanical), by = "env_id") %>% mutate(env_id = factor(env_id))
}
d <- join_detection(d, "region")
stopifnot(nrow(d) == 666, all(complete.cases(d[, c("det_natural", "det_human", "det_mechanical")])))
SETS <- list(total = "total_sal", total_diversity = c("total_sal", "shannon"),
             groups = c("src_natural", "src_human", "src_mechanical"),
             groups_diversity = c("src_natural", "src_human", "src_mechanical", "shannon"),
             groups_interaction = c("src_natural", "src_human", "src_mechanical", "src_natural:src_mechanical"),
             perceived_detected = c("src_natural", "src_human", "src_mechanical", "det_natural", "det_human", "det_mechanical"))
CONTRASTS <- tribble(~model, ~base, "total_diversity", "total", "groups_diversity", "groups", "groups_interaction", "groups", "perceived_detected", "groups")
variants <- list(primary = d, no_flagged = filter(d, !flagged), no_trained = filter(d, rater_type == "naive"),
                 without_E13_E28 = filter(d, !evaluator %in% c("E13", "E28")),
                 all_xl2 = da, full_recording_detection = join_detection(prepare(load_p35("wav")), "full"),
                 zone_adjusted = d, salience_curvature = d)
models <- list(); coefficients <- list(); increments <- list()
for (v in names(variants)) {
  x <- variants[[v]]
  physical <- if (v == "all_xl2") "xl2_laeq" else c("bk_laeq", "tonality")
  if (v == "zone_adjusted") physical <- c(physical, "zone")
  if (v == "salience_curvature") physical <- c(physical, "I((total_sal / 100)^2)")
  available <- if (v == "all_xl2") setdiff(names(SETS), "perceived_detected") else names(SETS)
  fits <- list()
  for (s in available) {
    rhs <- paste(c(physical, SETS[[s]]), collapse = " + ")
    fits[[s]] <- fit_backbone(rhs, "richness", x)
    m <- fits[[s]]
    msg <- m@optinfo$conv$lme4$messages
    stopifnot(!any(!grepl("boundary.*singular", msg)))
    models[[paste(v,s)]] <- tibble(variant=v, model=s, n=nobs(m), k=length(physical)+length(SETS[[s]]), r2m=r2_nakagawa(m)$R2_marginal)
    coefficients[[paste(v,s)]] <- as_tibble(summary(m)$coefficients, rownames="term") %>%
      mutate(variant=v, model=s, n=nobs(m), lo=Estimate-qt(0.975,df)*`Std. Error`, hi=Estimate+qt(0.975,df)*`Std. Error`)
    if (v=="primary" && s=="groups_interaction") {
      limits <- quantile(x$src_mechanical, c(0.05,0.95))
      slopes <- as_tibble(emmeans::emtrends(m, ~src_mechanical, var="src_natural",
                          at=list(src_mechanical=seq(limits[1],limits[2],length.out=101))))
      w(slopes, "natural_slope_by_mechanical.csv")
    }
  }
  inc <- bind_rows(lapply(seq_len(nrow(CONTRASTS)), function(i) {
    s <- CONTRASTS$model[i]; b <- CONTRASTS$base[i]
    if (!s %in% available) return(NULL)
    test <- anova(fits[[b]], fits[[s]], refit=TRUE)
    tibble(variant=v, model=s, base=b, n=nrow(x), chisq=test$Chisq[2], df=test$Df[2], p=test$`Pr(>Chisq)`[2])
  })) %>% mutate(q=p.adjust(p,"BH"), family_size=n())
  increments[[v]] <- inc
}
w(bind_rows(models), "models.csv"); w(bind_rows(coefficients), "coefficients.csv"); w(bind_rows(increments), "increments.csv")

# Fixed, outcome-independent source summaries. All regression fits use only the
# seven training sites; no model choice or tuning occurs inside this function.
cv <- function(x, s, physical=c("bk_laeq", "tonality")) {
  pred <- rep(NA_real_, nrow(x)); feats <- c(physical, SETS[[s]])
  for (site in unique(x$site)) {
    train <- x$site != site
    fit <- lm(reformulate(feats, "richness"), data=x[train,])
    pred[!train] <- predict(fit, newdata=x[!train,])
  }
  e <- tibble(env_id=x$env_id, site=x$site, observed=x$richness, predicted=pred) %>%
    group_by(env_id, site) %>% summarise(across(c(observed,predicted),mean), .groups="drop")
  list(score=tibble(model=s, n=nrow(x), n_stops=nrow(e), k=length(feats),
                    r2_record=1-sum((x$richness-pred)^2)/sum((x$richness-mean(x$richness))^2),
                    r2_env_mean=1-sum((e$observed-e$predicted)^2)/sum((e$observed-mean(e$observed))^2)),
       folds=e %>% group_by(site) %>% summarise(n_stops=n(), mse=mean((observed-predicted)^2), .groups="drop") %>% mutate(model=s))
}
scores <- list(); folds <- list()
for (v in names(variants)) {
  if (v %in% c("zone_adjusted", "salience_curvature")) next
  x <- variants[[v]]; physical <- if (v=="all_xl2") "xl2_laeq" else c("bk_laeq","tonality")
  for (s in if (v=="all_xl2") setdiff(names(SETS),"perceived_detected") else names(SETS)) {
    result <- cv(x,s,physical)
    scores[[paste(v,s)]] <- result$score %>% mutate(variant=v)
    folds[[paste(v,s)]] <- result$folds %>% mutate(variant=v)
  }
}
w(bind_rows(scores), "cv_scores.csv"); w(bind_rows(folds), "cv_site_errors.csv")
print(bind_rows(models) %>% filter(variant=="primary"))
print(bind_rows(increments), n=30)
print(bind_rows(scores) %>% select(variant, model, r2_env_mean) %>% pivot_wider(names_from=variant,values_from=r2_env_mean), width=180)
