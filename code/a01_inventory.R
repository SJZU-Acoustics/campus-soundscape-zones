# A01 — inventory of the sample: structure, outcome distributions, adjective
# structure, source-salience landscape by zone, acoustic and visual descriptives.
source(file.path("code", "helpers.R"))
OUT <- file.path("output", "analysis_01_inventory"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
d <- load_p35("all"); e <- load_p35_env(); p <- load_p35_people()
w <- function(x, f) write_csv(x, file.path(OUT, f))

# ---- 1. structure
w(d %>% count(zone, site, daytype) %>% pivot_wider(names_from = daytype, values_from = n, values_fill = 0) %>%
    mutate(total = weekday + weekend), "records_by_site_daytype.csv")
w(d %>% count(slot, daytype) %>% pivot_wider(names_from = daytype, values_from = n, values_fill = 0), "records_by_slot_daytype.csv")
w(d %>% group_by(session, date, daytype, slot_planned) %>% summarise(n = n(), n_wav = sum(wav_recorded), n_eval = n_distinct(evaluator), n_env = n_distinct(env_id), .groups = "drop"), "records_by_session.csv")
w(d %>% group_by(evaluator, cohort, rater_type) %>% summarise(n = n(), n_wav = sum(wav_recorded), n_sites = n_distinct(site), n_sessions = n_distinct(session), .groups = "drop"), "records_by_evaluator.csv")
cat("records", nrow(d), "| wav", sum(d$wav_recorded), "| envs", n_distinct(d$env_id), "| evaluators", n_distinct(d$evaluator), "\n")
cat("records per evaluator: median", median(table(d$evaluator)), "range", range(table(d$evaluator)), "\n")
cat("records per env:", paste(names(table(table(d$env_id))), table(table(d$env_id)), sep = "×", collapse = " "), "\n")

# ---- 2. outcomes
outs <- c("comfort", "richness", "rating_appropriateness", "rating_noisiness", "rating_auditory_env", "rating_visual_env")
desc <- d %>% summarise(across(all_of(outs), list(mean = ~mean(.x), sd = ~sd(.x), min = ~min(.x), max = ~max(.x)))) %>%
  pivot_longer(everything(), names_to = c("var", "stat"), names_pattern = "(.*)_(mean|sd|min|max)$") %>%
  pivot_wider(names_from = stat, values_from = value)
floorceil <- d %>% summarise(across(all_of(outs[3:6]), list(floor = ~mean(.x == 1), ceiling = ~mean(.x == 5)))) %>%
  pivot_longer(everything(), names_to = c("var", "stat"), names_pattern = "(.*)_(floor|ceiling)$") %>% pivot_wider(names_from = stat, values_from = value)
w(left_join(desc, floorceil, by = "var"), "outcome_descriptives.csv"); print(left_join(desc, floorceil, by = "var"), n = 10)
cm <- cor(d[, outs]); w(as.data.frame(cm) %>% mutate(var = rownames(cm)) %>% relocate(var), "outcome_correlations.csv")
cat("\nr(comfort, richness) =", round(cm["comfort", "richness"], 3), "| r(comfort, appropriateness) =", round(cm["comfort", "rating_appropriateness"], 3),
    "| r(comfort, noisiness) =", round(cm["comfort", "rating_noisiness"], 3), "| r(richness, noisiness) =", round(cm["richness", "rating_noisiness"], 3), "\n")

# ---- 3. adjective structure (circumplex check)
adj <- d %>% select(all_of(paste0("adj_", ADJECTIVES)))
pc <- prcomp(adj, scale. = TRUE)
ve <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)
load <- as.data.frame(round(pc$rotation[, 1:3], 3)) %>% mutate(adjective = ADJECTIVES) %>% relocate(adjective)
w(load, "adjective_pca_loadings.csv")
cat("\nadjective PCA variance %:", ve, "\n"); print(load)
cat("r(PC1, comfort) =", round(cor(pc$x[, 1], d$comfort), 3), " r(PC1, richness) =", round(cor(pc$x[, 1], d$richness), 3),
    " r(PC2, comfort) =", round(cor(pc$x[, 2], d$comfort), 3), " r(PC2, richness) =", round(cor(pc$x[, 2], d$richness), 3), "\n")
w(as.data.frame(round(cor(adj), 2)) %>% mutate(adjective = ADJECTIVES) %>% relocate(adjective), "adjective_correlations.csv")

# ---- 4. source-salience landscape
src_cols <- paste0("src_", SOURCES)
zone_src <- d %>% group_by(zone) %>% summarise(across(all_of(src_cols), mean), .groups = "drop")
w(zone_src, "source_salience_by_zone.csv")
prev <- d %>% summarise(across(all_of(src_cols), list(mean = ~mean(.x), prevalence = ~mean(.x > 0), share_ge5 = ~mean(.x >= 5)))) %>%
  pivot_longer(everything(), names_to = c("source", "stat"), names_pattern = "src_(.*)_(mean|prevalence|share_ge5)$") %>% pivot_wider(names_from = stat, values_from = value) %>% arrange(desc(mean))
w(prev, "source_prevalence.csv"); print(prev, n = 24)
top3 <- zone_src %>% pivot_longer(-zone, names_to = "source", values_to = "mean") %>% group_by(zone) %>% slice_max(mean, n = 3) %>%
  summarise(top3 = paste0(str_remove(source, "src_"), " (", round(mean, 1), ")", collapse = "; "))
print(top3)
grp <- d %>% group_by(zone) %>% summarise(natural = mean(src_natural), human = mean(src_human), mechanical = mean(src_mechanical), n_sources_heard = mean(rowSums(across(all_of(src_cols)) > 0)), .groups = "drop")
w(grp, "source_groups_by_zone.csv"); print(grp)
cat("mean number of sources heard (>0) per record:", round(mean(rowSums(d[, src_cols] > 0)), 2), "; other-sound named on", sum(d$other_sound_none == 0), "records\n")

# ---- 5. acoustic and visual descriptives by zone (env level)
ez <- e %>% mutate(zone = factor(ZONE_LEVELS[site], levels = ZONE_LEVELS)) %>% group_by(zone) %>%
  summarise(n_env = n(), n_wav = sum(wav_recorded), xl2_sd = sd(xl2_laeq), xl2_laeq = mean(xl2_laeq), bk_laeq = mean(bk_laeq, na.rm = TRUE),
            n5 = mean(n5, na.rm = TRUE), s = mean(s, na.rm = TRUE), naturalness = mean(naturalness), artificiality = mean(artificiality), vri = mean(vri), vdi = mean(vdi),
            vis_people = mean(vis_people), vis_vehicles = mean(vis_vehicles), .groups = "drop")
w(ez, "env_descriptives_by_zone.csv"); print(ez %>% mutate(across(where(is.numeric), ~round(.x, 2))), width = 200)
cat("\nenv-level XL2 LAeq: mean", round(mean(e$xl2_laeq), 1), "sd", round(sd(e$xl2_laeq), 1), "range", range(e$xl2_laeq), "\n")
cat("visual naturalness range", round(range(e$naturalness), 2), "; artificiality range", round(range(e$artificiality), 2), "; VRI table:", paste(names(table(e$vri)), table(e$vri), sep = "×", collapse = " "), "\n")
# person table summary
cat("\npeople: sex", paste(table(p$sex), collapse = "/"), "| age band", paste(names(table(p$age_band)), table(p$age_band), collapse = " / "), " | education", paste(table(p$education), collapse = "/"),
    "| trained", sum(p$rater_type == "trained"), "| WHO-5 mean", round(mean(p$who5_sum), 1), "low", sum(p$who5_low == 1), "| NS mean", round(mean(p$ns_sum), 1), "sd", round(sd(p$ns_sum), 1), "\n")
