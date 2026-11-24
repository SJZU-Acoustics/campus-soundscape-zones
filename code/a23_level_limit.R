# A23 — stop levels against a reference level of 55 dB(A), the daytime limit that the Chinese national standard
# GB 3096-2008 sets for Class 1 areas, by site, beside appropriateness. The comparison is descriptive and indicative:
# a stop measurement lasted about a minute, whereas the standard assesses long-term day and night levels.
source(file.path("code", "helpers.R"))
OUT <- file.path("output", "analysis_23_level_limit")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
LIMIT <- 55  # dB(A); every survey slot falls in the day period
d <- load_p35("all"); env <- load_p35_env(); sites <- load_p35_sites()
stopifnot(nrow(d) == 990, nrow(env) == 208, all(!is.na(env$xl2_laeq)))

stop_means <- d %>% group_by(env_id) %>%
  summarise(site = first(site), zone = first(zone), n_raters = n(),
            approp = mean(rating_appropriateness), comfort = mean(comfort),
            richness = mean(richness), .groups = "drop") %>%
  left_join(env %>% mutate(env_id = factor(env_id)) %>% select(env_id, xl2_laeq), by = "env_id") %>%
  mutate(above = xl2_laeq > LIMIT)
stopifnot(nrow(stop_means) == 208)

# stops above the reference level and mean ratings, by site (Table 1)
by_site <- stop_means %>% group_by(site, zone) %>%
  summarise(n_stops = n(), laeq_median = median(xl2_laeq), laeq_min = min(xl2_laeq), laeq_max = max(xl2_laeq),
            n_above = sum(above), share_above = mean(above),
            approp_mean = mean(approp), comfort_mean = mean(comfort), richness_mean = mean(richness), .groups = "drop") %>%
  left_join(sites %>% mutate(site = factor(site)) %>% select(site, site_name), by = "site") %>%
  relocate(site_name, .after = site)
write_csv(by_site, file.path(OUT, "limit_by_site.csv"))

overall <- stop_means %>% summarise(n_stops = n(), n_above = sum(above), share_above = mean(above),
                                    approp_above = mean(approp[above]), approp_below = mean(approp[!above]))
write_csv(overall, file.path(OUT, "limit_overall.csv"))

# ratings at stops above the reference level against ratings at the other stops
dd <- d %>% left_join(stop_means %>% select(env_id, xl2_laeq, above), by = "env_id")
m_all <- lmerTest::lmer(rating_appropriateness ~ above + (1 | evaluator) + (1 | env_id), data = dd)
contr <- broom.mixed::tidy(m_all, effects = "fixed", conf.int = TRUE) %>% filter(term == "aboveTRUE") %>%
  mutate(model = "appropriateness, unadjusted")
write_csv(contr, file.path(OUT, "limit_contrast.csv"))
