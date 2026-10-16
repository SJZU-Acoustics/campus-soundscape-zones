# A16b — Comfort and Richness of the four ICSAQ quadrants (raw means of the 826 classified ratings, and EMMs from the
# level + zone + quadrant model as in A16), used to compare the Comfort of the monotonous and chaotic quadrants.
source(file.path("code", "helpers.R"))
suppressPackageStartupMessages(library(emmeans))
OUT <- file.path("output", "analysis_16_quadrants")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
QL <- c("calm", "vibrant", "monotonous", "chaotic")
# Means that are exactly zero can come out as +/-1e-17 in floating point; values within 1e-9 of zero are set to zero first.
quad <- function(c, r) { c[abs(c) < 1e-9] <- 0; r[abs(r) < 1e-9] <- 0; factor(case_when(c == 0 | r == 0 ~ NA_character_, c > 0 & r < 0 ~ "calm", c > 0 & r > 0 ~ "vibrant", c < 0 & r < 0 ~ "monotonous", TRUE ~ "chaotic"), levels = QL) }
d <- load_p35("all") %>% mutate(quadrant = quad(comfort, richness)) %>% filter(!is.na(quadrant))
raw <- d %>% group_by(quadrant) %>% summarise(n = n(), comfort_mean = mean(comfort), comfort_sd = sd(comfort), richness_mean = mean(richness), richness_sd = sd(richness), .groups = "drop")
emm <- bind_rows(lapply(c("comfort", "richness"), function(y) { m <- fit_backbone("xl2_laeq_z + zone + quadrant", y, d)
  as_tibble(emmeans(m, "quadrant")) %>% transmute(quadrant, outcome = y, emmean, se = SE) }))
res <- raw %>% left_join(emm %>% filter(outcome == "comfort") %>% select(quadrant, comfort_emm = emmean, comfort_emm_se = se), by = "quadrant") %>%
  left_join(emm %>% filter(outcome == "richness") %>% select(quadrant, richness_emm = emmean, richness_emm_se = se), by = "quadrant")
write_csv(res, file.path(OUT, "comfort_by_quadrant.csv")); print(res %>% mutate(across(where(is.numeric), ~round(.x, 3))), width = 200)
