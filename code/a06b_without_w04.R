# A06b — the visual screen, audio-visual interactions and congruence re-run without session W04 (5 July, 77 records with
# recordings), whose 16 panoramas were taken on 6 and 7 July rather than during the walk.
source(file.path("code", "helpers.R"))
OUT <- file.path("output", "analysis_06_visual_to_perception"); w <- function(x, f) write_csv(x, file.path(OUT, f))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
VIS <- c("naturalness", "artificiality", "vdi", "vri", "vis_people", "vis_vehicles", "vis_sky", "vis_vegetation")
d0 <- load_p35("wav") %>% mutate(across(all_of(c(VIS, "src_natural", "src_human", "src_mechanical")), ~as.numeric(scale(.x)), .names = "{.col}_z"))
d <- d0 %>% filter(session != "W04"); cat("records:", nrow(d0), "->", nrow(d), "without W04\n")
OUTS <- c("comfort", "richness", "rating_noisiness", "rating_appropriateness")
scan <- expand_grid(outcome = OUTS, visual = VIS) %>% rowwise() %>% mutate(
  fit_b = list(fit_backbone(paste0("bk_laeq_z + ", visual, "_z"), outcome, d)), beta_bw = fixef(fit_b)[3], p_bw = summary(fit_b)$coefficients[3, 5],
  fit_w = list(fit_backbone(paste0("bk_laeq_z + zone + ", visual, "_z"), outcome, d)), beta_within = fixef(fit_w)[paste0(visual, "_z")], p_within = summary(fit_w)$coefficients[paste0(visual, "_z"), 5]) %>%
  ungroup() %>% select(-fit_b, -fit_w) %>% mutate(q_bw = p.adjust(p_bw, "BH"), q_within = p.adjust(p_within, "BH"), family_size = n())
full <- read_csv(file.path(OUT, "visual_screen.csv"), show_col_types = FALSE) %>% select(outcome, visual, beta_bw_all = beta_bw, q_bw_all = q_bw, beta_within_all = beta_within, q_within_all = q_within)
cmp <- scan %>% left_join(full, by = c("outcome", "visual")); w(cmp, "visual_screen_without_w04.csv")
print(cmp %>% mutate(across(c(beta_bw, beta_bw_all, beta_within, beta_within_all), ~round(.x, 3)), across(c(q_bw, q_bw_all, q_within, q_within_all), ~signif(.x, 2))) %>%
  select(outcome, visual, beta_bw_all, beta_bw, q_bw_all, q_bw, beta_within_all, beta_within, q_within_all, q_within) %>% arrange(outcome, q_bw_all), n = 32, width = 200)
pairs <- tribble(~heard, ~seen, "src_mechanical", "vis_vehicles", "src_human", "vis_people", "src_natural", "naturalness")
cong <- expand_grid(outcome = OUTS, pairs) %>% rowwise() %>% mutate(fit = list(fit_backbone(paste0("bk_laeq_z + zone + ", heard, "_z * ", seen, "_z"), outcome, d)),
  beta_int = fixef(fit)[paste0(heard, "_z:", seen, "_z")], p_int = summary(fit)$coefficients[paste0(heard, "_z:", seen, "_z"), 5]) %>% ungroup() %>% select(-fit) %>% mutate(q = p.adjust(p_int, "BH"), family_size = n())
w(cong, "congruence_without_w04.csv"); print(cong %>% mutate(beta_int = round(beta_int, 3), p_int = signif(p_int, 2), q = signif(q, 2)), n = 12, width = 150)
inter <- expand_grid(outcome = OUTS, visual = c("naturalness", "artificiality")) %>% rowwise() %>% mutate(fit = list(fit_backbone(paste0("bk_laeq_z * ", visual, "_z + zone"), outcome, d)),
  beta_int = fixef(fit)[paste0("bk_laeq_z:", visual, "_z")], p_int = summary(fit)$coefficients[paste0("bk_laeq_z:", visual, "_z"), 5]) %>% ungroup() %>% select(-fit) %>% mutate(q = p.adjust(p_int, "BH"))
w(inter, "audio_visual_interaction_without_w04.csv"); print(inter %>% mutate(beta_int = round(beta_int, 3), p_int = signif(p_int, 2), q = signif(q, 2)))
