# A06 — visual composition → the four outcomes, adjusting for level; audio-visual interaction and congruence.
# Identification: visual composition is ~83% a zone property (A03), so every effect is tested twice —
# between+within (no zone term) and within-zone (zone fixed factor).
source(file.path("code", "helpers.R"))
suppressPackageStartupMessages(library(performance))
OUT <- file.path("output", "analysis_06_visual_to_perception"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
w <- function(x, f) write_csv(x, file.path(OUT, f))
VIS <- c("naturalness", "artificiality", "vdi", "vri", "vis_people", "vis_vehicles", "vis_sky", "vis_vegetation")
d <- load_p35("wav") %>% mutate(across(all_of(c(VIS, "src_natural", "src_human", "src_mechanical")), ~as.numeric(scale(.x)), .names = "{.col}_z"))
OUTS <- c("comfort", "richness", "rating_noisiness", "rating_appropriateness")

# ---- 1. visual screen (family = 4 × 8 = 32), between+within and within-zone
scan <- expand_grid(outcome = OUTS, visual = VIS) %>% rowwise() %>% mutate(
  fit_b = list(fit_backbone(paste0("bk_laeq_z + ", visual, "_z"), outcome, d)),
  beta_bw = fixef(fit_b)[3], p_bw = summary(fit_b)$coefficients[3, 5],
  fit_w = list(fit_backbone(paste0("bk_laeq_z + zone + ", visual, "_z"), outcome, d)),
  beta_within = fixef(fit_w)[paste0(visual, "_z")], p_within = summary(fit_w)$coefficients[paste0(visual, "_z"), 5]) %>% ungroup() %>% select(-fit_b, -fit_w)
scan <- scan %>% mutate(q_bw = p.adjust(p_bw, "BH"), q_within = p.adjust(p_within, "BH"), family_size = n())
w(scan, "visual_screen.csv")
print(scan %>% mutate(across(c(beta_bw, beta_within), ~round(.x, 3)), across(c(p_bw, q_bw, p_within, q_within), ~signif(.x, 2))) %>% arrange(outcome, q_bw), n = 40, width = 160)
# within-zone variance share of the visual variables (how much is there to identify?)
vs <- load_p35_env() %>% mutate(zone = factor(ZONE_LEVELS[site], levels = ZONE_LEVELS))
share <- tibble(visual = VIS) %>% rowwise() %>% mutate(within_share = 1 - summary(lm(reformulate("zone", visual), data = vs))$r.squared) %>% ungroup()
w(share, "visual_within_zone_variance_share.csv"); print(share %>% mutate(within_share = round(within_share, 2)))

# ---- 2. audio-visual interaction: naturalness × level, artificiality × level (family = 4 × 2 = 8)
inter <- expand_grid(outcome = OUTS, visual = c("naturalness", "artificiality")) %>% rowwise() %>% mutate(
  fit = list(fit_backbone(paste0("bk_laeq_z * ", visual, "_z + zone"), outcome, d)),
  beta_int = fixef(fit)[paste0("bk_laeq_z:", visual, "_z")], p_int = summary(fit)$coefficients[paste0("bk_laeq_z:", visual, "_z"), 5]) %>% ungroup() %>% select(-fit) %>%
  mutate(q = p.adjust(p_int, "BH"), family_size = n())
w(inter, "audio_visual_interaction.csv"); print(inter %>% mutate(beta_int = round(beta_int, 3), p_int = signif(p_int, 2), q = signif(q, 2)))

# ---- 3. congruence: sources heard × matching visual class (family = 4 outcomes × 3 pairs = 12), within zone, level adjusted
pairs <- tribble(~heard, ~seen, "src_mechanical", "vis_vehicles", "src_human", "vis_people", "src_natural", "naturalness")
cong <- expand_grid(outcome = OUTS, pairs) %>% rowwise() %>% mutate(
  fit = list(fit_backbone(paste0("bk_laeq_z + zone + ", heard, "_z * ", seen, "_z"), outcome, d)),
  beta_heard = fixef(fit)[paste0(heard, "_z")], beta_seen = fixef(fit)[paste0(seen, "_z")],
  beta_int = fixef(fit)[paste0(heard, "_z:", seen, "_z")], p_int = summary(fit)$coefficients[paste0(heard, "_z:", seen, "_z"), 5]) %>% ungroup() %>% select(-fit) %>%
  mutate(q = p.adjust(p_int, "BH"), family_size = n())
w(cong, "congruence.csv"); print(cong %>% mutate(across(c(beta_heard, beta_seen, beta_int), ~round(.x, 3)), p_int = signif(p_int, 2), q = signif(q, 2)), n = 12, width = 150)

# ---- 4. the visual-environment rating itself: what drives it, and does it carry into comfort beyond sound?
mv <- fit_backbone("zone + naturalness_z + vdi_z", "rating_visual_env", d)
print(round(summary(mv)$coefficients[c("naturalness_z", "vdi_z"), ], 3))
mc <- fit_backbone("bk_laeq_z + zone + rating_visual_env", "comfort", d); mc0 <- fit_backbone("bk_laeq_z + zone", "comfort", d)
cat("comfort ~ level + zone + visual rating: visual-rating beta", round(fixef(mc)["rating_visual_env"], 3), "p", signif(summary(mc)$coefficients["rating_visual_env", 5], 2),
    "| R2m", round(r2_nakagawa(mc0)$R2_marginal, 3), "->", round(r2_nakagawa(mc)$R2_marginal, 3), "\n")
