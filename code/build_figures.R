# =============================================================================
# build_figures.R — Figures 1-6 of the paper.
# Reads the deposited tables and the analysis result tables in output/; the only
# model fitted here is the Figure 3a two-axis model, whose estimates are written
# to output/f3_richness_by_zone.csv for build_tables.R.
# Figure 1a draws the OpenStreetMap base layer in data/osm/ (see its README);
# without that file Figure 1 is skipped and Figures 2-6 are still built.
# =============================================================================
source(file.path("code", "style.R"))

sites <- load_p35_sites() %>% mutate(zone = factor(zone, levels = ZONE_LEVELS))
sess <- load_p35_sessions()
d_all <- load_p35("all")
zone_y <- function() scale_y_discrete(limits = rev(ZONE_LEVELS), labels = ZONE_LAB)

# =============================================================================
# Fig. 1 — sites and design: campus map (a) and session grid (b)
# =============================================================================
OSM_FILE <- file.path("data", "osm", "osm_campus_features.csv")
if (!file.exists(OSM_FILE)) {
  message("Figure 1 skipped: ", OSM_FILE, " not found (run python3 code/fetch_osm.py; see data/osm/README.md).")
} else {
osm <- readr::read_csv(OSM_FILE, show_col_types = FALSE) %>%
  mutate(gid = paste(osm_id, member_ref))
LAT0 <- mean(sites$lat_wgs84); ASP <- 1 / cos(LAT0 * pi / 180)
# Frame: the same box as the session grid in (b) (about 3:4), campus centred and
# set a little north so the legend strip above the arterial road stays clear of
# sites 1 and 2; the extra width shows the road and the blocks around the campus.
yr <- c(41.7354, 41.7473); xr <- 123.5107 + c(-1, 1) * (diff(yr) * ASP / 0.754) / 2
poly <- function(cls) filter(osm, class == cls)
f1a <- ggplot() +
  geom_polygon(data = poly("green"), aes(lon, lat, group = gid), fill = "#E3EFD9", colour = NA) +
  geom_polygon(data = poly("sport"), aes(lon, lat, group = gid), fill = "#F6E3CF", colour = NA) +
  geom_polygon(data = poly("water"), aes(lon, lat, group = gid), fill = "#CFE3F2", colour = NA) +
  geom_polygon(data = poly("building"), aes(lon, lat, group = gid), fill = "grey78", colour = NA) +
  geom_path(data = poly("road_minor"), aes(lon, lat, group = gid), colour = "white", linewidth = 0.45) +
  geom_path(data = poly("road_minor"), aes(lon, lat, group = gid), colour = "grey55", linewidth = 0.18) +
  geom_path(data = poly("road_major"), aes(lon, lat, group = gid), colour = "grey35", linewidth = 0.6) +
  geom_path(data = poly("campus"), aes(lon, lat, group = gid), colour = "grey20", linewidth = 0.35,
            linetype = "22") +
  geom_point(data = sites, aes(lon_wgs84, lat_wgs84, colour = zone, shape = zone), size = 2.6) +
  geom_text(data = sites, aes(lon_wgs84, lat_wgs84, label = site), family = "Helvetica",
            size = 2.6, nudge_x = 0.00055, nudge_y = 0.00030, hjust = 0, fontface = "bold") +
  # scale bar (200 m) and north arrow
  annotate("segment", x = xr[2] - 0.0009 - 200 / (111320 * cos(LAT0 * pi / 180)), xend = xr[2] - 0.0009,
           y = yr[1] + 0.0008, yend = yr[1] + 0.0008, linewidth = 0.8) +
  annotate("text", x = xr[2] - 0.0009 - 100 / (111320 * cos(LAT0 * pi / 180)), y = yr[1] + 0.0011,
           label = "200 m", size = 2.6, family = "Helvetica") +
  annotate("segment", x = xr[2] - 0.0010, xend = xr[2] - 0.0010, y = yr[2] - 0.0025, yend = yr[2] - 0.0009,
           linewidth = 0.5, arrow = arrow(length = unit(0.12, "cm"), type = "closed")) +
  annotate("text", x = xr[2] - 0.0010, y = yr[2] - 0.0028, label = "N", size = 2.8, family = "Helvetica") +
  scale_colour_manual(values = ZONE_COL, labels = paste(sites$site, ZONE_LAB[ZONE_LEVELS])) +
  scale_shape_manual(values = ZONE_SHAPE, labels = paste(sites$site, ZONE_LAB[ZONE_LEVELS])) +
  guides(colour = guide_legend(ncol = 2, byrow = TRUE), shape = guide_legend(ncol = 2, byrow = TRUE)) +
  coord_fixed(ratio = ASP, xlim = xr, ylim = yr, expand = FALSE) +
  theme_pub() +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.4),
        legend.position = "inside", legend.position.inside = c(0.006, 0.994),
        legend.justification.inside = c(0, 1), legend.text = element_text(size = 8),
        legend.key.size = unit(0.32, "cm"), legend.spacing.y = unit(0, "cm"),
        legend.spacing.x = unit(0.15, "cm"), legend.margin = margin(2, 4, 2, 2),
        legend.background = element_rect(fill = alpha("white", 0.88), colour = NA))

# ---- (b) session grid
SLOTS <- c("06-09", "09-12", "12-15", "15-18", "18-21")
DATES <- c("2026-07-02", "2026-07-04", "2026-07-05", "2026-07-06", "2026-07-07")
DLAB <- c("Thu\n2 Jul", "Sat\n4 Jul", "Sun\n5 Jul", "Mon\n6 Jul", "Tue\n7 Jul")
grid <- sess %>%
  mutate(date = as.character(date), slot = factor(slot_planned, levels = SLOTS),
         date_f = factor(date, levels = DATES, labels = DLAB),
         layer = ifelse(wav_recorded == 1, "Level, recording, panorama", "Level and panorama only"),
         lab = sprintf("%d\n%s", n_records, gsub("\\|", "+", cohort)))
f1b <- ggplot(grid, aes(date_f, slot)) +
  geom_tile(aes(fill = layer), colour = "white", linewidth = 0.8, width = 0.96, height = 0.92) +
  geom_text(aes(label = lab), size = 2.5, family = "Helvetica", lineheight = 0.9) +
  scale_fill_manual(values = c("Level, recording, panorama" = "#CFE3F2",
                               "Level and panorama only" = "grey88"),
                    breaks = c("Level, recording, panorama", "Level and panorama only")) +
  scale_x_discrete(drop = FALSE, expand = expansion(add = 0.55)) +
  scale_y_discrete(limits = rev(SLOTS), labels = function(x) gsub("-", "–", x)) +
  labs(x = "Survey day (July 2026)", y = "Time slot (h)") +
  theme_pub() +
  guides(fill = guide_legend(nrow = 2)) +
  theme(legend.position = "top", legend.location = "plot", legend.justification = "left",
        legend.direction = "horizontal", legend.key.size = unit(0.32, "cm"),
        legend.margin = margin(0, 0, 2, 0), axis.line = element_blank(),
        axis.ticks = element_blank())

f1 <- (f1a | f1b) + plot_layout(widths = c(1.1, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 10, face = "bold", family = "Helvetica"))
sv(f1, "fig1_design.png", W2, W2 * 0.50)
}

# =============================================================================
# Fig. 2 — zone profiles: outcome EMMs (a, b), source fingerprints (c), quadrant shares (d)
# =============================================================================
em <- OUTC("analysis_03_zone_profiles", "zone_emmeans.csv") %>%
  mutate(zone = factor(zone, levels = ZONE_LEVELS),
         outcome = factor(outcome, levels = names(OUT_LAB)))
emm_panel <- function(outs, xlab, xlim = NULL) {
  ggplot(filter(em, outcome %in% outs), aes(emmean, zone, colour = outcome)) +
    geom_vline(xintercept = if ("comfort" %in% outs) 0 else 3, colour = "grey60", linewidth = 0.3,
               linetype = "22") +
    geom_errorbar(aes(xmin = lower.CL, xmax = upper.CL), width = 0,
                   position = position_dodge(width = 0.55), linewidth = 0.4) +
    geom_point(position = position_dodge(width = 0.55), size = 1.9) +
    scale_colour_manual(values = OUT_COL[outs], labels = OUT_LAB[outs], breaks = outs) +
    scale_x_continuous(limits = xlim) +
    zone_y() + labs(x = xlab, y = NULL) + theme_pub() +
    theme(legend.position = "top", legend.location = "plot", legend.justification = "left",
          legend.direction = "horizontal", legend.key.size = unit(0.3, "cm"),
          legend.margin = margin(0, 0, 0, 0))
}
f2a <- emm_panel(c("comfort", "richness"), "Comfort or Richness (ICSAQ, −1 to 1)", c(-0.42, 0.42))
f2b <- emm_panel(c("rating_appropriateness", "rating_noisiness"), "Rating (1–5)", c(1.5, 4.1))

SRC_LAB <- c(wind = "wind", water = "water", insects = "insects", birds = "birds",
             animals = "animals", leaves = "leaves", conversation = "conversation",
             children = "children", playing = "playing", footsteps = "footsteps",
             running = "running", teaching = "teaching", broadcast = "broadcast",
             bike_voice = "bike voice", vehicles = "vehicles", vehicle_bumps = "vehicle jolts",
             rail = "rail", sweeping = "sweeping", mower = "mower", ac_units = "AC units",
             loading = "loading", luggage = "luggage", gates = "gates", construction = "construction")
fp <- OUTC("analysis_03_zone_profiles", "source_fingerprint_long.csv") %>%
  mutate(zone = factor(zone, levels = ZONE_LEVELS), source = factor(source, levels = SOURCES))
fp_lim <- ceiling(max(abs(fp$excess)) * 2) / 2
cat("F2c: excess salience range", round(range(fp$excess), 2), "-> limits ±", fp_lim, "\n")
f2c <- ggplot(fp, aes(source, zone, fill = excess)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  annotate("segment", x = c(6.5, 13.5), xend = c(6.5, 13.5), y = 0.4, yend = 8.6,
           linewidth = 0.5, colour = "black") +
  scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00", midpoint = 0,
                       limits = c(-fp_lim, fp_lim), breaks = c(-3, 0, 3),
                       name = "Excess salience\n(points, 0–10 scale)") +
  scale_x_discrete(labels = SRC_LAB, expand = c(0, 0)) +
  zone_y() + labs(x = NULL, y = NULL) + theme_pub(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
        axis.line = element_blank(), axis.ticks = element_blank(),
        legend.position = "right", legend.title = element_text(size = 7.5),
        legend.key.height = unit(0.5, "cm"), legend.key.width = unit(0.28, "cm"),
        legend.text = element_text(size = 7.5))

qr <- OUTC("analysis_16_quadrants", "quadrant_by_recording.csv") %>%
  mutate(zone = factor(zone, levels = ZONE_LEVELS), quadrant = factor(quadrant, levels = QUAD_LEVELS)) %>%
  count(zone, quadrant) %>% group_by(zone) %>% mutate(pct = 100 * n / sum(n)) %>% ungroup()
f2d <- ggplot(qr, aes(pct, zone, fill = quadrant)) +
  geom_col(width = 0.72, colour = NA) +
  scale_fill_manual(values = QUAD_FILL, breaks = QUAD_LEVELS, labels = QUAD_LAB) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.02)), breaks = seq(0, 100, 25)) +
  zone_y() + labs(x = "Share of stops (%)", y = NULL) + theme_pub() +
  guides(fill = guide_legend(nrow = 2)) +
  theme(legend.position = "top", legend.location = "plot", legend.justification = "left",
        legend.direction = "horizontal", legend.key.size = unit(0.3, "cm"),
        legend.margin = margin(0, 0, 0, 0))

f2 <- (f2a | f2b) / ((f2c | f2d) + plot_layout(widths = c(1.45, 1))) + plot_layout(heights = c(1, 1.1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 10, face = "bold", family = "Helvetica"))
sv(f2, "fig2_zones.png", W2, W2 * 0.86)

# =============================================================================
# Fig. 4 — the two axes: comfort x richness at recording level (a), the three
# headline sources' pulls (b), the Richness interaction (c, d)
# =============================================================================
qr_env <- OUTC("analysis_16_quadrants", "quadrant_by_recording.csv") %>%
  mutate(zone = factor(zone, levels = ZONE_LEVELS))
r_env <- cor(qr_env$comfort, qr_env$richness)
cat(sprintf("F3a: %d recordings; recording-level r(comfort, richness) = %.3f\n", nrow(qr_env), r_env))
f3a <- ggplot(qr_env, aes(richness, comfort, colour = zone, shape = zone)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3, linetype = "22") +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3, linetype = "22") +
  geom_point(size = 1.8, alpha = 0.9) +
  annotate("text", x = c(-0.68, 0.68, -0.68, 0.68), y = c(0.64, 0.64, -0.64, -0.64),
           label = unname(QUAD_LAB[QUAD_LEVELS]), colour = "grey45",
           hjust = c(0, 1, 0, 1), size = 2.8, family = "Helvetica", fontface = "italic") +
  scale_colour_manual(values = ZONE_COL, labels = ZONE_LAB) +
  scale_shape_manual(values = ZONE_SHAPE, labels = ZONE_LAB) +
  scale_x_continuous(limits = c(-0.7, 0.7), breaks = seq(-0.6, 0.6, 0.3)) +
  scale_y_continuous(limits = c(-0.7, 0.7), breaks = seq(-0.6, 0.6, 0.3)) +
  coord_fixed() +
  labs(x = "Richness (ICSAQ, −1 to 1)", y = "Comfort (ICSAQ, −1 to 1)") +
  theme_pub() +
  theme(legend.position = "right", legend.text = element_text(size = 8),
        legend.key.size = unit(0.3, "cm"), legend.margin = margin(0, 0, 0, 2))

ss <- OUTC("analysis_05_sources_to_perception", "single_source_screen.csv") %>%
  filter(source %in% c("birds", "vehicles", "conversation"), outcome %in% c("comfort", "richness")) %>%
  mutate(source = factor(source, levels = rev(c("birds", "conversation", "vehicles"))),
         outcome = factor(outcome, levels = c("comfort", "richness")),
         lo = beta_per_point - 1.96 * se, hi = beta_per_point + 1.96 * se)
f3b <- ggplot(ss, aes(beta_per_point, source, colour = outcome)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3, linetype = "22") +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, position = position_dodge(width = 0.5),
                 linewidth = 0.4) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  scale_colour_manual(values = OUT_COL, labels = OUT_LAB) +
  scale_x_continuous(limits = c(-0.025, 0.045), breaks = c(-0.02, 0, 0.02, 0.04)) +
  scale_y_discrete(labels = c(birds = "Birds", conversation = "Conversation", vehicles = "Vehicles")) +
  labs(x = "Change per salience point\n(ICSAQ units)", y = NULL) + theme_pub() +
  theme(legend.position = "top", legend.justification = "right",
        legend.key.size = unit(0.3, "cm"))
slopes <- OUTC("analysis_20_storyline_checks", "natural_slope_by_mechanical.csv")
f3c <- ggplot(slopes, aes(src_mechanical, src_natural.trend)) +
  geom_hline(yintercept=0, colour="grey60", linewidth=0.3, linetype="22") +
  geom_ribbon(aes(ymin=lower.CL,ymax=upper.CL), fill=OUT_COL[["richness"]], alpha=0.20) +
  geom_line(colour=OUT_COL[["richness"]], linewidth=0.65) +
  scale_x_continuous(breaks=0:3) +
  labs(x="Mechanical-source salience (0–10)",
       y="Richness change per point\nof natural-source salience") + theme_pub()
MIX_LAB <- c(total="Total salience (3)", total_diversity="Total + diversity (4)",
             groups="Three source groups (5)", groups_diversity="Groups + diversity (6)",
             groups_interaction="Groups + interaction (6)", perceived_detected="Perceived + detected (8)",
             pca3="Three source PCs (5)")
mix <- OUTC("analysis_20_storyline_checks", "cv_scores.csv") %>% filter(variant=="primary") %>% select(model,r2_env_mean)
mix <- bind_rows(mix, OUTC("analysis_19_source_structure", "cv_r2_by_set.csv") %>%
                  filter(outcome=="richness",scheme=="site",set=="pca3") %>% transmute(model=set,r2_env_mean)) %>%
  mutate(model=factor(model,levels=rev(names(MIX_LAB))))
f3d <- ggplot(mix,aes(r2_env_mean,model)) +
  geom_vline(xintercept=0, colour="grey60",linewidth=0.3,linetype="22") +
  geom_segment(aes(x=0,xend=r2_env_mean,yend=model),colour="grey75",linewidth=0.4) +
  geom_point(colour=OUT_COL[["richness"]],size=2) +
  geom_text(aes(label=sprintf("%.2f",r2_env_mean)),hjust=-0.55,size=2.7,family="Helvetica") +
  scale_x_continuous(limits=c(0,0.25),breaks=c(0,0.1,0.2)) +
  scale_y_discrete(labels=MIX_LAB) +
  labs(x=expression(atop("Richness prediction"~italic(R)^2,"at held-out sites (stop means)")),y=NULL) + theme_pub()
f3 <- ((f3a | f3b) / (f3c | f3d)) + plot_layout(widths = c(1.25, 1), heights=c(1.2,1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 10, face = "bold", family = "Helvetica"))
sv(f3, "fig4_axes.png", W2, W2 * 0.93)

# =============================================================================
# Fig. 5 — the five-feature model: forests (a, b), predictor ladder (c); the
# held-out comparison of feature sets is Fig. 6b
# =============================================================================
hv <- OUTC("analysis_12_robustness", "headline_model_variants.csv") %>%
  filter(variant %in% c("full_666", "all_990_xl2"), term %in% FEAT_LEVELS) %>%
  transmute(term = factor(term, levels = rev(FEAT_LEVELS)),
            outcome = factor(outcome, levels = names(OUT_LAB)),
            est = Estimate, lo = Estimate - 1.96 * `Std. Error`, hi = Estimate + 1.96 * `Std. Error`,
            sample = ifelse(variant == "full_666", "666 ratings (30-s region)", "990 ratings (XL2 level)"))
forest_panel <- function(outs, xlab, xlim) {
  dd <- filter(hv, outcome %in% outs)
  ggplot(dd, aes(est, term, colour = outcome, shape = sample)) +
    geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3, linetype = "22") +
    geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, show.legend = FALSE,
                   position = position_dodge(width = 0.6), linewidth = 0.4) +
    geom_point(position = position_dodge(width = 0.6), size = 1.9, fill = "white", show.legend = TRUE) +
    OUT_SCALE +
    scale_shape_manual(values = c("666 ratings (30-s region)" = 16, "990 ratings (XL2 level)" = 21)) +
    scale_x_continuous(limits = xlim) +
    scale_y_discrete(labels = FEAT_LAB) +
    labs(x = xlab, y = NULL) + theme_pub() +
    guides(colour = "none", shape = guide_legend(order = 2)) +
    theme(axis.text.y = element_markdown(size = 9, colour = "black"))
}
# one colour scale for every panel, so that the extracted legend lists all four outcomes
OUT_SCALE <- scale_colour_manual(values = OUT_COL, labels = OUT_LAB, breaks = names(OUT_COL),
                                 limits = names(OUT_COL), drop = FALSE)
f4a <- forest_panel(c("comfort", "richness"), "Change per SD (ICSAQ units)", c(-0.15, 0.16))
f4b <- forest_panel(c("rating_appropriateness", "rating_noisiness"), "Change per SD (rating units)",
                    c(-0.36, 0.42))

lad <- OUTC("analysis_09_variance_partitioning", "predictor_ladder_variance.csv") %>%
  filter(step != "none") %>%
  mutate(step = factor(step, levels = c("level", "level_tonality", "plus_sources", "plus_zone"),
                       labels = c("Level", "Tonality", "Sources", "Zone")),
         outcome = factor(outcome, levels = names(OUT_LAB)))
f4c <- ggplot(lad, aes(step, env_var_explained, colour = outcome, group = outcome)) +
  geom_line(linewidth = 0.5, show.legend = FALSE) + geom_point(size = 1.9) +
  OUT_SCALE +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 0.35)) +
  labs(x = "Predictors added cumulatively", y = "Reduction in between-stop\nvariance (share)") + theme_pub() +
  guides(colour = guide_legend(order = 1))

# one legend for the whole figure, drawn in the empty fourth cell: the outcome colours from c
# and the sample shapes from a, extracted with cowplot (patchwork's guide_area does not
# combine with the freed panel c)
leg_src <- f4a + guides(colour = guide_legend(order = 1), shape = guide_legend(order = 2)) +
  theme(legend.position = "right", legend.direction = "vertical", legend.box = "vertical",
        legend.box.just = "left", legend.key.size = unit(0.3, "cm"), legend.spacing.y = unit(0.15, "cm"),
        legend.text = element_text(size = 9))
f4_leg <- cowplot::get_plot_component(leg_src, "guide-box-right", return_all = FALSE)
nl <- function(p) p + theme(legend.position = "none")
f4 <- (nl(f4a) | nl(f4b)) / (free(nl(f4c), type = "label") | wrap_elements(full = f4_leg)) +
  plot_annotation(tag_levels = list(c("a", "b", "c", ""))) &
  theme(plot.tag = element_text(size = 10, face = "bold", family = "Helvetica"))
sv(f4, "fig5_model.png", W2, W2 * 0.80)

# =============================================================================
# Fig. 3 — the place judgement: richness -> appropriateness by zone with the rater's own
# and the co-walkers' Richness (a; A10, A21), quadrant EMMs (b)
# =============================================================================
mz <- fit_backbone("xl2_laeq_z + zone * richness + comfort", "rating_appropriateness", d_all)
m0 <- fit_backbone("xl2_laeq_z + zone + richness + comfort", "rating_appropriateness", d_all)
lrt <- anova(m0, mz, refit = TRUE)
tr <- as_tibble(emtrends(mz, "zone", var = "richness")) %>%
  transmute(zone, slope = richness.trend, se = SE, lo = lower.CL, hi = upper.CL)
f5_out <- tr %>% mutate(comfort_beta = lme4::fixef(m0)["comfort"],
                        comfort_se = summary(m0)$coefficients["comfort", "Std. Error"],
                        richness_common = lme4::fixef(m0)["richness"],
                        richness_common_se = summary(m0)$coefficients["richness", "Std. Error"],
                        lrt_chisq = lrt$Chisq[2], lrt_df = lrt$Df[2], lrt_p = lrt$`Pr(>Chisq)`[2],
                        n_records = nrow(d_all))
readr::write_csv(f5_out, file.path(OUTP, "f3_richness_by_zone.csv"))
cat("F3a two-axis model: comfort beta", round(f5_out$comfort_beta[1], 3), "; common richness beta",
    round(f5_out$richness_common[1], 3), "; zone x richness LRT chisq", round(lrt$Chisq[2], 1),
    "df", lrt$Df[2], "p", signif(lrt$`Pr(>Chisq)`[2], 2), "\n")
print(tr %>% mutate(across(where(is.numeric), ~round(.x, 3))))
cwsl <- OUTC("analysis_21_cowalker_salience", "zone_slopes_cowalker.csv") %>% filter(slope == "richness")
own_a21 <- OUTC("analysis_21_cowalker_salience", "zone_slopes_own.csv") %>% filter(slope == "richness")
stopifnot(max(abs(own_a21$trend[match(as.character(tr$zone), own_a21$zone)] - tr$slope)) < 1e-6)
SRC_LAB <- c(own = "Rater's own Richness", cowalker = "Co-walkers' mean Richness")
trb <- bind_rows(tr %>% mutate(src = "own"),
                 cwsl %>% transmute(zone, slope = trend, se = SE, lo = lower.CL, hi = upper.CL, src = "cowalker")) %>%
  mutate(zone = factor(zone, levels = ZONE_LEVELS), src = factor(src, levels = names(SRC_LAB)))
f5a <- ggplot(trb, aes(slope, zone, shape = src)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3, linetype = "22") +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, linewidth = 0.4, colour = "#009E73",
                position = position_dodge(width = 0.6)) +
  geom_point(size = 2.1, colour = "#009E73", fill = "white", position = position_dodge(width = 0.6)) +
  scale_shape_manual(values = c(own = 16, cowalker = 21), labels = SRC_LAB, breaks = names(SRC_LAB)) +
  scale_x_continuous(limits = c(-1, 2), breaks = seq(-1, 2, 0.5)) +
  zone_y() + labs(x = "Appropriateness per unit Richness (rating units)", y = NULL) + theme_pub() +
  theme(legend.position = "top", legend.location = "plot", legend.justification = "left",
        legend.direction = "horizontal", legend.key.size = unit(0.3, "cm"), legend.text = element_text(size = 8),
        legend.margin = margin(0, 0, 0, 0))

qe <- OUTC("analysis_16_quadrants", "appropriateness_emmeans_quadrant_by_zone.csv") %>%
  mutate(zone = factor(zone, levels = ZONE_LEVELS), quadrant = factor(quadrant, levels = QUAD_LEVELS))
f5b <- ggplot(qe, aes(emmean, zone, colour = quadrant)) +
  geom_vline(xintercept = 3, colour = "grey60", linewidth = 0.3, linetype = "22") +
  geom_errorbar(aes(xmin = lower.CL, xmax = upper.CL), width = 0,
                 position = position_dodge(width = 0.7), linewidth = 0.4) +
  geom_point(position = position_dodge(width = 0.7), size = 1.7) +
  scale_colour_manual(values = QUAD_COL, breaks = QUAD_LEVELS, labels = QUAD_LAB) +
  scale_x_continuous(limits = c(1, 5), breaks = 1:5) +
  zone_y() + labs(x = "Appropriateness (1–5)", y = NULL) + theme_pub() +
  theme(legend.position = "top", legend.location = "plot", legend.justification = "left",
        legend.direction = "horizontal", legend.key.size = unit(0.3, "cm"),
        legend.margin = margin(0, 0, 0, 0))
f5 <- (f5a | f5b) + plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 10, face = "bold", family = "Helvetica"))
sv(f5, "fig3_fit.png", W2, W2 * 0.50)

# =============================================================================
# Fig. 6 — what measurements and other listeners recover at held-out sites:
# the predictor hierarchy (a; A21) and further feature sets (b; A11)
# =============================================================================
R2_HELD <- "Prediction *R*<sup>2</sup> at held-out sites<br>(stop means)"
H6 <- c(level = "Level", level_tonality = "Level + tonality", detected = "+ detected sources",
        cowalker = "+ co-walkers' sources", own = "+ rater's own sources")
hh <- OUTC("analysis_21_cowalker_salience", "heldout_hierarchy.csv") %>%
  filter(features %in% names(H6)) %>%
  mutate(features = factor(features, levels = rev(names(H6))), outcome = factor(outcome, levels = names(OUT_LAB)))
a11 <- OUTC("analysis_11_feature_families", "cv_r2_by_family.csv") %>% filter(scheme == "site")
chk <- hh %>% filter(features %in% c("level", "level_tonality", "own")) %>%
  mutate(fam = c(level = "level_only", level_tonality = "level_tonality", own = "acoustic_sourcegroups")[as.character(features)]) %>%
  left_join(a11 %>% select(outcome, fam = family, a11 = r2_env_mean), by = c("outcome", "fam"))
stopifnot(max(abs(chk$r2_stop_mean - chk$a11)) < 1e-9)
f6a <- ggplot(hh, aes(r2_stop_mean, features, colour = outcome)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3, linetype = "22") +
  geom_point(position = position_dodge(width = 0.6), size = 1.9) +
  OUT_SCALE + scale_x_continuous(limits = c(-0.05, 0.62), breaks = seq(0, 0.6, 0.2)) +
  scale_y_discrete(labels = H6) +
  labs(x = R2_HELD, y = NULL) + theme_pub() + theme(axis.title.x = element_markdown(lineheight = 1.15))
FAM6 <- c(acoustic = "Full acoustic set (7)", visual = "Visual classes (6)", source_groups = "Three source groups",
          sources_24 = "24 single sources", all = "All features (12)", acoustic_sourcegroups = "Five-feature model")
cv6 <- a11 %>% filter(family %in% names(FAM6)) %>%
  mutate(family = factor(family, levels = rev(names(FAM6))), outcome = factor(outcome, levels = names(OUT_LAB)),
         clipped = r2_env_mean < -0.3, x = pmax(r2_env_mean, -0.3))
f6b <- ggplot(cv6, aes(x, family, colour = outcome)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3, linetype = "22") +
  geom_point(data = filter(cv6, !clipped), position = position_dodge(width = 0.6), size = 1.9) +
  geom_text(data = filter(cv6, clipped), aes(label = "<"), position = position_dodge(width = 0.6),
            size = 3, family = "Helvetica", show.legend = FALSE) +
  OUT_SCALE + scale_x_continuous(limits = c(-0.3, 0.62), breaks = seq(-0.2, 0.6, 0.2)) +
  scale_y_discrete(labels = FAM6) +
  labs(x = R2_HELD, y = NULL) + theme_pub() + theme(axis.title.x = element_markdown(lineheight = 1.15))
f6 <- (f6a | f6b) + plot_layout(guides = "collect") + plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 10, face = "bold", family = "Helvetica"),
        legend.position = "bottom", legend.direction = "horizontal", legend.key.size = unit(0.3, "cm"),
        legend.margin = margin(0, 0, 0, 0), legend.text = element_text(size = 8))
sv(f6, "fig6_measurement.png", W2, W2 * 0.46)
cat("done\n")
