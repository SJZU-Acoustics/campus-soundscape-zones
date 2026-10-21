# =============================================================================
# build_tables.R — LaTeX fragments of Tables 1-2 and Supplementary Tables S1-S34.
# Each fragment is a bare booktabs tabular; captions and table environments are
# in the paper and its Supplementary Information.
# Sources: the deposited tables (Table 1) and the analysis result tables in output/.
# =============================================================================
source(file.path("code", "style.R"))

fnum <- function(x, d = 3) ifelse(is.na(x), "--", sub("^-", "$-$", sprintf(paste0("%.", d, "f"), x)))
fp <- function(p) ifelse(is.na(p), "--", ifelse(p < 0.001, "$<$0.001", sprintf("%.3f", p)))
fci <- function(lo, hi, d = 2) paste0(fnum(lo, d), ", ", fnum(hi, d))
bse <- function(b, se, d = 3) paste0(fnum(b, d), " (", fnum(se, d), ")")
esc <- function(x) gsub("%", "\\\\%", gsub("&", "\\\\&", x))
frag <- function(rows, header, align, name, extra_head = NULL) {
  lines <- c(paste0("\\begin{tabular}{", align, "}"), "\\toprule", extra_head,
             paste0(header, " \\\\"), "\\midrule", rows, "\\bottomrule", "\\end{tabular}")
  dest <- if (startsWith(name, "si_")) SITAB else TAB
  writeLines(lines, file.path(dest, name)); cat("wrote", name, "(", length(rows), "rows )\n")
}
row <- function(...) paste(c(...), collapse = " & ")
mid <- "\\midrule"
sub_head <- function(txt, ncol) paste0("\\multicolumn{", ncol, "}{l}{\\emph{", txt, "}} \\\\")
OUTN <- c(comfort = "Comfort", richness = "Richness", rating_appropriateness = "Appropriateness",
          rating_noisiness = "Noisiness", rating_auditory_env = "Auditory environment",
          rating_visual_env = "Visual environment")
OUT4 <- names(OUT_LAB)
ZN <- function(z) ZONE_LAB[as.character(z)]
TERMN <- c(level_z = "$L_{\\mathrm{Aeq}}$", bk_laeq_z = "$L_{\\mathrm{Aeq}}$ (30-s)",
           xl2_laeq_z = "$L_{\\mathrm{Aeq}}$ (XL2)", tonality_z = "Tonality",
           src_natural_z = "Natural sound", src_human_z = "Human sound",
           src_mechanical_z = "Mechanical sound", `(Intercept)` = "Intercept")
IND <- c(bk_laeq = "$L_{\\mathrm{Aeq}}$", bk_l5 = "$L_{\\mathrm{A5}}$", bk_l95 = "$L_{\\mathrm{A95}}$",
         bk_dl5_95 = "$L_{\\mathrm{A5}}-L_{\\mathrm{A95}}$", n = "$N$", n5 = "$N_{5}$",
         n5_n95 = "$N_{5}/N_{95}$", s = "$S$", r = "$R$", f = "$F$", tonality = "Tonality",
         xl2_laeq = "$L_{\\mathrm{Aeq}}$ (XL2)", xl2_dl5_95 = "$L_{\\mathrm{AF5}}-L_{\\mathrm{AF95}}$ (XL2)")
SRCN <- c(wind = "Wind", water = "Water", insects = "Insects", birds = "Birds", animals = "Animals",
          leaves = "Leaves", conversation = "Conversation", children = "Children", playing = "Playing",
          footsteps = "Footsteps", running = "Running", teaching = "Teaching", broadcast = "Broadcast",
          bike_voice = "Shared-bike voice prompts", vehicles = "Vehicles", vehicle_bumps = "Vehicle jolts",
          rail = "Rail", sweeping = "Sweeping", mower = "Mower", ac_units = "Air-conditioning units",
          loading = "Loading", luggage = "Luggage trolleys", gates = "Gate barriers",
          construction = "Construction")
GRPN <- c(natural = "Natural", human = "Human", mechanical = "Mechanical")

d_all <- load_p35("all"); sites <- load_p35_sites()
env <- load_p35_env() %>% mutate(site = as.integer(as.character(site))) %>%
  left_join(sites %>% select(site, zone), by = "site")
stopifnot(!any(is.na(env$zone)))

# =============================================================================
# Table 1 (main) — sample and sites
# =============================================================================
env_z <- env %>% mutate(zone = factor(zone, levels = ZONE_LEVELS)) %>% group_by(zone) %>%
  summarise(n_env = n(), n_wav = sum(wav_recorded == 1),
            xl2_m = mean(xl2_laeq), xl2_sd = sd(xl2_laeq),
            bk_m = mean(bk_laeq, na.rm = TRUE), bk_sd = sd(bk_laeq, na.rm = TRUE),
            ton_m = mean(tonality, na.rm = TRUE), ton_sd = sd(tonality, na.rm = TRUE),
            nat_m = mean(naturalness), nat_sd = sd(naturalness), .groups = "drop")
rec_z <- d_all %>% group_by(zone) %>%
  summarise(n_rec = n(), n_eval = n_distinct(evaluator), s_nat = mean(src_natural),
            s_hum = mean(src_human), s_mech = mean(src_mechanical), .groups = "drop")
chk <- OUTC("analysis_01_inventory", "env_descriptives_by_zone.csv")
stopifnot(all(abs(env_z$xl2_m - chk$xl2_laeq[match(env_z$zone, chk$zone)]) < 1e-6),
          all(abs(env_z$bk_m - chk$bk_laeq[match(env_z$zone, chk$zone)]) < 1e-6),
          all(env_z$n_wav == chk$n_wav[match(env_z$zone, chk$zone)]))
chk2 <- OUTC("analysis_01_inventory", "source_groups_by_zone.csv")
stopifnot(all(abs(rec_z$s_nat - chk2$natural[match(rec_z$zone, chk2$zone)]) < 1e-6))
t1 <- env_z %>% left_join(rec_z, by = "zone") %>% left_join(sites %>% select(zone, site), by = "zone")
tot <- tibble(zone = "All", site = "", n_rec = nrow(d_all), n_env = nrow(env),
              n_wav = sum(env$wav_recorded == 1),
              xl2_m = mean(env$xl2_laeq), xl2_sd = sd(env$xl2_laeq),
              bk_m = mean(env$bk_laeq, na.rm = TRUE), bk_sd = sd(env$bk_laeq, na.rm = TRUE),
              ton_m = mean(env$tonality, na.rm = TRUE), ton_sd = sd(env$tonality, na.rm = TRUE),
              nat_m = mean(env$naturalness), nat_sd = sd(env$naturalness),
              s_nat = mean(d_all$src_natural), s_hum = mean(d_all$src_human), s_mech = mean(d_all$src_mechanical))
t1r <- function(r, lab) row(lab,
                             paste0(fnum(r$xl2_m, 1), " (", fnum(r$xl2_sd, 1), ")"),
                             paste0(fnum(r$bk_m, 1), " (", fnum(r$bk_sd, 1), ")"),
                             paste0(fnum(r$ton_m, 2), " (", fnum(r$ton_sd, 2), ")"),
                             fnum(r$nat_m, 2), fnum(r$s_nat, 2), fnum(r$s_hum, 2), fnum(r$s_mech, 2))
rows <- c(sapply(seq_len(nrow(t1)), function(i) {
  r <- t1[i, ]; t1r(r, paste0(r$site, " ", ZN(r$zone)))
}) %>% paste0(" \\\\"), mid, paste0(t1r(tot, "All sites"), " \\\\"))
frag(rows, row("Site",
               "\\shortstack{$L_{\\mathrm{Aeq}}$ XL2\\\\(dB(A))}", "\\shortstack{$L_{\\mathrm{Aeq}}$ 30-s\\\\(dB(A))}",
               "Tonality", "\\shortstack{Natural-\\\\ness}", "Natural", "Human", "Mechanical"),
     "lccccccc", "tab1_sites.tex",
     extra_head = "& & & & & \\multicolumn{3}{c}{Source salience (0--10)} \\\\ \\cmidrule(lr){6-8}")

# =============================================================================
# Table 2 (main) — the five-feature model (666 records), the 990-record
# replication and the natural-mechanical exchange rate
# =============================================================================
hv <- OUTC("analysis_12_robustness", "headline_model_variants.csv")
lad <- OUTC("analysis_09_variance_partitioning", "predictor_ladder_variance.csv") %>% filter(step == "plus_sources")
full <- filter(hv, variant == "full_666"); rep990 <- filter(hv, variant == "all_990_xl2")
stopifnot(all(abs(full$r2m[match(OUT4, full$outcome)] - lad$r2m[match(OUT4, lad$outcome)]) < 1e-6))
cell <- function(df, term, o) { r <- df[df$term == term & df$outcome == o, ]; c(bse(r$Estimate, r$`Std. Error`), fp(r$`Pr(>|t|)`)) }
t2rows <- c(sub_head("Five-feature model, 666 ratings with the 30-s region", 9),
            sapply(FEAT_LEVELS, function(tm) paste0(row(TERMN[tm], unlist(lapply(OUT4, function(o) cell(full, tm, o)))), " \\\\")),
            paste0(row("$R^2_{\\mathrm{m}}$ / $R^2_{\\mathrm{c}}$",
                       unlist(lapply(OUT4, function(o) c(paste0(fnum(lad$r2m[lad$outcome == o], 2), " / ",
                                                                 fnum(lad$r2c[lad$outcome == o], 2)), "")))), " \\\\"),
            mid, sub_head("Full-sample analysis using XL2 level (990 questionnaires; no tonality)", 9),
            sapply(c("level_z", "src_natural_z", "src_human_z", "src_mechanical_z"), function(tm)
              paste0(row(TERMN[tm], unlist(lapply(OUT4, function(o) cell(rep990, tm, o)))), " \\\\")))
frag(t2rows, row("Term", "$\\beta$ (SE)", "$p$", "$\\beta$ (SE)", "$p$", "$\\beta$ (SE)", "$p$", "$\\beta$ (SE)", "$p$"),
     "lcccccccc", "tab2_model.tex",
     extra_head = "& \\multicolumn{2}{c}{Comfort} & \\multicolumn{2}{c}{Richness} & \\multicolumn{2}{c}{Appropriateness} & \\multicolumn{2}{c}{Noisiness} \\\\ \\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}\\cmidrule(lr){8-9}")

# exchange-rate block (A13; per salience point, LAeq-adjusted additive model)
ai <- OUTC("analysis_13_compensation", "additive_and_interaction.csv")
bt <- OUTC("analysis_13_compensation", "offset_ratio_bootstrap.csv")
xr <- c("comfort", "rating_appropriateness")   # the two bootstrapped ratios (A13)
xrows <- sapply(xr, function(o) {
  a <- ai[ai$outcome == o, ]; b <- bt[bt$outcome == o, ]
  paste0(row(OUTN[o], bse(a$b_nat, NA, 3) %>% sub(" \\(--\\)", "", .), fp(a$p_nat),
             fnum(a$b_mech, 3), fp(a$p_mech), fnum(a$offset_ratio, 2),
             if (nrow(b)) fci(b$ci_lo, b$ci_hi, 2) else "--", fnum(a$b_int, 3), fp(a$p_int)), " \\\\")
})
# the same model with the co-walkers' salience (A21: mean of the other evaluators at the stop)
cw <- OUTC("analysis_21_cowalker_salience", "exchange_ratio_cowalker.csv")
cwrows <- sapply(xr, function(o) { a <- cw[cw$outcome == o, ]
  paste0(row(OUTN[o], fnum(a$b_nat, 3), fp(a$p_nat), fnum(a$b_mech, 3), fp(a$p_mech), fnum(a$ratio, 2),
             fci(a$ci_lo, a$ci_hi, 2), fnum(a$b_int, 3), fp(a$p_int)), " \\\\") })
xrows <- c(sub_head("Salience reported by the rater", 9), xrows, sub_head("Salience reported by the co-walkers at the stop", 9), cwrows)
frag(xrows, row("Outcome", "Natural $\\beta$", "$p$", "Mechanical $\\beta$", "$p$", "Ratio", "95\\% CI", "Interaction", "$p$"),
     "lcccccccc", "tab2b_exchange.tex")

# =============================================================================
# Zone level targets (A14b zone-specific slopes, 30-s layer); SI since the 16 October 2026 reframe
# =============================================================================
zs <- OUTC("analysis_14_neutral_levels", "neutral_levels_zone_slopes.csv") %>%
  mutate(zone = factor(zone, levels = ZONE_LEVELS)) %>% arrange(outcome, zone)
t3cell <- function(r) {
  sig <- r$slope_hi < 0
  slope <- paste0(fnum(10 * r$slope_per_dB, 2), " (", fci(10 * r$slope_lo, 10 * r$slope_hi, 2), ")")
  if (!sig) return(c(slope, "--$^{a}$", "--"))
  inr <- r$within_observed_range
  neut <- paste0(fnum(r$neutral_level, 1), " (", fci(r$ci_lo, r$ci_hi, 1), ")", if (!inr) "$^{b}$" else "")
  head <- fnum(r$headroom_dB, 1)
  c(slope, neut, head)
}
t3block <- function(o, lab) c(sub_head(lab, 5), sapply(ZONE_LEVELS, function(z) {
  r <- zs[zs$outcome == o & zs$zone == z, ]
  paste0(row(ZN(z), paste0(fnum(r$median_level, 1), " (", fnum(r$lo, 1), "--", fnum(r$hi, 1), ")"), t3cell(r)), " \\\\")
}))
t3rows <- c(t3block("comfort", "Comfort (midpoint $=$ 0)"), mid, t3block("rating_appropriateness", "Appropriateness (midpoint $=$ 3)"))
frag(t3rows, row("Site", "\\shortstack{Observed level\\\\(dB(A))}", "\\shortstack{Slope per 10 dB\\\\(95\\% CI)}", "\\shortstack{Midpoint-crossing level\\\\(95\\% CI)}", "\\shortstack{Difference from median\\\\(dB)}"),
     "lcccc", "si_zone_targets.tex")

# =============================================================================
# SI fragments (named by content; numbered in the SI tex in citation order)
# =============================================================================
# --- A01: adjective PCA loadings; outcome correlations; source prevalence
pca <- OUTC("analysis_01_inventory", "adjective_pca_loadings.csv")
frag(paste0(apply(pca, 1, function(r) row(tools::toTitleCase(r[["adjective"]]), fnum(as.numeric(r[["PC1"]]), 2),
                                            fnum(as.numeric(r[["PC2"]]), 2), fnum(as.numeric(r[["PC3"]]), 2))), " \\\\"),
     row("Adjective", "PC1", "PC2", "PC3"), "lccc", "si_adjective_pca.tex")
oc <- OUTC("analysis_01_inventory", "outcome_correlations.csv") %>% filter(var %in% OUT4) %>% select(var, all_of(OUT4))
frag(paste0(sapply(seq_len(nrow(oc)), function(i) row(OUTN[oc$var[i]], sapply(oc$var, function(v) fnum(oc[[v]][i], 2)))), " \\\\"),
     row("", OUTN[oc$var]), paste0("l", strrep("c", nrow(oc))), "si_outcome_correlations.tex")
sp <- OUTC("analysis_01_inventory", "source_prevalence.csv") %>% mutate(source = factor(source, levels = SOURCES)) %>% arrange(source)
frag(paste0(apply(sp, 1, function(r) row(SRCN[r[["source"]]], fnum(as.numeric(r[["mean"]]), 2), fnum(100 * as.numeric(r[["prevalence"]]), 1),
                                           fnum(100 * as.numeric(r[["share_ge5"]]), 1))), " \\\\"),
     row("Source", "Mean salience (0--10)", "Heard (\\%)", "Salience $\\geq$ 5 (\\%)"), "lccc", "si_source_prevalence.tex")

# --- A02: predictive comparison of the two windows; indicator correlations
wp <- OUTC("analysis_02_acoustic_windows", "window_predictive_comparison.csv")
frag(paste0(apply(wp, 1, function(r) row(OUTN[r[["outcome"]]], bse(as.numeric(r[["beta_bk"]]), as.numeric(r[["se_bk"]])), fp(as.numeric(r[["p_bk"]])), fnum(as.numeric(r[["r2m_bk"]]), 3),
                                          bse(as.numeric(r[["beta_xl2"]]), as.numeric(r[["se_xl2"]])), fp(as.numeric(r[["p_xl2"]])), fnum(as.numeric(r[["r2m_xl2"]]), 3),
                                          fnum(as.numeric(r[["both_beta_bk"]]), 3), fp(as.numeric(r[["both_p_bk"]])), fnum(as.numeric(r[["both_beta_xl2"]]), 3), fp(as.numeric(r[["both_p_xl2"]])))), " \\\\"),
     row("Outcome", "$\\beta$ (SE)", "$p$", "$R^2_{\\mathrm{m}}$", "$\\beta$ (SE)", "$p$", "$R^2_{\\mathrm{m}}$", "$\\beta$ 30-s", "$p$", "$\\beta$ XL2", "$p$"),
     "lcccccccccc", "si_window_predictive.tex",
     extra_head = "& \\multicolumn{3}{c}{30-s region} & \\multicolumn{3}{c}{XL2 full measurement} & \\multicolumn{4}{c}{Both in one model} \\\\ \\cmidrule(lr){2-4}\\cmidrule(lr){5-7}\\cmidrule(lr){8-11}")
ic <- OUTC("analysis_02_acoustic_windows", "indicator_correlations.csv")
frag(paste0(sapply(seq_len(nrow(ic)), function(i) row(IND[ic$indicator[i]], sapply(ic$indicator, function(v) fnum(ic[[v]][i], 2)))), " \\\\"),
     row("", IND[ic$indicator]), paste0("l", strrep("c", nrow(ic))), "si_indicator_correlations.tex")

# --- A04: single-indicator screen; beyond-level screen
si4 <- OUTC("analysis_04_acoustic_to_perception", "single_indicator_screen.csv") %>% mutate(outcome = factor(outcome, levels = OUT4)) %>% arrange(outcome)
frag(paste0(sapply(seq_len(nrow(si4)), function(i) { r <- si4[i, ]; row(OUTN[as.character(r$outcome)], IND[r$indicator], bse(r$beta, r$se), fp(r$p), fnum(r$r2m, 3), fp(r$q)) }), " \\\\"),
     row("Outcome", "Indicator", "$\\beta$ per SD (SE)", "$p$", "$R^2_{\\mathrm{m}}$", "$q$"), "llcccc", "si_single_indicator.tex")
bl <- OUTC("analysis_04_acoustic_to_perception", "beyond_level_screen.csv") %>% mutate(outcome = factor(outcome, levels = OUT4)) %>% arrange(outcome)
frag(paste0(sapply(seq_len(nrow(bl)), function(i) { r <- bl[i, ]; row(OUTN[as.character(r$outcome)], IND[r$indicator], fnum(r$beta_ind, 3), fp(r$p_ind), fnum(r$beta_laeq, 3), fnum(r$delta_r2m, 3), fp(r$q)) }), " \\\\"),
     row("Outcome", "Indicator", "$\\beta$ indicator", "$p$", "$\\beta$ $L_{\\mathrm{Aeq}}$", "$\\Delta R^2_{\\mathrm{m}}$", "$q$"), "llccccc", "si_beyond_level.tex")
# joint acoustic model (A04): level + tonality + sharpness + level range — the lower block of the same SI table
jm <- OUTC("analysis_04_acoustic_to_perception", "multi_indicator_models.csv") %>% filter(term != "(Intercept)")
JTN <- c(bk_laeq_z = "$L_{\\mathrm{Aeq}}$", tonality_z = "Tonality", s_z = "$S$", bk_dl5_95_z = "$L_{\\mathrm{A5}}-L_{\\mathrm{A95}}$")
jrows <- c(sapply(names(JTN), function(tm) paste0(row(JTN[tm], unlist(lapply(OUT4, function(o) { r <- jm[jm$term == tm & jm$outcome == o, ]; c(bse(r$Estimate, r$`Std. Error`), fp(r$`Pr(>|t|)`)) }))), " \\\\")),
           paste0(row("$R^2_{\\mathrm{m}}$ / $R^2_{\\mathrm{c}}$", unlist(lapply(OUT4, function(o) { r <- jm[jm$term == "bk_laeq_z" & jm$outcome == o, ]; c(paste0(fnum(r$r2m, 2), " / ", fnum(r$r2c, 2)), "") }))), " \\\\"))
frag(jrows, row("Term", rep(c("$\\beta$ (SE)", "$p$"), 4)), "lcccccccc", "si_joint_acoustic.tex",
     extra_head = "& \\multicolumn{2}{c}{Comfort} & \\multicolumn{2}{c}{Richness} & \\multicolumn{2}{c}{Appropriateness} & \\multicolumn{2}{c}{Noisiness} \\\\ \\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}\\cmidrule(lr){8-9}")

# --- A05: single-source screen (24 x 4, beta per point with q)
ss <- OUTC("analysis_05_sources_to_perception", "single_source_screen.csv")
ssw <- ss %>% mutate(cellv = paste0(fnum(beta_per_point, 3), " (", fp(q), ")")) %>% select(source, outcome, cellv) %>%
  tidyr::pivot_wider(names_from = outcome, values_from = cellv) %>% mutate(source = factor(source, levels = SOURCES)) %>% arrange(source)
grp_of <- setNames(rep(names(SOURCE_GROUPS), lengths(SOURCE_GROUPS)), unlist(SOURCE_GROUPS))
ssrows <- unlist(lapply(names(SOURCE_GROUPS), function(g) c(sub_head(paste(GRPN[g], "sources"), 5),
  paste0(sapply(SOURCE_GROUPS[[g]], function(s) { r <- ssw[ssw$source == s, ]; row(SRCN[s], r$comfort, r$richness, r$rating_appropriateness, r$rating_noisiness) }), " \\\\"))))
frag(ssrows, row("Source", "Comfort", "Richness", "Appropriateness", "Noisiness"), "lcccc", "si_single_source.tex")

# --- A06: visual screen (+ without 5 July), within-zone share, AV interaction, congruence
vs <- OUTC("analysis_06_visual_to_perception", "visual_screen.csv"); vs2 <- OUTC("analysis_06_visual_to_perception", "visual_screen_without_w04.csv")
VISN <- c(naturalness = "Naturalness", artificiality = "Artificiality", vdi = "VDI", vri = "VRI",
          vis_people = "People share", vis_vehicles = "Vehicle share", vis_sky = "Sky share", vis_vegetation = "Vegetation share")
vs <- vs %>% left_join(vs2 %>% select(outcome, visual, beta_bw2 = beta_bw, q_bw2 = q_bw), by = c("outcome", "visual")) %>%
  mutate(outcome = factor(outcome, levels = OUT4)) %>% arrange(outcome)
frag(paste0(sapply(seq_len(nrow(vs)), function(i) { r <- vs[i, ]; row(OUTN[as.character(r$outcome)], VISN[r$visual], fnum(r$beta_bw, 3), fp(r$q_bw), fnum(r$beta_within, 3), fp(r$q_within), fnum(r$beta_bw2, 3), fp(r$q_bw2)) }), " \\\\"),
     row("Outcome", "Visual variable", "$\\beta$", "$q$", "$\\beta$", "$q$", "$\\beta$", "$q$"), "llcccccc", "si_visual_screen.tex",
     extra_head = "& & \\multicolumn{2}{c}{Between zones} & \\multicolumn{2}{c}{Within zones} & \\multicolumn{2}{c}{Between, without 5 July} \\\\ \\cmidrule(lr){3-4}\\cmidrule(lr){5-6}\\cmidrule(lr){7-8}")
avi <- OUTC("analysis_06_visual_to_perception", "audio_visual_interaction.csv"); cg <- OUTC("analysis_06_visual_to_perception", "congruence.csv")
frag(paste0(sapply(seq_len(nrow(avi)), function(i) { r <- avi[i, ]; row(OUTN[r$outcome], VISN[r$visual], fnum(r$beta_int, 3), fp(r$p_int), fp(r$q)) }), " \\\\"),
     row("Outcome", "Visual variable", "$\\beta$ (visual $\\times$ $L_{\\mathrm{Aeq}}$)", "$p$", "$q$"), "llccc", "si_av_interaction.tex")
HS <- c(src_mechanical = "Mechanical sound", src_human = "Human sound", src_natural = "Natural sound")
SN <- c(vis_vehicles = "Vehicle share", vis_people = "People share", naturalness = "Naturalness", vis_vegetation = "Vegetation share")
frag(paste0(sapply(seq_len(nrow(cg)), function(i) { r <- cg[i, ]; row(OUTN[r$outcome], HS[r$heard], SN[r$seen], fnum(r$beta_heard, 3), fnum(r$beta_seen, 3), fnum(r$beta_int, 3), fp(r$p_int), fp(r$q)) }), " \\\\"),
     row("Outcome", "Heard", "Seen", "$\\beta$ heard", "$\\beta$ seen", "$\\beta$ heard $\\times$ seen", "$p$", "$q$"), "lllccccc", "si_congruence.tex")

# --- A07: day type and slot (ratings; environment)
rd <- OUTC("analysis_07_daytype_slot", "ratings_daytype_slot.csv")
frag(paste0(sapply(seq_len(nrow(rd)), function(i) { r <- rd[i, ]; row(OUTN[r$outcome], fnum(r$weekday, 2), fnum(r$weekend, 2), fp(r$p_daytype), fp(r$q_daytype), fp(r$p_slot), fp(r$q_slot), fp(r$p_daytype_sessionRE)) }), " \\\\"),
     row("Outcome", "Weekday", "Weekend", "$p$", "$q$", "$p$", "$q$", "$p$"), "lccccccc", "si_daytype_ratings.tex",
     extra_head = "& \\multicolumn{2}{c}{EMM} & \\multicolumn{2}{c}{Day type} & \\multicolumn{2}{c}{Slot} & Session RE \\\\ \\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}")
ed <- OUTC("analysis_07_daytype_slot", "env_level_daytype_slot.csv")
EDN <- c(xl2_laeq = "$L_{\\mathrm{Aeq}}$, XL2 (dB(A))", bk_laeq = "$L_{\\mathrm{Aeq}}$, 30-s (dB(A))", n_sources = "Sources heard (count)",
         src_natural = "Natural salience", src_human = "Human salience", src_mechanical = "Mechanical salience", tonality = "Tonality", naturalness = "Naturalness")
frag(paste0(sapply(seq_len(nrow(ed)), function(i) { r <- ed[i, ]; row(ifelse(is.na(EDN[r$variable]), r$variable, EDN[r$variable]), r$n, fnum(r$weekday, 2), fnum(r$weekend, 2), fp(r$p_daytype), fp(r$q_daytype), esc(r$slot_emm), fp(r$p_slot), fp(r$q_slot)) }), " \\\\"),
     row("Variable", "Stops", "Weekday", "Weekend", "$p$", "$q$", "Slot EMMs", "$p$", "$q$"), "lcccccp{4.2cm}cc", "si_daytype_environment.tex")

# --- A08: person-level main effects and moderation; random slopes
PVN <- c(ns_z = "Noise sensitivity (per SD)", who5_z = "WHO-5 (per SD)", rater_type = "Trained rater", sex = "Sex (female)", education = "Education (masters)")
pm <- OUTC("analysis_08_person_level", "person_main_effects.csv") %>% mutate(outcome = factor(outcome, levels = OUT4)) %>% arrange(outcome)
frag(paste0(sapply(seq_len(nrow(pm)), function(i) { r <- pm[i, ]; row(OUTN[as.character(r$outcome)], PVN[r$pv], fnum(r$beta, 3), fp(r$p), fnum(r$beta_990, 3), fp(r$p_990), fp(r$q)) }), " \\\\"),
     row("Outcome", "Person variable", "$\\beta$ (666)", "$p$", "$\\beta$ (990)", "$p$", "$q$"), "llccccc", "si_person_main.tex")
pmod <- OUTC("analysis_08_person_level", "person_moderation_of_level_slope.csv") %>% mutate(outcome = factor(outcome, levels = OUT4)) %>% arrange(outcome)
frag(paste0(sapply(seq_len(nrow(pmod)), function(i) { r <- pmod[i, ]; row(OUTN[as.character(r$outcome)], PVN[r$pv], fnum(r$beta_int, 3), fp(r$p_int), fnum(r$beta_int_990, 3), fp(r$p_int_990), fp(r$q)) }), " \\\\"),
     row("Outcome", "Person variable", "$\\beta$ $\\times L_{\\mathrm{Aeq}}$ (666)", "$p$", "$\\beta$ $\\times L_{\\mathrm{Aeq}}$ (990)", "$p$", "$q$"), "llccccc", "si_person_moderation.tex")
# one source for every column: the ML random-slope fit and its within-evaluator permutation null, all four outcomes (A09b)
rs <- OUTC("analysis_09_variance_partitioning", "random_slope_permutation_all.csv") %>% mutate(outcome = factor(outcome, levels = OUT4)) %>% arrange(outcome)
frag(paste0(sapply(seq_len(nrow(rs)), function(i) { r <- rs[i, ]; row(OUTN[as.character(r$outcome)], fnum(r$mean_slope, 3), fnum(r$sd_slope, 3),
  paste0(fnum(r$blup_slope_min, 3), " to ", fnum(r$blup_slope_max, 3)), paste0(r$n_positive_slopes, " of ", r$n_evaluators), fnum(r$lrt_observed, 1), fnum(r$null_q95, 1), fnum(r$null_q99, 1), fp(r$p_permutation)) }), " \\\\"),
     row("Outcome", "Mean", "SD", "Range", "Positive", "LRT", "Null 95th", "Null 99th", "$p_{\\mathrm{perm}}$"), "lcccccccc", "si_random_slopes.tex")

# --- A09: variance partition (990) + recording-mean reliability
vp <- OUTC("analysis_09_variance_partitioning", "variance_partition_990.csv") %>% select(outcome, grp, share) %>% tidyr::pivot_wider(names_from = grp, values_from = share)
rl <- OUTC("analysis_09_variance_partitioning", "recording_mean_reliability.csv")
vp <- vp %>% left_join(rl, by = "outcome") %>% mutate(outcome = factor(outcome, levels = OUT4)) %>% arrange(outcome)
frag(paste0(sapply(seq_len(nrow(vp)), function(i) { r <- vp[i, ]; row(OUTN[as.character(r$outcome)], fnum(100 * r$site, 1), fnum(100 * r$env_id, 1), fnum(100 * r$evaluator, 1), fnum(100 * r$session, 1), fnum(100 * r$Residual, 1),
  fnum(r$icc1_place, 2), fnum(r$icc1_evaluator, 2), fnum(r$k_raters_per_recording, 2), fnum(r$icc_k_place, 2)) }), " \\\\"),
     row("Outcome", "Site", "Stop", "Evaluator", "Session", "Residual", "ICC(1) place", "ICC(1) evaluator", "$k$", "ICC($k$) place"), "lccccccccc", "si_variance_partition.tex",
     extra_head = "& \\multicolumn{5}{c}{Variance share (\\%)} & \\multicolumn{4}{c}{Reliability of the stop mean} \\\\ \\cmidrule(lr){2-6}\\cmidrule(lr){7-10}")
ladf <- OUTC("analysis_09_variance_partitioning", "predictor_ladder_variance.csv") %>% mutate(outcome = factor(outcome, levels = OUT4)) %>% arrange(outcome)
STEPN <- c(none = "Null", level = "$L_{\\mathrm{Aeq}}$", level_tonality = "+ tonality", plus_sources = "+ three source groups", plus_zone = "+ zone")
frag(paste0(sapply(seq_len(nrow(ladf)), function(i) { r <- ladf[i, ]; row(OUTN[as.character(r$outcome)], STEPN[r$step], fnum(r$var_evaluator, 4), fnum(r$var_env, 4), fnum(r$var_resid, 4), fnum(r$r2m, 3), fnum(r$r2c, 3), fnum(100 * r$env_var_explained, 1)) }), " \\\\"),
     row("Outcome", "Step", "Var(evaluator)", "Var(stop)", "Var(residual)", "$R^2_{\\mathrm{m}}$", "$R^2_{\\mathrm{c}}$", "Stop variance explained (\\%)"), "llcccccc", "si_predictor_ladder.tex")

# --- A10: two-axis appropriateness slopes by site
f5 <- readr::read_csv(file.path(OUTP, "f3_richness_by_zone.csv"), show_col_types = FALSE)
f5c <- OUTC("analysis_10_appropriateness", "two_axis_zone_by_comfort_slopes.csv")   # A10b: the zone × Comfort companion
stopifnot(setequal(f5c$zone, f5$zone))
frag(paste0(sapply(seq_len(nrow(f5)), function(i) { r <- f5[i, ]; k <- f5c[f5c$zone == r$zone, ]; row(paste0(sites$site[sites$zone == r$zone], " ", ZN(r$zone)), bse(r$slope, r$se), fci(r$lo, r$hi, 2), bse(k$slope, k$se), fci(k$lo, k$hi, 2)) }), " \\\\"),
     row("Site", "Richness slope (SE)", "95\\% CI", "Comfort slope (SE)", "95\\% CI"), "lcccc", "si_richness_slopes_two_axis.tex",
     extra_head = "& \\multicolumn{2}{c}{Zone $\\times$ Richness model (Fig.~3a)} & \\multicolumn{2}{c}{Zone $\\times$ Comfort model} \\\\ \\cmidrule(lr){2-3}\\cmidrule(lr){4-5}")

# --- A11: cross-validated R2 by family (recording means; site- and evaluator-held-out)
cv <- OUTC("analysis_11_feature_families", "cv_r2_by_family.csv")
FAMN <- c(level_only = "$L_{\\mathrm{Aeq}}$ only", level_tonality = "$L_{\\mathrm{Aeq}}$ + tonality", acoustic = "Full acoustic set", source_groups = "Three source groups",
          sources_24 = "24 single sources", n_sources = "Number of sources", visual = "Visual classes", acoustic_visual = "Acoustic + visual",
          acoustic_sourcegroups = "Five-feature model", all = "All twelve features")
cvw <- cv %>% mutate(cellv = fnum(r2_env_mean, 2)) %>% select(outcome, family, scheme, k, cellv) %>% tidyr::pivot_wider(names_from = c(outcome, scheme), values_from = cellv)
cvw <- cvw[match(names(FAMN), cvw$family), ]
frag(paste0(sapply(seq_len(nrow(cvw)), function(i) { r <- cvw[i, ]; row(FAMN[r$family], r$k, unlist(lapply(OUT4, function(o) c(r[[paste0(o, "_site")]], r[[paste0(o, "_evaluator")]])))) }), " \\\\"),
     row("Feature set", "$k$", rep(c("Site", "Evaluator"), 4)), "lccccccccc", "si_cv_by_family.tex",
     extra_head = "& & \\multicolumn{2}{c}{Comfort} & \\multicolumn{2}{c}{Richness} & \\multicolumn{2}{c}{Appropriateness} & \\multicolumn{2}{c}{Noisiness} \\\\ \\cmidrule(lr){3-4}\\cmidrule(lr){5-6}\\cmidrule(lr){7-8}\\cmidrule(lr){9-10}")

# --- A12: robustness variants (five-feature model), zone omnibus variants, CLMMs
VARN <- c(full_666 = "Primary (666)", no_flagged = "Without flagged rows", no_trained = "Without trained raters", session_RE = "+ session random effect",
          xl2_on_666 = "XL2 level (666)", all_990_xl2 = "All 990, XL2 level", no_flagged_990 = "990 without flagged rows")
hvw <- hv %>% filter(term %in% FEAT_LEVELS) %>% mutate(cellv = paste0(fnum(Estimate, 3), " (", fp(`Pr(>|t|)`), ")")) %>% select(outcome, term, variant, cellv) %>%
  tidyr::pivot_wider(names_from = variant, values_from = cellv)
rrows <- unlist(lapply(OUT4, function(o) c(sub_head(OUTN[o], 8), paste0(sapply(FEAT_LEVELS, function(tm) { r <- hvw[hvw$outcome == o & hvw$term == tm, ]
  row(TERMN[tm], sapply(names(VARN), function(v) ifelse(is.na(r[[v]]), "--", r[[v]]))) }), " \\\\"))))
frag(rrows, row("Term", paste0("V", seq_along(VARN))), "lccccccc", "si_robustness_variants.tex")
cl <- bind_rows(OUTC("analysis_12_robustness", "clmm_rating_appropriateness.csv"), OUTC("analysis_12_robustness", "clmm_rating_noisiness.csv"))
frag(paste0(sapply(seq_len(nrow(cl)), function(i) { r <- cl[i, ]; row(OUTN[r$outcome], TERMN[r$term], bse(r$Estimate, r$`Std. Error`), fnum(r$`z value`, 2), fp(r$`Pr(>|z|)`)) }), " \\\\"),
     row("Outcome", "Term", "Log-odds (SE)", "$z$", "$p$"), "llccc", "si_clmm.tex")

# --- A13: prediction grid
pg <- OUTC("analysis_13_compensation", "prediction_grid.csv")
frag(paste0(sapply(seq_len(nrow(pg)), function(i) { r <- pg[i, ]; row(r$src_natural, r$src_mechanical, fnum(r$comfort, 3), fnum(r$appropriateness, 2)) }), " \\\\"),
     row("Natural salience", "Mechanical salience", "Predicted comfort", "Predicted appropriateness"), "cccc", "si_prediction_grid.tex")

# --- A15: richness indices
ri <- OUTC("analysis_15_richness_diversity", "single_index_screen.csv")
RIN <- c(n_sources = "Number of sources heard", n_salient = "Number of salient sources ($\\geq$ 5)", total_sal = "Total salience", max_sal = "Maximum salience",
         shannon = "Shannon diversity", evenness = "Evenness", dominance = "Dominance", simpson = "Simpson diversity")
frag(paste0(sapply(seq_len(nrow(ri)), function(i) { r <- ri[i, ]; row(OUTN[r$outcome], ifelse(is.na(RIN[r$index]), r$index, RIN[r$index]), fnum(r$beta, 3), fp(r$p), fnum(r$r2m, 3), fnum(r$aic, 1), fp(r$q)) }), " \\\\"),
     row("Outcome", "Index", "$\\beta$ per SD", "$p$", "$R^2_{\\mathrm{m}}$", "AIC", "$q$"), "llccccc", "si_richness_indices.tex")
rc <- OUTC("analysis_20_storyline_checks", "models.csv") %>% filter(variant=="primary")
rci <- OUTC("analysis_20_storyline_checks", "increments.csv")
rcc <- OUTC("analysis_20_storyline_checks", "cv_scores.csv")
RCN <- c(total="Total salience", total_diversity="Total + diversity", groups="Three source groups",
         groups_diversity="Groups + diversity", groups_interaction="Groups + interaction", perceived_detected="Perceived + detected groups")
rc <- rc %>% left_join(filter(rci,variant=="primary") %>% select(model,chisq,df,p,q),by="model") %>%
  left_join(filter(rcc,variant=="primary") %>% select(model,r2_env_mean),by="model")
frag(paste0(sapply(seq_len(nrow(rc)),function(i) { r<-rc[i,]; row(RCN[r$model],r$k,fnum(r$r2m,3),fnum(r$r2_env_mean,3),
       ifelse(is.na(r$chisq),"--",paste0(fnum(r$chisq,1)," (",r$df,")")),fp(r$p),fp(r$q)) })," \\\\"),
     row("Feature set","$k$","$R^2_{\\mathrm{m}}$","CV $R^2$","LRT (d.f.)","$p$","$q$"),"lcccccc","si_richness_combinations.tex")
rv <- c(primary="Primary", no_flagged="Without quality flags", no_trained="Without trained raters", without_E13_E28="Without two endpoint raters",
        all_xl2="All 990, XL2 level",zone_adjusted="Zone adjusted",salience_curvature="Salience curvature adjusted")
rsens <- OUTC("analysis_20_storyline_checks","coefficients.csv") %>% filter(model=="groups_interaction",term=="src_natural:src_mechanical",variant %in% names(rv)) %>%
  left_join(filter(rci,model=="groups_interaction") %>% select(variant,q),by="variant")
cvwide <- rcc %>% filter(model %in% c("groups","groups_interaction")) %>% select(variant,model,r2_env_mean) %>% pivot_wider(names_from=model,values_from=r2_env_mean)
rsens <- rsens %>% left_join(cvwide,by="variant")
frag(paste0(sapply(names(rv),function(v) {r<-filter(rsens,variant==v); row(rv[v],r$n,paste0(fnum(r$Estimate,3)," (",fci(r$lo,r$hi,3),")"),fp(r$q),fnum(r$groups,3),fnum(r$groups_interaction,3))})," \\\\"),
     row("Variant","Ratings","Interaction (95\\% CI)","$q$","CV additive","CV interaction"),"llcccc","si_richness_sensitivity.tex")

# --- A16: quadrant contrasts within each site
cm <- OUTC("analysis_16_quadrants", "chaotic_minus_monotonous_by_zone.csv"); vc <- OUTC("analysis_16_quadrants", "vibrant_minus_calm_by_zone.csv")
frag(paste0(sapply(ZONE_LEVELS, function(z) { a <- cm[cm$zone == z, ]; b <- vc[vc$zone == z, ]
  row(ZN(z), bse(a$estimate, a$SE, 2), fp(a$p.value), fp(a$p_holm), bse(b$estimate, b$SE, 2), fp(b$p.value), fp(b$p_holm)) }), " \\\\"),
     row("Site", "Estimate (SE)", "$p$", "$p$ (Holm)", "Estimate (SE)", "$p$", "$p$ (Holm)"), "lcccccc", "si_quadrant_contrasts_by_zone.tex",
     extra_head = "& \\multicolumn{3}{c}{Low Comfort: high $-$ low Richness} & \\multicolumn{3}{c}{High Comfort: high $-$ low Richness} \\\\ \\cmidrule(lr){2-4}\\cmidrule(lr){5-7}")

# --- A17: ecoacoustic index screen (primary + full window), increments
es <- OUTC("analysis_17_ecoacoustic_indices", "single_index_screen_primary_with_fullwindow.csv") %>% mutate(outcome = factor(outcome, levels = OUT4)) %>% arrange(outcome, index)
ECN <- c(ACI = "ACI", ADI = "ADI", ADI_10k = "ADI (10 kHz)", AEI = "AEI", AEI_10k = "AEI (10 kHz)", BI = "BI", BI_2_8k = "BI (2--8 kHz)", H_total = "$H$ (total entropy)",
         Hf = "$H_f$", Ht = "$H_t$", NDSI = "NDSI", NDSI_kasten = "NDSI (Kasten bands)")
frag(paste0(sapply(seq_len(nrow(es)), function(i) { r <- es[i, ]; row(OUTN[as.character(r$outcome)], ECN[r$index], bse(r$beta, r$se), fp(r$p), fp(r$q), fnum(r$r2m, 3), fnum(r$beta_full, 3), fp(r$p_full)) }), " \\\\"),
     row("Outcome", "Index", "$\\beta$ per SD (SE)", "$p$", "$q$", "$R^2_{\\mathrm{m}}$", "$\\beta$ full window", "$p$"), "llcccccc", "si_ecoacoustic_screen.tex",
     extra_head = "& & \\multicolumn{4}{c}{30-s region} & \\multicolumn{2}{c}{Full recording} \\\\ \\cmidrule(lr){3-6}\\cmidrule(lr){7-8}")
inc <- OUTC("analysis_17_ecoacoustic_indices", "increment_beyond_five_feature.csv") %>% mutate(outcome = factor(outcome, levels = OUT4)) %>% arrange(outcome, index)
frag(paste0(sapply(seq_len(nrow(inc)), function(i) { r <- inc[i, ]; row(OUTN[as.character(r$outcome)], ECN[r$index], fnum(r$beta, 3), fp(r$p), fp(r$q), fnum(r$r2m_base, 3), fnum(r$r2m_with, 3), fnum(r$d_aic, 1)) }), " \\\\"),
     row("Outcome", "Index", "$\\beta$ per SD", "$p$", "$q$", "$R^2_{\\mathrm{m}}$ base", "$R^2_{\\mathrm{m}}$ with index", "$\\Delta$AIC"), "llcccccc", "si_ecoacoustic_increment.tex")
# secondary ecoacoustic screen (A17): the 54 further scikit-maad indices, family 216
sec_raw <- OUTC("analysis_17_ecoacoustic_indices", "single_index_screen_secondary.csv")
stopifnot(all(sec_raw$family_size == 216), n_distinct(sec_raw$index) == 54)
sec <- sec_raw %>% mutate(cellv = paste0(fnum(beta, 3), " (", fp(q), ")")) %>% select(outcome, index, cellv) %>% tidyr::pivot_wider(names_from = outcome, values_from = cellv)
sec <- sec[match(unique(sec_raw$index), sec$index), ]
frag(paste0(sapply(seq_len(nrow(sec)), function(i) { r <- sec[i, ]; row(gsub("_", "\\\\_", r$index), sapply(OUT4, function(o) r[[o]])) }), " \\\\"),
     row("Index", OUTN[OUT4]), "lcccc", "si_ecoacoustic_secondary.tex")

# --- A18: detected vs perceived agreement (clip-level max probability), detected-group models, head-to-head
agd <- OUTC("analysis_18_detected_sources", "agreement_detected_vs_perceived.csv") %>% filter(score == "det")
agd <- agd %>% mutate(ord = match(source, c(SOURCES, "natural", "human", "mechanical"))) %>% arrange(ord)
frag(paste0(sapply(seq_len(nrow(agd)), function(i) { r <- agd[i, ]; nm <- ifelse(r$level == "group", paste(GRPN[r$source], "(group mean)"), SRCN[r$source])
  row(nm, fnum(r$perceived_mean, 2), fnum(100 * r$perceived_prevalence, 0), fnum(r$detected_mean, 3), fnum(r$rho, 2), fp(r$p), fp(r$q)) }), " \\\\"),
     row("Source", "Perceived", "Heard (\\%)", "Detected", "$\\rho$", "$p$", "$q$"), "lcccccc", "si_detected_agreement.tex")
dg <- OUTC("analysis_18_detected_sources", "detected_groups_models.csv") %>% mutate(outcome = factor(outcome, levels = OUT4)) %>% arrange(outcome)
DTN <- c(det_natural_z = "Detected natural", det_human_z = "Detected human", det_mechanical_z = "Detected mechanical")
frag(paste0(sapply(seq_len(nrow(dg)), function(i) { r <- dg[i, ]; row(OUTN[as.character(r$outcome)], DTN[r$term], bse(r$Estimate, r$`Std. Error`), fp(r$`Pr(>|t|)`), fp(r$q), fnum(r$r2m_level_only, 3), fnum(r$r2m_detected, 3), fnum(r$r2m_level_tonality_perceived, 3), fnum(r$r2m_both, 3)) }), " \\\\"),
     row("Outcome", "Term", "$\\beta$ per SD (SE)", "$p$", "$q$", "Level only", "+ detected", "+ perceived", "Both"), "llccccccc", "si_detected_models.tex",
     extra_head = "& & & & & \\multicolumn{4}{c}{$R^2_{\\mathrm{m}}$} \\\\ \\cmidrule(lr){6-9}")
hh <- OUTC("analysis_18_detected_sources", "head_to_head_detected_vs_perceived.csv") %>% mutate(outcome = factor(outcome, levels = OUT4)) %>% arrange(outcome)
HHN <- c(DTN, src_natural_z = "Perceived natural", src_human_z = "Perceived human", src_mechanical_z = "Perceived mechanical")
frag(paste0(sapply(seq_len(nrow(hh)), function(i) { r <- hh[i, ]; row(OUTN[as.character(r$outcome)], HHN[r$term], bse(r$Estimate, r$`Std. Error`), fp(r$`Pr(>|t|)`)) }), " \\\\"),
     row("Outcome", "Term", "$\\beta$ per SD (SE)", "$p$"), "llcc", "si_detected_head_to_head.tex")

# --- A19: cross-validated R2 for four ways of summarising the sources
cs <- OUTC("analysis_19_source_structure", "cv_r2_by_set.csv") %>% mutate(cellv = fnum(r2_env_mean, 2)) %>% select(outcome, set, scheme, cellv) %>% tidyr::pivot_wider(names_from = c(outcome, scheme), values_from = cellv)
SETN <- c(apriori3 = "A priori three groups", ward3 = "Ward three clusters", pca3 = "Three principal components", sources24 = "24 single sources")
cs <- cs[match(names(SETN), cs$set), ]
frag(paste0(sapply(seq_len(nrow(cs)), function(i) { r <- cs[i, ]; row(SETN[r$set], unlist(lapply(OUT4, function(o) c(r[[paste0(o, "_site")]], r[[paste0(o, "_evaluator")]])))) }), " \\\\"),
     row("Source set", rep(c("Site", "Evaluator"), 4)), "lcccccccc", "si_source_sets_cv.tex",
     extra_head = "& \\multicolumn{2}{c}{Comfort} & \\multicolumn{2}{c}{Richness} & \\multicolumn{2}{c}{Appropriateness} & \\multicolumn{2}{c}{Noisiness} \\\\ \\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}\\cmidrule(lr){8-9}")

# --- A03: zone omnibus and pairwise (Holm) for the four outcomes
zom <- OUTC("analysis_03_zone_profiles", "zone_omnibus.csv")
frag(paste0(sapply(seq_len(nrow(zom)), function(i) { r <- zom[i, ]; row(OUTN[r$outcome], fnum(r$chisq, 1), r$df, fp(r$p_lrt), fnum(r$r2_marginal, 3)) }), " \\\\"),
     row("Outcome", "$\\chi^2$", "df", "$p$ (LRT)", "$R^2_{\\mathrm{m}}$"), "lcccc", "si_zone_omnibus.tex")
emz <- OUTC("analysis_03_zone_profiles", "zone_emmeans.csv") %>% mutate(cellv = paste0(fnum(emmean, 2), " (", fci(lower.CL, upper.CL, 2), ")")) %>% select(zone, outcome, cellv) %>%
  tidyr::pivot_wider(names_from = outcome, values_from = cellv)
frag(paste0(sapply(ZONE_LEVELS, function(z) { r <- emz[emz$zone == z, ]; row(ZN(z), sapply(OUT4, function(o) r[[o]])) }), " \\\\"),
     row("Site", OUTN[OUT4]), "lcccc", "si_zone_emmeans.tex")
cat("all tables written\n")

# Established headline sensitivity checks (A20; no predictor search).
CVN <- c(primary="Primary", no_flagged="Without flags", no_trained="Without trained raters", without_E13_E28="Without two endpoint raters")
ce <- OUTC("analysis_20_storyline_checks", "core_exchange_sensitivity.csv")
crows <- unlist(lapply(c("wav","all"), function(ly) c(sub_head(ifelse(ly=="wav","30-s level","XL2 level, all dates"),5),
  paste0(sapply(names(CVN), function(v) {
    rr <- filter(ce, layer==ly, variant==v); a <- filter(rr,term=="src_natural"); b <- filter(rr,term=="src_mechanical")
    row(CVN[v],a$n,paste0(fnum(a$Estimate,3)," (",fci(a$lo,a$hi,3),")"),paste0(fnum(b$Estimate,3)," (",fci(b$lo,b$hi,3),")"),fnum(a$ratio,2))
  })," \\\\"))))
frag(crows,row("Variant","Ratings","Natural (95\\% CI)","Mechanical (95\\% CI)","Ratio"),"llccc","si_exchange_sensitivity.tex")
cc <- OUTC("analysis_20_storyline_checks", "core_context_sensitivity.csv")
crows <- unlist(lapply(c("wav","all"), function(ly) c(sub_head(ifelse(ly=="wav","30-s level","XL2 level, all dates"),6),
  paste0(sapply(names(CVN), function(v) {
    rr <- filter(cc,layer==ly,variant==v); a <- filter(rr,axis=="richness"); b <- filter(rr,axis=="comfort")
    row(CVN[v],a$n,fnum(a$chisq,1),fp(a$p),fnum(b$chisq,1),fp(b$p))
  })," \\\\"))))
frag(crows,row("Variant","Ratings","Richness $\\chi^2$","$p$","Comfort $\\chi^2$","$p$"),"llcccc","si_context_sensitivity.tex")

# Count versus total salience, on A15's common combination sample.
ct <- OUTC("analysis_15_richness_diversity", "richness_combination_coefficients.csv") %>%
  filter((model=="count" & term=="n_sources_z") | (model=="count_total" & term %in% c("n_sources_z","total_sal_z")))
frag(paste0(sapply(seq_len(nrow(ct)),function(i) {r<-ct[i,];row(ifelse(r$model=="count","Count alone","Count + total"),
      ifelse(r$term=="n_sources_z","Source count","Total salience"),paste0(fnum(r$Estimate,3)," (",fci(r$lo,r$hi,3),")"),fp(r$`Pr(>|t|)`))})," \\\\"),
     row("Richness model","Term","$\\beta$ per SD (95\\% CI)","$p$"),"llcc","si_richness_count_total.tex")

# --- A21: co-walkers' salience and the site-held-out hierarchy (added 16 October 2026)
A21 <- "analysis_21_cowalker_salience"
ff <- OUTC(A21, "five_feature_own_vs_cowalker.csv")
GR <- c(natural = "Natural sound", human = "Human sound", mechanical = "Mechanical sound")
c21 <- function(o, g, src, layer) {
  tm <- paste0(if (src == "own") "src_" else "cw_", g, "_z")
  r <- ff[ff$outcome == o & ff$term == tm & ff$salience == src & ff$layer == layer, ]; stopifnot(nrow(r) == 1); r }
cwm <- unlist(lapply(OUT4, function(o) c(sub_head(OUTN[o], 7), sapply(names(GR), function(g) {
  a <- c21(o, g, "own", "wav"); b <- c21(o, g, "cowalker", "wav"); x <- c21(o, g, "cowalker", "xl2_990")
  paste0(row(GR[g], bse(a$Estimate, a$`Std. Error`), fp(a$`Pr(>|t|)`), bse(b$Estimate, b$`Std. Error`), fp(b$q),
             bse(x$Estimate, x$`Std. Error`), fp(x$`Pr(>|t|)`)), " \\\\") }))))
frag(cwm, row("Source group", "Own $\\beta$ (SE)", "$p$", "Co-walkers' $\\beta$ (SE)", "$q$", "Co-walkers', 990 $\\beta$ (SE)", "$p$"),
     "lcccccc", "si_cowalker_models.tex")
ri <- OUTC(A21, "richness_interaction_own_vs_cowalker.csv")
SAL <- c(own = "Rater's own", cowalker = "Co-walkers'")
frag(sapply(c("own", "cowalker"), function(s) { r <- ri[ri$salience == s, ]
  paste0(row(SAL[s], fnum(r$b_int, 4), fci(r$lo, r$hi, 3), fnum(r$chisq, 1), fp(r$p), fnum(r$slope_p05, 3), fnum(r$slope_p95, 3)), " \\\\") }),
  row("Salience", "Interaction", "95\\% CI", "$\\chi^2$ (1 d.f.)", "$p$", "Natural slope at P5", "Natural slope at P95"),
  "lcccccc", "si_cowalker_richness.tex")
za <- OUTC(A21, "appropriateness_zone_own_vs_cowalker.csv")
zsl <- bind_rows(OUTC(A21, "zone_slopes_own.csv"), OUTC(A21, "zone_slopes_cowalker.csv"))
zcell <- function(sl, s, z) { r <- zsl[zsl$slope == sl & zsl$salience == s & zsl$zone == z, ]; stopifnot(nrow(r) == 1)
  paste0(fnum(r$trend, 2), " (", fci(r$lower.CL, r$upper.CL, 2), ")") }
zrow <- function(lab, f) paste0(row(lab, f("own"), f("cowalker")), " \\\\")
arows <- c(zrow("Comfort", function(s) bse(za$b_comfort[za$salience == s], za$se_comfort[za$salience == s], 2)),
           zrow("Richness", function(s) bse(za$b_richness[za$salience == s], za$se_richness[za$salience == s], 2)),
           zrow("Zone $\\times$ Comfort, $\\chi^2$ (7 d.f.), $p$", function(s) paste0(fnum(za$chisq_zone_comfort[za$salience == s], 1), ", ", fp(za$p_zone_comfort[za$salience == s]))),
           zrow("Zone $\\times$ Richness, $\\chi^2$ (7 d.f.), $p$", function(s) paste0(fnum(za$chisq_zone_richness[za$salience == s], 1), ", ", fp(za$p_zone_richness[za$salience == s]))),
           mid, sub_head("Richness slope by site (95\\% CI)", 3), sapply(ZONE_LEVELS, function(z) zrow(ZN(z), function(s) zcell("richness", s, z))),
           mid, sub_head("Comfort slope by site (95\\% CI)", 3), sapply(ZONE_LEVELS, function(z) zrow(ZN(z), function(s) zcell("comfort", s, z))))
frag(arows, row("Term", "Rater's own Comfort and Richness", "Co-walkers' mean Comfort and Richness"), "lcc", "si_cowalker_appropriateness.tex")
qd <- OUTC(A21, "quadrant_pairwise_own_vs_cowalker.csv")
qrows <- sapply(unique(qd$contrast), function(k) {
  q <- function(s) { r <- qd[qd$contrast == k & qd$salience == s, ]; c(bse(r$estimate, r$SE, 2), fp(r$p.value)) }
  display_k <- k
  for (key in names(QUAD_LAB)) display_k <- gsub(key, paste0("(", gsub("−", "$-$", gsub(" ", ", ", QUAD_LAB[[key]], fixed = TRUE), fixed = TRUE), ")"), display_k, fixed = TRUE)
  paste0(row(sub(" - ", " $-$ ", display_k), q("own"), q("cowalker")), " \\\\") })
frag(qrows, row("Contrast", "Own quadrants, estimate (SE)", "Holm $p$", "Co-walkers' quadrants, estimate (SE)", "Holm $p$"), "lcccc", "si_cowalker_quadrants.tex")
h <- OUTC(A21, "heldout_hierarchy.csv")
HL <- c(level = "$L_{\\mathrm{Aeq}}$", level_tonality = "$L_{\\mathrm{Aeq}}$ and tonality",
        detected = "\\quad + detected source groups", cowalker = "\\quad + co-walkers' source groups",
        own = "\\quad + rater's own source groups", cowalker_interaction = "\\quad + co-walkers' groups and natural $\\times$ mechanical",
        own_interaction = "\\quad + own groups and natural $\\times$ mechanical")
hrows <- sapply(names(HL), function(k) paste0(row(HL[k], sapply(OUT4, function(o) { v <- h$r2_stop_mean[h$features == k & h$outcome == o]
  if (length(v)) fnum(v, 3) else "--" })), " \\\\"))
frag(hrows, row("Predictors", "Comfort", "Richness", "Appropriateness", "Noisiness"), "lcccc", "si_heldout_hierarchy.tex")
