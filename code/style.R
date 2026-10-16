# =============================================================================
# style.R — shared plot style and data entry for the display items.
# Style: Helvetica, 8-10 pt at print size, left/bottom spines only, outward
# ticks, no grid, Okabe-Ito colours, PNG 600 dpi.
# Data: the deposited tables via code/helpers.R, and the result tables that the
# analysis modules write to output/<analysis>/.
# =============================================================================

source(file.path("code", "helpers.R"))
suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(ggtext); library(emmeans)
})

FIG <- file.path("output", "figures"); TAB <- file.path("output", "tables"); OUTP <- "output"
# Supplementary table bodies go to output/tables/si/; floats and captions are in
# the paper's Supplementary Information.
SITAB <- file.path(TAB, "si")
for (d in c(FIG, TAB, SITAB, OUTP)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
OUTC <- function(a, f) readr::read_csv(file.path("output", a, f),
                                       show_col_types = FALSE, guess_max = 5000)

mm2in <- function(mm) mm / 25.4
W1 <- mm2in(85); W2 <- mm2in(178)

# Okabe-Ito, eight colours (the eight zones need all eight); the six-colour
# house subset first, then yellow and black.
okabe_ito <- c("#0072B2", "#E69F00", "#009E73", "#D55E00", "#56B4E9", "#CC79A7",
               "#F0E442", "#000000")

# Zone order and labels (site names as in data/processed/p35_sites.csv).
ZONE_LAB <- c(boundary = "Entrance plaza", activity = "Tieshi Square",
              residential = "Dormitory 4", water = "Central waterway",
              sport = "Sports ground", green = "Paddy field",
              landscape = "Bawang Academy", teaching = "Teaching courtyard")
ZONE_FUNC <- c(boundary = "boundary", activity = "activity", residential = "residential",
               water = "waterside", sport = "sports", green = "green",
               landscape = "landscape", teaching = "teaching")
ZONE_COL <- setNames(okabe_ito[c(4, 2, 6, 1, 3, 7, 5, 8)], ZONE_LEVELS)
ZONE_SHAPE <- setNames(c(16, 17, 15, 18, 16, 17, 15, 18), ZONE_LEVELS)

OUT_LAB <- c(comfort = "Comfort", richness = "Richness",
             rating_appropriateness = "Appropriateness", rating_noisiness = "Noisiness")
OUT_COL <- c(comfort = "#0072B2", richness = "#E69F00",
             rating_appropriateness = "#009E73", rating_noisiness = "#D55E00")
GROUP_COL <- c(natural = "#009E73", human = "#E69F00", mechanical = "#D55E00")
QUAD_LEVELS <- c("calm", "vibrant", "monotonous", "chaotic")
QUAD_COL <- c(calm = "#56B4E9", vibrant = "#009E73", monotonous = "#CC79A7", chaotic = "#D55E00")
# Lighter tints of the quadrant colours for filled areas (F2d stacked bars); the
# saturated QUAD_COL stays on points and lines (F5b) so the hue mapping is shared.
lighten_hex <- function(col, f) {
  m <- grDevices::col2rgb(col)
  grDevices::rgb(t(m + (255 - m) * f), maxColorValue = 255)
}
QUAD_FILL <- setNames(lighten_hex(QUAD_COL, 0.42), names(QUAD_COL))

theme_pub <- function(base_size = 9, axis_title_size = 10) {
  theme_classic(base_size = base_size, base_family = "Helvetica") %+replace%
    theme(axis.line = element_line(colour = "black", linewidth = 0.4),
          axis.ticks = element_line(colour = "black", linewidth = 0.35),
          axis.ticks.length = unit(0.10, "cm"),
          axis.title = element_text(size = axis_title_size),
          axis.title.y = element_text(size = axis_title_size, angle = 90,
                                      margin = margin(r = 4)),
          axis.title.x = element_text(size = axis_title_size, margin = margin(t = 3)),
          axis.text = element_text(size = base_size, colour = "black"),
          legend.text = element_text(size = base_size),
          legend.title = element_blank(),
          legend.key = element_blank(),
          legend.key.size = unit(0.35, "cm"),
          legend.background = element_blank(),
          legend.margin = margin(0, 0, 0, 0),
          panel.grid = element_blank(),
          plot.title = element_blank(),
          plot.background = element_rect(fill = "white", colour = NA),
          panel.background = element_rect(fill = "white", colour = NA),
          plot.tag = element_text(size = 10, face = "bold", family = "Helvetica"),
          plot.margin = margin(4, 6, 4, 4))
}
sv <- function(p, name, w, h) {
  ggsave(file.path(FIG, name), p, width = w, height = h, dpi = 600, bg = "white",
         device = ragg::agg_png)
  cat("wrote", name, "\n")
}

# Plotmath labels (ISO 80000-8: italic quantity letter, roman subscript).
LAEQ <- 'italic(L)[plain("Aeq")]'
lab_laeq <- function(unit = "dB(A)", prefix = "") {
  parse(text = paste0(if (nzchar(prefix)) paste0('"', prefix, '"*') else "",
                      LAEQ, '*" (', unit, ')"'))
}

# Five-feature terms
FEAT_LEVELS <- c("level_z", "tonality_z", "src_natural_z", "src_human_z", "src_mechanical_z")
FEAT_LAB <- c(level_z = "*L*<sub>Aeq</sub>", tonality_z = "Tonality",
              src_natural_z = "Natural sound", src_human_z = "Human sound",
              src_mechanical_z = "Mechanical sound")
FEAT_LAB_PLAIN <- c(level_z = "LAeq", tonality_z = "Tonality",
                    src_natural_z = "Natural", src_human_z = "Human",
                    src_mechanical_z = "Mechanical")
