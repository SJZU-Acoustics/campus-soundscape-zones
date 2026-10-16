# helpers.R — shared data entry, variable sets and model helpers.
# Every analysis module source()s this file. The tables come from the deposited
# workbook through code/load_data.R; analysis modules never re-filter or
# re-derive design columns.
#
# Samples:
#   "all"  990 records — XL2 levels, questionnaire, visual (2 July has no B&K block)
#   "wav"  666 records — the B&K 30-s block present (primary acoustic layer)

suppressPackageStartupMessages({
  library(tidyverse)
  library(lmerTest)
  library(broom)
  library(broom.mixed)
})

# Pin marginal-mean/trend inference instead of inheriting package/session defaults.
# A03 explicitly uses Satterthwaite; other analyses use Kenward-Roger.
# Required dependencies and unlimited size thresholds prevent silent fallback.
stopifnot(requireNamespace("pbkrtest", quietly = TRUE))
emmeans::emm_options(lmer.df = "kenward-roger", disable.pbkrtest = FALSE,
                     pbkrtest.limit = Inf, disable.lmerTest = FALSE,
                     lmerTest.limit = Inf)

source(file.path("code", "load_data.R"))

ZONE_LEVELS <- c("boundary", "activity", "residential", "water", "sport", "green", "landscape", "teaching")
SLOT_LEVELS <- c("06-09", "09-12", "12-15", "15-18", "18-21")

OUTCOMES <- tribble(
  ~var,                     ~label,                        ~role,
  "comfort",                "Comfort (ICSAQ)",             "focal",
  "richness",               "Richness (ICSAQ)",            "focal",
  "rating_appropriateness", "Sound-place appropriateness", "focal",
  "rating_noisiness",       "Perceived noisiness",         "comparison",
  "rating_auditory_env",    "Auditory environment",        "comparison",
  "rating_visual_env",      "Visual environment",          "comparison",
)

ADJECTIVES <- c("comfortable", "chaotic", "interesting", "boring",
                "rich", "monotonous", "lively", "plain")

SOURCES <- c("wind", "water", "insects", "birds", "animals", "leaves",
             "conversation", "children", "playing", "footsteps", "running",
             "teaching", "broadcast", "bike_voice", "vehicles", "vehicle_bumps",
             "rail", "sweeping", "mower", "ac_units", "loading", "luggage",
             "gates", "construction")
SOURCE_GROUPS <- list(
  natural    = c("wind", "water", "insects", "birds", "animals", "leaves"),
  human      = c("conversation", "children", "playing", "footsteps", "running",
                 "teaching", "broadcast"),
  mechanical = c("bike_voice", "vehicles", "vehicle_bumps", "rail", "sweeping",
                 "mower", "ac_units", "loading", "luggage", "gates", "construction")
)

ACOUSTIC <- tribble(
  ~var,        ~label,                          ~family,    ~layer,
  "bk_laeq",   "LAeq (30-s region)",            "level",    "bk",
  "bk_l5",     "L5 (30-s)",                     "level",    "bk",
  "bk_l95",    "L95 (30-s)",                    "level",    "bk",
  "bk_dl5_95", "L5-L95 (30-s)",                 "dynamics", "bk",
  "n",         "Loudness N (whole region)",     "loudness", "bk",
  "n5",        "Loudness N5",                   "loudness", "bk",
  "n5_n95",    "N5/N95",                        "dynamics", "bk",
  "s",         "Sharpness S",                   "spectral", "bk",
  "r",         "Roughness R",                   "temporal", "bk",
  "f",         "Fluctuation strength F",        "temporal", "bk",
  "tonality",  "Tonality",                      "spectral", "bk",
  "xl2_laeq",  "LAeq (XL2, full measurement)",  "level",    "xl2",
  "xl2_dl5_95","LAF5-LAF95 (XL2)",              "dynamics", "xl2",
)

VISUAL <- c("vis_sky", "vis_vegetation", "vis_water", "vis_ground", "vis_buildings",
            "vis_pavement", "vis_people", "vis_vehicles", "vis_other",
            "naturalness", "artificiality", "vri", "vdi", "vei", "vdo")

read_proc <- function(name) read_frozen_csv(name, guess_max = 2000)

load_p35_people <- function() {
  read_proc("p35_people.csv") %>%
    mutate(sex = factor(sex, levels = c(1, 2), labels = c("male", "female")),
           education = factor(education, levels = c(5, 6), labels = c("undergraduate", "masters")),
           cohort = factor(cohort, levels = c("B-1", "B-2", "3", "4", "5", "R")),
           rater_type = factor(rater_type, levels = c("naive", "trained")))
}

load_p35_env <- function() {
  read_proc("p35_env.csv") %>%
    mutate(site = factor(site), daytype = factor(daytype, levels = c("weekday", "weekend")))
}

load_p35_sessions <- function() read_proc("p35_sessions.csv")
load_p35_sites <- function() read_proc("p35_sites.csv")

load_p35 <- function(sample = c("all", "wav")) {
  sample <- match.arg(sample)
  people <- load_p35_people() %>%
    select(evaluator, sex, age_band, hukou_liaoning, education, cohort, rater_type,
           who5_sum, who5_low, ns_sum)
  d <- read_proc("p35_records.csv") %>%
    left_join(people, by = "evaluator") %>%
    mutate(
      evaluator = factor(evaluator),
      env_id = factor(env_id),
      session = factor(session),
      site = factor(site),
      zone = factor(zone, levels = ZONE_LEVELS),
      daytype = factor(daytype, levels = c("weekday", "weekend")),
      slot = factor(slot, levels = SLOT_LEVELS),
      wav_recorded = as.integer(wav_recorded),
      flagged = flags != "" & !is.na(flags),
      src_natural = rowMeans(across(all_of(paste0("src_", SOURCE_GROUPS$natural)))),
      src_human = rowMeans(across(all_of(paste0("src_", SOURCE_GROUPS$human)))),
      src_mechanical = rowMeans(across(all_of(paste0("src_", SOURCE_GROUPS$mechanical))))
    )
  # z-scores for the acoustic indicators, computed within the layer's own sample
  for (v in ACOUSTIC$var) d[[paste0(v, "_z")]] <- as.numeric(scale(d[[v]]))
  if (sample == "wav") d <- filter(d, wav_recorded == 1)
  d
}

# BH-FDR over a tidy coefficient table: adds q and the stated family size
bh_family <- function(tab, p_col = "p.value") {
  tab %>% mutate(q = p.adjust(.data[[p_col]], method = "BH"), family_size = n())
}

# Resampling can produce difficult variance-component fits. Keep the original
# optimiser where it converges; retry only numerical failures, without dropping
# a bootstrap/permutation draw. A zero variance is a valid boundary fit.
checked_resample_lmer <- function(...) {
  args <- list(...)
  bad <- function(m) {
    msg <- m@optinfo$conv$lme4$messages
    any(!grepl("boundary.*singular", msg)) || any(m@optinfo$conv$opt != 0)
  }
  m <- suppressWarnings(suppressMessages(do.call(lme4::lmer, args)))
  if (!bad(m)) return(m)
  original <- paste(m@optinfo$conv$lme4$messages, collapse = "; ")
  old_ll <- as.numeric(logLik(m))
  for (optimizer in c("bobyqa", "Nelder_Mead")) {
    args$control <- lme4::lmerControl(optimizer = optimizer, optCtrl = list(maxfun = 200000))
    refit <- suppressWarnings(suppressMessages(do.call(lme4::lmer, args)))
    if (!bad(refit)) {
      cat("Resolved resampling convergence with", optimizer, "; logLik change",
          signif(as.numeric(logLik(refit)) - old_ll, 4), "; initial:", original, "\n")
      return(refit)
    }
  }
  stop("Unresolved numerical convergence in a resampling draw: ", original)
}

# Resample complete stops in draw order, giving each repeated draw its own
# random-effect identifier. This is equivalent to binding one mutated frame per
# sampled stop, but avoids thousands of small data-frame constructions.
resample_stops <- function(data, draws) {
  indices <- split(seq_len(nrow(data)), as.character(data$env_id))
  selected <- indices[as.character(draws)]
  stopifnot(all(lengths(selected) > 0))
  out <- data[unlist(selected, use.names = FALSE), ]
  out$env_id <- rep(paste0(as.character(draws), "_", seq_along(draws)), lengths(selected))
  out
}

# Default backbone: crossed random intercepts for evaluator and recording.
fit_backbone <- function(formula_rhs, outcome, data) {
  f <- as.formula(paste0(outcome, " ~ ", formula_rhs, " + (1 | evaluator) + (1 | env_id)"))
  lmerTest::lmer(f, data = data, REML = TRUE)
}
