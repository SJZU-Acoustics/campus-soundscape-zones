# A10b — does the Comfort slope on appropriateness differ between zones? Symmetric with the zone × Richness test of
# Figure 3a: appropriateness ~ XL2 level + zone + Comfort + Richness against the same model with zone × Comfort
# (ML likelihood-ratio test, 7 d.f.); per-site Comfort slopes with 95% CIs from the interaction model.
# A single test, reported with its own p.
source(file.path("code", "helpers.R"))
suppressPackageStartupMessages(library(emmeans))
OUT <- file.path("output", "analysis_10_appropriateness")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
d <- load_p35("all")
m0 <- fit_backbone("xl2_laeq_z + zone + comfort + richness", "rating_appropriateness", d)
m1 <- fit_backbone("xl2_laeq_z + zone * comfort + richness", "rating_appropriateness", d)
lrt <- anova(m0, m1, refit = TRUE)
# the same test on top of the zone × Richness model of Fig. 5a (companion)
mr0 <- fit_backbone("xl2_laeq_z + zone * richness + comfort", "rating_appropriateness", d)
mr1 <- fit_backbone("xl2_laeq_z + zone * richness + zone * comfort", "rating_appropriateness", d)
lrt2 <- anova(mr0, mr1, refit = TRUE)
co <- summary(m0)$coefficients["comfort", ]
out_lrt <- tibble(test = c("zone x Comfort | zone + Comfort + Richness (primary)", "zone x Comfort | zone x Richness + Comfort (companion)"),
                  chisq = c(lrt$Chisq[2], lrt2$Chisq[2]), df = c(lrt$Df[2], lrt2$Df[2]), p = c(lrt$`Pr(>Chisq)`[2], lrt2$`Pr(>Chisq)`[2]),
                  common_slope = co[["Estimate"]], common_se = co[["Std. Error"]])
write_csv(out_lrt, file.path(OUT, "two_axis_zone_by_comfort_lrt.csv")); print(out_lrt)
tr <- as_tibble(emtrends(m1, "zone", var = "comfort")) %>% transmute(zone, slope = comfort.trend, se = SE, lo = lower.CL, hi = upper.CL)
write_csv(tr, file.path(OUT, "two_axis_zone_by_comfort_slopes.csv")); print(tr %>% mutate(across(where(is.numeric), ~round(.x, 3))))
