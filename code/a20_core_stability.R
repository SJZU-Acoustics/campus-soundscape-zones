# A20 — fixed sensitivity checks of the headline results: the natural–mechanical exchange ratio for Comfort and the
# zone dependence of the Comfort and Richness slopes on appropriateness, on both level layers, with and without
# flagged records, trained raters and the two endpoint raters.
source(file.path("code", "helpers.R"))
OUT <- file.path("output", "analysis_20_storyline_checks")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
da <- load_p35("all"); dw <- load_p35("wav")
subsets <- function(d) list(primary=d, no_flagged=filter(d,!flagged),
  no_trained=filter(d,rater_type=="naive"),
  without_E13_E28=filter(d,!evaluator %in% c("E13","E28")))
exchange <- list(); context <- list(); slopes <- list()
for (layer in c("wav","all")) {
  variants <- subsets(if (layer=="wav") dw else da)
  level <- if (layer=="wav") "bk_laeq" else "xl2_laeq"
  for (v in names(variants)) {
    x <- variants[[v]]
    m <- fit_backbone(paste(level,"+ src_natural + src_mechanical + src_human"),"comfort",x)
    z <- as_tibble(summary(m)$coefficients,rownames="term") %>%
      filter(term %in% c("src_natural","src_mechanical")) %>%
      mutate(layer=layer,variant=v,n=nobs(m),lo=Estimate-qt(.975,df)*`Std. Error`,hi=Estimate+qt(.975,df)*`Std. Error`,
        ratio=-fixef(m)["src_mechanical"]/fixef(m)["src_natural"])
    exchange[[paste(layer,v)]] <- z
    m0 <- fit_backbone(paste(level,"+ zone + comfort + richness"),"rating_appropriateness",x)
    mr <- fit_backbone(paste(level,"+ zone * richness + comfort"),"rating_appropriateness",x)
    mc <- fit_backbone(paste(level,"+ zone * comfort + richness"),"rating_appropriateness",x)
    for (axis in c("richness","comfort")) {
      fit <- if (axis=="richness") mr else mc
      test <- anova(m0,fit,refit=TRUE)
      context[[paste(layer,v,axis)]] <- tibble(layer=layer,variant=v,axis=axis,n=nobs(fit),chisq=test$Chisq[2],df=test$Df[2],p=test$`Pr(>Chisq)`[2])
      tr <- as_tibble(emmeans::emtrends(fit,~zone,var=axis))
      names(tr)[names(tr)==paste0(axis,".trend")] <- "slope"
      slopes[[paste(layer,v,axis)]] <- tr %>% mutate(layer=layer,variant=v,axis=axis,n=nobs(fit))
    }
  }
}
write_csv(bind_rows(exchange),file.path(OUT,"core_exchange_sensitivity.csv"))
write_csv(bind_rows(context),file.path(OUT,"core_context_sensitivity.csv"))
write_csv(bind_rows(slopes),file.path(OUT,"core_zone_slopes_sensitivity.csv"))
print(bind_rows(exchange) %>% filter(term=="src_natural") %>% select(layer,variant,n,ratio))
print(bind_rows(context),n=16)
