# Produces: output/tables/primary_bayesian.csv; output/intermediate/diagnostics/primary_bayesian_ppc.csv; output/intermediate/models/primary_m3_*.rds
# Modelo bayesiano logit-normal (M3, prior meia-Normal(0, 1) nos DPs) no conjunto principal sem sobreposição
# (PRIMARY_POPULATIONS, code/helpers/setup.R), nas versões observada e melhor estimativa dos autores.
# Mesma especificação de 03i_fit_bayesian_corrected.R; diagnósticos e checagem preditiva.
suppressPackageStartupMessages({ library(brms); library(dplyr); library(posterior) })
source(file.path("code", "helpers", "setup.R"))
options(mc.cores = 2)
SEED <- 20261019
agec_all <- read.csv(file.path(DERIVED, "age_strata_corrected_with_epifloripa.csv"))   # as oito populações
agec7 <- subset(agec_all, study %in% PRIMARY_POPULATIONS)   # análise principal sem sobreposição (code/helpers/setup.R)
best  <- read.csv(file.path(DERIVED, "author_best_estimates_s6.csv"))
cen <- read.csv(file.path(OCT_IBGE, "ibge_census_brazil_2022_agegroups.csv"))
bands <- c("60-64", "65-69", "70-74", "75-79", "80-84", "85+")
pop6 <- c(cen$population[match(bands[1:5], cen$age_group)], sum(cen$population[cen$age_group %in% c("85-89", "90+")]))
author_best_rows <- function(d) {
  for (s in c("Bottino_2008", "Lopes_2012", "EpiFloripa_2024")) {
    f <- best$factor[best$study == s]; i <- d$study == s & !is.na(d$cases)
    pp <- pmin(pmax(d$cases[i], 0.5) / d$n[i] * f, 0.95); d$yi[i] <- qlogis(pp); d$vi[i] <- 1 / (d$n[i] * pp * (1 - pp))
  }
  d
}
versions <- list(observed = agec7, author_best = author_best_rows(agec7))
pri <- c(prior(normal(-2.5, 2), class = "b"), prior(normal(0, 1), class = "sd"))
f3 <- bf(yi | se(sei) ~ 0 + band + (1 | study) + (1 | study:row))
q <- function(x) quantile(x, c(.5, .025, .975))
out <- list(); ppc <- list()
for (v in names(versions)) {
  d <- versions[[v]]; d$band <- factor(d$band, levels = bands); d$sei <- sqrt(d$vi)
  m <- brm(f3, data = d, family = gaussian(), prior = pri, chains = 4, iter = 4000, warmup = 2000, cores = 2, seed = SEED,
           control = list(adapt_delta = 0.99), backend = "rstan", file = file.path(OUT_MOD, paste0("primary_m3_", v)), refresh = 0)
  bv <- grep("^b_band", variables(m), value = TRUE); dr <- as_draws_df(m); b <- as.matrix(dr[, bv])
  ss <- dr$sd_study__Intercept; sr <- dr[["sd_study:row__Intercept"]]; set.seed(SEED); z1 <- rnorm(200); z2 <- rnorm(200)
  pa <- sapply(seq_len(nrow(b)), function(i) sum(pop6 * colMeans(plogis(outer(z1 * ss[i] + z2 * sr[i], b[i, ], "+")))))
  s <- summarise_draws(as_draws_df(m), "rhat", "ess_bulk", "ess_tail"); s <- s[!grepl("^lp__|^lprior|^r_|^z_", s$variable), ]
  dg <- data.frame(max_rhat = max(s$rhat, na.rm = TRUE), min_ess_bulk = min(s$ess_bulk, na.rm = TRUE), min_ess_tail = min(s$ess_tail, na.rm = TRUE),
                   divergent = sum(subset(nuts_params(m), Parameter == "divergent__")$Value))
  out[[v]] <- rbind(data.frame(version = v, estimand = "typical study", t(q(as.vector(plogis(b) %*% pop6))), dg),
                    data.frame(version = v, estimand = "population-averaged", t(q(pa)), dg))
  pp <- posterior_predict(m)
  ppc[[v]] <- data.frame(version = v, ppp_mean = mean(rowMeans(pp) > mean(d$yi)), ppp_sd = mean(apply(pp, 1, sd) > sd(d$yi)),
                         coverage95 = mean(d$yi >= apply(pp, 2, quantile, .025) & d$yi <= apply(pp, 2, quantile, .975)))
}
out <- bind_rows(out); names(out)[3:5] <- c("cases_median", "cases_lb", "cases_ub"); ppc <- bind_rows(ppc)
print(out, digits = 4); print(ppc, digits = 3)
write.csv(out, file.path(OUT_TAB, "primary_bayesian.csv"), row.names = FALSE)
write.csv(ppc, file.path(OUT_DIAG, "primary_bayesian_ppc.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(OUT_DIAG, "sessioninfo_16_fit_bayesian_primary_set.txt"))
