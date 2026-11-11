# Produces: output/tables/primary_loo_prior_sensitivity.csv; output/intermediate/models/primary_m3_author_best_*.rds
# Sensibilidade da carga de 2022 no conjunto principal (PRIMARY_POPULATIONS, code/helpers/setup.R), com as
# estimativas dos autores (mesma regra de 03g/03j/03l):
#   (1) deixar uma população de fora de cada vez, no modelo etário frequentista (03j) e no bayesiano (03l);
#   (2) priors alternativos para os DPs dos efeitos aleatórios no modelo bayesiano: meia-Normal(0; 0,5) e
#       meia-t(3; 0; 2,5), como no conjunto R4 de 03i, ao lado da meia-Normal(0; 1) da análise principal.
# Estimando: população típica (efeitos zerados) e média populacional, aplicados ao Censo 2022.
suppressPackageStartupMessages({ library(metafor); library(brms); library(dplyr); library(posterior); library(MASS, exclude = "select") })
source(file.path("code", "helpers", "setup.R"))
options(mc.cores = 2)
SEED <- 20261019; set.seed(SEED)
agec_all <- read.csv(file.path(DERIVED, "age_strata_corrected_with_epifloripa.csv"))
best <- read.csv(file.path(DERIVED, "author_best_estimates_s6.csv"))
cen <- read.csv(file.path(OCT_IBGE, "ibge_census_brazil_2022_agegroups.csv"))
bands <- c("60-64", "65-69", "70-74", "75-79", "80-84", "85+")
pop6 <- c(cen$population[match(bands[1:5], cen$age_group)], sum(cen$population[cen$age_group %in% c("85-89", "90+")]))
author_best_rows <- function(d) {   # mesma regra de 03g (efeito de desenho 1)
  for (s in c("Bottino_2008", "Lopes_2012", "EpiFloripa_2024")) {
    f <- best$factor[best$study == s]; i <- d$study == s & !is.na(d$cases)
    pp <- pmin(pmax(d$cases[i], 0.5) / d$n[i] * f, 0.95); d$yi[i] <- qlogis(pp); d$vi[i] <- 1 / (d$n[i] * pp * (1 - pp))
  }
  d
}
d0 <- author_best_rows(subset(agec_all, study %in% PRIMARY_POPULATIONS))
d0$band <- factor(d0$band, levels = bands); d0$sei <- sqrt(d0$vi)
q <- function(x) quantile(x, c(.5, .025, .975), names = FALSE)
label <- c(Herrera_2002 = "Catanduva", Lopes_2012 = "Ribeirão Preto", Scazufca_2008 = "São Paulo SPAH",
           Cesar_2015 = "Tremembé", EpiFloripa_2024 = "Florianópolis")
omit <- c("none", PRIMARY_POPULATIONS)

# (1a) frequentista
fit_age <- function(d) rma.mv(yi, vi, mods = ~ 0 + band, random = ~ 1 | study / row, data = d, method = "REML")
freq <- bind_rows(lapply(omit, function(o) {
  d <- subset(d0, study != o); m <- fit_age(d); b <- coef(m); s2 <- m$sigma2
  set.seed(SEED); dr <- mvrnorm(20000, b, vcov(m)); z <- matrix(rnorm(400 * 2), ncol = 2)
  pa <- sapply(1:6, function(j) sapply(1:2000, function(i) mean(plogis(dr[i, j] + z[, 1] * sqrt(s2[1]) + z[, 2] * sqrt(s2[2])))))
  rbind(data.frame(model = "frequentist", analysis = "leave one out", omitted = o, prior_sd = "NA", estimand = "typical study",
                   t(q(plogis(dr) %*% pop6)), k = m$s.nlevels[1]),
        data.frame(model = "frequentist", analysis = "leave one out", omitted = o, prior_sd = "NA", estimand = "population-averaged",
                   t(q(pa %*% pop6)), k = m$s.nlevels[1]))
}))

# (1b, 2) bayesiano: mesma especificação de 03l
f3 <- bf(yi | se(sei) ~ 0 + band + (1 | study) + (1 | study:row))
priors <- list("half-Normal(0, 1)" = c(prior(normal(-2.5, 2), class = "b"), prior(normal(0, 1), class = "sd")),
               "half-Normal(0, 0.5)" = c(prior(normal(-2.5, 2), class = "b"), prior(normal(0, 0.5), class = "sd")),
               "half-t(3, 0, 2.5)" = c(prior(normal(-2.5, 2), class = "b"), prior(student_t(3, 0, 2.5), class = "sd")))
tag <- c("half-Normal(0, 1)" = "hn1", "half-Normal(0, 0.5)" = "hn05", "half-t(3, 0, 2.5)" = "ht25")
summ <- function(m, analysis, o, p) {
  bv <- grep("^b_band", variables(m), value = TRUE); dr <- as_draws_df(m); b <- as.matrix(dr[, bv])
  ss <- dr$sd_study__Intercept; sr <- dr[["sd_study:row__Intercept"]]; set.seed(SEED); z1 <- rnorm(200); z2 <- rnorm(200)
  pa <- sapply(seq_len(nrow(b)), function(i) sum(pop6 * colMeans(plogis(outer(z1 * ss[i] + z2 * sr[i], b[i, ], "+")))))
  s <- summarise_draws(dr, "rhat", "ess_bulk"); s <- s[!grepl("^lp__|^lprior|^r_|^z_", s$variable), ]
  dg <- data.frame(max_rhat = max(s$rhat, na.rm = TRUE), min_ess_bulk = min(s$ess_bulk, na.rm = TRUE),
                   divergent = sum(subset(nuts_params(m), Parameter == "divergent__")$Value))
  rbind(data.frame(model = "Bayesian", analysis = analysis, omitted = o, prior_sd = p, estimand = "typical study", t(q(as.vector(plogis(b) %*% pop6))), dg),
        data.frame(model = "Bayesian", analysis = analysis, omitted = o, prior_sd = p, estimand = "population-averaged", t(q(pa)), dg))
}
bay <- list()
for (p in names(priors)) {
  m <- brm(f3, data = d0, family = gaussian(), prior = priors[[p]], chains = 4, iter = 4000, warmup = 2000, cores = 2, seed = SEED,
           control = list(adapt_delta = 0.99), backend = "rstan", refresh = 0,
           file = file.path(OUT_MOD, if (p == "half-Normal(0, 1)") "primary_m3_author_best" else paste0("primary_m3_author_best_", tag[[p]])))
  bay[[length(bay) + 1]] <- summ(m, "prior", "none", p)
  if (p == "half-Normal(0, 1)") for (o in PRIMARY_POPULATIONS) {
    mo <- brm(f3, data = subset(d0, study != o), family = gaussian(), prior = priors[[p]], chains = 4, iter = 4000, warmup = 2000,
              cores = 2, seed = SEED, control = list(adapt_delta = 0.99), backend = "rstan", refresh = 0,
              file = file.path(OUT_MOD, paste0("primary_m3_author_best_loo_", o)))
    bay[[length(bay) + 1]] <- summ(mo, "leave one out", o, p)
  }
}
bay <- bind_rows(bay)
out <- bind_rows(freq, bay); names(out)[6:8] <- c("cases_median", "cases_lb", "cases_ub")
out$omitted_label <- ifelse(out$omitted == "none", "none", label[out$omitted])
print(out, digits = 4)
# a linha sem omissão do frequentista tem de reproduzir a análise principal (primary_burden_estimands.csv)
ref <- read.csv(file.path(OUT_TAB, "primary_burden_estimands.csv"))
r0 <- ref$cases_median[ref$version == "author_best" & ref$estimand == "typical study"]
stopifnot(abs(out$cases_median[out$model == "frequentist" & out$omitted == "none" & out$estimand == "typical study"] / r0 - 1) < 0.01)
write.csv(out, file.path(OUT_TAB, "primary_loo_prior_sensitivity.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(OUT_DIAG, "sessioninfo_03r_loo_prior_sensitivity.txt"))
