# Produces: output/tables/bayesian_corrected_models.csv, bayesian_leave_one_study_out.csv; output/intermediate/diagnostics/bayesian_posterior_predictive_checks.csv; output/intermediate/models/corrected_*.rds
# Modelos bayesianos sobre os dados corrigidos, comparados entre especificações.
#
# M2: binomial + spline (como o manuscrito), contagens corrigidas (85+ aberto, ponto médio 89),
#     só estudos com contagens (o ELSI não entra num binomial).
# M3: normal na escala logito com erro-padrão conhecido (inclui o ELSI ponderado), faixa como fator,
#     interceptos aleatórios de estudo e de faixa-dentro-de-estudo. Análogo bayesiano de 03a_fit_corrected_frequentist.R.
# M3 com três priores para os DPs (meia-Normal 0,5; meia-Normal 1; meia-t3 2,5) e para as médias.
# Diagnósticos: R-hat, ESS, divergências; checagem preditiva; influência (deixa um estudo de fora em M3).
# Estimandos de carga: estudo típico (efeitos = 0) e média populacional (integra os efeitos).
suppressPackageStartupMessages({ library(brms); library(dplyr); library(posterior) })
source(file.path("code", "helpers", "setup.R"))
options(mc.cores = 2)
SEED <- 20261019

agec <- read.csv(file.path(DERIVED, "age_strata_corrected.csv"))
cen  <- read.csv(file.path(OCT_IBGE, "ibge_census_brazil_2022_agegroups.csv"))
bands <- c("60-64", "65-69", "70-74", "75-79", "80-84", "85+")
pop6 <- c(cen$population[match(bands[1:5], cen$age_group)], sum(cen$population[cen$age_group %in% c("85-89", "90+")]))
mid6 <- c(62, 67, 72, 77, 82, 89)
agec$band <- factor(agec$band, levels = bands); agec$sei <- sqrt(agec$vi)
agec$age_mid <- mid6[as.integer(agec$band)]

diag_of <- function(m) {
  s <- summarise_draws(as_draws_df(m), "rhat", "ess_bulk", "ess_tail")
  s <- s[!grepl("^lp__|^lprior|^r_|^zs_|^z_", s$variable), ]
  data.frame(max_rhat = max(s$rhat, na.rm = TRUE), min_ess_bulk = min(s$ess_bulk, na.rm = TRUE),
             min_ess_tail = min(s$ess_tail, na.rm = TRUE), divergent = sum(subset(nuts_params(m), Parameter == "divergent__")$Value))
}
fit <- function(formula, data, family, prior, file) brm(formula, data = data, family = family, prior = prior, chains = 4,
  iter = 4000, warmup = 2000, cores = 2, seed = SEED, control = list(adapt_delta = 0.99), backend = "rstan",
  file = file.path(OUT_MOD, file), refresh = 0)

## M2 binomial spline nos dados corrigidos
cnt <- subset(agec, !is.na(cases))
m2 <- fit(bf(cases | trials(n) ~ s(age_mid, bs = "tp", k = 5) + (1 | study)), cnt, binomial(),
          c(prior(normal(0, 5), class = "Intercept"), prior(normal(0, 5), class = "b"), prior(normal(0, 2.5), class = "sd")),
          "corrected_m2_binomial_spline")
g2 <- data.frame(age_mid = mid6, n = 1L, study = "new")
e2 <- posterior_epred(m2, newdata = g2, re_formula = NA)

## M3 normal-logito com SE conhecido, três priores
pri3 <- list(
  hn05 = c(prior(normal(-2.5, 2), class = "b"), prior(normal(0, 0.5), class = "sd")),
  hn1  = c(prior(normal(-2.5, 2), class = "b"), prior(normal(0, 1), class = "sd")),
  ht25 = c(prior(normal(-2.5, 2), class = "b"), prior(student_t(3, 0, 2.5), class = "sd")))
f3 <- bf(yi | se(sei) ~ 0 + band + (1 | study) + (1 | study:row))
m3 <- lapply(names(pri3), function(p) fit(f3, agec, gaussian(), pri3[[p]], paste0("corrected_m3_", p)))
names(m3) <- names(pri3)

# brms renomeia os níveis do fator: "60-64" -> "b_band60M64", "85+" -> "b_band85P"; a ordem segue os níveis
band_vars <- function(m) { v <- grep("^b_band", variables(m), value = TRUE); stopifnot(length(v) == 6); v }
burden_typical <- function(m) {
  b <- as_draws_matrix(m, variable = band_vars(m))
  as.vector(plogis(b) %*% pop6)
}
burden_popavg <- function(m, nodes = 200) {
  dr <- as_draws_df(m); b <- as.matrix(dr[, band_vars(m)])
  ss <- dr$sd_study__Intercept; sr <- dr[["sd_study:row__Intercept"]]
  z1 <- rnorm(nodes); z2 <- rnorm(nodes)
  sapply(seq_len(nrow(b)), function(i) sum(pop6 * colMeans(plogis(outer(z1 * ss[i] + z2 * sr[i], b[i, ], "+")))))
}
q <- function(x) quantile(x, c(.5, .025, .975))
set.seed(SEED)
rows <- list(
  data.frame(model = "M2 binomial spline, counts only (no ELSI)", estimand = "typical study", t(q(as.vector(e2 %*% pop6))), diag_of(m2)))
for (p in names(m3)) {
  rows[[length(rows) + 1]] <- data.frame(model = paste("M3 logit-normal, prior SD", p), estimand = "typical study", t(q(burden_typical(m3[[p]]))), diag_of(m3[[p]]))
  rows[[length(rows) + 1]] <- data.frame(model = paste("M3 logit-normal, prior SD", p), estimand = "population-averaged", t(q(burden_popavg(m3[[p]]))), diag_of(m3[[p]]))
}
res <- bind_rows(rows); names(res)[3:5] <- c("cases_median", "cases_lb", "cases_ub")
write.csv(res, file.path(OUT_TAB, "bayesian_corrected_models.csv"), row.names = FALSE)
print(res, digits = 4)

## Checagem preditiva (estatística: média e DP dos yi; e cobertura por linha)
pp <- posterior_predict(m3$hn1)
ppc <- data.frame(stat = c("mean yi", "sd yi"), observed = c(mean(agec$yi), sd(agec$yi)),
                  ppp = c(mean(rowMeans(pp) > mean(agec$yi)), mean(apply(pp, 1, sd) > sd(agec$yi))))
cover <- mean(agec$yi >= apply(pp, 2, quantile, .025) & agec$yi <= apply(pp, 2, quantile, .975))
ppc <- rbind(ppc, data.frame(stat = "row-level 95% predictive coverage", observed = cover, ppp = NA))
pp2 <- posterior_predict(m2); cover2 <- mean(cnt$cases >= apply(pp2, 2, quantile, .025) & cnt$cases <= apply(pp2, 2, quantile, .975))
ppc <- rbind(ppc, data.frame(stat = "M2 row-level 95% predictive coverage (counts)", observed = cover2, ppp = NA))
write.csv(ppc, file.path(OUT_DIAG, "bayesian_posterior_predictive_checks.csv"), row.names = FALSE)
print(ppc)

## Influência: deixa um estudo de fora (M3, prior hn1)
infl <- bind_rows(lapply(unique(agec$study), function(s) {
  d <- subset(agec, study != s); d$row <- seq_len(nrow(d))
  m <- update(m3$hn1, newdata = d, refresh = 0, seed = SEED, cores = 2)
  data.frame(omitted = s, t(q(burden_typical(m))), diag_of(m))
}))
names(infl)[2:4] <- c("cases_median", "cases_lb", "cases_ub")
write.csv(infl, file.path(OUT_TAB, "bayesian_leave_one_study_out.csv"), row.names = FALSE)
print(infl, digits = 4)
writeLines(capture.output(sessionInfo()), file.path(OUT_DIAG, "sessioninfo_13_fit_bayesian_corrected.txt"))
