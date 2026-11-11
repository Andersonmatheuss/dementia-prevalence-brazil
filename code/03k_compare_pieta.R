# Produces: output/tables/pieta_out_of_sample.csv
# Comparação fora da amostra com o estudo Pietà (Caeté, MG; 75+ anos), achado na busca.
# O Pietà não entra no agrupamento: a prevalência global só foi publicada em resumos de congresso,
# que o critério da revisão exclui (Table S16, categoria 1). Aqui ele serve de teste externo da curva
# etária do conjunto principal (PRIMARY_POPULATIONS, code/helpers/setup.R): prevalência bruta esperada para a estrutura etária do Pietà, contra a
# observada.
# Fontes (conferidas no texto):
#   estrutura etária: Caramelli et al. 2011, Arq Neuropsiquiatr 69:579 (PMID 21877023), Discussão:
#     75-79 47,4%; 80-84 29,6%; 85-89 14,1%; 90+ 8,9%
#   demência: 174/639 (Alzheimers Dement 2013;9(4S):P4-152, doi 10.1016/j.jalz.2013.05.1543)
#             176/639 (Alzheimers Dement 2009;5(4S):P3-163, doi 10.1016/j.jalz.2009.04.1037)
suppressPackageStartupMessages({ library(metafor); library(dplyr); library(MASS, exclude = "select") })
source(file.path("code", "helpers", "setup.R"))
set.seed(20261019)
agec_all <- read.csv(file.path(DERIVED, "age_strata_corrected_with_epifloripa.csv"))   # as oito populações
agec7 <- subset(agec_all, study %in% PRIMARY_POPULATIONS)   # análise principal sem sobreposição (code/helpers/setup.R)
best  <- read.csv(file.path(DERIVED, "author_best_estimates_s6.csv"))
bands <- c("60-64", "65-69", "70-74", "75-79", "80-84", "85+")
fit_age <- function(d) { d$band <- factor(d$band, levels = bands); rma.mv(yi, vi, mods = ~ 0 + band, random = ~ 1 | study / row, data = d, method = "REML") }
author_best_rows <- function(d) {   # mesma regra de 03g e 03j (efeito de desenho 1)
  for (s in c("Bottino_2008", "Lopes_2012", "EpiFloripa_2024")) {
    f <- best$factor[best$study == s]; i <- d$study == s & !is.na(d$cases)
    pp <- pmin(pmax(d$cases[i], 0.5) / d$n[i] * f, 0.95); d$yi[i] <- qlogis(pp); d$vi[i] <- 1 / (d$n[i] * pp * (1 - pp))
  }
  d
}
versions <- list(observed = agec7, author_best = author_best_rows(agec7))

mix <- c(0, 0, 0, 0.474, 0.296, 0.141 + 0.089); stopifnot(abs(sum(mix) - 1) < 1e-9)
obs <- data.frame(source = c("abstract 2013 (174/639)", "abstract 2009 (176/639)"), x = c(174, 176), n = 639)
q <- function(x) quantile(x, c(.5, .025, .975), names = FALSE)
out <- list()
for (v in names(versions)) {
  m <- fit_age(versions[[v]]); b <- coef(m); s2 <- m$sigma2
  dr <- mvrnorm(20000, b, vcov(m))
  typical <- plogis(dr) %*% mix                                        # efeitos de estudo = 0
  newpop  <- plogis(dr + rnorm(20000, 0, sqrt(s2[1])) +                 # nova população: sorteia o efeito do
                    matrix(rnorm(20000 * 6, 0, sqrt(s2[2])), ncol = 6)) %*% mix   # estudo e o de faixa
  out[[v]] <- data.frame(version = v, estimand = c("typical population", "new population (prediction)"),
                         rbind(q(typical), q(newpop)))
}
out <- bind_rows(out); names(out)[3:5] <- c("expected", "lb", "ub")
ci <- t(mapply(function(x, n) binom.test(x, n)$conf.int, obs$x, obs$n))
obs <- cbind(obs, observed = obs$x / obs$n, obs_lb = ci[, 1], obs_ub = ci[, 2])
res <- merge(out, obs[1, ], by = NULL)
res$observed_inside_interval <- res$observed >= res$lb & res$observed <= res$ub
res$ratio_observed_to_expected <- res$observed / res$expected
print(res, digits = 3); print(obs, digits = 3)
write.csv(res, file.path(OUT_TAB, "pieta_out_of_sample.csv"), row.names = FALSE)
