# Produces: output/tables/pba_ppv_unexamined.csv
# Análise probabilística de viés (PBA) dos inquéritos em duas fases do conjunto principal (Herrera, Lopes,
# EpiFloripa) com um parâmetro a mais: a razão r entre o VPP dos positivos no rastreio NÃO examinados e o VPP
# dos examinados. A PBA de 03j supõe r = 1. Aqui:
#   p = s+ x [e x VPP + (1 - e) x min(r x VPP, 1)] + (1 - s+) x q,  e = examinados / positivos,
# com s+, VPP e q sorteados como em 03j (Jeffreys; q por agrupamento parcial dos quatro inquéritos que
# examinaram negativos). Cenários: r fixo em 0,5, 0,75, 1, 1,25 e 1,5, e r ~ Uniforme(0,5; 1,5), comum aos
# três inquéritos em cada iteração. Mesma regra de 03g/03j para levar a correção às faixas etárias.
suppressPackageStartupMessages({ library(metafor); library(dplyr); library(MASS, exclude = "select") })
source(file.path("code", "helpers", "setup.R"))
set.seed(20261021)
agec_all <- read.csv(file.path(DERIVED, "age_strata_corrected_with_epifloripa.csv"))
cen <- read.csv(file.path(OCT_IBGE, "ibge_census_brazil_2022_agegroups.csv"))
bands <- c("60-64", "65-69", "70-74", "75-79", "80-84", "85+")
pop6 <- c(cen$population[match(bands[1:5], cen$age_group)], sum(cen$population[cen$age_group %in% c("85-89", "90+")]))
fit_age2 <- function(d) { d$band <- factor(d$band, levels = bands); rma.mv(yi, vi, mods = ~ 0 + band, random = ~ 1 | study / row, data = d, method = "REML") }
q3 <- function(x) quantile(x, c(.5, .025, .975), names = FALSE)
apply_factor <- function(d, s, f) {   # mesma regra de 03g/03j: razão nas faixas, teto 0,95
  j <- d$study == s & !is.na(d$cases); p <- pmin(pmax(d$cases[j], 0.5) / d$n[j] * f, 0.95)
  d$yi[j] <- qlogis(p); d$vi[j] <- 1 / (d$n[j] * p * (1 - p)); d
}
d0 <- subset(agec_all, study %in% PRIMARY_POPULATIONS)
vest <- read.csv(VERIFIED_ESTIMATES)
two <- c("Herrera_2002", "Bottino_2008", "Lopes_2012", "Ribeiro_2013", "EpiFloripa_2024")
e <- vest[match(two, vest$study), ]
ph <- data.frame(study = two, n = c(1656, 1563, 1145, 736, 1184), pos = e$screen_positive, evaluated = e$positives_examined,
                 cases_pos = e$cases_in_positives, neg_exam = e$negatives_examined,
                 cases_neg = ifelse(is.na(e$cases_in_negatives), 0, e$cases_in_negatives), observed = e$cases_total / e$n_total)
two_primary <- intersect(ph$study, PRIMARY_POPULATIONS)
N <- 2000
ng <- subset(ph, neg_exam > 0)
gm <- rma.glmm(measure = "PLO", xi = cases_neg, ni = neg_exam, data = ng, method = "ML")
mu <- as.numeric(coef(gm)); se_mu <- as.numeric(gm$se); tau <- sqrt(gm$tau2)
draw_q <- function(r) {   # agrupamento parcial, como em 03j
  th <- rnorm(50 * N, rnorm(50 * N, mu, se_mu), tau)
  if (r$neg_exam == 0) return(plogis(th[seq_len(N)]))
  plogis(sample(th, N, replace = TRUE, prob = dbinom(r$cases_neg, r$neg_exam, plogis(th))))
}
scen <- list("r = 0.50" = function() rep(0.5, N), "r = 0.75" = function() rep(0.75, N), "r = 1.00" = function() rep(1, N),
             "r = 1.25" = function() rep(1.25, N), "r = 1.50" = function() rep(1.5, N), "r ~ Uniform(0.5, 1.5)" = function() runif(N, 0.5, 1.5))
out <- list()
for (sc in names(scen)) {
  rr <- scen[[sc]]()
  pars <- bind_rows(lapply(two_primary, function(s) {
    r <- ph[ph$study == s, ]; ex <- r$evaluated / r$pos
    sp <- rbeta(N, r$pos + 0.5, r$n - r$pos + 0.5); pv <- rbeta(N, r$cases_pos + 0.5, r$evaluated - r$cases_pos + 0.5); qq <- draw_q(r)
    data.frame(study = s, iter = seq_len(N), p_corr = sp * (ex * pv + (1 - ex) * pmin(rr * pv, 1)) + (1 - sp) * qq)
  }))
  pars$factor <- pars$p_corr / ph$observed[match(pars$study, ph$study)]
  res <- t(sapply(seq_len(N), function(i) {
    d <- d0
    for (s in two_primary) d <- apply_factor(d, s, pars$factor[pars$study == s & pars$iter == i])
    m <- fit_age2(d); b <- coef(m)
    c(bias_only = sum(plogis(b) * pop6), total = sum(plogis(mvrnorm(1, b, vcov(m))) * pop6))
  }))
  sp <- pars %>% group_by(study) %>% summarise(p = median(p_corr), .groups = "drop")
  out[[sc]] <- data.frame(scenario = sc, t(q3(res[, "total"])), share_of_variance_from_screening = var(res[, "bias_only"]) / var(res[, "total"]),
                          corrected_catanduva = sp$p[sp$study == "Herrera_2002"], corrected_ribeirao = sp$p[sp$study == "Lopes_2012"],
                          corrected_florianopolis = sp$p[sp$study == "EpiFloripa_2024"])
}
out <- bind_rows(out); names(out)[2:4] <- c("cases_median", "cases_lb", "cases_ub")
print(out, digits = 4)
# r = 1 é a PBA de 03j (variante de agrupamento parcial), dentro do erro de Monte Carlo
ref <- read.csv(file.path(OUT_TAB, "pba_screening_correction.csv"))
r1 <- ref$cases_median[ref$variant == "partial_pooling" & grepl("^total", ref$quantity)]
stopifnot(abs(out$cases_median[out$scenario == "r = 1.00"] / r1 - 1) < 0.03)
write.csv(out, file.path(OUT_TAB, "pba_ppv_unexamined.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(OUT_DIAG, "sessioninfo_03s_pba_ppv_unexamined.txt"))
