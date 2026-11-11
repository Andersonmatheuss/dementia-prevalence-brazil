# Produces: output/tables/primary_age_strata.csv, primary_burden_estimands.csv, primary_external_comparison_elsi.csv, primary_voi.csv, primary_burden_sensitivity_sets.csv, primary_elsi_algorithms_65plus.csv, pba_screening_correction.csv, pba_study_corrections.csv, uncertainty_decomposition_2022.csv, age_scale_correction_strata.csv, age_scale_correction_burden.csv
# Análises secundárias no conjunto principal sem sobreposição (PRIMARY_POPULATIONS, code/helpers/setup.R), nas versões
# "observada" e "melhor estimativa dos autores" (03g_fit_author_best_estimates.R): IC por faixa, três estimandos da carga,
# comparação fora da amostra com o ELSI e VOI (redução esperada de variância).
suppressPackageStartupMessages({ library(metafor); library(dplyr); library(MASS, exclude = "select") })
source(file.path("code", "helpers", "setup.R"))
set.seed(20261019)
agec_all <- read.csv(file.path(DERIVED, "age_strata_corrected_with_epifloripa.csv"))   # as oito populações
agec7 <- subset(agec_all, study %in% PRIMARY_POPULATIONS)   # análise principal sem sobreposição (code/helpers/setup.R)
best  <- read.csv(file.path(DERIVED, "author_best_estimates_s6.csv"))
cen <- read.csv(file.path(OCT_IBGE, "ibge_census_brazil_2022_agegroups.csv"))
bands <- c("60-64", "65-69", "70-74", "75-79", "80-84", "85+")
pop6 <- c(cen$population[match(bands[1:5], cen$age_group)], sum(cen$population[cen$age_group %in% c("85-89", "90+")]))
fit_age <- function(d) { d$band <- factor(d$band, levels = bands); rma.mv(yi, vi, mods = ~ 0 + band, random = ~ 1 | study / row, data = d, method = "REML") }

author_best_rows <- function(d) {   # mesma regra de 03g (efeito de desenho 1)
  for (s in c("Bottino_2008", "Lopes_2012", "EpiFloripa_2024")) {
    f <- best$factor[best$study == s]; i <- d$study == s & !is.na(d$cases)
    pp <- pmin(pmax(d$cases[i], 0.5) / d$n[i] * f, 0.95); d$yi[i] <- qlogis(pp); d$vi[i] <- 1 / (d$n[i] * pp * (1 - pp))
  }
  d
}
versions <- list(observed = agec7, author_best = author_best_rows(agec7))
# Mesmos sorteios da curva principal que 03g grava e que 03m, 03n e 03p usam, para que a carga de 2022
# seja um único número em todas as tabelas (o modelo de 03g com efeito de desenho 1 é o mesmo de author_best_rows).
common_dr <- qlogis(as.matrix(read.csv(file.path(DERIVED, "primary_age_schedule_draws.csv"), check.names = FALSE))[, bands])
drs <- list()

q <- function(x) quantile(x, c(.5, .025, .975))
strata <- list(); est <- list(); oos <- list(); voi <- list()
for (v in names(versions)) {
  d <- versions[[v]]; m <- fit_age(d); b <- coef(m); V <- vcov(m); s2 <- m$sigma2
  strata[[v]] <- data.frame(version = v, band = bands, prev = plogis(b), lb = plogis(b - 1.96 * sqrt(diag(V))), ub = plogis(b + 1.96 * sqrt(diag(V))))
  dr <- mvrnorm(20000, b, V)   # sorteado sempre, para não mudar a sequência aleatória dos passos seguintes
  if (v == "author_best") dr <- common_dr
  drs[[v]] <- dr
  u <- rnorm(20000, 0, sqrt(s2[1])); w <- matrix(rnorm(20000 * 6, 0, sqrt(s2[2])), ncol = 6)
  z <- matrix(rnorm(400 * 2), ncol = 2)
  pa <- sapply(1:6, function(j) sapply(1:2000, function(i) mean(plogis(dr[i, j] + z[, 1] * sqrt(s2[1]) + z[, 2] * sqrt(s2[2])))))
  est[[v]] <- data.frame(version = v, estimand = c("typical study", "population-averaged", "new population (prediction)"),
                         rbind(q(plogis(dr) %*% pop6), q(pa %*% pop6), q(plogis(dr + u + w) %*% pop6)), check.names = FALSE)
  # fora da amostra: o ELSI não está no conjunto principal; suas faixas vêm do arquivo das oito populações
  mw <- m; el <- subset(agec_all, study == "Bertola_2023")
  pr <- predict(mw, newmods = diag(6), transf = transf.ilogit); ey <- el$yi[match(bands, el$band)]
  oos[[v]] <- data.frame(version = v, band = bands, pooled_without_elsi = pr$pred, pi_lb = pr$pi.lb, pi_ub = pr$pi.ub,
                         elsi = plogis(ey), inside_pi = plogis(ey) >= pr$pi.lb & plogis(ey) <= pr$pi.ub,
                         ratio = pr$pred / plogis(ey), p_diff = 2 * pnorm(-abs((coef(mw) - ey) / sqrt(diag(vcov(mw)) + el$vi[match(bands, el$band)]))))
  # VOI (delta)
  se_b <- function(mm) { p <- plogis(coef(mm)); g <- pop6 * p * (1 - p); sqrt(as.numeric(t(g) %*% vcov(mm) %*% g)) }
  se0 <- se_b(m); share <- pop6 / sum(pop6)
  for (k in c(1, 3, 5)) for (n in c(1000, 5000)) {
    if (k > 1 && n == 5000) next
    sims <- replicate(200, {
      add <- bind_rows(lapply(seq_len(k), function(j) {
        nn <- pmax(round(n * share), 5); x <- rbinom(6, nn, plogis(b + rnorm(1, 0, sqrt(s2[1])) + rnorm(6, 0, sqrt(s2[2]))))
        es <- escalc("PLO", xi = x, ni = nn, add = 0.5, to = "only0"); data.frame(study = paste0("new", j), band = bands, yi = es$yi, vi = es$vi)
      }))
      dd <- bind_rows(d[, c("study", "band", "yi", "vi")], add); dd$row <- seq_len(nrow(dd)); tryCatch(se_b(fit_age(dd)), error = function(e) NA)
    })
    voi[[length(voi) + 1]] <- data.frame(version = v, new_surveys = k, n_each = n, current_se = se0, expected_se = mean(sims, na.rm = TRUE),
                                         variance_reduction_pct = 100 * (1 - mean(sims^2, na.rm = TRUE) / se0^2))
  }
}
strata <- bind_rows(strata); est <- bind_rows(est); names(est)[3:5] <- c("cases_median", "cases_lb", "cases_ub")
oos <- bind_rows(oos); voi <- bind_rows(voi)
# Sensibilidades fixadas com a regra de sobreposição (SENSITIVITY_SETS): carga 2022, estimando típico,
# nas duas versões. Semente própria, depois das tabelas acima, para não mudar os sorteios delas.
set.seed(20261019)
versions_all <- list(observed = agec_all, author_best = author_best_rows(agec_all))
sens <- bind_rows(lapply(names(SENSITIVITY_SETS), function(s) bind_rows(lapply(names(versions_all), function(v) {
  m <- fit_age(subset(versions_all[[v]], study %in% SENSITIVITY_SETS[[s]])); dr <- mvrnorm(20000, coef(m), vcov(m))
  data.frame(set = s, version = v, k = m$s.nlevels[1], t(q(plogis(dr) %*% pop6)), check.names = FALSE)
}))))
names(sens)[4:6] <- c("cases_median", "cases_lb", "cases_ub")
# Mesmo inquérito, outro algoritmo. Aliberti et al. 2026 (Alzheimers Dement 22:e71430, PMID 42056069,
# texto integral PMC13128339) classificaram o ELSI-Brazil 2015-16, 65+, com um algoritmo harmonizado com o HRS:
# 8,2% ponderada pelo desenho amostral; 8,0% (IC 95% 6,9-9,0) padronizada por idade e sexo. Para comparar na mesma
# faixa etária, as taxas ponderadas de Bertola 2023 e a curva principal vão para a população de 65+ do Censo 2022.
# Semente própria, depois das tabelas acima, para não mudar os sorteios delas.
set.seed(20261020)
w65 <- pop6[-1] / sum(pop6[-1])
el <- subset(agec_all, study == "Bertola_2023"); el <- el[match(bands, el$band), ]
de <- sapply(2:6, function(j) rnorm(20000, el$yi[j], sqrt(el$vi[j])))
p65 <- list(data.frame(source = "ELSI-Brazil, Bertola 2023 algorithm (weighted age-specific rates)", t(q(plogis(de) %*% w65)), check.names = FALSE))
for (v in names(versions)) {
  dr <- drs[[v]]   # os mesmos sorteios da carga de 2022
  p65[[v]] <- data.frame(source = paste("Primary age schedule,", c(observed = "observed counts", author_best = "authors' estimates")[[v]]),
                         t(q(plogis(dr[, -1]) %*% w65)), check.names = FALSE)
}
p65 <- bind_rows(p65); names(p65)[2:4] <- c("prev", "lb", "ub"); p65$age_structure <- "Census 2022, 65+"
p65 <- bind_rows(p65, data.frame(source = rep("ELSI-Brazil, harmonised algorithm (Aliberti 2026), as published", 2),
                                 prev = c(0.082, 0.080), lb = c(NA, 0.069), ub = c(NA, 0.090),
                                 age_structure = c("ELSI-Brazil 2015-16, survey-weighted", "age- and sex-standardised (marginal standardisation)")))
print(strata, digits = 3); print(est, digits = 4); print(oos, digits = 3); print(voi, digits = 3); print(sens, digits = 4); print(p65, digits = 3)
write.csv(strata, file.path(OUT_TAB, "primary_age_strata.csv"), row.names = FALSE)
write.csv(est, file.path(OUT_TAB, "primary_burden_estimands.csv"), row.names = FALSE)
write.csv(oos, file.path(OUT_TAB, "primary_external_comparison_elsi.csv"), row.names = FALSE)
write.csv(voi, file.path(OUT_TAB, "primary_voi.csv"), row.names = FALSE)
write.csv(sens, file.path(OUT_TAB, "primary_burden_sensitivity_sets.csv"), row.names = FALSE)
write.csv(p65, file.path(OUT_TAB, "primary_elsi_algorithms_65plus.csv"), row.names = FALSE)
# --- Incerteza da correção do rastreio no conjunto principal ---
# Não substitui a análise principal (estimativa de cada autor, 03g e as tabelas acima). Especificação fixada antes de rodar:
# (1) Análise probabilística de viés (PBA, erro total) dos inquéritos em duas fases do conjunto principal (Herrera,
#     Lopes, EpiFloripa). Em cada iteração sorteia, por estudo, s+ (positivos no rastreio), VPP (casos entre os
#     positivos examinados) e q = 1 - VPN (casos entre os negativos), com Beta de Jeffreys nos dois primeiros;
#     corrige o estudo por p = s+ x VPP + (1 - s+) x q, leva o fator às faixas como em 03g, reajusta a curva
#     etária e sorteia o erro aleatório dos coeficientes. Os estudos de uma fase entram como na análise principal.
#     q: agrupamento parcial (Bayes empírico) dos quatro inquéritos que examinaram negativos (Bottino 6/157, Lopes
#     1/30, Correa Ribeiro 0/44, EpiFloripa 0/16; modelo logit de efeitos aleatórios, rma.glmm); o Herrera, sem
#     negativos examinados, recebe a distribuição preditiva. Variante sem agrupamento: Jeffreys pelos negativos de
#     cada estudo; Herrera uniforme em [0; 0,05] (a faixa de VPN dos cenários de 03f_add_epifloripa.R).
# (2) Decomposição aproximada da variância da carga 2022: verificação (parâmetros do rastreio), variação entre
#     populações e amostragem (modelo com e sem os componentes de variância, método delta).
# (3) Escala da correção por idade (Lopes e EpiFloripa, os estudos do conjunto principal com fator do autor):
#     razão constante (análise principal), deslocamento constante no logito (odds) ou diferença absoluta constante,
#     com o mesmo total corrigido em cada estudo.
# Bloco no fim do script, depois de gravadas as tabelas acima, para não mudar nenhuma delas; roda antes
# da Figura 2 (05a_plot_national_estimates.py), que lê a PBA. Semente própria.
set.seed(20261020)
fit_age2 <- function(d, ...) { d$band <- factor(d$band, levels = bands); rma.mv(yi, vi, mods = ~ 0 + band, random = ~ 1 | study / row, data = d, method = "REML", ...) }
q3 <- function(x) quantile(x, c(.5, .025, .975), names = FALSE)
apply_factor <- function(d, s, f) {   # mesma regra de 03g: razão nas faixas, teto 0,95
  j <- d$study == s & !is.na(d$cases); p <- pmin(pmax(d$cases[j], 0.5) / d$n[j] * f, 0.95)
  d$yi[j] <- qlogis(p); d$vi[j] <- 1 / (d$n[j] * p * (1 - p)); d
}
d0 <- subset(agec_all, study %in% PRIMARY_POPULATIONS)

# --- quantidades do rastreio conferidas (mesma regra e mesmos n de primeira fase de 03f) ---
vest <- read.csv(VERIFIED_ESTIMATES)
two <- c("Herrera_2002", "Bottino_2008", "Lopes_2012", "Ribeiro_2013", "EpiFloripa_2024")
e <- vest[match(two, vest$study), ]
ph <- data.frame(study = two, n = c(1656, 1563, 1145, 736, 1184), pos = e$screen_positive, evaluated = e$positives_examined,
                 cases_pos = e$cases_in_positives, neg_exam = e$negatives_examined,
                 cases_neg = ifelse(is.na(e$cases_in_negatives), 0, e$cases_in_negatives), observed = e$cases_total / e$n_total)
stopifnot(ph$neg_exam[ph$study == "Herrera_2002"] == 0, all(ph$neg_exam[ph$study != "Herrera_2002"] > 0))
two_primary <- intersect(ph$study, PRIMARY_POPULATIONS)

# --- (1) PBA ---
N <- 2000
ng <- subset(ph, neg_exam > 0)
gm <- rma.glmm(measure = "PLO", xi = cases_neg, ni = neg_exam, data = ng, method = "ML")
mu <- as.numeric(coef(gm)); se_mu <- as.numeric(gm$se); tau <- sqrt(gm$tau2)
draw_q <- list(
  partial_pooling = function(r) {
    th <- rnorm(50 * N, rnorm(50 * N, mu, se_mu), tau)   # preditiva no logito, com a incerteza de mu
    if (r$neg_exam == 0) return(plogis(th[seq_len(N)]))
    plogis(sample(th, N, replace = TRUE, prob = dbinom(r$cases_neg, r$neg_exam, plogis(th))))
  },
  no_pooling = function(r) if (r$neg_exam == 0) runif(N, 0, 0.05) else rbeta(N, r$cases_neg + 0.5, r$neg_exam - r$cases_neg + 0.5)
)
pba <- list(); pars_all <- list()
for (v in names(draw_q)) {
  pars <- bind_rows(lapply(two_primary, function(s) {
    r <- ph[ph$study == s, ]
    sp <- rbeta(N, r$pos + 0.5, r$n - r$pos + 0.5); pv <- rbeta(N, r$cases_pos + 0.5, r$evaluated - r$cases_pos + 0.5); qq <- draw_q[[v]](r)
    data.frame(variant = v, study = s, iter = seq_len(N), s_pos = sp, ppv = pv, q = qq, p_corr = sp * pv + (1 - sp) * qq)
  }))
  pars$factor <- pars$p_corr / ph$observed[match(pars$study, ph$study)]
  res <- t(sapply(seq_len(N), function(i) {
    d <- d0
    for (s in two_primary) d <- apply_factor(d, s, pars$factor[pars$study == s & pars$iter == i])
    m <- fit_age2(d); b <- coef(m)
    c(bias_only = sum(plogis(b) * pop6), total = sum(plogis(mvrnorm(1, b, vcov(m))) * pop6))
  }))
  # erro aleatório só, com os parâmetros do rastreio na média
  dm <- d0
  for (s in two_primary) dm <- apply_factor(dm, s, mean(pars$factor[pars$study == s]))
  mc <- fit_age2(dm); conv <- as.vector(plogis(mvrnorm(20000, coef(mc), vcov(mc))) %*% pop6)
  pba[[v]] <- data.frame(variant = v, quantity = c("total error (verification and random)", "verification parameters only",
                                                   "random error only (verification parameters at their means)"),
                         rbind(q3(res[, "total"]), q3(res[, "bias_only"]), q3(conv)),
                         share_of_variance_from_verification = c(var(res[, "bias_only"]) / var(res[, "total"]), NA, NA))
  pars_all[[v]] <- pars
  if (v == "partial_pooling") { res_main <- res; mc_main <- mc; dm_main <- dm }
}
pba <- bind_rows(pba); names(pba)[3:5] <- c("cases_median", "cases_lb", "cases_ub")
pars_all <- bind_rows(pars_all)
study_pba <- pars_all %>% group_by(variant, study) %>%
  summarise(across(c(s_pos, ppv, q, p_corr), list(median = median, lb = ~ quantile(.x, .025), ub = ~ quantile(.x, .975))), .groups = "drop")

# --- (2) decomposição da variância (variante principal da PBA) ---
g <- function(m) { p <- plogis(coef(m)); pop6 * p * (1 - p) }
v_re <- as.numeric(t(g(mc_main)) %*% vcov(mc_main) %*% g(mc_main))
mfe <- fit_age2(dm_main, sigma2 = c(0, 0)); v_fe <- as.numeric(t(g(mfe)) %*% vcov(mfe) %*% g(mfe))
sh_ver <- var(res_main[, "bias_only"]) / var(res_main[, "total"])
decomp <- data.frame(component = c("verification (screening parameters)", "between-population variation", "sampling"),
                     share_of_variance = c(sh_ver, (1 - sh_ver) * (1 - v_fe / v_re), (1 - sh_ver) * v_fe / v_re))

# --- (3) escala da correção por idade (estimativas dos autores; mesmo total por estudo) ---
scaled <- function(scale) {
  d <- d0
  for (s in c("Lopes_2012", "EpiFloripa_2024")) {
    f <- best$factor[best$study == s]; j <- d$study == s & !is.na(d$cases)
    nb <- d$n[j]; pb <- pmax(d$cases[j], 0.5) / nb; target <- sum(nb * pb * f) / sum(nb)
    p <- switch(scale,
                ratio = pmin(pb * f, 0.95),
                odds = { dl <- uniroot(function(x) sum(nb * plogis(qlogis(pb) + x)) / sum(nb) - target, c(-5, 5), tol = 1e-10)$root; plogis(qlogis(pb) + dl) },
                difference = pb + (target - sum(nb * pb) / sum(nb)))
    d$yi[j] <- qlogis(p); d$vi[j] <- 1 / (nb * p * (1 - p))
  }
  m <- fit_age2(d); dr <- mvrnorm(20000, coef(m), vcov(m))
  if (scale == "ratio") dr <- common_dr   # razão constante = análise principal: mesmos sorteios
  cs <- as.vector(plogis(dr) %*% pop6)
  p65 <- as.vector(plogis(dr[, -1]) %*% (pop6[-1] / sum(pop6[-1])))
  list(strata = data.frame(scale = scale, band = bands, prev = plogis(coef(m)), lb = plogis(coef(m) - 1.96 * sqrt(diag(vcov(m)))),
                           ub = plogis(coef(m) + 1.96 * sqrt(diag(vcov(m))))),
       burden = data.frame(scale = scale, t(q3(cs)), t(q3(p65))))
}
sc <- lapply(c("ratio", "odds", "difference"), scaled)
age_scale <- bind_rows(lapply(sc, `[[`, "strata"))
el <- subset(agec_all, study == "Bertola_2023"); age_scale$elsi <- plogis(el$yi[match(age_scale$band, el$band)])
age_scale$ratio_to_elsi <- age_scale$prev / age_scale$elsi
age_burden <- bind_rows(lapply(sc, `[[`, "burden"))
names(age_burden)[2:7] <- c("cases_median", "cases_lb", "cases_ub", "prev65_median", "prev65_lb", "prev65_ub")

print(pba, digits = 4); print(as.data.frame(study_pba), digits = 3); print(decomp, digits = 3)
print(age_scale, digits = 3); print(age_burden, digits = 4)
cat(sprintf("q agrupado: mu = %.3f (logito), tau2 = %.4f, q mediano = %.4f\n", mu, gm$tau2, plogis(mu)))

# a escala de razão tem de reproduzir a análise principal (primary_burden_estimands.csv), dentro do erro de Monte Carlo
ref <- est$cases_median[est$version == "author_best" & est$estimand == "typical study"]
stopifnot("Escala de razão = análise principal" = abs(age_burden$cases_median[age_burden$scale == "ratio"] / ref - 1) < 0.01,
          "PBA sem NA" = !anyNA(pba$cases_median), "PBA com N iterações" = nrow(res_main) == N)

write.csv(pba, file.path(OUT_TAB, "pba_screening_correction.csv"), row.names = FALSE)
write.csv(study_pba, file.path(OUT_TAB, "pba_study_corrections.csv"), row.names = FALSE)
write.csv(decomp, file.path(OUT_TAB, "uncertainty_decomposition_2022.csv"), row.names = FALSE)
write.csv(age_scale, file.path(OUT_TAB, "age_scale_correction_strata.csv"), row.names = FALSE)
write.csv(age_burden, file.path(OUT_TAB, "age_scale_correction_burden.csv"), row.names = FALSE)

writeLines(capture.output(sessionInfo()), file.path(OUT_DIAG, "sessioninfo_14_fit_primary_set.txt"))
