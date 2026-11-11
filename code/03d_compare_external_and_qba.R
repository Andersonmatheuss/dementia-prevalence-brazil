# Produces: output/tables/external_comparison_elsi_by_age.csv, triangulation_national_estimates.csv, qba_two_phase_screen_sensitivity.csv, qba_elsi_algorithm_misclassification.csv
# Comparação externa com o ELSI-Brazil e análise quantitativa de viés de verificação.
#
# (1) Fora da amostra: modelo por idade SEM o Bertola/ELSI, comparado com o ELSI ponderado por faixa.
#     Só agregados publicados (Bertola 2023, Tabela 2); nenhum microdado do ELSI é usado.
# (2) Triangulação das estimativas nacionais disponíveis, todas levadas à população do Censo 2022.
# (3) QBA por cenários — NENHUM cenário é apresentado como prevalência verdadeira:
#     (a) estudos em duas fases sem correção (Herrera, Bottino, Lopes, Correa Ribeiro): só os positivos no
#         rastreio passaram por avaliação clínica; prevalência corrigida = observada / Se_rastreio,
#         supondo especificidade clínica ~1 e seguimento completo da fase 2. Âncora empírica: os próprios
#         autores corrigiram Bottino (6,8 -> 12,9) e Lopes (5,9 -> 12,5), razão ~1,9-2,1.
#     (b) classificação algorítmica do ELSI: p_corr = (p_obs - (1 - Sp)) / (Se - (1 - Sp)).
suppressPackageStartupMessages({ library(metafor); library(dplyr); library(MASS, exclude = "select") })
source(file.path("code", "helpers", "setup.R"))
set.seed(20261019)

agec <- read.csv(file.path(DERIVED, "age_strata_corrected.csv"))
cen  <- read.csv(file.path(OCT_IBGE, "ibge_census_brazil_2022_agegroups.csv"))
bands <- c("60-64", "65-69", "70-74", "75-79", "80-84", "85+")
pop6 <- c(cen$population[match(bands[1:5], cen$age_group)], sum(cen$population[cen$age_group %in% c("85-89", "90+")]))
pop7 <- setNames(cen$population, cen$age_group)
fit_age <- function(d) { d$band <- factor(d$band, levels = bands); rma.mv(yi, vi, mods = ~ 0 + band, random = ~ 1 | study / row, data = d, method = "REML") }
burden <- function(m, n = 20000) { dr <- plogis(mvrnorm(n, coef(m), vcov(m))); cs <- as.vector(dr %*% pop6); quantile(cs, c(.5, .025, .975)) }

## (1) Fora da amostra -------------------------------------------------------------------------
m_wo <- fit_age(subset(agec, study != "Bertola_2023"))
els  <- subset(agec, study == "Bertola_2023")
els$band <- factor(els$band, levels = bands)
pr <- predict(m_wo, newmods = diag(6), transf = transf.ilogit)   # pi = nova população (inclui sigma2 de estudo + estrato)
oos <- data.frame(band = bands, pooled_without_elsi = pr$pred, ci_lb = pr$ci.lb, ci_ub = pr$ci.ub,
                  pi_lb = pr$pi.lb, pi_ub = pr$pi.ub, elsi_weighted = plogis(els$yi[match(bands, els$band)]))
oos$elsi_inside_pi <- oos$elsi_weighted >= oos$pi_lb & oos$elsi_weighted <= oos$pi_ub
oos$ratio_pooled_to_elsi <- oos$pooled_without_elsi / oos$elsi_weighted
z <- (coef(m_wo) - els$yi[match(bands, els$band)]) / sqrt(diag(vcov(m_wo)) + els$vi[match(bands, els$band)])
oos$z_logit_difference <- z; oos$p_difference <- 2 * pnorm(-abs(z))
write.csv(oos, file.path(OUT_TAB, "external_comparison_elsi_by_age.csv"), row.names = FALSE)
print(oos, digits = 3)

## (2) Triangulação (população 60+ do Censo 2022 = 32.113.490) -------------------------------
elsi7 <- c("60-64" = .032, "65-69" = .027, "70-74" = .039, "75-79" = .087, "80-84" = .126, "85-89" = .213, "90+" = .428)
tri <- data.frame(
  source = c("ELSI-Brazil 2015-16, weighted age-specific (Bertola 2023) x Census 2022",
             "ELSI-Brazil, authors' own national estimate (2019 population)",
             "Delphi consensus 8.5% (Aliberti 2025), 2019 population (as published)",
             "Delphi consensus 8.5% applied to Census 2022 population 60+",
             "Manuscript, age-specific frequentist (original, 2022)",
             "Manuscript, Bayesian spline (original, 2022)",
             "Audit: corrected pooled age schedule, typical-study estimand (2022)",
             "Audit: corrected pooled age schedule without ELSI (2022)"),
  cases = c(sum(elsi7 * pop7[names(elsi7)]), 1757480, 2460000, 0.085 * sum(pop7), 2267772, 2419715,
            burden(fit_age(agec))[1], burden(m_wo)[1]),
  verification = c("computed here from published weighted rates", "PubMed abstract PMID 36682021", "PubMed abstract PMID 39938925",
                   "computed here", "manuscript line 62", "manuscript line 62", "script 05/08", "script 08"))
write.csv(tri, file.path(OUT_TAB, "triangulation_national_estimates.csv"), row.names = FALSE)
print(tri, digits = 7)

## (3a) QBA: sensibilidade do rastreio nos estudos em duas fases --------------------------------
two_phase <- c("Herrera_2002", "Bottino_2008", "Lopes_2012", "Ribeiro_2013")
qba_a <- bind_rows(lapply(c(1, 0.9, 0.8, 0.7, 0.6, 0.5), function(se) {
  d <- agec; i <- d$study %in% two_phase & !is.na(d$cases)
  p <- pmin((d$cases[i] / d$n[i]) / se, 0.99)
  d$yi[i] <- qlogis(p); d$vi[i] <- 1 / (d$n[i] * p * (1 - p))   # variância do n observado (otimista)
  b <- burden(fit_age(d))
  data.frame(screen_sensitivity = se, cases_median = b[1], cases_lb = b[2], cases_ub = b[3])
}))
write.csv(qba_a, file.path(OUT_TAB, "qba_two_phase_screen_sensitivity.csv"), row.names = FALSE)
print(qba_a, digits = 7)

## (3b) QBA: classificação algorítmica do ELSI ----------------------------------------------------
grid <- expand.grid(Se = c(0.6, 0.7, 0.8, 0.9, 1.0), Sp = c(0.96, 0.97, 0.98, 0.99, 1.0))
qba_b <- bind_rows(lapply(seq_len(nrow(grid)), function(j) {
  se <- grid$Se[j]; sp <- grid$Sp[j]
  pc <- (elsi7 - (1 - sp)) / (se - (1 - sp))
  data.frame(Se = se, Sp = sp, any_negative_stratum = any(pc < 0),
             cases_2022 = if (any(pc < 0)) NA else sum(pmin(pc, 1) * pop7[names(elsi7)]),
             corrected_prev_60plus = if (any(pc < 0)) NA else sum(pmin(pc, 1) * pop7[names(elsi7)]) / sum(pop7))
}))
write.csv(qba_b, file.path(OUT_TAB, "qba_elsi_algorithm_misclassification.csv"), row.names = FALSE)
print(qba_b, digits = 4)
writeLines(capture.output(sessionInfo()), file.path(OUT_DIAG, "sessioninfo_08_compare_external_and_qba.txt"))
