# Produces: output/intermediate/primary_age_schedule_draws.csv, output/tables/author_best_pooled.csv, author_best_age_strata.csv, author_best_burden_projections.csv, author_best_projections_regions.csv, author_best_sensitivity_sets.csv; output/intermediate/author_best_estimates_s6.csv
# Agrupamento e carga no conjunto principal sem sobreposição (PRIMARY_POPULATIONS), com as
# sensibilidades de SENSITIVITY_SETS; `best` continua com as oito populações (conjunto S6).
#
# Regra: usar o estimando de prevalência populacional que cada estudo publicou —
#   ELSI: ponderada (5,8%; 4,7-7,2);  Bottino: corrigida 12,9%;  Lopes: corrigida 12,5%;
#   EpiFloripa: corrigida 9,2% (7,45-10,95);  demais (Herrera, César, Scazufca, Correa Ribeiro): observada,
#   porque não publicaram correção (Herrera não avaliou negativos; César e Scazufca são de fase única).
# Variância: IC publicado quando existe; para Bottino e Lopes (sem IC no resumo), binomial no n observado,
# o que SUBESTIMA a incerteza das estimativas corrigidas — inflada por efeito de desenho 2 como sensibilidade.
# Idade: as contagens por faixa de Bottino, Lopes e EpiFloripa são escaladas pelo fator global (corrigida/observada).
suppressPackageStartupMessages({ library(metafor); library(dplyr); library(MASS, exclude = "select") })
source(file.path("code", "helpers", "setup.R"))
set.seed(20261019)
agec7 <- read.csv(file.path(DERIVED, "age_strata_corrected_with_epifloripa.csv"))
cen <- read.csv(file.path(OCT_IBGE, "ibge_census_brazil_2022_agegroups.csv"))
p24 <- read.csv(file.path(DERIVED, "ibge_projection_2024_60plus.csv"))
bands <- c("60-64", "65-69", "70-74", "75-79", "80-84", "85+")
pop6 <- c(cen$population[match(bands[1:5], cen$age_group)], sum(cen$population[cen$age_group %in% c("85-89", "90+")]))
fit_age <- function(d) { d$band <- factor(d$band, levels = bands); rma.mv(yi, vi, mods = ~ 0 + band, random = ~ 1 | study / row, data = d, method = "REML") }
vci <- function(lb, ub) ((qlogis(ub) - qlogis(lb)) / 3.92)^2

# Estimativa populacional publicada por cada autor, da extração conferida (VERIFIED_ESTIMATES):
# ponderada (ELSI 5,8 [4,7-7,2]; Correa Ribeiro 16,9 [14,4-19,8]) ou corrigida pela acurácia do rastreio
# (Bottino 12,9, sem IC publicado; Lopes 12,5 [10,6-14,4]; EpiFloripa 9,2 [7,45-10,95]); nos demais,
# casos/avaliados. O fator que leva as faixas etárias à estimativa do autor usa o total da MESMA tabela
# das faixas (Bottino 107/1.477; Lopes 67/1.052; EpiFloripa 47/1.043): com o bruto geral (107/1.563;
# 68/1.145) a correção ficava 6-7% maior que a publicada. No Correa Ribeiro o fator (16,9/16,8 = 1,004)
# é desprezível e as faixas ficam com as contagens.
est <- read.csv(VERIFIED_ESTIMATES)
ids <- c("Bertola_2023", "Bottino_2008", "Lopes_2012", "EpiFloripa_2024", "Herrera_2002", "Cesar_2015", "Scazufca_2008", "Ribeiro_2013")
e <- est[match(ids, est$study), ]
has_author <- !is.na(e$author_estimate_pct)
best <- data.frame(
  study = ids,
  prev  = ifelse(has_author, e$author_estimate_pct / 100, e$cases_total / e$n_total),
  n     = e$n_total,
  lb    = ifelse(has_author, e$author_lb / 100, NA), ub = ifelse(has_author, e$author_ub / 100, NA),
  observed = ifelse(ids %in% c("Bottino_2008", "Lopes_2012", "EpiFloripa_2024"), e$age_table_cases / e$age_table_n, NA))
best$factor <- ifelse(is.na(best$observed), 1, best$prev / best$observed)
best$primary <- best$study %in% PRIMARY_POPULATIONS   # análise principal sem sobreposição (code/helpers/setup.R)
run <- function(deff, ids = PRIMARY_POPULATIONS) {
  keep <- best$study %in% ids; b <- best[keep, ]; ha <- has_author[keep]; b$yi <- qlogis(b$prev)
  # IC publicado quando existe; sem IC (Bottino 12,9) variância binomial com efeito de desenho; contagens: binomial
  b$vi <- ifelse(!is.na(b$lb), vci(b$lb, b$ub), deff / (b$n * b$prev * (1 - b$prev)))
  b$vi[!ha] <- 1 / (b$n[!ha] * b$prev[!ha] * (1 - b$prev[!ha]))
  m <- rma(yi, vi, data = b, method = "REML", test = "knha"); p <- predict(m, transf = transf.ilogit)
  d <- subset(agec7, study %in% ids)
  for (s in intersect(c("Bottino_2008", "Lopes_2012", "EpiFloripa_2024"), ids)) {
    f <- best$factor[best$study == s]; i <- d$study == s & !is.na(d$cases)
    pp <- pmin(pmax(d$cases[i], 0.5) / d$n[i] * f, 0.95); d$yi[i] <- qlogis(pp); d$vi[i] <- deff / (d$n[i] * pp * (1 - pp))
  }
  ma <- fit_age(d); dr <- plogis(mvrnorm(20000, coef(ma), vcov(ma)))
  # Sorteios da curva etária principal (efeito de desenho 1, conjunto principal), reaproveitados por
  # 03m (UFs), 03n (cenários de tendência) e 03p (escolaridade) para que somem com a carga nacional.
  if (deff == 1 && identical(ids, PRIMARY_POPULATIONS)) {
    colnames(dr) <- bands; write.csv(dr, file.path(DERIVED, "primary_age_schedule_draws.csv"), row.names = FALSE) }
  tot <- function(popv) quantile(as.vector(dr %*% popv), c(.5, .025, .975))
  to6 <- function(yr, loc = "Brasil") { x <- subset(p24, local == loc & year == yr); c(x$population[match(bands[1:5], x$age)], sum(x$population[x$age %in% c("85-89", "90+")])) }
  yrs <- c(2022L, PROJECTION_YEARS)   # 2022 = Censo; demais anos = revisão 2024 do IBGE
  regions <- bind_rows(lapply(c("Norte", "Nordeste", "Sudeste", "Sul", "Centro-Oeste"), function(loc) bind_rows(lapply(yrs, function(y) {
    pv <- to6(y, loc); data.frame(deff = deff, local = loc, year = y, pop60 = sum(pv), t(tot(pv)), check.names = FALSE) }))))
  list(pooled = data.frame(deff = deff, k = m$k, prev = p$pred, hk_lb = p$ci.lb, hk_ub = p$ci.ub, pi_lb = p$pi.lb, pi_ub = p$pi.ub, tau2 = m$tau2, I2 = m$I2),
       strata = data.frame(deff = deff, band = bands, prev = plogis(coef(ma))),
       burden = data.frame(deff = deff, year = yrs, population = c(sum(pop6), sapply(PROJECTION_YEARS, function(y) sum(to6(y)))),
                           rbind(tot(pop6), t(sapply(PROJECTION_YEARS, function(y) tot(to6(y))))), check.names = FALSE),
       regions = regions)
}
r1 <- run(1); r2 <- run(2)
pooled <- rbind(r1$pooled, r2$pooled); strata <- rbind(r1$strata, r2$strata); burden <- rbind(r1$burden, r2$burden)
names(burden)[4:6] <- c("cases_median", "cases_lb", "cases_ub")
# Decomposição simétrica do aumento desde 2022 em tamanho da população 60+ e estrutura etária (medianas)
dec_size <- function(C0, P0, C1, P1) (P1 - P0) * (C0 / P0 + C1 / P1) / 2
burden <- burden |> group_by(deff) |> mutate(
  growth_size_component = ifelse(year == 2022, NA, dec_size(cases_median[year == 2022], population[year == 2022], cases_median, population)),
  growth_age_structure_component = ifelse(year == 2022, NA, cases_median - cases_median[year == 2022] - growth_size_component)) |> ungroup() |> as.data.frame()
regions <- r1$regions; names(regions)[5:7] <- c("cases_median", "cases_lb", "cases_ub")
# Sensibilidades fixadas com a regra de sobreposição (SENSITIVITY_SETS), efeito de desenho 1
sens <- bind_rows(lapply(names(SENSITIVITY_SETS), function(s) {
  r <- run(1, SENSITIVITY_SETS[[s]]); b22 <- r$burden[r$burden$year == 2022, ]
  data.frame(set = s, r$pooled[, c("k", "prev", "hk_lb", "hk_ub", "pi_lb", "pi_ub", "tau2")],
             cases_2022_median = b22[["50%"]], cases_2022_lb = b22[["2.5%"]], cases_2022_ub = b22[["97.5%"]])
}))
print(best); print(pooled, digits = 3); print(strata, digits = 3); print(burden, digits = 4); print(sens, digits = 3)
write.csv(best, file.path(DERIVED, "author_best_estimates_s6.csv"), row.names = FALSE)
write.csv(pooled, file.path(OUT_TAB, "author_best_pooled.csv"), row.names = FALSE)
write.csv(strata, file.path(OUT_TAB, "author_best_age_strata.csv"), row.names = FALSE)
write.csv(burden, file.path(OUT_TAB, "author_best_burden_projections.csv"), row.names = FALSE)
write.csv(regions, file.path(OUT_TAB, "author_best_projections_regions.csv"), row.names = FALSE)
write.csv(sens, file.path(OUT_TAB, "author_best_sensitivity_sets.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(OUT_DIAG, "sessioninfo_11_fit_author_best_estimates.txt"))
