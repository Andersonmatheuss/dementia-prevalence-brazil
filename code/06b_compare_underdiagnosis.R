# 22 — Quanto o denominador muda a estimativa publicada de subdiagnóstico (20/10/2026, post hoc, pedido do PI).
# Miguel et al. 2026 (Age Ageing 55:afag163, doi:10.1093/ageing/afag163, PMID 42248800; artigo e material
# suplementar lidos em 20/10/2026): subdiagnóstico = 1 - diagnosticados / esperados.
#   Diagnosticados: usuários de anticolinesterásico no SUS divididos, em sequência, por 0,70 (obtêm pelo SUS),
#   0,80 (elegíveis que usam), 0,80 (Alzheimer leve ou moderado) e 0,60 (Alzheimer entre as demências).
#   Esperados: 8,5% (Delphi) da população de 60+ da projeção do IBGE (revisão 2024).
#   Tabelas S8-S10 do suplemento (60+, total): usuários, população e esperados de 2019, 2022 e 2024, abaixo.
# A reconstrução reproduz os três valores publicados (81,6%, 84,9% e 84,3%). Os diagnosticados de 2022 não mudam
# com o denominador; o subdiagnóstico de cada estimativa nacional de 2022 é 1 - diagnosticados / casos dessa
# estimativa (as nossas, no Censo 2022). Determinístico: não há sorteio.
suppressPackageStartupMessages(library(dplyr))
source(file.path("code", "helpers", "setup.R"))

premises <- 0.70 * 0.80 * 0.80 * 0.60
pub <- data.frame(year = c(2019, 2022, 2024), users = c(121778, 109649, 122323), population = c(28950761, 31864499, 34169617),
                  expected = c(2460815, 2708482, 2904417), published = c(0.816, 0.849, 0.843))   # Tabelas S8-S10 e texto
pub$underdiagnosis <- 1 - pub$users / premises / pub$expected
stopifnot("A reconstrução reproduz os valores publicados" = all(abs(round(pub$underdiagnosis, 3) - pub$published) < 1e-9))
diagnosed <- pub$users[pub$year == 2022] / premises

pop60 <- read.csv(file.path(OUT_TAB, "author_best_burden_projections.csv")) %>% filter(deff == 1, year == 2022) %>% pull(population)
est <- read.csv(file.path(OUT_TAB, "primary_burden_estimands.csv"))
pba <- read.csv(file.path(OUT_TAB, "pba_screening_correction.csv"))
tri <- read.csv(file.path(OUT_TAB, "triangulation_national_estimates.csv"))
pick <- function(d, ...) { r <- filter(d, ...); stopifnot(nrow(r) == 1); r }
src <- bind_rows(
  data.frame(denominator = "Delphi 8.5% on the 2022 IBGE projection (as published)", population = pub$population[pub$year == 2022],
             cases = pub$expected[pub$year == 2022], lb = NA, ub = NA),
  data.frame(denominator = "ELSI-Brazil weighted rates x Census 2022", population = pop60,
             cases = pick(tri, grepl("weighted age-specific", source))$cases, lb = NA, ub = NA),
  with(pick(est, version == "observed", estimand == "typical study"),
       data.frame(denominator = "This analysis, observed counts", population = pop60, cases = cases_median, lb = cases_lb, ub = cases_ub)),
  with(pick(est, version == "author_best", estimand == "typical study"),
       data.frame(denominator = "This analysis, primary (authors' estimates)", population = pop60, cases = cases_median, lb = cases_lb, ub = cases_ub)),
  with(pick(pba, variant == "partial_pooling", quantity == "total error (verification and random)"),
       data.frame(denominator = "This analysis, probabilistic bias analysis", population = pop60, cases = cases_median, lb = cases_lb, ub = cases_ub))
)
und <- function(cases) 1 - diagnosed / cases
out <- src %>% mutate(prevalence = cases / population, underdiagnosis = und(cases), underdiagnosis_lb = und(lb), underdiagnosis_ub = und(ub),
                      diagnosed_2022 = diagnosed)
print(pub, digits = 4); print(out, digits = 4)
write.csv(out, file.path(OUT_TAB, "underdiagnosis_by_denominator_2022.csv"), row.names = FALSE)
