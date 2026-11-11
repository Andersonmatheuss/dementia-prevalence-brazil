# Produces: output/tables/projection_trend_scenarios.csv
# Projeções 2030-2070 sob três cenários de tendência da prevalência por idade, a partir da curva da análise
# principal (sorteios de 03g) aplicada a 2022:
#   constante  — pressuposto das projeções principais;
#   queda      — prevalência por idade cai 13% por década (razão 0,87 a cada 10 anos), a magnitude da queda
#                da incidência em sete coortes da Europa e dos EUA, 1988-2015 (Wolters et al., Neurology
#                2020;95:e519-31, doi:10.1212/WNL.0000000000010022), aplicada aqui à prevalência como cenário;
#   aumento    — simétrico (razão 1,13 por década), em linha com o aumento da prevalência padronizada em três
#                de cinco sítios latino-americanos do 10/66 entre 2003-2006 e 2016-2020 (Llibre-Guerra et al.,
#                JAMA Neurol 2026;83:906-14, doi:10.1001/jamaneurol.2026.2305).
# Os fatores são cenários, não estimativas de tendência no Brasil. A incerteza é a da curva etária; o fator
# entra sem incerteza própria.
suppressPackageStartupMessages({ library(dplyr) })
source(file.path("code", "helpers", "setup.R"))

bands <- c("60-64", "65-69", "70-74", "75-79", "80-84", "85+")
dr <- as.matrix(read.csv(file.path(DERIVED, "primary_age_schedule_draws.csv"), check.names = FALSE))[, bands]
p24 <- read.csv(file.path(DERIVED, "ibge_projection_2024_60plus.csv"))
cen <- read.csv(file.path(OCT_IBGE, "ibge_census_brazil_2022_agegroups.csv"))
pop_census <- c(cen$population[match(bands[1:5], cen$age_group)], sum(cen$population[cen$age_group %in% c("85-89", "90+")]))
to6 <- function(yr, loc) { x <- subset(p24, local == loc & year == yr)
  c(x$population[match(bands[1:5], x$age)], sum(x$population[x$age %in% c("85-89", "90+")])) }
scen <- c("constant" = 1, "decline 13% per decade" = 0.87, "increase 13% per decade" = 1.13)
q <- function(v) c(median(v), quantile(v, c(.025, .975), names = FALSE))

res <- bind_rows(lapply(c("Brasil", "Norte", "Nordeste", "Sudeste", "Sul", "Centro-Oeste"), function(loc) {
  bind_rows(lapply(names(scen), function(sn) bind_rows(lapply(c(2022L, PROJECTION_YEARS), function(y) {
    # 2022: Censo para o Brasil (como na carga nacional), ano 2022 da projeção para as regiões (como em 03g)
    p <- if (y == 2022L && loc == "Brasil") pop_census else to6(y, loc)
    f <- scen[[sn]]^((y - 2022) / 10)
    cs <- q(as.vector(dr %*% p) * f)
    data.frame(local = loc, scenario = sn, prevalence_ratio_vs_2022 = f, year = y, pop60 = sum(p),
               cases_median = cs[1], cases_lb = cs[2], cases_ub = cs[3])
  }))))
}))
# Conferência: o cenário constante reproduz a carga nacional de 03g
ab <- read.csv(file.path(OUT_TAB, "author_best_burden_projections.csv")) |> filter(deff == 1)
k <- res |> filter(local == "Brasil", scenario == "constant") |> inner_join(ab, by = "year")
stopifnot(nrow(k) == 6, max(abs(k$cases_median.x / k$cases_median.y - 1)) < 1e-8)
print(res |> filter(local == "Brasil", year %in% c(2022L, 2030L, 2050L, 2070L)), digits = 3)
write.csv(res, file.path(OUT_TAB, "projection_trend_scenarios.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(OUT_DIAG, "sessioninfo_03n_prevalence_trend_scenarios.txt"))
