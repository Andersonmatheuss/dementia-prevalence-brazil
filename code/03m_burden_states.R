# Produces: output/intermediate/ibge_projection_2024_60plus_states.csv; output/tables/burden_states.csv
# Carga esperada por Unidade da Federação (27 UFs), 2022 e 2030-2070: a curva etária da análise principal
# (sorteios de 03g, efeito de desenho 1) aplicada à estrutura etária de cada UF na revisão 2024 das projeções
# do IBGE. Como nas macrorregiões, é contagem esperada sob uma curva nacional comum: as diferenças entre UFs
# vêm só da estrutura etária, não de diferenças de risco. 2022 = ano 2022 da projeção (como em 03g para as
# regiões), para que UFs, regiões e anos usem a mesma fonte.
suppressPackageStartupMessages({ library(readxl); library(dplyr); library(tidyr) })
source(file.path("code", "helpers", "setup.R"))

bands <- c("60-64", "65-69", "70-74", "75-79", "80-84", "85+")
dr <- as.matrix(read.csv(file.path(DERIVED, "primary_age_schedule_draws.csv"), check.names = FALSE))[, bands]
yrs <- c(2022L, PROJECTION_YEARS)

region_of <- c(RO = "Norte", AC = "Norte", AM = "Norte", RR = "Norte", PA = "Norte", AP = "Norte", TO = "Norte",
               MA = "Nordeste", PI = "Nordeste", CE = "Nordeste", RN = "Nordeste", PB = "Nordeste", PE = "Nordeste",
               AL = "Nordeste", SE = "Nordeste", BA = "Nordeste", MG = "Sudeste", ES = "Sudeste", RJ = "Sudeste",
               SP = "Sudeste", PR = "Sul", SC = "Sul", RS = "Sul", MS = "Centro-Oeste", MT = "Centro-Oeste",
               GO = "Centro-Oeste", DF = "Centro-Oeste")

x <- suppressMessages(read_excel(file.path(RAW_IBGE_2024, "projecoes_2024_tab2_grupo_quinquenal.xlsx"), sheet = 1, skip = 5))
names(x)[1:5] <- c("age", "sex", "code", "abbr", "local")
x$age <- trimws(x$age)
# Rio Grande do Sul vem sem sigla na planilha (só o nome): completada pelo nome
x$abbr[x$local == "Rio Grande do Sul"] <- "RS"
st <- x |> filter(sex == "Ambos", abbr %in% names(region_of), age %in% c(bands[1:5], "85-89", "90+")) |>
  select(age, code, abbr, local, all_of(as.character(yrs))) |>
  pivot_longer(all_of(as.character(yrs)), names_to = "year", values_to = "population") |>
  mutate(year = as.integer(year), population = as.numeric(population),
         band = ifelse(age %in% c("85-89", "90+"), "85+", age)) |>
  group_by(code, abbr, local, year, band) |> summarise(population = sum(population), .groups = "drop")
stopifnot(n_distinct(st$abbr) == 27, !anyNA(st$population))
write.csv(st, file.path(DERIVED, "ibge_projection_2024_60plus_states.csv"), row.names = FALSE)

q <- function(v) c(median = median(v), lb = unname(quantile(v, .025)), ub = unname(quantile(v, .975)))
res <- bind_rows(lapply(unique(st$abbr), function(a) {
  s <- st[st$abbr == a, ]
  pv <- function(y) { z <- s[s$year == y, ]; z$population[match(bands, z$band)] }
  c22 <- as.vector(dr %*% pv(2022L))
  bind_rows(lapply(yrs, function(y) {
    p <- pv(y); cs <- as.vector(dr %*% p); ch <- cs / c22 - 1
    data.frame(code = s$code[1], abbr = a, state = s$local[1], region = region_of[[a]], year = y,
               pop60 = sum(p), share80plus = sum(p[5:6]) / sum(p),
               cases_median = q(cs)[[1]], cases_lb = q(cs)[[2]], cases_ub = q(cs)[[3]],
               prev_median = q(cs / sum(p))[[1]], prev_lb = q(cs / sum(p))[[2]], prev_ub = q(cs / sum(p))[[3]],
               change_vs_2022_median = if (y == 2022L) NA else q(ch)[[1]],
               change_vs_2022_lb = if (y == 2022L) NA else q(ch)[[2]],
               change_vs_2022_ub = if (y == 2022L) NA else q(ch)[[3]])
  }))
}))
res <- res |> arrange(year, code)

# Conferência: a soma das UFs reproduz a população 60+ das macrorregiões (03b) e a carga regional de 03g
reg_pop <- read.csv(file.path(DERIVED, "ibge_projection_2024_60plus.csv")) |> filter(local != "Brasil") |>
  group_by(local, year) |> summarise(pop = sum(population), .groups = "drop")
chk <- res |> group_by(region, year) |> summarise(pop = sum(pop60), .groups = "drop") |>
  inner_join(reg_pop, by = c("region" = "local", "year"))
stopifnot(nrow(chk) == 5 * length(yrs), all(abs(chk$pop.x / chk$pop.y - 1) < 1e-4))
reg <- read.csv(file.path(OUT_TAB, "author_best_projections_regions.csv"))
chk2 <- bind_rows(lapply(split(st, list(st$year, region_of[st$abbr]), drop = TRUE), function(z) {
  pv <- tapply(z$population, z$band, sum)[bands]
  data.frame(region = region_of[[z$abbr[1]]], year = z$year[1], cases = median(as.vector(dr %*% pv)))
})) |> inner_join(reg, by = c("region" = "local", "year"))
cat("Máxima diferença relativa entre UFs somadas e a carga regional de 03g:",
    signif(max(abs(chk2$cases / chk2$cases_median - 1)), 3), "\n")
stopifnot(max(abs(chk2$cases / chk2$cases_median - 1)) < 1e-4)

print(res |> filter(year %in% c(2022L, 2050L)) |> select(abbr, year, pop60, cases_median, prev_median, change_vs_2022_median), digits = 3)
write.csv(res, file.path(OUT_TAB, "burden_states.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(OUT_DIAG, "sessioninfo_03m_burden_states.txt"))
