# Produces: output/intermediate/ibge_projection_2024_60plus.csv; output/tables/projections_2018_vs_2024.csv
# Projeções 2030 a 2070 com a revisão 2024 das projeções do IBGE (base Censo 2022). Este script prepara os
# denominadores e compara as edições 2018 e 2024; a carga com a curva da análise principal, nacional e por
# região, sai de 03g_fit_author_best_estimates.R.
#
# O manuscrito cita "Population Projections ... 2000-2070 (2024)", mas os denominadores usados
# (ibge_projection_*_2030_2050_agegroups.csv) são da edição 2018 (SIDRA 7358), anterior ao Censo.
# Fonte nova: ftp.ibge.gov.br/Projecao_da_Populacao/Projecao_da_Populacao_2024/
#   projecoes_2024_tab2_grupo_quinquenal.xlsx (SHA-256 no log).
# Prevalência por idade mantida constante (mesmo pressuposto do manuscrito).
suppressPackageStartupMessages({ library(readxl); library(dplyr); library(tidyr); library(metafor); library(MASS, exclude = "select") })
source(file.path("code", "helpers", "setup.R"))
set.seed(20261019)

f24 <- file.path(RAW_IBGE_2024, "projecoes_2024_tab2_grupo_quinquenal.xlsx")
cat("SHA-256 da fonte 2024:", system2("shasum", c("-a", "256", shQuote(f24)), stdout = TRUE), "\n")
x <- suppressMessages(read_excel(f24, sheet = 1, skip = 5))
names(x)[1:5] <- c("age", "sex", "code", "abbr", "local")
p24 <- x |> filter(sex == "Ambos", local %in% c("Brasil", "Norte", "Nordeste", "Sudeste", "Sul", "Centro-Oeste"),
                  age %in% c("60-64", "65-69", "70-74", "75-79", "80-84", "85-89", "90+")) |>
  select(age, local, all_of(as.character(c(2022L, PROJECTION_YEARS)))) |>
  pivot_longer(all_of(as.character(c(2022L, PROJECTION_YEARS))), names_to = "year", values_to = "population") |>
  mutate(year = as.integer(year), population = as.numeric(population))
write.csv(p24, file.path(DERIVED, "ibge_projection_2024_60plus.csv"), row.names = FALSE)

cen  <- read.csv(file.path(OCT_IBGE, "ibge_census_brazil_2022_agegroups.csv"))
p18  <- read.csv(file.path(OCT_IBGE, "ibge_projection_brazil_2030_2050_agegroups.csv"))
den <- bind_rows(
  transmute(cen, source = "Census 2022", year = 2022L, age = age_group, population),
  transmute(p18, source = "IBGE projection, 2018 edition", year, age = age_group, population),
  p24 |> filter(local == "Brasil") |> transmute(source = "IBGE projection, 2024 revision", year, age, population))
tot_den <- den |> group_by(source, year) |> summarise(pop60 = sum(population), .groups = "drop")
print(tot_den)

# Prevalências por idade: (a) original do manuscrito (7 faixas); (b) corrigida (6 faixas, 85+ aberta)
orig <- read.csv(file.path(OUT_TAB, "age_specific_reproduced.csv"))
prev_orig <- setNames(orig$est, orig$age)
agec <- read.csv(file.path(DERIVED, "age_strata_corrected.csv"))
bands <- c("60-64", "65-69", "70-74", "75-79", "80-84", "85+")
agec$band <- factor(agec$band, levels = bands)
mc <- rma.mv(yi, vi, mods = ~ 0 + band, random = ~ 1 | study / row, data = agec, method = "REML")
draws <- plogis(mvrnorm(20000, coef(mc), vcov(mc)))
colnames(draws) <- bands

to6 <- function(d) {   # agrega 85-89 + 90+ em 85+
  d |> mutate(band = ifelse(age %in% c("85-89", "90+"), "85+", age)) |> group_by(band) |>
    summarise(population = sum(population), .groups = "drop") |> arrange(match(band, bands))
}
res <- list()
for (src in unique(den$source)) for (yr in unique(den$year[den$source == src])) {
  d <- den |> filter(source == src, year == yr)
  orig_cases <- sum(prev_orig[d$age] * d$population)
  d6 <- to6(d); cs <- as.vector(draws[, d6$band] %*% d6$population)
  res[[length(res) + 1]] <- data.frame(source = src, year = yr, pop60 = sum(d$population),
    cases_original_schedule = round(orig_cases),
    cases_corrected_median = round(median(cs)), cases_corrected_lb = round(quantile(cs, 0.025)),
    cases_corrected_ub = round(quantile(cs, 0.975)))
}
res <- bind_rows(res)

# Decomposição simétrica (tamanho x estrutura etária), base Censo 2022, com a revisão 2024
dec <- function(C0, P0, C1, P1) { size <- (P1 - P0) * (C0 / P0 + C1 / P1) / 2; c(size = size, structure = (C1 - C0) - size) }
b <- res |> filter(source == "Census 2022")
for (yr in PROJECTION_YEARS) {
  t <- res |> filter(source == "IBGE projection, 2024 revision", year == yr)
  dd <- dec(b$cases_corrected_median, b$pop60, t$cases_corrected_median, t$pop60)
  res$growth_size_component[res$source == t$source & res$year == yr] <- round(dd["size"])
  res$growth_age_structure_component[res$source == t$source & res$year == yr] <- round(dd["structure"])
}
write.csv(res, file.path(OUT_TAB, "projections_2018_vs_2024.csv"), row.names = FALSE)
print(res)

# As contagens regionais saem de 03g_fit_author_best_estimates.R, com a curva da análise principal
# (author_best_projections_regions.csv).
writeLines(capture.output(sessionInfo()), file.path(OUT_DIAG, "sessioninfo_06_project_ibge_2024.txt"))
