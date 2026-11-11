# Verificações dos resultados existentes, sem reajustar modelos.
source(file.path("code", "helpers", "setup.R"))

read_table <- function(filename, required) {
  file <- file.path(OUT_TAB, filename)
  stopifnot("CSV de teste ausente" = file.exists(file))
  data <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop(filename, ": colunas ausentes: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  data
}

one_row <- function(data, condition, label) {
  index <- which(!is.na(condition) & condition)
  if (length(index) != 1L) {
    stop(label, ": esperada exatamente uma linha; encontradas ", length(index), call. = FALSE)
  }
  data[index, , drop = FALSE]
}

assert_close <- function(actual, expected, tolerance, label) {
  stopifnot("Valor deve ser numérico, único e finito" =
              is.numeric(actual) && length(actual) == 1L && is.finite(actual))
  if (abs(actual - expected) > tolerance) {
    stop(sprintf("%s: obtido %.10g; esperado %.10g; tolerância absoluta %.10g.",
                 label, actual, expected, tolerance), call. = FALSE)
  }
  stopifnot(abs(actual - expected) <= tolerance)
}

replication <- read_table("replication_frequentist.csv", "status")
stopifnot("Status inesperado" = all(replication$status %in% c("reproduced", "new")),
          "Esperadas 72 linhas reproduzidas" = sum(replication$status == "reproduced") == 72L)
cat(sprintf("Replicação frequentista: %d reproduced; %d new.\n",
            sum(replication$status == "reproduced"), sum(replication$status == "new")))

discrepancies <- read_table("raw_chain_discrepancies.csv", character())
stopifnot("A cadeia bruta contém discrepâncias" = nrow(discrepancies) == 0L)
cat("Cadeia bruta: 0 discrepâncias.\n")

overall <- read_table("corrected_overall_by_set.csv", c("set", "prev", "k"))
s0 <- one_row(overall, overall$set == "S0", "Conjunto S0")
s4 <- one_row(overall, overall$set == "S4", "Conjunto S4")
assert_close(s0$prev, 0.0854, 0.0005, "S0 prev")
assert_close(s4$prev, 0.0831, 0.0005, "S4 prev")
stopifnot("S4 deve conter 7 estudos" = s4$k == 7L)

epifloripa <- read_table("corrected_overall_with_epifloripa.csv", c("set", "prev", "k"))
s6 <- one_row(epifloripa, startsWith(epifloripa$set, "S6"), "Conjunto S6")
assert_close(s6$prev, 0.0773, 0.0005, "S6 prev")
stopifnot("S6 deve conter 8 estudos" = s6$k == 8L)

pooled <- read_table("author_best_pooled.csv", c("deff", "prev"))
best <- one_row(pooled, pooled$deff == 1, "Estimativa dos autores, deff = 1")
# 19/10/2026: conjunto principal sem sobreposição (k = 5); antes 0,1001 com as oito populações
stopifnot("Conjunto principal com 5 populações" = best$k == 5L)
assert_close(best$prev, 0.0947, 0.0005, "Author-best prev")

burden <- read_table("author_best_burden_projections.csv", c("deff", "year", "population", "cases_median"))
# projeções até 2070 (revisão 2024 do IBGE); população 60+ de 2070 conferida na planilha do IBGE
stopifnot("Anos 2022 e 2030 a 2070" = all(c(2022, 2030, 2040, 2050, 2060, 2070) %in% burden$year[burden$deff == 1]))
assert_close(burden$population[burden$deff == 1 & burden$year == 2070], 75292150, 0.5, "População 60+ em 2070")
best_2022 <- one_row(burden, burden$deff == 1 & burden$year == 2022, "Carga 2022, deff = 1")
# 19/10/2026: 2.750.451 antes da extração conferida; 2.739.180 com as oito populações; 2.658.583 sem sobreposição
assert_close(best_2022$cases_median, 2658583, 0.02 * 2658583, "Author-best cases 2022")

projections <- read_table("projections_2018_vs_2024.csv", c("source", "year", "cases_original_schedule"))
projection_2030 <- one_row(projections,
  projections$source == "IBGE projection, 2024 revision" & projections$year == 2030,
  "IBGE revisão 2024, ano 2030")
assert_close(projection_2030$cases_original_schedule, 3007273, 1, "IBGE cases 2030")

# 19/10/2026, extração conferida no texto integral: o VPP passou a usar só os casos entre os positivos
# examinados (Bottino 101/164; Lopes 67/130). Antes: Lopes 0,1234 e Bottino 0,1296 (VPP com todos os casos).
qba <- read_table("qba_two_phase_documented_study_level.csv", c("study", "npv", "adjusted"))
lopes <- one_row(qba, qba$study == "Lopes_2012" & abs(qba$npv - 0.97) < 1e-10, "Lopes NPV 0.97")
bottino <- one_row(qba, qba$study == "Bottino_2008" & abs(qba$npv - 0.97) < 1e-10, "Bottino NPV 0.97")
assert_close(lopes$adjusted, 0.1220, 0.0005, "Lopes adjusted")
assert_close(bottino$adjusted, 0.1237, 0.0005, "Bottino adjusted")
stopifnot("Correa Ribeiro deve estar na análise de viés" = any(qba$study == "Ribeiro_2013"))

# Correções da extração conferida (fulltext_strata.csv) que precisam chegar ao modelo por idade
age <- read.csv(file.path(DERIVED, "age_strata_corrected.csv"))
rib85 <- one_row(age, age$study == "Ribeiro_2013" & age$band == "85+", "Ribeiro 85+")
stopifnot("Ribeiro 85+ deve ser 63/129 (85-89 + >90)" = rib85$cases == 63 && rib85$n == 129)
lop60 <- one_row(age, age$study == "Lopes_2012" & age$band == "60-64", "Lopes 60-64")
stopifnot("Lopes 60-64 deve ter n = 265 (Tabela 3)" = lop60$n == 265)
best_est <- read.csv(file.path(DERIVED, "author_best_estimates_s6.csv"))
assert_close(best_est$observed[best_est$study == "Bottino_2008"], 107 / 1477, 1e-9, "Bottino observed (age-table total)")
assert_close(best_est$observed[best_est$study == "Lopes_2012"], 67 / 1052, 1e-9, "Lopes observed (age-table total)")
# Análise principal sem sobreposição (00_setup.R, 19/10/2026): sem ELSI, Bottino nem Correa Ribeiro
sets <- read.csv(file.path(dir_proc, "analysis_sets.csv"))
stopifnot("Conjunto principal = 5 populações" = sum(sets$in_primary) == 5L,
          "ELSI, Bottino e Correa Ribeiro fora da principal" =
            !any(sets$in_primary[sets$study %in% c("Bertola_2023", "Bottino_2008", "Ribeiro_2013")]))
sens <- read_table("primary_burden_sensitivity_sets.csv", c("set", "version", "k", "cases_median"))
stopifnot("Sensibilidades: 5 conjuntos x 2 versões" = nrow(sens) == 10L)
all8 <- one_row(sens, sens$set == "all eight populations" & sens$version == "author_best", "Oito populações")
# as oito populações reproduzem o conjunto principal anterior (2.741.478 casos em 2022)
assert_close(all8$cases_median, 2741478, 0.02 * 2741478, "Author-best cases 2022, all eight populations")
# 65+ no Censo 2022: taxas ponderadas de Bertola contra o algoritmo harmonizado de Aliberti 2026 (20/10/2026)
p65 <- read_table("primary_elsi_algorithms_65plus.csv", c("source", "prev", "lb", "ub", "age_structure"))
stopifnot("Comparação 65+: 3 linhas calculadas e 2 publicadas" = nrow(p65) == 5L)
assert_close(p65$prev[grepl("^ELSI-Brazil, Bertola", p65$source)], 0.0776, 0.0005, "ELSI-Brazil (Bertola) 65+ on Census 2022")
cat("Referências numéricas: todas dentro das tolerâncias.\n")
cat("ALL TESTS PASSED\n")
