# Produces: output/intermediate/census2022_literacy_60plus.csv; output/tables/education_standardised_burden_2022.csv, education_standardisation_inputs.csv
# Carga de 2022 padronizada pela escolaridade (alfabetização), nacional, por macrorregião e por UF.
# Os cinco estudos do conjunto principal são do Sudeste e do Sul; o país e cada UF têm proporções diferentes de
# pessoas idosas não alfabetizadas. A curva etária principal (sorteios de 03g) é reponderada por
#   p[a, UF] = p[a] * {1 + i[a, UF] (RR - 1)} / {1 + i[a, BR] (RR - 1)} * {1 + i[BR] (RR - 1)} / {1 + i[S] (RR - 1)}
# onde i = proporção de não alfabetizados (por faixa a, na UF e no Brasil; i[BR] e i[S] aos 60+, no Brasil e
# nas amostras dos estudos) e RR = razão de prevalências não alfabetizados vs alguma escolaridade, agrupada nos
# estudos do conjunto principal que a relatam. O RR é bruto (não ajustado por idade) e parte do gradiente pode
# vir de viés dos instrumentos cognitivos; o resultado é um limite superior do efeito da escolaridade, e a carga
# sem padronização continua sendo a principal. Só 2022: a escolaridade das coortes futuras muda.
# Denominadores: Censo 2022, SIDRA tabela 9542 (pessoas de 15 anos ou mais, total e alfabetizadas, por grupos de
# idade; universo), baixada uma vez e guardada em output/intermediate (como o mapa de 05d). Não há 65-69 e 70-74
# separados na tabela: as faixas 65-69 e 70-74 usam a proporção de 65-74, e 80-84 e 85+ a de 80+.
suppressPackageStartupMessages({ library(dplyr); library(metafor); library(jsonlite) })
source(file.path("code", "helpers", "setup.R"))
set.seed(20261019)

bands <- c("60-64", "65-69", "70-74", "75-79", "80-84", "85+")
dr <- as.matrix(read.csv(file.path(DERIVED, "primary_age_schedule_draws.csv"), check.names = FALSE))[, bands]
nd <- nrow(dr)

# --- Censo 2022: alfabetização por idade, Brasil e UFs (cache) ---------------------------------------------
cache <- file.path(DERIVED, "census2022_literacy_60plus.csv")
if (!file.exists(cache)) {
  url <- paste0("https://apisidra.ibge.gov.br/values/t/9542/n1/all/n3/all/v/950/p/2022/c59/1023,1024/c2/6794/c86/95251/",
                "c287/6613,6614,6615,6616,6617,3000,9486,113623")
  j <- fromJSON(url); hdr <- unlist(j[1, ]); j <- j[-1, ]
  col_of <- function(codes) names(j)[sapply(names(j), function(nm) all(codes %in% j[[nm]]))][1]
  c_lit <- col_of(c("1023", "1024")); c_age <- col_of(c("6613", "113623"))
  raw <- data.frame(geo_code = as.integer(j$D1C), geo = j$D1N, literacy = ifelse(j[[c_lit]] == "1023", "literate", "illiterate"),
                    age_code = j[[c_age]], persons = as.numeric(j$V))
  write.csv(raw, cache, row.names = FALSE)
}
lit <- read.csv(cache)
stopifnot(n_distinct(lit$geo_code) == 28, !anyNA(lit$persons))
# Faixas derivadas: 60-64 (idades simples), 65-74 (65+ menos 75+), 75-79 (75+ menos 80+), 80+
g <- lit |> group_by(geo_code, geo, literacy) |> summarise(
  `60-64` = sum(persons[age_code %in% 6613:6617]),
  `65-74` = persons[age_code == 3000] - persons[age_code == 9486],
  `75-79` = persons[age_code == 9486] - persons[age_code == 113623],
  `80+`   = persons[age_code == 113623], .groups = "drop")
lb <- c("60-64", "65-74", "75-79", "80+")
ill <- g |> filter(literacy == "illiterate"); tot <- g |> group_by(geo_code, geo) |> summarise(across(all_of(lb), sum), .groups = "drop")
ill <- ill[match(tot$geo_code, ill$geo_code), ]
share <- as.matrix(ill[, lb]) / as.matrix(tot[, lb]); rownames(share) <- tot$geo_code
pop60_census <- rowSums(as.matrix(tot[, lb])); names(pop60_census) <- tot$geo_code
i60 <- rowSums(as.matrix(ill[, lb])) / pop60_census
cat("População 60+ do Censo (tabela 9542), Brasil:", pop60_census["1"], "\n")
to6 <- function(v) c(v[1], v[2], v[2], v[3], v[4], v[4])   # 4 faixas do Censo -> 6 faixas da curva

# --- Estudos: proporção de não alfabetizados e razão de prevalências ----------------------------------------
st <- read.csv(VERIFIED_STRATA) |> filter(study %in% PRIMARY_POPULATIONS, domain == "education", !is.na(n))
sx <- st |> group_by(study, category) |> summarise(n = sum(n), cases = sum(cases), .groups = "drop")
w <- data.frame(study = unique(sx$study))
w$n_none <- sx$n[sx$category == "none"][match(w$study, sx$study[sx$category == "none"])]
w$c_none <- sx$cases[sx$category == "none"][match(w$study, sx$study[sx$category == "none"])]
w$n_any  <- sx$n[sx$category == "any"][match(w$study, sx$study[sx$category == "any"])]
w$c_any  <- sx$cases[sx$category == "any"][match(w$study, sx$study[sx$category == "any"])]
es <- escalc("RR", ai = c_none, n1i = n_none, ci = c_any, n2i = n_any, data = w)
m <- rma(yi, vi, data = es, method = "REML", test = "knha")
iS <- sum(w$n_none) / sum(w$n_none + w$n_any)
rr_draw <- exp(rnorm(nd, coef(m)[[1]], m$se))   # incerteza do RR, independente da curva etária
cat(sprintf("RR agrupado (k = %d): %.2f (%.2f a %.2f); não alfabetizados nas amostras: %.1f%%; Brasil 60+: %.1f%%\n",
            m$k, exp(coef(m)), exp(m$ci.lb), exp(m$ci.ub), 100 * iS, 100 * i60["1"]))

# --- Carga padronizada -------------------------------------------------------------------------------------
states <- read.csv(file.path(DERIVED, "ibge_projection_2024_60plus_states.csv")) |> filter(year == 2022L)
region_of <- c(`1` = "Norte", `2` = "Nordeste", `3` = "Sudeste", `4` = "Sul", `5` = "Centro-Oeste")
q <- function(v) c(median(v), quantile(v, c(.025, .975), names = FALSE))
cen <- read.csv(file.path(OCT_IBGE, "ibge_census_brazil_2022_agegroups.csv"))
pop_br <- c(cen$population[match(bands[1:5], cen$age_group)], sum(cen$population[cen$age_group %in% c("85-89", "90+")]))
adj_draws <- function(geo) {   # matriz nd x 6 de prevalências padronizadas
  iu <- to6(share[as.character(geo), ]); ib <- to6(share["1", ])
  f_nat <- (1 + i60["1"] * (rr_draw - 1)) / (1 + iS * (rr_draw - 1))
  f_age <- sapply(seq_along(bands), function(k) (1 + iu[k] * (rr_draw - 1)) / (1 + ib[k] * (rr_draw - 1)))
  dr * f_age * f_nat
}
row_out <- function(level, code, name, pop6) {
  un <- as.vector(dr %*% pop6); ad <- as.vector(adj_draws(code) %*% pop6); r <- ad / un
  data.frame(level = level, code = code, name = name, pop60 = sum(pop6),
             illiterate_60plus = if (as.character(code) %in% names(i60)) i60[[as.character(code)]] else NA,
             cases_unadjusted_median = q(un)[1], cases_unadjusted_lb = q(un)[2], cases_unadjusted_ub = q(un)[3],
             cases_education_median = q(ad)[1], cases_education_lb = q(ad)[2], cases_education_ub = q(ad)[3],
             prev_unadjusted = q(un)[1] / sum(pop6), prev_education = q(ad)[1] / sum(pop6),
             ratio_median = q(r)[1], ratio_lb = q(r)[2], ratio_ub = q(r)[3])
}
uf <- bind_rows(lapply(unique(states$code), function(cd) {
  s <- states[states$code == cd, ]; row_out("state", cd, s$abbr[1], s$population[match(bands, s$band)]) }))
# Regiões: soma dos sorteios das UFs (cada UF com a sua proporção de não alfabetizados)
reg <- bind_rows(lapply(names(region_of), function(rc) {
  cds <- unique(states$code[states$code %/% 10 == as.integer(rc)])
  pops <- lapply(cds, function(cd) { s <- states[states$code == cd, ]; s$population[match(bands, s$band)] })
  un <- Reduce(`+`, lapply(pops, function(p) as.vector(dr %*% p)))
  ad <- Reduce(`+`, Map(function(cd, p) as.vector(adj_draws(cd) %*% p), cds, pops)); r <- ad / un
  P <- sum(unlist(pops))
  data.frame(level = "region", code = as.integer(rc), name = region_of[[rc]], pop60 = P, illiterate_60plus = NA,
             cases_unadjusted_median = q(un)[1], cases_unadjusted_lb = q(un)[2], cases_unadjusted_ub = q(un)[3],
             cases_education_median = q(ad)[1], cases_education_lb = q(ad)[2], cases_education_ub = q(ad)[3],
             prev_unadjusted = q(un)[1] / P, prev_education = q(ad)[1] / P, ratio_median = q(r)[1], ratio_lb = q(r)[2], ratio_ub = q(r)[3])
}))
nat <- row_out("country", 1L, "Brasil", pop_br)
out <- bind_rows(nat, reg, uf)
inputs <- data.frame(quantity = c("pooled prevalence ratio, illiterate vs any schooling (crude)", "RR 95% CI lower", "RR 95% CI upper",
                                  "k studies", "tau2", "illiterate share, pooled study samples", "illiterate share, Brazil 60+ (Census 2022)",
                                  "Census 2022 population 60+, table 9542"),
                     value = c(exp(coef(m)[[1]]), exp(m$ci.lb), exp(m$ci.ub), m$k, m$tau2, iS, i60[["1"]], pop60_census[["1"]]))
print(es[, c("study", "n_none", "c_none", "n_any", "c_any")]); print(inputs, digits = 4)
print(out[, c("level", "name", "illiterate_60plus", "cases_unadjusted_median", "cases_education_median", "ratio_median")], digits = 3)
write.csv(out, file.path(OUT_TAB, "education_standardised_burden_2022.csv"), row.names = FALSE)
write.csv(inputs, file.path(OUT_TAB, "education_standardisation_inputs.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(OUT_DIAG, "sessioninfo_03p_education_standardisation.txt"))
