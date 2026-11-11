# 03 — Reprodução independente dos resultados frequentistas do manuscrito de 18/10/2026.
#
# O script que gerou esses números não existe em disco (procurado no Pro, no T7 e no
# GitHub em 19/10/2026). Aqui a análise é reimplementada a partir das ENTRADAS por
# estudo salvas em out/2026 (meta_ready*.csv), seguindo o que o Methods descreve:
# logito, variância binomial, efeitos aleatórios REML (metafor), GLMM como sensibilidade.
# Cada número é comparado com o valor impresso no manuscrito (linha do texto convertido).
#
# Saídas: output/tables/replication_frequentist.csv (reportado x reproduzido)
#         output/tables/age_specific_reproduced.csv
suppressPackageStartupMessages({ library(metafor); library(dplyr); library(tidyr) })
source(file.path("code", "helpers", "setup.R"))

prim <- read.csv(file.path(OCT_OUT, "study_level_data", "meta_ready_primary_community.csv"))
allr <- read.csv(file.path(OCT_OUT, "study_level_data", "meta_ready.csv"))
flag <- read.csv(file.path(OCT_OUT, "study_level_data", "study_analysis_flags.csv"))
jbi  <- read.csv(file.path(dir_proc, "risk_of_bias.csv"))
pop  <- read.csv(file.path(OCT_IBGE, "ibge_census_brazil_2022_agegroups.csv"))
popr <- read.csv(file.path(OCT_IBGE, "ibge_census_regions_2022_agegroups.csv"))
proj <- read.csv(file.path(OCT_IBGE, "ibge_projection_brazil_2030_2050_agegroups.csv"))

# Recalcula yi/vi das contagens (não confia nos yi/vi salvos) e confere com os salvos.
recalc <- function(d) {
  cnt <- !is.na(d$cases) & !is.na(d$total)
  es  <- escalc(measure = "PLO", xi = d$cases[cnt], ni = d$total[cnt])
  stopifnot(max(abs(es$yi - d$yi[cnt])) < 1e-10, max(abs(es$vi - d$vi[cnt])) < 1e-10)
  d
}
prim <- recalc(prim); allr <- recalc(allr)

rows <- list()
add <- function(item, reported, reproduced, ms_line, note = "") {
  rows[[length(rows) + 1]] <<- data.frame(item = item, reported = reported,
    reproduced = reproduced, abs_diff = abs(reported - reproduced), ms_line = ms_line, note = note)
}

## 1. Análise primária (k = 11) ---------------------------------------------------
ov <- subset(prim, outcome == "dementia" & domain == "overall")
m0 <- fit_re(ov); p0 <- pred_p(m0)
add("primary k", 11, m0$k, 62)
add("primary N", 18642, sum(ov$total), 62)
add("primary prevalence %", 8.5, round(100 * p0["est"], 1), 62)
add("primary CI lower %", 4.9, round(100 * p0["lb"], 1), 62)
add("primary CI upper %", 14.5, round(100 * p0["ub"], 1), 62)
add("primary PI lower %", 1.2, round(100 * p0["pi_lb"], 1), 62)
add("primary PI upper %", 42.2, round(100 * p0["pi_ub"], 1), 62)
add("primary I2 %", 99.1, round(m0$I2, 1), 62)
add("primary tau2", 1.01, round(m0$tau2, 2), 62)

# GLMM binomial-normal (ML), sensibilidade
g0 <- rma.glmm(measure = "PLO", xi = cases, ni = total, data = ov)
add("GLMM prevalence %", 8.5, round(100 * transf.ilogit(g0$beta[1]), 1), 68)

# Influência e leave-one-out
inf <- influence(m0)
infl_names <- ov$Study[inf$is.infl]
add("influential study is Magalhaes (1=yes)", 1, as.numeric(identical(infl_names, "Magalhaes_2008")), 68)
loo <- leave1out(m0, transf = transf.ilogit)
add("excluding influential study %", 6.9, round(100 * loo$estimate[ov$Study == "Magalhaes_2008"], 1), 68)

# Egger (sei como preditor)
eg <- regtest(m0, predictor = "sei")
add("Egger p", 0.148, round(eg$pval, 3), 42)

## 2. Sensibilidades -----------------------------------------------------------------
sens <- function(d) round(100 * pred_p(fit_re(d))["est"], 1)
add("DSM-IV only %", 6.1, sens(subset(ov, DiagnosticGroup == "DSM-IV")), 68)
add("65+ only %", 6.4, sens(subset(ov, AgeLower >= 65)), 68)
ov_all <- subset(allr, outcome == "dementia" & domain == "overall" &
                   Study %in% flag$Study[flag$included_all_settings_overall_sensitivity])
add("all-settings k", 14, nrow(ov_all), 68)
add("all-settings %", 10.1, sens(ov_all), 68)
rob <- setNames(jbi$risk_of_bias, jbi$Study)
add("low-risk only %", 10.8, sens(subset(ov, rob[Study] == "low")), 68)
add("low/moderate risk %", 9.2, sens(subset(ov, rob[Study] %in% c("low", "moderate"))), 68)

## 3. Metarregressão -----------------------------------------------------------------
ov$year_c <- ov$Year - mean(ov$Year)
ov$diag_other <- as.numeric(ov$DiagnosticGroup != "DSM-IV")
mr <- rma(yi, vi, mods = ~ year_c + diag_other, data = ov, method = "REML")
add("MR primary beta diagnostic", 1.78, round(coef(mr)["diag_other"], 2), 68)
add("MR primary beta year", -0.05, round(coef(mr)["year_c"], 2), 68)
add("MR primary pseudo-R2 %", 31.2, round(100 * (m0$tau2 - mr$tau2) / m0$tau2, 1), 68)
add("MR primary residual tau2", 0.69, round(mr$tau2, 2), 68)
add("MR primary QM", 6.45, round(mr$QM, 2), 68)
ov_all$year_c <- ov_all$Year - mean(ov_all$Year)
ov_all$diag_other <- as.numeric(ov_all$DiagnosticGroup != "DSM-IV")
ov_all$hc <- as.numeric(ov_all$RecruitmentSetting == "health care-based")
mra <- rma(yi, vi, mods = ~ year_c + diag_other + hc, data = ov_all, method = "REML")
add("MR all-settings beta diagnostic", 1.66, round(coef(mra)["diag_other"], 2), 68)
add("MR all-settings beta health care", 0.56, round(coef(mra)["hc"], 2), 68)
# Sensibilidade nova (não está no manuscrito): metarregressão sem Magalhães
mr_nomag <- rma(yi, vi, mods = ~ year_c + diag_other, data = subset(ov, Study != "Magalhaes_2008"), method = "REML")

## 4. Prevalência por idade e carga nacional -------------------------------------------
ages <- c("60-64", "65-69", "70-74", "75-79", "80-84", "85-89", "90+")
age_d <- subset(prim, outcome == "dementia" & domain == "age" & subgroup %in% ages)
age_tab <- bind_rows(lapply(ages, function(a) {
  d <- subset(age_d, subgroup == a); m <- fit_re(d); p <- pred_p(m)
  data.frame(age = a, k = m$k, n_with_counts = sum(d$total, na.rm = TRUE),
             studies = paste(d$Study, collapse = ";"), est = p["est"], lb = p["lb"], ub = p["ub"],
             pi_lb = p["pi_lb"], pi_ub = p["pi_ub"], tau2 = m$tau2)
}))
age_tab <- left_join(age_tab, transmute(pop, age = age_group, population), by = "age") |>
  mutate(cases = est * population, cases_lb = lb * population, cases_ub = ub * population)
write.csv(age_tab, file.path(OUT_TAB, "age_specific_reproduced.csv"), row.names = FALSE)

rep_age <- data.frame(age = ages, k = c(4, 8, 8, 8, 9, 8, 4), N = c(789, 3157, 2503, 1625, 1093, 609, 58),
                      prev = c(3.1, 2.6, 4.4, 8.9, 16.5, 25.6, 51.5))
for (i in seq_along(ages)) {
  add(paste("age", ages[i], "k"), rep_age$k[i], age_tab$k[i], 208)
  add(paste("age", ages[i], "N"), rep_age$N[i], age_tab$n_with_counts[i], 208)
  add(paste("age", ages[i], "prevalence %"), rep_age$prev[i], round(100 * age_tab$est[i], 1), 208)
}
add("sum of age-stratum N", 9834, sum(age_tab$n_with_counts), 208, "N do estrato soma 9.834, nao 18.642")
nat <- sum(age_tab$cases)
add("national burden 2022", 2267772, round(nat), 62)
add("national burden lower (sum of stratum lower bounds)", 1450858, round(sum(age_tab$cases_lb)), 62,
    "o 'IC' do manuscrito e a soma dos limites inferiores dos estratos")
add("national burden upper (sum of stratum upper bounds)", 3428787, round(sum(age_tab$cases_ub)), 62)
add("implied national prevalence %", 7.1, round(100 * nat / sum(pop$population), 1), 62)
add("share of burden aged 80+ %", 50.8, round(100 * sum(age_tab$cases[5:7]) / nat, 1), 64)
add("share of burden aged 90+ %", NA, round(100 * age_tab$cases[7] / nat, 1), 208, "nao reportado; base: 58 pessoas")

## 5. Regiões, projeções, decomposição --------------------------------------------------
reg <- popr |> left_join(transmute(age_tab, age_group = age, est), by = "age_group") |>
  group_by(region) |> summarise(pop = sum(population), cases = sum(est * population)) |>
  mutate(prev = 100 * cases / pop)
add("Southeast 2022 cases", 1056094, round(reg$cases[reg$region == "Sudeste"]), 66)
add("Northeast 2022 cases", 592932, round(reg$cases[reg$region == "Nordeste"]), 66)
pj <- proj |> left_join(transmute(age_tab, age_group = age, est), by = "age_group") |>
  group_by(year) |> summarise(pop = sum(population), cases = sum(est * population))
add("projected 2030 cases", 3254496, round(pj$cases[pj$year == 2030]), 66)
add("projected 2050 cases", 6401585, round(pj$cases[pj$year == 2050]), 66)
# Decomposição simétrica (média das taxas brutas): componente de tamanho da população
dec <- function(y) {
  P0 <- sum(pop$population); C0 <- nat; P1 <- pj$pop[pj$year == y]; C1 <- pj$cases[pj$year == y]
  size <- (P1 - P0) * (C0 / P0 + C1 / P1) / 2
  c(size = size, structure = (C1 - C0) - size)
}
add("2030 growth component", 740089, round(dec(2030)["size"]), 66)
add("2050 ageing component", 1278312, round(dec(2050)["structure"]), 66)

## 6. Subgrupos -----------------------------------------------------------------------------
sg <- function(dom, lev) round(100 * pred_p(fit_re(subset(prim, outcome == "dementia" & domain == dom & subgroup == lev)))["est"], 1)
add("women %", 11.4, sg("sex", "Female"), 64)
add("men %", 9.2, sg("sex", "Male"), 64)
add("no formal education %", 24.6, sg("education", "Illiterate"), 64)
add("tertiary education %", 3.1, sg("education", "College"), 64)
add("class A %", 6.7, sg("social_class", "ClasseA"), 64)
add("class E %", 20.4, sg("social_class", "ClasseE"), 64)
add("pre-COVID %", 8.7, sg("covid_period", "pre_covid"), 64)
add("post-COVID %", 6.9, sg("covid_period", "post_covid"), 64)
# Os QM da Table S5 só fecham com modelo multinível com intercepto aleatório de ESTUDO
# (conferido em 19/10/2026: REML univariado, EE e DL não fecham). Esse modelo não deixa
# heterogeneidade residual do contraste entre estratos; a variante estudo/estrato a inclui.
qm <- function(dom, levs = NULL, nested = FALSE) {
  d <- subset(prim, outcome == "dementia" & domain == dom)
  if (!is.null(levs)) d <- subset(d, subgroup %in% levs)
  d$obs <- seq_len(nrow(d))
  rnd <- if (nested) ~ 1 | Study / obs else ~ 1 | Study
  rma.mv(yi, vi, mods = ~ factor(subgroup), random = rnd, data = d, method = "REML")
}
tab_s5 <- list(age = list(ages, 662.92), sex = list(NULL, 17.65), education = list(NULL, 320.87),
               social_class = list(NULL, 53.26), covid_period = list(NULL, 0.05))
for (dom in names(tab_s5)) {
  add(paste("QM", dom), tab_s5[[dom]][[2]], round(qm(dom, tab_s5[[dom]][[1]])$QM, 2), 0, "Table S5; rma.mv random=~1|Study")
  fn <- qm(dom, tab_s5[[dom]][[1]], nested = TRUE)
  add(paste("QM", dom, "with stratum residual (new)"), NA, round(fn$QM, 2), 0,
      sprintf("nao reportado; p = %.2g", fn$QMp))
}

res <- bind_rows(rows)
res$status <- ifelse(is.na(res$reported), "new",
                     ifelse(res$abs_diff <= 0.051 | (res$reported > 1000 & res$abs_diff <= 2), "reproduced", "DIFFERENT"))
write.csv(res, file.path(OUT_TAB, "replication_frequentist.csv"), row.names = FALSE)
print(res[, c("item", "reported", "reproduced", "status")], row.names = FALSE)
cat("\nMetarregressao sem Magalhaes (sensibilidade nova):\n"); print(coef(summary(mr_nomag)))
cat("\nEstudos por estrato etario:\n"); print(age_tab[, c("age", "k", "n_with_counts", "studies")], row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(OUT_DIAG, "sessioninfo_03_reproduce_frequentist.txt"))
