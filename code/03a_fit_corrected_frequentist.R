# Produces: output/tables/corrected_overall_by_set.csv, corrected_metaregression_diagnostic.csv, corrected_age_strata.csv, corrected_national_burden_2022.csv; output/intermediate/age_strata_corrected.csv
# Reanálise frequentista com as correções verificadas na fonte.
#
# NÃO é análise pré-especificada: são correções de erro de extração/elegibilidade e
# sensibilidades post hoc, cada uma isolando um problema documentado em
# data/raw/extraction/source_checks.csv.
#
# Conjuntos (cada um parte do anterior):
#   S0 original k = 11 (reproduzido em 02b_reproduce_frequentist.R)
#   S1 sem populações duplicadas: Tatsch 2006 (= amostra do Bottino 2008) e Nitrini 2005 (= coorte do Herrera 2002)
#   S2 sem Ramos-Cerqueira 2005 (busca de casos por agentes de saúde: 2% é piso, não prevalência)
#   S3 Bertola 2023 com a prevalência ponderada publicada, 5,8% (4,7-7,2), em vez da bruta 364/5249
#   S4 sem Magalhães 2008 (CAMCOG < 80 sem ajuste em população 74,7% analfabeta + demência "possível")
#   S5 S4 sem Correa Ribeiro 2013 (clientes de plano privado): só amostras da população geral
#   S4adj S4 com as prevalências corrigidas pelos próprios autores nos estudos em duas fases
#        (Bottino 12,9%; Lopes 12,5%). A variância dessas estimativas não está no resumo: cenário
#        com variância binomial do n observado inflada por efeito de desenho 1 e 2.
suppressPackageStartupMessages({ library(metafor); library(dplyr); library(MASS, exclude = "select") })
source(file.path("code", "helpers", "setup.R"))
set.seed(20261019)

prim <- read.csv(file.path(OCT_OUT, "study_level_data", "meta_ready_primary_community.csv"))
pop  <- read.csv(file.path(OCT_IBGE, "ibge_census_brazil_2022_agegroups.csv"))
logit <- function(p) log(p / (1 - p))
vi_from_ci <- function(lb, ub) ((logit(ub) - logit(lb)) / (2 * qnorm(0.975)))^2

ov <- subset(prim, outcome == "dementia" & domain == "overall")[, c("Study", "total", "cases", "yi", "vi", "DiagnosticGroup")]
ov$prev <- ov$cases / ov$total
set_bertola_weighted <- function(d) {
  i <- d$Study == "Bertola_2023"
  d$yi[i] <- logit(0.058); d$vi[i] <- vi_from_ci(0.047, 0.072); d$prev[i] <- 0.058; d
}
set_two_phase_adjusted <- function(d, deff) {
  adj <- c(Bottino_2008 = 0.129, Lopes_2012 = 0.125)
  for (s in names(adj)) {
    i <- d$Study == s; p <- adj[[s]]
    d$yi[i] <- logit(p); d$vi[i] <- deff / (d$total[i] * p * (1 - p)); d$prev[i] <- p
  }
  d
}
dup <- c("Tatsch_2006", "Nitrini_2005")
# S0 reproduz a primeira versão com a planilha original; de S1 em diante, totais conferidos no texto integral
# (Magalhães: 231 casos, soma de sexo e de escolaridade na Tabela 2; a planilha tinha 229, soma das faixas de idade).
vs <- read.csv(VERIFIED_STRATA)
ov_v <- ov
for (s in intersect(ov_v$Study, vs$study)) {
  r <- subset(vs, study == s & domain == "overall" & category == "all")
  i <- ov_v$Study == s
  if (nrow(r) == 1 && (r$n != ov_v$total[i] || r$cases != ov_v$cases[i])) {
    es <- escalc("PLO", xi = r$cases, ni = r$n)
    ov_v$total[i] <- r$n; ov_v$cases[i] <- r$cases; ov_v$yi[i] <- es$yi; ov_v$vi[i] <- es$vi; ov_v$prev[i] <- r$cases / r$n
  }
}
sets <- list(
  S0 = ov,
  S1 = subset(ov_v, !Study %in% dup),
  S2 = subset(ov_v, !Study %in% c(dup, "RamosCerqueira_2005")),
  S3 = set_bertola_weighted(subset(ov_v, !Study %in% c(dup, "RamosCerqueira_2005"))),
  S4 = set_bertola_weighted(subset(ov_v, !Study %in% c(dup, "RamosCerqueira_2005", "Magalhaes_2008")))
)
sets$S5 <- subset(sets$S4, Study != "Ribeiro_2013")
sets$S4adj_deff1 <- set_two_phase_adjusted(sets$S4, 1)
sets$S4adj_deff2 <- set_two_phase_adjusted(sets$S4, 2)

summ <- function(name, d) {
  m  <- rma(yi, vi, data = d, method = "REML")
  mh <- rma(yi, vi, data = d, method = "REML", test = "knha")
  p  <- predict(m, transf = transf.ilogit); ph <- predict(mh, transf = transf.ilogit)
  data.frame(set = name, k = m$k, studies = paste(d$Study, collapse = ";"),
             prev = p$pred, ci_lb = p$ci.lb, ci_ub = p$ci.ub,
             hk_ci_lb = ph$ci.lb, hk_ci_ub = ph$ci.ub, pi_lb = ph$pi.lb, pi_ub = ph$pi.ub,
             tau2 = m$tau2, I2 = m$I2)
}
overall <- bind_rows(Map(summ, names(sets), sets))
write.csv(overall, file.path(OUT_TAB, "corrected_overall_by_set.csv"), row.names = FALSE)
print(overall[, c("set", "k", "prev", "ci_lb", "ci_ub", "hk_ci_lb", "hk_ci_ub", "pi_lb", "pi_ub", "tau2")], digits = 3)

# Metarregressão do critério diagnóstico nos conjuntos corrigidos
mr <- bind_rows(lapply(c("S0", "S1", "S2", "S3", "S4"), function(s) {
  d <- sets[[s]]; d$other <- as.numeric(d$DiagnosticGroup != "DSM-IV")
  if (length(unique(d$other)) < 2) return(NULL)
  f <- rma(yi, vi, mods = ~ other, data = d, method = "REML", test = "knha")
  data.frame(set = s, k = f$k, n_other = sum(d$other), beta_other = coef(f)[2],
             ci_lb = f$ci.lb[2], ci_ub = f$ci.ub[2], p_knha = f$pval[2])
}))
write.csv(mr, file.path(OUT_TAB, "corrected_metaregression_diagnostic.csv"), row.names = FALSE)
print(mr, digits = 3)

## Prevalência por idade com o estrato aberto 85+ --------------------------------------------
# Fonte: extração conferida no texto integral (VERIFIED_STRATA). Regra: casos e n de cada faixa vêm da
# MESMA tabela do artigo. Nos estudos em duas fases essa tabela exclui os positivos não examinados
# (Bottino 1.477 de 1.563; Lopes 1.052 de 1.145), e o fator de correção de 03g_fit_author_best_estimates.R usa o mesmo total.
# Correções em relação à planilha de out/2026: Correa Ribeiro 85+ = 85-89 + >90 (63/129; a planilha
# perdeu a linha >90 e tinha 32/79); Lopes 60-64 = 265 (Tabela 3; a planilha usava 269, da Tabela 1).
keep_age <- c("Bottino_2008", "Cesar_2015", "Herrera_2002", "Lopes_2012", "Ribeiro_2013", "Scazufca_2008", "Bertola_2023")
va <- subset(vs, domain == "age" & study %in% keep_age)
age_long <- list()
w85 <- pop$population[pop$age_group == "85-89"]; w90 <- pop$population[pop$age_group == "90+"]
for (s in keep_age) {
  d <- subset(va, study == s)
  if (s == "Bertola_2023") {
    # Bertola só publica percentuais ponderados (Tabela 2): IC na escala logito
    for (g in c("60-64", "65-69", "70-74", "75-79", "80-84")) {
      r <- d[d$category == g, ]
      age_long[[length(age_long) + 1]] <- data.frame(study = s, band = g, yi = logit(r$prev_pct / 100),
        vi = vi_from_ci(r$prev_lb / 100, r$prev_ub / 100), cases = NA, n = NA, note = "ponderado (Tabela 2)")
    }
    # 85+ aproximado com pesos do Censo 2022 (85-89 : 90+)
    r85 <- d[d$category == "85-89", ]; r90 <- d[d$category == "90+", ]
    p85_ <- r85$prev_pct / 100; p90_ <- r90$prev_pct / 100
    v85 <- vi_from_ci(r85$prev_lb / 100, r85$prev_ub / 100); v90 <- vi_from_ci(r90$prev_lb / 100, r90$prev_ub / 100)
    p85 <- (w85 * p85_ + w90 * p90_) / (w85 + w90)
    se  <- sqrt((w85 / (w85 + w90))^2 * v85 * (p85_ * (1 - p85_))^2 +
                (w90 / (w85 + w90))^2 * v90 * (p90_ * (1 - p90_))^2)   # delta na escala de proporção
    age_long[[length(age_long) + 1]] <- data.frame(study = s, band = "85+", yi = logit(p85),
      vi = (se / (p85 * (1 - p85)))^2, cases = NA, n = NA, note = "85-89 e 90+ ponderados combinados com pesos do Censo")
    next
  }
  for (g in c("60-64", "65-69", "70-74", "75-79", "80-84")) {
    r <- d[d$category == g, ]
    if (nrow(r)) {
      es <- escalc("PLO", xi = r$cases, ni = r$n)
      age_long[[length(age_long) + 1]] <- data.frame(study = s, band = g, yi = es$yi, vi = es$vi, cases = r$cases, n = r$n,
        note = paste0(r$source, if (r$category_source != g) paste0("; faixa ", r$category_source) else ""))
    }
  }
  r85p <- d[d$category == "85+", ]; r8589 <- d[d$category %in% c("85-89", "90+"), ]
  stopifnot("85+ ambíguo" = xor(nrow(r85p) == 1, nrow(r8589) == 2))
  x <- if (nrow(r85p)) r85p$cases else sum(r8589$cases); n <- if (nrow(r85p)) r85p$n else sum(r8589$n)
  es <- escalc("PLO", xi = x, ni = n)
  age_long[[length(age_long) + 1]] <- data.frame(study = s, band = "85+", yi = es$yi, vi = es$vi, cases = x, n = n,
    note = if (nrow(r85p)) paste0(r85p$source, "; faixa ", r85p$category_source) else paste0(r8589$source[1], "; 85-89 + 90+ somados"))
}
age_long <- bind_rows(age_long); age_long$row <- seq_len(nrow(age_long))
# conferência: nos estudos com contagens, as faixas somam o total da tabela de idade
chk <- aggregate(cbind(cases, n) ~ study, data = subset(age_long, !is.na(n)), FUN = sum)
print(chk)
write.csv(age_long, file.path(DERIVED, "age_strata_corrected.csv"), row.names = FALSE)

bands <- c("60-64", "65-69", "70-74", "75-79", "80-84", "85+")
popb <- c(pop$population[match(bands[1:5], pop$age_group)], w85 + w90)

burden_from <- function(dat, label) {
  dat$band <- factor(dat$band, levels = bands)
  m <- rma.mv(yi, vi, mods = ~ 0 + band, random = ~ 1 | study / row, data = dat, method = "REML")
  b <- coef(m); V <- vcov(m)
  draws <- mvrnorm(20000, b, V)                        # incerteza conjunta das médias por estrato
  cases_mean <- plogis(draws) %*% popb                  # estimando: prevalência do estudo típico
  s2 <- m$sigma2                                        # [1] entre estudos, [2] entre estratos no estudo
  u <- rnorm(20000, 0, sqrt(s2[1])); v <- matrix(rnorm(20000 * 6, 0, sqrt(s2[2])), ncol = 6)
  cases_newpop <- plogis(draws + u + v) %*% popb        # previsão para uma nova população
  # média populacional (integra o efeito aleatório): Monte Carlo com 400 nós por draw
  z <- matrix(rnorm(400 * 2, 0, 1), ncol = 2)
  pa <- sapply(seq_len(6), function(j) {
    sapply(seq_len(2000), function(i) mean(plogis(draws[i, j] + z[, 1] * sqrt(s2[1]) + z[, 2] * sqrt(s2[2]))))
  })
  cases_popavg <- pa %*% popb
  q <- function(x) quantile(x, c(0.5, 0.025, 0.975))
  list(model = m,
       strata = data.frame(set = label, band = bands, k = as.vector(table(dat$band)), prev = plogis(b),
                           lb = plogis(b - 1.96 * sqrt(diag(V))), ub = plogis(b + 1.96 * sqrt(diag(V))), population = popb),
       total = data.frame(set = label, estimand = c("typical-study prevalence (RE = 0)", "new population (prediction)",
                                                    "population-averaged over study effects"),
                          rbind(q(cases_mean), q(cases_newpop), q(cases_popavg)),
                          sigma2_study = s2[1], sigma2_stratum = s2[2], check.names = FALSE))
}
b_all <- burden_from(age_long, "corrected age set (S4 studies)")
b_s5  <- burden_from(subset(age_long, study != "Ribeiro_2013"), "general-population only (S5)")
strata_orig <- bind_rows(b_all$strata, b_s5$strata)
tot <- bind_rows(b_all$total, b_s5$total)
names(tot)[3:5] <- c("cases_median", "cases_lb", "cases_ub")
write.csv(strata_orig, file.path(OUT_TAB, "corrected_age_strata.csv"), row.names = FALSE)
write.csv(tot, file.path(OUT_TAB, "corrected_national_burden_2022.csv"), row.names = FALSE)
print(strata_orig, digits = 3); print(tot, digits = 4)
writeLines(capture.output(sessionInfo()), file.path(OUT_DIAG, "sessioninfo_05_fit_corrected_frequentist.txt"))
