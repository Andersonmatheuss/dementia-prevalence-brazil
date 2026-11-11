# Produces: output/tables/qba_two_phase_documented_study_level.csv, qba_two_phase_documented_burden.csv
# QBA dos estudos em duas fases com as quantidades DOCUMENTADAS nas fontes (em vez de um cenário genérico).
#
# Fonte das quantidades: Chaimowicz F, Burdorf A. PLoS One 2015;10(7):e0131979 (PMC4488471, texto integral),
# que cita os artigos primários:
#   Herrera 2002: n = 1.656 rastreados; 234 positivos (14,1%); 220 avaliados; 118 casos -> VPP = 118/220.
#                 Nenhum negativo avaliado: VPN desconhecido.
#   Bottino 2008: n = 1.563; 250 positivos (16,0%); 164 avaliados; 107 casos -> VPP = 107/164.
#                 Estimativa corrigida dos autores = 12,9% -> VPN implícito ~ 0,97.
#   Lopes 2012:   n = 1.145; 217 positivos (18,9%); 130 avaliados; 68 casos -> VPP = 0,52;
#                 7,8% dos negativos avaliados: 3,3% com demência -> VPN = 0,967 (o 12,5% dos autores).
# Correção: p = s+ x VPP + (1 - s+) x (1 - VPN)   [fórmula citada por Chaimowicz; assume que os positivos não
# avaliados se parecem com os avaliados e ignora a não resposta da fase 1].
# Cenários de VPN: 1,00 (sem falso-negativo: só corrige os positivos não avaliados), 0,99, 0,97 (âncora
# empírica de Lopes e Bottino), 0,95. Correa Ribeiro (duas etapas) fica sem correção: faltam os dados de fase.
# Na carga, o fator de correção de cada estudo (corrigida/observada) escala suas contagens por idade:
# aproximação declarada (supõe correção proporcional igual em todas as idades).
suppressPackageStartupMessages({ library(metafor); library(dplyr); library(MASS, exclude = "select") })
source(file.path("code", "helpers", "setup.R"))
set.seed(20261019)

# Quantidades do rastreio conferidas no texto integral (VERIFIED_ESTIMATES). VPP = casos entre os
# positivos examinados / positivos examinados (no Bottino, 6 dos 107 casos e, no Lopes, 1 dos 68 vieram
# dos negativos examinados e não entram no VPP). Correa Ribeiro: duas fases sem exame da maioria dos
# negativos (Figura 1: 736 rastreados, 271 positivos, 221 examinados, 115 casos; os 683 analisados
# excluem os positivos não examinados).
est <- read.csv(VERIFIED_ESTIMATES)
two <- c("Herrera_2002", "Bottino_2008", "Lopes_2012", "Ribeiro_2013")
e <- est[match(two, est$study), ]
ph <- data.frame(study = two, n = ifelse(e$study == "Ribeiro_2013", 736, e$n_total), pos = e$screen_positive,
                 evaluated = e$positives_examined, cases = e$cases_in_positives, cases_total = e$cases_total,
                 n_observed_denominator = e$n_total)
ph$s_pos <- ph$pos / ph$n; ph$ppv <- ph$cases / ph$evaluated; ph$observed <- ph$cases_total / ph$n_observed_denominator
npv_grid <- c(1.00, 0.99, 0.97, 0.95)
adj <- bind_rows(lapply(npv_grid, function(npv) transform(ph, npv = npv, adjusted = s_pos * ppv + (1 - s_pos) * (1 - npv))))
adj$factor <- adj$adjusted / adj$observed
write.csv(adj, file.path(OUT_TAB, "qba_two_phase_documented_study_level.csv"), row.names = FALSE)
print(adj[, c("study", "npv", "observed", "adjusted", "factor")], digits = 3)
# checagem: com VPN 0,967 o Lopes deve dar ~12,5% (valor dos autores)
chk <- with(ph[ph$study == "Lopes_2012", ], s_pos * ppv + (1 - s_pos) * (1 - 0.967)); cat("Lopes, VPN 0,967:", round(100 * chk, 1), "%\n")

agec <- read.csv(file.path(DERIVED, "age_strata_corrected.csv"))
cen  <- read.csv(file.path(OCT_IBGE, "ibge_census_brazil_2022_agegroups.csv"))
bands <- c("60-64", "65-69", "70-74", "75-79", "80-84", "85+")
pop6 <- c(cen$population[match(bands[1:5], cen$age_group)], sum(cen$population[cen$age_group %in% c("85-89", "90+")]))
fit_age <- function(d) { d$band <- factor(d$band, levels = bands); rma.mv(yi, vi, mods = ~ 0 + band, random = ~ 1 | study / row, data = d, method = "REML") }

prim <- read.csv(file.path(OCT_OUT, "study_level_data", "meta_ready_primary_community.csv"))
s4 <- c("Bertola_2023", "Bottino_2008", "Cesar_2015", "Herrera_2002", "Lopes_2012", "Ribeiro_2013", "Scazufca_2008")
ov <- subset(prim, outcome == "dementia" & domain == "overall" & Study %in% s4)
ov$yi[ov$Study == "Bertola_2023"] <- qlogis(0.058); ov$vi[ov$Study == "Bertola_2023"] <- ((qlogis(0.072) - qlogis(0.047)) / 3.92)^2

res <- bind_rows(lapply(c(NA, npv_grid), function(npv) {
  d <- agec; o <- ov
  if (!is.na(npv)) for (s in ph$study) {
    f <- adj$factor[adj$study == s & adj$npv == npv]
    i <- d$study == s & !is.na(d$cases)
    p <- pmin(d$cases[i] / d$n[i] * f, 0.95); d$yi[i] <- qlogis(p); d$vi[i] <- 1 / (d$n[i] * p * (1 - p))
    po <- adj$adjusted[adj$study == s & adj$npv == npv]; j <- o$Study == s
    o$yi[j] <- qlogis(po); o$vi[j] <- 1 / (o$total[j] * po * (1 - po))
  }
  m <- fit_age(d); dr <- plogis(mvrnorm(20000, coef(m), vcov(m))); cs <- as.vector(dr %*% pop6)
  mo <- rma(yi, vi, data = o, method = "REML", test = "knha"); po <- predict(mo, transf = transf.ilogit)
  data.frame(scenario = if (is.na(npv)) "observed (no correction)" else sprintf("two-phase corrected, NPV = %.2f", npv),
             pooled_prev_S4 = po$pred, pooled_lb = po$ci.lb, pooled_ub = po$ci.ub,
             cases_2022_median = median(cs), cases_2022_lb = quantile(cs, .025), cases_2022_ub = quantile(cs, .975))
}))
write.csv(res, file.path(OUT_TAB, "qba_two_phase_documented_burden.csv"), row.names = FALSE)
print(res, digits = 4)
writeLines(capture.output(sessionInfo()), file.path(OUT_DIAG, "sessioninfo_09_qba_two_phase_documented.txt"))
