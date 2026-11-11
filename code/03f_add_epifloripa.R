# Produces: output/tables/corrected_overall_with_epifloripa.csv, corrected_age_strata_with_epifloripa.csv, qba_two_phase_with_epifloripa.csv; output/intermediate/age_strata_corrected_with_epifloripa.csv; data/processed/analysis_sets.csv
# Inclusão do EpiFloripa 2024 (estudo elegível ausente do manuscrito) no conjunto corrigido.
#
# Fonte: Lopes MA, Stoppe Junior A, Correa Neto Y, D'Orsi E. Aging Med (Milton) 2024;7(5):571-9.
# doi:10.1002/agm2.12362; PMID 39507227; PMC11535170 (texto integral).
#   Coorte EpiFloripa (>= 60 desde 2009-10), amostra por conglomerados; rastreio na 2a onda, diagnóstico
#   2017-2018 (NIA, consenso de dois médicos). 1.184 avaliados de 1.488 elegíveis; 243 positivos no rastreio,
#   102 (41,9%) avaliados; 16 negativos avaliados; VPP 46%, VPN 100%. 47 casos.
#   Bruta 4,5% (IC95% 3,24-5,76; = 47/1.043 classificados); estimada 9,2% (7,45-10,95).
#   Tabela 1 (casos/classificados): 65-69 0/214; 70-74 5/294; 75-79 12/222; 80-84 8/195; >=85 22/118.
# O estudo cumpre os critérios do próprio protocolo (população geral, critério reconhecido, confirmação clínica).
# Entra no conjunto corrigido com a prevalência OBSERVADA, como Bottino e Lopes;
# a correção pela acurácia do rastreio entra no QBA (fator documentado = s+ x VPP / observada).
# Análise principal sem sobreposição (PRIMARY_POPULATIONS, code/helpers/setup.R): este script grava os
# papéis em data/processed/analysis_sets.csv e roda o QBA no conjunto principal e nas oito populações.
suppressPackageStartupMessages({ library(metafor); library(dplyr); library(MASS, exclude = "select") })
source(file.path("code", "helpers", "setup.R"))
set.seed(20261019)

prim <- read.csv(file.path(OCT_OUT, "study_level_data", "meta_ready_primary_community.csv"))
agec <- read.csv(file.path(DERIVED, "age_strata_corrected.csv"))
cen  <- read.csv(file.path(OCT_IBGE, "ibge_census_brazil_2022_agegroups.csv"))
bands <- c("60-64", "65-69", "70-74", "75-79", "80-84", "85+")
pop6 <- c(cen$population[match(bands[1:5], cen$age_group)], sum(cen$population[cen$age_group %in% c("85-89", "90+")]))
fit_age <- function(d) { d$band <- factor(d$band, levels = bands); rma.mv(yi, vi, mods = ~ 0 + band, random = ~ 1 | study / row, data = d, method = "REML") }
burden <- function(m) { cs <- as.vector(plogis(mvrnorm(20000, coef(m), vcov(m))) %*% pop6); quantile(cs, c(.5, .025, .975)) }

# --- EpiFloripa: global e por idade (contagens observadas; correção de continuidade só no zero) ---
# faixas do EpiFloripa da extração conferida (Tabela 1, colunas Total: sem + com demência)
vs <- read.csv(VERIFIED_STRATA)
epi <- subset(vs, study == "EpiFloripa_2024" & domain == "age")
epi <- epi[match(bands[2:6], epi$category), ]
stopifnot(!anyNA(epi$n))
epi_x <- epi$cases; epi_n <- epi$n
es <- escalc("PLO", xi = epi_x, ni = epi_n, add = 0.5, to = "only0")
epi_age <- data.frame(study = "EpiFloripa_2024", band = bands[2:6], yi = es$yi, vi = es$vi, cases = epi_x, n = epi_n,
                      note = "Lopes 2024 Table 1; 65-69 com zero casos (correcao 0,5)")
agec7 <- bind_rows(agec, epi_age); agec7$row <- seq_len(nrow(agec7))
write.csv(agec7, file.path(DERIVED, "age_strata_corrected_with_epifloripa.csv"), row.names = FALSE)
epi_ov <- escalc("PLO", xi = 47, ni = 1043)

s4 <- c("Bertola_2023", "Bottino_2008", "Cesar_2015", "Herrera_2002", "Lopes_2012", "Ribeiro_2013", "Scazufca_2008")
ov <- subset(prim, outcome == "dementia" & domain == "overall" & Study %in% s4)[, c("Study", "total", "cases", "yi", "vi")]
ov$yi[ov$Study == "Bertola_2023"] <- qlogis(0.058); ov$vi[ov$Study == "Bertola_2023"] <- ((qlogis(0.072) - qlogis(0.047)) / 3.92)^2
ov6 <- rbind(ov, data.frame(Study = "EpiFloripa_2024", total = 1043, cases = 47, yi = epi_ov$yi, vi = epi_ov$vi))

pool <- function(d, label) {
  m <- rma(yi, vi, data = d, method = "REML", test = "knha"); p <- predict(m, transf = transf.ilogit)
  data.frame(set = label, k = m$k, prev = p$pred, hk_lb = p$ci.lb, hk_ub = p$ci.ub, pi_lb = p$pi.lb, pi_ub = p$pi.ub, tau2 = m$tau2, I2 = m$I2)
}
overall <- bind_rows(pool(ov, "S4 (7 studies)"), pool(ov6, "S6 = S4 + EpiFloripa (8 studies)"),
                     pool(subset(ov6, Study != "Ribeiro_2013"), "S7 = S6 without health-plan sample (7 studies)"),
                     pool(subset(ov6, Study %in% PRIMARY_POPULATIONS), "Primary, no overlap (observed counts)"),
                     bind_rows(lapply(names(SENSITIVITY_SETS), function(s)
                       pool(subset(ov6, Study %in% SENSITIVITY_SETS[[s]]), paste0("Sensitivity: ", s, " (observed counts)")))))
print(overall, digits = 3)

# Papel de cada população no conjunto sem sobreposição (lido por 04a_build_supplement_tables.R e 06a_check_reproduction.R)
sets <- data.frame(study = POPULATIONS_ALL, report_id = unname(STUDY_REPORT), role = unname(STUDY_ROLE),
                   in_primary = POPULATIONS_ALL %in% PRIMARY_POPULATIONS)
for (s in names(SENSITIVITY_SETS)) sets[[paste0("in_", gsub("[^a-z]+", "_", tolower(s)))]] <- POPULATIONS_ALL %in% SENSITIVITY_SETS[[s]]
write.csv(sets, file.path(dir_proc, "analysis_sets.csv"), row.names = FALSE)

# --- QBA documentado incluindo o EpiFloripa (VPN = 1 publicado; s+ = 243/1184; VPP = 47/102) ---
# Quantidades do rastreio conferidas (VERIFIED_ESTIMATES); VPP só com os casos entre os positivos
# examinados. n = amostra da primeira fase (Correa Ribeiro: 736 rastreados, Figura 1; EpiFloripa: 1.184
# avaliados); o denominador observado é o de cada artigo (Ribeiro 683 e EpiFloripa 1.043 excluem os
# positivos não examinados).
est <- read.csv(VERIFIED_ESTIMATES)
two <- c("Herrera_2002", "Bottino_2008", "Lopes_2012", "Ribeiro_2013", "EpiFloripa_2024")
e <- est[match(two, est$study), ]
ph <- data.frame(study = two, n = c(1656, 1563, 1145, 736, 1184), pos = e$screen_positive, evaluated = e$positives_examined,
                 cases = e$cases_in_positives, cases_total = e$cases_total, n_observed_denominator = e$n_total)
stopifnot(identical(ph$n[1:3], as.numeric(e$n_total[1:3])))
ph$s_pos <- ph$pos / ph$n; ph$ppv <- ph$cases / ph$evaluated; ph$observed <- ph$cases_total / ph$n_observed_denominator
# Conjunto principal (sem sobreposição: Herrera, Lopes e EpiFloripa são os de duas fases) e, como
# sensibilidade, as oito populações.
qba_set <- function(ids, set_label) bind_rows(lapply(c(NA, 1, 0.99, 0.97, 0.95), function(npv) {
  d <- subset(agec7, study %in% ids); o <- subset(ov6, Study %in% ids)
  if (!is.na(npv)) for (s in intersect(ph$study, ids)) {
    r <- ph[ph$study == s, ]; npv_s <- if (s == "EpiFloripa_2024") max(npv, 0.97) else npv   # EpiFloripa: VPN publicado = 1
    p_adj <- r$s_pos * r$ppv + (1 - r$s_pos) * (1 - npv_s); f <- p_adj / r$observed
    i <- d$study == s & !is.na(d$cases); p <- pmin(pmax(d$cases[i], 0.5) / d$n[i] * f, 0.95)
    d$yi[i] <- qlogis(p); d$vi[i] <- 1 / (d$n[i] * p * (1 - p))
    j <- o$Study == s; o$yi[j] <- qlogis(p_adj); o$vi[j] <- 1 / (o$total[j] * p_adj * (1 - p_adj))
  }
  b <- burden(fit_age(d)); po <- pool(o, "")
  data.frame(set = set_label, k = po$k,
             scenario = if (is.na(npv)) "observed (no correction)" else sprintf("two-phase corrected, NPV = %.2f", npv),
             pooled_prev = po$prev, pooled_lb = po$hk_lb, pooled_ub = po$hk_ub,
             cases_2022_median = b[1], cases_2022_lb = b[2], cases_2022_ub = b[3])
}))
scen <- bind_rows(qba_set(PRIMARY_POPULATIONS, "primary"), qba_set(POPULATIONS_ALL, "all eight populations"))
age_fit <- fit_age(agec7)
strata <- data.frame(band = bands, prev = plogis(coef(age_fit)), lb = plogis(coef(age_fit) - 1.96 * sqrt(diag(vcov(age_fit)))),
                     ub = plogis(coef(age_fit) + 1.96 * sqrt(diag(vcov(age_fit)))))
print(strata, digits = 3); print(scen, digits = 4)
write.csv(overall, file.path(OUT_TAB, "corrected_overall_with_epifloripa.csv"), row.names = FALSE)
write.csv(strata, file.path(OUT_TAB, "corrected_age_strata_with_epifloripa.csv"), row.names = FALSE)
write.csv(scen, file.path(OUT_TAB, "qba_two_phase_with_epifloripa.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(OUT_DIAG, "sessioninfo_10_add_epifloripa.txt"))
