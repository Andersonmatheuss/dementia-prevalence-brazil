# Produces: output/tables/corrected_subgroups.csv, corrected_subgroup_tests.csv
# Subgrupos (sexo, escolaridade, classe): agrupamentos da primeira versão contra a extração conferida.
#
# "original k=11": os mesmos agrupamentos da primeira versão, com a planilha de out/2026 (inclui
# Magalhães, Tatsch, Nitrini e Ramos-Cerqueira; rótulos de escolaridade não harmonizados).
# "verified (primary, no overlap)": as populações da análise principal sem sobreposição (PRIMARY_POPULATIONS,
# code/helpers/setup.R), com os estratos conferidos no
# texto integral (VERIFIED_STRATA). Escolaridade em dois níveis, nenhuma (analfabeto/0 anos) contra
# alguma, a única harmonização segura entre estudos com faixas diferentes (1-3, 1-4, 1-5...); classe
# só com ABIPEME (Lopes e César no conjunto principal; Bottino em sensibilidade). A renda do Correa Ribeiro não é
# classe e fica de fora.
# Teste de diferença entre estratos: modelo multinível com estudo e estrato-no-estudo.
suppressPackageStartupMessages({ library(metafor); library(dplyr) })
source(file.path("code", "helpers", "setup.R"))

prim <- read.csv(file.path(OCT_OUT, "study_level_data", "meta_ready_primary_community.csv"))
vs <- read.csv(VERIFIED_STRATA)
primary_set <- PRIMARY_POPULATIONS

# primeira versão, com os rótulos de origem
orig <- subset(prim, outcome == "dementia" & domain %in% c("sex", "education", "social_class"))
orig$domain[orig$domain == "social_class"] <- "class"
orig <- orig[, c("Study", "domain", "subgroup", "yi", "vi")]

# extração conferida: soma as faixas de escolaridade "alguma"; contagens com correção de 0,5 só no zero
ver <- subset(vs, study %in% primary_set & domain %in% c("sex", "education", "class") & !is.na(n))
ver <- aggregate(cbind(n, cases) ~ study + domain + category, data = ver, FUN = sum)
es <- escalc("PLO", xi = ver$cases, ni = ver$n, add = 0.5, to = "only0")
ver <- data.frame(Study = ver$study, domain = ver$domain, subgroup = ver$category, yi = es$yi, vi = es$vi,
                  n = ver$n, cases = ver$cases)

one <- function(d, set) {
  bind_rows(lapply(split(d, d$subgroup), function(x) {
    if (nrow(x) < 2) return(data.frame(set = set, subgroup = x$subgroup[1], k = nrow(x), studies = paste(x$Study, collapse = ";"),
                                        prev = plogis(x$yi[1]), lb = NA, ub = NA))
    m <- rma(yi, vi, data = x, method = "REML"); p <- predict(m, transf = transf.ilogit)
    data.frame(set = set, subgroup = x$subgroup[1], k = m$k, studies = paste(x$Study, collapse = ";"), prev = p$pred, lb = p$ci.lb, ub = p$ci.ub)
  }))
}
out <- list(); tests <- list()
for (dm in c("sex", "education", "class")) {
  sets <- list("original k=11" = subset(orig, domain == dm), "verified (primary, no overlap)" = subset(ver, domain == dm))
  out[[dm]] <- bind_rows(Map(one, sets, names(sets))) |> mutate(domain = dm)
  for (nm in names(sets)) {
    d <- sets[[nm]]; d$obs <- seq_len(nrow(d))
    f <- rma.mv(yi, vi, mods = ~ factor(subgroup), random = ~ 1 | Study / obs, data = d, method = "REML")
    tests[[length(tests) + 1]] <- data.frame(domain = dm, set = nm, k_studies = length(unique(d$Study)), QM = f$QM, df = f$QMdf[1], p = f$QMp)
  }
}
res <- bind_rows(out); tst <- bind_rows(tests)
write.csv(res, file.path(OUT_TAB, "corrected_subgroups.csv"), row.names = FALSE)
write.csv(tst, file.path(OUT_TAB, "corrected_subgroup_tests.csv"), row.names = FALSE)
print(res[, c("domain", "set", "subgroup", "k", "prev", "lb", "ub")], digits = 3, row.names = FALSE)
print(tst, digits = 3)
# Peso do Magalhães no estrato "sem escolaridade" da primeira versão
il <- subset(prim, outcome == "dementia" & domain == "education" & subgroup == "Illiterate")
print(il[, c("Study", "total", "cases")])
writeLines(capture.output(sessionInfo()), file.path(OUT_DIAG, "sessioninfo_07_fit_corrected_subgroups.txt"))
