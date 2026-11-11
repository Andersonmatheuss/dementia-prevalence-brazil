# Produces: output/tables/reconciliation_bridge_2022.csv
# Dados da figura-ponte da reconciliação (2022, população 60+ do Censo): parte das taxas ponderadas do
# ELSI-Brazil e soma, passo a passo, o efeito de cada escolha de verificação de casos nos inquéritos em duas
# fases do conjunto principal. Os degraus vêm todos da mesma série de cenários (qba_two_phase_with_epifloripa,
# conjunto principal), para que somem; as estimativas de referência (estimativas dos autores, análise
# probabilística de viés, Delphi) entram como marcadores, fora da soma. Nenhum valor é a prevalência verdadeira.
suppressPackageStartupMessages({ library(dplyr) })
source(file.path("code", "helpers", "setup.R"))

tri <- read.csv(file.path(OUT_TAB, "triangulation_national_estimates.csv"))
qba <- read.csv(file.path(OUT_TAB, "qba_two_phase_with_epifloripa.csv")) |> filter(set == "primary")
est <- read.csv(file.path(OUT_TAB, "primary_burden_estimands.csv"))
pba <- read.csv(file.path(OUT_TAB, "pba_screening_correction.csv"))
get_tri <- function(pat) { v <- tri$cases[grepl(pat, tri$source, fixed = TRUE)]; stopifnot(length(v) == 1); v }
get_qba <- function(s) { r <- qba[qba$scenario == s, ]; stopifnot(nrow(r) == 1); r }

elsi <- get_tri("ELSI-Brazil 2015-16, weighted age-specific (Bertola 2023) x Census 2022")
steps <- data.frame(
  order = 1:6,
  type = "step",
  label = c("ELSI-Brazil weighted rates", "Meta-analysis, observed counts",
            "Unexamined screen-positives corrected (NPV = 1.00)", "NPV = 0.99", "NPV = 0.97", "NPV = 0.95"),
  cases_median = c(elsi, sapply(c("observed (no correction)", "two-phase corrected, NPV = 1.00",
                                  "two-phase corrected, NPV = 0.99", "two-phase corrected, NPV = 0.97",
                                  "two-phase corrected, NPV = 0.95"), function(s) get_qba(s)$cases_2022_median)),
  cases_lb = c(NA, sapply(c("observed (no correction)", "two-phase corrected, NPV = 1.00", "two-phase corrected, NPV = 0.99",
                            "two-phase corrected, NPV = 0.97", "two-phase corrected, NPV = 0.95"), function(s) get_qba(s)$cases_2022_lb)),
  cases_ub = c(NA, sapply(c("observed (no correction)", "two-phase corrected, NPV = 1.00", "two-phase corrected, NPV = 0.99",
                            "two-phase corrected, NPV = 0.97", "two-phase corrected, NPV = 0.95"), function(s) get_qba(s)$cases_2022_ub)),
  source = c("triangulation_national_estimates.csv", rep("qba_two_phase_with_epifloripa.csv (primary)", 5)))
steps$increment <- c(NA, diff(steps$cases_median))

ab <- est[est$version == "author_best" & est$estimand == "typical study", ]
pb <- pba[pba$variant == "partial_pooling" & pba$quantity == "total error (verification and random)", ]
refs <- data.frame(
  order = 7:9, type = "reference",
  label = c("Primary analysis (authors' population estimates)", "Probabilistic bias analysis", "Delphi consensus 8.5%"),
  cases_median = c(ab$cases_median, pb$cases_median, get_tri("Delphi consensus 8.5% applied to Census 2022 population 60+")),
  cases_lb = c(ab$cases_lb, pb$cases_lb, NA), cases_ub = c(ab$cases_ub, pb$cases_ub, NA),
  source = c("primary_burden_estimands.csv", "pba_screening_correction.csv", "triangulation_national_estimates.csv"),
  increment = NA)
out <- rbind(steps, refs)
print(out, digits = 3)
write.csv(out, file.path(OUT_TAB, "reconciliation_bridge_2022.csv"), row.names = FALSE)
