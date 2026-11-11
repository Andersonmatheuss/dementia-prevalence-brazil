# Produces: output/tables/voi_expected_variance_reduction.csv, voi_uncertainty_components.csv
# Valor de nova informação como REDUÇÃO ESPERADA DE VARIÂNCIA (pré-posterior), não EVSI.
#
# EVSI formal exige um modelo de decisão com função de perda (ex.: dimensionar serviços); não há
# decisão especificada, então EVSI não é definido aqui. O que se calcula: quanto a incerteza da
# carga nacional 2022 (estimando = prevalência do estudo típico, modelo multinível de 03a_fit_corrected_frequentist.R)
# cairia, em média, com cada desenho de estudo novo, simulando dados do próprio modelo ajustado.
# E, à parte, quanto da incerteza total vem da verificação diagnóstica (cenários de 03d_compare_external_and_qba.R),
# que nenhum novo inquérito com o mesmo instrumento reduz; só um estudo de validação reduz.
suppressPackageStartupMessages({ library(metafor); library(dplyr) })
source(file.path("code", "helpers", "setup.R"))
set.seed(20261019)
NSIM <- 200

agec <- read.csv(file.path(DERIVED, "age_strata_corrected.csv"))
cen  <- read.csv(file.path(OCT_IBGE, "ibge_census_brazil_2022_agegroups.csv"))
bands <- c("60-64", "65-69", "70-74", "75-79", "80-84", "85+")
pop6 <- c(cen$population[match(bands[1:5], cen$age_group)], sum(cen$population[cen$age_group %in% c("85-89", "90+")]))
fit_age <- function(d) { d$band <- factor(d$band, levels = bands); rma.mv(yi, vi, mods = ~ 0 + band, random = ~ 1 | study / row, data = d, method = "REML") }
se_burden <- function(m) {   # delta: d(casos)/d(beta_j) = pop_j * p_j * (1 - p_j)
  p <- plogis(coef(m)); g <- pop6 * p * (1 - p); sqrt(as.numeric(t(g) %*% vcov(m) %*% g))
}
m0 <- fit_age(agec); se0 <- se_burden(m0); s2 <- m0$sigma2; b0 <- coef(m0)
cat("SE atual da carga (delta):", round(se0), "\n")

age_share_census <- pop6 / sum(pop6)
designs <- list(
  "1 survey, n = 1,000, census age mix"           = list(n = 1000, w = age_share_census, k = 1),
  "1 survey, n = 5,000, census age mix"           = list(n = 5000, w = age_share_census, k = 1),
  "1 survey, n = 1,000, 85+ oversampled to 30%"   = list(n = 1000, w = c(age_share_census[1:5] / sum(age_share_census[1:5]) * 0.7, 0.3), k = 1),
  "3 surveys, n = 1,000 each, census age mix"     = list(n = 1000, w = age_share_census, k = 3),
  "5 surveys, n = 1,000 each, census age mix"     = list(n = 1000, w = age_share_census, k = 5)
)
voi <- bind_rows(lapply(names(designs), function(nm) {
  ds <- designs[[nm]]
  se_new <- replicate(NSIM, {
    add <- bind_rows(lapply(seq_len(ds$k), function(j) {
      u <- rnorm(1, 0, sqrt(s2[1])); v <- rnorm(6, 0, sqrt(s2[2]))
      n <- pmax(round(ds$n * ds$w), 5); x <- rbinom(6, n, plogis(b0 + u + v))
      es <- escalc("PLO", xi = x, ni = n, add = 0.5, to = "only0")
      data.frame(study = paste0("new", j), band = bands, yi = es$yi, vi = es$vi)
    }))
    d <- bind_rows(agec[, c("study", "band", "yi", "vi")], add); d$row <- seq_len(nrow(d))
    tryCatch(se_burden(fit_age(d)), error = function(e) NA)
  })
  data.frame(design = nm, current_se = se0, expected_se_after = mean(se_new, na.rm = TRUE),
             expected_variance_reduction_pct = 100 * (1 - mean(se_new^2, na.rm = TRUE) / se0^2),
             failed_fits = sum(is.na(se_new)))
}))
print(voi, digits = 4)

# Parcela da incerteza devida à verificação diagnóstica (cenário de sensibilidade do rastreio 0,5-1,0,
# tratado como uniforme): variância entre cenários vs variância amostral/heterogeneidade
qa <- read.csv(file.path(OUT_TAB, "qba_two_phase_screen_sensitivity.csv"))
fit_lin <- approxfun(qa$screen_sensitivity, qa$cases_median)
var_ascert <- var(fit_lin(runif(1e5, 0.5, 1)))
share <- data.frame(component = c("sampling + between-study heterogeneity (delta SE^2)", "screen-sensitivity scenarios Se ~ U(0.5, 1)"),
                    variance = c(se0^2, var_ascert))
share$share_pct <- 100 * share$variance / sum(share$variance)
print(share, digits = 4)
# Estudo de validação que fixa a sensibilidade do rastreio em +/- 0,05 em torno de 0,75
var_after_validation <- var(fit_lin(runif(1e5, 0.70, 0.80)))
share_val <- data.frame(component = "ascertainment variance after a validation study (Se ~ U(0.70, 0.80))",
                        variance = var_after_validation, share_pct = NA)
out <- list(voi = voi, share = bind_rows(share, share_val))
write.csv(voi, file.path(OUT_TAB, "voi_expected_variance_reduction.csv"), row.names = FALSE)
write.csv(out$share, file.path(OUT_TAB, "voi_uncertainty_components.csv"), row.names = FALSE)
print(out$share, digits = 4)
writeLines(capture.output(sessionInfo()), file.path(OUT_DIAG, "sessioninfo_12_value_of_information.txt"))
