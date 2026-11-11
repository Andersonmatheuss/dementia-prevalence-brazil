# 04 — Reimplementação do modelo bayesiano descrito no manuscrito de 18/10/2026.
#
# O script original de out/2026 não existe; o de mar/2026 (GitHub, 04_bayesian_hierarchical.R)
# usa outra especificação (k = 6, priores t de Student, critério diagnóstico como covariável e
# pós-estratificação no nível "CAMDEX"). Aqui segue-se o texto do Methods e da Table S8:
#   cases | trials(n) ~ s(age_mid, bs = "tp", k = 5) + (1 | study)
#   Normal(0, 5) para intercepto e efeitos fixos; meia-Normal(0, 2,5) para o DP de estudo;
#   prior padrão do brms para o spline; 4 cadeias x 4000 iterações (2000 de aquecimento);
#   adapt_delta = 0,99. Dados: contagens por faixa de 5 anos dos estudos da análise primária
#   (o Bertola só tem percentuais e não entra num binomial).
# Duas leituras de "marginal over the study-level random intercept" são testadas:
#   (a) condicional no efeito aleatório igual a zero (re_formula = NA, como no código de março);
#   (b) integrada sobre a distribuição do efeito aleatório (novo estudo, gaussiano).
# Semente: a do código de março (20260328); o manuscrito não informa semente.
suppressPackageStartupMessages({ library(brms); library(dplyr); library(posterior) })
source(file.path("code", "helpers", "setup.R"))
options(mc.cores = 2)
SEED <- 20260328

prim <- read.csv(file.path(OCT_OUT, "study_level_data", "meta_ready_primary_community.csv"))
pop  <- read.csv(file.path(OCT_IBGE, "ibge_census_brazil_2022_agegroups.csv"))
mid  <- c("60-64" = 62, "65-69" = 67, "70-74" = 72, "75-79" = 77, "80-84" = 82, "85-89" = 87, "90+" = 92)

dat <- prim |>
  filter(outcome == "dementia", domain == "age", subgroup %in% names(mid), !is.na(cases)) |>
  transmute(study = Study, age_band = subgroup, age_mid = unname(mid[subgroup]),
            cases = as.integer(cases), n = as.integer(total))
cat("linhas:", nrow(dat), "| estudos:", n_distinct(dat$study), "| N:", sum(dat$n), "\n")
write.csv(dat, file.path(DERIVED, "bayes_input_reproduction.csv"), row.names = FALSE)

pri <- c(prior(normal(0, 5), class = "Intercept"),
         prior(normal(0, 5), class = "b"),
         prior(normal(0, 2.5), class = "sd"))   # truncado em 0 pelo brms: meia-Normal

fit_one <- function(formula, file) {
  brm(formula, data = dat, family = binomial(), prior = pri, chains = 4, iter = 4000, warmup = 2000,
      cores = 2, seed = SEED, control = list(adapt_delta = 0.99), backend = "rstan",
      file = file.path(OUT_MOD, file), refresh = 0)
}
m_spl <- fit_one(bf(cases | trials(n) ~ s(age_mid, bs = "tp", k = 5) + (1 | study)), "repro_spline")
m_lin <- fit_one(bf(cases | trials(n) ~ age_mid + (1 | study)), "repro_linear")

grid <- data.frame(age_band = names(mid), age_mid = unname(mid), n = 1L, study = "new_study") |>
  left_join(transmute(pop, age_band = age_group, population), by = "age_band")

poststrat <- function(m, how) {
  ep <- if (how == "conditional_re0") {
    posterior_epred(m, newdata = grid, re_formula = NA)
  } else {
    set.seed(SEED)
    posterior_epred(m, newdata = grid, allow_new_levels = TRUE, sample_new_levels = "gaussian")
  }
  cases <- as.vector(ep %*% grid$population)
  list(age = data.frame(age_band = grid$age_band, median = apply(ep, 2, median),
                        lb = apply(ep, 2, quantile, 0.025), ub = apply(ep, 2, quantile, 0.975)),
       total = quantile(cases, c(0.5, 0.025, 0.975)))
}

diag_of <- function(m) {
  s <- summarise_draws(as_draws_df(m), "rhat", "ess_bulk", "ess_tail")
  s <- s[!grepl("^lp__|^lprior|^r_study|^zs_", s$variable), ]
  data.frame(max_rhat = max(s$rhat, na.rm = TRUE), min_ess_bulk = min(s$ess_bulk, na.rm = TRUE),
             min_ess_tail = min(s$ess_tail, na.rm = TRUE),
             divergent = sum(subset(nuts_params(m), Parameter == "divergent__")$Value))
}

out <- list()
for (mm in list(list("spline", m_spl), list("linear", m_lin))) {
  for (how in c("conditional_re0", "marginal_new_study")) {
    ps <- poststrat(mm[[2]], how)
    out[[length(out) + 1]] <- data.frame(model = mm[[1]], prediction = how,
      cases_median = ps$total[1], cases_lb = ps$total[2], cases_ub = ps$total[3],
      prev_median = ps$total[1] / sum(grid$population), diag_of(mm[[2]]))
    write.csv(ps$age, file.path(OUT_TAB, sprintf("bayes_repro_age_%s_%s.csv", mm[[1]], how)), row.names = FALSE)
  }
}
res <- bind_rows(out)
reported <- data.frame(model = c("spline", "linear"), reported_cases = c(2419715, 2342112),
                       reported_lb = c(1233646, 1281791), reported_ub = c(4316011, 4197731))
res <- left_join(res, reported, by = "model")
write.csv(res, file.path(OUT_TAB, "replication_bayesian.csv"), row.names = FALSE)
print(res)
print(read.csv(file.path(OUT_TAB, "bayes_repro_age_spline_conditional_re0.csv")))
writeLines(capture.output(sessionInfo()), file.path(OUT_DIAG, "sessioninfo_04_reproduce_bayesian.txt"))
