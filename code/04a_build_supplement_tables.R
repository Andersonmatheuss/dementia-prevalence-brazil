# Gera apenas o suplemento; não altera os CSVs nem reajusta modelos.
source(file.path("code", "helpers", "setup.R"))

read_source <- function(filename, directory = OUT_TAB, required = character()) {
  file <- file.path(directory, filename)
  if (!file.exists(file)) stop("CSV ausente: ", file, call. = FALSE)
  data <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8")
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop(filename, ": colunas ausentes: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (!nrow(data) || anyDuplicated(names(data))) {
    stop("CSV vazio ou com cabeçalho duplicado: ", filename, call. = FALSE)
  }
  data
}

fmt_number <- function(x, digits = 3L, trim = TRUE) {
  result <- formatC(x, format = "f", digits = digits, decimal.mark = ".", big.mark = "")
  if (trim && digits > 0L) result <- sub("[.]$", "", sub("0+$", "", result))
  result[is.na(x)] <- "NA"
  result
}

fmt_count <- function(x) {
  result <- formatC(x, format = "f", digits = 0L, big.mark = ",", decimal.mark = ".")
  result[is.na(x)] <- "NA"
  result
}

fmt_percent <- function(x, scale = 100) {
  result <- paste0(fmt_number(x * scale, 1L, trim = FALSE), "%")
  result[is.na(x)] <- "NA"
  result
}

fmt_p <- function(x) {
  result <- formatC(x, format = "g", digits = 3L, flag = "#", decimal.mark = ".")
  result[is.na(x)] <- "NA"
  result
}

# Mantém as qualificações textuais da extração.
fmt_percent_text <- function(x) {
  vapply(as.character(x), function(value) {
    if (is.na(value)) return("NA")
    value <- trimws(value)
    match <- regexpr("^[0-9]+(?:[.][0-9]+)?", value, perl = TRUE)
    if (match[1] < 0L) return(value)
    prefix <- regmatches(value, match)
    suffix <- substring(value, nchar(prefix) + 1L)
    suffix <- sub("^%", "", suffix)
    # Só o intervalo imediatamente após a prevalência é reformatado.
    interval <- regexec("^ *\\(([0-9.]+)[-–]([0-9.]+)\\)", suffix, perl = TRUE)
    parts <- regmatches(suffix, interval)[[1]]
    if (length(parts)) {
      bounds <- paste0(" (", fmt_number(as.numeric(parts[2]), 1L, FALSE), "–",
                       fmt_number(as.numeric(parts[3]), 1L, FALSE), ")")
      suffix <- paste0(bounds, substring(suffix, nchar(parts[1]) + 1L))
    }
    paste0(fmt_percent(as.numeric(prefix), scale = 1), suffix)
  }, character(1), USE.NAMES = FALSE)
}

fmt_count_text <- function(x) {
  vapply(as.character(x), function(value) {
    if (is.na(value)) return("NA")
    matches <- gregexpr("[0-9]+", value, perl = TRUE)
    tokens <- regmatches(value, matches)[[1]]
    if (length(tokens)) regmatches(value, matches) <- list(fmt_count(as.numeric(tokens)))
    value
  }, character(1), USE.NAMES = FALSE)
}

interval_spec <- function(lower, upper, label, kind = "percent") {
  list(lower = lower, upper = upper, label = label, kind = kind)
}

format_table <- function(data, percent = character(), intervals = list(),
                         counts = character(), pvalues = character()) {
  required <- c(percent, counts, pvalues,
                unlist(lapply(intervals, function(spec) c(spec$lower, spec$upper))))
  missing <- setdiff(required, names(data))
  if (length(missing)) stop("Colunas de formatação ausentes: ", paste(missing, collapse = ", "))
  count_names <- c("n", "n_analysed", "n_with_counts", "pos", "evaluated", "population", "pop60",
                   "cases", "reported_cases", "growth_size_component", "growth_age_structure_component",
                   "failed_fits", "divergent", counts)
  p_names <- c("p", "p_knha", "p_difference", "ppp", pvalues)
  output <- data
  for (column in names(data)) {
    value <- data[[column]]
    output[[column]] <- if (column %in% percent) {
      fmt_percent(value)
    } else if (grepl("_pct$", column) || column == "I2") {
      fmt_percent_text(value)
    } else if (column %in% p_names) {
      fmt_p(value)
    } else if (column %in% count_names || startsWith(column, "cases_")) {
      if (is.numeric(value)) fmt_count(value) else fmt_count_text(value)
    } else if (column %in% c("tau2", "sigma2_study", "sigma2_stratum")) {
      fmt_number(value, 3L, trim = FALSE)
    } else if (column == "max_rhat") {
      fmt_number(value, 4L, trim = FALSE)
    } else if (column %in% c("current_se", "expected_se_after", "variance")) {
      fmt_count(value)
    } else if (is.numeric(value)) {
      fmt_number(value)
    } else {
      as.character(value)
    }
  }
  # Junta limites sem alterar a escala da estimativa.
  for (spec in intervals) {
    lower <- data[[spec$lower]]
    upper <- data[[spec$upper]]
    formatter <- switch(spec$kind,
      percent = function(x) fmt_number(100 * x, 1L, trim = FALSE),
      count = fmt_count,
      number = function(x) fmt_number(x, 3L, trim = FALSE),
      stop("Escala de intervalo desconhecida: ", spec$kind))
    value <- paste0("(", formatter(lower), "–", formatter(upper), ")")
    value[is.na(lower) | is.na(upper)] <- "NA"
    output[[spec$lower]] <- value
    names(output)[names(output) == spec$lower] <- spec$label
    output[[spec$upper]] <- NULL
  }
  output
}

escape_cell <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- "NA"
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub("|", "&#124;", x, fixed = TRUE)
  trimws(gsub("[\r\n\t]+", " ", x))
}

markdown_table <- function(data) {
  stopifnot(ncol(data) > 0L, nrow(data) > 0L)
  row <- function(cells) paste0("| ", paste(escape_cell(cells), collapse = " | "), " |")
  c(row(names(data)), row(rep("---", ncol(data))),
    vapply(seq_len(nrow(data)), function(i) row(unlist(data[i, ], use.names = FALSE)), character(1)))
}

document <- c("# Supplementary tables", "",
  "Percentages are rounded to one decimal; percentage interval bounds are on the 0–100 scale. CI denotes confidence interval, CrI credible interval, and PI prediction interval. NA and source annotations are retained; unavailable intervals are not inferred.", "")
table_count <- 0L
add_heading <- function(number, title) {
  document <<- c(document, paste0("## eTable A", number, ". ", title), "")
}
add_table <- function(data, filename, directory = "output/tables", note = "", caption = NULL) {
  if (!is.null(caption)) document <<- c(document, paste0("**", caption, "**"), "")
  source <- if (nzchar(directory)) paste0(directory, "/", filename) else filename
  source_note <- paste0("Source: `", source, "`.", if (nzchar(note)) paste0(" ", note) else "")
  document <<- c(document, markdown_table(data), "", source_note, "")
  table_count <<- table_count + 1L
}

add_heading(1, "Study-level extraction and role in the analysis")
columns <- c("report_id", "independent_population_id", "region", "ascertainment_type", "n_analysed", "cases",
             "author_reported_prevalence_pct", "author_adjusted_or_weighted_pct", "role_reason")
data <- read_source("extraction.csv", dir_proc, columns)
# Papel na análise principal sem sobreposição, do analysis_sets.csv.
sets <- read.csv(file.path(dir_proc, "analysis_sets.csv"), stringsAsFactors = FALSE)
data$role_in_this_analysis <- sets$role[match(data$report_id, sets$report_id)]
data$role_in_this_analysis[is.na(data$role_in_this_analysis)] <- "not pooled"
# Texto de exibição: tira do role_reason as referências à versão anterior da revisão (o CSV fica como está).
clean_reason <- c(
  "; the review's own Table S16 lists it as excluded" = "",
  " for selection bias" = "",
  "review's own Table S16 lists it as excluded" = "",
  "meets the review's own screening-only exclusion criterion; Table S16 lists it as excluded" = "classification by screening scores alone",
  "eligible under the review's own criteria; published Oct 2024 before the search end date; " = "",
  "use weighted estimate; reclassify COVID period (data 2015-2016); authors' own burden 1,757,480 in 2019 not cited" = "weighted estimate used; data collected in 2015-2016",
  " (review rule, Table S16 category 1)" = "",
  " (scripts/36)" = "")
for (k in names(clean_reason)) data$role_reason <- gsub(k, clean_reason[[k]], data$role_reason, fixed = TRUE)
data$role_reason <- sub(";\\s*$", "", data$role_reason)
columns <- append(columns, "role_in_this_analysis", after = match("role_reason", columns) - 1L)
add_table(format_table(data[, columns]), "extraction.csv", directory = "data/processed",
          note = "role_in_this_analysis comes from `data/processed/analysis_sets.csv` (primary analysis without overlapping samples).")

add_heading(2, "Source verification and adjudication")
columns <- c("study_id", "pmid", "doi", "issue", "verdict")
data <- read_source("source_checks.csv", RAW_EXTRACTION, columns)
add_table(format_table(data[, columns]), "source_checks.csv", directory = "data/raw/extraction")

add_heading(3, "Replication of frequentist results")
columns <- c("item", "reported", "reproduced", "status")
data <- read_source("replication_frequentist.csv", required = columns)[, columns]
for (column in c("reported", "reproduced")) {
  data[[column]] <- vapply(seq_len(nrow(data)), function(i) {
    value <- data[[column]][i]
    item <- data$item[i]
    if (grepl("%", item, fixed = TRUE)) return(fmt_percent(value, scale = 1))
    if (grepl("tau2", item, fixed = TRUE)) return(fmt_number(value, 3L, trim = FALSE))
    if (grepl("(^| )p$", item)) return(fmt_p(value))
    if (grepl("burden|cases|component|(^| )[Nk]$", item)) return(fmt_count(value))
    fmt_number(value)
  }, character(1))
}
add_table(data, "replication_frequentist.csv")

add_heading(4, "Replication of Bayesian burden estimates")
data <- read_source("replication_bayesian.csv")
add_table(format_table(data, percent = "prev_median", intervals = list(
  interval_spec("cases_lb", "cases_ub", "95% posterior interval (cases)", "count"),
  interval_spec("reported_lb", "reported_ub", "Reported 95% interval (cases)", "count"))),
  "replication_bayesian.csv", note = "Conditional estimates have credible intervals; marginal new-study estimates include predictive heterogeneity.")

add_heading(5, "Pooled prevalence across report sets")
data <- read_source("corrected_overall_by_set.csv")
add_table(format_table(data, percent = "prev", intervals = list(
  interval_spec("ci_lb", "ci_ub", "95% CI (%)"),
  interval_spec("hk_ci_lb", "hk_ci_ub", "95% HK CI (%)"),
  interval_spec("pi_lb", "pi_ub", "95% PI (%)"))), "corrected_overall_by_set.csv",
  note = "HK denotes Hartung–Knapp; tau2 is on the logit scale.", caption = "Report sets S0 to S5")
data <- read_source("corrected_overall_with_epifloripa.csv")
add_table(format_table(data, percent = "prev", intervals = list(
  interval_spec("hk_lb", "hk_ub", "95% HK CI (%)"),
  interval_spec("pi_lb", "pi_ub", "95% PI (%)"))), "corrected_overall_with_epifloripa.csv",
  note = "tau2 is on the logit scale.", caption = "Sets including EpiFloripa")

add_heading(6, "Diagnostic ascertainment meta-regression")
data <- read_source("corrected_metaregression_diagnostic.csv")
add_table(format_table(data, intervals = list(interval_spec("ci_lb", "ci_ub", "95% CI (logit)", "number"))),
  "corrected_metaregression_diagnostic.csv", note = "beta_other and its confidence limits are on the logit scale.")

add_heading(7, "Age-specific prevalence under alternative analyses")
data <- read_source("age_specific_reproduced.csv")
add_table(format_table(data, percent = "est", intervals = list(
  interval_spec("lb", "ub", "95% CI (%)"), interval_spec("pi_lb", "pi_ub", "95% PI (%)"),
  interval_spec("cases_lb", "cases_ub", "95% CI (cases)", "count"))), "age_specific_reproduced.csv",
  note = "Original age bands retain separate 85–89 and 90+ groups.", caption = "Original analysis")
data <- read_source("corrected_age_strata.csv")
add_table(format_table(data, percent = "prev", intervals = list(interval_spec("lb", "ub", "95% CI (%)"))),
  "corrected_age_strata.csv", caption = "Report sets S4 and S5")
data <- read_source("corrected_age_strata_with_epifloripa.csv")
add_table(format_table(data, percent = "prev", intervals = list(interval_spec("lb", "ub", "95% CI (%)"))),
  "corrected_age_strata_with_epifloripa.csv", caption = "Analysis including EpiFloripa")
data <- read_source("author_best_age_strata.csv", required = c("deff", "band", "prev"))
data <- data[!is.na(data$deff) & data$deff == 1, , drop = FALSE]
add_table(format_table(data, percent = "prev"), "author_best_age_strata.csv",
  note = "Restricted to deff = 1; the CSV contains no interval bounds. The corrected analyses use an open 85+ group.",
  caption = "Author-best estimates, design effect 1")

add_heading(8, "National dementia burden in 2022 by report set")
data <- read_source("corrected_national_burden_2022.csv")
add_table(format_table(data, intervals = list(interval_spec("cases_lb", "cases_ub", "95% simulation interval (cases)", "count"))),
  "corrected_national_burden_2022.csv", note = "Limits are simulation percentiles; intervals for the new-population estimand are predictive.")

add_heading(9, "National burden using author-best estimates, 2022 to 2070")
data <- read_source("author_best_burden_projections.csv")
add_table(format_table(data, intervals = list(interval_spec("cases_lb", "cases_ub", "95% simulation interval (cases)", "count"))),
  "author_best_burden_projections.csv", note = "Both design-effect scenarios are retained; limits are simulation percentiles. 2022 uses the Census; later years the 2024 revision of the IBGE projections. Growth components: symmetric decomposition of the change since 2022 into population size and age structure (medians).")

add_heading(10, "National and regional burden under IBGE population projections")
data <- read_source("projections_2018_vs_2024.csv")
add_table(format_table(data, intervals = list(interval_spec("cases_corrected_lb", "cases_corrected_ub", "95% simulation interval (cases)", "count"))),
  "projections_2018_vs_2024.csv", caption = "Comparison of population projection revisions")
data <- read_source("author_best_projections_regions.csv")
add_table(format_table(data, intervals = list(interval_spec("cases_lb", "cases_ub", "95% simulation interval (cases)", "count"))),
  "author_best_projections_regions.csv", caption = "Regional expected counts under the national primary age schedule (2024 revision of the IBGE projections)",
  note = "Expected counts under one national age schedule; they carry no information about regional differences in risk.")

add_heading(11, "Documented two-phase ascertainment bias analyses")
data <- read_source("qba_two_phase_documented_study_level.csv")
add_table(format_table(data, percent = c("s_pos", "ppv", "observed", "npv", "adjusted")),
  "qba_two_phase_documented_study_level.csv", caption = "Study-level correction scenarios")
data <- read_source("qba_two_phase_with_epifloripa.csv")
add_table(format_table(data, percent = "pooled_prev", intervals = list(
  interval_spec("pooled_lb", "pooled_ub", "95% HK CI (%)"),
  interval_spec("cases_2022_lb", "cases_2022_ub", "95% simulation interval (cases)", "count"))),
  "qba_two_phase_with_epifloripa.csv", caption = "Pooled scenarios including EpiFloripa")

add_heading(12, "Sensitivity to ELSI algorithm misclassification")
data <- read_source("qba_elsi_algorithm_misclassification.csv")
add_table(format_table(data, percent = c("Se", "Sp", "corrected_prev_60plus")),
  "qba_elsi_algorithm_misclassification.csv", note = "NA results remain missing when the source flags a negative corrected stratum.")

add_heading(13, "Age-specific external comparison with weighted ELSI estimates")
data <- read_source("external_comparison_elsi_by_age.csv")
add_table(format_table(data, percent = c("pooled_without_elsi", "elsi_weighted"), intervals = list(
  interval_spec("ci_lb", "ci_ub", "95% CI (%)"), interval_spec("pi_lb", "pi_ub", "95% PI (%)"))),
  "external_comparison_elsi_by_age.csv")

add_heading(14, "Triangulation of national dementia estimates")
add_table(format_table(read_source("triangulation_national_estimates.csv")), "triangulation_national_estimates.csv")

add_heading(15, "Expected value of information and uncertainty components")
add_table(format_table(read_source("voi_expected_variance_reduction.csv")), "voi_expected_variance_reduction.csv",
  note = "Standard errors are expressed in cases.", caption = "Expected variance reduction by survey design")
add_table(format_table(read_source("voi_uncertainty_components.csv")), "voi_uncertainty_components.csv",
  note = "Variance is expressed in squared cases; missing shares are retained.", caption = "Components of uncertainty")

add_heading(16, "Bayesian models and diagnostics")
for (filename in c("bayesian_corrected_models.csv", "bayesian_leave_one_study_out.csv")) {
  data <- read_source(filename)
  add_table(format_table(data, intervals = list(interval_spec("cases_lb", "cases_ub", "95% CrI (cases)", "count"))),
    filename, caption = if (filename == "bayesian_corrected_models.csv") "Model and estimand sensitivity" else "Leave-one-study-out analysis")
}
data <- read_source("bayesian_posterior_predictive_checks.csv", OUT_DIAG, c("stat", "observed", "ppp"))
display <- format_table(data)
coverage <- grepl("coverage", data$stat, fixed = TRUE)
display$observed[coverage] <- fmt_percent(data$observed[coverage])
add_table(display, "bayesian_posterior_predictive_checks.csv", directory = "output/intermediate/diagnostics",
  note = "Predictive coverage is a percentage; mean yi and sd yi use the logit scale; ppp is a posterior predictive tail probability.",
  caption = "Posterior predictive checks")

add_heading(17, "Subgroup prevalence and moderator tests")
data <- read_source("corrected_subgroups.csv")
add_table(format_table(data, percent = "prev", intervals = list(interval_spec("lb", "ub", "95% CI (%)"))),
  "corrected_subgroups.csv", caption = "Subgroup estimates")
add_table(format_table(read_source("corrected_subgroup_tests.csv")), "corrected_subgroup_tests.csv",
  caption = "Moderator tests")

add_heading(18, "Full-text extraction: verification against the source")
columns <- c("report_id", "field", "current_value", "verified_value", "source", "status", "analysis_impact")
data <- read_source("fulltext_adjudication.csv", RAW_EXTRACTION, columns)
add_table(format_table(data[, columns]), "fulltext_adjudication.csv", directory = "data/raw/extraction")

add_heading(19, "Risk of bias (JBI checklist for prevalence studies)")
jbi_long <- read_source("fulltext_jbi.csv", RAW_EXTRACTION, c("report_id", "item", "final"))
jbi_long <- jbi_long[!is.na(jbi_long$final) & nzchar(jbi_long$final), ]
items <- paste0("Q", 1:9)
wide <- with(jbi_long, tapply(final, list(report_id, item), function(x) x[1]))[, items]
n_no <- rowSums(wide == "No", na.rm = TRUE)
n_unclear <- rowSums(wide == "Unclear", na.rm = TRUE)
jbi <- data.frame(Study = rownames(wide), wide, check.names = FALSE, row.names = NULL)
names(jbi)[names(jbi) %in% items] <- paste0("JBI_", items)
jbi$risk_of_bias <- ifelse(n_no == 0 & n_unclear <= 1, "low",
                           ifelse(n_no <= 2 & n_unclear <= 2, "moderate", "high"))
rule <- "Low: 0 No and <=1 Unclear; moderate: 1-2 No and <=2 Unclear; high: >=3 No or >2 Unclear."
columns <- c("Study", paste0("JBI_Q", 1:9), "risk_of_bias")
add_table(format_table(jbi[, columns]), "fulltext_jbi.csv", directory = "data/raw/extraction",
  note = paste0("One row per report assessed; reports excluded from synthesis were not assessed. Q1-Q9 are the JBI items. Classification rule: ", rule))

# Valida cada bloco antes e depois da gravação.
validate_markdown <- function(lines) {
  headings <- grep("^## eTable A[0-9]+[.] ", lines, value = TRUE)
  numbers <- as.integer(sub("^## eTable A([0-9]+)[.].*$", "\\1", headings))
  stopifnot("Numeração de eTables inválida" = identical(numbers, seq_len(19L)))
  indices <- which(startsWith(lines, "|"))
  blocks <- split(indices, cumsum(c(TRUE, diff(indices) != 1L)))
  stopifnot("Número inesperado de tabelas" = length(blocks) == table_count,
            "Número inesperado de fontes" = sum(startsWith(lines, "Source: `")) == table_count)
  for (block in blocks) {
    rows <- lines[block]
    pipe_counts <- nchar(rows) - nchar(gsub("|", "", rows, fixed = TRUE))
    stopifnot("Tabela markdown irregular" = length(rows) >= 3L &&
                all(pipe_counts == pipe_counts[1]) && all(endsWith(rows, "|")),
              "Separador markdown inválido" = grepl("^\\|(?: --- \\|)+$", rows[2], perl = TRUE))
  }
  invisible(length(blocks))
}

validate_markdown(document)
output_file <- file.path("supplementary", "supplement_tables.md")
writeLines(enc2utf8(document), output_file, useBytes = TRUE)
written <- readLines(output_file, encoding = "UTF-8", warn = FALSE)
stopifnot("Arquivo gravado difere do conteúdo validado" = identical(written, document))
validate_markdown(written)
cat(sprintf("Gerado %s: 19 eTables, %d tabelas markdown válidas.\n", output_file, table_count))
