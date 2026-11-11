# =============================================================================
# Script 01: Data Harmonization and Preparation
# =============================================================================
#
# Title:   Prevalence of Dementia in Brazil: A Systematic Review and Meta-Analysis
# Journal: [Submission target]
# Authors: [Author list]
#
# Purpose:
#   Loads and harmonizes study-level data from 15 Brazilian epidemiological
#   studies (2002–2025), converts heterogeneous reporting formats (raw counts
#   and published prevalence with 95% CI) into a unified meta-analysis ready
#   dataset (meta_ready), and defines shared label maps and helper functions
#   used by all downstream scripts.
#
# Inputs:
#   data/raw/PrevDem.xlsx  — (data embedded directly in this script)
#
# Outputs (objects created in R environment):
#   meta_ready            — harmonized data frame with effect sizes (yi, vi)
#   dementia_overall_meta — subset for all-cause dementia, overall domain
#   lancet_*              — display label maps (Lancet/JAMA editorial style)
#   lbl()                 — label lookup helper
#   run_meta()            — random-effects meta-analysis wrapper
#   extract_meta_summary()— model summary extractor
#
# Requirements:
#   R >= 4.3.0
#   Packages: dplyr, tidyr, purrr, stringr, tibble, metafor
#
# Execution: source("code/01_data_harmonization.R")
#            OR called automatically by run_all.R
# =============================================================================

# -- Packages ------------------------------------------------------------------
packages <- c("dplyr", "tidyr", "purrr", "stringr", "tibble", "metafor")
new_pkgs <- packages[!packages %in% installed.packages()[, "Package"]]
if (length(new_pkgs)) install.packages(new_pkgs)

library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(tibble)
library(metafor)

options(stringsAsFactors = FALSE)

# =============================================================================
# 1. Study-level data
# =============================================================================
# 15 Brazilian community- and clinic-based studies on dementia prevalence.
# Data are embedded here to ensure full reproducibility without external
# dependency on the Excel workbook. The source Excel file (PrevDem.xlsx)
# is provided in data/raw/ for audit purposes.

study_names <- c(
  "Bertola_2023","Cesar_2015","FerreiraFilho_2020","Jacinto_2011","Lopes_2012",
  "Ribeiro_2013","Tatsch_2006","Scazufca_2008","RamosCerqueira_2005","Magalhaes_2008",
  "Bottino_2008","Herrera_2002","Nitrini_2005","Bezerra_2024","Luna_2025"
)

df <- data.frame(
  Study = study_names,
  Year = c(2023,2015,2020,2011,2012,2013,2006,2008,2005,2008,2008,2002,2005,2024,2024),
  Design = c("CS","CS","CS","CS","CS","CS","CS","CS","CS","CS","CS","CS","Prospective cohort","CS","Retrospective analysis of medical records"),
  Region = c("Brazil","Sudeste","Sul","Sudeste","Sudeste","Sudeste","Sudeste","Sudeste","Sudeste","Nordeste","Sudeste","Sudeste","Sudeste","Norte","Nordeste"),
  Coleta = c("Domiciliar","Domiciliar","Atenção Primária","Ambulatório hospitalar","Domiciliar","Domiciliar","Domiciliar","Domiciliar","Domiciliar","Domiciliar","Domiciliar","Domiciliar","Domiciliar","Aldeia Indigena (Domiciliar)","Ambulatório hospitalar"),
  Zona = c("U/R","U/R","U","U","U","U","U","U","U","R","U","U","U","U","U"),
  DiagnosticCriteria = c("NA","NIA-AA","NA","DSM IV","DSM IV","DSM IV","DSM IV","DSM IV","DSM IV","CAMDEX","DSM IV","DSM IV","DSM IV","NA","NIA-AA"),

  N_total    = c("5249","630","299","248","1,145","683","1,563","2072","2,222","466","1563","1656","1393","141","145"),
  N_dementia = c("364","110","38","21","68","115","64","105","45","229","107","118","113","18","65"),
  N_CIND     = c("422","135","104","22","16","NA","25","NA","NA","NA","25","NA","NA","16","12"),

  White     = c("NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","1089","NA","NA","NA","24"),
  dem_white = c("NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","67","NA","NA","NA","NA"),
  Brown     = c("NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","222","NA","NA","NA","107"),
  dem_Brown = c("NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","22","NA","NA","NA","NA"),
  Black     = c("NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","109","NA","NA","NA","2"),
  dem_Black = c("NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","13","NA","NA","NA","NA"),
  Asian     = c("NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","53","NA","NA","NA","1"),
  dem_Asian = c("NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","3","NA","NA","NA","NA"),

  total60_64 = c("NA","152","NA","NA","269","NA","NA","NA","NA","NA","368","NA","NA","NA","NA"),
  dem60_64   = c("3.2 (2.1–4.7)","8","NA","NA","5","NA","NA","NA","NA","NA","9","NA","NA","NA","NA"),

  total65_69 = c("NA","152","NA","NA","267","79","NA","889","836","NA","320","614","NA","NA","NA"),
  dem65_69   = c("2.7 (1.7–4.3)","16","NA","NA","6","3","NA","20","1","NA","13","10","NA","NA","NA"),

  total60_69 = c("NA","304","164","NA","536","79","NA","889","836","225","688","614","NA","42","NA"),
  dem60_69   = c("NA","24","11","NA","11","3","NA","20","1","99","22","10","NA","4","NA"),

  total70_74 = c("NA","117","NA","NA","238","163","NA","556","648","NA","311","470","NA","NA","NA"),
  dem70_74   = c("3.9 (2.7–5.5)","18","NA","NA","15","8","NA","11","8","NA","22","15","NA","NA","NA"),

  total75_79 = c("NA","101","NA","NA","154","158","NA","346","348","NA","252","266","NA","NA","NA"),
  dem75_79   = c("8.7 (6.5–11.4)","27","NA","NA","15","14","NA","27","9","NA","24","21","NA","NA","NA"),

  total70_79 = c("NA","218","94","NA","392","321","NA","902","996","148","563","736","NA","25","NA"),
  dem70_79   = c("NA","45","10","NA","30","22","NA","38","17","78","46","36","NA","6","NA"),

  total80_84 = c("NA","63","NA","NA","84","154","NA","169","224","73","128","198","NA","NA","NA"),
  dem80_84   = c("12.6 (9.2–17.0)","19","NA","NA","13","27","NA","23","7","39","17","30","NA","NA","NA"),

  total85_89 = c("NA","45","NA","NA","27","79","NA","112","166","NA","72","108","NA","NA","NA"),
  dem85_89   = c("21.3 (16.5–27.1)","22","NA","NA","5","32","NA","24","20","NA","11","42","NA","NA","NA"),

  total80_plus = c("NA","108","41","NA","102","233","NA","281","390","88","236","306","NA","NA","NA"),
  dem80_plus   = c("NA","41","17","NA","26","59","NA","47","27","52","39","72","NA","NA","NA"),

  total90_plus = c("NA","NA","NA","NA","17","NA","NA","NA","NA","15","26","NA","NA","NA","NA"),
  dem90_plus   = c("42.8 (32.8–53.4)","NA","NA","NA","8","NA","NA","NA","NA","13","11","NA","NA","NA","NA"),

  totalF = c("3144","397","202","75","726","484","115","1255","1291","259","1011","979","806","82","NA"),
  demF   = c("249","62","27","18","43","90","46","68","29","130","74","83","80","9","NA"),

  totalM = c("2105","233","97","21","419","199","48","817","931","207","466","677","587","59","NA"),
  demM   = c("115","48","11","3","25","25","14","37","16","101","33","35","33","9","NA"),

  totalClasseA = c("NA","26","NA","NA","117","NA","NA","NA","NA","NA","70","NA","NA","NA","NA"),
  demClasseA   = c("NA","3","NA","NA","8","NA","NA","NA","NA","NA","2","NA","NA","NA","NA"),
  totalClasseB = c("NA","184","NA","NA","364","121","NA","NA","NA","NA","465","NA","NA","NA","NA"),
  demClasseB   = c("NA","17","NA","NA","14","13","NA","NA","NA","NA","24","NA","NA","NA","NA"),
  totalClasseC = c("NA","315","NA","NA","317","202","NA","NA","NA","NA","492","NA","NA","NA","NA"),
  demClasseC   = c("NA","66","NA","NA","19","30","NA","NA","NA","NA","25","NA","NA","NA","NA"),
  totalClasseD = c("NA","104","NA","NA","186","225","NA","NA","NA","NA","298","NA","NA","NA","NA"),
  demClasseD   = c("NA","23","NA","NA","14","36","NA","NA","NA","NA","32","NA","NA","NA","NA"),
  totalClasseE = c("NA","1","NA","NA","61","101","NA","NA","NA","NA","133","NA","NA","NA","NA"),
  demClasseE   = c("NA","1","NA","NA","9","28","NA","NA","NA","NA","21","NA","NA","NA","NA"),

  totallessprimary = c("1229","NA","51","NA","NA","108","NA","NA","NA","NA","NA","NA","NA","NA","NA"),
  demlessprimary   = c("70","NA","16","NA","NA","34","NA","NA","NA","NA","NA","NA","NA","NA","NA"),
  totalprimary     = c("1561","NA","68","NA","449","139","NA","NA","NA","NA","NA","590","NA","60","NA"),
  demprimary       = c("73","NA","5","NA","27","30","NA","NA","NA","NA","NA","26","NA","9","NA"),
  totalmiddle      = c("NA","83","NA","NA","100","213","NA","NA","NA","NA","184","356","NA","32","NA"),
  demmiddle        = c("NA","8","NA","NA","5","20","NA","NA","NA","NA","5","18","NA","2","NA"),
  totalhigh        = c("592","51","132","NA","197","195","NA","NA","NA","NA","145","NA","NA","26","NA"),
  demhigh          = c("11","2","12","NA","7","15","NA","NA","NA","NA","6","NA","NA","1","NA"),
  totalcollege     = c("287","61","48","NA","203","NA","NA","NA","NA","NA","291","NA","NA","2","NA"),
  demcollege       = c("5","4","5","NA","6","NA","NA","NA","NA","NA","9","NA","NA","0","NA"),
  totaliletrado    = c("1132","87","56","NA","105","25","NA","NA","NA","348","230","567","NA","21","NA"),
  demilitrado      = c("190","20","15","NA","19","13","NA","NA","NA","177","43","69","NA","6","NA"),

  stringsAsFactors = FALSE,
  check.names = FALSE
)

# COVID-19 pandemic period: pre-pandemic = data year < 2020
df$covid_period <- ifelse(df$Year < 2020, "pre_covid", "post_covid")

# =============================================================================
# 2. Harmonization helper functions
# =============================================================================

clean_num <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub(",", "", x)
  x <- gsub("\\.(?=\\d{3}(\\D|$))", "", x, perl = TRUE)
  x[x %in% c("", "NA")] <- NA
  suppressWarnings(as.numeric(x))
}

is_prev_ci <- function(x) {
  x <- trimws(as.character(x))
  !is.na(x) & x != "" & str_detect(x, "\\(")
}

is_count <- function(x) {
  x <- trimws(as.character(x))
  !is.na(x) & x != "" & !str_detect(x, "\\(") & x != "NA"
}

extract_prev_ci <- function(x) {
  x <- trimws(as.character(x))
  if (is.na(x) || x == "" || !str_detect(x, "\\(")) return(c(NA_real_, NA_real_, NA_real_))
  x <- str_replace_all(x, "\u2013", "-")
  nums <- str_extract_all(x, "\\d+\\.?\\d*")[[1]]
  if (length(nums) < 3) return(c(NA_real_, NA_real_, NA_real_))
  c(as.numeric(nums[1]), as.numeric(nums[2]), as.numeric(nums[3]))
}

calc_prev_from_counts <- function(x, n) {
  if (is.na(x) || is.na(n) || n <= 0 || x < 0 || x > n) return(c(NA_real_, NA_real_, NA_real_))
  bt <- prop.test(x = x, n = n, correct = FALSE)
  c(x / n, bt$conf.int[1], bt$conf.int[2])
}

logit <- function(p) log(p / (1 - p))

# =============================================================================
# 3. Analysis map — outcome × domain × subgroup combinations
# =============================================================================

analysis_map <- tribble(
  ~outcome,    ~domain,         ~subgroup,       ~total_col,          ~event_col,
  "dementia",  "overall",       "overall",       "N_total",           "N_dementia",
  "CIND",      "overall",       "overall",       "N_total",           "N_CIND",
  "dementia",  "race",          "White",         "White",             "dem_white",
  "dementia",  "race",          "Brown",         "Brown",             "dem_Brown",
  "dementia",  "race",          "Black",         "Black",             "dem_Black",
  "dementia",  "race",          "Asian",         "Asian",             "dem_Asian",
  "dementia",  "age",           "60-64",         "total60_64",        "dem60_64",
  "dementia",  "age",           "65-69",         "total65_69",        "dem65_69",
  "dementia",  "age",           "60-69",         "total60_69",        "dem60_69",
  "dementia",  "age",           "70-74",         "total70_74",        "dem70_74",
  "dementia",  "age",           "75-79",         "total75_79",        "dem75_79",
  "dementia",  "age",           "70-79",         "total70_79",        "dem70_79",
  "dementia",  "age",           "80-84",         "total80_84",        "dem80_84",
  "dementia",  "age",           "85-89",         "total85_89",        "dem85_89",
  "dementia",  "age",           "80+",           "total80_plus",      "dem80_plus",
  "dementia",  "age",           "90+",           "total90_plus",      "dem90_plus",
  "dementia",  "sex",           "Female",        "totalF",            "demF",
  "dementia",  "sex",           "Male",          "totalM",            "demM",
  "dementia",  "social_class",  "ClasseA",       "totalClasseA",      "demClasseA",
  "dementia",  "social_class",  "ClasseB",       "totalClasseB",      "demClasseB",
  "dementia",  "social_class",  "ClasseC",       "totalClasseC",      "demClasseC",
  "dementia",  "social_class",  "ClasseD",       "totalClasseD",      "demClasseD",
  "dementia",  "social_class",  "ClasseE",       "totalClasseE",      "demClasseE",
  "dementia",  "education",     "LessPrimary",   "totallessprimary",  "demlessprimary",
  "dementia",  "education",     "Primary",       "totalprimary",      "demprimary",
  "dementia",  "education",     "Middle",        "totalmiddle",       "demmiddle",
  "dementia",  "education",     "High",          "totalhigh",         "demhigh",
  "dementia",  "education",     "College",       "totalcollege",      "demcollege",
  "dementia",  "education",     "Illiterate",    "totaliletrado",     "demilitrado"
)

analysis_map$filter_col <- NA_character_
analysis_map$filter_val <- NA_character_

analysis_map <- bind_rows(
  analysis_map,
  tribble(
    ~outcome,    ~domain,         ~subgroup,     ~total_col,  ~event_col,    ~filter_col,    ~filter_val,
    "dementia",  "covid_period",  "pre_covid",   "N_total",   "N_dementia",  "covid_period", "pre_covid",
    "dementia",  "covid_period",  "post_covid",  "N_total",   "N_dementia",  "covid_period", "post_covid"
  )
)

# =============================================================================
# 4. Harmonize all prevalence data to a unified long-format data frame
# =============================================================================

harmonized_prev <- map_dfr(seq_len(nrow(analysis_map)), function(i) {
  outcome_i    <- analysis_map$outcome[i]
  domain_i     <- analysis_map$domain[i]
  subgroup_i   <- analysis_map$subgroup[i]
  total_i      <- analysis_map$total_col[i]
  event_i      <- analysis_map$event_col[i]
  filter_col_i <- analysis_map$filter_col[i]
  filter_val_i <- analysis_map$filter_val[i]

  df_i <- if (!is.na(filter_col_i)) {
    df %>% filter(.data[[filter_col_i]] == filter_val_i)
  } else {
    df
  }

  tmp <- df_i %>%
    transmute(
      Study, Year, Design, Region, Coleta, Zona, DiagnosticCriteria,
      outcome   = outcome_i,
      domain    = domain_i,
      subgroup  = subgroup_i,
      total_raw = .data[[total_i]],
      event_raw = .data[[event_i]]
    ) %>%
    filter(!(is.na(total_raw) & is.na(event_raw)))

  from_counts <- tmp %>%
    filter(is_count(total_raw), is_count(event_raw)) %>%
    transmute(
      Study, Year, Design, Region, Coleta, Zona, DiagnosticCriteria,
      outcome, domain, subgroup, source = "count",
      total = clean_num(total_raw),
      cases = clean_num(event_raw)
    ) %>%
    filter(!is.na(total), !is.na(cases), total > 0, cases >= 0, cases <= total) %>%
    rowwise() %>%
    mutate(vals = list(calc_prev_from_counts(cases, total))) %>%
    mutate(prev = vals[[1]][1], ci.lb = vals[[1]][2], ci.ub = vals[[1]][3]) %>%
    ungroup() %>%
    select(-vals)

  from_prev <- tmp %>%
    filter(is_prev_ci(event_raw)) %>%
    mutate(vals = map(event_raw, extract_prev_ci)) %>%
    mutate(
      total  = clean_num(total_raw),
      cases  = NA_real_,
      prev   = map_dbl(vals, 1) / 100,
      ci.lb  = map_dbl(vals, 2) / 100,
      ci.ub  = map_dbl(vals, 3) / 100,
      source = "reported_prev"
    ) %>%
    select(Study, Year, Design, Region, Coleta, Zona, DiagnosticCriteria,
           outcome, domain, subgroup, source, total, cases, prev, ci.lb, ci.ub)

  bind_rows(from_counts, from_prev)
})

# Compute logit-scale effect sizes (yi) and variances (vi) for meta-analysis
prepare_meta_data <- function(dat) {
  dat_count <- dat %>%
    filter(!is.na(cases), !is.na(total), total > 0, cases >= 0, cases <= total) %>%
    mutate(cases = as.numeric(cases), total = as.numeric(total))

  dat_count_es <- if (nrow(dat_count) > 0) {
    escalc(measure = "PLO", xi = cases, ni = total, data = dat_count) %>%
      as.data.frame() %>%
      select(Study, Year, Design, Region, Coleta, Zona, DiagnosticCriteria,
             outcome, domain, subgroup, source, total, cases, prev, ci.lb, ci.ub, yi, vi)
  } else {
    data.frame()
  }

  dat_prev <- dat %>%
    filter(
      is.na(cases),
      !is.na(prev), !is.na(ci.lb), !is.na(ci.ub),
      prev > 0, prev < 1, ci.lb > 0, ci.ub < 1,
      ci.lb < prev, prev < ci.ub
    ) %>%
    mutate(
      yi = logit(prev),
      vi = ((logit(ci.ub) - logit(ci.lb)) / (2 * 1.96))^2
    ) %>%
    filter(is.finite(yi), is.finite(vi), !is.na(vi), vi > 0) %>%
    select(Study, Year, Design, Region, Coleta, Zona, DiagnosticCriteria,
           outcome, domain, subgroup, source, total, cases, prev, ci.lb, ci.ub, yi, vi)

  bind_rows(dat_count_es, dat_prev)
}

meta_ready <- prepare_meta_data(harmonized_prev)

# Convenience subset: overall all-cause dementia (used by multiple scripts)
dementia_overall_meta <- meta_ready %>%
  filter(outcome == "dementia", domain == "overall", subgroup == "overall")

# =============================================================================
# 5. Editorial label maps (Lancet / JAMA style)
# =============================================================================

lancet_outcome <- c(
  "dementia" = "All-cause dementia",
  "CIND"     = "Cognitive impairment, no dementia (CIND)"
)

lancet_domain <- c(
  "overall"      = "Overall",
  "age"          = "Age group (years)",
  "sex"          = "Sex",
  "race"         = "Race/skin colour (self-reported, IBGE)",
  "social_class" = "Socioeconomic class (ABEP criteria)",
  "education"    = "Educational attainment",
  "covid_period" = "COVID-19 pandemic period"
)

lancet_subgroup <- c(
  "overall"     = "Overall",
  "60-64"       = "60\u201364 years",
  "65-69"       = "65\u201369 years",
  "60-69"       = "60\u201369 years",
  "70-74"       = "70\u201374 years",
  "75-79"       = "75\u201379 years",
  "70-79"       = "70\u201379 years",
  "80-84"       = "80\u201384 years",
  "85-89"       = "85\u201389 years",
  "80+"         = "\u226580 years",
  "90+"         = "\u226590 years",
  "Female"      = "Women",
  "Male"        = "Men",
  "White"       = "White (branca)",
  "Brown"       = "Brown (parda)",
  "Black"       = "Black (preta)",
  "Asian"       = "Asian (amarela)",
  "ClasseA"     = "Economic class A (highest)",
  "ClasseB"     = "Economic class B",
  "ClasseC"     = "Economic class C",
  "ClasseD"     = "Economic class D",
  "ClasseE"     = "Economic class E (lowest)",
  "LessPrimary" = "Less than primary (ISCED 0\u20131)",
  "Primary"     = "Primary education (ISCED 1)",
  "Middle"      = "Lower secondary education (ISCED 2)",
  "High"        = "Upper secondary education (ISCED 3\u20134)",
  "College"     = "Tertiary education (ISCED 5\u20138)",
  "Illiterate"  = "No formal education / illiterate",
  "pre_covid"   = "Pre-pandemic (\u22642019)",
  "post_covid"  = "Pandemic / post-pandemic (\u22652020)",
  # Grouped analysis labels
  "Sudeste"             = "South-east Brazil",
  "Nao_Sudeste"         = "Other regions",
  "Urbano"              = "Urban",
  "Rural_Misto"         = "Rural / mixed (urban\u2013rural)",
  "DSM_IV"              = "DSM-IV",
  "Outros"              = "Other criteria (NIA-AA, CAMDEX, or not reported)",
  "Baixa_Escolaridade"  = "Low educational attainment (illiterate, less than primary, or primary)",
  "Maior_Escolaridade"  = "Higher educational attainment (lower secondary and above)"
)

lancet_moderator <- c(
  "Design"             = "Study design",
  "Region"             = "Geographic region",
  "Coleta"             = "Data collection setting",
  "Zona"               = "Area type (urban/rural)",
  "DiagnosticCriteria" = "Diagnostic criteria"
)

lancet_grouping <- c(
  "Sudeste_vs_Nao_Sudeste"      = "South-east vs other regions",
  "Urbano_vs_Rural_Misto"       = "Urban vs rural/mixed setting",
  "DSM_IV_vs_Outros"            = "DSM-IV vs other diagnostic criteria",
  "Baixa_vs_Maior_Escolaridade" = "Low vs higher educational attainment"
)

lbl <- function(x, map) ifelse(x %in% names(map), unname(map[x]), x)

# =============================================================================
# 6. Shared analysis helper functions
# =============================================================================

# Fit random-effects model (REML) for a given data subset
run_meta <- function(dat) {
  if (nrow(dat) == 0) return(NULL)
  tryCatch(rma.uni(yi = yi, vi = vi, data = dat, method = "REML"), error = function(e) NULL)
}

# Extract tidy summary from a fitted rma model
extract_meta_summary <- function(model, label) {
  if (is.null(model)) return(NULL)
  pr <- predict(model, transf = transf.ilogit)
  tibble(
    subgroup    = label,
    k           = model$k,
    pooled_prev = pr$pred,
    ci.lb       = pr$ci.lb,
    ci.ub       = pr$ci.ub,
    pi.lb       = pr$pi.lb,
    pi.ub       = pr$pi.ub,
    tau2        = model$tau2,
    I2          = model$I2,
    H2          = model$H2,
    QE          = model$QE,
    QEp         = model$QEp
  )
}

# Subgroup moderator test (between-group heterogeneity, Q_M statistic)
run_grouped_subgroup <- function(data, varname, labels = NULL) {
  dat <- data %>% filter(!is.na(.data[[varname]]), .data[[varname]] != "")
  dat[[varname]] <- if (!is.null(labels)) factor(dat[[varname]], levels = labels) else factor(dat[[varname]])
  fits <- split(dat, dat[[varname]])
  bind_rows(lapply(names(fits), function(nm) extract_meta_summary(fits[[nm]] |> run_meta(), nm)))
}

run_grouped_moderator_test <- function(data, varname) {
  dat <- data %>% filter(!is.na(.data[[varname]]), .data[[varname]] != "")
  if (nrow(dat) < 2) return(NULL)
  dat[[varname]] <- factor(dat[[varname]])
  fit <- tryCatch(
    rma.uni(yi = yi, vi = vi, mods = as.formula(paste0("~ `", varname, "` - 1")),
            data = dat, method = "REML"),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)
  tibble(moderator = varname, QM = fit$QM, QM_df = fit$m, QM_pval = fit$QMp)
}

cat("Script 01 complete: meta_ready (", nrow(meta_ready), "rows,",
    n_distinct(meta_ready$Study), "studies) created.\n")
