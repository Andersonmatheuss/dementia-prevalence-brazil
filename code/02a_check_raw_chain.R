# 02 — Confere a cadeia dos dados brutos: PrevDem.xlsx (planilha original, T7) -> dados embutidos no
# script 01 do GitHub (mar/2026) -> entradas de out/2026 (meta_ready.csv). Compara célula a célula
# as contagens que entram nas análises. Só leitura; grava o relatório de divergências.
suppressPackageStartupMessages({ library(readxl); library(dplyr); library(tidyr) })
source(file.path("code", "helpers", "setup.R"))
num <- function(x) suppressWarnings(as.numeric(gsub(",", "", trimws(as.character(x)))))

x <- suppressMessages(read_excel(file.path(RAW, "Prevalence_Dementia_23_03", "PrevDem.xlsx"), sheet = "Sheet2"))
names(x)[1] <- "Study"
# script 01 do GitHub: executa só a construção do data.frame df (até a linha 'stringsAsFactors')
src <- readLines(file.path(GH_DIR, "01_data_harmonization.R"))
a <- grep("^study_names <- c\\(", src); b <- grep("^\\s*check.names = FALSE", src)[1] + 1
env <- new.env(); eval(parse(text = src[a:b]), envir = env); gh <- env$df
cols_x  <- c("N total", "N dementia", "total60-64", "dem60-64", "total65-69", "dem65-69", "total70-74", "dem70-74", "total75-79", "dem75-79",
             "total80-84", "dem80-84", "total85-89", "dem85-89", "total80+", "dem80+", "total90+", "dem90+", "totalF", "demF")
cols_gh <- c("N_total", "N_dementia", "total60_64", "dem60_64", "total65_69", "dem65_69", "total70_74", "dem70_74", "total75_79", "dem75_79",
             "total80_84", "dem80_84", "total85_89", "dem85_89", "total80_plus", "dem80_plus", "total90_plus", "dem90_plus", "totalF", "demF")
stopifnot(nrow(x) == nrow(gh))
cmp <- bind_rows(lapply(seq_along(cols_x), function(j) data.frame(row = seq_len(nrow(gh)), study_gh = gh$Study, study_xlsx = x$Study,
  field = cols_gh[j], xlsx = as.character(x[[cols_x[j]]]), github = as.character(gh[[cols_gh[j]]]))))
cmp$equal <- mapply(function(a, b) { na <- num(a); nb <- num(b); (is.na(na) && is.na(nb) && (toupper(trimws(a)) == toupper(trimws(b)) | (is.na(a) | a %in% c("NA","")) & (is.na(b) | b %in% c("NA","")))) || isTRUE(all.equal(na, nb)) }, cmp$xlsx, cmp$github)
cat("células comparadas:", nrow(cmp), "| divergentes:", sum(!cmp$equal), "\n")
print(subset(cmp, !equal), row.names = FALSE)
# GitHub (mar) x entradas de out/2026: total e casos globais por estudo
oct <- read.csv(file.path(OCT_OUT, "study_level_data", "meta_ready.csv")) |> filter(outcome == "dementia", domain == "overall")
key <- function(s) sub("_\\d{4}$", "", s)
m <- merge(data.frame(k = key(gh$Study), gh_total = num(gh$N_total), gh_cases = num(gh$N_dementia)),
           data.frame(k = key(oct$Study), oct_total = oct$total, oct_cases = oct$cases), all = TRUE)
print(m); write.csv(subset(cmp, !equal), file.path(OUT_TAB, "raw_chain_discrepancies.csv"), row.names = FALSE)
