# Produces: output/tables/primary_age_band_contributions.csv (Supplementary Table S6b)
# Number of studies, participants and observed cases contributing to each age band of the primary age model
# (five populations without overlapping samples).
suppressPackageStartupMessages(library(dplyr))
source(file.path("code", "helpers", "setup.R"))
bands <- c("60-64", "65-69", "70-74", "75-79", "80-84", "85+")
agec <- read.csv(file.path(DERIVED, "age_strata_corrected_with_epifloripa.csv"))
prim <- subset(agec, study %in% PRIMARY_POPULATIONS & !is.na(n))
tab <- prim %>%
  mutate(band = factor(band, levels = bands)) %>%
  group_by(band) %>%
  summarise(k_studies = n_distinct(study), studies = paste(sort(unique(study)), collapse = ";"),
            participants = sum(n), cases_observed = sum(cases, na.rm = TRUE), .groups = "drop")
stopifnot(nrow(tab) == length(bands), all(tab$k_studies >= 1))
print(as.data.frame(tab))
write.csv(tab, file.path(OUT_TAB, "primary_age_band_contributions.csv"), row.names = FALSE)
