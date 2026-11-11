# Produces: output/figures/fig2_national.{pdf,png}, fig3_age.{pdf,png}, figS1_subgroups.{pdf,png}, figS2_regions.{pdf,png}, figS3_maps.{pdf,png} (Figures 2-3 and S1-S3)
# Figuras do manuscrito com tema único (painéis da versão de abril, dados do conjunto principal corrigido).
# Entradas: output/tables/*.csv, output/figures/national_estimates_data.csv (code/05a),
# output/intermediate/models/primary_m3_author_best.rds, Censo 2022 por faixa (data/raw/Dementia_prevalence/ibge_sidra_inputs)
# e estratos verificados (data/raw/extraction/fulltext_strata.csv). O mapa vem do geobr (read_region, 2020), baixado na execução.

suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(dplyr); library(tidyr); library(sf); library(geobr); library(brms)
})
set.seed(20261108)
args <- commandArgs(FALSE)
root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", args, value = TRUE))), ".."))
TAB <- file.path(root, "output", "tables"); FIG <- file.path(root, "output", "figures")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

# ---- tema ------------------------------------------------------------------------------------
INK <- "#1F2933"; SOFT <- "#5B6770"; RULE <- "#C9D1D9"; GRIDC <- "#ECEFF2"
PRIMARY <- "#1F5A8C"; OBS <- "#C0703C"; BAYES <- "#3E8E7E"; GREY <- "#8A949E"
REG_COL <- c(North = "#2F6F8F", Northeast = "#C0703C", `Central-West` = "#7A9E4D", Southeast = "#6C5B9E", South = "#B04A5A")
REG_ORDER <- names(REG_COL)
REG_PT <- c(Norte = "North", Nordeste = "Northeast", `Centro-Oeste` = "Central-West", Sudeste = "Southeast", Sul = "South")
FONT <- "Helvetica"

th <- theme_minimal(base_size = 10, base_family = FONT) +
  theme(text = element_text(colour = INK),
        plot.title = element_text(face = "bold", size = 11, hjust = 0, margin = margin(b = 2)),
        plot.subtitle = element_text(size = 9, colour = SOFT, hjust = 0, margin = margin(b = 8)),
        plot.title.position = "plot", plot.caption.position = "plot",
        plot.caption = element_text(size = 7.5, colour = SOFT, hjust = 0, margin = margin(t = 8)),
        axis.title = element_text(size = 9, colour = SOFT), axis.text = element_text(size = 8.5, colour = INK),
        panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(colour = GRIDC, linewidth = 0.4),
        axis.line.x = element_line(colour = RULE, linewidth = 0.5), axis.ticks.x = element_line(colour = RULE, linewidth = 0.4),
        legend.position = "top", legend.justification = "left", legend.title = element_blank(),
        legend.text = element_text(size = 8.5, colour = INK), legend.key.height = unit(3, "mm"),
        strip.text = element_text(face = "bold", size = 9, hjust = 0, colour = INK),
        plot.margin = margin(8, 10, 8, 8))
pctl <- function(x) paste0(x, "%")
save2 <- function(p, stem, w, h) {
  ggsave(file.path(FIG, paste0(stem, ".pdf")), p, width = w, height = h, units = "mm", device = cairo_pdf, bg = "white")
  ggsave(file.path(FIG, paste0(stem, ".png")), p, width = w, height = h, units = "mm", dpi = 300, bg = "white")
}

# ---- dados -----------------------------------------------------------------------------------
bands <- c("60-64", "65-69", "70-74", "75-79", "80-84", "85+")
blab <- c("60–64", "65–69", "70–74", "75–79", "80–84", "≥85")
age <- read.csv(file.path(TAB, "primary_age_strata.csv")) |> mutate(band = factor(band, bands, blab))
cen <- read.csv(file.path(root, "data/raw/Dementia_prevalence/ibge_sidra_inputs/ibge_census_brazil_2022_agegroups.csv")) |>
  filter(sex == "Total") |> mutate(band = ifelse(age_group %in% c("85-89", "90+"), "85+", age_group)) |>
  group_by(band) |> summarise(pop = sum(population), .groups = "drop") |> mutate(band = factor(band, bands, blab))
m <- readRDS(file.path(root, "output/intermediate/models/primary_m3_author_best.rds"))
dr <- as_draws_df(m)
bcols <- paste0("b_band", c("60M64", "65M69", "70M74", "75M79", "80M84", "85P"))
bay <- data.frame(band = factor(blab, blab), t(sapply(bcols, function(cn) quantile(plogis(dr[[cn]]), c(0.5, 0.025, 0.975)))))
names(bay)[2:4] <- c("prev", "lb", "ub")
regall <- read.csv(file.path(TAB, "author_best_projections_regions.csv")) |> filter(deff == 1) |>
  mutate(region = factor(REG_PT[local], levels = REG_ORDER), prev = cases_median / pop60 * 100,
         lb = cases_lb / pop60 * 100, ub = cases_ub / pop60 * 100)
nat <- read.csv(file.path(TAB, "author_best_burden_projections.csv")) |> filter(deff == 1)
est <- read.csv(file.path(TAB, "primary_burden_estimands.csv")) |> filter(version == "author_best", estimand == "typical study")
sg <- read.csv(file.path(TAB, "corrected_subgroups.csv")) |> filter(set == "verified (primary, no overlap)")
tst <- read.csv(file.path(TAB, "corrected_subgroup_tests.csv")) |> filter(set == "verified (primary, no overlap)")
strata <- read.csv(file.path(root, "data/raw/extraction/fulltext_strata.csv"))
prim <- c("Herrera_2002", "Lopes_2012", "Scazufca_2008", "Cesar_2015", "EpiFloripa_2024")
nstr <- strata |> filter(study %in% prim) |> group_by(domain, category) |> summarise(n = sum(n), .groups = "drop")

# ---- Figura 2: estimativas nacionais ----------------------------------------------------------
ne <- read.csv(file.path(FIG, "national_estimates_data.csv"))
grp <- unique(ne$group)
rows <- do.call(rbind, lapply(grp, function(g) rbind(
  data.frame(group = g, label = g, kind = "head", estimate_millions = NA, lower_millions = NA, upper_millions = NA),
  ne[ne$group == g, ])))
rows$y <- -seq_len(nrow(rows))
rows$txt <- with(rows, ifelse(kind == "head", "",
  ifelse(is.na(estimate_millions), sprintf("%.2f–%.2f", lower_millions, upper_millions),
  ifelse(is.na(lower_millions), sprintf("%.2f", estimate_millions),
         sprintf("%.2f (%.2f–%.2f)", estimate_millions, lower_millions, upper_millions)))))
kcol <- c(published = GREY, census = INK, main = PRIMARY, scenario = OBS, range = OBS)
pts <- filter(rows, kind != "head")
band_obs <- ne$estimate_millions[ne$label == "Observed counts"]
band_pri <- ne$estimate_millions[ne$label == "Authors' population estimates (primary)"]
p2 <- ggplot(pts) +
  annotate("rect", xmin = band_obs, xmax = band_pri, ymin = -Inf, ymax = Inf, fill = "#F4EFE8") +
  geom_segment(data = filter(pts, kind != "range"), aes(x = lower_millions, xend = upper_millions, y = y, yend = y, colour = kind),
               linewidth = 0.8, na.rm = TRUE) +
  geom_segment(data = filter(pts, kind == "range"), aes(x = lower_millions, xend = upper_millions, y = y, yend = y),
               colour = OBS, alpha = 0.35, linewidth = 3) +
  geom_point(aes(estimate_millions, y, colour = kind, shape = kind), size = 2.3, fill = "white", stroke = 1, na.rm = TRUE) +
  scale_colour_manual(values = kcol, guide = "none") +
  scale_shape_manual(values = c(published = 21, census = 16, main = 16, scenario = 16, range = 16), guide = "none") +
  scale_x_continuous(limits = c(1, 5), breaks = seq(1, 5, 0.5), expand = c(0, 0),
                     sec.axis = dup_axis(name = NULL, breaks = NULL)) +
  scale_y_continuous(breaks = rows$y, labels = rows$label, limits = c(min(rows$y) - 0.6, -0.4), expand = c(0, 0),
                     sec.axis = dup_axis(breaks = rows$y, labels = rows$txt, name = NULL)) +
  labs(x = "People with dementia (millions)", y = NULL,
       title = "Estimates of the number of people with dementia in Brazil",
       subtitle = "Shaded band: difference between the observed-count and the screening-corrected primary estimate, 2022 Census") +
  th + theme(axis.text.y.left = element_text(hjust = 0, size = 8.5, colour = INK,
                                             face = ifelse(rows$kind == "head", "bold", "plain")),
             axis.text.y.right = element_text(hjust = 1, size = 8, colour = SOFT),
             panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(colour = GRIDC, linewidth = 0.4))
save2(p2, "fig2_national", 180, 120)

# ---- Figura 3: idade -------------------------------------------------------------------------
ab <- age |> filter(version == "author_best"); ob <- age |> filter(version == "observed")
dd <- bind_rows(mutate(ab, series = "Authors' population estimates (primary)"), mutate(ob, series = "Observed counts"),
                mutate(bay, series = "Bayesian model (primary)")) |>
  mutate(series = factor(series, c("Authors' population estimates (primary)", "Bayesian model (primary)", "Observed counts")),
         x = as.integer(band), prev = prev * 100, lb = lb * 100, ub = ub * 100)
scol <- c(`Authors' population estimates (primary)` = PRIMARY, `Bayesian model (primary)` = BAYES, `Observed counts` = OBS)
dodge <- c(`Authors' population estimates (primary)` = -0.12, `Bayesian model (primary)` = 0, `Observed counts` = 0.12)
dd$xd <- dd$x + dodge[as.character(dd$series)]
p3a <- ggplot(dd, aes(xd, prev, colour = series)) +
  geom_ribbon(data = filter(dd, series == "Authors' population estimates (primary)"), aes(x = x, ymin = lb, ymax = ub),
              fill = PRIMARY, alpha = 0.10, colour = NA) +
  geom_line(aes(group = series), linewidth = 0.7) +
  geom_linerange(aes(ymin = lb, ymax = ub), linewidth = 0.6) +
  geom_point(size = 2, shape = 21, fill = "white", stroke = 0.9) +
  scale_colour_manual(values = scol) + scale_x_continuous(breaks = 1:6, labels = blab) +
  scale_y_continuous(labels = pctl, limits = c(0, 60), breaks = seq(0, 60, 10), expand = c(0, 0)) +
  labs(title = "A  Prevalence by age", x = "Age group (years)", y = "Prevalence (%, 95% CI)") + th
cases <- ab |> left_join(cen, by = "band") |> mutate(c = prev * pop / 1e6, lo = lb * pop / 1e6, hi = ub * pop / 1e6,
                                                     share = c / sum(c) * 100)
p3b <- ggplot(cases, aes(band, c)) +
  geom_col(fill = PRIMARY, width = 0.62, alpha = 0.9) +
  geom_linerange(aes(ymin = lo, ymax = hi), colour = INK, linewidth = 0.5) +
  geom_text(aes(y = hi, label = sprintf("%.2fM\n%.0f%%", c, share)), vjust = -0.25, size = 2.7, colour = INK, family = FONT, lineheight = 0.9) +
  scale_y_continuous(limits = c(0, 1.45), breaks = seq(0, 1.4, 0.2), expand = c(0, 0)) +
  labs(title = "B  Expected cases by age group, 2022", x = "Age group (years)",
       y = "People with dementia (millions, 95% CI)") + th
save2(p3a + p3b + plot_layout(widths = c(1.15, 1)) + plot_annotation(theme = theme(legend.position = "top")), "fig3_age", 180, 100)

# ---- Figura S1: subgrupos (dot-and-interval) -------------------------------------------------
lab_sg <- c(none = "No formal education", any = "Any formal education", female = "Women", male = "Men",
            A = "Class A", B = "Class B", C = "Class C", D = "Class D", E = "Class E")
dom_lab <- c(education = "A  Educational attainment", sex = "B  Sex", class = "C  Socioeconomic class")
s1 <- sg |> mutate(dom = factor(dom_lab[domain], dom_lab), lab = lab_sg[subgroup],
                   n = nstr$n[match(paste(domain, subgroup), paste(nstr$domain, nstr$category))],
                   prev = prev * 100, lb = lb * 100, ub = ub * 100)
s1$lab <- factor(s1$lab, rev(lab_sg))
pv <- tst |> mutate(dom = factor(dom_lab[domain], dom_lab), ptxt = ifelse(p < 0.001, "P < 0.001", sprintf("P = %.2f", p)),
                    k = sapply(domain, function(d) max(s1$k[s1$domain == d])))
s1 <- left_join(s1, select(pv, dom, ptxt), by = "dom")
pS1 <- ggplot(s1, aes(prev, lab)) +
  geom_segment(aes(x = lb, xend = ub, yend = lab), colour = PRIMARY, linewidth = 0.9, alpha = 0.55) +
  geom_point(size = 2.6, colour = PRIMARY) +
  geom_text(aes(x = 122, label = sprintf("%.1f%% (%.1f–%.1f)", prev, lb, ub)), hjust = 1, size = 2.8, colour = SOFT, family = FONT) +
  facet_grid(dom ~ ., scales = "free_y", space = "free_y", switch = "y") +
  geom_text(data = pv, aes(x = 122, y = Inf, label = sprintf("k = %d · %s", k, ptxt)), inherit.aes = FALSE,
            hjust = 1, vjust = 1.4, size = 2.7, colour = INK, fontface = "italic", family = FONT) +
  scale_x_continuous(labels = pctl, limits = c(0, 122), breaks = seq(0, 90, 10), expand = c(0, 0)) +
  labs(x = "Pooled prevalence (%, 95% CI), observed counts", y = NULL,
       title = "Dementia prevalence by educational attainment, sex and socioeconomic class",
       subtitle = "Primary set; k, number of studies; P, mixed-effects omnibus test") +
  th + theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0, face = "bold", hjust = 0, vjust = 1),
             panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(colour = GRIDC, linewidth = 0.4),
             panel.spacing.y = unit(4, "mm"))
save2(pS1, "figS1_subgroups", 180, 110)

# ---- Figura S2: regiões (área, prevalência, aumento) -----------------------------------------
pA <- ggplot(regall, aes(year, cases_median / 1e6, fill = region)) +
  geom_area(colour = "white", linewidth = 0.3, alpha = 0.92) +
  scale_fill_manual(values = REG_COL) +
  scale_x_continuous(breaks = seq(2030, 2070, 10) |> c(2022) |> sort(), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 10, 2), expand = expansion(mult = c(0, 0.03))) +
  labs(title = "A  Expected people with dementia by macroregion, 2022–2070", x = NULL, y = "Millions") + th
d2 <- regall |> filter(year %in% c(2022, 2050))
pB <- ggplot(d2, aes(prev, region)) +
  geom_segment(data = pivot_wider(select(d2, region, year, prev), names_from = year, values_from = prev, names_prefix = "y"),
               aes(x = y2022, xend = y2050, y = region, yend = region), inherit.aes = FALSE, colour = RULE, linewidth = 1.6) +
  geom_point(aes(colour = factor(year)), size = 3) +
  geom_text(aes(label = sprintf("%.1f", prev), colour = factor(year)), vjust = -1.1, size = 2.7, family = FONT, show.legend = FALSE) +
  scale_colour_manual(values = c(`2022` = GREY, `2050` = PRIMARY)) +
  scale_y_discrete(limits = rev(REG_ORDER)) +
  scale_x_continuous(labels = pctl, limits = c(7, 12), breaks = 7:12) +
  labs(title = "B  Prevalence, 2022 and 2050", x = "Prevalence among people aged ≥60 (%)", y = NULL) +
  th + theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(colour = GRIDC, linewidth = 0.4))
d3 <- regall |> group_by(region) |> mutate(chg = (cases_median / cases_median[year == 2022] - 1) * 100) |>
  filter(year %in% c(2030, 2050)) |> ungroup()
pC <- ggplot(d3, aes(chg, factor(region, rev(REG_ORDER)), fill = factor(year))) +
  geom_col(position = position_dodge(width = 0.72), width = 0.66) +
  geom_text(aes(label = sprintf("+%.0f%%", chg)), position = position_dodge(width = 0.72), hjust = -0.15, size = 2.6,
            colour = INK, family = FONT) +
  scale_fill_manual(values = c(`2030` = "#9DB8D0", `2050` = PRIMARY), breaks = c("2030", "2050")) +
  scale_x_continuous(labels = pctl, limits = c(0, 300), breaks = seq(0, 250, 50), expand = c(0, 0)) +
  labs(title = "C  Increase in cases vs 2022", x = "Increase (%)", y = NULL) +
  th + theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(colour = GRIDC, linewidth = 0.4))
save2(pA / (pB | pC) + plot_layout(heights = c(1, 1)), "figS2_regions", 180, 170)

# ---- Figura S3: mapas coropléticos -----------------------------------------------------------
reg <- regall |> filter(year %in% c(2022, 2030, 2050))
base <- reg |> filter(year == 2022) |> select(region, base = cases_median)
reg <- reg |> left_join(base, by = "region") |> mutate(chg = (cases_median / base - 1) * 100)
br <- read_region(year = 2020, showProgress = FALSE) |> mutate(region = REG_PT[name_region]) |> st_transform(5880)
cent <- suppressWarnings(st_point_on_surface(br)); xy <- cbind(st_drop_geometry(cent["region"]), st_coordinates(cent))
nudge <- data.frame(region = REG_ORDER, dx = c(0, 80e3, -60e3, 380e3, -420e3), dy = c(80e3, 0, 40e3, -120e3, -60e3))
pal <- c("#EEF3F7", "#C9DCEA", "#93BCD8", "#5A97C1", "#2F6FA3", "#16466F")
fillsc <- scale_fill_gradientn(colours = pal, limits = c(7, 11.5), breaks = 7:11, labels = pctl, name = "Expected prevalence, ≥60 years",
                               guide = guide_colourbar(barwidth = unit(60, "mm"), barheight = unit(2.4, "mm"), title.position = "top"))
one <- function(yr, letter) {
  d <- filter(reg, year == yr) |> mutate(region = as.character(region))
  g <- left_join(br, d, by = "region")
  lab <- left_join(xy, d, by = "region") |> left_join(nudge, by = "region") |>
    mutate(txt = if (yr == 2022) sprintf("%s\n%.2fM · %.1f%%", region, cases_median / 1e6, prev)
                 else sprintf("%s\n%.2fM · %.1f%%\n+%.0f%%", region, cases_median / 1e6, prev, chg))
  n <- filter(nat, year == yr); ncas <- if (yr == 2022) est$cases_median else n$cases_median
  sub <- sprintf("Brazil: %.2f million · %.1f%%%s", ncas / 1e6, ncas / n$population * 100,
                 if (yr == 2022) "" else sprintf(" · +%.0f%% vs 2022", (n$cases_median / nat$cases_median[nat$year == 2022] - 1) * 100))
  ggplot(g) + geom_sf(aes(fill = prev), colour = "white", linewidth = 0.5) +
    geom_label(data = lab, aes(X + dx, Y + dy, label = txt), size = 2.3, family = FONT, colour = INK, lineheight = 0.95,
               fill = alpha("white", 0.85), linewidth = 0, label.padding = unit(0.12, "lines")) +
    fillsc + labs(title = sprintf("%s  %d", letter, yr), subtitle = sub, x = NULL, y = NULL) + coord_sf(datum = NA) +
    th + theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_blank(),
               panel.grid.major.y = element_blank(), legend.position = "bottom", legend.title = element_text(size = 8.5, colour = SOFT))
}
ord <- reg |> filter(year == 2050) |> arrange(cases_median) |> pull(region) |> as.character()
bars <- reg |> mutate(region = factor(as.character(region), levels = ord))
pBar <- ggplot(bars, aes(cases_median / 1e6, region, fill = prev, group = factor(year))) +
  geom_col(position = position_dodge(width = 0.8), width = 0.74, colour = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%d  %.2fM", year, cases_median / 1e6)), position = position_dodge(width = 0.8),
            hjust = -0.08, size = 2.3, colour = INK, family = FONT) +
  fillsc + guides(fill = "none") + scale_x_continuous(limits = c(0, 4.2), breaks = 0:3, expand = c(0, 0)) +
  labs(title = "D  Expected people with dementia by macroregion", subtitle = "Bars coloured by expected prevalence in each year",
       x = "Millions", y = NULL) +
  th + theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(colour = GRIDC, linewidth = 0.4))
maps <- (one(2022, "A") | one(2030, "B") | one(2050, "C"))
pS3 <- (maps / pBar) + plot_layout(heights = c(1.25, 1), guides = "collect") & theme(legend.position = "bottom")
pS3 <- pS3 + plot_annotation(caption = "Map labels: expected people with dementia (millions) · expected prevalence · change in cases since 2022. National age schedule of the primary analysis applied to the 2024 IBGE projections;\nmodel-based, not region-specific estimates.",
                             theme = theme(plot.caption = element_text(size = 7.5, colour = SOFT, hjust = 0, family = FONT)))
save2(pS3, "figS3_maps", 250, 200)
cat("ok\n")
