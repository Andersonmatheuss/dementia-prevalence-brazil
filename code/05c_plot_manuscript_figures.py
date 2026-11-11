#!/usr/bin/env python3
# Produces: output/figures/figure2_age.{pdf,png}, figure3_subgroups.{pdf,png}, figure4_regions.{pdf,png} (Figures 2-4)
"""Figuras 2 a 4 do manuscrito, no formato da versão de abril, com o conjunto principal corrigido.

Figura 2: prevalência por faixa etária (estimativas dos autores e contagens observadas) e carga por faixa em 2022.
Figura 3: prevalência por escolaridade, sexo e classe socioeconômica (conjunto principal verificado, contagens observadas).
Figura 4: prevalência e carga esperadas por macrorregião, 2022 a 2050 (projeções IBGE, revisão 2024).

Lê só CSVs de output/tables/ e a contagem do Censo 2022 por faixa etária (data/raw, SIDRA 9514).
Uso (o code/run_all.sh chama): ~/.venvs/plots/bin/python -I code/05c_plot_manuscript_figures.py
"""
import csv
import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAB = ROOT / "output" / "tables"
FIG = ROOT / "output" / "figures"
FIG.mkdir(parents=True, exist_ok=True)
CENSUS = ROOT / "data" / "raw" / "Dementia_prevalence" / "ibge_sidra_inputs" / "ibge_census_brazil_2022_agegroups.csv"

spec = importlib.util.spec_from_file_location("figure_style", Path(__file__).resolve().parent / "helpers" / "figure_style.py")
fs = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fs)

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

matplotlib.rcParams["pdf.fonttype"] = 42
BANDS = ["60-64", "65-69", "70-74", "75-79", "80-84", "85+"]
LABEL = {"60-64": "60–64", "65-69": "65–69", "70-74": "70–74", "75-79": "75–79", "80-84": "80–84", "85+": "≥85"}
REG_PT = {"Norte": "North", "Nordeste": "Northeast", "Sudeste": "Southeast", "Sul": "South", "Centro-Oeste": "Central-West"}


def rows(name):
    return list(csv.DictReader(open(TAB / name, encoding="utf-8")))


def save(fig, stem):
    for ext in ("pdf", "png"):
        fig.savefig(FIG / f"{stem}.{ext}", dpi=600, metadata={"CreationDate": None} if ext == "pdf" else None)
    plt.close(fig)


def census_bands():
    n = {}
    for r in csv.DictReader(open(CENSUS, encoding="utf-8")):
        if r["sex"] == "Total":
            b = "85+" if r["age_group"] in ("85-89", "90+") else r["age_group"]
            n[b] = n.get(b, 0) + int(r["population"])
    return n


def figure2(base):
    st = {(r["version"], r["band"]): r for r in rows("primary_age_strata.csv")}
    n = census_bands()
    fig, (a, b) = plt.subplots(1, 2, figsize=(180 * fs.MM, 75 * fs.MM), constrained_layout=True)
    x = np.arange(len(BANDS))
    for off, (v, lab, col) in zip((-0.12, 0.12), [("author_best", "Authors' population estimates", fs.SERIES[0]),
                                                   ("observed", "Observed counts", fs.SERIES[1])]):
        p = np.array([float(st[(v, k)]["prev"]) for k in BANDS]) * 100
        lo = np.array([float(st[(v, k)]["lb"]) for k in BANDS]) * 100
        hi = np.array([float(st[(v, k)]["ub"]) for k in BANDS]) * 100
        a.errorbar(x + off, p, yerr=[p - lo, hi - p], fmt="o", ms=3.5, color=col, capsize=0, lw=1, label=lab)
    a.set_xticks(x, [LABEL[k] for k in BANDS]); a.set_xlabel("Age group (years)"); a.set_ylabel("Prevalence (%)")
    a.legend(frameon=False, loc="upper left"); fs.grid_y_only(a); fs.hide_spines(a); fs.panel_label(a, "A", base)
    p = np.array([float(st[("author_best", k)]["prev"]) * n[k] for k in BANDS]) / 1e3
    lo = np.array([float(st[("author_best", k)]["lb"]) * n[k] for k in BANDS]) / 1e3
    hi = np.array([float(st[("author_best", k)]["ub"]) * n[k] for k in BANDS]) / 1e3
    b.bar(x, p, color=fs.SERIES[0], width=0.65)
    b.errorbar(x, p, yerr=[p - lo, hi - p], fmt="none", ecolor=fs.INK, lw=0.8)
    b.set_xticks(x, [LABEL[k] for k in BANDS]); b.set_xlabel("Age group (years)")
    b.set_ylabel("People with dementia, 2022 (thousands)")
    fs.grid_y_only(b); fs.hide_spines(b); fs.panel_label(b, "B", base)
    save(fig, "figure2_age")


def figure3(base):
    sg = [r for r in rows("corrected_subgroups.csv") if r["set"] == "verified (primary, no overlap)"]
    tests = {r["domain"]: r for r in rows("corrected_subgroup_tests.csv") if r["set"] == "verified (primary, no overlap)"}
    panels = [("education", [("none", "No formal education"), ("any", "Any formal education")], "Educational attainment"),
              ("sex", [("female", "Women"), ("male", "Men")], "Sex"),
              ("class", [(c, f"Class {c}") for c in "ABCDE"], "Socioeconomic class")]
    fig, axes = plt.subplots(1, 3, figsize=(180 * fs.MM, 70 * fs.MM), constrained_layout=True)
    for ax, (dom, cats, title), letter in zip(axes, panels, "ABC"):
        d = {r["subgroup"]: r for r in sg if (r["domain"] == dom or (dom == "sex" and r["domain"] in ("", "sex")))}
        p = np.array([float(d[c]["prev"]) for c, _ in cats]) * 100
        lo = np.array([float(d[c]["lb"]) for c, _ in cats]) * 100
        hi = np.array([float(d[c]["ub"]) for c, _ in cats]) * 100
        k = d[cats[0][0]]["k"]
        x = np.arange(len(cats))
        ax.errorbar(x, p, yerr=[p - lo, hi - p], fmt="o", ms=4, color=fs.SERIES[0], lw=1, capsize=0)
        ax.set_xticks(x, [l for _, l in cats], rotation=0 if len(cats) < 4 else 0)
        t = tests[dom]
        pv = float(t["p"]); ptxt = "P < .001" if pv < 0.001 else f"P = {pv:.2f}".replace("0.", ".")
        ax.set_title(f"{title} (k = {k}; {ptxt})", fontsize=base, loc="left")
        ax.set_ylabel("Prevalence (%)" if letter == "A" else "")
        fs.grid_y_only(ax); fs.hide_spines(ax); fs.panel_label(ax, letter, base)
    save(fig, "figure3_subgroups")


def figure4(base):
    r = [x for x in rows("author_best_projections_regions.csv") if x["deff"] == "1"]
    d = {(REG_PT[x["local"]], int(x["year"])): x for x in r}
    regs = fs.REGION_ORDER
    fig, axes = plt.subplots(2, 2, figsize=(180 * fs.MM, 130 * fs.MM), constrained_layout=True)
    for ax, yr, letter in [(axes[0, 0], 2022, "A"), (axes[0, 1], 2050, "B")]:
        p = np.array([float(d[(g, yr)]["cases_median"]) / float(d[(g, yr)]["pop60"]) for g in regs]) * 100
        lo = np.array([float(d[(g, yr)]["cases_lb"]) / float(d[(g, yr)]["pop60"]) for g in regs]) * 100
        hi = np.array([float(d[(g, yr)]["cases_ub"]) / float(d[(g, yr)]["pop60"]) for g in regs]) * 100
        x = np.arange(len(regs))
        ax.bar(x, p, color=[fs.REGION[g] for g in regs], width=0.65)
        ax.errorbar(x, p, yerr=[p - lo, hi - p], fmt="none", ecolor=fs.INK, lw=0.8)
        ax.set_xticks(x, regs); ax.set_ylabel(f"Expected prevalence, {yr} (%)")
        fs.grid_y_only(ax); fs.hide_spines(ax); fs.panel_label(ax, letter, base)
    ax = axes[1, 0]
    yrs = [2022, 2030, 2050]
    w = 0.26
    x = np.arange(len(regs))
    for i, yr in enumerate(yrs):
        v = np.array([float(d[(g, yr)]["cases_median"]) for g in regs]) / 1e3
        ax.bar(x + (i - 1) * w, v, width=w, color=fs.SEQUENTIAL[2 + 2 * i], label=str(yr))
    ax.set_xticks(x, regs); ax.set_ylabel("People with dementia (thousands)")
    ax.legend(frameon=False, loc="upper left"); fs.grid_y_only(ax); fs.hide_spines(ax); fs.panel_label(ax, "C", base)
    ax = axes[1, 1]
    for i, yr in enumerate([2030, 2050]):
        v = np.array([(float(d[(g, yr)]["cases_median"]) / float(d[(g, 2022)]["cases_median"]) - 1) * 100 for g in regs])
        ax.bar(x + (i - 0.5) * 0.36, v, width=0.36, color=fs.SEQUENTIAL[4 + 2 * i], label=str(yr))
    ax.set_xticks(x, regs); ax.set_ylabel("Increase relative to 2022 (%)")
    ax.legend(frameon=False, loc="upper center", ncol=2); fs.grid_y_only(ax); fs.hide_spines(ax); fs.panel_label(ax, "D", base)
    save(fig, "figure4_regions")


if __name__ == "__main__":
    base = fs.apply_style(180)
    figure2(base)
    figure3(base)
    figure4(base)
