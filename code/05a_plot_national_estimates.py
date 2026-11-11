#!/usr/bin/env python3
"""Figura 2: estimativas nacionais de pessoas com demência no Brasil, lado a lado.

Lê só CSVs gerados pelos scripts R (output/tables/) e grava PDF + PNG na largura final
(180 mm), com o estilo canônico do grupo (figure_style.py do estudo dementia-death, lido no
lugar, sem cópia). Grava também o CSV com exatamente o que foi desenhado.

Entradas fora dos CSVs, publicadas sem taxa por idade, o que impede levá-las ao Censo 2022:
- os 2,2 milhões da classificação por aprendizado de máquina (Cleret de Langavant 2020,
  PMID 32885026, resumo conferido em 19/10/2026), para o ELSI 2016;
- o GBD 2019 para o Brasil, todas as idades: 1.849.981 (II 95% 1.630.321 a 2.097.631) em 2019
  (GBD 2019 Dementia Forecasting Collaborators, Lancet Public Health 2022;7:e105-25, Tabela,
  p. e113; PDF conferido em 20/10/2026).

Uso (o code/run_all.sh chama): ~/.venvs/plots/bin/python -I code/05a_plot_national_estimates.py
"""
import csv
import hashlib
import importlib.util
import os
from pathlib import Path

AUDIT = Path(__file__).resolve().parents[1]          # raiz do projeto
TAB = AUDIT / "output" / "tables"
FIG = AUDIT / "output" / "figures"
FIG.mkdir(parents=True, exist_ok=True)

STYLE = Path(__file__).resolve().parent / "helpers" / "figure_style.py"
spec = importlib.util.spec_from_file_location("figure_style", STYLE)
fs = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fs)

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

matplotlib.rcParams["pdf.fonttype"] = 42   # texto pesquisável e editável no PDF


def rows(name):
    return list(csv.DictReader(open(TAB / name, encoding="utf-8")))


def pick(table, **kw):
    out = [r for r in table if all(r[k] == v for k, v in kw.items())]
    assert len(out) == 1, (kw, len(out))
    return out[0]


M = 1e6
tri = {r["source"]: float(r["cases"]) for r in rows("triangulation_national_estimates.csv")}
est = rows("primary_burden_estimands.csv")
bay = rows("primary_bayesian.csv")
# análise de viés no conjunto principal sem sobreposição (a outra linha do CSV é a das oito populações)
qba = {r["scenario"]: r for r in rows("qba_two_phase_with_epifloripa.csv") if r["set"] == "primary"}
sens = rows("primary_burden_sensitivity_sets.csv")
pba = rows("pba_screening_correction.csv")   # script 22
elsi_grid = [float(r["cases_2022"]) for r in rows("qba_elsi_algorithm_misclassification.csv") if r["cases_2022"] != "NA"]


def ci(r, a, b, c):
    return float(r[a]) / M, float(r[b]) / M, float(r[c]) / M


# (grupo, rótulo, estimativa, li, ls, tipo); tipo: published | census | first | main | scenario | range
data = [
    ("As published", "ELSI-Brazil, authors' estimate (2019 population)",
     tri["ELSI-Brazil, authors' own national estimate (2019 population)"] / M, None, None, "published"),
    ("As published", "Machine-learning classification of ELSI-Brazil (2015–16)", 2.2, None, None, "published"),
    ("As published", "Delphi consensus (2019 population)",
     tri["Delphi consensus 8.5% (Aliberti 2025), 2019 population (as published)"] / M, None, None, "published"),
    ("As published", "GBD 2019, all ages (95% UI)", 1849981 / M, 1630321 / M, 2097631 / M, "published"),
    ("Applied to the 2022 Census", "ELSI-Brazil weighted age-specific rates",
     tri["ELSI-Brazil 2015-16, weighted age-specific (Bertola 2023) x Census 2022"] / M, None, None, "census"),
    ("Applied to the 2022 Census", "Delphi prevalence of 8.5%",
     tri["Delphi consensus 8.5% applied to Census 2022 population 60+"] / M, None, None, "census"),
    ("This analysis, 2022 Census", "Observed counts",
     *ci(pick(est, version="observed", estimand="typical study"), "cases_median", "cases_lb", "cases_ub"), "main"),
    ("This analysis, 2022 Census", "Authors' population estimates (primary)",
     *ci(pick(est, version="author_best", estimand="typical study"), "cases_median", "cases_lb", "cases_ub"), "main"),
    ("This analysis, 2022 Census", "Primary, population-averaged",
     *ci(pick(est, version="author_best", estimand="population-averaged"), "cases_median", "cases_lb", "cases_ub"), "main"),
    ("This analysis, 2022 Census", "All eight populations (sensitivity)",
     *ci(pick(sens, set="all eight populations", version="author_best"), "cases_median", "cases_lb", "cases_ub"), "main"),
    ("This analysis, 2022 Census", "Primary, Bayesian (95% CrI)",
     *ci(pick(bay, version="author_best", estimand="typical study"), "cases_median", "cases_lb", "cases_ub"), "main"),
    ("Bias-analysis scenarios, 2022 Census", "Two-phase surveys corrected, NPV 1.00",
     *ci(qba["two-phase corrected, NPV = 1.00"], "cases_2022_median", "cases_2022_lb", "cases_2022_ub"), "scenario"),
    ("Bias-analysis scenarios, 2022 Census", "Two-phase surveys corrected, NPV 0.97",
     *ci(qba["two-phase corrected, NPV = 0.97"], "cases_2022_median", "cases_2022_lb", "cases_2022_ub"), "scenario"),
    ("Bias-analysis scenarios, 2022 Census", "Two-phase surveys corrected, NPV 0.95",
     *ci(qba["two-phase corrected, NPV = 0.95"], "cases_2022_median", "cases_2022_lb", "cases_2022_ub"), "scenario"),
    ("Bias-analysis scenarios, 2022 Census", "Probabilistic bias analysis of two-phase surveys",
     *ci(pick(pba, variant="partial_pooling", quantity="total error (verification and random)"), "cases_median", "cases_lb", "cases_ub"),
     "scenario"),
    ("Bias-analysis scenarios, 2022 Census", "ELSI-Brazil algorithm misclassification (range)",
     None, min(elsi_grid) / M, max(elsi_grid) / M, "range"),
]

# Linhas: cabeçalho de grupo + itens, de cima para baixo
layout, y, last = [], 0.0, None
for g, lab, x, lo, hi, kind in data:
    if g != last:
        if last is not None:
            y += 0.55
        layout.append(("header", g, y))
        y += 1.0
        last = g
    layout.append(("item", (lab, x, lo, hi, kind), y))
    y += 1.0
n_rows = y

base = fs.apply_style(180)
W = 180 * fs.MM
H = (n_rows * 0.165 + 0.55)          # ~4,2 mm por linha + margens
fig = plt.figure(figsize=(W, H))
# Coluna de rótulos à esquerda (texto), painel no meio, coluna de valores à direita
ax = fig.add_axes([0.385, 0.42 / H, 0.43, 1 - 0.55 / H])
ax.set_ylim(n_rows - 0.4, -0.6)
# Eixo pelo dado: até o meio milhão acima do maior limite, para nenhum intervalo ser cortado
x_hi = max(v for _, _, x, lo, hi, _ in data for v in (x, hi) if v is not None)
x_max = (int(x_hi * 2) + 1) / 2
ticks = [1.0 + 0.5 * i for i in range(int((x_max - 1.0) * 2) + 1)]
ax.set_xlim(1.0, x_max + 0.1)
ax.set_xticks(ticks)
ax.set_xticklabels([f"{t:.1f}" for t in ticks])
ax.set_yticks([])
ax.spines["left"].set_visible(False)
fs.grid_x_only(ax)
ax.set_xlabel("People with dementia (millions)")

colour = {"published": fs.INK_SOFT, "census": fs.INK_SOFT, "first": fs.INK_SOFT,
          "main": fs.NAVY, "scenario": fs.CORAL, "range": fs.CORAL}
drawn = []
x_lab = -0.88   # em fração do eixo: início da coluna de rótulos
x_val = 1.03
for kind_row, payload, yy in layout:
    if kind_row == "header":
        ax.text(x_lab, yy, payload, transform=ax.get_yaxis_transform(), ha="left", va="center",
                fontsize=base, fontweight="bold", color=fs.INK, clip_on=False)
        continue
    lab, x, lo, hi, kind = payload
    c = colour[kind]
    ax.text(x_lab + 0.03, yy, lab, transform=ax.get_yaxis_transform(), ha="left", va="center",
            fontsize=base - 0.5, color=fs.INK, clip_on=False)
    if kind == "range":
        ax.plot([lo, hi], [yy, yy], color=c, lw=3.2, alpha=0.35, solid_capstyle="butt", zorder=2)
        txt = f"{lo:.2f} to {hi:.2f}"
    else:
        if lo is not None:
            ax.plot([lo, hi], [yy, yy], color=c, lw=0.9, zorder=2)
        face = "white" if kind == "published" else c
        ax.plot([x], [yy], marker="o", ms=4.2, mfc=face, mec=c, mew=0.9, ls="none", zorder=3)
        txt = f"{x:.2f}" + (f" ({lo:.2f} to {hi:.2f})" if lo is not None else "")
    ax.text(x_val, yy, txt, transform=ax.get_yaxis_transform(), ha="left", va="center",
            fontsize=base - 0.5, color=fs.INK, clip_on=False)
    drawn.append({"group": [g for g, l, *_ in data if l == lab][0], "label": lab, "kind": kind,
                  "estimate_millions": "" if x is None else round(x, 4),
                  "lower_millions": "" if lo is None else round(lo, 4),
                  "upper_millions": "" if hi is None else round(hi, 4)})

ax.text(x_val, -0.95, "Millions (95% interval)", transform=ax.get_yaxis_transform(), ha="left",
        va="center", fontsize=base - 0.5, color=fs.INK_SOFT, clip_on=False)

for ext in ("pdf", "png"):
    # sem data de criacao no PDF, para que a reexecucao gere o mesmo arquivo byte a byte
    fig.savefig(FIG / f"national_estimates.{ext}", dpi=600, metadata={"CreationDate": None} if ext == "pdf" else None)
plt.close(fig)

with open(FIG / "national_estimates_data.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, list(drawn[0].keys()))
    w.writeheader()
    w.writerows(drawn)
print("figure_style.py sha256:", hashlib.sha256(STYLE.read_bytes()).hexdigest()[:16], "| linhas:", len(drawn),
      "| altura (pol):", round(H, 2))
