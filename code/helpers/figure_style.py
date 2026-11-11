"""Sistema de design das figuras, no padrão de periódico biomédico de alto impacto.

Três decisões estruturam tudo o que segue.

**1. A figura não carrega título.** Em Lancet, NEJM e JAMA o título e a nota
metodológica vivem na *legenda*, composta pela revista junto ao texto — nunca
desenhados dentro da arte. Título embutido é a marca mais reconhecível de figura
amadora, e além disso rouba área do painel e duplica o que a legenda dirá. As
legendas são geradas em ``output/figures/figure_captions.md``, prontas para
colar no manuscrito.

**2. O cromo é recessivo.** Sem grade vertical; grade horizontal em fio de
cabelo, atrás das marcas; duas espinhas apenas; marcas de escala curtas e para
fora. O leitor deve ver os dados, não o gráfico.

**3. A cor é verificada, não estimada.** A paleta é conferida com
``scripts/validate_palette.js`` do guia de visualização: banda de luminosidade,
piso de croma, separação sob deuteranopia/protanopia/tritanopia, piso de visão
normal e contraste contra a superfície.

Histórico, porque a versão anterior deste módulo afirmava validação que NÃO se
sustentava. O conjunto ``#1F4E79 / #E76F51 / #2A9D8F / #8C6D1F / #3F444A``
**reprovou** em três das cinco verificações:

* o grafite ``#3F444A`` tinha croma 0,012 — lia-se como cinza, não como cor de
  série, exatamente o defeito que motivara aposentar o ``#6C757D``;
* grafite contra azul-marinho dava ΔE 8,2 na **visão normal**, abaixo do piso de
  15: eram indistinguíveis mesmo para quem enxerga todas as cores. E esse par
  carregava Sudeste contra Centro-Oeste, que é justamente o contraste que o
  leitor precisa fazer nas figuras 7 e 10;
* o próprio azul-marinho ``#1F4E79`` ficava fora da banda de luminosidade
  (0,414) e abaixo do piso de croma (0,088).

A paleta atual é o subconjunto de Okabe-Ito **sem o rosa-arroxeado**, e passa em
todas as verificações com os cinco pares. Restam dois avisos de contraste
(laranja 2,19:1 e azul-céu 2,25:1), que o guia admite mediante codificação
secundária — atendida aqui por rótulo direto ou legenda com texto em todas as
figuras que usam essas duas cores, além do CSV que acompanha cada figura.
"""

from __future__ import annotations

import matplotlib
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Paleta — validada; ver o cabeçalho do módulo
# ---------------------------------------------------------------------------
INK = "#1A1A1A"          # texto principal
INK_SOFT = "#5A6068"     # texto secundário, anotações
RULE = "#3F444A"         # espinhas e marcas de escala
GRID = "#E7E9EC"         # grade horizontal, um passo da superfície
SURFACE = "#FFFFFF"
ABSENT = "#D9DEE3"       # dado ausente, explicitamente cinzento

# Sequência categórica FIXA. Nunca reciclada: uma 6ª série vira "outras",
# múltiplos pequenos ou codificação composta.
#
# Okabe-Ito sem o rosa-arroxeado. Conferida com --pairs all: banda de
# luminosidade, croma, separação sob deuteranopia/protanopia/tritanopia e piso
# de visão normal, todas aprovadas. NÃO altere um valor isolado: o piso vale
# para TODOS os pares, e trocar uma cor pode reprovar um par que não a envolve
# diretamente. Reexecute o validador a cada mudança.
SERIES = ["#0072B2", "#D55E00", "#009E73", "#E69F00", "#56B4E9"]
NAVY, CORAL, TEAL, OCHRE, SKY = SERIES
CHARCOAL = SKY          # nome antigo, mantido para não quebrar importações

# Rampa SEQUENCIAL de um hue só, claro→escuro. Proporções não têm ponto neutro
# natural, então uma escala divergente as descreveria mal: ela afirmaria um
# meio-termo que o dado não tem.
# Ancorada no azul de série, para que heatmap e mapa não falem um azul diferente
# do das linhas. Luminosidade estritamente decrescente — condição da rampa
# sequencial, verificada em ``tests/test_audit_regressions.py``.
SEQUENTIAL = ["#F1F6FA", "#D2E4F1", "#A8CBE4", "#77AED4", "#3F8FC0", "#0072B2", "#004E7A"]

# Papéis semânticos fixos, idênticos em todas as figuras
ROLE = {
    "specific": NAVY,
    "unspecified": CORAL,
    "any_mention": NAVY,
    "underlying": TEAL,
    "alzheimer": NAVY,
    "vascular": TEAL,
    "other_specified": OCHRE,
    "mixed": CHARCOAL,
}
REGION = {
    "North": TEAL,
    "Northeast": CORAL,
    "Southeast": NAVY,
    "South": OCHRE,
    "Central-West": CHARCOAL,
}
# Ordem canônica do IBGE, e a ordem em que o leitor brasileiro espera as regiões.
# Ordenar pelo nome em inglês produz Central-West, North, Northeast, South,
# Southeast — uma sequência que não corresponde a nada, nem geográfico nem
# estatístico, e que separa Sul de Sudeste. Toda figura agrupada por região usa
# esta lista; ``REGION`` já está declarada nesta ordem, mas ordenações derivadas
# de ``sorted()`` sobre o nome precisam desta chave explícita.
REGION_ORDER = list(REGION)

MM = 1 / 25.4
WIDTHS_MM = [90, 180]

# Corpo tipográfico por largura. A 90 mm a figura ocupa uma coluna e é reduzida
# na impressão, então o corpo sobe para permanecer acima do mínimo legível de
# ~6 pt; a 180 mm ela atravessa a página e o corpo cai.
FONT_PT = {90: 7.0, 180: 7.5}


def apply_style(width_mm: int = 180) -> float:
    """Aplica o estilo e devolve o corpo base em pontos."""
    base = FONT_PT.get(width_mm, 7.5)
    plt.rcParams.update({
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "savefig.transparent": False,
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"],
        "font.size": base,
        "text.color": INK,
        "axes.labelcolor": INK,
        "axes.labelsize": base,
        "axes.labelpad": 3.0,
        "axes.titlesize": base,
        "axes.edgecolor": RULE,
        "axes.linewidth": 0.6,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "axes.axisbelow": True,          # a grade fica ATRÁS das marcas
        "grid.color": GRID,
        "grid.linewidth": 0.4,
        "grid.linestyle": "-",           # nunca tracejada: lê-se como limiar
        "xtick.color": RULE,
        "ytick.color": RULE,
        "xtick.labelcolor": INK,
        "ytick.labelcolor": INK,
        "xtick.labelsize": base - 0.5,
        "ytick.labelsize": base - 0.5,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.major.size": 2.5,
        "ytick.major.size": 2.5,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.major.pad": 2.0,
        "ytick.major.pad": 2.0,
        "legend.frameon": False,
        "legend.fontsize": base - 0.5,
        "legend.handlelength": 1.4,
        "legend.handletextpad": 0.5,
        "legend.borderaxespad": 0.0,
        "legend.labelspacing": 0.35,
        "lines.linewidth": 1.5,
        "lines.solid_capstyle": "round",
        "patch.linewidth": 0.0,
    })
    return base


def grid_x_only(ax) -> None:
    """Grade só no eixo x (para barras horizontais)."""
    ax.grid(True, axis="x")
    ax.grid(False, axis="y")


def grid_y_only(ax) -> None:
    """Grade só no eixo y — o padrão para séries temporais."""
    ax.grid(True, axis="y")
    ax.grid(False, axis="x")


def no_grid(ax) -> None:
    ax.grid(False)


def panel_label(ax, letter: str, base: float, dx: float = -0.10, dy: float = 1.03) -> None:
    """Rótulo de painel (A, B, C…) no canto superior esquerdo, fora da área de dados."""
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=base + 1.5,
            fontweight="bold", va="bottom", ha="left", color=INK)


def direct_label(ax, x, y, text: str, colour: str, base: float,
                 dx: float = 0.006, va: str = "center") -> None:
    """Rótulo junto ao fim da série, dispensando caixa de legenda.

    O texto vai na cor da série apenas aqui, onde ele *é* a legenda; nas demais
    posições o texto usa tinta neutra.
    """
    # Peso NORMAL. "semibold" não existe em Helvetica Neue e o matplotlib cai
    # silenciosamente em 700, isto é, negrito — verificado no log da renderização
    # ("Failed to find font weight semibold, now using 700"). Rótulo em negrito
    # espalhado pela área de dados não é registro de periódico: a cor já
    # identifica a série, e o peso só acrescenta ruído.
    ax.annotate(text, xy=(x, y), xytext=(dx, 0), textcoords="offset fontsize",
                fontsize=base - 0.5, color=colour, va=va, ha="left")


def hide_spines(ax, keep: tuple[str, ...] = ("left", "bottom")) -> None:
    for side in ["top", "right", "left", "bottom"]:
        ax.spines[side].set_visible(side in keep)


def tidy_axis_limits(ax, values, pad_frac: float = 0.08, axis: str = "y",
                     floor: float | None = None, ceil: float | None = None) -> None:
    """Ajusta o limite aos dados, com folga — nunca 0-100 por reflexo.

    Uma proporção que varia de 63 % a 84 % desenhada num eixo de 0 a 100 gasta
    três quintos do painel em branco e achata justamente a variação que a figura
    existe para mostrar.
    """
    lo, hi = float(min(values)), float(max(values))
    span = hi - lo or max(abs(hi), 1.0)
    lo, hi = lo - pad_frac * span, hi + pad_frac * span
    if floor is not None:
        lo = max(lo, floor)
    if ceil is not None:
        hi = min(hi, ceil)
    (ax.set_ylim if axis == "y" else ax.set_xlim)(lo, hi)


def sequential_cmap(name: str = "study_sequential"):
    return matplotlib.colors.LinearSegmentedColormap.from_list(name, SEQUENTIAL)


def year_axis(ax, years, n_ticks: int = 6, right_pad_years: float = 0.0) -> None:
    """Marcas de ano SOMENTE dentro do intervalo observado.

    Estender o limite do eixo para acomodar rótulos junto à série faz o
    localizador automático gerar marcas além do último ano com dado — um eixo
    que chega a 2028 num estudo que termina em 2023 afirma uma cobertura que
    não existe. Aqui o limite cresce, mas as marcas param no dado.
    """
    y0, y1 = int(min(years)), int(max(years))
    span = y1 - y0
    step = max(1, int(round(span / max(n_ticks - 1, 1))))
    for nice in (1, 2, 5, 10, 20, 25):        # passos legíveis
        if step <= nice:
            step = nice
            break
    ticks = list(range(y0, y1 + 1, step))
    if ticks[-1] != y1 and (y1 - ticks[-1]) >= step / 2:
        ticks.append(y1)
    ax.set_xticks(ticks)
    ax.set_xlim(y0 - 0.4, y1 + right_pad_years + 0.4)


def _relative_luminance(rgb) -> float:
    """Luminância relativa da WCAG 2.1, com a linearização sRGB."""
    def channel(c: float) -> float:
        return c / 12.92 if c <= 0.03928 else ((c + 0.055) / 1.055) ** 2.4
    r, g, b = (channel(float(c)) for c in rgb[:3])
    return 0.2126 * r + 0.7152 * g + 0.0722 * b


def contrast_ratio(a, b) -> float:
    """Razão de contraste da WCAG entre duas cores. 4,5:1 é o piso de texto pequeno."""
    la, lb = _relative_luminance(a), _relative_luminance(b)
    hi, lo = max(la, lb), min(la, lb)
    return (hi + 0.05) / (lo + 0.05)


def text_on(face) -> str:
    """Tinta sobre um preenchimento: a que der MAIOR contraste, medido.

    Antes isto era um limiar fixo de luminância em 0,62, e o número estava
    errado: no passo de 76 % da rampa sequencial ele já mandava texto BRANCO,
    que ali rende 2,51:1 — abaixo do piso de 4,5:1 da WCAG para texto pequeno —
    enquanto a tinta escura no mesmo tom renderia 7,25:1. Um limiar de
    luminância não sabe disso porque não é a luminância que se lê, é o contraste
    entre as duas cores.

    Comparar as duas razões elimina o número mágico e se autocorrige se a rampa
    mudar. Não garante 4,5:1: perto do ponto de cruzamento nenhuma das duas
    tintas alcança o piso, e é a rampa que precisa ceder, não a regra. Coberto
    por ``tests/test_audit_regressions.py``.
    """
    white, ink = (1.0, 1.0, 1.0), (0.102, 0.102, 0.102)
    return "#FFFFFF" if contrast_ratio(face, white) > contrast_ratio(face, ink) else INK


def shade_period(ax, x0: float, x1: float, label: str | None, base: float,
                 y: float | None = None) -> None:
    """Faixa de período (p.ex. 2020-2022), discreta e rotulada no topo."""
    ax.axvspan(x0, x1, color=GRID, alpha=0.75, zorder=0, linewidth=0)
    if label:
        top = ax.get_ylim()[1] if y is None else y
        ax.text((x0 + x1) / 2, top, label, ha="center", va="top",
                fontsize=base - 1.5, color=INK_SOFT)
