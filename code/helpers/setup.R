# 00_setup.R — semente, caminhos e funcoes comuns do projeto.
# O run_all.R o executa primeiro; cada script o carrega de novo para tambem rodar sozinho
# (Rscript code/NNx_nome.R). Caminhos SEMPRE relativos ao projeto (here).

set.seed(20261019)  # semente do projeto; script estocastico fixa a sua antes de sortear


# Caminhos canonicos. Nomes em ingles porque o alvo sao journals
# internacionais: o repositorio vira material suplementar do artigo, e figura
# e tabela sobem com o nome que ja aparece no texto.
dir_raw  <- file.path("data", "raw")        
dir_proc <- file.path("data", "processed")
dir_fig  <- file.path("output", "figures")
dir_tab  <- file.path("output", "tables")

# Nomes usados pelos scripts da reanalise (vieram da auditoria de 19/10/2026).
RAW            <- dir_raw
OCT_OUT        <- file.path(RAW, "Dementia_prevalence", "dementia_outputs")   # saidas de out/2026 (sem o script)
OCT_IBGE       <- file.path(RAW, "Dementia_prevalence", "ibge_sidra_inputs")  # Censo 2022 e projecao 2018 (SIDRA)
GH_DIR         <- file.path(RAW, "dementia-prevalence-brazil-addc1b1")        # codigo de mar/2026 (GitHub)
RAW_EXTRACTION <- file.path(RAW, "extraction")             # extracao e conferencia contra as fontes
# Extracao conferida no texto integral (19/10/2026): estratos com tabela e pagina de cada valor, e os
# numeros por estudo (estimativas publicadas, rastreio). A analise corrigida le daqui; a reproducao
# da primeira versao (scripts 03 e 04) continua lendo a planilha de out/2026 (OCT_OUT).
VERIFIED_STRATA    <- file.path(RAW_EXTRACTION, "fulltext_strata.csv")
VERIFIED_ESTIMATES <- file.path(RAW_EXTRACTION, "fulltext_study_estimates.csv")
RAW_IBGE_2024  <- file.path(RAW, "ibge_projection_2024")   # projecoes do IBGE, revisao 2024
DERIVED  <- file.path("output", "intermediate")
OUT_TAB  <- dir_tab
OUT_FIG  <- dir_fig
OUT_MOD  <- file.path("output", "intermediate", "models")       # fora do git: o brmsfit guarda a tabela em fit$data
OUT_DIAG <- file.path("output", "intermediate", "diagnostics")
for (d in c(DERIVED, OUT_TAB, OUT_FIG, OUT_MOD, OUT_DIAG)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# O brms reaproveita o ajuste gravado em OUT_MOD; sem isto ("never", o padrao) ele o reaproveita
# mesmo com dados novos, e uma correcao de extracao nao chegaria ao modelo bayesiano.
options(brms.file_refit = "on_change")

# Conjuntos de analise (decisao do PI em 19/10/2026: "a analise principal deve ir sem ter overllaping").
# Regras: o ELSI-Brazil cobre o territorio de todos os estudos locais e e o comparador externo, entao fica
# so como comparador; em Sao Paulo, onde Bottino (PROTER, distritos e datas nao informados) e Scazufca
# (SPAH, Butanta, 2003-2005) podem ter participantes em comum, fica a amostra com base e periodo
# documentados (Scazufca); Correa Ribeiro (clientes de plano de saude) nao e amostra da populacao geral.
# Sensibilidades fixadas junto com a regra, antes de qualquer resultado com o conjunto novo.
POPULATIONS_ALL <- c("Bertola_2023", "Bottino_2008", "Cesar_2015", "EpiFloripa_2024", "Herrera_2002",
                     "Lopes_2012", "Ribeiro_2013", "Scazufca_2008")
PRIMARY_POPULATIONS <- c("Cesar_2015", "EpiFloripa_2024", "Herrera_2002", "Lopes_2012", "Scazufca_2008")
SENSITIVITY_SETS <- list(
  "all eight populations"       = POPULATIONS_ALL,
  "primary + Bottino"           = c(PRIMARY_POPULATIONS, "Bottino_2008"),
  "primary + Correa Ribeiro"    = c(PRIMARY_POPULATIONS, "Ribeiro_2013"),
  "primary + ELSI-Brazil"       = c(PRIMARY_POPULATIONS, "Bertola_2023"),
  "Bottino instead of Scazufca" = c(setdiff(PRIMARY_POPULATIONS, "Scazufca_2008"), "Bottino_2008"))
# Papel de cada populacao e o report_id correspondente em data/raw/extraction/extraction_master.csv
STUDY_ROLE <- c(Bertola_2023 = "external comparator", Bottino_2008 = "sensitivity", Ribeiro_2013 = "sensitivity",
                setNames(rep("primary", length(PRIMARY_POPULATIONS)), PRIMARY_POPULATIONS))[POPULATIONS_ALL]
STUDY_REPORT <- c(Bertola_2023 = "Bertola_2023", Bottino_2008 = "Bottino_2008", Cesar_2015 = "Cesar_2016",
                  EpiFloripa_2024 = "EpiFloripa_2024", Herrera_2002 = "Herrera_2002", Lopes_2012 = "Lopes_2012",
                  Ribeiro_2013 = "CorreaRibeiro_2013", Scazufca_2008 = "Scazufca_2008")[POPULATIONS_ALL]
stopifnot(all(PRIMARY_POPULATIONS %in% POPULATIONS_ALL), !any(c("Bertola_2023", "Bottino_2008", "Ribeiro_2013") %in% PRIMARY_POPULATIONS))

# Anos das projecoes (revisao 2024 do IBGE, que vai ate 2070); 2070 a pedido do PI em 19/10/2026.
PROJECTION_YEARS <- c(2030L, 2040L, 2050L, 2060L, 2070L)

# Agrupamento que fecha com o manuscrito da primeira versao: rma REML no logito, variancia binomial.
fit_re <- function(d) metafor::rma(yi = yi, vi = vi, data = d, method = "REML")
pred_p <- function(m) {
  p <- predict(m, transf = metafor::transf.ilogit)
  c(est = p$pred, lb = p$ci.lb, ub = p$ci.ub, pi_lb = p$pi.lb, pi_ub = p$pi.ub)
}

# Grafico: as Figuras 1 e 2 saem do matplotlib (code/05a e 05b), com o estilo
# canonico do grupo (figure_style.py do estudo dementia-death, lido no lugar).

# --- Reprodutibilidade -------------------------------------------------------
# Na PRIMEIRA vez em cada projeto, rodar uma vez no console:
#   renv::init()          # cria renv.lock com as versoes exatas
#   renv::snapshot()      # apos instalar/atualizar qualquer pacote
# O cache global (~/.cache/R/renv) evita rebaixar pacote ja usado noutro projeto.
#
# snapshot() so registra o que aparece em library() no codigo do projeto.
# type = "all" trava a biblioteca inteira: o lock passa de 300 KB e amarra o
# projeto a pacote que ele nao usa. Use so se for esse o objetivo.
