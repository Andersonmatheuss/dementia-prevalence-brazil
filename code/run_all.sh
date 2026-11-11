#!/usr/bin/env bash
# Reproduces every table, figure and supplementary table of the manuscript from data/, in one pass.
#
# The brms fits (code/02c, 03i, 03l) reuse the cached models in output/intermediate/models/ while data
# and model code are unchanged (options(brms.file_refit = "on_change") in code/helpers/setup.R).
#
# Usage (from any directory):
#   code/run_all.sh            one pass
#   code/run_all.sh --check    two passes; every output must be byte-identical between them
# Environment:
#   PY=/path/to/python         Python interpreter (default: python3); needs matplotlib (code/requirements.txt)
set -euo pipefail
cd "$(dirname "$0")/.."
export PY=${PY:-python3}
LOGS=$(mktemp -d)

# Order matters: every step reads only files written by the steps above it (or by data/).
STEPS=(
  code/02a_check_raw_chain.R
  code/02b_reproduce_frequentist.R
  code/02c_reproduce_bayesian.R
  code/03a_fit_corrected_frequentist.R
  code/03b_project_ibge_2024.R
  code/03c_fit_corrected_subgroups.R
  code/03d_compare_external_and_qba.R
  code/03e_qba_two_phase_documented.R
  code/03f_add_epifloripa.R
  code/03g_fit_author_best_estimates.R
  code/03h_value_of_information.R
  code/03i_fit_bayesian_corrected.R
  code/03j_fit_primary_set.R
  code/03k_compare_pieta.R
  code/03l_fit_bayesian_primary_set.R
  code/03m_burden_states.R
  code/03n_prevalence_trend_scenarios.R
  code/03o_reconciliation_bridge.R
  code/03p_education_standardisation.R
  code/03q_age_band_contributions.R
  code/03r_loo_prior_sensitivity_primary.R
  code/03s_pba_ppv_unexamined.R
  code/04a_build_supplement_tables.R
  code/05a_plot_national_estimates.py
  code/05c_plot_manuscript_figures.py
  code/05d_plot_figures_manuscript_layout.R
  code/06a_check_reproduction.R
  code/06b_compare_underdiagnosis.R
)

run_step () {
  local s=$1 log="$LOGS/$(basename "$1").log" cmd
  case $s in
    *.py) cmd=("$PY" -I "$s") ;;
    *.R)  cmd=(Rscript "$s") ;;
  esac
  echo "  $(date +%H:%M:%S) $s"
  if ! "${cmd[@]}" > "$log" 2>&1; then
    echo "FAILED: $s (log: $log)"; tail -20 "$log"; exit 1
  fi
}

run_chain () {
  mkdir -p output/intermediate output/figures output/tables supplementary
  for s in "${STEPS[@]}"; do run_step "$s"; done
}

# Files the chain writes. The cached brms fits (output/intermediate/models) are inputs, not outputs.
list_outputs () {
  {
    find output data/processed -type f ! -path 'output/intermediate/models/*' ! -name '.*'
    ls supplementary/supplement_tables.md
  } | sort
}

hash_outputs () { list_outputs | while read -r f; do shasum -a 256 "$f"; done; }

if [[ "${1:-}" != "--check" ]]; then
  echo "run_all: one pass"
  run_chain
  echo "done: $(list_outputs | wc -l | tr -d ' ') outputs (logs in $LOGS)"
  exit 0
fi

# --check: the chain and its inputs must not change between the two passes.
guard () { find code data/raw output/intermediate/models -type f \
             -exec shasum -a 256 {} + | sort | shasum -a 256; }
before=$(guard)
echo "PASS 1"; run_chain; hash_outputs > "$LOGS/pass1.sha"
echo "PASS 2"; run_chain; hash_outputs > "$LOGS/pass2.sha"
if [[ "$(guard)" != "$before" ]]; then
  echo "ABORTED: code or inputs changed during the check; run it again without editing."; exit 2
fi
if diff "$LOGS/pass1.sha" "$LOGS/pass2.sha"; then
  echo "CHECK: all $(wc -l < "$LOGS/pass1.sha" | tr -d ' ') outputs identical between the two passes."
else
  echo "CHECK FAILED: outputs above differ between the two passes."; exit 1
fi
