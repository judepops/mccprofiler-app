#!/usr/bin/env bash
# Re-run the whole audit, in dependency order, after the feature space changes.
#
#   conda activate cd4env
#   bash audit/scripts/run_all.sh [outdir]
#
# Everything here reads the BUILT STORE, so rebuild it first:
#   conda activate cd4env && python backend/scripts/build_store.py
#   conda activate mccapp && python backend/scripts/verify_store.py
#
# verify_store.py will FAIL after a feature change, and that is correct: its
# EXPECTED constants encode the current numbers. Update them deliberately,
# never silently, and only once you have read why each one moved.
set -uo pipefail
OUT="${1:-audit/scripts/output}"
mkdir -p "$OUT"
cd "$(dirname "$0")/../.." || exit 1

run () {                       # run <script> <description>
  local f="audit/scripts/$1.py" log="$OUT/$1.txt"
  printf '\n=== %-38s %s\n' "$1" "$2"
  if python "$f" > "$log" 2>&1; then
    printf '    ok, %s lines -> %s\n' "$(wc -l < "$log")" "$log"
  else
    printf '    FAILED, see %s\n' "$log"; tail -3 "$log" | sed 's/^/      /'
  fi
}

echo "audit re-run, $(date -u +%Y-%m-%dT%H:%M)"
echo "store built: $(python -c "import json;print(json.load(open('store/manifest.json'))['built'][:19])" 2>/dev/null || echo unknown)"

echo; echo "--- 1. the substrate: what is in it and how far to trust it ---"
run diagnose_feature_novelty        "do gene properties explain the features?"
run diagnose_dimension_trust        "which components rest on reproducible features"
run diagnose_pc_names               "how much of each axis name is true"

echo; echo "--- 2. Aim 2: is the structure categorical? ---"
run experiment_cluster_search       "16 conditions, four cluster tests each"
run experiment_rotate_axes          "would rotation make the axes nameable"

echo; echo "--- 3. Aim 3: external sets, in confound order ---"
run diagnose_external_structure_corrected  "displacement, magnitude-corrected"
run diagnose_external_density_stratified   "+ gene density"
run diagnose_chromosome_confound           "+ chromosome"
run diagnose_external_power                "minimum detectable effect"
run diagnose_external_axis_rank            "where each set leans"
run diagnose_definition_type               "does definition type predict displacement"

echo; echo "--- 4. predictive claims, which need the length control ---"
run diagnose_loeuf_confounds        "length, expression, CpG vs the shape substrate"

echo; echo "--- 5. checks that were once ad hoc ---"
run diagnose_adhoc_checks           "coverage, taxonomy robustness, DBD split, chr19"

echo; echo "--- 6. figures ---"
run make_audit_figures              "six composites into audit/figures/"

echo; echo "done. outputs in $OUT"
echo "Now re-read audit/CURRENT_FINDINGS.md against these: every component count,"
echo "variance share and effect size in it will have moved."
