#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: $0 ANALYSIS_EXECUTABLE PYTHIA_FILE_LIST [RUN_NAME]" >&2
  echo "Optional environment: ROOT_BIN, MAX_NOON_ENTRIES, NOON_SEED" >&2
  exit 2
}

[[ $# -eq 2 || $# -eq 3 ]] || usage
analysis_bin="$1"
file_list="$2"
run_name="${3:-run_$(date +%Y%m%d_%H%M%S)}"
max_noon_entries="${MAX_NOON_ENTRIES:--1}"
noon_seed="${NOON_SEED:-12345}"

package_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
if [[ -n "${ROOT_BIN:-}" ]]; then
  root_bin="$ROOT_BIN"
else
  root_bin="$(command -v root || true)"
fi

[[ -x "$analysis_bin" ]] || { echo "Analysis executable not found: $analysis_bin" >&2; exit 1; }
[[ -f "$file_list" ]] || { echo "PYTHIA file list not found: $file_list" >&2; exit 1; }
[[ -n "$root_bin" && -x "$root_bin" ]] || { echo "ROOT not found; activate ROOT or set ROOT_BIN" >&2; exit 1; }

output_dir="$package_dir/outputs/$run_name"
[[ ! -e "$output_dir" ]] || { echo "Output already exists: $output_dir" >&2; exit 1; }
mkdir -p -- "$output_dir"

weights="$package_dir/inputs/reweighting_factors_ep10x100_test17-18_updated_analysis_macro.txt"
flux="$package_dir/inputs/noon_flux.root"
pythia_output="$output_dir/pythia_weighted.root"
pythia_qa="$output_dir/pythia_qa"
noon_output="$output_dir/n00n_afterburner"
noon_qa="$output_dir/n00n_qa"

cpp_quote() {
  local value="$1"
  [[ "$value" != *$'\n'* && "$value" != *$'\r'* ]] || return 1
  value="${value//\\/\\\\}"
  value="${value//\"/\\\"}"
  printf '"%s"' "$value"
}

echo "[1/5] PYTHIA analysis with n00n flux weights and STAR-like response"
"$analysis_bin" py8 "$file_list" "$pythia_output" "$weights"

echo "[2/5] PYTHIA/detector QA"
pythia_call="$package_dir/pythia/plot_pythia_post_generation_essentials_macro2.C+($(cpp_quote "$pythia_output"),$(cpp_quote "$flux"),$(cpp_quote "$pythia_qa"))"
"$root_bin" -l -b -q "$pythia_call"

echo "[3/5] n00n neutron afterburner"
mkdir -p -- "$noon_output"
(
  cd -- "$package_dir/n00n"
  "$root_bin" -l -b -q 'NeutronGenerator.cxx+' \
    "run_pythia_noon.C+($(cpp_quote "$pythia_output"),$max_noon_entries,$(cpp_quote "$noon_output"),$noon_seed,0.,5.,0,100000)"
)

neutron_output="$noon_output/pythia_noon.root"
echo "[4/5] n00n neutron QA"
(
  cd -- "$package_dir/n00n"
  "$root_bin" -l -b -q \
    "plot_pythia_noon_qa.C+($(cpp_quote "$neutron_output"),$(cpp_quote "$noon_qa"),$(cpp_quote "$flux"))"
)

echo "[5/5] Persisted-output validation"
"$root_bin" -l -b -q \
  "$package_dir/scripts/validate_chain.C+($(cpp_quote "$pythia_output"),$(cpp_quote "$neutron_output"))"

echo "Chain completed successfully."
echo "PYTHIA output: $pythia_output"
echo "Neutron output: $neutron_output"
echo "QA directories: $pythia_qa and $noon_qa"

