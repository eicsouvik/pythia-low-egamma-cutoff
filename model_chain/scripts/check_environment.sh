#!/usr/bin/env bash
set -euo pipefail

package_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"

fail=0
if [[ -n "${ROOT_BIN:-}" && -x "$ROOT_BIN" ]]; then
  root_bin="$ROOT_BIN"
else
  root_bin="$(command -v root || true)"
fi
if [[ -z "$root_bin" ]]; then
  echo "MISSING: ROOT executable (activate ROOT before running the chain)" >&2
  fail=1
else
  echo "ROOT: $root_bin"
  root_config="$(dirname -- "$root_bin")/root-config"
  [[ -x "$root_config" ]] && "$root_config" --version || true
fi

for path in \
  "$package_dir/inputs/noon_flux.root" \
  "$package_dir/inputs/reweighting_factors_ep10x100_test17-18_updated_analysis_macro.txt" \
  "$package_dir/n00n/NeutronGenerator.cxx" \
  "$package_dir/n00n/Data/ENDF_Au197.txt"; do
  if [[ ! -f "$path" ]]; then
    echo "MISSING: $path" >&2
    fail=1
  else
    echo "FOUND: $path"
  fi
done

if [[ $fail -ne 0 ]]; then
  exit 1
fi
echo "ROOT/n00n runtime inputs are available."
echo "The PYTHIA stage additionally requires the STAR software environment,"
echo "a compiled analysis executable, and a generator ROOT-file list."

