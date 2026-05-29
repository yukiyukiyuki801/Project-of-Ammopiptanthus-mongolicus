#!/usr/bin/env bash
set -euo pipefail

AMMO_DIR="${AMMO_DIR:-/home/liao/DataVol1/personal2/sun/ammo}"
RONA_DIR="${RONA_DIR:-$AMMO_DIR/rona}"
PYRONA="${PYRONA:-pyRona.py}"

cd "$RONA_DIR"

required=(
  "ammo_current.txt"
  "associations.csv"
  "ammo_outliers_3overlap.recode.lfmm_imputed.lfmm"
  "var_name.txt"
)

for f in "${required[@]}"; do
  [[ -s "$f" ]] || { echo "Missing required file: $RONA_DIR/$f" >&2; exit 1; }
done

run_pyrona() {
  local scenario="$1"
  local period="$2"
  local future="$3"
  local pdf_out="$4"
  local txt_out="$5"

  [[ -s "$future" ]] || { echo "Missing future climate file: $RONA_DIR/$future" >&2; exit 1; }

  "$PYRONA" lfmm \
    -pc ammo_current.txt \
    -fc "$future" \
    -out "$pdf_out" \
    -P 1 \
    -assoc associations.csv \
    -geno ammo_outliers_3overlap.recode.lfmm_imputed.lfmm \
    -covar_names var_name.txt \
    -covars 5 \
    > "$txt_out"
}

run_pyrona "ssp126" "2040" "ammo_rona_ssp126_2040-2070.txt" "ammo_lfmm2_rona_126_2040.pdf" "ssp126_2040_rona.txt"
run_pyrona "ssp126" "2070" "ammo_rona_ssp126_2070-2100.txt" "ammo_lfmm2_rona_126_2070.pdf" "ssp126_2070_rona.txt"
run_pyrona "ssp585" "2040" "ammo_rona_ssp585_2040-2070.txt" "ammo_lfmm2_rona_585_2040.pdf" "ssp585_2040_rona.txt"
run_pyrona "ssp585" "2070" "ammo_rona_ssp585_2070-2100.txt" "ammo_lfmm2_rona_585_2070.pdf" "ssp585_2070_rona.txt"

echo "pyRona analyses completed in $RONA_DIR"
