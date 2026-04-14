#!/usr/bin/env bash

set -euo pipefail

# ============================================================
# PREP_LIGANDS CLI TOOL
# ============================================================

# -----------------------------
# Default values
# -----------------------------
PH=7.4
CPU=""
CPU_PERCENT=0.8
OUTDIR="prep_out"

# -----------------------------
# Argument parsing
# -----------------------------
usage() {
    echo "Usage: prep_ligands --input FILE [--ph FLOAT] [--cpu INT | --cpu_percent FLOAT] --outdir DIR"
    exit 1
}

while [[ $# -gt 0 ]]; do
    case $1 in
        --input) INPUT="$2"; shift 2 ;;
        --ph) PH="$2"; shift 2 ;;
        --cpu) CPU="$2"; shift 2 ;;
        --cpu_percent) CPU_PERCENT="$2"; shift 2 ;;
        --outdir) OUTDIR="$2"; shift 2 ;;
        *) usage ;;
    esac
done

[[ -z "${INPUT:-}" ]] && usage

# -----------------------------
# Directory structure
# -----------------------------
SMI_DIR="$OUTDIR/smiles"
PROT_DIR="$OUTDIR/protonation"
CONF_DIR="$OUTDIR/3D"
PDBQT_DIR="$OUTDIR/pdbqt"
LOG_DIR="$OUTDIR/logs"

mkdir -p "$SMI_DIR" "$PROT_DIR" "$CONF_DIR" "$PDBQT_DIR" "$LOG_DIR"

# -----------------------------
# CPU handling
# -----------------------------
TOTAL_CPU=$(nproc)

if [[ -n "$CPU" ]]; then
    JOBS="$CPU"
else
    JOBS=$(python3 - <<EOF
import math
print(max(1, int($TOTAL_CPU * $CPU_PERCENT)))
EOF
)
fi

echo "[INFO] CPUs: $TOTAL_CPU | Using: $JOBS"

# ============================================================
# STEP 1: Input processing
# ============================================================

echo "[INFO] Step 1: Processing input..."

python3 py/input_handler.py \
    --input "$INPUT" \
    --out_smi "$SMI_DIR/master.smi" \
    --id_map "$OUTDIR/id_map.csv" \
    > "$LOG_DIR/input.log" 2>&1

# ============================================================
# STEP 2: Protonation
# ============================================================

echo "[INFO] Step 2: Protonation..."

bash lib/protonate.sh "$SMI_DIR/master.smi" "$PROT_DIR" "$PH" "$JOBS"

# ============================================================
# STEP 3: 3D generation
# ============================================================

echo "[INFO] Step 3: 3D conformers..."

find "$PROT_DIR/protonated_smi" -name "*.smi" | \
parallel -j "$JOBS" '
python3 py/generate_3d.py --input {} --outdir '"$CONF_DIR"' > '"$LOG_DIR"'/3D_{/}.log 2>&1 || true
'

# ============================================================
# STEP 4: PDBQT conversion
# ============================================================

echo "[INFO] Step 4: PDBQT conversion..."

find "$CONF_DIR" -name "*_3D.sdf" | \
parallel -j "$JOBS" '
base=$(basename {} _3D.sdf)
mk_prepare_ligand.py -i {} -o '"$PDBQT_DIR"'/$base.pdbqt \
> '"$LOG_DIR"'/meeko_$base.log 2>&1
'

echo "[SUCCESS] prep_ligands complete"