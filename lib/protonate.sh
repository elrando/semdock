#!/usr/bin/env bash
set -euo pipefail

INPUT_SMI="$1"
OUT_BASE="$2"
PH="$3"
JOBS="$4"

# -----------------------------
# DIRECTORIES
# -----------------------------
SPLIT_DIR="$OUT_BASE/chunks"
OUT_DIR="$OUT_BASE/protonated_chunks"
LOG_DIR="$OUT_BASE/logs"
FINAL_FILE="$OUT_BASE/protonated_master.smi"
NEXT_INPUT="$OUT_BASE/protonated_smi"

rm -rf "$NEXT_INPUT"
mkdir -p "$SPLIT_DIR" "$OUT_DIR" "$LOG_DIR" "$NEXT_INPUT"

echo "[INFO] Protonation started"
echo "[INFO] Input: $INPUT_SMI"
echo "[INFO] pH: $PH"
echo "[INFO] Jobs: $JOBS"

# -----------------------------
# CHECK INPUT
# -----------------------------
if [[ ! -s "$INPUT_SMI" ]]; then
    echo "[ERROR] Empty or missing input SMILES"
    exit 1
fi

# -----------------------------
# SPLIT INPUT
# -----------------------------
split -n l/$JOBS "$INPUT_SMI" "$SPLIT_DIR/chunk_"

echo "[INFO] Total ligands:"
wc -l "$INPUT_SMI"

# -----------------------------
# PROTONATION (PARALLEL)
# -----------------------------
echo "[INFO] Running Dimorphite-DL..."

find "$SPLIT_DIR" -name "chunk_*" | parallel -j "$JOBS" '
chunk={}
base=$(basename "$chunk")

outfile='"$OUT_DIR"'/${base}_pH.smi
logfile='"$LOG_DIR"'/${base}.log
tmpfile='"$OUT_DIR"'/${base}_tmp.smi

dimorphite_dl "$chunk" \
    --ph_min '"$PH"' \
    --ph_max '"$PH"' \
    --max_variants 1 \
    > "$tmpfile" 2>> "$logfile"

paste -d " " "$tmpfile" <(awk "{print \$2}" "$chunk") > "$outfile"
rm -f "$tmpfile"
'

# -----------------------------
# MERGE
# -----------------------------
cat "$OUT_DIR"/*_pH.smi | sed '/^$/d' > "$FINAL_FILE"

echo "[INFO] Protonated molecules:"
wc -l "$FINAL_FILE"

# -----------------------------
# FIXED SPLITTING (IMPORTANT FIX)
# -----------------------------
echo "[INFO] Creating per-ligand SMILES..."

while read -r smi name; do
    [[ -z "$smi" || -z "$name" ]] && continue
    echo "$smi $name" > "$NEXT_INPUT/${name}.smi"
done < "$FINAL_FILE"

echo "[INFO] Generated files:"
ls "$NEXT_INPUT" | wc -l

echo "[SUCCESS] Protonation complete"
