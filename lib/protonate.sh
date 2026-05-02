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

TOTAL=$(wc -l < "$INPUT_SMI")
echo "[INFO] Total ligands: $TOTAL"

# -----------------------------
# SPLIT INPUT
# -----------------------------
split -n l/$JOBS "$INPUT_SMI" "$SPLIT_DIR/chunk_"

# -----------------------------
# PROTONATION (SAFE VERSION)
# -----------------------------
echo "[INFO] Running Dimorphite-DL..."

find "$SPLIT_DIR" -name "chunk_*" | parallel -j "$JOBS" '
chunk={}
base=$(basename "$chunk")

outfile='"$OUT_DIR"'/${base}_pH.smi
logfile='"$LOG_DIR"'/${base}.log

# Temporary files
tagged_in='"$OUT_DIR"'/${base}_tagged.smi
prot_out='"$OUT_DIR"'/${base}_prot_raw.smi

# -----------------------------------
# STEP 1: Create tagged input
# Format: SMILES<TAB>ID
# -----------------------------------
awk "{print \$1\"\t\"\$2}" "$chunk" > "$tagged_in"

# -----------------------------------
# STEP 2: Protonate (SMILES only)
# -----------------------------------
cut -f1 "$tagged_in" | dimorphite_dl \
    --ph_min '"$PH"' \
    --ph_max '"$PH"' \
    --max_variants 1 \
    > "$prot_out" 2>> "$logfile"

# -----------------------------------
# STEP 3: Safe recombination
# -----------------------------------
paste "$prot_out" <(cut -f2 "$tagged_in") > "$outfile"

# Cleanup temp
rm -f "$tagged_in" "$prot_out"
'

# -----------------------------
# MERGE OUTPUT
# -----------------------------
cat "$OUT_DIR"/*_pH.smi | sed '/^$/d' > "$FINAL_FILE"

if [[ ! -s "$FINAL_FILE" ]]; then
    echo "[ERROR] Protonation produced no output"
    exit 1
fi

PROT_COUNT=$(wc -l < "$FINAL_FILE")
echo "[INFO] Protonated molecules: $PROT_COUNT"

# -----------------------------
# CREATE PER-LIGAND FILES
# -----------------------------
echo "[INFO] Creating per-ligand SMILES..."

FAIL_SPLIT=0

while read -r smi name; do
    if [[ -z "$smi" || -z "$name" ]]; then
        ((FAIL_SPLIT++))
        continue
    fi
    echo "$smi $name" > "$NEXT_INPUT/${name}.smi"
done < "$FINAL_FILE"

GENERATED=$(ls "$NEXT_INPUT" | wc -l)

echo "[INFO] Generated per-ligand files: $GENERATED"

# -----------------------------
# FINAL CHECKS
# -----------------------------
if [[ "$GENERATED" -eq 0 ]]; then
    echo "[ERROR] No valid protonated ligands generated"
    exit 1
fi

if [[ "$GENERATED" -lt "$PROT_COUNT" ]]; then
    echo "[WARNING] Some ligands failed during splitting: $FAIL_SPLIT"
fi

echo "[SUCCESS] Protonation complete"
