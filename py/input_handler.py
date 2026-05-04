#!/usr/bin/env python3

import argparse
import csv
import re
import sys
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

# ============================================================
# CLI
# ============================================================

parser = argparse.ArgumentParser(
    description="Universal molecule input handler (robust, HTVS-ready)"
)
parser.add_argument("--input", required=True)
parser.add_argument("--out_smi", required=True)
parser.add_argument("--id_map", required=True)
args = parser.parse_args()

ext = args.input.split(".")[-1].lower()

count, fail = 0, 0

# ============================================================
# Helpers
# ============================================================

def detect_source(raw_id):
    if not raw_id:
        return ("Custom", "")

    raw_id = str(raw_id).strip()

    if re.match(r'^ZINC', raw_id, re.I):
        return ("ZINC", raw_id)

    if re.match(r'^CHEMBL', raw_id, re.I):
        return ("ChEMBL", raw_id)

    if raw_id.isdigit():
        return ("PubChem", raw_id)

    return ("Custom", raw_id)


def clean_mol(mol):
    """Standard cleanup pipeline"""
    mol = rdMolStandardize.Cleanup(mol)
    mol = rdMolStandardize.FragmentParent(mol)
    Chem.SanitizeMol(mol)
    return mol


def detect_columns(headers):
    """Improved CSV column detection"""

    smiles_col = None
    id_col = None

    id_priority = [
        "compound_cid",
        "compound cid",
        "cid",
        "pubchem_cid",
        "pubchem cid",
        "pubchem id",
        "pubchem"
        "zinc_id",
        "zinc id",
        "zinc",
        "chembl_id",
        "chembl id",
        "chembl",
        "id"
    ]

    smiles_priority = [
        "smiles",
        "canonical_smiles",
        "canonical smiles",
        "smi"
    ]

    lower_map = {h.lower(): h for h in headers}


     # SMILES detection
    for key in smiles_priority:
        if key in lower_map:
            smiles_col = lower_map[key]
            break

    # ID detection
    for key in id_priority:
        if key in lower_map:
            id_col = lower_map[key]
            break

    # Fallback fuzzy smiles
    if smiles_col is None:
        for h in headers:
            if "smile" in h.lower():
                smiles_col = h
                break

    # Fallback ID = first column
    if id_col is None:
        id_col = headers[0]

    return smiles_col, id_col


# ============================================================
# Processing
# ============================================================

try:
    with open(args.out_smi, "w") as out_smi, open(args.id_map, "w") as map_file:

        map_file.write("Ligand,Source,SourceID\n")

        # ------------------------------------------------------------
        # SDF
        # ------------------------------------------------------------
        if ext == "sdf":

            suppl = Chem.SDMolSupplier(args.input, sanitize=True, removeHs=False)

            for i, mol in enumerate(suppl):

                if mol is None:
                    fail += 1
                    continue

                try:
                    mol = clean_mol(mol)

                    smi = Chem.MolToSmiles(mol, isomericSmiles=True)

                    raw_id = mol.GetProp("_Name") if mol.HasProp("_Name") else f"mol_{i+1}"
                    source, sid = detect_source(raw_id)

                    count += 1
                    name = f"mol_{count}"

                    out_smi.write(f"{smi} {name}\n")
                    map_file.write(f"{name},{source},{sid}\n")

                except Exception as e:
                    print(f"[ERROR SDF] {e}", file=sys.stderr)
                    fail += 1

        # ------------------------------------------------------------
        # SMI
        # ------------------------------------------------------------
        elif ext == "smi":

            with open(args.input) as f:
                for i, line in enumerate(f):

                    parts = line.strip().split()
                    if not parts:
                        continue

                    smi = parts[0]
                    raw_id = parts[1] if len(parts) > 1 else f"mol_{i+1}"

                    try:
                        mol = Chem.MolFromSmiles(smi)

                        if mol is None:
                            fail += 1
                            continue

                        mol = clean_mol(mol)

                        smi = Chem.MolToSmiles(mol, isomericSmiles=True)

                        source, sid = detect_source(raw_id)

                        count += 1
                        name = f"mol_{count}"

                        out_smi.write(f"{smi} {name}\n")
                        map_file.write(f"{name},{source},{sid}\n")

                    except Exception as e:
                        print(f"[ERROR SMI] {e}", file=sys.stderr)
                        fail += 1

        # ------------------------------------------------------------
        # CSV
        # ------------------------------------------------------------
        elif ext == "csv":

            with open(args.input) as f:
                sample = f.read(5000)

            try:
                dialect = csv.Sniffer().sniff(sample)
            except:
                dialect = csv.excel

            with open(args.input) as f:
                reader = csv.DictReader(f, dialect=dialect)

                headers = reader.fieldnames

                if not headers:
                    raise ValueError("CSV has no headers")

                smiles_col, id_col = detect_columns(headers)

                if smiles_col is None:
                    raise ValueError("No SMILES column detected")

                for i, row in enumerate(reader):

                    smi = str(row.get(smiles_col, "")).strip()
                    raw_id = row.get(id_col, f"mol_{i+1}")

                    if not smi:
                        continue

                    try:
                        mol = Chem.MolFromSmiles(smi)

                        if mol is None:
                            fail += 1
                            continue

                        mol = clean_mol(mol)

                        smi = Chem.MolToSmiles(mol, isomericSmiles=True)

                        source, sid = detect_source(raw_id)

                        count += 1
                        name = f"mol_{count}"

                        out_smi.write(f"{smi} {name}\n")
                        map_file.write(f"{name},{source},{sid}\n")

                    except Exception as e:
                        print(f"[ERROR CSV] {e}", file=sys.stderr)
                        fail += 1

        else:
            raise ValueError(f"Unsupported file format: {ext}")

# ============================================================
# FINAL
# ============================================================

except Exception as e:
    print(f"[FATAL ERROR] {e}", file=sys.stderr)
    sys.exit(1)

print("================================")
print(f"Processed: {count}")
print(f"Failed: {fail}")
print("================================")

if count == 0:
    sys.exit(2)
