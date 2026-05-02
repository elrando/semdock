#!/usr/bin/env python3

import argparse
import os
import sys
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

# ============================================================
# CLI
# ============================================================
parser = argparse.ArgumentParser(description="Generate 3D conformers from SMILES")
parser.add_argument("--input", required=True, help="Input .smi file (single ligand)")
parser.add_argument("--outdir", required=True, help="Output directory for SDF")

parser.add_argument("--nconf", type=int, default=12)
parser.add_argument("--keep", type=int, default=1)
parser.add_argument("--maxiters", type=int, default=500)
parser.add_argument("--seed", type=int, default=42)

args = parser.parse_args()

# ============================================================
# Validate input (IMPORTANT FIX)
# ============================================================
if not os.path.exists(args.input):
    print(f"[ERROR] Input file not found: {args.input}")
    sys.exit(1)

if os.path.getsize(args.input) == 0:
    print(f"[ERROR] Input file is empty: {args.input}")
    sys.exit(1)

# ============================================================
# Setup output
# ============================================================
os.makedirs(args.outdir, exist_ok=True)

infile = args.input
base = os.path.basename(infile).replace(".smi", "")
outfile = os.path.join(args.outdir, f"{base}_3D.sdf")

# ============================================================
# Main execution
# ============================================================
try:
    print("===================================")
    print("Ligand:", infile)

    # -----------------------------
    # READ SMILES
    # -----------------------------
    with open(infile) as f:
        line = f.readline().strip()

    if not line:
        raise ValueError("Empty SMILES file")

    parts = line.split()
    smiles = parts[0]
    name = parts[1] if len(parts) > 1 else base

    # -----------------------------
    # MOLECULE CREATION
    # -----------------------------
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")

    mol.SetProp("_Name", name)

    # IMPORTANT: hydrogens added AFTER protonation step
    mol = Chem.AddHs(mol)

    print("Atoms:", mol.GetNumAtoms())
    print("Rotatable bonds:", Descriptors.NumRotatableBonds(mol))

    # ============================================================
    # CONFORMER GENERATION
    # ============================================================
    params = AllChem.ETKDGv3()
    params.randomSeed = args.seed
    params.pruneRmsThresh = 0.5
    params.useSmallRingTorsions = True
    params.numThreads = 1

    conf_ids = AllChem.EmbedMultipleConfs(
        mol,
        numConfs=args.nconf,
        params=params
    )

    print("Conformers generated:", len(conf_ids))

    if len(conf_ids) == 0:
        raise ValueError("No conformers generated")

    energies = []

    # ============================================================
    # ENERGY OPTIMIZATION
    # ============================================================
    if AllChem.MMFFHasAllMoleculeParams(mol):
        print("Force field: MMFF")

        results = AllChem.MMFFOptimizeMoleculeConfs(
            mol,
            maxIters=args.maxiters
        )

        for cid, (not_conv, energy) in enumerate(results):
            energies.append((cid, energy, not_conv))

    else:
        print("Force field: UFF")

        for cid in conf_ids:
            not_conv = AllChem.UFFOptimizeMolecule(
                mol,
                confId=cid,
                maxIters=args.maxiters
            )
            ff = AllChem.UFFGetMoleculeForceField(mol, confId=cid)
            energy = ff.CalcEnergy()
            energies.append((cid, energy, not_conv))

    # ============================================================
    # SORT CONFORMERS
    # ============================================================
    energies.sort(key=lambda x: x[1])

    print("Conformer energies:")
    for cid, energy, nc in energies:
        status = "OK" if nc == 0 else "NOT_CONVERGED"
        print(f"conf {cid} : {energy:.4f} ({status})")

    # ============================================================
    # SELECT BEST CONFORMERS
    # ============================================================
    converged = [e for e in energies if e[2] == 0]

    if converged:
        print("Using only converged conformers")
        selected = converged
    else:
        print("WARNING: No converged conformers, using all")
        selected = energies

    # ============================================================
    # WRITE OUTPUT
    # ============================================================
    writer = Chem.SDWriter(outfile)

    written = 0

    for cid, energy, nc in selected[:args.keep]:

        mol.SetProp("Energy", str(energy))
        mol.SetProp("Converged", "yes" if nc == 0 else "no")
        mol.SetProp(
            "Method",
            "MMFF" if AllChem.MMFFHasAllMoleculeParams(mol) else "UFF"
        )

        writer.write(mol, confId=cid)
        written += 1

    writer.close()

    if written == 0:
        raise ValueError("No conformers written to output")

    print("SUCCESS")
    sys.exit(0)

# ============================================================
# ERROR HANDLING
# ============================================================
except Exception as e:
    print("FAILED:", str(e))
    sys.exit(1)
