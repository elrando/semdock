<p align="center">
  <picture>
   <source media="(prefers-color-scheme: dark)" srcset="https://github.com/elrando/semdock/blob/main/asset/banner.jpg">
    <source media="(prefers-color-scheme: light)" srcset="https://github.com/elrando/semdock/blob/main/asset/banner.jpg">
    <img src="https://github.com/elrando/semdock/blob/main/asset/banner.jpg" alt="semdock banner" width="850">
  </picture>
</p>

# semdock [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**semdock** is a CPU-parallelized molecular docking workflow for large-scale virtual screening using **AutoDock Vina**.  
The pipeline automates ligand preparation from multiple chemical formats, receptor preparation, docking execution with job resumption, and post-docking analysis with optional pose similarity scoring.

The workflow consists of four command-line tools:

- **sem_prep** – Converts raw chemical datasets (`.sdf`, `.smi`, `.csv`) into docking-ready PDBQT ligands.
- **jamreceptor** – Prepares receptor structures and generates docking grid configuration files.
- **sem_dock** – Runs CPU-parallel AutoDock Vina docking with resume support.
- **sem_filter** – Analyzes docking results, ranks ligands, filters hits, and optionally computes pose similarity scores (SimScore).

---



# Workflow overview

```text
Input molecules
   ↓
sem_prep
   ↓
Docking-ready ligands (.pdbqt)
   ↓
jamreceptor
   ↓
Prepared receptor + conf.txt
   ↓
sem_dock
   ↓
Docking outputs
   ↓
sem_filter
   ↓
Ranked + filtered results
```

---

# Features

## sem_prep
- Accepts:
  - `.sdf`
  - `.smi`
  - `.csv`
- Ligand standardization
- Salt removal
- Protonation using **Dimorphite-DL**
- 3D conformer generation using **RDKit**
- Geometry optimization
- PDBQT conversion using **Meeko**
- Source tracking:
  - PubChem
  - ZINC
  - ChEMBL
  - Custom IDs

Outputs:
- standardized ligands
- protonated SMILES
- 3D SDF
- PDBQT files
- `id_map.csv`

---

## jamreceptor
Prepares receptor structures for docking.

Features:
- Converts receptor PDB → PDBQT
- Optional binding pocket detection using **fpocket**
- Generates docking box coordinates
- Writes `conf.txt`

Outputs:
- `receptor.pdbqt`
- `conf.txt`

---

## sem_dock
Runs AutoDock Vina docking in parallel.

Features:
- CPU parallelization
- Resume interrupted jobs
- Reads docking parameters from `conf.txt`
- Supports user-defined:
  - exhaustiveness
  - num_modes
  - energy_range

Outputs:
- docked PDBQT files
- docking logs

---

## sem_filter
Post-docking analysis and filtering.

Features:
- Generates ranked results for **all ligands**
- Threshold-based filtering
- Top-N filtering
- Optional descriptor calculation
- Optional pose similarity scoring (`--simscore`)

SimScore:
- Calculated from Vina pose RMSD clustering
- Reported only when `--simscore` flag is used
- **Not used for ranking**
- If only one docking mode exists:
  - `SimScore = 0`

Outputs:
- `Results.csv`
- `Filtered_Results.csv`
- `Filtered_Features.csv`

---

# System setup & Installation

**semdock** is designed for Linux-based systems and is intended to run inside a **Conda environment** for simplified dependency management and reproducibility.

These instructions assume a clean Ubuntu/Linux system.

---

## 1. Update system and install essential packages

Open a terminal and run:

```bash
sudo apt update && sudo apt upgrade -y
sudo apt install -y build-essential nano wget curl git cmake libboost-all-dev
```

> **Note:** You may be prompted to enter your superuser password.

---

## 2. Install Miniconda and create environment

Conda provides isolated software environments for reproducible installation.

### Download and install Miniconda

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc
```

### Create and activate semdock environment

```bash
conda create --name semdock python=3.10 pip
conda activate semdock
conda config --env --add channels conda-forge
```

---

## 3. Install software dependencies

Install required packages:

```bash
conda install -c conda-forge rdkit meeko parallel pymol-open-source
pip install dimorphite_dl
```

---

## 4. Install AutoDockTools (MGLTools)

AutoDockTools is required by `jamreceptor`.

### Download and extract

```bash
mkdir ~/programs
cd ~/programs
wget https://ccsb.scripps.edu/mgltools/download/491/mgltools_Linux-x86_64_1.5.7.tar.gz
tar -zxf mgltools_Linux-x86_64_1.5.7.tar.gz
```

### Install

```bash
cd mgltools_x86_64Linux2_1.5.7
./install.sh
```

### Add to PATH

```bash
echo 'export PATH=$HOME/programs/mgltools_x86_64Linux2_1.5.7/bin:$PATH' >> ~/.bashrc
source ~/.bashrc
```

---

## 5. Install fpocket

fpocket is used by `jamreceptor` for binding pocket detection.

```bash
cd ~/programs
git clone https://github.com/Discngine/fpocket.git
cd fpocket
make
sudo make install
```

---

## 6. Install AutoDock Vina

Clone and build AutoDock Vina:

```bash
cd ~/programs
git clone https://github.com/ccsb-scripps/AutoDock-Vina
cd AutoDock-Vina
mkdir build
cd build
cmake ..
make -j$(nproc)
```

Add Vina to PATH:

```bash
echo 'export PATH=$HOME/programs/AutoDock-Vina/build:$PATH' >> ~/.bashrc
source ~/.bashrc
```

---

## 7. Install semdock

Clone repository:

```bash
cd ~/programs
git clone https://github.com/elrando/semdock.git
cd semdock
```

Make scripts executable:

```bash
chmod +x sem_* jamreceptor fail_check
chmod +x py/*.py
chmod +x lib/*.sh
```

Add semdock to PATH:

```bash
echo 'export PATH=$HOME/programs/semdock:$PATH' >> ~/.bashrc
source ~/.bashrc
```

Reactivate environment:

```bash
conda activate semdock
```

This enables execution of:

- `sem_prep`
- `jamreceptor`
- `sem_dock`
- `sem_filter`
- `fail_check`

from any terminal window.

---

# Usage

The semdock workflow consists of four sequential steps:

1. Ligand preparation (`sem_prep`)
2. Receptor preparation (`jamreceptor`)
3. Molecular docking (`sem_dock`)
4. Result filtering (`sem_filter`)

---

## 1. Prepare ligands

Convert raw compound libraries into docking-ready PDBQT files.

Supported input formats:

* `.sdf`
* `.smi`
* `.csv`

Example:

```bash
conda activate semdock

mkdir my_project
cd my_project

sem_prep --input compounds.csv --ph 7.4 --cpu 16 --out_dir library
```

Alternative CPU allocation:

```bash
sem_prep --input compounds.csv --ph 7.4 --cpu_fraction 0.8 --out_dir library
```

Output:

```text
library/
├── smiles/
├── protonation/
├── 3D/
├── pdbqt/
├── logs/
└── id_map.csv
```

---

## 2. Prepare receptor

Prepare receptor structure and docking grid.

Create receptor directory:

```bash
mkdir receptor_A
cd receptor_A
```

Place receptor PDB file inside this directory, then run:

```bash
jamreceptor
```

The script will prompt for:

* receptor PDB filename
* chain selection
* binding pocket selection
* grid padding

Output:

```text
receptor.pdbqt
conf.txt
grid_box.py
```

Optional grid visualization:

```bash
pymol grid_box.py
```

---

## 3. Run docking

Dock prepared ligands using AutoDock Vina.

Run inside receptor directory:

```bash
sem_dock \
    --config conf.txt \
    --lig_dir ../library \
    --cpu 16 \
    --out_dir .
```

Alternative CPU allocation:

```bash
sem_dock \
    --config conf.txt \
    --lig_dir ../library \
    --cpu_fraction 0.8 \
    --out_dir .
```

Optional docking parameters:

```bash
sem_dock \
    --config conf.txt \
    --lig_dir ../library \
    --cpu 16 \
    --exhaustiveness 16 \
    --max_poses 20 \
    --energy_range 5 \
    --out_dir .
```

Resume interrupted jobs:

Re-run the same command. Previously completed ligands are skipped automatically.

Docking output:

```text
vina_out/
vina_logs/
ligands.txt
failed_docking.txt
```

---

## 4. Filter results

Analyze docking results and rank compounds.

Basic threshold filtering:

```bash
sem_filter \
    --dock_dir . \
    --prep_dir ../library \
    --threshold -7.0
```

Top-N filtering:

```bash
sem_filter \
    --dock_dir . \
    --prep_dir ../library \
    --top 100
```

Enable SimScore calculation:

```bash
sem_filter \
    --dock_dir . \
    --prep_dir ../library \
    --top 100 \
    --simscore
```

Output:

```text
analysis/
├── Results.csv
├── Filtered_Results.csv
└── Filtered_Features.csv
```

Notes:

* `SimScore` is calculated only when `--simscore` is provided.
* `SimScore` is reported only and is **not used for ranking**.
* Ligands with a single docking mode receive `SimScore = 0`.

---

## 5. Check failed ligand preparation

Use `fail_check` to identify ligands that failed during preparation.

```bash
fail_check --prep_dir library
```

Detects failures including:

* missing 3D conformers
* missing SDF files
* missing PDBQT files
* incomplete ligand processing

```
```


# Limitations

- Rigid receptor docking only
- Single protonation state per ligand
- Single conformer per ligand
- Binding affinities are approximate screening scores
- Receptor quality strongly affects docking reliability
- SimScore reflects pose clustering only and does not predict binding validity

---

# Citation

If you use semdock in research, please cite:
- Emon, S. (2026). Scalable Virtual Screening Pipeline: CPU‑Parallelized Docking and Robust Ligand Preparation from Raw Chemical Data (Version 1). Zenodo. https://doi.org/10.5281/zenodo.19998830
- Barbosa Pereira, P. J., Ripoll-Rozada, J., Macedo-Ribeiro, S., & Manso, J. A. (2025). Protocol for an automated virtual screening pipeline including library generation and docking evaluation. STAR Protocols, 6(4), 104161. https://doi.org/10.1016/j.xpro.2025.104161
- Manso, J.A. (2025). A collection of bash scripts for a virtual screening pipeline: from compound library generation to docking score evaluation. Zenodo, https://doi.org/10.5281/zenodo.15577778.
- Eberhardt, J., Santos-Martins, D., Tillack, A. F., & Forli, S. (2021). AutoDock Vina 1.2.0: New Docking Methods, Expanded Force Field, and Python Bindings. Journal of Chemical Information and Modeling, 61(8), 3891–3898. https://doi.org/10.1021/acs.jcim.1c00203
- Morris, G.M., Huey, R., Lindstrom, W., Sanner, M.F., Belew, R.K., Goodsell, D.S., and Olson, A.J. (2009). AutoDock4 and AutoDockTools4: Automated Docking with Selective Receptor Flexibility. J Comput Chem 30, 2785–2791. https://doi.org/10.1002/jcc.21256.
- Trott, O., and Olson, A.J. (2010). AutoDock Vina: Improving the speed and accuracy of docking with a new scoring function, efficient optimization, and multithreading. Journal of Computational Chemistry 31, 455–461. https://doi.org/10.1002/jcc.21334.
- Le Guilloux, V., Schmidtke, P., and Tuffery, P. (2009). Fpocket: An open source platform for ligand pocket detection. BMC Bioinformatics 10, 168. https://doi.org/10.1186/1471-2105-10-168.
- Irwin, J.J., and Shoichet, B.K. (2005). ZINC − A Free Database of Commercially Available Compounds for Virtual Screening. J. Chem. Inf. Model. 45, 177–182. https://doi.org/10.1021/ci049714+.
- Sterling, T., and Irwin, J.J. (2015). ZINC 15 – Ligand Discovery for Everyone. J. Chem. Inf. Model. 55, 2324-2337. https://doi.org/10.1021/acs.jcim.5b00559.
- Irwin et al. (2020). ZINC20—A Free Ultralarge-Scale Chemical Database for Ligand Discovery. J. Chem. Inf. Model. 60, 6065–6073. https://doi.org/10.1021/acs.jcim.0c00675.
---

# License

This project is licensed under the MIT License - see the LICENSE file for details.

---

# Contact

**Shahariar Emon**  
Email: smemon275@gmail.com  
GitHub: https://github.com/elrando
