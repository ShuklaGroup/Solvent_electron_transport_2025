# Solvent Electron Transport in Peptides (2025)

![Representative Figure](figures/fig1.png)

This repository supports the analysis of how different solvents impact electron transport in short peptides (4-mer and 5-mer). It contains scripts for simulation and analysis, as well as precomputed features extracted from molecular dynamics (MD) simulations. The goal is to understand solvent effects on hydrogen bonding, secondary structure, and conduction-relevant features.

---

## 🔍 What’s in this repository?

### ✅ `analysis_scripts/`

This folder includes high-level Python scripts to analyze features extracted from MD simulations.

- `4mer_main_analysis.py`: Analyzes 4-mer peptides across different solvents. It processes the features and hydrogen bonding data to evaluate conduction-related properties.
- `5mer_main_analysis.py`: Similar to the 4-mer script but focused on 5-mer peptide systems, including more complex hydrogen bonding patterns.

---

### 📦 `feature_pkls/`

This folder contains preprocessed features saved as `.pkl` files (Python pickles), ready for machine learning or further analysis. The features are divided by peptide length (4-mer vs 5-mer) and solvent conditions.

#### For 4-mers:
- Files like `feats_acn.pkl`, `feats_water.pkl`, etc., contain general features for different solvents.
- Files starting with `feats_HO_*.pkl` contain **hydrogen bond occupancy** (HO) features, e.g., `feats_HO_tfe.pkl` corresponds to 4-mer hydrogen bonding in TFE.

#### For 5-mers:
- Files like `feat_H_14_acn.pkl`, `feat_H_25_tfe.pkl`, etc., represent **specific hydrogen bonds**, such as:
  - `H_14`: hydrogen bond between residues 1 and 4
  - `H_15`: between residues 1 and 5
  - `H_25`: between residues 2 and 5
- These are provided across solvents (acetonitrile, glycerol, TFE).
- General feature files (e.g., `feats_acn.pkl`) are also provided for the 5-mers.

---

### 🧪 `simulation_scripts/`

This folder contains scripts for setting up and running MD simulations of peptides in various solvents.

- `minimization_NVT.py`: Minimizes the system energy under NVT ensemble.
- `equilibration.py`: Equilibrates the system before production.
- `NPT.py`: Runs simulations under constant pressure and temperature.
- `production.py`: Carries out the final production run for analysis.

---

### 🖼️ `figures/`

Contains figures used for visualization or in publications.  
- `fig1.png`: A representative figure showing key insights or results (e.g., hydrogen bonding trends across solvents).

---

### 📜 Other Files

- `README.md`: You’re reading it.
- `LICENSE`: The license governing the usage and distribution of this repository.

---

## 🧠 Context & Purpose

This project aims to correlate **solvent environment** with **electron transport properties** in peptides, focusing on how solvent-dependent hydrogen bonding and conformational states affect conductance. It is part of an effort to design solvent-optimized peptide-based materials for nanoscale electronic applications.

---

## 🚀 Getting Started

1. Run the MD simulations using the scripts in `simulation_scripts/` (or use existing data).
2. Load the feature `.pkl` files from `feature_pkls/` for downstream machine learning or statistical analysis.
3. Use `analysis_scripts/` to reproduce key figures or extract insights.
4. Visualize outputs using `figures/` or generate new ones.

---

## 🧬 Citation

If you use this repository or its contents in your work, please cite the associated publication (to be added once published).

---

Feel free to open issues or contact the maintainer for questions or suggestions!
