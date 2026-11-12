# NIST-OrgThermo-ML Database

[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![Python](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/)
[![Data Size](https://img.shields.io/badge/molecules-11,860-green.svg)](data/ml_dataset.csv)
[![DOI](https://img.shields.io/badge/DOI-pending-orange.svg)](https://zenodo.org/)

A curated thermochemistry database of **standard enthalpies of formation (ΔfH°)** for **11,860 organic molecules**, combining experimental data from the **NIST Chemistry WebBook** and computational data from the **QM9 dataset** (B3LYP/6-31G(2df,p)).

---

## 📊 Database Statistics

| Metric | Value |
|--------|-------|
| **Unique molecular structures** | 11,860 |
| **Molecules with gas-phase ΔfH°** | 11,803 |
| **Experimental gas-phase values (NIST)** | 166 |
| **Computed gas-phase values (QM9)** | 11,637 |
| **Total experimental observations** | 302 (gas/liquid/solid) |
| **Total computational observations** | 12,000 |
| **Experimental gold standard (Review method)** | 5 |
| **NIST ↔ QM9 matched by InChIKey** | 88 |
| **QM9 vs Exp MAE** | 80.12 kJ/mol |
| **QM9 vs Exp R²** | 0.9624 |

---

## 🚀 Quickstart

```python
import pandas as pd

# ML-ready dataset (one row per molecule with features + target)
ml_df = pd.read_csv('data/ml_dataset.csv')
print(f"Molecules: {len(ml_df)}")
print(f"Experimental targets: {(ml_df.target_label == 'experimental').sum()}")
print(f"Computed targets: {(ml_df.target_label == 'computed').sum()}")

# Structure registry
struct_df = pd.read_csv('data/structure_registry.csv')

# Experimental observations (NIST)
exp_obs = pd.read_csv('data/experimental_observations.csv')

# Computational observations (QM9)
comp_obs = pd.read_csv('data/computational_observations.csv')

# Quality control report
import json
with open('data/qc_report.json') as f:
    qc = json.load(f)
```

### Load with helper script
```python
from scripts.load_database import load_nist_qm9_database

db = load_nist_qm9_database('data/')
# Returns dict with: ml_dataset, structure_registry, experimental_observations, 
# computational_observations, experimental_gold, computational_qm9, 
# combined_master, qc_report, source_statistics
```

---

## 📁 Data Files

| File | Rows | Description |
|------|------|-------------|
| `data/structure_registry.csv` | 11,860 | Unique molecules with identifiers (InChIKey, SMILES, CAS, formula, MW) |
| `data/experimental_observations.csv` | 302 | NIST ΔfH° observations (gas/liquid/solid) with uncertainties, methods, refs |
| `data/computational_observations.csv` | 12,000 | QM9 ΔfH°(gas) at B3LYP/6-31G(2df,p) with all QM9 properties |
| `data/experimental_gold.csv` | 5 | High-quality NIST "Review" method values with uncertainties |
| `data/computational_qm9.csv` | 12,000 | Full QM9 records (U₀, H, G, ZPVE, dipole, HOMO/LUMO, etc.) |
| `data/combined_master.csv` | 12,302 | All observations merged with structure info |
| `data/ml_dataset.csv` | 11,860 | **ML-ready**: features (atom counts, DoU, RDKit descriptors) + target |
| `data/qc_report.json` | — | Quality metrics, validation results, duplicate checks |
| `data/source_statistics.csv` | — | Per-source summary statistics |

See [DATA_DICTIONARY.md](docs/DATA_DICTIONARY.md) for complete column documentation.

---

## 🔬 Methodology

### Data Sources
- **NIST Chemistry WebBook** (SRD 69): Evaluated experimental thermochemistry, accessed Nov 2026
- **QM9 Dataset** (Ramakrishnan et al., *Sci. Data* 2014): 134k molecules at B3LYP/6-31G(2df,p)

### Thermochemical Cycle (QM9 → ΔfH°)
```
ΔfH°(mol, 298 K) = Σ nᵢ × ΔfH°_exp(atomᵢ) 
                   - [Σ nᵢ × U₀_lit(atomᵢ) - U₀_QM9(mol)] 
                   + [H_QM9(298) - U₀_QM9(mol)]
```
- Atomic U₀ from literature at **same DFT level** (B3LYP/6-31G(2df,p))
- Experimental atomic ΔfH°: H=218.0, C=716.7, N=472.7, O=249.2, F=79.4 kJ/mol
- 1 Hartree = 2625.49962 kJ/mol

### Deduplication
Structures unified by **InChIKey** (preferred) → InChI → SMILES → CAS → name+formula.  
**Zero duplicate InChIKeys** in final registry.

See [METHODOLOGY.md](docs/METHODOLOGY.md) for full details.

---

## ✅ Quality Control

| Check | Result |
|-------|--------|
| Duplicate InChIKeys | 0 |
| Duplicate CAS | 11,622 (QM9 lacks CAS) |
| Invalid SMILES (RDKit) | 0 |
| Missing InChIKey | 21 |
| Missing SMILES | 21 |
| QM9 vs Exp (88 matched) MAE | 80.12 kJ/mol |
| QM9 vs Exp R² | 0.9624 |

### Known Limitations
1. **QM9 accuracy**: B3LYP/6-31G(2df,p) has systematic error (~80 kJ/mol MAE vs experiment)
2. **Element coverage**: QM9 limited to C, H, O, N, F; NIST has broader elements
3. **Phase data**: Liquid/solid ΔfH° sparse in both sources
4. **Temperature**: All values at 298.15 K; no temperature dependence
5. **Subset**: QM9 processed 12,000/133,885 molecules (9%) for tractability

---

## 📦 Installation

```bash
pip install -r requirements.txt
# Requires: pandas, numpy, rdkit-pypi
```

---

## 🔬 Reproducibility

```bash
# 1. Scrape NIST (resumable, incremental)
python run_large_scrape.py

# 2. Process QM9
python process_qm9.py --max-molecules 12000

# 3. Build final database
python build_final_database.py
```

Full environment and pipeline details in [REPRODUCIBILITY.md](docs/REPRODUCIBILITY.md).

---

## 📖 Citation

If you use this database, please cite:

```bibtex
@dataset{nist_orgthermo_ml_2026,
  title = {NIST-OrgThermo-ML: Curated Thermochemistry Database for Organic Molecules from NIST and QM9},
  author = {Anwesh Pandey},
  year = {2026},
  publisher = {GitHub},
  url = {https://github.com/anweshpandey/thermo-nist-qm9-database},
  doi = {pending},
  note = {11,860 molecules; 166 experimental + 11,637 computed ΔfH°(gas) at 298.15 K}
}
```

Also see [CITATION.cff](CITATION.cff) for machine-readable citation metadata.

### Primary References
1. **NIST Chemistry WebBook**, NIST Standard Reference Database 69, https://webbook.nist.gov
2. Ramakrishnan, R. et al. "Quantum chemistry structures and properties of 134 kilo molecules." *Scientific Data* **1**, 140022 (2014). DOI: [10.1038/sdata.2014.22](https://doi.org/10.1038/sdata.2014.22)
3. Feller, D. et al. "Thermochemistry of the QM9 dataset." *J. Chem. Phys.* **140**, 18A505 (2014).

---

## 📄 License

This database is released under **Creative Commons Attribution 4.0 International (CC BY 4.0)**.

- **You may**: Share, copy, redistribute, adapt, and build upon the data for any purpose, even commercially.
- **You must**: Give appropriate credit, provide a link to the license, and indicate if changes were made.

See [LICENSE](LICENSE) for full legal text.

---

## 🤝 Contributing

Issues and PRs welcome for:
- Data corrections/updates
- Additional validation checks
- Improved documentation
- Extended processing scripts

---

## 📬 Contact

- **Issues**: [GitHub Issues](https://github.com/anweshpandey/thermo-nist-qm9-database/issues)

---

*Generated: 2026-11-12 | NIST-OrgThermo-ML v1.0*
