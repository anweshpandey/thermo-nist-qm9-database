# NIST-OrgThermo-ML Database - Data Dictionary

## Overview
This database contains thermochemical data (standard enthalpies of formation) for organic molecules from experimental (NIST Chemistry WebBook) and computational (QM9 DFT) sources.

## Files

### 1. structure_registry.csv
Unique molecular structures with identifiers.

| Column | Type | Description |
|--------|------|-------------|
| molecule_id | int | Unique integer identifier |
| source | str | Source database (NIST_WebBook, QM9) |
| source_id | str | Source-specific ID (CAS, gdb_X) |
| name | str | Chemical name |
| formula | str | Molecular formula (Hill order) |
| cas_number | str | CAS Registry Number |
| inchi | str | IUPAC Standard InChI |
| inchikey | str | IUPAC Standard InChIKey (27 chars) |
| smiles | str | Canonical SMILES |
| molecular_weight | float | Molecular weight (g/mol) |
| molecular_charge | int | Formal charge |
| isotope_info | str | Isotope labeling info (empty if natural) |
| structure_key | str | Deduplication key (InChIKey preferred) |

### 2. experimental_observations.csv
Experimental thermochemical observations from NIST.

| Column | Type | Description |
|--------|------|-------------|
| observation_id | str | Unique observation ID |
| molecule_id | int | Foreign key to structure_registry |
| source_database | str | NIST_Chemistry_WebBook |
| target_property | str | delta_f_H_gas, delta_f_H_liquid, delta_f_H_solid |
| target_value | float | Enthalpy value (kJ/mol) |
| target_unit | str | kJ/mol |
| phase | str | gas, liquid, solid |
| temperature_K | float | Reference temperature (typically 298.15) |
| pressure | str | Reference pressure (1 bar) |
| uncertainty | float | Standard uncertainty (kJ/mol) |
| method | str | NIST method code (Review, Ccb, AVG, etc.) |
| reference | str | Literature citation |
| source_record_url | str | NIST WebBook URL |
| quality_flag | str | experimental, experimental_evaluated |
| label_origin | str | experimental |

### 3. computational_observations.csv
Computational thermochemical observations from QM9.

| Column | Type | Description |
|--------|------|-------------|
| observation_id | str | Unique observation ID |
| molecule_id | int | Foreign key to structure_registry |
| source_database | str | QM9 |
| target_property | str | delta_f_H_gas |
| target_value | float | Computed ΔfH° (kJ/mol) |
| target_unit | str | kJ/mol |
| phase | str | gas |
| temperature_K | float | 298.15 |
| pressure | str | 1 bar |
| uncertainty | float | None (no uncertainty for DFT) |
| method | str | B3LYP/6-31G(2df,p)//B3LYP/6-31G(2df,p) |
| reference | str | QM9 dataset (Ramakrishnan et al., 2014) |
| source_record_url | str | Figshare URL |
| quality_flag | str | computational_DFT_B3LYP_6-31G2dfp |
| label_origin | str | computed |

### 4. experimental_gold.csv
High-quality experimental subset (Review method, with uncertainty).

### 5. computational_qm9.csv
Full QM9 computational records with all properties.

### 6. combined_master.csv
All observations merged with structure info.

### 7. ml_dataset.csv
ML-ready dataset: one row per molecule with features and target.

| Column | Type | Description |
|--------|------|-------------|
| molecule_id | int | Structure registry ID |
| name | str | Chemical name |
| formula | str | Molecular formula |
| inchikey | str | InChIKey |
| smiles | str | SMILES |
| molecular_weight | float | MW (g/mol) |
| n_C, n_H, n_O, n_N | int | Atom counts |
| n_F, n_Cl, n_Br, n_I | int | Halogen counts |
| n_halogen | int | Total halogens |
| DoU | float | Degree of unsaturation |
| delta_f_H_gas_kJ_mol | float | Target value |
| target_source | str | NIST_WebBook or QM9 |
| target_label | str | experimental or computed |
| rdkit_* | float | RDKit descriptors |

### 8. qc_report.json
Quality control metrics.

### 9. source_statistics.csv
Per-source statistics.

## Thermochemical Cycle (QM9)

ΔfH°(298 K) = Σ nᵢ × ΔfH°_exp(atomᵢ, 298 K) - E_atom(0 K) + H_corr(0→298 K)

Where:
- E_atom(0 K) = Σ nᵢ × U₀_lit(atomᵢ) - U₀_QM9(mol)
- H_corr(0→298 K) = H_QM9(298) - U₀_QM9(mol)
- U₀_lit from B3LYP/6-31G(2df,p) literature values

## Quality Flags

- experimental: Standard experimental value
- experimental_evaluated: NIST Review method (evaluated/averaged)
- computational_DFT_B3LYP_6-31G2dfp: B3LYP/6-31G(2df,p) DFT

## Label Origins

- experimental: Measured/evaluated experimental value
- computed: Derived from quantum chemical calculations

## Validation

QM9 vs Experimental (88 matched by InChIKey):
- MAE: 80.12 kJ/mol
- RMSE: 87.83 kJ/mol
- R²: 0.9624

## References

1. NIST Chemistry WebBook, NIST Standard Reference Database Number 69
2. Ramakrishnan et al., "Quantum chemistry structures and properties of 134 kilo molecules", Scientific Data 1, 140022 (2014)
