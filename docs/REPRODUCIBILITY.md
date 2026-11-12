# Reproducibility Report

## Environment
- Python: 3.10+
- Key packages: pandas, numpy, rdkit, tqdm
- RDKit version: 2026.09.1

## Data Acquisition

### NIST Chemistry WebBook
- Source: https://webbook.nist.gov/cgi/cbook.cgi
- Scraper: nist_thermo package (custom)
- Date accessed: 2026-11-11
- Rate limit: 0.5-1.0 seconds between requests
- Formulas searched: Comprehensive organic classes (alkanes, alkenes, aromatics, alcohols, ethers, carbonyls, acids, esters, amines, nitriles, amides, heterocycles, halogenated)

### QM9 Dataset
- Source: https://figshare.com/articles/dataset/QM9/9668873
- Original publication: Ramakrishnan et al., Scientific Data 1, 140022 (2014)
- Level of theory: B3LYP/6-31G(2df,p)
- Properties: U₀, H(298), G(298), ZPVE, dipole, polarizability, HOMO/LUMO, etc.
- Molecules: ~134k organic molecules (C, H, O, N, F)

## Processing Pipeline

1. **NIST Scraping** (`run_large_scrape.py`, `run_10k_scrape.py`)
   - Search by molecular formula
   - Extract all isomers per formula
   - Parse gas/liquid/solid ΔfH°, entropy, Cp, phase change data
   - Deduplicate by CAS number

2. **QM9 Processing** (`process_qm9.py`)
   - Load SDF (3D structures) and CSV (properties)
   - Match by mol_id (gdb_X)
   - Compute ΔfH°(298 K) via thermochemical cycle:
     - Literature atomic U₀ at B3LYP/6-31G(2df,p)
     - Experimental atomic ΔfH° at 298 K
     - QM9 molecular U₀ and H(298)
   - Generate SMILES, InChI, InChIKey via RDKit
   - Compute molecular descriptors

3. **Database Integration** (`build_final_database.py`)
   - Unify structures by InChIKey
   - Separate experimental vs computational observations
   - Create structure registry with molecule_id
   - Build ML dataset with features

## Thermochemical Cycle Details

### Atomic Reference Values (B3LYP/6-31G(2df,p))
| Element | U₀ (Hartree) | ΔfH°_exp (kJ/mol) |
|---------|--------------|-------------------|
| H       | -0.500273    | 218.0             |
| C       | -37.846772   | 716.7             |
| N       | -54.583861   | 472.7             |
| O       | -75.064579   | 249.2             |
| F       | -99.718730   | 79.4              |

### Conversion
1 Hartree = 2625.49962 kJ/mol

### Formula
ΔfH°(mol, 298) = Σ nᵢ × ΔfH°_exp(atomᵢ) - [Σ nᵢ × U₀_lit(atomᵢ) - U₀_QM9(mol)] + [H_QM9(298) - U₀_QM9(mol)]

## Quality Control

- Structure deduplication by InChIKey
- Formula validation via RDKit
- Charge neutrality check
- Experimental vs computational comparison (88 matched)
- Duplicate CAS/InChIKey detection
- Missing identifier reporting

## Known Limitations

1. **QM9 accuracy**: B3LYP/6-31G(2df,p) has systematic errors (MAE ~80 kJ/mol vs experiment)
2. **NIST coverage**: Limited to molecules in WebBook (~hundreds with ΔfH°)
3. **Phase data**: Liquid/solid ΔfH° sparse in both sources
4. **Elements**: QM9 limited to C, H, O, N, F; NIST has broader elemental coverage
5. **Temperature**: All values at 298.15 K; no temperature dependence

## Running the Pipeline

```bash
# 1. Scrape NIST (resumable)
python run_large_scrape.py

# 2. Process QM9
python process_qm9.py --max-molecules 12000

# 3. Build final database
python build_final_database.py
```

## Output Statistics (Current Run)

- Unique molecular structures: 11860
- Experimental observations: 302
- Computational observations: 12000
- Experimental gold standard: 5
- ML dataset size: 11860 molecules
- ML targets available: 11803

Generated: 2026-11-11T15:14:45.020605
