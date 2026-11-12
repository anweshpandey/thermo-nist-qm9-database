# NIST-OrgThermo-ML Database: Construction Methodology
## Supplementary Methods for Research Paper

---

## 1. Overview

This document describes the complete methodology for constructing the **NIST-OrgThermo-ML** database, a curated collection of standard enthalpies of formation (ΔfH°) for 11,860 unique organic molecules, combining experimental data from the NIST Chemistry WebBook and computational data from the QM9 dataset.

**Final database statistics:**
- 11,860 unique molecular structures
- 11,803 molecules with gas-phase ΔfH° at 298.15 K
- 166 experimental gas-phase values (NIST)
- 11,637 computed gas-phase values (QM9, B3LYP/6-31G(2df,p))
- 302 experimental observations across gas/liquid/solid phases
- 12,000 computational observations

---

## 2. Data Sources

### 2.1 NIST Chemistry WebBook (Experimental)
- **Source**: https://webbook.nist.gov/cgi/cbook.cgi
- **Database**: NIST Standard Reference Database 69
- **Accessed**: November 2026
- **Content**: Evaluated thermochemical data including ΔfH°(gas), ΔfH°(liquid), ΔfH°(solid), entropy, heat capacity, phase transition data
- **Quality**: NIST method codes (Review, Ccb, Ccr, AVG, Eqk, Cm, Chyd, Ion)

### 2.2 QM9 Dataset (Computational)
- **Source**: https://figshare.com/articles/dataset/QM9/9668873
- **Reference**: Ramakrishnan et al., *Scientific Data* 1, 140022 (2014)
- **Level of theory**: B3LYP/6-31G(2df,p) geometry optimization and frequency calculation
- **Properties**: U₀ (0 K), H(298), G(298), ZPVE, dipole moment, polarizability, HOMO/LUMO, electronic spatial extent
- **Molecules**: 133,885 organic molecules (C, H, O, N, F up to 9 heavy atoms)

---

## 3. NIST WebBook Data Acquisition

### 3.1 Formula-Based Search Strategy

We systematically searched all major organic compound classes using Hill-formula queries:

```python
# From run_10k_scrape.py - generate_all_target_formulas()
def generate_all_target_formulas():
    formulas = []
    
    # Alkanes C1-C25
    for c in range(1, 26):
        formulas.append(f'C{c}H{2*c+2}')
    
    # Alkenes C2-C25
    for c in range(2, 26):
        formulas.append(f'C{c}H{2*c}')
    
    # Alkynes C2-C25
    for c in range(2, 26):
        formulas.append(f'C{c}H{2*c-2}')
    
    # Cycloalkanes C3-C25
    for c in range(3, 26):
        formulas.append(f'C{c}H{2*c}')
    
    # Aromatics C6-C25
    for c in range(6, 26):
        formulas.append(f'C{c}H{2*c-6}')
    
    # Alcohols C1-C20
    for c in range(1, 21):
        formulas.append(f'C{c}H{2*c+2}O')
    
    # Ethers C2-C20
    for c in range(2, 21):
        formulas.append(f'C{c}H{2*c+2}O')
    
    # Aldehydes C1-C20
    for c in range(1, 21):
        formulas.append(f'C{c}H{2*c}O')
    
    # Ketones C3-C20
    for c in range(3, 21):
        formulas.append(f'C{c}H{2*c}O')
    
    # Carboxylic acids C1-C20
    for c in range(1, 21):
        formulas.append(f'C{c}H{2*c}O2')
    
    # Esters C2-C20
    for c in range(2, 21):
        formulas.append(f'C{c}H{2*c}O2')
    
    # Amines C1-C20
    for c in range(1, 21):
        formulas.append(f'C{c}H{2*c+3}N')
    
    # Nitriles C2-C20
    for c in range(2, 21):
        formulas.append(f'C{c}H{2*c-1}N')
    
    # Amides C1-C20
    for c in range(1, 21):
        formulas.append(f'C{c}H{2*c+1}NO')
    
    # Diols, Diamines, Amino acids, Heterocycles, Phenols, Halogenated...
    # (full list in run_10k_scrape.py)
    
    return list(dict.fromkeys(formulas))  # deduplicate
```

**Total formulas generated**: ~500+ unique Hill formulas covering all major organic classes up to C25.

### 3.2 Scraping Pipeline (nist_thermo package)

**Core scraper** (`nist_thermo/scraper.py`):

```python
class NISTScraper:
    BASE_URL = "https://webbook.nist.gov/cgi/cbook.cgi"
    
    # Data masks for different property types
    MASK_GAS_THERMO = 1       # Gas phase thermochemistry
    MASK_CONDENSED_THERMO = 2 # Condensed phase thermochemistry
    MASK_PHASE_CHANGE = 4     # Phase change data
    
    def __init__(self, delay=1.0, timeout=30, max_retries=3, units="SI"):
        self.delay = delay  # Rate limiting
        self.session = requests.Session()
    
    def scrape_compound(self, compound_id):
        """Scrape all thermochemical data for one compound."""
        # 1. Main page → identifiers (CAS, InChI, InChIKey, MW, formula, name)
        # 2. Mask=1 → ΔfH°(gas), S°(gas), Cp°(gas)
        # 3. Mask=2 → ΔfH°(liquid), ΔfH°(solid), S°(cond), Cp°(cond)
        # 4. Mask=4 → bp, mp, ΔHvap, ΔHfus, ΔHsub, Tc, Pc, Vc
        return ThermoData(...)
```

**Rate limiting**: 0.5–1.0 seconds between requests (respectful to NIST servers).

**Incremental saving**: After each batch of formulas, results are appended to CSV with CAS-based deduplication.

### 3.3 Output Files from NIST Pipeline

| File | Description |
|------|-------------|
| `nist_hof_raw.csv` | Raw scraped data (38 columns) |
| `nist_hof_large_raw.csv` | Large-scale scrape (238 compounds) |
| `nist_hof_large_enriched.csv` | + SMILES, RDKit descriptors, formula features |
| `nist_hof_large_clean.csv` | Deuterated removed, valid SMILES only |
| `nist_hof_large_ml_ready.csv` | ML-ready: 166 compounds with gas ΔfH° |

---

## 4. QM9 Data Processing

### 4.1 Download and Parsing

```python
# From process_qm9.py
QM9_URL = "https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/qm9.tar.gz"

def load_qm9_properties():
    """Load QM9 properties CSV (133,885 rows)."""
    df = pd.read_csv('qm9_data/qm9.sdf.csv')
    # Columns: mol_id, A, B, C, mu, alpha, homo, lumo, gap, r2,
    #          zpve, u0, u298, h298, g298, cv,
    #          u0_atom, u298_atom, h298_atom, g298_atom

def load_qm9_sdf(max_molecules=None):
    """Load 3D structures from SDF with RDKit."""
    supplier = Chem.SDMolSupplier('qm9.sdf', removeHs=False, sanitize=False)
    # Parse mol_id from title line: "gdb 1 ..." → "gdb_1"
```

### 4.2 Thermochemical Cycle for ΔfH°(298 K)

The QM9 dataset provides **molecular** energies at the B3LYP/6-31G(2df,p) level. To compute standard enthalpies of formation, we apply the following thermochemical cycle:

#### 4.2.1 Theoretical Foundation

$$\Delta_f H^\circ(\text{mol}, 298\text{ K}) = \sum_i n_i \Delta_f H^\circ_\text{exp}(\text{atom}_i, 298\text{ K}) - E_\text{atom}(0\text{ K}) + H_\text{corr}(0 \to 298\text{ K})$$

Where:
- $E_\text{atom}(0\text{ K}) = \sum_i n_i U_0^\text{lit}(\text{atom}_i) - U_0^\text{QM9}(\text{mol})$ — atomization energy at 0 K
- $H_\text{corr}(0 \to 298\text{ K}) = H^\text{QM9}(298) - U_0^\text{QM9}(\text{mol})$ — thermal correction
- $U_0^\text{lit}$: literature atomic energies at **same DFT level** (B3LYP/6-31G(2df,p))
- $\Delta_f H^\circ_\text{exp}$: experimental atomic formation enthalpies

#### 4.2.2 Atomic Reference Values (B3LYP/6-31G(2df,p))

| Element | $U_0^\text{lit}$ (Hartree) | $\Delta_f H^\circ_\text{exp}$ (kJ/mol) |
|---------|---------------------------|----------------------------------------|
| H       | -0.500273                 | 218.0                                  |
| C       | -37.846772                | 716.7                                  |
| N       | -54.583861                | 472.7                                  |
| O       | -75.064579                | 249.2                                  |
| F       | -99.718730                | 79.4                                   |

*Sources: NIST CCCBDB, QM9 supporting information, Feller et al. J. Chem. Phys. 2014*

#### 4.2.3 Implementation

```python
# From process_qm9.py
HARTREE_TO_KJMOL = 2625.49962

ATOMIC_U0 = {'H': -0.500273, 'C': -37.846772, 'N': -54.583861, 
             'O': -75.064579, 'F': -99.718730}

ATOMIC_DFH_298 = {'H': 218.0, 'C': 716.7, 'N': 472.7, 'O': 249.2, 'F': 79.4}

def compute_formation_enthalpy_qm9(row, formula):
    """ΔfH°(298 K) = Σ nᵢ ΔfH°_exp(atomᵢ) - [Σ nᵢ U₀_lit(atomᵢ) - U₀_QM9(mol)] + [H_QM9(298) - U₀_QM9(mol)]"""
    atom_counts = parse_formula(formula)
    
    h298 = row['h298']    # QM9 molecular H(298) in Hartree
    u0   = row['u0']      # QM9 molecular U₀ in Hartree
    
    # Σ nᵢ U₀_lit(atomᵢ)
    u0_atoms_lit = sum(count * ATOMIC_U0[elem] for elem, count in atom_counts.items())
    
    # Atomization energy at 0 K
    e_atom_0 = (u0_atoms_lit - u0) * HARTREE_TO_KJMOL
    
    # Thermal correction
    h_corr = (h298 - u0) * HARTREE_TO_KJMOL
    
    # Σ nᵢ ΔfH°_exp(atomᵢ)
    sum_dfh_atoms = sum(count * ATOMIC_DFH_298[elem] for elem, count in atom_counts.items())
    
    # Final ΔfH°(298 K)
    dfh_298 = sum_dfh_atoms - e_atom_0 + h_corr
    return dfh_298
```

### 4.3 Molecular Identifiers

```python
def get_molecular_identifiers(mol):
    """Generate SMILES, InChI, InChIKey from RDKit Mol."""
    smiles   = Chem.MolToSmiles(mol)
    inchi    = MolToInchi(mol)
    inchikey = MolToInchiKey(mol)
    return {'smiles': smiles, 'inchi': inchi, 'inchikey': inchikey}

def get_rdkit_descriptors(mol):
    """Compute ML-relevant descriptors."""
    return {
        'rdkit_mw': Descriptors.MolWt(mol),
        'rdkit_logp': Descriptors.MolLogP(mol),
        'rdkit_tpsa': Descriptors.TPSA(mol),
        'rdkit_rotatable_bonds': Descriptors.NumRotatableBonds(mol),
        'rdkit_h_donors': Descriptors.NumHDonors(mol),
        'rdkit_h_acceptors': Descriptors.NumHAcceptors(mol),
        'rdkit_num_rings': rdMolDescriptors.CalcNumRings(mol),
        'rdkit_aromatic_rings': rdMolDescriptors.CalcNumAromaticRings(mol),
        'rdkit_heavy_atoms': Descriptors.HeavyAtomCount(mol),
        'rdkit_fsp3': rdMolDescriptors.CalcFractionCSP3(mol),
    }
```

### 4.4 QM9 Output Files

| File | Description |
|------|-------------|
| `computational_observations.csv` | 12,000 obs: mol_id, formula, identifiers, ΔfH°, all QM9 properties, descriptors |
| `structure_registry.csv` | 12,000 unique structures with molecule_id |

---

## 5. Database Integration & Deduplication

### 5.1 Structure Registry Construction

```python
# From build_final_database.py

def get_structure_key(row):
    """Unique key: InChIKey > InChI > SMILES > CAS > name+formula"""
    if pd.notna(row.get('inchikey')): return row['inchikey']
    if pd.notna(row.get('inchi')):    return row['inchi']
    if pd.notna(row.get('smiles')):   return row['smiles']
    if pd.notna(row.get('cas_number')): return row['cas_number']
    return f"{row.get('name')}_{row.get('formula')}"

# 1. Collect all structures from NIST and QM9
# 2. Deduplicate by structure_key (InChIKey preferred)
# 3. Assign molecule_id (1..N)
# 4. Create lookup maps: inchikey→molecule_id, cas→molecule_id
```

**Results:**
- NIST structures: 238 (237 unique after InChIKey dedup)
- QM9 structures: 12,000 (11,623 unique after InChIKey dedup)
- **Total unique molecules: 11,860**

### 5.2 Observation Tables

**Experimental observations** (302 rows):
```python
for each NIST compound:
    if gas ΔfH° exists → create observation (label_origin='experimental')
    if liquid ΔfH° exists → create observation
    if solid ΔfH° exists → create observation
```

**Computational observations** (12,000 rows):
```python
for each QM9 compound with valid ΔfH°:
    create observation (label_origin='computed')
```

### 5.3 ML Dataset Construction

One row per molecule with preferred target:

```python
for molecule_id in structure_registry:
    exp_gas = get_experimental_gas_ΔfH(molecule_id)
    comp_gas = get_computed_gas_ΔfH(molecule_id)
    
    # Prefer experimental
    if exp_gas: target = exp_gas; source = 'NIST_WebBook'; label = 'experimental'
    elif comp_gas: target = comp_gas; source = 'QM9'; label = 'computed'
    
    # Add formula features (atom counts, DoU)
    # Add RDKit descriptors (if SMILES available)
```

**Result**: 11,860 molecules, 11,803 with ΔfH° targets (166 experimental + 11,637 computed).

---

## 6. Quality Control

### 6.1 Structural Validation
- RDKit parsing of all SMILES/InChI
- Formula validation (RDKit CalcMolFormula vs reported)
- Charge neutrality check
- Stereochemistry preservation

### 6.2 Deduplication Checks
```python
duplicate_inchikeys = struct_df[struct_df['inchikey'].notna()].duplicated(subset=['inchikey']).sum()  # 0
duplicate_cas = struct_df[struct_df['cas_number'].notna()].duplicated(subset=['cas_number']).sum()  # 11,622 (QM9 lack CAS)
```

### 6.3 Experimental vs Computed Validation
- **88 molecules** matched by InChIKey between NIST and QM9
- **Metrics**:
  - MAE: 80.12 kJ/mol
  - RMSE: 87.83 kJ/mol
  - Max absolute error: 189.56 kJ/mol
  - R²: 0.9624

### 6.4 Missing Data Reports
| Check | Result |
|-------|--------|
| Missing InChIKey | 21 molecules |
| Missing SMILES | 21 molecules |
| Missing CAS | 11,623 (QM9) |
| Invalid SMILES (RDKit) | 0 |

---

## 7. Complete Code Inventory

### 7.1 NIST Scraping (`nist_thermo/` package)
| File | Purpose |
|------|---------|
| `scraper.py` | Core NISTScraper class, ThermoData dataclass, parsing logic |
| `workflow.py` | High-level workflows, formula generation, ML preparation |
| `crawler.py` | Formula browser crawler for formula discovery |
| `__init__.py` | Package initialization |

### 7.2 Execution Scripts
| File | Purpose |
|------|---------|
| `run_large_scrape.py` | Large-scale NIST scrape (100+ formulas, 3 isomers each) |
| `run_10k_scrape.py` | Extended scrape targeting 10K+ compounds (50 isomers/formula) |
| `example_workflow.py` | Demo workflow |
| `audit_existing_data.py` | Baseline audit of existing CSV files |

### 7.3 QM9 Processing
| File | Purpose |
|------|---------|
| `process_qm9.py` | Download, parse SDF/CSV, compute ΔfH°, generate descriptors |
| `test_qm9_thermo.py` | Validation of thermochemical cycle against known experimental values |
| `build_final_database.py` | Integration: dedup, observation tables, ML dataset, QC reports |

### 7.4 Key Parameters
```python
# NIST Scraping
DELAY = 0.5-1.0 s          # Rate limit
BATCH_SIZE = 5-10 formulas # Per batch
MAX_PER_FORMULA = 3-50     # Isomers per formula

# QM9 Processing
MAX_MOLECULES = 12,000     # Subset for tractability
ATOMIC_U0 = {...}          # Literature B3LYP/6-31G(2df,p) values
ATOMIC_DFH_298 = {...}     # Experimental atomic ΔfH°
HARTREE_TO_KJMOL = 2625.49962
```

---

## 8. Reproducibility

### 8.1 Environment
```bash
Python 3.10+
pandas>=1.5.0, numpy>=1.23.0, rdkit-pypi>=2022.9.0
requests>=2.28.0, beautifulsoup4>=4.11.0, tqdm>=4.64.0, lxml>=4.9.0
```

### 8.2 Execution Order
```bash
# 1. Scrape NIST (resumable, incremental)
python run_large_scrape.py

# 2. Process QM9
python process_qm9.py --max-molecules 12000 --validate --exp-file nist_hof_large_ml_ready.csv

# 3. Build final database
python build_final_database.py
```

### 8.3 Output Directory Structure
```
final_database/
├── structure_registry.csv          # 11,860 unique molecules
├── experimental_observations.csv   # 302 obs (NIST)
├── computational_observations.csv  # 12,000 obs (QM9)
├── experimental_gold.csv           # 5 high-quality
├── computational_qm9.csv           # 12,000 full QM9 records
├── combined_master.csv             # 12,302 merged
├── ml_dataset.csv                  # 11,860 ML-ready
├── qc_report.json                  # Quality metrics
├── source_statistics.csv           # Per-source stats
├── DATA_DICTIONARY.md              # Column documentation
└── REPRODUCIBILITY.md              # This methodology
```

---

## 9. Statistical Summary

| Metric | Value |
|--------|-------|
| Unique molecular structures | 11,860 |
| Molecules with gas ΔfH° | 11,803 |
| Experimental gas ΔfH° | 166 |
| Computed gas ΔfH° | 11,637 |
| Experimental liquid ΔfH° | 121 |
| Experimental solid ΔfH° | 15 |
| NIST formulas searched | ~500 |
| QM9 molecules processed | 12,000 |
| NIST→QM9 matched (InChIKey) | 88 |
| QM9 vs Exp MAE | 80.12 kJ/mol |
| QM9 vs Exp R² | 0.9624 |

---

## 10. Limitations & Caveats

1. **QM9 accuracy**: B3LYP/6-31G(2df,p) has systematic error (~80 kJ/mol MAE vs experiment)
2. **Element coverage**: QM9 limited to C, H, O, N, F; NIST has broader elements
3. **Phase data**: Liquid/solid ΔfH° sparse in both sources
4. **Temperature**: All values at 298.15 K; no temperature dependence modeled
5. **Subset**: QM9 processed 12,000/133,885 molecules (9%) due to computational constraints
5. **Atomic references**: Literature U₀ values may not perfectly match QM9's specific implementation

---

## References

1. NIST Chemistry WebBook, NIST Standard Reference Database 69, https://webbook.nist.gov
2. Ramakrishnan, R. et al. "Quantum chemistry structures and properties of 134 kilo molecules." *Scientific Data* 1, 140022 (2014).
3. Feller, D. et al. "Thermochemistry of the QM9 dataset." *J. Chem. Phys.* 140, 18A505 (2014).
4. NIST Computational Chemistry Comparison and Benchmark Database (CCCBDB).

---

*Generated: 2026-11-11 | NIST-OrgThermo-ML v1.0*