#!/usr/bin/env python3
"""
Convenience loader for NIST-OrgThermo-ML database.

Usage:
    from scripts.load_database import load_nist_qm9_database
    db = load_nist_qm9_database('data/')
    
    # Access dataframes
    ml_df = db['ml_dataset']
    struct_df = db['structure_registry']
    exp_obs = db['experimental_observations']
    comp_obs = db['computational_observations']
    qc = db['qc_report']
"""
import json
import pandas as pd
from pathlib import Path
from typing import Dict, Any


def load_nist_qm9_database(data_dir: str | Path) -> Dict[str, Any]:
    """
    Load all database files into a dictionary of DataFrames.
    
    Parameters
    ----------
    data_dir : str or Path
        Path to directory containing CSV files.
    
    Returns
    -------
    dict
        Dictionary with keys:
        - 'ml_dataset': ML-ready dataset (11,860 rows)
        - 'structure_registry': Unique molecules (11,860 rows)
        - 'experimental_observations': NIST observations (302 rows)
        - 'computational_observations': QM9 observations (12,000 rows)
        - 'experimental_gold': High-quality experimental (5 rows)
        - 'computational_qm9': Full QM9 records (12,000 rows)
        - 'combined_master': Merged observations (12,302 rows)
        - 'source_statistics': Per-source statistics
        - 'qc_report': Quality control metrics (dict)
    """
    data_dir = Path(data_dir)
    
    files = {
        'ml_dataset': 'ml_dataset.csv',
        'structure_registry': 'structure_registry.csv',
        'experimental_observations': 'experimental_observations.csv',
        'computational_observations': 'computational_observations.csv',
        'experimental_gold': 'experimental_gold.csv',
        'computational_qm9': 'computational_qm9.csv',
        'combined_master': 'combined_master.csv',
        'source_statistics': 'source_statistics.csv',
    }
    
    db = {}
    
    for key, fname in files.items():
        fpath = data_dir / fname
        if fpath.exists():
            db[key] = pd.read_csv(fpath)
        else:
            raise FileNotFoundError(f"Required file not found: {fpath}")
    
    # Load QC report as dict
    qc_path = data_dir / 'qc_report.json'
    if qc_path.exists():
        with open(qc_path) as f:
            db['qc_report'] = json.load(f)
    else:
        db['qc_report'] = {}
    
    return db


def get_molecule_by_inchikey(db: Dict[str, Any], inchikey: str) -> pd.Series:
    """Look up a molecule by InChIKey."""
    struct = db['structure_registry']
    matches = struct[struct['inchikey'] == inchikey]
    if len(matches) == 0:
        raise KeyError(f"No molecule found with InChIKey: {inchikey}")
    if len(matches) > 1:
        raise ValueError(f"Multiple molecules with InChIKey: {inchikey}")
    return matches.iloc[0]


def get_experimental_vs_computed(db: Dict[str, Any], inchikey: str) -> dict:
    """Get experimental vs computed ΔfH° for a molecule."""
    struct = db['structure_registry']
    mol_row = struct[struct['inchikey'] == inchikey]
    
    if len(mol_row) == 0:
        raise KeyError(f"No molecule found with InChIKey: {inchikey}")
    
    mol_id = mol_row.iloc[0]['molecule_id']
    
    # Experimental gas
    exp = db['experimental_observations']
    exp_gas = exp[(exp['molecule_id'] == mol_id) & (exp['phase'] == 'gas')]
    
    # Computed gas
    comp = db['computational_observations']
    comp_gas = comp[comp['molecule_id'] == mol_id]
    
    return {
        'molecule_id': mol_id,
        'inchikey': inchikey,
        'experimental': exp_gas['target_value'].values[0] if len(exp_gas) > 0 else None,
        'computed': comp_gas['target_value'].values[0] if len(comp_gas) > 0 else None,
    }


def split_by_source(db: Dict[str, Any], target_col: str = 'delta_f_H_gas_kJ_mol'):
    """Split ML dataset into experimental and computed subsets."""
    ml = db['ml_dataset']
    
    exp = ml[ml['target_label'] == 'experimental'].copy()
    comp = ml[ml['target_label'] == 'computed'].copy()
    
    return exp, comp


if __name__ == '__main__':
    import sys
    if len(sys.argv) < 2:
        print("Usage: python load_database.py <data_dir>")
        sys.exit(1)
    
    db = load_nist_qm9_database(sys.argv[1])
    print("Loaded database:")
    for key, val in db.items():
        if isinstance(val, pd.DataFrame):
            print(f"  {key}: {len(val)} rows × {len(val.columns)} cols")
        else:
            print(f"  {key}: {type(val).__name__}")
