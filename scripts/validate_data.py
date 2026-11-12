#!/usr/bin/env python3
"""
Data validation script for NIST-OrgThermo-ML database.
Run: python scripts/validate_data.py data/
"""
import sys
import json
import pandas as pd
import numpy as np
from pathlib import Path


def validate_structure_registry(df: pd.DataFrame) -> list:
    """Validate structure registry."""
    errors = []
    
    # Expected columns
    expected_cols = [
        'molecule_id', 'source', 'source_id', 'name', 'formula',
        'cas_number', 'inchi', 'inchikey', 'smiles', 'molecular_weight',
        'molecular_charge', 'isotope_info', 'structure_key'
    ]
    missing = set(expected_cols) - set(df.columns)
    if missing:
        errors.append(f"Missing columns: {missing}")
    
    # Row count
    if len(df) != 11860:
        errors.append(f"Row count: {len(df)} != 11860")
    
    # Unique molecule_ids
    if df['molecule_id'].nunique() != len(df):
        errors.append("Duplicate molecule_id values")
    
    # Unique structure_keys (InChIKey preferred)
    inchikey_col = df['inchikey'].dropna()
    if inchikey_col.duplicated().any():
        errors.append(f"Duplicate InChIKeys: {inchikey_col.duplicated().sum()}")
    
    # Source distribution
    source_counts = df['source'].value_counts().to_dict()
    expected_sources = {'QM9': 11623, 'NIST_WebBook': 237}
    for src, exp in expected_sources.items():
        actual = source_counts.get(src, 0)
        if actual != exp:
            errors.append(f"Source {src}: {actual} != {exp}")
    
    return errors


def validate_experimental_observations(df: pd.DataFrame) -> list:
    """Validate experimental observations."""
    errors = []
    
    expected_cols = [
        'observation_id', 'molecule_id', 'source_database', 'target_property',
        'target_value', 'target_unit', 'phase', 'temperature_K', 'pressure',
        'uncertainty', 'method', 'reference', 'source_record_url',
        'quality_flag', 'label_origin'
    ]
    missing = set(expected_cols) - set(df.columns)
    if missing:
        errors.append(f"Missing columns: {missing}")
    
    if len(df) != 302:
        errors.append(f"Row count: {len(df)} != 302")
    
    # Gas phase count
    gas_count = (df['phase'] == 'gas').sum()
    if gas_count != 166:
        errors.append(f"Gas phase observations: {gas_count} != 166")
    
    # Foreign key integrity
    return errors


def validate_computational_observations(df: pd.DataFrame) -> list:
    """Validate computational observations."""
    errors = []
    
    if len(df) != 12000:
        errors.append(f"Row count: {len(df)} != 12000")
    
    # All gas phase
    if not (df['phase'] == 'gas').all():
        errors.append("Not all computational observations are gas phase")
    
    # Method consistency
    method = df['method'].iloc[0]
    if method != 'B3LYP/6-31G(2df,p)//B3LYP/6-31G(2df,p)':
        errors.append(f"Unexpected method: {method}")
    
    return errors


def validate_ml_dataset(df: pd.DataFrame) -> list:
    """Validate ML dataset."""
    errors = []
    
    if len(df) != 11860:
        errors.append(f"Row count: {len(df)} != 11860")
    
    # Target source distribution
    exp_count = (df['target_label'] == 'experimental').sum()
    comp_count = (df['target_label'] == 'computed').sum()
    if exp_count != 166:
        errors.append(f"Experimental targets: {exp_count} != 166")
    if comp_count != 11637:
        errors.append(f"Computed targets: {comp_count} != 11637")
    
    # Target value range
    target_min = df['delta_f_H_gas_kJ_mol'].min()
    target_max = df['delta_f_H_gas_kJ_mol'].max()
    if target_min < -1000 or target_max > 1000:
        errors.append(f"Target range suspicious: [{target_min}, {target_max}]")
    
    return errors


def validate_qc_report(qc: dict) -> list:
    """Validate QC report matches data."""
    errors = []
    
    # Structure registry
    sr = qc.get('structure_registry', {})
    if sr.get('total_unique_molecules') != 11860:
        errors.append(f"QC total_unique_molecules: {sr.get('total_unique_molecules')} != 11860")
    
    # Validation metrics
    val = qc.get('validation', {}).get('qm9_vs_experimental', {})
    if abs(val.get('mean_absolute_error_kJ_mol', 0) - 80.12) > 0.01:
        errors.append(f"QC MAE: {val.get('mean_absolute_error_kJ_mol')} != 80.12")
    if abs(val.get('r_squared', 0) - 0.9624) > 0.0001:
        errors.append(f"QC R²: {val.get('r_squared')} != 0.9624")
    
    return errors


def validate_foreign_keys(data_dir: Path) -> list:
    """Validate foreign key relationships."""
    errors = []
    
    struct = pd.read_csv(data_dir / 'structure_registry.csv')
    valid_ids = set(struct['molecule_id'])
    
    # Check experimental observations
    exp = pd.read_csv(data_dir / 'experimental_observations.csv')
    invalid = set(exp['molecule_id']) - valid_ids
    if invalid:
        errors.append(f"Experimental observations: {len(invalid)} invalid molecule_ids")
    
    # Check computational observations
    comp = pd.read_csv(data_dir / 'computational_observations.csv')
    invalid = set(comp['molecule_id']) - valid_ids
    if invalid:
        errors.append(f"Computational observations: {len(invalid)} invalid molecule_ids")
    
    # Check ML dataset
    ml = pd.read_csv(data_dir / 'ml_dataset.csv')
    invalid = set(ml['molecule_id']) - valid_ids
    if invalid:
        errors.append(f"ML dataset: {len(invalid)} invalid molecule_ids")
    
    return errors


def main():
    if len(sys.argv) < 2:
        print("Usage: python validate_data.py <data_directory>")
        sys.exit(1)
    
    data_dir = Path(sys.argv[1])
    if not data_dir.exists():
        print(f"Data directory not found: {data_dir}")
        sys.exit(1)
    
    print(f"Validating data in: {data_dir}")
    print("=" * 60)
    
    all_errors = []
    
    # Load and validate each file
    validators = [
        ("structure_registry.csv", validate_structure_registry),
        ("experimental_observations.csv", validate_experimental_observations),
        ("computational_observations.csv", validate_computational_observations),
        ("ml_dataset.csv", validate_ml_dataset),
    ]
    
    for fname, validator in validators:
        fpath = data_dir / fname
        if not fpath.exists():
            all_errors.append(f"Missing file: {fname}")
            continue
        
        df = pd.read_csv(fpath)
        errors = validator(df)
        if errors:
            print(f"\n❌ {fname}:")
            for e in errors:
                print(f"  - {e}")
            all_errors.extend([f"{fname}: {e}" for e in errors])
        else:
            print(f"✓ {fname}: {len(df)} rows")
    
    # Validate QC report
    qc_path = data_dir / 'qc_report.json'
    if qc_path.exists():
        with open(qc_path) as f:
            qc = json.load(f)
        errors = validate_qc_report(qc)
        if errors:
            print(f"\n❌ qc_report.json:")
            for e in errors:
                print(f"  - {e}")
            all_errors.extend([f"qc_report.json: {e}" for e in errors])
        else:
            print(f"✓ qc_report.json: valid")
    else:
        all_errors.append("Missing qc_report.json")
    
    # Validate foreign keys
    fk_errors = validate_foreign_keys(data_dir)
    if fk_errors:
        print(f"\n❌ Foreign key validation:")
        for e in fk_errors:
            print(f"  - {e}")
        all_errors.extend(fk_errors)
    else:
        print(f"✓ Foreign keys: all valid")
    
    print("\n" + "=" * 60)
    if all_errors:
        print(f"❌ VALIDATION FAILED: {len(all_errors)} error(s)")
        sys.exit(1)
    else:
        print("✅ ALL VALIDATIONS PASSED")
        sys.exit(0)


if __name__ == '__main__':
    main()
