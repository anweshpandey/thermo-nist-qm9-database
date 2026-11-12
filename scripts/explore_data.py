#!/usr/bin/env python3
"""
Quick exploration and visualization scripts for NIST-OrgThermo-ML database.

Run: python scripts/explore_data.py data/
"""
import sys
import pandas as pd
import numpy as np
from pathlib import Path

try:
    import matplotlib.pyplot as plt
    import seaborn as sns
    HAS_PLOTTING = True
except ImportError:
    HAS_PLOTTING = False
    print("Note: matplotlib/seaborn not installed. Install with: pip install matplotlib seaborn")


def load_db(data_dir):
    """Load all data."""
    data_dir = Path(data_dir)
    return {
        'ml': pd.read_csv(data_dir / 'ml_dataset.csv'),
        'struct': pd.read_csv(data_dir / 'structure_registry.csv'),
        'exp_obs': pd.read_csv(data_dir / 'experimental_observations.csv'),
        'comp_obs': pd.read_csv(data_dir / 'computational_observations.csv'),
        'qc': pd.read_json(data_dir / 'qc_report.json'),
    }


def print_summary(db):
    """Print database summary."""
    ml = db['ml']
    struct = db['struct']
    exp_obs = db['exp_obs']
    comp_obs = db['comp_obs']
    
    print("=" * 60)
    print("NIST-OrgThermo-ML Database Summary")
    print("=" * 60)
    print(f"\nStructure Registry: {len(struct)} unique molecules")
    print(f"  Sources: {struct['source'].value_counts().to_dict()}")
    print(f"  With InChIKey: {struct['inchikey'].notna().sum()}")
    print(f"  With SMILES: {struct['smiles'].notna().sum()}")
    print(f"  With CAS: {struct['cas_number'].notna().sum()}")
    
    print(f"\nML Dataset: {len(ml)} molecules")
    print(f"  Experimental targets: {(ml['target_label'] == 'experimental').sum()}")
    print(f"  Computed targets: {(ml['target_label'] == 'computed').sum()}")
    print(f"  Target range: [{ml['delta_f_H_gas_kJ_mol'].min():.1f}, {ml['delta_f_H_gas_kJ_mol'].max():.1f}] kJ/mol")
    print(f"  Target mean: {ml['delta_f_H_gas_kJ_mol'].mean():.1f} kJ/mol")
    print(f"  Target std: {ml['delta_f_H_gas_kJ_mol'].std():.1f} kJ/mol")
    
    print(f"\nExperimental Observations: {len(exp_obs)}")
    print(f"  Phases: {exp_obs['phase'].value_counts().to_dict()}")
    print(f"  Methods: {exp_obs['method'].value_counts().to_dict()}")
    print(f"  With uncertainty: {exp_obs['uncertainty'].notna().sum()}")
    
    print(f"\nComputational Observations: {len(comp_obs)}")
    print(f"  Method: {comp_obs['method'].iloc[0]}")
    print(f"  Target range: [{comp_obs['target_value'].min():.1f}, {comp_obs['target_value'].max():.1f}] kJ/mol")


def analyze_elemental_composition(db):
    """Analyze elemental composition."""
    ml = db['ml']
    
    print("\n" + "=" * 60)
    print("Elemental Composition Analysis")
    print("=" * 60)
    
    for elem in ['C', 'H', 'O', 'N', 'F', 'Cl', 'Br', 'I']:
        col = f'n_{elem}'
        if col in ml.columns:
            non_zero = (ml[col] > 0).sum()
            max_val = ml[col].max()
            print(f"  {elem}: {non_zero} molecules ({non_zero/len(ml)*100:.1f}%), max count: {max_val}")


def analyze_molecular_weight(db):
    """Analyze molecular weight distribution."""
    ml = db['ml']
    
    print("\n" + "=" * 60)
    print("Molecular Weight Distribution")
    print("=" * 60)
    
    mw = ml['molecular_weight'].dropna()
    print(f"  Count: {len(mw)}")
    print(f"  Range: [{mw.min():.1f}, {mw.max():.1f}] g/mol")
    print(f"  Mean: {mw.mean():.1f} g/mol")
    print(f"  Median: {mw.median():.1f} g/mol")
    print(f"  Std: {mw.std():.1f} g/mol")
    
    # By source
    for label in ['experimental', 'computed']:
        subset = ml[ml['target_label'] == label]['molecular_weight'].dropna()
        print(f"  {label}: n={len(subset)}, mean={subset.mean():.1f}, median={subset.median():.1f}")


def compare_exp_vs_comp(db):
    """Compare experimental vs computed for matched molecules."""
    struct = db['struct']
    exp_obs = db['exp_obs']
    comp_obs = db['comp_obs']
    
    print("\n" + "=" * 60)
    print("Experimental vs Computed Comparison")
    print("=" * 60)
    
    # Find matched by InChIKey
    exp_struct = struct[struct['source'] == 'NIST_WebBook'][['molecule_id', 'inchikey']].dropna()
    comp_struct = struct[struct['source'] == 'QM9'][['molecule_id', 'inchikey']].dropna()
    
    matched = pd.merge(exp_struct, comp_struct, on='inchikey', suffixes=('_exp', '_comp'))
    print(f"  Matched by InChIKey: {len(matched)} molecules")
    
    if len(matched) > 0:
        # Get values
        exp_vals = exp_obs[exp_obs['molecule_id'].isin(matched['molecule_id_exp']) & 
                          (exp_obs['phase'] == 'gas')].set_index('molecule_id')['target_value']
        comp_vals = comp_obs[comp_obs['molecule_id'].isin(matched['molecule_id_comp'])].set_index('molecule_id')['target_value']
        
        # Align by InChIKey
        comparison = []
        for _, row in matched.iterrows():
            exp_v = exp_vals.get(row['molecule_id_exp'])
            comp_v = comp_vals.get(row['molecule_id_comp'])
            if exp_v is not None and comp_v is not None:
                comparison.append({'inchikey': row['inchikey'], 'exp': exp_v, 'comp': comp_v})
        
        if comparison:
            comp_df = pd.DataFrame(comparison)
            comp_df['diff'] = comp_df['comp'] - comp_df['exp']
            comp_df['abs_diff'] = comp_df['diff'].abs()
            
            print(f"  Valid pairs: {len(comp_df)}")
            print(f"  MAE: {comp_df['abs_diff'].mean():.2f} kJ/mol")
            print(f"  RMSE: {np.sqrt((comp_df['diff']**2).mean()):.2f} kJ/mol")
            print(f"  Max abs error: {comp_df['abs_diff'].max():.2f} kJ/mol")
            print(f"  R²: {np.corrcoef(comp_df['exp'], comp_df['comp'])[0,1]**2:.4f}")
            print(f"  Mean exp: {comp_df['exp'].mean():.2f}")
            print(f"  Mean comp: {comp_df['comp'].mean():.2f}")


def plot_distributions(db, output_dir=None):
    """Generate distribution plots."""
    if not HAS_PLOTTING:
        print("\nSkipping plots (matplotlib/seaborn not available)")
        return
    
    ml = db['ml']
    exp_obs = db['exp_obs']
    comp_obs = db['comp_obs']
    
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    
    # 1. Target distribution by source
    ax = axes[0, 0]
    for label, color in [('experimental', 'blue'), ('computed', 'orange')]:
        subset = ml[ml['target_label'] == label]['delta_f_H_gas_kJ_mol']
        ax.hist(subset, bins=50, alpha=0.6, label=label, color=color, density=True)
    ax.set_xlabel('ΔfH° (kJ/mol)')
    ax.set_ylabel('Density')
    ax.set_title('Target Distribution by Source')
    ax.legend()
    
    # 2. Molecular weight distribution
    ax = axes[0, 1]
    for label, color in [('experimental', 'blue'), ('computed', 'orange')]:
        subset = ml[ml['target_label'] == label]['molecular_weight']
        ax.hist(subset, bins=50, alpha=0.6, label=label, color=color, density=True)
    ax.set_xlabel('Molecular Weight (g/mol)')
    ax.set_ylabel('Density')
    ax.set_title('Molecular Weight Distribution')
    ax.legend()
    
    # 3. Degree of Unsaturation
    ax = axes[0, 2]
    for label, color in [('experimental', 'blue'), ('computed', 'orange')]:
        subset = ml[ml['target_label'] == label]['DoU']
        ax.hist(subset, bins=20, alpha=0.6, label=label, color=color, density=True)
    ax.set_xlabel('Degree of Unsaturation')
    ax.set_ylabel('Density')
    ax.set_title('DoU Distribution')
    ax.legend()
    
    # 4. Experimental method distribution
    ax = axes[1, 0]
    method_counts = exp_obs['method'].value_counts()
    method_counts.plot(kind='bar', ax=ax, color='steelblue')
    ax.set_xlabel('NIST Method')
    ax.set_ylabel('Count')
    ax.set_title('Experimental Methods')
    ax.tick_params(axis='x', rotation=45)
    
    # 5. Experimental phase distribution
    ax = axes[1, 1]
    phase_counts = exp_obs['phase'].value_counts()
    phase_counts.plot(kind='bar', ax=ax, color='coral')
    ax.set_xlabel('Phase')
    ax.set_ylabel('Count')
    ax.set_title('Experimental Phases')
    
    # 6. Target vs Molecular Weight scatter
    ax = axes[1, 2]
    exp_ml = ml[ml['target_label'] == 'experimental']
    comp_ml = ml[ml['target_label'] == 'computed']
    ax.scatter(exp_ml['molecular_weight'], exp_ml['delta_f_H_gas_kJ_mol'], 
               alpha=0.6, label='Experimental', color='blue', s=20)
    ax.scatter(comp_ml['molecular_weight'], comp_ml['delta_f_H_gas_kJ_mol'], 
               alpha=0.3, label='Computed', color='orange', s=5)
    ax.set_xlabel('Molecular Weight (g/mol)')
    ax.set_ylabel('ΔfH° (kJ/mol)')
    ax.set_title('Target vs Molecular Weight')
    ax.legend()
    
    plt.tight_layout()
    
    if output_dir:
        out_path = Path(output_dir) / 'database_exploration.png'
        plt.savefig(out_path, dpi=150, bbox_inches='tight')
        print(f"\nPlot saved to: {out_path}")
    else:
        plt.show()


def main():
    if len(sys.argv) < 2:
        print("Usage: python explore_data.py <data_directory> [output_dir]")
        sys.exit(1)
    
    data_dir = sys.argv[1]
    output_dir = sys.argv[2] if len(sys.argv) > 2 else None
    
    db = load_db(data_dir)
    print_summary(db)
    analyze_elemental_composition(db)
    analyze_molecular_weight(db)
    compare_exp_vs_comp(db)
    
    if HAS_PLOTTING:
        plot_distributions(db, output_dir)
    
    print("\n✅ Exploration complete")


if __name__ == '__main__':
    main()
