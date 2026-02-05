"""
Torsion-Score Analysis Script

Extracts base-pair scores from motif_basepairs.csv and links them with
backbone torsion angles from data/torsions/ for correlation analysis.

Usage:
    python torsion_scores_analysis.py [--em N] [--xray N] [--output FILE]

Options:
    --em N       Number of EM base-pairs to sample (default: 500)
    --xray N     Number of X-ray base-pairs to sample (default: 500)
    --output     Output CSV file (default: torsion_scores.csv)
"""

import argparse
import json
import os
import pandas as pd
import numpy as np
from pathlib import Path


# Torsion angles to extract
TORSION_ANGLES = ['alpha', 'beta', 'gamma', 'delta', 'epsilon', 'zeta', 'chi',
                  'eta', 'v1', 'v2', 'v3', 'v4']


def load_torsions(pdb_id: str, torsion_dir: str = "data/torsions") -> dict:
    """Load torsion data for a PDB structure."""
    torsion_file = Path(torsion_dir) / f"{pdb_id}.json"
    if torsion_file.exists():
        with open(torsion_file) as f:
            return json.load(f)
    return {}


def get_sugar_pucker(delta: float) -> str:
    """Classify sugar pucker based on delta torsion angle."""
    if delta is None:
        return "unknown"
    if 55 < delta < 95:
        return "C3'-endo"
    elif 130 < delta < 165:
        return "C2'-endo"
    else:
        return "intermediate"


def get_chi_conformation(chi: float) -> str:
    """Classify glycosidic bond conformation based on chi angle."""
    if chi is None:
        return "unknown"
    # anti: chi ~ -160° (range roughly -180 to -60)
    # syn: chi ~ 60° (range roughly 30 to 90)
    if -180 <= chi <= -60 or 120 <= chi <= 180:
        return "anti"
    elif 30 <= chi <= 90:
        return "syn"
    else:
        return "intermediate"


def analyze_basepairs(input_csv: str, n_em: int, n_xray: int,
                      torsion_dir: str, output_csv: str):
    """
    Main analysis function.

    Args:
        input_csv: Path to motif_basepairs.csv
        n_em: Number of EM base-pairs to sample
        n_xray: Number of X-ray base-pairs to sample
        torsion_dir: Directory containing torsion JSON files
        output_csv: Output CSV file path
    """
    print(f"Loading base-pair data from {input_csv}...")
    df = pd.read_csv(input_csv, low_memory=False)

    print(f"Total base-pairs in file: {len(df):,}")

    # Get available torsion files
    available_pdb_ids = set(
        f.replace('.json', '')
        for f in os.listdir(torsion_dir)
        if f.endswith('.json')
    )
    print(f"Structures with torsion data: {len(available_pdb_ids):,}")

    # Filter to structures that have torsion data
    df = df[df['pdb_id'].isin(available_pdb_ids)]
    print(f"Base-pairs with torsion data available: {len(df):,}")

    # Split by method (handle different naming conventions)
    em_df = df[df['method'].str.upper() == 'EM']
    xray_df = df[df['method'].str.upper().str.contains('X-RAY|XRAY', na=False)]

    print(f"\nEM base-pairs: {len(em_df):,}")
    print(f"X-ray base-pairs: {len(xray_df):,}")

    # Sample
    n_em_actual = min(n_em, len(em_df))
    n_xray_actual = min(n_xray, len(xray_df))

    print(f"\nSampling {n_em_actual} EM and {n_xray_actual} X-ray base-pairs...")

    em_sample = em_df.sample(n=n_em_actual, random_state=42) if n_em_actual > 0 else pd.DataFrame()
    xray_sample = xray_df.sample(n=n_xray_actual, random_state=42) if n_xray_actual > 0 else pd.DataFrame()

    sampled = pd.concat([em_sample, xray_sample], ignore_index=True)
    print(f"Total sampled: {len(sampled):,}")

    # Cache for torsion data
    torsion_cache = {}

    # Build output records
    records = []
    missing_torsions = 0

    for idx, row in sampled.iterrows():
        pdb_id = row['pdb_id']

        # Load torsions (with caching)
        if pdb_id not in torsion_cache:
            torsion_cache[pdb_id] = load_torsions(pdb_id, torsion_dir)
        torsions = torsion_cache[pdb_id]

        # Get residue identifiers
        res1 = row['res1']  # e.g., "Y-G-2-"
        res2 = row['res2']  # e.g., "Y-U-12-"

        # Get torsions for each residue
        t1 = torsions.get(res1, {})
        t2 = torsions.get(res2, {})

        if not t1 and not t2:
            missing_torsions += 1
            continue

        # Build record with base-pair info
        record = {
            'pdb_id': pdb_id,
            'method': row['method'],
            'resolution': row.get('resolution', None),
            'deposition_year': row.get('deposition_year', None),
            'motif_type': row.get('motif_type', None),
            'base_pair': row['base_pair'],
            'res1': res1,
            'res2': res2,
            'lw_notation': row['lw_notation'],
            'bp_type': row['bp_type'],
            'bp_stability_class': row.get('bp_stability_class', None),
            'is_interchain': row.get('is_interchain', False),

            # Score and penalties
            'basepair_score': row['basepair_score'],
            'isPoor': row.get('isPoor', False),
            'geometry_penalty': row.get('geometry_penalty', 0),
            'hbond_penalty': row.get('hbond_penalty', 0),

            # Geometry parameters (from scoring)
            'shear': row.get('shear', None),
            'stretch': row.get('stretch', None),
            'stagger': row.get('stagger', None),
            'buckle': row.get('buckle', None),
            'propeller': row.get('propeller', None),
            'opening': row.get('opening', None),

            # H-bond info
            'hbond_score': row.get('hbond_score', None),
            'number_of_hbonds': row.get('number_of_hbonds', None),

            # Issues
            'issues': row.get('issues', None),
        }

        # Add torsion angles for residue 1
        for angle in TORSION_ANGLES:
            record[f'{angle}_1'] = t1.get(angle, None)

        # Add torsion angles for residue 2
        for angle in TORSION_ANGLES:
            record[f'{angle}_2'] = t2.get(angle, None)

        # Compute differences and averages for key angles
        for angle in ['delta', 'zeta', 'epsilon', 'gamma', 'chi']:
            v1 = t1.get(angle)
            v2 = t2.get(angle)
            if v1 is not None and v2 is not None:
                record[f'{angle}_diff'] = abs(v1 - v2)
                record[f'{angle}_avg'] = (v1 + v2) / 2
            else:
                record[f'{angle}_diff'] = None
                record[f'{angle}_avg'] = None

        # Sugar pucker classification
        record['pucker_1'] = get_sugar_pucker(t1.get('delta'))
        record['pucker_2'] = get_sugar_pucker(t2.get('delta'))
        record['pucker_match'] = record['pucker_1'] == record['pucker_2']

        # Chi conformation
        record['chi_conf_1'] = get_chi_conformation(t1.get('chi'))
        record['chi_conf_2'] = get_chi_conformation(t2.get('chi'))
        record['chi_conf_match'] = record['chi_conf_1'] == record['chi_conf_2']

        records.append(record)

    print(f"\nRecords with torsion data: {len(records):,}")
    print(f"Missing torsion data: {missing_torsions:,}")

    # Create output DataFrame
    output_df = pd.DataFrame(records)

    # Save to CSV
    output_df.to_csv(output_csv, index=False)
    print(f"\nSaved to {output_csv}")

    # Print summary statistics
    print_summary(output_df)

    return output_df


def print_summary(df: pd.DataFrame):
    """Print summary statistics of the analysis."""
    print("\n" + "=" * 70)
    print("SUMMARY STATISTICS")
    print("=" * 70)

    print(f"\nBy Method:")
    print(df.groupby('method')['basepair_score'].agg(['count', 'mean', 'std']).round(2))

    print(f"\nBy Sugar Pucker Match:")
    print(df.groupby('pucker_match')['basepair_score'].agg(['count', 'mean', 'std']).round(2))

    print(f"\nBy Pucker Type (res1):")
    print(df.groupby('pucker_1')['basepair_score'].agg(['count', 'mean', 'std']).round(2))

    print(f"\nCorrelations with basepair_score:")
    numeric_cols = ['delta_diff', 'zeta_diff', 'epsilon_diff', 'gamma_diff', 'chi_diff',
                    'delta_avg', 'chi_avg']
    for col in numeric_cols:
        if col in df.columns:
            valid = df[[col, 'basepair_score']].dropna()
            if len(valid) > 10:
                corr = valid[col].corr(valid['basepair_score'])
                print(f"  {col}: r = {corr:+.4f}")

    # Mean ± SD for torsion differences
    print(f"\nTorsion Difference Statistics (Mean ± SD):")
    for angle in ['delta', 'zeta', 'epsilon', 'gamma', 'chi']:
        col = f'{angle}_diff'
        if col in df.columns:
            valid = df[col].dropna()
            if len(valid) > 0:
                mean = valid.mean()
                std = valid.std()
                print(f"  {angle}: {mean:.1f}° ± {std:.1f}° (n={len(valid)})")
                print(f"    Thresholds: +1SD = {mean + std:.1f}°, +2SD = {mean + 2*std:.1f}°")


def main():
    parser = argparse.ArgumentParser(
        description="Analyze relationship between backbone torsions and base-pair scores"
    )
    parser.add_argument(
        '--input', '-i',
        default='motif_basepairs.csv',
        help='Input CSV file with base-pair scores (default: motif_basepairs.csv)'
    )
    parser.add_argument(
        '--em',
        type=int,
        default=500,
        help='Number of EM base-pairs to sample (default: 500)'
    )
    parser.add_argument(
        '--xray',
        type=int,
        default=500,
        help='Number of X-ray base-pairs to sample (default: 500)'
    )
    parser.add_argument(
        '--torsion-dir',
        default='data/torsions',
        help='Directory containing torsion JSON files (default: data/torsions)'
    )
    parser.add_argument(
        '--output', '-o',
        default='torsion_scores.csv',
        help='Output CSV file (default: torsion_scores.csv)'
    )
    parser.add_argument(
        '--all',
        action='store_true',
        help='Process all base-pairs (ignore --em and --xray limits)'
    )

    args = parser.parse_args()

    if args.all:
        args.em = float('inf')
        args.xray = float('inf')

    analyze_basepairs(
        input_csv=args.input,
        n_em=args.em,
        n_xray=args.xray,
        torsion_dir=args.torsion_dir,
        output_csv=args.output
    )


if __name__ == "__main__":
    main()
