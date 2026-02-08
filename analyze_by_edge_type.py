#!/usr/bin/env python3
"""
Analyze geometric and H-bond parameters by edge type (Leontis-Westhof notation).

This script efficiently processes data to compare parameter distributions
across different edge types (cWW, tWW, cWH, etc.) without overloading memory.
"""

#python analyze_by_edge_type.py

import json
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from collections import defaultdict
import sys
from config import Config
from g_quads import g_quads


def normalize_bp_type(bp_type: str) -> str:
    """
    Normalize bp_type so A-U == U-A, G-C == C-G, etc.
    """
    if not bp_type or "-" not in bp_type:
        return bp_type or "unknown"
    parts = bp_type.split("-")
    if len(parts) != 2:
        return bp_type
    return "-".join(sorted(parts))


def is_adjacent_pair(res_1: str, res_2: str) -> bool:
    """
    Return True if the two residues are adjacent in sequence (same chain, residue numbers differ by 1).
    Excludes base pairs between consecutive residues (e.g. 74-75) which are often stacking/artifacts.
    """
    if not res_1 or not res_2:
        return False
    parts1 = res_1.split('-')
    parts2 = res_2.split('-')
    if len(parts1) < 3 or len(parts2) < 3:
        return False
    try:
        chain1, num1 = parts1[0], int(parts1[2])
        chain2, num2 = parts2[0], int(parts2[2])
    except (ValueError, IndexError):
        return False
    if chain1 != chain2:
        return False
    return abs(num1 - num2) == 1


def load_basepair_sample(basepair_dir: Path, max_files: int = None, sample_fraction: float = None,
                         include_pdb_ids: Optional[set] = None) -> List[Dict]:
    """
    Load base pair data efficiently.

    Args:
        basepair_dir: Directory containing base pair JSON files
        max_files: Maximum number of files to process (None = process all)
        sample_fraction: Fraction of files to sample (None = use max_files)
        include_pdb_ids: Set of PDB IDs (uppercase) to include. If provided, only these are processed.

    Returns:
        List of base pair dictionaries
    """
    json_files = list(basepair_dir.glob("*.json"))

    # Filter to only allowed PDB IDs if provided
    if include_pdb_ids:
        before = len(json_files)
        json_files = [p for p in json_files if p.stem.upper() in include_pdb_ids]
        print(f"Filtered to {len(json_files)} base pair files from allowlist ({before} total available)")
    
    if sample_fraction:
        n_files = max(100, int(len(json_files) * sample_fraction))
        json_files = np.random.choice(json_files, size=min(n_files, len(json_files)), replace=False)
    elif max_files is None:
        # Process all files
        pass  # Use all files
    else:
        json_files = json_files[:max_files]
    
    print(f"Loading base pairs from {len(json_files)} files...")
    
    all_basepairs = []
    files_processed = 0
    
    for bp_file in json_files:
        try:
            with open(bp_file, 'r') as f:
                data = json.load(f)
            
            if isinstance(data, list):
                base_pairs = data
            elif isinstance(data, dict):
                base_pairs = data.get('base_pairs', [])
            else:
                continue
            
            # Normalize bp_type for symmetry
            for bp in base_pairs:
                bp["bp_type_norm"] = normalize_bp_type(bp.get("bp_type", "unknown"))
            all_basepairs.extend(base_pairs)
            files_processed += 1
            
            if files_processed % 500 == 0:
                print(f"  Processed {files_processed}/{len(json_files)} files, collected {len(all_basepairs):,} base pairs...")
        
        except Exception as e:
            continue
    
    print(f"✓ Loaded {len(all_basepairs):,} base pairs from {files_processed} files")
    return all_basepairs


def create_basepair_lookup(basepairs: List[Dict]) -> Dict:
    """
    Create a lookup dictionary: (res_1, res_2) -> base pair info including lw.
    """
    lookup = {}
    for bp in basepairs:
        res_1 = bp.get('res_1', '')
        res_2 = bp.get('res_2', '')
        if res_1 and res_2:
            # Create canonical key (sorted)
            key = tuple(sorted([res_1, res_2]))
            lookup[key] = {
                'lw': bp.get('lw', 'unknown'),
                'bp_type': normalize_bp_type(bp.get('bp_type', 'unknown')),
                'shear': bp.get('shear'),
                'stretch': bp.get('stretch'),
                'stagger': bp.get('stagger'),
                'buckle': bp.get('buckle'),
                'propeller': bp.get('propeller'),
                'opening': bp.get('opening'),
                'hbond_score': bp.get('hbond_score')
            }
    return lookup


def load_hbond_sample_with_lw(hbond_dir: Path, bp_lookup: Dict, max_files: int = None,
                                sample_fraction: float = None, max_hbonds_per_file: int = None,
                                include_pdb_ids: Optional[set] = None) -> pd.DataFrame:
    """
    Load H-bond data and link with base pair edge types.
    Only loads H-bonds that belong to known base pairs.

    Args:
        hbond_dir: Directory containing H-bond CSV files
        bp_lookup: Base pair lookup dictionary
        max_files: Maximum number of files to process (None = process all)
        sample_fraction: Fraction of files to sample
        max_hbonds_per_file: Maximum H-bonds to load per file (for memory efficiency)
        include_pdb_ids: Set of PDB IDs (uppercase) to include. If provided, only these are processed.

    Returns:
        DataFrame with H-bond data and edge type information
    """
    csv_files = list(hbond_dir.glob("*.csv"))

    # Filter to only allowed PDB IDs if provided
    if include_pdb_ids:
        before = len(csv_files)
        csv_files = [p for p in csv_files if p.stem.upper() in include_pdb_ids]
        print(f"Filtered to {len(csv_files)} H-bond files from allowlist ({before} total available)")
    
    if sample_fraction:
        n_files = max(100, int(len(csv_files) * sample_fraction))
        csv_files = np.random.choice(csv_files, size=min(n_files, len(csv_files)), replace=False)
    elif max_files is None:
        # Process all files
        pass  # Use all files
    else:
        csv_files = csv_files[:max_files]
    
    print(f"\nLoading H-bonds from {len(csv_files)} files...")
    
    all_hbonds = []
    files_processed = 0
    total_collected = 0
    
    for csv_file in csv_files:
        try:
            df = pd.read_csv(csv_file)
            
            # Filter to RNA-RNA only
            df = df[(df['res_type_1'] == 'RNA') & (df['res_type_2'] == 'RNA')]
            
            if len(df) == 0:
                continue
            
            # Sample if too many (only if max_hbonds_per_file is set)
            if max_hbonds_per_file and len(df) > max_hbonds_per_file:
                df = df.sample(n=max_hbonds_per_file, random_state=42)
            
            # Track statistics before filtering
            total_before = len(df)
            
            # Add edge type information
            df['lw'] = 'unknown'
            df['bp_type'] = 'unknown'
            
            matched_count = 0
            for idx, row in df.iterrows():
                pair_key = tuple(sorted([row['res_1'], row['res_2']]))
                if pair_key in bp_lookup:
                    df.at[idx, 'lw'] = bp_lookup[pair_key]['lw']
                    df.at[idx, 'bp_type'] = bp_lookup[pair_key]['bp_type']
                    matched_count += 1
            
            # Only keep H-bonds with known edge types
            df = df[df['lw'] != 'unknown']
            
            # Debug info (only print occasionally to avoid spam)
            if files_processed % 1000 == 0 and files_processed > 0:
                print(f"    File {files_processed}: {total_before} total H-bonds, {matched_count} matched to base pairs, {len(df)} kept")
            
            if len(df) > 0:
                all_hbonds.append(df)
                total_collected += len(df)
                files_processed += 1
                
                if files_processed % 500 == 0:
                    print(f"  Processed {files_processed}/{len(csv_files)} files, collected {total_collected:,} H-bonds with edge types...")
        
        except Exception as e:
            continue
    
    if not all_hbonds:
        print("Warning: No H-bonds with edge type information found")
        return pd.DataFrame()
    
    combined_df = pd.concat(all_hbonds, ignore_index=True)
    print(f"✓ Loaded {len(combined_df):,} H-bonds with edge type information")
    
    return combined_df


def group_geometry_by_bp_and_edge(basepairs: List[Dict], min_edge_count: int = 1000) -> Dict:
    grouped = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    counts = defaultdict(lambda: defaultdict(int))
    totals = defaultdict(int)

    for bp in basepairs:
        bp_type = bp.get("bp_type_norm") or normalize_bp_type(bp.get("bp_type", "unknown"))
        lw = bp.get('lw', 'unknown')
        if not bp_type or not lw or lw == 'unknown':
            continue
        totals[bp_type] += 1
        counts[bp_type][lw] += 1
        for param in ['shear', 'stretch', 'stagger', 'buckle', 'propeller', 'opening', 'hbond_score']:
            value = bp.get(param)
            if value is not None:
                try:
                    grouped[bp_type][lw][param].append(float(value))
                except (ValueError, TypeError):
                    pass

    results = {}
    for bp_type, edges in grouped.items():
        results[bp_type] = {}
        total_bp = totals[bp_type]
        other_bucket = defaultdict(list)
        for lw, params in edges.items():
            edge_count = counts[bp_type][lw]
            target = params if (edge_count >= min_edge_count or total_bp < min_edge_count) else other_bucket
            target_params = {}
            for param, values in params.items():
                if len(values) > 0:
                    # Only use ideal_value=0.0 for geometry params, not hbond_score
                    ideal = 0.0 if param != 'hbond_score' else None
                    target_params[param] = calculate_statistics_by_group(values, f"{bp_type}_{lw}_{param}", ideal_value=ideal)
            if target is params:
                results[bp_type][lw] = target_params
        # add other if populated
        if other_bucket:
            other_stats = {}
            for param, values in other_bucket.items():
                if len(values) > 0:
                    ideal = 0.0 if param != 'hbond_score' else None
                    other_stats[param] = calculate_statistics_by_group(values, f"{bp_type}_OTHER_{param}", ideal_value=ideal)
            if other_stats:
                results[bp_type]['_OTHER'] = other_stats
    return results


def group_hbonds_by_bp_and_edge(df: pd.DataFrame, min_edge_count: int = 1000) -> Dict:
    if len(df) == 0:
        return {}
    df = df.copy()
    df['bp_type_norm'] = df['bp_type'].apply(normalize_bp_type)
    df['dihedral_norm'] = ((df['dihedral_angle'] + 180) % 360) - 180
    df['angle_min'] = df[['angle_1', 'angle_2']].min(axis=1)

    def separate_dihedral_groups(dihedral_values: List[float]) -> Tuple[List[float], List[float]]:
        cis_values = []
        trans_values = []
        for val in dihedral_values:
            if pd.isna(val) or not np.isfinite(val):
                continue
            val_float = float(val)
            if -50.0 <= val_float <= 50.0:
                cis_values.append(val_float)
            elif abs(val_float) >= 140.0:
                trans_values.append(val_float)
        return cis_values, trans_values

    grouped = {}
    for bp_type, bp_df in df.groupby('bp_type_norm'):
        total_bp = len(bp_df)
        edges_result = {}
        edge_counts = bp_df['lw'].value_counts().to_dict()
        other_frames = []

        for lw, lw_df in bp_df.groupby('lw'):
            if lw == 'unknown' or pd.isna(lw):
                continue
            edge_count = edge_counts.get(lw, 0)
            target_df = lw_df if (edge_count >= min_edge_count or total_bp < min_edge_count) else None
            if target_df is None:
                other_frames.append(lw_df)
                continue

            all_dihedrals = target_df['dihedral_norm'].tolist()
            cis_vals, trans_vals = separate_dihedral_groups(all_dihedrals)

            edge_stats = {
                'distance': calculate_statistics_by_group(target_df['distance'].tolist(), f"{bp_type}_{lw}_distance"),
                'angle_1': calculate_statistics_by_group(target_df['angle_1'].tolist(), f"{bp_type}_{lw}_angle1"),
                'angle_2': calculate_statistics_by_group(target_df['angle_2'].tolist(), f"{bp_type}_{lw}_angle2"),
                'angle_min': calculate_statistics_by_group(target_df['angle_min'].tolist(), f"{bp_type}_{lw}_angle_min"),
                'dihedral': calculate_statistics_by_group(all_dihedrals, f"{bp_type}_{lw}_dihedral"),
                'quality': calculate_statistics_by_group(target_df['score'].tolist(), f"{bp_type}_{lw}_quality")
            }
            edge_stats['dihedral_cis'] = calculate_statistics_by_group(cis_vals, f"{bp_type}_{lw}_dihedral_cis") if cis_vals else {"error": "No CIS dihedral data", "count": 0}
            edge_stats['dihedral_trans'] = calculate_statistics_by_group(trans_vals, f"{bp_type}_{lw}_dihedral_trans") if trans_vals else {"error": "No TRANS dihedral data", "count": 0}

            edges_result[lw] = edge_stats

        # other bucket
        if other_frames:
            other_df = pd.concat(other_frames, ignore_index=True)
            all_dihedrals = other_df['dihedral_norm'].tolist()
            cis_vals, trans_vals = separate_dihedral_groups(all_dihedrals)
            edge_stats = {
                'distance': calculate_statistics_by_group(other_df['distance'].tolist(), f"{bp_type}_OTHER_distance"),
                'angle_1': calculate_statistics_by_group(other_df['angle_1'].tolist(), f"{bp_type}_OTHER_angle1"),
                'angle_2': calculate_statistics_by_group(other_df['angle_2'].tolist(), f"{bp_type}_OTHER_angle2"),
                'angle_min': calculate_statistics_by_group(other_df['angle_min'].tolist(), f"{bp_type}_OTHER_angle_min"),
                'dihedral': calculate_statistics_by_group(all_dihedrals, f"{bp_type}_OTHER_dihedral"),
                'quality': calculate_statistics_by_group(other_df['score'].tolist(), f"{bp_type}_OTHER_quality")
            }
            edge_stats['dihedral_cis'] = calculate_statistics_by_group(cis_vals, f"{bp_type}_OTHER_dihedral_cis") if cis_vals else {"error": "No CIS dihedral data", "count": 0}
            edge_stats['dihedral_trans'] = calculate_statistics_by_group(trans_vals, f"{bp_type}_OTHER_dihedral_trans") if trans_vals else {"error": "No TRANS dihedral data", "count": 0}
            edges_result['_OTHER'] = edge_stats

        grouped[bp_type] = edges_result

    return grouped


def calculate_statistics_by_group(values: List[float], group_name: str, ideal_value: float = None) -> Dict:
    """Calculate statistics for a group of values.
    
    Args:
        values: List of values to analyze
        group_name: Name of the group (for reference)
        ideal_value: Ideal value for this parameter (None = no ideal, use mean-based)
    
    Returns:
        Dictionary with statistics including both mean-based and ideal-based calculations
    """
    if not values or len(values) == 0:
        return {"error": "No data available", "count": 0}
    
    values_array = np.array(values)
    values_array = values_array[~np.isnan(values_array)]
    values_array = values_array[np.isfinite(values_array)]
    
    if len(values_array) == 0:
        return {"error": "No valid data", "count": 0}
    
    mean_val = float(np.mean(values_array))
    std_val = float(np.std(values_array))
    
    result = {
        "count": len(values_array),
        "mean": mean_val,
        "std": std_val,
        "median": float(np.median(values_array)),
        "min": float(np.min(values_array)),
        "max": float(np.max(values_array)),
        "p2.5": float(np.percentile(values_array, 2.5)),
        "p5": float(np.percentile(values_array, 5)),
        "p25": float(np.percentile(values_array, 25)),
        "p50": float(np.percentile(values_array, 50)),
        "p75": float(np.percentile(values_array, 75)),
        "p95": float(np.percentile(values_array, 95)),
        "p97.5": float(np.percentile(values_array, 97.5)),
        # Mean-based deviations
        "mean_minus_05sd": float(mean_val - 0.5 * std_val),
        "mean_plus_05sd": float(mean_val + 0.5 * std_val),
        "mean_minus_1sd": float(mean_val - 1.0 * std_val),
        "mean_plus_1sd": float(mean_val + 1.0 * std_val)
    }
    
    # Add ideal-based calculations if ideal_value is provided
    if ideal_value is not None:
        # For geometry, we want deviation from ideal (not from mean):
        # std_from_ideal = RMS from the ideal value = sqrt(mean((x - ideal)^2))
        # This properly captures how far values sit from the target of 0.
        ideal_rms = float(np.sqrt(np.mean((values_array - ideal_value) ** 2)))
        
        result["ideal"] = ideal_value
        result["sd_from_ideal"] = ideal_rms
        result["ideal_minus_05sd"] = float(ideal_value - 0.5 * ideal_rms)
        result["ideal_plus_05sd"] = float(ideal_value + 0.5 * ideal_rms)
        result["ideal_minus_1sd"] = float(ideal_value - 1.0 * ideal_rms)
        result["ideal_plus_1sd"] = float(ideal_value + 1.0 * ideal_rms)
        result["ideal_minus_15sd"] = float(ideal_value - 1.5 * ideal_rms)
        result["ideal_plus_15sd"] = float(ideal_value + 1.5 * ideal_rms)
    
    return result



def analyze_geometry_by_edge_type(basepairs: List[Dict]) -> Dict:
    """Analyze geometry parameters grouped by bp_type and edge type."""
    print("\nAnalyzing geometry parameters by base pair type and edge type...")
    return group_geometry_by_bp_and_edge(basepairs, min_edge_count=1000)


def analyze_hbonds_by_edge_type(df: pd.DataFrame) -> Dict:
    """Analyze H-bond parameters grouped by bp_type and edge type."""
    print("\nAnalyzing H-bond parameters by base pair type and edge type...")
    return group_hbonds_by_bp_and_edge(df, min_edge_count=1000)


def print_summary(geometry_results: Dict, hbond_results: Dict):
    """Print summary of findings."""
    print("\n" + "="*80)
    print("SUMMARY BY EDGE TYPE")
    print("="*80)
    
    # Count base pairs by edge type
    print("\nBase Pair Counts by Edge Type:")
    lw_counts = {}
    for lw, params in geometry_results.items():
        if lw.startswith('_'):
            continue
        if 'shear' in params:
            lw_counts[lw] = params['shear'].get('count', 0)
    
    for lw, count in sorted(lw_counts.items(), key=lambda x: x[1], reverse=True):
        print(f"  {lw:6s}: {count:8,} base pairs")
    
    # Canonical vs non-canonical comparison
    if '_CANONICAL' in geometry_results and '_NONCANONICAL' in geometry_results:
        print("\n" + "-"*80)
        print("CANONICAL vs NON-CANONICAL COMPARISON (Geometry)")
        print("-"*80)
        
        for param in ['shear', 'stretch', 'stagger', 'buckle', 'propeller', 'opening', 'hbond_score']:
            if param in geometry_results['_CANONICAL'] and param in geometry_results['_NONCANONICAL']:
                can = geometry_results['_CANONICAL'][param]
                non = geometry_results['_NONCANONICAL'][param]
                
                if 'error' not in can and 'error' not in non:
                    print(f"\n{param.upper()}:")
                    print(f"  Canonical:    mean={can['mean']:7.3f}, std={can['std']:7.3f}")
                    print(f"                Mean ±0.5 SD: [{can['mean_minus_05sd']:7.3f}, {can['mean_plus_05sd']:7.3f}]")
                    print(f"                Mean ±1.0 SD: [{can['mean_minus_1sd']:7.3f}, {can['mean_plus_1sd']:7.3f}]")
                    if 'ideal' in can:
                        print(f"                Ideal (0) ±0.5 SD: [{can['ideal_minus_05sd']:7.3f}, {can['ideal_plus_05sd']:7.3f}]")
                        print(f"                Ideal (0) ±1.0 SD: [{can['ideal_minus_1sd']:7.3f}, {can['ideal_plus_1sd']:7.3f}] (n={can['count']:,})")
                    else:
                        print(f"                (n={can['count']:,})")
                    print(f"  Non-canonical: mean={non['mean']:7.3f}, std={non['std']:7.3f}")
                    print(f"                Mean ±0.5 SD: [{non['mean_minus_05sd']:7.3f}, {non['mean_plus_05sd']:7.3f}]")
                    print(f"                Mean ±1.0 SD: [{non['mean_minus_1sd']:7.3f}, {non['mean_plus_1sd']:7.3f}]")
                    if 'ideal' in non:
                        print(f"                Ideal (0) ±0.5 SD: [{non['ideal_minus_05sd']:7.3f}, {non['ideal_plus_05sd']:7.3f}]")
                        print(f"                Ideal (0) ±1.0 SD: [{non['ideal_minus_1sd']:7.3f}, {non['ideal_plus_1sd']:7.3f}] (n={non['count']:,})")
                    else:
                        print(f"                (n={non['count']:,})")
                    print(f"  Difference:    std_diff={abs(can['std'] - non['std']):7.3f}")
    
    # H-bond comparison
    if '_CANONICAL' in hbond_results and '_NONCANONICAL' in hbond_results:
        print("\n" + "-"*80)
        print("CANONICAL vs NON-CANONICAL COMPARISON (H-bonds)")
        print("-"*80)
        
        for param in ['distance', 'angle_min', 'quality']:
            if param in hbond_results['_CANONICAL'] and param in hbond_results['_NONCANONICAL']:
                can = hbond_results['_CANONICAL'][param]
                non = hbond_results['_NONCANONICAL'][param]
                
                if 'error' not in can and 'error' not in non:
                    print(f"\n{param.upper()}:")
                    print(f"  Canonical:    mean={can['mean']:7.3f}, std={can['std']:7.3f}")
                    print(f"                Mean ±0.5 SD: [{can['mean_minus_05sd']:7.3f}, {can['mean_plus_05sd']:7.3f}]")
                    print(f"                Mean ±1.0 SD: [{can['mean_minus_1sd']:7.3f}, {can['mean_plus_1sd']:7.3f}] (n={can['count']:,})")
                    print(f"  Non-canonical: mean={non['mean']:7.3f}, std={non['std']:7.3f}")
                    print(f"                Mean ±0.5 SD: [{non['mean_minus_05sd']:7.3f}, {non['mean_plus_05sd']:7.3f}]")
                    print(f"                Mean ±1.0 SD: [{non['mean_minus_1sd']:7.3f}, {non['mean_plus_1sd']:7.3f}] (n={non['count']:,})")
                    print(f"  Difference:    std_diff={abs(can['std'] - non['std']):7.3f}")


def load_allowed_pdb_ids(csv_path: Path, exclude_pdb_ids: set) -> set:
    """
    Load allowed PDB IDs from uniqueRNAs.csv and remove excluded ones (e.g., G-quads).

    Args:
        csv_path: Path to uniqueRNAs.csv (must have 'unique_pdb_id' column)
        exclude_pdb_ids: Set of PDB IDs (uppercase) to exclude

    Returns:
        Set of allowed PDB IDs (uppercase)
    """
    if not csv_path.exists():
        print(f"Error: Unique RNAs file not found: {csv_path}")
        sys.exit(1)

    df = pd.read_csv(csv_path)
    all_ids = {pdb_id.strip().upper() for pdb_id in df['unique_pdb_id'].dropna()}
    allowed = all_ids - exclude_pdb_ids
    print(f"Loaded {len(all_ids)} unique PDB IDs from {csv_path}")
    print(f"Excluded {len(all_ids) - len(allowed)} G-quad PDB IDs, {len(allowed)} remaining")
    return allowed


def main():
    """Main analysis function."""
    config = Config()
    basepair_dir = Path(config.BASEPAIR_DIR)
    hbond_dir = Path(config.HBOND_DIR)
    output_file = Path('edge_type_analysis.json')

    # Build allowlist: only PDB IDs in uniqueRNAs.csv, minus G-quads
    exclude_pdb_ids = {p.upper() for p in g_quads}
    unique_rnas_path = Path("data/uniqueRNAs.csv")
    allowed_pdb_ids = load_allowed_pdb_ids(unique_rnas_path, exclude_pdb_ids)

    if not basepair_dir.exists():
        print(f"Error: Base pair directory not found: {basepair_dir}")
        sys.exit(1)

    if not hbond_dir.exists():
        print(f"Error: H-bond directory not found: {hbond_dir}")
        sys.exit(1)

    print("="*80)
    print("EDGE TYPE ANALYSIS")
    print("="*80)
    print("\nThis script will:")
    print("  1. Load a sample of base pair data")
    print("  2. Load matching H-bond data")
    print("  3. Analyze parameters by edge type (lw notation)")
    print("  4. Compare canonical vs non-canonical pairs")
    print("\nUsing efficient sampling to avoid memory issues...")
    
    # Load base pairs (process all files for complete analysis)
    print("\nProcessing ALL base pair files (this may take a few minutes)...")
    basepairs = load_basepair_sample(basepair_dir, max_files=None, include_pdb_ids=allowed_pdb_ids)

    # Exclude adjacent base pairs (same chain, residue numbers differ by 1)
    before_adj = len(basepairs)
    basepairs = [bp for bp in basepairs if not is_adjacent_pair(bp.get('res_1', ''), bp.get('res_2', ''))]
    excluded_adj = before_adj - len(basepairs)
    if excluded_adj > 0:
        print(f"✓ Excluded {excluded_adj:,} adjacent base pairs (residue diff=1), {len(basepairs):,} remaining")

    if len(basepairs) == 0:
        print("Error: No base pair data loaded")
        sys.exit(1)
    
    # Create lookup for linking H-bonds to base pairs
    print("\nCreating base pair lookup...")
    bp_lookup = create_basepair_lookup(basepairs)
    print(f"✓ Created lookup for {len(bp_lookup):,} base pairs")
    
    # Load H-bonds with edge type information (process all files, no sampling limit)
    print("\nProcessing ALL H-bond files (this may take a few minutes)...")
    print("Note: Loading all H-bonds (no per-file limit) to get complete dataset")
    hbond_df = load_hbond_sample_with_lw(
        hbond_dir,
        bp_lookup,
        max_files=None,
        max_hbonds_per_file=None,
        include_pdb_ids=allowed_pdb_ids
    )
    
    # Analyze geometry by edge type
    geometry_results = analyze_geometry_by_edge_type(basepairs)
    
    # Analyze H-bonds by edge type
    hbond_results = {}
    if len(hbond_df) > 0:
        hbond_results = analyze_hbonds_by_edge_type(hbond_df)
    
    # Print summary
    print_summary(geometry_results, hbond_results)
    
    # Compile results
    results = {
        'summary': {
            'total_basepairs': len(basepairs),
            'total_hbonds': len(hbond_df),
            'analysis_date': pd.Timestamp.now().isoformat(),
            'note': 'Sampled data for efficiency'
        },
        'geometry_by_edge_type': geometry_results,
        'hbond_by_edge_type': hbond_results
    }
    
    # Save to JSON
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    print("\n" + "="*80)
    print(f"✓ Analysis complete! Results saved to: {output_file}")
    print("="*80)
    print("\nKey findings:")
    print("  - _CANONICAL and _NONCANONICAL groups show aggregate statistics")
    print("  - Compare std values to see if edge-specific thresholds are needed")
    print("  - If std differs significantly, consider edge-specific thresholds")
    print("="*80)


if __name__ == '__main__':
    main()

