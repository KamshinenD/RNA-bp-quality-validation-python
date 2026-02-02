#!/usr/bin/env python3
"""
Analyze distributions from RNA structure data to recommend threshold values.

This script analyzes base pair geometry and H-bond parameters from all available
data files and outputs statistical distributions and recommended thresholds.
"""

import json
import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Optional
import sys


def load_all_basepairs(basepair_dir: Path) -> List[Dict]:
    """Load all base pair JSON files and collect geometry parameters."""
    all_basepairs = []
    basepair_files = list(basepair_dir.glob("*.json"))
    
    print(f"Loading {len(basepair_files)} base pair files...")
    
    for bp_file in basepair_files:
        try:
            with open(bp_file, 'r') as f:
                data = json.load(f)
            
            # Handle both list and dict formats
            if isinstance(data, list):
                base_pairs = data
            elif isinstance(data, dict):
                base_pairs = data.get('base_pairs', [])
            else:
                continue
            
            all_basepairs.extend(base_pairs)
            
        except Exception as e:
            print(f"Warning: Error loading {bp_file.name}: {e}")
            continue
    
    print(f"✓ Loaded {len(all_basepairs)} total base pairs")
    return all_basepairs


# H-bond CSV analysis excluded as requested
# def load_all_hbonds(hbond_dir: Path) -> pd.DataFrame:
#     """Load all H-bond CSV files and combine into single DataFrame."""
#     ...


def calculate_statistics(values: List[float], param_name: str) -> Dict:
    """Calculate comprehensive statistics for a parameter."""
    if not values or len(values) == 0:
        return {"error": "No data available"}
    
    values_array = np.array(values)
    values_array = values_array[~np.isnan(values_array)]  # Remove NaN
    
    if len(values_array) == 0:
        return {"error": "No valid data after removing NaN"}
    
    stats = {
        "count": len(values_array),
        "mean": float(np.mean(values_array)),
        "std": float(np.std(values_array)),
        "min": float(np.min(values_array)),
        "max": float(np.max(values_array)),
        "median": float(np.median(values_array)),
        "percentiles": {
            "p5": float(np.percentile(values_array, 5)),
            "p10": float(np.percentile(values_array, 10)),
            "p25": float(np.percentile(values_array, 25)),
            "p50": float(np.percentile(values_array, 50)),
            "p75": float(np.percentile(values_array, 75)),
            "p90": float(np.percentile(values_array, 90)),
            "p95": float(np.percentile(values_array, 95)),
            "p99": float(np.percentile(values_array, 99))
        }
    }
    
    # Recommended thresholds based on percentiles
    # For ranges (min/max), use 5th and 95th percentiles
    # For single thresholds, use mean ± 2 std or 95th percentile
    stats["recommended_thresholds"] = {
        "p5_p95_range": {
            "min": stats["percentiles"]["p5"],
            "max": stats["percentiles"]["p95"]
        },
        "p10_p90_range": {
            "min": stats["percentiles"]["p10"],
            "max": stats["percentiles"]["p90"]
        },
        "mean_plus_minus_2std": {
            "min": stats["mean"] - 2 * stats["std"],
            "max": stats["mean"] + 2 * stats["std"]
        },
        "mean_plus_minus_1std": {
            "min": stats["mean"] - stats["std"],
            "max": stats["mean"] + stats["std"]
        }
    }
    
    return stats


def analyze_geometry_parameters(basepairs: List[Dict]) -> Dict:
    """Analyze geometry parameters from base pairs."""
    print("\nAnalyzing geometry parameters...")
    
    geometry_params = {
        'shear': [],
        'stretch': [],
        'stagger': [],
        'buckle': [],
        'propeller': [],
        'opening': []
    }
    
    for bp in basepairs:
        for param in geometry_params.keys():
            if param in bp and bp[param] is not None:
                try:
                    geometry_params[param].append(float(bp[param]))
                except (ValueError, TypeError):
                    pass
    
    results = {}
    for param, values in geometry_params.items():
        print(f"  {param}: {len(values)} values")
        results[param] = calculate_statistics(values, param)
    
    return results


# H-bond CSV analysis excluded as requested
# def analyze_hbond_parameters(hbonds_df: pd.DataFrame) -> Dict:
#     """Analyze H-bond parameters from H-bond data."""
#     ...


def analyze_dssr_hbond_score(basepairs: List[Dict]) -> Dict:
    """Analyze DSSR H-bond score from base pairs."""
    print("\nAnalyzing DSSR H-bond scores...")
    
    hbond_scores = []
    for bp in basepairs:
        if 'hbond_score' in bp and bp['hbond_score'] is not None:
            try:
                score = float(bp['hbond_score'])
                if score > 0:  # Only include pairs with H-bonds
                    hbond_scores.append(score)
            except (ValueError, TypeError):
                pass
    
    print(f"  hbond_score: {len(hbond_scores)} values (score > 0)")
    return calculate_statistics(hbond_scores, 'hbond_score')


def generate_recommendations(geometry_stats: Dict, dssr_stats: Optional[Dict] = None) -> Dict:
    """Generate recommended threshold values based on statistics.
    
    For parameters that should be around 0 (shear, stretch, stagger, buckle, propeller, opening),
    we use standard deviation-based thresholds (mean ± N*std) which is more scientifically
    appropriate than percentiles for normally distributed parameters centered at 0.
    """
    recommendations = {}
    
    print("\nGenerating threshold recommendations...")
    print("  Using standard deviation-based thresholds for parameters centered around 0")
    
    # Geometry thresholds - using standard deviation from 0 (ideal value)
    recommendations['geometry'] = {}
    
    # Standard deviation multipliers for different strictness levels
    # 1 std = ~68% of data, 2 std = ~95%, 3 std = ~99.7%
    std_multiplier = 1.0  # Use 1 std for ~68% coverage (default recommendation)
    
    for param, stats in geometry_stats.items():
        if 'error' not in stats:
            mean = stats['mean']  # Keep for reference, but don't use in calculation
            std = stats['std']
            
            # All geometry parameters should be around 0 (ideal value)
            # So we measure deviation from 0, not from the observed mean
            if param in ['shear', 'stagger']:
                # These are typically absolute values (max thresholds)
                # Use N*std as max deviation from 0
                max_value = std_multiplier * std
                recommendations['geometry'][f'{param.upper()}_MAX'] = {
                    'value': max_value,
                    'rationale': f'{std_multiplier}*std ({std:.3f}) = {max_value:.3f} (deviation from ideal 0)',
                    'ideal_value': 0.0,
                    'observed_mean': mean,
                    'std': std,
                    'std_multiplier': std_multiplier,
                    'current_value': None
                }
            elif param == 'stretch':
                # Range parameter - use 0 ± N*std
                min_value = -std_multiplier * std
                max_value = std_multiplier * std
                recommendations['geometry']['STRETCH_MIN'] = {
                    'value': min_value,
                    'rationale': f'-{std_multiplier}*std ({std:.3f}) = {min_value:.3f} (deviation from ideal 0)',
                    'ideal_value': 0.0,
                    'observed_mean': mean,
                    'std': std,
                    'std_multiplier': std_multiplier,
                    'current_value': None
                }
                recommendations['geometry']['STRETCH_MAX'] = {
                    'value': max_value,
                    'rationale': f'+{std_multiplier}*std ({std:.3f}) = {max_value:.3f} (deviation from ideal 0)',
                    'ideal_value': 0.0,
                    'observed_mean': mean,
                    'std': std,
                    'std_multiplier': std_multiplier,
                    'current_value': None
                }
            elif param in ['buckle', 'opening']:
                # Absolute max thresholds - use N*std as max deviation from 0
                max_value = std_multiplier * std
                recommendations['geometry'][f'{param.upper()}_MAX'] = {
                    'value': max_value,
                    'rationale': f'{std_multiplier}*std ({std:.3f}) = {max_value:.3f} (absolute deviation from ideal 0)',
                    'ideal_value': 0.0,
                    'observed_mean': mean,
                    'std': std,
                    'std_multiplier': std_multiplier,
                    'current_value': None
                }
            elif param == 'propeller':
                # Range parameter (can be negative) - use 0 ± N*std
                min_value = -std_multiplier * std
                max_value = std_multiplier * std
                recommendations['geometry']['PROPELLER_MIN'] = {
                    'value': min_value,
                    'rationale': f'-{std_multiplier}*std ({std:.3f}) = {min_value:.3f} (deviation from ideal 0)',
                    'ideal_value': 0.0,
                    'observed_mean': mean,
                    'std': std,
                    'std_multiplier': std_multiplier,
                    'current_value': None
                }
                recommendations['geometry']['PROPELLER_MAX'] = {
                    'value': max_value,
                    'rationale': f'+{std_multiplier}*std ({std:.3f}) = {max_value:.3f} (deviation from ideal 0)',
                    'ideal_value': 0.0,
                    'observed_mean': mean,
                    'std': std,
                    'std_multiplier': std_multiplier,
                    'current_value': None
                }
    
    # Add alternative recommendations with different std multipliers for comparison
    recommendations['geometry_alternatives'] = {
        'note': 'Alternative thresholds using different standard deviation multipliers (deviation from ideal 0)',
        '1_std': {},  # ~68% of data
        '2_std': {},  # ~95% of data (default recommendation)
        '3_std': {}   # ~99.7% of data
    }
    
    for multiplier, label in [(1.0, '1_std'), (2.0, '2_std'), (3.0, '3_std')]:
        for param, stats in geometry_stats.items():
            if 'error' not in stats:
                std = stats['std']
                
                if param in ['shear', 'stagger']:
                    # Max deviation from 0
                    max_value = multiplier * std
                    recommendations['geometry_alternatives'][label][f'{param.upper()}_MAX'] = max_value
                elif param == 'stretch':
                    # Range: 0 ± N*std
                    recommendations['geometry_alternatives'][label]['STRETCH_MIN'] = -multiplier * std
                    recommendations['geometry_alternatives'][label]['STRETCH_MAX'] = multiplier * std
                elif param in ['buckle', 'opening']:
                    # Absolute max deviation from 0
                    recommendations['geometry_alternatives'][label][f'{param.upper()}_MAX'] = multiplier * std
                elif param == 'propeller':
                    # Range: 0 ± N*std
                    recommendations['geometry_alternatives'][label]['PROPELLER_MIN'] = -multiplier * std
                    recommendations['geometry_alternatives'][label]['PROPELLER_MAX'] = multiplier * std
    
    # DSSR H-bond score threshold (from base pair data, not H-bond CSV)
    # For H-bond score, we still use percentile since it's not centered around 0
    if dssr_stats and 'error' not in dssr_stats:
        recommendations['dssr'] = {}
        recommendations['dssr']['HBOND_SCORE_MIN'] = {
            'value': dssr_stats['percentiles']['p5'],
            'rationale': '5th percentile from DSSR H-bond score distribution (pairs with H-bonds)',
            'mean': dssr_stats['mean'],
            'std': dssr_stats['std'],
            'current_value': None
        }
        # Also provide std-based alternative
        mean = dssr_stats['mean']
        std = dssr_stats['std']
        recommendations['dssr']['HBOND_SCORE_MIN_std_based'] = {
            'value': mean - 2.0 * std,  # Mean - 2*std
            'rationale': f'Mean ({mean:.3f}) - 2*std ({std:.3f}) = {mean - 2.0*std:.3f}',
            'mean': mean,
            'std': std,
            'std_multiplier': 2.0,
            'current_value': None
        }
    
    return recommendations


def main():
    """Main analysis function."""
    # Set up paths
    basepair_dir = Path('./data/basepairs')
    output_file = Path('threshold_analysis.json')
    
    if not basepair_dir.exists():
        print(f"Error: Base pair directory not found: {basepair_dir}")
        sys.exit(1)
    
    print("=" * 60)
    print("RNA Structure Threshold Analysis")
    print("(H-bond CSV analysis excluded)")
    print("=" * 60)
    
    # Load data
    basepairs = load_all_basepairs(basepair_dir)
    
    if len(basepairs) == 0:
        print("Error: No base pair data loaded")
        sys.exit(1)
    
    # Analyze parameters (geometry and DSSR score only)
    geometry_stats = analyze_geometry_parameters(basepairs)
    dssr_stats = analyze_dssr_hbond_score(basepairs)
    
    # Generate recommendations
    recommendations = generate_recommendations(geometry_stats, dssr_stats)
    
    # Compile results
    results = {
        'summary': {
            'total_basepairs': len(basepairs),
            'analysis_date': pd.Timestamp.now().isoformat(),
            'note': 'H-bond CSV analysis excluded as requested'
        },
        'geometry_statistics': geometry_stats,
        'dssr_hbond_score_statistics': dssr_stats,
        'recommended_thresholds': recommendations
    }
    
    # Save to JSON
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"\n{'=' * 60}")
    print(f"✓ Analysis complete! Results saved to: {output_file}")
    print(f"{'=' * 60}")
    
    # Print summary of recommendations
    print("\nRecommended Thresholds Summary:")
    print("-" * 60)
    for category, params in recommendations.items():
        print(f"\n{category.upper()}:")
        for param_name, param_data in params.items():
            if 'value' in param_data:
                print(f"  {param_name}: {param_data['value']:.3f}")
                print(f"    ({param_data['rationale']})")
    
    print(f"\n{'=' * 60}")
    print("For detailed statistics and distributions, see:", output_file)
    print(f"{'=' * 60}")


if __name__ == '__main__':
    main()

