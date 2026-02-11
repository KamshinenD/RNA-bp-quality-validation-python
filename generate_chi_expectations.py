#!/usr/bin/env python3
"""
Generate expected chi glycosidic bond conformations per (edge_type, bp_type).

Analyzes all torsion + basepair data to determine whether each (edge, bp_type)
combination expects anti, syn, or both conformations. Output is saved to
data/chi_expectations.json and consumed by scorer2.py at runtime.

Usage:
    python generate_chi_expectations.py
"""

import json
from collections import defaultdict
from pathlib import Path


# Chi conformation boundaries (same as torsion_scores_analysis.py)
def classify_chi(chi: float) -> str:
    """Classify chi angle as anti, syn, or intermediate."""
    if chi is None:
        return "unknown"
    if -180 <= chi <= -60 or 120 <= chi <= 180:
        return "anti"
    elif 0 <= chi <= 90:
        return "syn"
    else:
        return "intermediate"


# Canonical base mapping (modified bases -> parent)
MODIFIED_BASE_MAP = {
    'DA': 'A', 'DC': 'C', 'DG': 'G', 'DT': 'U', 'DU': 'U',
    'PSU': 'U', 'H2U': 'U', '5MU': 'U', '5MC': 'C', '1MA': 'A',
    '2MG': 'G', '7MG': 'G', 'M2G': 'G', 'OMG': 'G', 'YYG': 'G',
    'OMC': 'C', 'OMU': 'U', '4SU': 'U', 'I': 'G',
}


def map_to_canonical(bp_type: str) -> str:
    """Map a bp_type like 'DC-DG' to canonical 'C-G'."""
    parts = bp_type.split('-')
    if len(parts) != 2:
        return bp_type
    mapped = []
    for p in parts:
        if p in MODIFIED_BASE_MAP:
            mapped.append(MODIFIED_BASE_MAP[p])
        elif p in ('A', 'C', 'G', 'U'):
            mapped.append(p)
        else:
            mapped.append(MODIFIED_BASE_MAP.get(p, p))
    return '-'.join(mapped)


def main():
    torsion_dir = Path('data/torsions')
    basepair_dir = Path('data/basepairs')
    output_file = Path('data/chi_expectations.json')

    if not torsion_dir.exists():
        print(f"Error: {torsion_dir} not found")
        return
    if not basepair_dir.exists():
        print(f"Error: {basepair_dir} not found")
        return

    # Collect chi conformations per (edge_type, bp_type, residue_position)
    # Structure: counts[edge_type][bp_type] = {'anti': N, 'syn': N, 'intermediate': N, 'count': N}
    counts = defaultdict(lambda: defaultdict(lambda: {'anti': 0, 'syn': 0, 'intermediate': 0, 'count': 0}))

    torsion_files = sorted(torsion_dir.glob('*.json'))
    print(f"Processing {len(torsion_files)} torsion files...")

    processed = 0
    for tf in torsion_files:
        pdb_id = tf.stem
        bp_file = basepair_dir / f'{pdb_id}.json'
        if not bp_file.exists():
            continue

        try:
            with open(tf) as f:
                torsions = json.load(f)
            with open(bp_file) as f:
                basepairs = json.load(f)
        except (json.JSONDecodeError, IOError):
            continue

        for bp in basepairs:
            edge_type = bp.get('lw', '') or '_OTHER'
            bp_type_raw = bp.get('bp_type', '')
            bp_type = map_to_canonical(bp_type_raw)

            # Only consider canonical base pair types
            parts = bp_type.split('-')
            if len(parts) != 2 or any(p not in ('A', 'C', 'G', 'U') for p in parts):
                continue

            res_1 = bp.get('res_1', '')
            res_2 = bp.get('res_2', '')

            # Check chi for both residues
            for res_id in [res_1, res_2]:
                t = torsions.get(res_id, {})
                chi = t.get('chi')
                if chi is None:
                    continue

                conf = classify_chi(chi)
                if conf == 'unknown':
                    continue

                bucket = counts[edge_type][bp_type]
                bucket[conf] += 1
                bucket['count'] += 1

        processed += 1
        if processed % 500 == 0:
            print(f"  Processed {processed}/{len(torsion_files)} structures...")

    print(f"Processed {processed} structures total.")

    # Build expectations: for each (edge_type, bp_type) with >= 100 observations
    MIN_OBSERVATIONS = 100
    expectations = {}

    # Also compute global stats for _OTHER fallback
    global_counts = {'anti': 0, 'syn': 0, 'intermediate': 0, 'count': 0}

    for edge_type, bp_types in sorted(counts.items()):
        edge_dict = {}
        for bp_type, bucket in sorted(bp_types.items()):
            total = bucket['count']
            if total < MIN_OBSERVATIONS:
                continue

            anti_pct = bucket['anti'] / total * 100
            syn_pct = bucket['syn'] / total * 100
            inter_pct = bucket['intermediate'] / total * 100

            # Determine expected conformation
            if anti_pct >= 95:
                expected = 'anti'
            elif syn_pct >= 95:
                expected = 'syn'
            elif anti_pct >= 5 and syn_pct >= 5:
                expected = 'both'
            elif anti_pct >= syn_pct:
                expected = 'anti'
            else:
                expected = 'syn'

            edge_dict[bp_type] = {
                'expected': expected,
                'anti_pct': round(anti_pct, 1),
                'syn_pct': round(syn_pct, 1),
                'intermediate_pct': round(inter_pct, 1),
                'count': total,
            }

            # Accumulate global stats
            for k in ['anti', 'syn', 'intermediate', 'count']:
                global_counts[k] += bucket[k]

        if edge_dict:
            expectations[edge_type] = edge_dict

    # Add global fallback _OTHER
    if global_counts['count'] > 0:
        total = global_counts['count']
        anti_pct = global_counts['anti'] / total * 100
        syn_pct = global_counts['syn'] / total * 100
        expectations['_OTHER'] = {
            'expected': 'anti',  # anti is overwhelmingly dominant globally
            'anti_pct': round(anti_pct, 1),
            'syn_pct': round(syn_pct, 1),
            'count': total,
        }

    # Save
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, 'w') as f:
        json.dump(expectations, f, indent=2)

    print(f"\nSaved chi expectations to {output_file}")
    print(f"Edge types: {len(expectations) - 1} + _OTHER fallback")

    # Print summary
    for edge_type in sorted(expectations.keys()):
        if edge_type == '_OTHER':
            continue
        bp_types = expectations[edge_type]
        print(f"\n  {edge_type}: {len(bp_types)} bp types")
        for bp_type, info in sorted(bp_types.items()):
            print(f"    {bp_type}: expected={info['expected']} "
                  f"(anti={info['anti_pct']}%, syn={info['syn_pct']}%, n={info['count']})")


if __name__ == '__main__':
    main()
