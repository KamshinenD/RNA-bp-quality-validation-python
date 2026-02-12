#!/usr/bin/env python3
"""
Generate expected chi glycosidic bond conformations per (edge_type, bp_type).

Analyzes all torsion + basepair data to determine whether each (edge, bp_type)
combination expects anti, syn, or both conformations. Output is saved to
data/chi_expectations.json and consumed by scorer2.py at runtime.

Uncommon edge types (containing '.' or '--' or empty) are folded into _OTHER.

Usage:
    python generate_chi_expectations.py
"""

import csv
import json
from collections import defaultdict
from pathlib import Path

from g_quads import g_quads


# Chi conformation boundaries (must match scorer2._classify_chi)
# Chi = base/sugar orientation (O4′-C1′-N1-C2 pyrimidines, O4′-C1′-N9-C4 purines)
# Ranges match X3DNA-DSSR: anti = [-180,-90]∪[90,180], syn = (-90,90)
# https://x3dna.org/highlights/the-chi-x-torsion-angle-characterizes-base-sugar-relative-orientation
def classify_chi(chi: float) -> str:
    """Classify chi angle as anti or syn."""
    if chi is None:
        return "unknown"
    if chi <= -90 or chi >= 90:
        return "anti"
    else:
        return "syn"


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


def is_main_edge(edge_type: str) -> bool:
    """Return True for well-defined 3-char Leontis-Westhof edges (e.g. cWW, tSH)."""
    if not edge_type or edge_type == '--':
        return False
    if '.' in edge_type:
        return False
    return True


def format_json_compact(expectations: dict) -> str:
    """Format JSON so each bp_type entry is on a single line."""
    lines = ['{']
    edge_keys = list(expectations.keys())
    for ei, edge_type in enumerate(edge_keys):
        val = expectations[edge_type]
        edge_comma = ',' if ei < len(edge_keys) - 1 else ''

        if isinstance(val, dict) and 'expected' in val:
            # Top-level _OTHER is itself a single entry
            lines.append(f'  "{edge_type}": {json.dumps(val)}{edge_comma}')
        else:
            # Edge type with bp_type children
            lines.append(f'  "{edge_type}": {{')
            bp_keys = list(val.keys())
            for bi, bp_type in enumerate(bp_keys):
                bp_comma = ',' if bi < len(bp_keys) - 1 else ''
                lines.append(f'    "{bp_type}": {json.dumps(val[bp_type])}{bp_comma}')
            lines.append(f'  }}{edge_comma}')
    lines.append('}')
    return '\n'.join(lines)


def main():
    torsion_dir = Path('data/torsions')
    basepair_dir = Path('data/basepairs')
    unique_csv = Path('data/uniqueRNAS.csv')
    output_file = Path('data/chi_expectations.json')

    for p in [torsion_dir, basepair_dir, unique_csv]:
        if not p.exists():
            print(f"Error: {p} not found")
            return

    # Load unique PDB IDs and exclude g-quadruplexes
    g_quad_set = {p.upper() for p in g_quads}
    unique_pdb_ids = set()
    with open(unique_csv) as f:
        for row in csv.DictReader(f):
            pdb_id = row['unique_pdb_id'].strip().upper()
            if pdb_id not in g_quad_set:
                unique_pdb_ids.add(pdb_id)
    print(f"Using {len(unique_pdb_ids)} unique PDB IDs (excluded {len(g_quad_set)} g-quadruplexes)")

    # Collect chi conformations per (edge_type, bp_type)
    # Structure: counts[edge_type][bp_type] = {'anti': N, 'syn': N, 'intermediate': N, 'count': N}
    counts = defaultdict(lambda: defaultdict(lambda: {'anti': 0, 'syn': 0, 'intermediate': 0, 'count': 0}))

    torsion_files = sorted(tf for tf in torsion_dir.glob('*.json') if tf.stem.upper() in unique_pdb_ids)
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
            raw_edge = bp.get('lw', '') or ''

            # Route uncommon edges into _OTHER
            edge_type = raw_edge if is_main_edge(raw_edge) else '_OTHER'

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

    # Build expectations: for each (edge_type, bp_type) with >= 1000 observations
    MIN_OBSERVATIONS = 1000
    expectations = {}

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

        if edge_dict:
            expectations[edge_type] = edge_dict

    # Save with compact formatting (one bp per line)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, 'w') as f:
        f.write(format_json_compact(expectations))

    print(f"\nSaved chi expectations to {output_file}")
    edge_count = len([k for k in expectations if k != '_OTHER'])
    print(f"Edge types: {edge_count} named + _OTHER fallback")

    # Print summary
    for edge_type in sorted(expectations.keys()):
        bp_types = expectations[edge_type]
        if isinstance(bp_types, dict) and 'expected' not in bp_types:
            print(f"\n  {edge_type}: {len(bp_types)} bp types")
            for bp_type, info in sorted(bp_types.items()):
                print(f"    {bp_type}: expected={info['expected']} "
                      f"(anti={info['anti_pct']}%, syn={info['syn_pct']}%, n={info['count']})")
        else:
            print(f"\n  {edge_type}: expected={bp_types.get('expected','?')} (n={bp_types.get('count','?')})")


if __name__ == '__main__':
    main()
