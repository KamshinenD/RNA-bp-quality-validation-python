#!/usr/bin/env python3
"""
Pre-compute and cache full structure scores for all unique PDB IDs.
This avoids recomputing full structure scores for each motif.
"""

import json
import sys
from pathlib import Path
import glob
import subprocess

def find_unique_pdb_ids(motifs_dir='motifs'):
    """Find all unique PDB IDs from motif CIF files."""
    motif_files = glob.glob(f"{motifs_dir}/*.cif")
    pdb_ids = set()
    
    for motif_file in motif_files:
        # Extract PDB ID from motif filename (e.g., HAIRPIN-10-...-7KGB-1.cif -> 7KGB)
        parts = Path(motif_file).stem.split('-')
        for part in reversed(parts):
            if len(part) == 4 and part[0].isdigit() and part[1:].isalnum():
                pdb_ids.add(part)
                break
    
    return sorted(pdb_ids)

def cache_full_structure_score(pdb_id, cache_dir='full_structure_cache'):
    """Compute and cache full structure score for a PDB ID."""
    cache_file = Path(cache_dir) / f"{pdb_id}.json"
    
    # Skip if already cached
    if cache_file.exists():
        try:
            with open(cache_file, 'r') as f:
                data = json.load(f)
                if 'full_structure_score' in data:
                    return True
        except:
            pass
    
    # Run app.py to get full structure score
    try:
        result = subprocess.run(
            [sys.executable, 'app.py', '--pdb_id', pdb_id],
            capture_output=True,
            text=True,
            timeout=300  # 5 minute timeout
        )
        
        if result.returncode == 0:
            # Read the report.json
            report_file = Path('report.json')
            if report_file.exists():
                with open(report_file, 'r') as f:
                    data = json.load(f)
                    full_score = data.get('overall_score')
                    
                    if full_score is not None:
                        # Save to cache
                        Path(cache_dir).mkdir(exist_ok=True)
                        cache_data = {
                            'pdb_id': pdb_id,
                            'full_structure_score': full_score,
                            'full_structure_grade': data.get('grade', 'N/A'),
                            'total_base_pairs': data.get('total_base_pairs', 0),
                            'num_nucleotides': data.get('num_nucleotides', 0)
                        }
                        with open(cache_file, 'w') as f:
                            json.dump(cache_data, f, indent=2)
                        return True
    except subprocess.TimeoutExpired:
        print(f"  ⚠ Timeout for {pdb_id}")
    except Exception as e:
        print(f"  ✗ Error for {pdb_id}: {e}")
    
    return False

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Cache full structure scores")
    parser.add_argument('--motifs-dir', default='motifs', help='Motifs directory')
    parser.add_argument('--cache-dir', default='full_structure_cache', help='Cache directory')
    parser.add_argument('--pdb-id', help='Cache specific PDB ID only')
    
    args = parser.parse_args()
    
    if args.pdb_id:
        # Cache single PDB ID
        print(f"Caching full structure score for {args.pdb_id}...")
        success = cache_full_structure_score(args.pdb_id, args.cache_dir)
        if success:
            print(f"✓ Cached {args.pdb_id}")
        else:
            print(f"✗ Failed to cache {args.pdb_id}")
            sys.exit(1)
    else:
        # Find and cache all unique PDB IDs
        print("Finding unique PDB IDs from motifs...")
        pdb_ids = find_unique_pdb_ids(args.motifs_dir)
        print(f"Found {len(pdb_ids)} unique PDB IDs")
        
        Path(args.cache_dir).mkdir(exist_ok=True)
        
        cached = 0
        failed = 0
        
        for i, pdb_id in enumerate(pdb_ids, 1):
            print(f"[{i}/{len(pdb_ids)}] Processing {pdb_id}...", end=' ')
            if cache_full_structure_score(pdb_id, args.cache_dir):
                cached += 1
                print("✓")
            else:
                failed += 1
                print("✗")
        
        print(f"\n{'='*60}")
        print(f"Cached: {cached}/{len(pdb_ids)}")
        print(f"Failed: {failed}/{len(pdb_ids)}")
        print(f"{'='*60}")

if __name__ == '__main__':
    main()

