#!/usr/bin/env python3
"""
Cache all unique full RNA structure scores.
Can be run locally or on cluster.
Skips already cached PDB IDs (resumable).
"""

import sys
import json
import glob
from pathlib import Path
import subprocess
from collections import Counter

def find_unique_pdb_ids(motifs_dir='unique_motifs'):
    """Find all unique PDB IDs from motif CIF files."""
    motif_files = glob.glob(f"{motifs_dir}/*.cif")
    pdb_ids = set()
    
    print(f"Scanning {len(motif_files)} motif files...")
    
    for motif_file in motif_files:
        # Extract PDB ID from motif filename (e.g., HAIRPIN-10-...-7KGB-1.cif -> 7KGB)
        parts = Path(motif_file).stem.split('-')
        for part in reversed(parts):
            if len(part) == 4 and part[0].isdigit() and part[1:].isalnum():
                pdb_ids.add(part)
                break
    
    return sorted(pdb_ids)

def cache_full_structure_score(pdb_id, cache_dir='full_structure_cache', timeout=300):
    """Compute and cache full structure score for a PDB ID."""
    cache_file = Path(cache_dir) / f"{pdb_id}.json"
    
    # Skip if already cached
    if cache_file.exists():
        try:
            with open(cache_file, 'r') as f:
                data = json.load(f)
                if 'full_structure_score' in data and data['full_structure_score'] is not None:
                    return True, "already_cached"
        except:
            # Corrupted cache file, re-cache
            pass
    
    # Check if data files exist
    bp_file = Path(f'data/basepairs/{pdb_id}.json')
    hb_file = Path(f'data/hbonds/{pdb_id}.csv')
    
    if not bp_file.exists():
        return False, "no_basepair_file"
    
    if not hb_file.exists():
        return False, "no_hbond_file"
    
    # Run app.py to get full structure score
    try:
        result = subprocess.run(
            [sys.executable, 'app.py', '--pdb_id', pdb_id],
            capture_output=True,
            text=True,
            timeout=timeout
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
                            'total_base_pairs': data.get('total_base_pairs', 0),
                            'num_nucleotides': data.get('num_nucleotides', 0)
                        }
                        with open(cache_file, 'w') as f:
                            json.dump(cache_data, f, indent=2)
                        return True, "cached"
                    else:
                        return False, "no_score_in_report"
            else:
                return False, "no_report_file"
        else:
            # Check if it's because of no base pairs
            if "No base pairs found" in result.stdout or "No base pairs found" in result.stderr:
                return False, "no_base_pairs"
            return False, f"exit_code_{result.returncode}"
            
    except subprocess.TimeoutExpired:
        return False, "timeout"
    except Exception as e:
        return False, f"error_{str(e)[:50]}"

def main():
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Cache all unique full RNA structure scores from motifs"
    )
    parser.add_argument(
        '--motifs-dir',
        default='unique_motifs',
        help='Motifs directory (default: unique_motifs)'
    )
    parser.add_argument(
        '--cache-dir',
        default='full_structure_cache',
        help='Cache directory (default: full_structure_cache)'
    )
    parser.add_argument(
        '--timeout',
        type=int,
        default=300,
        help='Timeout per PDB ID in seconds (default: 300)'
    )
    parser.add_argument(
        '--skip-existing',
        action='store_true',
        default=True,
        help='Skip already cached PDB IDs (default: True)'
    )
    
    args = parser.parse_args()
    
    print("="*80)
    print("CACHING ALL UNIQUE RNA STRUCTURE SCORES")
    print("="*80)
    
    # Find unique PDB IDs
    print(f"\nFinding unique PDB IDs from {args.motifs_dir}...")
    pdb_ids = find_unique_pdb_ids(args.motifs_dir)
    print(f"Found {len(pdb_ids)} unique PDB IDs")
    
    # Check existing cache
    cache_dir = Path(args.cache_dir)
    cache_dir.mkdir(exist_ok=True)
    
    existing_cache = set()
    if args.skip_existing:
        for cache_file in cache_dir.glob("*.json"):
            try:
                with open(cache_file, 'r') as f:
                    data = json.load(f)
                    if 'full_structure_score' in data and data['full_structure_score'] is not None:
                        existing_cache.add(cache_file.stem)
            except:
                pass
        
        if existing_cache:
            print(f"Found {len(existing_cache)} already cached PDB IDs")
            pdb_ids_to_cache = [pdb for pdb in pdb_ids if pdb not in existing_cache]
            print(f"Need to cache: {len(pdb_ids_to_cache)} PDB IDs")
        else:
            pdb_ids_to_cache = pdb_ids
    else:
        pdb_ids_to_cache = pdb_ids
    
    if not pdb_ids_to_cache:
        print("\n✅ All PDB IDs already cached!")
        return 0
    
    # Cache all PDB IDs
    print(f"\n{'='*80}")
    print(f"Caching {len(pdb_ids_to_cache)} PDB IDs...")
    print(f"{'='*80}\n")
    
    results = Counter()
    failed_pdb_ids = []
    
    for i, pdb_id in enumerate(pdb_ids_to_cache, 1):
        print(f"[{i}/{len(pdb_ids_to_cache)}] {pdb_id}...", end=' ', flush=True)
        
        success, reason = cache_full_structure_score(pdb_id, args.cache_dir, args.timeout)
        
        if success:
            if reason == "already_cached":
                print("✓ (already cached)")
            else:
                print("✓")
            results[reason] += 1
        else:
            print(f"✗ ({reason})")
            results[reason] += 1
            failed_pdb_ids.append((pdb_id, reason))
        
        # Progress update every 10
        if i % 10 == 0:
            cached_so_far = results["cached"] + results["already_cached"]
            print(f"  Progress: {cached_so_far}/{len(pdb_ids_to_cache)} cached")
    
    # Summary
    print(f"\n{'='*80}")
    print("SUMMARY")
    print(f"{'='*80}")
    print(f"Total unique PDB IDs: {len(pdb_ids)}")
    print(f"Already cached: {results['already_cached']}")
    print(f"Newly cached: {results['cached']}")
    print(f"Failed: {sum(v for k, v in results.items() if k not in ['cached', 'already_cached'])}")
    
    if results:
        print(f"\nBreakdown:")
        for reason, count in sorted(results.items(), key=lambda x: -x[1]):
            print(f"  {reason}: {count}")
    
    if failed_pdb_ids:
        print(f"\n{'='*80}")
        print(f"FAILED PDB IDs ({len(failed_pdb_ids)}):")
        print(f"{'='*80}")
        
        # Group by reason
        by_reason = {}
        for pdb_id, reason in failed_pdb_ids:
            if reason not in by_reason:
                by_reason[reason] = []
            by_reason[reason].append(pdb_id)
        
        for reason, pdb_list in sorted(by_reason.items()):
            print(f"\n{reason} ({len(pdb_list)}):")
            for pdb_id in pdb_list[:10]:
                print(f"  {pdb_id}")
            if len(pdb_list) > 10:
                print(f"  ... and {len(pdb_list) - 10} more")
        
        # Save failed list
        failed_file = Path('failed_pdb_ids.txt')
        with open(failed_file, 'w') as f:
            for pdb_id, _ in failed_pdb_ids:
                f.write(f"{pdb_id}\n")
        print(f"\nFailed PDB IDs saved to: {failed_file}")
    
    total_cached = results['cached'] + results['already_cached'] + len(existing_cache)
    print(f"\n{'='*80}")
    print(f"✅ Total cached: {total_cached}/{len(pdb_ids)}")
    print(f"{'='*80}")
    
    return 0 if total_cached == len(pdb_ids) else 1

if __name__ == '__main__':
    sys.exit(main())

