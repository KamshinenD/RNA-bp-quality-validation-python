#!/usr/bin/env python3
"""
FAST batch processing - skips network calls for speed.

This version:
- Skips nucleotide count downloads
- Skips validation metrics API calls
- Only does the essential scoring
- 10-30x faster than regular app.py
"""

import os
import sys
import subprocess
import time
import csv
import json
from pathlib import Path
from typing import Set, List, Optional, Dict

# Import the scoring components directly
sys.path.insert(0, str(Path(__file__).parent))
from config import Config
from utils.data_loader import DataLoader
from scorer2 import Scorer
from utils.report_generator import ReportGenerator


def get_processed_pdb_ids(csv_file: str) -> Set[str]:
    """Get set of PDB IDs already in scores_summary.csv."""
    processed = set()
    if os.path.exists(csv_file):
        try:
            with open(csv_file, 'r') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    pdb_id = row.get('PDB_ID', '').strip().upper()
                    if pdb_id:
                        processed.add(pdb_id)
        except Exception as e:
            print(f"Warning: Could not read {csv_file}: {e}")
    return processed


def get_all_pdb_ids(basepairs_dir: str) -> List[str]:
    """Extract all PDB IDs from basepairs JSON files."""
    pdb_ids = []
    basepairs_path = Path(basepairs_dir)
    
    if not basepairs_path.exists():
        print(f"Error: Directory {basepairs_dir} does not exist")
        return []
    
    json_files = list(basepairs_path.glob('*.json'))
    for json_file in json_files:
        pdb_id = json_file.stem.upper()
        if pdb_id:
            pdb_ids.append(pdb_id)
    
    return sorted(set(pdb_ids))


def load_metadata_cache(cache_file: str = 'metadata_cache.json') -> dict:
    """Load cached metadata."""
    if Path(cache_file).exists():
        try:
            with open(cache_file, 'r') as f:
                return json.load(f)
        except Exception:
            return {}
    return {}


def get_nucleotide_count(pdb_id: str, cache: dict, data_loader) -> int:
    """Get nucleotide count from cache, or download if not cached."""
    nuc_counts = cache.get('nucleotide_counts', {})
    if pdb_id in nuc_counts:
        return nuc_counts[pdb_id]
    
    # Not in cache - download it
    print(f"    (Downloading nucleotide count for {pdb_id}...)")
    try:
        count = data_loader.get_nucleotide_count(pdb_id)
        # Update cache
        cache.setdefault('nucleotide_counts', {})[pdb_id] = count
        return count
    except Exception:
        return 0


def get_validation_metrics(pdb_id: str, cache: dict) -> Optional[dict]:
    """Get validation metrics from cache, or download if not cached."""
    val_metrics = cache.get('validation_metrics', {})
    if pdb_id in val_metrics and val_metrics[pdb_id] is not None:
        return val_metrics[pdb_id]
    
    # Not in cache - download it
    print(f"    (Downloading validation metrics for {pdb_id}...)")
    try:
        metrics = data_loader.get_validation_metrics(pdb_id)
        # Update cache
        cache.setdefault('validation_metrics', {})[pdb_id] = metrics
        return metrics
    except Exception:
        return None


def process_single_rna_fast(pdb_id: str, config, data_loader, scorer, report_gen, cache: dict) -> bool:
    """
    Process a single RNA structure using cached metadata when possible.
    
    Returns:
        True if successful, False otherwise
    """
    try:
        # Load data (local files only - fast)
        basepair_data = data_loader.load_basepairs(pdb_id)
        hbond_data = data_loader.load_hbonds(pdb_id)
        torsion_data = data_loader.load_torsions(pdb_id, quiet=True)

        if basepair_data is None or hbond_data is None:
            print(f"  ✗ Could not load data for {pdb_id}")
            return False

        if len(basepair_data) == 0:
            print(f"  ✗ No base pairs found for {pdb_id}")
            return False

        # Get nucleotide count from cache (or download if not cached)
        num_nucleotides = get_nucleotide_count(pdb_id, cache, data_loader)

        # Get validation metrics from cache (or download if not cached)
        validation_metrics = get_validation_metrics(pdb_id, cache)

        # Score the structure
        result = scorer.score_structure(basepair_data, hbond_data, torsion_data=torsion_data)
        
        # Convert result to dictionary for CSV export
        result_dict = scorer.export_to_dict(result)
        result_dict['pdb_id'] = pdb_id
        result_dict['analysis_type'] = 'baseline'
        result_dict['num_nucleotides'] = num_nucleotides
        
        # Save/update CSV summary (with validation metrics)
        report_gen.save_score_summary_csv(
            result_dict, 
            'scores_summary.csv',
            hbond_data=hbond_data,
            validation_metrics=validation_metrics
        )
        
        return True
        
    except Exception as e:
        print(f"  ✗ Error processing {pdb_id}: {e}")
        return False


def main():
    """Main execution function."""
    # Configuration
    BASEPAIRS_DIR = 'data/basepairs'
    SCORES_CSV = 'scores_summary.csv'
    DELAY_BETWEEN_RUNS = 0.5  # Reduced delay since no network calls
    
    print("="*80)
    print("FAST BATCH RNA SCORING PROCESSOR")
    print("(Uses cached metadata when available)")
    print("="*80)
    
    # Load metadata cache
    cache = load_metadata_cache()
    cache_file = 'metadata_cache.json'
    
    if cache:
        print(f"\nLoaded metadata cache:")
        print(f"  Cached nucleotide counts: {len(cache.get('nucleotide_counts', {}))}")
        print(f"  Cached validation metrics: {len(cache.get('validation_metrics', {}))}")
    else:
        print(f"\nNo metadata cache found. Will download as needed.")
        print(f"  (Run cache_metadata.py first to pre-download all metadata)")
    
    # Initialize components
    config = Config()
    data_loader = DataLoader(config)
    scorer = Scorer(config)
    report_gen = ReportGenerator(config)
    
    # Get all PDB IDs
    print(f"\nScanning {BASEPAIRS_DIR} for RNA structures...")
    all_pdb_ids = get_all_pdb_ids(BASEPAIRS_DIR)
    
    if not all_pdb_ids:
        print("No PDB IDs found!")
        return
    
    print(f"Found {len(all_pdb_ids)} unique PDB IDs")
    
    # Get already processed PDB IDs
    print(f"\nChecking {SCORES_CSV} for already processed structures...")
    processed = get_processed_pdb_ids(SCORES_CSV)
    print(f"Found {len(processed)} already processed structures")
    
    # Filter out already processed
    to_process = [pdb_id for pdb_id in all_pdb_ids if pdb_id not in processed]
    
    if not to_process:
        print("\n✓ All structures have already been processed!")
        return
    
    print(f"\n{len(to_process)} structures remaining to process")
    
    # Estimate time
    # If cached: ~1-2 seconds per structure
    # If not cached: ~10-30 seconds per structure (network calls)
    cached_count = sum(1 for pdb_id in to_process 
                      if pdb_id in cache.get('nucleotide_counts', {}) 
                      and pdb_id in cache.get('validation_metrics', {}))
    uncached_count = len(to_process) - cached_count
    
    if cached_count == len(to_process):
        estimated_seconds = len(to_process) * 1.5
        print(f"Estimated time: ~{estimated_seconds/3600:.1f} hours")
        print(f"(All metadata cached - ~1.5 seconds per structure)")
    else:
        cached_time = cached_count * 1.5
        uncached_time = uncached_count * 20  # ~20 seconds for network calls
        estimated_seconds = cached_time + uncached_time
        print(f"Estimated time: ~{estimated_seconds/3600:.1f} hours")
        print(f"  ({cached_count} cached: ~1.5s each, {uncached_count} need download: ~20s each)")
        print(f"  (Run cache_metadata.py first to speed this up!)")
    
    # Ask for confirmation
    response = input(f"\nProceed with FAST processing of {len(to_process)} structures? (yes/no): ")
    if response.lower() not in ['yes', 'y']:
        print("Cancelled.")
        return
    
    # Process each PDB ID
    successful = 0
    failed = 0
    failed_pdb_ids = []
    
    start_time = time.time()
    
    for i, pdb_id in enumerate(to_process, 1):
        print(f"\n[{i}/{len(to_process)}] {pdb_id}...", end=' ', flush=True)
        
        if process_single_rna_fast(pdb_id, config, data_loader, scorer, report_gen, cache):
            successful += 1
            print("✓")
        else:
            failed += 1
            failed_pdb_ids.append(pdb_id)
            print("✗")
        
        # Small delay to prevent overheating (much shorter since no network)
        if i < len(to_process) and DELAY_BETWEEN_RUNS > 0:
            time.sleep(DELAY_BETWEEN_RUNS)
        
        # Progress update every 100 structures
        if i % 100 == 0:
            elapsed = time.time() - start_time
            rate = i / elapsed
            remaining = len(to_process) - i
            eta_seconds = remaining / rate
            eta_minutes = eta_seconds / 60
            print(f"\n  Progress: {i}/{len(to_process)} | "
                  f"Success: {successful} | Failed: {failed} | "
                  f"ETA: {eta_minutes:.1f} minutes")
    
    # Save updated cache
    if cache:
        with open(cache_file, 'w') as f:
            json.dump(cache, f, indent=2)
        print(f"\nUpdated metadata cache saved to: {cache_file}")
    
    # Summary
    elapsed = time.time() - start_time
    print(f"\n{'='*80}")
    print("PROCESSING COMPLETE")
    print(f"{'='*80}")
    print(f"Total processed: {len(to_process)}")
    print(f"Successful: {successful}")
    print(f"Failed: {failed}")
    print(f"Time elapsed: {elapsed/60:.1f} minutes ({elapsed/3600:.2f} hours)")
    print(f"Average time per structure: {elapsed/len(to_process):.2f} seconds")
    
    if failed_pdb_ids:
        print(f"\nFailed PDB IDs ({len(failed_pdb_ids)}):")
        for pdb_id in failed_pdb_ids[:20]:
            print(f"  - {pdb_id}")
        if len(failed_pdb_ids) > 20:
            print(f"  ... and {len(failed_pdb_ids) - 20} more")


if __name__ == '__main__':
    main()

