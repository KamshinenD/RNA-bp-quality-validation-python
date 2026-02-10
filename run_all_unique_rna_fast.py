#!/usr/bin/env python3
"""
FAST batch processing for UNIQUE RNAs only - skips network calls for speed.

This version:
- Processes only unique_pdb_id's from uniqueRNAs.csv
- Skips nucleotide count downloads (uses cache)
- Skips validation metrics API calls (uses cache)
- Only does the essential scoring
- Outputs to scores_summary_unique.csv
"""

import os
import sys
import subprocess
import time
import csv
import json
import pandas as pd
from pathlib import Path
from typing import Set, List, Optional, Dict

# Import the scoring components directly
sys.path.insert(0, str(Path(__file__).parent))
from config import Config
from utils.data_loader import DataLoader
from scorer2 import Scorer
from utils.report_generator import ReportGenerator


def get_unique_pdb_ids(unique_csv: str = 'uniqueRNAs.csv') -> List[str]:
    """Get list of unique PDB IDs from uniqueRNAs.csv."""
    unique_ids = []
    
    if not os.path.exists(unique_csv):
        print(f"Error: {unique_csv} not found!")
        return []
    
    try:
        df = pd.read_csv(unique_csv)
        if 'unique_pdb_id' not in df.columns:
            print(f"Error: 'unique_pdb_id' column not found in {unique_csv}")
            return []
        
        # Get unique PDB IDs and convert to uppercase
        unique_ids = df['unique_pdb_id'].dropna().str.strip().str.upper().unique().tolist()
        
        return sorted(unique_ids)
    
    except Exception as e:
        print(f"Error reading {unique_csv}: {e}")
        return []


def get_processed_pdb_ids(csv_file: str) -> Set[str]:
    """Get set of PDB IDs already in scores_summary_unique.csv."""
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


def get_validation_metrics(pdb_id: str, cache: dict, data_loader) -> Optional[dict]:
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
        validation_metrics = get_validation_metrics(pdb_id, cache, data_loader)

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
            'scores_summary_unique.csv',  # Changed output file
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
    UNIQUE_CSV = 'data/uniqueRNAs.csv'
    SCORES_CSV = 'scores_summary_unique.csv'  # Changed output file
    DELAY_BETWEEN_RUNS = 0.5
    
    print("="*80)
    print("FAST BATCH RNA SCORING PROCESSOR - UNIQUE RNAs ONLY")
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
    
    # Get unique PDB IDs from uniqueRNAs.csv
    print(f"\nReading unique PDB IDs from {UNIQUE_CSV}...")
    all_pdb_ids = get_unique_pdb_ids(UNIQUE_CSV)
    
    if not all_pdb_ids:
        print("No unique PDB IDs found!")
        return
    
    print(f"Found {len(all_pdb_ids)} unique PDB IDs")
    
    # Verify that basepair files exist
    basepairs_path = Path(BASEPAIRS_DIR)
    available_ids = []
    missing_ids = []
    
    for pdb_id in all_pdb_ids:
        json_file = basepairs_path / f"{pdb_id.lower()}.json"
        if json_file.exists():
            available_ids.append(pdb_id)
        else:
            missing_ids.append(pdb_id)
    
    if missing_ids:
        print(f"\nWarning: {len(missing_ids)} PDB IDs don't have basepair files:")
        for pdb_id in missing_ids[:10]:
            print(f"  - {pdb_id}")
        if len(missing_ids) > 10:
            print(f"  ... and {len(missing_ids) - 10} more")
    
    print(f"\n{len(available_ids)} structures have basepair files available")
    
    # Get already processed PDB IDs
    print(f"\nChecking {SCORES_CSV} for already processed structures...")
    processed = get_processed_pdb_ids(SCORES_CSV)
    print(f"Found {len(processed)} already processed structures")
    
    # Filter out already processed
    to_process = [pdb_id for pdb_id in available_ids if pdb_id not in processed]
    
    if not to_process:
        print("\n✓ All unique structures have already been processed!")
        return
    
    print(f"\n{len(to_process)} unique structures remaining to process")
    
    # Estimate time
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
    response = input(f"\nProceed with FAST processing of {len(to_process)} unique structures? (yes/no): ")
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
        
        # Small delay to prevent overheating
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
    if len(to_process) > 0:
        print(f"Average time per structure: {elapsed/len(to_process):.2f} seconds")
    
    if failed_pdb_ids:
        print(f"\nFailed PDB IDs ({len(failed_pdb_ids)}):")
        for pdb_id in failed_pdb_ids[:20]:
            print(f"  - {pdb_id}")
        if len(failed_pdb_ids) > 20:
            print(f"  ... and {len(failed_pdb_ids) - 20} more")
    
    print(f"\nResults saved to: {SCORES_CSV}")


if __name__ == '__main__':
    main()