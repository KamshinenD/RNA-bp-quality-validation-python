#!/usr/bin/env python3
"""
FAST parallel metadata caching - downloads multiple RNAs simultaneously.

Uses multiprocessing to download nucleotide counts and validation metrics
in parallel, achieving 10-20x speedup.
"""

import json
import time
import requests
from pathlib import Path
from typing import Dict, Optional, Tuple
import sys
from multiprocessing import Pool, Manager
from functools import partial

# Import data loader to reuse its methods
sys.path.insert(0, str(Path(__file__).parent))
from config import Config
from utils.data_loader import DataLoader


def get_all_pdb_ids(basepairs_dir: str) -> list:
    """Get all PDB IDs from basepairs directory."""
    basepairs_path = Path(basepairs_dir)
    if not basepairs_path.exists():
        return []
    
    json_files = list(basepairs_path.glob('*.json'))
    return sorted(set([f.stem.upper() for f in json_files]))


def download_nucleotide_count(args: Tuple[str, int]) -> Tuple[str, int]:
    """Download nucleotide count for a single PDB ID (worker function)."""
    pdb_id, worker_id = args
    import tempfile
    import os
    import shutil
    
    # Create a unique temporary file for this worker process
    # Use process ID + worker ID to ensure uniqueness
    import os as os_module
    process_id = os_module.getpid()
    temp_cif = f"temp_structure_{process_id}_{worker_id}.cif"
    
    try:
        # Create a temporary data loader for this worker
        config = Config()
        data_loader = DataLoader(config)
        
        # Temporarily override get_nucleotide_count to use our unique temp file
        original_method = data_loader.get_nucleotide_count
        
        def get_nuc_count_with_unique_file(pdb_id):
            """Wrapper that uses unique temp file."""
            # Download CIF to unique temp file
            if not data_loader.download_cif(pdb_id, temp_cif):
                return 0
            
            # Count nucleotides
            count = data_loader.count_nucleotides_from_cif(temp_cif)
            
            # Clean up temp file
            try:
                if os.path.exists(temp_cif):
                    os.unlink(temp_cif)
            except:
                pass
            
            return count
        
        # Use our wrapper
        count = get_nuc_count_with_unique_file(pdb_id)
        return (pdb_id, count)
        
    except Exception as e:
        # Clean up on error
        try:
            if os.path.exists(temp_cif):
                os.unlink(temp_cif)
        except:
            pass
        return (pdb_id, 0)


def download_validation_metrics(pdb_id: str) -> Tuple[str, Optional[Dict]]:
    """Download validation metrics for a single PDB ID (worker function).
    
    Uses the SAME extraction logic as data_loader.get_validation_metrics()
    to ensure consistency.
    """
    try:
        # Import data loader to use its extraction method
        config = Config()
        data_loader = DataLoader(config)
        
        # Use the same method as data_loader
        metrics = data_loader.get_validation_metrics(pdb_id)
        return (pdb_id, metrics)
            
    except Exception as e:
        return (pdb_id, None)


def process_batch_parallel(items, worker_func, num_workers=20, desc="Processing"):
    """Process items in parallel batches."""
    results = {}
    total = len(items)
    
    print(f"\n{desc} {total} items with {num_workers} parallel workers...")
    start_time = time.time()
    
    with Pool(processes=num_workers) as pool:
        # Process in batches to show progress
        batch_size = 100
        for i in range(0, total, batch_size):
            batch = items[i:i+batch_size]
            batch_results = pool.map(worker_func, batch)
            
            # Update results
            for pdb_id, result in batch_results:
                results[pdb_id] = result
            
            # Show progress
            completed = min(i + batch_size, total)
            elapsed = time.time() - start_time
            rate = completed / elapsed if elapsed > 0 else 0
            eta = (total - completed) / rate if rate > 0 else 0
            
            print(f"  Progress: {completed}/{total} ({completed*100/total:.1f}%) | "
                  f"Rate: {rate:.1f}/sec | ETA: {eta/60:.1f} min")
    
    elapsed = time.time() - start_time
    print(f"  Completed in {elapsed/60:.1f} minutes ({elapsed:.1f} seconds)")
    
    return results


def main():
    """Main function to cache all metadata in parallel."""
    CACHE_FILE = 'metadata_cache.json'
    BASEPAIRS_DIR = 'data/basepairs'
    NUM_WORKERS = 20  # Number of parallel downloads
    
    print("="*80)
    print("FAST PARALLEL METADATA CACHING")
    print(f"(Using {NUM_WORKERS} parallel workers)")
    print("="*80)
    
    # Load existing cache
    cache = {}
    if Path(CACHE_FILE).exists():
        try:
            with open(CACHE_FILE, 'r') as f:
                cache = json.load(f)
            print(f"\nLoaded existing cache: {CACHE_FILE}")
            print(f"  Cached nucleotide counts: {len(cache.get('nucleotide_counts', {}))}")
            print(f"  Cached validation metrics: {len(cache.get('validation_metrics', {}))}")
        except Exception as e:
            print(f"Warning: Could not load cache: {e}")
            cache = {}
    
    # Get all PDB IDs
    print(f"\nScanning {BASEPAIRS_DIR}...")
    all_pdb_ids = get_all_pdb_ids(BASEPAIRS_DIR)
    print(f"Found {len(all_pdb_ids)} PDB IDs")
    
    # Check what needs to be cached
    nuc_counts = cache.get('nucleotide_counts', {})
    val_metrics = cache.get('validation_metrics', {})
    
    need_nuc = [pdb_id for pdb_id in all_pdb_ids if pdb_id not in nuc_counts]
    need_val = [pdb_id for pdb_id in all_pdb_ids if pdb_id not in val_metrics or val_metrics[pdb_id] is None]
    
    print(f"\nNucleotide counts needed: {len(need_nuc)}")
    print(f"Validation metrics needed: {len(need_val)}")
    
    if not need_nuc and not need_val:
        print("\n✓ All metadata already cached!")
        return
    
    # Estimate time
    total_needed = max(len(need_nuc), len(need_val))
    # With NUM_WORKERS parallel, each download takes ~3 seconds
    # So we process NUM_WORKERS every ~3 seconds
    estimated_seconds = (total_needed / NUM_WORKERS) * 3
    estimated_minutes = estimated_seconds / 60
    
    print(f"\nEstimated time: ~{estimated_minutes:.1f} minutes")
    print(f"(With {NUM_WORKERS} parallel workers)")
    
    # Ask for confirmation
    response = input(f"\nProceed with parallel caching? (yes/no): ")
    if response.lower() not in ['yes', 'y']:
        print("Cancelled.")
        return
    
    # Cache nucleotide counts in parallel
    if need_nuc:
        print(f"\n{'='*80}")
        print(f"Caching nucleotide counts ({len(need_nuc)} remaining)...")
        print(f"{'='*80}")
        
        # Prepare arguments (pdb_id, worker_id for each worker)
        # Each worker gets a unique ID to avoid temp file conflicts
        args_list = [(pdb_id, i % NUM_WORKERS) for i, pdb_id in enumerate(need_nuc)]
        
        # Process in parallel
        results = process_batch_parallel(
            args_list,
            download_nucleotide_count,
            num_workers=NUM_WORKERS,
            desc="Downloading nucleotide counts"
        )
        
        # Update cache
        cache.setdefault('nucleotide_counts', {}).update(results)
        
        # Save cache
        with open(CACHE_FILE, 'w') as f:
            json.dump(cache, f, indent=2)
        print(f"  Cache saved: {len(results)} nucleotide counts")
    
    # Cache validation metrics in parallel
    if need_val:
        print(f"\n{'='*80}")
        print(f"Caching validation metrics ({len(need_val)} remaining)...")
        print(f"{'='*80}")
        
        # Process in parallel
        results = process_batch_parallel(
            need_val,
            download_validation_metrics,
            num_workers=NUM_WORKERS,
            desc="Downloading validation metrics"
        )
        
        # Update cache
        cache.setdefault('validation_metrics', {}).update(results)
        
        # Save cache
        with open(CACHE_FILE, 'w') as f:
            json.dump(cache, f, indent=2)
        print(f"  Cache saved: {len(results)} validation metrics")
    
    # Final save
    with open(CACHE_FILE, 'w') as f:
        json.dump(cache, f, indent=2)
    
    print(f"\n{'='*80}")
    print("CACHING COMPLETE")
    print(f"{'='*80}")
    print(f"Cache saved to: {CACHE_FILE}")
    print(f"  Nucleotide counts: {len(cache.get('nucleotide_counts', {}))}")
    print(f"  Validation metrics: {len(cache.get('validation_metrics', {}))}")
    print(f"\nYou can now run batch processing with cached data!")


if __name__ == '__main__':
    main()
