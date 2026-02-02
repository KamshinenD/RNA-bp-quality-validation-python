#!/usr/bin/env python3
"""
Pre-download and cache validation metrics and nucleotide counts for all RNAs.

This allows batch processing to use cached data instead of making network calls.
Run this once before batch processing to speed things up.
"""

import json
import time
import requests
from pathlib import Path
from typing import Dict, Optional
import sys

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


def cache_nucleotide_count(pdb_id: str, data_loader: DataLoader, cache: Dict) -> int:
    """Get nucleotide count, using cache if available."""
    if pdb_id in cache.get('nucleotide_counts', {}):
        return cache['nucleotide_counts'][pdb_id]
    
    print(f"  Downloading nucleotide count for {pdb_id}...", end=' ', flush=True)
    try:
        count = data_loader.get_nucleotide_count(pdb_id)
        cache.setdefault('nucleotide_counts', {})[pdb_id] = count
        print(f"✓ ({count})")
        return count
    except Exception as e:
        print(f"✗ ({e})")
        cache.setdefault('nucleotide_counts', {})[pdb_id] = 0
        return 0


def cache_validation_metrics(pdb_id: str, cache: Dict) -> Optional[Dict]:
    """Get validation metrics, using cache if available."""
    if pdb_id in cache.get('validation_metrics', {}):
        return cache['validation_metrics'][pdb_id]
    
    print(f"  Downloading validation metrics for {pdb_id}...", end=' ', flush=True)
    try:
        url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id.upper()}"
        response = requests.get(url, timeout=30)
        
        if response.status_code == 200:
            data = response.json()
            metrics = {}
            
            # Extract relevant fields
            if 'rcsb_entry_info' in data:
                info = data['rcsb_entry_info']
                metrics['experimental_method'] = info.get('experimental_method', '')
                metrics['deposition_date'] = info.get('deposition_date', '')
            
            if 'exptl' in data and data['exptl']:
                exptl = data['exptl'][0]
                metrics['structure_determination_method'] = exptl.get('method', '')
            
            if 'rcsb_entry_container_identifiers' in data:
                identifiers = data['rcsb_entry_container_identifiers']
                metrics['em_resolution'] = identifiers.get('em_resolution', '')
                metrics['em_diffraction_resolution'] = identifiers.get('em_diffraction_resolution', '')
            
            if 'refine' in data and data['refine']:
                refine = data['refine'][0]
                metrics['refinement_resolution'] = refine.get('ls_d_res_high', '')
                metrics['r_free'] = refine.get('ls_R_free_R_factor', '')
                metrics['r_work'] = refine.get('ls_R_factor_R_work', '')
            
            if 'pdbx_vrpt_summary' in data:
                summary = data['pdbx_vrpt_summary']
                if summary:
                    metrics['clashscore'] = summary[0].get('clashscore', '')
                    metrics['rnasuiteness'] = summary[0].get('rnasuiteness', '')
                    metrics['angles_rmsz'] = summary[0].get('angles_rmsz', '')
                    metrics['bonds_rmsz'] = summary[0].get('bonds_rmsz', '')
                    metrics['percent_ramachandran_outliers'] = summary[0].get('percent_ramachandran_outliers', '')
                    metrics['percent_rotamer_outliers'] = summary[0].get('percent_rotamer_outliers', '')
            
            cache.setdefault('validation_metrics', {})[pdb_id] = metrics
            print("✓")
            return metrics
        else:
            print(f"✗ (HTTP {response.status_code})")
            cache.setdefault('validation_metrics', {})[pdb_id] = None
            return None
            
    except Exception as e:
        print(f"✗ ({e})")
        cache.setdefault('validation_metrics', {})[pdb_id] = None
        return None


def main():
    """Main function to cache all metadata."""
    CACHE_FILE = 'metadata_cache.json'
    BASEPAIRS_DIR = 'data/basepairs'
    DELAY = 0.5  # Delay between API calls to be respectful
    
    print("="*80)
    print("METADATA CACHING FOR BATCH PROCESSING")
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
    
    # Initialize data loader
    config = Config()
    data_loader = DataLoader(config)
    
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
    estimated_minutes = total_needed * (DELAY + 2) / 60  # ~2 seconds per call + delay
    print(f"\nEstimated time: ~{estimated_minutes:.1f} minutes")
    
    # Ask for confirmation
    response = input(f"\nProceed with caching metadata? (yes/no): ")
    if response.lower() not in ['yes', 'y']:
        print("Cancelled.")
        return
    
    # Cache nucleotide counts
    if need_nuc:
        print(f"\n{'='*80}")
        print(f"Caching nucleotide counts ({len(need_nuc)} remaining)...")
        print(f"{'='*80}")
        
        for i, pdb_id in enumerate(need_nuc, 1):
            print(f"[{i}/{len(need_nuc)}] {pdb_id}:", end=' ')
            cache_nucleotide_count(pdb_id, data_loader, cache)
            
            # Save cache periodically
            if i % 50 == 0:
                with open(CACHE_FILE, 'w') as f:
                    json.dump(cache, f, indent=2)
                print(f"  (Cache saved - {i}/{len(need_nuc)} complete)")
            
            time.sleep(DELAY)
    
    # Cache validation metrics
    if need_val:
        print(f"\n{'='*80}")
        print(f"Caching validation metrics ({len(need_val)} remaining)...")
        print(f"{'='*80}")
        
        for i, pdb_id in enumerate(need_val, 1):
            print(f"[{i}/{len(need_val)}] {pdb_id}:", end=' ')
            cache_validation_metrics(pdb_id, cache)
            
            # Save cache periodically
            if i % 50 == 0:
                with open(CACHE_FILE, 'w') as f:
                    json.dump(cache, f, indent=2)
                print(f"  (Cache saved - {i}/{len(need_val)} complete)")
            
            time.sleep(DELAY)
    
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

