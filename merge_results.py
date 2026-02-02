#!/usr/bin/env python3
"""Merge individual CSV results into a single summary file."""

import csv
from pathlib import Path
from collections import Counter
import sys


def merge_csv_files(results_dir="results_individual", output_file="scores_summary_merged.csv"):
    """
    Merge all individual CSV files into one summary.
    
    Args:
        results_dir: Directory containing individual CSV files
        output_file: Output merged CSV file
    """
    results_path = Path(results_dir)
    
    if not results_path.exists():
        print(f"Error: Directory {results_dir} does not exist!")
        return False
    
    # Find all CSV files
    csv_files = list(results_path.glob("*.csv"))
    
    if not csv_files:
        print(f"Error: No CSV files found in {results_dir}")
        return False
    
    print("="*70)
    print("MERGING INDIVIDUAL RESULTS")
    print("="*70)
    print(f"Found {len(csv_files)} individual result files")
    print(f"Output file: {output_file}\n")
    
    # Read all data
    all_rows = []
    header = None
    failed_files = []
    
    for i, csv_file in enumerate(sorted(csv_files), 1):
        try:
            with open(csv_file, 'r') as f:
                reader = csv.DictReader(f)
                
                # Get header from first file
                if header is None:
                    header = reader.fieldnames
                
                # Read all rows from this file
                for row in reader:
                    all_rows.append(row)
            
            # Progress indicator
            if i % 100 == 0:
                print(f"  Processed {i}/{len(csv_files)} files...")
                
        except Exception as e:
            print(f"  Warning: Failed to read {csv_file.name}: {e}")
            failed_files.append(csv_file.name)
    
    if not all_rows:
        print("\nError: No data rows found!")
        return False
    
    # Remove duplicates (keep last occurrence)
    seen_pdbs = {}
    for row in all_rows:
        pdb_id = row.get('PDB_ID', '')
        seen_pdbs[pdb_id] = row
    
    unique_rows = list(seen_pdbs.values())
    duplicates_removed = len(all_rows) - len(unique_rows)
    
    # Sort by PDB_ID
    unique_rows.sort(key=lambda x: x.get('PDB_ID', ''))
    
    # Write merged file
    try:
        with open(output_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=header)
            writer.writeheader()
            writer.writerows(unique_rows)
        
        print(f"\n{'='*70}")
        print("✓ MERGE COMPLETE!")
        print("="*70)
        print(f"Total structures: {len(unique_rows)}")
        if duplicates_removed > 0:
            print(f"Duplicates removed: {duplicates_removed}")
        if failed_files:
            print(f"Failed files: {len(failed_files)}")
        print(f"Output file: {output_file}\n")
        
        # Show statistics
        print_statistics(unique_rows)
        
        if failed_files:
            print("\nFailed files:")
            for fname in failed_files[:10]:
                print(f"  - {fname}")
            if len(failed_files) > 10:
                print(f"  ... and {len(failed_files) - 10} more")
        
        print("="*70 + "\n")
        return True
        
    except Exception as e:
        print(f"\nError writing output file: {e}")
        return False


def print_statistics(rows):
    """Print summary statistics."""
    if not rows:
        return
    
    print("Summary Statistics:")
    print("-"*70)
    
    # Score distribution
    grades = {'EXCELLENT': 0, 'GOOD': 0, 'FAIR': 0, 'POOR': 0}
    total_score = 0
    total_bps = 0
    
    for row in rows:
        try:
            score = float(row.get('Overall_Score', 0))
            total_score += score
            
            bps = int(row.get('Total_Base_Pairs', 0))
            total_bps += bps
            
            # Assign grade
            if score >= 85:
                grades['EXCELLENT'] += 1
            elif score >= 70:
                grades['GOOD'] += 1
            elif score >= 50:
                grades['FAIR'] += 1
            else:
                grades['POOR'] += 1
        except:
            pass
    
    total = len(rows)
    avg_score = total_score / total if total > 0 else 0
    avg_bps = total_bps / total if total > 0 else 0
    
    print(f"Average score: {avg_score:.1f}/100")
    print(f"Average base pairs per structure: {avg_bps:.1f}")
    print(f"\nGrade Distribution:")
    
    for grade in ['EXCELLENT', 'GOOD', 'FAIR', 'POOR']:
        count = grades[grade]
        pct = (count / total * 100) if total > 0 else 0
        print(f"  {grade:12s}: {count:5d} ({pct:5.1f}%)")


def main():
    """Main entry point."""
    import argparse
    
    parser = argparse.ArgumentParser(description="Merge individual CSV results")
    parser.add_argument(
        '--input-dir',
        default='results_individual',
        help='Directory containing individual CSV files'
    )
    parser.add_argument(
        '--output',
        default='scores_summary_merged.csv',
        help='Output merged CSV file'
    )
    
    args = parser.parse_args()
    
    success = merge_csv_files(args.input_dir, args.output)
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()

