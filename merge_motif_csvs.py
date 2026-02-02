#!/usr/bin/env python3
"""
Merge all individual motif CSV files into a single final CSV.
This script combines all CSV files from motif_csvs/ directory into scores_motifs_summary.csv
"""

import csv
import sys
from pathlib import Path
from collections import OrderedDict

def merge_motif_csvs(csv_dir='motif_csvs', output_file='scores_motifs_summary.csv'):
    """
    Merge all individual motif CSV files into a single CSV.
    
    Args:
        csv_dir: Directory containing individual motif CSV files
        output_file: Path to output merged CSV file
    """
    csv_dir_path = Path(csv_dir)
    
    if not csv_dir_path.exists():
        print(f"Error: Directory '{csv_dir}' does not exist!")
        sys.exit(1)
    
    # Find all CSV files
    csv_files = sorted(csv_dir_path.glob("*.csv"))
    
    if not csv_files:
        print(f"Error: No CSV files found in '{csv_dir}'!")
        sys.exit(1)
    
    print(f"Found {len(csv_files)} CSV files to merge")
    print(f"Output file: {output_file}")
    print("-" * 80)
    
    # Collect all rows and fieldnames
    all_rows = []
    fieldnames_set = set()
    
    # Read all CSV files
    for csv_file in csv_files:
        try:
            with open(csv_file, 'r', newline='') as f:
                reader = csv.DictReader(f)
                # Collect fieldnames
                if reader.fieldnames:
                    fieldnames_set.update(reader.fieldnames)
                # Read all rows
                for row in reader:
                    all_rows.append(row)
        except Exception as e:
            print(f"Warning: Error reading {csv_file.name}: {e}")
            continue
    
    if not all_rows:
        print("Error: No data rows found in any CSV files!")
        sys.exit(1)
    
    # Use OrderedDict to preserve column order (use first file's order if available)
    # Otherwise, sort alphabetically
    if csv_files:
        try:
            with open(csv_files[0], 'r', newline='') as f:
                reader = csv.DictReader(f)
                if reader.fieldnames:
                    fieldnames = list(OrderedDict.fromkeys(reader.fieldnames))
                    # Add any missing fieldnames
                    for fn in sorted(fieldnames_set):
                        if fn not in fieldnames:
                            fieldnames.append(fn)
                else:
                    fieldnames = sorted(fieldnames_set)
        except:
            fieldnames = sorted(fieldnames_set)
    else:
        fieldnames = sorted(fieldnames_set)
    
    # Write merged CSV
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    try:
        with open(output_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            
            # Write all rows
            for row in all_rows:
                # Ensure all fieldnames are present in row
                complete_row = {fn: row.get(fn, '') for fn in fieldnames}
                writer.writerow(complete_row)
        
        print(f"\n{'='*80}")
        print(f"SUCCESS: Merged {len(csv_files)} CSV files")
        print(f"Total rows: {len(all_rows)}")
        print(f"Output: {output_file}")
        print(f"{'='*80}")
        return True
        
    except Exception as e:
        print(f"\nError writing merged CSV: {e}")
        import traceback
        print(traceback.format_exc())
        return False

if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Merge all individual motif CSV files into a single CSV"
    )
    parser.add_argument(
        '--csv-dir',
        default='motif_csvs',
        help='Directory containing individual motif CSV files (default: motif_csvs)'
    )
    parser.add_argument(
        '--output',
        default='scores_motifs_summary.csv',
        help='Output merged CSV file (default: scores_motifs_summary.csv)'
    )
    
    args = parser.parse_args()
    
    success = merge_motif_csvs(args.csv_dir, args.output)
    sys.exit(0 if success else 1)

