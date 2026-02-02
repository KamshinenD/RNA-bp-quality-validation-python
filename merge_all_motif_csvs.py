#!/usr/bin/env python3
"""
Merge all motif CSV files from cluster jobs into one master file.
Removes duplicates and ensures data integrity.
"""

import csv
from pathlib import Path
import argparse

def merge_csvs(input_dir='motif_results_individual', output_file='all_motifs_scored.csv'):
    """
    Merge all CSV files from motif tasks into a single file.
    Preserves existing data in output_file if it exists.
    Removes duplicates based on Motif_Name.
    
    Args:
        input_dir: Directory containing individual task CSV files
        output_file: Output merged CSV file (will preserve existing data if it exists)
    """
    input_path = Path(input_dir)
    
    if not input_path.exists():
        print(f"Error: Directory '{input_dir}' not found!")
        return
    
    # Find all CSV files
    csv_files = sorted(input_path.glob('*.csv'))
    
    if not csv_files:
        print(f"No CSV files found in '{input_dir}/'")
        return
    
    print("=" * 60)
    print("Merging All Motif Scores")
    print("=" * 60)
    print(f"Found {len(csv_files)} CSV files to merge")
    
    # Read existing merged file if it exists
    existing_motifs = set()
    existing_rows = []
    fieldnames = None
    
    if Path(output_file).exists():
        print(f"\nReading existing merged file: {output_file}")
        try:
            with open(output_file, 'r', newline='') as f:
                reader = csv.DictReader(f)
                if fieldnames is None:
                    fieldnames = reader.fieldnames
                for row in reader:
                    motif_name = row.get('Motif_Name', '').strip()
                    if motif_name:
                        existing_motifs.add(motif_name)
                        existing_rows.append(row)
            print(f"  ✓ Found {len(existing_rows)} existing motifs in merged file")
        except Exception as e:
            print(f"  ⚠ Warning: Could not read existing file: {e}")
    
    # Read all rows from new task CSV files
    all_rows = existing_rows.copy()  # Start with existing rows
    seen_motifs = existing_motifs.copy()  # Start with existing motifs
    duplicates = 0
    
    for csv_file in csv_files:
        try:
            with open(csv_file, 'r', newline='') as f:
                reader = csv.DictReader(f)
                
                if fieldnames is None:
                    fieldnames = reader.fieldnames
                
                for row in reader:
                    motif_name = row.get('Motif_Name', '').strip()
                    if motif_name:
                        # Check for duplicates (skip if already in existing file or already seen)
                        if motif_name in seen_motifs:
                            duplicates += 1
                            continue
                        seen_motifs.add(motif_name)
                        all_rows.append(row)
                        
            count = sum(1 for _ in open(csv_file)) - 1  # Subtract header
            print(f"  ✓ {csv_file.name}: {count} rows")
            
        except Exception as e:
            print(f"  ✗ Error reading {csv_file.name}: {e}")
    
    if not all_rows:
        print("No data found in CSV files!")
        return
    
    # Write merged file
    try:
        with open(output_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(all_rows)
        
        print(f"\n{'='*60}")
        print("MERGE COMPLETE")
        print(f"{'='*60}")
        print(f"Existing motifs: {len(existing_rows)}")
        print(f"New motifs added: {len(all_rows) - len(existing_rows)}")
        print(f"Total unique motifs: {len(all_rows)}")
        if duplicates > 0:
            print(f"Duplicates removed: {duplicates}")
        print(f"Output file: {output_file}")
        print(f"{'='*60}")
        
    except Exception as e:
        print(f"Error writing merged file: {e}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge all motif CSV files from cluster jobs")
    
    parser.add_argument(
        '--input-dir',
        type=str,
        default='motif_results_individual',
        help='Directory containing individual CSV files (default: motif_results_individual)'
    )
    
    parser.add_argument(
        '--output',
        type=str,
        default='all_motifs_scored.csv',
        help='Output merged CSV file (default: all_motifs_scored.csv)'
    )
    
    args = parser.parse_args()
    
    merge_csvs(args.input_dir, args.output)

