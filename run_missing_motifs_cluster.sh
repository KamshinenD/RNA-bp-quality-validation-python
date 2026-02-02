#!/bin/bash
#SBATCH --job-name=missing_motifs
#SBATCH --output=logs/missing_motifs_%A_%a.out
#SBATCH --error=logs/missing_motifs_%A_%a.err
#SBATCH --time=08:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --array=0-599

# ==============================================================================
# Process Missing Motifs - SLURM Array Job
# ==============================================================================
# Processes missing motifs by splitting work across array tasks.
# Each task writes to its own CSV to avoid conflicts, then merge afterward.
#
# Usage:
#   # First, generate the missing motifs list
#   python process_missing_motifs.py --find-only
#
#   # Then submit the job
#   sbatch run_missing_motifs_cluster.sh
#
#   # After completion, merge results
#   python merge_missing_motifs.py
# ==============================================================================

# Create directories
mkdir -p logs
mkdir -p missing_motifs_results

# Configuration
MISSING_LIST_FILE="all_missing_motifs_list.txt"
TASK_ID=$SLURM_ARRAY_TASK_ID
NUM_TASKS=$SLURM_ARRAY_TASK_COUNT
OUTPUT_CSV="missing_motifs_results/missing_task_${TASK_ID}.csv"

echo "========================================"
echo "Array Task: $TASK_ID / $((NUM_TASKS - 1))"
echo "========================================"

# Check if missing motifs list exists
if [ ! -f "$MISSING_LIST_FILE" ]; then
    echo "ERROR: Missing motifs list file '$MISSING_LIST_FILE' not found!"
    echo "Run: python process_all_missing_motifs.py"
    exit 1
fi

# Count total missing motifs
TOTAL_MISSING=$(wc -l < "$MISSING_LIST_FILE")
echo "Total missing motifs: $TOTAL_MISSING"

# Calculate which motifs this task should process
MOTIFS_PER_TASK=$(( (TOTAL_MISSING + NUM_TASKS - 1) / NUM_TASKS ))
START_LINE=$((TASK_ID * MOTIFS_PER_TASK + 1))
END_LINE=$((START_LINE + MOTIFS_PER_TASK - 1))

# Ensure we don't go beyond the last line
if [ $END_LINE -gt $TOTAL_MISSING ]; then
    END_LINE=$TOTAL_MISSING
fi

# Skip if no motifs to process
if [ $START_LINE -gt $TOTAL_MISSING ]; then
    echo "No motifs to process for task $TASK_ID"
    exit 0
fi

MOTIFS_TO_PROCESS=$((END_LINE - START_LINE + 1))
echo "This task will process: motifs from line $START_LINE to $END_LINE ($MOTIFS_TO_PROCESS motifs)"
echo ""

SUCCESS=0
FAILED=0

# Process motifs assigned to this task
LINE_NUM=0
while IFS= read -r line; do
    ((LINE_NUM++))
    
    # Skip lines before our start
    if [ $LINE_NUM -lt $START_LINE ]; then
        continue
    fi
    
    # Stop after our end
    if [ $LINE_NUM -gt $END_LINE ]; then
        break
    fi
    
    # Parse line: motif_name|cif_file_path
    MOTIF_NAME=$(echo "$line" | cut -d'|' -f1)
    CIF_FILE=$(echo "$line" | cut -d'|' -f2)
    
    PROGRESS=$((LINE_NUM - START_LINE + 1))
    PCT=$(echo "scale=1; $PROGRESS * 100 / $MOTIFS_TO_PROCESS" | bc)
    
    echo "[$PROGRESS/$MOTIFS_TO_PROCESS] ($PCT%) $MOTIF_NAME"
    
    # Extract PDB ID (4 chars: 1 digit + 3 alphanumeric)
    PDB_ID=$(echo "$MOTIF_NAME" | grep -oE '[0-9][A-Z0-9]{3}' | tail -1)
    
    if [ -z "$PDB_ID" ]; then
        echo "  ✗ Could not extract PDB ID"
        ((FAILED++))
        continue
    fi
    
    # Parse CIF file for chain and residue range
    if [ ! -f "$CIF_FILE" ]; then
        echo "  ✗ CIF file not found: $CIF_FILE"
        ((FAILED++))
        continue
    fi
    
    CHAIN=$(grep "^ATOM" "$CIF_FILE" | head -1 | awk '{print $6}')
    RES_NUMS=$(grep "^ATOM" "$CIF_FILE" | awk '{print $5}' | sort -n | uniq)
    START_RES=$(echo "$RES_NUMS" | head -1)
    END_RES=$(echo "$RES_NUMS" | tail -1)
    
    if [ -z "$CHAIN" ] || [ -z "$START_RES" ] || [ -z "$END_RES" ]; then
        echo "  ✗ Could not parse CIF file"
        ((FAILED++))
        continue
    fi
    
    # Extract exact residue IDs from CIF file (format: CHAIN-BASE-RESNUM-)
    # This gives us the exact list of residues in the motif
    MOTIF_RESIDUES=$(grep "^ATOM" "$CIF_FILE" | awk -v chain="$CHAIN" '{
        res_num = $5
        base = $4
        # Format: CHAIN-BASE-RESNUM-
        printf "%s-%s-%s-,", chain, base, res_num
    }' | sort -u | tr '\n' ',' | sed 's/,$//')
    
    echo "  PDB: $PDB_ID, Chain: $CHAIN, Range: $START_RES-$END_RES"
    echo "  Motif residues: $(echo "$MOTIF_RESIDUES" | tr ',' '\n' | wc -l | tr -d ' ') residues from CIF"
    
    # Create unique report file for this motif to avoid conflicts
    MOTIF_REPORT="motif_report_${MOTIF_NAME//[^a-zA-Z0-9]/_}.json"
    
    # Run motif scoring with exact residue list from CIF file
    python app.py --pdb_id "$PDB_ID" --motif "$START_RES" "$END_RES" --chain "$CHAIN" --motif-residues "$MOTIF_RESIDUES" 2>&1 | tail -5
    
    exit_code=$?
    
    # Immediately copy motif_report.json to unique file to prevent overwriting
    if [ -f "motif_report.json" ]; then
        cp "motif_report.json" "$MOTIF_REPORT"
    fi
    
    if [ $exit_code -eq 0 ] || [ $exit_code -eq 1 ]; then
        # Extract data from the unique motif_report file
        if [ -f "$MOTIF_REPORT" ]; then
            # Check if file is empty or too small to be valid JSON
            if [ ! -s "$MOTIF_REPORT" ]; then
                echo "  ✗ Empty report file for $MOTIF_NAME"
                rm -f "$MOTIF_REPORT"
                ((FAILED++))
            else
                python -c "
import json
import csv
from pathlib import Path
import sys

motif_name = '$MOTIF_NAME'
csv_file = '$OUTPUT_CSV'
report_file = '$MOTIF_REPORT'

try:
    with open(report_file, 'r') as f:
        content = f.read().strip()
        if not content:
            print(f'  ✗ Empty JSON file for {motif_name}')
            sys.exit(1)
        data = json.loads(content)
except json.JSONDecodeError as e:
    print(f'  ✗ Invalid JSON in report file for {motif_name}: {e}')
    sys.exit(1)
except Exception as e:
    print(f'  ✗ Error reading report file for {motif_name}: {e}')
    sys.exit(1)
    
# VERIFY: Check that the motif name matches the PDB ID in the report
report_pdb = data.get('pdb_id', '').upper()
expected_pdb = '$PDB_ID'.upper()
if report_pdb != expected_pdb:
    print(f'  ⚠ ERROR: PDB mismatch! Expected {expected_pdb}, got {report_pdb}')
    print(f'  This indicates a race condition - wrong motif was scored!')
    print(f'  Skipping this entry to prevent data corruption.')
    exit(1)

# Verify chain and residue range match
report_chain = str(data.get('motif_chain', ''))
expected_chain = '$CHAIN'
if report_chain != expected_chain:
    print(f'  ⚠ WARNING: Chain mismatch! Expected {expected_chain}, got {report_chain}')

report_range = data.get('motif_range', '')
expected_range = '$START_RES-$END_RES'
if report_range != expected_range:
    print(f'  ⚠ WARNING: Range mismatch! Expected {expected_range}, got {report_range}')

# Calculate pairing percentage
motif_length = data.get('motif_num_nucleotides', 0)
num_paired = data.get('num_paired_nucleotides', 0)
pairing_pct = round((num_paired / motif_length * 100), 1) if motif_length > 0 else 0

# Build row
motif_data = {
    'Motif_Name': motif_name,
    'PDB_ID': data.get('pdb_id', 'N/A'),
    'Chain': data.get('motif_chain', 'N/A'),
    'Residue_Range': data.get('motif_range', 'N/A'),
    'Motif_Score': data.get('motif_score', 0),
    'Full_Structure_Score': data.get('full_structure_score', 0),
    'Score_Difference': data.get('score_difference', 0),
    'Motif_Length': data.get('motif_num_nucleotides', 0),
    'Full_Structure_Length': data.get('full_structure_num_nucleotides', 0),
    'Total_Base_Pairs': data.get('total_base_pairs', 0),
    'Num_Problematic_BPs': data.get('num_problematic_bps', 0),
    'Num_Paired_Nucleotides': num_paired,
    'Pairing_Percentage': pairing_pct,
    'Geom_Misaligned': data.get('geometry_issues', {}).get('misaligned', 0),
    'Geom_Twisted': data.get('geometry_issues', {}).get('twisted', 0),
    'Geom_NonCoplanar': data.get('geometry_issues', {}).get('non_coplanar', 0),
    'Geom_PoorHBond': data.get('geometry_issues', {}).get('poor_hbond', 0),
    'Geom_ZeroHBond': data.get('geometry_issues', {}).get('zero_hbond', 0),
    'HBond_BadDistance': data.get('hbond_issues', {}).get('bad_distance', 0),
    'HBond_BadAngles': data.get('hbond_issues', {}).get('bad_angles', 0),
    'HBond_BadDihedral': data.get('hbond_issues', {}).get('bad_dihedral', 0),
    'HBond_WeakQuality': data.get('hbond_issues', {}).get('weak_quality', 0),
    'HBond_IncorrectCount': data.get('hbond_issues', {}).get('incorrect_count', 0),
}

# Dominant issues
all_issues = {}
for k, v in data.get('geometry_issues', {}).items():
    if v > 0:
        all_issues[f'geom_{k}'] = v
for k, v in data.get('hbond_issues', {}).items():
    if v > 0:
        all_issues[f'hbond_{k}'] = v
sorted_issues = sorted(all_issues.items(), key=lambda x: x[1], reverse=True)
motif_data['Dominant_Issues'] = '; '.join([f'{issue}({count})' for issue, count in sorted_issues]) or 'none'

# Write to CSV
fieldnames = list(motif_data.keys())
file_exists = Path(csv_file).exists()

with open(csv_file, 'a', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    if not file_exists:
        writer.writeheader()
    writer.writerow(motif_data)
"
                if [ $? -eq 0 ]; then
                    echo "  ✓ Success"
                    # Clean up the temporary report file
                    rm -f "$MOTIF_REPORT"
                    ((SUCCESS++))
                else
                    echo "  ✗ Data validation failed or JSON error"
                    rm -f "$MOTIF_REPORT"
                    ((FAILED++))
                fi
            fi
        else
            echo "  ✗ No motif report generated for $MOTIF_NAME"
            ((FAILED++))
        fi
    else
        echo "  ✗ Scoring failed (exit code: $exit_code)"
        ((FAILED++))
    fi
    
    # Clean up any leftover motif_report.json to prevent contamination
    rm -f motif_report.json report.json
done < "$MISSING_LIST_FILE"

echo ""
echo "========================================"
echo "Task $TASK_ID Summary"
echo "========================================"
echo "Processed: $((SUCCESS + FAILED)) motifs"
echo "Success:   $SUCCESS"
echo "Failed:    $FAILED"
echo "Output:    $OUTPUT_CSV"
echo "========================================"

exit 0


