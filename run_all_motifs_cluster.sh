#!/bin/bash
#SBATCH --job-name=motifs_fixed
#SBATCH --output=logs/motifs_%A_%a.out
#SBATCH --error=logs/motifs_%A_%a.err
#SBATCH --time=05:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --array=0-999

mkdir -p logs motif_results_individual full_structure_cache

TASK_ID=$SLURM_ARRAY_TASK_ID
NUM_TASKS=$SLURM_ARRAY_TASK_COUNT

# CRITICAL FIX: Create motif list file ONCE, use across all tasks
MOTIF_LIST_FILE="all_motifs_list.txt"

# Only task 0 creates the list (or use one pre-created on login node)
if [ ! -f "$MOTIF_LIST_FILE" ]; then
    if [ $TASK_ID -eq 0 ]; then
        echo "Task 0: Creating master motif list..."
        find motifs -name "*.cif" -type f > "$MOTIF_LIST_FILE.tmp"
        mv "$MOTIF_LIST_FILE.tmp" "$MOTIF_LIST_FILE"
        echo "Task 0: Created list with $(wc -l < "$MOTIF_LIST_FILE") motifs"
    else
        # Wait for task 0 to create the file
        echo "Task $TASK_ID: Waiting for motif list..."
        for i in {1..60}; do
            if [ -f "$MOTIF_LIST_FILE" ]; then
                break
            fi
            sleep 2
        done
        
        if [ ! -f "$MOTIF_LIST_FILE" ]; then
            echo "ERROR: Motif list file not created!"
            exit 1
        fi
    fi
fi

# Read motif list
mapfile -t MOTIF_FILES < "$MOTIF_LIST_FILE"
TOTAL_MOTIFS=${#MOTIF_FILES[@]}

echo "========================================"
echo "Task $TASK_ID: Found $TOTAL_MOTIFS motifs"
echo "========================================"

# Calculate work distribution
FILES_PER_TASK=$(( (TOTAL_MOTIFS + NUM_TASKS - 1) / NUM_TASKS ))
START_INDEX=$((TASK_ID * FILES_PER_TASK))
END_INDEX=$((START_INDEX + FILES_PER_TASK - 1))

if [ $END_INDEX -ge $TOTAL_MOTIFS ]; then
    END_INDEX=$((TOTAL_MOTIFS - 1))
fi

if [ $START_INDEX -ge $TOTAL_MOTIFS ]; then
    echo "No work for task $TASK_ID"
    exit 0
fi

MOTIFS_TO_PROCESS=$((END_INDEX - START_INDEX + 1))
echo "Processing motifs $START_INDEX to $END_INDEX ($MOTIFS_TO_PROCESS motifs)"

OUTPUT_CSV="motif_results_individual/motifs_task_${TASK_ID}.csv"

SUCCESS=0
FAILED=0
SKIPPED=0

for ((i=START_INDEX; i<=END_INDEX; i++)); do
    MOTIF_FILE="${MOTIF_FILES[$i]}"
    MOTIF_NAME=$(basename "$MOTIF_FILE" .cif)

    if [ -f "all_motifs_scored.csv" ] && grep -q "^${MOTIF_NAME}," "all_motifs_scored.csv" 2>/dev/null; then
        ((SKIPPED++))
        continue
    fi
    
    # Skip if already in this task's CSV
    if [ -f "$OUTPUT_CSV" ] && grep -q "^${MOTIF_NAME}," "$OUTPUT_CSV" 2>/dev/null; then
        continue
    fi
    
    PROGRESS=$((i - START_INDEX + 1))
    echo "[$PROGRESS/$MOTIFS_TO_PROCESS] $MOTIF_NAME"
    
    # Extract PDB ID
    PDB_ID=$(echo "$MOTIF_NAME" | grep -oE '[0-9][A-Z0-9]{3}' | tail -1)
    
    if [ -z "$PDB_ID" ]; then
        echo "  ✗ Could not extract PDB ID"
        ((FAILED++))
        continue
    fi
    
    # Parse CIF
    CHAIN=$(grep "^ATOM" "$MOTIF_FILE" | head -1 | awk '{print $6}')
    RES_NUMS=$(grep "^ATOM" "$MOTIF_FILE" | awk '{print $5}' | sort -n | uniq)
    START_RES=$(echo "$RES_NUMS" | head -1)
    END_RES=$(echo "$RES_NUMS" | tail -1)
    
    if [ -z "$CHAIN" ] || [ -z "$START_RES" ] || [ -z "$END_RES" ]; then
        echo "  ✗ Could not parse CIF"
        ((FAILED++))
        continue
    fi
    
    MOTIF_RESIDUES=$(grep "^ATOM" "$MOTIF_FILE" | awk -v chain="$CHAIN" '{
        printf "%s-%s-%s-\n", chain, $4, $5
    }' | sort -u | tr '\n' ',' | sed 's/,$//')
    
    # Run scoring
    python app.py --pdb_id "$PDB_ID" --motif "$START_RES" "$END_RES" \
        --chain "$CHAIN" --motif-residues "$MOTIF_RESIDUES" 2>&1 | tail -3
    
    EXIT_CODE=$?
    
    if [ -f "motif_report.json" ]; then
        python3 << EOF
import json, csv
from pathlib import Path

try:
    with open('motif_report.json', 'r') as f:
        data = json.load(f)
    
    motif_length = data.get('motif_num_nucleotides', 0)
    num_paired = data.get('num_paired_nucleotides', 0)
    pairing_pct = round((num_paired / motif_length * 100), 1) if motif_length > 0 else 0
    
    row = {
        'Motif_Name': '$MOTIF_NAME',
        'PDB_ID': data.get('pdb_id', 'N/A'),
        'Chain': data.get('motif_chain', 'N/A'),
        'Residue_Range': data.get('motif_range', 'N/A'),
        'Motif_Score': data.get('motif_score', 0),
        'Full_Structure_Score': data.get('full_structure_score', 0),
        'Score_Difference': data.get('score_difference', 0),
        'Motif_Length': motif_length,
        'Num_Paired_Nucleotides': num_paired,
        'Pairing_Percentage': pairing_pct,
        'Total_Base_Pairs': data.get('total_base_pairs', 0),
    }
    
    csv_file = '$OUTPUT_CSV'
    file_exists = Path(csv_file).exists()
    
    with open(csv_file, 'a', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=row.keys())
        if not file_exists:
            writer.writeheader()
        writer.writerow(row)
    print("  ✓")
except Exception as e:
    print(f"  ✗ {e}")
    exit(1)
EOF
        
        if [ $? -eq 0 ]; then
            ((SUCCESS++))
        else
            ((FAILED++))
        fi
        rm -f motif_report.json
    else
        echo "  ✗ No report"
        ((FAILED++))
    fi
    
    rm -f report.json motif_report.json
done

echo "Task $TASK_ID: Success=$SUCCESS, Failed=$FAILED, Skipped=$SKIPPED"