#!/bin/bash
#SBATCH --job-name=cache_scores
#SBATCH --output=logs/cache_scores_%A_%a.out
#SBATCH --error=logs/cache_scores_%A_%a.err
#SBATCH --time=04:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --array=0-99

# ==============================================================================
# Cache Full Structure Scores - SLURM Array Job
# ==============================================================================
# Pre-computes full structure scores for all unique PDB IDs.
# This avoids recomputing full structure scores for each motif.
#
# Usage:
#   sbatch cache_full_scores_cluster.sh
#
# After completion, run the main motif scoring job.
# ==============================================================================

# Create directories
mkdir -p logs
mkdir -p full_structure_cache

# Configuration
MOTIFS_DIR="motifs"
TASK_ID=$SLURM_ARRAY_TASK_ID
NUM_TASKS=$SLURM_ARRAY_TASK_COUNT

echo "========================================"
echo "Array Task: $TASK_ID / $((NUM_TASKS - 1))"
echo "========================================"

# Find all unique PDB IDs from motif files
echo "Finding unique PDB IDs from motifs..."
UNIQUE_PDB_FILE="unique_pdb_ids.txt"

if [ ! -f "$UNIQUE_PDB_FILE" ]; then
    echo "Generating unique PDB ID list..."
    find "$MOTIFS_DIR" -name "*.cif" -type f | while read cif_file; do
        motif_name=$(basename "$cif_file" .cif)
        # Extract PDB ID (4 chars: 1 digit + 3 alphanumeric)
        pdb_id=$(echo "$motif_name" | grep -oE '[0-9][A-Z0-9]{3}' | tail -1)
        if [ -n "$pdb_id" ]; then
            echo "$pdb_id"
        fi
    done | sort -u > "$UNIQUE_PDB_FILE"
    echo "Generated $UNIQUE_PDB_FILE"
fi

TOTAL_PDB_IDS=$(wc -l < "$UNIQUE_PDB_FILE")
echo "Total unique PDB IDs: $TOTAL_PDB_IDS"

# Calculate which PDB IDs this task should process
PDB_PER_TASK=$(( (TOTAL_PDB_IDS + NUM_TASKS - 1) / NUM_TASKS ))
START_LINE=$((TASK_ID * PDB_PER_TASK + 1))
END_LINE=$((START_LINE + PDB_PER_TASK - 1))

# Ensure we don't go beyond the last line
if [ $END_LINE -gt $TOTAL_PDB_IDS ]; then
    END_LINE=$TOTAL_PDB_IDS
fi

# Skip if no PDB IDs to process
if [ $START_LINE -gt $TOTAL_PDB_IDS ]; then
    echo "No PDB IDs to process for task $TASK_ID"
    exit 0
fi

PDB_IDS_TO_PROCESS=$((END_LINE - START_LINE + 1))
echo "This task will process: PDB IDs from line $START_LINE to $END_LINE ($PDB_IDS_TO_PROCESS PDB IDs)"
echo ""

SUCCESS=0
FAILED=0
SKIPPED=0

# Process PDB IDs assigned to this task
LINE_NUM=0
while IFS= read -r PDB_ID; do
    ((LINE_NUM++))
    
    # Skip lines before our start
    if [ $LINE_NUM -lt $START_LINE ]; then
        continue
    fi
    
    # Stop after our end
    if [ $LINE_NUM -gt $END_LINE ]; then
        break
    fi
    
    PROGRESS=$((LINE_NUM - START_LINE + 1))
    PCT=$(echo "scale=1; $PROGRESS * 100 / $PDB_IDS_TO_PROCESS" | bc)
    
    echo "[$PROGRESS/$PDB_IDS_TO_PROCESS] ($PCT%) $PDB_ID"
    
    # Check if already cached
    CACHE_FILE="full_structure_cache/${PDB_ID}.json"
    if [ -f "$CACHE_FILE" ]; then
        # Verify cache is valid
        if python3 -c "import json; json.load(open('$CACHE_FILE'))" 2>/dev/null; then
            echo "  ✓ Already cached"
            ((SKIPPED++))
            continue
        else
            echo "  ⚠ Invalid cache, recomputing..."
            rm -f "$CACHE_FILE"
        fi
    fi
    
    # Run app.py to compute full structure score
    python app.py --pdb_id "$PDB_ID" 2>&1 | tail -3
    
    exit_code=$?
    
    if [ $exit_code -eq 0 ]; then
        # Check if report.json was created and has valid data
        if [ -f "report.json" ]; then
            # Extract score and save to cache
            python3 -c "
import json
from pathlib import Path

try:
    with open('report.json', 'r') as f:
        data = json.load(f)
        full_score = data.get('overall_score')
        
        if full_score is not None:
            cache_data = {
                'pdb_id': '$PDB_ID',
                'full_structure_score': full_score,
                'full_structure_grade': data.get('grade', 'N/A'),
                'total_base_pairs': data.get('total_base_pairs', 0),
                'num_nucleotides': data.get('num_nucleotides', 0)
            }
            Path('full_structure_cache').mkdir(exist_ok=True)
            with open('full_structure_cache/$PDB_ID.json', 'w') as f:
                json.dump(cache_data, f, indent=2)
            print(f'  ✓ Cached: {full_score}/100')
        else:
            print('  ✗ No score in report')
            exit(1)
except Exception as e:
    print(f'  ✗ Error: {e}')
    exit(1)
"
            if [ $? -eq 0 ]; then
                ((SUCCESS++))
            else
                ((FAILED++))
            fi
        else
            echo "  ✗ No report.json generated"
            ((FAILED++))
        fi
    else
        echo "  ✗ Failed (exit code: $exit_code)"
        ((FAILED++))
    fi
    
    # Clean up
    rm -f report.json
    
done < "$UNIQUE_PDB_FILE"

echo ""
echo "========================================"
echo "Task $TASK_ID Summary"
echo "========================================"
echo "Processed: $((SUCCESS + FAILED + SKIPPED)) PDB IDs"
echo "Cached:    $SUCCESS"
echo "Skipped:   $SKIPPED (already cached)"
echo "Failed:    $FAILED"
echo "========================================"

exit 0

