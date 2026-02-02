#!/bin/bash
#SBATCH --job-name=cache_full_scores
#SBATCH --output=logs/cache_full_%A_%a.out
#SBATCH --error=logs/cache_full_%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --array=0-99

# Create logs directory
mkdir -p logs
mkdir -p full_structure_cache

# Get number of array tasks from SLURM
NUM_TASKS=${SLURM_ARRAY_TASK_COUNT:-100}

# Find all unique PDB IDs from motif files
echo "Finding unique PDB IDs from motifs..."
UNIQUE_PDB_IDS=($(python3 << 'PYEOF'
import glob
from pathlib import Path

motifs_dir = "motifs"
motif_files = glob.glob(f"{motifs_dir}/*.cif")
pdb_ids = set()

for motif_file in motif_files:
    parts = Path(motif_file).stem.split('-')
    for part in reversed(parts):
        if len(part) == 4 and part[0].isdigit() and part[1:].isalnum():
            pdb_ids.add(part)
            break

for pdb_id in sorted(pdb_ids):
    print(pdb_id)
PYEOF
))

TOTAL_PDB_IDS=${#UNIQUE_PDB_IDS[@]}

echo "========================================"
echo "Array Task: $SLURM_ARRAY_TASK_ID"
echo "Total unique PDB IDs: $TOTAL_PDB_IDS"
echo "Number of array tasks: $NUM_TASKS"
echo "========================================"

if [ $TOTAL_PDB_IDS -eq 0 ]; then
    echo "ERROR: No PDB IDs found!"
    exit 1
fi

# Calculate PDB IDs per task based on ACTUAL number of tasks
PDB_IDS_PER_TASK=$(( (TOTAL_PDB_IDS + NUM_TASKS - 1) / NUM_TASKS ))

echo "PDB IDs per task: $PDB_IDS_PER_TASK"

# Calculate which PDB IDs this task should process
START_INDEX=$((SLURM_ARRAY_TASK_ID * PDB_IDS_PER_TASK))
END_INDEX=$((START_INDEX + PDB_IDS_PER_TASK - 1))

# Ensure we don't go beyond the last PDB ID
if [ $END_INDEX -ge $TOTAL_PDB_IDS ]; then
    END_INDEX=$((TOTAL_PDB_IDS - 1))
fi

# Skip if no PDB IDs to process
if [ $START_INDEX -ge $TOTAL_PDB_IDS ]; then
    echo "No PDB IDs to process for task $SLURM_ARRAY_TASK_ID"
    exit 0
fi

SUCCESS=0
FAILED=0

echo "Processing PDB IDs $START_INDEX to $END_INDEX ($((END_INDEX - START_INDEX + 1)) PDB IDs)"

# Process all PDB IDs assigned to this task
for ((i=START_INDEX; i<=END_INDEX; i++)); do
    PDB_ID="${UNIQUE_PDB_IDS[$i]}"
    
    echo "Processing: $PDB_ID (PDB ID $((i + 1))/$TOTAL_PDB_IDS)"
    
    # Check if already cached
    CACHE_FILE="full_structure_cache/${PDB_ID}.json"
    if [ -f "$CACHE_FILE" ]; then
        echo "  → Already cached, skipping"
        ((SUCCESS++))
        continue
    fi
    
    # Cache full structure score
    python3 cache_full_structure_scores.py --pdb-id "$PDB_ID" --cache-dir full_structure_cache
    
    exit_code=$?
    if [ $exit_code -eq 0 ]; then
        echo "  ✓ Successfully cached $PDB_ID"
        ((SUCCESS++))
    else
        echo "  ✗ Failed to cache $PDB_ID (exit code: $exit_code)"
        ((FAILED++))
    fi
    echo "----------------------------------------"
done

echo "Task $SLURM_ARRAY_TASK_ID completed: $SUCCESS cached, $FAILED failed"
exit 0