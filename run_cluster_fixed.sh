#!/bin/bash
#SBATCH --job-name=rna_quality
#SBATCH --output=logs/rna_quality_%A_%a.out
#SBATCH --error=logs/rna_quality_%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --array=0-999

# Create logs directory
mkdir -p logs

# Get list of all PDB files - FIXED VERSION
BASEPAIRS_DIR="data/basepairs"

# Use find to reliably get all JSON files regardless of count
echo "Finding all JSON files in $BASEPAIRS_DIR..."
PDB_FILES=($(find "$BASEPAIRS_DIR" -maxdepth 1 -name "*.json" -type f | sort))
TOTAL_FILES=${#PDB_FILES[@]}

# Extract just filenames without path
for i in "${!PDB_FILES[@]}"; do
    PDB_FILES[$i]=$(basename "${PDB_FILES[$i]}")
done

echo "========================================"
echo "Array Task: $SLURM_ARRAY_TASK_ID"
echo "Total files found: $TOTAL_FILES"
echo "========================================"

if [ $TOTAL_FILES -eq 0 ]; then
    echo "ERROR: No JSON files found!"
    exit 1
fi

# Calculate files per task
FILES_PER_TASK=$(( (TOTAL_FILES + 999) / 1000 ))

# Calculate which files this task should process
START_INDEX=$((SLURM_ARRAY_TASK_ID * FILES_PER_TASK))
END_INDEX=$((START_INDEX + FILES_PER_TASK - 1))

# Ensure we don't go beyond the last file
if [ $END_INDEX -ge $TOTAL_FILES ]; then
    END_INDEX=$((TOTAL_FILES - 1))
fi

# Skip if no files to process (for the last few tasks if file count doesn't divide evenly)
if [ $START_INDEX -ge $TOTAL_FILES ]; then
    echo "No files to process for task $SLURM_ARRAY_TASK_ID"
    exit 0
fi

SUCCESS=0
WARNING=0
FAILED=0

echo "Processing files $START_INDEX to $END_INDEX ($((END_INDEX - START_INDEX + 1)) files)"

# Process all files assigned to this task
for ((i=START_INDEX; i<=END_INDEX; i++)); do
    PDB_FILE="${PDB_FILES[$i]}"
    PDB_ID="${PDB_FILE%.json}"
    FULL_PATH="${BASEPAIRS_DIR}/${PDB_FILE}"
    
    echo "Processing: $PDB_ID (file $((i + 1))/$TOTAL_FILES)"
    
    python app.py --pdb_id "$PDB_ID" --csv "scores_summary.csv"
    
    exit_code=$?
    if [ $exit_code -eq 0 ]; then
        echo "Successfully processed $PDB_ID"
        ((SUCCESS++))
    elif [ $exit_code -eq 1 ]; then
        echo "Processed $PDB_ID (poor quality)"
        ((WARNING++))
    else
        echo "Failed to process $PDB_ID (exit code: $exit_code)"
        ((FAILED++))
    fi
    echo "----------------------------------------"
done

echo "Task $SLURM_ARRAY_TASK_ID completed: $SUCCESS success, $WARNING warnings, $FAILED failed"
exit 0

