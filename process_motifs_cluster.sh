#!/bin/bash
#SBATCH --job-name=process_motifs
#SBATCH --output=logs/process_motifs_%A_%a.out
#SBATCH --error=logs/process_motifs_%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --array=0-999

# Create necessary directories
mkdir -p logs
mkdir -p reports
mkdir -p motif_csvs

# Get list of all motif CIF files
MOTIFS_DIR="unique_motifs"

echo "Finding all motif CIF files in $MOTIFS_DIR..."
MOTIF_FILES=($(find "$MOTIFS_DIR" -maxdepth 1 -name "*.cif" -type f | sort))
TOTAL_MOTIFS=${#MOTIF_FILES[@]}

echo "========================================"
echo "Array Task: $SLURM_ARRAY_TASK_ID"
echo "Total motifs found: $TOTAL_MOTIFS"
echo "========================================"

if [ $TOTAL_MOTIFS -eq 0 ]; then
    echo "ERROR: No motif CIF files found!"
    exit 1
fi

# Calculate motifs per task
MOTIFS_PER_TASK=$(( (TOTAL_MOTIFS + 999) / 1000 ))

# Calculate which motifs this task should process
START_INDEX=$((SLURM_ARRAY_TASK_ID * MOTIFS_PER_TASK))
END_INDEX=$((START_INDEX + MOTIFS_PER_TASK - 1))

# Ensure we don't go beyond the last motif
if [ $END_INDEX -ge $TOTAL_MOTIFS ]; then
    END_INDEX=$((TOTAL_MOTIFS - 1))
fi

# Skip if no motifs to process
if [ $START_INDEX -ge $TOTAL_MOTIFS ]; then
    echo "No motifs to process for task $SLURM_ARRAY_TASK_ID"
    exit 0
fi

SUCCESS=0
FAILED=0

echo "Processing motifs $START_INDEX to $END_INDEX ($((END_INDEX - START_INDEX + 1)) motifs)"

# Process all motifs assigned to this task
for ((i=START_INDEX; i<=END_INDEX; i++)); do
    MOTIF_FILE="${MOTIF_FILES[$i]}"
    MOTIF_NAME=$(basename "$MOTIF_FILE" .cif)
    
    echo "Processing: $MOTIF_NAME (motif $((i + 1))/$TOTAL_MOTIFS)"
    
    # Check if already processed
    REPORT_FILE="reports/${MOTIF_NAME}.json"
    CSV_FILE="motif_csvs/${MOTIF_NAME}.csv"
    
    if [ -f "$REPORT_FILE" ] && [ -f "$CSV_FILE" ]; then
        echo "  → Already processed, skipping"
        ((SUCCESS++))
        continue
    fi
    
    # Process motif (each writes to its own files - no race conditions)
    python3 app.py --motif-name "$MOTIF_NAME" \
        --output-dir reports \
        --csv-dir motif_csvs
    
    exit_code=$?
    if [ $exit_code -eq 0 ]; then
        echo "  ✓ Successfully processed $MOTIF_NAME"
        ((SUCCESS++))
    else
        echo "  ✗ Failed to process $MOTIF_NAME (exit code: $exit_code)"
        ((FAILED++))
    fi
    echo "----------------------------------------"
done

echo "Task $SLURM_ARRAY_TASK_ID completed: $SUCCESS success, $FAILED failed"
exit 0

