#!/bin/bash
# Batch process all RNA structures
# This is a simpler shell script alternative to run_all_rnas.py

# Configuration
BASEPAIRS_DIR="data/basepairs"
SCORES_CSV="scores_summary.csv"
DELAY=2  # seconds between runs
BATCH_SIZE=50  # longer break every N structures

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo "=========================================="
echo "BATCH RNA SCORING PROCESSOR (Shell Script)"
echo "=========================================="

# Get all PDB IDs
echo ""
echo "Scanning $BASEPAIRS_DIR for RNA structures..."
PDB_IDS=$(ls "$BASEPAIRS_DIR"/*.json 2>/dev/null | sed 's/.*\///; s/\.json$//' | tr '[:lower:]' '[:upper:]' | sort -u)

TOTAL=$(echo "$PDB_IDS" | wc -l | tr -d ' ')
echo "Found $TOTAL unique PDB IDs"

# Get already processed PDB IDs
if [ -f "$SCORES_CSV" ]; then
    PROCESSED=$(tail -n +2 "$SCORES_CSV" | cut -d',' -f1 | tr '[:lower:]' '[:upper:]' | sort -u | wc -l | tr -d ' ')
    echo "Found $PROCESSED already processed structures"
else
    PROCESSED=0
    echo "scores_summary.csv not found (will create new)"
fi

# Count remaining
REMAINING=$((TOTAL - PROCESSED))
echo ""
echo "Structures remaining to process: $REMAINING"

if [ "$REMAINING" -eq 0 ]; then
    echo -e "${GREEN}All structures have already been processed!${NC}"
    exit 0
fi

# Estimate time
ESTIMATED_MIN=$((REMAINING * (DELAY + 10) / 60))
echo "Estimated time: ~$ESTIMATED_MIN minutes"
echo "(Assuming ~10 seconds per structure + ${DELAY}s delay)"

# Ask for confirmation
echo ""
read -p "Proceed with processing $REMAINING structures? (yes/no): " response
if [ "$response" != "yes" ] && [ "$response" != "y" ]; then
    echo "Cancelled."
    exit 0
fi

# Process each PDB ID
SUCCESSFUL=0
FAILED=0
COUNT=0

for PDB_ID in $PDB_IDS; do
    # Check if already processed (simple check)
    if [ -f "$SCORES_CSV" ] && grep -qi "^$PDB_ID," "$SCORES_CSV"; then
        continue
    fi
    
    COUNT=$((COUNT + 1))
    echo ""
    echo "=========================================="
    echo "[$COUNT/$REMAINING] Processing: $PDB_ID"
    echo "=========================================="
    
    # Run scoring
    if python3 app.py --pdb_id "$PDB_ID" > /dev/null 2>&1; then
        echo -e "${GREEN}✓ Successfully processed $PDB_ID${NC}"
        SUCCESSFUL=$((SUCCESSFUL + 1))
    else
        echo -e "${RED}✗ Error processing $PDB_ID${NC}"
        FAILED=$((FAILED + 1))
    fi
    
    # Rate limiting
    if [ $COUNT -lt $REMAINING ]; then
        echo "Waiting $DELAY seconds before next run..."
        sleep $DELAY
    fi
    
    # Longer break every BATCH_SIZE structures
    if [ $((COUNT % BATCH_SIZE)) -eq 0 ] && [ $COUNT -lt $REMAINING ]; then
        echo ""
        echo "=========================================="
        echo "Processed $COUNT structures. Taking a longer break..."
        echo "Progress: $SUCCESSFUL successful, $FAILED failed"
        echo "=========================================="
        sleep 10
    fi
done

# Summary
echo ""
echo "=========================================="
echo "PROCESSING COMPLETE"
echo "=========================================="
echo "Total processed: $COUNT"
echo "Successful: $SUCCESSFUL"
echo "Failed: $FAILED"

