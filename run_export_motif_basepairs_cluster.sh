#!/bin/bash
#SBATCH --job-name=export_motif_bp
#SBATCH --output=logs/export_motif_bp_%A_%a.out
#SBATCH --error=logs/export_motif_bp_%A_%a.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --array=0-999

# Simple submission:
#   sbatch run_export_motif_basepairs_cluster.sh
#
# Each array task processes a chunk of PDB IDs sequentially and writes shards under motif_base_pair/<PDB_ID>/.
# After all tasks complete, merge shards:
#   python3 export_motif_basepairs.py --merge-shards --shard-dir motif_base_pair --output motif_basepairs.csv

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT_DIR"
mkdir -p logs

IFS=$'\n' read -r -d '' -a PDB_IDS < <(ls motifs/*.cif 2>/dev/null | sed -E 's/.*-([0-9A-Z]{4})-[^-]*$/\1/' | tr a-z A-Z | sort -u && printf '\0')
TOTAL=${#PDB_IDS[@]}

if [ "$TOTAL" -eq 0 ]; then
  echo "No motif CIF files found in motifs/; aborting."
  exit 1
fi

TASK_ID=${SLURM_ARRAY_TASK_ID:-0}

# Fixed one-PDB-per-task; excess array indices exit cleanly.
CHUNK_SIZE=1

START=$(( TASK_ID * CHUNK_SIZE ))
END=$(( START + CHUNK_SIZE - 1 ))

if [ "$START" -ge "$TOTAL" ]; then
  echo "Task ID $TASK_ID start index $START beyond total $TOTAL; exiting."
  exit 0
fi
if [ "$END" -ge "$TOTAL" ]; then
  END=$(( TOTAL - 1 ))
fi

echo "Array task $TASK_ID processing PDB indices [$START, $END] (CHUNK_SIZE=$CHUNK_SIZE of TOTAL=$TOTAL)"

for ((i=START; i<=END; i++)); do
  PDB_ID=${PDB_IDS[$i]}
  echo "  -> PDB_ID=$PDB_ID"
  python3 export_motif_basepairs.py \
    --motifs-dir motifs \
    --shard-by-pdb \
    --shard-dir motif_base_pair \
    --pdb "$PDB_ID"
done

echo "Task $TASK_ID done."

