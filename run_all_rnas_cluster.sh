#!/bin/bash
#SBATCH --job-name=cache_and_score_rnas
#SBATCH --output=logs/cache_and_score_%j.out
#SBATCH --error=logs/cache_and_score_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

# Simple pipeline for the cluster:
# 1) Cache nucleotide counts and validation metrics
# 2) Run fast scoring on all RNAs (uses the cache to skip network calls)

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT_DIR"
mkdir -p logs

echo "==> [$(date)] Step 1: cache_metadata_parallel.py"
python3 cache_metadata_parallel.py

echo "==> [$(date)] Step 2: run_all_rnas_fast.py"
python3 run_all_rnas_fast.py

echo "==> [$(date)] Done."

