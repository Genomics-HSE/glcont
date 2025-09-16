#!/usr/bin/env bash
# Usage: JOBS=12 ./run_pairs_parallel.sh pairs.samples.txt
# Requires: GNU parallel

set -euo pipefail

command -v parallel >/dev/null 2>&1 || {
  echo "[ERROR] GNU parallel is required. Install it (e.g., 'sudo apt-get install parallel') and retry." >&2
  exit 1
}

PAIRS_TXT="$1"
JOBS="${JOBS:-8}"

folder=$(basename "$PAIRS_TXT" .samples.txt)
mkdir -p "$folder"

# Skip header row if present
tail -n +2 "$PAIRS_TXT" | while IFS=, read -r loc1 sample1 loc2 sample2; do
  # Build inputs
  bam1="${loc1}/mtdna_bams/${sample1}.final.bam"
  bam2="${loc2}/mtdna_bams/${sample2}.final.bam"

#   pair_dir="${folder}/${sample1}_${sample2}"
#   mkdir -p "$pair_dir"

  echo "[INFO] Running simulations for ${loc1}_${sample1} and ${loc2}_${sample2} ..."
#   echo "       BAM1: $bam1"
#   echo "       BAM2: $bam2"
#   echo "       OUT : $pair_dir"

  # Run all (h, p, v) in parallel
  parallel -j "$JOBS" --halt soon,fail=1 '
    h={1}; p={2}; v={3};
    echo "  [RUN] h=$h p=$p v=$v"
    python3 make_simulations.py \
      --bam1 "'"$bam1"'" \
      --bam2 "'"$bam2"'" \
      --out "'"$folder"'"/'"$sample1"'_'"$sample2"'_cov${h}_p${p}_v${v}.bam \
      --h "$h" \
      --p "$p" \
      --seed "$v"
  ' ::: 5 10 30 ::: 50 60 70 80 90 95 ::: $(seq 0 2)

  echo "[OK] Finished ${sample1}_${sample2}"
done
