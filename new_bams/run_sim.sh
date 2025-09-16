#!/usr/bin/env bash
# Usage: ./run_pairs.sh pairs.csv

set -euo pipefail

PAIRS_TXT="$1"
folder=$(basename "$PAIRS_TXT" .samples.txt )
mkdir -p "$folder"
# Skip header and process each line
tail "$PAIRS_TXT" | while IFS=, read -r loc1 sample1 loc2 sample2
do
    
    bam1="${loc1}/mtdna_bams/${sample1}.final.bam"
    bam2="${loc2}/mtdna_bams/${sample2}.final.bam"

    pair_dir="${sample1}_${sample2}"
    mkdir -p "$folder"
    mkdir -p "$folder/$pair_dir"
    echo "$bam1 $bam2"

    echo "[INFO] Running simulation for ${loc1}_${sample1} and ${loc2}_${sample2} ..."
    for h in 5 10 30; do
        for p in 50 60 70 80 90 95; do
            for v in {0..14}; do
                echo "  Coverage depth: $h, Endogenous proportion: $p%, version: $v"
                python3 make_simulations.py \
                    --bam1 "$bam1" \
                    --bam2 "$bam2" \
                    --out "${folder}/${sample1}_${sample2}_cov${h}_p${p}_v${v}.bam" \
                    --h "$h" \
                    --p "$p" \
                    --seed "$v"
            done
        done
    done
done