#!/usr/bin/env bash
# Usage: ./calc_stats.sh pairs.csv

set -euo pipefail

PAIRS_TXT="$1"
folder=$(basename "$PAIRS_TXT" .samples.txt )

tail "$PAIRS_TXT" | while IFS=, read -r loc1 sample1 loc2 sample2
pair_dir="${sample1}_${sample2}"
for h in 5 10 30; do
    for p in 50 60 70 80 90 95; do
        for v in {0..14}; do
            for model in 0 1; do
                for consensus in all/consensuses/${sample2}.final.fa mpileup; do
                    bam=${folder}/${pair_dir}/${sample1}_${sample2}_cov${h}_p${p}_v${v}.bam
                    P=$(python3 ../glcont.py --ref ../rCRS.fa --bam "$bam" --consensus "$consensus" --model $model)
                    if [[ "$consensus" != "mpileup" ]]; then
                        consensus="True"
                    fi
                    echo "$sample1,$sample2,$h,$p,$v,$model,$consensus,$P"
                done
        done
    done
done
