#!/usr/bin/env bash
# Usage: ./calc_stats.sh pairs.csv

set -euo pipefail

PAIRS_TXT="$1"
folder=$(basename "$PAIRS_TXT" .samples.txt )
tail "$PAIRS_TXT" | while IFS=, read -r loc1 sample1 loc2 sample2
do
    contaminant_file="all_consensuses/${sample2}.final.fa"
    for h in 5 10 30; do
        for p in 50 60 70 80 90 95; do
            for v in {0..2}; do
                for model in 0 1; do
                    for consensus in all_consensuses/${sample1}.final.fa mpileup; do
                        bam=${folder}/${sample1}_${sample2}_cov${h}_p${p}_v${v}.bam
                        
                        sh glcont.sh ../rCRS.fa "$bam" -c "$contaminant_file" -s "$consensus" -m $model
                        exit 0
                        # if [[ "$consensus" != "mpileup" ]]; then
    #                         consensus="True"
    #                     fi
    #                     echo "$sample1,$sample2,$h,$p,$v,$model,$consensus,$P"
                    done
                done
            done
        done
    done
done
