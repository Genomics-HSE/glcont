#!/usr/bin/env bash
# Usage: ./calc_stats.sh pairs.csv

set -euo pipefail

PAIRS_TXT="$1"
folder=$(basename "$PAIRS_TXT" .samples.txt )

# Create a function for the inner loop
process_combination() {
    local sample1="$1"
    local sample2="$2"
    local h="$3"
    local p="$4"
    local v="$5"
    local model="$6"
    local consensus="$7"
    local folder="$8"
    
    local contaminant_file="new_bam_consensuses/${sample2}.final.fa"
    local bam="${folder}/${sample1}_${sample2}_cov${h}_p${p}_v${v}.bam"
    
    sh glcont.sh rCRS.fa "$bam" -c "$contaminant_file" -s "$consensus" -m "$model"
}

export -f process_combination

# Generate all parameter combinations and feed to parallel
tail "$PAIRS_TXT" | while IFS=, read -r loc1 sample1 loc2 sample2; do
    for h in 5 10 30; do
        for p in 50 60 70 80 90 95; do
            for v in {0..2}; do
                for model in 0 1; do
                    for consensus in new_bam_consensuses/${sample1}.final.fa mpileup; do
                        echo "$sample1 $sample2 $h $p $v $model $consensus $folder"
                    done
                done
            done
        done
    done
done | parallel -j $(nproc) --colsep ' ' '
    process_combination {1} {2} {3} {4} {5} {6} {7} {8}
'