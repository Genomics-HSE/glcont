#!/usr/bin/env bash
# Usage: ./calc_stats.sh pairs.csv

# Create a function for the inner loop
process_combination() {
    local sample1="$1"
    local sample2="$2"
    local h="$3"
    local p="$4"
    local v="$5"
    local model="$6"
    local consensus="$7"
    
    local contaminant_file="real_bam/${sample2}.fa"
    local bam="real_bam/${sample1}_${sample2}/${sample1}_${sample2}_cov${h}_p${p}_v${v}.bam"
    
    sh glcont.sh rCRS.fa "$bam" -c "$contaminant_file" -s "$consensus" -m "$model"
}

export -f process_combination

# Generate all parameter combinations and feed to parallel
sample1="in1"
sample2="in2"
for h in 30; do
    for p in 50 60 70 80 90 95; do
        for v in {0..9}; do
            for model in 0 1; do
                for consensus in real_bam/${sample1}.fa mpileup; do
                    echo "$sample1 $sample2 $h $p $v $model $consensus"
                done
            done
        done
    done
done | parallel -j $(nproc) --colsep ' ' '
    process_combination {1} {2} {3} {4} {5} {6} {7}
'