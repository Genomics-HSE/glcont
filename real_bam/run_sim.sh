#!/usr/bin/env bash
    
bam1="in1.bam"
bam2="in2.bam"
sample1="in1"
sample2="in2"
pair_dir="in1_in2"
    mkdir -p "$pair_dir"
    echo "$bam1 $bam2"
    for h in 30; do
        for p in 50 60 70 80 90 95; do
            for v in {0..2}; do
                echo "  Coverage depth: $h, Endogenous proportion: $p%, version: $v"
                python3 make_simulations.py \
                    --bam1 "$bam1" \
                    --bam2 "$bam2" \
                    --out "${pair_dir}/${sample1}_${sample2}_cov${h}_p${p}_v${v}.bam" \
                    --h "$h" \
                    --p "$p" \
                    --seed "$v"
            done
        done
    done
done