#!/usr/bin/env bash

bam1="in1.bam"
bam2="in2.bam"
sample1="in1"
sample2="in2"
pair_dir="in1_in2"

mkdir -p "$pair_dir"

# Функция для запуска симуляции
run_simulation() {
    local h=$1
    local p=$2
    local v=$3
    
    echo "Coverage depth: $h, Endogenous proportion: $p%, version: $v"
    python3 make_simulations.py \
        --bam1 "$bam1" \
        --bam2 "$bam2" \
        --out "${pair_dir}/${sample1}_${sample2}_cov${h}_p${p}_v${v}.bam" \
        --h "$h" \
        --p "$p" \
        --seed "$v"
}

export -f run_simulation
export bam1 bam2 sample1 sample2 pair_dir

# Запускаем параллельно
parallel --jobs 16 run_simulation ::: 30 ::: 50 60 70 80 90 95 ::: 0 1 2 3 4 5 6 7 8 9