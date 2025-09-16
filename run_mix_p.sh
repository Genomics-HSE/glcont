#!/bin/bash

export -f
parallel --jobs 8 '
    i={1}
    j={2}
    first_bam=in$i
    second_bam=in$j
    mkdir -p simulations/${first_bam}_${second_bam}
    python3 real_bam_sim.py \
        --bam1 real_bam/${first_bam}.bam \
        --bam2 real_bam/${second_bam}.bam \
        --out simulations/${first_bam}_${second_bam}/${first_bam}_${second_bam}_cov{3}_{4}_v{5}.bam \
        --p {4} --seed {5} --h {3}
' ::: {1..10} ::: {2..10} ::: 5 10 30 ::: 50 60 70 80 90 95 ::: {0..14}
