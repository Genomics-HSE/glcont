for i in {1..10} do
    for j in {$i+1..10} do
        first_bam=$i
        second_bam=$i
        mkdir -p simulations/${first_bam}_${second_bam}
        for p in 50 60 70 80 90 95
        do
        for h in 5 10 30
        do
            for v in {0..14}
            do
                python3 real_bam_sim.py --bam1 real_bam/${first_bam}.bam --bam2 real_bam/${second_bam}.bam  --out simulations/${first_bam}_${second_bam}/${first_bam}_${second_bam}_cov${h}_${p}_v${v}.bam --p ${p} --seed ${v} --h ${h}
            done
        done
        done
    done
done