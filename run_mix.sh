for p in 50 55 60 65 70 75 80 85 90 95
do
for h in 5 10 30
do
    for v in {1..15}
    do
        python3 real_bam_sim.py --bam1 real_bam/in1.bam --bam2 real_bam/in2.bam  --out in1_in2_cov${h}_${p}_v${v}.bam --p ${p} --seed ${v} --h ${h}
        mv in1_in2_cov${h}_${p}_v${v}.bam real_bam_mix/in1_in2_cov${h}_${p}_v${v}.bam
    done
done
done