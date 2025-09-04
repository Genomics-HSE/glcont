import numpy as np
from gl_cont import run_glcont
from contamsim import generate2

ref_fname = 'refchrm.fa'
bam_fname = 'simulated.bam'
genomes_fname = 'contaminants.fa'
chrom = 'chrM'
nIter = 50000
n_point = 50
enodo_props = [0.6, 0.8]
it = 0
y = np.array([np.zeros((n_point, 3)) for j in range(2)])
while it < n_point:
    for i in range(2):
        enodo_prop = enodo_props[i]
        contam_prop = (1-enodo_prop) / 2
        hap1, hap2, hap3 = generate2(enodo_prop, contam_prop, contam_prop)
        if len(set([hap1, hap2, hap3])) != 3:
            continue
        y[i][it] = np.mean(
            run_glcont(
                ref_fname, 
                genomes_fname, 
                bam_fname, 
                n_iterations=nIter, 
                chrom=chrom, 
                model=1)[5000::4], axis=0)
    it+=1
np.save('results_ternary', y)