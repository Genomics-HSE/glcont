import numpy as np
from gl_cont import run_glcont
from contamsim import generate

from mt_tree_analisys import *
ref_fname = 'refchrm.fa'
bam_fname = 'simulated.bam'
genomes_fname = 'contaminants.fa'
nIter = 10000
chrom = 'chrM'
nPoints = 50
y = np.zeros((nPoints, 4))
d = np.zeros((nPoints, 4))
endo_proportions = [0.5, 0.75, 0.9, 0.95]
tree = import_tree('mttree.json')

for i in range(nPoints):
    for j in range(4):
        endo_proportion = endo_proportions[j]
        hap1, hap2 = generate(endo_proportion)
        genomic_distance = calc_genomic_distance(hap1, hap2, tree)
        P = run_glcont(ref_fname, genomes_fname,  f'simulated.bam', n_iterations=nIter, chrom='chrM', model=1)
        y[i, j] = np.mean(P[:,0])
        print(P)
        d[i, j] = genomic_distance
        
np.save('errors', np.stack((y, d), axis=0))
print(y)