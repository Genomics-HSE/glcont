import sys
import os

# Get the absolute path to the directory containing the module
# module_path = os.path.abspath('/Users/nikita/Desktop/HSE/genomic/glcont') 

# Add the path to sys.path
# sys.path.append(module_path)

from contamsim import generate_certain
import numpy as np
ref_fname = 'rCRS.fa'
bam_fname = 'simulated.bam'
genomes_fname = 'contaminants.fa'
nIter = 20000
second_haplogroup = 'I1c1a'
chrom = 'chrM'
nPoints = 15
endo_proportions = [0.5, 0.6, 0.7, 0.8, 0.9, 0.95]
for k, depth in enumerate([5, 10, 30]):
    for i in range(nPoints):
        for j, endo_proportion in enumerate(endo_proportions):
            generate_certain(endo_proportion, depth, 'fasta/X2b5.fasta', f'fasta/{second_haplogroup}.fasta', f'simulations/X2b5_{second_haplogroup}/X2b5_{second_haplogroup}_cov_{depth}_{int(100*endo_proportion)}_v{i}.bam')
            # os.replace(f'/Users/nikita/Desktop/HSE/genomic/glcont/X2b5_X2b6_cov_{depth}_{int(100*endo_proportion)}_v{i}.bam', f'X2b5_X2b6_cov_{depth}_{int(100*endo_proportion)}_v{i}.bam')
            
            
            