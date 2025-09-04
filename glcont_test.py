import numpy as np
from gl_cont import run_glcont
from contamsim import generate, make_fastq

from mt_tree_analisys import *
ref_fname = 'refchrm.fa'
bam_fname = 'simulated.bam'
genomes_fname = 'contaminants.fa'
nIter = 50000
chrom = 'chrM'
p = 0.8

make_fastq(['fasta/F4.fasta', 'fasta/H1ap.fasta'], [p, 1-p], 0.001, 0,0,30,'simulated.bam')

run_glcont(ref_fname, genomes_fname,  f'simulated.bam', n_iterations=nIter, chrom='chrM', model=0, true_endo='X2b5.fasta')