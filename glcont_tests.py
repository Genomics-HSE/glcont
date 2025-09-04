import numpy as np
from gl_cont import run_glcont
from contamsim import generate, make_fastq
import csv

import numpy as np
import scipy.stats as st

# def confidence_interval(data, confidence=0.95):
    # """
    # Calculates the confidence interval for a given dataset.
# 
    # Args:
        # data (list or numpy.ndarray): The input data.
        # confidence (float): The desired confidence level (default is 0.95).
# 
    # Returns:
        # tuple: A tuple containing the lower and upper bounds of the confidence interval.
    # """
    # a = 1.0 * np.array(data)
    # n = len(a)
    # m, se = np.mean(a), st.sem(a)
    # h = se * st.t.ppf((1 + confidence) / 2., n-1)
    # return m-h, m+h


def x2b5_x2b6_check(true_cons = None):
    ref_fname = 'rCRS.fa'
    genomes_fname = 'simulations/X2b5_X2b6_simulations/contaminants.fa'
    
    iterations = 20000
    proportions = [0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.99]
    for v in [1,2,3]:
        for cov in [5, 10 ,30]:
            results = []
            results.append(['proportion, endogenous'])
            for proportion in proportions:
                bam_fname = f'simulations/X2b5_X2b6_simulations/X2b5_X2b6_cov_{cov}_{int(100*proportion)}_v{v}.bam'
                res = run_glcont(ref_fname, genomes_fname,  bam_fname, n_iterations=iterations, chrom='chrM', model=1, true_endo=true_cons)
                res = res[2000::2] # отбраковка 10%
                endo = np.mean(res[:, 0])
                sd = np.std(res[:, 0])
                ci = [np.quantile(res[:, 0], 0.025), np.quantile(res[:, 0], 0.975)]
                results.append([proportion, endo, sd, ci[0], ci[1]])
            if true_cons!=None:
                with open(f'simulations/X2b5_X2b6_simulations/result(fixcons)_cov_{cov}_v{v}.tsv', 'w', newline='') as tsvfile:
                    writer = csv.writer(tsvfile, delimiter='\t')
                    writer.writerows(results)
            else:
                with open(f'simulations/X2b5_X2b6_simulations/result_cov_{cov}_v{v}.tsv', 'w', newline='') as tsvfile:
                    writer = csv.writer(tsvfile, delimiter='\t')
                    writer.writerows(results) 
                
                
def x2b5_ici1a_check():
    ref_fname = 'refchrm.fa'
    genomes_fname = 'contaminants.fa'
    
    iterations = 10000
    proportions = [0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.99]
    # depth = 10
    for v in [2,3]:
        for cov in [5]:
            results = []
            results.append(['proportion, endogenous'])
            for proportion in proportions:
        
                bam_fname = f'X2b5_X2b6_simulations/X2b5_X2b6_cov_{cov}_{int(100*proportion)}_v{v}.bam'
                res = run_glcont(ref_fname, genomes_fname,  bam_fname, n_iterations=iterations, chrom='chrM', model=1)
                res = res[5000::4] # отбраковка 10%
                endo = np.mean(res[:, 0])
                sd = np.std(res[:, 0])
                ci = [np.quantile(res[:, 0], 0.025), np.quantile(res[:, 0], 0.975)]
                results.append([proportion, endo, sd, ci[0], ci[1]])
            with open(f'X2b5_X2b6_fixcont_contamix/result_cov_{cov}_v{v}.tsv', 'w', newline='') as tsvfile:
                writer = csv.writer(tsvfile, delimiter='\t')
                writer.writerows(results)
    
def check_1_5():
    ref_fname = 'refchrm.fa'
    genomes_fname = '5.fa'
    
    iterations = 10000
    proportions = [0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.99]
    # depth = 10
    for v in [1]:
        for cov in [10, 30]:
            results = []
            results.append(['proportion, endogenous'])
            for proportion in proportions:
        
                bam_fname = f'1-5_cov{cov}_{int(100*proportion)}_v{v}.bam'
                res = run_glcont(ref_fname, genomes_fname,  bam_fname, n_iterations=iterations, chrom='chrM', model=0)
                res = res[5000::4] # отбраковка 10%
                endo = np.mean(res[:, 0])
                sd = np.std(res[:, 0])
                ci = [np.quantile(res[:, 0], 0.025), np.quantile(res[:, 0], 0.975)]
                results.append([proportion, endo, sd, ci[0], ci[1]])
            with open(f'real_simulations/result(fixcons)_cov_{cov}_v{v}.tsv', 'w', newline='') as tsvfile:
                writer = csv.writer(tsvfile, delimiter='\t')
                writer.writerows(results)
    
x2b5_x2b6_check(true_cons=None)
# check_1_5()
# x2b5_ici1a_check()