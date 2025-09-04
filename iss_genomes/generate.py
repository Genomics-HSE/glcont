import os
import numpy as np

proportions = [0.5, 0.6, 0.7, 0.8, 0.9, 0.95]
cov = [5, 10, 30]
num_repeats = 15

first_haplogroup = 'X2b5'
second_haplogroup = 'I1c1a'

def generate_data(first_haplogroup, second_haplogroup, proportions, cov, num_repeats):
    L = 16569  # Length of the mitochondrial genome
    for i in proportions:
        for c in cov:
            n1_reads = int(L)
            for v in range(num_repeats):
                os.system('')