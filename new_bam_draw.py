import numpy
import matplotlib.pyplot as p_list
import sys
import pandas as pd
pair = sys.argv[1]
import os
from tqdm import tqdm
samples1 = []
samples2 = []
consensuses = []
models = []
covs = []
vs = []
means = []
ps = []
with open(pair) as f:
    lines = f.readlines()
    for line in tqdm(lines):
        loc1, sample1, loc2, sample2 = line.strip().split(',')
        for consensus in ['mpileup', 'true']:
            for p in [50, 60, 70, 80, 90, 95]:
                for cov in [5, 10, 30]:
                    for model in [0, 1]:
                        for v in range(3):
                            file_path=f'Results/{sample1}_{sample2}_cov{cov}_p{p}_v{v}_{consensus}_model{model}/{sample1}_{sample2}_cov{cov}_p{p}_v{v}_chrM_ra.final.results'
                            if os.path.exists(file_path):
                                with open(file_path) as f1:
                                    f1.readline()
                                    mean = float(f1.readline().split(',')[0])
                                samples1.append(sample1)
                                samples2.append(sample2)
                                covs.append(cov)
                                vs.append(v)
                                consensuses.append(consensus)
                                models.append(model)
                                ps.append(p)
                                means.append(mean)
                            
print(1)
df = pd.DataFrame({
    'sample1' : samples1,
    'sample2' : samples2,
    'model' : models,
    'proportion' : ps,
    'version' : vs,
    'coverage' : covs,
    "consensus" : consensuses,
    'mean' : means 
})     
print(2)
df.to_csv('1.txt')                 
                            
                    
                    