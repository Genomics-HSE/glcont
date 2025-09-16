#%%
import numpy as np
import matplotlib.pyplot as plt
import Levenshtein

#%%
consensuses = []
for i in range(1, 11):
    consensuses.append(open(f'real_bam/in{i}.fa').read().replace('\n','')[5:])
# %%
M = np.zeros((10,10))
for i in range(10):
    for j in range(10):
        if i == j:
            M[i,j] = 0
        else:
            M[i,j] = Levenshtein.distance(consensuses[i], consensuses[j])
# %%
print(M)
# %%
np.savetxt('real_bam/consensus_distance_matrix.txt', M, fmt='%d',delimiter='\t')
