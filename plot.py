#%%
import numpy as np
import matplotlib.pyplot as plt
#%%
def get_glcont(second_haplogroup, consensus, model = 0):
    y = np.load(f'X2b5_{second_haplogroup}_model_{model}_consensus_{consensus}_average.npy')
    return y

endo_proportions = [0.5, 0.6, 0.7, 0.8, 0.9, 0.95]#%%

# %%
haplogroups = ['X2b6', 'I1c1a']
endo_proportions = [0.5, 0.6, 0.7, 0.8, 0.9, 0.95]
depths = [5, 10, 30]
models = [0, 1]
def plot_compare(depth, second_haplogroup):
    plt.rcParams.update({'axes.titlesize': 'small'})
    depths = {5:0, 10:1, 30:2}
    plt.figure(figsize=(30, 10))
    ax = plt.subplot(1,4,1)
    for i, model in enumerate(models):
        for j, cons in enumerate(['mpileup']):
            y = get_glcont(second_haplogroup, cons, model)
            mean = np.mean(y[depths[depth]], axis=0)
            std = np.std(y[depths[depth]], axis=0)
            if 2*i+j == 0:
                plt.subplot(1, 4, 2*i+j+1)
            else:
                plt.subplot(1, 4, 2*i+j+1, sharey=ax)
            plt.plot([0.5, 0.95],[0.5 ,0.95])
            plt.plot(endo_proportions, mean)
            # plt.xticks(endo_proportions[::2])
            plt.fill_between(endo_proportions, mean-std, mean + std, alpha=0.1)
            plt.title(f'glcont_model{model}\nconsensus_{cons}')
            plt.xlabel("True")
            plt.ylabel('Predicted')
    plt.show()
    
# for d in [5, 10, 30]:
#     for h in ['X2b6', 'I1c1a']:
#         plot_compare(d, h)
            
    # %%
second_haplogroup = 'I1c1a'
depth = 5
plot_compare(depth, second_haplogroup)    
# %%
