#%%
import numpy as np
import matplotlib.pyplot as plt
#%%
def get_glcont(second_haplogroup, consensus, model = 0):
    y = np.load(f'glcont/X2b5_{second_haplogroup}_model_{model}_consensus_{consensus}_average.npy')
    return y

endo_proportions = [0.5, 0.6, 0.7, 0.8, 0.9, 0.95]#%%
def get_contamix(second_haplogroup, consensus = 0):
    cons = '' if consensus == 0 else '_cons'
    y_cont = np.zeros((3, 15, 6))
    for i, cov in enumerate([5, 10, 30]):
        for v in range(15):
            for k, prop in enumerate(endo_proportions):
                # print(i, v, k)
                with open(f'X2b5_{second_haplogroup}{cons}/X2b5_{second_haplogroup}_cov_{cov}_{int(100*prop)}_v{v}_mt.summary.txt') as f:
                    y_cont[i, v, k] = float(
                        f.read().splitlines()[10].split()[2])
    return y_cont
# %%
# y = get_glcont('X2b6', False, 0)
# y_cont = get_contamix('X2b6', 0)
# %%
haplogroups = ['X2b6', 'I1c1a']
consensuses = ['true', 'contamix', 'samtools', 'mpileup']
consensuses = ['true', 'mpileup']
endo_proportions = [0.5, 0.6, 0.7, 0.8, 0.9, 0.95]
depths = [5, 10, 30]
models = [0, 1]
def plot_compare(depth, second_haplogroup):
    plt.rcParams.update({'axes.titlesize': 'small'})
    depths = {5:0, 10:1, 30:2}
    plt.figure(figsize=(20, 8))
    ax = plt.subplot(2,len(models)*len(consensuses),1)
    for i, model in enumerate(models):
        for j, cons in enumerate(consensuses):
            y = get_glcont(second_haplogroup, cons, model)
            mean = np.mean(y[depths[depth]], axis=0)
            std = np.std(y[depths[depth]], axis=0)
            if 2*i+j == 0:
                plt.subplot(2, len(models)*len(consensuses), len(models)*i+j+1)
            else:
                plt.subplot(2, len(models)*len(consensuses), len(models)*i+j+1, sharey=ax)
            plt.plot([0.5, 0.95],[0.5 ,0.95])
            plt.plot(endo_proportions, mean)
            plt.xticks(endo_proportions[::2])
            plt.fill_between(endo_proportions, mean-std, mean + std, alpha=0.1)
            plt.title(f'glcont_model{model}\nconsensus_{cons}')
            plt.xlabel("True")
            plt.ylabel('Predicted')
    
    for j, cons in enumerate(['old_cons', 'true_cons']):
        y = get_contamix(second_haplogroup, j)
        mean = np.mean(y[depths[depth]], axis=0)
        std = np.std(y[depths[depth]], axis=0)
        plt.subplot(2, len(models)*len(consensuses), len(models)*len(consensuses)+j+1, sharey=ax)
        plt.plot([0.5, 0.95],[0.5 ,0.95])
        plt.plot(endo_proportions, mean)
        plt.xticks(endo_proportions[::2])
        plt.fill_between(endo_proportions, mean-std, mean + std, alpha=0.1)
        plt.title(f'contamix\nconsensus_{cons}')
        plt.xlabel("True")
        plt.ylabel('Predicted')
    plt.subplots_adjust(hspace=1, wspace=0.5)
    plt.suptitle(f"Haplogroup:{second_haplogroup}, depth:{depth}", fontsize=14)
    plt.savefig(f"modification_depths{depth}_{second_haplogroup}.pdf")
    # plt.show()
for d in [5, 10, 30]:
    for h in ['X2b6', 'I1c1a']:
        plot_compare(d, h)
            
    
    
