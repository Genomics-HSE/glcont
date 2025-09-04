#!/usr/bin/env python3
# coding: utf-8
# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def legend_without_duplicate_labels(figure):
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    figure.figlegend(by_label.values(), by_label.keys(), loc='lower left')
plt.rcParams['text.usetex'] = True
arr = np.load('errors.npy')
endo_proportions = [0.5, 0.75, 0.9, 0.95]
# %%
for i in range(4):
    endo_proportion = endo_proportions[i]
    plt.subplot(2, 2, i+1)
    plt.title(f'$p_{{endo}} = {endo_proportion}$')
    plt.scatter(arr[1, :, i], np.abs(arr[0, :, i]), s=0.7)
    plt.plot([0,100], [endo_proportion, endo_proportion], '--', c='r', alpha=0.2)
    plt.xlabel('Distance')
    plt.ylabel('Estimated endogeneous proportion')
plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, 
                    top=0.9, wspace=0.4,hspace=0.4)
# %%
plt.boxplot((np.abs(arr[0, :, 0]), np.abs(arr[0, :, 1]), np.abs(arr[0, :, 2]),np.abs(arr[0, :, 3])), labels=[r'$p_{endo} = 0.5$', r'$p_{endo} = 0.75$', r'$p_{endo} = 0.9$', r'$p_{endo} = 0.95$'])
plt.ylabel('Estimated endogeneous proportion')
plt.xlabel('Real endogeneous proportion')
plt.show()


# %%
