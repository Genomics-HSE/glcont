#!/usr/bin/env python3
# coding: utf-8
# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import ternary

def legend_without_duplicate_labels(figure):
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    figure.figlegend(by_label.values(), by_label.keys(), loc='lower left')
plt.rcParams['text.usetex'] = True
arr = np.load('results_ternary.npy')
# print(arr4)
# %%

points_arr = [arr[0], arr[1]]
p = [0.6,0.8]
# 
# ternary.
plt.figure(figsize=(20, 10))
for i in range(2):
    
    ax = plt.subplot(1,2, i+1)
    figure, tax = ternary.figure(ax, scale=1)
    tax.set_title(f"Glcont two contaminants simulation with $p= ({p[i]}, {(1-p[i])/2:.2f}, {(1-p[i])/2:.2f}$) ")

    # Example data: list of (A, B, C) triples where A + B + C = 1
    points = points_arr[i]

    # Plot the points
    tax.scatter(points, marker='o', alpha=0.5, s=2, color='blue', label="Samples")
    tax.scatter([np.mean(points, axis=0)],marker='x', label="Mean prediction")
    tax.line([p[i], 1-p[i],0], [p[i], 0, 1-p[i]], color='green', linestyle=":", label="True proportion")
    # Beautify
    tax.boundary()
    tax.gridlines(multiple=0.1, color="gray")
    fontsize = 10
    tax.left_axis_label("Contaminant 2", offset=0.14, fontsize=fontsize)
    tax.right_axis_label("Contaminant 1", offset=0.14, fontsize=fontsize)
    tax.bottom_axis_label("Endogenous", offset=-0.04, fontsize=fontsize)
    # tax.legend()

    tax.ticks(axis='lbr', multiple=0.2, tick_formats="%.1f", fontsize=fontsize)
    tax.clear_matplotlib_ticks()

    # Show the plot
legend_without_duplicate_labels(plt)
plt.tight_layout()
tax.show()
# %%