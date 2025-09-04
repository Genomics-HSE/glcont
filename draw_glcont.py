#!/usr/bin/env python3
# coding: utf-8
# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import csv
from pprint import pprint
def draw_cons_vs_fix():
    for cov in [5]:
        for v in [1, 2,3]:
            results = []
            with open(f"simulations/X2b5_X2b6_simulations/result(fixcons)_cov_{cov}_v{v}.tsv") as fd:
                rd = csv.reader(fd, delimiter="\t", quotechar='"')
                next(rd)
                for row in rd:
                    results.append([float(x) for x in row])
                results = np.array(results)
                print(results)
            true_contamination = results[:, 0]
            # print('True contamination')
            # pprint(true_contamination)
            glcont_contamination_fix = results[:,1]
            # print('glcont_contamination')
            # pprint(glcont_contamination)

            si_left_glcont_fix = results[:,3]
            si_right_glcont_fix = results[:,4]
            results = []
            with open(f"simulations/X2b5_X2b6_simulations/result_cov_{cov}_v{v}.tsv") as fd:
                rd = csv.reader(fd, delimiter="\t", quotechar='"')
                next(rd)
                for row in rd:
                    results.append([float(x) for x in row])
                results = np.array(results)
                print(results)
            true_contamination = results[:, 0]
            glcont_contamination = results[:,1]

            si_left_glcont = results[:,3]
            si_right_glcont = results[:,4]
            # print(np.array([si_left, si_right]))

            contamix_contamination, si_left_contamix, si_right_contamix = [], [], []
            for i in [50, 55, 60, 65, 70, 75,80 , 85, 90, 95, 99]:
                with open(f'contamix/calc_consensus/X2b5_X2b6_cov_{cov}_{i}_v{v}_mt.summary.txt') as f:
                    lines = f.read().splitlines()
                    contamix_contamination.append(float(lines[10].split()[2]))
                    si_left_contamix.append(float(lines[9].split()[0]))
                    si_right_contamix.append(float(lines[9].split()[1]))
            
            contamix_contamination_fix, si_left_contamix_fix, si_right_contamix_fix = [], [], []
            for i in [50, 55, 60, 65, 70, 75,80 , 85, 90, 95, 99]:
                with open(f'contamix/fix_consensus/X2b5_X2b6_cov_{cov}_{i}_v{v}_mt.summary.txt') as f:
                    lines = f.read().splitlines()
                    contamix_contamination_fix.append(float(lines[10].split()[2]))
                    si_left_contamix_fix.append(float(lines[9].split()[0]))
                    si_right_contamix_fix.append(float(lines[9].split()[1]))        

            plt.figure(figsize=(20,20))
            plt.subplot(2,2,2)
            plt.plot([0.5, 1],[0.5, 1], label='True contamination', color='gray')
            sns.lineplot(x=true_contamination, y=glcont_contamination_fix, label="glcont")
            plt.fill_between(x=true_contamination, y1=si_left_glcont_fix, y2=si_right_glcont_fix, label='glcont ci', alpha=0.1)
            plt.legend()
            plt.xlabel("Simulated contamination proportion")
            plt.ylabel("Predicted contamination proportion")
            plt.title(f'coverage: {cov}, version: {v}, True_cons')
            
            plt.subplot(2,2,1)
            plt.plot([0.5, 1],[0.5, 1], label='True contamination', color='gray')
            sns.lineplot(x=true_contamination, y=glcont_contamination, label="glcont")
            plt.fill_between(x=true_contamination, y1=si_left_glcont, y2=si_right_glcont, label='glcont ci', alpha=0.1)
            plt.legend()
            plt.xlabel("Simulated contamination proportion")
            plt.ylabel("Predicted contamination proportion")
            plt.title(f'coverage: {cov}, version: {v}, Calc_cons')
            
            plt.subplot(2,2,3)
            plt.plot([0.5, 1],[0.5, 1], label='True contamination', color='gray')
            sns.lineplot(x=true_contamination, y=contamix_contamination, label="contamix")
            plt.fill_between(x=true_contamination, y1=si_left_contamix, y2=si_right_contamix, label='contamix ci', alpha=0.1)
            plt.legend()
            plt.xlabel("Simulated contamination proportion")
            plt.ylabel("Predicted contamination proportion")
            plt.title(f'coverage: {cov}, version: {v} Calc_cons')
            
            plt.subplot(2,2,4)
            plt.plot([0.5, 1],[0.5, 1], label='True contamination', color='gray')
            sns.lineplot(x=true_contamination, y=contamix_contamination_fix, label="contamix")
            plt.fill_between(x=true_contamination, y1=si_left_contamix_fix, y2=si_right_contamix_fix, label='contamix ci', alpha=0.1)
            plt.legend()
            plt.xlabel("Simulated contamination proportion")
            plt.ylabel("Predicted contamination proportion")
            plt.title(f'coverage: {cov}, version: {v}, fix_cons')
            
            plt.show()
        # %%
draw_cons_vs_fix()

# %%
