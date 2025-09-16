import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import numpy as np
from gl_cont import run_glcont
# from contamsim import generate_X2b5_X2b6
# from mt_tree_analisys import *
import multiprocessing as mp

ref_fname = '../rCRS.fa'
nIter = 15000
chrom = 'chrM'
repeats = 3
endo_proportions = [0.5, 0.6, 0.7, 0.8, 0.9, 0.95]
# tree = import_tree('mttree.json')
def calc_contamination(folder, sam1, sam2, model, consensus, v):
    y_part = np.zeros((3, 6))  # Для каждого v создаем массив для 3 значений cov и 6 prop
    for k, cov in enumerate([5, 10, 30]):
        for j, prop in enumerate(endo_proportions):
            bam_fname = f'{folder}/{sam1}.final_{sam2}.final/{sam1}.final_{sam2}.final_p{int(100*prop)}_cov{cov}_v{v}.bam'
            contaminant = f'all_consensuses/{sam2}.final.fa'
            P = run_glcont(ref_fname, contaminant,  
                          bam_fname, 
                          n_iterations=nIter,
                          chrom='chrM',
                          model=model,
                          consensus=consensus)
            y_part[k, j] = np.mean(P[:,0])
            # print({y_part[k, j]}')
    return y_part

def process_repeat(args):
    folder, sam1, sam2, model, consensus, v = args
    return calc_contamination(folder, sam1, sam2, model, consensus, v)


    
def calc_average_q(folder, sam1, sam2, consensus_file, contaminant):
    with mp.Pool() as pool:
        for model in [0, 1]:
            for consensus in [f'{consensus_file}', 'mpileup']:
                # Подготовка аргументов
                args_list = [(folder,  sam1, sam2, model, consensus, v) for v in range(repeats)]
                
                # Запускаем параллельное выполнение
                results = pool.map(process_repeat, args_list)
                
                # Собираем результаты
                y = np.stack(results, axis=1)  # Ось 1 соответствует repeats
                cons = consensus if consensus != f'{consensus_file}' else 'true'
                np.save(f'{sam1}_{sam2}_model_{model}_consensus_{cons}_average', y)
                
if __name__ == '__main__':
    # calc_average_q('X2b6')  # Запуск функции для второго гаплогруппы
    import sys
    folder = sys.argv[1]  # Первый аргумент - имя первого bam файла без расширения
    samples_file = f'{folder}.txt.samples.txt'
    with open(samples_file, 'r') as f:
        lines = f.readlines()
    for line in lines:
        line = line.strip()
        if not line:
            continue
        loc1, sam1, loc2, sam2 = line.split(',')
        consensus_file = f'all_consensuses/{sam1}.final.fa'
        contaminant = f'all_consensuses/{sam2}.final.fa'
        calc_average_q(folder, sam1, sam2, consensus_file, contaminant)