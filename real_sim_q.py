
import numpy as np
from gl_cont import run_glcont
from contamsim import generate_X2b5_X2b6
from mt_tree_analisys import *
import multiprocessing as mp

ref_fname = 'rCRS.fa'
bam_fname = 'real_bam/in3.bam'
genomes_fname = 'contaminants.fa'
nIter = 15000
chrom = 'chrM'
repeats = 15
endo_proportions = [0.5, 0.6, 0.7, 0.8, 0.9, 0.95]
tree = import_tree('mttree.json')
first_bam = 'in2'
second_bam = 'in3'
def calc_contamination(v, model, second_bam, consensus='samtools'):
    y_part = np.zeros((3, 6))  # Для каждого v создаем массив для 3 значений cov и 6 prop
    for k, cov in enumerate([5, 10, 30]):
    # for k, cov in enumerate([10, 30]):
        for j, prop in enumerate(endo_proportions):
            P = run_glcont(ref_fname, f'real_bam/{second_bam}.fa',  
                          f'simulations/{first_bam}_{second_bam}/{first_bam}_{second_bam}_cov{cov}_{int(100*prop)}_v{v}.bam', 
                          n_iterations=nIter, chrom='chrM', model=model, consensus=consensus)
            y_part[k, j] = np.mean(P[:,0])
            print(f'simulations/{first_bam}_{second_bam}/{first_bam}_{second_bam}_cov{cov}_{int(100*prop)}_v{v}.bam : {y_part[k, j]}')
    return y_part

def process_repeat(args):
    v, model, second_bam, consensus = args
    return calc_contamination(v, model, second_bam, consensus)


    
def calc_average_q(second_bam):
    with mp.Pool() as pool:
        for model in [0, 1]:
            # for consensus in ['samtools', 'bcftools', 'contamix', 'X2b5.fasta']:
            for consensus in [f'real_bam/{first_bam}.fa', 'mpileup']:
                # Подготовка аргументов
                args_list = [(v, model, second_bam, consensus) for v in range(repeats)]
                
                # Запускаем параллельное выполнение
                results = pool.map(process_repeat, args_list)
                
                # Собираем результаты
                y = np.stack(results, axis=1)  # Ось 1 соответствует repeats
                cons = consensus if consensus != f'real_bam/{first_bam}.fa' else 'true'
                np.save(f'{first_bam}_{second_bam}_model_{model}_consensus_{cons}_average', y)
                
if __name__ == '__main__':
    # calc_average_q('X2b6')  # Запуск функции для второго гаплогруппы
    calc_average_q('in3')  # Запуск функции для I1c1a