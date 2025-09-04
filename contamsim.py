#!/usr/bin/env python
# coding: utf-8

#simple contamination read create

import numpy as np
import os
import sys
import numpy as np
import argparse
import random


folder_path = "fasta"
files = [f for f in os.listdir(folder_path) if os.path.isfile(os.path.join(folder_path, f))]
if not len(sys.argv)> 1:
    fnames = ['fasta/U5a2e.fasta', 'fasta/X2b5.fasta']
    proportions = [0.6, 0.4]
    err_base = 0.001
    del_base = 0
    ins_base = 0.0
    coverage = 30
    output = 'simulated_data.bam'
else:
    parser = argparse.ArgumentParser(description='Call out my name. Adepti Xiao. I will be here')
    parser.add_argument('--fnames', nargs='*', help='genomes filenames')
    parser.add_argument('--proportions', nargs='*',type=float , help='proportions of contaminants')
    parser.add_argument('--coverage', default=20, type=float, help='coverage')
    parser.add_argument('--err_base',default=0.01, help='coverage')
    parser.add_argument('--del_base',default=0., help='coverage')
    parser.add_argument('--ins_base',default=0., help='coverage')
    parser.add_argument('--output',default='simulated_data.bam',type=str, help='name of final file')
    args = parser.parse_args()
    fnames = args.fnames
    proportions = args.proportions
    err_base = args.del_base
    del_base = args.del_base
    ins_base = args.ins_base
    coverage = args.coverage
    output = args.output

def make_fastq(genomes_fnames, proportion, err_base, del_base, ins_base, coverage, output):
    
    assert np.sum(proportion) == 1
    
    cat_str = 'samtools cat'
    os.system('rm contaminants.fa')
    os.system('touch contaminants.fa')
    for i, genome_file in enumerate(genomes_fnames):
        os.system(f'cp "{genome_file}" ./genome_{i}.fa')
        if i != 0:
            os.system(f'cat "{genome_file}" >> contaminants.fa')
            
        with open(f'genome_{i}.fa', "r") as f:
            lines = f.readlines()

# Change the first line (be sure to preserve a newline if you want one)
        lines[0] = '>chrM\n'

# Write modified lines back to the file
        with open(f'genome_{i}.fa', "w") as f:
            f.writelines(lines)
            
        os.system(f'''sed -i -e "1 s/.*/>chrM/" genome_{i}.fa''')
        fname =  f'genome_{i}.fa'
        command_line = \
        f'simlord\
        -rr {fname}\
        -pi {ins_base}\
        -pd {del_base}\
        -ps {err_base}\
        -fl 100\
        -c {proportion[i] * coverage}\
        genome_{i}'
        os.system(command_line)
        os.system(f'samtools view -b genome_{i}.sam > genome_{i}.bam')
        os.system(f'rm genome_{i}.sam genome_{i}.fa genome_{i}.fastq')
        cat_str += f' genome_{i}.bam'
    
    os.system(f'{cat_str} | samtools sort > {output}') 
    for i, genome_file in enumerate(genomes_fnames):
        os.system(f'rm genome_{i}.bam')
    print('finish')

def main():
    fnames = ['fasta/X2b5.fasta', 'fasta/I1c1a.fasta'] #['fasta/X2b5.fasta', 'fasta/X2b6.fasta']# ['X2b5.fasta', 'I1c1a.fasta']  #
    for i in [0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.99]:
        for cov in [5, 10, 30]:
            for v in [1,2,3]:
    # proportions=[0.9, 0.1]
    # fnames = ['fasta/F1c1a1.fasta', 'fasta/J1c12b.fasta']
        # make_fastq(fnames, proportions, err_base, del_base, ins_base, coverage, f'simulated2.bam')
                make_fastq(fnames, [i, 1-i], err_base, del_base, ins_base, cov, f'X2b5_I1c1a_cov_{cov}_{int(100*i)}_v{v}.bam')

    # endo_proportions = [0.5, 0.6, 0.75, 0.9]
    # folder_path = "fasta"
    # files = [f for f in os.listdir(folder_path) if os.path.isfile(os.path.join(folder_path, f))]

    # print(files)

    # with open('fasta/F1c1a1.fasta') as f:
    #     line1 = f.readlines()[1].strip().replace(' ','')
    #     print(len(line1))   
    # with open('contaminants.fa') as f:
    #     line2 = f.readlines()[1].strip().replace(' ','')
    #     print(len(line2))
    # with open('simulated2_chrM.fa') as f:
    #     line3 = f.read().replace('\n','').split('>chrM')[1]
    #     # print(line3)
    #     print(len(line3))

    # with open('simulated2_chrM_genomes.fa') as f:
    #     lines = f.read().replace('\n','').split('>chrM')[1].split('>X2b6')
    #     line4 = lines[0].strip().replace(' ','')
    #     line5 = lines[1].strip().replace(' ','')
    #     print(len(line4), len(line5))
    #     # print('-' in line5)

    # with open('simulated2_chrM_aligned.fa') as f:
    #     lines = f.read().replace('\n','').split('>chrM')[1].split('>X2b6')
    #     line4 = lines[0]
    #     line5 = lines[1]
    #     print(len(line4), len(line5))
    #     print('-' in line4)


    # with open('simulated2_chrM_genomes1.fa') as f:
    #     lines = f.read().replace('\n','').split('>chrM')
    #     line4 = lines[1]
    #     line5 = lines[2]
        
    # len(line4)
    # for i in range(len(line4)):
    #     if line4[i] != line5[i]:
    #         print(i, line4[i], line5[i])
    # len(line3)
    # c = 0
    # for i in range(len(line4)):
    #     if line3[i] != line4[i] and line3[i]!='n' and line4[i]!='n':
    #         c += 1
    # c
    # from Levenshtein import distance
    # distance(line3, line4)
    # distance(line4, line5)
    # len(line4)
    # len(line5)
    
# def make_mix(name1, name2, p=0.8):
    

    
def generate(d):
    
    while True:
        random_file1 = random.choice(files)
        random_file2 = random.choice(files)
        print(random_file1, random_file2)
        if random_file1 != random_file2 and '+' not in random_file1 and '+' not in random_file2 and random_file1[-5:]=='fasta' and random_file2[-5:]=='fasta':
            break
    fnames = [f'fasta/{random_file1}', f'fasta/{random_file2}']
    proportions=[d, 1-d]
    make_fastq(fnames, proportions, 0.001, 0, 0, 30, f'simulated.bam')
    return random_file1[:-6], random_file2[:-6]



def generate2(p1, p2, p3):
    
    while True:
        random_file1 = random.choice(files)
        random_file2 = random.choice(files)
        random_file3 = random.choice(files)
        print(random_file1, random_file2)
        if len(set([random_file3, random_file1, random_file2])) == 3 and '+' not in random_file1 and '+' not in random_file2 and '+' not in random_file3 and random_file1[-5:]=='fasta' and random_file2[-5:]=='fasta' and random_file3[-5:]=='fasta' :
            break
    fnames = [f'fasta/{random_file1}', f'fasta/{random_file2}', f'fasta/{random_file3}']
    proportions=[p1, p2, p3]
    make_fastq(fnames, proportions, 0.001, 0, 0, 30, f'simulated.bam')
    return random_file1[:-6], random_file2[:-6], random_file3[:-6]
    
# generate2(0.6, 0.2, 0.2)
# generate(0.5)


    

if __name__ == '__main__':
    main()
    
    
def generate_X2b5_X2b6(d, depth, name = "simulated.bam"):
    fnames = [f'fasta/X2b5.fasta', f'fasta/X2b6.fasta']
    proportions=[d, 1-d]
    make_fastq(fnames, proportions, 0.01, 0, 0, depth, name)
    
def generate_certain(d, depth, name1, name2, name = "simulated.bam"):
    fnames = [name1, name2]
    proportions=[d, 1-d]
    make_fastq(fnames, proportions, 0.01, 0, 0, depth, name)