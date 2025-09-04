#!/usr/bin/env python3
# coding: utf-8

import os
import numpy as np

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