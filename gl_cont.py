#!/usr/bin/env python3
# coding: utf-8


import sys
import argparse
import os
from collections import Counter
import pysam
import numpy as np
from tqdm import tqdm, trange
from scipy.special import binom
import scipy.stats as st
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from multiprocess import Pool
from MN import *
import subprocess
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import arviz
from statsmodels.tsa.stattools import acf as autocorr
from contamsim import generate
import logging
import subprocess


def command_run(command, show=False):
    output = subprocess.getoutput(command)
    if show ==True:
        print(output)
    logging.warning(f"Command:{command}")
    logging.warning(output)
# logging.basicConfig(format="{levelname}:{name}:{message}", style="{")

def get_mean_base_calling_error(bam_fname):
    bam = pysam.AlignmentFile(bam_fname, "rb")
    N_bases = 0
    Sum_of_errors = 0
    for readId, read in enumerate(bam.fetch('chrM')):
        if not read.is_mapped or 'D' in read.cigarstring or 'I' in read.cigarstring:
            continue
        
        q_scores = read.query_qualities
        for q in q_scores:
            Sum_of_errors += 10**(- q/10)
            N_bases += 1
    return Sum_of_errors/N_bases        
        
        

def get_base_err1(bam_fname, ref, aln_pos, same_set):
    bam = pysam.AlignmentFile(bam_fname, "rb")
    correct = 0
    incorrect = 0
    for readId, read in enumerate(bam.fetch('chrM')):
        
        if not read.is_mapped or 'D' in read.cigarstring or 'I' in read.cigarstring:
            continue
        
        seq = read.query_sequence
        pos = read.reference_start
        
        if read.cigartuples[0][0] == 4: #read is soft clipped
            left_trim = read.cigartuples[0][1]
            seq = seq[left_trim:]
                                
        
        offset = 0
        debug_str = ''

        for k in range(len(seq)):
            if aln_pos[pos+k] in same_set:
                if seq[k].upper() == ref[aln_pos[pos+k]]:
                    correct+=1
                else:
                    # print(pos, k, readId)
                    incorrect += 1
    return correct, incorrect, incorrect/(correct + incorrect)
                
        



def preprocess(ref_fname, contaminants_fname, bam_fname, chrom='chrM', consensus='samtools'):
    filepath = bam_fname
    basename = os.path.basename(bam_fname)[:-4]
    pysam.index(filepath)
    if not os.path.exists(f'Results/{basename}'):
        os.mkdir(f'Results/{basename}')
    base = f'Results/{basename}/{basename}_{chrom}'
    logging.basicConfig(
     filename=f"{base}.log",
     encoding="utf-8",
     filemode="w",
     format="{asctime} - {levelname} - {message}",
     style="{",
     datefmt="%Y-%m-%d %H:%M",
     level='DEBUG'
     )

    
    command_run(f"samtools view {filepath} {chrom} -o {base}.bam")
    
    if consensus == 'samtools':
        command_run(f'samtools consensus -m simple --min-MQ 13 --min-BQ 20 -c 0.51 -q -d 1 -H 0.001 --show-ins no  -o {base}.fa {base}.bam')
    elif consensus == "bcftools":
        command_run(f'bcftools mpileup -d 2000 -m 3 -C50 -q 30 -EQ 20 -f {ref_fname} {base}.bam | bcftools call -Oz -m --ploidy 1 > {base}.vcf.gz')
        command_run(f'bcftools index {base}.vcf.gz')
        command_run(f'bcftools consensus -f {ref_fname} {base}.vcf.gz > {base}.fa')
    elif consensus == 'contamix':
        command_run(f'bcftools mpileup -d 2000 -m 3 -C50 -q 30 -EQ 20 -f {ref_fname} {base}.bam | bcftools call -Oz -m --ploidy 1 > {base}.vcf.gz')
        command_run(f'perl CnsMaj3_1.pl -i {base}.vcf -o {base}.fa -l 16569 -cov 1 -diff 0.5 -idiff 0.5 -h chrM  -callindels no > {base}.cns')
    elif consensus == 'mpileup':
        command_run(f'samtools mpileup -f {ref_fname} {base}.bam > {base}.mpileup')
        command_run(f'python3 mpileup_reference.py --ref {ref_fname} --mpileup {base}.mpileup -o {base}.fa')
    else:
        print('consensus')
        command_run(f'cat {consensus} > {base}.fa') 
    command_run(f'cat {base}.fa > {base}_{consensus}.fa') 
        
    # command_run(f'python3 cns.py -i {base}.vcf -o {base}.fa -l 16569 --cov 1 --diff 0.5 --idiff 0.5 --header chrM  --callindels yes > {base}.cns', show=True)
    
    logging.warning('#CONSENSUS IS READY')
    
    command_run(f'cat {base}.fa {contaminants_fname} > {base}_genomes.fa');
    # gather consensus and possible contaminants together
    command_run(f'mafft {base}_genomes.fa >  {base}_aligned.fa') # do multiple alignment
    aligned_genomes = f'{base}_aligned.fa'
    logging.warning("#ALL GENOMES ARE READY")
    command_run(f'bwa index -a bwtsw {base}.fa') #indexing consensus
    command_run(f'samtools faidx {base}.fa')
    command_run(f'rm {base}.dict')
    command_run(f'samtools fastq {bam_fname} > {base}.fq')
    command_run(f"bwa mem {base}.fa {base}.fq | samtools sort -O BAM -o {base}_ra.sort.bam")    
    
    
    # total_reads = int(pysam.view('-c', bam_fname))
    # need_reads = 64000
    # proportion = need_reads / total_reads #for tests only!
    # print(proportion)

    # -s {123+proportion}
    
    # os.system(f'samtools view  -O BAM {base}_ra.sort.bam | samtools sort > {base}_ra.sort.rmdup.bam')
    
    pysam.index(f'{base}_ra.sort.bam');
    command_run(f'samtools calmd -Erb {base}_ra.sort.bam {base}.fa > {base}_ra.final.bam');
    bam_final = f'{base}_ra.final.bam'
    command_run(f'samtools index {base}_ra.final.bam')
    command_run(f'rm {base}_ra.sai')
    logging.warning("#BAM FILE IS READY")
    return bam_final, aligned_genomes


def make_genomes_arr(genomes_fname):
    genomes = list()
    for record in SeqIO.parse(genomes_fname, "fasta"):
        genomes.append(str(record.seq))
    genomes_arr = np.array([list(x) for x in genomes], dtype = 'S1')
    return genomes_arr


def get_same(genomes_arr):
    same_positions = []
    for i in range(genomes_arr.shape[1]):
        if len(np.unique(genomes_arr[:,i])) == 1:
            same_positions.append(i)
    return set(same_positions)


def get_base_err(bam_fname, same_dict, chrom='chrM'):
    samfile = pysam.AlignmentFile(bam_fname, "rb" )
    same_positions = list(same_dict.keys())
    correct = 0
    total = 0
    samfile = pysam.AlignmentFile(bam_fname, "rb" )
    for pileupcolumn in tqdm(samfile.pileup(chrom)):
        pos = pileupcolumn.pos
        if pos not in same_positions:
            continue

        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
                total += 1
                # query position is None if is_del or is_refskip is set.
                nbase =  pileupread.alignment.query_sequence[pileupread.query_position]
                if nbase == same[pos].decode('ascii').upper():
                    correct += 1
    base_err = 1 - correct/total
    samfile.close()
    return base_err


def get_aln_pos(reference):
    aln_coor = []
    for i in range(len(reference)):
        if reference[i] != '-':
            aln_coor.append(i)
            
    return np.asarray(aln_coor)


def do_mcmc(MC, base_err, n_iterations = 50000, output_file='', n_threads=8, model=0, show_each=10, show_time=False):
    p_list = list()
    N_ess = list()
    if output_file != '':
        res = open(output_file,'w')
    num_reads, num_genomes  = MC.shape
    logging.info(f'there is {MC.shape[0]} reads and {MC.shape[1]} genomes')
    p = np.random.dirichlet([1]*num_genomes)
    # p = np.array([1,0.0])
    # base_err=0.0000
    if show_time:
        iterations = trange(n_iterations, leave=True)
    else: iterations = range(n_iterations)
    for i in iterations:
        
        func = lambda x: get_Zi(MC, p, base_err, x)
        
        # Z = np.array(pool.map_async(func, range(num_reads)).get())
        Z = np.array([func(s) for s in range(num_reads) ], dtype=int)
        eta = get_eta(Z, num_genomes)
        if model == 0:
            p0 = np.random.beta(0.1 + eta[0],0.1+num_reads-eta[0])
            p_other = np.random.dirichlet(0.1+ eta[1:])
            p_other *= (1.-p0)

            p[0] = p0
            p[1:] = p_other
        else:
            p = np.random.dirichlet(1+eta)
        p_list.append(p)
        # N_ess.append(arviz.ess(np.array(p_list), method='tail'))
        if show_time:
            iterations.set_description(f'contamination is {1-p[0]}')
    return p_list

def run_glcont(ref_fname, genomes_fname, bam_fname, n_iterations=10000, chrom='chrM', model = 1, consensus='samtools', show_time=False):
    if chrom not in ['chrM', 'chrY']:
        raise Exception('Wrong chromosome')
    bam, genomes = preprocess(ref_fname, genomes_fname, bam_fname, consensus=consensus)
    genomes_arr = make_genomes_arr(genomes)
    same = get_same(genomes_arr)
    genomes0 = (''.join( np.array(genomes_arr, dtype = str)[0])).upper()
    aln_coords = get_aln_pos(genomes0)
    M, N, base_err = get_MN(genomes_arr, bam, aln_coords, same, chrom=chrom)
    # base_err = 0.0001
    # base_err = get_base_err1(bam, genomes0, aln_coords, same)[2]
    base_err = get_mean_base_calling_error(bam)
    # base_err = 0.000
    np.savetxt('M.txt', M)
    np.savetxt('N.txt', N)
    
    
    idx = np.where((M[:,0]!=M[:,1]))[0] # Так можно?
    M = M[idx]
    N = N[idx]
    # print(M[227])
    # 
    # idx = np.where((M[:,0]==1))[0]
    # M = M[idx]
    # N = N[idx]
    idx = np.where((M[:,0]>M[:,1]))[0]
    # print(len(idx))
    # print(len(idx)/len(M))
    # print(M.shape)
    # return 0
    
    logging.info(f"M<N: {(M<N).sum()}")
    # base_err=0.0001
    # print(base_err)
    MC = get_mc(M, N, base_err)
    # print(M)
    # print(MC)
    
    idx = np.where((MC[:,0]!=MC[:,1]))[0] # Так можно?
    MC = MC[idx]
    logging.info(f'{base_err=}')
    # print(f'{base_err=}')
    
    P = do_mcmc(MC, base_err,  n_iterations, n_threads=1, model=model, show_time=show_time)
    endo_p = np.array(P)[:,0]
    with open(f'{bam[:-4]}.results','w') as f:
        print(f'Mean,q-5,q-95', file=f)
        print(f'{np.mean(endo_p)},{np.quantile(endo_p, 0.05)}, {np.quantile(endo_p, 0.95)}', file=f)
    return np.array(P)


def neff(arr):
    n = len(arr)
    acf = autocorr(arr, nlags=n, fft=True)
    sums = 0
    for k in range(1, len(acf)):
        sums = sums + (n-k)*acf[k]/n

# else:
    # m
def main():
    if not len(sys.argv) > 1:
        ref_fname = 'rCRS.fa'
        bam_fname =  f'simulations/X2b5_I1c1a/X2b5_I1c1a_cov_30_80_v10.bam'
        genomes_fname = 'fasta/I1c1a.fasta'
        nIter = 10000
        chrom = 'chrM'
    else:
        parser = argparse.ArgumentParser(description='Supply reference fasta and bam file')
        parser.add_argument('ref',
                            help='reference fasta')
        parser.add_argument('bam',
                            help='bam file')
        parser.add_argument('contaminants',
                            help='fasta with all possible contaminants', default='contaminants.fa')
        parser.add_argument('nIter',
                            help='number of iteration mcmc', default=10000)
        parser.add_argument('chrom',
                            help='number of iteration mcmc', default='chrM')
        args = parser.parse_args()
        ref_fname = args.ref
        bam_fname = args.bam
        genomes_fname = args.contaminants
        nIter = int(args.nIter)
        chrom = args.chrom
    print("filename = ", bam_fname)
    P = run_glcont(ref_fname, genomes_fname, bam_fname, n_iterations=10000, chrom='chrM', model=0, consensus='mpileup', show_time=False)
    print("Proportion is: ", np.mean(P[1000::5,0]))
if __name__ == '__main__':
    # P = run_glcont(ref_fname, genomes_fname, bam_fname, n_iterations=10000, chrom='chrM', model=0, true_endo=None, show_time=False)
    # generate(0.55)
    # print(get_mean_base_calling_error('simulations/X2b5_I1c1a/X2b5_I1c1a_cov_30_80_v0.bam'))
    main()
