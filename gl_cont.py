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
from MN import *
from Bio import SeqIO
# from contamsim import generate
import logging



def get_mean_base_calling_error(bam_fname):
    bam = pysam.AlignmentFile(bam_fname, "rb")
    N_bases = 0
    Sum_of_errors = 0
    for readId, read in enumerate(bam.fetch('chrM')):
<<<<<<< HEAD
        if  read.is_unmapped or 'D' in read.cigarstring or 'I' in read.cigarstring:
=======
        if  read.unmapped or 'D' in read.cigarstring or 'I' in read.cigarstring:
>>>>>>> 6c3c4c5 (pysam update to support old versions)
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
        
<<<<<<< HEAD
        if not read.is_unmapped or 'D' in read.cigarstring or 'I' in read.cigarstring:
=======
        if read.is_unmapped or 'D' in read.cigarstring or 'I' in read.cigarstring:
>>>>>>> 6c3c4c5 (pysam update to support old versions)
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



def run_glcont(bam_fname, genomes_fname, n_iterations=10000, output="Results", chrom='chrM', model = 1, consensus='samtools', show_time=False):
    if chrom not in ['chrM', 'chrY']:
        raise Exception('Wrong chromosome')
    print(bam_fname)
    genomes_arr = make_genomes_arr(genomes_fname)
    same = get_same(genomes_arr)
    genomes0 = (''.join( np.array(genomes_arr, dtype = str)[0])).upper()
    aln_coords = get_aln_pos(genomes0)
    M, N, base_err = get_MN(genomes_arr, bam_fname, aln_coords, same, chrom=chrom)
    base_err = get_mean_base_calling_error(bam_fname)
    basename = os.path.basename(f"{bam_fname}").split('.')[0]
    # print(f'!{basename}')
    # base = f'Results/{basename}/'
    # print(f'!{base}')
    np.savetxt(f'{output}/M.txt', M)
    np.savetxt(f'{output}/N.txt', N)


    idx = np.where((M[:,0]!=M[:,1]))[0] # Так можно?
    M = M[idx]
    N = N[idx]


    logging.info(f"M<N: {(M<N).sum()}")

    MC = get_mc(M, N, base_err)


    idx = np.where((MC[:,0]!=MC[:,1]))[0] # Так можно?
    MC = MC[idx]
    # logging.info(f'{base_err=}')


    P = do_mcmc(MC, base_err,  n_iterations, n_threads=1, model=model, show_time=show_time)
    endo_p = np.array(P)[:,0]
    with open(f'{bam_fname[:-4]}.results','w') as f:
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
    parser = argparse.ArgumentParser(description='Supply reference fasta and bam file')
    parser.add_argument('--bam',
                        help='bam file')
    parser.add_argument('--genomes',
                        help='fasta consensus + contaminants', default='contaminants.fa')
    parser.add_argument('--nIter',
                        help='number of iteration mcmc', default=5000)
    parser.add_argument('--chrom',
                        help='number of iteration mcmc', default='chrM')
    parser.add_argument('--model',
                        help='model 0 or 1', default=0)
    parser.add_argument('--burn_in',
                        help='burn_in period', default=5)
    parser.add_argument('--output',
                         help='where to put files', default='Results')
    args = parser.parse_args()
    # ref_fname = args.ref
    bam_fname = args.bam
    genomes_fname = args.genomes
    chrom = args.chrom
    nIter = int(args.nIter)
    chrom = args.chrom
    model = int(args.model)
    burn_in = int(args.burn_in)
    output = args.output
    P = run_glcont(bam_fname, genomes_fname, output=output, n_iterations=nIter, chrom=chrom, model=model, consensus='mpileup', show_time=False)
    endo_prop_mean = np.mean(P[int(nIter/10)::burn_in,0])
    print(endo_prop_mean)
    
if __name__ == '__main__':
    # P = run_glcont(ref_fname, genomes_fname, bam_fname, n_iterations=10000, chrom='chrM', model=0, true_endo=None, show_time=False)
    # generate(0.55)
    # print(get_mean_base_calling_error('simulations/X2b5_I1c1a/X2b5_I1c1a_cov_30_80_v0.bam'))
    main()
