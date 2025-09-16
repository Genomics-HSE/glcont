#!/usr/bin/env python3

"""
# José Víctor Moreno Mayar <morenomayar@gmail.com> #
Переписанны на Python 3
оригинальный скрипт cns_maj.pl
"""

import sys
import re
import argparse
from typing import List, Dict, Tuple

"""
n means covered but tied under mindiff
ATGC=high cov, high Q, maj
atgc=high cov, low Q
N=low/no cov
n=no maj/alt low Q
if cov=4, cov 4 IS taken into account
ONLY take into account HIGHQ and HIGHCOV INDELs
"""

def print_usage():
    print("""
USAGE: {} <arguments>

Arguments:
    -i, --input       STR     vcf obtained with samtools mpileup -g and bcftools -cg
    -o, --output      STR     Output fasta file
    -l, --length      INT     Length of the reference that was used to map the reads (not multifasta)
    --cov             INT     Minimum depth of coverage required to call a cns base
    --diff            FLOAT   Minimum fraction of concordant read bases to call a cns base
    --idiff           FLOAT   Minimum fraction of concordant read bases to call a cns indel
    -h, --header      STR     Fasta header for the output file
    --callindels      STR     If set to "yes", indels are called. Other values turn this option off.

Sample run: ./cns_maj.py -i sample.vcf -o sample.fasta -l 16569 --cov 5 --diff 0.7 --idiff 0.7 -h sample

Note: All of the arguments are required.
""".format(sys.argv[0]))
    sys.exit(0)

def main():
    parser = argparse.ArgumentParser(description='Consensus sequence caller')
    parser.add_argument('-i', '--input', required=True, help='Input VCF file')
    parser.add_argument('-o', '--output', required=True, help='Output FASTA file')
    parser.add_argument('-l', '--length', type=int, required=True, help='Reference length')
    parser.add_argument('--cov', type=int, required=True, help='Minimum coverage')
    parser.add_argument('--diff', type=float, required=True, help='Minimum SNP fraction difference')
    parser.add_argument('--idiff', type=float, required=True, help='Minimum INDEL fraction difference')
    parser.add_argument('--header', required=True, help='FASTA header')
    parser.add_argument('--callindels', required=True, choices=['yes', 'no'], help='Call indels')
    
    args = parser.parse_args()
    
    min_cov = args.cov
    min_diff = args.diff
    indel_min_diff = args.idiff
    genome_size = args.length
    call_indels = args.callindels == 'yes'
    
    print(f"Min_depth={min_cov}")
    print(f"Min_diff={min_diff}")
    print(f"Indel_min_diff={indel_min_diff}")
    print(f"Initial_size={genome_size}")
    print("Used the following SNPs:")
    
    cns = ['N'] * genome_size
    indels = []
    het_indels = []
    
    with open(args.input, 'r') as vcf:
        for line in vcf:
            if line.startswith('#'):
                continue
                
            is_indel = False
            is_het = False
            cns_base = 'N'
            diff = 0.0
            depth = 0
            high_quality = False
            
            fields = line.strip().split('\t')
            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4]
            
            if ',' in alt:
                is_het = True
                
            info = fields[7]
            info_fields = info.split(';')
            
            dp4_match = None
            dp_match = None
            for field in info_fields:
                if field.startswith('DP4='):
                    dp4_match = field
                elif field.startswith('DP='):
                    dp_match = field
            
            ref_alleles = 0
            alt_alleles = 0
            total_depth = 0
            
            if dp4_match:
                high_quality = True
                dp4_values = list(map(int, dp4_match[4:].split(',')))
                ref_alleles = dp4_values[0] + dp4_values[1]
                alt_alleles = dp4_values[2] + dp4_values[3]
                total_depth = ref_alleles + alt_alleles
            elif dp_match:
                total_depth = int(dp_match[3:])
                ref_alleles = total_depth
            
            if info_fields[0] == "INDEL":
                is_indel = True
            
            # Handle SNPs
            if not is_het and not is_indel:
                if total_depth >= min_cov and high_quality:
                    if ref_alleles / total_depth > min_diff:
                        diff = ref_alleles / total_depth
                        cns_base = ref
                    else:
                        if alt_alleles / total_depth > min_diff and alt != '.':
                            diff = alt_alleles / total_depth
                            cns_base = alt
                            print(f"{pos}\t{ref}\t{alt}\t{total_depth}\t{ref_alleles}\t{alt_alleles}\t{diff}\t{cns_base}")
                        else:
                            cns_base = 'n'
                else:
                    if total_depth >= min_cov and not high_quality:
                        cns_base = ref.lower()
                
                cns[pos-1] = cns_base
            
            # Handle multiallelic SNPs
            if is_het and not is_indel:
                if ',' in alt:
                    alt = alt.split(',')[0]
                
                if total_depth >= min_cov and high_quality:
                    if ref_alleles / total_depth > min_diff:
                        diff = ref_alleles / total_depth
                        cns_base = ref
                    else:
                        if alt_alleles / total_depth > min_diff and alt != '.':
                            diff = alt_alleles / total_depth
                            cns_base = alt
                            print(f"{pos}\t{ref}\t{alt}\t{total_depth}\t{ref_alleles}\t{alt_alleles}\t{diff}\t{cns_base}\tMultiAllelic")
                        else:
                            cns_base = 'n'
                else:
                    if total_depth >= min_cov and not high_quality:
                        cns_base = ref.lower()
                
                cns[pos-1] = cns_base
            
            # Collect INDELs for later processing
            if is_indel and call_indels:
                indels.append(line)
    
    original_cns = cns.copy()
    
    if call_indels:
        print("Used the following INDELs:")
        shift = 0
        subst_shift = 0
        
        for indel in indels:
            depth = 0
            is_het = False
            ref_selected = False
            alt_selected = False
            high_quality = False
            
            fields = indel.strip().split('\t')
            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4]
            
            if ',' in alt:
                is_het = True
                
            info = fields[7]
            info_fields = info.split(';')
            
            dp4_match = None
            dp_match = None
            for field in info_fields:
                if field.startswith('DP4='):
                    dp4_match = field
                elif field.startswith('DP='):
                    dp_match = field
            
            ref_alleles = 0
            alt_alleles = 0
            total_depth = 0
            
            if dp4_match:
                high_quality = True
                dp4_values = list(map(int, dp4_match[4:].split(',')))
                ref_alleles = dp4_values[0] + dp4_values[1]
                alt_alleles = dp4_values[2] + dp4_values[3]
                total_depth = ref_alleles + alt_alleles
            elif dp_match:
                total_depth = int(dp_match[3:])
                ref_alleles = total_depth
            
            # Handle non-het INDELs
            if not is_het and alt != '.':
                if total_depth >= min_cov and high_quality:
                    if ref_alleles / total_depth > indel_min_diff:
                        diff = ref_alleles / total_depth
                        ref_selected = True
                    else:
                        if alt_alleles / total_depth > indel_min_diff:
                            diff = alt_alleles / total_depth
                            alt_selected = True
                else:
                    ref_selected = True
                
                if alt_selected:
                    ref_len = len(ref)
                    alt_len = len(alt)
                    shift = ref_len - alt_len
                    new_size = genome_size - shift
                    modified_cns = []
                    
                    # Copy everything before the INDEL
                    for i in range(pos - 1 + subst_shift):
                        modified_cns.append(cns[i])
                    
                    # Insert the ALT sequence
                    for base in alt:
                        modified_cns.append(base)
                    
                    # Copy everything after the INDEL with shift adjustment
                    for i in range(pos - 1 + subst_shift + ref_len, genome_size):
                        modified_cns.append(original_cns[i - subst_shift])
                    
                    print(f"{pos}\t{ref}\t{alt}\t{total_depth}\t{ref_alleles}\t{alt_alleles}\t{diff}\t{genome_size}\t{new_size}\t{ref_len}\t{alt_len}\t{alt}")
                    
                    subst_shift -= shift
                    cns = modified_cns
                    genome_size = new_size
            
            # Handle het INDELs
            if is_het and alt != '.':
                if ',' in alt:
                    alt = alt.split(',')[0]
                
                if total_depth >= min_cov and high_quality:
                    if ref_alleles / total_depth > indel_min_diff:
                        diff = ref_alleles / total_depth
                        ref_selected = True
                    else:
                        if alt_alleles / total_depth > indel_min_diff:
                            diff = alt_alleles / total_depth
                            alt_selected = True
                else:
                    ref_selected = True
                
                if alt_selected:
                    ref_len = len(ref)
                    alt_len = len(alt)
                    shift = ref_len - alt_len
                    new_size = genome_size - shift
                    modified_cns = []
                    
                    # Copy everything before the INDEL
                    for i in range(pos - 1 + subst_shift):
                        modified_cns.append(cns[i])
                    
                    # Insert the ALT sequence
                    for base in alt:
                        modified_cns.append(base)
                    
                    # Copy everything after the INDEL with shift adjustment
                    for i in range(pos - 1 + subst_shift + ref_len, genome_size):
                        modified_cns.append(original_cns[i - subst_shift])
                    
                    print(f"{pos}\t{ref}\t{alt}\t{total_depth}\t{ref_alleles}\t{alt_alleles}\t{diff}\t{genome_size}\t{new_size}\t{ref_len}\t{alt_len}\t{alt}\tMultiallelic")
                    
                    subst_shift -= shift
                    cns = modified_cns
                    genome_size = new_size
    
    print(f"The final size is: {genome_size}")
    
    with open(args.output, 'w') as out:
        out.write(f">{args.header}\n")
        out.write(''.join(cns))
        out.write("\n")

if __name__ == "__main__":
    main()