from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pysam
import numpy as np

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


def get_base_err(bam_fname, same_dict):
    samfile = pysam.AlignmentFile(bam_fname, "rb" )
    same_positions = list(same_dict.keys())
    correct = 0
    total = 0
    samfile = pysam.AlignmentFile(bam_fname, "rb" )
    for pileupcolumn in tqdm(samfile.pileup("chrM")):
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


def get_num_indels(bam_fname, trunc = 7):
    samfile = pysam.AlignmentFile(bam_fname, "rb" )
    num_reads = 0
    for read in samfile.fetch('chrM'):
        if not read.is_mapped or read.pos < trunc:
            continue
        if "I" in read.cigarstring:
            num_reads += 1
    samfile.close()
    return num_reads


def get_cigar_string(bam_fname):
    ''''
    This function calculate mapped reads
    '''
    samfile = pysam.AlignmentFile(bam_fname, "rb" )
    
    for read in samfile.fetch('chrM'):
        print(read.cigartuples)
        samfile.close()
    
    
def calculate_likelihood(probs, mc):
    probs = np.asarray(probs)
    num_reads, num_genomes = mc.shape
    log_l = 0
    for i in range(num_reads):
        log_l += np.log((probs*MC[i,:]).sum())
    return log_l


def get_probs(mc, p):
    num_reads, num_genomes = mc.shape
    p = np.asarray(p)
    
    probs = np.zeros_like(mc)
    # probs = np.zeros(num_genomes, dtype = float)
    for i in range(num_reads):
        s = 0
        for j in range(num_genomes):
            probs[i, j] = mc[i, j] * p[j]
            s += probs[i, j]
        for j in range(num_genomes):
            probs[i, j] = probs[i, j] / s
    return probs


# def consensus_caller(ref_fname, bam_fname):
#     base = bam_fname[:-4]
#     os.system(f"samtools view {bam_fname} chrM -o {base+'_mt.bam'}")
#     base = base + '_mt'
#     os.system(f'samtools consensus -o {base}_st.fa {bam_fname}') #st means samtools
#     os.system(f'bwa index -a bwtsw {base}_st.fa') #indexing consensus
#     os.system(f'samtools faidx {base}_st.fa')
#     os.system(f'rm {base}.dict')
#     os.system(f'picard CreateSequenceDictionary R={base}.fa O={base}.dict')
#     os.system(f'samtools fastq {bam_fname} > {base}.fq')
#     os.system(f'bwa aln -l 1000 -t 10 {base}_st.fa {base}.fq > {base}_ra.sai')
#     os.system(f"bwa samse -r '@RG\\tID:{base}\\tLB:{base}_L1\\tPL:ILLUMINA\\tSM:{base}' {base}_st.fa {base}_ra.sai {base}.fq |  samtools sort -O BAM -o {base}_ra.sort.bam")
#     os.system(f'samtools index {base}_ra.sort.bam')
#     consensus = bam2consensus(f'{base}_st.fa', f'{base}_ra.sort.bam')
#     consensus_fa = '>chrM\n'+''.join(consensus) +'\n'
#     with open(f'{base}.fa', 'w') as new_genomes:
#         new_genomes.write(consensus_fa)
#     return consensus
    
    