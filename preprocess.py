import os
import pysam



def preprocess(ref_fname, genomes_fname, bam_fname):
    
    base = bam_fname[:-4]
    
    pysam.index(bam_fname);
    
    os.system(f"samtools view {bam_fname} chrM -o {base+'_mt.bam'}")
    base = base + '_mt'
    
    print('#EXTRACTING MTDNA OK')
    
    # Один из вариантов получения консенсуса, работает не очень хорошо.
    '''consensus = bam2consensus(ref_fname, bam_fname)
    consensus_fa = '>chrM\n'+''.join(consensus) +'\n'
    with open(f'{base}.fa', 'w') as new_genomes:
        new_genomes.write(consensus_fa)'''
    
    
    # os.system(f'samtools consensus -o {base}.fa {bam_fname}')
    
    #os.system(f'sh gatkconsensus.sh {ref_fname} {base}.bam  {base}1.fa')
    #os.system(f'''awk -i '/^>/{{print ">chrM"; next}}{{print}}' {base}1.fa > {base}.fa''')
    #os.system(f'rm {base}1.fa')
    
    # os.system(f'cp genome_0.fa {base}.fa')
    #-C50
    os.system(f'bcftools mpileup  -d 2000 -m 3  -q 30 -EQ 20 -f {ref_fname} {base}.bam | bcftools call -m --ploidy 1 > {base}.vcf')

    os.system(f'perl CnsMaj3_1.pl -i {base}.vcf -o {base}.fa -l 16569 -cov 1 -diff 0.1 -idiff 0.1 -h {base} -callindels no > {base}.cns')
    
    # return 0
    
    print('#CONSENSUS IS READY')
    
    os.system(f'cat {base}.fa {genomes_fname} > {base}_genomes.fa'); # gather consensus and possible contaminants together
    os.system(f'mafft {base}_genomes.fa >  {base}_aligned.fa') # do multiple alignment
    aligned_genomes = f'{base}_aligned.fa'
    print("#ALL GENOMES ARE READY")
    # new_cons = list(SeqIO.parse(f'{base}_aligned.fa', "fasta"))[0]
    # SeqIO.write(new_cons, f'{base}.real.fa', "fasta") #gives you reference after it have been realigned with MAFFT
    os.system(f'bwa index -a bwtsw {base}.fa') #indexing consensus
    os.system(f'samtools faidx {base}.fa')
    os.system(f'rm {base}.dict')
    # os.system(f'picard CreateSequenceDictionary R={base}.fa O={base}.dict')
    os.system(f'samtools fastq {bam_fname} > {base}.fq')
    # os.system(f'bwa aln -l 1000 -t 10 {base}.fa {base}.fq > {base}_ra.sai')
    print(123)
    #-r '@RG\\tID:{base}\\tLB:{base}_L1\\tPL:ILLUMINA\\tSM:{base}'
    # os.system(f"bwa samse  {base}.fa {base}_ra.sai {base}.fq | samtools sort -O BAM -o {base}_ra.sort.bam")
    os.system(f"bwa mem {base}.fa {base}.fq | samtools sort -O BAM -o {base}_ra.sort.bam")
    # os.system(f'picard MarkDuplicates I={base}_ra.sort.bam O={base}_ra.sort.rmdup.bam METRICS_FILE=metrics.txt TMP_DIR=temp REMOVE_DUPLICATES=true ASSUME_SORTED=true ') #VALIDATION_STRINGENCY=LENIENT
    
    total_reads = int(pysam.view('-c', bam_fname))
    
    need_reads = 64000
    proportion = need_reads / total_reads
    # print(proportion)

    # -s {123+proportion}
    os.system(f'samtools view  -O BAM {base}_ra.sort.bam | samtools sort > {base}_ra.sort.rmdup.bam')
    
    pysam.index(f'{base}_ra.sort.rmdup.bam');
    os.system(f'samtools calmd -Erb {base}_ra.sort.rmdup.bam {base}.fa > {base}_ra.final.bam 2>/dev/null');
    bam_final = f'{base}_ra.final.bam'
    os.system(f'samtools index {base}_ra.final.bam')
    os.system(f'rm {base}_ra.sai')
    print("#BAM FILE IS READY")
    return bam_final, aligned_genomes