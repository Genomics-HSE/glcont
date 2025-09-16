#!/bin/bash



preprocess() {
    local ref_fname="$1"
    local contaminants_fname="$2"
    local bam_fname="$3"
    local chrom="${4:-chrM}"
    local consensus="${5:-samtools}"
    local nIter="${5:-5000}"
    local burn_in="${5:-5}"
    
    local filepath="$bam_fname"
    local basename=$(basename "$bam_fname" .bam)
    
    # Create directory if it doesn't exist
    if [ ! -d "Results/$basename" ]; then
        mkdir -p "Results/$basename"
    fi
    
    local base="Results/$basename/${basename}_${chrom}"
    
    # Create log file
    > "${base}.log"
    
    # Index the BAM file
    samtools index "$filepath"
    
    # Extract chromosome-specific BAM
    samtools view "$filepath" "$chrom" -o "${base}.bam"
    
    # Generate consensus based on method
    case "$consensus" in
        "samtools")
            samtools consensus -m simple --min-MQ 13 --min-BQ 20 -c 0.51 -q -d 1 -H 0.001 --show-ins no -o "${base}.fa" "${base}.bam"
            ;;
        "bcftools")
            bcftools mpileup -d 2000 -m 3 -C50 -q 30 -EQ 20 -f "$ref_fname" "${base}.bam" | bcftools call -Oz -m --ploidy 1 > "${base}.vcf.gz"
            bcftools index "${base}.vcf.gz"
            bcftools consensus -f "$ref_fname" "${base}.vcf.gz" > "${base}.fa"
            ;;
        "contamix")
            bcftools mpileup -d 2000 -m 3 -C50 -q 30 -EQ 20 -f "$ref_fname" "${base}.bam" | bcftools call -Oz -m --ploidy 1 > "${base}.vcf"
            perl CnsMaj3_1.pl -i "${base}.vcf" -o "${base}.fa" -l 16569 -cov 1 -diff 0.5 -idiff 0.5 -h chrM -callindels no > "${base}.cns"
            ;;
        "mpileup")
            samtools mpileup -f "$ref_fname" "${base}.bam" > "${base}.mpileup"
            python3 mpileup_reference.py --ref "$ref_fname" --mpileup "${base}.mpileup" -o "${base}.fa"
            ;;
        *)
            # If consensus is a file, copy it
            if [ -f "$consensus" ]; then
                cat "$consensus" > "${base}.fa"
            else
                echo "Error: Invalid consensus method or file not found: $consensus"
                exit 1
            fi
            ;;
    esac
    
    # Copy consensus file with consensus method in name
    cp "${base}.fa" "${base}_${consensus}.fa"
    
    echo "#CONSENSUS IS READY"
    
    # Combine consensus and contaminants
    cat "${base}.fa" "$contaminants_fname" > "${base}_genomes.fa"
    
    # Perform multiple alignment
    mafft "${base}_genomes.fa" > "${base}_aligned.fa"
    
    local aligned_genomes="${base}_aligned.fa"
    echo "#ALL GENOMES ARE READY"
    
    # Index consensus
    bwa index -a bwtsw "${base}.fa"
    samtools faidx "${base}.fa"
    
    # Remove any existing dict file
    rm -f "${base}.dict"
    
    # Convert BAM to FASTQ
    samtools fastq "$bam_fname" > "${base}.fq"
    
    # Align reads to consensus
    bwa mem "${base}.fa" "${base}.fq" | samtools sort -O BAM -o "${base}_ra.sort.bam"
    
    # Index and process final BAM
    samtools index "${base}_ra.sort.bam"
    samtools calmd -Erb "${base}_ra.sort.bam" "${base}.fa" > "${base}_ra.final.bam"
    
    local bam_final="${base}_ra.final.bam"
    samtools index "$bam_final"
    
    # Clean up temporary files
    rm -f "${base}_ra.sai"
    
    echo "#BAM FILE IS READY"
    # python3 gl_cont.py --help
   python3 gl_cont.py  --bam ${bam_final} --genomes ${aligned_genomes} --output "Results/$basename" # --nIter ${nIter} --chrom ${chrom} --model ${model} --burn_in ${burn_in}
    
    echo "$bam_final,$aligned_genomes"

}


#!/bin/bash

# Default values
CONTAMINANTS="contaminants.fa"
NITER=5000
CHROM="chrM"
MODEL=0
CONSENSUS=mpileup
BURN_IN=5

# Function to display usage
usage() {
    echo "Usage: $0 [options] <ref> <bam>"
    echo "Supply reference fasta and bam file"
    echo ""
    echo "Positional arguments:"
    echo "  ref          reference fasta"
    echo "  bam          bam file"
    echo ""
    echo "Optional arguments:"
    echo "  -c, --contaminants FILE   fasta with all possible contaminants (default: $CONTAMINANTS)"
    echo "  -n, --niter INTEGER       number of iteration mcmc (default: $NITER)"
    echo "  -C, --chrom CHROM         chromosome (default: $CHROM)"
    echo "  -m, --model INTEGER       model 0 or 1 (default: $MODEL)"
    echo "  -s, --consensus METHOD    consensus: contamix, samtools, mpileup or provided (default: $CONSENSUS)"
    echo "  -b, --burn-in INTEGER     burn_in period (default: $BURN_IN)"
    echo "  -h, --help                show this help message and exit"
    exit 1
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -c|--contaminants)
            CONTAMINANTS="$2"
            shift 2
            ;;
        -n|--niter)
            NITER="$2"
            shift 2
            ;;
        -C|--chrom)
            CHROM="$2"
            shift 2
            ;;
        -m|--model)
            MODEL="$2"
            shift 2
            ;;
        -s|--consensus)
            CONSENSUS="$2"
            shift 2
            ;;
        -b|--burn-in)
            BURN_IN="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        -*)
            echo "Unknown option: $1"
            usage
            ;;
        *)
            # Positional arguments
            if [[ -z "$REF" ]]; then
                REF="$1"
            elif [[ -z "$BAM" ]]; then
                BAM="$1"
            else
                echo "Too many positional arguments"
                usage
            fi
            shift
            ;;
    esac
done

# Check required arguments
if [[ -z "$REF" ]] || [[ -z "$BAM" ]]; then
    echo "Error: ref and bam arguments are required"
    usage
fi

# Assign to variables with similar names to Python version
ref_fname="$REF"
bam_fname="$BAM"
contaminants_fname="$CONTAMINANTS"
chrom="$CHROM"
nIter="$NITER"
model="$MODEL"
burn_in="$BURN_IN"
consensus="$CONSENSUS"

# Convert numeric values to integers (shell variables are strings by default)
nIter=$((nIter))
model=$((model))
burn_in=$((burn_in))

# Display the parsed values (for debugging)
echo "Reference: $ref_fname"
echo "BAM file: $bam_fname"
echo "Contaminants: $contaminants_fname"
echo "Chromosome: $chrom"
echo "Iterations: $nIter"
echo "Model: $model"
echo "Consensus method: $consensus"
echo "Burn-in: $burn_in"


preprocess ${ref_fname} ${contaminants_fname} ${bam_fname} ${chrom} $consensus $nIter

# Example of how to use these variables with the preprocess function from previous example
# preprocess "$ref_fname" "$contaminants_fname" "$bam_fname" "$chrom" "$consensus"



# Example usage:
# preprocess "reference.fa" "contaminants.fa" "input.bam" "chrM" "samtools"