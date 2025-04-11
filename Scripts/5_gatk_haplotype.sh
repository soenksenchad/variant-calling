#!/bin/bash
# File: Scripts/5_gatk_haplotype.sh
#===============================================================================
# GATK HaplotypeCaller Script
# Runs GATK HaplotypeCaller on processed BAM files
#===============================================================================

# Source the conda environment
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate BFFcent

# Get the sample name
filename=$(basename "${INPUT_R1}")
sample_name="${filename%_1.fastq.gz}"
input_bam="${SAMTOOLS_DIR}/${sample_name}_samtools_pipe.bam"

echo "Starting GATK HaplotypeCaller for sample: ${sample_name}"
echo "Input BAM: ${input_bam}"

# Check input files exist
if [ ! -f "${input_bam}" ]; then
    echo "Error: Input BAM file not found: ${input_bam}"
    exit 1
fi

if [ ! -f "${input_bam}.bai" ]; then
    echo "Error: BAM index not found: ${input_bam}.bai"
    exit 1
fi

# Check for required reference files
dict_file="${REF_GENOME%.*}.dict"
if [ ! -f "${dict_file}" ]; then
    echo "Error: Reference dictionary not found: ${dict_file}"
    exit 1
fi

if [ ! -f "${REF_GENOME}.fai" ]; then
    echo "Error: Reference index not found: ${REF_GENOME}.fai"
    exit 1
fi

# Run GATK HaplotypeCaller
echo "Running GATK HaplotypeCaller..."
gatk --java-options "-Xmx110g -Djava.io.tmpdir=$SLURM_SCRATCH" HaplotypeCaller \
    -R "${REF_GENOME}" \
    -I "${input_bam}" \
    -O "${GATK_DIR}/${sample_name}.g.vcf.gz" \
    --native-pair-hmm-threads 28 \
    -ERC GVCF \
    -L 3  # Process only chromosome 1 (Comment must be on its own line)

if [ $? -ne 0 ]; then
    echo "Error: GATK HaplotypeCaller failed"
    exit 1
fi

echo "GATK HaplotypeCaller complete for sample: ${sample_name}"
