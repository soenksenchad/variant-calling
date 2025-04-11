#!/bin/bash
# File: Scripts/3_alignment.sh
#===============================================================================
# Alignment Script
# Aligns trimmed reads to reference genome using BWA-MEM2
#===============================================================================

# Source the conda environment
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate BFFcent

# Get the sample name
filename=$(basename "${INPUT_R1}")
sample_name="${filename%_1.fastq.gz}"

echo "Starting alignment for sample: ${sample_name}"

# Define input files
TRIMMED_R1="${FASTP_DIR}/${sample_name}_1_trimmed.fastq.gz"
TRIMMED_R2="${FASTP_DIR}/${sample_name}_2_trimmed.fastq.gz"

# Check if trimmed files exist
if [ ! -f "${TRIMMED_R1}" ] || [ ! -f "${TRIMMED_R2}" ]; then
    echo "Error: Trimmed files not found"
    echo "Expected R1: ${TRIMMED_R1}"
    echo "Expected R2: ${TRIMMED_R2}"
    exit 1
fi

# Function to extract read group information
get_read_group() {
    local fastq=$1
    local sample_name=$2
    
    # Extract header information
    local header=$(zcat "${fastq}" | head -n 1)
    local flowcell=$(echo "${header}" | awk -F: '{print $3}')
    local lane=$(echo "${header}" | awk -F: '{print $4}')
    local barcode=$(echo "${header}" | awk -F: '{print $10}')
    
    # Construct read group string
    echo "@RG\tID:${sample_name}.${flowcell}.${lane}\tSM:${sample_name}\tLB:Lib2\tPU:${flowcell}.${lane}.${barcode}\tPL:Illumina"
}

# Get read group information
RG_STRING=$(get_read_group "${TRIMMED_R1}" "${sample_name}")

# Run BWA-MEM2 alignment
echo "Running BWA-MEM2 alignment..."
bwa-mem2 mem -t 24 \
    -R "${RG_STRING}" \
    "${REF_GENOME}" \
    "${TRIMMED_R1}" \
    "${TRIMMED_R2}" \
    > "${BWAMEM2_DIR}/${sample_name}.sam"

if [ $? -ne 0 ]; then
    echo "Error: BWA-MEM2 alignment failed"
    exit 1
fi

echo "Alignment complete for sample: ${sample_name}"
