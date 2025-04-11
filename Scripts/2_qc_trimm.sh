#!/bin/bash
# File: Scripts/2_qc_trimm.sh
#===============================================================================
# Quality Control and Trimming Script
# Runs FastQC on raw reads and performs adapter trimming with FastP
#===============================================================================

# Source the conda environment
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate BFFcent

# Extract filename and sample name
filename=$(basename "${INPUT_R1}")
sample_name="${filename%_1.fastq.gz}"

echo "Processing sample: ${sample_name}"
echo "Input R1: ${INPUT_R1}"
echo "Input R2: ${INPUT_R2}"

# Create output directories if they don't exist
mkdir -p "${PRETRIM_FASTQC_DIR}" "${FASTP_DIR}" "${POSTTRIM_FASTQC_DIR}"

# Run pre-trim FastQC
echo "Running pre-trim FastQC..."
fastqc -o "${PRETRIM_FASTQC_DIR}" -t 4 "${INPUT_R1}" "${INPUT_R2}"
if [ $? -ne 0 ]; then
    echo "Error: Pre-trim FastQC failed"
    exit 1
fi

# Run FastP
echo "Running FastP trimming..."
fastp \
    -i "${INPUT_R1}" \
    -I "${INPUT_R2}" \
    -o "${FASTP_DIR}/${sample_name}_1_trimmed.fastq.gz" \
    -O "${FASTP_DIR}/${sample_name}_2_trimmed.fastq.gz" \
    --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
    --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    --detect_adapter_for_pe \
    --cut_right \
    --cut_right_window_size 4 \
    --cut_right_mean_quality 20 \
    --thread 4 \
    --json "${FASTP_DIR}/${sample_name}.json" \
    --html "${FASTP_DIR}/${sample_name}.html"

if [ $? -ne 0 ]; then
    echo "Error: FastP trimming failed"
    exit 1
fi

# Run post-trim FastQC
echo "Running post-trim FastQC..."
fastqc -o "${POSTTRIM_FASTQC_DIR}" -t 4 \
    "${FASTP_DIR}/${sample_name}_1_trimmed.fastq.gz" \
    "${FASTP_DIR}/${sample_name}_2_trimmed.fastq.gz"

if [ $? -ne 0 ]; then
    echo "Error: Post-trim FastQC failed"
    exit 1
fi

echo "QC and trimming complete for sample: ${sample_name}"
