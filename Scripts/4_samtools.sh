#!/bin/bash
# File: Scripts/4_samtools.sh
#===============================================================================
# SAMtools Processing Script
# Processes SAM files through SAMtools pipeline
#===============================================================================

# Source the conda environment
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate BFFcent

# Get the sample name
filename=$(basename "${INPUT_R1}")
sample_name="${filename%_1.fastq.gz}"
input_sam="${BWAMEM2_DIR}/${sample_name}.sam"

echo "Starting SAMtools processing for sample: ${sample_name}"

# Convert SAM to BAM and filter unmapped reads
echo "Converting SAM to BAM..."
samtools view -@ 8 -u \
    -o "${SAMTOOLS_DIR}/${sample_name}_justmap.bam" \
    "${input_sam}"

if [ $? -ne 0 ]; then
    echo "Error: SAM to BAM conversion failed"
    exit 1
fi

# Process BAM file through samtools pipeline
echo "Running SAMtools pipeline..."
samtools collate -@ 8 -O -u \
    "${SAMTOOLS_DIR}/${sample_name}_justmap.bam" - | \
samtools fixmate -@ 8 -m -u - - | \
samtools sort -@ 8 -T /scratch/alpine/soenksen@colostate.edu/tmp/chrm3 -O bam - | \
samtools markdup -@ 8 \
    -f "${STATS_DIR}/${sample_name}_stats_file.txt" \
    -S -d 2500 \
    --mode t \
    --include-fails \
    - "${SAMTOOLS_DIR}/${sample_name}_samtools_pipe.bam"

if [ $? -ne 0 ]; then
    echo "Error: SAMtools pipeline failed"
    exit 1
fi

# Index the final BAM file
echo "Indexing BAM file..."
samtools index -@ 8 "${SAMTOOLS_DIR}/${sample_name}_samtools_pipe.bam"

if [ $? -ne 0 ]; then
    echo "Error: BAM indexing failed"
    exit 1
fi

# Clean up intermediate files
rm -f "${SAMTOOLS_DIR}/${sample_name}_justmap.bam"

echo "SAMtools processing complete for sample: ${sample_name}"
