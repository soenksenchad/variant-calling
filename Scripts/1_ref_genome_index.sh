#!/bin/bash
# File: Scripts/1_ref_genome_index.sh
#===============================================================================
# Reference Genome Indexing Script
# Creates BWA-MEM2 index and required GATK files for reference genome
#===============================================================================

# Source the conda environment
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate BFFcent

echo "Starting reference genome preparation..."
echo "Reference genome: ${REF_GENOME}"

# Check if BWA-MEM2 index files exist
echo "Checking for existing BWA-MEM2 indices..."
missing_index=false
for ext in ".amb" ".ann" ".bwt.2bit.64" ".pac"; do
    if [ ! -f "${REF_GENOME}${ext}" ]; then
        missing_index=true
        break
    fi
done

if [ "$missing_index" = true ]; then
    echo "Creating BWA-MEM2 index..."
    bwa-mem2 index "${REF_GENOME}"
    if [ $? -ne 0 ]; then
        echo "Error: BWA-MEM2 indexing failed"
        exit 1
    fi
else
    echo "BWA-MEM2 index already exists"
fi

# Create GATK sequence dictionary
dict_file="${REF_GENOME%.*}.dict"
echo "Checking for sequence dictionary..."
if [ ! -s "${dict_file}" ]; then
    echo "Creating sequence dictionary..."
    rm -f "${dict_file}"  # Remove if exists but empty
    gatk CreateSequenceDictionary \
        -R "${REF_GENOME}" \
        -O "${dict_file}"
    
    if [ ! -s "${dict_file}" ]; then
        echo "Error: Failed to create sequence dictionary"
        exit 1
    fi
else
    echo "Sequence dictionary already exists"
fi

# Create reference FASTA index
echo "Checking for FASTA index..."
if [ ! -f "${REF_GENOME}.fai" ]; then
    echo "Creating FASTA index..."
    samtools faidx "${REF_GENOME}"
    if [ $? -ne 0 ]; then
        echo "Error: Failed to create FASTA index"
        exit 1
    fi
else
    echo "FASTA index already exists"
fi

echo "Reference genome preparation complete"
