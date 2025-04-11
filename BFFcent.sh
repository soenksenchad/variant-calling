#!/bin/bash

#===============================================================================
# BFFcent Pipeline Controller
# Description: Main controller script for genomic analysis pipeline
# Usage: ./BFFcent.sh <ref_genome>
#===============================================================================

set -e  # Exit on error

#-------------------------------------------------------------------------------
# Configuration
#-------------------------------------------------------------------------------

# Check command line arguments
if [ $# -ne 1 ]; then
    echo "Usage: $0 <ref_genome>"
    echo "Example: $0 /path/to/reference.fasta"
    exit 1
fi

# Base directories
PROJECT_DIR=$(pwd)
INPUT_DIR="${PROJECT_DIR}/input"
RESULTS_DIR="${PROJECT_DIR}/results"
SCRIPTS_DIR="${PROJECT_DIR}/Scripts"
LOGS_DIR="${RESULTS_DIR}/logs"

# Input
REF_GENOME=$(realpath "$1")

# Output directories
PRETRIM_FASTQC_DIR="${RESULTS_DIR}/pretrim_fastqc"
FASTP_DIR="${RESULTS_DIR}/fastp"
POSTTRIM_FASTQC_DIR="${RESULTS_DIR}/posttrim_fastqc"
BWAMEM2_DIR="${RESULTS_DIR}/bwamem2"
SAMTOOLS_DIR="${RESULTS_DIR}/samtools"
STATS_DIR="${RESULTS_DIR}/stats"
GATK_DIR="${RESULTS_DIR}/GATK"
GENDB_DIR="${RESULTS_DIR}/GenDB"
VCF_DIR="${RESULTS_DIR}/VCF"
FILTER_DIR="${RESULTS_DIR}/VCF/filterSNPs"
PLINK_DIR="${RESULTS_DIR}/plink2"

# SLURM parameters
PARTITION="amilan"
NODES="1"

#-------------------------------------------------------------------------------
# Validation
#-------------------------------------------------------------------------------

# Check if reference genome exists
if [ ! -f "${REF_GENOME}" ]; then
    echo "Error: Reference genome file not found: ${REF_GENOME}"
    exit 1
fi

# Check if input directory exists and contains fastq files
if [ ! -d "${INPUT_DIR}" ]; then
    echo "Error: Input directory not found: ${INPUT_DIR}"
    exit 1
fi

# Check for fastq files
FASTQ_COUNT=$(ls "${INPUT_DIR}"/*_1.fastq.gz 2>/dev/null | wc -l)
if [ "${FASTQ_COUNT}" -eq 0 ]; then
    echo "Error: No FASTQ files found in ${INPUT_DIR}"
    exit 1
fi

#-------------------------------------------------------------------------------
# Directory Setup
#-------------------------------------------------------------------------------

# Create all necessary directories
echo "Creating output directories..."
for dir in "${LOGS_DIR}" "${PRETRIM_FASTQC_DIR}" "${FASTP_DIR}" "${POSTTRIM_FASTQC_DIR}" \
    "${BWAMEM2_DIR}" "${SAMTOOLS_DIR}" "${STATS_DIR}" "${GATK_DIR}" "${GENDB_DIR}" \
    "${VCF_DIR}" "${FILTER_DIR}" "${PLINK_DIR}"; do
    mkdir -p "${dir}"
    if [ ! -d "${dir}" ]; then
        echo "Error: Failed to create directory: ${dir}"
        exit 1
    fi
done

#-------------------------------------------------------------------------------
# Helper Functions
#-------------------------------------------------------------------------------

# Submit a SLURM job and return the job ID
submit_job() {
    local script=$1
    local job_name=$2
    local dependencies=$3
    local resources=$4

    # Verify script exists
    if [ ! -f "${script}" ]; then
        echo "Error: Script not found: ${script}"
        exit 1
    fi

    # Export variables for the job
    export PROJECT_DIR INPUT_DIR RESULTS_DIR REF_GENOME
    export PRETRIM_FASTQC_DIR FASTP_DIR POSTTRIM_FASTQC_DIR
    export BWAMEM2_DIR SAMTOOLS_DIR STATS_DIR
    export GATK_DIR GENDB_DIR VCF_DIR FILTER_DIR PLINK_DIR

    # Build sbatch command
    local sbatch_cmd="sbatch --parsable"
    sbatch_cmd+=" --job-name=${job_name}"
    sbatch_cmd+=" --output=${LOGS_DIR}/${job_name}_%j.out"
    sbatch_cmd+=" --error=${LOGS_DIR}/${job_name}_%j.err"
    sbatch_cmd+=" --partition=${PARTITION}"
    sbatch_cmd+=" --nodes=${NODES}"

    # Add dependencies if specified
    if [ -n "${dependencies}" ]; then
        sbatch_cmd+=" --dependency=${dependencies}"
    fi

    # Add additional resources if specified
    if [ -n "${resources}" ]; then
        sbatch_cmd+=" ${resources}"
    fi

    # Add script to command
    sbatch_cmd+=" ${script}"

    # Submit job and capture job ID
    local job_id
    job_id=$(eval "${sbatch_cmd}")

    # Check if job submission was successful
    if [ $? -eq 0 ]; then
        echo "${job_id}"
    else
        echo "Error: Failed to submit job ${job_name}"
        exit 1
    fi
}

#-------------------------------------------------------------------------------
# Pipeline Execution
#-------------------------------------------------------------------------------

echo "Starting BFFcent pipeline..."
echo "Input directory: ${INPUT_DIR}"
echo "Reference genome: ${REF_GENOME}"
echo "Results directory: ${RESULTS_DIR}"

# Initialize array for sample job IDs
declare -a sample_job_ids

# Step 1: Reference Genome Indexing
echo "1. Submitting reference genome indexing job..."
index_job_id=$(submit_job "${SCRIPTS_DIR}/1_ref_genome_index.sh" \
    "index_ref" "" "--cpus-per-task=16")

# Process each sample
for fastq1 in "${INPUT_DIR}"/*_1.fastq.gz; do
    # Extract sample name without path
    sample_name=$(basename "${fastq1}" "_1.fastq.gz")
    fastq2="${INPUT_DIR}/${sample_name}_2.fastq.gz"
    
    if [ ! -f "${fastq2}" ]; then
        echo "Warning: Missing R2 file for ${sample_name}"
        continue
    fi
    
    echo "Processing sample: ${sample_name}"
    
    # Export sample-specific variables
    export INPUT_R1="${fastq1}"
    export INPUT_R2="${fastq2}"
    
    # Step 2: QC and Trimming
    echo "Running QC and trimming for ${sample_name}..."
    qc_job_id=$(submit_job "${SCRIPTS_DIR}/2_qc_trimm.sh" \
        "qc_${sample_name}" "afterok:${index_job_id}" \
        "--cpus-per-task=4 --mem=16G")
    
    if [ -z "${qc_job_id}" ]; then
        echo "Error: Failed to submit QC job for ${sample_name}"
        continue
    fi
    
    # Step 3: Alignment
    echo "Running alignment for ${sample_name}..."
    align_job_id=$(submit_job "${SCRIPTS_DIR}/3_alignment.sh" \
        "align_${sample_name}" "afterok:${qc_job_id}" \
        "--cpus-per-task=24")
    
    if [ -z "${align_job_id}" ]; then
        echo "Error: Failed to submit alignment job for ${sample_name}"
        continue
    fi
    
    # Step 4: Samtools Processing
    echo "Running samtools processing for ${sample_name}..."
    samtools_job_id=$(submit_job "${SCRIPTS_DIR}/4_samtools.sh" \
        "samtools_${sample_name}" "afterok:${align_job_id}" \
        "--cpus-per-task=24")
    
    if [ -z "${samtools_job_id}" ]; then
        echo "Error: Failed to submit samtools job for ${sample_name}"
        continue
    fi
    
    # Step 5: GATK HaplotypeCaller
    echo "Running GATK HaplotypeCaller for ${sample_name}..."
    gatk_job_id=$(submit_job "${SCRIPTS_DIR}/5_gatk_haplotype.sh" \
        "gatk_${sample_name}" "afterok:${samtools_job_id}" \
        "--cpus-per-task=28 --time=8:00:00")

    if [ -z "${gatk_job_id}" ]; then
        echo "Error: Failed to submit GATK job for ${sample_name}"
        continue
    fi
    
    # Store the final job ID for this sample
    sample_job_ids+=("${gatk_job_id}")
    echo "Successfully submitted all jobs for ${sample_name}"
done

# Only proceed with joint calling if we have submitted jobs successfully
if [ ${#sample_job_ids[@]} -gt 0 ]; then
    # Create dependency string for GATK joint calling
    dependency_string="afterok"
    for job_id in "${sample_job_ids[@]}"; do
        dependency_string+=":${job_id}"
    done

    # Step 6: GATK Pipeline
    echo "Submitting GATK joint calling pipeline..."
    gatk_pipe_job_id=$(submit_job "${SCRIPTS_DIR}/6_gatk_pipe.sh" \
        "gatk_pipe" "${dependency_string}" \
        "--cpus-per-task=8")

    if [ -z "${gatk_pipe_job_id}" ]; then
        echo "Error: Failed to submit GATK joint calling job"
        exit 1
    fi

    # Step 7: PLINK2 Analysis
    echo "Submitting PLINK2 analysis..."
    plink_job_id=$(submit_job "${SCRIPTS_DIR}/7_plink2.sh" \
        "plink" "afterok:${gatk_pipe_job_id}" \
        "--cpus-per-task=8")

    if [ -z "${plink_job_id}" ]; then
        echo "Error: Failed to submit PLINK2 analysis job"
        exit 1
    fi
else
    echo "Error: No sample processing jobs were successfully submitted"
    exit 1
fi

echo "All jobs submitted successfully!"
echo "Monitor progress in: ${LOGS_DIR}"
echo "Results will be in: ${RESULTS_DIR}"
