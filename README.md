# nf_variant_caller: Maize Variant Calling Pipeline

A Nextflow pipeline for variant discovery and population genomics analysis of maize hapoloid crossing experiment.

## Overview

This pipeline performs the following steps:
- Quality control and trimming of Illumina paired-end reads (FASTP)
- Read alignment to reference genome (BWA-MEM2)
- SAM/BAM processing and duplicate marking (Samtools)
- Variant calling (GATK HaplotypeCaller)
- Joint genotyping and variant filtering (GATK)

## Requirements
- Nextflow 21.10.0 or later
- Conda or Mamba (recommended for faster environment creation)
- SLURM workload manager (optional, for HPC execution)

## Prerequisites
- Nextflow (>= 21.10.0)
- SLURM scheduler
- Reference genome and index files
- Conda with a variant-calling environment (containing bwa-mem2, samtools, gatk, fastp)

## Reference Genome Preparation

**Important:** Before running the pipeline, you must manually index your reference genome and place it in the `reference` directory.

### Setting Up The Reference Directory

1. Place your reference genome file in the `reference` directory:
```bash
# Create the reference directory if it doesn't exist
mkdir -p reference

# Copy your reference genome to the reference directory
cp /path/to/your/reference.fa reference/
```

2. Index the reference genome file:
```bash
# Create a SLURM job script for indexing
cat > index_genome.sh << 'EOF'
#!/bin/bash
#SBATCH --partition=amem
#SBATCH --qos=mem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=600G
#SBATCH --time=24:00:00

# Path to reference genome in the reference directory
REF_GENOME="./reference/your_reference.fa"

# Load required modules or activate conda environment
module load bwa-mem2 samtools gatk
# OR
# source activate variant-calling

# Create BWA-MEM2 index
echo "Creating BWA-MEM2 index..."
bwa-mem2 index $REF_GENOME

# Create SAMtools FASTA index
echo "Creating SAMtools FASTA index..."
samtools faidx $REF_GENOME

# Create GATK sequence dictionary
echo "Creating GATK sequence dictionary..."
gatk CreateSequenceDictionary -R $REF_GENOME
EOF

# Submit the job
sbatch index_genome.sh
```

The pipeline requires all these index files to be present in the `reference` directory alongside the reference FASTA file:
- BWA-MEM2 index files (*.amb, *.ann, *.bwt, *.pac, *.sa, *.0123)
- SAMtools FASTA index (.fai)
- GATK sequence dictionary (.dict)

### Alternative: Using An Existing Pre-Indexed Reference Genome

If you already have an indexed reference genome, you can create symbolic links to your existing files:

```bash
# Create the reference directory
mkdir -p reference

# Create symbolic links to reference genome and index files
ln -s /path/to/your/indexed/reference.fa reference/
ln -s /path/to/your/indexed/reference.fa.* reference/
ln -s /path/to/your/indexed/reference.dict reference/
```

Make sure all necessary index files are linked or copied to the reference directory.

## Usage

### 1. Installation

Clone this repository:
```bash
git clone https://github.com/yourusername/variant_calling.git
cd variant_calling
```

### 2. Conda Environment

The pipeline will automatically create and manage the required Conda environment using the provided `environment.yml` file. The environment will be created during the first run and reused in subsequent runs, so you don't need to manually create it.

If you want to make changes to the environment, simply edit the `environment.yml` file before running the pipeline.

### 3. Sample file preparation

Prepare your sample file (`samples.txt`) with tab-separated values:
```
sample_id    read1    read2
sample1    /path/to/sample1_R1.fastq.gz    /path/to/sample1_R2.fastq.gz
sample2    /path/to/sample2_R1.fastq.gz    /path/to/sample2_R2.fastq.gz
```

If all of your files are in a directory you can use the script create_sample_list_from_dir.sh to automate the creation of the samples.txt file

```
./modules/local/create_sample_list_from_dir.sh /path/to/fastq_files/dir/
```

### 4. Reference genome

The pipeline is configured to use the maize reference genome. The path is set in the `nextflow.config` file. Before running, make sure to update this path to point to your copy of the reference genome:

```groovy
// In nextflow.config
params {
    reference_genome = "/path/to/maize_reference_genome.fasta"  // Update this path
    // ...
}
```

If needed, you can temporarily override this setting at runtime:
```bash
nextflow run main.nf --reference_genome /alternate/path/to/genome.fasta
```

### 5. Running the pipeline

Basic usage:
```bash
nextflow run main.nf \
  --input samples.txt \
  --outdir results
```

Using SLURM with resume capability:
```bash
nextflow run main.nf \
  -profile slurm \
  -resume \
  --outdir results
```

#### Additional Parameters

```bash
# With SLURM account
nextflow run main.nf \
  --input samples.txt \
  --outdir results \
  --slurm_account your_account \
  --slurm_partition partition_name \
  --slurm_qos qos_name

# With custom resources
nextflow run main.nf \
  --input samples.txt \
  --outdir results \
  --fastp_partition short-cpu \
  --alignment_partition long-highmem
```

## Pipeline Parameters

### Main Parameters
- `--input`: Path to samples.txt file (default: 'samples.txt')
- `--outdir`: Output directory (default: './results')
- `--temp_dir`: Directory for temporary files (default: './temp')

### SLURM Parameters
- `--slurm_account`: SLURM account (optional)
- `--slurm_partition`: Default SLURM partition (default: 'day-long-cpu')
- `--slurm_qos`: Default SLURM QOS (optional)

### Process-specific Parameters
Each process can have a custom partition and QOS defined:
- `--fastp_partition`, `--fastp_qos`
- `--alignment_partition`, `--alignment_qos`
- etc.

## Output Structure

```
results/
├── ref_index/          # Indexed reference genome files
├── trimming/           # Trimmed FASTQ files and QC reports
├── alignment/          # Aligned BAM files
├── readgroups/         # Read group information and summary
├── samtools/           # Processed BAM files
├── gatk_haplotype/     # Individual GVCF files
├── gatk_pipe/          # Joint-called VCF files
└── pipeline_info/      # Execution reports and logs
```

## Read Group Extraction

The pipeline automatically extracts read group information from FASTQ headers to ensure proper sample tracking throughout the analysis. This information is used for the alignment step and is also saved in tabular format for reference.

The read group information includes:
- Sample ID
- Read Group ID (constructed from sample name, library, flowcell, and lane)
- Sample Name (SM tag)
- Library (LB tag)
- Platform Unit (PU tag, constructed from flowcell, lane, and barcode)
- Platform (PL tag, set to "Illumina")

A master table with all samples' read group information is generated in the `results/readgroups/` directory after pipeline execution. This file includes validation checks to ensure all necessary information was properly extracted.

## License

This pipeline is licensed under the MIT License.

## Contact

For questions or issues, please create a GitHub issue or contact the author.

## Troubleshooting

### Reference Genome Index Issues

If you encounter errors related to reference genome index files:

- **BWA-MEM2 Error: Unable to open genome.fa.bwt.2bit.64**: Ensure all BWA-MEM2 index files (`.amb`, `.ann`, `.bwt.2bit.64`, `.pac`, `.0123`) exist in the reference directory.
- **GATK Error**: Verify `.fai` and `.dict` files are present.
- Check `workdir_contents_before.txt` and `workdir_contents_after.txt` in the failed work directory for debugging information.

### Complete List of Required Index Files

The pipeline requires these index files alongside your reference FASTA:

| File Extension | Created By | Purpose |
|----------------|------------|---------|
| `.amb` | BWA-MEM2 index | BWA-MEM2 alignment |
| `.ann` | BWA-MEM2 index | BWA-MEM2 alignment |
| `.bwt` | BWA-MEM2 index | BWA-MEM2 alignment |
| `.bwt.2bit.64` | BWA-MEM2 index | BWA-MEM2 alignment |
| `.pac` | BWA-MEM2 index | BWA-MEM2 alignment |
| `.0123` | BWA-MEM2 index | BWA-MEM2 alignment |
| `.sa` | BWA-MEM2 index | BWA-MEM2 alignment |
| `.fai` | samtools faidx | GATK and samtools |
| `.dict` | GATK CreateSequenceDictionary | GATK tools |

### Common Errors and Solutions

- **Error: No files matched pattern cp "$REF_DIR"/"$REF_NAME"***: The reference directory doesn't contain the expected index files.
- **Error: Failed to create symlink for genome.fa.$ext**: The index file doesn't exist or cannot be accessed.
- **OOM (Out of Memory) Errors**: Increase memory allocation in nextflow.config for the affected process.

### Job Debugging

For failed processes, you can examine the work directory and logs:

```bash
# Find the work directory of a failed process
grep -A1 "Process.*terminated" .nextflow.log

# Examine output and error files
cd /path/to/work/directory
cat .command.out
cat .command.err
cat workdir_contents_before.txt
cat workdir_contents_after.txt
```
