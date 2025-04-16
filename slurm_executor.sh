#! /bin/bash

#SBATCH --job-name=variant_calling
#SBATCH --time=156:00:00
#SBATCH --partition=amilan
#SBATCH --cpus-per-task=1
#SBATCH --qos=long
#SBATCH --output=variant_calling_%j.out
#SBATCH --error=variant_calling_%j.err

source ~/.bashrc
nextflow run main.nf -profile slurm
