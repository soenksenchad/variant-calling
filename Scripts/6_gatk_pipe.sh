#!/bin/bash

# File: Scripts/6_gatk_pipe.sh
#===============================================================================
# GATK Joint Genotyping Script
# Performs joint genotyping on all gVCF files
#===============================================================================

# Source the conda environment
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate BFFcent

echo "Starting GATK joint genotyping pipeline"

# Step 1: GenomicsDBImport
echo "Running GenomicsDBImport..."
gvcf_params=""
for gvcf in "${GATK_DIR}"/*.g.vcf.gz; do
    gvcf_params="${gvcf_params} -V ${gvcf}"
done

if [ -d "${GENDB_DIR}" ]; then 
  rm -rf "${GENDB_DIR}"
fi

gatk --java-options "-Xmx12g -Djava.io.tmpdir=$SLURM_SCRATCH" GenomicsDBImport \
    ${gvcf_params} \
    --genomicsdb-workspace-path "${GENDB_DIR}" \
    --batch-size 50 \
    --interval-padding 100 \
    -L 3 \
    --genomicsdb-shared-posixfs-optimizations 

if [ $? -ne 0 ]; then
    echo "Error: GenomicsDBImport failed"
    exit 1
fi

# Step 2: GenotypeGVCFs
echo "Running GenotypeGVCFs..."
gatk --java-options "-Xmx12g -Djava.io.tmpdir=$SLURM_SCRATCH" GenotypeGVCFs \
    -R "${REF_GENOME}" \
    -V "gendb://${GENDB_DIR}" \
    -O "${VCF_DIR}/raw_variants.vcf.gz" \
    -L 3 

if [ $? -ne 0 ]; then
    echo "Error: GenotypeGVCFs failed"
    exit 1
fi

# Step 3: Select and filter SNPs
echo "Selecting and filtering SNPs..."
gatk --java-options "-Xmx12g -Djava.io.tmpdir=$SLURM_SCRATCH" SelectVariants \
    -R "${REF_GENOME}" \
    -V "${VCF_DIR}/raw_variants.vcf.gz" \
    --select-type-to-include SNP \
    -O "${FILTER_DIR}/raw_snps.vcf.gz"

if [ $? -ne 0 ]; then
    echo "Error: SNP selection failed"
    exit 1
fi

gatk --java-options "-Xmx12g -Djava.io.tmpdir=$SLURM_SCRATCH" VariantFiltration \
    -R "${REF_GENOME}" \
    -V "${FILTER_DIR}/raw_snps.vcf.gz" \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -filter "AF < 0.05" --filter-name "minMAF" \
    -filter "DP < 5 || DP > 100" --filter-name "Depth" \
    -O "${FILTER_DIR}/temp_filtered_snps.vcf.gz"

if [ $? -ne 0 ]; then
    echo "Error: SNP filtration failed"
    exit 1
fi

gatk SelectVariants \
  -V "${FILTER_DIR}/temp_filtered_snps.vcf.gz" \
  --select-type-to-include SNP \
  --restrict-alleles-to BIALLELIC \
  --exclude-filtered true \
  -O "${FILTER_DIR}/filtered_snps.vcf.gz"

echo "GATK joint genotyping pipeline complete"
