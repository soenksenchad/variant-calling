process GATK_JOINT_GENOTYPE {
    tag "joint_genotyping_${interval}"
    publishDir "${params.outdir}/gatk_pipe", mode: params.publish_dir_mode, pattern: "filtered_snps_${interval}.*"
    cpus 32
    memory '64 GB'
    time '24h'

    input:
    tuple val(interval), path(sample_map)
    tuple path(ref_genome), path(ref_genome_dict), path(ref_genome_fasta_index)

    output:
    tuple val(interval), path("filtered_snps_${interval}.vcf.gz"), emit: filtered_vcf_interval
    tuple val(interval), path("filtered_snps_${interval}.vcf.gz.tbi"), emit: filtered_vcf_interval_index

    script:
    """
    # Create temporary directory for GATK within task directory
    mkdir -p gatk_tmp
    gendb_path="gendb_${interval}"
    rm -rf \${gendb_path} # Clean up potential leftovers from previous runs
    mkdir -p \${gendb_path}

    # Step 1: Import GVCFs using sample map
    gatk --java-options "-Xmx${task.memory.toGiga() - 10}g -Djava.io.tmpdir=./gatk_tmp" GenomicsDBImport \\
        --genomicsdb-workspace-path \${gendb_path} \\
        --batch-size 50 \\
        --sample-name-map ${sample_map} \\
        -L ${interval} \\
        --interval-padding 100 \\
        --genomicsdb-shared-posixfs-optimizations \\
        --reader-threads ${task.cpus.intdiv(2) ?: 1}

    # Step 2: GenotypeGVCFs
    gatk --java-options "-Xmx${task.memory.toGiga() - 10}g -Djava.io.tmpdir=./gatk_tmp" GenotypeGVCFs \\
        -R "${ref_genome}" \\
        -V "gendb://\${gendb_path}" \\
        -O "raw_variants_${interval}.vcf.gz" \\
        -L ${interval} \\
        --reader-threads ${task.cpus.intdiv(2) ?: 1}

    # Step 3: Select and filter SNPs for the interval
    gatk --java-options "-Xmx${task.memory.toGiga() - 10}g -Djava.io.tmpdir=./gatk_tmp" SelectVariants \\
        -R "${ref_genome}" \\
        -V "raw_variants_${interval}.vcf.gz" \\
        --select-type-to-include SNP \\
        -O "raw_snps_${interval}.vcf.gz"

    gatk --java-options "-Xmx${task.memory.toGiga() - 10}g -Djava.io.tmpdir=./gatk_tmp" VariantFiltration \\
        -R "${ref_genome}" \\
        -V "raw_snps_${interval}.vcf.gz" \\
        -filter "QD < 2.0" --filter-name "QD2" \\
        -filter "QUAL < 30.0" --filter-name "QUAL30" \\
        -filter "SOR > 3.0" --filter-name "SOR3" \\
        -filter "FS > 60.0" --filter-name "FS60" \\
        -filter "MQ < 40.0" --filter-name "MQ40" \\
        -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \\
        -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \\
        -filter "AF < 0.05" --filter-name "minMAF" \\
        -filter "DP < 5 || DP > 100" --filter-name "Depth" \\
        -O "temp_filtered_snps_${interval}.vcf.gz"

    gatk SelectVariants \\
        -V "temp_filtered_snps_${interval}.vcf.gz" \\
        --select-type-to-include SNP \\
        --restrict-alleles-to BIALLELIC \\
        --exclude-filtered true \\
        -O "filtered_snps_${interval}.vcf.gz"
    """
}