process GATK_JOINT_GENOTYPE {
    tag "Interval $interval"
    publishDir "${params.outdir}/gatk_pipe", mode: params.publish_dir_mode
    cpus 32
    time '156h'

    input:
    tuple val(interval), path(sample_map), path(reference, stageAs: "genome.fa")

    output:
    tuple val(interval), path("${interval}.filtered.vcf.gz"), path("${interval}.filtered.vcf.gz.tbi"), emit: filtered_vcf_interval

    script:
    """
    # List current directory for debugging
    echo "Contents of working directory before symlinks:"
    ls -la > workdir_contents_before.txt
    
    # Define reference path variables
    REF_PATH="${params.reference_genome}"
    REF_DIR=\$(dirname "\$REF_PATH")
    REF_NAME=\$(basename "\$REF_PATH")
    REF_BASE=\$(echo "\$REF_NAME" | sed 's/\\.[^.]*\$//')
    
    # Create symlink for .fai index
    ln -s "\${REF_PATH}.fai" "genome.fa.fai"
    if [ ! -e "genome.fa.fai" ]; then
        echo "Error: genome.fa.fai not found. Check if \${REF_PATH}.fai exists."
        exit 1
    fi
    
    # Create symlink for .dict file
    ln -s "\${REF_DIR}/\${REF_BASE}.dict" "genome.dict"
    if [ ! -e "genome.dict" ]; then
        echo "Error: genome.dict not found. Check if \${REF_DIR}/\${REF_BASE}.dict exists."
        exit 1
    fi
    
    # List current directory after linking for debugging
    echo "Contents of working directory after symlinks:"
    ls -la > workdir_contents_after.txt
    
    # Create a unique workspace for this interval
    WORKSPACE="genomicsdb_${interval}"
    
    # Import GVCFs to GenomicsDB
    gatk --java-options "-Xmx${task.memory.toGiga()}g" GenomicsDBImport \\
        --genomicsdb-workspace-path \$WORKSPACE \\
        --batch-size 50 \\
        -L ${interval} \\
        --sample-name-map ${sample_map} \\
        --reader-threads ${Math.max(1, task.cpus/4 as int)} \\
        --max-num-intervals-to-import-in-parallel ${Math.max(1, task.cpus/8 as int)}

    # Run joint genotyping on the interval
    gatk --java-options "-Xmx${task.memory.toGiga()}g" GenotypeGVCFs \\
        -R genome.fa \\
        -V gendb://\$WORKSPACE \\
        -O ${interval}.vcf.gz \\
        -L ${interval} \\
        --tmp-dir=./tmp

    # Apply variant filtration
    gatk --java-options "-Xmx${task.memory.toGiga()}g" VariantFiltration \\
        -R genome.fa \\
        -V ${interval}.vcf.gz \\
        -O ${interval}.filtered.vcf.gz \\
        --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \\
        --filter-name "basic_snp_filter"
    """
}