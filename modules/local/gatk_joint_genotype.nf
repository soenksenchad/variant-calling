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
    // Create symlinks to all reference index files
    def ref_dir = new File(reference).getParent()
    def ref_name = reference.getName()
    def ref_base = ref_name.take(ref_name.lastIndexOf('.'))
    
    """
    # Create symlinks to all reference index files
    for idx_file in ${ref_dir}/${ref_name}.*; do
        ln -sf \$idx_file ./\$(basename \$idx_file)
    done
    
    # Special case for .dict file
    if [ -f "${ref_dir}/${ref_base}.dict" ]; then
        ln -sf "${ref_dir}/${ref_base}.dict" ./genome.dict
    fi
    
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