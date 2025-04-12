process GATK_HAPLOTYPECALLER {
    tag "$sample_id:$interval"
    publishDir "${params.outdir}/gatk_haplotype/${sample_id}", mode: params.publish_dir_mode
    cpus 28
    memory '64 GB'
    time '24h'

    input:
    tuple val(sample_id), path(bam_file), path(ref_genome), val(interval)

    output:
    tuple val(sample_id), val(interval), path("${sample_id}_${interval}.g.vcf.gz"), path("${sample_id}_${interval}.g.vcf.gz.tbi"), emit: gvcfs_interval

    script:
    """
    # Create temporary directory for GATK within task directory
    mkdir -p gatk_tmp
    
    gatk --java-options "-Xmx${task.memory.toGiga()}g -Djava.io.tmpdir=./gatk_tmp" HaplotypeCaller \\
        -R "${ref_genome}" \\
        -I "${bam_file}" \\
        -O "${sample_id}_${interval}.g.vcf.gz" \\
        --native-pair-hmm-threads ${task.cpus} \\
        -ERC GVCF \\
        -L ${interval}

    # Index the GVCF file
    gatk IndexFeatureFile -I "${sample_id}_${interval}.g.vcf.gz"
    """
}