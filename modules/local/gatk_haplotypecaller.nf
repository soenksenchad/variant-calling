process GATK_HAPLOTYPECALLER {
    tag "$meta.id - Interval $interval"
    publishDir "${params.outdir}/gatk_haplotype", mode: params.publish_dir_mode
    cpus 28
    memory '64 GB'
    time '24h'

    input:
    tuple val(meta), path(bam_file), path(reference), val(interval)

    output:
    tuple val(meta), val(interval), path("${meta.id}.${interval}.g.vcf.gz"), path("${meta.id}.${interval}.g.vcf.gz.tbi"), emit: gvcfs_interval

    script:
    """
    # Create temporary directory for GATK within task directory
    mkdir -p gatk_tmp
    
    gatk --java-options "-Xmx${task.memory.toGiga()}g -Djava.io.tmpdir=./gatk_tmp" HaplotypeCaller \\
        -R ${reference} \\
        -I ${bam_file} \\
        -O ${meta.id}.${interval}.g.vcf.gz \\
        -L ${interval} \\
        -ERC GVCF \\
        --native-pair-hmm-threads ${task.cpus} \\
        --smith-waterman JAVA

    # Index the GVCF file
    gatk IndexFeatureFile -I "${meta.id}.${interval}.g.vcf.gz"
    """
}