process GATK_HAPLOTYPECALLER {
    tag "$meta.id - Interval $interval"
    publishDir "${params.outdir}/gatk_haplotype", mode: params.publish_dir_mode
    cpus 28
    memory '64 GB'
    time '24h'

    input:
    tuple val(meta), path(bam_file), val(interval), path(reference, stageAs: "genome.fa")

    output:
    tuple val(meta), val(interval), path("${meta.id}.${interval}.g.vcf.gz"), path("${meta.id}.${interval}.g.vcf.gz.tbi"), emit: gvcfs_interval

    script:
    """
    # Copy all reference index files to work directory
    REF_PATH="${params.reference_genome}"
    REF_DIR=\$(dirname "\$REF_PATH")
    REF_NAME=\$(basename "\$REF_PATH")
    
    # Copy all index files with *.ext pattern
    cp "\$REF_DIR"/"\$REF_NAME".* ./
    
    # Copy dict file if it exists
    REF_BASE=\$(echo "\$REF_NAME" | sed 's/\\.[^.]*\$//')
    if [ -f "\$REF_DIR/\$REF_BASE.dict" ]; then
        cp "\$REF_DIR/\$REF_BASE.dict" ./genome.dict
    fi
    
    # Create temporary directory for GATK within task directory
    mkdir -p gatk_tmp
    
    gatk --java-options "-Xmx${task.memory.toGiga()}g -Djava.io.tmpdir=./gatk_tmp" HaplotypeCaller \\
        -R genome.fa \\
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