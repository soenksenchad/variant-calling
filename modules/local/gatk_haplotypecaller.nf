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