process BWAMEM2_ALIGN {
    tag "$meta.id"
    publishDir "${params.outdir}/alignment", mode: params.publish_dir_mode
    cpus 24
    memory '48 GB'
    time '24h'
    
    input:
    tuple val(meta), path(reads), path(reference, stageAs: "genome.fa")

    output:
    tuple val(meta), path("${meta.id}.aligned.bam"), emit: bam_files

    script:
    """
    # Align with bwa-mem2
    bwa-mem2 mem -t ${task.cpus} -M genome.fa ${reads[0]} ${reads[1]} | \
    samtools view -bS - > ${meta.id}.bam

    # Sort BAM file
    samtools sort -@ ${task.cpus} -o ${meta.id}.aligned.bam ${meta.id}.bam
    
    # Remove intermediate BAM
    rm ${meta.id}.bam
    """
}