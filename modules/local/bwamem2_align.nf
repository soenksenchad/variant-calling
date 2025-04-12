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
    
    # Align with bwa-mem2
    bwa-mem2 mem -t ${task.cpus} -M genome.fa ${reads[0]} ${reads[1]} | \
    samtools view -bS - > ${meta.id}.bam

    # Sort BAM file
    samtools sort -@ ${task.cpus} -o ${meta.id}.aligned.bam ${meta.id}.bam
    
    # Remove intermediate BAM
    rm ${meta.id}.bam
    """
}