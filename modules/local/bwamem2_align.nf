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
    // Create a symlink to all reference index files in current directory
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
    
    # Align with bwa-mem2
    bwa-mem2 mem -t ${task.cpus} -M genome.fa ${reads[0]} ${reads[1]} | \
    samtools view -bS - > ${meta.id}.bam

    # Sort BAM file
    samtools sort -@ ${task.cpus} -o ${meta.id}.aligned.bam ${meta.id}.bam
    
    # Remove intermediate BAM
    rm ${meta.id}.bam
    """
}