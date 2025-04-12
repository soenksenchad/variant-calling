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
    # List current directory for debugging
    echo "Contents of working directory before symlinks:"
    ls -la > workdir_contents_before.txt
    
    # Define reference path variables
    REF_PATH="${params.reference_genome}"
    REF_DIR=\$(dirname "\$REF_PATH")
    REF_NAME=\$(basename "\$REF_PATH")
    
    # Create symlinks for BWA-MEM2 index files
    for ext in amb ann bwt.2bit.64 pac 0123; do
        ln -s "\${REF_PATH}.\$ext" "genome.fa.\$ext"
        if [ ! -e "genome.fa.\$ext" ]; then
            echo "Error: Failed to create symlink for genome.fa.\$ext. Check if \${REF_PATH}.\$ext exists."
            exit 1
        fi
    done
    
    # Check for additional standard BWA-MEM2 files
    for ext in bwt sa; do
        if [ -e "\${REF_PATH}.\$ext" ]; then
            ln -s "\${REF_PATH}.\$ext" "genome.fa.\$ext"
            echo "Created symlink for genome.fa.\$ext"
        fi
    done
    
    # List current directory after linking for debugging
    echo "Contents of working directory after symlinks:"
    ls -la > workdir_contents_after.txt
    
    # Align with bwa-mem2
    bwa-mem2 mem -t ${task.cpus} -M genome.fa ${reads[0]} ${reads[1]} | \
    samtools view -bS - > ${meta.id}.bam

    # Sort BAM file
    samtools sort -@ ${task.cpus} -o ${meta.id}.aligned.bam ${meta.id}.bam
    
    # Remove intermediate BAM
    rm ${meta.id}.bam
    """
}