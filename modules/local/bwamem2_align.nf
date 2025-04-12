process BWAMEM2_ALIGN {
    tag "$meta.id"
    publishDir "${params.outdir}/alignment", mode: params.publish_dir_mode
    cpus 24
    memory '48 GB'
    time '24h'
    
    input:
    tuple val(meta), path(reads), path(reference, stageAs: "genome.fa"), path(rg_info)

    output:
    tuple val(meta), path("${meta.id}.aligned.bam"), emit: bam_files

    script:
    """
    #!/bin/bash
    set -e -o pipefail

    echo "Contents of working directory before symlinks:"
    ls -la > workdir_contents_before.txt
    
    # Use the interpolated Nextflow parameter
    REF_PATH="${params.reference_genome}"
    # Use shell commands to get dir/name
    REF_DIR=$(dirname "$REF_PATH")
    REF_NAME=$(basename "$REF_PATH")
    
    echo "Using reference path: $REF_PATH"
    echo "Reference Dir: $REF_DIR"
    echo "Reference Name: $REF_NAME"

    echo "Creating symlinks for BWA-MEM2 index files..."
    # Use shell variables for source path
    for ext in amb ann bwt.2bit.64 pac sa; do
        src_file="$REF_DIR/$REF_NAME.$ext"
        tgt_link="genome.fa.$ext"
        if [ -e "$src_file" ]; then
            ln -s "$src_file" "$tgt_link"
            if [ ! -e "$tgt_link" ]; then
                echo "Error: Failed to create symlink $tgt_link from $src_file"
                exit 1
            fi
        else
            echo "Error: Source index file $src_file does not exist."
            exit 1
        fi
    done
    
    # Check for additional optional BWA index file
    bwt_src="$REF_DIR/$REF_NAME.bwt"
    bwt_tgt="genome.fa.bwt"
    if [ -e "$bwt_src" ]; then
        ln -s "$bwt_src" "$bwt_tgt"
        echo "Created symlink for $bwt_tgt"
    fi

    echo "Contents of working directory after symlinks:"
    ls -la > workdir_contents_after.txt
    
    # Function to extract read group information
    get_read_group() {
        local fastq=$1
        local sample_name=$2
        
        # Extract header information
        local header=$(zcat "${fastq}" | head -n 1)
        local flowcell=$(echo "${header}" | awk -F: '{print $3}')
        local lane=$(echo "${header}" | awk -F: '{print $4}')
        local barcode=$(echo "${header}" | awk -F: '{print $10}')
        
        # Construct read group string
        echo "@RG\tID:${sample_name}.${flowcell}.${lane}\tSM:${sample_name}\tLB:Lib2\tPU:${flowcell}.${lane}.${barcode}\tPL:Illumina"
    }

    # Get read group information
    RG_STRING=$(get_read_group "${reads[0]}" "${meta.id}")
    echo "Constructed RG string: $RG_STRING"

    echo "Running BWA-MEM2 alignment..."
    # Pass the shell variable RG_STRING to -R
    bwa-mem2 mem -t ${task.cpus} -M -R "$RG_STRING" genome.fa ${reads[0]} ${reads[1]} | \
    samtools view -bS - > ${meta.id}.bam

    echo "Moving output BAM..."
    mv ${meta.id}.bam ${meta.id}.aligned.bam

    echo "Alignment finished."
    """
}