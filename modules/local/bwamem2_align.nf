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
    // Use a shell block with explicit newlines for clarity and correctness
    '''
    #!/bin/bash
    set -e -o pipefail

    echo "Contents of working directory before symlinks:"
    ls -la > workdir_contents_before.txt
    
    REF_PATH="${params.reference_genome}"
    REF_DIR=$(dirname "$REF_PATH")
    REF_NAME=$(basename "$REF_PATH")
    
    echo "Creating symlinks for BWA-MEM2 index files..."
    for ext in amb ann bwt.2bit.64 pac sa; do
        ln -s "${REF_PATH}.$ext" "genome.fa.$ext"
        if [ ! -e "genome.fa.$ext" ]; then
            echo "Error: Failed to create symlink for genome.fa.$ext. Check if ${REF_PATH}.$ext exists."
            exit 1
        fi
    done
    
    if [ -e "${REF_PATH}.bwt" ]; then
        ln -s "${REF_PATH}.bwt" "genome.fa.bwt"
        echo "Created symlink for genome.fa.bwt"
    fi
    
    echo "Contents of working directory after symlinks:"
    ls -la > workdir_contents_after.txt
    
    echo "Reading read group info from ${rg_info}..."
    RG_ID=$(awk -F'\t' '{print $1}' ${rg_info} | sed 's/ID://')
    RG_SM=$(awk -F'\t' '{print $2}' ${rg_info} | sed 's/SM://')
    RG_LB=$(awk -F'\t' '{print $3}' ${rg_info} | sed 's/LB://')
    RG_PU=$(awk -F'\t' '{print $4}' ${rg_info} | sed 's/PU://')
    RG_PL=$(awk -F'\t' '{print $5}' ${rg_info} | sed 's/PL://')

    # Check if variables were extracted
    if [ -z "$RG_ID" ] || [ -z "$RG_SM" ] || [ -z "$RG_LB" ] || [ -z "$RG_PU" ] || [ -z "$RG_PL" ]; then
        echo "Error: Failed to extract one or more read group fields from ${rg_info}"
        echo "Content of ${rg_info}:"
        cat ${rg_info}
        exit 1
    fi
    
    # Construct read group string for BWA-MEM2 -R option
    # Ensure tabs are literal tabs for the @RG header
    RG_STRING=$(printf '@RG\tID:%s\tSM:%s\tLB:%s\tPL:%s\tPU:%s' "$RG_ID" "$RG_SM" "$RG_LB" "$RG_PL" "$RG_PU")
    echo "Constructed RG string: $RG_STRING"

    echo "Running BWA-MEM2 alignment..."
    bwa-mem2 mem -t ${task.cpus} -M -R "$RG_STRING" genome.fa ${reads[0]} ${reads[1]} | \
    samtools view -bS - > ${meta.id}.bam

    echo "Moving output BAM..."
    mv ${meta.id}.bam ${meta.id}.aligned.bam

    echo "Alignment finished."
    '''
}