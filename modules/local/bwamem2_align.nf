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
    # List current directory for debugging
    echo "Contents of working directory before symlinks:"
    ls -la > workdir_contents_before.txt
    
    # Define reference path variables
    REF_PATH="${params.reference_genome}"
    REF_DIR=\$(dirname "\$REF_PATH")
    REF_NAME=\$(basename "\$REF_PATH")
    
    # Create symlinks for BWA-MEM2 index files
    for ext in amb ann bwt.2bit.64 pac sa; do
        ln -s "\${REF_PATH}.\$ext" "genome.fa.\$ext"
        if [ ! -e "genome.fa.\$ext" ]; then
            echo "Error: Failed to create symlink for genome.fa.\$ext. Check if \${REF_PATH}.\$ext exists."
            exit 1
        fi
    done
    
    # Check for additional optional BWA index file
    if [ -e "\${REF_PATH}.bwt" ]; then
        ln -s "\${REF_PATH}.bwt" "genome.fa.bwt"
        echo "Created symlink for genome.fa.bwt"
    fi
    
    # List current directory after linking for debugging
    echo "Contents of working directory after symlinks:"
    ls -la > workdir_contents_after.txt
    
    # Read read group info from file using corrected awk and sed
    RG_ID=\$(awk -F'\\t' '{print \$1}' ${rg_info} | sed 's/ID://')\n    RG_SM=\$(awk -F'\\t' '{print \$2}' ${rg_info} | sed 's/SM://')\n    RG_LB=\$(awk -F'\\t' '{print \$3}' ${rg_info} | sed 's/LB://')\n    RG_PU=\$(awk -F'\\t' '{print \$4}' ${rg_info} | sed 's/PU://')\n    RG_PL=\$(awk -F'\\t' '{print \$5}' ${rg_info} | sed 's/PL://')\n    \n    # Construct read group string, escaping special characters for bwa-mem2 -R\n    RG=\"@RG\\\\tID:\\${RG_ID}\\\\tSM:\\${RG_SM}\\\\tLB:\\${RG_LB}\\\\tPL:\\${RG_PL}\\\\tPU:\\${RG_PU}\"\n    \n    # Align with bwa-mem2 including read group\n    bwa-mem2 mem -t ${task.cpus} -M -R \"\\${RG}\" genome.fa ${reads[0]} ${reads[1]} | \\\n    samtools view -bS - > ${meta.id}.bam\n

    # Output the unsorted aligned BAM for SAMTOOLS_PROCESS
    mv ${meta.id}.bam ${meta.id}.aligned.bam
    """
}