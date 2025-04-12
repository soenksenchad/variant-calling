process CREATE_SAMPLE_MAP {
    tag "$interval"
    // Publish map files for potential debugging/review
    publishDir "${params.outdir}/sample_maps", mode: params.publish_dir_mode, pattern: "*.sample_map", saveAs: { filename -> filename }

    input:
    tuple val(interval), val(gvcf_sample_path_list) // e.g., ["1", [["sampleA", pathA], ["sampleB", pathB]]]

    output:
    tuple val(interval), path("${interval}.sample_map"), emit: sample_map

    script:
    // Use Groovy's collect and join for cleaner map file creation
    def sample_map_content = gvcf_sample_path_list.collect { sample_id, gvcf_path ->
        "${sample_id}\t${gvcf_path.toAbsolutePath()}" // Use absolute paths in map file
    }.join('\\n')

    """
    echo -e '${sample_map_content}' > "${interval}.sample_map"
    """
} 