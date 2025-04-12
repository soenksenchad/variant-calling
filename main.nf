#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Define Parameters
params.reads = "samples.txt" // Input samplesheet (sample,fastq_1,fastq_2)
params.reference_genome = "ref/genome.fasta" // Reference genome path (must already be indexed)
params.intervals_file = "intervals_chromo_only_nopos.list" // File listing intervals (1-10)
params.outdir = "./results"
params.publish_dir_mode = 'copy' // Or 'link', 'rellink', etc.

// Define Log Level
log.info """\
         V A R I A N T - C A L L I N G   P I P E L I N E
         ================================================
         reads              : ${params.reads}
         reference_genome   : ${params.reference_genome}
         intervals          : ${params.intervals_file}
         outdir             : ${params.outdir}
         ------------------------------------------------
         """
         .stripIndent()

// Include Modules
include { FASTP_TRIM } from './modules/local/fastp_trim'
include { BWAMEM2_ALIGN } from './modules/local/bwamem2_align'
include { SAMTOOLS_PROCESS } from './modules/local/samtools_process'
include { GATK_HAPLOTYPECALLER } from './modules/local/gatk_haplotypecaller'
include { CREATE_SAMPLE_MAP } from './modules/local/create_sample_map'
include { GATK_JOINT_GENOTYPE } from './modules/local/gatk_joint_genotype'

// Workflow Definition
workflow {

    // --- Input Channels ---
    ch_reads = Channel.fromPath(params.reads)
        .splitCsv(header: true, sep: '\t')
        .map { row -> 
            def meta = [:]
            meta.sample_id = row.sample_id
            
            [ meta, file(row.read1), file(row.read2) ]
        }
    ch_intervals = Channel.fromPath(params.intervals_file).splitText().map { it.trim() }.filter { it } // Emits '1', '2', ... '10'
    
    // Reference genome files (assumes pre-indexed)
    def ref_file = file(params.reference_genome)
    def ref_base = ref_file.getBaseName()
    def ref_dir = ref_file.getParent()
    def ref_name = ref_file.getName()
    
    // Expected index file paths
    def dict_file = file("${ref_dir}/${ref_base}.dict")
    def fai_file = file("${ref_dir}/${ref_name}.fai")
    
    // Check if index files exist
    if (!dict_file.exists()) {
        log.error "Reference dictionary file not found: ${dict_file}"
        log.error "Please index your reference genome before running the pipeline."
        log.error "See README.md for instructions on how to index your reference genome."
        exit 1
    }
    
    if (!fai_file.exists()) {
        log.error "Reference FASTA index file not found: ${fai_file}"
        log.error "Please index your reference genome before running the pipeline."
        log.error "See README.md for instructions on how to index your reference genome."
        exit 1
    }
    
    // Create reference channel
    ch_ref_indexed = Channel.of(tuple(ref_file, dict_file, fai_file))
    
    // --- Step 1: Read Trimming ---
    FASTP_TRIM(ch_reads)

    // --- Step 2: Alignment ---
    ch_align_input = FASTP_TRIM.out.trimmed_reads.combine(ch_ref_indexed) // Input: [sample_id, [r1,r2], [ref, dict, fai]]
    BWAMEM2_ALIGN(ch_align_input.map{ sample_id, reads, ref_files -> tuple(sample_id, reads, ref_files[0]) }) // Input: [sample_id, [r1,r2], ref_path]

    // --- Step 3: SAMtools Processing ---
    SAMTOOLS_PROCESS(BWAMEM2_ALIGN.out.bam_files) // Input: [sample_id, bam_path]

    // --- Step 4: Haplotype Calling (per sample, per interval) ---
    ch_haplotype_input = SAMTOOLS_PROCESS.out.processed_bam
                             .combine(ch_ref_indexed) // Combine with ref [sample_id, bam, [ref, dict, fai]]
                             .combine(ch_intervals)   // Combine with intervals [sample_id, bam, [ref, dict, fai], interval]
                             .map { sample_id, bam, ref_files, interval -> tuple(sample_id, bam, ref_files[0], interval) } // Reformat [sample_id, bam, ref, interval]

    GATK_HAPLOTYPECALLER(ch_haplotype_input) // Emits: [sample_id, interval, gvcf, tbi]

    // --- Step 5: Prepare for Joint Genotyping (Create Sample Map per Interval) ---
    ch_gvcf_info = GATK_HAPLOTYPECALLER.out.gvcfs_interval
        .map { sample_id, interval, gvcf_path, tbi_path ->
            // Check if the TBI file actually exists
            if (file(tbi_path).exists()) {
                return tuple(interval, tuple(sample_id, gvcf_path)) // Format for grouping: [interval, [sample_id, gvcf_path]]
            } else {
                log.warn "Missing index file ${tbi_path} for ${gvcf_path}. Skipping this GVCF for interval ${interval}."
                return null // Indicate failure
            }
        }
        .filter { it != null } // Remove entries where index was missing

    // Group by interval -> [interval, [[sample1, path1], [sample2, path2], ...]]
    ch_grouped_gvcf_info = ch_gvcf_info.groupTuple()

    CREATE_SAMPLE_MAP(ch_grouped_gvcf_info) // Emits: [interval, map_file_path]

    // --- Step 6: Joint Genotyping (per interval) ---
    // Pass the two required input channels as separate arguments
    GATK_JOINT_GENOTYPE(CREATE_SAMPLE_MAP.out.sample_map, ch_ref_indexed) // Emits: [interval, vcf_path, tbi_path]

    // --- Workflow Output ---
    ch_final_vcfs = GATK_JOINT_GENOTYPE.out.filtered_vcf_interval
    ch_final_vcfs.view { interval, vcf -> "Final VCF for interval ${interval}: ${vcf}" }

}

// Workflow Completion Message
workflow.onComplete {
    log.info ( workflow.success ? """
        Pipeline execution summary
        ---------------------------
        Completed successfully at: ${workflow.complete}
        Duration               : ${workflow.duration}
        CPU hours              : ${workflow.cpuHours}
        Memory Mb-hours        : ${workflow.memoryMbHours}
        Work directory         : ${workflow.workDir}
        Output directory       : ${params.outdir}
        ================================================
        Pipeline execution finished!
        """ : """
        Pipeline execution failed!
        ---------------------------
        Error summary: ${workflow.errorReport}
        Work directory: ${workflow.workDir}
        Check the log file for details: .nextflow.log
        ================================================
        """
    )
}
