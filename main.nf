#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Define Parameters
params.reads = "samples.txt" // Input samplesheet (sample,fastq_1,fastq_2)
params.ref_genome = "ref/genome.fasta.gz" // Reference genome
params.ref_genome_dict = "ref/genome.dict" // Reference dictionary (will be created if not exists)
params.ref_genome_fasta_index = "ref/genome.fasta.fai" // Reference FASTA index (will be created if not exists)
params.intervals_file = "intervals_chromo_only_nopos.list" // File listing intervals (1-10)
params.outdir = "./results"
params.publish_dir_mode = 'copy' // Or 'link', 'rellink', etc.

// Define Log Level
log.info """\
         V A R I A N T - C A L L I N G   P I P E L I N E
         ================================================
         reads         : ${params.reads}
         ref_genome    : ${params.ref_genome}
         intervals     : ${params.intervals_file}
         outdir        : ${params.outdir}
         ------------------------------------------------
         """
         .stripIndent()

// Include Modules
include { FASTP_TRIM } from './modules/local/fastp_trim'
include { REFERENCE_INDEX } from './modules/local/reference_index'
include { BWAMEM2_ALIGN } from './modules/local/bwamem2_align'
include { SAMTOOLS_PROCESS } from './modules/local/samtools_process'
include { GATK_HAPLOTYPECALLER } from './modules/local/gatk_haplotypecaller'
include { CREATE_SAMPLE_MAP } from './modules/local/create_sample_map'
include { GATK_JOINT_GENOTYPE } from './modules/local/gatk_joint_genotype'

// Workflow Definition
workflow {

    // --- Input Channels ---
    ch_reads = Channel.fromSamplesheet(params.reads) // Emits [ meta, reads_r1, reads_r2 ] -> use meta.id for sample_id
    ch_intervals = Channel.fromPath(params.intervals_file).splitText().map { it.trim() }.filter { it } // Emits '1', '2', ... '10'
    ch_ref_tuple = Channel.of(file(params.ref_genome), file(params.ref_genome_dict), file(params.ref_genome_fasta_index))

    // --- Step 1: Reference Indexing ---
    REFERENCE_INDEX(ch_ref_tuple)
    ch_ref_indexed = REFERENCE_INDEX.out.indexed_ref // Emits [ref, dict, fai]

    // --- Step 2: Read Trimming ---
    FASTP_TRIM(ch_reads.map{ meta, r1, r2 -> tuple(meta.id, [r1, r2]) }) // Input: [sample_id, [r1, r2]]

    // --- Step 3: Alignment ---
    ch_align_input = FASTP_TRIM.out.trimmed_reads.combine(ch_ref_indexed) // Input: [sample_id, [r1,r2], [ref, dict, fai]]
    BWAMEM2_ALIGN(ch_align_input.map{ sample_id, reads, ref_files -> tuple(sample_id, reads, ref_files[0]) }) // Input: [sample_id, [r1,r2], ref_path]

    // --- Step 4: SAMtools Processing ---
    SAMTOOLS_PROCESS(BWAMEM2_ALIGN.out.bam_files) // Input: [sample_id, bam_path]

    // --- Step 5: Haplotype Calling (per sample, per interval) ---
    ch_haplotype_input = SAMTOOLS_PROCESS.out.processed_bam
                             .combine(ch_ref_indexed) // Combine with ref [sample_id, bam, [ref, dict, fai]]
                             .combine(ch_intervals)   // Combine with intervals [sample_id, bam, [ref, dict, fai], interval]
                             .map { sample_id, bam, ref_files, interval -> tuple(sample_id, bam, ref_files[0], interval) } // Reformat [sample_id, bam, ref, interval]

    GATK_HAPLOTYPECALLER(ch_haplotype_input) // Emits: [sample_id, interval, gvcf, tbi]

    // --- Step 6: Prepare for Joint Genotyping (Create Sample Map per Interval) ---
    // Check for corresponding .tbi file existence before grouping
    ch_gvcf_info = GATK_HAPLOTYPECALLER.out.gvcfs_interval
        .map { sample_id, interval, gvcf_path, tbi_path ->
            // Check if the TBI file actually exists (belt and suspenders)
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

    // --- Step 7: Joint Genotyping (per interval) ---
    ch_joint_input = CREATE_SAMPLE_MAP.out.sample_map
                         .combine(ch_ref_indexed) // Combine with ref [interval, map_path, [ref, dict, fai]]

    GATK_JOINT_GENOTYPE(ch_joint_input) // Emits: [interval, vcf_path, tbi_path]

    // --- Workflow Output ---
    // Optional: collect outputs if needed, or rely on publishDir in processes
     ch_final_vcfs = GATK_JOINT_GENOTYPE.out.filtered_vcf_interval
     ch_final_vcfs.view { interval, vcf, tbi -> "Final VCF for interval ${interval}: ${vcf}" }

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
