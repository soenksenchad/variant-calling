#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Define Parameters
params.reads = "samples.txt" // Input samplesheet (sample,fastq_1,fastq_2)
params.reference_genome = "$projectDir/reference/Zmays_833_Zm-B73-REFERENCE-NAM-5.0.fa" // Reference genome in reference directory
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
            def meta = [id: row.sample_id]
            def read1 = file(row.read1, checkIfExists: true)
            def read2 = file(row.read2, checkIfExists: true)
            
            [ meta, [read1, read2] ]
        }
    
    ch_intervals = Channel.fromPath(params.intervals_file).splitText().map { it.trim() }.filter { it } // Emits '1', '2', ... '10'
    
    // Reference genome files (assumes pre-indexed)
    def ref_file = file(params.reference_genome, checkIfExists: true)
    
    // Create simple reference channel with just the reference file
    ch_ref_fasta = Channel.value(ref_file)
    
    // --- Step 1: Read Trimming ---
    FASTP_TRIM(ch_reads)

    // --- Step 2: Alignment ---
    ch_align_input = FASTP_TRIM.out.trimmed_reads.combine(ch_ref_fasta)
    BWAMEM2_ALIGN(ch_align_input)

    // --- Step 3: SAMtools Processing ---
    SAMTOOLS_PROCESS(BWAMEM2_ALIGN.out.bam_files)

    // --- Step 4: Haplotype Calling (per sample, per interval) ---
    ch_haplotype_input = SAMTOOLS_PROCESS.out.processed_bam
                             .combine(ch_intervals)   // [meta, bam, interval]
                             .combine(ch_ref_fasta)   // [meta, bam, interval, ref]

    GATK_HAPLOTYPECALLER(ch_haplotype_input)

    // --- Step 5: Prepare for Joint Genotyping (Create Sample Map per Interval) ---
    ch_gvcf_info = GATK_HAPLOTYPECALLER.out.gvcfs_interval
        .map { meta, interval, gvcf_path, tbi_path ->
            // Check if the TBI file actually exists
            if (file(tbi_path).exists()) {
                return tuple(interval, tuple(meta.id, gvcf_path)) // Format for grouping: [interval, [sample_id, gvcf_path]]
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
    ch_joint_genotype_input = CREATE_SAMPLE_MAP.out.sample_map
                                 .combine(ch_ref_fasta)  // [interval, map_file, ref_file]
    
    GATK_JOINT_GENOTYPE(ch_joint_genotype_input)

    // --- Workflow Output ---
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
