/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-cov main analysis workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { samplesheetToList } from 'plugin/nf-schema'

include { PREPARE_GENOME } from '../subworkflows/local/prepare_genome'
include { FASTQ_QC_TRIM } from '../subworkflows/local/fastq_qc_trim'
include { FASTQ_ALIGN } from '../subworkflows/local/fastq_align'
include { BAM_LIBRARY } from '../subworkflows/local/bam_library'
include { BAM_MERGE_LEVELS } from '../subworkflows/local/bam_merge_levels'
include { BAM_COVERAGE } from '../subworkflows/local/bam_coverage'
include { FASTQ_READ_LENGTH } from '../modules/local/fastq_read_length'
include { buildReads; assignRuns; validateSamplesheet } from '../subworkflows/local/utils_nfcov_pipeline'

workflow NFCOV {
    take:
    fasta         // path: genome fasta
    index         // path: prebuilt aligner index (optional)
    aligner       // string: 'bwa' | 'bwa-mem2' | 'bwa-mem3'
    catalog_gsize // val: map read-length -> gsize, or null
    explicit_egs  // val: int, or null
    read_length   // val: int, or null

    main:
    //
    // Read in the samplesheet, build per-run meta and validate
    //
    def samplesheet = samplesheetToList(params.input, "${projectDir}/assets/schema_input.json")
    validateSamplesheet(samplesheet)
    ch_reads = channel.fromList(assignRuns(buildReads(samplesheet)))

    //
    // Read QC and adapter trimming
    //
    FASTQ_QC_TRIM(ch_reads, params.skip_fastqc, params.skip_trimming)
    ch_processed_reads = FASTQ_QC_TRIM.out.reads

    //
    // Read length: explicit param, otherwise inferred from the first FASTQ when the
    // effective genome size still has to be resolved.
    //
    ch_read_length = channel.empty()
    if (read_length) {
        ch_read_length = channel.value(read_length as Integer)
    }
    else if (!explicit_egs) {
        FASTQ_READ_LENGTH(ch_reads.map { meta, fastqs -> [ meta, fastqs[0] ] }.first())
        ch_read_length = FASTQ_READ_LENGTH.out.read_length
            .map { _meta, txt -> txt.text.trim() as Integer }
            .collect()
            .map { lengths -> lengths.max() }
    }

    //
    // Reference preparation
    //
    PREPARE_GENOME(fasta, index, aligner, catalog_gsize, explicit_egs, ch_read_length)

    //
    // Alignment
    //
    FASTQ_ALIGN(ch_processed_reads, PREPARE_GENOME.out.index, PREPARE_GENOME.out.fasta, aligner)

    //
    // Per-library merge, deduplication (mark-only), clean filter, stats, Preseq
    //
    ch_ref = PREPARE_GENOME.out.fasta
        .combine(PREPARE_GENOME.out.fai.map { _meta, f -> f })
        .first()
    BAM_LIBRARY(
        FASTQ_ALIGN.out.bam,
        ch_ref,
        params.skip_preseq
    )

    //
    // Replicate merges: .mLb (tech reps) and .mRp (bio reps)
    //
    BAM_MERGE_LEVELS(BAM_LIBRARY.out.bam, ch_ref, params.save_library)

    //
    // Coverage bigWigs. One track per published level; RPGC is recomputed on each BAM.
    //
    ch_ref_cov = PREPARE_GENOME.out.egs
        .combine(PREPARE_GENOME.out.fasta)
        .combine(PREPARE_GENOME.out.fai.map { _meta, f -> f })
        .first()
    ch_cover_bams = channel.empty()
        .mix(BAM_MERGE_LEVELS.out.mLb)
        .mix(BAM_MERGE_LEVELS.out.mRp)
    if (params.save_library) {
        ch_cover_bams = ch_cover_bams.mix(BAM_LIBRARY.out.bam_bai)
    }
    BAM_COVERAGE(ch_cover_bams, ch_ref_cov)

    //
    // Collate and save software versions. Modules emit versions via the 'versions' topic.
    //
    def topic_versions = channel.topic('versions')
        .map { process, tool, version -> [ process[process.lastIndexOf(':') + 1..-1], "  ${tool}: ${version}" ] }
        .groupTuple(by: [0])
        .map { process, tool_versions ->
            "${process}:\n${tool_versions.unique().sort().join('\n')}"
        }

    topic_versions
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf-cov_software_mqc_versions.yml',
            sort: true,
            newLine: true
        )

    emit:
    clean_bam   = BAM_LIBRARY.out.bam          // channel: [ val(meta), path(*.bam) ]
    clean_bai   = BAM_LIBRARY.out.bai          // channel: [ val(meta), path(*.bai) ]
    mLb         = BAM_MERGE_LEVELS.out.mLb     // channel: [ val(meta), path(*.mLb.bam) ]
    mRp         = BAM_MERGE_LEVELS.out.mRp     // channel: [ val(meta), path(*.mRp.bam) ]
    bigwig      = BAM_COVERAGE.out.bigwig      // channel: [ val(meta), path(*.bigWig) ]
    dup_metrics = BAM_LIBRARY.out.dup_metrics  // channel: [ val(meta), path(*.metrics.txt) ]
    stats       = BAM_LIBRARY.out.stats
    flagstat    = BAM_LIBRARY.out.flagstat
    idxstats    = BAM_LIBRARY.out.idxstats
    preseq      = BAM_LIBRARY.out.preseq
    bam   = FASTQ_ALIGN.out.bam       // channel: [ val(meta), path(*.bam) ]
    reads = ch_processed_reads        // channel: [ val(meta), [ fastq_1(, fastq_2) ] ]
    fastqc_zip  = FASTQ_QC_TRIM.out.fastqc_zip
    fastqc_html = FASTQ_QC_TRIM.out.fastqc_html
    trim_log    = FASTQ_QC_TRIM.out.trim_log
    fasta = PREPARE_GENOME.out.fasta   // channel: path(genome.fa)
    fai   = PREPARE_GENOME.out.fai     // channel: [ val(meta), path(genome.fa.fai) ]
    index = PREPARE_GENOME.out.index   // channel: [ val(meta), path(index directory) ]
    egs   = PREPARE_GENOME.out.egs     // channel: val(effective genome size)
}
