/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-cov main analysis workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { samplesheetToList } from 'plugin/nf-schema'

include { PREPARE_GENOME } from '../subworkflows/local/prepare_genome'
include { FASTQ_QC_TRIM } from '../subworkflows/local/fastq_qc_trim'
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

    emit:
    reads = ch_processed_reads        // channel: [ val(meta), [ fastq_1(, fastq_2) ] ]
    fastqc_zip  = FASTQ_QC_TRIM.out.fastqc_zip
    fastqc_html = FASTQ_QC_TRIM.out.fastqc_html
    trim_log    = FASTQ_QC_TRIM.out.trim_log
    fasta = PREPARE_GENOME.out.fasta   // channel: path(genome.fa)
    fai   = PREPARE_GENOME.out.fai     // channel: [ val(meta), path(genome.fa.fai) ]
    index = PREPARE_GENOME.out.index   // channel: [ val(meta), path(index directory) ]
    egs   = PREPARE_GENOME.out.egs     // channel: val(effective genome size)
}
