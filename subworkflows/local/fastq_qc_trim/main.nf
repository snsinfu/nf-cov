/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Read QC and adapter trimming (FastQC + Trim Galore)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC     } from '../../../modules/nf-core/fastqc'
include { TRIMGALORE } from '../../../modules/nf-core/trimgalore'

workflow FASTQ_QC_TRIM {
    take:
    reads          // channel: [ val(meta), [ fastq_1(, fastq_2) ] ]
    skip_fastqc    // val: boolean
    skip_trimming  // val: boolean

    main:
    ch_trimmed    = reads
    ch_fastqc_zip = channel.empty()
    ch_fastqc_html = channel.empty()
    ch_trim_zip   = channel.empty()
    ch_trim_log   = channel.empty()

    if (!skip_trimming) {
        TRIMGALORE(reads)
        ch_trimmed  = TRIMGALORE.out.reads
        ch_trim_zip = TRIMGALORE.out.zip
        ch_trim_log = TRIMGALORE.out.log
    }

    if (!skip_fastqc) {
        FASTQC(ch_trimmed)
        ch_fastqc_zip  = FASTQC.out.zip
        ch_fastqc_html = FASTQC.out.html
    }

    emit:
    reads        = ch_trimmed       // channel: [ val(meta), [ fastq_1(, fastq_2) ] ]
    fastqc_zip   = ch_fastqc_zip    // channel: [ val(meta), path(*.zip) ]
    fastqc_html  = ch_fastqc_html   // channel: [ val(meta), path(*.html) ]
    trim_zip     = ch_trim_zip      // channel: [ val(meta), path(*.zip) ]
    trim_log     = ch_trim_log      // channel: [ val(meta), path(*report.txt) ]
}
