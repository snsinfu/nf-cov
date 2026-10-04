/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Per-library processing: merge runs, mark duplicates, clean filter, stats, Preseq
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SAMTOOLS_MERGE          } from '../../../modules/nf-core/samtools/merge'
include { SAMTOOLS_VIEW           } from '../../../modules/nf-core/samtools/view'
include { PICARD_MARKDUPLICATES   } from '../../../modules/nf-core/picard/markduplicates'
include { PRESEQ_LCEXTRAP         } from '../../../modules/nf-core/preseq/lcextrap'
include { BAM_STATS_SAMTOOLS      } from '../../nf-core/bam_stats_samtools'

workflow BAM_LIBRARY {
    take:
    ch_bam        // channel: [ val(meta run), path(*.bam) ]
    ch_ref        // channel: [ path(fasta), path(fai) ] (single element)
    skip_preseq   // val: boolean
    index_format  // val: 'bai' | 'csi' (BAM index format)

    main:
    //
    // The reference is provided as a single-element channel from the caller.
    //
    ch_fasta_fai = ch_ref.map { f, fi -> [ [:], f, fi ] }
    ch_merge_ref = ch_ref.map { f, fi -> [ [:], f, fi, [] ] }

    //
    // Merge all runs of a library
    //
    ch_lib_bams = ch_bam
        .map { meta, bam -> [ meta.library.toString(), meta, bam ] }
        .groupTuple(by: [0])
        .map { lib, metas, bams ->
            def m = [
                id            : lib,
                sample        : metas[0].sample,
                replicate     : metas[0].replicate,
                tech_replicate: metas[0].tech_replicate,
                library       : lib,
                single_end    : metas[0].single_end,
                level         : 'library'
            ]
            [ m, bams.flatten(), [] ]
        }

    SAMTOOLS_MERGE(ch_lib_bams, ch_merge_ref)
    ch_merged = SAMTOOLS_MERGE.out.bam

    //
    // Mark duplicates per library (mark-only; duplicates are removed by the clean filter)
    //
    PICARD_MARKDUPLICATES(ch_merged, ch_fasta_fai)
    ch_marked = PICARD_MARKDUPLICATES.out.bam

    //
    // Clean filter: primary, mapped, non-duplicate, non-supplementary reads
    //
    SAMTOOLS_VIEW(
        ch_marked.map { meta, bam -> [ meta, bam, [] ] },
        ch_fasta_fai,
        [ [], [] ],
        [ [], [] ],
        index_format
    )
    ch_clean_bam = SAMTOOLS_VIEW.out.bam
    //
    // samtools/view emits .bai and .csi as separate channels; only one is populated.
    //
    ch_clean_index = SAMTOOLS_VIEW.out.bai.mix(SAMTOOLS_VIEW.out.csi)
    ch_clean_bam_bai = ch_clean_bam.join(ch_clean_index, by: [0])

    //
    // Alignment stats
    //
    BAM_STATS_SAMTOOLS(ch_clean_bam_bai, ch_fasta_fai)

    //
    // Preseq library complexity (on the marked BAM, duplicates present)
    //
    ch_preseq = channel.empty()
    if (!skip_preseq) {
        PRESEQ_LCEXTRAP(ch_marked)
        ch_preseq = PRESEQ_LCEXTRAP.out.lc_extrap
    }

    emit:
    bam         = ch_clean_bam            // channel: [ val(meta), path(*.bam) ]
    bai         = ch_clean_index          // channel: [ val(meta), path(*.bai|*.csi) ]
    bam_bai     = ch_clean_bam_bai        // channel: [ val(meta), path(*.bam), path(*.bai|*.csi) ]
    marked_bam  = ch_marked               // channel: [ val(meta), path(*.bam) ]
    dup_metrics = PICARD_MARKDUPLICATES.out.metrics
    stats       = BAM_STATS_SAMTOOLS.out.stats
    flagstat    = BAM_STATS_SAMTOOLS.out.flagstat
    idxstats    = BAM_STATS_SAMTOOLS.out.idxstats
    preseq      = ch_preseq
}
