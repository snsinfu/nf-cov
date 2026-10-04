/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Prepare reference genome files (fai, aligner index, effective genome size)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BWA_INDEX          } from '../../../modules/nf-core/bwa/index'
include { BWAMEM2_INDEX      } from '../../../modules/nf-core/bwamem2/index'
include { BWAMEM3_INDEX      } from '../../../modules/nf-core/bwamem3/index'
include { SAMTOOLS_FAIDX     } from '../../../modules/nf-core/samtools/faidx'
include { KHMER_UNIQUEKMERS  } from '../../../modules/nf-core/khmer/uniquekmers'

workflow PREPARE_GENOME {
    take:
    fasta         // path: genome fasta
    index         // path: prebuilt index directory for the selected aligner (optional)
    aligner       // string: 'bwa' | 'bwa-mem2' | 'bwa-mem3'
    catalog_gsize // val: map read-length(string) -> effective genome size, or null
    explicit_egs  // val: int, or null
    read_length   // val: int, or null

    main:
    ch_fasta = channel.value(file(fasta, checkIfExists: true))

    //
    // FASTA index
    //
    SAMTOOLS_FAIDX(ch_fasta.map { item -> [ [:], item, [] ] }, false)
    ch_fai = SAMTOOLS_FAIDX.out.fai

    //
    // Build the aligner index if one was not supplied
    //
    if (!index) {
        if (aligner == 'bwa') {
            BWA_INDEX(ch_fasta.map { item -> [ [:], item ] })
            ch_index = BWA_INDEX.out.index
        }
        else if (aligner == 'bwa-mem2') {
            BWAMEM2_INDEX(ch_fasta.map { item -> [ [:], item ] })
            ch_index = BWAMEM2_INDEX.out.index
        }
        else if (aligner == 'bwa-mem3') {
            BWAMEM3_INDEX(ch_fasta.map { item -> [ [:], item ] })
            ch_index = BWAMEM3_INDEX.out.index
        }
        else {
            error("Invalid --aligner '${aligner}'. Use 'bwa', 'bwa-mem2' or 'bwa-mem3'.")
        }
    }
    else {
        ch_index = channel.value([ [:], file(index, checkIfExists: true) ])
    }

    //
    // Effective genome size for RPGC normalization.
    //   1. explicit --effective_genome_size
    //   2. catalog macs_gsize, keyed by read length (nearest key)
    //   3. khmer unique-kmers estimate (requires read length)
    //
    ch_egs = channel.empty()
    if (explicit_egs) {
        ch_egs = channel.value(explicit_egs as Long)
    }
    else {
        def chosen = null
        if (catalog_gsize && read_length) {
            def keys = catalog_gsize.keySet().collect { it.toString() as Integer }.sort()
            def best = null
            def best_dist = Integer.MAX_VALUE
            keys.each { k ->
                def d = Math.abs(k - (read_length as Integer))
                if (d < best_dist) {
                    best_dist = d
                    best = k
                }
            }
            if (best != null) {
                chosen = [key: best.toString(), value: catalog_gsize[best.toString()] as Long]
            }
        }

        if (chosen != null) {
            log.info "[nf-cov] Using catalog effective genome size ${chosen.value} (read length key '${chosen.key}')"
            ch_egs = channel.value(chosen.value)
        }
        else if (read_length) {
            KHMER_UNIQUEKMERS(ch_fasta.map { item -> [ [:], item ] }, read_length)
            ch_egs = KHMER_UNIQUEKMERS.out.kmers.map { _meta, kmers ->
                def txt = kmers.text.trim()
                txt ? txt.toLong() : 0L
            }
        }
        else {
            error("Effective genome size is required. Pass --effective_genome_size, use a --genome with a catalog 'macs_gsize', or pass --read_length so khmer can estimate it.")
        }
    }

    emit:
    fasta = ch_fasta    // channel: path(genome.fa)
    fai   = ch_fai      // channel: [ val(meta), path(genome.fa.fai) ]
    index = ch_index    // channel: [ val(meta), path(index directory) ]
    egs   = ch_egs      // channel: val(effective genome size)
}
