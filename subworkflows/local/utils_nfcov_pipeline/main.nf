/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-cov utility functions (samplesheet parsing, meta building, validation)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Build one reads item per samplesheet row: [ meta, [ fastq_1(, fastq_2) ] ]
//
def buildReads(samplesheet) {
    def reads = []
    samplesheet.each { item ->
        def meta = item[0]
        def fastq_1 = item[1]
        def fastq_2 = item[2]

        def tech = (meta.tech_replicate ?: 1) as Integer
        def replicate = meta.replicate as Integer
        def library = "${meta.id}_REP${replicate}_T${tech}".toString()
        def single_end = !fastq_2

        def m = [
            id            : meta.id,
            sample        : meta.id,
            replicate     : replicate,
            tech_replicate: tech,
            library       : library,
            single_end    : single_end,
            level         : 'library'
        ]
        def fastqs = single_end
            ? [ resolveInputPath(fastq_1) ]
            : [ resolveInputPath(fastq_1), resolveInputPath(fastq_2) ]
        reads << [ m, fastqs ]
    }
    return reads
}

//
// Samplesheet FASTQ paths are normally absolute; resolve relative paths against the
// pipeline root so committed test samplesheets work from any working directory.
//
def resolveInputPath(p) {
    def s = p.toString()
    if (s.startsWith('/') || s ==~ /^[A-Za-z]:[\\\/].*/) {
        return file(s)
    }
    return file("${projectDir}/${s}")
}

//
// Assign a run index and a unique read group to every sequencing run. Rows sharing
// (sample, replicate, tech_replicate) become multiple runs of the same library.
//
def assignRuns(reads) {
    def grouped = reads.groupBy { item -> item[0].library }
    def out = []
    grouped.each { _library, items ->
        items.eachWithIndex { item, idx ->
            def run = idx + 1
            def base = item[0]
            def read_group = "@RG\\tID:${base.library}.run${run}\\tSM:${base.sample}\\tLB:${base.library}\\tPL:ILLUMINA".toString()
            def meta = base + [ run: run, read_group: read_group ]
            out << [ meta, item[1] ]
        }
    }
    return out
}

//
// Validate the raw samplesheet list before any process runs.
//
def validateSamplesheet(samplesheet) {
    // Unique fastq_1
    def seen = [:] as Map
    samplesheet.each { _meta, fastq_1, _fastq_2 ->
        def f1 = fastq_1.toString()
        if (seen.containsKey(f1)) {
            error("Duplicate fastq_1 in samplesheet: ${f1}")
        }
        seen[f1] = true
    }

    // Consistent SE/PE layout within a library
    def lib_layout = [:] as Map
    samplesheet.each { meta, _fastq_1, fastq_2 ->
        def tech = (meta.tech_replicate ?: 1)
        def library = "${meta.id}_REP${meta.replicate}_T${tech}".toString()
        def layout = fastq_2 ? 'PE' : 'SE'
        if (lib_layout.containsKey(library) && lib_layout[library] != layout) {
            error("Library '${library}' mixes single-end and paired-end runs.")
        }
        lib_layout[library] = layout
    }

    // Consistent SE/PE layout within a sample (across biological replicates)
    def sample_layout = [:] as Map
    samplesheet.each { meta, _fastq_1, fastq_2 ->
        def layout = fastq_2 ? 'PE' : 'SE'
        if (sample_layout.containsKey(meta.id) && sample_layout[meta.id] != layout) {
            error("Sample '${meta.id}' mixes single-end and paired-end replicates.")
        }
        sample_layout[meta.id] = layout
    }
}

//
// Largest contig length recorded in a samtools faidx index. The .fai is plain text:
// NAME<TAB>LENGTH<TAB>OFFSET<TAB>LINEBASES<TAB>LINEWIDTH.
//
def maxContigLength(fai) {
    fai.readLines()
        .findAll { it?.trim() }
        .collect { it.split('\t')[1] as Long }
        .max() ?: 0L
}

//
// The BAI format stores positions as 2^29-byte (512 Mbp) blocks and cannot address a contig
// >= 2^29 bp. Fail early with an actionable message rather than letting samtools fail late.
//
def validateIndexFormat(Long maxLen, String format) {
    if (format == 'bai' && maxLen >= (1L << 29)) {
        error("Reference has a contig of ${maxLen} bp (>= 512 Mbp / 2^29), which the BAI " +
            "index format cannot address. Re-run with --bam_index_format csi.")
    }
    return true
}

