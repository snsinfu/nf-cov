# nf-cov

A minimal [nf-core](https://nf-co.re)-style Nextflow pipeline that produces
genome-wide coverage tracks (bigWig) from short-read FASTQ.


## Usage

Samplesheet:

```
sample,fastq_1,fastq_2,replicate,tech_replicate
WGS_HL60,DRR357125_1.fastq.gz,DRR357125_2.fastq.gz,1,1
```

Command line:

```bash
nextflow run snsinfu/nf-cov \
    -profile docker \
    --input samplesheet.csv \
    --genome GRCh38 \
    --outdir results
```

Flow:

```
samplesheet
  -> FastQC / Trim Galore
  -> bwa | bwa-mem2 | bwa-mem3 (per run)
  -> samtools merge (runs -> library)
  -> Picard MarkDuplicates (mark only)
  -> samtools view (clean filter)
  -> bamCoverage
  -> merge tech reps (.mLb) / bio reps (.mRp)
  -> bamCoverage
```


## Parameters

| Parameter                        | Default   | Description |
|----------------------------------|-----------|-------------|
| `--input`                        |           | Samplesheet CSV. |
| `--outdir`                       | `results` | Output directory. |
| `--genome`                       |           | iGenomes key (uses `conf/igenomes.config`). |
| `--igenomes_ignore`              | `false`   | Do not load the iGenomes config. |
| `--fasta`                        |           | Reference FASTA. |
| `--bwa_index`                    |           | Precomputed bwa index dir. |
| `--bwamem2_index`                |           | Precomputed bwa-mem2 index dir. |
| `--bwamem3_index`                |           | Precomputed bwa-mem3 index dir. |
| `--aligner`                      | `bwa`     | `bwa`, `bwa-mem2` or `bwa-mem3`. |
| `--effective_genome_size`        |           | RPGC effective genome size. |
| `--read_length`                  |           | Read length for catalog lookup / khmer (inferred when unset). |
| `--bin_size`                     | `1`       | bigWig bin size in bases. |
| `--normalization`                | `RPGC`    | `None`, `RPKM`, `CPM`, `BPM` or `RPGC`. |
| `--min_mapq`                     | `1`       | Minimum mapping quality for coverage. |
| `--fragment_size`                | `0`       | `0` = auto (PE estimate; SE no extension). |
| `--normalization_exclude_chroms` |           | Chromosomes excluded from normalization (default: catalog `mito_name`). |
| `--save_library`                 | `false`   | Publish per-tech-library bigWigs/BAMs. |
| `--skip_fastqc`                  | `false`   | Skip FastQC. |
| `--skip_trimming`                | `false`   | Skip TrimGalore. |
| `--skip_preseq`                  | `false`   | Skip Preseq. |

## Output

```
results/
├── <aligner>/
│   ├── merged_library/bigwig/<sample>_REP<bio>.mLb.bigWig
│   ├── merged_library/bam/<sample>_REP<bio>.mLb.bam
│   ├── merged_replicate/bigwig/<sample>.mRp.bigWig   # only for multi-replicate samples
│   ├── merged_replicate/bam/<sample>.mRp.bam
│   └── library/...                                   # only with --save_library
├── fastqc/
├── trimgalore/
├── samtools/{stats,flagstat,idxstats}/
├── picard/<library>.mkd.metrics.txt
├── preseq/
└── pipeline_info/
```

## License

MIT
