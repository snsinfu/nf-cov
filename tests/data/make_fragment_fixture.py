#!/usr/bin/env python3
"""Generate the mixed-fragment-length PE fixture sampleMix_{1,2}.fastq.gz.

Reads are simulated from tests/data/genome.fa at two known fragment lengths
(180 and 400 bp, read length 50) so the fragment-size filter can be tested by
removing one size and keeping the other. Run from anywhere:

    python3 tests/data/make_fragment_fixture.py
"""

import gzip
import os

HERE = os.path.dirname(os.path.abspath(__file__))
GENOME = os.path.join(HERE, "genome.fa")
READ_LEN = 50
FRAGMENTS = [180, 400]
PER_SIZE = 300
SPACING = 450  # > max fragment length, so fragments do not overlap
START = 100


def read_fasta(path):
    seqs = {}
    name, parts = None, []
    with open(path) as handle:
        for line in handle:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(parts)
                name, parts = line[1:].split()[0], []
            else:
                parts.append(line)
    if name is not None:
        seqs[name] = "".join(parts)
    return seqs


def revcomp(seq):
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def main():
    chrom = read_fasta(GENOME)["chr1"]
    qual = "I" * READ_LEN

    with gzip.open(os.path.join(HERE, "sampleMix_1.fastq.gz"), "wt") as r1, \
            gzip.open(os.path.join(HERE, "sampleMix_2.fastq.gz"), "wt") as r2:
        i = 0
        for size in FRAGMENTS:
            for _ in range(PER_SIZE):
                start = START + i * SPACING
                frag = chrom[start:start + size]
                if len(frag) != size:
                    raise RuntimeError("fixture runs past chr1 end")
                seq1 = frag[:READ_LEN]
                seq2 = revcomp(frag[size - READ_LEN:])
                name = "mix{}_{}".format(size, i)
                r1.write("@{}\n{}\n+\n{}\n".format(name, seq1, qual))
                r2.write("@{}\n{}\n+\n{}\n".format(name, seq2, qual))
                i += 1


if __name__ == "__main__":
    main()
