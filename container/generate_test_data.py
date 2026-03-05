#!/usr/bin/env python3
"""
Generate synthetic test data for the Psomagen Bulk RNA-Seq pipeline integration tests.

Creates:
  - genome.fa          : 10 kb synthetic genome (1 chromosome)
  - genes.gtf          : 3 genes (2 protein_coding, 1 lncRNA) with 2 exons each
  - sample1_1.fastq.gz : 200 paired-end reads, R1 (100 bp, fragment length 250)
  - sample1_2.fastq.gz : 200 paired-end reads, R2
  - samples.txt        : Sample manifest ("sample1")
  - NexteraPE-PE.fa    : Nextera adapter sequence for Trimmomatic
  - ref_human_geneid_genename_genebiotype.tsv : Gene reference for annotation

All data is deterministic (seed=42) for reproducibility.
Reads are extracted from exonic regions to guarantee HISAT2 alignment and nonzero counts.
"""

import gzip
import os
import random

random.seed(42)

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "test_data")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ── Parameters ──
GENOME_LEN = 10000
CHROM = "chr_test"
READ_LEN = 100
FRAGMENT_LEN = 250
NUM_READS = 200
QUALITY_CHAR = "I"  # Phred 40

# ── Gene definitions ──
# Each gene: (gene_id, gene_name, gene_type, exons=[(start, end), ...], strand)
GENES = [
    ("ENSG00000000001.1", "TEST_GENE_A", "protein_coding", [(500, 800), (900, 1200)], "+"),
    ("ENSG00000000002.1", "TEST_GENE_B", "protein_coding", [(2500, 2900), (3100, 3500)], "-"),
    ("ENSG00000000003.1", "TEST_LNCRNA", "lncRNA", [(5500, 5800), (6000, 6300)], "+"),
]


def random_seq(length):
    """Generate a random DNA sequence with ~40% GC content."""
    weights = [0.3, 0.2, 0.2, 0.3]  # A, C, G, T
    return "".join(random.choices("ACGT", weights=weights, k=length))


def reverse_complement(seq):
    comp = str.maketrans("ACGT", "TGCA")
    return seq.translate(comp)[::-1]


def write_genome():
    """Write a synthetic genome FASTA."""
    seq = random_seq(GENOME_LEN)
    path = os.path.join(OUTPUT_DIR, "genome.fa")
    with open(path, "w") as f:
        f.write(f">{CHROM}\n")
        for i in range(0, len(seq), 80):
            f.write(seq[i : i + 80] + "\n")
    return seq


def write_gtf():
    """Write a GTF annotation file with 3 genes, each with 2 exons."""
    path = os.path.join(OUTPUT_DIR, "genes.gtf")
    with open(path, "w") as f:
        for gene_id, gene_name, gene_type, exons, strand in GENES:
            gene_start = exons[0][0]
            gene_end = exons[-1][1]
            tx_id = gene_id.replace("ENSG", "ENST")

            attrs = (
                f'gene_id "{gene_id}"; transcript_id "{tx_id}"; '
                f'gene_type "{gene_type}"; gene_name "{gene_name}";'
            )

            f.write(f"{CHROM}\ttest\tgene\t{gene_start}\t{gene_end}\t.\t{strand}\t.\t{attrs}\n")
            f.write(f"{CHROM}\ttest\ttranscript\t{gene_start}\t{gene_end}\t.\t{strand}\t.\t{attrs}\n")

            for idx, (ex_start, ex_end) in enumerate(exons, 1):
                f.write(
                    f"{CHROM}\ttest\texon\t{ex_start}\t{ex_end}\t.\t{strand}\t.\t"
                    f'{attrs} exon_number {idx};\n'
                )


def write_reads(genome_seq):
    """Generate paired-end reads from exonic regions."""
    r1_path = os.path.join(OUTPUT_DIR, "sample1_1.fastq.gz")
    r2_path = os.path.join(OUTPUT_DIR, "sample1_2.fastq.gz")

    # Collect all exonic regions
    exonic_regions = []
    for _, _, _, exons, strand in GENES:
        for start, end in exons:
            if end - start >= READ_LEN:
                exonic_regions.append((start - 1, end, strand))  # 0-based

    r1_records = []
    r2_records = []

    for i in range(NUM_READS):
        region = random.choice(exonic_regions)
        reg_start, reg_end, strand = region

        # Pick a position within the exonic region
        max_start = reg_end - READ_LEN
        if max_start <= reg_start:
            pos = reg_start
        else:
            pos = random.randint(reg_start, max_start)

        frag_seq = genome_seq[pos : pos + FRAGMENT_LEN]
        if len(frag_seq) < READ_LEN * 2:
            frag_seq = genome_seq[pos : pos + READ_LEN * 2]

        r1_seq = frag_seq[:READ_LEN]
        r2_seq = reverse_complement(frag_seq[-READ_LEN:]) if len(frag_seq) >= READ_LEN * 2 else reverse_complement(frag_seq[READ_LEN:])

        if len(r1_seq) < READ_LEN:
            r1_seq += random_seq(READ_LEN - len(r1_seq))
        if len(r2_seq) < READ_LEN:
            r2_seq += random_seq(READ_LEN - len(r2_seq))

        qual = QUALITY_CHAR * READ_LEN
        r1_records.append(f"@read{i+1}/1\n{r1_seq}\n+\n{qual}\n")
        r2_records.append(f"@read{i+1}/2\n{r2_seq}\n+\n{qual}\n")

    with gzip.open(r1_path, "wt") as f:
        f.writelines(r1_records)
    with gzip.open(r2_path, "wt") as f:
        f.writelines(r2_records)


def write_samples():
    """Write the sample manifest."""
    path = os.path.join(OUTPUT_DIR, "samples.txt")
    with open(path, "w") as f:
        f.write("sample1\n")


def write_adapter():
    """Write the Nextera adapter file for Trimmomatic."""
    path = os.path.join(OUTPUT_DIR, "NexteraPE-PE.fa")
    with open(path, "w") as f:
        f.write(">Nextera_SE\n")
        f.write("CTGTCTCTTATACACATCT\n")


def write_ref_tsv():
    """Write the gene ID/name/biotype reference TSV."""
    path = os.path.join(OUTPUT_DIR, "ref_human_geneid_genename_genebiotype.tsv")
    with open(path, "w") as f:
        for gene_id, gene_name, gene_type, _, _ in GENES:
            f.write(f"{gene_id} {gene_name} {gene_type}\n")


if __name__ == "__main__":
    print("Generating synthetic test data...")
    genome = write_genome()
    print(f"  genome.fa          : {GENOME_LEN} bp, 1 chromosome ({CHROM})")
    write_gtf()
    print(f"  genes.gtf          : {len(GENES)} genes, {sum(len(g[3]) for g in GENES)} exons")
    write_reads(genome)
    print(f"  sample1_1.fastq.gz : {NUM_READS} PE reads, R1 ({READ_LEN} bp)")
    print(f"  sample1_2.fastq.gz : {NUM_READS} PE reads, R2 ({READ_LEN} bp)")
    write_samples()
    print(f"  samples.txt        : 1 sample")
    write_adapter()
    print(f"  NexteraPE-PE.fa    : Nextera adapter")
    write_ref_tsv()
    print(f"  ref_human_*.tsv    : {len(GENES)} genes")
    print("Done.")
