# Psomagen Bulk RNA-Seq Pipeline

Nextflow DSL2 pipeline for bulk RNA-seq data processed by Psomagen. Uses HISAT2 for alignment and includes Trimmomatic adapter/quality trimming for Nextera-prepared libraries.

This repo is analogous to [utdal/Bulk-RNA-Seq-Nextflow-Pipeline](https://github.com/utdal/Bulk-RNA-Seq-Nextflow-Pipeline) (STAR-based) and is consumed by the [HPC deployment framework](https://github.com/mwilde49/hpc).

## Pipeline Steps

| Step | Output Directory | Description |
|------|-----------------|-------------|
| 0 | `0_fastqc_output/` | FastQC + MultiQC quality reports (Phase 1) |
| 1 | `1_fastqc_output/` | FastQC on raw reads (Phase 2) |
| 2 | `2_trim_output/` | Trimmomatic adapter/quality trimming |
| 3 | `3_hisat2_mapping_output/` | HISAT2 splice-aware alignment |
| 3.1 | `3_1_map_metrics_output_qc/` | HISAT2 alignment summary CSV |
| 4 | `4_filter_output/` | sambamba dedup + bedtools blacklist filtering |
| 4.1 | `4_1_qualimap_filter_output_qc/` | Qualimap QC on filtered BAMs |
| 5 | `5_stringtie_counts_output/` | StringTie TPM quantification |
| 5.1 | `5_1_feature_counts_output/` | Rsubread featureCounts |
| 6 | `6_raw_counts_output/` | HTSeq raw counts matrix |
| 7 | `7_stats_output/` | Pipeline statistics and tool versions |

## Key Differences from UTDal Pipeline

- **Aligner**: HISAT2 (not STAR)
- **Trimming**: Trimmomatic with Nextera adapter removal (additional step)
- **Read suffixes**: `_1`/`_2` (Psomagen convention, not `_R1_001`/`_R2_001`)
- **Reference**: Gencode v48 (GRCh38 primary assembly)

## Usage

### Prerequisites

- Apptainer/Singularity container (see [Container](#container) below)
- HISAT2 index built from the reference genome
- Gencode v48 GTF annotation

### Running the Pipeline

1. Edit `pipeline.config` with your paths and parameters
2. List samples in `rna_seq_samples.txt` (one per line, no read suffix)
3. Run Phase 1 (QC only): set `run_fastqc = true` in config
4. Run Phase 2 (full pipeline): set `run_rna_pipeline = true` in config

```bash
nextflow run psomagen_bulk_rna_seq_pipeline.nf -c pipeline.config
```

## Container

The `container/` directory contains everything needed to build and test the Apptainer container.

### Building

```bash
cd container
./build.sh            # produces psomagen_v1.0.0.sif
./build.sh --force    # overwrite existing
```

### Testing

```bash
# Validate all tools are installed
./test_container.sh

# Run end-to-end integration test on synthetic data (~2 min)
python3 generate_test_data.py   # generate test data (if not present)
./test_integration.sh
```

### Included Tools

FastQC, MultiQC, HISAT2, Trimmomatic, Samtools, Sambamba, Bedtools, StringTie, HTSeq, Qualimap, R + Rsubread, Python + pandas, Nextflow, Java 17

## Reference Files

- `ref_human_geneid_genename_genebiotype.tsv` — 78,894 genes from Gencode v48 (gene_id, gene_name, gene_biotype)
- `NexteraPE-PE.fa` — Nextera adapter sequences for Trimmomatic
