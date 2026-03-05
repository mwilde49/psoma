# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Psomagen Bulk RNA-Seq Nextflow pipeline — analogous to [utdal/Bulk-RNA-Seq-Nextflow-Pipeline](https://github.com/utdal/Bulk-RNA-Seq-Nextflow-Pipeline) but using HISAT2 instead of STAR, with an added Trimmomatic trimming step for Nextera-prepared libraries. This repo combines the pipeline AND the Apptainer container (analogous to combining utdal + [mwilde49/bulkseq](https://github.com/mwilde49/bulkseq)). Consumed by [mwilde49/hpc](https://github.com/mwilde49/hpc).

## Architecture

- **`psomagen_bulk_rna_seq_pipeline.nf`** — Main Nextflow DSL2 workflow with 10 processes. HISAT2 mapping is inlined; other processes delegate to external shell scripts.
- **`pipeline.config`** — All parameters (paths, threads, Trimmomatic settings, phase control flags).
- **Shell scripts** (`*.sh`) — Called by Nextflow processes. Each handles one pipeline step.
- **`container/`** — Apptainer container definition, build script, and test suite.
  - `psomagen.def` — Container definition (Ubuntu 22.04 + Miniforge/mamba).
  - `build.sh` — Builds `psomagen_v<VERSION>.sif`.
  - `test_container.sh` — Validates all tools are installed (5 phases).
  - `test_integration.sh` — End-to-end 9-step test on synthetic data.
  - `generate_test_data.py` — Deterministic (seed=42) synthetic data generator.

## Key Commands

```bash
# Build container
cd container && ./build.sh

# Test container tools
./test_container.sh [path/to/psomagen_v*.sif]

# Generate synthetic test data
python3 generate_test_data.py

# Run integration tests
./test_integration.sh [path/to/psomagen_v*.sif]

# Run pipeline
nextflow run psomagen_bulk_rna_seq_pipeline.nf -c pipeline.config
```

## Important Conventions

- Output directories numbered 0–7 (one more than UTDal due to Trimmomatic step at position 2).
- Read suffixes: `_1`/`_2` (not `_R1_001`/`_R2_001`).
- BAM naming: `sample.sorted.bam` (HISAT2), not `sampleAligned.sortedByCoord.out.bam` (STAR).
- Reference: Gencode v48 (GRCh38). GTF uses `gene_type` attribute (not `gene_biotype`).
- Container uses bioconda `trimmomatic` wrapper (not `java -jar`).
- Strandedness: `--rf` (reverse-forward, matching dUTP/Nextera library prep).

## Ecosystem

Pipeline repo (this) → Container (.sif) → HPC framework (mwilde49/hpc) deploys via `pipeline_manager.sh`.
