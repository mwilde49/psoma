#!/bin/bash
set -euo pipefail

# ── End-to-end Nextflow pipeline test ──
# Runs the actual psomagen_bulk_rna_seq_pipeline.nf on synthetic test data
# inside the Apptainer container, validating the full Nextflow orchestration.
# Usage: ./test_nextflow_pipeline.sh [path/to/psomagen_v*.sif]

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
TEST_DATA="${SCRIPT_DIR}/test_data"

# Auto-detect newest .sif or use argument
if [[ -n "${1:-}" ]]; then
    SIF="$1"
else
    SIF="$(ls -t "${SCRIPT_DIR}"/psomagen_v*.sif 2>/dev/null | head -1)"
    if [[ -z "$SIF" ]]; then
        SIF="$(ls -t psomagen_v*.sif 2>/dev/null | head -1)"
    fi
fi

if [[ -z "$SIF" || ! -f "$SIF" ]]; then
    echo -e "${RED}ERROR: No .sif file found.${NC}"
    exit 1
fi

echo -e "${YELLOW}Nextflow pipeline test using: ${SIF}${NC}"
echo -e "${YELLOW}Test data: ${TEST_DATA}${NC}"
echo ""

# Verify test data exists
for f in genome.fa genes.gtf sample1_1.fastq.gz sample1_2.fastq.gz samples.txt NexteraPE-PE.fa ref_human_geneid_genename_genebiotype.tsv; do
    if [[ ! -f "${TEST_DATA}/${f}" ]]; then
        echo -e "${RED}ERROR: Missing test data file: ${TEST_DATA}/${f}${NC}"
        echo "Run: python3 generate_test_data.py"
        exit 1
    fi
done

# Create temp working directory (this becomes config_directory)
WORKDIR=$(mktemp -d)
trap "rm -rf $WORKDIR" EXIT

echo "Working directory: $WORKDIR"
echo ""

# ── Set up the working directory ──
# Copy all pipeline scripts
echo -n "  Setting up pipeline scripts... "
for f in \
    psomagen_bulk_rna_seq_pipeline.nf \
    trimmomatic_trimming.sh \
    dedup_and_filtering.sh \
    generate_fastqc_reports.sh \
    generate_map_metrics.py \
    generate_qualimap_report.sh \
    generate_raw_counts.sh \
    generate_stringtie_counts.sh \
    generate_stats.sh \
    fetch_genename_genebiotype_for_counts.py \
    rsubread_featurecount_script.R; do
    cp "${REPO_DIR}/${f}" "${WORKDIR}/"
done
echo -e "${GREEN}OK${NC}"

# Copy test reference files
echo -n "  Copying reference files... "
cp "${TEST_DATA}/genes.gtf" "${WORKDIR}/genes.gtf"
cp "${TEST_DATA}/NexteraPE-PE.fa" "${WORKDIR}/NexteraPE-PE.fa"
cp "${TEST_DATA}/ref_human_geneid_genename_genebiotype.tsv" "${WORKDIR}/ref_human_geneid_genename_genebiotype.tsv"
echo -e "${GREEN}OK${NC}"

# Create sample manifest
echo "sample1" > "${WORKDIR}/rna_seq_samples.txt"

# Create FASTQ directory with test reads
mkdir -p "${WORKDIR}/fastq"
cp "${TEST_DATA}/sample1_1.fastq.gz" "${WORKDIR}/fastq/"
cp "${TEST_DATA}/sample1_2.fastq.gz" "${WORKDIR}/fastq/"

# Create minimal BED files (single entry on a dummy chromosome — won't filter anything)
echo -e "chr_dummy\t0\t1" > "${WORKDIR}/blacklist.bed"
echo -e "chr_dummy\t0\t1" > "${WORKDIR}/filter.bed"

# Create a writable home for Nextflow (it needs ~/.nextflow)
NF_HOME="${WORKDIR}/.nf_home"
mkdir -p "${NF_HOME}"

# ── Build HISAT2 index from test genome ──
echo -n "  Building HISAT2 index... "
if apptainer exec --cleanenv --env PYTHONNOUSERSITE=1 \
    --bind "${TEST_DATA}:${TEST_DATA}" \
    --bind "${WORKDIR}:${WORKDIR}" \
    "$SIF" hisat2-build -p 2 \
    "${TEST_DATA}/genome.fa" \
    "${WORKDIR}/test_index" &>"${WORKDIR}/hisat2_build.log"; then
    echo -e "${GREEN}OK${NC}"
else
    echo -e "${RED}FAIL${NC}"
    echo "  Log: ${WORKDIR}/hisat2_build.log"
    tail -5 "${WORKDIR}/hisat2_build.log" 2>/dev/null || true
    exit 1
fi

# ── Generate Nextflow test config ──
cat > "${WORKDIR}/test_pipeline.config" << CFGEOF
params.proj_name = 'Integration-Test'
params.species = 'Human'

params.config_directory = '${WORKDIR}'
params.fastq_files = '${WORKDIR}/fastq'
params.samples_file = '${WORKDIR}/rna_seq_samples.txt'

params.fastqc_cores = 2

params.hisat2_index = '${WORKDIR}/test_index'
params.reference_gtf = '${WORKDIR}/genes.gtf'

params.exclude_bed_file_path = '${WORKDIR}/filter.bed'
params.blacklist_bed_file_path = '${WORKDIR}/blacklist.bed'

params.paired_end = true
params.read1_suffix = "_1"
params.read2_suffix = "_2"

params.illumina_clip_file = '${WORKDIR}/NexteraPE-PE.fa'
params.headcrop = 10
params.leading = 3
params.trailing = 3
params.slidingwindow = '4:15'
params.minlen = 36
params.illuminaclip_params = '2:30:10:5:true'

params.strand_st = '--rf'
params.strand_hts = 'reverse'
params.paired_hts = 'pos'

params.run_fastqc = false
params.run_rna_pipeline = true
CFGEOF

echo ""
echo -e "${YELLOW}Running Nextflow pipeline...${NC}"
echo ""

# ── Run the pipeline ──
if apptainer run --cleanenv \
    --env PYTHONNOUSERSITE=1 \
    --env HOME="${NF_HOME}" \
    --bind "${WORKDIR}:${WORKDIR}" \
    --bind "${TEST_DATA}:${TEST_DATA}" \
    --bind "${NF_HOME}:${NF_HOME}" \
    "$SIF" run "${WORKDIR}/psomagen_bulk_rna_seq_pipeline.nf" \
    -c "${WORKDIR}/test_pipeline.config" \
    &>"${WORKDIR}/nextflow_run.log"; then
    echo -e "${GREEN}Nextflow pipeline completed successfully${NC}"
else
    echo -e "${RED}Nextflow pipeline FAILED${NC}"
    echo ""
    echo "Last 30 lines of log:"
    tail -30 "${WORKDIR}/nextflow_run.log" 2>/dev/null || true
    echo ""
    echo "Full log: ${WORKDIR}/nextflow_run.log"
    # Don't exit yet — still check what was produced
fi

echo ""

# ── Validate outputs ──
echo -e "${YELLOW}Validating outputs...${NC}"
PASS=0
FAIL=0

check_output() {
    local desc="$1"
    local path="$2"
    echo -n "  ${desc}... "
    if [[ -e "$path" ]]; then
        if [[ -d "$path" ]]; then
            local count=$(find "$path" -type f | wc -l)
            echo -e "${GREEN}PASS${NC} (${count} files)"
        elif [[ -s "$path" ]]; then
            echo -e "${GREEN}PASS${NC} ($(wc -l < "$path") lines)"
        else
            echo -e "${RED}FAIL (empty)${NC}"
            FAIL=$((FAIL + 1))
            return
        fi
        PASS=$((PASS + 1))
    else
        echo -e "${RED}FAIL (missing)${NC}"
        FAIL=$((FAIL + 1))
    fi
}

check_output "Trimmomatic output"          "${WORKDIR}/2_trim_output/sample1_1.trim.fastq.gz"
check_output "Trimmomatic R2 output"       "${WORKDIR}/2_trim_output/sample1_2.trim.fastq.gz"
check_output "HISAT2 sorted BAM"           "${WORKDIR}/3_hisat2_mapping_output/sample1.sorted.bam"
check_output "HISAT2 BAM index"            "${WORKDIR}/3_hisat2_mapping_output/sample1.sorted.bam.bai"
check_output "HISAT2 summary report"       "${WORKDIR}/3_hisat2_mapping_output/log/sample1_h2report.txt"
check_output "Mapping metrics CSV"         "${WORKDIR}/3_1_map_metrics_output_qc/hisat2_summary.csv"
check_output "Dedup BAM"                   "${WORKDIR}/4_filter_output/sample1.dedup.bam"
check_output "Filtered BAM"                "${WORKDIR}/4_filter_output/sample1.filt.bam"
check_output "Filtered BAM index"          "${WORKDIR}/4_filter_output/sample1.filt.bam.bai"
check_output "Qualimap output dir"         "${WORKDIR}/4_1_qualimap_filter_output_qc"
check_output "StringTie assembled GTF"     "${WORKDIR}/5_stringtie_counts_output/sample1.assemb.gtf"
check_output "StringTie tab"               "${WORKDIR}/5_stringtie_counts_output/sample1.tab"
check_output "StringTie merged TPM"        "${WORKDIR}/5_stringtie_counts_output/genes.tpm.txt"
check_output "featureCounts CSV"           "${WORKDIR}/6_raw_counts_output/raw_feature_counts.csv"
check_output "HTSeq raw counts"            "${WORKDIR}/6_raw_counts_output/raw_htseq_counts.csv"
check_output "Pipeline stats log"          "$(ls ${WORKDIR}/7_pipeline_stats_*.log 2>/dev/null | head -1)"

echo ""
echo -e "${GREEN}PASSED: ${PASS}${NC}  ${RED}FAILED: ${FAIL}${NC}"

if [[ $FAIL -gt 0 ]]; then
    echo ""
    echo "Nextflow log: ${WORKDIR}/nextflow_run.log"
    exit 1
fi
