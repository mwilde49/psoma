#!/bin/bash
set -euo pipefail

# ── End-to-end integration test for the Psomagen pipeline container ──
# Runs 9 pipeline steps on synthetic test data (~2 min)
# Usage: ./test_integration.sh [path/to/psomagen_v*.sif]

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
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

echo -e "${YELLOW}Integration test using: ${SIF}${NC}"
echo -e "${YELLOW}Test data: ${TEST_DATA}${NC}"
echo ""

# Verify test data exists
for f in genome.fa genes.gtf sample1_1.fastq.gz sample1_2.fastq.gz samples.txt NexteraPE-PE.fa; do
    if [[ ! -f "${TEST_DATA}/${f}" ]]; then
        echo -e "${RED}ERROR: Missing test data file: ${TEST_DATA}/${f}${NC}"
        echo "Run: python generate_test_data.py to regenerate test data."
        exit 1
    fi
done

# Create temp working directory
WORKDIR=$(mktemp -d)
trap "rm -rf $WORKDIR" EXIT

echo "Working directory: $WORKDIR"
echo ""

PASS=0
FAIL=0

run_step() {
    local step_num="$1"
    local step_name="$2"
    shift 2
    echo -n "  Step ${step_num}: ${step_name}... "
    if apptainer exec --cleanenv --env PYTHONNOUSERSITE=1 \
        --bind "${TEST_DATA}:${TEST_DATA}" \
        --bind "${WORKDIR}:${WORKDIR}" \
        "$SIF" "$@" &>"${WORKDIR}/step${step_num}.log"; then
        echo -e "${GREEN}PASS${NC}"
        PASS=$((PASS + 1))
    else
        echo -e "${RED}FAIL${NC}"
        echo "  Log: ${WORKDIR}/step${step_num}.log"
        tail -5 "${WORKDIR}/step${step_num}.log" 2>/dev/null || true
        FAIL=$((FAIL + 1))
    fi
}

# ── Step 1: Build HISAT2 index from tiny genome ──
run_step 1 "HISAT2 genome index" \
    hisat2-build \
    -p 2 \
    "${TEST_DATA}/genome.fa" \
    "${WORKDIR}/test_index"

# ── Step 2: Trimmomatic trimming ──
run_step 2 "Trimmomatic trimming" \
    trimmomatic PE -phred33 -threads 2 \
    "${TEST_DATA}/sample1_1.fastq.gz" \
    "${TEST_DATA}/sample1_2.fastq.gz" \
    "${WORKDIR}/sample1_1.trim.fastq.gz" \
    "${WORKDIR}/sample1_1.trim_unpaired.fastq.gz" \
    "${WORKDIR}/sample1_2.trim.fastq.gz" \
    "${WORKDIR}/sample1_2.trim_unpaired.fastq.gz" \
    "ILLUMINACLIP:${TEST_DATA}/NexteraPE-PE.fa:2:30:10:5:true" \
    HEADCROP:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# ── Step 3: HISAT2 alignment ──
run_step 3 "HISAT2 alignment" \
    bash -c "hisat2 -p 2 -x ${WORKDIR}/test_index \
        -1 ${WORKDIR}/sample1_1.trim.fastq.gz \
        -2 ${WORKDIR}/sample1_2.trim.fastq.gz \
        -S ${WORKDIR}/sample1.sam \
        --summary-file ${WORKDIR}/sample1_h2report.txt && \
    samtools sort -@ 2 -o ${WORKDIR}/sample1.sorted.bam ${WORKDIR}/sample1.sam && \
    rm ${WORKDIR}/sample1.sam"

# ── Step 4: samtools index ──
run_step 4 "samtools index" \
    samtools index "${WORKDIR}/sample1.sorted.bam"

# ── Step 5: sambamba dedup ──
run_step 5 "sambamba markdup" \
    bash -c "sambamba markdup -r -t 2 ${WORKDIR}/sample1.sorted.bam ${WORKDIR}/sample1.dedup.bam && \
    samtools index ${WORKDIR}/sample1.dedup.bam"

# ── Step 6: StringTie ──
run_step 6 "StringTie quantification" \
    stringtie -p 2 -e -G "${TEST_DATA}/genes.gtf" --rf \
    "${WORKDIR}/sample1.dedup.bam" \
    -o "${WORKDIR}/sample1.assemb.gtf" \
    -A "${WORKDIR}/sample1.tab"

# ── Step 7: HTSeq count ──
run_step 7 "HTSeq count" \
    bash -c "samtools sort -o ${WORKDIR}/sample1.pos.bam ${WORKDIR}/sample1.dedup.bam && \
    samtools index ${WORKDIR}/sample1.pos.bam && \
    htseq-count -f bam -s reverse -r pos ${WORKDIR}/sample1.pos.bam ${TEST_DATA}/genes.gtf > ${WORKDIR}/sample1_counts.tsv"

# ── Step 8: featureCounts (Rsubread) ──
run_step 8 "featureCounts (Rsubread)" \
    Rscript -e "library(Rsubread); fc <- featureCounts(files='${WORKDIR}/sample1.dedup.bam', annot.ext='${TEST_DATA}/genes.gtf', isGTFAnnotationFile=TRUE, isPairedEnd=TRUE, nthreads=2); write.csv(fc\$counts, file='${WORKDIR}/feature_counts.csv')"

# ── Step 9: MultiQC ──
run_step 9 "MultiQC" \
    bash -c "cd ${WORKDIR} && multiqc ."

echo ""
echo -e "${GREEN}PASSED: ${PASS}${NC}  ${RED}FAILED: ${FAIL}${NC}"

if [[ $FAIL -gt 0 ]]; then
    exit 1
fi
