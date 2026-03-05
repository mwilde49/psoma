#!/bin/bash
set -euo pipefail

# ── Validate all tools are installed in the Psomagen container ──
# Usage: ./test_container.sh [path/to/psomagen_v*.sif]

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

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
    echo -e "${RED}ERROR: No .sif file found. Build first or pass path as argument.${NC}"
    exit 1
fi

echo -e "${YELLOW}Testing container: ${SIF}${NC}"
echo ""

PASS=0
FAIL=0

check_tool() {
    local name="$1"
    shift
    if apptainer exec "$SIF" "$@" &>/dev/null; then
        echo -e "  ${GREEN}PASS${NC}  $name"
        PASS=$((PASS + 1))
    else
        echo -e "  ${RED}FAIL${NC}  $name"
        FAIL=$((FAIL + 1))
    fi
}

# ── Phase 1: Check tools are on PATH ──
echo "Phase 1: Checking tools on PATH"
for tool in fastqc multiqc hisat2 trimmomatic samtools sambamba bedtools stringtie htseq-count qualimap R Rscript python nextflow java; do
    check_tool "$tool" which "$tool"
done
echo ""

# ── Phase 2: Check version commands ──
echo "Phase 2: Checking version output"
check_tool "fastqc --version"       fastqc --version
check_tool "multiqc --version"      multiqc --version
check_tool "hisat2 --version"       hisat2 --version
check_tool "trimmomatic -version"   trimmomatic -version
check_tool "samtools --version"     samtools --version
check_tool "sambamba --version"     sambamba --version
check_tool "bedtools --version"     bedtools --version
check_tool "stringtie --version"    stringtie --version
check_tool "htseq-count --version"  htseq-count --version
# Nextflow needs a writable home directory for its config
NF_TMP=$(mktemp -d)
echo -n "  "
if apptainer exec --env HOME="${NF_TMP}" "$SIF" nextflow -v &>/dev/null; then
    echo -e "${GREEN}PASS${NC}  nextflow -v"
    PASS=$((PASS + 1))
else
    echo -e "${RED}FAIL${NC}  nextflow -v"
    FAIL=$((FAIL + 1))
fi
rm -rf "${NF_TMP}"
check_tool "java -version"          java -version
echo ""

# ── Phase 3: Check R packages ──
echo "Phase 3: Checking R packages"
check_tool "Rsubread" Rscript -e 'library(Rsubread)'
echo ""

# ── Phase 4: Check Python packages ──
echo "Phase 4: Checking Python packages"
check_tool "pandas" python -c 'import pandas'
echo ""

# ── Phase 5: Smoke test ──
echo "Phase 5: Functional smoke test"
TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# Create a tiny FASTQ for FastQC smoke test
echo -e "@read1\nACGTACGT\n+\nIIIIIIII" > "${TMPDIR}/test.fastq"
check_tool "FastQC smoke test" fastqc "${TMPDIR}/test.fastq" --outdir="${TMPDIR}"

echo ""
echo -e "${GREEN}PASSED: ${PASS}${NC}  ${RED}FAILED: ${FAIL}${NC}"

if [[ $FAIL -gt 0 ]]; then
    exit 1
fi
