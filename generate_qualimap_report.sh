#!/bin/bash

set -e
set -o pipefail

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

trap 'echo -e "${RED}Error occurred on line $LINENO. Exiting...${NC}"' ERR

log() {
    echo -e "${YELLOW}[$(date +'%Y-%m-%d %H:%M:%S')]${NC} $1"
}

usage() {
    log "Usage: $0 <OUTPUT_DIRECTORY> <NUM_THREADS> <REFERENCE_GTF> <PAIRED_END (true/false)>"
    exit 1
}

if [ "$#" -ne 4 ]; then
    usage
fi

output_directory=$1
num_threads=$2
ref_genome_file=$3
paired_end=$4

memory_size="10G"

# INPUT & OUTPUT PATHS
filter_output_fpath="${output_directory}/4_filter_output"
qualimap_output_fpath="${output_directory}/4_1_qualimap_filter_output_qc"

mkdir -p "${qualimap_output_fpath}"
cd "${qualimap_output_fpath}"

# Check if reference GTF file exists
if [ ! -f "$ref_genome_file" ]; then
    echo -e "${RED}Error: Reference GTF file not found: $ref_genome_file${NC}"
    exit 1
fi

if [ ! -d "${filter_output_fpath}" ]; then
    echo -e "${RED}Error: Filtered files directory not found: $filter_output_fpath${NC}"
    exit 1
fi

for file in "$filter_output_fpath"/*.filt.bam; do
    log " "
    log "Running Qualimap on $file"

    if [ "$paired_end" = "true" ]; then
        log "Running Qualimap for paired-end data"
        log "qualimap rnaseq -pe -bam $file -outdir $qualimap_output_fpath -nt $num_threads -gtf $ref_genome_file --java-mem-size=$memory_size"
        qualimap rnaseq -pe -bam "$file" -outdir "$qualimap_output_fpath" -nt "$num_threads" -gtf "$ref_genome_file" --java-mem-size="$memory_size"
    else
        log "Running Qualimap for single-end data"
        log "qualimap rnaseq -bam $file -outdir $qualimap_output_fpath -nt $num_threads -gtf $ref_genome_file --java-mem-size=$memory_size"
        qualimap rnaseq -bam "$file" -outdir "$qualimap_output_fpath" -nt "$num_threads" -gtf "$ref_genome_file" --java-mem-size="$memory_size"
    fi
    log " "
done

log "Qualimap output directory: $qualimap_output_fpath"
log "Qualimap pipeline completed successfully."

exit 0
