#!/bin/bash

set -e
set -o pipefail

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

trap 'echo -e "${RED}Error occurred on line $LINENO. Exiting...${NC}"' ERR

# logging function
log() {
    echo -e "${YELLOW}[$(date +'%Y-%m-%d %H:%M:%S')]${NC} $1"
}

check_version() {
    local tool=$1
    local version_cmd=$2
    echo -e "\n$tool Version:" >> "$log_file"
    if command -v $tool &> /dev/null; then
        $version_cmd >> "$log_file" 2>&1
    else
        echo -e "${RED}$tool is not installed.${NC}" >> "$log_file"
    fi
}

usage() {
    log "Usage: $0 <CONFIG_DIRECTORY> <FASTQ_FILES> <PAIRED_END> <HEADCROP>"
    exit 1
}


if [ "$#" -ne 4 ]; then
    usage
fi

config_directory=$1
fastq_files=$2
paired_end=$3
headcrop=$4

mkdir -p "${config_directory}" || { echo -e "${RED}Failed to create ${config_directory}${NC}"; exit 1; }

fastqc_multiqc_output_dir="${config_directory}/1_fastqc_and_multiqc_reports"
trim_output_fpath="${config_directory}/2_trim_output"
map_output_fpath="${config_directory}/3_hisat2_mapping_output"
map_metrics_output_fpath="${config_directory}/3_1_map_metrics_output_qc"
filter_output_fpath="${config_directory}/4_filter_output"
qualimap_output_fpath="${config_directory}/4_1_qualimap_filter_output_qc"
counts_output_fpath="${config_directory}/5_stringtie_counts_output"
raw_counts_output_fpath="${config_directory}/6_raw_counts_output"

timestamp=$(date +"%Y%m%d_%H%M%S")
log_file="${config_directory}/7_pipeline_stats_${timestamp}.log"

{
  log "x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x"
  log "x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x"
  log "x - x - x - x - x - x - x - x - Psomagen Bulk RNA-seq Pipeline Snapshot - x - x - x - x - x - x - x - x - x - x -x"
  log "x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x"
  log "x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x"
  log " "
  log "Configuration Directory: ${config_directory}"
  log "Fastq files: ${fastq_files}"
  log "Paired End: ${paired_end}"
  log "Headcrop: ${headcrop}"
  log " "

  log "Pipeline Stages:"
  log "1. FastQC and MultiQC reports stored at: ${fastqc_multiqc_output_dir}"
  log "2. Trimmomatic trimming stored at: ${trim_output_fpath}"
  log "3. HISAT2 genome alignment stored at: ${map_output_fpath}"
  log "3_1. Mapping metrics stored at: ${map_metrics_output_fpath}"
  log "4. Post-alignment processing (Filtering, Deduplication) stored at: ${filter_output_fpath}"
  log "4_1. Qualimap QC stored at: ${qualimap_output_fpath}"
  log "5. StringTie counts stored at: ${counts_output_fpath}"
  log "6. Raw counts generated using HTSeq & Feature Counts stored at: ${raw_counts_output_fpath}"
  log " "

  log "Checking software versions..."

} > "$log_file"


check_version "nextflow" "nextflow -v"
check_version "conda" "conda --version"
check_version "trimmomatic" "trimmomatic -version"
check_version "fastqc" "fastqc --version"
check_version "multiqc" "multiqc --version"
check_version "hisat2" "hisat2 --version"
check_version "samtools" "samtools --version"
check_version "sambamba" "sambamba --version"
check_version "bedtools" "bedtools --version"
check_version "stringtie" "stringtie --version"
check_version "qualimap" "qualimap --version"
check_version "htseq-count" "pip show HTSeq"
check_version "R" "R --version"

{
  log "x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x"
  log "x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x"
  log "Pipeline stats saved to: $log_file"
} >> "$log_file"

log "Pipeline stats saved to: $log_file"

exit 0
