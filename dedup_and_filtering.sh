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
    log "Usage: $0 <OUTPUT_DIRECTORY> <NUM_THREADS> <BLACKLIST_BED_FILE> <FILTER_BED_FILE>"
    exit 1
}

if [ "$#" -ne 4 ]; then
    usage
fi

output_directory=$1
num_threads=$2
blist_bed_file=$3
exclude_bed_file=$4

if [ ! -d "$output_directory" ]; then
    echo -e "${RED}Error: Config directory not found: $output_directory${NC}"
    exit 1
fi

if [ ! -d "${output_directory}/3_hisat2_mapping_output" ]; then
    echo -e "${RED}Error: HISAT2 mapping output directory not found${NC}"
    exit 1
fi

if [ ! -f "$blist_bed_file" ] && [ ! -f "$exclude_bed_file" ]; then
    echo -e "${YELLOW}Warning: No BED files for filtering found. Proceeding without filtering.${NC}"
fi

# INPUT & OUTPUT PATHS
map_output_fpath="${output_directory}/3_hisat2_mapping_output"
filter_output_fpath="${output_directory}/4_filter_output"
filter_log_fpath="${filter_output_fpath}/log"

mkdir -p "$filter_output_fpath" "$filter_log_fpath"

log "Processing BAM files in ${map_output_fpath}"


for mapped_bam_file in "${map_output_fpath}"/*.sorted.bam; do
    if [[ -f "${mapped_bam_file}" ]]; then
        sample=$(basename "${mapped_bam_file}" .sorted.bam)
        log "Processing $sample"

        # Deduplication
        log "Deduplicating ${mapped_bam_file}"
        log "Command: sambamba markdup -r -t ${num_threads} ${mapped_bam_file} ${filter_output_fpath}/${sample}.dedup.bam"
        sambamba markdup -r -t ${num_threads} "${mapped_bam_file}" "${filter_output_fpath}/${sample}.dedup.bam" 2> "${filter_log_fpath}/${sample}_dedup_step.log"

        # Indexing the deduplicated sample
        log "Indexing ${filter_output_fpath}/${sample}.dedup.bam"
        log "Command: samtools index ${filter_output_fpath}/${sample}.dedup.bam"
        samtools index "${filter_output_fpath}/${sample}.dedup.bam" 2> "${filter_log_fpath}/${sample}_dedupindexing_step.log"

        # Filtering based on the presence of BED files
        log "Filtering ${sample}.dedup.bam"
        if [[ -f "${exclude_bed_file}" && -f "${blist_bed_file}" ]]; then
            log "Command: bedtools intersect -v -abam ${filter_output_fpath}/${sample}.dedup.bam -b ${exclude_bed_file} ${blist_bed_file} > ${filter_output_fpath}/${sample}.filt.bam"
            bedtools intersect -v -abam "${filter_output_fpath}/${sample}.dedup.bam" -b "${exclude_bed_file}" "${blist_bed_file}" > "${filter_output_fpath}/${sample}.filt.bam" 2> "${filter_log_fpath}/${sample}_filt_step.log"
        elif [[ -f "${exclude_bed_file}" ]]; then
            log "Command: bedtools intersect -v -abam ${filter_output_fpath}/${sample}.dedup.bam -b ${exclude_bed_file} > ${filter_output_fpath}/${sample}.filt.bam"
            bedtools intersect -v -abam "${filter_output_fpath}/${sample}.dedup.bam" -b "${exclude_bed_file}" > "${filter_output_fpath}/${sample}.filt.bam" 2> "${filter_log_fpath}/${sample}_filt_step.log"
        elif [[ -f "${blist_bed_file}" ]]; then
            log "Command: bedtools intersect -v -abam ${filter_output_fpath}/${sample}.dedup.bam -b ${blist_bed_file} > ${filter_output_fpath}/${sample}.filt.bam"
            bedtools intersect -v -abam "${filter_output_fpath}/${sample}.dedup.bam" -b "${blist_bed_file}" > "${filter_output_fpath}/${sample}.filt.bam" 2> "${filter_log_fpath}/${sample}_filt_step.log"
        else
            log "No BED files found for filtering. Skipping filtering step."
            cp "${filter_output_fpath}/${sample}.dedup.bam" "${filter_output_fpath}/${sample}.filt.bam"
        fi

        # Indexing the filtered sample
        log "Indexing ${filter_output_fpath}/${sample}.filt.bam"
        log "Command: samtools index ${filter_output_fpath}/${sample}.filt.bam"
        samtools index "${filter_output_fpath}/${sample}.filt.bam" 2> "${filter_log_fpath}/${sample}_filtindexing_step.log"
    else
        log "No file ${mapped_bam_file} found in ${map_output_fpath}"
    fi
done

log "Processing completed for all BAM files."

exit 0
