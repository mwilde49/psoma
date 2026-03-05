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
    log "Usage: $0 <CONFIG_DIR> <SAMPLES_FILE> <FASTQ_DIR> <PAIRED_END> <READ1_SUFFIX> <READ2_SUFFIX> <NUM_THREADS> <ILLUMINA_CLIP_FILE> <HEADCROP> <LEADING> <TRAILING> <SLIDINGWINDOW> <MINLEN> <ILLUMINACLIP_PARAMS>"
    exit 1
}

if [ "$#" -ne 14 ]; then
    usage
fi

config_directory=$1
samples_file=$2
fastq_dir=$3
paired_end=$4
read1_suffix=$5
read2_suffix=$6
num_threads=$7
illumina_clip_file=$8
headcrop=$9
leading=${10}
trailing=${11}
slidingwindow=${12}
minlen=${13}
illuminaclip_params=${14}

# OUTPUT PATHS
trim_output_fpath="${config_directory}/2_trim_output"
trim_log_fpath="${trim_output_fpath}/log"

mkdir -p "${trim_output_fpath}" "${trim_log_fpath}"

log "Starting Trimmomatic trimming..."
log "Samples file: ${samples_file}"
log "FASTQ directory: ${fastq_dir}"
log "Paired-end: ${paired_end}"
log "Threads: ${num_threads}"
log "Illumina clip file: ${illumina_clip_file}"
log "HEADCROP: ${headcrop}, LEADING: ${leading}, TRAILING: ${trailing}"
log "SLIDINGWINDOW: ${slidingwindow}, MINLEN: ${minlen}"
log "ILLUMINACLIP params: ${illuminaclip_params}"

for sample_name in $(cat "${samples_file}")
do
    log "Trimming sample: ${sample_name}"

    if [ "${paired_end}" = "true" ]; then
        log "Command: trimmomatic PE -phred33 -threads ${num_threads} ${fastq_dir}/${sample_name}${read1_suffix}.fastq.gz ${fastq_dir}/${sample_name}${read2_suffix}.fastq.gz ${trim_output_fpath}/${sample_name}${read1_suffix}.trim.fastq.gz ${trim_output_fpath}/${sample_name}${read1_suffix}.trim_unpaired.fastq.gz ${trim_output_fpath}/${sample_name}${read2_suffix}.trim.fastq.gz ${trim_output_fpath}/${sample_name}${read2_suffix}.trim_unpaired.fastq.gz ILLUMINACLIP:${illumina_clip_file}:${illuminaclip_params} HEADCROP:${headcrop} LEADING:${leading} TRAILING:${trailing} SLIDINGWINDOW:${slidingwindow} MINLEN:${minlen}"
        ( trimmomatic PE -phred33 -threads ${num_threads} \
            ${fastq_dir}/${sample_name}${read1_suffix}.fastq.gz \
            ${fastq_dir}/${sample_name}${read2_suffix}.fastq.gz \
            ${trim_output_fpath}/${sample_name}${read1_suffix}.trim.fastq.gz \
            ${trim_output_fpath}/${sample_name}${read1_suffix}.trim_unpaired.fastq.gz \
            ${trim_output_fpath}/${sample_name}${read2_suffix}.trim.fastq.gz \
            ${trim_output_fpath}/${sample_name}${read2_suffix}.trim_unpaired.fastq.gz \
            ILLUMINACLIP:${illumina_clip_file}:${illuminaclip_params} \
            HEADCROP:${headcrop} \
            LEADING:${leading} TRAILING:${trailing} \
            SLIDINGWINDOW:${slidingwindow} \
            MINLEN:${minlen} ) 2> ${trim_log_fpath}/${sample_name}_trim_step.log
    else
        log "Command: trimmomatic SE -phred33 -threads ${num_threads} ${fastq_dir}/${sample_name}${read1_suffix}.fastq.gz ${trim_output_fpath}/${sample_name}${read1_suffix}.trim.fastq.gz ILLUMINACLIP:${illumina_clip_file}:${illuminaclip_params} HEADCROP:${headcrop} LEADING:${leading} TRAILING:${trailing} SLIDINGWINDOW:${slidingwindow} MINLEN:${minlen}"
        ( trimmomatic SE -phred33 -threads ${num_threads} \
            ${fastq_dir}/${sample_name}${read1_suffix}.fastq.gz \
            ${trim_output_fpath}/${sample_name}${read1_suffix}.trim.fastq.gz \
            ILLUMINACLIP:${illumina_clip_file}:${illuminaclip_params} \
            HEADCROP:${headcrop} \
            LEADING:${leading} TRAILING:${trailing} \
            SLIDINGWINDOW:${slidingwindow} \
            MINLEN:${minlen} ) 2> ${trim_log_fpath}/${sample_name}_trim_step.log
    fi

    if [ "$?" -ne 0 ]; then
        echo -e "${RED}Trimmomatic failed for sample: ${sample_name}${NC}"
        exit 1
    fi
    log "Trimming complete for: ${sample_name}"
done

log "Trimmomatic trimming completed for all samples."

exit 0
