#!/bin/bash

set -e

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'  # No color

log() {
    echo -e "${YELLOW}[$(date +'%Y-%m-%d %H:%M:%S')]${NC} $1"
}


if [ "$#" -ne 3 ]; then
    echo -e "${RED}Error: Invalid number of arguments.${NC}"
    log "Usage: $0 <fastq_input_path> <num_threads> <output_dir>"
    exit 1
fi

# Input variables
fastq_input_path=$1
num_threads=$2
output_dir=$3/1_fastqc_and_multiqc_reports

# Check if the fastq input path exists
if [ ! -d "$fastq_input_path" ]; then
    echo -e "${RED}Error: Input directory $fastq_input_path does not exist.${NC}"
    exit 1
fi

if [ -z "$num_threads" ]; then
    num_threads=$(nproc)
    log "Number of threads not provided. Using default: $num_threads"
else
    log "Using $num_threads threads."
fi

if [ ! -d "$output_dir" ]; then
    log "Creating output directory: $output_dir"
    mkdir -p "$output_dir"
fi

cd $fastq_input_path

start_time=$(date +%s)

# Run FastQC on all fastq.gz files
log "Running FastQC on files in $fastq_input_path..."
for file in *.fastq.gz; do
    if [[ -f "$file" ]]; then
        log "Processing $file..."
        fastqc -t $num_threads "$file" --outdir="$output_dir"
    else
        echo -e "${RED}No fastq.gz files found in the directory.${NC}"
        exit 1
    fi
done

log "Successfully completed FastQC reports."

cd $output_dir

# Run MultiQC
log "Running MultiQC in $output_dir..."
multiqc .

log "Successfully completed MultiQC reports."

end_time=$(date +%s)
execution_time=$((end_time - start_time))

log "FastQC and MultiQC reports generated successfully."
log "Reports are stored at: $output_dir"
log "Total execution time: $(($execution_time / 60)) minutes and $(($execution_time % 60)) seconds."
