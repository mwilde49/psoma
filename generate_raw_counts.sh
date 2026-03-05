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
    log "Usage: $0 <CONFIG_DIRECTORY> <REFERENCE_GTF> <PAIRED_END> <STRAND_HTS> <PAIRED_HTS> <SPECIES>"
    exit 1
}

if [ "$#" -ne 6 ]; then
    usage
fi

config_directory=$1
reference_gtf=$2
paired_end=$3
strand_hts=$4
paired_hts=$5
species=$6

if [ ! -f "$reference_gtf" ]; then
    echo -e "${RED}Error: Reference GTF file not found: $reference_gtf${NC}"
    exit 1
fi

cd "$config_directory" || { echo -e "${RED}Error: Config directory does not exist${NC}"; exit 1; }

# INPUT & OUTPUT PATHS
filter_output_fpath="${config_directory}/4_filter_output"
raw_counts_output_fpath="${config_directory}/6_raw_counts_output"
raw_counts_log_fpath="${raw_counts_output_fpath}/log"
merged_param_string=""
header_param_string="Gene_ID,"

mkdir -p "$raw_counts_output_fpath" "$raw_counts_log_fpath"

log "Starting HTSeq processing..."

for filtered_file in "${filter_output_fpath}"/*.filt.bam; do
    if [[ -f "$filtered_file" ]]; then
        temp_file="${filtered_file##*/}"
        libr_name="${temp_file%.filt.bam}"

        # Sort and index BAM file
        log "Sorting and indexing BAM file: $libr_name"
        samtools sort -o "${raw_counts_output_fpath}/${libr_name}.sorted.filt.bam" -O bam "$filtered_file"
        samtools index "${raw_counts_output_fpath}/${libr_name}.sorted.filt.bam"

        # Add BAM file to merged string
        merged_param_string+="${raw_counts_output_fpath}/${libr_name}.sorted.filt.bam "
        header_param_string+="${libr_name},"
    else
        log "No BAM file found at $filtered_file, skipping..."
    fi
done


header_param_string="${header_param_string%,}"  # removes the trailing comma

log "Generating raw counts using HTSeq..."

if [ "$paired_end" = "true" ]; then
    log "Running HTSeq: ( ( htseq-count -f bam -s ${strand_hts} -r ${paired_hts} ${merged_param_string} ${reference_gtf} ) > ${raw_counts_output_fpath}/raw_htseq_counts.tsv ) 2> ${raw_counts_log_fpath}/count_HTS_step.log"
    ( ( htseq-count -f bam -s ${strand_hts} -r ${paired_hts} ${merged_param_string} ${reference_gtf} ) > ${raw_counts_output_fpath}/raw_htseq_counts.tsv ) 2> ${raw_counts_log_fpath}/count_HTS_step.log
else
    log "Running HTSeq: ( ( htseq-count -f bam -s ${strand_hts} ${merged_param_string} ${reference_gtf} ) > ${raw_counts_output_fpath}/raw_htseq_counts.tsv ) 2> ${raw_counts_log_fpath}/count_HTS_step.log"
    ( ( htseq-count -f bam -s ${strand_hts} ${merged_param_string} ${reference_gtf} ) > ${raw_counts_output_fpath}/raw_htseq_counts.tsv ) 2> ${raw_counts_log_fpath}/count_HTS_step.log
fi

log "Converting HTSeq output to CSV..."
cat "${raw_counts_output_fpath}/raw_htseq_counts.tsv" | tr "\t" "," > "${raw_counts_output_fpath}/raw_htseq_counts.csv"
echo "$header_param_string" | cat - "${raw_counts_output_fpath}/raw_htseq_counts.csv" > "${raw_counts_output_fpath}/temp.csv"
mv "${raw_counts_output_fpath}/temp.csv" "${raw_counts_output_fpath}/raw_htseq_counts.csv"

cd $raw_counts_output_fpath

if [[ "$species" =~ ^[Hh][Uu][Mm][Aa][Nn]$ ]]; then
    ref_gene_id_name_human_file="${config_directory}/ref_human_geneid_genename_genebiotype.tsv"
    python "${config_directory}/fetch_genename_genebiotype_for_counts.py" "${config_directory}" "${ref_gene_id_name_human_file}"
elif [[ "$species" =~ ^[Mm][Oo][Uu][Ss][Ee]$ ]]; then
    ref_gene_id_name_mouse_file="${config_directory}/ref_mouse_geneid_genename_genebiotype.tsv"
    python "${config_directory}/fetch_genename_genebiotype_for_counts.py" "${config_directory}" "${ref_gene_id_name_mouse_file}"
elif [[ "$species" =~ ^[Rr][Aa][Tt][Tt][Uu][Ss]$ ]]; then
    ref_gene_id_name_rattus_file="${config_directory}/ref_rattus_geneid_genename_genebiotype.tsv"
    python "${config_directory}/fetch_genename_genebiotype_for_counts.py" "${config_directory}" "${ref_gene_id_name_rattus_file}"
fi

log "HTSeq counts generated successfully."

log "Cleaning up sorted BAM files..."
rm -rf "${raw_counts_output_fpath}"/*.sorted.filt.bam "${raw_counts_output_fpath}"/*.sorted.filt.bam.bai

log "HTSeq pipeline completed successfully."

exit 0
