#!/bin/bash

set -e

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'  # No color

log() {
    echo -e "${YELLOW}[$(date +'%Y-%m-%d %H:%M:%S')]${NC} $1"
}

usage() {
    echo -e "${RED}Usage:${NC} $0 <OUTPUT_DIRECTORY> <CONFIG_DIRECTORY> <NUM_THREADS> <REFERENCE_GTF> <STRAND_ST> <SPECIES>"
    exit 1
}

if [ "$#" -ne 6 ]; then
    usage
fi

output_directory=$1
config_directory=$2
num_threads=$3
reference_gtf=$4
strand_st=$5
species=$6

cd $output_directory

# INPUT & OUTPUT PATHS
filter_output_fpath="${output_directory}/4_filter_output"
counts_output_fpath="${output_directory}/5_stringtie_counts_output"
counts_log_fpath="${counts_output_fpath}/log"

mkdir -p ${counts_output_fpath} ${counts_log_fpath}

log "Processing stringtie counts for files in ${filter_output_fpath}"
log " "
for filtered_file in "${filter_output_fpath}"/*.filt.bam; do
    if [[ -f "$filtered_file" ]]; then
        temp_file="${filtered_file##*/}"
        libr_name="${temp_file%.filt.bam}"
        log "${libr_name}"
        log "Running stringtie: ${filtered_file}"
        log "( stringtie -p ${num_threads} -e -v -G ${reference_gtf} ${strand_st} ${filtered_file} -o ${counts_output_fpath}/${libr_name}.assemb.gtf -A ${counts_output_fpath}/${libr_name}.tab ) 2> ${counts_log_fpath}/${libr_name}_stringtie_step.log"
        ( stringtie -p ${num_threads} -e -v -G ${reference_gtf} ${strand_st} ${filtered_file} -o ${counts_output_fpath}/${libr_name}.assemb.gtf -A ${counts_output_fpath}/${libr_name}.tab ) 2> ${counts_log_fpath}/${libr_name}_stringtie_step.log
    else
        log "No ${filtered_file}.filt.bam file found in $filter_output_fpath"
    fi

done


# merging the stringtie outputs
log "Merging stringtie outputs..."
mkdir -p "${output_directory}/stringtie"
cp -pr "${counts_output_fpath}"/*.tab "${output_directory}/stringtie"

ls ./stringtie/ > sample_list
FILE_NO=`cat sample_list | wc -l`
START_TIME=$(date +%s)

while read x
do
    cat stringtie/$x | awk '{ print $1"#"$2"#"$3"#"$4"#"$5"#"$6,$9;}' | sort -k1,1 > ./stringtie/$x.txt
done < sample_list

counter=0

while read x
do
    counter=$(($counter + 1))
    if [ $counter -eq 1 ]
    then
        cat ./stringtie/$x.txt > ./stringtie/tmp
    else
        join ./stringtie/tmp ./stringtie/$x.txt > ./stringtie/tmp2
        mv ./stringtie/tmp2 ./stringtie/tmp
    fi
done < sample_list

echo ensmbl_gene_id gene_name chr strand start stop `cat sample_list | tr '\n' " "` > ./stringtie/genes.tpm.txt
cat ./stringtie/tmp | tr "#" " " > ./stringtie/tmp2
cat ./stringtie/tmp2 | sort -k1,1 > ./stringtie/tmp
cat ./stringtie/tmp >> ./stringtie/genes.tpm.txt

log "Done merging output!"

FINISH_TIME=$(($(date +%s) - $START_TIME))
log "$FILE_NO stringtie jobs finished in $FINISH_TIME seconds!"

cp -pr "${output_directory}/stringtie/genes.tpm.txt" "${counts_output_fpath}"
rm -rf  "${output_directory}/stringtie" "${output_directory}/sample_list" "${counts_output_fpath}/tmp_*"

cd $counts_output_fpath

if [[ "$species" =~ ^[Hh][Uu][Mm][Aa][Nn]$ ]]; then
    ref_file="${config_directory}/ref_human_geneid_genename_genebiotype.tsv"
    awk 'NR==FNR {gene_biotype[$1]=$3; next} FNR==1 {print $0, "gene_biotype"; next} ($1 in gene_biotype) {print $0, gene_biotype[$1]}' $ref_file "${counts_output_fpath}/genes.tpm.txt" > "${counts_output_fpath}/genes_tpm.txt"
elif [[ "$species" =~ ^[Mm][Oo][Uu][Ss][Ee]$ ]]; then
    ref_file="${config_directory}/ref_mouse_geneid_genename_genebiotype.tsv"
    awk 'NR==FNR {gene_biotype[$1]=$3; next} FNR==1 {print $0, "gene_biotype"; next} ($1 in gene_biotype) {print $0, gene_biotype[$1]}' $ref_file "${counts_output_fpath}/genes.tpm.txt" > "${counts_output_fpath}/genes_tpm.txt"
elif [[ "$species" =~ ^[Rr][Aa][Tt][Tt][Uu][Ss]$ ]]; then
    ref_file="${config_directory}/ref_rattus_geneid_genename_genebiotype.tsv"
    awk 'NR==FNR {gene_biotype[$1]=$3; next} FNR==1 {print $0, "gene_biotype"; next} ($1 in gene_biotype) {print $0, gene_biotype[$1]}' $ref_file "${counts_output_fpath}/genes.tpm.txt" > "${counts_output_fpath}/genes_tpm.txt"
fi

log "Processing completed."

exit 0
