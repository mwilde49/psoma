#!/usr/bin/env nextflow

// Psomagen Bulk RNA-seq pipeline is a sequential pipeline that generates counts from Psomagen-processed fastq data
// Analogous to: utdal/Bulk-RNA-Seq-Nextflow-Pipeline
// Aligner: HISAT2 (replaces STAR)
// Trimming: Trimmomatic (adapter + quality trimming for Nextera adapters)
//
// Here are the steps involved in the Psomagen Bulk RNA-seq Analysis:
// 1. Quality Control - Generate FastQC and MultiQC reports
// 2. Trimming - Adapter and quality trimming using Trimmomatic
// 3. Genome Alignment - Map reads to the reference genome using HISAT2
// 3_1. Mapping Metrics - Generate mapping statistics and quality reports
// 4. Post-Alignment Processing - Filter, deduplicate, and index aligned reads
// 4_1. Quality Assessment - Generate Qualimap reports for aligned reads
// 5. Transcript Assembly and Quantification - Generate counts using StringTie
// 6. Raw Count Generation - Generate raw counts using HTSeq && Feature Counts using Rsubread's featureCounts


println "Psomagen Bulk RNA-seq Pipeline"
println "  "
println "Here are the steps involved in the Psomagen Bulk RNA-seq Analysis:"
println "1. Quality Control - Generate FastQC and MultiQC reports"
println "2. Trimming - Adapter and quality trimming using Trimmomatic"
println "3. Genome Alignment - Map reads to the reference genome using HISAT2"
println "3_1. Mapping Metrics - Generate mapping statistics and quality reports"
println "4. Post-Alignment Processing - Filter, deduplicate, and index aligned reads"
println "4_1. Quality Assessment - Generate Qualimap reports for aligned reads"
println "5. Transcript Assembly and Quantification - Generate counts using StringTie"
println "6. Raw Count Generation - Generate raw counts using HTSeq && Feature Counts using Rsubread's featureCounts"


process generate_fastqc_multiqc_reports {
    // Process to generate FastQC and MultiQC reports
    debug true
    errorStrategy 'terminate'

    input:
    val config_directory
    val fastq_files
    val fastqc_cores

    output:
    val "${config_directory}/1_fastqc_and_multiqc_reports", emit: fastqc_output

    script:
    """
    mkdir -p "${config_directory}/0_nextflow_logs"
    echo "[INFO] Starting FastQC and MultiQC at \$(date)" | tee -a ${config_directory}/0_nextflow_logs/fastqc_multiqc.log

    if [ ! -d "${fastq_files}" ]; then
        echo "[ERROR] FASTQ files directory not found: ${fastq_files}" | tee -a ${config_directory}/0_nextflow_logs/fastqc_multiqc.log
        exit 1
    fi

    echo 'bash ${config_directory}/generate_fastqc_reports.sh ${fastq_files} ${fastqc_cores} ${config_directory}'
    bash ${config_directory}/generate_fastqc_reports.sh ${fastq_files} ${fastqc_cores} ${config_directory}

    if [ "\$?" -ne 0 ]; then
        echo "[ERROR] FastQC failed on files in ${fastq_files}" | tee -a ${config_directory}/0_nextflow_logs/fastqc_multiqc.log
        exit 1
    fi
    echo "[INFO] FastQC completed successfully for all FASTQ files" | tee -a ${config_directory}/0_nextflow_logs/fastqc_multiqc.log
    """
}


process trimmomatic_trimming {
    // Process to perform Trimmomatic adapter and quality trimming
    debug true
    errorStrategy 'terminate'

    input:
    val config_directory
    val samples_file
    val fastq_files
    val is_paired_end
    val read1_suffix
    val read2_suffix
    val num_threads
    val illumina_clip_file
    val headcrop
    val leading
    val trailing
    val slidingwindow
    val minlen
    val illuminaclip_params

    output:
    val "${config_directory}/2_trim_output", emit: trimmed_files

    script:
    """
    mkdir -p "${config_directory}/0_nextflow_logs"
    echo "[INFO] Starting Trimmomatic trimming at \$(date)" | tee -a ${config_directory}/0_nextflow_logs/trimmomatic.log

    if [ ! -f "${samples_file}" ]; then
        echo "[ERROR] Samples file not found: ${samples_file}" | tee -a ${config_directory}/0_nextflow_logs/trimmomatic.log
        exit 1
    fi

    if [ ! -d "${fastq_files}" ]; then
        echo "[ERROR] FASTQ files directory not found: ${fastq_files}" | tee -a ${config_directory}/0_nextflow_logs/trimmomatic.log
        exit 1
    fi

    echo "[INFO] Running Trimmomatic with the following parameters:" | tee -a ${config_directory}/0_nextflow_logs/trimmomatic.log
    echo "  Samples file: ${samples_file}" | tee -a ${config_directory}/0_nextflow_logs/trimmomatic.log
    echo "  FASTQ directory: ${fastq_files}" | tee -a ${config_directory}/0_nextflow_logs/trimmomatic.log
    echo "  Paired-end: ${is_paired_end}" | tee -a ${config_directory}/0_nextflow_logs/trimmomatic.log
    echo "  Threads: ${num_threads}" | tee -a ${config_directory}/0_nextflow_logs/trimmomatic.log
    echo "  HEADCROP: ${headcrop}" | tee -a ${config_directory}/0_nextflow_logs/trimmomatic.log

    echo 'bash ${config_directory}/trimmomatic_trimming.sh ${config_directory} ${samples_file} ${fastq_files} ${is_paired_end} ${read1_suffix} ${read2_suffix} ${num_threads} ${illumina_clip_file} ${headcrop} ${leading} ${trailing} ${slidingwindow} ${minlen} ${illuminaclip_params}'
    bash ${config_directory}/trimmomatic_trimming.sh ${config_directory} ${samples_file} ${fastq_files} ${is_paired_end} ${read1_suffix} ${read2_suffix} ${num_threads} ${illumina_clip_file} ${headcrop} ${leading} ${trailing} ${slidingwindow} ${minlen} ${illuminaclip_params}

    if [ "\$?" -ne 0 ]; then
        echo "[ERROR] Trimmomatic trimming failed" | tee -a ${config_directory}/0_nextflow_logs/trimmomatic.log
        exit 1
    fi
    echo "[INFO] Trimmomatic trimming completed successfully at \$(date)" | tee -a ${config_directory}/0_nextflow_logs/trimmomatic.log
    """
}


process hisat2_mapping {
    // Process to perform HISAT2 alignment for RNA-seq data
    debug true
    errorStrategy 'terminate'

    input:
    val config_directory
    val samples_file
    val is_paired_end
    val read1_suffix
    val read2_suffix
    val num_threads
    val hisat2_index
    val trimmed_files

    output:
    val "${config_directory}/3_hisat2_mapping_output", emit: mapped_files

    script:
    def trim_output_fpath = "${config_directory}/2_trim_output"
    def map_output_fpath = "${config_directory}/3_hisat2_mapping_output"
    def map_log_fpath = "${map_output_fpath}/log"

    """
    mkdir -p "${config_directory}/0_nextflow_logs"
    echo "[INFO] Starting HISAT2 mapping at \$(date)" | tee -a ${config_directory}/0_nextflow_logs/hisat2_mapping.log

    if [ ! -f "${samples_file}" ]; then
        echo "[ERROR] Samples file not found: ${samples_file}" | tee -a ${config_directory}/0_nextflow_logs/hisat2_mapping.log
        exit 1
    fi

    if [ ! -d "${trim_output_fpath}" ]; then
        echo "[ERROR] Trimmed files directory not found: ${trim_output_fpath}" | tee -a ${config_directory}/0_nextflow_logs/hisat2_mapping.log
        exit 1
    fi

    mkdir -p ${map_output_fpath} ${map_log_fpath}

    for sample_name in \$(cat ${samples_file})
    do
        echo "Running HISAT2 for sample: \$sample_name"

        if [ "${is_paired_end}" = true ]; then
            hisat2 --summary-file ${map_log_fpath}/\${sample_name}_h2report.txt \
                -p ${num_threads} \
                -x ${hisat2_index} \
                -1 ${trim_output_fpath}/\${sample_name}${read1_suffix}.trim.fastq.gz \
                -2 ${trim_output_fpath}/\${sample_name}${read2_suffix}.trim.fastq.gz \
                -S ${map_output_fpath}/\${sample_name}.sam \
                2> ${map_log_fpath}/\${sample_name}_mapping_step.log
        else
            hisat2 --summary-file ${map_log_fpath}/\${sample_name}_h2report.txt \
                -p ${num_threads} \
                -x ${hisat2_index} \
                -U ${trim_output_fpath}/\${sample_name}${read1_suffix}.trim.fastq.gz \
                -S ${map_output_fpath}/\${sample_name}.sam \
                2> ${map_log_fpath}/\${sample_name}_mapping_step.log
        fi

        if [ "\$?" -ne 0 ]; then
            echo "[ERROR] HISAT2 mapping failed for sample: \${sample_name}" | tee -a ${config_directory}/0_nextflow_logs/hisat2_mapping.log
            exit 1
        fi

        echo "Sorting BAM for: \${sample_name}"
        samtools sort -@ ${num_threads} -o ${map_output_fpath}/\${sample_name}.sorted.bam ${map_output_fpath}/\${sample_name}.sam 2> ${map_log_fpath}/\${sample_name}_sorting_step.log
        rm ${map_output_fpath}/\${sample_name}.sam

        echo "Indexing BAM for: \${sample_name}"
        samtools index ${map_output_fpath}/\${sample_name}.sorted.bam 2> ${map_log_fpath}/\${sample_name}_indexing_step.log
    done

    echo "[INFO] HISAT2 mapping process completed at \$(date)" | tee -a ${config_directory}/0_nextflow_logs/hisat2_mapping.log
    """
}


process generate_mapping_metrics {
    // Process to generate mapping metrics and MultiQC report on the mapped files
    debug true
    errorStrategy 'terminate'

    input:
    val config_directory
    val mapped_files

    output:
    val "${config_directory}/3_1_map_metrics_output_qc", emit: map_metric_files

    script:
    """
    echo "[INFO] Starting mapping metrics generation at \$(date)" | tee -a ${config_directory}/0_nextflow_logs/mapping_metrics.log

    if [ ! -f "${config_directory}/generate_map_metrics.py" ]; then
        echo "[ERROR] Mapping metrics script not found: ${config_directory}/generate_map_metrics.py" | tee -a ${config_directory}/0_nextflow_logs/mapping_metrics.log
        exit 1
    fi

    if [ ! -d "${config_directory}/3_hisat2_mapping_output" ]; then
        echo "[ERROR] HISAT2 mapping output directory not found: ${config_directory}/3_hisat2_mapping_output" | tee -a ${config_directory}/0_nextflow_logs/mapping_metrics.log
        exit 1
    fi

    echo "[INFO] Running mapping metrics with the following parameters:" | tee -a ${config_directory}/0_nextflow_logs/mapping_metrics.log
    echo "  Config directory: ${config_directory}" | tee -a ${config_directory}/0_nextflow_logs/mapping_metrics.log
    echo "  Mapped files directory: ${config_directory}/3_hisat2_mapping_output" | tee -a ${config_directory}/0_nextflow_logs/mapping_metrics.log

    echo 'python "${config_directory}/generate_map_metrics.py" "${config_directory}" "${config_directory}/3_hisat2_mapping_output/log"' | tee -a ${config_directory}/0_nextflow_logs/mapping_metrics.log
    python "${config_directory}/generate_map_metrics.py" "${config_directory}" "${config_directory}/3_hisat2_mapping_output/log" | tee -a ${config_directory}/0_nextflow_logs/mapping_metrics.log

    cd "${config_directory}/3_hisat2_mapping_output"
    multiqc . | tee -a ${config_directory}/0_nextflow_logs/mapping_metrics.log
    mkdir -p "${config_directory}/3_1_map_metrics_output_qc"
    cp -pr multiqc_report.html multiqc_data "${config_directory}/3_1_map_metrics_output_qc" | tee -a ${config_directory}/0_nextflow_logs/mapping_metrics.log
    rm -r multiqc_report.html multiqc_data | tee -a ${config_directory}/0_nextflow_logs/mapping_metrics.log
    """
}


process filter_samples {
    // Process to filter and deduplicate aligned samples
    debug true
    errorStrategy 'terminate'

    input:
    val config_directory
    val num_threads
    val blist_bed_file
    val exclude_bed_file
    val mapped_files

    output:
    val "${config_directory}/4_filter_output", emit: filtered_files

    script:
    """
    echo "[INFO] Starting filtering and deduplication at \$(date)" | tee -a ${config_directory}/0_nextflow_logs/filtering.log

    if [ ! -f "${config_directory}/dedup_and_filtering.sh" ]; then
        echo "[ERROR] Deduplication and filtering script not found: ${config_directory}/dedup_and_filtering.sh" | tee -a ${config_directory}/0_nextflow_logs/filtering.log
        exit 1
    fi

    if [ ! -f "${blist_bed_file}" ] && [ ! -f "${exclude_bed_file}" ]; then
        echo "[ERROR] Both blacklist and exclude BED files are missing." | tee -a ${config_directory}/0_nextflow_logs/filtering.log
        exit 1
    fi

    echo "[INFO] Running deduplication and filtering with the following parameters:" | tee -a ${config_directory}/0_nextflow_logs/filtering.log
    echo "  Config directory: ${config_directory}" | tee -a ${config_directory}/0_nextflow_logs/filtering.log
    echo "  Number of threads: ${num_threads}" | tee -a ${config_directory}/0_nextflow_logs/filtering.log
    echo "  Blacklist BED file: ${blist_bed_file}" | tee -a ${config_directory}/0_nextflow_logs/filtering.log
    echo "  Exclude BED file: ${exclude_bed_file}" | tee -a ${config_directory}/0_nextflow_logs/filtering.log

    echo 'bash ${config_directory}/dedup_and_filtering.sh ${config_directory} ${num_threads} ${blist_bed_file} ${exclude_bed_file}' | tee -a ${config_directory}/0_nextflow_logs/filtering.log
    bash ${config_directory}/dedup_and_filtering.sh ${config_directory} ${num_threads} ${blist_bed_file} ${exclude_bed_file} | tee -a ${config_directory}/0_nextflow_logs/filtering.log
    if [ "\$?" -ne 0 ]; then
        echo "[ERROR] Filtering and deduplication failed" | tee -a ${config_directory}/0_nextflow_logs/filtering.log
        exit 1
    else
        echo "[INFO] Filtering and deduplication completed successfully at \$(date)" | tee -a ${config_directory}/0_nextflow_logs/filtering.log
    fi

    echo "[INFO] Filtering process completed at \$(date)" | tee -a ${config_directory}/0_nextflow_logs/filtering.log
    """
}


process generate_stringtie_counts {
    // Process to generate counts using StringTie
    debug true
    errorStrategy 'terminate'

    input:
    val config_directory
    val num_threads
    val reference_gtf
    val strand_st
    val species
    val filtered_files

    output:
    val "${config_directory}/5_stringtie_counts_output", emit: counts_files

    script:
    """
    echo "[INFO] Starting StringTie counts generation at \$(date)" | tee -a ${config_directory}/0_nextflow_logs/stringtie_counts.log

    if [ ! -f "${config_directory}/generate_stringtie_counts.sh" ]; then
        echo "[ERROR] StringTie counts script not found: ${config_directory}/generate_stringtie_counts.sh" | tee -a ${config_directory}/0_nextflow_logs/stringtie_counts.log
        exit 1
    fi

    if [ ! -f "${reference_gtf}" ]; then
        echo "[ERROR] Reference GTF file not found: ${reference_gtf}" | tee -a ${config_directory}/0_nextflow_logs/stringtie_counts.log
        exit 1
    fi

    echo "[INFO] Running StringTie with the following parameters:" | tee -a ${config_directory}/0_nextflow_logs/stringtie_counts.log
    echo "  Config directory: ${config_directory}" | tee -a ${config_directory}/0_nextflow_logs/stringtie_counts.log
    echo "  Number of threads: ${num_threads}" | tee -a ${config_directory}/0_nextflow_logs/stringtie_counts.log
    echo "  Reference GTF: ${reference_gtf}" | tee -a ${config_directory}/0_nextflow_logs/stringtie_counts.log
    echo "  Strand-specific: ${strand_st}" | tee -a ${config_directory}/0_nextflow_logs/stringtie_counts.log
    echo "  Species: ${species}" | tee -a ${config_directory}/0_nextflow_logs/stringtie_counts.log

    echo 'bash ${config_directory}/generate_stringtie_counts.sh ${config_directory} ${num_threads} ${reference_gtf} ${strand_st} ${species}' | tee -a ${config_directory}/0_nextflow_logs/stringtie_counts.log
    bash ${config_directory}/generate_stringtie_counts.sh ${config_directory} ${num_threads} ${reference_gtf} ${strand_st} ${species} | tee -a ${config_directory}/0_nextflow_logs/stringtie_counts.log

    if [ "\$?" -ne 0 ]; then
        echo "[ERROR] StringTie counts generation failed" | tee -a ${config_directory}/0_nextflow_logs/stringtie_counts.log
        exit 1
    else
        echo "[INFO] StringTie counts generated successfully at \$(date)" | tee -a ${config_directory}/0_nextflow_logs/stringtie_counts.log
    fi

    echo "[INFO] StringTie counts process completed at \$(date)" | tee -a ${config_directory}/0_nextflow_logs/stringtie_counts.log
    """
}


process generate_qualimap_reports {
    // Process to generate Qualimap reports for quality control of BAM files
    debug true
    errorStrategy 'terminate'

    input:
    val config_directory
    val reference_gtf
    val is_paired_end
    val num_threads
    val filtered_files

    output:
    val "${config_directory}/4_1_qualimap_filter_output_qc", emit: qualimap_files

    script:
    """
    echo "[INFO] Starting Qualimap report generation at \$(date)" | tee -a ${config_directory}/0_nextflow_logs/qualimap_report.log

    if [ ! -f "${config_directory}/generate_qualimap_report.sh" ]; then
        echo "[ERROR] Qualimap script not found: ${config_directory}/generate_qualimap_report.sh" | tee -a ${config_directory}/0_nextflow_logs/qualimap_report.log
        exit 1
    fi

    if [ ! -f "${reference_gtf}" ]; then
        echo "[ERROR] GTF file not found: ${reference_gtf}" | tee -a ${config_directory}/0_nextflow_logs/qualimap_report.log
        exit 1
    fi

    echo "[INFO] Running Qualimap with the following parameters:" | tee -a ${config_directory}/0_nextflow_logs/qualimap_report.log
    echo "  Config directory: ${config_directory}" | tee -a ${config_directory}/0_nextflow_logs/qualimap_report.log
    echo "  GTF file: ${reference_gtf}" | tee -a ${config_directory}/0_nextflow_logs/qualimap_report.log
    echo "  Paired-end: ${is_paired_end}" | tee -a ${config_directory}/0_nextflow_logs/qualimap_report.log
    echo "  Filtered BAM files: ${filtered_files}" | tee -a ${config_directory}/0_nextflow_logs/qualimap_report.log

    echo "bash ${config_directory}/generate_qualimap_report.sh ${config_directory} ${num_threads} ${reference_gtf} ${is_paired_end}" | tee -a ${config_directory}/0_nextflow_logs/qualimap_report.log
    bash "${config_directory}/generate_qualimap_report.sh" "${config_directory}" "${num_threads}" "${reference_gtf}" "${is_paired_end}" | tee -a ${config_directory}/0_nextflow_logs/qualimap_report.log

    if [ "\$?" -ne 0 ]; then
        echo "[ERROR] Qualimap report generation failed" | tee -a ${config_directory}/0_nextflow_logs/qualimap_report.log
        exit 1
    else
        echo "[INFO] Qualimap report generation completed successfully at \$(date)" | tee -a ${config_directory}/0_nextflow_logs/qualimap_report.log
    fi
    """
}


process generate_feature_counts {
    // Process to generate feature counts using Rsubread's featureCounts
    debug true
    errorStrategy 'terminate'

    input:
    val config_directory
    val num_threads
    val reference_gtf
    val is_paired_end
    val filtered_files

    output:
    val "${config_directory}/6_raw_counts_output", emit: feature_counts_files

    script:
    """
    echo "[INFO] Starting feature counts generation at \$(date)" | tee -a ${config_directory}/0_nextflow_logs/feature_counts.log

    if [ ! -f "${config_directory}/rsubread_featurecount_script.R" ]; then
        echo "[ERROR] Rsubread feature counts script not found: ${config_directory}/rsubread_featurecount_script.R" | tee -a ${config_directory}/0_nextflow_logs/feature_counts.log
        exit 1
    fi

    if [ ! -f "${reference_gtf}" ]; then
        echo "[ERROR] Reference GTF file not found: ${reference_gtf}" | tee -a ${config_directory}/0_nextflow_logs/feature_counts.log
        exit 1
    fi

    cd "${config_directory}"
    mkdir -p "${config_directory}/temp_feature_counts_dir"
    cp -r "${config_directory}/4_filter_output"/*.filt.bam "${config_directory}/temp_feature_counts_dir"

    echo "[INFO] Running featureCounts with the following parameters:" | tee -a ${config_directory}/0_nextflow_logs/feature_counts.log
    echo "  Config directory: ${config_directory}" | tee -a ${config_directory}/0_nextflow_logs/feature_counts.log
    echo "  Reference GTF: ${reference_gtf}" | tee -a ${config_directory}/0_nextflow_logs/feature_counts.log
    echo "  Paired-end: ${is_paired_end}" | tee -a ${config_directory}/0_nextflow_logs/feature_counts.log
    echo "  Filtered files directory: ${config_directory}/temp_feature_counts_dir" | tee -a ${config_directory}/0_nextflow_logs/feature_counts.log

    mkdir -p "${config_directory}/6_raw_counts_output"
    echo 'Rscript ${config_directory}/rsubread_featurecount_script.R "${config_directory}/temp_feature_counts_dir" ${reference_gtf} ${is_paired_end} ${num_threads} "${config_directory}/6_raw_counts_output/raw_feature_counts.csv"' | tee -a ${config_directory}/0_nextflow_logs/feature_counts.log
    Rscript ${config_directory}/rsubread_featurecount_script.R "${config_directory}/temp_feature_counts_dir" ${reference_gtf} ${is_paired_end} ${num_threads} "${config_directory}/6_raw_counts_output/raw_feature_counts.csv" | tee -a ${config_directory}/0_nextflow_logs/feature_counts.log
    rm -r "${config_directory}/temp_feature_counts_dir"

    if [ "\$?" -ne 0 ]; then
        echo "[ERROR] Feature counts generation failed" | tee -a ${config_directory}/0_nextflow_logs/feature_counts.log
        exit 1
    else
        echo "[INFO] Feature counts generated successfully at \$(date)" | tee -a ${config_directory}/0_nextflow_logs/feature_counts.log
    fi

    echo "[INFO] Feature counts process completed at \$(date)" | tee -a ${config_directory}/0_nextflow_logs/feature_counts.log
    """
}


process generate_raw_counts {
    // Process to generate raw counts using HTSeq-count
    debug true
    errorStrategy 'terminate'

    input:
    val config_directory
    val reference_gtf
    val is_paired_end
    val strand_hts
    val paired_hts
    val species
    val feature_counts_files

    output:
    val "${config_directory}/6_raw_counts_output", emit: raw_counts_files

    script:
    """
    echo "[INFO] Starting raw counts generation with HTSeq at \$(date)" | tee -a ${config_directory}/0_nextflow_logs/htseq_counts.log

    if [ ! -f "${config_directory}/generate_raw_counts.sh" ]; then
        echo "[ERROR] HTSeq-count script not found: ${config_directory}/generate_raw_counts.sh" | tee -a ${config_directory}/0_nextflow_logs/htseq_counts.log
        exit 1
    fi

    if [ ! -f "${reference_gtf}" ]; then
        echo "[ERROR] Reference GTF file not found: ${reference_gtf}" | tee -a ${config_directory}/0_nextflow_logs/htseq_counts.log
        exit 1
    fi

    echo "[INFO] Running HTSeq-count with the following parameters:" | tee -a ${config_directory}/0_nextflow_logs/htseq_counts.log
    echo "  Config directory: ${config_directory}" | tee -a ${config_directory}/0_nextflow_logs/htseq_counts.log
    echo "  Reference GTF: ${reference_gtf}" | tee -a ${config_directory}/0_nextflow_logs/htseq_counts.log
    echo "  Paired-end: ${is_paired_end}" | tee -a ${config_directory}/0_nextflow_logs/htseq_counts.log
    echo "  Strand: ${strand_hts}" | tee -a ${config_directory}/0_nextflow_logs/htseq_counts.log
    echo "  Paired HTSeq option: ${paired_hts}" | tee -a ${config_directory}/0_nextflow_logs/htseq_counts.log

    echo 'bash ${config_directory}/generate_raw_counts.sh ${config_directory} ${reference_gtf} ${is_paired_end} ${strand_hts} ${paired_hts} ${species}' | tee -a ${config_directory}/0_nextflow_logs/htseq_counts.log
    bash ${config_directory}/generate_raw_counts.sh ${config_directory} ${reference_gtf} ${is_paired_end} ${strand_hts} ${paired_hts} ${species} | tee -a ${config_directory}/0_nextflow_logs/htseq_counts.log

    if [ "\$?" -ne 0 ]; then
        echo "[ERROR] HTSeq-count generation failed" | tee -a ${config_directory}/0_nextflow_logs/htseq_counts.log
        exit 1
    else
        echo "[INFO] HTSeq-count generated successfully at \$(date)" | tee -a ${config_directory}/0_nextflow_logs/htseq_counts.log
    fi

    echo "[INFO] Raw counts process completed at \$(date)" | tee -a ${config_directory}/0_nextflow_logs/htseq_counts.log
    """
}


process generate_stats {
    // Process to generate stats about the pipeline
    debug true
    errorStrategy 'terminate'

    input:
    val config_directory
    val fastq_files
    val is_paired_end
    val headcrop
    val feature_counts_files

    script:
    """
    echo "[INFO] Starting stats generation at \$(date)" | tee -a ${config_directory}/0_nextflow_logs/7_pipeline_stats.log

    if [ ! -f "${config_directory}/generate_stats.sh" ]; then
        echo "[ERROR] generate_stats.sh script not found in ${config_directory}" | tee -a ${config_directory}/0_nextflow_logs/7_pipeline_stats.log
        exit 1
    fi

    echo "[INFO] Running stats with the following inputs:" | tee -a ${config_directory}/0_nextflow_logs/7_pipeline_stats.log
    echo "  Fastq files: ${fastq_files}" | tee -a ${config_directory}/0_nextflow_logs/7_pipeline_stats.log
    echo "  Paired end: ${is_paired_end}" | tee -a ${config_directory}/0_nextflow_logs/7_pipeline_stats.log
    echo "  Headcrop: ${headcrop}" | tee -a ${config_directory}/0_nextflow_logs/7_pipeline_stats.log

    echo "bash ${config_directory}/generate_stats.sh ${config_directory} ${fastq_files} ${is_paired_end} ${headcrop}" | tee -a ${config_directory}/0_nextflow_logs/7_pipeline_stats.log
    bash ${config_directory}/generate_stats.sh ${config_directory} ${fastq_files} ${is_paired_end} ${headcrop} | tee -a ${config_directory}/0_nextflow_logs/7_pipeline_stats.log

    if [ "\$?" -ne 0 ]; then
        echo "[ERROR] Stats generation failed" | tee -a ${config_directory}/0_nextflow_logs/7_pipeline_stats.log
        exit 1
    else
        echo "[INFO] Stats generation completed successfully at \$(date)" | tee -a ${config_directory}/0_nextflow_logs/7_pipeline_stats.log
    fi
    """
}


workflow {

    // input parameters
    config_directory = params.config_directory
    fastq_files = params.fastq_files
    samples_file = params.samples_file
    fastqc_cores = params.fastqc_cores
    is_paired_end = params.paired_end
    read1_suffix = params.read1_suffix
    read2_suffix = params.read2_suffix
    // hisat2 variables
    hisat2_index = params.hisat2_index
    // trimmomatic variables
    illumina_clip_file = params.illumina_clip_file
    headcrop = params.headcrop
    leading = params.leading
    trailing = params.trailing
    slidingwindow = params.slidingwindow
    minlen = params.minlen
    illuminaclip_params = params.illuminaclip_params
    // filtering variables
    blist_bed_file = params.blacklist_bed_file_path
    exclude_bed_file = params.exclude_bed_file_path
    // reference file and variables used for stringtie / raw-counts
    reference_gtf = params.reference_gtf
    strand_st = params.strand_st
    strand_hts = params.strand_hts
    paired_hts = params.paired_hts

    species = params.species

    // logging
    log.info " "
    log.info "Config directory: ${config_directory}"
    log.info "Fastq files directory: ${fastq_files}"
    log.info " "

    // Assertions to ensure that critical parameters are set
    assert config_directory != null : "ERROR: 'config_directory' must be defined!"
    assert fastq_files != null : "ERROR: 'fastq_files' must be defined!"
    assert samples_file != null : "ERROR: 'samples_file' must be defined!"
    assert fastqc_cores != null : "ERROR: 'fastqc_cores' must be defined!"
    assert is_paired_end != null : "ERROR: 'paired_end' must be defined!"
    assert read1_suffix != null : "ERROR: 'read1_suffix' must be defined!"

    // Running FastQC and MultiQC Reports
    if (params.run_fastqc) {
        log.info "Running FastQC and MultiQC"
        println "FastQC and MultiQC output directory: ${config_directory}/1_fastqc_and_multiqc_reports"
        println "No. of cores to be used for generating FastQC reports: ${fastqc_cores}"
        fastqc_output = generate_fastqc_multiqc_reports(config_directory, fastq_files, fastqc_cores)
    }

    if (params.run_rna_pipeline) {
        assert hisat2_index != null : "ERROR: 'hisat2_index' must be defined!"
        assert reference_gtf != null : "ERROR: 'reference_gtf' must be defined!"
        assert strand_st != null : "ERROR: 'strand_st' must be defined!"
        assert strand_hts != null : "ERROR: 'strand_hts' must be defined!"
        assert paired_hts != null : "ERROR: 'paired_hts' must be defined!"
        assert illumina_clip_file != null : "ERROR: 'illumina_clip_file' must be defined!"
    }

    // Running Psomagen RNA Pipeline
    if (params.run_rna_pipeline) {
        // Trimmomatic trimming
        log.info "Running Trimmomatic Trimming"
        println "Trimmomatic output directory: ${config_directory}/2_trim_output"
        trimmed_files = trimmomatic_trimming(
            config_directory, samples_file, fastq_files, is_paired_end, read1_suffix, read2_suffix,
            fastqc_cores, illumina_clip_file, headcrop, leading, trailing, slidingwindow, minlen, illuminaclip_params
        )

        // HISAT2 Mapping
        log.info "Running HISAT2 Mapping"
        println "HISAT2 mapping output directory: ${config_directory}/3_hisat2_mapping_output"
        mapped_files = hisat2_mapping(
            config_directory, samples_file, is_paired_end, read1_suffix, read2_suffix,
            fastqc_cores, hisat2_index, trimmed_files.trimmed_files
        )

        // Generating mapping metrics
        log.info "Generating Mapping Metrics"
        map_metric_files = generate_mapping_metrics(config_directory, mapped_files.mapped_files)

        // Filtering, deduplication and indexing
        log.info "Running Filtering and Deduplication"
        println "Filtering, deduplication and sorting output directory: ${config_directory}/4_filter_output"
        filtered_files = filter_samples(config_directory, fastqc_cores, blist_bed_file, exclude_bed_file, mapped_files.mapped_files)

        // Generate Qualimap reports
        log.info "Running Qualimap for quality control"
        println "Generating qualimap reports for filtered BAM files: ${config_directory}/4_1_qualimap_filter_output_qc"
        qualimap_files = generate_qualimap_reports(config_directory, reference_gtf, is_paired_end, fastqc_cores, filtered_files.filtered_files)

        // Generating stringtie and raw-counts
        log.info "Running StringTie and Raw Counts Generation"
        println "Generating StringTie and raw counts in: ${config_directory}/5_stringtie_counts_output and ${config_directory}/6_raw_counts_output"
        stringtie_files = generate_stringtie_counts(config_directory, fastqc_cores, reference_gtf, strand_st, species, filtered_files.filtered_files)
        feature_count_files = generate_feature_counts(config_directory, fastqc_cores, reference_gtf, is_paired_end, filtered_files.filtered_files)
        raw_count_files = generate_raw_counts(config_directory, reference_gtf, is_paired_end, strand_hts, paired_hts, species, feature_count_files.feature_counts_files)

        // Generate Stats
        log.info "Generating Pipeline Stats"
        println "Generating Stats for the current run: ${config_directory}/7_pipeline_stats_<TIMESTAMP>.log"
        generate_stats(config_directory, fastq_files, is_paired_end, headcrop, feature_count_files.feature_counts_files)

    }
}
