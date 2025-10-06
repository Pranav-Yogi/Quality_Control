#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ────────────────────────────────
// PARAMETERS
// ────────────────────────────────
params.raw_glob     = "{Path to your Raw data}/00_raw/*_R{1,2}.fastq.gz"                          // Path for Raw data
params.fastqc_out   = "{Path to your Deliverables directory}/01_QC/raw_read_qc/fastqc_reports"    // Path for QC report
params.fastqc_img   = "/home/cm/tools/fastqc_0.12.1.sif"                                          // Docker Image OR Singularity container path of fastqc
params.fastp_img    = "/home/cm/tools/fastp_0.23.4.sif"                                           // Docker Image OR Singularity container path of Fastp
params.min_length   = 50                                                                          // Minimum lenth for fastp trim 
params.deliverables = "{Path to your Deliverables directory}/02_Processedreads"                   // Path of Processed data

// ────────────────────────────────
// CHANNEL: PAIR R1/R2
// ────────────────────────────────
// CORRECTION: Removed the invalid 'flat' and 'fileExtension' parameters.
Channel
    .fromFilePairs(params.raw_glob)
    .set { read_pairs_ch }

// ────────────────────────────────
// PROCESS 1: FASTQC
// ────────────────────────────────
process FastQC_pair {
    tag { sampleId }

    cpus 3
    memory '4 GB'

    input:
    tuple val(sampleId), path(reads)

    publishDir "${params.fastqc_out}/${sampleId}", mode: 'move'

    output:
    tuple val(sampleId), path(reads)        // pass reads forward
    file "*.html"
    file "*.zip"

    // make sure output dir exists before running fastqc
    beforeScript "mkdir -p ${params.fastqc_out}/${sampleId}"

    script:
    """
    echo "Running FastQC for sample: ${sampleId}"
    singularity exec ${params.fastqc_img} fastqc -t ${task.cpus} ${reads.join(' ')} -o ./
    """
}

// ────────────────────────────────
// PROCESS 2: FASTP
// ────────────────────────────────
process Fastp_trimming {
    tag { sampleId }

    cpus 4
    memory '8 GB'

    input:
    tuple val(sampleId), path(reads)

    publishDir params.deliverables, mode: 'move'

    output:
    tuple val(sampleId), path("processed_data/${sampleId}_R1.fastp.fastq.gz"), path("processed_data/${sampleId}_R2.fastp.fastq.gz")
    file "fastp_reports/${sampleId}_fastp.html"
    file "fastp_reports/${sampleId}_fastp.json"
    file "logs/${sampleId}_fastp.log"
    file "logs/${sampleId}_processed_fastp.ok"

    script:
    """
    set -euo pipefail

    mkdir -p processed_data fastp_reports logs

    R1_out=processed_data/${sampleId}_R1.fastp.fastq.gz
    R2_out=processed_data/${sampleId}_R2.fastp.fastq.gz
    html_report=fastp_reports/${sampleId}_fastp.html
    json_report=fastp_reports/${sampleId}_fastp.json
    logpath=logs
    fastpok=\${logpath}/${sampleId}_processed_fastp.ok

    echo "[\$(date)] FASTP: ${sampleId} (R1=${reads[0]}, R2=${reads[1]})"

    /usr/bin/time -v singularity exec ${params.fastp_img} fastp \\
        --in1 ${reads[0]} \\
        --in2 ${reads[1]} \\
        --out1 \${R1_out} \\
        --out2 \${R2_out} \\
        --thread ${task.cpus} \\
        --length_required ${params.min_length} \\
        --correction \\
        --trim_poly_g \\
        --cut_front \\
        --cut_tail \\
        --qualified_quality_phred 30 \\
        --unqualified_percent_limit 30 \\
        --average_qual 30 \\
        --html \${html_report} \\
        --json \${json_report} \\
        --report_title ${sampleId}_fastp_report 2>&1 | tee \${logpath}/${sampleId}_fastp.log \\
        && touch \${fastpok}
    """
}

// ────────────────────────────────
// WORKFLOW
// ────────────────────────────────
workflow {
    FastQC_pair(read_pairs_ch) 
    Fastp_trimming(read_pairs_ch)
}
