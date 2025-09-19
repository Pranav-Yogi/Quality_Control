#!/usr/bin/bash

# Get user input for necessary paths and parameters
echo -n "Enter Deliverables Directory Path: "
read deliverables_dir

echo -n "Enter Rawdata Directory Path: "
read rwDir

# Get file extensions for R1 and R2 reads
for i in `ls ${rwDir} | grep _R1`;do echo $i;done
echo -n "Please look into the read names and provide the extension after sample name Example: _R1_001.fastq.gz: "
read "R1_extn"

for i in `ls ${rwDir} | grep _R2`;do echo $i;done
echo -n "Please look into the read names and provide the extension after sample name Example: _R2_001.fastq.gz: "
read "R2_extn"

# Get resource allocation details
echo -n "number of threads: "
read cpu

echo -n "Minimum length after trimming : "
read min_length

# Create the required directory structure
mkdir -p $deliverables_dir/00_logs
mkdir -p $deliverables_dir/01_QC/processed_read_qc/fastp_reports/
mkdir -p $deliverables_dir/02_processed_data/
mkdir -p $deliverables_dir/01_QC/processed_read_qc/fastqc_reports


# Define log path & adapter path & tool path
logpath="${deliverables_dir}/00_logs"
fastptool="{Path to your container}/fastp_0.23.4.sif"     # Provide path of your singularity fastp container
fastqctool="{Path to your container}/fastqc_0.12.1.sif"   # Provide path of your singularity fastqc container

# adapters="Path of adapter fasta if you have"

# Loop through each sample based on the R1 file extension
for i in $(ls $rwDir | grep ${R1_extn}); do
    sample=`basename $i ${R1_extn}`
    echo "----------------------------------------"
    echo "Processing sample: $sample"
    echo "----------------------------------------"

    # Define input and output file paths for the current sample
    R1_in=${rwDir}/${sample}${R1_extn}
    R2_in=${rwDir}/${sample}${R2_extn}

    R1_out=${deliverables_dir}/02_processed_data/${sample}_R1.fastp.fastq.gz
    R2_out=${deliverables_dir}/02_processed_data/${sample}_R2.fastp.fastq.gz

    html_report="${deliverables_dir}/01_QC/processed_read_qc/fastp_reports/${sample}_fastp.html"
    json_report="${deliverables_dir}/01_QC/processed_read_qc/fastp_reports/${sample}_fastp.json"

    # Define fastqc output directory
    mkdir -p ${deliverables_dir}/01_QC/processed_read_qc/fastqc_reports/$sample
    fastqc_out=${deliverables_dir}/01_QC/processed_read_qc/fastqc_reports/$sample

    # --- Step 1: Run fastp for trimming and quality filtering ---
    fastpok=${logpath}/${sample}_processed_fastp.ok
    if [ ! -f ${fastpok} ]; then
        echo "Running fastp for ${sample}..."
        /usr/bin/time -v ${fastptool} fastp \
            --in1 ${R1_in} \
            --in2 ${R2_in} \
            --out1 ${R1_out} \
            --out2 ${R2_out} \
            --thread ${cpu} \
            --length_required ${min_length} \
            --correction \
            --trim_poly_g \
            --cut_front \
            --cut_tail \
            --qualified_quality_phred 30 \
            --unqualified_percent_limit 30 \
            --average_qual 30 \
            --html ${html_report} \
            --json ${json_report} \
            --report_title ${sample}_fastp_report 2>&1 | tee ${logpath}/${sample}_fastp.log && touch ${logpath}/${sample}_processed_fastp.ok
    else
        echo "fastp already completed for ${sample}. Skipping."
    fi

    # --- Step 2: Run FastQC on the processed data from fastp ---
    # This step will only run if the fastp step was successful (the .ok file exists)
    if [ -f "${fastpok}" ]; then
        echo "Running FastQC on processed data for ${sample}..."
        /usr/bin/time -v ${fastqctool} fastqc \
            ${R1_out} \
            ${R2_out} \
            -o ${fastqc_out} \
            -t ${cpu} 2>&1 | tee ${logpath}/${sample}_processed_data_fastqc.log && touch ${logpath}/${sample}_processed_data_fastqc.ok
    else
        echo "fastp failed for ${sample}. Skipping FastQC."
    fi
done

echo "----------------------------------------"
echo "All samples processed."
echo "----------------------------------------"
