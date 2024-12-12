#!/bin/bash

# 2024-11-19
# Barbara Walkowiak bw18
# Script taken from https://github.com/BehjatiLab/PURPLE/blob/main/run_amber2.sh with my modifications

# This is with the option to add in multiple samples
# Adjust memory as required per sample - around 18GB
# This is for when you don't want to use GRIDDS in your output. 

# Load R and set the library path
module load R
export R_LIBS="/lustre/scratch126/casm/team274sb/ah39/Rlib" # from Angus 

# Define patient IDs (project 3434)
PATIENT_IDS=("PD62341aa" "PD62341ad" "PD62341h" "PD62341n" "PD62341q" "PD62341v" "PD63383t" "PD63383u" "PD63383w" "PD63383ae" "PD63383ak" "PD63383bb")

# Set common variables
AMBER_BASE_DIR="/lustre/scratch126/casm/team274sb/bw18/twins_sarcoma/purple/amber_output"
COBALT_BASE_DIR="/lustre/scratch126/casm/team274sb/bw18/twins_sarcoma/purple/cobalt_output"
GC_PROFILE="/lustre/scratch126/casm/team274sb/ah39/DNA_data/sim_links/scLeuk_3.7.24/GC_profile.1000bp.38.cnp"
REF_GENOME_VERSION="38"
REF_GENOME="/lustre/scratch124/casm/team78pipelines/canpipe/live/ref/human/GRCh38_full_analysis_set_plus_decoy_hla/genome.fa"
ENSEMBL_DATA_DIR="/lustre/scratch126/casm/team274sb/ah39/DNA_data/sim_links/scLeuk_3.7.24/ensembl_data"
GERMLINE_HOTSPOTS="/nfs/users/nfs_b/bw18/data/purple/KnownHotspots.germline.38.vcf.gz"
GERMLINE_DEL="/nfs/users/nfs_b/bw18/data/purple/cohort_germline_del_freq.38.csv"
PURPLE_JAR="/lustre/scratch125/casm/team274sb/ah39/code_packages/purple_v3.8.4.jar"
CIRCOS_BIN="/lustre/scratch125/casm/team274sb/ah39/code_packages/circos-0.69-9/bin/circos"

# Loop over patient IDs
for SAMPLE in "${PATIENT_IDS[@]}"; do
    # Construct paths based on current patient ID
    AMBER_DIR="${AMBER_BASE_DIR}/${SAMPLE}"
    COBALT_DIR="${COBALT_BASE_DIR}/${SAMPLE}"
    OUTPUT_DIR="/lustre/scratch126/casm/team274sb/bw18/twins_sarcoma/purple/purple_output/${SAMPLE}/without_GRIDDS_germline"

    # Create output directory if it doesn't exist
    mkdir -p $OUTPUT_DIR

    # Run PURPLE for current patient ID
    java -jar $PURPLE_JAR \
        -reference $SAMPLE \
        -amber $AMBER_DIR \
        -cobalt $COBALT_DIR \
        -gc_profile $GC_PROFILE \
        -ref_genome_version $REF_GENOME_VERSION \
        -ref_genome  $REF_GENOME \
        -ensembl_data_dir $ENSEMBL_DATA_DIR \
        -germline_hotspots $GERMLINE_HOTSPOTS \
        -germline_del_freq_file $GERMLINE_DEL \
        -output_dir $OUTPUT_DIR \

    echo "PURPLE job submitted for patient ID: $SAMPLE"

done

echo "All PURPLE jobs submitted successfully."