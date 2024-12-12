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
PATIENT_IDS=("PD62341aa" "PD62341ad" "PD62341h" "PD62341n" "PD62341q" "PD62341v" "PD63383t" "PD63383u" "PD63383w" "PD63383ae" "PD63383ak" "PD63383bb" "PD62341ae" "PD62341ag" "PD62341aj" "PD62341ak" "PD62341ap" "PD62341am" "PD62341b" "PD62341u" "PD63383ap" "PD63383aq")

# Set common variables
AMBER_BASE_DIR="/lustre/scratch126/casm/team274sb/bw18/twins_sarcoma/purple/amber_output"
COBALT_BASE_DIR="/lustre/scratch126/casm/team274sb/bw18/twins_sarcoma/purple/cobalt_output"
GC_PROFILE="/lustre/scratch126/casm/team274sb/ah39/DNA_data/sim_links/scLeuk_3.7.24/GC_profile.1000bp.38.cnp"
REF_GENOME_VERSION="38"
REF_GENOME="/lustre/scratch124/casm/team78pipelines/canpipe/live/ref/human/GRCh38_full_analysis_set_plus_decoy_hla/genome.fa"
ENSEMBL_DATA_DIR="/lustre/scratch126/casm/team274sb/ah39/DNA_data/sim_links/scLeuk_3.7.24/ensembl_data"
PURPLE_JAR="/lustre/scratch125/casm/team274sb/ah39/code_packages/purple_v3.8.4.jar"
CIRCOS_BIN="/lustre/scratch125/casm/team274sb/ah39/code_packages/circos-0.69-9/bin/circos"

# Loop over patient IDs
for TUMOR in "${PATIENT_IDS[@]}"; do
    # Construct paths based on current patient ID
    AMBER_DIR="${AMBER_BASE_DIR}/${TUMOR}"
    COBALT_DIR="${COBALT_BASE_DIR}/${TUMOR}"
    OUTPUT_DIR="/lustre/scratch126/casm/team274sb/bw18/twins_sarcoma/purple/purple_output/${TUMOR}/without_GRIDDS"

    # Create output directory if it doesn't exist
    mkdir -p $OUTPUT_DIR

    # Run PURPLE for current patient ID
    java -jar $PURPLE_JAR \
      -tumor $TUMOR \
      -amber $AMBER_DIR \
      -cobalt $COBALT_DIR \
      -gc_profile $GC_PROFILE \
      -ref_genome_version $REF_GENOME_VERSION \
      -ref_genome $REF_GENOME \
      -ensembl_data_dir $ENSEMBL_DATA_DIR \
      -output_dir $OUTPUT_DIR \
      -circos $CIRCOS_BIN

    echo "PURPLE job submitted for patient ID: $TUMOR"
done

echo "All PURPLE jobs submitted successfully."