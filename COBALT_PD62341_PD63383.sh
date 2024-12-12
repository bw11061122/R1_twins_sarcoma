#!/bin/bash

# 2024-11-19
# Barbara Walkowiak bw18
# Script taken from https://github.com/BehjatiLab/PURPLE/blob/main/run_cobalt2.sh with my modifications

#You need to remember to change the required memory per job
#Each sample takes around 12GB

# Load R and set the library path
module load R
export R_LIBS="/lustre/scratch125/casm/team274sb/ah39/Rlib"

# List of patient IDs
PATIENT_IDS=("PD62341aa" "PD62341ad" "PD62341h" "PD62341n" "PD62341q" "PD62341v" "PD63383t" "PD63383u" "PD63383w" "PD63383ae" "PD63383ak" "PD63383bb" "PD62341ae" "PD62341ag" "PD62341aj" "PD62341ak" "PD62341ap" "PD62341am" "PD62341b" "PD62341u" "PD63383ap" "PD63383aq")

# Base directories
BAM_BASE_DIR="/nfs/cancer_ref01/nst_links/live/3434"
OUTPUT_BASE_DIR="/lustre/scratch126/casm/team274sb/bw18/twins_sarcoma/purple/cobalt_output"
GC_PROFILE="/lustre/scratch126/casm/team274sb/ah39/DNA_data/sim_links/scLeuk_3.7.24/GC_profile.1000bp.38.cnp"
DIPLOID_BED="/lustre/scratch126/casm/team274sb/ah39/DNA_data/sim_links/scLeuk_3.7.24/DiploidRegions.38.bed.gz"
COBALT_JAR="/lustre/scratch125/casm/team274sb/ah39/code_packages/cobalt-1.14.1.jar"

# Set other variables
THREADS=2

# Loop through each patient ID
for TUMOR in "${PATIENT_IDS[@]}"; do
    # Define paths specific to each patient
    TUMOR_BAM="${BAM_BASE_DIR}/${TUMOR}/${TUMOR}.sample.dupmarked.bam"
    OUTPUT_DIR="${OUTPUT_BASE_DIR}/${TUMOR}"

    # Create output directory if it doesn't exist
    mkdir -p $OUTPUT_DIR

    # Run COBALT in tumor-only mode for each patient
    echo "Submitting job for patient: $TUMOR"
    java -Xmx28G -jar $COBALT_JAR \
        -tumor $TUMOR \
        -tumor_bam $TUMOR_BAM \
        -output_dir $OUTPUT_DIR \
        -threads $THREADS \
        -gc_profile $GC_PROFILE \
        -tumor_only_diploid_bed $DIPLOID_BED
done

# tumour_only_diploid_bed = Bed file of diploid regions of the genome