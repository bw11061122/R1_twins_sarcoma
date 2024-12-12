#!/bin/bash

# 2024-11-19
# Barbara Walkowiak bw18
# Script taken from https://github.com/BehjatiLab/PURPLE/blob/main/run_amber2.sh with my modifications

#This is with option to add in multiple samples 
#is cira 12GB per sample required - so 15GB is safe

# Load R and set the library path
module load R
export R_LIBS="/lustre/scratch125/casm/team274sb/ah39/Rlib"

# Define your patient IDs
PATIENT_IDS=("PD62341aa" "PD62341ad" "PD62341h" "PD62341n" "PD62341q" "PD62341v" "PD63383t" "PD63383u" "PD63383w" "PD63383ae" "PD63383ak" "PD63383bb" "PD62341ae" "PD62341ag" "PD62341aj" "PD62341ak" "PD62341ap" "PD62341am" "PD62341b" "PD62341u" "PD63383ap" "PD63383aq")

# Set common variables
BAM_BASE_DIR="/nfs/cancer_ref01/nst_links/live/3434"
OUTPUT_BASE_DIR="/lustre/scratch126/casm/team274sb/bw18/twins_sarcoma/purple/amber_output"
LOCI="/lustre/scratch126/casm/team274sb/ah39/DNA_data/sim_links/scLeuk_3.7.24/GermlineHetPon.38.vcf.gz"
AMBER_JAR="/lustre/scratch125/casm/team274sb/ah39/code_packages/amber-3.9.jar"
THREADS=2
REF_GENOME_VERSION="V38"

# Loop over patient IDs
for TUMOR in "${PATIENT_IDS[@]}"; do
    # Construct paths based on current patient ID
    TUMOR_BAM="${BAM_BASE_DIR}/${TUMOR}/${TUMOR}.sample.dupmarked.bam"
    OUTPUT_DIR="${OUTPUT_BASE_DIR}/${TUMOR}"

    # Create output directory if it doesn't exist
    mkdir -p $OUTPUT_DIR

    # Run AMBER for current patient ID
    java -Xmx32G -cp $AMBER_JAR com.hartwig.hmftools.amber.AmberApplication \
      -tumor $TUMOR \
      -tumor_bam $TUMOR_BAM \
      -output_dir $OUTPUT_DIR \
      -threads $THREADS \
      -loci $LOCI \
      -tumor_only_min_vaf 0.1 \
      -tumor_only_min_support 5 \
      -ref_genome_version $REF_GENOME_VERSION

    echo "AMBER job submitted for patient ID: $TUMOR"
done

# tumour_only_min_vaf = min VAF required to call it as present in the data 
# tumour_only_min_support = min nr of reads reporting the variant allele to call it as present in the data 

echo "All AMBER jobs submitted successfully."