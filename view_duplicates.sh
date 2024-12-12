#!/bin/bash

# 2024-11-22
# Barbara Walkowiak bw18

# I want to view duplicates in bam files for all 3434 samples using samtools 
# first cd /lustre/scratch126/casm/team274sb/bw18/twins_sarcoma/
# remember you have to load samtools if you want this script to run!

module load samtools-1.14/python-3.12.0

echo "try to do this for one sample" >> duplicates/dup_report.txt
samtools view /nfs/cancer_ref01/nst_links/live/3434/PD62341aa/PD62341aa.sample.dupmarked.bam | head >> duplicates/dup_report.sam
samtools flagstat /nfs/cancer_ref01/nst_links/live/3434/PD62341aa/PD62341aa.sample.dupmarked.bam | grep 'duplicates' >> duplicates/dup_report.txt
echo "done" >> duplicates/dup_report.txt

# for each sample, identify correct bam file and use samtools to view duplicates 
PATIENT_IDS=("PD62341aa" "PD62341ad" "PD62341h" "PD62341n" "PD62341q" "PD62341v" "PD63383t" "PD63383u" "PD63383w" "PD63383ae" "PD63383ak" "PD63383bb" "PD62341ae" "PD62341ag" "PD62341aj" "PD62341ak" "PD62341ap" "PD62341am" "PD62341b" "PD62341u" "PD63383ap" "PD63383aq")

for SAMPLE in "${PATIENT_IDS[@]}"; do
    BAM="/nfs/cancer_ref01/nst_links/live/3434/${SAMPLE}/${SAMPLE}.sample.dupmarked.bam"
    echo "Processing bam file $BAM" >> duplicates/dup_report.txt
    samtools view -f 0x400 "$BAM" | head >> duplicates/dup_report.txt
    samtools flagstat "$BAM" | grep 'duplicates' >> duplicates/dup_report.txt 
done 
