#!/bin/bash

# 2024-11-22
# Barbara Walkowiak bw18

# I want to view all reads that came from chr1 that span 
# first cd /lustre/scratch126/casm/team274sb/bw18/twins_sarcoma/
# remember you have to load samtools if you want this script to run!

module load samtools-1.14/python-3.12.0

# 1 I want to identify the number of duplicate reads for each file 
# for each sample, identify correct bam file and use samtools to view duplicates 
PATIENT_IDS=("PD62341aa" "PD62341ad" "PD62341h" "PD62341n" "PD62341q" "PD62341v" "PD63383t" "PD63383u" "PD63383w" "PD63383ae" "PD63383ak" "PD63383bb" "PD62341ae" "PD62341ag" "PD62341aj" "PD62341ak" "PD62341ap" "PD62341am" "PD62341b" "PD62341u" "PD63383ap" "PD63383aq")

for SAMPLE in "${PATIENT_IDS[@]}"; do
    BAM="/nfs/cancer_ref01/nst_links/live/3434/${SAMPLE}/${SAMPLE}.sample.dupmarked.bam"
    echo "Processing bam file $BAM" >> duplicates/dup_report.txt
    samtools flagstat "$BAM" | grep 'duplicates' >> duplicates/dup_report.txt 
done 

# 2 I want to get out all reads that span chr1_38827952 and write them to separate files so I can inspect them without staring at Jbrowse 
for SAMPLE in "${PATIENT_IDS[@]}"; do
    BAM="/nfs/cancer_ref01/nst_links/live/3434/${SAMPLE}/${SAMPLE}.sample.dupmarked.bam"
    echo "Processing bam file $BAM" >> duplicates/chr1_388_reads_all.sam
    echo "Processing bam file $BAM" >> duplicates/chr1_388_reads_duplicates.sam
    samtools view "$BAM" chr1:38827952-38827952 >> duplicates/chr1_388_reads_all.sam
    samtools view -f 1024 "$BAM" chr1:38827952-38827952 >> duplicates/chr1_388_reads_duplicates.sam 
done 
