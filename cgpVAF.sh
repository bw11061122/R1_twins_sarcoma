
# Barbara Walkowiak bw 18
# 2024-10-09
https://confluence.sanger.ac.uk/pages/viewpage.action?pageId=22710418

# Instructions on how to run the pileup (SNV) 
# The goal of the pileup is to create a matrix 
# Each row represents a mutation (identified in any of the samples)
# Each column is a sample that was examined in the study 
# Based on this, you can identify mutations present in a set of samples 
# Note that pileup is created from SNVs, exonerate is created from indels (the same logic)

# NOTE: your data has to be staged correctly 
https://confluence.sanger.ac.uk/display/CAS/Staging+files 
# use the stageBam.pl command to get your files on luster first 

# First we want to get the bed file with mutations to look at
# You can do this by: 
# moving to the directory where you have all CaVEMan files
# (if your files are gzipped, you need to gunzip them)
# creating sub files wth only mutations that passed all filters (grep PASS)
# from all files , cat -f 1,2,4,5 (chr, pos, ref, alt) and write to a single bed file 
# Done!
cd /nfs/cancer_ref01/nst_links/live/3434
ls
gunzip *CaVEMan.gzip
*.caveman | gunzip *.gz | grep 'PASS' | cat -f 1,2,4,5 > caveman_subs.bed

module load cgpVAFcommand # you need to load the module of interest to access the commands 

createVafCmd.pl -pid 3434 -o testdir 
# -pid is project ID
# -o is the output directory
# select caveman_c (caveman_java is an older implementation)
# this generates three files (run_bsub.sh, cpgVafChr.cmd and cpgVafConcat.cmd)

# check the command you have:
cgpVaf.pl -d tmp -o tmp \
-g /lustre/scratch119/casm/team78pipelines/canpipe/live/ref/human/GRCh37d5/genome.fa \
-hdr /lustre/scratch119/casm/team78pipelines/canpipe/live/ref/human/GRCh37d5/shared/ucscHiDepth_0.01_mrg1000_no_exon_coreChrs.bed.gz \
-a snp \
-b twins_unm_subs.bed
-nn PD23583b2 -tn PD62341  \
-pid 3434

# g: path to the reference human genome: make sure to use the right version (37 vs 38)
# hdr: path to high depth regions to exclude: need to be gzipped and tbi (indexing system)
# mq: mapping quality
# bq: base quality
# nn: normal sample (if your provide the bed file, it is okay if this is sth else to the in silico ref genome)
# tn: tumour sample (or normal sample in which to find mutations)

# Henry recommends running with different parameters (mq, bq) (you don't want to be too stringent)
# e.g., you don't want to label mutations that are present in all samples as 'absent' in some due to poor mapping quality

bash run_bsub.sh
# you want to provide the bed file (chr, position, ref, alt)
# the bed file has mutations you want to look at 
# if you have lots of samples, you may need to edit the run_bsub.sh file (e.g., change memory requirements)

# understanding the output
# you are effectively interested in:
# MRT - number of reads with the variant allele
# WRT - number of reads with the wild-type allele
# TOTAL DEPTH = MRT + WRT
# VAF = MRT / (MRT + WRT)
# in the pileup, you get several columns for each sample with different filters / relevant parameters 

# TRYING TO FIGURE THIS OUT 23/10/2024
log in to the farm

# create directories to store the output in 
cd /lustre/scratch126/casm/team274sb
mkdir -p bw18/twins_sarcoma/subs/subs_pileup_mq5_bq5_20241023

cd /nfs/cancer_ref01/nst_links/live/3434 # this is where the sample data is stored
gunzip *CaVEMan.gzip
*.caveman | grep 'PASS' | cat -f 1,2,4,5 > /lustre/scratch126/casm/team274sb/bw18/twins_sarcoma/subs/subs_pileup_mq5_bq5_20241023/twins_unm_subs.bed

cd /lustre/scratch126/casm/team274sb/bw18/twins_sarcoma/subs/subs_pileup_mq1_bq1_20241023
module load cgpVAFcommand
createVafCmd.pl -pid 3434 -o testdir
bash run_bsub.sh

# I cannot find the script in which the cpgVaf.pl command is present 
cgpVaf.pl -d tmp -o tmp \
-g /lustre/scratch119/casm/team78pipelines/canpipe/live/ref/human/GRCh37d5/genome.fa \
-hdr /lustre/scratch119/casm/team78pipelines/canpipe/live/ref/human/GRCh37d5/shared/ucscHiDepth_0.01_mrg1000_no_exon_coreChrs.bed.gz \
-a snp \
-b twins_unm_subs.bed
-nn PD23583b2 -tn PD62341  \
-pid 3434