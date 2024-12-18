# README

#########################################################################################
#########################################################################################

2024-10-11

# I downloaded all files from Henry Lee-Six from google drive (https://drive.google.com/drive/folders/1fcJSLfVH2AzwqqvmtmYHDhUED6Jw1_6Q)

# I first checked the number of lines in each file 
cd ~/Desktop/1SB/Analysis/Data # you are in the data folder 
wc -l *.tsv

# Output:
360941 PD62341v_PD62341q_snp_vaf.tsv
358938 PD63383w_PD63383t_snp_vaf.tsv
360933 PDv38is_wgs_PD62341v_snp_vaf.tsv
360933 PDv38is_wgs_PD63383w_snp_vaf.tsv

# Next I removed the header (in case the differences are due to differences in the header)
bash ../scripts/20241011_removeheader.sh

wc -l *_header.tsv

# Output:
360884 PD62341v_PD62341q_snp_vaf_removed_header.tsv
358883 PD63383w_PD63383t_snp_vaf_removed_header.tsv
360884 PDv38is_wgs_PD62341v_snp_vaf_removed_header.tsv
360884 PDv38is_wgs_PD63383w_snp_vaf_removed_header.tsv

# There is still sth odd about the 2nd file 
# NOTE: I redownloaded the files and repeated this in case sth went wrong with the download
# However, I got the same result again, so I think it must be sth wrong with the file itself

# I decided to read all files in R and merge them in R 
# I also wanted to check what was wrong with the second file 

# In R, there is an error when reading the second file (PD63383w_PD63383t_snp_vaf_removed_header.tsv):
# Warning message:
# In scan(file = file, what = what, sep = sep, quote = quote, dec = dec,  :embedded nul(s) found in input

# In this file, there is a row where the data is not read correctly (row 58041, see R script)
# Look at this row in the original data
# (NOTE: need to change row index because the original file also has a header)
sed -n '58040,58042p;58042q' PD63383w_PD63383t_snp_vaf_removed_header.tsv
sed -n '58097,58099p;58100q' PD63383w_PD63383t_snp_vaf.tsv

# I think the problem is likely that there is an entry '2PD63383w' appearing in a different column
# (It looks like one row is interpreted as following on for the previous row)
# I edited the file (moved the entry to a new line), however, this likely has happened because some data is missing
# Therefore, I need to check on Monday 14/10/2024 if something went wrong with the pileup such that the file is corrupt
# See R scripts (20241011_pileup_checks) for details on how the file was processed as of 2024/10/12

# I think there must be also some other issues with how this file was read
# Because some PD63383bb columns expected to be numeric contain non-numeric characters (?)

#########################################################################################
#########################################################################################
2024-10-16 

# OKAY so Henry re-run the pileup and now I have all 22 samples 
# which all should be complete (no issues with missing columns / rows etc)

run:
cd ~/Desktop/1SB/Analysis/Data # you are in the data folder 
bash ../scripts/20241011_removeheader.sh
wc -l *_header.tsv

output: 
362815 PD62341v_PD62341q_snp_vaf_removed_header.tsv
362815 PD63383w_PD63383t_snp_vaf_removed_header.tsv
362815 PDv38is_wgs_PD62341b_snp_vaf_removed_header.tsv
362815 PDv38is_wgs_PD63383w_snp_vaf_removed_header.tsv
1451260 total

#########################################################################################
#########################################################################################
2024-10-17 

# I want to move everything to Jupyterhub

# For now, I will move the merged pileup to my home directory on nfs to be able to work with it from where 
# I can actually access it and work with it 
cd Desktop/1SB/Analysis/Data
rsync -avhu pileup_merged_20241016.tsv bw18@gen3:3434_data/
rsync -avhu PD62341_PD63382_unmatched_indels_pass_or_fail.txt bw18@gen3:3434_data/
rsync -avhu PD62341_PD63382_unmatched_pass_indels.txt bw18@gen3:3434_data/

#########################################################################################
#########################################################################################
2024-10-19

# FILE ORGANIZATION
# Scripts:
20241011_pileup_checks: merge the pileups from HLS into one file
20241011_pileup_filters: create a list of mutations that pass filters and can be used for phylogeny
20241018_pileup_qcplots: output QC plots (basic VAF, coverage etc.) and misc checks 
20241018_pileup_mutTumour: output groups of muts that fall into categories useful for phylo reconstruction

#########################################################################################
#########################################################################################
2024-10-23

# I needed to get github set up and this was extremely nightmarish so currently I am using Github Desktop App
current organisation of the folders is:
1SB - everything that relates to your rotation
In this folder, I have:
- Data (all relevant data used)
- Scripts/R1 - R1 is a github repo where I track all relevant scripts
- Results (pdfs with results figures)
- Report (files relevant to my report + overleaf report format)
- Misc (everything else to keep things clean)
- We will need to open the Github Desktop app each time we want to commit stuff and add the changes with commit 

#########################################################################################
#########################################################################################
2024-10-24

# I am trying to run the low quality pileup (mq 5, bq 5)

cd /nfs/cancer_ref01/nst_links/live/3434
ls

# reading the GitHub page of CaVEMan, these are the output files you can see in 
PD62341aa.caveman_c.algbean			    
PD62341aa.caveman_c.annot.vcf.gz		    
PD62341aa.caveman_c.annot.vcf.gz.tbi		    
PD62341aa.caveman_c.cov_arr	- generated by the maximization + merge step 			    
PD62341aa.caveman_c.flag.vcf.gz			    
PD62341aa.caveman_c.flag.vcf.gz.tbi		    
PD62341aa.caveman_c.merged_split - output of the first step which creates a list of segments to be analysed 		    
PD62341aa.caveman_c.no.analysis.bed.gz	- exclude from analysis (likely poor coverage or mapping)	   
PD62341aa.caveman_c.no.analysis.bed.gz.tbi - exclude from analysis (likely poor coverage or mapping)	
PD62341aa.caveman_c.prob_arr - generated by the maximization + merge step 	    
PD62341aa.caveman_c.snps.vcf.gz	- generated by the expectation step: likely SNPs 		    
PD62341aa.caveman_c.snps.vcf.gz.tbi	- generated by the expectation step: likely SNPs 	
PD62341aa.caveman_c.vcf.gz - I presume these are the likely somatic mutations (?)			    
PD62341aa.caveman_c.vcf.gz.tbi - I presume these are the likely somatic mutations (?)

# tbi is an index file for tab-delimited files 

# now, for each folder you want to access the *.caveman_c.annot.vcf.gz file and do 
cd /nfs/cancer_ref01/nst_links/live/3434
for dir in */; do gunzip -c "$dir"/*.annot.vcf.gz > /lustre/scratch126/casm/team274sb/bw18/twins_sarcoma/subs/"$(basename "${dir"%/}")"_annot_vcf; done 
grep 'PASS' *_annot_vcf | cut -f 1,2,4,5 >> twins_subs.bed
sed 's/^[^:]*://g' twins_subs.bed > twins_unm_subs.bed
sort twins_unm_subs.bed | uniq > twins_unm_subs.bed 
# something went wrong given that I have ~3.5k unique subs, that cannot be true 

# okay I think I give up and will just copy the file from hl11 folder
cd /lustre/scratch126/casm/team274sb/bw18/twins_sarcoma/subs
cp ../../../hl11/twins_sarcoma/subs/subs_unmatched/twins_unm_subs.bed twins_unm_subs.bed

#########################################################################################
2024-10-25
# ok trying to understand how to actually do this 

ssh farm22
cd /lustre/scratch126/casm/team274sb/bw18
mkdir -p twins_sarcoma/subs/pileup_20241025
cp ../hl11/twins_sarcoma/subs/subs_unmatched/twins_unm_subs.bed twins_sarcoma/subs/twins_unm_subs.bed
cd twins_sarcoma/subs
source /software/CASM/etc/module.CGP
module load cgpVAFcommand
3 (choose caveman_c)
createVafCmd.pl -pid 3434 -o pileup_20241025 -g /lustre/scratch124/casm/team78pipelines/reference/human/GRCh38_full_analysis_set_plus_decoy_hla/genome.fa -hdr /lustre/scratch124/casm/team78pipelines/reference/human/GRCh38_full_analysis_set_plus_decoy_hla/HiDepth_v3.bed.gz -bo 1 -b twins_unm_subs.bed -mq 5 -bq 5

edit run_bsub.sh to increase memory requirements (because there are 22 samples and failed on standard)

#!/usr/bin/env bash
set -e;
set -o pipefail;
JID=`uuidgen`
bsub -oo cgpvaf%I.log -G team274-grp -q normal -J "$JID[1-96]%5" -n 1 -R 'select[mem>=4000] span[hosts=1] rusage[mem=4000]' -M2500 '/software/CASM/modules/installs/cgpVAFcommand/cgpVAFcommand-2.5.0/perl/bin/farm_idx_exec.pl cgpVafChr.cmd $LSB_JOBINDEX'
bsub -w "done($JID)||exit($JID)" -oo concat%I.log -q normal -G team274-grp -J 'batchjobs[1-4]%4' -n 1 -R 'select[mem>=1000] span[hosts=1] rusage[mem=1000]' -M1000 '/software/CASM/modules/installs/cgpVAFcommand/cgpVAFcommand-2.5.0/perl/bin/farm_idx_exec.pl cgpVafConcat.cmd $LSB_JOBINDEX'

bash run_bsub.sh

okay got this thing done!

cd /lustre/scratch126/casm/team274sb/bw18/twins_sarcoma/subs/pileup_20241025/output/PD62341v/snps
# get the tsv files onto your laptop 

rsync -avhu bw18@farm22:"/lustre/scratch126/casm/team274sb/bw18/twins_sarcoma/subs/pileup_20241025/output/PD62341v/snp/*.tsv" /Users/bw18/Desktop/1SB/Data/
rsync -avhu bw18@farm22:"/lustre/scratch126/casm/team274sb/bw18/twins_sarcoma/subs/pileup_20241025/output/PD63383w/snp/*.tsv" /Users/bw18/Desktop/1SB/Data/

# on your local: cd Desktop/1SB/Data/pileup_20241025
bash ../../Scripts/R1/20241011_removeheader.sh 

#########################################################################################
#########################################################################################
2024-10-28
# NOTE THAT I HAVE NOW DIFFERENT LISTS OF MUTATIONS
# 20241024 is 2,557 mutations 
# 20241028 is 742 mutations (added shearwater and low / high quality pileup)

#########################################################################################
2024-11-01
# update on scripts

current list of mutations: 413 (saved 30/10/2024)

# THESE ARE THE ONLY SCRIPTS THAT I AM CURRENTLY WORKING WITH

script to obtain this list of mutations: 20241011_pileup_filters

script to get QC plots for each sample: 20241018_pileup_qcplots

contamination script: 20241101_contamination

early developmental mutations script: 20241029_pileup_mutClasses_twins 

tumour evolution script: 20241026_pileup_mutClasses_tumour

#########################################################################################
#########################################################################################
2024-11-12

# get all indels identified in the whole genome sequencing data and write to a file to download 

ssh to the farm (ssh farm22)

# from the indel file for each sample, extract all introns:
# 1 remove the header
# 2 only allow indels which passed quality filters
# 3 require an annotated gene (require annotated protein)
# 4 remove indels in upstream, downstream and intronic variants 

ls /nfs/cancer_ref01/nst_links/live/3434/PD6*/

zcat /nfs/cancer_ref01/nst_links/live/3434/PD6*/*pindel.annot.vcf.gz | grep -v '#' | grep 'PASS' | grep 'VW=' | grep -v 'intron' | grep -v 'stream' | grep -v 'p.?' | wc -l

#########################################################################################
#########################################################################################
2024-11-13

# Telomere length 

scripts: 
bsub < scripts/20241113_telomere_length_bsub.sh

get scripts onto local 
rsync -avhu bw18@farm22:"scripts/*.sh" /Users/bw18/Desktop/1SB/Scripts/R1

get output onto local 
rsync -avhu bw18@farm22:"tel/*.csv" /Users/bw18/Desktop/1SB/Data

#########################################################################################
2024-11-14
# SUMMARY OF SCRIPTS 

1 Run pileup 

2 Filtering mutation set
20241011_pileup_filters.R > here you get 1,134 mutations which pass all quality filters
20241113_germlineCopyNumber.R > here, you further exclude germline mutations that are present in germline large structural variants which alter the copy number at a given region 
- don't exclude stuff that is not on a cluster that has wonky coverage and is good quality but not enriched in either twin 

3 Contamination estimation

4 Twin-twin transfusion estimation 

5 Phylogeny - assign all mutations to classes 

Misc
1 if you want to see the diff b/n 1002 and 1134 mutations, I made a quick script to look at what changed (20241114_compare_mutSets_1102_1134.R)

#########################################################################################
2024-11-19 

# get brass files (bedpe.gz)
rsync -avhu bw18@farm22:"/nfs/cancer_ref01/nst_links/live/3434/PD6*/*.brass.annot.bedpe.gz" /Users/bw18/Desktop/1SB/Data/

bash ../scripts/R1/20241119_removeheader_brass.sh

# run purple 

cd /lustre/scratch126/casm/team274sb/bw18/twins_sarcoma/scripts 
bsub COBALT_PD62341_PD63383.bsub
bsub AMBER_PD62341_PD63383.bsub

#########################################################################################
2024-11-20

# trying to figure out if I have PCR duplicates in my samples or not for chr1_388 mutation 

# How samtools finds duplicates in bam files 
# Duplicates are found by using the alignment data for each read (and its mate for paired reads). 
# Position and orientation (which strand it aligns against and in what direction) are used to to compare
# the reads against one another. If two (or more) reads have the same values 
# then the one with the highest base qualities is held to be the original and the others are the duplicates.

# Based on this, my understanding is that read pairs are considered duplicates if both mates are identical 
# You PCR the whole thing before you do the sequencing on the flowcell 

# want to look at bam files to see what was marked as duplicate 
cd /nfs/cancer_ref01/nst_links/live/3434/PD62341h 
samtools view -h -S -b -@ 6 -o

rsync -avhu bw18@farm22:"/lustre/scratch126/casm/team274sb/bw18/twins_sarcoma/purple/purple_output/*/without_GRIDDS/circos" /Users/bw18/Desktop/1SB
rsync -avhu bw18@farm22:"/lustre/scratch126/casm/team274sb/bw18/twins_sarcoma/purple/purple_output/*/without_GRIDDS/plot" /Users/bw18/Desktop/1SB
rsync -avhu bw18@farm22:"/lustre/scratch126/casm/team274sb/bw18/twins_sarcoma/purple/amber_output/*/*.baf.tsv.gz" /Users/bw18/Desktop/1SB

rsync -avhu bw18@farm22:"/nfs/cancer_ref01/nst_links/live/3434/*/*.dupmarked.bam.bas" /Users/bw18/Desktop/1SB

# get all .bam.bas files and write to /nfs
cat /nfs/cancer_ref01/nst_links/live/3434/*/*.sample.dupmarked.bam.bas > /nfs/users/nfs_b/bw18/cat_PD62341_PD63383.dupmarked.bam.bas

rsync -avhu bw18@farm22:"/lustre/scratch126/casm/team274sb/bw18/twins_sarcoma/duplicates/*" /Users/bw18/Desktop/1SB
# wrote a script to split this into files specific to each sample
# cd Desktop/1SB
bash Scripts/R1/split_chr1file.sh

#########################################################################################
# 2024-11-26
# I want to get the data from purple output 

rsync -avhu bw18@farm22:"/lustre/scratch126/casm/team274sb/bw18/twins_sarcoma/purple/cobalt_output/*/*.cobalt.ratio.tsv.gz" /Users/bw18/Desktop/1SB
rsync -avhu bw18@farm22:"/lustre/scratch126/casm/team274sb/bw18/twins_sarcoma/purple/purple_output/*/without_GRIDDS/*.purple.segment.tsv" /Users/bw18/Desktop/1SB
rsync -avhu bw18@farm22:"/lustre/scratch126/casm/team274sb/bw18/twins_sarcoma/purple/purple_output/*/without_GRIDDS/*.purple.purity.tsv" /Users/bw18/Desktop/1SB
rsync -avhu bw18@farm22:"/lustre/scratch126/casm/team274sb/bw18/twins_sarcoma/purple/purple_output/*/without_GRIDDS/*.purple.cnv.somatic.tsv" /Users/bw18/Desktop/1SB

#########################################################################################
#########################################################################################
# 2024-11-27
# I want to get the data from scRNA-seq data output 

# Data from Nathan 

# GOSH100 
CG_SB_NB13760628, CG_SB_NB13760629, CG_SB_NB13760630, CG_SB_NB14599068, CG_SB_NB14599069, CG_SB_NB14599070

# Single cell samples (fresh tissue)
CG_SB_NB13760628	GOSH-2023-100-C_1_G
CG_SB_NB13760629	GOSH-2023-100-C_2_G
CG_SB_NB13760630	GOSH-2023-100-C_3_G

# Nuclear samples (frozen tissue)
CG_SB_NB14599068	GOSH-2023-100-C-1a
CG_SB_NB14599069	GOSH-2023-100-C-1b
CG_SB_NB14599070	GOSH-2023-100-C-1c

# GOSH101 
CG_SB_NB13652544, CG_SB_NB13652545, CG_SB_NB13652546

# Single cell samples (fresh tissue)
CG_SB_NB13652544	GOSH-2023-101-C_1
CG_SB_NB13652545	GOSH-2023-101-C_2
CG_SB_NB13652546	GOSH-2023-101-C_3

# so to know what's where and to manage this data, I first made directories in Desktop/1SB/scRNAseq/
# each sample goes into a separate directory 
cd Desktop/1SB/scRNAseq/Data
mkdir -p CG_SB_NB13760628
mkdir -p CG_SB_NB13760629
mkdir -p CG_SB_NB13760630
mkdir -p CG_SB_NB14599068
mkdir -p CG_SB_NB14599069
mkdir -p CG_SB_NB14599070
mkdir -p CG_SB_NB13652544
mkdir -p CG_SB_NB13652545
mkdir -p CG_SB_NB13652546

# path to the data on the farm: /lustre/scratch126/casm/team274sb/project_folders/Sarcoma/sc_raw_data/nonRMS/
rsync -avhu bw18@farm22:"/lustre/scratch126/casm/team274sb/project_folders/Sarcoma/sc_raw_data/nonRMS/cellranger800_count_47451_CG_SB_NB13760628_GRCh38-1_2_0/" /Users/bw18/Desktop/1SB/scRNAseq/Data/CG_SB_NB13760628
rsync -avhu bw18@farm22:"/lustre/scratch126/casm/team274sb/project_folders/Sarcoma/sc_raw_data/nonRMS/cellranger800_count_47451_CG_SB_NB13760629_GRCh38-1_2_0/" /Users/bw18/Desktop/1SB/scRNAseq/Data/CG_SB_NB13760629
rsync -avhu bw18@farm22:"/lustre/scratch126/casm/team274sb/project_folders/Sarcoma/sc_raw_data/nonRMS/cellranger800_count_47451_CG_SB_NB13760630_GRCh38-1_2_0/" /Users/bw18/Desktop/1SB/scRNAseq/Data/CG_SB_NB13760630

rsync -avhu bw18@farm22:"/lustre/scratch126/casm/team274sb/project_folders/Sarcoma/sc_raw_data/nonRMS/cellranger800_count_48393_CG_SB_NB14599068_GRCh38-1_2_0/" /Users/bw18/Desktop/1SB/scRNAseq/Data/CG_SB_NB14599068
rsync -avhu bw18@farm22:"/lustre/scratch126/casm/team274sb/project_folders/Sarcoma/sc_raw_data/nonRMS/cellranger800_count_48393_CG_SB_NB14599069_GRCh38-1_2_0/" /Users/bw18/Desktop/1SB/scRNAseq/Data/CG_SB_NB14599069
rsync -avhu bw18@farm22:"/lustre/scratch126/casm/team274sb/project_folders/Sarcoma/sc_raw_data/nonRMS/cellranger800_count_48393_CG_SB_NB14599070_GRCh38-1_2_0/" /Users/bw18/Desktop/1SB/scRNAseq/Data/CG_SB_NB14599070

rsync -avhu bw18@farm22:"/lustre/scratch126/casm/team274sb/project_folders/Sarcoma/sc_raw_data/nonRMS/cellranger800_count_47089_CG_SB_NB13652544_GRCh38-1_2_0/" /Users/bw18/Desktop/1SB/scRNAseq/Data/CG_SB_NB13652544
rsync -avhu bw18@farm22:"/lustre/scratch126/casm/team274sb/project_folders/Sarcoma/sc_raw_data/nonRMS/cellranger800_count_47089_CG_SB_NB13652545_GRCh38-1_2_0/" /Users/bw18/Desktop/1SB/scRNAseq/Data/CG_SB_NB13652545
rsync -avhu bw18@farm22:"/lustre/scratch126/casm/team274sb/project_folders/Sarcoma/sc_raw_data/nonRMS/cellranger800_count_47089_CG_SB_NB13652546_GRCh38-1_2_0/" /Users/bw18/Desktop/1SB/scRNAseq/Data/CG_SB_NB13652546

# for each of those samples, there is a folder called filtered_feature_bc_matrix
# this contains 3 files (gzipped, can rsync and then get readable with gunzip *.gz from the directory)
# read below if you need a reminder of what those data files contain 

#########################################################################################
# VERY BASIC IDEA OF WHAT THE DATA IS HERE:
# I am using this GitHub to teach myself this: sounds like it could be very helpful
# https://crazyhottommy.github.io/scRNA-seq-workshop-Fall-2019/scRNAseq_workshop_0.html 

You will find the files Seurat works with in the filtered_feature_bc_matrix folder. 
There are 3 files in the folder:
# barcodes.tsv - this contains barcodes used in scRNA-seq (some ATCG sequences so I presume this is more for QC and I don't need to worry about this too much)
# features.tsv - this is a list of GENE IDs, genes and parameter you are looking at (Gene expression)
# matrix.mtx - this looks like it could be a count file??

# Basic checks 

# how many cells (barcodes)?
gzcat barcodes.tsv.gz | wc -l

# The `features.tsv.gz` contains the ENSEMBLE id and gene symbol
gzcat features.tsv.gz | head -5

# How many genes?
gzcat features.tsv.gz | wc -l # 33694

# matrix.mtx.gz is a sparse matrix which contains the non-zero counts
gzcat matrix.mtx.gz | head -10

# NOTE: Most of the entries in the final gene x cell count matrix are zeros. 
# Sparse matrix efficiently save the disk space by only recording the non-zero entries.

#########################################################################################
# Why do I have multiple files per sample? > this is count matrices for different lanes

# I am reading this here: https://bioinformatics-core-shared-training.github.io/UnivCambridge_ScRnaSeq_Nov2021/Markdowns/03_CellRanger.html 
# and it sounds like CellRanger can be run such that fastq from across the sample (ie for all lanes) is merged together
# therefore, I don't really get why this has not been done
# I thought of getting reports from sc_raw_data/cellranger_reports but there are no reports corresponding to my samples

#########################################################################################
# 2024-11-28

# I am trying to get the PURPLE output for those samples
# I ran PURPLE in tumour only mode but this did not work very well (I think)

# I discussed with Angus (slack) and he suggested running PURPLE in germline-only mode

# PURPLE command for germline_only mode 
from https://github.com/hartwigmedical/hmftools/blob/master/purple/README.md#germline-only

java -jar purple.jar \
   -reference COLO829R \
   -amber /path/COLO829/amber \
   -cobalt /path/COLO829/cobalt \
   -gc_profile /path/GC_profile.1000bp.37.cnp \
   -ref_genome_version 37 \
   -ref_genome /path/Homo_sapiens_assembly37.fasta \
   -ensembl_data_dir /path_to_ensembl_data_cache/ \
   -germline_hotspots /path/KnownHotspots.germline.37.vcf.gz \
   -germline_del_freq_file /path/cohort_germline_del_freq.csv \
   -germline_vcf /path/COLO829/COLO829.germline.vcf.gz \
   -run_drivers \
   -driver_gene_panel /path/DriverGenePanel.37.tsv \ 
   -output_dir /output/purple/ \

# I need the following files (extra to what I was using for tumour-only mode)
-germline_hotspots /path/KnownHotspots.germline.37.vcf.gz \
-germline_del_freq_file /path/cohort_germline_del_freq.csv \
-germline_vcf /path/COLO829/COLO829.germline.vcf.gz \ 
# this is going to possibly be a massive issue, since I do not know where the germline vcf is and how it relates to the files available in nfs/links/3434 for each sample 
# will try to run purple w/o this first and see what happens 

# I cannot find what the 'run drivers' command mean but I guess if it searches for some kind of genes I do not 
# care about this for the purpose of this analysis 
# I downloaded germline_hotspots and germline_del_freq_file from 
https://console.cloud.google.com/storage/browser/_details/hmf-public/HMFtools-Resources/dna_pipeline/v6_0/38/hmf_dna_pipeline_resources.38_v6.0.tar.gz;tab=live_object?inv=1&invt=AbiuIg 

hmf_public > HMFtools_resources > dna_pipeline > v6_0 > 38 > download the folder 
# I want to move germline hotspots and germline del freq file to the 

rsync -avhu ~/Desktop/1SB/PURPLE/purple_data/ bw18@farm22:/nfs/users/nfs_b/bw18/data

#########################################################################################
# 2024-12-01

# I am doing scRNA-seq now and trying to figure out how to do cell type annotation

# What kinds of cells can I expect in this dataset?

look at the following literature:
https://www.nature.com/articles/s43018-024-00743-y sarcoma microenvironment
# report the following cell types:
mast, fibroblast-like, epithelial-like, monocytes / macrophages, CD8+/CD4+ T cells, tumour cells, NK cells, B cell

https://www.nature.com/articles/s41388-024-03001-8?fromPaywallRec=false (undifferentiated pleomorphic sarcoma)
# monocyte, NK, T cell, mast cell, fibroblast, tumour cell, endothelial cell, pericyte
# what is a pericyte? cell that lines the wall of small blood vessels (capillaries)
# what is an endothelial cell? cell that lines all blood vessels (single layer)

# what else can we expect?
# any blood cells (basically, anything from the myeloid lineage)
# any cells that sit in blood vessels (because tumours do angiogenesis)
# tumour cells ofc 
# fibroblasts 
# here, can we also expect some nerve cells or microglia?

#########################################################################################
# 2024-12-04

# I now want to do single cell RNA-seq genotyping

# what I want to do is to look at all positions I am interested in
# I think I can start with the 255 mutations that I think are real and useful 
# I want to ask if they are present or absent in the cells that I have (which is ~7k per sample)

# I am creating a new script based on Nathan's scripts 
path to Nathan's scripts here:
/lustre/scratch126/casm/team274sb/na15/RMS/alleleintegrator/Oct10_2024/scripts/

Script based on his script is now on farm (/lustre/scratch126/bw18/twins_sarcoma/scripts)

cd Desktop/1SB/Data
rsync -avhu 20241204_mutations255.bed bw18@farm22:/lustre/scratch126/casm/team274sb/bw18/twins_sarcoma/scRNA_genotyping

# Excellent! AlleleIntegrator completed so I can now get the files from output and analyse
rsync -avhu bw18@farm22:"/lustre/scratch126/casm/team274sb/bw18/twins_sarcoma/scRNA_genotyping/output/" /Users/bw18/Desktop/1SB/scRNAseq/Data/

#########################################################################################
# 2024-12-10
# Mounting IGV on the farm (email with instructions from Nathan)

# NOTE: the guildelines in the email say REMOTE=casmins, this port has been changed at the beginning of December

vim ~/network_home
export MOUNTBASE=$HOME/volumes
export MOUNTNAME=bw18_network
export REMOTE=farm22
export REMOTELOC=/nfs/cancer_ref01/nst_links/live/3434/
export MUSER=bw18
export REMOTEPORT=farm22

vim ~/sshfs_config/cpg_mount.sh
#!/bin/bash
if [ $# -ne 1 ]; then
  (>&2 echo "ERROR: Please provide a fuse mount config file")
  exit 1
fi
if [ ! -e $1 ]; then
  (>&2 echo "ERROR: File does not exist: $1")
  exit 1
fi
set -e
source $1
  
PORT=22
if [ "x$REMOTEPORT" != "x" ]; then
  PORT=$REMOTEPORT
fi
 
set -u
FULLMOUNT=$MOUNTBASE/$MOUNTNAME
  
set +e
mkdir -p $FULLMOUNT
umount -f $FULLMOUNT >& /dev/null
set -e
sshfs -p $PORT -oauto_cache,cache=no,reconnect,ServerAliveInterval=15,ServerAliveCountMax=3,defer_permissions,follow_symlinks -o volname=$MOUNTNAME $MUSER@$REMOTE:$REMOTELOC $FULLMOUNT
echo "Mounted at: $FULLMOUNT"

alias mountIGV='sh ~/sshfs_config/cpg_mount.sh ~/sshfs_config/network_home'
alias umount='umount /Users/bw18/volumes/bw18_network'

# to mount a specific folder, you need to change the path to the thing on the farm and then 
# from your local do source ~/.bashrc > mountIGV > go to volumes > bw18_network > you can look at IGV now 

# 2024-12-11
# I want to look at the scRNA-seq bam files so mount IGV pointing to a different directory 

# bam files are here 
/lustre/scratch126/casm/team274sb/project_folders/Sarcoma/sc_bams/

export MOUNTBASE=$HOME/volumes
export MOUNTNAME=bw18_network
export REMOTE=farm22
export REMOTELOC=/lustre/scratch126/casm/team274sb/project_folders/Sarcoma/sc_bams/
export MUSER=bw18
export REMOTEPORT=farm22

#########################################################################################
# 2024-12-12
I am getting my scripts together so need to get the allele integrator script and outputs

# I decided to get all scripts from this rotation so I have record of everything 
rsync -avhu bw18@farm22:"/lustre/scratch126/casm/team274sb/bw18/twins_sarcoma/scripts/" /Users/bw18/Desktop/1SB/Scripts/R1/

# Get .err and .out from allele integrator 
664854 this was the jobID that worked and we used output from this 
(I have the error file as well but it's literally saying which packages were loaded in R)
rsync -avhu bw18@farm22:"/lustre/scratch126/casm/team274sb/bw18/twins_sarcoma/logs/*_664854.out" /Users/bw18/Desktop/1SB/Out/F7/

rsync -avhu bw18@farm22:"/nfs/users/nfs_b/bw18/3434_indels.txt" /Users/bw18/Desktop/1SB/Data/3434_indels_20241112.txt

#########################################################################################
# 2024-12-18
# Allele Integrator seems to be coming up with barcodes that don't exist in the Cell Ranger output so investigating what is going on 

# check barcodes in 
cd /lustre/scratch126/casm/team274sb/project_folders/Sarcoma/sc_raw_data/RMS/cellranger800_count_30855_CG_SB_NB8113359_GRCh38-1_2_0/filtered_feature_bc_matrix

rsync -avhu bw18@farm22:"/lustre/scratch126/casm/team274sb/project_folders/Sarcoma/sc_raw_data/RMS/cellranger800_count_30855_CG_SB_NB8113359_GRCh38-1_2_0/filtered_feature_bc_matrix/barcodes.tsv.gz" /Users/bw18/Desktop/1SB/Data/scRNAseq

rsync -avhu bw18@farm22:"/lustre/scratch126/casm/team274sb/project_folders/Sarcoma/sc_raw_data/nonRMS/cellranger800_count_47089_CG_SB_NB13652544_GRCh38-1_2_0/filtered_feature_bc_matrix/barcodes.tsv.gz" /Users/bw18/Desktop/1SB/Data/scRNAseq/barcodes_NB13652544.tsv.gz
rsync -avhu bw18@farm22:"/lustre/scratch126/casm/team274sb/project_folders/Sarcoma/sc_raw_data/nonRMS/cellranger800_count_47089_CG_SB_NB13652545_GRCh38-1_2_0/filtered_feature_bc_matrix/barcodes.tsv.gz" /Users/bw18/Desktop/1SB/Data/scRNAseq/barcodes_NB13652545.tsv.gz
rsync -avhu bw18@farm22:"/lustre/scratch126/casm/team274sb/project_folders/Sarcoma/sc_raw_data/nonRMS/cellranger800_count_47089_CG_SB_NB13652546_GRCh38-1_2_0/filtered_feature_bc_matrix/barcodes.tsv.gz" /Users/bw18/Desktop/1SB/Data/scRNAseq/barcodes_NB13652546.tsv.gz

rsync -avhu bw18@farm22:"/lustre/scratch126/casm/team274sb/project_folders/Sarcoma/sc_raw_data/nonRMS/cellranger800_count_47451_CG_SB_NB13760628_GRCh38-1_2_0/filtered_feature_bc_matrix/barcodes.tsv.gz" /Users/bw18/Desktop/1SB/Data/scRNAseq/barcodes_NB13760628.tsv.gz
rsync -avhu bw18@farm22:"/lustre/scratch126/casm/team274sb/project_folders/Sarcoma/sc_raw_data/nonRMS/cellranger800_count_47451_CG_SB_NB13760629_GRCh38-1_2_0/filtered_feature_bc_matrix/barcodes.tsv.gz" /Users/bw18/Desktop/1SB/Data/scRNAseq/barcodes_NB13760629.tsv.gz
rsync -avhu bw18@farm22:"/lustre/scratch126/casm/team274sb/project_folders/Sarcoma/sc_raw_data/nonRMS/cellranger800_count_47451_CG_SB_NB13760630_GRCh38-1_2_0/filtered_feature_bc_matrix/barcodes.tsv.gz" /Users/bw18/Desktop/1SB/Data/scRNAseq/barcodes_NB13760630.tsv.gz
