# Script to analyse the pileup (run 25/10/2024)
# 2024-10-25 
# Barbara Walkowiak bw18

# INPUT: pileup dataframes with mutation calls (from CaVEMan) for tumour and normal samples
# I ran the low quality pileup (-mq 5 -bq 5, see README 25/10/2024)
# obtained 4 separate .tsv files which I downloaded from the farm via rsync to /Desktop/1SB/Data/pileup_20241025

# In this script, I merge the file and prepare it for analysis (done in separate script)

# Load needed libraries
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(tidyr)
library(glue)

# COLUMNS (from the header)
##FORMAT-:DEP	Total reads covering this position (for subs del positions should be ignored)
##FORMAT-:FAZ	Reads presenting A for this position, forward strand
##FORMAT-:FCZ	Reads presenting C for this position, forward strand
##FORMAT-:FGZ	Reads presenting G for this position, forward strand
##FORMAT-:FTZ	Reads presenting T for this position, forward strand
##FORMAT-:MDR	Variant allele read directions 0=no reads; 1=Forward; 2=Reverse; 3=Forward + Reverse
##FORMAT-:MTR	Reads reporting the variant allele
##FORMAT-:OFS	Original filter status as defined in input vcf FILTER field
##FORMAT-:RAZ	Reads presenting A for this position, reverse strand
##FORMAT-:RCZ	Reads presenting C for this position, reverse strand
##FORMAT-:RGZ	Reads presenting G for this position, reverse strand
##FORMAT-:RTZ	Reads presenting T for this position, reverse strand
##FORMAT-:VAF	Variant Allele Fraction (excludes ambiguous reads if any)
##FORMAT-:WDR	Reference allele read directions 0=no reads; 1=Forward; 2=Reverse; 3=Forward + Reverse
##FORMAT-:WTR	Reads reporting the reference allele
##FILTER-:1	New filter status 1=not called in any sample
##FILTER-:2	New filter status 2=called in at least one sample
##FILTER-:3	New filter status 3=called + passed
##FILTER-:BD	Location from bed file
##ORIGINAL_FILTER-:CR	Position falls within a centromeric repeat using the supplied bed file
##ORIGINAL_FILTER-:DTH	Less than 1/3 mutant alleles were >= 25 base quality
##ORIGINAL_FILTER-:GI	Position falls within a germline indel using the supplied bed file
##ORIGINAL_FILTER-:HSD	Position falls within a high sequencing depth region using the supplied bed file
##ORIGINAL_FILTER-:MN	More than 0.05 of mutant alleles that were >= 15 base quality found in the matched normal
##ORIGINAL_FILTER-:MNP Tumour sample mutant allele proportion - normal sample mutant allele proportion < 0.2
##ORIGINAL_FILTER-:MQ	Mean mapping quality of the mutant allele reads was < 21
##ORIGINAL_FILTER-:PH	Mutant reads were on one strand (permitted proportion on other strand: 0.04), and mean mutant base quality was less than 21
##ORIGINAL_FILTER-:PT	Mutant alleles all on one direction of read (1rd allowed on opposite strand) and in second half of the read. Second half of read contains the motif GGC[AT]G in sequenced orientation and the mean base quality of all bases after the motif was less than 20
##ORIGINAL_FILTER-:RP	Coverage was less than 8 and no mutant alleles were found in the first 2/3 of a read (shifted 0.08 from the start and extended 0.08 more than 2/3 of the read length)
##ORIGINAL_FILTER-:SE	Coverage is >= 10 on each strand but mutant allele is only present on one strand
##ORIGINAL_FILTER-:SR	Position falls within a simple repeat using the supplied bed file
##ORIGINAL_FILTER-:VUM	Position has >= 3 mutant allele present in at least 1 percent unmatched normal samples in the unmatched VCF

###################################################################################################################################
###################################################################################################################################

setwd('/Users/bw18/Desktop/1SB/Data')
twins1 = data.table(read.csv('pileup_20241025/PD62341v_PD62341b_snp_vaf_removed_header.tsv', sep = '\t'))
twins2 = data.table(read.csv('pileup_20241025/PD63383w_PD63383t_snp_vaf_removed_header.tsv', sep = '\t'))
twins3 = data.table(read.csv('pileup_20241025/PD63383w_PD62341v_snp_vaf_removed_header.tsv', sep = '\t'))
twins4 = data.table(read.csv('pileup_20241025/PD62341v_PD63383w_snp_vaf_removed_header.tsv', sep = '\t'))

# Check the number of mutations in each of the pileups:
paste('Number of mutations in pileup PD62341v_PD62341b:', dim(twins1)[1]) # 362814
paste('Number of mutations in pileup PD63383w_PD63383t:', dim(twins2)[1]) # 362814
paste('Number of mutations in pileup PD63383w_PD62341v:', dim(twins3)[1]) # 362814
paste('Number of mutations in pileup PD62341v_PD63383w:', dim(twins4)[1]) # 362814

# For each dataframe, check which columns contain NA values
names(which(colSums(is.na(twins1)) > 0)) # only OFS columns 
names(which(colSums(is.na(twins2)) > 0)) # only OFS columns
names(which(colSums(is.na(twins3)) > 0)) # only OFS columns
names(which(colSums(is.na(twins4)) > 0)) # only OFS columns 

# Create a column to specify the mutation (position and DNA change)
twins1[, mut_ID := paste(Chrom, Pos, Ref, Alt, sep='_')]
twins2[, mut_ID := paste(Chrom, Pos, Ref, Alt, sep='_')]
twins3[, mut_ID := paste(Chrom, Pos, Ref, Alt, sep='_')]
twins4[, mut_ID := paste(Chrom, Pos, Ref, Alt, sep='_')]

######################################################################################################
###################################################################################################################################

# CREATE MERGED DATAFRAME (I WILL REDO THIS ONCE THE ISSUE WITH THE 2nd FILE IS RESOLVED)

# Info columns
cols_info = c('VariantID', 'Chrom', 'Pos', 'Ref', 'Alt', 'Qual', 'Filter', 
              'Gene', 'Transcript', 'RNA', 'CDS', 'Protein', 'Type', 'Effect')

# Table with information
df_info = twins1[, c(cols_info, 'mut_ID'), with=FALSE]

# Get sub dfs with only data and mutation ID 
twins1_sub = twins1[, c('X.Normal', cols_info) := NULL]
twins2_sub = twins2[, c('X.Normal', cols_info) := NULL]
twins3_sub = twins3[, c('X.Normal', cols_info) := NULL]
twins4_sub = twins4[, c('X.Normal', cols_info) := NULL]
# the column names are the same in PD62341 and PD63383 so I gather you can use those once 

# create a merged dataframe 
twins_m1 = merge(twins1_sub, twins2_sub)
twins_df = merge(twins_m1, df_info) # add the data 

# save the merged dataframe to a file
write.csv(twins_df, 'pileup_merged_20241025.tsv', row.names=FALSE)

###################################################################################################################################
###################################################################################################################################
