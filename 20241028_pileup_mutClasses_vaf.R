# Script to analyse the pileup (high quality, run 15/10/2024)
# Tumour specific mutations (classification by VAF)
# 2024-10-28
# Barbara Walkowiak bw18

# INPUT: 
# 1 merged pileup dataframes with mutation calls (from CaVEMan) for tumour and normal samples
# 2 list of mutations that passed required filters 

# Cleaned up script to identify mutations of interest 
# Script to look at tumour-specific mutations and tumour evolution 

###################################################################################################################################
# Load needed libraries
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(tidyr)
library(glue)
library(comprehenr) # list comprehension in R 
library(stringr)
library(beeswarm)
library(viridis)
library(grid)
library(pheatmap)

###################################################################################################################################
# INPUT FILES 
# Read the merged dataframe 
setwd('/Users/bw18/Desktop/1SB')
twins_dt = data.table(read.csv('Data/pileup_merged_20241016.tsv')) # import high quality pileup

# Drop columns with PD38is_wgs (used as reference)
twins_PDv38is = grep("PDv38is", names(twins_dt), value = TRUE)
twins_dt[, c(twins_PDv38is) := NULL]

# Filter to only include mutations retained for filtering 
muts = read.table('Data/mutations_include_20241028_3.txt') %>% unlist()
paste('Number of mutations that passed required filters:', length(muts)) # 255
twins_filtered_dt = twins_dt[mut_ID %in% muts]

# Add column to indicate chromosomes lost in the tumour
twins_filtered_dt[, loss := as.factor(fcase( 
  Chrom %in% c('chr1', 'chr18'), 'loss in tumour', # chr1 and chr18 segments lost in tumour samples
  !Chrom %in% c('chr1', 'chr18'), 'normal ploidy'
))]

###################################################################################################################################
# PLOT SETTINGS
# Specify colors for plotting 
col_tumour = '#ad0505'
col_normal = '#07a7d0'
col_PD62341 = "#a45f03"
col_PD63383 = "#d4b17e"
col_tumour_PD62341 = "#980505"
col_tumour_PD63383 = "#eb6767"
col_normal_PD62341 = "#0785a5"
col_normal_PD63383 = "#70c6db"
col_bar = '#8c08d3'

######################################################################################################
# SAMPLES
# create lists of possible samples of interest
samples_names = c("PD62341v", "PD62341q", "PD62341aa", "PD62341ad", "PD63383w", "PD63383t", "PD63383u", "PD63383ae", "PD63383ak", "PD63383bb", "PD62341ae", "PD62341ag", "PD62341aj", "PD62341ak", "PD62341am", "PD62341ap",  "PD63383ap", "PD63383aq", "PD62341b", "PD62341h", "PD62341n", "PD62341u")
samples_normal = c("PD62341v", "PD62341q", "PD62341aa", "PD62341ad", "PD63383w", "PD63383t", "PD63383u", "PD63383ae", "PD63383ak", "PD63383bb", "PD62341h", "PD62341n")
samples_tumour = c("PD62341ae", "PD62341ag", "PD62341aj", "PD62341ak", "PD62341am", "PD62341ap",  "PD63383ap", "PD63383aq", "PD62341b", "PD62341u")
samples_PD62341 = c("PD62341v", "PD62341q", "PD62341aa", "PD62341ad", "PD62341ae", "PD62341ag", "PD62341aj", "PD62341ak", "PD62341am", "PD62341ap", "PD62341b", "PD62341h", "PD62341n", "PD62341u")
samples_PD63383 = c("PD63383w", "PD63383t", "PD63383u", "PD63383ae", "PD63383ak", "PD63383ap", "PD63383aq", "PD63383bb")
samples_normal_PD62341 = c("PD62341v", "PD62341q", "PD62341aa", "PD62341ad", "PD62341h", "PD62341n")
samples_normal_PD63383 = c("PD63383w", "PD63383t", "PD63383u", "PD63383ae", "PD63383ak", "PD63383bb")
samples_tumour_PD62341 = c("PD62341ae", "PD62341ag", "PD62341aj", "PD62341ak", "PD62341am", "PD62341ap", "PD62341b", "PD62341u")
samples_tumour_PD63383 = c( "PD63383ap", "PD63383aq")

samples_mtr = grep("MTR", names(twins_dt), value = TRUE)
samples_dep = grep("DEP", names(twins_dt), value = TRUE)
samples_vaf = grep("VAF", names(twins_dt), value = TRUE)

samples_normal_mtr = paste(samples_normal, 'MTR', sep='_')
samples_normal_dep = paste(samples_normal, 'DEP', sep='_')
samples_normal_vaf = paste(samples_normal, 'VAF', sep='_')

samples_tumour_mtr = paste(samples_tumour, 'MTR', sep='_')
samples_tumour_dep = paste(samples_tumour, 'DEP', sep='_')
samples_tumour_vaf = paste(samples_tumour, 'VAF', sep='_')

samples_normal_PD62341_mtr = paste(samples_normal_PD62341, 'MTR', sep='_')
samples_normal_PD62341_dep = paste(samples_normal_PD62341, 'DEP', sep='_')
samples_normal_PD62341_vaf = paste(samples_normal_PD62341, 'VAF', sep='_')

samples_tumour_PD62341_mtr = paste(samples_tumour_PD62341, 'MTR', sep='_')
samples_tumour_PD62341_dep = paste(samples_tumour_PD62341, 'DEP', sep='_')
samples_tumour_PD62341_vaf = paste(samples_tumour_PD62341, 'VAF', sep='_')

samples_normal_PD63383_mtr = paste(samples_normal_PD63383, 'MTR', sep='_')
samples_normal_PD63383_dep = paste(samples_normal_PD63383, 'DEP', sep='_')
samples_normal_PD63383_vaf = paste(samples_normal_PD63383, 'VAF', sep='_')

samples_tumour_PD63383_mtr = paste(samples_tumour_PD63383, 'MTR', sep='_')
samples_tumour_PD63383_dep = paste(samples_tumour_PD63383, 'DEP', sep='_')
samples_tumour_PD63383_vaf = paste(samples_tumour_PD63383, 'VAF', sep='_')

samples_PD62341_mtr = paste(samples_PD62341, 'MTR', sep='_')
samples_PD62341_dep = paste(samples_PD62341, 'DEP', sep='_')
samples_PD62341_vaf = paste(samples_PD62341, 'VAF', sep='_')

samples_PD63383_mtr = paste(samples_PD63383, 'MTR', sep='_')
samples_PD63383_dep = paste(samples_PD63383, 'DEP', sep='_')
samples_PD63383_vaf = paste(samples_PD63383, 'VAF', sep='_')

######################################################################################################
# MTR, DEP, VAF DATA 

# subset MTR, DEP and VAF data
twins_filtered_mtr = twins_filtered_dt[,c('mut_ID', samples_mtr), with=FALSE]
twins_filtered_dep = twins_filtered_dt[,c('mut_ID', samples_dep), with=FALSE]
twins_filtered_vaf = twins_filtered_dt[,c('mut_ID', samples_vaf), with=FALSE]

# determine min, max and mean nr of variant reads / depth / vaf
twins_filtered_mtr[,mean_mtr := apply(.SD, 1, function(x) mean(x[x>0])), .SDcols = samples_mtr]
twins_filtered_mtr[,median_mtr := apply(.SD, 1, function(x) median(x[x>0])), .SDcols = samples_mtr]
twins_filtered_mtr[,min_mtr := apply(.SD, 1, function(x) min(x[x>0])), .SDcols = samples_mtr]
twins_filtered_mtr[,max_mtr := apply(.SD, 1, max), .SDcols = samples_mtr]

twins_filtered_dep[,mean_dep := apply(.SD, 1, function(x) mean(x[x>0])), .SDcols = samples_dep]
twins_filtered_dep[,median_dep := apply(.SD, 1, function(x) median(x[x>0])), .SDcols = samples_dep]
twins_filtered_dep[,min_dep := apply(.SD, 1, function(x) min(x[x>0])), .SDcols = samples_dep]
twins_filtered_dep[,max_dep := apply(.SD, 1, max), .SDcols = samples_dep]

twins_filtered_vaf[,mean_vaf := apply(.SD, 1, function(x) mean(x[x>0])), .SDcols = samples_vaf]
twins_filtered_vaf[,median_vaf := apply(.SD, 1, function(x) median(x[x>0])), .SDcols = samples_vaf]
twins_filtered_vaf[,min_vaf := apply(.SD, 1, function(x) min(x[x>0])), .SDcols = samples_vaf]
twins_filtered_vaf[,max_vaf := apply(.SD, 1, max), .SDcols = samples_vaf]

# determine this specifically in normal samples (in case of ploidy changes in the tumour)
twins_filtered_mtr[,mean_mtr_normal := apply(.SD, 1, function(x) mean(x[x>0])), .SDcols = samples_normal_mtr]
twins_filtered_mtr[,median_mtr_normal := apply(.SD, 1, function(x) median(x[x>0])), .SDcols = samples_normal_mtr]
twins_filtered_mtr[,min_mtr_normal := apply(.SD, 1, function(x) min(x[x>0])), .SDcols = samples_normal_mtr]
twins_filtered_mtr[,max_mtr_normal := apply(.SD, 1, max), .SDcols = samples_normal_mtr]

twins_filtered_dep[,mean_dep_normal := apply(.SD, 1, function(x) mean(x[x>0])), .SDcols = samples_normal_dep]
twins_filtered_dep[,median_dep_normal := apply(.SD, 1, function(x) median(x[x>0])), .SDcols = samples_normal_dep]
twins_filtered_dep[,min_dep_normal := apply(.SD, 1, function(x) min(x[x>0])), .SDcols = samples_normal_dep]
twins_filtered_dep[,max_dep_normal := apply(.SD, 1, max), .SDcols = samples_normal_dep]

twins_filtered_vaf[,mean_vaf_normal := apply(.SD, 1, function(x) mean(x[x>0])), .SDcols = samples_normal_vaf]
twins_filtered_vaf[,median_vaf_normal := apply(.SD, 1, function(x) median(x[x>0])), .SDcols = samples_normal_vaf]
twins_filtered_vaf[,min_vaf_normal := apply(.SD, 1, function(x) min(x[x>0])), .SDcols = samples_normal_vaf]
twins_filtered_vaf[,max_vaf_normal := apply(.SD, 1, max), .SDcols = samples_normal_vaf]

######################################################################################################
# PRESENCE / ABSENCE 

# Tumour-specific mutations 

# I will allow up to 3 normal samples due to the possibility of contamination
# I suspect that skin samples (aa, bb) and liver sample (PD62341h) are contaminated with reads from the tumour

# Mutations found in all tumour samples 
twins_filtered_vaf[sum_tumour == 10 & sum_normal == 0] # 0 

twins_filtered_vaf[sum_tumour == 10 & sum_normal == 1] # 6 
# "chr16_32621521_C_A" # mapping quality not excellent (double check)
# "chr17_78256883_C_T" # looks okay
# "chr22_17774668_G_A" # looks okay 
# "chr22_30553442_G_T" # looks okay
# chr13_60216504_A_C
# chr2_199565129_G_A  looks okay

twins_filtered_vaf[sum_tumour == 10 & sum_normal == 2] # 14
# chr10_100754461_C_T # looks okay
# chr11_134158143_C_T # looks okay 
# chr14_104122986_G_T
# "chr15_48353502_C_A" # looks okay
# "chr15_56691722_C_T" # looks okay
# "chr16_4003400_G_A" # looks okay
# "chr1_69984580_G_A" # can estimate contamination from this one!!
# "chr2_133311125_G_T" # looks okay
# chr3_49878199_C_A
# chr4_178731458_A_G
# "chr5_136349883_C_T" # looks okay
# "chr5_37790015_G_A" # looks okay
# "chr6_95827754_C_A" # mapping could be better but still decent
# "chr9_11364991_T_A" # looks okay

twins_filtered_vaf[sum_tumour == 10 & sum_normal == 3] # 6 
# it's all bb, h and aa (skin + liver)
# chr14_104500332_C_T # looks okay
# chr15_23462705_C_A # lots of mis-mapped reads 
# chr17_40061856_G_C # looks okay
# chr4_147908994_T_C # looks okay
# chrX_66719643_C_G # looks okay
# chrX_68803487_C_T # looks okay

twins_filtered_vaf[sum_tumour == 9 & sum_normal == 0] # 0 

twins_filtered_vaf[sum_tumour == 9 & sum_normal == 1] # 8
# chr12_61020950_C_A # looks okay
# chr19_3961910_G_A # looks okay
# chrX_46375010_C_T # looks okay
# chrX_56448796_C_A # looks okay

twins_filtered_vaf[sum_tumour == 9 & sum_normal == 2] # 19
# "chr12_10337296_G_T" # variable mapping (- ok verified this is poor quality (double check)
# "chr12_39248376_A_T" # looks okay
# chr13_100405282_G_T # looks okay
# chr1_119005653_G_C 
# chr1_60839899_G_A
# chr20_44114996_C_T
# "chr21_45729938_T_A" # looks okay
# "chr22_16694683_C_A" # looks okay
# "chr2_187893783_G_T" # looks okay
# chr3_153319439_A_G
# "chr3_94699090_G_T" # looks okay
# "chr4_180479634_G_A" # looks okay
# "chr5_112010886_G_A" # looks okay
# "chr7_12297469_T_G" # looks okay
# chr7_139659050_G_A
# chr8_54112197_T_C
# "chr8_54112201_C_A" # looks okay
# "chr9_94529993_C_T" # quite a lot of mis-mapped reads (double check)

# there are two mutations next to each other that show a very similar pattern and VAF 
# the other is: chr8_54112197_T_C, and is called as present in 9 tumour samples and 3 normal
# normal chr8_54112197_T_C: PD62341aa, h, PD63383 bb
# normal chr8_54112201_C_A: PD62341h, PD63383bb (3 reads in PD62341aa)

twins_filtered_vaf[sum_tumour == 9 & sum_normal == 3] # 15
# "chr15_76646556_T_C" # looks okay  
# "chr17_5029328_G_C" # looks okay 
# "chr2_124573116_C_T" # looks great  
# "chr2_1420635_A_T" # looks great # there is a beautiful mutation next to this one 
# chr4_180479634_G_A
# "chr4_46067117_A_G" # difficult region but looks okay
# "chr4_93081874_T_A" # looks great  
# "chr7_25564729_G_A" # looks okay
# "chr8_131571989_A_T" # looks okay
# chr8_46678612_T_C
# "chr8_97785584_T_C" # looks okay 
# "chr9_12824689_G_A" # looks okay 
# "chrX_124531709_C_T" # looks okay  
# "chrX_83934262_C_G" # looks okay

# very interesting mutation chr1_60839875_T_C 
# seems to be on the other chromosome 
# can compare reads with different muts to determine contamination in different tumour samples
# also, PD62341q for example has 5 of those reads that belong to the tumour chromosome 
# we don't know tho if the loss of chr1 is subclonal in some tumour samples - what about this?

twins_filtered_vaf[sum_tumour == 8 & sum_normal == 0] # empty  

twins_filtered_vaf[sum_tumour == 8 & sum_normal == 1] # 1
# "chr12_106354695_C_T" # mixed feelings - looks okay but quite a lot of insertions, double check

twins_filtered_vaf[sum_tumour == 8 & sum_normal == 2] # 2
# chr14_49354012_C_T
# chr2_78054577_G_A # looks okay
# chr2_57814739_A_C 

twins_filtered_vaf[sum_tumour == 7 & sum_normal == 0] # 0  

twins_filtered_vaf[sum_tumour == 7 & sum_normal == 1] # 2

twins_filtered_vaf[sum_tumour == 7 & sum_normal == 2] # 1

twins_filtered_vaf[sum_tumour == 7 & sum_normal == 3] # 1

# Not sure the below are worth it since we don't know if they tell us anything about the tumour phylogeny or not 
twins_filtered_vaf[sum_tumour == 6 & sum_normal == 0] # 0 

twins_filtered_vaf[sum_tumour == 6 & sum_normal == 1] # 4
# chr10_51734221_A_T # looks okay
# chr12_112652964_C_T # looks okay
# chr14_54164980_C_G 
# chr15_94939486_G_A # looks okay
# chr6_58221412_T_C # mutations on poorly mapped reads

twins_filtered_vaf[sum_tumour == 6 & sum_normal == 2] # 1
# chr2_158600954_G_A

twins_filtered_vaf[sum_tumour == 5 & sum_normal == 0] # 2
# chr16_17805832_C_T # looks okay 

twins_filtered_vaf[sum_tumour == 5 & sum_normal == 1] # 6
# "chr11_26060345_G_T" # looks okay
# chr12_112652964_C_T
# "chr13_103318073_T_A" # looks okay
# "chr17_70958446_C_T" # looks okay
# chr4_135461656_G_A
# "chr6_8310684_G_T" # looks okay

twins_filtered_vaf[sum_tumour == 4 & sum_normal == 0] # 1
# chr8_21444725_G_A

twins_filtered_vaf[sum_tumour == 4 & sum_normal == 1] # 3
# "chr12_114403258_G_A" # looks okay
# chr12_82832044_C_A
# "chr18_1379467_G_T"  # looks okay
# chr6_168731130_C_T

twins_filtered_vaf[sum_tumour == 3 & sum_normal == 0] # 4  
# chr13_28953051_G_T

twins_filtered_vaf[sum_tumour == 3 & sum_normal == 1] # 5
# "chr8_141860196_G_A" # looks okay

######################################################################################################
# Mutations specific to PD63383 tumour

twins_filtered_vaf[sum_normal == 0 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 0, mut_ID] # 10
# "chr10_100308403_C_T" # looks okay
# "chr11_80115755_C_T" # looks okay
# "chr19_43348626_C_A" # not sure, double check
# "chr1_51714096_G_C" # looks okay  
# "chr2_72939205_C_T" # looks okay 
# "chr5_157248612_A_G" # looks okay
# "chr7_13677323_A_G" # mapping not excellent, double check
# "chrX_72671786_A_T" # poor mapping + strand bias   

muts_PD63383_tumour1 = c("chr10_100308403_C_T", "chr11_80115755_C_T", "chr19_43348626_C_A", "chr1_51714096_G_C",
                         "chr2_72939205_C_T", "chr2_82147172_T_C", "chr5_157248612_A_G", "chr5_28014472_C_T",
                         "chr5_54017866_G_T", "chrX_48453603_G_A", "chrX_72671786_A_T", "chrX_77681162_G_A")
twins_info = twins_dt[, c('mut_ID', 'Gene', 'Transcript', 'RNA', 'CDS', 'Protein', 'Type', 'Effect'), with=FALSE]
twins_info[mut_ID %in% muts_PD63383_tumour1] # nothing in protein coding genes

twins_filtered_vaf[sum_normal_PD63383 == 1  & sum_normal_PD62341 == 0 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 0, mut_ID] # 15
# "chr10_31333522_G_A" # looks okay # bb
# "chr12_43132316_G_A" # looks okay but close to poly-T # bb
# "chr12_80762031_G_A" # looks okay # bb
# "chr1_160088990_G_A" # looks okay # bb
# "chr22_25281865_C_T" # looks okay # bb
# "chr3_137508691_C_A" # looks okay # bb
# chr4_179587218_G_A
# "chr6_86549064_C_T" # looks okay # bb
# "chr6_92179408_A_G" # looks okay # bb
# "chr7_148446413_G_A" # looks okay # bb
# "chr7_49732059_G_A" # looks okay # bb 
# "chr9_98275090_C_T" # looks okay # bb
# "chrX_117418608_G_A" # looks okay # bb
# "chrX_119915655_G_A" # looks okay # bb
# "chrX_36779694_C_T" # looks okay # bb
# "chrX_70352624_G_A" # looks okay # bb
# chrX_77681162_G_A # looks okay

# "chr7_49732059_G_A" note that there is a mutation next to this one 
# chr7_49732077_G_C which is present in more samples
# so can track subclonal evolution in a sense 
# good to check where chr7_49732077_G_C is present (germline?)
# I tried to check but not sure it was called anywhere 

# reads in normal samples are exclusively present in bb
# PD63383bb is skin, so I find it highly plausible that this is contamination with tumour cells
# especially that the kid had lesions in the skin - can I double check with Henry what he thinks?

muts_PD63383_tumour2 = c("chr10_31333522_G_A", "chr12_43132316_G_A", "chr12_80762031_G_A", "chr1_160088990_G_A",
                         "chr22_25281865_C_T", "chr3_137508691_C_A", "chr6_86549064_C_T", "chr6_92179408_A_G", 
                         "chr7_148446413_G_A", "chr7_49732059_G_A",  "chr9_98275090_C_T", "chrX_117418608_G_A",
                         "chrX_36779694_C_T", "chrX_70352624_G_A")
twins_info[mut_ID %in% muts_PD63383_tumour2] # nothing in protein coding genes

twins_filtered_vaf[sum_normal_PD63383 == 2  & sum_normal_PD62341 == 0 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 0, mut_ID] # 0
# empty

twins_filtered_vaf[sum_normal_PD63383 == 3  & sum_normal_PD62341 == 0 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 0, mut_ID] # 0
# empty 

twins_filtered_vaf[sum_normal_PD63383 == 0  & sum_normal_PD62341 == 1 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 0, mut_ID] # 0
# empty 

# how many mutations specific to PD63383 tumour have we identified?
muts_PD63383_tumour = c(muts_PD63383_tumour1, muts_PD63383_tumour2)
paste('Number of mutations private to PD63383 tumour:', length(muts_PD63383_tumour)) # 26 

######################################################################################################
# Mutations private to PD62341 tumour (subclonal)

# searching for mutations acquired after PD63383 tumour migration
# actually, good question: has the tumour migrated once only? or were there several transmissions?

# no muts present in 8/8 PD62341 and none PD63383 tumour
# actually not sure if it makes sense to split PD63383 / PD62341 - if there is nothing in PD63383 tumour, the normal samples couldn't have been contaminated
twins_filtered_vaf[sum_normal_PD62341 == 0 & sum_normal_PD63383 == 0 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 8] # 0
twins_filtered_vaf[sum_normal_PD62341 == 1 & sum_normal_PD63383 == 0 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 8] # 0
twins_filtered_vaf[sum_normal_PD62341 == 2 & sum_normal_PD63383 == 0 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 8] # 0

# no muts present in 7/8 PD62341 and none PD63383 tumour
twins_filtered_vaf[sum_normal == 0 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 7] # 0
twins_filtered_vaf[sum_normal == 1 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 7] # 0
twins_filtered_vaf[sum_normal == 2 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 7] # 0

# muts present in 6/8 PD62341 tumour and none PD63383 tumour (nothing real)
twins_filtered_vaf[sum_normal == 0 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 6] # 0
twins_filtered_vaf[sum_normal == 1 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 6, mut_ID] # 1
twins_filtered_vaf[sum_normal == 2 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 6, mut_ID] # 0

# muts present in 5/8 PD62341 and none PD63383 tumour (2 real ones!)
twins_filtered_vaf[sum_normal == 0 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 5, mut_ID] # 2
# "chr16_17805832_C_T" # looks okay 
# ae, ag, aj, b, u (interesting - all non-primary tumour so that would suggest the PD63383 emerged from primary)

twins_filtered_vaf[sum_normal == 0 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 5] # 0
# "chr16_17805832_C_T" # looks okay 

twins_filtered_vaf[sum_normal == 1 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 5] # 1
# "chr6_8310684_G_T" # looks okay!
# ae, ag, aj, b, u (again!) and present in PD62341h normal - likely contaminated

twins_filtered_vaf[sum_normal == 2 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 5, mut_ID] # 0
# empty

# muts present in 4/10 PD62341 and none PD63383 tumour
twins_filtered_vaf[sum_normal == 0 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 4, mut_ID] # 1

twins_filtered_vaf[sum_normal == 1 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 4, mut_ID] # 2
# "chr12_114403258_G_A" # looks okay
# ae, aj, b, u (2 reads in ag and aa, 1 in am and ap)
# "chr18_1379467_G_T" # looks okay 
# ae, ag, aj, u (3 reads in b, 2 in aa, 1 in v, am, also in h but likely contaminated)

twins_filtered_vaf[sum_normal == 2 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 4, mut_ID] # 0
# empty

# muts present in 3/10 PD62341 and none PD63383 tumour
twins_filtered_vaf[sum_normal == 0 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 3, mut_ID] # 3
# "chr13_28953051_G_T"

twins_filtered_vaf[sum_normal == 1 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 3, mut_ID] # 3

twins_filtered_vaf[sum_normal == 2 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 3, mut_ID] # 0

# muts present in 2/10 PD62341 and none PD63383 tumour
# not sure these are worth investigating in detail because ultimately each tumour sample would have been seeded by more than one cell probably
# so not very likely that you would identify anything that is truly clonal 
twins_filtered_vaf[sum_normal == 0 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 2, mut_ID] # 0
# "chr11_43029745_A_T" 
# "chr12_19555916_C_T" 
# "chr2_231591536_G_T"

twins_filtered_vaf[sum_normal == 1 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 2, mut_ID] # 1

twins_filtered_vaf[sum_normal == 2 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 2, mut_ID] # 2
# "chr9_39590690_G_A" # poor mapping 

######################################################################################################
# Mutations shared between tumour and one twin 

# Mutations present in ALL tumour samples and all / some normal PD62341 samples 
twins_filtered_vaf[sum_tumour == 10 & sum_normal_PD62341 == 6 & sum_normal_PD63383 == 0] # 1
# "chr14_105458006_C_A" our favourite one!

twins_filtered_vaf[sum_tumour == 10 & sum_normal_PD62341 == 5 & sum_normal_PD63383 == 0] # 0

twins_filtered_vaf[sum_tumour == 10 & sum_normal_PD62341 == 4 & sum_normal_PD63383 == 0] # 2
# "chr3_62055057_C_G" aa and h (possibly contamination), n and q (heart and pancreas), missing in v (?) and ad
# "chr3_62055077_G_C" aa and h (possibly contamination), n and q (heart and pancreas), missing in v (?) and ad

twins_filtered_vaf[sum_tumour == 10 & sum_normal_PD62341 == 3 & sum_normal_PD63383 == 0] # 0

twins_filtered_vaf[sum_tumour == 10 & sum_normal_PD62341 == 2 & sum_normal_PD63383 == 0] # 1
# chr6_95827754_C_A not the best mapping quality but fine, aa and h (possibly contamination) 
# chr10_100754461_C_T

twins_filtered_vaf[sum_tumour == 10 & sum_normal_PD62341 == 1 & sum_normal_PD63383 == 0] # 4
# "chr13_60216504_A_C"
# "chr16_32621521_C_A" # mapping quality not excellent (double check because it doesn't look that bad)
# "chr17_78256883_C_T" # looks okay
# "chr22_17774668_G_A" # looks okay 
# "chr22_30553442_G_T" # looks okay
# chr2_199565129_G_A

# Mutations present in 9 tumour samples and all / some PD62341 normal samples 
twins_filtered_vaf[sum_tumour == 9 & sum_normal_PD62341 == 6 & sum_normal_PD63383 == 0] # 0
twins_filtered_vaf[sum_tumour == 9 & sum_normal_PD62341 == 5 & sum_normal_PD63383 == 0] # 0
twins_filtered_vaf[sum_tumour == 9 & sum_normal_PD62341 == 4 & sum_normal_PD63383 == 0] # 0
twins_filtered_vaf[sum_tumour == 9 & sum_normal_PD62341 == 3 & sum_normal_PD63383 == 0] # 1
# chr8_46678612_T_C
# also there is a very cool mutation on the other chromosome (in all samples?) chr7_139659034_T_G
# missing in ak, which is technically tumour (primary)

twins_filtered_vaf[sum_tumour == 9 & sum_normal_PD62341 == 2 & sum_normal_PD63383 == 0] # 3
# chr20_44114996_C_T
# chr21_45729938_T_A # looks okay
# chr7_139659050_G_A
# all missing in ak 

twins_filtered_vaf[sum_tumour == 9 & sum_normal_PD62341 == 1 & sum_normal_PD63383 == 0] # 4
# mostly in h so all of these are likely contamination
# chr12_61020950_C_A # looks okay
# chr19_3961910_G_A # looks okay
# chrX_46375010_C_T # looks okay
# chrX_56448796_C_A # looks okay
# all missing in ak 

# Allow presence in 1 normal sample from PD63383 (because of contamination)
twins_filtered_vaf[sum_tumour == 10 & sum_normal_PD62341 == 6 & sum_normal_PD63383 == 1] # 1
# chr17_33422229_C_A

twins_filtered_vaf[sum_tumour == 10 & sum_normal_PD62341 == 5 & sum_normal_PD63383 == 1] # 3
# chr16_5479739_C_T # looks real (absent from PD62341v, present in PD63383bb)
# chr2_95662131_G_A # looks semi-okay (double check) (absent from PD62341v, present in PD63383bb)
# chr3_50106043_C_T # looks okay (absent from PD62341v, present in PD63383bb)

twins_filtered_vaf[sum_tumour == 10 & sum_normal_PD62341 == 4 & sum_normal_PD63383 == 1] # 0

twins_filtered_vaf[sum_tumour == 10 & sum_normal_PD62341 == 3 & sum_normal_PD63383 == 1] # 1
# chr5_44907911_G_A # looks okay (present in PD62341q, aa, h and PD63383bb)
# chrX_115495236_T_C

twins_filtered_vaf[sum_tumour == 10 & sum_normal_PD62341 == 2 & sum_normal_PD63383 == 1] # 9
# always PD62341aa, PD62341h and PD63383bb so very likely contamination 
# "chr14_104500332_C_T" # looks okay
# "chr15_23462705_C_A" # mis-mapping 
# "chr17_40061856_G_C" # looks okay 
# "chr4_147908994_T_C" # looks okay 
# "chrX_68803487_C_T" 
# "chrX_66719643_C_G" # looks okay

twins_filtered_vaf[sum_tumour == 10 & sum_normal_PD62341 == 1 & sum_normal_PD63383 == 1] # 12
# always PD62341aa or PD62341h and PD63383bb so very likely contamination 
# "chr11_134158143_C_T" 
# "chr14_104122986_G_T"
# "chr15_48353502_C_A" 
# "chr15_56691722_C_T" 
# "chr16_4003400_G_A" 
# "chr1_69984580_G_A" # looks like it's been deleted in tumour samples 
# "chr2_133311125_G_T" 
# "chr3_49878199_C_A"
# "chr4_178731458_A_G"
# "chr5_136349883_C_T" 
# "chr5_37790015_G_A" 
# "chr9_11364991_T_A"  

######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
# ANALYSIS OF THE SHARED MUTATIONS BETWEEN TUMOUR AND PD62341 NORMAL

# note that I am allowing presence in 1 sample from PD63383 because there could have been contamination
twins_filtered_vaf[sum_tumour >= 6 & sum_normal_PD62341 >= 4 & sum_normal_PD63383 <= 1] # 8 
# "chr14_105458006_C_A" # looks good 
# "chr16_5479739_C_T"  # looks good 
# "chr17_18402034_C_A" # mapping issues  
# "chr17_33422229_C_A" # 
# "chr2_95662131_G_A" # mapping issues  
# "chr3_50106043_C_T" # looks good 
# "chr3_62055057_C_G" # looks good   
# "chr3_62055077_G_C" # looks good  
# some do look real!

# which samples from PD62341 share these mutations (of those that look real)?
# "chr14_105458006_C_A" # v, q, aa, ad, n, h (>= 7 in all): clonal in ad?
# "chr16_5479739_C_T" # q, aa, ad, h, n (>= 12 in all, 3 in v): looks clonal pretty much
# "chr1_148424284_G_A" # v (4), q (9), h (7), n (5): not clonal
# "chr2_14504041_A_C" # v (7), q (9), aa (6), n (7) # generally low vaf
# "chr3_50106043_C_T" # q, aa, ad, h, n (>= 8 in all, 2 in v)
# "chr3_62055057_C_G"  # q, aa, h, n (>= 9, 2 in v, 1 in ad) 
# "chr3_62055077_G_C" # q, aa, h, n (>= 9 in all, 2, in v, 2 in ad)

twins_filtered_vaf[sum_tumour >= 6 & sum_normal_PD62341 <= 1 & sum_normal_PD63383 >= 4] # 2
# chr2_111459022_G_T
# chr4_190122142_C_A

######################################################################################################
######################################################################################################
# MUTATIONS SHARED B/N PD62341 NORMAL AND TUMOUR IN MORE DETAIL

# are there mutations shared with only a specific normal tissue of PD62341?
twins_filtered_vaf[sum_tumour >= 6 & sum_normal_PD62341 >= 1 & sum_normal_PD63383 == 0] # 21
# "chr10_100754461_C_T"
# "chr12_106354695_C_T"
# "chr13_60216504_A_C"
# "chr12_61020950_C_A"  # looks okay 
# "chr14_105458006_C_A" # looks okay
# "chr16_32621521_C_A"  # poor mapping 
# "chr17_78256883_C_T" # looks okay
# "chr19_3961910_G_A" # looks okay  
# "chr20_44114996_C_T"
# "chr21_45729938_T_A" # looks okay
# "chr22_17774668_G_A" # looks quite okay 
# "chr22_30553442_G_T" # generally looks okay
# "chr2_199565129_G_A" 
# "chr2_57814739_A_C" # poor mapping quality   
# "chr3_62055057_C_G" # looks okay  
# "chr3_62055077_G_C" # looks okay
# "chr6_95827754_C_A" # rather poor mapping quality  
# "chr7_139659050_G_A" # looks okay 
# "chr8_46678612_T_C" # looks okay
# "chrX_46375010_C_T" # looks okay  
# "chrX_56448796_C_A" # looks okay (but there may be some insertions)
