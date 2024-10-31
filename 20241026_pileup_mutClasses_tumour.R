# Script to analyse the pileup (high quality, run 15/10/2024)
# 2024-10-25
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
muts = read.table('Data/mutations_include_20241030_610.txt') %>% unlist()
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

# Add presence / absence based on MTR data
twins_filtered_mtr[, sum_tumour := rowSums(.SD>=4), .SDcols = samples_tumour_mtr]
twins_filtered_mtr[, sum_normal := rowSums(.SD>=4), .SDcols = samples_normal_mtr]
twins_filtered_mtr[, sum_PD62341 := rowSums(.SD>=4), .SDcols = samples_PD62341_mtr]
twins_filtered_mtr[, sum_PD63383 := rowSums(.SD>=4), .SDcols = samples_PD63383_mtr]
twins_filtered_mtr[, sum_tumour_PD62341 := rowSums(.SD>=4), .SDcols = samples_tumour_PD62341_mtr]
twins_filtered_mtr[, sum_tumour_PD63383 := rowSums(.SD>=4), .SDcols = samples_tumour_PD63383_mtr]
twins_filtered_mtr[, sum_normal_PD62341 := rowSums(.SD>=4), .SDcols = samples_normal_PD62341_mtr]
twins_filtered_mtr[, sum_normal_PD63383 := rowSums(.SD>=4), .SDcols = samples_normal_PD63383_mtr]

# Add presence / absence based on VAF data
twins_filtered_vaf[, sum_tumour := rowSums(.SD>=0.1), .SDcols = samples_tumour_vaf]
twins_filtered_vaf[, sum_normal := rowSums(.SD>=0.1), .SDcols = samples_normal_vaf]
twins_filtered_vaf[, sum_PD62341 := rowSums(.SD>=0.1), .SDcols = samples_PD62341_vaf]
twins_filtered_vaf[, sum_PD63383 := rowSums(.SD>=0.1), .SDcols = samples_PD63383_vaf]
twins_filtered_vaf[, sum_tumour_PD62341 := rowSums(.SD>=0.1), .SDcols = samples_tumour_PD62341_vaf]
twins_filtered_vaf[, sum_tumour_PD63383 := rowSums(.SD>=0.1), .SDcols = samples_tumour_PD63383_vaf]
twins_filtered_vaf[, sum_normal_PD62341 := rowSums(.SD>=0.1), .SDcols = samples_normal_PD62341_vaf]
twins_filtered_vaf[, sum_normal_PD63383 := rowSums(.SD>=0.1), .SDcols = samples_normal_PD63383_vaf]

######################################################################################################
######################################################################################################
# TUMOUR EVOLUTION

# Tumour-specific mutations 

# I will allow up to 3 normal samples due to the possibility of contamination
# I suspect that skin samples (aa, bb) and liver sample (PD62341h) are contaminated with reads from the tumour

# Mutations found in all tumour samples 
twins_filtered_mtr[sum_tumour == 10 & sum_normal == 0] # 1 
# "chr2_199565129_G_A" # looks okay

twins_filtered_mtr[sum_tumour == 10 & sum_normal == 1] # 4 
# "chr16_32621521_C_A" # mapping quality not excellent (double check)
# "chr17_78256883_C_T" # looks okay
# "chr22_17774668_G_A" # looks okay 
# "chr22_30553442_G_T" # looks okay

twins_filtered_mtr[sum_tumour == 10 & sum_normal == 2] # 11
# "chr13_60216504_A_C" # looks okay
# "chr15_48353502_C_A" # looks okay
# "chr15_56691722_C_T" # looks okay
# "chr16_4003400_G_A" # looks okay
# "chr1_69984580_G_A" # can estimate contamination from this one!!
# "chr2_133311125_G_T" # looks okay
# "chr5_136349883_C_T" # looks okay
# "chr5_37790015_G_A" # looks okay
# "chr6_95827754_C_A" # mapping could be better but still decent
# "chr9_11364991_T_A" # looks okay
# "chrX_68803487_C_T" # looks okay

twins_filtered_mtr[sum_tumour == 10 & sum_normal == 3] # 9 
# it's all bb, h and aa (skin + liver)
# chr10_100754461_C_T # looks okay
# chr11_134158143_C_T # looks okay
# chr14_104500332_C_T # looks okay
# chr15_23462705_C_A # lots of mis-mapped reads 
# chr17_40061856_G_C # looks okay
# chr4_130207996_C_T # looks okay
# chr4_147908994_T_C # looks okay
# chr4_178731458_A_G # looks okay
# chrX_66719643_C_G # looks okay

twins_filtered_mtr[sum_tumour == 9 & sum_normal == 0] # 2 
# chr12_106354695_C_T # looks okay, missing in PD62341ak
# chrX_115495236_T_C # looks okay to me (double check but I would believe it), missing in PD62341b

twins_filtered_mtr[sum_tumour == 9 & sum_normal == 1] # 7
# chr12_61020950_C_A # looks okay
# chr14_104122986_G_T # looks okay (nb some mis-mapping)
# chr7_115240063_T_C # looks okay
# chr7_120677593_C_T # looks okay
# chr8_46678612_T_C # looks okay
# chrX_46375010_C_T # looks okay
# chrX_56448796_C_A # looks okay

twins_filtered_mtr[sum_tumour == 9 & sum_normal == 2] # 13
# "chr12_10337296_G_T" # variable mapping (double check)
# "chr12_39248376_A_T" # looks okay
# "chr18_71485446_G_A" # looks like del of the other copy - can estimate contamination
# "chr21_45729938_T_A" # looks okay
# "chr22_16694683_C_A" # looks okay
# "chr2_119170376_C_T" # looks okay
# "chr2_187893783_G_T" # looks okay
# "chr3_94699090_G_T" # looks okay
# "chr4_180479634_G_A" # looks okay
# "chr5_112010886_G_A" # looks okay
# "chr7_12297469_T_G" # looks okay
# "chr8_54112201_C_A" # looks okay
# "chr9_94529993_C_T" # quite a lot of mis-mapped reads (double check)

# there are two mutations next to each other that show a very similar pattern and VAF 
# the other is: chr8_54112197_T_C, and is called as present in 9 tumour samples and 3 normal
# normal chr8_54112197_T_C: PD62341aa, h, PD63383 bb
# normal chr8_54112201_C_A: PD62341h, PD63383bb (3 reads in PD62341aa)

twins_filtered_mtr[sum_tumour == 9 & sum_normal == 3] # 21
# "chr13_100405282_G_T" # looks okay
# "chr15_76646556_T_C" # looks okay  
# "chr17_5029328_G_C" # looks okay 
# "chr1_119005653_G_C" # mis-mapped reads (double check cause it really doesn't look too bad) 
# "chr1_60839899_G_A" # looks great
# "chr2_124573116_C_T" # looks great  
# "chr2_1420635_A_T" # looks great # there is a beautiful mutation next to this one 
# "chr2_2023708_C_G"  # looks great  
# "chr2_232306029_G_C" # looks great
# "chr3_49878199_C_A" # looks great  
# "chr4_46067117_A_G" # difficult region but looks okay
# "chr4_93081874_T_A" # looks great  
# "chr5_119778240_G_C" # poor mapping quality (double check)
# "chr7_139659050_G_A" # looks okay
# "chr7_25564729_G_A" # looks okay
# "chr8_131571989_A_T" # looks okay
# "chr8_54112197_T_C" # looks okay
# "chr8_97785584_T_C" # looks okay 
# "chr9_12824689_G_A" # looks okay 
# "chrX_124531709_C_T" # looks okay  
# "chrX_83934262_C_G" # looks okay

# very interesting mutation chr1_60839875_T_C 
# seems to be on the other chromosome 
# can compare reads with different muts to determine contamination in different tumour samples
# also, PD62341q for example has 5 of those reads that belong to the tumour chromosome 
# we don't know tho if the loss of chr1 is subclonal in some tumour samples - what about this?

twins_filtered_mtr[sum_tumour == 8 & sum_normal == 0] # empty  

twins_filtered_mtr[sum_tumour == 8 & sum_normal == 1] # 1
# "chr2_57814739_A_C" # mixed feelings - looks okay but quite a lot of insertions, double check

twins_filtered_mtr[sum_tumour == 8 & sum_normal == 2] # 2
# chr2_78054577_G_A # looks okay
# chr3_153319439_A_G # looks okay

twins_filtered_mtr[sum_tumour == 7 & sum_normal == 0] # 0  
# chr7_152562130_G_A # looks good 

twins_filtered_mtr[sum_tumour == 7 & sum_normal == 1] # 2
# chr12_84029949_C_A # looks okay
# chr6_77631804_C_A # looks okay

twins_filtered_mtr[sum_tumour == 7 & sum_normal == 2] # 1
# chr2_158600954_G_A # looks okay

twins_filtered_mtr[sum_tumour == 7 & sum_normal == 3] # 1
# chr8_124634298_G_T # deletions, should have been removed?

# Not sure the below are worth it since we don't know if they tell us anything about the tumour phylogeny or not 
twins_filtered_mtr[sum_tumour == 6 & sum_normal == 0] # 0 

twins_filtered_mtr[sum_tumour == 6 & sum_normal == 1] # 3
# chr10_51734221_A_T # looks okay
# chr12_112652964_C_T # looks okay
# chr15_94939486_G_A # looks okay

twins_filtered_mtr[sum_tumour == 6 & sum_normal == 2] # 1
# chr18_64267234_G_A # looks okay, inspect further - is this on the deleted copy?

twins_filtered_mtr[sum_tumour == 5 & sum_normal == 0] # 3
# chr16_17805832_C_T # looks okay 
# chr7_153784080_C_G # don't think this is real
# chr8_21444725_G_A # looks okay

twins_filtered_mtr[sum_tumour == 5 & sum_normal == 1] # 7
# "chr11_26060345_G_T" # looks okay
# "chr12_129105882_A_G" # poor mapping 
# "chr13_103318073_T_A" # looks okay
# "chr14_54164980_C_G" # looks okay
# "chr17_70958446_C_T" # looks okay
# "chr6_168731130_C_T" # looks okay  
# "chr6_8310684_G_T" # looks okay

twins_filtered_mtr[sum_tumour == 4 & sum_normal == 0] # 1
# chr13_28953051_G_T # looks okay

twins_filtered_mtr[sum_tumour == 4 & sum_normal == 1] # 3
# "chr12_114403258_G_A" # looks okay
# "chr18_1379467_G_T"  # looks okay
# "chr4_135461656_G_A" # looks okay

twins_filtered_mtr[sum_tumour == 3 & sum_normal == 0] # 4  
# "chr11_43029745_A_T" # looks okay
# "chr12_19555916_C_T" # looks okay
# "chr12_82832044_C_A" # looks okay
# "chr2_231591536_G_T" # looks okay
# "chr5_141433408_A_T" # indels - why not thrown away?

twins_filtered_mtr[sum_tumour == 3 & sum_normal == 1] # 5
# "chr1_145588334_T_C" # poor mapping (not extremely bad)
# "chr1_148952577_G_A" # poor mapping 
# "chr4_179587218_G_A" # looks real!
# "chr8_141860196_G_A" # looks okay

######################################################################################################
# Mutations specific to PD63383 tumour

twins_filtered_mtr[sum_normal == 0 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 0] # 12
# "chr10_100308403_C_T" # looks okay
# "chr11_80115755_C_T" # looks okay
# "chr19_43348626_C_A" # not sure, double check
# "chr1_51714096_G_C" # looks okay  
# "chr2_72939205_C_T" # looks okay 
# "chr2_82147172_T_C" # looks okay 
# "chr5_157248612_A_G" # looks okay
# "chr5_54017866_G_T" # looks okay
# "chr5_54017866_G_T" # looks okay  
# "chr7_13677323_A_G" # mapping not excellent, double check
# "chrX_72671786_A_T" # poor mapping + strand bias   
# "chrX_77681162_G_A" # looks okay

muts_PD63383_tumour1 = c("chr10_100308403_C_T", "chr11_80115755_C_T", "chr19_43348626_C_A", "chr1_51714096_G_C",
                        "chr2_72939205_C_T", "chr2_82147172_T_C", "chr5_157248612_A_G", "chr5_28014472_C_T",
                        "chr5_54017866_G_T", "chrX_48453603_G_A", "chrX_72671786_A_T", "chrX_77681162_G_A")
twins_info = twins_dt[, c('mut_ID', 'Gene', 'Transcript', 'RNA', 'CDS', 'Protein', 'Type', 'Effect'), with=FALSE]
twins_info[mut_ID %in% muts_PD63383_tumour1] # nothing in protein coding genes

twins_filtered_mtr[sum_normal_PD63383 == 1  & sum_normal_PD62341 == 0 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 0] # 15
# "chr10_31333522_G_A" # looks okay # bb
# "chr12_43132316_G_A" # looks okay but close to poly-T # bb
# "chr12_80762031_G_A" # looks okay # bb
# "chr1_160088990_G_A" # looks okay # bb
# "chr22_25281865_C_T" # looks okay # bb
# "chr3_137508691_C_A" # looks okay # bb
# "chr6_86549064_C_T" # looks okay # bb
# "chr6_92179408_A_G" # looks okay # bb
# "chr7_148446413_G_A" # looks okay # bb
# "chr7_49732059_G_A" # looks okay # bb 
# "chr9_98275090_C_T" # looks okay # bb
# "chrX_117418608_G_A" # looks okay # bb
# "chrX_119915655_G_A" # looks okay # bb
# "chrX_36779694_C_T" # looks okay # bb
# "chrX_70352624_G_A" # looks okay # bb

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

twins_filtered_mtr[sum_normal_PD63383 == 2  & sum_normal_PD62341 == 0 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 0, mut_ID] # 0
# empty

twins_filtered_mtr[sum_normal_PD63383 == 3  & sum_normal_PD62341 == 0 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 0, mut_ID] # 0
# empty 

twins_filtered_mtr[sum_normal_PD63383 == 0  & sum_normal_PD62341 == 1 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 0, mut_ID] # 0
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
twins_filtered_mtr[sum_normal_PD62341 == 0 & sum_normal_PD63383 == 0 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 8] # 0
twins_filtered_mtr[sum_normal_PD62341 == 1 & sum_normal_PD63383 == 0 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 8] # 0
twins_filtered_mtr[sum_normal_PD62341 == 2 & sum_normal_PD63383 == 0 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 8] # 0

# no muts present in 7/8 PD62341 and none PD63383 tumour
twins_filtered_mtr[sum_normal == 0 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 7] # 0
twins_filtered_mtr[sum_normal == 1 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 7] # 0
twins_filtered_mtr[sum_normal == 2 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 7] # 0

# muts present in 6/8 PD62341 tumour and none PD63383 tumour (nothing real)
twins_filtered_mtr[sum_normal == 0 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 6] # 0
twins_filtered_mtr[sum_normal == 1 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 6] # 0
twins_filtered_mtr[sum_normal == 2 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 6] # 0

# muts present in 5/8 PD62341 and none PD63383 tumour (2 real ones!)
twins_filtered_mtr[sum_normal == 0 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 5] # 2
# "chr16_17805832_C_T" # looks okay 
# ae, ag, aj, b, u (interesting - all non-primary tumour so that would suggest the PD63383 emerged from primary)
# "chr7_153784080_C_G" # not believable, double check
twins_filtered_mtr[sum_normal == 1 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 5] # 1
# "chr6_8310684_G_T" # looks okay!
# ae, ag, aj, b, u (again!) and present in PD62341h normal - likely contaminated
twins_filtered_mtr[sum_normal == 2 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 5] # 0
# empty

# muts present in 4/10 PD62341 and none PD63383 tumour
twins_filtered_mtr[sum_normal == 0 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 4] # 1
# "chr13_28953051_G_T" # looks okay
# present in ae, ag, aj, u (3 reads in b and ap!)

twins_filtered_mtr[sum_normal == 1 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 4] # 2
# "chr12_114403258_G_A" # looks okay
# ae, aj, b, u (2 reads in ag and aa, 1 in am and ap)
# "chr18_1379467_G_T" # looks okay 
# ae, ag, aj, u (3 reads in b, 2 in aa, 1 in v, am, also in h but likely contaminated)

twins_filtered_mtr[sum_normal == 2 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 4] # 0

# muts present in 3/10 PD62341 and none PD63383 tumour
twins_filtered_mtr[sum_normal == 0 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 3] # 3
# "chr11_43029745_A_T" very low MTR numbers everywhere but looks okay on Jb
# "chr12_19555916_C_T" very low MTR numbers everywhere but looks okay on Jb
# "chr2_231591536_G_T" looks good (present in aj and ae, less in u)

twins_filtered_mtr[sum_normal == 1 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 3] # 3
# "chr1_248574696_C_G" looks okay

twins_filtered_mtr[sum_normal == 2 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 3, mut_ID] # 0
# "chr10_119182885_G_C" # poor mapping 

# muts present in 2/10 PD62341 and none PD63383 tumour
# not sure these are worth investigating in detail because ultimately each tumour sample would have been seeded by more than one cell probably
# so not very likely that you would identify anything that is truly clonal 
twins_filtered_mtr[sum_normal == 0 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 2] # 5
# "chr13_62110424_G_A" # looks good 
# "chr1_14644486_C_A" # looks good  
# "chr1_38311094_G_T # looks good 
# "chr1_83136026_A_T" # not great, double check 
# "chrX_9804641_G_A" # looks okay

twins_filtered_mtr[sum_normal == 1 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 2] # 4
# "chr11_18922568_G_T" # double check (repetitive region)
# "chr17_1199637_C_G" # poor mapping 
# "chr9_100061581_C_T" # maybe? weird region so may have thrown Illumina off (double check)
# "chr9_67625770_G_A" # some mapping issues 

twins_filtered_mtr[sum_normal == 2 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 2] # 2
# "chr19_50615732_C_T" # double check
# "chr9_39590690_G_A" # poor mapping 

######################################################################################################
# Mutations shared between tumour and one twin 

# Mutations present in ALL tumour samples and all / some normal PD62341 samples 
twins_filtered_mtr[sum_tumour == 10 & sum_normal_PD62341 == 6 & sum_normal_PD63383 == 0] # 1
# "chr14_105458006_C_A" our favourite one!

twins_filtered_mtr[sum_tumour == 10 & sum_normal_PD62341 == 5 & sum_normal_PD63383 == 0] # 0

twins_filtered_mtr[sum_tumour == 10 & sum_normal_PD62341 == 4 & sum_normal_PD63383 == 0] # 2
# "chr3_62055057_C_G" aa and h (possibly contamination), n and q (heart and pancreas), missing in v (?) and ad
# "chr3_62055077_G_C" aa and h (possibly contamination), n and q (heart and pancreas), missing in v (?) and ad

twins_filtered_mtr[sum_tumour == 10 & sum_normal_PD62341 == 3 & sum_normal_PD63383 == 0] # 0

twins_filtered_mtr[sum_tumour == 10 & sum_normal_PD62341 == 2 & sum_normal_PD63383 == 0] # 1
# chr6_95827754_C_A not the best mapping quality but fine, aa and h (possibly contamination) 

twins_filtered_mtr[sum_tumour == 10 & sum_normal_PD62341 == 1 & sum_normal_PD63383 == 0] # 4
# "chr16_32621521_C_A" # mapping quality not excellent (double check because it doesn't look that bad)
# "chr17_78256883_C_T" # looks okay
# "chr22_17774668_G_A" # looks okay 
# "chr22_30553442_G_T" # looks okay

# Mutations present in 9 tumour samples and all / some PD62341 normal samples 
twins_filtered_mtr[sum_tumour == 9 & sum_normal_PD62341 == 6 & sum_normal_PD63383 == 0] # 0
twins_filtered_mtr[sum_tumour == 9 & sum_normal_PD62341 == 5 & sum_normal_PD63383 == 0] # 0
twins_filtered_mtr[sum_tumour == 9 & sum_normal_PD62341 == 4 & sum_normal_PD63383 == 0] # 0
twins_filtered_mtr[sum_tumour == 9 & sum_normal_PD62341 == 3 & sum_normal_PD63383 == 0] # 1
# chr7_139659050_G_A looks okay
# also there is a very cool mutation on the other chromosome (in all samples?) chr7_139659034_T_G
# missing in ak, which is technically tumour (primary)

twins_filtered_mtr[sum_tumour == 9 & sum_normal_PD62341 == 2 & sum_normal_PD63383 == 0] # 2
# chr18_71485446_G_A # looks okay and on the segment lost in the other chromosome 
# chr21_45729938_T_A # looks okay
# chr2_119170376_C_T # very significant strand bias but technically fine? double check
# all missing in ak 

twins_filtered_mtr[sum_tumour == 9 & sum_normal_PD62341 == 1 & sum_normal_PD63383 == 0] # 5
# mostly in h so all of these are likely contamination
# chr12_61020950_C_A # looks okay
# chr7_120677593_C_T # looks okay
# chr8_46678612_T_C # looks okay
# chrX_46375010_C_T # looks okay
# chrX_56448796_C_A # looks okay
# all missing in ak 

# Allow presence in 1 normal sample from PD63383 (because of contamination)
twins_filtered_mtr[sum_tumour == 10 & sum_normal_PD62341 == 6 & sum_normal_PD63383 == 1] # 0
twins_filtered_mtr[sum_tumour == 10 & sum_normal_PD62341 == 5 & sum_normal_PD63383 == 1] # 3
# chr16_5479739_C_T # looks real (absent from PD62341v, present in PD63383bb)
# chr2_95662131_G_A # looks semi-okay (double check) (absent from PD62341v, present in PD63383bb)
# chr3_50106043_C_T # looks okay (absent from PD62341v, present in PD63383bb)

twins_filtered_mtr[sum_tumour == 10 & sum_normal_PD62341 == 4 & sum_normal_PD63383 == 1] # 0
twins_filtered_mtr[sum_tumour == 10 & sum_normal_PD62341 == 3 & sum_normal_PD63383 == 1] # 1
# chr5_44907911_G_A # looks okay (present in PD62341q, aa, h and PD63383bb)

twins_filtered_mtr[sum_tumour == 10 & sum_normal_PD62341 == 2 & sum_normal_PD63383 == 1] # 9
# always PD62341aa, PD62341h and PD63383bb so very likely contamination 
# "chr10_100754461_C_T" # looks okay
# "chr11_134158143_C_T" # looks okay
# "chr14_104500332_C_T" # looks okay
# "chr15_23462705_C_A" # mis-mapping 
# "chr17_40061856_G_C" # looks okay 
# "chr4_130207996_C_T" # looks okay 
# "chr4_147908994_T_C" # looks okay 
# "chr4_178731458_A_G" # looks okay
# "chrX_66719643_C_G" # looks okay

twins_filtered_mtr[sum_tumour == 10 & sum_normal_PD62341 == 1 & sum_normal_PD63383 == 1] # 10
# always PD62341aa or PD62341h and PD63383bb so very likely contamination 
# "chr13_60216504_A_C" 
# "chr15_48353502_C_A" 
# "chr15_56691722_C_T" 
# "chr16_4003400_G_A" 
# "chr1_69984580_G_A" # looks like it's been deleted in tumour samples 
# "chr2_133311125_G_T" 
# "chr5_136349883_C_T" 
# "chr5_37790015_G_A" 
# "chr9_11364991_T_A"  
# "chrX_68803487_C_T"

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
twins_filtered_mtr[sum_tumour >= 6 & sum_normal_PD62341 >= 4 & sum_normal_PD63383 <= 1] # 11 
# "chr14_105458006_C_A" # looks good 
# "chr16_5479739_C_T"  # looks good 
# "chr16_88009874_G_A"  # not sure 
# "chr20_44114996_C_T" # ok but check VAF 
# "chr2_95662131_G_A" # mapping issues  
# "chr3_50106043_C_T" # looks good 
# "chr3_62055057_C_G" # looks good   
# "chr3_62055077_G_C" # looks good  
# some do look real!

# which samples from PD62341 share these mutations (of those that look real)?
# "chr14_105458006_C_A" # v, q, aa, ad, n, h (>= 7 in all): clonal in ad?
# "chr16_5479739_C_T" # q, aa, ad, h, n (>= 12 in all, 3 in v): looks clonal pretty much
# "chr1_148424284_G_A" # v (4), q (9), h (7), n (5): not clonal
# "chr2_14504041_A_C" # v (7), q (9), aa (6), n (7) # generally low MTR
# "chr3_50106043_C_T" # q, aa, ad, h, n (>= 8 in all, 2 in v)
# "chr3_62055057_C_G"  # q, aa, h, n (>= 9, 2 in v, 1 in ad) 
# "chr3_62055077_G_C" # q, aa, h, n (>= 9 in all, 2, in v, 2 in ad)

twins_filtered_mtr[sum_tumour >= 6 & sum_normal_PD62341 <= 1 & sum_normal_PD63383 >= 4] # 0 

######################################################################################################
######################################################################################################
# MUTATIONS SHARED B/N PD62341 NORMAL AND TUMOUR IN MORE DETAIL

# are there mutations shared with only a specific normal tissue of PD62341?
twins_filtered_mtr[sum_tumour >= 6 & sum_normal_PD62341 >= 1 & sum_normal_PD63383 == 0] # 20
# "chr12_61020950_C_A"  # looks okay 
# "chr14_105458006_C_A" # looks okay
# "chr16_32621521_C_A"  # poor mapping 
# "chr17_78256883_C_T" # looks okay
# "chr18_64267234_G_A" # some mapping issues but generally seems okay
# "chr18_71485446_G_A" # looks okay 
# "chr21_45729938_T_A" # looks okay
# "chr22_17774668_G_A" # looks quite okay 
# "chr22_30553442_G_T" # generally looks okay
# "chr2_57814739_A_C" # poor mapping quality   
# "chr3_62055057_C_G" # looks okay  
# "chr3_62055077_G_C" # looks okay
# "chr6_58221412_T_C" # muts only on reads with poor mapping   
# "chr6_95827754_C_A" # rather poor mapping quality  
# "chr7_120677593_C_T" # insertions
# "chr7_139659050_G_A" # looks okay 
# "chr8_46678612_T_C" # looks okay
# "chrX_46375010_C_T" # looks okay  
# "chrX_56448796_C_A" # looks okay (but there may be some insertions)

# count samples with mutations (for all mutations)
colSums(twins_filtered_mtr[sum_tumour >= 6 & sum_normal_PD62341 >= 1 & sum_normal_PD63383 == 0,c(samples_normal_PD62341_mtr), with=FALSE] >= 4)
# v: 5, q: 9, aa: 12, ad: 6, h: 20, n: 5
colSums(twins_filtered_mtr[sum_tumour >= 6 & sum_normal_PD62341 >= 1 & sum_normal_PD63383 == 0,c(samples_normal_PD62341_mtr), with=FALSE] >= 10) # less likely to be contamination
# v: 0, q: 3, aa: 4, ad: 1, h: 4, n: 3

# run analysis for mutations with confirmed high quality 
muts_hq = c("chr12_61020950_C_A", "chr14_105458006_C_A", "chr17_78256883_C_T", "chr18_71485446_G_A",
            "chr19_3961910_G_A", "chr21_45729938_T_A", "chr22_17774668_G_A", "chr3_62055057_C_G",
            "chr3_62055077_G_C", "chr7_139659050_G_A", "chr8_46678612_T_C", "chrX_46375010_C_T", "chrX_56448796_C_A")
paste('Number of high-quality mutations shared between tumour and PD62341 normal:', length(muts_hq))

colSums(twins_filtered_mtr[mut_ID %in% muts_hq, c(samples_normal_PD62341_mtr), with=FALSE] >= 4)
# v: 1, q: 6, aa: 5, ad: 1, h: 12, n: 3
# the one that is not 'present' in h is there at 3 reads, and 5 reads in q (contamination?)
# h is at quite low depth there, so VAF is actually 0.1875 

colSums(twins_filtered_mtr[mut_ID %in% muts_hq, c(samples_normal_PD62341_mtr), with=FALSE] >= 10) # less likely to be contamination
# v: 0, q: 2, aa: 3, ad: 1, h: 3, n: 3

# we can also check VAF (> 0.1) to account for differences in coverage
colSums(twins_vaf[mut_ID %in% muts_hq, c(samples_normal_PD62341_vaf), with=FALSE] >= 0.1)
# v: 1, q: 6, aa: 4, ad: 1, h: 13, n: 4

# how does this compare to tumour?
colSums(twins_vaf[mut_ID %in% muts_hq, c(samples_tumour_vaf), with=FALSE] >= 0.1)
# all are 13 but PD62341ak (maybe contaminated?? some other cells? or duplicated chromosomes so lower VAF?? - but then this wouldn't drop it below 0.1)

# but then it is very rare that any of these mutations seem clonal
# I imagine this can be because of contamination with other cells (which is not known)
# for normal samples, bear in mind they are not monoclonal so really hard to tell
# also if there are really precursors you could have them at lower levels / freq so lower VAF
colSums(twins_vaf[mut_ID %in% muts_hq, c(samples_tumour_vaf), with=FALSE] >= 0.4)
colSums(twins_vaf[mut_ID %in% muts_hq, c(samples_normal_PD62341_vaf), with=FALSE] >= 0.4)
# 1 case in PD62341n and PD62341d 

# Examine VAF relationship 
mut_shared_tumour_normal = twins_vaf[mut_ID %in% muts_hq, c('mut_ID', samples_tumour_vaf, samples_normal_PD62341_vaf), with=FALSE]
mut_shared_melt = melt(mut_shared_tumour_normal, id.vars = 'mut_ID')
mut_shared_melt[, sample := tstrsplit(variable, '_', fixed=TRUE, keep = 1)]
mut_shared_melt[, status := as.factor(fcase( 
  sample %in% samples_normal, 'normal', # differs from 0.5 so biased and maybe we don't want it
  sample %in% samples_tumour, 'tumour'
))]
mut_shared_melt[, twin := as.factor(fcase( 
  sample %in% samples_PD62341, 'PD62341', # differs from 0.5 so biased and maybe we don't want it
  sample %in% samples_PD63383, 'PD63383'
))]
mut_shared_melt[, sample_type := as.factor(paste(status, twin, sep = '_'))]
mut_shared_melt[, mut_ID := factor(mut_ID, levels = c("chr14_105458006_C_A", "chr3_62055057_C_G", "chr3_62055077_G_C",
                                                      "chr8_46678612_T_C", "chr7_139659050_G_A", "chr18_71485446_G_A", 
                                                      "chr12_61020950_C_A", "chr17_78256883_C_T", "chr19_3961910_G_A", 
                                                      "chr21_45729938_T_A", "chr22_17774668_G_A", "chrX_46375010_C_T", "chrX_56448796_C_A"))] 

ggplot(mut_shared_melt, aes(x = value, y = mut_ID, col = status))+
  geom_point(size=2, position=position_dodge(.75))+
  theme_bw(base_size = 14)+
  scale_color_manual(values = c(col_normal, col_tumour))+
  labs(x = 'VAF', y = 'Mutation')+
  ggtitle('VAF in tumour vs normal')
ggsave('Results/20241020_p2_vaf_dist_for_each_shared_mut13.pdf', width=9, height=7)

# so would we conclude that bc most mutations are shared, the tumour arose in the liver?
# is VAF higher in the liver than other tissues? for mutations which are shared?
mut_shared_melt[, sample_name := factor(fcase(
  sample == 'PD62341h', 'liver',
  sample == 'PD62341q', 'pancreas',
  sample == 'PD62341n', 'heart',
  sample == 'PD62341v', 'spleen',
  sample == 'PD62341ad', 'cerebellum',
  sample == 'PD62341aa', 'skin'
))]
mut_shared_melt[, sample_name := factor(sample_name, 
                                        levels=c('liver', 'pancreas', 'heart', 'spleen', 'cerebellum', 'skin'))]

ggplot(mut_shared_melt[status=='normal'], aes(x = value, y = mut_ID, col = sample_name))+
  geom_point(size=3)+
  theme_bw(base_size = 14)+
  scale_color_manual(values = c('#af0101', '#cc7a3d', '#e6c994',
                                '#ffaa6a', '#74a892', '#004367'))+
  labs(x = 'VAF', y = 'Mutation', color = 'Sample')+
  ggtitle('VAF distribution in normal samples')
ggsave('Results/20241020_p2_vaf_dist_for_each_shared_mut13_normalsamples.pdf', width=10, height=7)

# can I see mutations only present in the liver of PD62341 and the tumour?
# there are 6: there is 1 present in PD62341q (5 reads, 3 reads in PD62341h so can be due to coverage differences)
twins_filtered_mtr[mut_ID %in% muts_hq & sum_tumour >= 6 & sum_normal_PD62341 == 1]
colSums(twins_filtered_mtr[mut_ID %in% muts_hq, c(samples_normal_PD62341_mtr), with=FALSE] >= 10) # less likely to be contamination

# note that a real mutation could be present in PD63383 due to contamination
twins_filtered_mtr[sum_tumour >= 6 & sum_normal_PD62341 >= 1 & sum_normal_PD63383 <= 1] # 108 

# better plots for VAF
mut_shared_melt[order(sample, value)]
mut_shared_melt[, sample := as.factor(sample)]

for (s in mut_shared_melt[sample_type=='normal_PD62341',sample] %>% unlist() %>% unique()){
  ggplot(mut_shared_melt[sample==s], aes(x=mut_ID, y=value))+
    geom_line(group = 1)+
    geom_point(size=1.5)+
    geom_area(aes(fill = mut_ID))+
    scale_fill_manual(values = col_normal_PD62341)+
    theme_bw(base_size = 10)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    labs(x = 'Mutation', y = 'VAF')+
    ggtitle(glue('VAF distribution in sample {s}'))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    ylim(c(0, 0.6))
  ggsave(glue('Results/20241020_p2_vaf_dist_13mut_normal_{s}.pdf'), width = 6, height = 2.5)
}

# look at each mutation one by one 
mut_hq13 = twins_vaf[mut_ID %in% muts_hq, c('mut_ID', samples_tumour_vaf, samples_normal_vaf), with=FALSE]
mut_hq13[, sum_all := rowSums(.SD >= 0.1), .SDcols = 2:23]
mut_hq13[, sum_normal := rowSums(.SD >= 0.1), .SDcols = samples_normal_vaf]
mut_hq13[, sum_tumour := rowSums(.SD >= 0.1), .SDcols = samples_tumour_vaf]
mut_hq13[, sum_normal_PD62341 := rowSums(.SD >= 0.1), .SDcols = samples_normal_PD62341_vaf]
mut_hq13[, sum_tumour_PD62341 := rowSums(.SD >= 0.1), .SDcols = samples_tumour_PD62341_vaf]
mut_hq13[, sum_normal_PD63383 := rowSums(.SD >= 0.1), .SDcols = samples_normal_PD63383_vaf]
mut_hq13[, sum_tumour_PD63383 := rowSums(.SD >= 0.1), .SDcols = samples_tumour_PD63383_vaf]
mut_hq13[order(sum_all)] # note that there are reads in PD63383bb but these are likely contamination (skin!!)

######################################################################################################
# Visualization: heatmaps  

# order the columns such that tumour and normal samples are separated and ordered by # of mutations present
mut_hq13_sub = mut_hq13[, 1:17]
setcolorder(mut_hq13_sub, c('mut_ID', 'PD63383aq_VAF', 'PD63383ap_VAF',
                            'PD62341ae_VAF', 'PD62341ag_VAF', 'PD62341aj_VAF', 'PD62341am_VAF', 'PD62341ap_VAF', 'PD62341b_VAF', 'PD62341u_VAF', 'PD62341ak_VAF',
                            'PD62341h_VAF', 'PD62341q_VAF', 'PD62341n_VAF', 'PD62341aa_VAF', 'PD62341ad_VAF', 'PD62341v_VAF'))
row_order = c(2, 8, 9, 11, 10, 4, 1, 3, 5, 12, 13)
mut_hq13_sub = mut_hq13_sub[row_order,]

# matrix with VAF data (not binary)
mat_vaf = as.matrix(mut_hq13_sub[,2:17])
rownames(mat_vaf) =  mut_hq13_sub[,1] %>% unlist()  
colnames(mat_vaf) = c('PD63383aq', 'PD63383ap',
                      'PD62341ae', 'PD62341ag', 'PD62341aj', 'PD62341am', 'PD62341ap', 'PD62341b', 'PD62341u', 'PD62341ak',
                      'PD62341h', 'PD62341q', 'PD62341n', 'PD62341aa', 'PD62341ad', 'PD62341v')

# annotation of column names
col_annotation = data.frame(Status = c(rep('tumour', 10), rep('normal', 6)))
rownames(col_annotation) = colnames(mat_vaf)
annotation_colors = list(
  Status = c(normal=col_normal, tumour=col_tumour))

pdf('Results/20241020_p2_heatmap_status_13mut_vaf.pdf')
pheatmap(mat_vaf,
         cellwidth=14, cellheight=14,
         annotation_col = col_annotation,
         annotation_colors = annotation_colors,
         main="Status of shared mutations", 
         legend = T,
         cluster_rows=F, cluster_cols = F, 
         show_rownames = T, show_colnames = T,
         fontsize=11, cexCol=6) 
dev.off()

# matrix with binary data (VAF > 0.1 - present, else - absent)
mut_hq13_sub_binary_vaf01 = data.table(mut_hq13_sub)

for (j in names(mut_hq13_sub_binary_vaf01)[2:17]){
  
  rows_to_0 = which(mut_hq13_sub_binary_vaf01[[j]] < 0.1)
  set(mut_hq13_sub_binary_vaf01, rows_to_0, j, 0)
  
  rows_to_1 = which(mut_hq13_sub_binary_vaf01[[j]] >= 0.1)
  set(mut_hq13_sub_binary_vaf01, rows_to_1, j, 1)
  
}

mat_vaf01 = as.matrix(mut_hq13_sub_binary_vaf01[,2:17])
rownames(mat_vaf01) =  mut_hq13_sub_binary_vaf01[,1] %>% unlist()  
colnames(mat_vaf01) = c('PD63383aq', 'PD63383ap',
                        'PD62341ae', 'PD62341ag', 'PD62341aj', 'PD62341am', 'PD62341ap', 'PD62341b', 'PD62341u', 'PD62341ak',
                        'PD62341h', 'PD62341q', 'PD62341n', 'PD62341aa', 'PD62341ad', 'PD62341v')

# annotation of column names
col_annotation = data.frame(Status = c(rep('tumour', 10), rep('normal', 6)))
rownames(col_annotation) = colnames(mat_vaf01)
annotation_colors = list(
  Status = c(normal=col_normal, tumour=col_tumour))

pdf('Results/20241020_p2_heatmap_status_13mut_vaf_binary.pdf')
pheatmap(mat_vaf01,
         cellwidth=14, cellheight=14,
         color = c('#597ad3', '#e76f04'),
         annotation_col = col_annotation,
         annotation_colors = annotation_colors,
         main="Status of shared mutations (VAF > 0.1)", 
         legend = F, legend_breaks = c(-0.5, 0.5, 1),
         cluster_rows=F, cluster_cols=F, 
         show_rownames = T, show_colnames = T,
         fontsize=11, cexCol=6) 
dev.off()

# repeat for MTR >= 4
mut_hq13_mtr = twins_mtr[mut_ID %in% muts_hq, c('mut_ID', samples_tumour_mtr, samples_normal_PD62341_mtr), with=FALSE]
setcolorder(mut_hq13_mtr, c('mut_ID', 'PD63383aq_MTR', 'PD63383ap_MTR',
                            'PD62341ae_MTR', 'PD62341ag_MTR', 'PD62341aj_MTR', 'PD62341am_MTR', 'PD62341ap_MTR', 'PD62341b_MTR', 'PD62341u_MTR', 'PD62341ak_MTR',
                            'PD62341h_MTR', 'PD62341q_MTR', 'PD62341n_MTR', 'PD62341aa_MTR', 'PD62341ad_MTR', 'PD62341v_MTR'))
row_order = c(2, 8, 9, 11, 10, 4, 1, 3, 5, 12, 13)
mut_hq13_mtr = mut_hq13_mtr[row_order,]

# matrix with VAF data (not binary)
mat_mtr = as.matrix(mut_hq13_mtr[,2:17])
rownames(mat_mtr) =  mut_hq13_mtr[,1] %>% unlist()  
colnames(mat_mtr) = c('PD63383aq', 'PD63383ap',
                      'PD62341ae', 'PD62341ag', 'PD62341aj', 'PD62341am', 'PD62341ap', 'PD62341b', 'PD62341u', 'PD62341ak',
                      'PD62341h', 'PD62341q', 'PD62341n', 'PD62341aa', 'PD62341ad', 'PD62341v')

# annotation of column names
col_annotation = data.frame(Status = c(rep('tumour', 10), rep('normal', 6)))
rownames(col_annotation) = colnames(mat_mtr)
annotation_colors = list(
  Status = c(normal=col_normal, tumour=col_tumour))

pdf('Results/20241020_p2_heatmap_status_13mut_mtr.pdf')
pheatmap(mat_mtr,
         cellwidth=14, cellheight=14,
         annotation_col = col_annotation,
         annotation_colors = annotation_colors,
         main="Status of shared mutations", 
         legend = T,
         cluster_rows=F, cluster_cols = F, 
         show_rownames = T, show_colnames = T,
         fontsize=11, cexCol=6) 
dev.off()

# matrix with binary data (MTR >= 4)
mut_hq13_binary_mtr4 = data.table(mut_hq13_mtr)

for (j in names(mut_hq13_binary_mtr4)[2:17]){
  
  rows_to_0 = which(mut_hq13_binary_mtr4[[j]] < 4)
  set(mut_hq13_binary_mtr4, rows_to_0, j, 0)
  
  rows_to_1 = which(mut_hq13_binary_mtr4[[j]] >= 4)
  set(mut_hq13_binary_mtr4, rows_to_1, j, 1)
  
}

mat_mtr4 = as.matrix(mut_hq13_binary_mtr4[,2:17])
rownames(mat_mtr4) =  mut_hq13_binary_mtr4[,1] %>% unlist()  
colnames(mat_mtr4) = c('PD63383aq', 'PD63383ap',
                       'PD62341ae', 'PD62341ag', 'PD62341aj', 'PD62341am', 'PD62341ap', 'PD62341b', 'PD62341u', 'PD62341ak',
                       'PD62341h', 'PD62341q', 'PD62341n', 'PD62341aa', 'PD62341ad', 'PD62341v')

# annotation of column names
col_annotation = data.frame(Status = c(rep('tumour', 10), rep('normal', 6)))
rownames(col_annotation) = colnames(mat_mtr4)
annotation_colors = list(
  Status = c(normal=col_normal, tumour=col_tumour))

pdf('Results/20241020_p2_heatmap_status_13mut_mtr_binary.pdf')
pheatmap(mat_mtr4,
         cellwidth=14, cellheight=14,
         color = c('#597ad3', '#e76f04'),
         annotation_col = col_annotation,
         annotation_colors = annotation_colors,
         main="Status of shared mutations (MTR >= 4)", 
         legend = F, legend_breaks = c(-0.5, 0.5, 1),
         cluster_rows=F, cluster_cols=F, 
         show_rownames = T, show_colnames = T,
         fontsize=11, cexCol=6) 
dev.off()

######################################################################################################
# select mutations which are only present in one twin (normal samples)
twins_filtered_mtr[sum_normal_PD62341 >= 4 & sum_normal_PD63383 <= 1] # 15
# "chr14_105458006_C_A" # looks okay 
# "chr16_5479739_C_T" # looks kind of okay, but suspicious nr of insertions
# "chr16_88009874_G_A" # looks kind of okay, some insertions but not many  
# "chr1_148424284_G_A" # looks good (not called)
# "chr20_44114996_C_T" # looks okay 
# "chr20_54999342_C_A" # looks okay
# "chr2_95662131_G_A" # I'd double check, some reads poor quality but fine otherwise  
# "chr3_50106043_C_T" # looks real 
# "chr3_62055057_C_G" # looks real 
# "chr3_62055077_G_C" # looks real 
# "chrX_106521135_C_A" # not entirely convinced but double check 

twins_filtered_mtr[sum_normal_PD62341 <= 1 & sum_normal_PD63383 >= 4, mut_ID] # 12
# "chr11_34011887_C_T" # looks real   
# "chr13_50815806_A_G" # looks okay but very low nr of MTR
# "chr1_196748545_G_A" # poor mapping quality
# "chr1_242992933_A_C" # poor mapping quality 
# "chr21_28681933_A_C"  # poor mapping quality 
# "chr21_40193588_G_A" # looks okay 
# "chr3_77633967_C_T" # looks real   
# "chr4_74625500_G_T" # looks okay  
# "chr4_75704880_G_A" # looks quite okay 
# "chr6_165179306_G_A" # looks real 
# "chr7_73831920_C_T" # looks real  
# "chrX_115066661_C_T" # looks real 

######################################################################################################
######################################################################################################
# MTR, DEP, VAF DATA TUMOUR (2024/10/29)

twins_filtered_dt[, sum_tumour_mtr := rowSums(.SD >= 4), .SDcols = samples_tumour_mtr]
twins_filtered_dt[, sum_tumour_vaf := rowSums(.SD >= 0.1), .SDcols = samples_tumour_vaf]

twins_tumour_dt = twins_filtered_dt[sum_tumour_mtr>0 & sum_tumour_vaf>0]
paste('Number of mutations present in the tumour:', dim(twins_tumour_dt)[1])

# subset MTR, DEP and VAF TUMOUR data
twins_tumour_mtr = twins_tumour_dt[,c('mut_ID', samples_tumour_mtr, 'sum_tumour_mtr', 'sum_tumour_vaf'), with=FALSE]
twins_tumour_dep = twins_tumour_dt[,c('mut_ID', samples_tumour_dep, 'sum_tumour_mtr', 'sum_tumour_vaf'), with=FALSE]
twins_tumour_vaf = twins_tumour_dt[,c('mut_ID', samples_tumour_vaf, 'sum_tumour_mtr', 'sum_tumour_vaf'), with=FALSE]

# determine min, max and mean nr of variant reads / depth / vaf (in all tumour samples)
twins_tumour_mtr[,mean_mtr := apply(.SD, 1, function(x) mean(x[x>0])), .SDcols = samples_tumour_mtr]
twins_tumour_mtr[,median_mtr := apply(.SD, 1, function(x) median(x[x>0])), .SDcols = samples_tumour_mtr]
twins_tumour_mtr[,min_mtr := apply(.SD, 1, function(x) min(x[x>0])), .SDcols = samples_tumour_mtr]
twins_tumour_mtr[,max_mtr := apply(.SD, 1, max), .SDcols = samples_tumour_mtr]

twins_tumour_dep[,mean_dep := apply(.SD, 1, function(x) mean(x[x>0])), .SDcols = samples_tumour_dep]
twins_tumour_dep[,median_dep := apply(.SD, 1, function(x) median(x[x>0])), .SDcols = samples_tumour_dep]
twins_tumour_dep[,min_dep := apply(.SD, 1, function(x) min(x[x>0])), .SDcols = samples_tumour_dep]
twins_tumour_dep[,max_dep := apply(.SD, 1, max), .SDcols = samples_tumour_dep]

twins_tumour_vaf[,mean_vaf := apply(.SD, 1, function(x) mean(x[x>0])), .SDcols = samples_tumour_vaf]
twins_tumour_vaf[,median_vaf := apply(.SD, 1, function(x) median(x[x>0])), .SDcols = samples_tumour_vaf]
twins_tumour_vaf[,min_vaf := apply(.SD, 1, function(x) min(x[x>0])), .SDcols = samples_tumour_vaf]
twins_tumour_vaf[,max_vaf := apply(.SD, 1, max), .SDcols = samples_tumour_vaf]

# Add presence / absence based on MTR data
twins_tumour_mtr[, sum_tumour := rowSums(.SD>=4), .SDcols = samples_tumour_mtr]
twins_tumour_mtr[, sum_tumour_PD62341 := rowSums(.SD>=4), .SDcols = samples_tumour_PD62341_mtr]
twins_tumour_mtr[, sum_tumour_PD63383 := rowSums(.SD>=4), .SDcols = samples_tumour_PD63383_mtr]

# Add presence / absence based on VAF data
twins_tumour_vaf[, sum_tumour := rowSums(.SD>=0.1), .SDcols = samples_tumour_vaf]
twins_tumour_vaf[, sum_tumour_PD62341 := rowSums(.SD>=0.1), .SDcols = samples_tumour_PD62341_vaf]
twins_tumour_vaf[, sum_tumour_PD63383 := rowSums(.SD>=0.1), .SDcols = samples_tumour_PD63383_vaf]

######################################################################################################
# Plot a histogram of all mutations across all tumour samples

# Mutations present in at least 1 tumour sample
pdf('Results/20241030_p1_hist_tumour_mtr.pdf')
hist(twins_tumour_mtr[, 2:11, with=FALSE] %>% unlist(), breaks=100, xlab = 'Number of variant reads', main = 'Distribution of variant reads in tumour samples')
dev.off()

pdf('Results/20241030_p1_hist_tumour_vaf.pdf')
hist(twins_tumour_vaf[, 2:11, with=FALSE] %>% unlist(), breaks=100, xlab = 'Variant allele frequency', main = 'Distribution of VAF in tumour samples')
dev.off()

pdf('Results/20241030_p1_hist_tumour_dep.pdf')
hist(twins_tumour_dep[, 2:11, with=FALSE] %>% unlist(), breaks=100, xlab = 'Total number of reads', main = 'Distribution of coverage in tumour samples')
dev.off()

# plot mean VAF
pdf('Results/20241030_p1_hist_tumour_vaf_mean.pdf')
hist(twins_tumour_vaf[,mean_vaf] %>% unlist(), breaks=100, xlab = 'Mean variant allele frequency', main = 'Distribution of mean VAF in tumour samples')
dev.off()

pdf('Results/20241030_p1_hist_tumour_vaf_median.pdf')
hist(twins_tumour_vaf[,median_vaf] %>% unlist(), breaks=100, xlab = 'Median variant allele frequency', main = 'Distribution of mean VAF in tumour samples')
dev.off()

# look only at mutations detected in all tumour samples
twins_tumour_mtr_all = twins_tumour_mtr[sum_tumour_mtr==10 & sum_tumour_vaf==10]
twins_tumour_vaf_all = twins_tumour_vaf[sum_tumour_mtr==10 & sum_tumour_vaf==10]

pdf('Results/20241030_p1_hist_tumour_mtr_mean_all.pdf')
hist(twins_tumour_mtr_all[,mean_mtr] %>% unlist(), breaks=100, xlab = 'Number of variant reads', main = 'Distribution of mean MTR in tumour samples')
dev.off()

pdf('Results/20241030_p1_hist_tumour_mtr_mean_all.pdf')
hist(twins_tumour_vaf_all[,mean_vaf] %>% unlist(), breaks=100, xlab = 'Variant allele frequency', main = 'Distribution of mean VAF in tumour samples')
dev.off()

pdf('Results/20241030_p1_hist_tumour_mtr_median_all.pdf')
hist(twins_tumour_vaf_all[,median_vaf] %>% unlist(), breaks=100, xlab = 'Variant allele frequency', main = 'Distribution of mean VAF in tumour samples')
dev.off()

median(twins_tumour_vaf_all[,median_vaf] %>% unlist()) # 0.3407
mean(twins_tumour_vaf_all[,median_vaf] %>% unlist()) # 0.3435

# run this for each sample separately to estimate contamination per-sample
for (sample in samples_tumour){
  
  s_mtr = paste0(sample, '_MTR')
  s_vaf = paste0(sample, '_VAF')
  
  sample_mtr = twins_tumour_mtr[, ..s_mtr] %>% unlist()
  sample_vaf = twins_tumour_vaf[, ..s_vaf] %>% unlist()
  sample_mtr_all = twins_tumour_mtr_all[, ..s_mtr] %>% unlist()
  sample_vaf_all = twins_tumour_vaf_all[, ..s_vaf] %>% unlist()
  
  pdf(glue('Results/20241030_p1_hist_tumour_mtr_{sample}.pdf'))
  hist(sample_mtr, xlab = 'Number of variant reads', main = glue('Distribution of variant reads in tumour sample: {sample}'))
  dev.off()
  pdf(glue('Results/20241030_p1_hist_tumour_vaf_{sample}.pdf'))
  hist(sample_vaf, breaks=100, xlab = 'Variant allele frequency', main = glue('Distribution of VAF in tumour sample: {sample}'))
  dev.off()
  pdf(glue('Results/20241030_p1_hist_tumour_mtr_all_{sample}.pdf'))
  hist(sample_mtr_all, xlab = 'Number of variant reads', main = glue('Distribution of variant reads in tumour sample: {sample}'))
  dev.off()
  pdf(glue('Results/20241030_p1_hist_tumour_vaf_all_{sample}.pdf'))
  hist(sample_vaf_all, breaks=100, xlab = 'Variant allele frequency', main = glue('Distribution of VAF in tumour sample: {sample}'))
  dev.off()
  
}

######################################################################################################
# Find and check mutations at specific VAF values
twins_tumour_vaf[mean_vaf > 0.95]


######################################################################################################
# plots modeled on 2015 PJC paper 

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

# OLD MUTATIONS 
# "chr7_152562130_G_A" # looks okay, missing in PD62341ak, PD62341ag, PD62341u



