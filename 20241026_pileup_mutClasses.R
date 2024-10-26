# Script to analyse the pileup (high quality, run 15/10/2024)
# 2024-10-25
# Barbara Walkowiak bw18

# INPUT: 
# 1 merged pileup dataframes with mutation calls (from CaVEMan) for tumour and normal samples
# 2 list of mutations that passed required filters 

# Cleaned up script to identify mutations of interest 

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
twins_dt = data.table(read.csv('Data/pileup_merged_20241016.tsv'))

# Drop columns with PD38is_wgs (used as reference)
twins_PDv38is = grep("PDv38is", names(twins_dt), value = TRUE)
twins_dt[, c(twins_PDv38is) := NULL]

# Filter to only include mutations retained for filtering 
muts = read.table('Data/mutations_include_20241024.txt') %>% unlist()
paste('Number of mutations that passed required filters:', length(muts)) # 1,966
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
# MUTATION CLASSES 

######################################################################################################
# EARLY DEVELOPMENT IN MONOZYGOTIC TWINS (normal samples)

######################################################################################################
# TUMOUR EVOLUTION

# Tumour-specific mutations 

# I will allow up to 3 normal samples due to the possibility of contamination
# I suspect that skin samples (aa, bb) and liver sample (PD62341h) are contaminated with reads from the tumour

# Mutations found in all tumour samples 
twins_filtered_mtr[sum_tumour == 10 & sum_normal == 0] # 1 
# "chr2_199565129_G_A" # looks okay

twins_filtered_mtr[sum_tumour == 10 & sum_normal == 1] # 4 
# "chr16_32621521_C_A" # mapping quality not excellent
# "chr17_78256883_C_T" # looks okay
# "chr22_17774668_G_A" # looks okay 
# "chr22_30553442_G_T" # looks okay

twins_filtered_mtr[sum_tumour == 10 & sum_normal == 2]
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

twins_filtered_mtr[sum_tumour == 10 & sum_normal == 3]
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
# "chrX_115495236_T_C" # looks okay to me (double check but I would believe it), missing in PD62341b

twins_filtered_mtr[sum_tumour == 9 & sum_normal == 1]
# chr12_61020950_C_A # looks okay
# chr14_104122986_G_T # looks okay (nb some mis-mapping)
# chr19_3961910_G_A # looks okay
# chr7_115240063_T_C # looks okay
# chr7_120677593_C_T # looks okay
# chr8_46678612_T_C # looks okay
# chrX_46375010_C_T # looks okay
# chrX_56448796_C_A # looks okay

twins_filtered_mtr[sum_tumour == 9 & sum_normal == 2]
# "chr12_10337296_G_T" # variable mapping - ok verified this is poor quality
# "chr12_39248376_A_T" # looks okay
# "chr18_71485446_G_A" # looks like del of the other copy - can estimate contamination
# "chr21_45729938_T_A" # looks okay
# "chr22_16694683_C_A" # looks okay
# "chr2_119170376_C_T" # really not sure, double check 
# "chr2_187893783_G_T" # looks okay
# "chr3_94699090_G_T" # looks okay
# "chr4_180479634_G_A" # looks okay
# "chr5_112010886_G_A" # looks okay
# "chr7_12297469_T_G" # looks okay
# "chr8_54112201_C_A" # looks okay
# "chr9_94529993_C_T" # quite a lot of mis-mapped reads 

# there are two mutations next to each other that show a very similar pattern and VAF 
# the other is: chr8_54112197_T_C, and is called as present in 9 tumour samples and 3 normal
# normal chr8_54112197_T_C: PD62341aa, h, PD63383 bb
# normal chr8_54112201_C_A: PD62341h, PD63383bb (3 reads in PD62341aa)

twins_filtered_mtr[sum_tumour == 9 & sum_normal == 3]
# "chr10_6288921_C_T" # poor mapping 
# "chr13_100405282_G_T" # looks okay
# "chr15_76646556_T_C" # looks okay  
# "chr17_5029328_G_C" # looks okay 
# "chr1_119005653_G_C" # mis-mapped reads  
# "chr1_60839899_G_A" # looks great
# very interesting mutation chr1_60839875_T_C 
# seems to be on the other chromosome 
# can compare reads with different muts to determine contamination in different tumour samples
# also, PD62341q for example has 5 of those reads that belong to the tumour chromosome 
# we don't know tho if the loss of chr1 is subclonal in some tumour samples - what about this?

# "chr2_124573116_C_T" # looks great  
# "chr2_1420635_A_T" # looks great 
# there is a beautiful mutation next to this one 
# chr2_1420659_C_A, present on the same chr in all samples

# "chr2_2023708_C_G"  # looks great  
# "chr2_232306029_G_C" # looks great
# "chr3_49878199_C_A" # looks great  
# "chr4_46067117_A_G" # difficult region but looks okay
# "chr4_93081874_T_A" # looks great  
# "chr5_119778240_G_C" # poor mapping quality 
# "chr7_139659050_G_A" # looks okay
# "chr7_25564729_G_A" # looks okay
# "chr8_131571989_A_T" # looks okay
# "chr8_54112197_T_C" # looks okay
# "chr8_97785584_T_C" # looks okay 
# "chr9_12824689_G_A" # looks okay 
# "chrX_124531709_C_T" # looks okay  
# "chrX_83934262_C_G" # looks okay

twins_filtered_mtr[sum_tumour == 8 & sum_normal == 0] # empty  

twins_filtered_mtr[sum_tumour == 8 & sum_normal == 1] # 1
# "chr2_57814739_A_C" # mixed feelings - looks okay but quite a lot of insertions, double check

twins_filtered_mtr[sum_tumour == 8 & sum_normal == 2]
# chr2_78054577_G_A # looks okay
# chr3_153319439_A_G # looks okay

twins_filtered_mtr[sum_tumour == 7 & sum_normal == 0] # 2  
# "chr14_45928737_G_T" # poor mapping
# "chr7_152562130_G_A" # looks okay, missing in PD62341ak, PD62341ag, PD62341u

twins_filtered_mtr[sum_tumour == 7 & sum_normal == 1] # 2
# chr12_84029949_C_A # looks okay
# chr6_77631804_C_A # looks okay

twins_filtered_mtr[sum_tumour == 7 & sum_normal == 2]
# chr2_158600954_G_A # looks okay
# chr7_143823602_A_G # poor mapping 

twins_filtered_mtr[sum_tumour == 7 & sum_normal == 3]
# chr12_131995737_T_C # poor mapping 
# chr12_87682552_T_C # poor mapping
# chr16_20542569_A_G # poor mapping
# chr19_88999_T_C # poor mapping
# chr2_24681477_T_A # poor mapping
# chr8_124634298_G_T # poor mapping
# chr8_98618014_G_A # poor mapping

# Not sure the below are worth it since we don't know if they tell us anything about the tumour phylogeny or not 
twins_filtered_mtr[sum_tumour == 6 & sum_normal == 0] # 1 
# chr14_45928734_G_A # poor mapping 

twins_filtered_mtr[sum_tumour == 6 & sum_normal == 1] # 5
# chr10_51734221_A_T # looks okay
# chr12_112652964_C_T # looks okay
# chr15_94939486_G_A # looks okay
# chr6_58221412_T_C # mutations on poorly mapped reads
# chr8_87226547_A_G # poor mapping 

twins_filtered_mtr[sum_tumour == 6 & sum_normal == 2]
# chr15_30811818_A_G # poor mapping 
# chr18_64267234_G_A # looks okay, inspect further - is this on the deleted copy?
# chr1_243074560_C_T # poor mapping
# chr2_113513495_A_G # poor mapping
# chr2_131290036_G_A # poor mapping 
# chr2_172537945_G_A # poor mapping
# chr3_195982104_C_T # poor mapping
# chr4_68662000_A_T # variable mapping, double check
# chr9_41808224_G_T # poor mapping 

twins_filtered_mtr[sum_tumour == 5 & sum_normal == 0] # 4
# chr16_17805832_C_T # looks okay 
# chr2_111603601_A_G # rather poor mapping quality 
# chr7_153784080_C_G # mutations on poorly mapped reads
# chr8_21444725_G_A # looks okay

twins_filtered_mtr[sum_tumour == 5 & sum_normal == 1] # 8
# "chr11_26060345_G_T" # looks okay
# "chr12_129105882_A_G" # poor mapping 
# "chr13_103318073_T_A" # looks okay
# "chr14_54164980_C_G" # looks okay
# "chr17_70958446_C_T" # looks okay
# "chr2_87466042_G_A" # poor mapping    
# "chr6_168731130_C_T" # looks okay  
# "chr6_8310684_G_T" # looks okay

twins_filtered_mtr[sum_tumour == 4 & sum_normal == 0] # 2
# chr13_28953051_G_T # looks okay
# chr15_28718667_A_G # poor mapping

twins_filtered_mtr[sum_tumour == 4 & sum_normal == 1] # 11
# chr12_114403258_G_A" # looks okay
# "chr15_23357707_C_T"  # poor mapping 
# "chr16_89042129_G_A"  # poor mapping 
# "chr18_1379467_G_T"  # looks okay
# "chr2_113447401_A_T" # poor mapping 
# "chr2_131534951_C_T" # poor mapping 
# "chr2_95748047_C_T" # very few reads, some mapping issues - mixed feelings
# "chr4_135461656_G_A" # looks okay
# "chr8_133213958_T_C" # deletions
# "chr8_33050650_A_C" # deletions  
# "chr9_6697738_A_T" # indels  

twins_filtered_mtr[sum_tumour == 3 & sum_normal == 0] # 16  
# "chr10_87249207_G_T" # poor mapping 
# "chr11_43029745_A_T" # looks okay
# "chr12_19555916_C_T" # looks okay
# "chr12_82832044_C_A" # looks okay
# "chr1_143987248_C_T" # poor mapping 
# "chr21_13120203_G_A" # poor mapping 
# "chr2_174270140_T_A" # poor mapping 
# "chr2_231591536_G_T" # looks okay
# "chr3_195504485_G_C" # poor mapping
# "chr5_141433408_A_T" # deletions
# "chr5_83262853_T_A" # insertions
# "chr6_91355049_C_T" # insertions
# "chr7_152796118_A_T" # variable mapping 
# "chr8_2342452_C_T" # poor mapping 
# "chr8_6177429_A_T" # poor mapping 
# "chr9_38999165_A_G" # poor mapping

twins_filtered_mtr[sum_tumour == 3 & sum_normal == 1] # 34
# "chr10_49979259_G_T" # poor mapping 
# "chr11_127782363_A_C" # poor mapping
# "chr11_4225536_T_C" # poor mapping  
# "chr11_4252756_G_T" # poor mapping 
# "chr11_67710577_G_T" # poor mapping 
# "chr12_75561694_T_C" # poor mapping 
# "chr15_30819247_G_A" # poor mapping
# "chr17_18716394_G_T" # poor mapping
# "chr18_15378233_T_C" # poor mapping
# "chr19_39186678_T_G" # poor mapping 
# "chr19_50086327_C_T" # poor mapping  
# "chr1_120476102_T_G" # poor mapping
# "chr1_145588334_T_C" # poor mapping
# "chr1_148952577_G_A" # poor mapping 
# "chr1_16649763_G_A" # poor mapping   
# "chr1_201459206_G_C" # poor mapping
# "chr1_248574696_C_G" # poor mapping
# "chr22_24260204_G_T" # poor mapping  
# "chr2_113565897_A_G" # poor mapping  
# "chr2_87257985_T_C" # poor mapping  
# "chr2_87470055_G_A" # poor mapping  
# "chr3_119005597_G_A" # poor mapping 
# "chr4_179587218_G_A" # looks real!
# "chr5_110777251_A_G" # poor mapping 
# "chr5_69476229_C_T" # poor mapping   
# "chr5_99523861_G_C" # poor mapping  
# "chr7_104559847_T_G" # poor mapping 
# "chr8_141860196_G_A" # looks okay
# "chr9_41083153_A_G" # poor mapping  
# "chr9_41948814_C_T" # poor mapping  
# "chr9_64485288_C_T" # poor mapping  
# "chrX_144081878_T_C" # poor mapping
# "chrX_97349036_G_A" # poor mapping

######################################################################################################
# Mutations specific to PD63383 tumour

twins_filtered_mtr[sum_normal == 0 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 0, mut_ID] # 14
# "chr10_100308403_C_T" # looks okay
# "chr11_80115755_C_T" # looks okay
# "chr19_43348626_C_A" # not sure, double check
# "chr1_51714096_G_C" # looks okay  
# "chr2_72939205_C_T" # looks okay 
# "chr2_82147172_T_C" # looks okay 
# "chr5_157248612_A_G" # looks okay
# "chr5_28014472_C_T" # looks okay   
# "chr5_54017866_G_T" # looks okay  
# "chr7_13677323_A_G" # mapping not excellent, double check
# "chr9_41622225_T_C" # poor mapping quality 
# "chrX_48453603_G_A" # looks believable but close to a region very prone to insertions - I think it's fine
# "chrX_72671786_A_T" # poor mapping + strand bias   
# "chrX_77681162_G_A" # looks okay

muts_PD63383_tumour1 = c("chr10_100308403_C_T", "chr11_80115755_C_T", "chr19_43348626_C_A", "chr1_51714096_G_C",
                        "chr2_72939205_C_T", "chr2_82147172_T_C", "chr5_157248612_A_G", "chr5_28014472_C_T",
                        "chr5_54017866_G_T", "chrX_48453603_G_A", "chrX_72671786_A_T", "chrX_77681162_G_A")
twins_info = twins_dt[, c('mut_ID', 'Gene', 'Transcript', 'RNA', 'CDS', 'Protein', 'Type', 'Effect'), with=FALSE]
twins_info[mut_ID %in% muts_PD63383_tumour1] # nothing in protein coding genes

twins_filtered_mtr[sum_normal_PD63383 == 1  & sum_normal_PD62341 == 0 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 0, mut_ID] # 16
# "chr10_31333522_G_A" # looks okay # bb
# "chr12_43132316_G_A" # looks okay but close to poly-T # bb
# "chr12_80762031_G_A" # looks okay # bb
# "chr1_160088990_G_A" # looks okay # bb
# "chr22_25281865_C_T" # looks okay # bb
# "chr2_130070401_C_G" # poor mapping 
# "chr3_137508691_C_A" # looks okay # bb
# "chr6_86549064_C_T" # looks okay # bb
# "chr6_92179408_A_G" # looks okay # bb
# "chr7_148446413_G_A" # looks okay # bb
# "chr7_49732059_G_A" # looks okay # bb
# note that there is a mutation next to this one 
# chr7_49732077_G_C which is present in more samples
# so can track subclonal evolution in a sense 
# good to check where chr7_49732077_G_C is present (germline?)
# I tried to check but not sure it was called anywhere 

# "chr9_98275090_C_T" # looks okay # bb
# "chrX_117418608_G_A" # looks okay # bb
# "chrX_119915655_G_A" # looks okay # bb
# "chrX_36779694_C_T" # looks okay # bb
# "chrX_70352624_G_A" # looks okay # bb

# reads in normal samples are exclusively present in bb
# PD63383bb is skin, so I find it highly plausible that this is contamination with tumour cells
# especially that the kid had lesions in the skin - can I double check with Henry what he thinks?

muts_PD63383_tumour2 = c("chr10_31333522_G_A", "chr12_43132316_G_A", "chr12_80762031_G_A", "chr1_160088990_G_A",
                            "chr22_25281865_C_T", "chr3_137508691_C_A", "chr6_86549064_C_T", "chr6_92179408_A_G", 
                            "chr7_148446413_G_A", "chr7_49732059_G_A",  "chr9_98275090_C_T", "chrX_117418608_G_A",
                            "chrX_36779694_C_T", "chrX_70352624_G_A")
twins_info[mut_ID %in% muts_PD63383_tumour2] # nothing in protein coding genes

twins_filtered_mtr[sum_normal_PD63383 == 2  & sum_normal_PD62341 == 0 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 0, mut_ID] # 1
# "chr10_73678114_G_A" poor mapping quality 

twins_filtered_mtr[sum_normal_PD63383 == 3  & sum_normal_PD62341 == 0 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 0, mut_ID] # 0
# empty 

twins_filtered_mtr[sum_normal_PD63383 == 0  & sum_normal_PD62341 == 1 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 0, mut_ID] # 0
# empty 

# how many mutations specific to PD63383 tumour have we identified?
muts_PD63383_tumour = c(muts_PD63383_tumour1, muts_PD63383_tumour2)
paste('Number of mutations private to PD63383 tumour:', length(muts_PD63383_tumour)) # 26 

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
# "chr16_32621521_C_A" # mapping quality not excellent
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

twins_filtered_mtr[sum_tumour == 9 & sum_normal_PD62341 == 2 & sum_normal_PD63383 == 0] # 3
# chr18_71485446_G_A # looks okay and on the segment lost in the other chromosome 
# chr21_45729938_T_A # looks okay
# chr2_119170376_C_T # poor quality (strand bias)
# all missing in ak 

twins_filtered_mtr[sum_tumour == 9 & sum_normal_PD62341 == 1 & sum_normal_PD63383 == 0] # 6
# mostly in h so all of these are likely contamination
# chr12_61020950_C_A # looks okay
# chr19_3961910_G_A # looks okay
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
twins_filtered_mtr[sum_tumour >= 6 & sum_normal_PD62341 >= 4 & sum_normal_PD63383 <= 1] # 17 
# "chr11_4314008_C_A" # mapping issues 
# "chr11_4314107_T_C" # mapping issues 
# "chr14_105458006_C_A" 
# "chr16_5479739_C_T"  
# "chr16_88009874_G_A"  # not sure 
# "chr17_18402034_C_A" # mapping issues  
# "chr17_72773680_C_T" # insertions
# "chr19_9123804_C_A" # insertions
# "chr1_148424284_G_A"  
# "chr20_44114996_C_T" # ok but check VAF 
# "chr2_14504041_A_C"   
# "chr2_95662131_G_A" # mapping issues  
# "chr3_50106043_C_T"   
# "chr3_58819009_G_A" # mapping issues   
# "chr3_62055057_C_G"   
# "chr3_62055077_G_C"  
# "chr7_40558111_A_C" # insertions 
# some do look real!

# which samples from PD62341 share these mutations (of those that look real)?
# "chr14_105458006_C_A" # v, q, aa, ad, n, h (>= 7 in all): clonal in ad?
# "chr16_5479739_C_T" # q, aa, ad, h, n (>= 12 in all, 3 in v): looks clonal pretty much
# "chr1_148424284_G_A" # v (4), q (9), h (7), n (5): not clonal
# "chr2_14504041_A_C" # v (7), q (9), aa (6), n (7) # generally low MTR
# "chr3_50106043_C_T" # q, aa, ad, h, n (>= 8 in all, 2 in v)
# "chr3_62055057_C_G"  # q, aa, h, n (>= 9, 2 in v, 1 in ad) 
# "chr3_62055077_G_C" # q, aa, h, n (>= 9 in all, 2, in v, 2 in ad)

twins_filtered_mtr[sum_tumour >= 6 & sum_normal_PD62341 <= 1 & sum_normal_PD63383 >= 4] # 7 
# "chr10_46137859_G_A" # mapping issues 
# "chr10_87249115_A_G" # mapping issues 
# "chr11_126712981_A_G" # insertions
# "chr13_18341947_A_G" # mapping issues
# "chr2_113598407_C_T" # mapping issues   
# "chr4_37131693_G_A" # insertions
# "chr8_6300151_C_A" # mapping issues
# nothing that looks real 

######################################################################################################
######################################################################################################
# MUTATIONS SHARED B/N PD62341 NORMAL AND TUMOUR IN MORE DETAIL

# are there mutations shared with only a specific normal tissue of PD62341?
twins_filtered_mtr[sum_tumour >= 6 & sum_normal_PD62341 >= 1 & sum_normal_PD63383 == 0] # 26

# "chr12_43772588_T_A" # insertions
# "chr12_61020950_C_A"  # looks okay 
# "chr14_105458006_C_A" # looks okay
# "chr16_32621521_C_A"  # poor mapping 
# "chr17_78256883_C_T" # looks okay
# "chr18_71485446_G_A" # looks okay 
# "chr19_245236_A_C" # poor mapping   
# "chr19_3961910_G_A" # looks okay  
# "chr21_45729938_T_A" # looks okay
# "chr22_17774668_G_A" # looks quite okay 
# "chr22_30553442_G_T" # generally looks okay
# "chr2_119170376_C_T" # insertions and strand bias 
# "chr2_131290036_G_A" # poor mapping quality 
# "chr2_57814739_A_C" # poor mapping quality   
# "chr3_195982104_C_T" # poor mapping quality  
# "chr3_58819009_G_A" # poor mapping quality    
# "chr3_62055057_C_G" # looks okay  
# "chr3_62055077_G_C" # looks okay
# "chr6_58221412_T_C" # muts only on reads with poor mapping   
# "chr6_95827754_C_A" # rather poor mapping quality  
# "chr7_120677593_C_T" # insertions
# "chr7_139659050_G_A" # looks okay 
# "chr7_40558111_A_C" # insertions  
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
twins_filtered_mtr[sum_normal_PD62341 >= 4 & sum_normal_PD63383 <= 1]
# "chr11_127782275_A_G" 
# "chr11_4313257_A_G"   
# "chr11_4314008_C_A"  
# "chr11_4314107_T_C"   
# "chr12_51026540_G_C"  
# "chr14_105458006_C_A" # looks okay 
# "chr15_22594356_G_A"  
# "chr15_28712470_T_C"  
# "chr16_5479739_C_T" # looks kind of okay, but suspicious nr of insertions
# "chr16_86987421_C_T" # most Ts would be due to insertions I believe 
# "chr16_88009874_G_A" # looks kind of okay, some insertions but not many  
# "chr17_18402034_C_A" 
# "chr17_72773680_C_T" # high insertion content 
# "chr19_13472747_T_A" # insertions 
# "chr19_468950_G_A"  # quite a lot of insertions  
# "chr19_9123804_C_A" # deletions 
# "chr1_148178303_A_G" # poor mapping quality  
# "chr1_148424284_G_A" # looks good 
# "chr1_202828783_A_G" # insertions 
# "chr20_44114996_C_T"  
# "chr20_54999342_C_A" # looks real 
# "chr20_54999345_T_C"  # looks real 
# "chr2_124594744_C_T" # poor mapping quality 
# "chr2_14504041_A_C" # poor mapping quality 
# "chr2_86720577_G_T" # poor quality 
# "chr2_95662131_G_A" # I'd double check, some reads poor quality but fine otherwise  
# "chr3_130995649_G_A" # insertions
# "chr3_50106043_C_T" # looks real 
# "chr3_58819009_G_A" # poor mapping quality 
# "chr3_62055057_C_G" # looks real 
# "chr3_62055077_G_C" # looks real 
# "chr5_70784531_C_A" # poor mapping quality   
# "chr5_71157247_C_T" # poor mapping quality 
# "chr7_40558111_A_C" # insertions 
# "chr7_94130613_A_G" # insertions  
# "chrX_106521135_C_A" # not entirely convinced but double check 
# "chrX_35808_T_A" # poor mapping quality 

twins_filtered_mtr[sum_normal_PD62341 <= 1 & sum_normal_PD63383 >= 4, mut_ID]
# "chr10_46137859_G_A" # poor mapping quality   
# "chr10_87249115_A_G" # poor mapping quality    
# "chr11_126712981_A_G" # insertions 
# "chr11_34011887_C_T" # looks real   
# "chr12_66882864_C_A" # insertions 
# "chr13_18341947_A_G" # poor mapping quality
# "chr13_50815806_A_G" # looks okay but very low nr of MTR
# "chr16_20481199_C_T" # poor mapping quality   
# "chr16_71358622_G_C" # insertions
# "chr17_18844573_T_C" # poor mapping quality 
# "chr19_53269593_G_C" # mapping issues  
# "chr1_1446554_A_G" # poor mapping quality  
# "chr1_196748545_G_A" # poor mapping quality
# "chr1_242992933_A_C" # poor mapping quality 
# "chr1_25399456_G_A" # poor mapping quality  
# "chr21_28681933_A_C"  # poor mapping quality 
# "chr21_40193588_G_A" # looks okay 
# "chr22_21781563_A_C" # insertions
# "chr2_113598407_C_T" # poor mapping quality 
# "chr3_77633967_C_T" # looks real   
# "chr4_37131693_G_A" # insertions 
# "chr4_74625500_G_T" # looks okay  
# "chr4_75704880_G_A" # looks quite okay 
# "chr5_23443714_T_G" # insertions  
# "chr6_165179306_G_A" # looks real 
# "chr7_152796704_A_G" # poor mapping quality   
# "chr7_73831920_C_T" # looks real  
# "chr8_6300151_C_A"  # poor mapping quality 
# "chrX_115066661_C_T" # looks real 

######################################################################################################
######################################################################################################
######################################################################################################
# Mutations that define each twin 

######################################################################################################

# private to PD63383
twins_filtered_mtr[sum_normal_PD62341 == 0 & sum_normal_PD63383 == 6] # 0
twins_filtered_mtr[sum_normal_PD62341 == 1 & sum_normal_PD63383 == 6] # 3
# all present in PD62341v (spleen?), all at VAF ~0.1?
# chr11_34011887_C_T looks okay (vaf ~0.25 in PD63383)
# chr4_75704880_G_A looks okay (vaf ~0.25 in PD63383)
# chr6_165179306_G_A looks okay (vaf ~0.25 in PD63383)

twins_filtered_mtr[sum_normal_PD62341 == 2 & sum_normal_PD63383 == 6]
# chr13_80917269_T_G # poor mapping 
# chr1_77763_G_A # poor mapping 
# chr3_105182212_T_C # poor mapping 
# chr3_165901319_C_A # poor mapping 
# chr6_55120674_A_G # not sure - checked again and I think this is rubbish 

twins_filtered_mtr[sum_normal_PD62341 == 0 & sum_normal_PD63383 == 5]
# chr13_50815806_A_G # looks okay (3 reads in PD63383u but VAF 0.15)

twins_filtered_mtr[sum_normal_PD62341 == 1 & sum_normal_PD63383 == 5]
# 1 in aa and the others in PD62341v, missing in PD63383u / PD63383ae / PD63383ak (dist)
# interestingly the 3 other mutations are absent from tumour samples: another evidence this is really PD62341 cells 
# these ones could be really cool to look at, talk to Henry about those!
# chr13_18341947_A_G # mis-mapping 
# chr17_18844573_T_C # mis-mapping 
# chr1_16849834_G_T # poor mapping 
# chr1_35202938_A_G # indels 
# chr4_74625500_G_T # oh that looks great! 
# chr7_73831920_C_T # looks okay

twins_filtered_mtr[sum_normal_PD62341 == 2 & sum_normal_PD63383 == 5]
# "chr13_45665901_G_A" # poor mapping  
# "chr14_103083068_T_C" # double check 
# "chr15_101386102_C_T" # deletions
# "chr15_22477455_A_T" # poor mapping
# "chr15_82429313_C_A" # poor mapping  
# "chr17_36270996_T_A" # poor mapping 
# "chr1_120463162_T_C" # poor mapping 
# "chr1_185892098_T_C" # poor mapping
# "chr1_223890733_G_A" # poor mapping 
# "chr1_84230759_C_T" # poor mapping  
# "chr2_24681471_A_G" # poor mapping  
# "chr2_5987933_C_T" # poor mapping 
# "chr2_72019504_G_A" # poor mapping 
# "chr3_195504100_C_T" # poor mapping  
# "chr5_69908084_C_T" # poor mapping  
# "chr7_125275065_A_C" # poor mapping
# "chr7_146018850_C_T" # poor mapping 
# "chr7_30415341_A_G" # poor mapping  
# "chr8_98742634_G_C" # poor mapping  
# "chr9_133191234_T_C" # poor mapping
# "chr9_94544944_G_C" # poor mapping  
# "chrX_135194514_G_C" # poor mapping / deletions

######################################################################################################

# private to PD62341 
twins_filtered_mtr[sum_normal_PD62341 == 6 & sum_normal_PD63383 == 0] # 1
# chr14_105458006_C_A favourite one 
# chr3_58819009_G_A poor mapping 

twins_filtered_mtr[sum_normal_PD62341 == 6 & sum_normal_PD63383 == 1] # 0
twins_filtered_mtr[sum_normal_PD62341 == 6 & sum_normal_PD63383 == 2] # 6
# chr13_19104709_A_G # poor mapping 
# chr15_49480646_T_A # poor mapping 
# chr17_33422229_C_A # looks good! 
# chr2_218846190_G_C # poor mapping
# chr3_195492439_C_T # poor mapping 
# chr9_40726116_G_A # poor mapping 

twins_filtered_mtr[sum_normal_PD62341 == 5 & sum_normal_PD63383 == 0] # 2
# chr12_51026540_G_C # deletions
# chr2_124594744_C_T # poor mapping 

twins_filtered_mtr[sum_normal_PD62341 == 5 & sum_normal_PD63383 == 1] # 7
# always missing in PD62341v and tend to be present in PD63383bb (contamination from the tumour?)
# chr16_5479739_C_T # looks real
# chr17_18402034_C_A # poor mapping
# chr19_468950_G_A # looks messy - double check with Henry
# chr1_148178303_A_G # poor mapping 
# chr2_113577037_G_A # poor mapping 
# chr2_95662131_G_A # looks okay?, double check
# chr3_50106043_C_T # looks okay?, double check 

twins_filtered_mtr[sum_normal_PD62341 == 5 & sum_normal_PD63383 == 2] # 30
# "chr10_46118239_C_A" # poor mapping 
# "chr13_38383597_T_C" # looks kind of okay? - I wouldn't include it
# "chr15_32446399_C_T" # poor mapping 
# "chr15_32446476_T_C" # poor mapping 
# "chr16_16235743_A_G" # poor mapping 
# "chr16_16237911_T_G" # poor mapping 
# "chr17_220430_C_T" # poor mapping    
# "chr1_143878377_T_C" # poor mapping 
# "chr1_144415293_C_T" # poor mapping 
# "chr1_25421823_C_T" # not convinced
# "chr1_47108559_G_A" # poor mapping  
# "chr22_10683960_T_A" # poor mapping 
# "chr22_22278002_C_G" # poor mapping 
# "chr2_111602759_A_G" # poor mapping 
# "chr2_19307044_C_A" # poor mapping 
# "chr2_231848552_A_G" # poor mapping
# "chr2_68517910_T_A" # poor mapping 
# "chr2_72068448_A_G" # poor mapping 
# "chr3_195991332_T_C" # poor mapping
# "chr5_125420017_A_T" # poor mapping
# "chr5_152309464_C_T" # poor mapping
# "chr5_167859432_G_A" # poor mapping
# "chr5_55443234_C_T" # poor mapping 
# "chr7_125701070_T_C" # looks somewhat plausible - I wouldn't include it
# "chr7_77001001_G_T" # poor mapping 
# "chr9_39073134_A_G" # poor mapping  
# "chr9_39815614_A_G" # poor mapping
# "chrX_122845887_T_C" # poor mapping
# "chrX_76327335_T_C" # poor mapping

######################################################################################################
######################################################################################################
# Can I make some nice plots with these 
muts_PD62341 = c("chr14_105458006_C_A", "chr17_33422229_C_A", "chr16_5479739_C_T",
                 "chr2_95662131_G_A", "chr3_50106043_C_T", "chr19_468950_G_A",
                 "chr13_38383597_T_C", "chr7_125701070_T_C") # note the last 2 ones I am not super convinced by 
mut_PD63383 = c("chr11_34011887_C_T", "chr4_75704880_G_A", "chr6_165179306_G_A",
                "chr13_50815806_A_G", "chr4_74625500_G_T", "chr7_73831920_C_T")
mut_early = c(mut_PD62341, mut_PD63383)

mut_PD62341_dt = twins_vaf[mut_ID %in% muts_PD62341, 1:23]
mut_PD62341_melt = melt(mut_PD62341_dt, id.vars = 'mut_ID')
mut_PD62341_melt[, sample := tstrsplit(variable, '_', fixed=TRUE, keep = 1)]
mut_PD62341_melt[, status := as.factor(fcase( 
  sample %in% samples_normal, 'normal', # differs from 0.5 so biased and maybe we don't want it
  sample %in% samples_tumour, 'tumour'
))]
mut_PD62341_melt[, twin := as.factor(fcase( 
  sample %in% samples_PD62341, 'PD62341', # differs from 0.5 so biased and maybe we don't want it
  sample %in% samples_PD63383, 'PD63383'
))]
mut_PD62341_melt[, sample_type := as.factor(paste(status, twin, sep = '_'))]
mut_PD62341_melt[, mut_ID := factor(mut_ID, levels = 
                                      c('chr16_5479739_C_T', 'chr17_33422229_C_A', 
                                        'chr14_105458006_C_A',  'chr3_50106043_C_T',
                                        'chr2_95662131_G_A', 'chr7_125701070_T_C',
                                        'chr13_38383597_T_C', 'chr19_468950_G_A'))]

# plot option 1: we can make VAF plots for different samples
for (sample_name in mut_PD62341_melt[,sample] %>% unlist() %>% unique()){
  ggplot(mut_PD62341_melt[sample==sample_name], aes(x=mut_ID, y=value))+
    geom_line(group = 1)+
    geom_point(size=1.5)+
    geom_area(aes(fill = mut_ID))+
    scale_fill_manual(values = col_normal_PD62341)+
    theme_bw(base_size = 10)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    labs(x = 'Mutation', y = 'VAF')+
    ggtitle(glue('VAF distribution in sample {sample_name}'))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    ylim(c(0, 0.7))
  ggsave(glue('Results/20241024_p3_vaf_dist_PD62341_mut8_{sample_name}.pdf'), width=9, height=7)
}

# can we aggregate samples by sample type and plot the mean VAF (+ standard deviation)
# we should use a different color for each category 
mut_PD62341_agg = data.table(aggregate(value~mut_ID+sample_type, mut_PD62341_melt, FUN=mean))
mut_PD62341_agg[, mut_ID := factor(mut_ID, levels = 
                                     c('chr16_5479739_C_T', 'chr17_33422229_C_A', 
                                       'chr14_105458006_C_A',  'chr3_50106043_C_T',
                                       'chr2_95662131_G_A', 'chr7_125701070_T_C',
                                       'chr13_38383597_T_C', 'chr19_468950_G_A'))]

ggplot(mut_PD62341_agg, aes(x=mut_ID, y=value, col = sample_type))+
  geom_line(group = 1)+
  geom_point(size=1.5)+
  geom_area(aes(fill = mut_ID))+
  scale_color_manual(values = c(col_normal_PD62341, col_normal_PD63383, col_tumour_PD62341, col_tumour_PD63383))+
  theme_bw(base_size = 12)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = 'Mutation', y = 'VAF')+
  facet_wrap(~sample_type)+
  theme(legend.position = "none")+
  ggtitle(glue('VAF distribution for twin-specific mutations'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylim(c(0, 0.7))
ggsave(glue('Results/20241024_p3_vaf_dist_PD62341_mut8_aggregated.pdf'), width=9, height=7)

# what if I do this without 
samples_contaminated = c('PD62341v', 'PD62341aa', 'PD62341h', 'PD63383bb')
mut_PD62341_agg = data.table(aggregate(value~mut_ID+sample_type, mut_PD62341_melt[!sample %in% samples_contaminated], FUN=mean))
mut_PD62341_agg[, mut_ID := factor(mut_ID, levels = 
                                     c('chr16_5479739_C_T', 'chr17_33422229_C_A', 
                                       'chr14_105458006_C_A',  'chr3_50106043_C_T',
                                       'chr2_95662131_G_A', 'chr7_125701070_T_C',
                                       'chr13_38383597_T_C', 'chr19_468950_G_A'))]

ggplot(mut_PD62341_agg, aes(x=mut_ID, y=value, col = sample_type))+
  geom_line(group = 1)+
  geom_point(size=1.5)+
  geom_area(aes(fill = mut_ID))+
  scale_color_manual(values = c(col_normal_PD62341, col_normal_PD63383, col_tumour_PD62341, col_tumour_PD63383))+
  theme_bw(base_size = 12)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = 'Mutation', y = 'VAF')+
  facet_wrap(~sample_type, ncol = 1)+
  theme(legend.position = "none")+
  ggtitle(glue('VAF distribution for twin-specific mutations\n(excluded contaminated samples)'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylim(c(0, 0.7))
ggsave(glue('Results/20241024_p3_vaf_dist_PD62341_mut8_aggregated_rmcontamination.pdf'), width=5, height=8)

ggplot(mut_PD62341_agg, aes(x=value, y=mut_ID, col = sample_type))+
  geom_point(size=3)+
  scale_color_manual(values = c(col_normal_PD62341, col_normal_PD63383, col_tumour_PD62341, col_tumour_PD63383))+
  theme_bw(base_size = 12)+
  theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_line( size=.1, color="black"))+
  labs(x = 'VAF', y = 'Mutation')+
  ggtitle(glue('VAF distribution for PD62341 mutations\n(excluded contaminated samples)'))
ggsave(glue('Results/20241024_p3_vaf_dist_PD62341_mut8_aggregated_rmcontamination_reorderxy.pdf'), width=6.5, height=8)

# PD63383 
muts_PD63383 = c("chr11_34011887_C_T", "chr4_75704880_G_A", "chr6_165179306_G_A",
                 "chr13_50815806_A_G", "chr4_74625500_G_T", "chr7_73831920_C_T")


mut_PD63383_dt = twins_vaf[mut_ID %in% muts_PD63383, 1:23]
mut_PD63383_melt = melt(mut_PD63383_dt, id.vars = 'mut_ID')
mut_PD63383_melt[, sample := tstrsplit(variable, '_', fixed=TRUE, keep = 1)]
mut_PD63383_melt[, status := as.factor(fcase( 
  sample %in% samples_normal, 'normal', # differs from 0.5 so biased and maybe we don't want it
  sample %in% samples_tumour, 'tumour'))]
mut_PD63383_melt[, twin := as.factor(fcase( 
  sample %in% samples_PD62341, 'PD62341', # differs from 0.5 so biased and maybe we don't want it
  sample %in% samples_PD63383, 'PD63383'
))]
mut_PD63383_melt[, sample_type := as.factor(paste(status, twin, sep = '_'))]
mut_PD63383_melt[, mut_ID := factor(mut_ID, levels = 
                                      c('chr16_5479739_C_T', 'chr17_33422229_C_A', 
                                        'chr14_105458006_C_A',  'chr3_50106043_C_T',
                                        'chr2_95662131_G_A', 'chr7_125701070_T_C',
                                        'chr13_38383597_T_C', 'chr19_468950_G_A'))]

# plot option 1: we can make VAF plots for different samples
for (sample_name in mut_PD63383_melt[,sample] %>% unlist() %>% unique()){
  ggplot(mut_PD63383_melt[sample==sample_name], aes(x=mut_ID, y=value))+
    geom_line(group = 1)+
    geom_point(size=1.5)+
    geom_area(aes(fill = mut_ID))+
    scale_fill_manual(values = col_normal_PD63383)+
    theme_bw(base_size = 10)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    labs(x = 'Mutation', y = 'VAF')+
    ggtitle(glue('VAF distribution in sample {sample_name}'))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    ylim(c(0, 0.7))
  ggsave(glue('Results/20241024_p3_vaf_dist_PD63383_mut8_{sample_name}.pdf'), width=9, height=7)
}

# can we aggregate samples by sample type and plot the mean VAF (+ standard deviation)
# we should use a different color for each category 
mut_PD63383_agg = data.table(aggregate(value~mut_ID+sample_type, mut_PD63383_melt, FUN=mean))
mut_PD63383_agg[, mut_ID := factor(mut_ID, levels = 
                                     c('chr11_34011887_C_T', 'chr4_75704880_G_A', 
                                       'chr4_74625500_G_T',  'chr6_165179306_G_A',
                                       'chr13_50815806_A_G', 'chr7_73831920_C_T'))]
mut_PD63383_agg[, sample_type := factor(sample_type, levels = 
                                          c('normal_PD63383', 'normal_PD62341',
                                            'tumour_PD63383', 'tumour_PD62341'))]

ggplot(mut_PD63383_agg, aes(x=mut_ID, y=value, col = sample_type))+
  geom_line(group = 1)+
  geom_point(size=1.5)+
  geom_area(aes(fill = mut_ID))+
  scale_color_manual(values = c(col_normal_PD62341, col_normal_PD63383, col_tumour_PD62341, col_tumour_PD63383))+
  theme_bw(base_size = 12)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = 'Mutation', y = 'VAF')+
  facet_wrap(~sample_type)+
  theme(legend.position = "none")+
  ggtitle(glue('VAF distribution for PD63383 specific mutations'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylim(c(0, 0.7))
ggsave(glue('Results/20241024_p3_vaf_dist_PD63383_mut6_aggregated.pdf'), width=9, height=7)

# what if I do this without contaminating samples
samples_contaminated = c('PD62341v', 'PD62341aa', 'PD62341h', 'PD63383bb')
mut_PD63383_agg = data.table(aggregate(value~mut_ID+sample_type, mut_PD63383_melt[!sample %in% samples_contaminated], FUN=mean))
mut_PD63383_agg[, mut_ID := factor(mut_ID, levels = 
                                     c('chr11_34011887_C_T', 'chr4_75704880_G_A', 
                                       'chr4_74625500_G_T',  'chr6_165179306_G_A',
                                       'chr13_50815806_A_G', 'chr7_73831920_C_T'))]
mut_PD63383_agg[, sample_type := factor(sample_type, levels = 
                                          c('normal_PD63383', 'normal_PD62341',
                                            'tumour_PD63383', 'tumour_PD62341'))]

ggplot(mut_PD63383_agg, aes(x=mut_ID, y=value, col = sample_type))+
  geom_line(group = 1)+
  geom_point(size=1.5)+
  geom_area(aes(fill = mut_ID))+
  scale_color_manual(values = c(col_normal_PD62341, col_normal_PD63383, col_tumour_PD62341, col_tumour_PD63383))+
  theme_bw(base_size = 12)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = 'Mutation', y = 'VAF')+
  facet_wrap(~sample_type, ncol = 1)+
  theme(legend.position = "none")+
  ggtitle(glue('VAF distribution for PD63383 specific mutations\n(excluded contaminated samples)'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylim(c(0, 0.7))
ggsave(glue('Results/20241024_p3_vaf_dist_PD63383_mut6_aggregated_rmcontamination.pdf'), width=5, height=8)

ggplot(mut_PD63383_agg, aes(x=value, y=mut_ID, col = sample_type))+
  geom_point(size=3)+
  scale_color_manual(values = c(col_normal_PD62341, col_normal_PD63383, col_tumour_PD62341, col_tumour_PD63383))+
  theme_bw(base_size = 12)+
  theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_line( size=.1, color="black"))+
  labs(x = 'VAF', y = 'Mutation')+
  ggtitle(glue('VAF distribution for PD63383 mutations\n(excluded contaminated samples)'))
ggsave(glue('Results/20241024_p3_vaf_dist_PD63383_mut8_aggregated_rmcontamination_reorderxy.pdf'), width=6.5, height=8)

# okay so these mutations are certainly real, now it would be good to get a distribution of those for each
mut_PD63383_melt[, mut_ID := factor(mut_ID, levels = 
                                      c('chr11_34011887_C_T', 'chr4_75704880_G_A', 
                                        'chr4_74625500_G_T',  'chr6_165179306_G_A',
                                        'chr13_50815806_A_G', 'chr7_73831920_C_T'))]
mut_PD63383_melt[, sample_type := factor(sample_type, levels = 
                                           c('normal_PD63383', 'normal_PD62341',
                                             'tumour_PD63383', 'tumour_PD62341'))]

ggplot(mut_PD63383_melt, aes(x=mut_ID, y=value, col = sample_type))+
  geom_point(size=2)+
  scale_color_manual(values = c(col_normal_PD63383, col_normal_PD62341, col_tumour_PD63383, col_tumour_PD62341))+
  theme_bw(base_size = 12)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = 'Mutation', y = 'VAF')+
  ggtitle(glue('VAF distribution of PD63383-restricted mutations'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylim(c(0, 0.4))
ggsave(glue('Results/20241024_p3_vaf_dist_PD63383_mut6_dist_in_samples.pdf'), width=9, height=7)

# Mutation selection part 2
# PD62341 
# I would say with the three last mutations (chr7, chr13, chr19)
# it doesn't look like they are specific to PD62341
# they are at reasonable VAF in both twins - what does this tell us?
# let's do this again for the 5 mutations we are confident about
muts_PD62341 = c("chr14_105458006_C_A", "chr17_33422229_C_A", "chr16_5479739_C_T",
                 "chr2_95662131_G_A", "chr3_50106043_C_T") # note the last 2 ones I am not super convinced by 

mut_PD62341_dt = twins_vaf[mut_ID %in% muts_PD62341, 1:23]
mut_PD62341_melt = melt(mut_PD62341_dt, id.vars = 'mut_ID')
mut_PD62341_melt[, sample := tstrsplit(variable, '_', fixed=TRUE, keep = 1)]
mut_PD62341_melt[, status := as.factor(fcase( 
  sample %in% samples_normal, 'normal', # differs from 0.5 so biased and maybe we don't want it
  sample %in% samples_tumour, 'tumour'
))]
mut_PD62341_melt[, twin := as.factor(fcase( 
  sample %in% samples_PD62341, 'PD62341', # differs from 0.5 so biased and maybe we don't want it
  sample %in% samples_PD63383, 'PD63383'
))]
mut_PD62341_melt[, sample_type := as.factor(paste(status, twin, sep = '_'))]
mut_PD62341_melt[, mut_ID := factor(mut_ID, levels = 
                                      c('chr16_5479739_C_T', 'chr17_33422229_C_A', 
                                        'chr14_105458006_C_A',  'chr3_50106043_C_T',
                                        'chr2_95662131_G_A'))]

# plot option 1: we can make VAF plots for different samples
for (sample_name in mut_PD62341_melt[,sample] %>% unlist() %>% unique()){
  ggplot(mut_PD62341_melt[sample==sample_name], aes(x=mut_ID, y=value))+
    geom_line(group = 1)+
    geom_point(size=1.5)+
    geom_area(aes(fill = mut_ID))+
    scale_fill_manual(values = col_normal_PD62341)+
    theme_bw(base_size = 10)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    labs(x = 'Mutation', y = 'VAF')+
    ggtitle(glue('VAF distribution in sample {sample_name}'))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    ylim(c(0, 0.7))
  ggsave(glue('Results/20241024_p3_vaf_dist_PD62341_mut5_{sample_name}.pdf'), width=9, height=7)
}

# can we aggregate samples by sample type and plot the mean VAF (+ standard deviation)
# we should use a different color for each category 
mut_PD62341_agg = data.table(aggregate(value~mut_ID+sample_type, mut_PD62341_melt, FUN=mean))
mut_PD62341_agg[, mut_ID := factor(mut_ID, levels = 
                                     c('chr16_5479739_C_T', 'chr17_33422229_C_A', 
                                       'chr14_105458006_C_A',  'chr3_50106043_C_T',
                                       'chr2_95662131_G_A'))]

ggplot(mut_PD62341_agg, aes(x=mut_ID, y=value, col = sample_type))+
  geom_line(group = 1)+
  geom_point(size=1.5)+
  geom_area(aes(fill = mut_ID))+
  scale_color_manual(values = c(col_normal_PD62341, col_normal_PD63383, col_tumour_PD62341, col_tumour_PD63383))+
  theme_bw(base_size = 12)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = 'Mutation', y = 'VAF')+
  facet_wrap(~sample_type)+
  theme(legend.position = "none")+
  ggtitle(glue('VAF distribution for twin-specific mutations'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylim(c(0, 0.7))
ggsave(glue('Results/20241024_p3_vaf_dist_PD62341_mut5_aggregated.pdf'), width=9, height=7)

# what if I do this without 
samples_contaminated = c('PD62341v', 'PD62341aa', 'PD62341h', 'PD63383bb')
mut_PD62341_agg = data.table(aggregate(value~mut_ID+sample_type, mut_PD62341_melt[!sample %in% samples_contaminated], FUN=mean))
mut_PD62341_agg[, mut_ID := factor(mut_ID, levels = 
                                     c('chr16_5479739_C_T', 'chr17_33422229_C_A', 
                                       'chr14_105458006_C_A',  'chr3_50106043_C_T',
                                       'chr2_95662131_G_A'))]

ggplot(mut_PD62341_agg, aes(x=mut_ID, y=value, col = sample_type))+
  geom_line(group = 1)+
  geom_point(size=1.5)+
  geom_area(aes(fill = mut_ID))+
  scale_color_manual(values = c(col_normal_PD62341, col_normal_PD63383, col_tumour_PD62341, col_tumour_PD63383))+
  theme_bw(base_size = 12)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = 'Mutation', y = 'VAF')+
  facet_wrap(~sample_type, ncol = 1)+
  theme(legend.position = "none")+
  ggtitle(glue('VAF distribution for twin-specific mutations\n(excluded contaminated samples)'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylim(c(0, 0.7))
ggsave(glue('Results/20241024_p3_vaf_dist_PD62341_mut5_aggregated_rmcontamination.pdf'), width=5, height=8)

ggplot(mut_PD62341_agg, aes(x=value, y=mut_ID, col = sample_type))+
  geom_point(size=3)+
  scale_color_manual(values = c(col_normal_PD62341, col_normal_PD63383, col_tumour_PD62341, col_tumour_PD63383))+
  theme_bw(base_size = 12)+
  theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_line( size=.1, color="black"))+
  labs(x = 'VAF', y = 'Mutation')+
  ggtitle(glue('VAF distribution for PD62341 mutations\n(excluded contaminated samples)'))
ggsave(glue('Results/20241024_p3_vaf_dist_PD62341_mut5_aggregated_rmcontamination_reorderxy.pdf'), width=6.5, height=8)

mut_PD62341_melt = data.table(aggregate(value~mut_ID+sample_type, mut_PD62341_melt[!sample %in% samples_contaminated], FUN=mean))
mut_PD62341_melt[, mut_ID := factor(mut_ID, levels = 
                                      c('chr16_5479739_C_T', 'chr17_33422229_C_A', 
                                        'chr14_105458006_C_A',  'chr3_50106043_C_T',
                                        'chr2_95662131_G_A'))]

ggplot(mut_PD62341_melt, aes(x=mut_ID, y=value, col = sample_type))+
  geom_point(size=2)+
  scale_color_manual(values = c(col_normal_PD62341, col_normal_PD63383, col_tumour_PD62341, col_tumour_PD63383))+
  theme_bw(base_size = 12)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = 'Mutation', y = 'VAF')+
  ggtitle(glue('VAF distribution of PD62341-restricted mutations'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylim(c(0, 0.7))
ggsave(glue('Results/20241024_p3_vaf_dist_PD62341_mut5_dist_in_samples.pdf'), width=9, height=7)

ggplot(mut_PD62341_melt, aes(x=value, y=mut_ID, col = sample_type))+
  geom_point(size=2)+
  scale_color_manual(values = c(col_normal_PD62341, col_normal_PD63383, col_tumour_PD62341, col_tumour_PD63383))+
  theme_bw(base_size = 12)+
  theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_line( size=.1, color="black"))+
  labs(x = 'VAF', y = 'Mutation')+
  ggtitle(glue('VAF distribution of PD62341-restricted mutations'))+
  xlim(c(0, 0.7))
ggsave(glue('Results/20241024_p3_vaf_dist_PD62341_mut5_dist_in_samples_reorderxy.pdf'), width=9, height=7)

######################################################################################################
# TRACKING TUMOUR EVOLUTION
# look at mutations which are present in a subset of tumour samples
twins_filtered_mtr[sum_normal == 0 & sum_tumour >= 2, mut_ID] #  mutations are only informative if shared 
# 88 such mutations 

# mutations exclusive to PD63383 tumour samples (NB allowing some normal bc contamination)

# mutations exclusive to PD62341 tumour samples
# that would be evo after the PD63383 clone spread and was isolated
# that actually makes sense given PD63383 is probably subclonal from some PD62341

# no muts present in 8/10 PD62341 and none PD63383 tumour
twins_filtered_mtr[sum_normal == 0 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 8] # 0
twins_filtered_mtr[sum_normal == 1 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 8] # 0
twins_filtered_mtr[sum_normal == 2 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 8] # 0

# no muts present in 7/10 PD62341 and none PD63383 tumour
twins_filtered_mtr[sum_normal == 0 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 7] # 0
twins_filtered_mtr[sum_normal == 1 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 7] # 0
twins_filtered_mtr[sum_normal == 2 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 7] # 0

# muts present in 6/10 PD62341 tumour and none PD63383 tumour (nothing real)
twins_filtered_mtr[sum_normal == 0 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 6] # 0

twins_filtered_mtr[sum_normal == 1 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 6, mut_ID] # 2
# chr6_58221412_T_C # muts only present on badly mapping reads
# chr8_87226547_A_G # poor mapping

twins_filtered_mtr[sum_normal == 2 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 6, mut_ID] # 1
# chr4_68662000_A_T # not convincing 

# muts present in 5/10 PD62341 and none PD63383 tumour (3 real ones!)
twins_filtered_mtr[sum_normal == 0 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 5, mut_ID] # 3
# "chr13_48918872_G_T" # poor quality (deletions)
# "chr16_17805832_C_T" # looks okay 
# ae, ag, aj, b, u (interesting - all non-primary tumour so that would suggest the PD63383 emerged from primary)
# "chr7_153784080_C_G" # poor mapping

twins_filtered_mtr[sum_normal == 1 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 5, mut_ID] # 1
# "chr6_8310684_G_T" # looks okay!
# ae, ag, aj, b, u (again!) and present in PD62341h normal - likely contaminated

twins_filtered_mtr[sum_normal == 2 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 5, mut_ID] # 5
# "chr15_22512512_G_A" # poor mapping 
# "chr17_22060965_A_G" # poor mapping
# "chr19_56105967_G_A" # poor mapping
# "chr22_12193936_C_T" # poor mapping
# "chr4_1896684_T_C" # insertions   

# muts present in 4/10 PD62341 and none PD63383 tumour
twins_filtered_mtr[sum_normal == 0 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 4, mut_ID] # 2
# "chr13_28953051_G_T" # looks okay
# present in ae, ag, aj, u (3 reads in b and ap!)
# "chr5_95743104_A_T" # insertions

twins_filtered_mtr[sum_normal == 1 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 4, mut_ID] # 4
# "chr12_114403258_G_A" # looks okay
# ae, aj, b, u (2 reads in ag and aa, 1 in am and ap)
# "chr18_1379467_G_T" # looks okay 
# ae, ag, aj, u (3 reads in b, 2 in aa, 1 in v, am, also in h but likely contaminated)

# "chr15_23357707_C_T" # poor mapping
# "chr2_95748047_C_T" # not believable 

twins_filtered_mtr[sum_normal == 2 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 4, mut_ID] # 10
# "chr14_19666039_G_A" # poor mapping 
# "chr16_69734015_A_G" # insertions
# "chr1_83135148_G_C" # poor mapping
# "chr3_198134013_T_C" # poor mapping
# "chr6_84206157_C_T" # poor mapping 
# "chr6_93467795_A_T" # insertions
# "chr8_7611749_C_A" # poor mapping  
# "chr9_138176658_A_G" # poor mapping 
# "chr9_490527_T_C" # insertions 
# "chrX_114602547_A_G" # poor mapping 

# muts present in 3/10 PD62341 and none PD63383 tumour
twins_filtered_mtr[sum_normal == 0 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 3, mut_ID] # 12
# "chr10_87249207_G_T" poor mapping
# "chr11_43029745_A_T" very low MTR numbers everywhere but looks okay on Jb
# "chr12_19555916_C_T" very low MTR numbers everywhere but looks okay on Jb
# "chr20_16480818_T_A" insertions 
# "chr2_174270140_T_A" poor mapping 
# "chr2_231591536_G_T" looks good (present in aj and ae, less in u)
# "chr3_195504485_G_C" poor mapping 
# "chr5_83262853_T_A" insertions 
# "chr6_91355049_C_T" insertions
# "chr7_152796118_A_T" poor mapping 
# "chr7_41392814_T_A" insertions  
# "chr9_38999165_A_G" poor mapping

twins_filtered_mtr[sum_normal == 1 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 3, mut_ID] # 15
# "chr10_49979259_G_T" #rather poor mapping 
# "chr11_4225536_T_C"  # poor mapping 
# "chr11_67710577_G_T" # rather poor mapping
# "chr12_53568402_T_A" # insertions 
# "chr15_30819247_G_A" # insertions
# "chr17_18716394_G_T" # poor mapping  
# "chr19_39186678_T_G" # insertions 
# "chr1_16649763_G_A" # poor mapping  
# "chr2_113565897_A_G" # poor mapping 
# "chr2_87257985_T_C" # poor mapping 
# "chr2_87470055_G_A" # quite poor mapping 
# "chr5_69476229_C_T"  # looks okay (present in ap, aa, am but low levels in other samples too)
# "chr9_64485288_C_T" # poor mapping
# "chrX_144081878_T_C" # poor mapping 
# "chrX_97349036_G_A" # insertions

twins_filtered_mtr[sum_normal == 2 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 3, mut_ID] # 15
# "chr10_119182885_G_C" poor mapping 
# "chr11_4223668_G_C" poor mapping   
# "chr12_10969656_C_T" deletions
# "chr15_75272871_A_G"  poor mapping 
# "chr15_82479918_T_C" poor mapping  
# "chr16_88294968_G_T" mixed feelings about this one 
# "chr17_72783209_A_C" insertions
# "chr1_119314101_A_C" looks okay but seems to be present in all samples on Jbrowse 
# "chr2_113524448_G_T" poor mapping
# "chr3_132470491_G_C" poor mapping  
# "chr4_88172593_A_T" insertions   
# "chr5_75041702_C_A" insertions
# "chr7_152800119_T_C" poor mapping   
# "chr8_91155897_T_C" deletions   
# "chrX_49347094_T_G" poor mapping 
# really nothing shared :((( 

# muts present in 2/10 PD62341 and none PD63383 tumour
twins_filtered_mtr[sum_normal == 0 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 2, mut_ID] # 32
twins_filtered_mtr[sum_normal == 1 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 2, mut_ID] # 24
twins_filtered_mtr[sum_normal == 2 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 2, mut_ID] # 26

######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
# MUTATIONS PRIVATE TO ONE TWIN (REGARDLESS OF THEIR TUMOUR STATUS OR NOT)
# Note that given how early those twins split, tumours will likely have some of these
# This really is somewhat disappointing 

twins_filtered_mtr[sum_normal_PD62341 == 6 & sum_normal_PD63383 == 0, mut_ID]
# "chr14_105458006_C_A" # our favorite mutation!
# "chr3_58819009_G_A" # poor mapping 

twins_filtered_mtr[sum_normal_PD62341 == 5 & sum_normal_PD63383 == 0, mut_ID]
# "chr12_51026540_G_C" # deletions
# "chr2_124594744_C_T" # poor mapping 
# "chr7_40558111_A_C" # insertions  
# "chr7_94130613_A_G" # insertions

twins_filtered_mtr[sum_normal_PD62341 == 4 & sum_normal_PD63383 == 0, mut_ID]
# "chr11_4313257_A_G"  # poor mapping
# "chr19_13472747_T_A" # insertions 
# "chr20_54999342_C_A" # looks okay?
# "chr3_130995649_G_A" # insertions
# "chr3_62055057_C_G" # looks okay 
# "chr3_62055077_G_C" # looks okay

twins_filtered_mtr[sum_normal_PD62341 == 0 & sum_normal_PD63383 == 6, mut_ID] # 0
twins_filtered_mtr[sum_normal_PD62341 == 0 & sum_normal_PD63383 == 5, mut_ID] # 0
# "chr13_50815806_A_G" # looks real but quite low VAF
# "chr4_37131693_G_A" # deletions, really not useful

twins_filtered_mtr[sum_normal_PD62341 == 0 & sum_normal_PD63383 == 4, mut_ID] # 0
# "chr10_87249115_A_G" # poor mapping 
# "chr1_196748545_G_A" # poor mapping on most reads
# "chr1_242992933_A_C" # poor mapping
# "chr21_28681933_A_C" # poor mapping 
# "chr7_152796704_A_G" # poor mapping   

######################################################################################################
######################################################################################################
# HOW ABOUT WE LOOK FOR MUTATIONS ID AS GERMLINE IN ONE TWIN BUT NOT OTHER
# just look at normal samples as tumour is shared (but allow muts to be present in the tumour)

######################################################################################################
# Mutations private to PD62341
twins_dt_germlinePD62341 = twins_dt_filters[f7_likelyGermline_PD62341==1 & f7_likelyGermline_PD63383==0] # 569
twins_dt_germlinePD62341 = twins_dt_germlinePD62341[, c('mut_ID', samples_normal_mtr, samples_normal_vaf, samples_normal_dep), with=F]

# add useful columns
twins_dt_germlinePD62341[,mean_vaf_PD62341 := apply(.SD, 1, function(x) mean(x[x>0])), .SDcols = samples_normal_PD62341_vaf]
twins_dt_germlinePD62341[,median_vaf_PD62341 := apply(.SD, 1, function(x) median(x[x>0])), .SDcols = samples_normal_PD62341_vaf]
twins_dt_germlinePD62341[,min_vaf_PD62341 := apply(.SD, 1, function(x) min(x[x>0])), .SDcols = samples_normal_PD62341_vaf]
twins_dt_germlinePD62341[,max_vaf_PD62341 := apply(.SD, 1, max), .SDcols = samples_normal_PD62341_vaf]

twins_dt_germlinePD62341[,mean_vaf_PD63383 := apply(.SD, 1, function(x) mean(x[x>0])), .SDcols = samples_normal_PD63383_vaf]
twins_dt_germlinePD62341[,median_vaf_PD63383 := apply(.SD, 1, function(x) median(x[x>0])), .SDcols = samples_normal_PD63383_vaf]
twins_dt_germlinePD62341[,min_vaf_PD63383 := apply(.SD, 1, function(x) min(x[x>0])), .SDcols = samples_normal_PD63383_vaf]
twins_dt_germlinePD62341[,max_vaf_PD63383 := apply(.SD, 1, max), .SDcols = samples_normal_PD63383_vaf]

twins_dt_germlinePD62341[,mean_mtr_PD62341 := apply(.SD, 1, function(x) mean(x[x>0])), .SDcols = samples_normal_PD62341_mtr]
twins_dt_germlinePD62341[,median_mtr_PD62341 := apply(.SD, 1, function(x) median(x[x>0])), .SDcols = samples_normal_PD62341_mtr]
twins_dt_germlinePD62341[,min_mtr_PD62341 := apply(.SD, 1, function(x) min(x[x>0])), .SDcols = samples_normal_PD62341_mtr]
twins_dt_germlinePD62341[,max_mtr_PD62341 := apply(.SD, 1, max), .SDcols = samples_normal_PD62341_mtr]

twins_dt_germlinePD62341[,mean_mtr_PD63383 := apply(.SD, 1, function(x) mean(x[x>0])), .SDcols = samples_normal_PD63383_mtr]
twins_dt_germlinePD62341[,median_mtr_PD63383 := apply(.SD, 1, function(x) median(x[x>0])), .SDcols = samples_normal_PD63383_mtr]
twins_dt_germlinePD62341[,min_mtr_PD63383 := apply(.SD, 1, function(x) min(x[x>0])), .SDcols = samples_normal_PD63383_mtr]
twins_dt_germlinePD62341[,max_mtr_PD63383 := apply(.SD, 1, max), .SDcols = samples_normal_PD63383_mtr]

twins_dt_germlinePD62341[,sum_vaf_PD62341 := rowSums(.SD>=0.1), .SDcols = samples_normal_PD62341_vaf]
twins_dt_germlinePD62341[,sum_vaf_PD63383 := rowSums(.SD>=0.1), .SDcols = samples_normal_PD63383_vaf]

twins_dt_germlinePD62341[,sum_mtr_PD62341 := rowSums(.SD>=4), .SDcols = samples_normal_PD62341_mtr]
twins_dt_germlinePD62341[,sum_mtr_PD63383 := rowSums(.SD>=4), .SDcols = samples_normal_PD63383_mtr]

# can I find mutations that are absent from PD63383
twins_dt_germlinePD62341[sum_mtr_PD62341 == 6 & sum_mtr_PD63383 == 0]
# chr14 is great, chr3 maps very poorly
twins_dt_germlinePD62341[sum_mtr_PD62341 == 5 & sum_mtr_PD63383 == 0]
# chr2_124594744_C_T poor mapping
twins_dt_germlinePD62341[sum_mtr_PD62341 == 4 & sum_mtr_PD63383 == 0]
# chr8_12200512_C_T poor mapping 
twins_dt_germlinePD62341[sum_mtr_PD62341 == 3 & sum_mtr_PD63383 == 0]
# chr15_30241671_T_G poor mapping
# chr16_22554658_G_A poor mapping
# chr8_8059732_G_A poor mapping
# chr8_8159193_G_C poor mapping 
# chr9_41310691_A_G poor mapping

######################################################################################################
# Mutations private to PD63383

twins_dt_germlinePD63383 = twins_dt_filters[f6_likelyGermline_PD62341==0 & f6_likelyGermline_PD63383==1] # 765
twins_dt_germlinePD63383 = twins_dt_germlinePD63383[, c('mut_ID', samples_normal_mtr, samples_normal_vaf, samples_normal_dep), with=F]

twins_dt_germlinePD63383[,mean_vaf_PD62341 := apply(.SD, 1, function(x) mean(x[x>0])), .SDcols = samples_normal_PD62341_vaf]
twins_dt_germlinePD63383[,median_vaf_PD62341 := apply(.SD, 1, function(x) median(x[x>0])), .SDcols = samples_normal_PD62341_vaf]
twins_dt_germlinePD63383[,min_vaf_PD62341 := apply(.SD, 1, function(x) min(x[x>0])), .SDcols = samples_normal_PD62341_vaf]
twins_dt_germlinePD63383[,max_vaf_PD62341 := apply(.SD, 1, max), .SDcols = samples_normal_PD62341_vaf]

twins_dt_germlinePD63383[,mean_vaf_PD63383 := apply(.SD, 1, function(x) mean(x[x>0])), .SDcols = samples_normal_PD63383_vaf]
twins_dt_germlinePD63383[,median_vaf_PD63383 := apply(.SD, 1, function(x) median(x[x>0])), .SDcols = samples_normal_PD63383_vaf]
twins_dt_germlinePD63383[,min_vaf_PD63383 := apply(.SD, 1, function(x) min(x[x>0])), .SDcols = samples_normal_PD63383_vaf]
twins_dt_germlinePD63383[,max_vaf_PD63383 := apply(.SD, 1, max), .SDcols = samples_normal_PD63383_vaf]

twins_dt_germlinePD63383[,mean_mtr_PD62341 := apply(.SD, 1, function(x) mean(x[x>0])), .SDcols = samples_normal_PD62341_mtr]
twins_dt_germlinePD63383[,median_mtr_PD62341 := apply(.SD, 1, function(x) median(x[x>0])), .SDcols = samples_normal_PD62341_mtr]
twins_dt_germlinePD63383[,min_mtr_PD62341 := apply(.SD, 1, function(x) min(x[x>0])), .SDcols = samples_normal_PD62341_mtr]
twins_dt_germlinePD63383[,max_mtr_PD62341 := apply(.SD, 1, max), .SDcols = samples_normal_PD62341_mtr]

twins_dt_germlinePD63383[,mean_mtr_PD63383 := apply(.SD, 1, function(x) mean(x[x>0])), .SDcols = samples_normal_PD63383_mtr]
twins_dt_germlinePD63383[,median_mtr_PD63383 := apply(.SD, 1, function(x) median(x[x>0])), .SDcols = samples_normal_PD63383_mtr]
twins_dt_germlinePD63383[,min_mtr_PD63383 := apply(.SD, 1, function(x) min(x[x>0])), .SDcols = samples_normal_PD63383_mtr]
twins_dt_germlinePD63383[,max_mtr_PD63383 := apply(.SD, 1, max), .SDcols = samples_normal_PD63383_mtr]

twins_dt_germlinePD63383[,sum_vaf_PD62341 := rowSums(.SD>=0.1), .SDcols = samples_normal_PD62341_vaf]
twins_dt_germlinePD63383[,sum_vaf_PD63383 := rowSums(.SD>=0.1), .SDcols = samples_normal_PD63383_vaf]

twins_dt_germlinePD63383[,sum_mtr_PD62341 := rowSums(.SD>=4), .SDcols = samples_normal_PD62341_mtr]
twins_dt_germlinePD63383[,sum_mtr_PD63383 := rowSums(.SD>=4), .SDcols = samples_normal_PD63383_mtr]

twins_dt_germlinePD63383[sum_mtr_PD62341 == 0 & sum_mtr_PD63383 == 6] # 0
twins_dt_germlinePD63383[sum_mtr_PD62341 == 0 & sum_mtr_PD63383 == 5] # 0 
twins_dt_germlinePD63383[sum_mtr_PD62341 == 0 & sum_mtr_PD63383 == 4] # 0
twins_dt_germlinePD63383[sum_mtr_PD62341 == 0 & sum_mtr_PD63383 == 3] # 5
# chr10_87241684_G_T # poor mapping 
# chr11_121233628_C_T # deletions
# chr15_22581513_C_G # poor mapping 
# chr22_15685254_C_T # poor mapping
# chr7_94210766_T_A # deletions 

# am I inadvertently filtering out these mutations or what is going so horribly wrong?
# maybe I should run this on all mutations (search for ones not present in one twin) and see how it goes?
twins_normal_mtr = twins_dt[, c('mut_ID', samples_normal_mtr), with=F]
twins_normal_mtr[, sum_normal_PD62341 := rowSums(.SD>=4), .SDcols = samples_normal_PD62341_mtr]
twins_normal_mtr[, sum_normal_PD63383 := rowSums(.SD>=4), .SDcols = samples_normal_PD63383_mtr]

twins_normal_mtr[sum_normal_PD62341 == 6 & sum_normal_PD63383 == 0] # only chr3 and chr14
twins_normal_mtr[sum_normal_PD62341 == 5 & sum_normal_PD63383 == 0]
# none of these are real 
# "chr12_51026540_G_C" 
# "chr2_124594744_C_T" 
# "chr5_92878529_T_C"  
# "chr7_40558111_A_C" 
# "chr7_94130613_A_G"  
# "chrX_140013223_A_C"

twins_normal_mtr[sum_normal_PD62341 == 4 & sum_normal_PD63383 == 0]
# "chr11_4313257_A_G"   
# "chr12_112622862_G_A" 
# "chr12_22869396_A_T"  
# "chr14_23581034_C_T" 
# "chr15_20314435_C_T"  
# "chr15_20626012_C_T"  
# "chr15_28360258_G_A"  
# "chr15_30593218_G_A" 
# "chr17_20755253_T_C"  
# "chr19_13472747_T_A"  
# "chr1_148532730_C_G"  
# "chr1_149377844_T_G" 
# "chr1_37002676_T_A"   
# "chr20_54999342_C_A"  
# "chr20_62417660_C_T"  
# "chr3_102628981_A_C" 
# "chr3_130995649_G_A"  
# "chr3_62055057_C_G"   
# "chr3_62055077_G_C"   
# "chr5_160709961_T_G" 
# "chr6_93632726_C_T"   
# "chr8_12200512_C_T"   
# "chr9_91173137_T_C"   
# "chrX_49318267_G_A"  
# "chrX_50061924_G_A"   
# "chrX_89802106_C_T"

twins_normal_mtr[sum_normal_PD62341 == 0 & sum_normal_PD63383 == 6] # 0
twins_normal_mtr[sum_normal_PD62341 == 0 & sum_normal_PD63383 == 5] # 4
# "chr13_50815806_A_G" # looks quite real but low VAF
# "chr20_1613209_A_G"  # very odd mapping
# "chr3_188006515_A_T" # poor mapping
# "chr4_37131693_G_A" # insertions

twins_normal_mtr[sum_normal_PD62341 == 0 & sum_normal_PD63383 == 4] 
# "chr10_87249115_A_G" 
# "chr11_90003912_A_T" 
# "chr1_196748545_G_A" 
# "chr1_242992933_A_C"
# "chr1_47504134_T_C"  
# "chr1_82552867_T_C"  
# "chr21_28681933_A_C" # poor mapping but not terrible so maybe double check
# "chr7_152796704_A_G" # poor mapping but not terrible so maybe double check



