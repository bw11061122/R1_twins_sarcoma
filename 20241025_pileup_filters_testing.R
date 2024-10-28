# Script to analyse the pileup (run 15/10/2024)
# 2024-10-11 - 2024-10-14
# Barbara Walkowiak bw18

# INPUT: merged pileup dataframes with mutation calls (from CaVEMan) for tumour and normal samples
# dataframe created in the script: 20241015_pileup_checks.R (pileup run 14/10/2024-15/10/2024)

# The goal of this script is to test the specificity and sensitivity of different thresholds applied for filtering 

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
library(VGAM)

###################################################################################################################################
# Read the merged dataframe 
setwd('/Users/bw18/Desktop/1SB')
twins_dt = data.table(read.csv('Data/pileup_merged_20241016.tsv'))
twins_lq = data.table(read.csv('Data/pileup_merged_20241025.tsv'))

# Drop columns with PD38is_wgs (used as reference)
twins_PDv38is = grep("PDv38is", names(twins_dt), value = TRUE)
twins_dt[, c(twins_PDv38is) := NULL]

###################################################################################################################################
# Specify settings for plotting 
col_tumour = '#ad0505'
col_normal = '#07a7d0'
col_PD62341 = "#d37e03"
col_PD63383 = "#d4b17e"
col_tumour_PD62341 = "#980505"
col_tumour_PD63383 = "#eb6767"
col_normal_PD62341 = "#0785a5"
col_normal_PD63383 = "#70c6db"
col_bar = '#8c08d3'

######################################################################################################
# Samples 
samples = colnames(twins_dt[,c(seq(2, 330, 15))])
samples = sapply(strsplit(samples,"_"), `[`, 1)
paste('Number of samples analysed:', length(samples)) # remove the PD38is_wgs sample: not from the samples of interest

# Create lists of possible samples of interest
samples_names = c("PD62341v", "PD62341q", "PD62341aa", "PD62341ad", "PD63383w", "PD63383t", "PD63383u", "PD63383ae", "PD63383ak", "PD63383bb","PD62341ae", "PD62341ag", "PD62341aj", "PD62341ak", "PD62341am", "PD62341ap", "PD63383ap", "PD63383aq", "PD62341b", "PD62341h", "PD62341n", "PD62341u")
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
# Create relevant subsetted dataframes 
twins_mtr = twins_dt[,c('mut_ID', samples_mtr), with=FALSE]
twins_dep = twins_dt[,c('mut_ID', samples_dep), with=FALSE]
twins_vaf = twins_dt[,c('mut_ID', samples_vaf), with=FALSE]

# determine min, max and mean nr of variant reads / depth / vaf
twins_mtr[,mean_mtr := apply(.SD, 1, function(x) mean(x[x>0])), .SDcols = samples_mtr]
twins_mtr[,median_mtr := apply(.SD, 1, function(x) median(x[x>0])), .SDcols = samples_mtr]
twins_mtr[,min_mtr := apply(.SD, 1, function(x) min(x[x>0])), .SDcols = samples_mtr]
twins_mtr[,max_mtr := apply(.SD, 1, max), .SDcols = samples_mtr]

twins_dep[,mean_dep := apply(.SD, 1, function(x) mean(x[x>0])), .SDcols = samples_dep]
twins_dep[,median_dep := apply(.SD, 1, function(x) median(x[x>0])), .SDcols = samples_dep]
twins_dep[,min_dep := apply(.SD, 1, function(x) min(x[x>0])), .SDcols = samples_dep]
twins_dep[,max_dep := apply(.SD, 1, max), .SDcols = samples_dep]

twins_vaf[,mean_vaf := apply(.SD, 1, function(x) mean(x[x>0])), .SDcols = samples_vaf]
twins_vaf[,median_vaf := apply(.SD, 1, function(x) median(x[x>0])), .SDcols = samples_vaf]
twins_vaf[,min_vaf := apply(.SD, 1, function(x) min(x[x>0])), .SDcols = samples_vaf]
twins_vaf[,max_vaf := apply(.SD, 1, max), .SDcols = samples_vaf]

# determine this specifically in normal samples (in case of ploidy changes in the tumour)
twins_mtr[,mean_mtr_normal := apply(.SD, 1, function(x) mean(x[x>0])), .SDcols = samples_normal_mtr]
twins_mtr[,median_mtr_normal := apply(.SD, 1, function(x) median(x[x>0])), .SDcols = samples_normal_mtr]
twins_mtr[,min_mtr_normal := apply(.SD, 1, function(x) min(x[x>0])), .SDcols = samples_normal_mtr]
twins_mtr[,max_mtr_normal := apply(.SD, 1, max), .SDcols = samples_normal_mtr]

twins_dep[,mean_dep_normal := apply(.SD, 1, function(x) mean(x[x>0])), .SDcols = samples_normal_dep]
twins_dep[,median_dep_normal := apply(.SD, 1, function(x) median(x[x>0])), .SDcols = samples_normal_dep]
twins_dep[,min_dep_normal := apply(.SD, 1, function(x) min(x[x>0])), .SDcols = samples_normal_dep]
twins_dep[,max_dep_normal := apply(.SD, 1, max), .SDcols = samples_normal_dep]

twins_vaf[,mean_vaf_normal := apply(.SD, 1, function(x) mean(x[x>0])), .SDcols = samples_normal_vaf]
twins_vaf[,median_vaf_normal := apply(.SD, 1, function(x) median(x[x>0])), .SDcols = samples_normal_vaf]
twins_vaf[,min_vaf_normal := apply(.SD, 1, function(x) min(x[x>0])), .SDcols = samples_normal_vaf]
twins_vaf[,max_vaf_normal := apply(.SD, 1, max), .SDcols = samples_normal_vaf]

######################################################################################################
# Define list of poorly and well mutations to check specificity and sensitivity of filters (based on Jbrowse)
poorly_mapped_muts = c("chr2_111603601_A_G", "chr7_153784080_C_G", "chr10_49979259_G_T", "chr11_127782363_A_C", "chr11_4225536_T_C",  
                       "chr11_4252756_G_T", "chr11_67710577_G_T", "chr12_75561694_T_C", "chr15_30819247_G_A", "chr17_18716394_G_T", 
                       "chr18_15378233_T_C", "chr19_39186678_T_G", "chr19_50086327_C_T", "chr1_120476102_T_G", "chr1_145588334_T_C",
                       "chr1_148952577_G_A", "chr1_16649763_G_A", "chr1_201459206_G_C", "chr1_248574696_C_G", "chr22_24260204_G_T",
                       "chr2_113565897_A_G", "chr2_87257985_T_C", "chr2_87470055_G_A", "chr3_119005597_G_A", "chr9_41083153_A_G",  
                       "chr9_41948814_C_T", "chr9_64485288_C_T", "chr1_243074560_C_T", "chr2_113513495_A_G", "chr2_131290036_G_A", 
                       "chr2_172537945_G_A", "chr3_195982104_C_T", "chr10_46137859_G_A", "chr10_87249115_A_G", "chr11_126712981_A_G",
                       "chr13_18341947_A_G", "chr2_113598407_C_T", "chr4_37131693_G_A", "chr8_6300151_C_A", "chr2_131290036_G_A", 
                       "chr2_57814739_A_C", "chr3_195982104_C_T", "chr3_58819009_G_A", "chr17_18844573_T_C", "chr19_53269593_G_C",  
                       "chr1_1446554_A_G", "chr1_196748545_G_A", "chr1_242992933_A_C", "chr1_25399456_G_A", "chr21_28681933_A_C", 
                       "chr13_80917269_T_G", "chr1_77763_G_A", "chr3_105182212_T_C", "chr3_165901319_C_A", "chr13_18341947_A_G", 
                       "chr17_18844573_T_C", "chr1_16849834_G_T", "chr15_22477455_A_T", "chr15_82429313_C_A", "chr17_36270996_T_A", 
                       "chr1_120463162_T_C", "chr1_185892098_T_C", "chr1_223890733_G_A", "chr1_84230759_C_T", "chr2_24681471_A_G",  
                       "chr2_5987933_C_T", "chr2_72019504_G_A", "chr3_195504100_C_T", "chr5_69908084_C_T", "chr7_125275065_A_C", 
                       "chr7_146018850_C_T", "chr7_30415341_A_G", "chr8_98742634_G_C", "chr9_133191234_T_C", "chr9_94544944_G_C",   
                       "chr15_32446399_C_T", "chr15_32446476_T_C", "chr16_16235743_A_G", "chr16_16237911_T_G", "chr17_220430_C_T",     
                       "chr1_143878377_T_C", "chr1_144415293_C_T", "chr1_47108559_G_A", "chr22_10683960_T_A", "chr22_22278002_C_G",  
                       "chr2_111602759_A_G", "chr2_19307044_C_A", "chr2_231848552_A_G", "chr2_68517910_T_A", "chr2_72068448_A_G",  
                       "chr3_195991332_T_C", "chr5_125420017_A_T", "chr5_152309464_C_T", "chr5_167859432_G_A", "chr5_55443234_C_T",  
                       "chr15_22512512_G_A", "chr17_22060965_A_G", "chr19_56105967_G_A", "chr22_12193936_C_T", "chr10_119182885_G_C",
                       "chr14_19666039_G_A", "chr16_69734015_A_G", "chr1_83135148_G_C", "chr3_198134013_T_C", "chr6_84206157_C_T", 
                       "chr10_119182885_G_C", "chr11_4223668_G_C", "chr12_10969656_C_T" ,"chr15_75272871_A_G", "chr15_82479918_T_C",   
                       "chr17_72783209_A_C", "chr2_113524448_G_T", "chr3_132470491_G_C", "chr4_88172593_A_T", "chr5_75041702_C_A", 
                       "chr7_152800119_T_C", "chr8_91155897_T_C", "chrX_49347094_T_G", "chr15_28398565_A_G", "chr15_32482084_G_A", 
                       "chr16_14864167_G_C", "chr17_36230372_A_G", "chr3_75122630_G_T", "chr4_26808742_A_G", "chr4_55624665_T_A")  
paste("Number of poorly mapped mutations:", length(poorly_mapped_muts))

well_mapped_muts = c("chr10_72139813_C_T", "chr17_33422229_C_A", "chr15_68707110_G_A", "chr3_189240074_C_A", "chr3_95470911_T_A",
                     "chr5_58182780_T_A", "chr8_91252357_A_C", "chrX_64561555_C_T", "chr12_112652964_C_T", "chrX_77681162_G_A",
                     "chr3_62055077_G_C",  "chr12_61020950_C_A", "chr19_3961910_G_A", "chr10_51734221_A_T", "chr3_62055057_C_G",
                     "chr17_78256883_C_T", "chr22_17774668_G_A", "chr22_30553442_G_T", "chr7_115240063_T_C",  "chr6_160324993_G_A",
                     "chr7_120677593_C_T", "chr8_46678612_T_C", "chrX_46375010_C_T", "chrX_56448796_C_A", "chr13_103318073_T_A",
                     "chr14_54164980_C_G", "chr17_70958446_C_T", "chr13_60216504_A_C", "chr15_48353502_C_A", "chr15_56691722_C_T",
                     "chr16_4003400_G_A", "chr2_187893783_G_T", "chr3_94699090_G_T", "chr4_180479634_G_A", "chr5_112010886_G_A",
                     "chr7_12297469_T_G", "chr8_54112201_C_A", "chr11_134158143_C_T", "chr14_104500332_C_T", "chrX_115495236_T_C",
                     "chr17_40061856_G_C", "chr4_130207996_C_T", "chr4_147908994_T_C", "chr4_178731458_A_G",  "chrX_66719643_C_G", 
                     "chr2_133311125_G_T", "chr5_136349883_C_T", "chr5_37790015_G_A", "chr9_11364991_T_A", "chrX_68803487_C_T",
                     "chr11_34011887_C_T", "chr4_75704880_G_A", "chr6_165179306_G_A", "chr1_51714096_G_C", "chr2_72939205_C_T", 
                     "chr2_82147172_T_C", "chr5_157248612_A_G", "chr5_28014472_C_T", "chr5_54017866_G_T", "chr12_80762031_G_A",
                     "chr1_160088990_G_A", "chr22_25281865_C_T", "chr6_86549064_C_T", "chr6_92179408_A_G", "chr7_148446413_G_A",
                     "chr7_49732059_G_A", "chr9_98275090_C_T", "chrX_117418608_G_A",  "chrX_119915655_G_A", "chrX_36779694_C_T", 
                     "chrX_70352624_G_A", "chr12_39248376_A_T", "chr10_100754461_C_T", "chr13_100405282_G_T", "chr12_84029949_C_A",
                     "chr13_50815806_A_G", "chr1_119314101_A_C", "chr12_19555916_C_T", "chr2_231591536_G_T", "chr16_17805832_C_T")
paste("Number of well mapped mutations:", length(well_mapped_muts))

######################################################################################################
# Filters with arbitrary thresholds 1: coverage 

pdf('Results/20241026_p4_checking_muts_coverage.pdf')
par(mfrow = c(1,2))
hist(twins_dep[mut_ID %in% poorly_mapped_muts, 2:23] %>% unlist(), xlab = 'Coverage', main = 'Coverage of poorly mapped mutations')
hist(twins_dep[mut_ID %in% well_mapped_muts, 2:23] %>% unlist(), xlab = 'Coverage', main = 'Coverage of well mapped mutations')
dev.off()

# test filter over coverage across all samples 
paste('Poorly mapped mutations with min coverage below 15:', dim(twins_dep[mut_ID %in% poorly_mapped_muts & min_dep < 15])[1])
paste('Poorly mapped mutations with max coverage above 60:', dim(twins_dep[mut_ID %in% poorly_mapped_muts & max_dep > 60])[1]) 
paste('Poorly mapped mutations with median coverage below 20:', dim(twins_dep[mut_ID %in% poorly_mapped_muts & median_dep < 20])[1]) 
paste('Poorly mapped mutations with median coverage above 50:', dim(twins_dep[mut_ID %in% poorly_mapped_muts & median_dep > 50])[1]) 

paste('Well mapped mutations with min coverage below 15:', dim(twins_dep[mut_ID %in% well_mapped_muts & min_dep < 15])[1]) 
paste('Well mapped mutations with max coverage above 60:', dim(twins_dep[mut_ID %in% well_mapped_muts & max_dep > 60])[1]) 
paste('Well mapped mutations with median coverage below 20:', dim(twins_dep[mut_ID %in% well_mapped_muts & median_dep < 20])[1]) 
paste('Well mapped mutations with median coverage above 50:', dim(twins_dep[mut_ID %in% well_mapped_muts & median_dep > 50])[1]) 

# test filter over coverage across normal samples 
paste('Poorly mapped mutations with min coverage below 15 (normal samples):', dim(twins_dep[mut_ID %in% poorly_mapped_muts & min_dep_normal < 15])[1]) 
paste('Poorly mapped mutations with max coverage above 60 (normal samples):', dim(twins_dep[mut_ID %in% poorly_mapped_muts & max_dep_normal > 60])[1]) 
paste('Poorly mapped mutations with median coverage below 20 (normal samples):', dim(twins_dep[mut_ID %in% poorly_mapped_muts & median_dep_normal < 20])[1]) 
paste('Poorly mapped mutations with median coverage above 50 (normal samples):', dim(twins_dep[mut_ID %in% poorly_mapped_muts & median_dep_normal > 50])[1]) 

paste('Well mapped mutations with min coverage below 15 (normal samples):', dim(twins_dep[mut_ID %in% well_mapped_muts & min_dep_normal < 15])[1]) 
paste('Well mapped mutations with max coverage above 60 (normal samples):', dim(twins_dep[mut_ID %in% well_mapped_muts & max_dep_normal > 60])[1]) 
paste('Well mapped mutations with median coverage below 20 (normal samples):', dim(twins_dep[mut_ID %in% well_mapped_muts & median_dep_normal < 20])[1]) 
paste('Well mapped mutations with median coverage above 50 (normal samples):', dim(twins_dep[mut_ID %in% well_mapped_muts & median_dep_normal > 50])[1]) 

# Can we test for different thresholds
min_median = c(10, 12, 14, 16, 18, 20)
max_median = c(50, 55, 60, 65, 70, 80)

for (min in min_median){
  
  nr_poor = dim(twins_dep[mut_ID %in% poorly_mapped_muts & median_dep < min])[1]
  nr_well = dim(twins_dep[mut_ID %in% well_mapped_muts & median_dep < min])[1]
  nr_poor_normal = dim(twins_dep[mut_ID %in% poorly_mapped_muts & median_dep_normal < min])[1]
  nr_well_normal = dim(twins_dep[mut_ID %in% well_mapped_muts & median_dep_normal < min])[1]
  print(paste('Threshold:', min))
  print(paste('Poorly mapped mutations excluded:', nr_poor / length(poorly_mapped_muts)))
  print(paste('Well mapped mutations excluded:', nr_well / length(well_mapped_muts)))
  print(paste('Threshold normal:', min))
  print(paste('Poorly mapped mutations excluded:', nr_poor_normal / length(poorly_mapped_muts)))
  print(paste('Well mapped mutations excluded:', nr_well_normal / length(well_mapped_muts)))
  
}

for (max in max_median){
  
  nr_poor = dim(twins_dep[mut_ID %in% poorly_mapped_muts & median_dep > max])[1]
  nr_well = dim(twins_dep[mut_ID %in% well_mapped_muts & median_dep > max])[1]
  nr_poor_normal = dim(twins_dep[mut_ID %in% poorly_mapped_muts & median_dep_normal > max])[1]
  nr_well_normal = dim(twins_dep[mut_ID %in% well_mapped_muts & median_dep_normal > max])[1]
  print(paste('Threshold:', max))
  print(paste('Poorly mapped mutations excluded:', nr_poor / length(poorly_mapped_muts)))
  print(paste('Well mapped mutations excluded:', nr_well / length(well_mapped_muts)))
  print(paste('Threshold normal:', max))
  print(paste('Poorly mapped mutations excluded:', nr_poor_normal / length(poorly_mapped_muts)))
  print(paste('Well mapped mutations excluded:', nr_well_normal / length(well_mapped_muts)))
  
}

# you may want to use median depth of 15 reads across all samples

min_all = c(5, 8, 10, 12, 15)
max_all = c(70, 75, 80, 90, 100)

for (min in min_all){
  
  nr_poor = dim(twins_dep[mut_ID %in% poorly_mapped_muts & min_dep < min])[1]
  nr_well = dim(twins_dep[mut_ID %in% well_mapped_muts & min_dep < min])[1]
  nr_poor_normal = dim(twins_dep[mut_ID %in% poorly_mapped_muts & min_dep_normal < min])[1]
  nr_well_normal = dim(twins_dep[mut_ID %in% well_mapped_muts & min_dep_normal < min])[1]
  print(paste('Threshold:', min))
  print(paste('Poorly mapped mutations excluded:', nr_poor / length(poorly_mapped_muts)))
  print(paste('Well mapped mutations excluded:', nr_well / length(well_mapped_muts)))
  print(paste('Threshold normal:', min))
  print(paste('Poorly mapped mutations excluded:', nr_poor_normal / length(poorly_mapped_muts)))
  print(paste('Well mapped mutations excluded:', nr_well_normal / length(well_mapped_muts)))
  
}

for (max in max_all){
  
  nr_poor = dim(twins_dep[mut_ID %in% poorly_mapped_muts & max_dep > max])[1]
  nr_well = dim(twins_dep[mut_ID %in% well_mapped_muts & max_dep > max])[1]
  nr_poor_normal = dim(twins_dep[mut_ID %in% poorly_mapped_muts & max_dep_normal > max])[1]
  nr_well_normal = dim(twins_dep[mut_ID %in% well_mapped_muts & max_dep_normal > max])[1]
  print(paste('Threshold:', max))
  print(paste('Poorly mapped mutations excluded:', nr_poor / length(poorly_mapped_muts)))
  print(paste('Well mapped mutations excluded:', nr_well / length(well_mapped_muts)))
  print(paste('Threshold normal:', max))
  print(paste('Poorly mapped mutations excluded:', nr_poor_normal / length(poorly_mapped_muts)))
  print(paste('Well mapped mutations excluded:', nr_well_normal / length(well_mapped_muts)))
  
}

# I think based on this let's filter out mutations with 
# 1 median depth in normal cells < 15
# 2 median depth in normal cells > 60 (this is quite rare anyway)

######################################################################################################
# Filters with arbitrary thresholds 2: beta-binomial (rho - overdispersion parameter)

twins_normal = twins_dt[, c('mut_ID', samples_normal_mtr, samples_normal_dep), with=FALSE]
twins_normal[, sum_MTR_PD62341 := rowSums(.SD), .SDcols = samples_normal_PD62341_mtr]
twins_normal[, sum_DEP_PD62341 := rowSums(.SD), .SDcols = samples_normal_PD62341_dep]
twins_normal[, sum_MTR_PD63383 := rowSums(.SD), .SDcols = samples_normal_PD63383_mtr]
twins_normal[, sum_DEP_PD63383 := rowSums(.SD), .SDcols = samples_normal_PD63383_dep]
twins_normal[, sum_MTR_all := sum_MTR_PD62341 + sum_MTR_PD63383]
twins_normal[, sum_DEP_all := sum_DEP_PD62341 + sum_DEP_PD63383]

# I am running this filter only on the data from normal mutations (12 samples)
germline_mutations = read.table('Data/mutations_likely_germline_20241026.txt', header = F) %>% unlist()
twins_normal_filtered = twins_normal[!mut_ID %in% germline_mutations]
paste('Number of mutations which are likely not germline:', dim(twins_normal_filtered)[1]) # 11,987

# initialize possible values of rho to screen 
rhovec = 10^seq(-6,-0.05,by=0.01) # values to search through (grid search)
rhos = c()

# for each row in the dataframe, identify the vector of MTRs and vector of DEPs
for (i in 1:dim(twins_normal_filtered)[1]){
  
  mtrs = as.numeric(twins_normal_filtered[i, samples_normal_mtr, with=FALSE] %>% unlist())
  deps = as.numeric(twins_normal_filtered[i, samples_normal_dep, with=FALSE] %>% unlist())
  mu=sum(mtrs)/sum(deps)
  ll = sapply(rhovec, function(rhoj) sum(dbetabinom(x=mtrs, size=deps, rho=rhoj, prob=mu, log=T)))
  rhos[i] = rhovec[ll==max(ll)][1]
  
}

# merge rhos with the data table 
twins_normal_filtered = cbind(twins_normal_filtered, data.table(rhos))

thresholds = c(0.0005, 0.001, 0.005, 0.01, 0.025, 0.05, 0.075, 0.1)

for (t in thresholds){
  
  print(paste('Threshold:', t))
  muts_betabinomial = twins_normal_filtered[rhos < t, mut_ID] %>% unlist()
  print(paste('Number of mutations excluded at this rho threshold:', length(muts_betabinomial)))
  print(paste('Poorly mapped mutations excluded:', sum(poorly_mapped_muts %in% muts_betabinomial) / length(poorly_mapped_muts)))
  print(paste('Well mapped mutations excluded:', sum(well_mapped_muts %in% muts_betabinomial) / length(well_mapped_muts)))
  
}


muts_betabinomial = twins_normal_filtered[rhos < 0.001, mut_ID] %>% unlist()
to_vec(for (mut in well_mapped_muts) if (mut %in% muts_betabinomial) mut)
# some of these are mutations exclusive to PD63383 tumour
# I think we need to be mindful of these and make sure we are not removing those 
# maybe require this threshold unless sum_PD63383_tumour == 2 (as this is one case where this could be an issue)

######################################################################################################
# Filters with arbitrary thresholds 3: ratio of reads in low- vs high-quality pileup 

# read in files with high and low quality pileup 
twins_hq = data.table(read.csv('Data/pileup_merged_20241016.tsv'))
twins_lq = data.table(read.csv('Data/pileup_merged_20241025.tsv'))

twins_PDv38is = grep("PDv38is", names(twins_hq), value = TRUE)
twins_hq[, c(twins_PDv38is) := NULL]

twins_mtr_hq = twins_hq[,c('mut_ID', samples_mtr), with=FALSE]
twins_dep_hq = twins_hq[,c('mut_ID', samples_dep), with=FALSE]
twins_vaf_hq = twins_hq[,c('mut_ID', samples_vaf), with=FALSE]
twins_mtr_lq = twins_lq[,c('mut_ID', samples_mtr), with=FALSE]
twins_dep_lq = twins_lq[,c('mut_ID', samples_dep), with=FALSE]
twins_vaf_lq = twins_lq[,c('mut_ID', samples_vaf), with=FALSE]

twins_mtr_hq[,sum_mtr_mapped := rowSums(.SD), .SDcols = samples_mtr]
twins_mtr_lq[,sum_mtr_mapped := rowSums(.SD), .SDcols = samples_mtr]

twins_dep_hq[,sum_dep_mapped := rowSums(.SD), .SDcols = samples_dep]
twins_dep_lq[,sum_dep_mapped := rowSums(.SD), .SDcols = samples_dep]

muts_hq = merge(twins_mtr_hq[, c('mut_ID','sum_mtr_mapped'), with=FALSE], twins_dep_hq[, c('mut_ID', 'sum_dep_mapped'), with=FALSE])
muts_lq = merge(twins_mtr_lq[, c('mut_ID','sum_mtr_mapped'), with=FALSE], twins_dep_lq[, c('mut_ID', 'sum_dep_mapped'), with=FALSE])

setnames(muts_hq, c('sum_mtr_mapped', 'sum_dep_mapped'), c('sum_mtr_hq', 'sum_dep_hq'))
setnames(muts_lq, c('sum_mtr_mapped', 'sum_dep_mapped'), c('sum_mtr_lq', 'sum_dep_lq'))

muts_hl = merge(muts_hq, muts_lq)
muts_hl[, ratio_mtr := sum_mtr_hq / sum_mtr_lq]
muts_hl[, ratio_dep := sum_dep_hq / sum_dep_lq]
muts_hl[, ratio_vaf := (sum_mtr_hq / sum_dep_hq) / (sum_mtr_lq / sum_dep_lq)]

# test different ratios
ratios = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)

for (ratio in ratios){
  
  muts_ratio = muts_hl[ratio_mtr < ratio, mut_ID] %>% unlist()
  print(paste('Number of mutations excluded at HQ/LQ ratio:', ratio, ':', length(muts_ratio)))
  print(paste('Poorly mapped mutations excluded:', sum(poorly_mapped_muts %in% muts_ratio) / length(poorly_mapped_muts)))
  print(paste('Well mapped mutations excluded:', sum(well_mapped_muts %in% muts_ratio) / length(well_mapped_muts)))
  
}

muts_ratio07 = muts_hl[ratio_mtr < 0.7, mut_ID] %>% unlist()
to_vec(for (mut in well_mapped_muts) if (mut %in% muts_ratio07) mut)
# does Jbrowse show you the poorly mapped reads? because this looks really clean actually 
# either way, it does look like distribution wise this mutation is not very helpful 
# therefore, I think this threshold (maybe 0.6 to be sure) is okay
