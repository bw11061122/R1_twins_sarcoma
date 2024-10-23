# Script to analyse the pileup (run 15/10/2024)
# 2024-10-11 - 2024-10-14
# Barbara Walkowiak bw18

# INPUT: merged pileup dataframes with mutation calls (from CaVEMan) for tumour and normal samples
# with filters 
# NOTE that for the moment, I am investigating if there is sth wrong in how the pileup was run
# I may alter the input df once everything is resolved

# NOTE ON ADDING COLUMNS WITH FILTERS 
# I want to add a column which says 0 if sample is fine and 1 if sample is wrong 

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
# Read the merged dataframe 
setwd('/Users/bw18/Desktop/1SB')
twins_dt = data.table(read.csv('Data/pileup_merged_20241016.tsv'))

# Drop columns with PD38is_wgs (used as reference)
twins_PDv38is = grep("PDv38is", names(twins_dt), value = TRUE)
twins_dt[, c(twins_PDv38is) := NULL]

# Filter to only include mutations retained for filtering 
muts = read.table('Data/mutations_include_20241023.txt') %>% unlist()
paste('Number of mutations that passed required filters:', length(muts)) # 1,966
twins_filtered_dt = twins_dt[mut_ID %in% muts]

###################################################################################################################################
# Specify settings for plotting 
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
# Add column to indicate chromosomes lost in the tumour
twins_filtered_dt[, loss := as.factor(fcase( 
  Chrom %in% c('chr1', 'chr18'), 'loss in tumour', # chr1 and chr18 segments lost in tumour samples
  !Chrom %in% c('chr1', 'chr18'), 'normal ploidy'
))]

######################################################################################################
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

# Note: it's find that you see some columns where there is the max VAF in normal samples is 0 
# can i see these mutations
twins_filtered_vaf[, sum_vaf_normal0 := rowSums(.SD==0), .SDcols = samples_normal_vaf]
twins_filtered_vaf[sum_vaf_normal0==12, mut_ID]
# all look good 
# "chr10_72139813_C_T" 
# "chr15_68707110_G_A" 
# "chr3_189240074_C_A" 
# "chr3_95470911_T_A"  
# "chr5_28014472_C_T"  
# "chr5_58182780_T_A"  
# "chr8_91252357_A_C"  
# "chrX_64561555_C_T" 

######################################################################################################
# Examine interesting classes of mutations

twins_filtered_mtr[, sum_tumour := rowSums(.SD>=4), .SDcols = samples_tumour_mtr]
twins_filtered_mtr[, sum_normal := rowSums(.SD>=4), .SDcols = samples_normal_mtr]
twins_filtered_mtr[, sum_PD62341 := rowSums(.SD>=4), .SDcols = samples_PD62341_mtr]
twins_filtered_mtr[, sum_PD63383 := rowSums(.SD>=4), .SDcols = samples_PD63383_mtr]
twins_filtered_mtr[, sum_tumour_PD62341 := rowSums(.SD>=4), .SDcols = samples_tumour_PD62341_mtr]
twins_filtered_mtr[, sum_tumour_PD63383 := rowSums(.SD>=4), .SDcols = samples_tumour_PD63383_mtr]
twins_filtered_mtr[, sum_normal_PD62341 := rowSums(.SD>=4), .SDcols = samples_normal_PD62341_mtr]
twins_filtered_mtr[, sum_normal_PD63383 := rowSums(.SD>=4), .SDcols = samples_normal_PD63383_mtr]

# select mutations which are only / mostly present in the tumour

twins_filtered_mtr[sum_tumour == 10 & sum_normal == 0] # 1 
# "chr2_199565129_G_A" # looks okay

twins_filtered_mtr[sum_tumour == 9 & sum_normal == 0] # 2 
# chr12_106354695_C_T # looks okay, missing in PD62341ak
# "chrX_115495236_T_C" # looks okay to me (double check but I would believe it), missing in PD62341b

twins_filtered_mtr[sum_tumour == 8 & sum_normal == 0] # empty  

twins_filtered_mtr[sum_tumour == 7 & sum_normal == 0] # 2  
# "chr14_45928737_G_T" # poor mapping
# "chr7_152562130_G_A" # looks okay, missing in PD62341ak, PD62341ag, PD62341u

twins_filtered_mtr[sum_tumour == 10 & sum_normal == 1] # 4 
# "chr16_32621521_C_A" # mapping quality not excellent
# "chr17_78256883_C_T"
# "chr22_17774668_G_A"  
# "chr22_30553442_G_T" 

twins_filtered_mtr[sum_tumour >= 6 & sum_normal <= 1] # 26 
twins_filtered_mtr[sum_tumour >= 6 & sum_normal <= 1, mut_ID] 

# "chr10_51734221_A_T"  
# "chr12_106354695_C_T" 
# "chr12_112652964_C_T" 
# "chr12_61020950_C_A" 
# "chr12_84029949_C_A"  
# "chr14_104122986_G_T" 
# "chr14_45928734_G_A" # mapping quality issues 
# "chr14_45928737_G_T" # mapping quality issues 
# "chr15_94939486_G_A"  
# "chr16_32621521_C_A" # mapping quality not excellent
# "chr17_78256883_C_T"  
# "chr19_3961910_G_A"  
# "chr22_17774668_G_A"  
# "chr22_30553442_G_T"  
# "chr2_199565129_G_A"  
# "chr2_57814739_A_C" # fair share of insertions 
# "chr6_58221412_T_C" # mapping issues   
# "chr6_77631804_C_A"   
# "chr7_115240063_T_C"  
# "chr7_120677593_C_T" 
# "chr7_152562130_G_A"  
# "chr8_46678612_T_C"   
# "chr8_87226547_A_G" # mapping issues   
# "chrX_115495236_T_C" 
# "chrX_46375010_C_T"   
# "chrX_56448796_C_A" 

# NOTE there is a mutation X_56448822_T_C which looks real and must have arisen earlier
# present in all samples and not just the tumour + C_A is in some reads with T_C but not all

######################################################################################################
# select mutations which are only present in the tumour and either twin

twins_filtered_mtr[sum_tumour >= 6 & sum_normal_PD62341 >= 4 & sum_normal_PD63383 == 0] # 5 (3 look real)
# "chr14_105458006_C_A" 
# "chr3_58819009_G_A" # mapping issues   
# "chr3_62055057_C_G"   
# "chr3_62055077_G_C" 
# "chr7_40558111_A_C" # insertions 

twins_filtered_mtr[sum_tumour >= 6 & sum_normal_PD62341 >= 1 & sum_normal_PD63383 == 1] # 108
# may allow one to be present in PD63383 bc of contamination 

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
# TRACKING TUMOUR EVOLUTION
# look at mutations which are present in a subset of tumour samples
twins_filtered_mtr[sum_normal == 0 & sum_tumour >= 2, mut_ID] #  mutations are only informative if shared 
# 88 such mutations 

# mutations exclusive to PD63383 tumour samples (NB allowing some normal bc contamination)
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
# "chrX_48453603_G_A" # looks believable but close to a region very prone to insertions
# "chrX_72671786_A_T" # poor mapping + strand bias   
# "chrX_77681162_G_A" # looks okay

PD63383_tumour_private1 = twins_filtered_mtr[sum_normal == 0 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 0, mut_ID] %>% unlist()
twins_info = twins_dt[, c('mut_ID', 'Gene', 'Transcript', 'RNA', 'CDS', 'Protein', 'Type', 'Effect'), with=FALSE]
twins_info[mut_ID %in% PD63383_tumour_private1] # nothing in protein coding genes

twins_filtered_mtr[sum_normal == 1 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 0, mut_ID] # 16
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

# "chr9_98275090_C_T" # looks okay # bb
# "chrX_117418608_G_A" # looks okay # bb
# "chrX_119915655_G_A" # looks okay # bb
# "chrX_36779694_C_T" # looks okay # bb
# "chrX_70352624_G_A" # looks okay # bb

# reads in normal samples are exclusively present in bb
# PD63383bb is skin, so I find it highly plausible that this is contamination with tumour cells
# especially that the kid had lesions in the skin - can I double check with Henry what he thinks?

PD63383_tumour_private2 = twins_filtered_mtr[sum_normal == 1 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 0, mut_ID] %>% unlist()
twins_info[mut_ID %in% PD63383_tumour_private2] # nothing in protein coding genes
# one mutation in POTEF (downstream) - doesn't sound like it would affect biology 

twins_filtered_mtr[sum_normal == 2 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 0, mut_ID] # 1
# "chr10_73678114_G_A" poor mapping quality 

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

twins_dt_germlinePD63383 = twins_dt_filters[f7_likelyGermline_PD62341==0 & f7_likelyGermline_PD63383==1] # 765
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

######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
# Another hypothesis that could explain mutation chr14_105458006_C_A 
# is if twins split later BUT cells were anatomically segregated such that all cells 
# in one twin but not another inherited a mutation that was very early on 
# maybe it is good to check whether there are mutations shared by some tissues in PD62341 and PD63383 
# let's also require max mtr >= 10
# tbf I am getting convinced that the VAST majority of those are artefacts 

twins_normal_mtr[,mean_mtr := apply(.SD, 1, function(x) mean(x[x>0])), .SDcols = samples_normal_mtr]
twins_normal_mtr[,median_mtr := apply(.SD, 1, function(x) median(x[x>0])), .SDcols = samples_normal_mtr]
twins_normal_mtr[,min_mtr := apply(.SD, 1, function(x) min(x[x>0])), .SDcols = samples_normal_mtr]
twins_normal_mtr[,max_mtr := apply(.SD, 1, max), .SDcols = samples_normal_mtr]
twins_normal_mtr[,sum_mtr0 := rowSums(.SD==0), .SDcols = samples_normal_mtr]

twins_normal_mtr[sum_normal_PD62341 == 2 & sum_normal_PD63383 == 2 & max_mtr >= 10] # 11
# "chr15_28398565_A_G" # poor mapping 
# "chr15_32482084_G_A" # poor mapping 
# "chr16_14864167_G_C" # poor mapping 
# "chr17_36230372_A_G" # poor mapping 
# "chr3_75122630_G_T" # poor mapping  
# "chr4_26808742_A_G" # insertions 
# "chr4_55624665_T_A" # poor mapping 
# "chr5_12769974_A_G" # poor mapping 
# "chr8_143864098_G_A" # poor mapping 
# "chr8_7423511_C_T" # poor mapping   
# "chr9_110300113_C_G" # very low coverage in several samples

twins_normal_mtr[sum_normal_PD62341 == 3 & sum_normal_PD63383 == 3 & max_mtr >= 10 & sum_mtr0 >= 1] # 6
# mean MTR of those is all ~4-5

twins_normal_mtr[sum_normal_PD62341 == 4 & sum_normal_PD63383 == 4 & max_mtr >= 10 & sum_mtr0 >= 1] # 10
# mean MTR of those is all ~4-6

twins_normal_mtr[sum_normal_PD62341 == 5 & sum_normal_PD63383 == 5 & max_mtr >= 10 & sum_mtr0 >= 1] # 16
# mean MTR of those is all ~5-6, in the case it is > 10, the mapping quality is poor

######################################################################################################
######################################################################################################
### RANDOM
twins_vaf[mut_ID =='chr14_105458006_C_A']
# VAF in PD62341v == 0.1892
# would be good to know the aggregate VAF in normal tissues
# my exact binomial classified it as likely germline in PD62341 so maybe this is also a filter to look at 


######################################################################################################
######################################################################################################

# OTHER MUTATIONS (OLD FILTERING)
# Checking some mutations on JBrowse:
# "chr10_100754461_C_T" - looks great 
# "chr13_100405282_G_T"- looks great 
# "chr6_160324993_G_A"
# "chr2_95914901_C_A" = a bit questionable but I would take it 
# "chr1_1433044_C_T" - I'd buy it 
# "chr12_39248376_A_T" - looks great 

# not looking very good:
# "chr10_45856845_G_A"
# "chr11_4313654_C_T"
# "chr11_4344722_C_A" - not terrible but quite poor mapping quality 
# "chr21_45822642_T_G" - phasing looks good but v questionable mapping
# "chr1_16891590_G_A" - sometimes questionable mapping - but I would believe this is real
# "chr1_111340150_T_G" - clear issues with mapping 
# "chr7_154076175_C_T" - issues with mapping
# "chr4_15679_A_T" - huge issues with mapping 
# "chr22_20717458_T_C" - definitely has some mapping issues 
# "chr1_109412643_A_C" - insertions in most reads 
# "chr18_41634003_C_G" - mapping issues 
# "chr14_82508185_C_A" - next to an insertion 
# "chr9_66076217_A_G" - looks terrible 
# "chr7_154076175_C_T" - poor mapping quality 
# "chr19_11361902_T_C" - very close to insertions so quite questionable 
# "chr10_30842788_C_A" - close to insertions 
# "chr14_19129291_C_T" - poor mapping quality 
# "chr1_248453436_T_C" - poor mapping quality 

# I think sample q is either contaminated or common origin 
# chr16_32621521_C_A mapping quality could be better
# chr22_16694683_C_A mapping quality is maybe not excellent but at the same time maybe it's okay - can check 
# chr2_57814739_A_C this is a bit too close to quite many insertions for comfort 


