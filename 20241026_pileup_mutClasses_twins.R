# Script to analyse the pileup (high quality, run 15/10/2024)
# 2024-10-25
# Barbara Walkowiak bw18

# INPUT: 
# 1 merged pileup dataframes with mutation calls (from CaVEMan) for tumour and normal samples
# 2 list of mutations that passed required filters 

# Cleaned up script to identify mutations of interest 
# Script to look at twin-specific mutations and early development

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
col_normal_PD62341 = "#058064"
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
# Mutations that define each twin 

######################################################################################################
# Mutations private to PD63383

twins_filtered_mtr[sum_normal_PD62341 == 0 & sum_normal_PD63383 == 6] # 0
twins_filtered_mtr[sum_normal_PD62341 == 1 & sum_normal_PD63383 == 6] # 3
# all present in PD62341v (spleen?), all at VAF ~0.1?
# chr11_34011887_C_T looks okay (vaf ~0.25 in PD63383)
# chr4_75704880_G_A looks okay (vaf ~0.25 in PD63383)
# chr6_165179306_G_A looks okay (vaf ~0.25 in PD63383)

twins_filtered_mtr[sum_normal_PD62341 == 2 & sum_normal_PD63383 == 6] # 5
# chr13_80917269_T_G # poor mapping 
# chr1_77763_G_A # poor mapping 
# chr3_105182212_T_C # poor mapping 
# chr3_165901319_C_A # poor mapping 
# chr6_55120674_A_G # not sure - checked again and I think this is rubbish 

twins_filtered_mtr[sum_normal_PD62341 == 0 & sum_normal_PD63383 == 5] # 1
# chr13_50815806_A_G # looks okay (3 reads in PD63383u but VAF 0.15)

twins_filtered_mtr[sum_normal_PD62341 == 1 & sum_normal_PD63383 == 5] # 6
# 1 in aa and the others in PD62341v, missing in PD63383u / PD63383ae / PD63383ak (dist)
# interestingly the 3 other mutations are absent from tumour samples: another evidence this is really PD62341 cells 
# these ones could be really cool to look at, talk to Henry about those!
# chr13_18341947_A_G # mis-mapping 
# chr17_18844573_T_C # mis-mapping 
# chr1_16849834_G_T # poor mapping 
# chr1_35202938_A_G # indels 
# chr4_74625500_G_T # oh that looks great! 
# chr7_73831920_C_T # looks okay

twins_filtered_mtr[sum_normal_PD62341 == 2 & sum_normal_PD63383 == 5] # 22
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

# Create a list that contains mutations private to PD63383 
muts_PD63383 = c("chr11_34011887_C_T", "chr4_75704880_G_A", "chr6_165179306_G_A",
                        "chr4_74625500_G_T", "chr7_73831920_C_T", "chr13_50815806_A_G")
paste('Number of mutations specific to PD63383:', length(muts_PD63383))

######################################################################################################
# Mutations private to PD62341 
twins_filtered_mtr[sum_normal_PD62341 == 6 & sum_normal_PD63383 == 0] # 2
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
# chr19_468950_G_A # looks messy - double check with Henry (not taking it)
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

muts_PD62341 = c("chr14_105458006_C_A", "chr17_33422229_C_A", "chr16_5479739_C_T", "chr2_95662131_G_A", "chr3_50106043_C_T") 
paste('Number of mutations specific to PD62341:', length(muts_PD62341))

muts_early = c(muts_PD62341, muts_PD63383)
paste('Number of informative early developmental mutations identified:', length(muts_early))

######################################################################################################
######################################################################################################
# Data visualization  

######################################################################################################
# Mutations specific to PD62341 
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

# VAF plots for different samples 
for (sample_name in mut_PD62341_melt[,sample] %>% unlist() %>% unique()){
  ggplot(mut_PD62341_melt[sample==sample_name], aes(x=mut_ID, y=value))+
    geom_line(group = 1)+
    geom_point(size=1.5)+
    geom_area(aes(fill = mut_ID))+
    scale_fill_manual(values = col_normal_PD62341)+
    theme_bw(base_size = 10)+
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
ggsave(glue('Results/20241024_p3_vaf_dist_PD62341_mut5_aggregated.pdf'), width=9, height=7)

# VAF distribution of each mutation in each sample 
ggplot(mut_PD62341_melt, aes(x=value, y=mut_ID, colour=sample_type))+
  geom_point(size=3.5)+
  scale_color_manual(values = c(col_normal_PD62341, col_normal_PD63383, col_tumour_PD62341, col_tumour_PD63383), name = 'Sample')+
  theme_bw(base_size = 18)+
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())+
  labs(x = 'VAF', y = 'Mutation')+
  ggtitle(glue('PD62341-specific mutations'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlim(c(0, 0.7))+
  geom_vline(xintercept = 0.5, color = "black", size=0.7)
ggsave(glue('Results/20241024_p3_vaf_dist_PD62341_mut5_samples.pdf'), width=9, height=7)

# only normal samples considered 
ggplot(mut_PD62341_melt[sample_type %in% c('normal_PD62341', 'normal_PD63383')], aes(x=value, y=mut_ID, colour=sample_type))+
  geom_point(size=3.5)+
  scale_color_manual(values = c(col_normal_PD62341, col_normal_PD63383), name = 'Sample')+
  theme_bw(base_size = 18)+
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())+
  labs(x = 'VAF', y = 'Mutation')+
  ggtitle(glue('PD62341-specific mutations'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlim(c(0, 0.7))+
  geom_vline(xintercept = 0.5, color = "black", size=0.7)
ggsave(glue('Results/20241024_p3_vaf_dist_PD62341_mut5_samples_normal.pdf'), width=9, height=7)

# are these mutations classified as germline in PD62341?
muts_germline_PD62341 = twins_dt_filters[f6_likelyGermline_PD62341==1, mut_ID] %>% unlist()
paste('Fraction of PD62341-specific mutations also identified as PD62341 likely germline:', sum(muts_PD62341 %in% muts_germline_PD62341)/length(muts_PD62341))
# all identified as likely germline - which may suggest that 

# maybe heatmap is a good idea in the end 
mut_PD62341_dt_normal = mut_PD62341_dt[, c('mut_ID', samples_normal_vaf), with=FALSE]
setcolorder(mut_PD62341_dt_normal, c('mut_ID', "PD62341aa_VAF",  "PD62341ad_VAF",  "PD62341n_VAF", "PD62341q_VAF", "PD62341h_VAF",  "PD62341v_VAF" ,
                              "PD63383w_VAF",  "PD63383t_VAF",  "PD63383u_VAF",  "PD63383ae_VAF", "PD63383ak_VAF", "PD63383bb_VAF"))
# matrix with VAF data (not binary)
mat_vaf = as.matrix(mut_PD62341_dt_normal[,2:13])
rownames(mat_vaf) =  mut_PD62341_dt_normal[,1] %>% unlist()  
colnames(mat_vaf) = c("PD62341aa",  "PD62341ad",  "PD62341n", "PD62341q", "PD62341h",  "PD62341v",
                      "PD63383w",  "PD63383t",  "PD63383u",  "PD63383ae", "PD63383ak", "PD63383bb")

# annotation of column names
col_annotation = data.frame(Twin = c(rep('PD62341', 6), rep('PD63383', 6)))
rownames(col_annotation) = colnames(mat_vaf)
annotation_colors = list(Twin = c(PD62341=col_normal_PD62341, PD63383=col_normal_PD63383))

pdf('Results/20241024_p4_heatmap_status_PD62341_mut5.pdf')
pheatmap(mat_vaf,
         cellwidth=14, cellheight=14,
         annotation_col = col_annotation,
         annotation_colors = annotation_colors,
         main="Status of PD62341 mutations", 
         legend = F,
         cluster_rows=F, cluster_cols = F, 
         show_rownames = T, show_colnames = T,
         fontsize=11, cexCol=6) 
dev.off()

######################################################################################################
# Mutations specific to PD63383 
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
                                      c("chr11_34011887_C_T", "chr4_75704880_G_A", "chr6_165179306_G_A",
                                        "chr4_74625500_G_T", "chr7_73831920_C_T", "chr13_50815806_A_G"))]

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
  ggsave(glue('Results/20241024_p4_vaf_dist_PD63383_mut6_{sample_name}.pdf'), width=9, height=7)
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
ggsave(glue('Results/20241024_p4_vaf_dist_PD63383_mut6_aggregated.pdf'), width=9, height=7)

# VAF distribution of each mutation in each sample 
ggplot(mut_PD63383_melt, aes(x=value, y=mut_ID, colour=sample_type))+
  geom_point(size=3.5)+
  scale_color_manual(values = c(col_normal_PD62341, col_normal_PD63383, col_tumour_PD62341, col_tumour_PD63383), name = 'Sample')+
  theme_bw(base_size = 18)+
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())+
  labs(x = 'VAF', y = 'Mutation')+
  ggtitle(glue('PD63383-specific mutations'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlim(c(0, 0.7))+
  geom_vline(xintercept = 0.25, color = "black", size=0.7)+
  geom_vline(xintercept = 0.125, color = "grey", size=0.7)
ggsave(glue('Results/20241024_p3_vaf_dist_PD63383_mut6_samples.pdf'), width=9, height=7)

# only normal samples considered 
ggplot(mut_PD63383_melt[sample_type %in% c('normal_PD62341', 'normal_PD63383')], aes(x=value, y=mut_ID, colour=sample_type))+
  geom_point(size=3.5)+
  scale_color_manual(values = c(col_normal_PD62341, col_normal_PD63383), name = 'Sample')+
  theme_bw(base_size = 18)+
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())+
  labs(x = 'VAF', y = 'Mutation')+
  ggtitle(glue('PD63383-specific mutations'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlim(c(0, 0.7))+
  geom_vline(xintercept = 0.25, color = "black", size=0.7)+
  geom_vline(xintercept = 0.125, color = "grey", size=0.7)
ggsave(glue('Results/20241024_p3_vaf_dist_PD63383_mut6_samples_normal.pdf'), width=9, height=7)

# are these mutations classified as germline in PD62341?
muts_germline_PD63383 = twins_dt_filters[f6_likelyGermline_PD63383==1, mut_ID] %>% unlist()
paste('Fraction of PD63383-specific mutations also identified as PD63383 likely germline:', sum(muts_PD63383 %in% muts_germline_PD63383)/length(muts_PD63383))
# none identified as likely germline - which may suggest that the true VAF is more likely 0.25 or 0.125

# maybe heatmap is a good idea in the end 
mut_PD63383_dt_normal = mut_PD63383_dt[, c('mut_ID', samples_normal_vaf), with=FALSE]
setcolorder(mut_PD63383_dt_normal, c('mut_ID', "PD62341aa_VAF",  "PD62341ad_VAF",  "PD62341n_VAF", "PD62341q_VAF", "PD62341h_VAF",  "PD62341v_VAF" ,
                                     "PD63383w_VAF",  "PD63383t_VAF",  "PD63383u_VAF",  "PD63383ae_VAF", "PD63383ak_VAF", "PD63383bb_VAF"))
# matrix with VAF data (not binary)
mat_vaf = as.matrix(mut_PD63383_dt_normal[,2:13])
rownames(mat_vaf) =  mut_PD63383_dt_normal[,1] %>% unlist()  
colnames(mat_vaf) = c("PD62341aa",  "PD62341ad",  "PD62341n", "PD62341q", "PD62341h",  "PD62341v",
                      "PD63383w",  "PD63383t",  "PD63383u",  "PD63383ae", "PD63383ak", "PD63383bb")

# annotation of column names
col_annotation = data.frame(Twin = c(rep('PD62341', 6), rep('PD63383', 6)))
rownames(col_annotation) = colnames(mat_vaf)
annotation_colors = list(Twin = c(PD62341=col_normal_PD62341, PD63383=col_normal_PD63383))

pdf('Results/20241024_p4_heatmap_status_PD63383_mut6.pdf')
pheatmap(mat_vaf,
         cellwidth=14, cellheight=14,
         annotation_col = col_annotation,
         annotation_colors = annotation_colors,
         main="Status of PD63383 mutations", 
         legend = F,
         cluster_rows=F, cluster_cols = F, 
         show_rownames = T, show_colnames = T,
         fontsize=11, cexCol=6) 
dev.off()

######################################################################################################
# Reconstructing the phylogeny 

# Mutations specific to PD62341 are present at VAF 0.5 
# That means that they are present in all cells of PD62341 (on one copy of the chr)
# They are absent from PD63383 



######################################################################################################
######################################################################################################
# Can I do a PCA on filtered mutations and see if samples separate?

library(gridExtra)
library(ggplot2)
library(dplyr)
library(tibble)
library(ggrepel)
library("FactoMineR") ## PCA
library("factoextra") ## PCA 

twins_filtered_vaf_mat = twins_filtered_vaf[, c(1:23), with=FALSE]
mut_mat = data.table(t(twins_filtered_vaf_mat[, c(2:23)]))
mut_mat[, sample_name := colnames(twins_filtered_vaf_mat[,2:23])]
mut_mat[, sample := tstrsplit(sample_name, '_', fixed=TRUE, keep = 1)]
mut_mat[, status := as.factor(fcase( 
  sample %in% samples_normal, 'normal', # differs from 0.5 so biased and maybe we don't want it
  sample %in% samples_tumour, 'tumour'
))]
mut_mat[, twin := as.factor(fcase( 
  sample %in% samples_PD62341, 'PD62341', # differs from 0.5 so biased and maybe we don't want it
  sample %in% samples_PD63383, 'PD63383'
))]
mut_mat[, sample_type := as.factor(paste(status, twin, sep = '_'))]

res.pca <- PCA(mut_mat[,1:418],  graph = FALSE) # PCA on JUST numerical
var <- get_pca_var(res.pca)
pdf('Results/20241028_p5_pca2.pdf')
PCA_figure <- fviz_pca_ind(res.pca,
                              geom = 'point',
                              legend.title = "Sample",
                              pointsize = 2, pointshape = 21, 
                              mean.point = FALSE,
                              fill.ind = as.factor(mut_mat[, sample_type]),
                              habillage = as.factor(mut_mat[, sample_type]),
                              title="PCA: 418 filtered mutations",
                              xlim = c(-30, 30),
                              ylim = c(-30, 30),
                              palette = c(col_normal_PD62341, col_normal_PD63383, col_tumour_PD62341, col_tumour_PD63383))
PCA_figure
dev.off()


