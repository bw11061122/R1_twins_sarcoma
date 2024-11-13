# Script to analyse the pileup (run 15/10/2024)
# 2024-10-11 - 2024-10-14
# Barbara Walkowiak bw18

# INPUT: merged pileup dataframes with mutation calls (from CaVEMan) for tumour and normal samples
# dataframe created in the script: 20241015_pileup_checks.R (pileup run 14/10/2024-15/10/2024)

# OUTPUT:
# 1 list of mutations which pass required filters and will be used for phylogeny reconstruction
# 2 dataframe with pass/fail status for each mutation for each filter 

# NOTE ON ADDING COLUMNS WITH FILTERS 
# 0 indicates that filter is passed; 1 indicates that filter is failed (sample does not meet the criteria and should be excluded) 

###################################################################################################################################
# LIBRARIES

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
library(grid)
library(Biostrings)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)

###################################################################################################################################
# INPUT FILES 
# Read the merged dataframe 
setwd('/Users/bw18/Desktop/1SB')
twins_dt = data.table(read.csv('Data/pileup_merged_20241016.tsv')) # import high quality pileup

# Drop columns with PD38is_wgs (used as reference)
twins_PDv38is = grep("PDv38is", names(twins_dt), value = TRUE)
twins_dt[, c(twins_PDv38is) := NULL]

# Create a dataframe that only includes mutations retained post filtering  
muts = read.table('Data/mutations_include_20241106_1002.txt') %>% unlist()
paste('Number of mutations that passed required filters:', length(muts)) # 1002
twins_filtered_dt = twins_dt[mut_ID %in% muts]

###################################################################################################################################
# PLOT SETTINGS

# Specify settings for plotting 
col_tumour = '#ad0505'
col_normal = '#07a7d0'
col_PD62341 = "#8909c1"
col_PD63383 = "#bca4f6"
col_tumour_PD62341 = "#980505"
col_tumour_PD63383 = "#eb6767"
col_normal_PD62341 = "#0785a5"
col_normal_PD63383 = "#70c6db"
col_bar = '#e87811'

######################################################################################################
# SAMPLES

# Examine which samples are being analysed 
samples = grep("PD6", names(twins_dt), value = TRUE)
samples = sapply(strsplit(samples,"_"), `[`, 1)
paste('Number of samples analysed:', length(samples)) # remove the PD38is_wgs sample: not from the samples of interest

# create lists of possible samples of interest
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

# Print the number of unique mutations identified in at least one sample (ie how many mutations do I have?)
paste('Number of unique mutations identified across samples:', dim(twins_dt)[1]) # 362,814
paste('Number of unique mutations identified across samples (check):', length(twins_dt[, mut_ID] %>% unlist() %>% unique())) # 362,814

######################################################################################################
# Accounting for twin-twin transfusion (spleen samples) - I want to have quantitative estimates of how much transfer there is 

muts_PD63383 = c(muts_normal_only_val, 'chr1_38827952_C_A')
muts_PD62341 = c("chr14_105458006_C_A", "chr17_33422229_C_A", "chr15_49480646_T_A",
                 "chr16_5479739_C_T", "chr2_95662131_G_A", "chr3_50106043_C_T",
                 "chr3_62055057_C_G", "chr3_62055077_G_C", "chr20_44114996_C_T") 
mut_early = c(muts_PD62341, muts_PD63383)

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
                                      c('chr16_5479739_C_T','chr15_49480646_T_A',
                                        'chr17_33422229_C_A','chr14_105458006_C_A',  
                                        'chr3_50106043_C_T', 'chr2_95662131_G_A',
                                        "chr3_62055057_C_G", "chr3_62055077_G_C", 
                                        "chr20_44114996_C_T"))]

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

# if this occurs, the spleen in PD62341 should be more similar to PD63383 than other PD62341 tissues
# at the same time, we expect spleen in PD63383 to be more similar to PD62341 than other PD63383 tissues 
# NB this is complicated because the tumour (which arose in PD62341) contaminated the skin sample (PD63383bb)

# check PD62341v and PD63383w in PD63383-specific mutations
mut_PD63383_melt[, sample_type := as.factor(paste(
  status, twin, sep = '_'))]

mut_PD63383_melt[, sample_type2 := as.factor(fcase(
  sample == 'PD62341v', 'normal_PD62341_spleen',
  sample == 'PD63383w', 'normal_PD63383_spleen',
  !sample %in% c('PD62341v', 'PD63383w'), paste(status, twin, sep = '_')
))]

mut_PD63383_melt[, sample_type3 := as.factor(fcase(
  sample == 'PD62341v', 'normal_PD62341_spleen',
  sample == 'PD63383w', 'normal_PD63383_spleen',
  sample == 'PD63383bb', 'normal_PD63383_skin',
  !sample %in% c('PD62341v', 'PD63383w', 'PD63383bb'), paste(status, twin, sep = '_')
))]

mut_PD62341_melt[, sample_type2 := as.factor(fcase(
  sample == 'PD62341v', 'normal_PD62341_spleen',
  sample == 'PD63383w', 'normal_PD63383_spleen',
  !sample %in% c('PD62341v', 'PD63383w'), paste(status, twin, sep = '_')
))]

mut_PD62341_melt[, sample_type3 := as.factor(fcase(
  sample == 'PD62341v', 'normal_PD62341_spleen',
  sample == 'PD63383w', 'normal_PD63383_spleen',
  sample == 'PD63383bb', 'normal_PD63383_skin',
  !sample %in% c('PD62341v', 'PD63383w', 'PD63383bb'), paste(status, twin, sep = '_')
))]

ggplot(mut_PD63383_melt[status=='normal'], aes(x=mut_ID, y=value, color=sample_type2))+
  geom_point(size=2.5, position = position_jitterdodge(0.4))+
  scale_color_manual(values = c(col_PD62341, '#C42F0F', col_PD63383, '#F99A49'))+
  theme_classic(base_size = 14)+
  labs(x = 'Mutation', y = 'VAF', col = 'Sample category')+
  ggtitle(glue('PD63383-specific mutations'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylim(c(0, 0.7))
ggsave(glue('Results/20241109_p3_vaf_dist_PD63383_muts_samples_normal_labelspleen.pdf'), width=7, height=4.5)

ggplot(mut_PD62341_melt[status=='normal'], aes(x=mut_ID, y=value, color=sample_type2))+
  geom_point(size=2.5, position = position_jitterdodge(0.4))+
  scale_color_manual(values = c(col_PD62341, '#C42F0F', col_PD63383, '#F99A49'))+
  theme_classic(base_size = 14)+
  labs(x = 'Mutation', y = 'VAF', col = 'Sample category')+
  ggtitle(glue('PD62341-specific mutations'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylim(c(0, 0.7))
ggsave(glue('Results/20241109_p3_vaf_dist_PD62341_muts_samples_normal_labelspleen.pdf'), width=7, height=4.5)

# Is there a way I can quantify this 
# Maybe first, compare mean VAF for 5 non-spleen samples and spleen 
means_PD62341muts_spleen_PD62341 = mut_PD62341_melt[sample == 'PD62341v', c('mut_ID', 'value'), with=FALSE]
means_PD62341muts_nonspleen_PD62341 = data.table(mut_PD62341_melt[sample %in% c('PD62341q','PD62341h', 'PD62341n', 'PD62341aa', 'PD62341ad'), mean(value), by = 'mut_ID'])
means_PD62341muts_spleen_PD63383 = mut_PD62341_melt[sample == 'PD63383w', c('mut_ID', 'value'), with=FALSE]
means_PD62341muts_nonspleen_PD63383 = data.table(mut_PD62341_melt[sample %in% c('PD63383t', 'PD63383u', 'PD63383ak','PD63383ae'), mean(value), by = 'mut_ID']) # ignore skin as this is contaminated
means_PD62341muts = merge(merge(means_PD62341muts_spleen_PD62341, means_PD62341muts_nonspleen_PD62341, by = 'mut_ID'),
                          merge(means_PD62341muts_spleen_PD63383, means_PD62341muts_nonspleen_PD63383, by = 'mut_ID'), by = 'mut_ID')
setnames(means_PD62341muts, c('value.x', 'V1.x', 'value.y', 'V1.y'), c('PD62341_spleen', 'PD62341_nonspleen', 'PD63383_spleen', 'PD63383_nonspleen'))

# I would do the calculations for mutations where VAF in PD63383 non-spleen is 0 for now
# VAF_PD62341_spleen = VAF_PD62341_spleen * fraction_PD62341_spleen + VAF_PD63383 * fraction_PD63383 
# VAF_PD63383_spleen = VAF_PD62341 * fraction_PD62341 + VAF_PD63383_spleen * fraction_PD63383_spleen
means_PD62341muts_clean = means_PD62341muts[PD63383_nonspleen==0]
means_PD62341muts_clean[, PD62341_spleen_fPD62341 := PD62341_spleen / PD62341_nonspleen]
means_PD62341muts_clean[, PD63383_spleen_fPD62341 := PD63383_spleen / PD62341_nonspleen]

ggplot(means_PD62341muts_clean, aes(x=PD62341_nonspleen, y=PD62341_spleen))+
  geom_point(size=2.5)+
  theme_classic(base_size = 14)+
  labs(x = 'mean VAF PD62341 (non-spleen)', y = 'VAF PD62341 (spleen)')+
  ggtitle(glue('PD62341-specific mutations'))+
  coord_equal(ratio = 1)+
  xlim(c(0, 0.5))+
  ylim(c(0, 0.5))
ggsave(glue('Results/20241109_spleen_vs_nonspleen_PD62341.pdf'), width=4, height=4)

# do the same for PD63383 specific mutations
means_PD63383muts_spleen_PD62341 = mut_PD63383_melt[sample == 'PD62341v', c('mut_ID', 'value'), with=FALSE]
means_PD63383muts_nonspleen_PD62341 = data.table(mut_PD63383_melt[sample %in% c('PD62341q','PD62341h', 'PD62341n', 'PD62341aa', 'PD62341ad'), mean(value), by = 'mut_ID'])
means_PD63383muts_spleen_PD63383 = mut_PD63383_melt[sample == 'PD63383w', c('mut_ID', 'value'), with=FALSE]
means_PD63383muts_nonspleen_PD63383 = data.table(mut_PD63383_melt[sample %in% c('PD63383t', 'PD63383u', 'PD63383ak','PD63383ae'), mean(value), by = 'mut_ID']) # ignore skin as this is contaminated
means_PD63383muts = merge(merge(means_PD63383muts_spleen_PD62341, means_PD63383muts_nonspleen_PD62341, by = 'mut_ID'),
                          merge(means_PD63383muts_spleen_PD63383, means_PD63383muts_nonspleen_PD63383, by = 'mut_ID'), by = 'mut_ID')
setnames(means_PD63383muts, c('value.x', 'V1.x', 'value.y', 'V1.y'), c('PD62341_spleen', 'PD62341_nonspleen', 'PD63383_spleen', 'PD63383_nonspleen'))
means_PD63383muts[, PD62341_spleen_fPD63383 := PD62341_spleen / PD63383_nonspleen] # contamination
means_PD63383muts[, PD63383_spleen_fPD63383 := PD63383_spleen / PD63383_nonspleen] # purity 

mut_PD62341_melt[, sample := factor(sample, levels = 
                                      c('PD62341ad', 'PD62341n', 'PD62341q', 'PD62341h', 'PD62341aa', 'PD62341v',
                                        'PD63383w', 'PD63383u', 'PD63383t', 'PD63383ae', 'PD63383ak', 'PD63383bb',
                                        'PD62341b', 'PD62341u', 'PD62341ae', 'PD62341ag', 'PD62341aj', 'PD62341ak', 'PD62341ap', 'PD62341am',
                                        'PD63383ap', 'PD63383aq'))]
mut_PD63383_melt[, sample := factor(sample, levels = 
                                      c('PD62341ad', 'PD62341n', 'PD62341q', 'PD62341h', 'PD62341aa', 'PD62341v',
                                        'PD63383w', 'PD63383u', 'PD63383t', 'PD63383ae', 'PD63383ak', 'PD63383bb',
                                        'PD62341b', 'PD62341u', 'PD62341ae', 'PD62341ag', 'PD62341aj', 'PD62341ak', 'PD62341ap', 'PD62341am',
                                        'PD63383ap', 'PD63383aq'))]


# plot VAF for each mutation across samples so we can see samples separately
for (mut in muts_PD62341){
  dt = mut_PD62341_melt[mut_ID == mut]
  ggplot(dt %>% arrange(sample_type), aes(x=sample, y=value, col = sample_type))+
    geom_point(size=2.5)+
    theme_classic(base_size = 14)+
    labs(x = 'sample', y = 'VAF')+
    scale_color_manual(values = c(col_PD62341, col_PD63383, col_tumour, col_tumour_PD62341))+
    ggtitle(glue('{mut}'))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(glue('Results/20241109_PD62341_spec_by_mut_{mut}.pdf'), width=6, height=3.5)
}

for (mut in muts_PD63383){
  dt = mut_PD63383_melt[mut_ID == mut]
  ggplot(dt %>% arrange(sample_type), aes(x=sample, y=value, col = sample_type))+
    geom_point(size=2.5)+
    theme_classic(base_size = 14)+
    labs(x = 'sample', y = 'VAF')+
    scale_color_manual(values = c(col_PD62341, col_PD63383, col_tumour, col_tumour_PD62341))+
    ggtitle(glue('{mut}'))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(glue('Results/20241109_PD63383_spec_by_mut_{mut}.pdf'), width=6, height=3.5)
}
