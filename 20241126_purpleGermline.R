###################################################################################################################################
# SCRIPT 2

# Script to remove germline mutations in copy number-altered regions
# 2024-11-26
# Barbara Walkowiak bw18

# INPUT: 
# 1 PURPE dt (.cnv.segment) for each sample 

# OUTPUT:
# 1 list of mutations to exclude based on presence in copy number germline variants 

# Script to identify regions of copy number alterations (deletions or duplications) in the germline
# Identify mutations in those regions where the VAF is 0.3 and so those can be explained as a germline copy number alteration 

# UPDATE 04/12/2024: 
# After discussion with Henry, we decided to scrap PURPLE for now and use something else to get rid of the germline mutations
# on copy number-altered regions, as this is too painful and looks like it's making our life harder not easier 

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
library(beeswarm)
library(viridis)
library(grid)
library(pheatmap)
library(ggrepel)

###################################################################################################################################
# INPUT FILES 

# Read the merged dataframe 
setwd('/Users/bw18/Desktop/1SB')
twins_dt = data.table(read.csv('Data/twins_dt_20241114_1134.csv')) # import high quality pileup with added filters

# Drop columns with PD38is_wgs (used as reference)
twins_PDv38is = grep("PDv38is", names(twins_dt), value = TRUE)
twins_dt[, c(twins_PDv38is) := NULL]

# Create a dataframe that only includes mutations retained post filtering  
muts_included = read.table('Data/mutations_include_20241114_1134.txt') %>% unlist()
paste('Number of mutations that passed required filters:', length(muts_included)) # 1134

muts_germline = read.table('Data/mutations_putativeGermline_20241114.txt') %>% unlist()
paste('Number of mutations identified as putative germline:', length(muts_germline)) # 333,031

###################################################################################################################################
# PLOT SETTINGS

# Specify colors for plotting 
col_tumour = '#ad0505'
col_normal = '#07a7d0'
col_PD62341 = "#0ac368"
col_PD63383 = "#a249e8"
col_tumour_PD62341 = "#099272"
col_tumour_PD63383 = "#6F09D4"
col_normal_PD62341 = "#71D99B"
col_normal_PD63383 = "#C99DF6"

col_background = '#c8c8c8'
col_germline = '#d3b28b'
col_other = '#25b6db'
col_normal_all = '#5265c9'
col_normal_all_removed = '#b9dff5'
col_clusters = '#bf068a'
col_coverage = '#bf7a06'
col_clusters_coverage = '#810cb9'
col_excluded = '#a80505'

######################################################################################################
# SAMPLES

# Create lists of possible samples of interest
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
# AGGREGATED MTR, DEP and VAF (by sample category) 

# Aggregated VAF in tumour samples 
twins_dt[, dep_all_normal_PD62341 := rowSums(.SD), .SDcols = samples_normal_PD62341_dep]
twins_dt[, dep_all_normal_PD63383 := rowSums(.SD), .SDcols = samples_normal_PD63383_dep]
twins_dt[, dep_all_tumour_PD62341 := rowSums(.SD), .SDcols = samples_tumour_PD62341_dep]
twins_dt[, dep_all_tumour_PD63383 := rowSums(.SD), .SDcols = samples_tumour_PD63383_dep]
twins_dt[, dep_all_normal := rowSums(.SD), .SDcols = samples_normal_dep]
twins_dt[, dep_all_tumour := rowSums(.SD), .SDcols = samples_tumour_dep]
twins_dt[, dep_all := rowSums(.SD), .SDcols = samples_dep]

twins_dt[, mtr_all_normal_PD62341 := rowSums(.SD), .SDcols = samples_normal_PD62341_mtr]
twins_dt[, mtr_all_normal_PD63383 := rowSums(.SD), .SDcols = samples_normal_PD63383_mtr]
twins_dt[, mtr_all_tumour_PD62341 := rowSums(.SD), .SDcols = samples_tumour_PD62341_mtr]
twins_dt[, mtr_all_tumour_PD63383 := rowSums(.SD), .SDcols = samples_tumour_PD63383_mtr]
twins_dt[, mtr_all_normal := rowSums(.SD), .SDcols = samples_normal_mtr]
twins_dt[, mtr_all_tumour := rowSums(.SD), .SDcols = samples_tumour_mtr]
twins_dt[, mtr_all := rowSums(.SD), .SDcols = samples_mtr]

twins_dt[, vaf_all_normal_PD62341 := mtr_all_normal_PD62341 / dep_all_normal_PD62341]
twins_dt[, vaf_all_normal_PD63383 := mtr_all_normal_PD63383 / dep_all_normal_PD63383]
twins_dt[, vaf_all_tumour_PD62341 := mtr_all_tumour_PD62341 / dep_all_tumour_PD62341]
twins_dt[, vaf_all_tumour_PD63383 := mtr_all_tumour_PD63383 / dep_all_tumour_PD63383]
twins_dt[, vaf_all_normal := mtr_all_normal / dep_all_normal]
twins_dt[, vaf_all_tumour := mtr_all_tumour / dep_all_tumour]
twins_dt[, vaf_all := mtr_all / dep_all]

######################################################################################################
# Load purple data 

# Load purple.segment data 
purple_segment = data.table()
for (sample in samples_normal){
  purple_dt = data.table(read.table(paste0('PURPLE/', sample, '.purple.segment.tsv'), sep='\t', fill = TRUE, header=TRUE))
  purple_dt[, sample := sample]
  purple_segment = rbind(purple_segment, purple_dt)
}

# Number of segments identified across each sample
table(purple_segment[, sample]) # not the same number, but reasonable 
hist(purple_segment[, fittedTumorCopyNumber]) # only 0 values here?

# Load purple.cnv.somatic data
# The copy number file TUMOR.purple.cnv.somatic.tsv contains the copy number profile of all (contiguous) segments of the tumor sample
# copyNumber = fitted absolute copy number of segment adjusted for purity and ploidy
purple = data.table()
for (sample in samples_normal){
  purple_dt = data.table(read.table(paste0('PURPLE/', sample, '.purple.cnv.somatic.tsv'), sep='\t', fill = TRUE, header=TRUE))
  purple_dt[, sample := sample]
  purple = rbind(purple, purple_dt)
}

paste('Total number of copy number regions identified across all normal samples:', dim(purple)[1])

# Basic checks 
# What are the copy number values reported in PURPLE?
paste('Minimum copy number value identified:', min(purple[, copyNumber])) # -7.7
paste('Maximum copy number value identified:', max(purple[, copyNumber])) # 26.9
# This seems wrong to me. It isn't possible that the copy number is below 0. 

# Show the distribution of copy number calls across samples 
pdf('Results/20241127_hist_purple_copynumber.pdf')
hist(purple[, copyNumber] %>% unlist(), xlab = 'Copy Number', breaks=50)
dev.off()

pdf('Results/20241127_hist_purple_copynumber_from0to3.pdf')
hist(purple[copyNumber > 0 & copyNumber < 3, copyNumber] %>% unlist(), xlab = 'Copy Number', breaks=50)
abline(v = 2.1, col = 'purple', lwd = 2.5)
dev.off()

# checked with Nathan: use 2.5 as a threshold for CN gain and 1.5 for CN loss 
# Number of copy number region GAINS in each sample
table(purple[copyNumber > 2.5, sample])

# Which regions are called as copy number germline in each sample?

# Samples with a single copy number germline change
purple[copyNumber > 2.5 & sample=='PD62341n'] # chr4 124138001 124341000
purple[copyNumber > 2.5 & sample=='PD62341q'] # chr4 124145001 124341000
purple[copyNumber > 2.5 & sample=='PD62341v'] # chr4 124138001 124342000
purple[copyNumber > 2.5 & sample=='PD63383ak'] # chr4 124138001 124342000
purple[copyNumber > 2.5 & sample=='PD63383u'] # chr4 124138001 124342000

purple[copyNumber > 2.5 & sample=='PD63383w'] # 2 copy number changes 
# chr4 124138001 124341000
# chr14 106075001 106349000
# this is not called in PURPLE in other samples but it is very clear that this region has a higher coverage
# note that PURPLE also calls a loss in this sample next to this gain 

# other clean samples
purple[copyNumber > 2.5 & sample=='PD63383t'] # 6
# chr12  31849001  31911000 # clear coverage change across all normal samples 
# chr4 124138001 124341000
# chr5  69561001  71157000 # clear coverage change across all normal samples 
# chr7  65774001  65824000 # clear coverage change across all normal samples 
# chr7 100959001 100976000 # clear coverage change across all normal samples
# chr9  33794001  33805000 # doesn't look like a difference in coverage to me

purple[copyNumber > 2.5 & sample=='PD62341ad']
# chr16  32004001  34210000 # coverage change across all samples 
# chr16  34210001  34917000 # coverage change across all samples
# chr16  34917001  34957000 # coverage change across all samples
# chr4 124138001 124341000

purple[copyNumber > 2.5 & sample=='PD62341aa'] # 72 regions identified 

# samples we think are contaminated with the tumour 
purple[copyNumber > 2.5 & sample=='PD62341aa'] # 87 regions 
purple[copyNumber > 2.5 & sample=='PD62341h'] # 8 regions
# copy nr on chromosome 10, 38375501 to 46332000 
# copy nr on chromosome 16, 32012001 to 35016000
# copy nr on chromosome 4
purple[copyNumber > 2.5 & sample=='PD63383bb'] # 5 regions  
# copy nr on chromosome 3 10311001  11176000
# copy nr on chromosome 3  50111001  50660000
# copy nr on chromosome 6 33658001  34148000
# copy nr on chromosome X 58075001  58295000

# copy number loss
table(purple[copyNumber < 1.5, sample])
# PD62341aa - 48
# PD62341h = 2
# PD63383ae - 212
# PD63383bb - 13
# PD63383t - 102
# PD63383w - 1
# PD63383ak, u - 0, PD62341ad, n, v, q - 0

purple[copyNumber > 0 & copyNumber < 1.5 & sample=='PD62341aa'] # 26 regions  
# 13 regions on chromosome 10, all ~lost - likely tumour infiltration
# 5 regions on chromosome 18 - likely tumour infiltration
# 3 regions lost on chromosome X, 2 chr5, 1 chr10, 1 chr16, 1 chr17 - tumour?
purple[copyNumber > 0 & copyNumber < 1.5 & sample=='PD62341h'] # 2 regions  
# maps to chr1 and chr18 - likely tumour infiltration
purple[copyNumber > 0 & copyNumber < 1.5 & sample=='PD63383bb'] # 10 regions  
# 5 losses on chr1, 1 on chr18, 1 on chr17 and chr20, 2 on X
# I could believe this is tumour infiltration 

# Samples which I am quite sure are clean and not tumour infiltrated 
purple[copyNumber > 0 & copyNumber < 1.5 & sample=='PD63383w'] # 1 region  
# chr14 105865501-106075000
purple[copyNumber > 0 & copyNumber < 1.5 & sample=='PD63383ae'] # 180 regions  
# cannot explain what happened with this sample 
# reported losses across most chromosomes 
purple[copyNumber > 0 & copyNumber < 1.5 & sample=='PD63383t'] # 93 regions  
# reported losses across most chromosomes

# Identify regions considered as copy number germline in all samples 
# Use threshold of copy number 2.5 to remove mutations present on regions of 3 copies 
germline_cn = c()
for (chr in twins_coord[, Chrom] %>% unlist() %>% unique()){
  tc = twins_coord[Chrom==chr] # extract the data for each chromosome 
  purple_dt = purple[chromosome==chr & copyNumber >= 2.5]
  starts = purple_dt[, start] %>% unlist() %>% unique()
  ends = purple_dt[, end] %>% unlist() %>% unique()
  res = tc[purple_dt, on = .(Pos >= start, Pos <= end), nomatch = 0]
  germline_cn = c(germline_cn, res[, mut_ID] %>% unlist())
}

# Add column to indicate mutation being present on a germline copy number region 
germline_cn_list = germline_cn %>% unique() %>% unlist()
twins_dt[, germlineCN := fcase(
  mut_ID %in% germline_cn_list, 'copy number change (PURPLE)',
  mut_ID %in% muts_cov, 'excluded based on coverage',
  !mut_ID %in% c(germline_cn_list, muts_cov), 'no copy number change')]

twins_deps = twins_dt[Chrom!='chrY', c('Chrom', 'Pos', 'germlineCN', samples_dep), with=FALSE]
twins_deps_melt = melt(twins_deps, id.vars = c('Chrom', 'Pos', 'germlineCN'))
twins_deps_melt[, sample := tstrsplit(variable, '_', fixed=TRUE, keep=1)]
twins_deps_melt[, Chrom := factor(Chrom, levels = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
                                                    "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
                                                    "chr20", "chr21", "chr22", "chrX"))]
twins_deps_melt[, germlineCN := factor(germlineCN, levels = c('no copy number change', 'excluded based on coverage', 'copy number change (PURPLE)'))]

# plot different classes of mutations ~ coverage for each sample
for (s in samples_normal){
  ggplot(twins_deps_melt[sample==s] %>% arrange(germlineCN), aes(x = Pos, y = value, color = germlineCN, alpha= germlineCN, order = germlineCN))+
    geom_point(size=1.5)+
    scale_color_manual(values = c('grey', 'goldenrod', 'purple'))+
    scale_alpha_manual(values = c(0.1, 0.6, 0.8), guide = 'none')+
    theme_classic(base_size = 12)+
    facet_wrap(~Chrom,  scales = "free_x")+
    labs(x = 'Genomic position', y = 'Coverage ({s})', color = 'germline CN?')+
    guides(color = guide_legend(override.aes = list(alpha = 1) ) )+
    ggtitle(glue('{s}'))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5)) # smaller x axis ticks 
  ggsave(glue('Results/20241114_p2_dep_across_genome_{s}_allmuts.pdf'), height=10, width=15)
}

# plot again, but exclude mutations filtered out based on coverage
for (s in samples_normal){
  ggplot(twins_deps_melt[sample==s & germlineCN != 'excluded based on coverage'] %>% arrange(germlineCN), aes(x = Pos, y = value, color = germlineCN, alpha= germlineCN, order = germlineCN))+
    geom_point(size=1.5)+
    guides(color = guide_legend(override.aes = list(alpha = 1) ) )+
    scale_color_manual(values = c('grey', 'purple'))+
    scale_alpha_manual(values = c(0.1, 0.8), guide = 'none')+
    theme_classic(base_size = 12)+
    facet_wrap(~Chrom)+
    labs(x = 'Genomic position', y = 'Coverage ({s})', color = 'germline CN?')+
    ggtitle(glue('{s}'))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5)) # remove x axis ticks 
  ggsave(glue('Results/20241114_p2_dep_across_genome_{s}_allmuts_exclCov.pdf'), height=10, width=15)
}

# would it make sense to merge those regions or perhaps first plot them across the genome?
chromosome = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
          "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", 
          "chr20", "chr21", "chr22", "chrX")
# Source: https://www.ncbi.nlm.nih.gov/grc/human/data (Hg38), accessed 14/11/2024
chrom_length = as.numeric(c(248956422, 242193529, 198295559, 190214555,
                            181538259, 170805979, 159345973, 145138636,
                            138394717, 133797422, 135086622, 133275309,
                            114364328, 107043718, 101991189, 90338345,
                            83257441, 80373285, 58617616, 64444167,
                            46709983, 50818468, 156040895))
chrom_dt = data.table(cbind(chromosome, chrom_length))

purple = merge(purple, chrom_dt, by = 'chromosome')
purple[, chromosome := factor(chromosome, levels = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
                                             "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", 
                                             "chr20", "chr21", "chr22", "chrX"))]


######################################################################################################
# Filtering based on PURPLE outputs 
# Identify mutations close to indels present on each chromosome
# Use threshold of copy number 2.9 to remove mutations present on regions of 3 copies 
germline_cn = c()
for (chr in twins_coord[, Chrom] %>% unlist() %>% unique()){
  tc = twins_coord[Chrom==chr] # extract the data for each chromosome 
  purple_dt = purple[chromosome==chr & copyNumber >= 2.5]
  starts = purple_dt[, start] %>% unlist() %>% unique()
  ends = purple_dt[, end] %>% unlist() %>% unique()
  res = tc[purple_dt, on = .(Pos >= start, Pos <= end), nomatch = 0]
  germline_cn = c(germline_cn, res[, mut_ID] %>% unlist())
}




