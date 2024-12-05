###################################################################################################################################
# SCRIPT 9

# Script to analyse the scRNA-seq genotyping allele integrator output 
# 2024-12-04
# Barbara Walkowiak bw18

# INPUT: 
# 1 .tsv dataframes generated from alleleIntegrator
# allele integrator ran 04/12/2024; jobID 664854
# downloaded files from the farm via rsync -avhu bw18@farm22:"/lustre/scratch126/casm/team274sb/bw18/twins_sarcoma/scRNA_genotyping/output/" /Users/bw18/Desktop/1SB/scRNAseq/Data/

# OUTPUT:

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
library(RColorBrewer)
library(Biostrings)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)

###################################################################################################################################
# INPUT FILES 

setwd('/Users/bw18/Desktop/1SB')

# load files from all samples (I )

samples_scRNAseq = c('NB13652544', 'NB13652545', 'NB13652546',
                     'NB13760628', 'NB13760629', 'NB13760630',
                     'NB14599068', 'NB14599069', 'NB14599070')

ai_counts_dt = data.table()
for (sample in samples_scRNAseq){
  ai_counts = fread(paste0('scRNAseq/Data/', sample, '_somaticVar_scRNA_alleleCounts.tsv'), sep='\t', fill = TRUE, header=TRUE)
  ai_counts[, sample_ID := sample]
  ai_counts_dt = rbind(ai_counts_dt, ai_counts)
}

# add source of the sample
ai_counts_dt[, sample := factor(fcase(
  sample_ID %in% c('NB13652544', 'NB13652545', 'NB13652546'), 'PD62341',
  sample_ID %in% c('NB13760628', 'NB13760629', 'NB13760630'), 'PD63383',
  sample_ID %in% c('NB14599068', 'NB14599069', 'NB14599070'), 'PD63383_nuclear'
))]

# load list of final mutations used for the phylogeny
muts_dt = data.table(read.csv('Data/20241114_599muts_QCJbrowse.csv', header = T))
muts_final = muts_dt[Jbrowse.quality == 'Y', mut_ID] %>% unlist()
paste('Number of mutations that passed QC:', length(muts_final)) # 255 

muts_final_dt = data.table(muts_final)
setnames(muts_final_dt, 'muts_final', 'mut_ID')
muts_final_dt[, mut_ID0 := substr(mut_ID, 1, nchar(mut_ID)-4)]

# basic checks on the data
paste('Number of cells with detected read spanning mutation position:', length(ai_counts_dt[, barcode] %>% unlist() %>% unique()))
# 716 

# what exactly is the Tot column?
hist(ai_counts_dt[, Tot], xlab = 'Number of reads', main = 'Number of mutation-spanning reads in a cell')

# identify mutation ID 
ai_counts_dt[, Chrom := paste0('chr', chr)] # set the same chromosome name as mut_ID 
ai_counts_dt[, mut_ID0 := paste0(Chrom, '_', pos)]
ai_counts_dt = merge(ai_counts_dt, muts_final_dt, by = 'mut_ID0')

# number of mutations identified
paste('Number of mutation-spanning positions identified:', length(ai_counts_dt[, mut_ID] %>% unlist() %>% unique()))
# 93 so this is actually quite a lot 

# For each read, determine if it is mutant or wt
ai_counts_dt[, Ref := tstrsplit(mut_ID, '_', fixed=TRUE, keep=3)]
ai_counts_dt[, Alt := tstrsplit(mut_ID, '_', fixed=TRUE, keep=4)]
mut_cols = c('A', 'C', 'G', 'T')
ai_counts_dt[, Read := apply(ai_counts_dt[, ..mut_cols], 1, function(row) {
  read = mut_cols[which(row > 0)]
  paste(read, collapse = ', ')
})]
ai_counts_dt[, status := factor(fcase(
  Read == Ref, 'wt',
  Read == Alt, 'mutant'
))]

table(ai_counts_dt[, status])
# 25 mutant reads
# 730 wt reads

# what mutations are those?
# read in the dataframe with mutation group assignment
muts_assignment = fread('Data/20241205_muts_classes_255.csv', sep = ',', header = TRUE)
ai_counts_dt = merge(ai_counts_dt, muts_assignment, by = 'mut_ID')
table(ai_counts_dt[, mut_class])

#in which clusters are there cells with those reads? what about cells with mutant reads specifically?




