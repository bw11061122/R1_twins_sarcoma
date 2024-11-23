###################################################################################################################################
# SCRIPT 6

# Script to analyse the pileup (high quality, run 15/10/2024)
# 2024-11-22
# Barbara Walkowiak bw18

# INPUT: 
# 1 sam file with reads spanning chr138827952_C_A

# OUTPUT:
# 1 try to figure out if any reads are possible PCR duplicates 

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
list_sam_PD62341 = paste0('Data/chr1_388_PD62341/', list.files(path = "Data/chr1_388_PD62341/", pattern = "*.txt"))
sam_PD62341 = do.call(bind_rows, lapply(list_sam_PD62341, function(x) read.table(x, sep='\t', fill = TRUE, header=FALSE)))
sam_PD62341 = sam_PD62341[, c(1:10)]
list_sam_PD63383 = paste0('Data/chr1_388_PD63383/', list.files(path = "Data/chr1_388_PD63383/", pattern = "*.txt"))
sam_PD63383 = do.call(bind_rows, lapply(list_sam_PD63383, function(x) read.table(x, sep='\t', fill = TRUE, header=FALSE)))
sam_PD63383 = sam_PD63383[, c(1:10)]

PD62341aa = read.table('Data/chr1_388_reads_PD62341aa.txt', fill=TRUE)
PD62341ak = read.table('Data/chr1_388_reads_PD62341ak.txt', fill=TRUE)
PD62341aa = PD62341aa[, c(1:10)]
PD62341ak = PD62341ak[, c(1:10)]

sam = rbind(sam_PD62341, PD62341aa, PD62341ak, sam_PD63383)

