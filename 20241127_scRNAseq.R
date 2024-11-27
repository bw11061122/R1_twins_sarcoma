###################################################################################################################################
# SCRIPT 8

# Script to remove germline mutations in copy number-altered regions
# 2024-11-27
# Barbara Walkowiak bw18

# INPUT: 
# 1 PILEUP dt (.cnv.segment) for each sample 

# OUTPUT:
# 1 list of mutations to exclude based on presence in copy number germline variants 

# Script to identify regions of copy number alterations (deletions or duplications) in the germline
# Identify mutations in those regions where the VAF is 0.3 and so those can be explained as a germline copy number alteration 

###################################################################################################################################
# LIBRARIES
library(Seurat)
library(ggplot2)
library(kableExtra)

###################################################################################################################################
# INPUT DATA 

setwd('/Users/bw18/Desktop/1SB')

experiment_name = "Covid Example"
dataset_loc <- "./"
ids <- c("PBMC2", "PBMC3", "T021PBMC", "T022PBMC")

d10x.metrics <- lapply(ids, function(i){
  # remove _Counts is if names don't include them
  metrics <- read.csv(file.path(dataset_loc,paste0(i,"_Counts/outs"),"metrics_summary.csv"), colClasses = "character")
})
experiment.metrics <- do.call("rbind", d10x.metrics)
rownames(experiment.metrics) <- ids

sequencing_metrics <- data.frame(t(experiment.metrics[,c(4:17,1,18,2,3,19,20)]))

row.names(sequencing_metrics) <- gsub("\\."," ", rownames(sequencing_metrics))