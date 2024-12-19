SCRIPTS 
Barbara Walkowiak bw18

1 MAIN SCRIPTS 
- F1 20241011_pileup_filters.R - Filtering script (up until 225 muts)
- F2 20241101_tumourInfiltration.R - Estimation of tumour infiltration for each sample
- F3 20241119_phylogenyFinalSet.R - Reconstruction of the phylogeny 
- F4 20241113_twinTransfusion.R - Twin-twin transfusion estimates 
- F5 20241204_ABC_grid.R - ABC simulations 
- F6 20241127_scRNAseq.R - scRNA-seq analysis (basic QC and clustering)
- F7 20241204_scRNAseq_genotyping.r - scRNA-seq genotyping (Allele Integrator)

2 OTHER IMPORTANT SCRIPTS
- 20241011_removeheader.sh – short bash script to remove headers from out files from the pileup 
- 20241015_pileup_checks.R – check that everything went okay with the 2024/10/15 pileup and merge into a .tsv file (pileup_merged_20241026.tsv)
- 20241025_cpgVAF.txt – running low quality pileup on the farm (describes how I did this) (I didn’t run the high quality pileup script)
- 20241025_pileup_checks.R – check that everything went okay with the 2024/10/25 pileup and merge into a .tsv file (pileup_merged_20241025.tsv)
- 20241027_driver_analysis.R (searched for mutations (subs, indels, chromosomal rearrangements) in the dataset to see if any would be likely additional drivers; we did not find anything convincing) 
- 20241113_telomere_length_bsub.sh – get telomere length files (submission)
- 20241113_telomere_length.sh – get telomere length files (script)
- 20241113_twinTransfusionTelomeres.R – analyse telomere length files (we did this to see if there’s any evidence of variable proliferation between PD62341 / PD63383 cells, especially in the spleen sample)
- 20241119_removeheader_brass.sh – short bash script to remove header from BRASS output files (.bedpe)
- 20241122_chr1_duplicates.R (tried to look at SAM flags of all reads spanning the position of chr1;38827952 C>A mutation to see if higher rate of duplicates could explain the inflated VAF)
- 20241123_graphics.R (pie charts plots for phylogenetic trees with shaded VAF)
- 20241127_ABC.R – this was my attempt to do the simulations by simulating cell division and embryo splitting on a tree, moved to a grid-based approach  
- AlleleIntegrator_scRNAgenotyping.r – script to run Allele Integrator on the farm 
- AlleleIntegrator_scRNAgenotyping.bsub – AlleIe Integrator submission to the farm
- chr1_388_cutFile.sh – creates a shorter version of the file with chr1;38827952 C>A reads for each sample (_short.txt file) 
- split_chr1file.sh – split chr1;38827952 C>A reads by sample 
- view_dups_chr1_388.bsub – get all duplicate reads spanning the mutation chr1;38827952 C>A using samtools view (submission)
- view_dups_chr1_388.sh – get all duplicate reads spanning the mutation chr1;38827952 C>A using samtools view (script)
- README1.txt – I keep a readme file for running stuff on the farm and transferring scripts to have a record of what I was doing, when and why. Does not contain any actual analysis

3 PACKAGES VERSIONS
R version 4.4.1 (2024-06-14)
Seurat_5.1.0                      
SeuratObject_5.0.2                                   
abc_2.2.1                                                              
plyr_1.8.9                        
ggrepel_0.9.6                    
pheatmap_1.0.12                   
BSgenome.Hsapiens.UCSC.hg38_1.4.5 
BSgenome_1.72.0                  
rtracklayer_1.64.0                
BiocIO_1.14.0                     
GenomicRanges_1.56.2             
Biostrings_2.72.1                 
GenomeInfoDb_1.40.1               
XVector_0.44.0                   
IRanges_2.38.1                    
S4Vectors_0.42.1                  
BiocGenerics_0.50.0              
glue_1.8.0                       
purrr_1.0.2                       
tidyr_1.3.1                      
tibble_3.2.1                      
tidyverse_2.0.0                   
ggplot2_3.5.1                    
dplyr_1.1.4                       
data.table_1.16.2                




