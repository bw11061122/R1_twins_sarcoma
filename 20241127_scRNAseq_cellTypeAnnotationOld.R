###################################################################################################################################
# Cell type annotation 

# One possible option is to look at expression of known markers

# I found https://www-nature-com.ezp.lib.cam.ac.uk/articles/s41388-024-03001-8 (March 2024, Oncogene)
# which reports the scRNA-seq of undifferentiated pleomorphic sarcoma (which is kind of what I have; pleomorphic means occurring in different forms)
# they identify: monocytes (macrophages), tumour cells, fibroblast, T cell, endothelial, mast, NK, pericyte cells
# perhaps do this for tumour from each twin separately and then can compare clusters 

# run pre-processing on aggregated tumour from one twin
PD62341.d = sapply(PD62341_ids, function(i){
  d10x = Read10X(file.path(data.loc, "Data", i, "filtered_feature_bc_matrix"))
  colnames(d10x) = paste(sapply(strsplit(colnames(d10x),split="-"),'[[',1L),i,sep="-")
  d10x
})
PD62341.data = do.call("cbind", PD62341.d) # that just treats this as having > 20k cells - is this the right way to do this?

PD62341t = CreateSeuratObject(counts = PD62341.data, project = "PD62341.sarcoma", names.field = 1, names.delim = "_", meta.data = NULL)
PD62341t[["percent.mt"]] = PercentageFeatureSet(PD62341t, pattern = "^MT-")
PD62341t@meta.data$CellID = rownames(PD62341t@meta.data) # add the cellID to the metadata
meta = PD62341t@meta.data # create a dt with all the metadata of interest, including cell ID 
keeprows = meta[meta$nCount_RNA > 1000 & meta$nCount_RNA < 30000,] # require number of UMIs to be higher than minimum (1000 UMIs)
keeprows = keeprows[keeprows$nFeature_RNA > 300 & meta$nFeature_RNA < 1000,] # require number of features to be higher than minimum (300 genes)
keeprows = keeprows[keeprows$percent.mt <= 10,] # exclude cells with mitochondrial content above a threshold (10%)
cellids_keep = rownames(keeprows) # which CellIDs to keep 
PD62341t_filtered = subset(PD62341t, subset = CellID %in% cellids_keep) # select cells of interest
keepfeatures=rownames(PD62341t_filtered)[which(!rownames(PD62341t_filtered) %in% genes_to_exclude)]
PD62341t_filtered = subset(PD62341t_filtered, features=keepfeatures)
PD62341t_filtered = NormalizeData(object = PD62341t_filtered, normalization.method = "LogNormalize", scale.factor = 10000) # default 1e4
PD62341t_filtered = FindVariableFeatures(PD62341t_filtered, selection.method = "vst", nfeatures = 2000) # default 2k 
PD62341t_filtered = ScaleData(PD62341t_filtered, vars.to.regress = "percent.mt") # Data scaling performed by default on the 2,000 variable features, some tutorials recommend regressing out mt content  
PD62341t_filtered = RunPCA(PD62341t_filtered, npcs = 50, features=VariableFeatures(object = PD62341t_filtered), verbose = F) # npcs = number of PCs to compute and store (50 by default)
PD62341t_filtered = FindNeighbors(PD62341t_filtered, dims = 1:50, verbose = T)
PD62341t_filtered = FindClusters(PD62341t_filtered, resolution = 0.4, verbose = T)
PD62341t_filtered = RunUMAP(PD62341t_filtered, dims = 1:50, verbose = T, min.dist = 0.6, n.neighbors = 40)

pdf(glue('Figures/F6/20241208_umap_PD62341t_rmGenes_aggTumour.pdf'), height = 4, width = 4)
DimPlot(PD62341t_filtered, reduction = "umap", label = TRUE)
dev.off() 

# first, start by looking at expression of each of the possible marker genes 
marker_genes = c("PTPRC", "CD68", "CD14", "CD163", "AIF1", # monocyte 
                 "COL11A1", "COL3A1", # fibroblasts 
                 "CD8A", "CD4", "CD3D", # T cells
                 "GNLY", "NKG7", # NK cells
                 "PLVAP", "PECAM1", "VWF", # endothelial
                 "ACTA2", "NOTCH3", "RGS5", # pericyte 
                 "ZNF341", "MN1", # tumour cells 
                 "TMEM119", "C1Q", "PROS1", "IBA1")  

pdf(glue('Figures/F6/20241208_featurePlot_marker_genes_tPD62341_agg.pdf'), width=12, height=16)
FeaturePlot(PD62341t_filtered, features = marker_genes, order=TRUE)
dev.off()

###################################################################################################################################
# Cell type annotation (SEURAT TUTORIAL) - I TRIED BUT DOESN'T SEEM TO WORK WITH THE FETAL DATASET

# another approach is to maybe try and do data integration with Seurat 
# library(devtools)
# devtools::install_github('satijalab/seurat-data')
library(SeuratData)
AvailableData() # what kinds of datasets would be useful for this sample? some sarcoma? any fetal dataset?

# install a reference package and load the data from the reference 
InstallData("panc8")
panc8 = LoadData("panc8")
pancreas.ref = subset(panc8, tech %in% c("celseq2", "smartseq2")) # select scRNA-seq generated with Cel-seq and SMART-seq
pancreas.ref[["RNA"]] = split(pancreas.ref[["RNA"]], f = pancreas.ref$tech)

# data pre-processing (no integration)
pancreas.ref = NormalizeData(pancreas.ref)
pancreas.ref = FindVariableFeatures(pancreas.ref)
pancreas.ref = ScaleData(pancreas.ref)
pancreas.ref = RunPCA(pancreas.ref)
pancreas.ref = FindNeighbors(pancreas.ref, dims = 1:30)
pancreas.ref = FindClusters(pancreas.ref)
pancreas.ref = RunUMAP(pancreas.ref, dims = 1:30)
DimPlot(pancreas.ref, group.by = c("tech", "celltype"))

# integrate datasets from different sequencing tech into a shared reference
pancreas.ref = IntegrateLayers(object = pancreas.ref, method = CCAIntegration, orig.reduction = "pca",
                               new.reduction = "integrated.cca", verbose = FALSE)
pancreas.ref = FindNeighbors(pancreas.ref, reduction = "integrated.cca", dims = 1:30)
pancreas.ref = FindClusters(pancreas.ref)
pancreas.ref = RunUMAP(pancreas.ref, reduction = "integrated.cca", dims = 1:30)
DimPlot(pancreas.ref, group.by = c("tech", "celltype"))

# create a query dataset (two other technologies)
pancreas.query = subset(panc8, tech %in% c("fluidigmc1", "celseq"))
pancreas.query = NormalizeData(pancreas.query)
pancreas.anchors = FindTransferAnchors(reference = pancreas.ref, query = pancreas.query, dims = 1:30,
                                       reference.reduction = "pca")
predictions = TransferData(anchorset = pancreas.anchors, refdata = pancreas.ref$celltype, dims = 1:30)
pancreas.query = AddMetaData(pancreas.query, metadata = predictions)
pancreas.query$prediction.match = pancreas.query$predicted.id == pancreas.query$celltype
table(pancreas.query$prediction.match)

# use reference dataset to assign cell type to the query
pancreas.ref = RunUMAP(pancreas.ref, dims = 1:30, reduction = "integrated.cca", return.model = TRUE)
pancreas.query = MapQuery(anchorset = pancreas.anchors, reference = pancreas.ref, query = pancreas.query,
                          refdata = list(celltype = "celltype"), reference.reduction = "pca", reduction.model = "umap")
p1 = DimPlot(pancreas.ref, reduction = "umap", group.by = "celltype", label = TRUE, label.size = 3,
             repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
p2 = DimPlot(pancreas.query, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
             label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
p1 + p2

# maybe we could try with fetal data as reference? could identify immune cells / fibroblasts / endothelial
remotes::install_github('satijalab/azimuth', ref = 'master')
InstallData("fetusref")
fetal = LoadData("fetusref", "azimuth") 
# I am getting an error where fetusref is installed but cannot be loaded 
# attempted solving by installing Azimuth package (also Satija lab) but doesn't seem like I can get all the packages reasonably easily so I think I'm going to drop this for now
fetal.ref = subset(fetal) # select scRNA-seq generated with Cel-seq and SMART-seq
fetal.ref[["RNA"]] = split(fetal.ref[["RNA"]], f = fetal.ref$tech)

# data pre-processing (no integration)
fetal.ref = NormalizeData(fetal.ref)
fetal.ref = FindVariableFeatures(fetal.ref)
fetal.ref = ScaleData(fetal.ref)
fetal.ref = RunPCA(fetal.ref)
fetal.ref = FindNeighbors(fetal.ref, dims = 1:30)
fetal.ref = FindClusters(fetal.ref)
fetal.ref = RunUMAP(fetal.ref, dims = 1:30)
DimPlot(fetal.ref, group.by = c("tech", "celltype"))

# integrate datasets from different sequencing tech into a shared reference
fetal.ref = IntegrateLayers(object = fetal.ref, method = CCAIntegration, orig.reduction = "pca",
                            new.reduction = "integrated.cca", verbose = FALSE)
fetal.ref = FindNeighbors(fetal.ref, reduction = "integrated.cca", dims = 1:30)
fetal.ref = FindClusters(fetal.ref)
fetal.ref = RunUMAP(fetal.ref, reduction = "integrated.cca", dims = 1:30)
DimPlot(fetal.ref, group.by = c("tech", "celltype"))


# Cell type annotation based on markers of each cluster 

# Let's agree that we are using clusters with the following parameters:
# N_PCs = 70, resolution = 0.2, min_dist = 0.3, n_neighbors = 30
PD62341t = CreateSeuratObject(counts = PD62341.data, project = "PD62341.sarcoma", names.field = 1, names.delim = "_", meta.data = NULL)
PD62341t[["percent.mt"]] = PercentageFeatureSet(PD62341t, pattern = "^MT-")
PD62341t@meta.data$CellID = rownames(PD62341t@meta.data) # add the cellID to the metadata
meta = PD62341t@meta.data # create a dt with all the metadata of interest, including cell ID 
keeprows = meta[meta$nCount_RNA > 1000,] # require number of UMIs to be higher than minimum (1000 UMIs)
keeprows = keeprows[keeprows$nFeature_RNA > 300,] # require number of features to be higher than minimum (300 genes)
keeprows = keeprows[keeprows$percent.mt <= 10,] # exclude cells with mitochondrial content above a threshold (10%)
cellids_keep = rownames(keeprows) # which CellIDs to keep 
PD62341t_filtered = subset(PD62341t, subset = CellID %in% cellids_keep) # select cells of interest
keepfeatures=rownames(PD62341t_filtered)[which(!rownames(PD62341t_filtered) %in% genes_to_exclude)]
PD62341t_filtered = subset(PD62341t_filtered, features=keepfeatures)
PD62341t_filtered = NormalizeData(object = PD62341t_filtered, normalization.method = "LogNormalize", scale.factor = 10000) # default 1e4
PD62341t_filtered = FindVariableFeatures(PD62341t_filtered, selection.method = "vst", nfeatures = 2000) # default 2k 
PD62341t_filtered = ScaleData(PD62341t_filtered, vars.to.regress = "percent.mt") # Data scaling performed by default on the 2,000 variable features, some tutorials recommend regressing out mt content  
PD62341t_filtered = RunPCA(PD62341t_filtered, npcs = 70, features=VariableFeatures(object = PD62341t_filtered), verbose = F) # npcs = number of PCs to compute and store (50 by default)
PD62341t_filtered = FindNeighbors(PD62341t_filtered, dims = 1:70, verbose = T)
PD62341t_filtered = FindClusters(PD62341t_filtered, resolution = 0.2, verbose = T)
PD62341t_filtered = RunUMAP(PD62341t_filtered, dims = 1:70, verbose = T, min.dist = 0.3, n.neighbors = 30)
DimPlot(PD62341t_filtered, reduction = "umap", label = TRUE) # okay so we have 14 clusters for the moment 

pdf('Figures/F6/20241208_umap_PD62341t_14clusters.pdf')
DimPlot(PD62341t_filtered, reduction = "umap", label = TRUE) # okay so we have 14 clusters for the moment 
dev.off()

# identify top 5 genes for each cluster 
PD62341t_filtered.markers = FindAllMarkers(PD62341t_filtered, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top4 = PD62341t_filtered.markers %>% group_by(cluster) %>% top_n(n = 4, wt = avg_log2FC)
DotPlot(PD62341t_filtered, features = top4$gene) + coord_flip()

top10 = PD62341t_filtered.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

pdf('Figures/F6/20241208_dotplot_PD62341t_14clusters.pdf', height=16, width=10)
DotPlot(PD62341t_filtered, features = top4$gene) + coord_flip()
dev.off()

# also do a feature plot as this will help identify the easiest clusters (tumour, fibroblasts)
pdf(glue('Figures/F6/20241208_featurePlot_knownMarkers_tPD62341_agg_14clusters.pdf'), width=12, height=16)
FeaturePlot(PD62341t_filtered, features = marker_genes, order=TRUE)
dev.off()

# checking genes for each cluster 
# 0 - looks like something from the bone marrow. could be macrophages? weirdly also has some brain-related stuff. but also, very low PTPRC expression?
# 1 - tumour 
# 2 - endothelial 
# 3 - looks like some cells from the bone marrow 
# 4 - tumour 
# 5 - looks like meyloid cells (+ immune subset)
# 6 - fibroblasts?
# 7 - express PTPRC and some immune-enriched stuff, so likely another myeloid set
# 8 - fibroblasts?
# 9 - tumour 
# 10 - fibroblasts?
# 11 - quite unclear 
# 12 - quite unclear 
# 13 - looks endothelial

# I am not sure what clusters 0 and 3 are 
# Find all markers for those two clusters
cluster0.markers <- FindMarkers(PD62341t_filtered, ident.1 = 0)
# based on the stuff that's up (HBA1, HBA2, HGB1, ALAS2, TMCC2) - feels like this has to be erythrocytes
cluster3.markers <- FindMarkers(PD62341t_filtered, ident.1 = 3)
# also feels like this has to be erythroid cells 

# can we use markers for blood cells?
# PBMC markers from https://satijalab.org/seurat/articles/pbmc3k_tutorial.html 
pbmc_markers = c("CCR7", # naive T cell (CD4)
                 "CD14", "LYZ", # monocyte 
                 "S100A4", # memory T cell (CD4)
                 "CD8A", # CD8+ T cell 
                 "CD14", "CD68", # myeloid-like stuff I guess
                 "FCGR3A", "MS4A7", # monocyte 
                 "GNLY", "NKG7", # NK
                 "FCER1A", "CST3", # DC
                 "CD34", "RUNDC3A", "TMEM119",
                 "CD45",
                 "HBA1", "HBA2", "HBB", "ALAS2", "AHSP") # hemoglobin 

pdf(glue('Figures/F6/20241208_featurePlot_pbcmMarkers_tPD62341_agg.pdf'), width=12, height=12)
FeaturePlot(PD62341t_filtered, features = pbmc_markers, order=TRUE)
dev.off()

# Another way of doing this is what I found here: https://bioconductor.org/packages/release/bioc/vignettes/SingleR/inst/doc/SingleR.html

###################################################################################################################################
# Logistic regression on HCA 

# the first thing I can do is probably check the expression of some adrenal markers in my dataset
# markers from https://pmc.ncbi.nlm.nih.gov/articles/PMC10443814/#sd 
adrenal_markers = c("NR5A1", "CYP11A1", "MC2R",
                    "STAR", "CLRN1", "MIR202HG", "FAM166B", "CYP11B1", "MT3", 
                    "CCN3", "SULT2A1") # CCN3 is a DZ (definitive zone) marker, SULT2A1 is an FZ (fetal zone) marker 

pdf(glue('Figures/F6/20241208_featurePlot_adrenalMarkers_tPD62341_agg.pdf'), width=12, height=8)
FeaturePlot(PD62341t_filtered, features = adrenal_markers, order=TRUE)
dev.off()

# generally I feel like I am getting too many clusters
# tumour clearly expresses some adrenal markers (like CYP11A1 and SULT2A1)
# however, the 0-3 cluster also seems to express some of those 

# Logistic regression on HCA 
# the way to do this is described in Matthew Young's paper here:
# https://www.science.org/doi/10.1126/science.aat1699

# how this is done in the paper:
# To measure the similarity of tumor cells to either normal mature kidney cells or fetal cells,
# we trained a logistic regression model using elastic net regularization on the cellular identities
# defined by the clusters in normal mature and fetal epithelial and vascular compartments. In
# training this model we set alpha=0.99 to produce strong regularization but to prevent the
# exclusion of strongly co-linear genes (would be good for you to check what they mean by regularization)

# To obtain regression coefficients specific to each cluster in our training data we fit a series
# of N binomial logistic regression models, where N is the number of clusters in the training data
# (i.e., one-versus-rest binomial logistic regression). To prevent the observed frequencies of cells
# (which we do not expect to accurately reflect the true abundances in situ) from biasing the
# regression coefficients we use an offset for each model given by log(f / (1-f))
# where f is the fraction of cells in the cluster being trained.

# In each case, we performed 10-fold cross validation and selected the regularization co-
# efficient, lambda, to be as large as possible (i.e., as few non-zero coefficients as possible) such
# that the cross validated accuracy was within 1 standard deviation of the minimum.
# These models were then used to calculate a predicted similarity to each of the fetal and
# normal clusters for each cell in the tumor map. In calculating the predicted values, an offset of 0
# was used. Softmax normalization was not used to allow for the possibility that tumor cells do
# not resemble any of the fetal or normal cells in the training set. Predicted logits were then
# averaged within each tumor cluster and converted to probabilities for visualization.

# required libraries
library(glmnet)
library(ggplot2)
library(cowplot)
library(foreach)
library(doMC)
BiocManager::install(c("monocle"))
BiocManager::install(c('DelayedArray', 'DelayedMatrixStats', 'org.Hs.eg.db', 'org.Mm.eg.db'))
devtools::install_github("cole-trapnell-lab/garnett")

# Run logistic regression analysis 
# code from https://github.com/constantAmateur/scKidneyTumors/blob/master/postPipeline.R

# define required functions (as per https://github.com/constantAmateur/scKidneyTumors/blob/master/similarity.R)
# get offset (what does this actually mean?)
getPopulationOffset = function(y){
  if(!is.factor(y)) # convert the input to factor level 
    y=factor(y)
  if(length(levels(y))!=2)
    stop("y must be a two-level factor")
  off = sum(y==levels(y)[2])/length(y)
  off = log(off/(1-off))
  return(rep(off,length(y)))
}

# Do the OvR fit for every variable.  This just does a simple CV selection of regularisation amount.  Far from ideal, but should be good enough for the main conclusions.
multinomialFitCV = function(x,y,nParallel=1,...){
  fits = list()
  if(nParallel>1)
    registerDoMC(cores=nParallel)
  #Do them in order of size
  marks = names(sort(table(as.character(y))))
  for(mark in marks){
    message(sprintf("Fitting model for variable %s",mark))
    fac = factor(y==mark)
    #The two main modes of failure are too few positives and errors constructing lambda.  These should be handled semi-gracefully
    fits[[mark]] = tryCatch(
      cv.glmnet(x,fac,offset=getPopulationOffset(fac),family='binomial',intercept=FALSE,alpha=0.99,nfolds=10,type.measure='class',parallel=nParallel>1,...),
      error = function(e) {
        tryCatch(
          cv.glmnet(x,fac,offset=getPopulationOffset(fac),family='binomial',intercept=FALSE,alpha=0.99,nfolds=10,type.measure='class',parallel=nParallel>1,lambda=exp(seq(-10,-3,length.out=100)),...),
          error = function(e) {
            warning(sprintf("Could not fit model for variable %s",mark))
            return(NULL)
          })
      })
  }
  return(fits)
}

# Load training data
loadTrainingData = function(dirs){
  mDat = NULL
  toc = NULL
  for(dir in dirs){
    metadata = readRDS(file.path(dir,'metadata.RDS'))
    ttoc = readRDS(file.path(dir,'tableOfCounts.RDS'))
    rownames(metadata$meta.data) = paste0(metadata$meta.data$Sanger_study_ID,'___',gsub('^[0-9]+_','',rownames(metadata$meta.data)))
    colnames(ttoc) = rownames(metadata$meta.data)
    if(is.null(mDat)){
      mDat = metadata$meta.data
    }else{
      mDat = rbind(mDat,metadata$meta.data)
    }
    if(is.null(toc)){
      toc = ttoc
    }else{
      toc = cbind(toc,ttoc)
    }
  }
  return(list(toc=toc,mDat=mDat))
}

# where do I get the data from??
# I am considering possibly this https://descartes.brotmanbaty.org/bbi/human-gene-expression-during-development/
# (you can try and get the .Rds data files for the different organs that they have)

adrenal_rds = readRDS('Data/adrenal_classifier.Rds')


# read the data 
plotDir='/Results'
tfs = read.table('~/Projects/Common/Homo_sapiens_transcription_factors_gene_list.txt',sep='\t',header=TRUE)
rtoc = readRDS('~/scratch/KidneySC/globalRawTableOfCounts.RDS')

# exclude the following genes
excludeGenes='^(EIF[0-9]|RPL[0-9]|RPS[0-9]|RPN1|POLR[0-9]|SNX[0-9]|HSP[AB][0-9]|H1FX|H2AF[VXYZ]|PRKA|NDUF[ABCSV]|ATP[0-9]|PSM[ABCDEFG][0-9]|UBA[0-9]|UBE[0-9]|USP[0-9]|TXN)'
coreExcludeGenes = unique(c(grep('\\.[0-9]+_',rownames(rtoc),value=TRUE), # poorly characterized (what does this mean??)
                            grep('MALAT1',rownames(rtoc),value=TRUE), # contamination - why are we excluding this? how do we know this is contamination
                            grep('^HB[BGMQDAZE][12]?_',rownames(rtoc),value=TRUE), # contamination (apparently you should filter out blood markers)
                            grep('^MT-',rownames(rtoc),value=TRUE), # mt genes  
                            grep(excludeGenes,rownames(rtoc),value=TRUE) # housekeeping genes
))


# Training models 

###################
# Foetal Clusters
trainDat = loadTrainingData(file.path(plotDir,'barcodeLevelGroup_foetalEpitheliumAndVascularV3/batchCorrected/'))
trainDat$mDat$Trainer = trainDat$mDat$res.1

#Add in the MAST cells
extras = sim$loadTrainingData(file.path(plotDir,'barcodeLevelGroup_tumourImmune/batchCorrected/'))
#Get the cells and their class
cells = c(rownames(trainDat$mDat),rownames(extras$mDat[extras$mDat$res.1=='16',]))
classes = c(trainDat$mDat$res.1,rep('MAST',sum(extras$mDat$res.1=='16')))
#Any genes to exclude go here
excludeGenes=coreExcludeGenes
#Any genes that we wish to give extra weight should go here
includeGenes=c()

#Get and normalise the data
dat = Seurat::LogNormalize(rtoc[,cells])
dat = dat[(rownames(dat) %in% includeGenes) | (Matrix::rowSums(dat>0)>3 & !(rownames(dat)%in%excludeGenes)),]
dat = t(dat)
fitFoetalClusters = sim$multinomialFitCV(dat,classes,nParallel=10)
#Get the genes that are going to be informative
saveRDS(fitFoetalClusters,file.path(plotDir,'lrFoetalClustersV4.RDS'))

###################
# Normal Clusters
trainDat = loadTrainingData(file.path(plotDir,'barcodeLevelGroup_normalEpitheliumAndVascularWithoutProximalTubularV2/batchCorrected/'))
#Label and merge clusters

#annotMap = c(EN17='G',EN0='PT1',EN7='PT1',EN6='PT2',EN5='PT2',EN4='PT2',EN9='PT2',EN3='PT2',EN18='PT3',EN19='D',EN21='C1',EN22='C2',EN20='P',EN15='U1/2',EN25='U3/4',EN23='M',EN24='F',EN16='GE',EN12='DV',EN10='AV')
#Drop the ambiguous clusters and merge the ones that should really be the same
#trainDat$mDat = trainDat$mDat[paste0('EN',trainDat$mDat$res.1) %in% names(annotMap),]
#trainDat$mDat$Trainer = annotMap[paste0('EN',trainDat$mDat$res.1)]
trainDat$mDat$Trainer = trainDat$mDat$res.1

#Add in the MAST cells
extras = sim$loadTrainingData(file.path(plotDir,'barcodeLevelGroup_tumourImmune/batchCorrected/'))

#Get the cells and their class
cells = c(rownames(trainDat$mDat),rownames(extras$mDat[extras$mDat$res.1=='16',]))
classes = c(trainDat$mDat$res.1,rep('MAST',sum(extras$mDat$res.1=='16')))

#Any genes to exclude go here
excludeGenes=coreExcludeGenes

#Any genes that we wish to give extra weight should go here
includeGenes=c()

#Get and normalise the data
dat = Seurat::LogNormalize(rtoc[,cells])
dat = dat[(rownames(dat) %in% includeGenes) | (Matrix::rowSums(dat>0)>3 & !(rownames(dat)%in%excludeGenes)),]
dat = t(dat)
fitNormalClusters = multinomialFitCV(dat,classes,nParallel=10)

#Get the genes that are going to be informative
saveRDS(fitNormalClusters,file.path(plotDir,'lrNormalClustersWithControl.RDS'))

###################################################################################################################################
# Myeloid cells - what's up with those? 
# where are the myeloid cells and what can I say about the expression of CD14 and CD68?

# How do I get reads from those cells to do the genotyping?
