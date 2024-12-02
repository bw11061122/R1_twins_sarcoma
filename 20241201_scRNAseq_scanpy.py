###################################################################################################################################
# SCRIPT 8

# Script to remove germline mutations in copy number-altered regions
# 2024-11-27
# Barbara Walkowiak bw18

# INPUT: 
# 1 scRNAseq CellRanger output 

# OUTPUT: 
# 1 basic scRNA-seq analysis (just teaching myself how to do this)

# Note: I am getting started by following this tutorial: 
# https://scanpy.readthedocs.io/en/stable/tutorials/basics/clustering.html

# %%
###################################################################################################################################

# Scanpy scRNA-seq analysis 

###################################################################################################################################
# LIBRARIES
import scanpy as sc
import anndata as ad

###################################################################################################################################
# INPUT DATA 

samples = {
    "s1d1": "s1d1_filtered_feature_bc_matrix.h5",
    "s1d3": "s1d3_filtered_feature_bc_matrix.h5",
}
adatas = {}

for sample_id, filename in samples.items():
    path = EXAMPLE_DATA.fetch(filename)
    sample_adata = sc.read_10x_h5(path)
    sample_adata.var_names_make_unique()
    adatas[sample_id] = sample_adata

adata = ad.concat(adatas, label="sample")
adata.obs_names_make_unique()
print(adata.obs["sample"].value_counts())
adata

###################################################################################################################################
# BASIC QC (NR COUNTS, NR FEATURES, MT GENES)

# mitochondrial genes, "MT-" for human, "Mt-" for mouse
adata.var["mt"] = adata.var_names.str.startswith("MT-")
# ribosomal genes
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
# hemoglobin genes
adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")

sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True
)

sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.4,
    multi_panel=True,
)

sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")

sc.pp.filter_cells(adata, min_genes=100)
sc.pp.filter_genes(adata, min_cells=3)

# %%
###################################################################################################################################
# DOUBLET DETECTION 
sc.pp.scrublet(adata, batch_key="sample")

# %%
###################################################################################################################################
# NORMALIZATION 
# Saving count data
adata.layers["counts"] = adata.X.copy()
# Normalizing to median total counts
sc.pp.normalize_total(adata)
# Logarithmize the data
sc.pp.log1p(adata)

# %%
###################################################################################################################################
# FEATURE SELECTION
# reproduces methods from Seurat / CellRanger / Seurat v3
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="sample")
sc.pl.highly_variable_genes(adata)

# %%
###################################################################################################################################
# DIMENSIONALITY REDUCTION

# PCA 
sc.tl.pca(adata)
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True)
sc.pl.pca(
    adata,
    color=["sample", "sample", "pct_counts_mt", "pct_counts_mt"],
    dimensions=[(0, 1), (2, 3), (0, 1), (2, 3)],
    ncols=2,
    size=2,
)

# Compute the neighbouhood graph 
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.pl.umap(
    adata,
    color="sample",
    # Setting a smaller point size to get prevent overlap
    size=2,
)

# %%
###################################################################################################################################
# CLUSTERING

# Using the igraph implementation and a fixed number of iterations can be significantly faster, especially for larger datasets
sc.tl.leiden(adata, flavor="igraph", n_iterations=2)
sc.pl.umap(adata, color=["leiden"])

# %%
###################################################################################################################################
# CLUSTERING to re-assess QC metrics 

sc.pl.umap(
    adata,
    color=["leiden", "predicted_doublet", "doublet_score"],
    # increase horizontal space between panels
    wspace=0.5,
    size=3,
)

sc.pl.umap(
    adata,
    color=["leiden", "log1p_total_counts", "pct_counts_mt", "log1p_n_genes_by_counts"],
    wspace=0.5,
    ncols=2,
)

# %%
###################################################################################################################################
# CLUSTERING for cell type annotation 

# test different resolution parameters 
for res in [0.02, 0.5, 2.0]:
    sc.tl.leiden(
        adata, key_added=f"leiden_res_{res:4.2f}", resolution=res, flavor="igraph"
    )

sc.pl.umap(
    adata,
    color=["leiden_res_0.02", "leiden_res_0.50", "leiden_res_2.00"],
    legend_loc="on data",
)

# %%
###################################################################################################################################
# MANUAL CELL TYPE ANNOTATION WITH EXPECTED MARKERS 

# Which cell types do I actually expect in this data?
# I found https://www-nature-com.ezp.lib.cam.ac.uk/articles/s41388-024-03001-8 (March 2024, Oncogene)
# which reports the scRNA-seq of undifferentiated pleomorphic sarcoma (which is kind of what I have; pleomorphic means occurring in different forms)
# they identify: monocytes (macrophages), tumour cells, fibroblast, T cell, endothelial, mast, NK, pericyte cells

# Define a set of marker genes that you expect in the dataset
# I am doing this based on this paper, which is okay for the moment but I would need to do more reading about this 
# how many genes should you use per cell cluster actually??
marker_genes = {
    "Monocyte": ["PTPRC", "CD68", "CD163", "AIF1"],
    "Fibroblast": ["COL11A1", "COL3A1"], # collagens so this makes sense  
    "T cell": ["CD8A", "CD4", "CD3D"], # CD8, CD4 also makes sense 
    "NK": ["GNLY", "NKG7", "CD247", "FCER1G", "TYROBP", "KLRG1", "FCGR3A"], # can leave suggested markers
    "Endothelial": ["PLVAP", "PECAM1", "VWF"],
    "Mast": ["KIT", "TPSAB1", "TPSAB2"],
    "Pericyte": ["ACTA2", "NOTCH3", "RGS5"],
    "Tumour": ["ZNF341", "MN1"], # genes identified in the likely driver fusion, not sure what else one would add
}

sc.pl.dotplot(adata, marker_genes, groupby="leiden_res_0.02", standard_scale="var")

# Coarse cell types 
adata.obs["cell_type_lvl1"] = adata.obs["leiden_res_0.02"].map(
    {
        "0": "Lymphocytes",
        "1": "Monocytes",
        "2": "Erythroid",
        "3": "B Cells",
    }
)
sc.pl.dotplot(adata, marker_genes, groupby="leiden_res_0.50", standard_scale="var")

# %%
###################################################################################################################################
# IDENTIFY CLUSTER MARKERS 

# Obtain cluster-specific differentially expressed genes
sc.tl.rank_genes_groups(adata, groupby="leiden_res_0.50", method="wilcoxon")

sc.pl.rank_genes_groups_dotplot(
    adata, groupby="leiden_res_0.50", standard_scale="var", n_genes=5
)

sc.get.rank_genes_groups_df(adata, group="7").head(5)

dc_cluster_genes = sc.get.rank_genes_groups_df(adata, group="7").head(5)["names"]
sc.pl.umap(
    adata,
    color=[*dc_cluster_genes, "leiden_res_0.50"],
    legend_loc="on data",
    frameon=False,
    ncols=3,
)
