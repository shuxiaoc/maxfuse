module unload python/python-3.6.2
module load python/python-3.8.2
conda activate scglue2
R

library(anndata)
#### greenleaf data ####
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
library(dplyr)
library(ggplot2)


setwd("/home/mnt/nzh/nzhanglab/project/shuang/scATAC/data/10x_RNA_ATAC_GreenleafCortical")

load("Writeup14n_10x_greenleaf.RData")

greenleaf.rna = greenleaf@assays$RNA@counts
greenleaf.atac = greenleaf@assays$ATAC@counts

greenleaf.obj.rna <- CreateSeuratObject(
  counts = greenleaf.rna,
  assay = "RNA"
)

greenleaf.obj.atac <- CreateSeuratObject(
  counts = greenleaf.atac,
  assay = "RNA"
)

greenleaf$celltype -> greenleaf_celltype

greenleaf.obj.rna$celltype <- greenleaf_celltype
greenleaf.obj.atac$celltype <- greenleaf_celltype

greenleaf.obj.rna@meta.data$domain <- "scRNA-seq"
greenleaf.obj.atac@meta.data$domain <- "scATAC-seq"


setwd("/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/scglue/greenleaf")
library(SeuratDisk)
SaveH5Seurat(greenleaf.obj.rna, filename = "greenleaf_RNA_v2.h5Seurat")
Convert("greenleaf_RNA_v2.h5Seurat", dest = "h5ad")
SaveH5Seurat(greenleaf.obj.atac, filename = "greenleaf_ATAC_v2.h5Seurat")
Convert("greenleaf_ATAC_v2.h5Seurat", dest = "h5ad")
