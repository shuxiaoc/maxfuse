module unload python/python-3.6.2
module load python/python-3.8.2
conda activate scglue2
R

library(anndata)
#### pbmc data ####
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
library(dplyr)
library(ggplot2)


setwd("/home/mnt/nzh/nzhanglab/project/SingleCellAlignment/data/10x_RNA_ATAC_PBMC")
load("pbmc_chromvar_annotated.rda")
pbmc.rna = pbmc@assays[["RNA"]]@counts
pbmc.obj.rna <- CreateSeuratObject(
  counts = pbmc.rna,
  assay = "RNA"
)

pbmc.atac = pbmc@assays[["ATAC"]]@counts
pbmc.obj.atac <- CreateSeuratObject(
  counts = pbmc.atac,
  assay = "RNA"
)

pbmc$citeseq.celltype -> pbmc_celltype

pbmc.obj.rna$celltype <- pbmc_celltype
pbmc.obj.atac$celltype <- pbmc_celltype

pbmc.obj.rna@meta.data$domain <- "scRNA-seq"
pbmc.obj.atac@meta.data$domain <- "scATAC-seq"

setwd("/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/scglue/pbmc")
library(SeuratDisk)
SaveH5Seurat(pbmc.obj.rna, filename = "pbmc_RNA_v2.h5Seurat")
Convert("pbmc_RNA_v2.h5Seurat", dest = "h5ad")
SaveH5Seurat(pbmc.obj.atac, filename = "pbmc_ATAC_v2.h5Seurat")
Convert("pbmc_ATAC_v2.h5Seurat", dest = "h5ad")
