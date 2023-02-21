### prepare data for scglue ###
#setwd("/Users/sijia_work/Documents/nancy_projects/scATAC_RNA/glue/")

module unload python/python-3.6.2
module load python/python-3.8.2
conda activate scglue2
R
library(BiocGenerics, lib.loc="/home/mnt/nzh/nzhanglab/project/shuang/miniconda3/envs/scglue2/envs/MAESTRO/lib/R/library")
library(S4Vectors, lib.loc="/home/mnt/nzh/nzhanglab/project/shuang/miniconda3/envs/scglue2/envs/MAESTRO/lib/R/library")
library(IRanges, lib.loc="/home/mnt/nzh/nzhanglab/project/shuang/miniconda3/envs/scglue2/envs/MAESTRO/lib/R/library")
library(GenomeInfoDb, lib.loc="/home/mnt/nzh/nzhanglab/project/shuang/miniconda3/envs/scglue2/envs/MAESTRO/lib/R/library")

library(anndata)

library(Seurat)
library(Signac,lib.loc="/home/mnt/nzh/nzhanglab/project/shuang/miniconda3/envs/scglue2/envs/MAESTRO/lib/R/library")
library(tables)
library(reticulate)

library(SingleCellExperiment)
library(DropletUtils)
library(scater)
library(ggplot2)

#### e18 mouse data ####
setwd("/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/data/10x_RNA_ATAC_EmbryonicMouseBrain/")
#setwd("/Users/sijia_work/Dropbox/SingleCellAlignment/data/10x_RNA_ATAC_EmbryonicMouseBrain/")

e18=readRDS("e18.4.20210917.rds")
table(e18$celltype)

e18mouseRNA = e18@assays[["RNA"]]@counts
e18mouseATAC = e18@assays[["ATAC"]]@counts

e18.obj.rna <- CreateSeuratObject(
  counts = e18mouseRNA,
  assay = "RNA"
)
e18.obj.atac <- CreateSeuratObject(
  counts = e18mouseATAC,
  assay = "RNA"
)

e18.obj.rna$celltype <- e18$celltype
e18.obj.atac$celltype <- e18$celltype

e18.obj.rna@meta.data$domain <- "scRNA-seq"
e18.obj.atac@meta.data$domain <- "scATAC-seq"


setwd("/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/scglue/e18mouse")
library(SeuratDisk)
SaveH5Seurat(e18.obj.rna, filename = "e18mouse_RNA_v2.h5Seurat")
Convert("e18mouse_RNA_v2.h5Seurat", dest = "h5ad")
SaveH5Seurat(e18.obj.atac, filename = "e18mouse_ATAC_v2.h5Seurat")
Convert("e18mouse_ATAC_v2.h5Seurat", dest = "h5ad")
