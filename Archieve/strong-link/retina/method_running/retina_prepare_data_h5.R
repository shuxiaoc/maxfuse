module unload python/python-3.6.2
module load python/python-3.8.2
conda activate /home/mnt/nzh/nzhanglab/project/shuang/miniconda3/envs/scglue2
R

#### retina data ####
library(Seurat)
library(Signac)
library(GenomeInfoDb)
library(dplyr)
library(ggplot2)


setwd("/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/data/Retina/")
retina<-readRDS("data/retina_peak.rds")

meta <- read.csv("data/meta20k.csv")
colnames(meta)[1]<-c("barcode")
retina.rna = retina@assays$RNA@counts

meta1 <- meta[match(colnames(retina.rna), meta$barcode),]

meta1$annotation ->retina_celltype

meta_subset <- read.csv("data/meta_20k.csv")
colnames(meta_subset)[1]<-c("barcode")
subset_retina.rna <- retina.rna[,meta_subset$barcode]
meta_subset$annotation ->retina_celltype


retina.obj.rna <- CreateSeuratObject(
  counts = subset_retina.rna,
  assay = "RNA"
)
retina.obj.rna$celltype <- retina_celltype


retina.atac = retina@assays$peak@counts
subset_retina.atac <- retina.atac[,meta_subset$barcode]


retina.obj.atac <- CreateSeuratObject(
  counts = subset_retina.atac,
  assay = "RNA"
)

retina.obj.atac$celltype <- retina_celltype

retina.obj.rna@meta.data$domain <- "scRNA-seq"
retina.obj.atac@meta.data$domain <- "scATAC-seq"


setwd("/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/scglue/retina")
library(SeuratDisk)
SaveH5Seurat(retina.obj.rna, filename = "retina_RNA.h5Seurat")
Convert("retina_RNA.h5Seurat", dest = "h5ad")
SaveH5Seurat(retina.obj.atac, filename = "retina_ATAC.h5Seurat")
Convert("retina_ATAC.h5Seurat", dest = "h5ad")
