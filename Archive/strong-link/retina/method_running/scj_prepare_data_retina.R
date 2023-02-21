qlogin -now no
module unload python/python-3.6.2
module load python/python-3.8.2
conda activate /home/mnt/nzh/nzhanglab/project/shuang/miniconda3/envs/scglue2
R


### prepare scJoint input of h5 file ###

library(Seurat)
library(Signac)
library(tables)
library(reticulate)

library(SingleCellExperiment)
library(DropletUtils)
library(scater)
library(ggplot2)
setwd("/home/mnt/nzh/nzhanglab/project/shuang/scATAC/data/retina")

retina<-readRDS("data/retina_seurat.rds")

meta <- read.csv("data/meta.csv")
colnames(meta)[1]<-c("barcode")
retina.rna = retina@assays$RNA@counts
retina.atac = retina@assays$ATAC@counts

meta_subset <- read.csv("data/meta_20k.csv")
colnames(meta_subset)[1]<-c("barcode")
subset_retina.rna <- retina.rna[,meta_subset$barcode]
subset_retina.atac <- retina.atac[,meta_subset$barcode]
meta_subset$annotation ->retina_celltype

row.names(meta_subset) <-meta_subset$barcode
retina.activity=subset_retina.atac
retina.obj.rna <- CreateSeuratObject(
  counts = subset_retina.rna,
  meta.data=meta_subset,
  assay = "RNA"
)

# Only keep common genes between two dataset
common_genes <- intersect(rownames(retina.rna),
                          rownames(retina.atac))
length(common_genes) #21369

retina.obj.activity <- CreateSeuratObject(
  counts = retina.activity,
  assay = "RNA",
  meta.data=meta_subset
)

### create logcounts ###
activity.sce <- as.SingleCellExperiment(retina.obj.activity)
rna.sce <- as.SingleCellExperiment(retina.obj.rna)

# Extract the logcounts data from sce object
exprs_atac <- logcounts(activity.sce[common_genes, ])
exprs_rna <- logcounts(rna.sce[common_genes, ])

source("/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/scJoint/data_to_h5.R")
#source("/Users/sijia_work/Documents/nancy_projects/scATAC_RNA/maestro/scJoint/data_to_h5.R")
write_h5_scJoint(exprs_list = list(rna = exprs_rna,
                                   atac = exprs_atac),
                 h5file_list = c("/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/scJoint/retina/exprs_retina_rna.h5",
                                 "/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/scJoint/retina/exprs_retina_atac.h5"))

write_csv_scJoint(cellType_list =  list(rownames(meta_subset)),
                  csv_list = c("/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/scJoint/retina/cellname_cellType_retina.csv"))
write_csv_scJoint(cellType_list =  list(retina_celltype),
                  csv_list = c("/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/scJoint/retina/cellType_retina.csv"))


write_h5_scJoint(exprs_list = list(rna = exprs_rna,
                                   atac = exprs_atac),
                 h5file_list = c("/Users/sijia_work/Documents/nancy_projects/scATAC_RNA/maestro/scJoint/retina/exprs_10xretina_rna.h5",
                                 "/Users/sijia_work/Documents/nancy_projects/scATAC_RNA/maestro/scJoint/retina/exprs_10xretina_atac.h5"))
write_csv_scJoint(cellType_list =  list(retina_celltype),
                  csv_list = c("/Users/sijia_work/Documents/nancy_projects/scATAC_RNA/maestro/scJoint/retina/cellType_10xretina.csv"))



### final output ###
predict_label <- read.csv("/Users/sijia_work/Documents/nancy_projects/scATAC_RNA/maestro/scJoint/retina/barcode_labels/retina_predictlabel.csv")
rna_pred <- read.csv("/Users/sijia_work/Documents/nancy_projects/scATAC_RNA/maestro/scJoint/retina/barcode_labels/exprs_retina_rna_predictions.txt",header=F,sep=" ")
atac_pred_label <- read.csv("/Users/sijia_work/Documents/nancy_projects/scATAC_RNA/maestro/scJoint/retina/barcode_labels/exprs_retina_atac_knn_predictions.txt",header=F)
atac_pred <- read.csv("/Users/sijia_work/Documents/nancy_projects/scATAC_RNA/maestro/scJoint/retina/barcode_labels/exprs_retina_atac_predictions.txt",header=F,sep=" ")
idx_label <- read.csv("/Users/sijia_work/Documents/nancy_projects/scATAC_RNA/maestro/scJoint/retina/barcode_labels/label_to_idx.txt",header=F,sep=" ")

atac_embed <-read.csv("/Users/sijia_work/Documents/nancy_projects/scATAC_RNA/maestro/scJoint/retina/barcode_labels/exprs_retina_atac_embeddings.txt",header=F,sep=" ")
rna_embed <- read.csv("/Users/sijia_work/Documents/nancy_projects/scATAC_RNA/maestro/scJoint/retina/barcode_labels/exprs_retina_rna_embeddings.txt",header=F,sep=" ")
