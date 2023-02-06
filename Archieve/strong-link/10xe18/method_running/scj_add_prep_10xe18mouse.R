qlogin -now no
module unload python/python-3.6.2
module load python/python-3.8.2
conda activate MAESTRO
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
setwd("/home/mnt/nzh/nzhanglab/project/shuang/scATAC/data/10x_RNA_ATAC_EmbryonicMouseBrain/")

e18=readRDS("e18.4.20210917.rds")

table(e18$celltype)

e18mouseRNA = e18@assays[["RNA"]]@counts

e18mouseATAC = e18@assays[["ATAC"]]@counts

DefaultAssay(e18)='ATAC'
Annotation(e18)
frags=UpdatePath(Fragments(e18)[[1]], new.path = 'e18_mouse_brain_fresh_5k_atac_fragments.tsv.gz')
Fragments(e18)=NULL
e18=SetAssayData(e18, slot = "fragments", new.data = frags)

e18gene.activities <- GeneActivity(e18)

library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
library(dplyr)
library(ggplot2)


e18[["ACTIVITY"]] <- CreateAssayObject(counts = e18gene.activities)
DefaultAssay(e18) <- "ACTIVITY"
# SCTransform normalization and PCA dimensional reduction on gene activity
e18<- SCTransform(e18, assay="ACTIVITY", verbose = FALSE, new.assay.name = 'SCT.ACTIVITY') %>% RunPCA(verbose=F,  reduction.name = 'pca.activity') %>% RunUMAP(verbose=F, dims = 1:20, reduction='pca.activity', reduction.name='umap.activity')


e18$celltype -> e18_celltype

e18.obj.rna <- CreateSeuratObject(
  counts = e18mouseRNA,
  assay = "RNA"
)

# Only keep common genes between two dataset
common_genes <- intersect(rownames(e18.obj.rna),
                          rownames(e18gene.activities))
length(common_genes)

e18.obj.activity <- CreateSeuratObject(
  counts = e18gene.activities,
  assay = "RNA"
)

### create logcounts ###
activity.sce <- as.SingleCellExperiment(e18.obj.activity)
rna.sce <- as.SingleCellExperiment(e18.obj.rna)

# Extract the logcounts data from sce object
exprs_atac <- logcounts(activity.sce[common_genes, ])
exprs_rna <- logcounts(rna.sce[common_genes, ])

source("/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/scJoint/data_to_h5.R")
#source("/Users/sijia_work/Documents/nancy_projects/scATAC_RNA/maestro/scJoint/data_to_h5.R")
write_h5_scJoint(exprs_list = list(rna = exprs_rna,
                                   atac = exprs_atac),
                 h5file_list = c("/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/scJoint/e18mouse/exprs_10xe18_rna.h5",
                                 "/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/scJoint/e18mouse/exprs_10xe18_atac.h5"))

write_csv_scJoint(cellType_list =  list(names(e18_celltype)),
                  csv_list = c("/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/scJoint/e18mouse/cellname_cellType_10xe18.csv"))
write_csv_scJoint(cellType_list =  list(e18_celltype),
                  csv_list = c("/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/scJoint/e18mouse/cellType_10xe18.csv"))



### final output ###
e18_predict_label <- read.csv("/Users/sijia_work/Documents/nancy_projects/scATAC_RNA/maestro/scJoint/e18mouse/barcode_labels/e18_predictlabel_celltype.csv")
e18_rna_pred <- read.csv("/Users/sijia_work/Documents/nancy_projects/scATAC_RNA/maestro/scJoint/e18mouse/barcode_labels/exprs_10xe18_rna_predictions.txt",header=F,sep=" ")
e18_atac_pred_label <- read.csv("/Users/sijia_work/Documents/nancy_projects/scATAC_RNA/maestro/scJoint/e18mouse/barcode_labels/exprs_10xe18_atac_knn_predictions.txt",header=F)
e18_atac_pred <- read.csv("/Users/sijia_work/Documents/nancy_projects/scATAC_RNA/maestro/scJoint/e18mouse/barcode_labels/exprs_10xe18_atac_predictions.txt",header=F,sep=" ")
e18_idx_label <- read.csv("/Users/sijia_work/Documents/nancy_projects/scATAC_RNA/maestro/scJoint/e18mouse/barcode_labels/label_to_idx.txt",header=F,sep=" ")

e18_atac_embed <-read.csv("/Users/sijia_work/Documents/nancy_projects/scATAC_RNA/maestro/scJoint/e18mouse/barcode_labels/exprs_10xe18_atac_embeddings.txt",header=F,sep=" ")
e18_rna_embed <- read.csv("/Users/sijia_work/Documents/nancy_projects/scATAC_RNA/maestro/scJoint/e18mouse/barcode_labels/exprs_10xe18_rna_embeddings.txt",header=F,sep=" ")
