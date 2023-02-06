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
setwd("/home/mnt/nzh/nzhanglab/project/SingleCellAlignment/data/10x_RNA_ATAC_PBMC")
#setwd("/Users/sijia_work/Dropbox/SingleCellAlignment/data/10x_RNA_ATAC_PBMC/")

load("pbmc_chromvar_annotated.rda")

pbmc.rna = pbmc@assays[["RNA"]]@counts

library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
library(dplyr)
library(ggplot2)

DefaultAssay(pbmc)='ATAC'
Annotation(pbmc)
frags=UpdatePath(Fragments(pbmc)[[1]], new.path = 'pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz')
Fragments(pbmc)=NULL
pbmc=SetAssayData(pbmc, slot = "fragments", new.data = frags)

gene.activities <- GeneActivity(pbmc)
pbmc[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(pbmc) <- "ACTIVITY"
# SCTransform normalization and PCA dimensional reduction on gene activity
pbmc<- SCTransform(pbmc, assay="ACTIVITY", verbose = FALSE, new.assay.name = 'SCT.ACTIVITY') %>% RunPCA(verbose=F,  reduction.name = 'pca.activity') %>% RunUMAP(verbose=F, dims = 1:20, reduction='pca.activity', reduction.name='umap.activity')
save(pbmc, file='0418_pbmc_activity.rda')


pbmc$citeseq.celltype -> pbmc_celltype

pbmc.obj.rna <- CreateSeuratObject(
  counts = pbmc.rna,
  assay = "RNA"
)

# Only keep common genes between two dataset
common_genes <- intersect(rownames(pbmc.rna),
                          rownames(gene.activities))
length(common_genes)

pbmc.obj.activity <- CreateSeuratObject(
  counts = gene.activities,
  assay = "RNA"
)

### create logcounts ###
activity.sce <- as.SingleCellExperiment(pbmc.obj.activity)
rna.sce <- as.SingleCellExperiment(pbmc.obj.rna)

# Extract the logcounts data from sce object
exprs_atac <- logcounts(activity.sce[common_genes, ])
exprs_rna <- logcounts(rna.sce[common_genes, ])

source("/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/scJoint/data_to_h5.R")
#source("/Users/sijia_work/Documents/nancy_projects/scATAC_RNA/maestro/scJoint/data_to_h5.R")
write_h5_scJoint(exprs_list = list(rna = exprs_rna,
                                   atac = exprs_atac),
                 h5file_list = c("/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/scJoint/pbmc/exprs_10xPBMC_rna.h5",
                                 "/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/scJoint/pbmc/exprs_10xPBMC_atac.h5"))

write_csv_scJoint(cellType_list =  list(names(pbmc_celltype)),
                  csv_list = c("/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/scJoint/pbmc/cellname_cellType_10xPBMC.csv"))


write_h5_scJoint(exprs_list = list(rna = exprs_rna,
                                   atac = exprs_atac),
                 h5file_list = c("/Users/sijia_work/Documents/nancy_projects/scATAC_RNA/maestro/scJoint/pbmc/exprs_10xPBMC_rna.h5",
                                 "/Users/sijia_work/Documents/nancy_projects/scATAC_RNA/maestro/scJoint/pbmc/exprs_10xPBMC_atac.h5"))
write_csv_scJoint(cellType_list =  list(pbmc_celltype),
                  csv_list = c("/Users/sijia_work/Documents/nancy_projects/scATAC_RNA/maestro/scJoint/pbmc/cellType_10xPBMC.csv"))



### final output ###
predict_label <- read.csv("/Users/sijia_work/Documents/nancy_projects/scATAC_RNA/maestro/scJoint/pbmc/barcode_labels/pbmc_predictlabel.csv")
rna_pred <- read.csv("/Users/sijia_work/Documents/nancy_projects/scATAC_RNA/maestro/scJoint/pbmc/barcode_labels/exprs_10xPBMC_rna_predictions.txt",header=F,sep=" ")
atac_pred_label <- read.csv("/Users/sijia_work/Documents/nancy_projects/scATAC_RNA/maestro/scJoint/pbmc/barcode_labels/exprs_10xPBMC_atac_knn_predictions.txt",header=F)
atac_pred <- read.csv("/Users/sijia_work/Documents/nancy_projects/scATAC_RNA/maestro/scJoint/pbmc/barcode_labels/exprs_10xPBMC_atac_predictions.txt",header=F,sep=" ")
idx_label <- read.csv("/Users/sijia_work/Documents/nancy_projects/scATAC_RNA/maestro/scJoint/pbmc/barcode_labels/label_to_idx.txt",header=F,sep=" ")

atac_embed <-read.csv("/Users/sijia_work/Documents/nancy_projects/scATAC_RNA/maestro/scJoint/pbmc/barcode_labels/exprs_10xPBMC_atac_embeddings.txt",header=F,sep=" ")
rna_embed <- read.csv("/Users/sijia_work/Documents/nancy_projects/scATAC_RNA/maestro/scJoint/pbmc/barcode_labels/exprs_10xPBMC_rna_embeddings.txt",header=F,sep=" ")
