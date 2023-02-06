qlogin -now no
module unload python/python-3.6.2
module load python/python-3.8.2
conda activate scjoint
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
setwd("/home/mnt/nzh/nzhanglab/project/shuang/scATAC/data/10x_RNA_ATAC_GreenleafCortical")

load("Writeup14n_10x_greenleaf.RData")

greenleafRNA = greenleaf@assays$RNA@counts
all_data <-greenleaf
Seurat::DefaultAssay(all_data) <- "ATAC"
gene_activities <- Signac::GeneActivity(all_data)
greenleafATAC = gene_activities


library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
library(dplyr)
library(ggplot2)

#greenleaf<- SCTransform(greenleaf, assay="geneActivity", verbose = FALSE, new.assay.name = 'SCT.ACTIVITY') %>% RunPCA(verbose=F,  reduction.name = 'pca.activity') %>% RunUMAP(verbose=F, dims = 1:20, reduction='pca.activity', reduction.name='umap.activity')

greenleaf$celltype -> greenleaf_celltype

greenleaf.obj.rna <- CreateSeuratObject(
  counts = greenleafRNA,
  assay = "RNA"
)

# Only keep common genes between two dataset
rownames(greenleafATAC) <-gsub("ATAC-","",rownames(greenleafATAC),fixed=T)

common_genes <- intersect(rownames(greenleafRNA),
                          rownames(greenleafATAC))
length(common_genes)

greenleaf.obj.activity <- CreateSeuratObject(
  counts = greenleafATAC,
  assay = "RNA"
)

### create logcounts ###
activity.sce <- as.SingleCellExperiment(greenleaf.obj.activity)
rna.sce <- as.SingleCellExperiment(greenleaf.obj.rna)

# Extract the logcounts data from sce object
exprs_atac <- logcounts(activity.sce[common_genes, ])
exprs_rna <- logcounts(rna.sce[common_genes, ])

source("/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/scJoint/data_to_h5.R")
#source("/Users/sijia_work/Documents/nancy_projects/scATAC_RNA/maestro/scJoint/data_to_h5.R")
write_h5_scJoint(exprs_list = list(rna = exprs_rna,
                                   atac = exprs_atac),
                 h5file_list = c("/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/scJoint/greenleaf/exprs_10xgreenleaf_rna.h5",
                                 "/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/scJoint/greenleaf/exprs_10xgreenleaf_atac.h5"))

write_csv_scJoint(cellType_list =  list(names(greenleaf_celltype)),
                  csv_list = c("/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/scJoint/greenleaf/cellname_cellType_10xgreenleaf.csv"))
write_csv_scJoint(cellType_list =  list(greenleaf_celltype),
                  csv_list = c("/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/scJoint/greenleaf/cellType_10xgreenleaf.csv"))


write_h5_scJoint(exprs_list = list(rna = exprs_rna,
                                   atac = exprs_atac),
                 h5file_list = c("/Users/sijia_work/Documents/nancy_projects/scATAC_RNA/maestro/scJoint/greenleaf/exprs_10xgreenleaf_rna.h5",
                                 "/Users/sijia_work/Documents/nancy_projects/scATAC_RNA/maestro/scJoint/greenleaf/exprs_10xgreenleaf_atac.h5"))
write_csv_scJoint(cellType_list =  list(greenleaf_celltype),
                  csv_list = c("/Users/sijia_work/Documents/nancy_projects/scATAC_RNA/maestro/scJoint/greenleaf/cellType_10xgreenleaf.csv"))

save(gene_activities,file="/home/mnt/nzh/nzhanglab/project/shuang/scATAC/data/10x_RNA_ATAC_GreenleafCortical/seurat_geneactivity.RData")

### final output ###
predict_label <- read.csv("/Users/sijia_work/Documents/nancy_projects/scATAC_RNA/maestro/scJoint/greenleaf/barcode_labels/greenleaf_predictlabel.csv")
rna_pred <- read.csv("/Users/sijia_work/Documents/nancy_projects/scATAC_RNA/maestro/scJoint/greenleaf/barcode_labels/exprs_10xgreenleaf_rna_predictions.txt",header=F,sep=" ")
atac_pred_label <- read.csv("/Users/sijia_work/Documents/nancy_projects/scATAC_RNA/maestro/scJoint/greenleaf/barcode_labels/exprs_10xgreenleaf_atac_knn_predictions.txt",header=F)
atac_pred <- read.csv("/Users/sijia_work/Documents/nancy_projects/scATAC_RNA/maestro/scJoint/greenleaf/barcode_labels/exprs_10xgreenleaf_atac_predictions.txt",header=F,sep=" ")
idx_label <- read.csv("/Users/sijia_work/Documents/nancy_projects/scATAC_RNA/maestro/scJoint/greenleaf/barcode_labels/label_to_idx.txt",header=F,sep=" ")

atac_embed <-read.csv("/Users/sijia_work/Documents/nancy_projects/scATAC_RNA/maestro/scJoint/greenleaf/barcode_labels/exprs_10xgreenleaf_atac_embeddings.txt",header=F,sep=" ")
rna_embed <- read.csv("/Users/sijia_work/Documents/nancy_projects/scATAC_RNA/maestro/scJoint/greenleaf/barcode_labels/exprs_10xgreenleaf_rna_embeddings.txt",header=F,sep=" ")
