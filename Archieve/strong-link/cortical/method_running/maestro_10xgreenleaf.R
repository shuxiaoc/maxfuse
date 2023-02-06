#### 10x greenleaf maestro gene activity score calculation ####
qlogin -now no
module unload python/python-3.6.2
module load python/python-3.8.2
conda activate MAESTRO
R

### MAESTRO GEX on greenleaf data ###

library(MAESTRO)
library(Seurat)
library(Signac)
library(tables)
library(reticulate)
library(SummarizedExperiment)
use_python("/home/stat/shuang91/miniconda3/envs/MAESTRO/bin/python3.8",required = TRUE)
setwd("/home/mnt/nzh/nzhanglab/project/shuang/scATAC/data/10x_RNA_ATAC_GreenleafCortical/")
load("Writeup14n_10x_greenleaf.RData")
greenleafRNA = greenleaf@assays$RNA@counts
#greenleafATAC = greenleaf@assays$ATAC@counts
# greenleaf$celltype

all_data <-greenleaf
Seurat::DefaultAssay(all_data) <- "ATAC"
gene_activities <- Signac::GeneActivity(all_data)

greenleaf.rna <- CreateSeuratObject(counts = greenleafRNA, project = "greenleaf.rna", min.cells = 10, min.features = 200, assay = 'RNA')

greenleaf_rna = greenleaf.rna@assays[["RNA"]]@counts

greenleaf.RNA.res <- RNARunSeurat(inputMat = greenleaf_rna,
                               project = "greenleaf",
                               orig.ident = NULL,
                               min.c = 10,
                               min.g = 200,
                               variable.genes = 2000,
                               organism = "GRCh38",
                               dims.use = 1:15,
                               cluster.res = 0.6,
                               only.pos = FALSE,
                               genes.test.use = "presto",
                               genes.cutoff = 1e-05,
                               genes.pct = 0.1,
                               genes.logfc = 0.25,
                               outdir = "/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/MAESTRO/greenleaf/gex/R_ver/")

data(human.immune.CIBERSORT)
greenleaf.RNA.res$RNA <- RNAAnnotateCelltype(RNA = greenleaf.RNA.res$RNA,
                                          gene = greenleaf.RNA.res$genes,
                                          signatures = "human.immune.CIBERSORT",
                                          min.score = 0.05,
                                          outdir = "/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/MAESTRO/greenleaf/gex/R_ver/")

saveRDS(greenleaf.RNA.res, "/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/MAESTRO/greenleaf/gex/R_ver/greenleaf_res.rds")


###  MAESTRO ATAC on greenleaf data ###

greenleaf_atac = gene_activities

greenleaf.ATAC.res <- ATACRunSeurat(inputMat = greenleaf_atac,
                                 project = "greenleaf_atac",
                                 min.c = 50,
                                 min.p = 500,
                                 method = "LSI",
                                 dims.use = 1:30,
                                 cluster.res = 0.6,
                                 only.pos = TRUE,
                                 peaks.test.use = "presto",
                                 peaks.cutoff = 1e-05,
                                 peaks.pct = 0.1,
                                 peaks.logfc = 0.2,
                                 outdir ="/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/MAESTRO/greenleaf/atac/R_ver/")

### manually calculate gene activity score ####

setwd("/home/mnt/nzh/nzhanglab/project/shuang/scATAC/data/10x_RNA_ATAC_GreenleafCortical/")

greenleafRNA = greenleaf@assays$RNA@counts
greenleafATAC = greenleaf@assays$ATAC@counts


greenleaf.atac.mat <-greenleafATAC
row.names(greenleaf.atac.mat) <-gsub("-","_",row.names(greenleaf.atac.mat),fixed=T )
greenleaf.atac.mat <- Matrix::Matrix(greenleaf.atac.mat, sparse = TRUE)
row.names(greenleaf.atac.mat) <-gsub("-","_",row.names(greenleaf.atac.mat),fixed=T)
peaks_list = rownames(greenleaf.atac.mat)

if(class(greenleaf.atac.mat) != "dgCMatrix"){
  greenleaf.atac.mat = as(as.matrix(greenleaf.atac.mat), "dgCMatrix")
}
data(GRCh38.refgenes.genescore.adjusted)
refgenes.genescore = GRCh38.refgenes.genescore.adjusted
source_python(paste(system.file(package="MAESTRO"), "ATACCalculateGenescore.py", sep="/"))

rp_result = calculate_RP_score(cell_peaks = greenleaf.atac.mat, peaks_list = peaks_list, gene_bed_df = refgenes.genescore,
                               genes_list = NULL, decay = 1000, model = "Enhanced")
greenleaf_rpmatrix <- rp_result[[1]]


### use python functions to convert the matrix to R readable dataframes ###
library(reticulate)
use_python("/home/stat/shuang91/miniconda3/envs/MAESTRO/bin/python3.8",required = TRUE)
scipy <- import("scipy", convert = FALSE)
scipy$io$mmwrite("/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/MAESTRO/greenleaf/atac/R_ver/greenleaf_rpmatrix.mtx",greenleaf_rpmatrix)
greenleafmatrix <-Matrix::readMM("/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/MAESTRO/greenleaf/atac/R_ver/greenleaf_rpmatrix.mtx")
row.names(greenleafmatrix) <- rp_result[[2]]
colnames(greenleafmatrix) <- colnames(greenleaf.atac.mat)
greenleafmatrix ->maestro_greenleaf_score
saveRDS(maestro_greenleaf_score,"/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/MAESTRO/greenleaf/atac/R_ver/MAESTRO_greenleaf_rpmatrix.RDS")


### continue pipeline ###
greenleaf.gene <- maestro_greenleaf_score
greenleaf.ATAC.res$ATAC <- ATACAttachGenescore(ATAC = greenleaf.ATAC.res$ATAC, RPmatrix = greenleaf.gene)
head(greenleaf.ATAC.res$ATAC@meta.data)

data(human.immune.CIBERSORT)
greenleaf.ATAC.res$ATAC <- ATACAnnotateCelltype(ATAC = greenleaf.ATAC.res$ATAC,
                                             signatures = human.immune.CIBERSORT,
                                             min.score = 0.1,
                                             genes.test.use = "presto",
                                             genes.cutoff = 1E-5,
                                             outdir = "/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/MAESTRO/greenleaf/atac/R_ver/")

greenleaf.ATAC.res$ATAC <- ATACAnnotateChromatinAccessibility(ATAC = greenleaf.ATAC.res$ATAC,
                                                           peaks = greenleaf.ATAC.res$peaks,
                                                           project = "atac_greenleaf",
                                                           giggle.path = "/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/MAESTRO/annotation/giggle.all",
                                                           organism = "GRCh38")

saveRDS(greenleaf.ATAC.res, "/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/MAESTRO/greenleaf/atac/R_ver/greenleaf_atac.rds")

### MAESTRO integration ###
MAESTRO integrate-init --rna-object /home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/MAESTRO/greenleaf/gex/R_ver/greenleaf_res.rds --atac-object /home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/MAESTRO/greenleaf/atac/R_ver/greenleaf_atac.rds --directory /home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/MAESTRO/greenleaf/ --outprefix integrate
snakemake -np
snakemake -j 16
