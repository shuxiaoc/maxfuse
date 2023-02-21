qlogin -now no
module unload python/python-3.6.2
module load python/python-3.8.2
conda activate MAESTRO
R

### MAESTRO GEX on pbmc data ###

library(MAESTRO)
library(Seurat)
library(Signac)
library(tables)
library(reticulate)
use_python("/home/stat/shuang91/miniconda3/envs/MAESTRO/bin/python3.8",required = TRUE)
setwd("/home/mnt/nzh/nzhanglab/project/SingleCellAlignment/data/10x_RNA_ATAC_PBMC")
load("pbmc_chromvar_annotated.rda")

pbmc.rna = pbmc@assays[["RNA"]]@counts

pbmc.atac = pbmc@assays[["ATAC"]]@counts
pbmc$citeseq.celltype -> pbmc_celltype

pbmc.RNA.res <- RNARunSeurat(inputMat = pbmc.rna,
                                 project = "pbmc",
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
                                 outdir = "/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/MAESTRO/pbmc/gex/R_ver/")

data(human.immune.CIBERSORT)
pbmc.RNA.res$RNA <- RNAAnnotateCelltype(RNA = pbmc.RNA.res$RNA,
                                            gene = pbmc.RNA.res$genes,
                                            signatures = "human.immune.CIBERSORT",
                                            min.score = 0.05,
                                            outdir = "/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/MAESTRO/pbmc/gex/R_ver/")

saveRDS(pbmc.RNA.res, "/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/MAESTRO/pbmc/gex/R_ver/pbmc_res.rds")


###  MAESTRO ATAC on pbmc data ###
pbmc.ATAC.res <- ATACRunSeurat(inputMat = pbmc.atac,
                               project = "pbmc_atac",
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
                               outdir ="/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/MAESTRO/pbmc/atac/R_ver/")

### manually calculate gene activity score ####
pbmc.atac.mat <-as.matrix(pbmc.atac)
row.names(pbmc.atac.mat) <-gsub("-","_",row.names(pbmc.atac.mat),fixed=T )
pbmc.atac.mat <- Matrix::Matrix(pbmc.atac.mat, sparse = TRUE)
row.names(pbmc.atac.mat) <-gsub("-","_",row.names(pbmc.atac.mat),fixed=T)
peaks_list = rownames(pbmc.atac.mat)

if(class(pbmc.atac.mat) != "dgCMatrix"){
  pbmc.atac.mat = as(as.matrix(pbmc.atac.mat), "dgCMatrix")
}
data(GRCh38.refgenes.genescore.adjusted)
refgenes.genescore = GRCh38.refgenes.genescore.adjusted
source_python(paste(system.file(package="MAESTRO"), "ATACCalculateGenescore.py", sep="/"))

rp_result = calculate_RP_score(cell_peaks = pbmc.atac.mat, peaks_list = peaks_list, gene_bed_df = refgenes.genescore,
                               genes_list = NULL, decay = 1000, model = "Enhanced")
pbmc_rpmatrix <- rp_result[[1]]


### use python functions to convert the matrix to R readable dataframes ###
library(reticulate)
use_python("/home/stat/shuang91/miniconda3/envs/MAESTRO/bin/python3.8",required = TRUE)
scipy <- import("scipy", convert = FALSE)
scipy$io$mmwrite("/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/MAESTRO/pbmc/atac/R_ver/pbmc_rpmatrix.mtx",pbmc_rpmatrix)
pbmcmatrix <-Matrix::readMM("/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/MAESTRO/pbmc/atac/R_ver/pbmc_rpmatrix.mtx")
row.names(pbmcmatrix) <- rp_result[[2]]
colnames(pbmcmatrix) <- colnames(pbmc.atac.mat)
pbmcmatrix ->maestro_pbmc_score
saveRDS(maestro_pbmc_score,"/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/MAESTRO/pbmc/atac/R_ver/MAESTRO_pbmc_rpmatrix.RDS")


### continue pipeline ###
pbmc.gene <- maestro_pbmc_score
pbmc.ATAC.res$ATAC <- ATACAttachGenescore(ATAC = pbmc.ATAC.res$ATAC, RPmatrix = pbmc.gene)
head(pbmc.ATAC.res$ATAC@meta.data)

data(human.immune.CIBERSORT)
pbmc.ATAC.res$ATAC <- ATACAnnotateCelltype(ATAC = pbmc.ATAC.res$ATAC,
                                           signatures = human.immune.CIBERSORT,
                                           min.score = 0.1,
                                           genes.test.use = "presto",
                                           genes.cutoff = 1E-5,
                                           outdir = "/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/MAESTRO/pbmc/atac/R_ver/")

pbmc.ATAC.res$ATAC <- ATACAnnotateChromatinAccessibility(ATAC = pbmc.ATAC.res$ATAC,
                                                         peaks = pbmc.ATAC.res$peaks,
                                                         project = "atac_pbmc",
                                                         giggle.path = "/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/MAESTRO/annotation/giggle.all",
                                                         organism = "GRCh38")

saveRDS(pbmc.ATAC.res, "/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/MAESTRO/pbmc/atac/R_ver/pbmc_atac.rds")


### data integration ###
cd /home/stat/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/MAESTRO/pbmc

MAESTRO integrate-init --rna-object /home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/MAESTRO/pbmc/gex/R_ver/pbmc_res.rds --atac-object /home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/MAESTRO/pbmc/atac/R_ver/pbmc_atac.rds --directory /home/stat/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/MAESTRO/pbmc/integration --outprefix R_pipeline_MAESTRO_integration
cd integration/
snakemake -np
snakemake -j 16
