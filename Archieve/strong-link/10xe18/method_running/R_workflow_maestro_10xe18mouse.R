qlogin -now no
module unload python/python-3.6.2
module load python/python-3.8.2
conda activate MAESTRO
R

### MAESTRO GEX on e18 data ###
library(MAESTRO)
library(Seurat)
library(Signac)
library(reticulate)
use_python("/home/stat/shuang91/miniconda3/envs/MAESTRO/bin/python3.8",required = TRUE)

setwd("/home/mnt/nzh/nzhanglab/project/SingleCellAlignment/data/10x_RNA_ATAC_EmbryonicMouseBrain")
e18=readRDS("e18.4.20210917.rds")

e18mouseRNA = e18@assays[["RNA"]]@counts

e18mouseATAC = e18@assays[["ATAC"]]@counts

e18$celltype -> e18_celltype


e18mouse.RNA.res <- RNARunSeurat(inputMat = e18mouseRNA,
                             project = "e18mouse",
                             orig.ident = NULL,
                             min.c = 10,
                             min.g = 200,
                             variable.genes = 2000,
                             organism = "GRCm38",
                             dims.use = 1:15,
                             cluster.res = 0.6,
                             only.pos = FALSE,
                             genes.test.use = "presto",
                             genes.cutoff = 1e-05,
                             genes.pct = 0.1,
                             genes.logfc = 0.25,
                             outdir = "/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/MAESTRO/e18mouse/gex/R_ver/")

data(mouse.brain.ALLEN)
e18mouse.RNA.res$RNA <- RNAAnnotateCelltype(RNA = e18mouse.RNA.res$RNA,
                                        gene = e18mouse.RNA.res$genes,
                                        signatures = "mouse.brain.ALLEN",
                                        min.score = 0.05,
                                        outdir = "/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/MAESTRO/e18mouse/gex/R_ver/")

saveRDS(e18mouse.RNA.res, "/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/MAESTRO/e18mouse/gex/R_ver/e18m_res.rds")

### MAESTRO ATAC on e18mouse data ###

###  MAESTRO ATAC on e18mouse data ###
e18mouse.ATAC.res <- ATACRunSeurat(inputMat = e18mouseATAC,
                                   project = "e18mouse_atac",
                                   method = "LSI",
                                   dims.use = 1:30,
                                   cluster.res = 0.6,
                                   only.pos = TRUE,
                                   peaks.test.use = "presto",
                                   peaks.cutoff = 1e-05,
                                   peaks.pct = 0.1,
                                   peaks.logfc = 0.2,
                                   outdir ="/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/MAESTRO/e18mouse/atac/R_ver/")
dim(e18mouse.ATAC.res$ATAC)
### manually calculate gene activity score ####
e18mouse.atac.mat <-as.matrix(e18mouseATAC)
row.names(e18mouse.atac.mat) <-gsub("-","_",row.names(e18mouse.atac.mat),fixed=T )
e18mouse.atac.mat <- Matrix::Matrix(e18mouse.atac.mat, sparse = TRUE)
peaks_list = rownames(e18mouse.atac.mat)

if(class(e18mouse.atac.mat) != "dgCMatrix"){
  e18mouse.atac.mat = as(as.matrix(e18mouse.atac.mat), "dgCMatrix")
}
data(GRCm38.refgenes.genescore.adjusted)
refgenes.genescore = GRCm38.refgenes.genescore.adjusted
source_python(paste(system.file(package="MAESTRO"), "ATACCalculateGenescore.py", sep="/"))

rp_result = calculate_RP_score(cell_peaks = e18mouse.atac.mat, peaks_list = peaks_list, gene_bed_df = refgenes.genescore,
                               genes_list = NULL, decay = 1000, model = "Enhanced")
e18mouse_rpmatrix <- rp_result[[1]]


### use python functions to convert the matrix to R readable dataframes ###
library(reticulate)
use_python("/home/stat/shuang91/miniconda3/envs/MAESTRO/bin/python3.8",required = TRUE)
scipy <- import("scipy", convert = FALSE)
scipy$io$mmwrite("/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/MAESTRO/e18mouse/atac/R_ver/e18mouse_rpmatrix.mtx",e18mouse_rpmatrix)
e18mousematrix <-Matrix::readMM("/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/MAESTRO/e18mouse/atac/R_ver/e18mouse_rpmatrix.mtx")
row.names(e18mousematrix) <- rp_result[[2]]
colnames(e18mousematrix) <- colnames(e18mouse.atac.mat)
e18mousematrix ->maestro_e18mouse_score
saveRDS(maestro_e18mouse_score,"/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/MAESTRO/e18mouse/atac/R_ver/MAESTRO_e18mouse_rpmatrix.RDS")


### continue pipeline ###
e18mouse.gene <- maestro_e18mouse_score
e18mouse.ATAC.res$ATAC <- ATACAttachGenescore(ATAC = e18mouse.ATAC.res$ATAC, RPmatrix = e18mouse.gene)
head(e18mouse.ATAC.res$ATAC@meta.data)

data(mouse.brain.ALLEN)
e18mouse.ATAC.res$ATAC <- ATACAnnotateCelltype(ATAC = e18mouse.ATAC.res$ATAC,
                                               signatures = mouse.brain.ALLEN,
                                               min.score = 0.1,
                                               genes.test.use = "presto",
                                               genes.cutoff = 1E-5,
                                               outdir = "/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/MAESTRO/e18mouse/atac/R_ver/")

e18mouse.ATAC.res$ATAC <- ATACAnnotateChromatinAccessibility(ATAC = e18mouse.ATAC.res$ATAC,
                                                             peaks = e18mouse.ATAC.res$peaks,
                                                             project = "atac_e18mouse",
                                                             giggle.path = "/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/MAESTRO/annotation/giggle.all",
                                                             organism = "GRCm38")

saveRDS(e18mouse.ATAC.res, "/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/MAESTRO/e18mouse/atac/R_ver/e18mouse_atac.rds")

### data integation ###
cd /home/stat/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/MAESTRO/e18mouse

MAESTRO integrate-init --rna-object /home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/MAESTRO/e18mouse/gex/R_ver/e18m_res.rds --atac-object /home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/MAESTRO/e18mouse/atac/R_ver/e18mouse_atac.rds --directory /home/stat/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/MAESTRO/e18mouse/integration --outprefix R_pipeline_MAESTRO_integration
cd integration/
snakemake -np
snakemake -j 16

R
library(MAESTRO)
Seurat1 = readRDS("/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/MAESTRO/e18mouse/atac/R_ver/e18mouse_atac.rds")
Seurat2 = readRDS("/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/MAESTRO/e18mouse/gex/R_ver/e18m_res.rds")

Combine.res = Incorporate(ATAC = Seurat1$ATAC, RNA = Seurat2$RNA,assembly="GRCm38", RPmatrix = NULL, project = "R_ver_integrate", dims.use = 1:30, RNA.res = 0.6, ATAC.res = 0.6)

Seurat.combined = Combine.res$CombinedObj
saveRDS(Seurat.combined, paste0("R_ver_integrate", "_integrate_Object.rds"))
