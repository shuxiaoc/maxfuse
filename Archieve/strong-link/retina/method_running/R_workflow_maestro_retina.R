#### 10x retina maestro gene activity score calculation ####
qlogin -now no
module unload python/python-3.6.2
module load python/python-3.8.2
mamba activate MAESTRO
R

### MAESTRO GEX on retina data ###

library(MAESTRO)
library(Seurat)
library(Signac)
library(tables)
library(reticulate)
library(SummarizedExperiment)
use_python("/home/mnt/nzh/nzhanglab/project/shuang/miniconda3/envs/scglue2/envs/MAESTRO/bin/python3.8",required = TRUE)
setwd("/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/data/Retina/data")
retina<-readRDS("retina_peak.rds")

retinaRNA <- retina@assays$RNA@counts
retinaATAC <- retina@assays$ATAC@counts


meta_subset <- read.csv("meta_20k.csv")
colnames(meta_subset)[1]<-c("barcode")
subset_retina.rna <- retinaRNA[,meta_subset$barcode]
subset_retina.atac <- retinaATAC[,meta_subset$barcode]
meta_subset$annotation ->retina_celltype

retinaRNA = retina@assays$RNA@counts
#retinaATAC = retina@assays$ATAC@counts
# retina$celltype
meta <- read.csv("data/meta.csv")
colnames(meta)[1]<-c("barcode")
meta1 <- meta[match(colnames(retina.rna), meta$barcode),]
meta1$annotation ->retina_celltype
row.names(meta1)<-meta1[,1]

retina.RNA.res <- RNARunSeurat(inputMat = subset_retina.rna,
                                  project = "retina",
                                  orig.ident = NULL,
                                  mito = FALSE,
                                  min.c = 0,
                                  min.g = 0,
                                  variable.genes = 2000,
                                  organism = "GRCh38",
                                  dims.use = 1:15,
                                  cluster.res = 0.6,
                                  only.pos = FALSE,
                                  genes.test.use = "presto",
                                  genes.cutoff = 1e-05,
                                  genes.pct = 0.1,
                                  genes.logfc = 0.25,
                                  outdir = "/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/MAESTRO/retina/gex/R_ver/")


data(human.immune.CIBERSORT)
retina.RNA.res$RNA <- RNAAnnotateCelltype(RNA = retina.RNA.res$RNA,
                                             gene = retina.RNA.res$genes,
                                             signatures = "human.immune.CIBERSORT",
                                             min.score = 0.05,
                                             outdir = "/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/MAESTRO/retina/gex/R_ver/")

saveRDS(retina.RNA.res, "/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/MAESTRO/retina/gex/R_ver/retina_res.rds")


###  MAESTRO ATAC on retina data ###

retina.atac <- CreateSeuratObject(counts = retina_atac, project = "retina.atac",  assay = 'RNA',meta.data=meta1)


retina.ATAC.res <- ATACRunSeurat(inputMat = subset_retina.atac,
                                    project = "retina_atac",
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
                                    outdir ="/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/MAESTRO/retina/atac/R_ver/")

### manually calculate gene activity score ####

setwd("/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/data/Retina/data")
retina<-readRDS("retina_peak.rds")

retinaRNA = retina@assays$RNA@counts
retinaATAC = retina@assays$peak@counts
meta_subset <- read.csv("meta_20k.csv")
colnames(meta_subset)[1]<-c("barcode")
subset_retina.rna <- retinaRNA[,meta_subset$barcode]
subset_retina.atac <- retinaATAC[,meta_subset$barcode]


retina.atac.mat <-subset_retina.atac
row.names(retina.atac.mat) <-gsub("-","_",row.names(retina.atac.mat),fixed=T )
retina.atac.mat <- Matrix::Matrix(retina.atac.mat, sparse = TRUE)
row.names(retina.atac.mat) <-gsub("-","_",row.names(retina.atac.mat),fixed=T)
peaks_list = rownames(retina.atac.mat)

if(class(retina.atac.mat) != "dgCMatrix"){
  retina.atac.mat = as(as.matrix(retina.atac.mat), "dgCMatrix")
}
data(GRCh38.refgenes.genescore.adjusted)
refgenes.genescore = GRCh38.refgenes.genescore.adjusted
source_python(paste(system.file(package="MAESTRO"), "ATACCalculateGenescore.py", sep="/"))

rp_result = calculate_RP_score(cell_peaks = retina.atac.mat, peaks_list = peaks_list, gene_bed_df = refgenes.genescore,
                               genes_list = NULL, decay = 1000, model = "Enhanced")
retina_rpmatrix <- rp_result[[1]]


### use python functions to convert the matrix to R readable dataframes ###
library(reticulate)
use_python("/home/mnt/nzh/nzhanglab/project/shuang/miniconda3/envs/scglue2/envs/MAESTRO/bin/python3.8",required = TRUE)
scipy <- import("scipy", convert = FALSE)
scipy$io$mmwrite("/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/MAESTRO/retina/atac/R_ver/retina_rpmatrix.mtx",retina_rpmatrix)
retinamatrix <-Matrix::readMM("/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/MAESTRO/retina/atac/R_ver/retina_rpmatrix.mtx")
row.names(retinamatrix) <- rp_result[[2]]
colnames(retinamatrix) <- colnames(retina.atac.mat)
retinamatrix ->maestro_retina_score
saveRDS(maestro_retina_score,"/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/MAESTRO/retina/atac/R_ver/MAESTRO_retina_rpmatrix.RDS")


### continue pipeline ###
retina.gene <- maestro_retina_score
retina.ATAC.res$ATAC <- ATACAttachGenescore(ATAC = retina.ATAC.res$ATAC, RPmatrix = retina.gene)
head(retina.ATAC.res$ATAC@meta.data)

data(human.immune.CIBERSORT)
retina.ATAC.res$ATAC <- ATACAnnotateCelltype(ATAC = retina.ATAC.res$ATAC,
                                                signatures = human.immune.CIBERSORT,
                                                min.score = 0.1,
                                                genes.test.use = "presto",
                                                genes.cutoff = 1E-5,
                                                outdir = "/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/MAESTRO/retina/atac/R_ver/")

retina.ATAC.res$ATAC<- AddMetaData(
  object = retina.ATAC.res$ATAC,
  metadata = meta1,
)


retina.ATAC.res$ATAC <- ATACAnnotateChromatinAccessibility(ATAC = retina.ATAC.res$ATAC,
                                                              peaks = retina.ATAC.res$peaks,
                                                              project = "retina_atac",
                                                              giggle.path = "/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/MAESTRO/annotation/giggle.all",
                                                              organism = "GRCh38")

saveRDS(retina.ATAC.res, "/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/MAESTRO/retina/atac/R_ver/retina_atac.rds")

### MAESTRO integration ###
MAESTRO integrate-init --rna-object /home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/MAESTRO/retina/gex/R_ver/retina_res.rds --atac-object /home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/MAESTRO/retina/atac/R_ver/retina_atac.rds --directory /home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/MAESTRO/retina/ --outprefix integrate
snakemake -np
snakemake -j 16

### doesn't work so I tried line by line to resolve ###
library(MAESTRO)
retina.RNA.res <- readRDS('/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/MAESTRO/retina/gex/R_ver/retina_res.rds')
retina.ATAC.res <- readRDS('/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/MAESTRO/retina/atac/R_ver/retina_atac.rds')

length(intersect(colnames(retina.RNA.res$RNA),colnames(retina.ATAC.res$ATAC)))
rna_input <-RenameCells(retina.RNA.res[["RNA"]], new.names = paste0("RNA_", colnames(x = retina.RNA.res[["RNA"]])))
retina.coembedded.res <- Incorporate(RNA = rna_input,
                                     ATAC = retina.ATAC.res$ATAC,
                                     project = "retina_integrated",
                                     method = "MAESTRO")

retina_embeddings <- retina.coembedded.res$CombinedObj@reductions$pca@cell.embeddings

retina_maestro <- readRDS("retina_integrate_Object.rds")

retina_embeddings <- retina_maestro@reductions$pca@cell.embeddings


saveRDS(retina_embeddings, "/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/MAESTRO/retina/retina_embeddings.rds")
saveRDS(retina.coembedded.res, "/home/mnt/nzh/nzhanglab/project/shuang/scATAC/comparison_methods/MAESTRO/retina/retina_integrated.rds")

