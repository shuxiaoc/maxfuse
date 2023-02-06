# full data set run for related spatial analysis
# for harmony

library(Seurat)
library(harmony)
library(Matrix)
library(matrixStats)

root_dir = '/tonsil_v2/'
out_dir = '/tonsil_v2/match/match_output/full/'
out_indx = 15

##
out_dir =paste0(out_dir,"hm/")
dir.create(out_dir)

rna = readMM(paste0(root_dir,"/RNA/tonsil_rna_0510.txt"))
protein = read.csv(paste0(root_dir,"/Codex/FCS_output_DeepCell_extOnly/formatch_clusters_x28_y715V2.csv"))

meta_rna = read.csv(paste0(root_dir,"/RNA/tonsil_rna_0510_meta.csv"))

names(protein)[names(protein) == 'collagen.IV'] <- 'collagen IV'
names(protein)[names(protein) == 'HLA.DR'] <- 'HLA DR'

rna_names = read.csv("/tonsil_v2/RNA/tonsil_rna_0510_names.csv") # rna names always the same
colnames(rna) = rna_names$names

#### for bsc
rownames(rna) = paste0("rna", c(1:nrow(rna)))
rownames(protein) = paste0("pro", c(1:nrow(protein)))

# change name
correspondence = read.csv('/tonsil_v2/match/protein_rna_name_conversionV11.csv')
correspondence = correspondence[!apply(correspondence == "", 1, all),]
rna_list = c()
protein_list = c()
for (j in c(1:dim(correspondence)[1])){
  protein_n = as.character(correspondence[j,1])
  rna_n = as.character(correspondence[j,2])
  if (grepl("Ignore", rna_n, fixed = TRUE)){
    next
  }
  rna_n = strsplit(rna_n, '/')[[1]]
  for(r in rna_n){
    if (r %in% rna_names$names){
      rna_list = c(rna_list, r)
      protein_list = c(protein_list, protein_n)
    }
  }
}
# get clean shared features
rna.shared = as.matrix(rna[,rna_list[protein_list %in% colnames(protein)]]) # protein object
protein.shared = as.matrix(protein[,protein_list[protein_list %in% colnames(protein)]]) # rna object
colnames(protein.shared) = rna_list[protein_list %in% colnames(protein)] # make sure feature names same

# copy sp filtering
rna.shared.sub = rna.shared[,colSds(rna.shared)>0.5]
protein.shared.sub = protein.shared[,colSds(protein.shared)>0.1]
rownames(rna.shared.sub) = as.character(c(1:nrow(rna.shared.sub)))
rownames(protein.shared.sub) = as.character(c(1:nrow(protein.shared.sub)))

# then we construct the seurat objects
x_obj=CreateSeuratObject(counts=t(rna.shared.sub),assay="x")
x_obj <- NormalizeData(x_obj)
#x_obj <- FindVariableFeatures(x_obj, selection.method = "vst", nfeatures = 3000) # no need to select variable genes in this case
x_obj <- ScaleData(x_obj, features = rownames(x_obj))
# add suerat object datay
y_obj=CreateSeuratObject(counts=t(protein.shared.sub),assay="y")
y_obj <- NormalizeData(y_obj)
y_obj <- ScaleData(y_obj, features = rownames(y_obj))
list_modality=list(x_obj,y_obj)
features=intersect(colnames(rna.shared.sub),colnames(protein.shared.sub))
# run harmony in seurat, need to make a new seurat object
xy_obj = CreateSeuratObject(counts=cbind(t(rna.shared.sub[,features]), t(protein.shared.sub[,features])))
#xy_obj = SetAssayData(xy_obj, slot = "scale.data", cbind(x_obj@assays$x@scale.data[features,], y_obj@assays$y@scale.data[features,])) # takes very long
xy_obj <- NormalizeData(xy_obj)
xy_obj = ScaleData(xy_obj)
xy_obj = RunPCA(xy_obj, features = rownames(xy_obj), npcs = 15, verbose = FALSE)
xy_obj@meta.data$orig = c(rep("x",dim(rna.shared.sub)[1]), rep("x",dim(protein.shared.sub)[1]))
# cbind together, scale within modality is better
xy_obj <- xy_obj %>% RunHarmony("orig")
embedding = Embeddings(xy_obj, 'harmony')[,c(1:out_indx)]
name_1 = "full_embed_x0.csv"
name_2 = "full_embed_y0.csv"
# does not directly produce matching info, produce later using knn with embeddning distance matrix
write.csv(embedding[c(1:ncol(x_obj)),c(1:out_indx)], paste0(out_dir,name_1),
        row.names=FALSE) # need to decide output pca cell
write.csv(embedding[c((ncol(x_obj) + 1):(ncol(x_obj) + ncol(y_obj))),c(1:out_indx)],
        paste0(out_dir,name_2), row.names=FALSE) # need to decide
write.csv(data.frame(method = "hm"), paste0(out_dir,"metrics.csv"), row.names=FALSE) 
