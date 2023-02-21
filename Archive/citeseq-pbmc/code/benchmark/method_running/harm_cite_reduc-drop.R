#harmony for umap viz and cf matrix, drop antibody version
library(Seurat)
library(Matrix)
library(matrixStats)
library(harmony)
# read in files
out_root = "/bench_test3/output/reduction-drop/"
in_root = "/bench_test3/input/reduction/"

out_dir =paste0(out_root,"hm/")
in_dir = in_root
out_indx = 15

# read
dropped_pro = read.csv("/bench_test3/input/rank30.csv")
target = as.character(dropped_pro$target)

rna = readMM(paste0(in_dir,"rna.txt"))
protein = read.csv(paste0(in_dir,"pro.csv"))
protein = protein[,target]

meta = read.csv(paste0(in_dir,"meta.csv"))
rna_names = read.csv("/bench_test3/input/citeseq_rna_names.csv") # rna names always the same
colnames(rna) = rna_names$names
# change name
correspondence = read.csv('protein_rna_name_conversionV11.csv')
correspondence = correspondence[!apply(correspondence == "", 1, all),]
rna_list = c()
protein_list = c()
for (i in c(1:dim(correspondence)[1])){
  protein_n = as.character(correspondence[i,1])
  rna_n = as.character(correspondence[i,2])
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
# change name end
# first filtering step should be same as in sp
rna.shared = as.matrix(rna[,rna_list[protein_list %in% colnames(protein)]]) # protein object
protein.shared = as.matrix(protein[,protein_list[protein_list %in% colnames(protein)]]) # rna object
colnames(protein.shared) = rna_list[protein_list %in% colnames(protein)] # make sure feature names same
dim(protein.shared)
dim(rna.shared)
# copy sp filtering
rna.shared.sub = rna.shared[,colSds(rna.shared)>0.5]
protein.shared.sub = protein.shared[,colSds(protein.shared)>0.05]
rownames(rna.shared.sub) = as.character(c(1:nrow(rna.shared.sub)))
rownames(protein.shared.sub) = as.character(c(1:nrow(protein.shared.sub)))
#
dim(protein.shared.sub)
dim(rna.shared.sub)
#

# then we construct the seurat objects
x_obj=CreateSeuratObject(counts=t(rna.shared.sub),assay="x")
x_obj <- NormalizeData(x_obj)
x_obj <- FindVariableFeatures(x_obj, selection.method = "vst", nfeatures = 3000)
x_obj <- ScaleData(x_obj, features = rownames(x_obj))
# add suerat object datay
y_obj=CreateSeuratObject(counts=t(protein.shared.sub),assay="y")
y_obj <- NormalizeData(y_obj)
y_obj <- ScaleData(y_obj, features = rownames(y_obj))
#list_modality=list(x_obj,y_obj)
# get shared clean features
features=intersect(colnames(rna.shared.sub),colnames(protein.shared.sub))
# run harmony in seurat, need to make a new seurat object
xy_obj = CreateSeuratObject(counts=cbind(t(rna.shared.sub[,features]), t(protein.shared.sub[,features])))
#xy_obj = SetAssayData(xy_obj, slot = "scale.data", cbind(x_obj@assays$x@scale.data[features,], y_obj@assays$y@scale.data[features,])) # takes very long
xy_obj <- NormalizeData(xy_obj)
xy_obj <- ScaleData(xy_obj, features = rownames(xy_obj))
xy_obj <- RunPCA(xy_obj, features = rownames(xy_obj), npc = length(features))
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
