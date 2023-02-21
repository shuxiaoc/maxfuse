#harmony benchmark
library(Seurat)
library(Matrix)
library(matrixStats)
library(harmony)
# read in files
out_root = "/asap/output/"
in_root = "/asap/data/"
out_indx = 15

out_dir =paste0(out_root,"hm/")
in_dir = in_root
dir.create(out_root)
dir.create(out_dir)
# read

protein = read.csv(paste0(in_dir,"adt_pbmc.csv"))
protein = protein[,-which(names(protein) %in% c("X","barcode","CD4.1",'CD8a','CD11b.1'))]# not used channels
colnames(protein) = gsub('\\.','-', colnames(protein))
colnames(protein) = gsub('-$','', colnames(protein))

meta = read.csv(paste0(in_dir,"asap_pbmc_meta.csv"))

atacactivity = readMM(paste0(in_dir,"genescore_pbmc.txt"))
atacactivity = as.matrix(atacactivity)
gas_names = read.csv(paste0(in_dir ,'genescore_names_pbmc.csv'))
colnames(atacactivity) = gas_names$names

## remove
atacactivity = atacactivity[meta$human_ann != "dirt",]
protein = protein[meta$human_ann != "dirt",]
meta = meta[meta$human_ann != "dirt",]
##

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
    if (r %in% gas_names$names){
      rna_list = c(rna_list, r)
      protein_list = c(protein_list, protein_n)
    }
  }
}

act.shared = as.matrix(atacactivity[,rna_list[protein_list %in% colnames(protein)]]) # protein object
protein.shared = as.matrix(protein[,protein_list[protein_list %in% colnames(protein)]]) # rna object
colnames(protein.shared) = rna_list[protein_list %in% colnames(protein)] # make sure feature names same

# copy sp filtering to produce better output
act.shared.sub = act.shared[,colSds(act.shared)>0.5]
protein.shared.sub = protein.shared[,colSds(protein.shared)>0.1]
rownames(act.shared.sub) = paste0("d1",as.character(c(1:nrow(act.shared.sub))))
rownames(protein.shared.sub) = paste0("d2",as.character(c(1:nrow(protein.shared.sub))))
# then we construct the seurat objects
x_obj=CreateSeuratObject(counts=t(act.shared.sub),assay="x")
#x_obj <- NormalizeData(x_obj)
#x_obj <- FindVariableFeatures(x_obj, selection.method = "vst", nfeatures = 3000)
x_obj <- ScaleData(x_obj, features = rownames(x_obj))
# add suerat object datay
y_obj=CreateSeuratObject(counts=t(protein.shared.sub),assay="y")
y_obj <- NormalizeData(y_obj)
y_obj <- ScaleData(y_obj, features = rownames(y_obj))
#list_modality=list(x_obj,y_obj)
# get shared clean features
features=intersect(colnames(act.shared.sub),colnames(protein.shared.sub))
# run harmony in seurat, need to make a new seurat object
xy_obj = CreateSeuratObject(counts=cbind(t(act.shared.sub[,features]), t(protein.shared.sub[,features])))
xy_obj = SetAssayData(xy_obj, slot = "scale.data", cbind(x_obj@assays$x@scale.data[features,], y_obj@assays$y@scale.data[features,])) # takes very long
xy_obj = RunPCA(xy_obj, features = rownames(xy_obj), npcs = out_indx, verbose = FALSE)
xy_obj@meta.data$orig = c(rep("x",dim(act.shared.sub)[1]), rep("x",dim(protein.shared.sub)[1]))
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
