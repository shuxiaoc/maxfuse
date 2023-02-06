#seurat benchmark
library(Seurat)
library(Matrix)
library(matrixStats)
# read in files
out_root = "/abseq/output/"
in_root = "/abseq/data_prep/"
out_indx = 15

out_dir =paste0(out_root,"sr/")
in_dir = in_root
dir.create(out_root)
dir.create(out_dir)
# read
rna = readMM(paste0(in_dir,"abseqwta_rna.txt"))
rna = as.matrix(rna)
protein = read.csv(paste0(in_dir,"abseqwta_pro.csv"))
protein = protein[,-1]# remove the row index
meta = read.csv(paste0(in_dir,"abseqwta_meta.csv"))
rna_names = read.csv(paste0(in_dir,"abseqwta_rna_names.csv")) # rna names always the same
colnames(rna) = rna_names$names

rna = rna[meta$hm_annotate != "dirty",]
protein = protein[meta$hm_annotate != "dirty",]
meta = meta[meta$hm_annotate != "dirty",]

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
# copy sp filtering
rna.shared.sub = rna.shared[,colSds(rna.shared)>0.2]
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
# get transfer anchor
features=intersect(rownames(x_obj),rownames(y_obj))
pre.anchors <- FindTransferAnchors(reference = x_obj, query = y_obj, 
                                   dims = 1:20, features = features)
predictions <- TransferData(anchorset = pre.anchors, refdata = colnames(x_obj), 
                            dims = 1:20)
full_df = data.frame(idx2 = c(1:length(predictions$predicted.id)) -1,  idx1 = as.integer(predictions$predicted.id) -1,
                     score = predictions$prediction.score.max) # mind the r index difference
# get integration embedding
print("starting seurat integration")
Int.anchors <- FindIntegrationAnchors(object.list = list_modality,
                                      dims = 1:20, anchor.features =features, k.filter = 10)
xy_int <- IntegrateData(anchorset = Int.anchors, dims = 1:20, k.weight = 10)
# 
DefaultAssay(xy_int) <- "integrated"
xy_int <- ScaleData(xy_int, verbose = FALSE)
xy_int <- RunPCA(xy_int, npcs = out_indx, verbose = FALSE) # index of pca, 15 as fusion
embedding = xy_int@reductions$pca@cell.embeddings
name_1 = "full_embed_x0.csv"
name_2 = "full_embed_y0.csv"
#pathout = out_dir
write.csv(embedding[c(1:ncol(x_obj)),c(1:out_indx)], paste0(out_dir,name_1), row.names=FALSE) # need to decide output pca cell
write.csv(embedding[c((ncol(x_obj) + 1):(ncol(x_obj) + ncol(y_obj))),c(1:out_indx)],
          paste0(out_dir,name_2), row.names=FALSE) # need to decide
write.csv(full_df, paste0(out_dir,"full_idx.csv"), row.names=FALSE) # need to decide
write.csv(data.frame(method = "sr"), paste0(out_dir,"metrics.csv"), row.names=FALSE) 
