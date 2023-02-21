# full data set run for related spatial analysis
# for bindsc
# too many cells BindSC can not handle
# split data into two as input, then combine the results later

library(Seurat)
library(bindSC)
library(Matrix)
library(matrixStats)

root_dir = '/tonsil_v2/'
out_dir = '/tonsil_v2/match/match_output/'

rna = readMM(paste0(root_dir,"/RNA/tonsil_rna_0510.txt"))
protein = read.csv(paste0(root_dir,"/Codex/FCS_output_DeepCell_extOnly/formatch_clusters_x28_y715V2.csv"))

meta_rna = read.csv(paste0(root_dir,"/RNA/tonsil_rna_0510_meta.csv"))

names(protein)[names(protein) == 'collagen.IV'] <- 'collagen IV' # quick name change
names(protein)[names(protein) == 'HLA.DR'] <- 'HLA DR'

rna_names = read.csv("//tonsil_v2/RNA/tonsil_rna_0510_names.csv") # rna names always the same
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
rna.shared.sub = rna.shared[,colSds(rna.shared)>0.5 & colSds(protein.shared)>0.1]
protein.shared.sub = protein.shared[,colSds(rna.shared)>0.5 & colSds(protein.shared)>0.1]

# get cluster for bindsc x, using all x features
x_obj=CreateSeuratObject(counts=t(protein.shared.sub),assay="x")
x_obj <- NormalizeData(x_obj)
x_obj <- ScaleData(x_obj, features = rownames(x_obj))
x_obj <- RunPCA(x_obj, features = rownames(x_obj))
x_obj <- FindNeighbors(x_obj, dims = 1:15)
x_obj <- FindClusters(x_obj, resolution = 1)
x_cluster = as.factor(paste0('x_',as.character(Idents(x_obj))))

# get cluster for bindsc y, using all y features (variable)
y_obj=CreateSeuratObject(counts=t(rna),assay="y")
y_obj <- NormalizeData(y_obj)
y_obj <- ScaleData(y_obj, features = rownames(y_obj))
y_obj <- FindVariableFeatures(y_obj, nfeatures = 2000)
y_obj <- RunPCA(y_obj, features = VariableFeatures(object = y_obj))
y_obj <- FindNeighbors(y_obj, dims = 1:15)
y_obj <- FindClusters(y_obj, resolution = 1)
y_cluster = as.factor(paste0('y_',as.character(Idents(y_obj))))

y_input_features = VariableFeatures(object = y_obj)

## for Z0
z_obj=CreateSeuratObject(counts=t(rna.shared.sub),assay="z")
z_obj <- NormalizeData(z_obj)

## now gather all the actual inputs
x_input = x_obj@assays$x@data
y_input = as.matrix(as.data.frame(y_obj@assays$y@data[y_input_features,]))
z0_input = z_obj@assays$z@data


###### problem: bindsc can not do on all cells, too large
###### to help it lets cut the codex data in half and do it
cdx_cells = dim(x_input)[2]
x_input_p1 = x_input[,c(1:89459)]
x_input_p2 = x_input[,c(89460:178919)]


# start bindsc part1
res <- BiCCA( X = x_input_p1 ,
              Y =  y_input, 
              Z0 = z0_input, 
              X.clst = x_cluster[c(1:89459)],
              Y.clst = y_cluster,
              alpha = 0.1, 
              lambda = 0.7,
              K = 15,
              temp.path  = "out",
              num.iteration = 50,
              tolerance = 0.01,
              save = TRUE,
              parameter.optimize = FALSE,
              block.size = 0)

name_1 = "BSC_full_embed_x0_finV_p1.csv"
name_2 = "BSC_full_embed_y0_finV_p1.csv"
pathout = out_dir

out_indx = 15
write.csv(data.frame(res$r)[,c(1:out_indx)], paste0(out_dir,name_1), row.names=FALSE) # rna embed
write.csv(data.frame(res$u)[,c(1:out_indx)], paste0(out_dir,name_2), row.names=FALSE) # pro embed


# start bindsc part1
res <- BiCCA( X = x_input_p2 ,
              Y =  y_input, 
              Z0 = z0_input, 
              X.clst = x_cluster[c(89460:178919)],
              Y.clst = y_cluster,
              alpha = 0.1, 
              lambda = 0.7,
              K = 15,
              temp.path  = "out",
              num.iteration = 50,
              tolerance = 0.01,
              save = TRUE,
              parameter.optimize = FALSE,
              block.size = 0)

name_1 = "BSC_full_embed_x0_finV_p2.csv"
name_2 = "BSC_full_embed_y0_finV_p2.csv"
pathout = out_dir

out_indx = 15
write.csv(data.frame(res$r)[,c(1:out_indx)], paste0(out_dir,name_1), row.names=FALSE) # rna embed
write.csv(data.frame(res$u)[,c(1:out_indx)], paste0(out_dir,name_2), row.names=FALSE) # pro embed


