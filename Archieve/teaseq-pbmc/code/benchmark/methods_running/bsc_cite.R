#bindsc benchmark
library(Seurat)
library(Matrix)
library(matrixStats)
# read in files
out_root = "/ICICLE/output/"
in_root = "/ICICLE/data/"
out_indx = 15

out_dir =paste0(out_root,"bsc/")
in_dir = in_root
dir.create(out_root)
dir.create(out_dir)

protein = read.csv(paste0(in_dir,"adt.csv"))
colnames(protein) = gsub('\\.','-', colnames(protein)) # change name formatting
colnames(protein) = gsub('-$','', colnames(protein))
protein$cell_barcode <- NULL
protein$total <- NULL

meta = read.csv(paste0(in_dir,"atac_meta.csv"))

atacactivity = readMM(paste0(in_dir,"genescore_tea.txt"))
atacactivity = as.matrix(atacactivity)
gas_names = read.csv(paste0(in_dir ,'genescore_names_tea.csv'))
colnames(atacactivity) = gas_names$names

#### for bsc
rownames(atacactivity) = paste0("act", c(1:nrow(atacactivity)))
rownames(protein) = paste0("pro", c(1:nrow(protein)))

# change name
correspondence = read.csv('conversion_v12.csv')# updated to v12 as there are some differently named proteins
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

# copy sp filtering
act.shared.sub = act.shared[,colSds(act.shared)>0.36 & colSds(protein.shared)>3.6]
protein.shared.sub = protein.shared[,colSds(act.shared)>0.36 & colSds(protein.shared)>3.6]
rownames(act.shared.sub) = as.character(c(1:nrow(act.shared.sub)))
rownames(protein.shared.sub) = as.character(c(1:nrow(protein.shared.sub)))

# get cluster for bindsc x, using all x features
xc_obj=CreateSeuratObject(counts=t(protein),assay="x")
xc_obj <- NormalizeData(xc_obj)
xc_obj <- ScaleData(xc_obj, features = rownames(xc_obj))
xc_obj <- RunPCA(xc_obj, features = rownames(xc_obj))
xc_obj <- FindNeighbors(xc_obj, dims = 1:15)
xc_obj <- FindClusters(xc_obj, resolution = 1)
x_cluster = as.factor(paste0('x_',as.character(Idents(xc_obj))))

# get cluster for bindsc x, using all x features
x_obj=CreateSeuratObject(counts=t(protein.shared.sub),assay="x")
x_obj <- NormalizeData(x_obj)
x_obj <- ScaleData(x_obj, features = rownames(x_obj))# not used

# get cluster for bindsc y, using all y features (variable)
y_obj=CreateSeuratObject(counts=t(atacactivity),assay="y")
#y_obj <- NormalizeData(y_obj)
y_obj <- ScaleData(y_obj, features = rownames(y_obj))
y_obj <- FindVariableFeatures(y_obj, nfeatures = 3000)
y_obj <- RunPCA(y_obj, features = VariableFeatures(object = y_obj))
y_obj <- FindNeighbors(y_obj, dims = 1:25) # based on elbow plot
y_obj <- FindClusters(y_obj, resolution = 1)
y_cluster = as.factor(paste0('y_',as.character(Idents(y_obj))))

y_input_features = VariableFeatures(object = y_obj)

## for Z0
z_obj=CreateSeuratObject(counts=t(act.shared.sub),assay="z")
#z_obj <- NormalizeData(z_obj)
z_obj <- ScaleData(z_obj)

## now gather all the actual inputs
x_input = x_obj@assays$x@data
y_input = as.matrix(as.data.frame(y_obj@assays$y@data[y_input_features,]))
z0_input = z_obj@assays$z@data


library(bindSC)
# start bindsc
res <- BiCCA( X = x_input ,
              Y =  y_input, 
              Z0 = z0_input, 
              X.clst = x_cluster,
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

name_1 = "full_embed_x0.csv"
name_2 = "full_embed_y0.csv"
pathout = out_dir
write.csv(data.frame(res$r)[,c(1:out_indx)], paste0(out_dir,name_1), row.names=FALSE) # rna embed
write.csv(data.frame(res$u)[,c(1:out_indx)], paste0(out_dir,name_2), row.names=FALSE) # pro embed
write.csv(data.frame(method = "bsc"), paste0(out_dir,"metrics.csv"), row.names=FALSE) 
