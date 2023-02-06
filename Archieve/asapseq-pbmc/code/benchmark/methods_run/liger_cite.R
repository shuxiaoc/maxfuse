#liger benchmark
library(rliger)
library(Matrix)
library(matrixStats)
# read in files
out_root = "/asap/output/"
in_root = "/asap/data/"
out_indx = 15

out_dir =paste0(out_root,"lg/")
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
# then we construct the liger objects
ligerobj=createLiger( list(x = t(act.shared.sub), y = t(protein.shared.sub)), remove.missing = FALSE)
###Start integration
features=intersect(colnames(act.shared.sub),colnames(protein.shared.sub)) # shared features accross datasets with good quality
# default preprocessing
ligerobj <- rliger::normalize(ligerobj, remove.missing = FALSE)
# do not need to select genes
#ligerobj <- selectGenes(ifnb_liger, var.thresh = 0, alpha.thresh=1)
ligerobj@var.genes=features #  just use all
ligerobj <- scaleNotCenter(ligerobj, remove.missing = FALSE)
ligerobj <- optimizeALS(ligerobj, k = 20,remove.missing = TRUE)
ligerobj <- quantile_norm(ligerobj)
embedding = ligerobj@H.norm[,c(1:out_indx)]
name_1 = "full_embed_x0.csv"
name_2 = "full_embed_y0.csv"
# no avaliable matching information from liger thus not saved out
# will use knn to serach matching on embedding in downstreatm analysis
# check what cell is filtered out 
`%notin%` <- Negate(`%in%`)
filtered = 
  c(rownames(act.shared.sub), rownames(protein.shared.sub))[c(rownames(act.shared.sub), rownames(protein.shared.sub)) %notin% rownames(ligerobj@H.norm)]
filtered_id = as.integer(gsub("d1", "", filtered)) # no cells filtered

write.csv(embedding[c(1:4360),], # no cells filltered out dump way to save out
          paste0(out_dir,name_1), row.names=FALSE) # need to decide output pca cell
write.csv(embedding[c(4361:8720),],
          paste0(out_dir,name_2), row.names=FALSE) # need to decide
write.csv(data.frame(method = "lg"), paste0(out_dir,"metrics.csv"), row.names=FALSE) 

