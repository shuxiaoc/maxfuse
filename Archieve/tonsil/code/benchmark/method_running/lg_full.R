# full data set run for related spatial analysis
# for liger

library(rliger)
library(Matrix)
library(matrixStats)

root_dir = '/tonsil_v2/'
out_dir = '/tonsil_v2/match/match_output/full/'
out_indx = 15
`%notin%` <- Negate(`%in%`)

out_dir =paste0(out_dir,"lg/")
dir.create(out_dir)

rna = readMM(paste0(root_dir,"/RNA/tonsil_rna_0510.txt"))
protein = read.csv(paste0(root_dir,"/Codex/FCS_output_DeepCell_extOnly/formatch_clusters_x28_y715V2.csv"))

meta_rna = read.csv(paste0(root_dir,"/RNA/tonsil_rna_0510_meta.csv"))

names(protein)[names(protein) == 'collagen.IV'] <- 'collagen IV'
names(protein)[names(protein) == 'HLA.DR'] <- 'HLA DR'

rna_names = read.csv("/tonsil_v2/RNA/tonsil_rna_0510_names.csv") # rna names always the same
colnames(rna) = rna_names$names

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
rownames(rna.shared.sub) = paste0("d1",as.character(c(1:nrow(rna.shared.sub))))
rownames(protein.shared.sub) = paste0("d2",as.character(c(1:nrow(protein.shared.sub))))

# then we construct the liger objects
ligerobj=createLiger( list(x = t(rna.shared.sub), y = t(protein.shared.sub)), remove.missing = FALSE)
###Start integration
features=intersect(colnames(rna.shared.sub),colnames(protein.shared.sub)) # shared features accross datasets with good quality
# default preprocessing
ligerobj <- rliger::normalize(ligerobj, remove.missing = FALSE)
# do not need to select genes
#ligerobj <- selectGenes(ligerobj, var.thres= 0,unshared = TRUE,
#                        unshared.datasets = list(2), unshared.thresh= 0, alpha.thresh = 1) # unimf version of liger
ligerobj@var.genes=features #  only used for length 
ligerobj <- scaleNotCenter(ligerobj, remove.missing = FALSE)
ligerobj <- optimizeALS(ligerobj, k = 20,remove.missing = FALSE)
#ligerobj <- optimizeALS(ligerobj, use.unshared = TRUE, k = 20,remove.missing = FALSE)
ligerobj <- quantile_norm(ligerobj)
embedding = ligerobj@H.norm[,c(1:out_indx)]

name_1 = "full_embed_x0.csv"
name_2 = "full_embed_y0.csv"
# no avaliable matching information from liger thus not saved out
# will use knn to serach matching on embedding in downstreatm analysis


# before proceed, make sure what cells got deleted
a1 = rownames(rna.shared.sub)[rownames(rna.shared.sub) %notin% rownames(embedding)]
b1  = length(a1) # 39 rna cells got removed

a2 = rownames(protein.shared.sub)[rownames(protein.shared.sub) %notin% rownames(embedding)]
b2  = length(a2) # 8 cdx cells got removed

rn = nrow(rna.shared.sub)
pn = nrow(protein.shared.sub)

write.csv(embedding[c(1:(rn - b1)),c(1:out_indx)],
          paste0(out_dir,name_1), row.names=FALSE) # need to decide output pca cell
write.csv(embedding[c((rn + 1 - b1):(rn + pn - b1 - b2)),c(1:out_indx)],
          paste0(out_dir,name_2), row.names=FALSE) # need to decide
write.csv(data.frame(method = "lg"), paste0(out_dir,"metrics.csv"), row.names=FALSE) 

## get ids of missing cells

a1_s = as.integer(gsub("d1", "", a1))
a2_s = as.integer(gsub("d2", "", a2))

write.csv(data.frame(id = a1_s),paste0(out_dir,'d1_id.csv'))
write.csv(data.frame(id = a2_s),paste0(out_dir,'d2_id.csv'))
