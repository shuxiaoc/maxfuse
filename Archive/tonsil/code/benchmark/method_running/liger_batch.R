#liger benchmark on 10k30k 5 batch cells, result used to produce matching accu and slt ari repeats
library(rliger)
library(Matrix)
library(matrixStats)
# read in files
out_root = "/tonsil_v2/match/bench_out/"
in_root = "/tonsil_v2/match/bench_input/"
batch = 5
out_indx = 15
temp_l1 = list()
temp_l2 = list()
`%notin%` <- Negate(`%in%`)

for(i in c(1:5)){
  batch_name = paste0("b",as.character(i),"/")
  out_dir =paste0(out_root,batch_name,"lg/")
  in_dir = paste0(in_root,batch_name)
  dir.create(paste0(out_root,batch_name))
  dir.create(out_dir)
  # read
  rna = readMM(paste0(in_dir,"rna.txt"))
  protein = read.csv(paste0(in_dir,"pro.csv"))
  meta_rna = read.csv(paste0(in_dir,"meta_rna.csv"))
  meta_pro = read.csv(paste0(in_dir,"meta_pro.csv"))
  
  # note this version caused name different, correct back
  names(protein)[names(protein) == 'collagen.IV'] <- 'collagen IV'
  names(protein)[names(protein) == 'HLA.DR'] <- 'HLA DR'
  
  rna_names = read.csv("/tonsil_v2/RNA/tonsil_rna_0510_names.csv") # rna names always the same
  colnames(rna) = rna_names$names
  # change name
  correspondence = read.csv('/tonsil_v2/match/protein_rna_name_conversionV11.csv')
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
  # copy sp filtering to produce better output
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
  #ligerobj <- selectGenes(ifnb_liger, var.thresh = 0, alpha.thresh=1)
  ligerobj@var.genes=features #  just use all
  ligerobj <- scaleNotCenter(ligerobj, remove.missing = FALSE)
  ligerobj <- optimizeALS(ligerobj, k = 20,remove.missing = FALSE)
  ligerobj <- quantile_norm(ligerobj)
  embedding = ligerobj@H.norm[,c(1:out_indx)]
  
  name_1 = "full_embed_x0.csv"
  name_2 = "full_embed_y0.csv"
  # no avaliable matching information from liger thus not saved out
  # will use knn to serach matching on embedding in downstreatm analysis
  
  # detect cells that was removed during liger step
  a1 = rownames(rna.shared.sub)[rownames(rna.shared.sub) %notin% rownames(embedding)]
  b1  = length(a1)
  #temp_l1[[i]] = a1
  
  a2 = rownames(protein.shared.sub)[rownames(protein.shared.sub) %notin% rownames(embedding)]
  b2  = length(a2)
  #temp_l2[[i]] = a2
  
  rn = nrow(rna.shared.sub)
  pn = nrow(protein.shared.sub)
  
  write.csv(embedding[c(1:(rn - b1)),c(1:out_indx)],
            paste0(out_dir,name_1), row.names=FALSE) # need to decide output pca cell
  write.csv(embedding[c((rn + 1 - b1):(rn + pn - b1 - b2)),c(1:out_indx)],
            paste0(out_dir,name_2), row.names=FALSE) # need to decide
  write.csv(data.frame(method = "lg"), paste0(out_dir,"metrics.csv"), row.names=FALSE) 
  
  # get updated orig
  ox = read.csv(paste0(in_dir,"orig_x.csv"))
  oy = read.csv(paste0(in_dir,"orig_y.csv"))
  a1_s = as.integer(gsub("d1", "", a1))
  a2_s = as.integer(gsub("d2", "", a2))
  
  if ( identical(a1_s, integer(0)) ){
    write.csv(ox, paste0(out_dir,'orig_x.csv'))
    write.csv(oy[-a2_s,], paste0(out_dir,'orig_y.csv'))
    write.csv(data.frame(id = a2_s),paste0(out_dir,'d2_id.csv'))
  }
  if ( identical(a2_s, integer(0)) ){
    write.csv(ox[-a1_s,], paste0(out_dir,'orig_x.csv'))
    write.csv(oy, paste0(out_dir,'orig_y.csv'))
    write.csv(data.frame(id = a1_s),paste0(out_dir,'d1_id.csv'))
  }
  else{
    write.csv(ox[-a1_s,], paste0(out_dir,'orig_x.csv'))
    write.csv(oy[-a2_s,], paste0(out_dir,'orig_y.csv'))
    write.csv(data.frame(id = a1_s),paste0(out_dir,'d1_id.csv'))
    write.csv(data.frame(id = a2_s),paste0(out_dir,'d2_id.csv'))
  }
  
}