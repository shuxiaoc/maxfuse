# bindsc benchmark, full antibody panel
library(rliger)
library(Matrix)
library(matrixStats)
# read in files
out_root = "/bench_test3/output/"
in_root = "/bench_test3/input/"
batch = 5
out_indx = 15

for(i in c(1:batch)){
  batch_name = paste0("b",as.character(i),"/")
  out_dir =paste0(out_root,batch_name,"lg/")
  in_dir = paste0(in_root,batch_name)
  dir.create(paste0(out_root,batch_name))
  dir.create(out_dir)
  # read
  rna = readMM(paste0(in_dir,"rna.txt"))
  protein = read.csv(paste0(in_dir,"pro.csv"))
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
  if (dim(embedding)[1] != 20000) {
    break
  }
  name_1 = "full_embed_x0.csv"
  name_2 = "full_embed_y0.csv"
  # no avaliable matching information from liger thus not saved out
  # will use knn to serach matching on embedding in downstreatm analysis
  write.csv(embedding[c(1:nrow(rna.shared.sub)),c(1:out_indx)],
            paste0(out_dir,name_1), row.names=FALSE) # need to decide output pca cell
  write.csv(embedding[c((nrow(rna.shared.sub) + 1):(nrow(rna.shared.sub) + nrow(protein.shared.sub))),c(1:out_indx)],
            paste0(out_dir,name_2), row.names=FALSE) # need to decide
  write.csv(data.frame(method = "lg"), paste0(out_dir,"metrics.csv"), row.names=FALSE) 
}