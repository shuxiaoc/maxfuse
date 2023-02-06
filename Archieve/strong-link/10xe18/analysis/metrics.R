require(tidyverse)
# get_ovlp = function(x, k) {
#   x1 = x[1:k]
#   x2 = x[(k+1):length(x)]
#   return(length(intersect(x1, x2)) / k)
# }
# 
# get_local_struct_metric_each = function(X_full, X_reduced, k=100, eps=0.1, verbose=FALSE) {
#   browser()
#   print("Constructing k-NN graph for full data...")
#   nn_full = RANN::nn2(data=X_full, k=k, eps=eps)
#   print("Constructing k-NN graph for reduced data...")
#   nn_reduced = RANN::nn2(data=X_reduced, k=k, eps=eps)
#   print("Computing the metric...")
#   scores = apply(cbind(nn_full$nn.idx, nn_reduced$nn.idx), MARGIN=2, get_ovlp, k=k)
#   return(mean(scores))
# }
clust_labels_to_int = function(clust_labels) {
  uniq_labels = unique(clust_labels)
  clust_labels_int = rep(NA, length(clust_labels))
  
  for (ii in 1:length(clust_labels)) {
    clust_labels_int[ii] = match(clust_labels[ii], uniq_labels)
  }
  
  return(clust_labels_int)
}


compute_cor = function(x, method) {
  # Assume x is a vector of length 2n
  # x1 = the first half, and x2 = the second half.
  # This function returns the Kendall's tau between x1 and x2
  n = length(x)/2
  x1 = x[1:n]
  x2 = x[(n+1):length(x)]
  return(cor(x1, x2, method=method))
}

get_pairwise_euclid_dist = function(X) {
  X = as.matrix(X)
  
  get_norm_sq = function(x) {
    return(sum(x^2))
  }
  n = nrow(X)
  
  norm_sq_mat = apply(X, MARGIN=1, get_norm_sq)
  norm_sq_mat = matrix(norm_sq_mat, nrow=n, ncol=n, byrow=FALSE)
  dist_mat = norm_sq_mat + t(norm_sq_mat) - 2*X%*%t(X)
  dist_mat = sqrt(make_positive(dist_mat))
  
  return(dist_mat)
}

make_positive = function(X) {
  return(X - min(X))
}



get_struc_alignment_metric = function(X_full, X_reduced, convert_to_rank=TRUE, verbose=TRUE) {
  if (verbose) {
    print("Computing similarity matrix for full data...")
  }
  # sim_full = cor(t(X_full), method='pearson')
  sim_full = -get_pairwise_euclid_dist(X_full) 
  if (verbose) {
    print("Computing similarity matrix for reduced data...")
  }
  # sim_reduced = cor(t(X_reduced), method='pearson')
  sim_reduced = -get_pairwise_euclid_dist(X_reduced) 
  
  if (convert_to_rank) {
    if (verbose) {
      print("Converting to rank...")
    }
    sim_full = t(apply(sim_full, 1, rank, na.last=TRUE))
    sim_reduced = t(apply(sim_reduced, 1, rank, na.last=TRUE))
  }
  
  if (verbose) {
    print("Computing the metric...")
  }
  scores = apply(cbind(sim_full, sim_reduced), MARGIN=1, FUN=compute_cor, method='pearson')
  
  return(mean(scores))
}


sam = function(orig_fname='data/feature_deletion/orig', 
               embed_fname='data/feature_deletion/embedding', 
               n_idx=20, data_idx='x') {
  ### data_idx is either 'x' or 'y'.
  ### Say data_idx == 'x'.
  ### The original x data (i.e., cells, features, and a column 'label' of cluster labels) 
  ### is under orig_fname_x.csv.
  ### The embeddings are in embed_fname_xi.csv where
  ### 0 <= i <= n_idx.
  ### Calculate the return sam_x
  sam_vec = rep(0, n_idx+1)
  orig_data = read_csv(paste0(orig_fname, '_', data_idx, '.csv'), col_types = cols())
  orig_data = as.matrix(orig_data %>% select(-label))
  #print(dim(orig_data)) # add debug 0717 liger
  
  for (i in 0:n_idx) {
    #print(paste0("Now at iteration ", i))
    # read in x
    reduced_data = read_csv(paste0(embed_fname, '_', data_idx, i, '.csv'), col_types = cols())
    #print(dim(reduced_data)) # add debug 0717 liger
    sam_vec[i+1] = get_struc_alignment_metric(orig_data, reduced_data, verbose=FALSE)
    #print(sam_vec)
  }
  
  return(sam_vec)
}

slt = function(orig_fname='data/feature_deletion/orig', 
               embed_fname='data/feature_deletion/embedding', 
               n_idx=20) {
  ### The original x data (i.e., cells, features, and a column 'label' of cluster labels) 
  ### is under orig_fname_x.csv.
  ### The embeddings are in embed_fname_xi.csv where
  ### 0 <= i <= n_idx.
  ### The y data and embeddings are stored similarly.
  ### For each 0 <= i <= n_idx, calculate the following quantities:
  ###   1. slt_mix: one minus normalized Silhouette width with label being dataset index (x or y),
  ###               this is a measure of mixing
  ###   2. slt_clust: normalized Silhouette width with label being cell type clusters,
  ###                 this is a measure of how well the embeddings are if used to do clustering 
  ###   3. slt_f1: 2 * slt_mix * slt_clust / (slt_mix + slt_clust) --- F1 score of slt_mix and slt_clust
  ###
  ### Each one of slt_mix, slt_clust, slt_f1 is a length (n_idx+1) vector.
  ### Return cbind(slt_mix, slt_clust, slt_f1)
  
  SLT_MIN = -1
  SLT_MAX = 1
  slt_mix = rep(0, n_idx+1)
  slt_clust = rep(0, n_idx+1)
  slt_f1 = rep(0, n_idx+1)
  
  # cell type cluster labels
  x_label_vec = read_csv(paste0(orig_fname, '_x.csv'), col_types = cols()) %>% select(label) %>% pull()
  y_label_vec = read_csv(paste0(orig_fname, '_y.csv'), col_types = cols()) %>% select(label) %>% pull()
  clust_label = c(x_label_vec, y_label_vec)
  clust_label = clust_labels_to_int(clust_label)
  
  
  # dataset identifiers 
  #data_label = c(rep(1, length(x_label_vec)), rep(2, length(y_label_vec)))
  data_label = c(rep(1, length(x_label_vec)), rep(2, length(y_label_vec)))
  data_label = clust_labels_to_int(data_label)
  print(paste0("class datalabel is ", class(data_label)))
  print(paste0("length datalabel is ", length(data_label)))
  
  for (i in 0:n_idx) {
    reduced_data_x = as.matrix(read_csv(paste0(embed_fname, '_x', i, '.csv'), col_types = cols()))
    reduced_data_y = as.matrix(read_csv(paste0(embed_fname, '_y', i, '.csv'), col_types = cols()))
    reduced_data = rbind(reduced_data_x, reduced_data_y)
    dist_mat = get_pairwise_euclid_dist(reduced_data) 
    # calculate metrics
    curr_slt_mix = 1 - (cluster::silhouette(data_label, dmatrix=dist_mat)[, 3] - SLT_MIN) / (SLT_MAX - SLT_MIN)
    curr_slt_clust = (cluster::silhouette(clust_label, dmatrix=dist_mat)[, 3] - SLT_MIN) / (SLT_MAX - SLT_MIN)
    curr_slt_f1 = 2 * curr_slt_mix * curr_slt_clust / (curr_slt_mix + curr_slt_clust)
    slt_mix[i+1] = mean(curr_slt_mix)
    slt_clust[i+1] = mean(curr_slt_clust)
    slt_f1[i+1] = mean(curr_slt_f1)
  }
  
  return(cbind(slt_mix, slt_clust, slt_f1))
}


ari = function(orig_fname='data/feature_deletion/orig', 
               embed_fname='data/feature_deletion/embedding', 
               n_idx=20, kmeans_nsim=10) {
  ### The original x data (i.e., cells, features, and a column 'label' of cluster labels) 
  ### is under orig_fname_x.csv.
  ### The embeddings are in embed_fname_xi.csv where
  ### 0 <= i <= n_idx.
  ### The y data and embeddings are stored similarly.
  ### For each 0 <= i <= n_idx, calculate the following quantities:
  ###   1. ari_mix: one minus normalized adjusted random index between the k-means clustering result (with two clusters)
  ###               and the dataset index (x and y). This is a measure of mixing.
  ###   2. ari_clust: normalized adjusted random index between the k-means clustering result (with ground truth number of clusters) 
  ###                 and the ground truth cell type cluster labels.
  ###                 This is a measure of how well the embeddings are if used to do clustering 
  ###   3. ari_f1: 2 * ari * slt_clust / (ari_mix + ari_clust) --- F1 score of ari_mix and ari_clust
  ###
  ### Since the k-means implementation in R is highly dependent on random initializations, we do k-means for kmeans_nsim times
  ### and average the result
  ### Each one of ari_mix, ari_clust, ari_f1 is a length (n_idx+1) vector.
  ### Return cbind(ari_mix, ari_clust, ari_f1)
  
  ARI_MIN = -1
  ARI_MAX = 1
  ari_mix = rep(0, n_idx+1)
  ari_clust = rep(0, n_idx+1)
  ari_f1 = rep(0, n_idx+1)
  
  # cell type cluster labels
  x_label_vec = read_csv(paste0(orig_fname, '_x.csv'), col_types = cols()) %>% select(label) %>% pull()
  y_label_vec = read_csv(paste0(orig_fname, '_y.csv'), col_types = cols()) %>% select(label) %>% pull()
  clust_label = c(x_label_vec, y_label_vec)
  clust_label = clust_labels_to_int(clust_label)
  n_clusters = length(unique(clust_label))
  
  # dataset identifiers 
  data_label = c(rep(1, length(x_label_vec)), rep(2, length(y_label_vec)))
  
  for (i in 0:n_idx) {
    reduced_data_x = as.matrix(read_csv(paste0(embed_fname, '_x', i, '.csv'), col_types = cols()))
    reduced_data_y = as.matrix(read_csv(paste0(embed_fname, '_y', i, '.csv'), col_types = cols()))
    reduced_data = rbind(reduced_data_x, reduced_data_y)
    curr_avg_ari_mix = 0
    curr_avg_ari_clust = 0
    curr_avg_ari_f1 = 0
    for (j in 1:kmeans_nsim) {
      # calculate metrics
      est_data_label = kmeans(reduced_data, centers=2)$cluster
      est_clust_label = kmeans(reduced_data, centers=n_clusters)$cluster
      curr_ari_mix = 1 - (mclust::adjustedRandIndex(est_data_label, data_label) - ARI_MIN)  / (ARI_MAX - ARI_MIN)
      curr_ari_clust = (mclust::adjustedRandIndex(est_clust_label, clust_label) - ARI_MIN)  / (ARI_MAX - ARI_MIN)
      curr_ari_f1 = 2 * curr_ari_mix * curr_ari_clust / (curr_ari_mix + curr_ari_clust)
      # average over different kmeans runs
      curr_avg_ari_mix = curr_avg_ari_mix + curr_ari_mix / kmeans_nsim
      curr_avg_ari_clust = curr_avg_ari_clust + curr_ari_clust / kmeans_nsim
      curr_avg_ari_f1 = curr_avg_ari_f1 + curr_ari_f1 / kmeans_nsim
    }
    # save the results
    ari_mix[i+1] = curr_avg_ari_mix
    ari_clust[i+1] = curr_avg_ari_clust
    ari_f1[i+1] = curr_avg_ari_f1 
  }
  
  return(cbind(ari_mix, ari_clust, ari_f1))
}


kbet = function(orig_fname='data/feature_deletion/orig', 
                embed_fname='data/feature_deletion/embedding', 
                n_idx=20, kmeans_nsim=10) {
  ### The original x data (i.e., cells, features, and a column 'label' of cluster labels) 
  ### is under orig_fname_x.csv.
  ### The embeddings are in embed_fname_xi.csv where
  ### 0 <= i <= n_idx.
  ### The y data and embeddings are stored similarly.
  ### For each 0 <= i <= n_idx, calculate the kBET test statistic and return it.
  ### We try different 'testSize' parameter in {10, 15, 20}, and we average all the results to get the final result.
  
  kbet_vec = rep(0, n_idx+1)
  
  # cell type cluster labels
  x_label_vec = read_csv(paste0(orig_fname, '_x.csv'), col_types = cols()) %>% select(label) %>% pull()
  y_label_vec = read_csv(paste0(orig_fname, '_y.csv'), col_types = cols()) %>% select(label) %>% pull()
  clust_label = c(x_label_vec, y_label_vec)
  unique_labels = unique(clust_label)
  
  # dataset identifiers 
  data_label = c(rep(1, length(x_label_vec)), rep(2, length(y_label_vec)))
  
  # hyperparams
  sizes = c(10, 15, 20)
  
  for (i in 0:n_idx) {
    reduced_data_x = as.matrix(read_csv(paste0(embed_fname, '_x', i, '.csv'), col_types = cols()))
    reduced_data_y = as.matrix(read_csv(paste0(embed_fname, '_y', i, '.csv'), col_types = cols()))
    reduced_data = rbind(reduced_data_x, reduced_data_y)
    
    curr_kbet = c() 
    # for (curr_label in unique_labels) {
      # mask = (clust_label == curr_label)
      # curr_reduced_data = reduced_data[mask, ]
      # curr_data_label = data_label[mask]
      
      # if (length(curr_data_label) < 5) { 
        # too few samples, do not count
        # next
      # }
      for (curr_size in sizes) {
        # calculate kBET test statistic
        curr_kbet_trial = 1 - as.numeric(
          kBET::kBET(
            df = reduced_data, batch = data_label, do.pca = F, 
            addTest = F, testSize = curr_size, n_repeat = 10, plot = F, verbose = F
          )$summary[1, 2]
        )
        # print(curr_kbet_trial)
        # accumulate the result if not NA
        curr_kbet = c(curr_kbet, curr_kbet_trial)
      }
    # }
    kbet_vec[i+1] = mean(curr_kbet, na.rm=TRUE)
  }
  
  return(kbet_vec)
}


lisi = function(orig_fname='data/feature_deletion/orig', 
                embed_fname='data/feature_deletion/embedding', 
                n_idx=20) {
  ### The original x data (i.e., cells, features, and a column 'label' of cluster labels) 
  ### is under orig_fname_x.csv.
  ### The embeddings are in embed_fname_xi.csv where
  ### 0 <= i <= n_idx.
  ### The y data and embeddings are stored similarly.
  ### For each 0 <= i <= n_idx, calculate the following quantities:
  ###   1. lisi_mix: LISI (NOT NORMALIED) where the label is given by the dataset indices (x or y). This is a measure of mixing.
  ###   2. lisi_clust: 1 - LISI (NOT NORMALIZED) where the label is given by the 
  ###                  ground truth cell type cluster labels.
  ###                 This is a measure of how well the embeddings are if used to do clustering 
  ### Note that there is no a priori bound on the range of LISI, so there is no canonical way to normalized LISI.
  ### Thus, we do not calculate the "F1" score.
  ###
  ### Each one of lisi_mix, lisi_clust is a length (n_idx+1) vector.
  ### Return cbind(lisi_mix, lisi_clust)
  
  lisi_mix = rep(0, n_idx+1)
  lisi_clust = rep(0, n_idx+1)
  
  # cell type cluster labels
  x_label_vec = read_csv(paste0(orig_fname, '_x.csv'), col_types = cols()) %>% select(label) %>% pull()
  y_label_vec = read_csv(paste0(orig_fname, '_y.csv'), col_types = cols()) %>% select(label) %>% pull()
  clust_label = c(x_label_vec, y_label_vec)
  clust_label = clust_labels_to_int(clust_label)
  
  # dataset identifiers 
  data_label = c(rep(1, length(x_label_vec)), rep(2, length(y_label_vec)))
  
  labels = as.data.frame(cbind(data_label, clust_label))
  colnames(labels) = c('data_label', 'clust_label')
  
  for (i in 0:n_idx) {
    reduced_data_x = as.matrix(read_csv(paste0(embed_fname, '_x', i, '.csv'), col_types = cols()))
    reduced_data_y = as.matrix(read_csv(paste0(embed_fname, '_y', i, '.csv'), col_types = cols()))
    reduced_data = as.data.frame(rbind(reduced_data_x, reduced_data_y))
    
    # calculate LISI
    lisi_res <- lisi::compute_lisi(reduced_data, labels, c('data_label', 'clust_label'), perplexity = 40)
    # save the results
    lisi_mix[i+1] = median(lisi_res$data_label)
    lisi_clust[i+1] = median(1 - lisi_res$clust_label)
  }
  
  return(cbind(lisi_mix, lisi_clust))
}



mix = function(orig_fname='data/feature_deletion/orig', 
               embed_fname='data/feature_deletion/embedding', 
               n_idx=20) {
  ### The original x and y data (i.e., cells, features, and a column 'label' of cluster labels) 
  ### are under orig_fname_x.csv, orig_fname_y.csv, respectively.
  ### The embeddings are in embed_fname_xi.csv and embed_fname_yi.csv where
  ### 0 <= i <= n_idx.
  ### Calculate the return avg_mix. 
  
  mix_vec = rep(0, n_idx+1)
  x_labels = read_csv(paste0(orig_fname, '_x.csv'), col_types = cols()) %>% select(label) %>% pull()
  y_labels = read_csv(paste0(orig_fname, '_y.csv'), col_types = cols()) %>% select(label) %>% pull()
  uniq_labels = unique(c(x_labels, y_labels))
  
  for (i in 0:n_idx) {
    # print(paste0("Now at iteration ", i))
    # read in embeddings
    x_embed = read_csv(paste0(embed_fname, '_x', i, '.csv'), col_types = cols())
    y_embed = read_csv(paste0(embed_fname, '_y', i, '.csv'), col_types = cols())
    for (k in 1:length(uniq_labels)) {
      # calculate mixing for each cluster
      cluster = uniq_labels[k]
      p = ncol(x_embed)
      stat_dist = rep(NA, p)
      for (j in 1:p) {
        # do two sample testing for each column of the embeddings 
        # use 1 - dist to make sure the larger, the better
        stat_dist[j] = 1 - twosamples::ks_stat(x_embed[x_labels==cluster, j]%>%pull(), y_embed[y_labels==cluster, j]%>%pull())
      }
      # the median of the test statistics for all columns is taken as the mixing metric for this cluster 
      mix_vec[i+1] = (mix_vec[i+1] + median(stat_dist)) / length(uniq_labels)
    }
  }
  return(mix_vec)
}
