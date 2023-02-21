#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

metrics_fname = args[1]
orig_fname = args[2]
embed_fname = args[3]
n_idx = as.integer(args[4]) 

# Compute the following metrics:
# - sam_x: structure alignment metric for x data (the larger, the better)
# - sam_y: structure alignment metric for y data (the larger, the better)
# - slt_mix: mixing via Silhouette width (the larger, the better)
# - slt_clust: quality of embeddings for clustering via Silhouette width (the larger, the better)
# - slt_f1: an integrated metric using both slt_mix and slt_clust (the larger, the better)
# - ari_mix: mixing via adjusted random index (the larger, the better)
# - ari_clust: quality of embeddings for clustering via adjusted random index (the larger, the better)
# - lisi_mix: mixing via Local Inverse Simpson’s Index (LISI) (the larger, the better)
# - lisi_clust: quality of embeddings for clustering via LISI (the larger, the better)
# - kbet: mixing via k-nearest neighbour batch effect test (kBET) (the larger, the better)
# - avg_mix: mixing metric via two sample test, averaged over all clusters (the larger, the better)
setwd("./")
source("metrics.R")

# load existing metrics
metrics = read_csv(metrics_fname, col_types=cols())


# calculate structure alignment metrics
print(paste0(format(Sys.Date(), "%c"), ': calculating structure alignment metrics...'))
sam_x = sam(orig_fname=orig_fname, embed_fname=embed_fname,
            n_idx=n_idx, data_idx='x')
sam_y= sam(orig_fname=orig_fname, embed_fname=embed_fname,
           n_idx=n_idx, data_idx='y')
#print(sam_x)
#print(sam_y)
metrics = metrics %>% add_column(sam_x=sam_x) %>% add_column(sam_y=sam_y)
#print(metrics)
# calculate Silhouette width
print(paste0(format(Sys.Date(), "%c"), ': calculating Silhouette width...'))
slt_res = slt(orig_fname=orig_fname, embed_fname=embed_fname, n_idx=n_idx)
#print(slt_res)
metrics = metrics %>% add_column(slt_mix=slt_res[, 1]) %>% add_column(slt_clust=slt_res[, 2]) %>% add_column(slt_f1=slt_res[, 3])
#print(metrics)
# calculate ARI
print(paste0(format(Sys.Date(), "%c"), ': calculating adjusted random index...'))
ari_res = ari(orig_fname=orig_fname, embed_fname=embed_fname, n_idx=n_idx)
metrics = metrics %>% add_column(ari_mix=ari_res[, 1]) %>% add_column(ari_clust=ari_res[, 2]) %>% add_column(ari_f1=ari_res[, 3])

# calculate LISI
print(paste0(format(Sys.Date(), "%c"), ': calculating Local Inverse Simpson’s Index...'))
lisi_res = lisi(orig_fname=orig_fname, embed_fname=embed_fname, n_idx=n_idx)
metrics = metrics %>% add_column(lisi_mix=lisi_res[, 1]) %>% add_column(lisi_clust=lisi_res[, 2])

# calculate mixing averaged over clusters
print(paste0(format(Sys.Date(), "%c"), ': calculating mixing quality...'))
avg_mix = mix(orig_fname=orig_fname, embed_fname=embed_fname, 
              n_idx=n_idx) 
metrics = metrics %>% add_column(avg_mix=avg_mix)

# save metrics, because the calculation of kBET is substantially slower.
write_csv(metrics, metrics_fname)
#print(paste0(format(Sys.Date(), "%c"), ': nearly done...'))

#### not calculating kBet here because too slow for this stage
# calculate kBET
#print(paste0(format(Sys.Date(), "%c"), ': calculating kBET...'))
#kbet_res = kbet(orig_fname=orig_fname, embed_fname=embed_fname, n_idx=n_idx)
#metrics = metrics %>% add_column(kBET=kbet_res)

#write_csv(metrics, metrics_fname)
#print(paste0(format(Sys.Date(), "%c"), ': done!'))