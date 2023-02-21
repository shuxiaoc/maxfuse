# no condo env requirement
# used to calc slt and ari for all methods

# for mf
/usr/bin/Rscript calculate_metrics.R '/abseq/output/mf/metrics.csv' '/abseq/data_prep/orig' '/abseq/output/mf/full_embed' 0 &

# for sr
/usr/bin/Rscript calculate_metrics.R '/abseq/output/sr/metrics.csv' '/abseq/data_prep/orig' '/abseq/output/sr/full_embed' 0 &

# for lg
/usr/bin/Rscript calculate_metrics.R '/abseq/output/lgunimf/metrics.csv' '/abseq/data_prep/orig' '/abseq/output/lgunimf/full_embed' 0 &

# for hm
/usr/bin/Rscript calculate_metrics.R '/abseq/output/hm/metrics.csv' '/abseq/data_prep/orig' '/abseq/output/hm/full_embed' 0 &

# for bsc
/usr/bin/Rscript calculate_metrics.R '/abseq/output/bsc/metrics.csv' '/abseq/data_prep/orig' '/abseq/output/bsc/full_embed' 0