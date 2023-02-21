# code to calc slt and ari f1

# for mf
/usr/bin/Rscript calculate_metrics.R '/asap/output/mf/metrics.csv' '/asap/data/orig' '/asap/output/mf/full_embed' 0 &

# for sr
/usr/bin/Rscript calculate_metrics.R '/asap/output/sr/metrics.csv' '/asap/data/orig' '/asap/output/sr/full_embed' 0 &

# for lg
/usr/bin/Rscript calculate_metrics.R '/asap/output/lg/metrics.csv' '/asap/data/orig' '/asap/output/lg/full_embed' 0 &

# for hm
/usr/bin/Rscript calculate_metrics.R '/asap/output/hm/metrics.csv' '/asap/data/orig' '/asap/output/hm/full_embed' 0 &

# for bsc
/usr/bin/Rscript calculate_metrics.R '/asap/output/bsc/metrics.csv' '/asap/data/orig' '/asap/output/bsc/full_embed' 0