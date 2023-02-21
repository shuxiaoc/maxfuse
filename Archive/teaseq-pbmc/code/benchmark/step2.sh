# no condo env requirement
# for mf
/usr/bin/Rscript calculate_metrics.R '/ICICLE/output/mf/metrics.csv' '/ICICLE/data/orig' '/ICICLE/output/mf/full_embed' 0 &

# for sr
/usr/bin/Rscript calculate_metrics.R '/ICICLE/output/sr/metrics.csv' '/ICICLE/data/orig' '/ICICLE/output/sr/full_embed' 0 &

# for lg
/usr/bin/Rscript calculate_metrics.R '/ICICLE/output/lg/metrics.csv' '/ICICLE/data/orig_lg' '/ICICLE/output/lg/full_embed' 0 &

# for hm
/usr/bin/Rscript calculate_metrics.R '/ICICLE/output/hm/metrics.csv' '/ICICLE/data/orig' '/ICICLE/output/hm/full_embed' 0 &

# for bsc
/usr/bin/Rscript calculate_metrics.R '/ICICLE/output/bsc/metrics.csv' '/ICICLE/data/orig' '//ICICLE/output/bsc/full_embed' 0