# quick code to calc slt ari f1 scores for all methods
# for mf
/usr/bin/Rscript calculate_metrics.R '/bench_test4/output/mf/metrics.csv' '/bench_test4/input/orig' '/bench_test4/output/mf/full_embed' 0 &

# for sr
/usr/bin/Rscript calculate_metrics.R '/bench_test4/output/sr/metrics.csv' '/bench_test4/input/orig' '/bench_test4/output/sr/full_embed' 0 &

# for lg
/usr/bin/Rscript calculate_metrics.R '/bench_test4/output/lgunimf/metrics.csv' '/bench_test4/input/orig' '/bench_test4/output/lgunimf/full_embed' 0 &

# for hm
/usr/bin/Rscript calculate_metrics.R '/bench_test4/output/hm/metrics.csv' '/bench_test4/input/orig' '/bench_test4/output/hm/full_embed' 0

# for bsc
/usr/bin/Rscript calculate_metrics.R '/bench_test4/output/bsc/metrics.csv' '/bench_test4/input/orig' '/bench_test4/output/bsc/full_embed' 0