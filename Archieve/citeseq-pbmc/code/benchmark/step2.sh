# script to produce ARI f1 score and SLT f1 score
# only done on the full panel, not done one the dropping verions

# b1-5 for mf
/usr/bin/Rscript calculate_metrics.R '//bench_test3/output/b1/mf/metrics.csv' '/bench_test3/input/b1/orig' '/bench_test3/output/b1/mf/full_embed' 0 &
/usr/bin/Rscript calculate_metrics.R '/bench_test3/output/b2/mf/metrics.csv' '/bench_test3/input/b2/orig' '/bench_test3/output/b2/mf/full_embed' 0 &
/usr/bin/Rscript calculate_metrics.R '/bench_test3/output/b3/mf/metrics.csv' '/bench_test3/input/b3/orig' '/bench_test3/output/b3/mf/full_embed' 0 &
/usr/bin/Rscript calculate_metrics.R '/bench_test3/output/b4/mf/metrics.csv' '/bench_test3/input/b4/orig' '/bench_test3/output/b4/mf/full_embed' 0 &
/usr/bin/Rscript calculate_metrics.R '/bench_test3/output/b5/mf/metrics.csv' '/bench_test3/input/b5/orig' '/bench_test3/output/b5/mf/full_embed' 0 &
wait
# b1-5 for sr
/usr/bin/Rscript calculate_metrics.R '/bench_test3/output/b1/sr/metrics.csv' '/bench_test3/input/b1/orig' '/bench_test3/output/b1/sr/full_embed' 0 &
/usr/bin/Rscript calculate_metrics.R '/bench_test3/output/b2/sr/metrics.csv' '/bench_test3/input/b2/orig' '/bench_test3/output/b2/sr/full_embed' 0 &
/usr/bin/Rscript calculate_metrics.R '/bench_test3/output/b3/sr/metrics.csv' '/bench_test3/input/b3/orig' '/bench_test3/output/b3/sr/full_embed' 0 &
/usr/bin/Rscript calculate_metrics.R '/bench_test3/output/b4/sr/metrics.csv' '/bench_test3/input/b4/orig' '/bench_test3/output/b4/sr/full_embed' 0 &
/usr/bin/Rscript calculate_metrics.R '/bench_test3/output/b5/sr/metrics.csv' '/bench_test3/input/b5/orig' '/bench_test3/output/b5/sr/full_embed' 0 &
wait
# b1-5 for lg
/usr/bin/Rscript calculate_metrics.R '/bench_test3/output/b1/lgunimf/metrics.csv' '/bench_test3/input/b1/orig' '/bench_test3/output/b1/lgunimf/full_embed' 0 &
/usr/bin/Rscript calculate_metrics.R '/bench_test3/output/b2/lgunimf/metrics.csv' '/bench_test3/input/b2/orig' '/bench_test3/output/b2/lgunimf/full_embed' 0 &
/usr/bin/Rscript calculate_metrics.R '/bench_test3/output/b3/lgunimf/metrics.csv' '/bench_test3/input/b3/orig' '/bench_test3/output/b3/lgunimf/full_embed' 0 &
/usr/bin/Rscript calculate_metrics.R '/bench_test3/output/b4/lgunimf/metrics.csv' '/bench_test3/input/b4/orig' '/bench_test3/output/b4/lgunimf/full_embed' 0 &
/usr/bin/Rscript calculate_metrics.R '/bench_test3/output/b5/lgunimf/metrics.csv' '/bench_test3/input/b5/orig' '/bench_test3/output/b5/lgunimf/full_embed' 0 &
wait
# b1-5 for hm
/usr/bin/Rscript calculate_metrics.R '/bench_test3/output/b1/hm/metrics.csv' '/bench_test3/input/b1/orig' '/bench_test3/output/b1/hm/full_embed' 0 &
/usr/bin/Rscript calculate_metrics.R '/bench_test3/output/b2/hm/metrics.csv' '/bench_test3/input/b2/orig' '/bench_test3/output/b2/hm/full_embed' 0 &
/usr/bin/Rscript calculate_metrics.R '/bench_test3/output/b3/hm/metrics.csv' '/bench_test3/input/b3/orig' '/bench_test3/output/b3/hm/full_embed' 0 &
/usr/bin/Rscript calculate_metrics.R '/bench_test3/output/b4/hm/metrics.csv' '/bench_test3/input/b4/orig' '/bench_test3/output/b4/hm/full_embed' 0 &
/usr/bin/Rscript calculate_metrics.R '/bench_test3/output/b5/hm/metrics.csv' '/bench_test3/input/b5/orig' '/bench_test3/output/b5/hm/full_embed' 0
wait
# b1-5 for bsc
/usr/bin/Rscript calculate_metrics.R '/bench_test3/output/b1/bsc/metrics.csv' '/bench_test3/input/b1/orig' '/bench_test3/output/b1/bsc/full_embed' 0 &
/usr/bin/Rscript calculate_metrics.R '/bench_test3/output/b2/bsc/metrics.csv' '/bench_test3/input/b2/orig' '/bench_test3/output/b2/bsc/full_embed' 0 &
/usr/bin/Rscript calculate_metrics.R '/bench_test3/output/b3/bsc/metrics.csv' '/bench_test3/input/b3/orig' '/bench_test3/output/b3/bsc/full_embed' 0 &
/usr/bin/Rscript calculate_metrics.R '/bench_test3/output/b4/bsc/metrics.csv' '/bench_test3/input/b4/orig' '/bench_test3/output/b4/bsc/full_embed' 0 &
/usr/bin/Rscript calculate_metrics.R '/bench_test3/output/b5/bsc/metrics.csv' '/bench_test3/input/b5/orig' '/bench_test3/output/b5/bsc/full_embed' 0