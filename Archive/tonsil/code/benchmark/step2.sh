# no condo env requirement
# b1-5 for mf
/usr/bin/Rscript calculate_metrics.R '/tonsil_v2/match/bench_out/b1/mf/metrics.csv' '/tonsil_v2/match/bench_input/b1/orig' '/tonsil_v2/match/bench_out/b1/mf/full_embed' 0 &
/usr/bin/Rscript calculate_metrics.R '/tonsil_v2/match/bench_out/b2/mf/metrics.csv' '/tonsil_v2/match/bench_input/b2/orig' '/tonsil_v2/match/bench_out/b2/mf/full_embed' 0 &
/usr/bin/Rscript calculate_metrics.R '/tonsil_v2/match/bench_out/b3/mf/metrics.csv' '/tonsil_v2/match/bench_input/b3/orig' '/tonsil_v2/match/bench_out/b3/mf/full_embed' 0 &
/usr/bin/Rscript calculate_metrics.R '/tonsil_v2/match/bench_out/b4/mf/metrics.csv' '/tonsil_v2/match/bench_input/b4/orig' '/tonsil_v2/match/bench_out/b4/mf/full_embed' 0 &
/usr/bin/Rscript calculate_metrics.R '/tonsil_v2/match/bench_out/b5/mf/metrics.csv' '/tonsil_v2/match/bench_input/b5/orig' '/tonsil_v2/match/bench_out/b5/mf/full_embed' 0 &
wait
# b1-5 for sr
/usr/bin/Rscript calculate_metrics.R '/tonsil_v2/match/bench_out/b1/sr/metrics.csv' '/tonsil_v2/match/bench_input/b1/orig' '/tonsil_v2/match/bench_out/b1/sr/full_embed' 0 &
/usr/bin/Rscript calculate_metrics.R '/tonsil_v2/match/bench_out/b2/sr/metrics.csv' '/tonsil_v2/match/bench_input/b2/orig' '/tonsil_v2/match/bench_out/b2/sr/full_embed' 0 &
/usr/bin/Rscript calculate_metrics.R '/tonsil_v2/match/bench_out/b3/sr/metrics.csv' '/tonsil_v2/match/bench_input/b3/orig' '/tonsil_v2/match/bench_out/b3/sr/full_embed' 0 &
/usr/bin/Rscript calculate_metrics.R '/tonsil_v2/match/bench_out/b4/sr/metrics.csv' '/tonsil_v2/match/bench_input/b4/orig' '/tonsil_v2/match/bench_out/b4/sr/full_embed' 0 &
/usr/bin/Rscript calculate_metrics.R '/tonsil_v2/match/bench_out/b5/sr/metrics.csv' '/tonsil_v2/match/bench_input/b5/orig' '/tonsil_v2/match/bench_out/b5/sr/full_embed' 0 &
wait
# b1-5 for lg
/usr/bin/Rscript calculate_metrics.R '/tonsil_v2/match/bench_out/b1/lg/metrics.csv' '/tonsil_v2/match/bench_out/b1/lg/orig' '/tonsil_v2/match/bench_out/b1/lg/full_embed' 0 &
/usr/bin/Rscript calculate_metrics.R '/tonsil_v2/match/bench_out/b2/lg/metrics.csv' '/tonsil_v2/match/bench_out/b2/lg/orig' '/tonsil_v2/match/bench_out/b2/lg/full_embed' 0 &
/usr/bin/Rscript calculate_metrics.R '/tonsil_v2/match/bench_out/b3/lg/metrics.csv' '/tonsil_v2/match/bench_out/b3/lg/orig' '/tonsil_v2/match/bench_out/b3/lg/full_embed' 0 &
/usr/bin/Rscript calculate_metrics.R '/tonsil_v2/match/bench_out/b4/lg/metrics.csv' '/tonsil_v2/match/bench_out/b4/lg/orig' '/tonsil_v2/match/bench_out/b4/lg/full_embed' 0 &
/usr/bin/Rscript calculate_metrics.R '/tonsil_v2/match/bench_out/b5/lg/metrics.csv' '/tonsil_v2/match/bench_out/b5/lg/orig' '/tonsil_v2/match/bench_out/b5/lg/full_embed' 0 &
wait
# b1-5 for hm
/usr/bin/Rscript calculate_metrics.R '/tonsil_v2/match/bench_out/b1/hm/metrics.csv' '/tonsil_v2/match/bench_input/b1/orig' '/tonsil_v2/match/bench_out/b1/hm/full_embed' 0 &
/usr/bin/Rscript calculate_metrics.R '/tonsil_v2/match/bench_out/b2/hm/metrics.csv' '/tonsil_v2/match/bench_input/b2/orig' '/tonsil_v2/match/bench_out/b2/hm/full_embed' 0 &
/usr/bin/Rscript calculate_metrics.R '/tonsil_v2/match/bench_out/b3/hm/metrics.csv' '/tonsil_v2/match/bench_input/b3/orig' '/tonsil_v2/match/bench_out/b3/hm/full_embed' 0 &
/usr/bin/Rscript calculate_metrics.R '/tonsil_v2/match/bench_out/b4/hm/metrics.csv' '/tonsil_v2/match/bench_input/b4/orig' '/tonsil_v2/match/bench_out/b4/hm/full_embed' 0 &
/usr/bin/Rscript calculate_metrics.R '/tonsil_v2/match/bench_out/b5/hm/metrics.csv' '/tonsil_v2/match/bench_input/b5/orig' '/tonsil_v2/match/bench_out/b5/hm/full_embed' 0 &
wait
# b1-5 for bsc
/usr/bin/Rscript calculate_metrics.R '/tonsil_v2/match/bench_out/b1/bsc/metrics.csv' '/tonsil_v2/match/bench_input/b1/orig' '/tonsil_v2/match/bench_out/b1/bsc/full_embed' 0 &
/usr/bin/Rscript calculate_metrics.R '/tonsil_v2/match/bench_out/b2/bsc/metrics.csv' '/tonsil_v2/match/bench_input/b2/orig' '/tonsil_v2/match/bench_out/b2/bsc/full_embed' 0 &
/usr/bin/Rscript calculate_metrics.R '/tonsil_v2/match/bench_out/b3/bsc/metrics.csv' '/tonsil_v2/match/bench_input/b3/orig' '/tonsil_v2/match/bench_out/b3/bsc/full_embed' 0 &
/usr/bin/Rscript calculate_metrics.R '/tonsil_v2/match/bench_out/b4/bsc/metrics.csv' '/tonsil_v2/match/bench_input/b4/orig' '/tonsil_v2/match/bench_out/b4/bsc/full_embed' 0 &
/usr/bin/Rscript calculate_metrics.R '/tonsil_v2/match/bench_out/b5/bsc/metrics.csv' '/tonsil_v2/match/bench_input/b5/orig' '/tonsil_v2/match/bench_out/b5/bsc/full_embed' 0
# b1-5 for lg
/usr/bin/Rscript calculate_metrics.R '/tonsil_v2/match/bench_out/b1/lg/metrics.csv' '/tonsil_v2/match/bench_input/b1/orig' '/tonsil_v2/match/bench_out/b1/lgunimf/full_embed' 0 &
/usr/bin/Rscript calculate_metrics.R '/tonsil_v2/match/bench_out/b2/lg/metrics.csv' '/tonsil_v2/match/bench_input/b2/orig' '/tonsil_v2/match/bench_out/b2/lgunimf/full_embed' 0 &
/usr/bin/Rscript calculate_metrics.R '/tonsil_v2/match/bench_out/b3/lg/metrics.csv' '/tonsil_v2/match/bench_input/b3/orig' '/tonsil_v2/match/bench_out/b3/lgunimf/full_embed' 0 &
/usr/bin/Rscript calculate_metrics.R '/tonsil_v2/match/bench_out/b4/lg/metrics.csv' '/tonsil_v2/match/bench_input/b4/orig' '/tonsil_v2/match/bench_out/b4/lgunimf/full_embed' 0 &
/usr/bin/Rscript calculate_metrics.R '/tonsil_v2/match/bench_out/b5/lg/metrics.csv' '/tonsil_v2/match/bench_input/b5/orig' '/tonsil_v2/match/bench_out/b5/lgunimf/full_embed' 0 
