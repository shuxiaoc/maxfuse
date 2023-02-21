# no condo env requirement
# calculate ari and slt f1 scores

# for mf
/usr/bin/Rscript calculate_metrics.R '/home/bkzhu/super_mario/atac_bench_nrz/retina/mf/metrics.csv' '/home/bkzhu/super_mario/atac_bench_nrz/retina/data/orig' '/home/bkzhu/super_mario/atac_bench_nrz/retina/mf/full_embed' 0 &

# for scjoint
/usr/bin/Rscript calculate_metrics.R '/home/bkzhu/super_mario/atac_bench_nrz/retina/scj/metrics.csv' '/home/bkzhu/super_mario/atac_bench_nrz/retina/data/orig' '/home/bkzhu/super_mario/atac_bench_nrz/retina/scj/full_embed' 0 &

# for maestro
/usr/bin/Rscript calculate_metrics.R '/home/bkzhu/super_mario/atac_bench_nrz/retina/ms/metrics.csv' '/home/bkzhu/super_mario/atac_bench_nrz/retina/data/orig' '/home/bkzhu/super_mario/atac_bench_nrz/retina/ms/full_embed' 0 &

# for glue
/usr/bin/Rscript calculate_metrics.R '/home/bkzhu/super_mario/atac_bench_nrz/retina/glue/metrics.csv' '/home/bkzhu/super_mario/atac_bench_nrz/retina/data/orig' '/home/bkzhu/super_mario/atac_bench_nrz/retina/glue/full_embed' 0