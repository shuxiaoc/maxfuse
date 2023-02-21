# no condo env requirement
# calculate ari and slt f1 scores

# for mf
/usr/bin/Rscript calculate_metrics.R '/home/bkzhu/super_mario/atac_bench_nrz/10x_e18/mf/metrics.csv' '/home/bkzhu/super_mario/atac_bench_nrz/10x_e18/data/orig' '/home/bkzhu/super_mario/atac_bench_nrz/10x_e18/mf/full_embed' 0 &

# for scjoint
/usr/bin/Rscript calculate_metrics.R '/home/bkzhu/super_mario/atac_bench_nrz/10x_e18/scjoint/metrics.csv' '/home/bkzhu/super_mario/atac_bench_nrz/10x_e18/data/orig' '/home/bkzhu/super_mario/atac_bench_nrz/10x_e18/scjoint/raw_labels/full_embed' 0 &

# for maestro
/usr/bin/Rscript calculate_metrics.R '/home/bkzhu/super_mario/atac_bench_nrz/10x_e18/maestro/metrics.csv' '/home/bkzhu/super_mario/atac_bench_nrz/10x_e18/data/orig' '/home/bkzhu/super_mario/atac_bench_nrz/10x_e18/maestro/full_embed' 0 &

# for glue
/usr/bin/Rscript calculate_metrics.R '/home/bkzhu/super_mario/atac_bench_nrz/10x_e18/glue/metrics.csv' '/home/bkzhu/super_mario/atac_bench_nrz/10x_e18/data/orig' '/home/bkzhu/super_mario/atac_bench_nrz/10x_e18/glue/full_embed' 0