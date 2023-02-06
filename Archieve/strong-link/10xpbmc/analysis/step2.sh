# no condo env requirement
# calculate ari and slt f1 scores

# for mf
/usr/bin/Rscript calculate_metrics.R '/home/bkzhu/super_mario/atac_bench_nrz/10xpbmc/mf/metrics.csv' '/home/bkzhu/super_mario/atac_bench_nrz/10xpbmc/data/orig' '/home/bkzhu/super_mario/atac_bench_nrz/10xpbmc/mf/full_embed' 0 &

# for scjoint
/usr/bin/Rscript calculate_metrics.R '/home/bkzhu/super_mario/atac_bench_nrz/10xpbmc/scJoint/metrics.csv' '/home/bkzhu/super_mario/atac_bench_nrz/10xpbmc/data/orig' '/home/bkzhu/super_mario/atac_bench_nrz/10xpbmc/scJoint/raw_labels/full_embed' 0 &

# for maestro
/usr/bin/Rscript calculate_metrics.R '/home/bkzhu/super_mario/atac_bench_nrz/10xpbmc/maestro/metrics.csv' '/home/bkzhu/super_mario/atac_bench_nrz/10xpbmc/data/orig' '/home/bkzhu/super_mario/atac_bench_nrz/10xpbmc/maestro/full_embed' 0 &

# for glue
/usr/bin/Rscript calculate_metrics.R '/home/bkzhu/super_mario/atac_bench_nrz/10xpbmc/glue/metrics.csv' '/home/bkzhu/super_mario/atac_bench_nrz/10xpbmc/data/orig' '/home/bkzhu/super_mario/atac_bench_nrz/10xpbmc/glue/full_embed' 0