# no condo env requirement
# calculate ari and slt f1 scores

# for mf
/usr/bin/Rscript calculate_metrics.R '/home/bkzhu/super_mario/atac_bench_nrz/greanleaf_cortical/mf/metrics.csv' '/home/bkzhu/super_mario/atac_bench_nrz/greanleaf_cortical/data/orig' '/home/bkzhu/super_mario/atac_bench_nrz/greanleaf_cortical/mf/full_embed' 0 &

# for scjoint
/usr/bin/Rscript calculate_metrics.R '/home/bkzhu/super_mario/atac_bench_nrz/greanleaf_cortical/scjoint/metrics.csv' '/home/bkzhu/super_mario/atac_bench_nrz/greanleaf_cortical/data/orig' '/home/bkzhu/super_mario/atac_bench_nrz/greanleaf_cortical/scjoint/raw_labels/full_embed' 0 &

# for maestro
/usr/bin/Rscript calculate_metrics.R '/home/bkzhu/super_mario/atac_bench_nrz/greanleaf_cortical/maestro/metrics.csv' '/home/bkzhu/super_mario/atac_bench_nrz/greanleaf_cortical/data/orig' '/home/bkzhu/super_mario/atac_bench_nrz/greanleaf_cortical/maestro/full_embed' 0 &

# for glue
/usr/bin/Rscript calculate_metrics.R '/home/bkzhu/super_mario/atac_bench_nrz/greanleaf_cortical/glue/metrics.csv' '/home/bkzhu/super_mario/atac_bench_nrz/greanleaf_cortical/data/orig' '/home/bkzhu/super_mario/atac_bench_nrz/greanleaf_cortical/glue/full_embed' 0