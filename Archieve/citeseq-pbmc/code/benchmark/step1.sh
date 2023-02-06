## script to run all methods at the same time:
## step one of benchmarking, full antibody panel version
python maxfuse_cite.py &
/usr/bin/Rscript seurat_cite.R &
/usr/bin/Rscript liger_cite.R &
/usr/bin/Rscript harm_cite.R &
/usr/bin/Rscript bsc_cite.R