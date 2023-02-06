## script to run all methods at the same time:
## step one of benchmarking, dropping antibody panel version
python maxfuse_cite-drop.py &
/usr/bin/Rscript seurat_cite_drop.R &
/usr/bin/Rscript liger_cite_drop.R &
/usr/bin/Rscript harm_cite_drop.R &
/usr/bin/Rscript bsc_cite_drop.R
