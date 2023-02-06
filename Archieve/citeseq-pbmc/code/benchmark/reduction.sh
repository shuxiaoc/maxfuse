## run all methods, result used for umap viz and confuse matrix plotting
python cite_mf_reduction.py &
/usr/bin/Rscript seurat_cite_reduc.R &
/usr/bin/Rscript liger_cite_reduction.R &
/usr/bin/Rscript harm_cite_reduc.R &
/usr/bin/Rscript bsc_cite_reduction.R
