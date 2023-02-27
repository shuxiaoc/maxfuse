# MaxFuse: MAtching X-modality via FUzzy Smoothed Embedding


<img src="https://github.com/shuxiaoc/maxfuse/blob/main/media/ai_generated_icon.png" width="100" height="100">

## Description

MaxFuse is a Python package for integrating single-cell datasets from different modalities with no overlapping features and/or under low signal-to-noise ratio regimes. For most single-cell cross modality integration methods, the feasibility of cross-modal integration relies on the existence of highly correlated, a priori 'linked' features.  When such linked features are few or uninformative, a scenario that we call 'weak linkage', existing methods fail.  We developed MaxFuse, a cross-modal data integration method that, through iterative co-embedding, data smoothing, and cell matching, leverages all information in each modality to obtain high-quality integration. A prototypical example of weak linkage is the integration of **spatial proteomic data** with **single-cell sequencing data**. For details, please refer to the [paper](https://www.biorxiv.org/content/10.1101/2023.01.12.523851).

This work has been led by Shuxiao Chen from [Zongming Lab](http://www-stat.wharton.upenn.edu/~zongming/) @Upenn and Bokai Zhu from [Nolan lab](https://web.stanford.edu/group/nolan/) @Stanford.

<img src="https://github.com/shuxiaoc/maxfuse/blob/main/media/fig1.png" width="800" height="280">

## Installation
MaxFuse is hosted on `pypi` and can be installed via `pip`. We recommend working with a fresh virtual environment. In the following example we use conda.

```
conda create -n maxfuse python=3.8
conda activate maxfuse
python -m pip install maxfuse
```

## Vignettes

<!-- linke to sphinx ? -->

<!-- two ipynb link: -->
Example1: Protein -- RNA test run on ground-truth CITE-seq [here](https://github.com/shuxiaoc/maxfuse/blob/main/docs/citeseq_pbmc_evaluate.ipynb).

Example2: Protein -- RNA test run on tissue [here](https://github.com/shuxiaoc/maxfuse/blob/main/docs/tonsil_codex_rnaseq.ipynb).

Note in cases when integrating single cell data across **protein** and **RNA** modalities, many times the nomenclature of features are different (e.g., mRNA ```ITGAM``` could be named as ```CD11b-1``` when used as antibody). We gathered a [.csv](https://github.com/shuxiaoc/maxfuse/blob/main/docs/protein_gene_conversion.csv) file that covers many of such naming conversions and used during the ```MaxFuse``` process. Of course, this is not a complete conversion, and users should manually add in new naming conversions if they were not included in this .csv file. 

## Code archive

The analysis presented in the [manuscript](https://www.biorxiv.org/content/10.1101/2023.01.12.523851) was also deposited in this GitHub repository, under this [folder](https://github.com/shuxiaoc/maxfuse/tree/main/Archive). Note in the manuscript we used a development version of ```MaxFuse``` with slightly different grammar and can also be found there. If you require additional information on the analysis/data, please contact Zongming Ma (zongming@wharton.upenn.edu).

## License

```MaxFuse``` is under the [Academic Software License Agreement](https://github.com/shuxiaoc/maxfuse/blob/main/LICENSE), please use accordingly.
