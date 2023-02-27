# this code is for MaxFuse matching and integration (for colon cells)

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay
import umap
import anndata as ad
import scanpy as sc

import sys
sys.path.append("../MaxFuse_devo/09302022V/") # MaxFuse devo code version, for reprodcution reasons
import match
from scipy.io import mmread
sys.executable

################# file reading

## some pre formated hubmap data, refer to script in https://github.com/shuxiaoc/maxfuse/tree/main/Archive/hubmap/code/preparation
## for details

correspondence = pd.read_csv('/home/bkzhu/MaxFuse/production/hubmap/protein_rna_name_conversionV7.csv')
rna = mmread("/home/bkzhu/MaxFuse/production/hubmap/match/data/CL_immune_rna1X.txt")
rna_stro  = mmread("/home/bkzhu/MaxFuse/production/hubmap/match/data/CL_stroma_rna1X.txt")
rna_epi  = mmread("/home/bkzhu/MaxFuse/production/hubmap/match/data/CL_epi_rna1X.txt")

protein = pd.read_csv("/home/bkzhu/MaxFuse/hubmap_phase2/codex_cleaned/CL_immune.csv")
protein_stro = pd.read_csv("/home/bkzhu/MaxFuse/hubmap_phase2/codex_cleaned/CL_stromal.csv")
protein_epi = pd.read_csv("/home/bkzhu/MaxFuse/hubmap_phase2/codex_cleaned/CL_epi_codex.csv")
protein_epi=protein_epi.rename(columns = {'Cell.Type2':'cluster.term'}) # make sure column name same

### read in pre annotated labels
protein_all1 = protein.append(protein_stro)
protein_all = protein_all1.append(protein_epi)

labels_rna = pd.read_csv('/home/bkzhu/MaxFuse/production/hubmap/match/data/CL_immune_rna_meta1X.csv')
labels_rna_stro = pd.read_csv('/home/bkzhu/MaxFuse/production/hubmap/match/data/CL_stroma_rna_meta1X.csv')
labels_rna_epi = pd.read_csv('/home/bkzhu/MaxFuse/production/hubmap/match/data/CL_epi_rna_meta1X.csv')

labels_rna_all1 = labels_rna.append(labels_rna_stro)
labels_rna_all = labels_rna_all1.append(labels_rna_epi)

rna_names = pd.read_csv('/home/bkzhu/MaxFuse/production/hubmap/match/data/SB_immune_rna_names.csv')
rna = pd.DataFrame(rna.todense())
rna.columns =rna_names['names']

rna_stro = pd.DataFrame(rna_stro.todense())
rna_stro.columns =rna_names['names']

rna_epi = pd.DataFrame(rna_epi.todense())
rna_epi.columns =rna_names['names']

rna_all1 = rna.append(rna_stro)
rna_all = rna_all1.append(rna_epi)

rna_all = rna_all.drop(drop_rna, axis=0)
protein_all = protein_all.drop(drop_pro, axis=0)

labels_protein = protein_all['cluster.term']
labels_rna_all = labels_rna_all.drop(drop_rna, axis=0)

# rename t cells
labels_protein = labels_protein.str.replace('CD8\+','T cells')
labels_protein = labels_protein.str.replace('CD4\+ T cell','T cells')
labels_rna_all['CellType2'] = labels_rna_all['CellType2'].str.replace('CD8\+ T cell','T cells')
labels_rna_all['CellType2'] = labels_rna_all['CellType2'].str.replace('CD4\+ T cell','T cells')

# clusters to smooth
labels2 = labels_protein.to_numpy()
labels1 = labels_rna_all['CellType2'].to_numpy()
labels2_s = protein_all['seurat_clusters'].to_numpy()
labels1_s = labels_rna_all['seurat_clusters'].to_numpy()

# make adata
rna_adata = ad.AnnData(
    rna_all.to_numpy(), dtype=np.float32
)
rna_adata.var_names = rna_names['names']
rna_adata.obs_names = labels_rna_all['Unnamed: 0']

protein_adata = ad.AnnData(
    protein_all[protein_all.columns[2:49]].to_numpy(), dtype=np.float32 # some markers not full
    # eg olfm4 fap cd25, coiiiv, ck7
)
protein_adata.var_names = protein_all[protein_all.columns[2:49]].columns

# name conversion
rna_protein_correspondence = []
for i in range(correspondence.shape[0]):
    curr_protein_name, curr_rna_names = correspondence.iloc[i]
    if curr_protein_name not in protein_adata.var_names:
        continue
    if curr_rna_names.find('Ignore') != -1:
        continue
    curr_rna_names = curr_rna_names.split('/')
    for r in curr_rna_names:
        if r in rna_adata.var_names:
            rna_protein_correspondence.append([r, curr_protein_name])
rna_protein_correspondence = np.array(rna_protein_correspondence)

rna_shared = rna_adata[:, rna_protein_correspondence[:, 0]].X
protein_shared = protein_adata[:, rna_protein_correspondence[:, 1]].X

## preprocessing
# remove static features
mask = ((rna_shared.std(axis=0) > 1e-5) & (protein_shared.std(axis=0) > 1e-5))
rna_shared = rna_shared[:, mask]
protein_shared = protein_shared[:, mask]

# normalize shared RNA counts
rna_shared = ad.AnnData(rna_shared)
sc.pp.normalize_total(rna_shared)
sc.pp.log1p(rna_shared)
sc.pp.scale(rna_shared)
rna_shared = rna_shared.X
# extract shared protein, no need to normalize for this data since prenormed by John already
protein_shared = protein_shared.copy()

# normalize RNA counts
sc.pp.normalize_total(rna_adata)
sc.pp.log1p(rna_adata)
sc.pp.highly_variable_genes(rna_adata)
# only retain highly variable genes
rna_adata = rna_adata[:, rna_adata.var.highly_variable].copy()
sc.pp.scale(rna_adata)

# double check if static
rna_active = rna_adata.X
protein_active = protein_adata.X
rna_active = rna_active[:, rna_active.std(axis=0) > 1e-5]
protein_active = protein_active[:, protein_active.std(axis=0) > 1e-5]




###
######################## start maxfuse process

spm = match.SuperMario(
    shared_arr1=rna_shared,
    shared_arr2=protein_shared,
    active_arr1=rna_active,
    active_arr2=protein_active,
    method='centroid_shrinkage',
    labels1=labels1_s.astype(str), # if None, then use scanpy clustering pipeline
    labels2=labels2_s.astype(str)
)

spm.split_into_batches(
    max_outward_size=5000,
    matching_ratio=6,
    metacell_size=2,
    verbose=True,
    seed = 42
)

spm.construct_graphs(
    n_neighbors1=15,
    n_neighbors2=15,
    svd_components1=30,
    svd_components2=20,
    resolution1=1,
    resolution2=1,
    randomized_svd=False, 
    svd_runs=1,
    resolution_tol=0.1,
    leiden_runs=1,
    leiden_seed=None,
    verbose=True
)

spm.find_initial_pivots(
    wt1=0.3, wt2=0.3,
    svd_components1=10, svd_components2=15,
    randomized_svd=False, svd_runs=1,
    verbose=True
)

spm.refine_pivots(
    wt1=0.3, wt2=0.3,
    svd_components1=200, svd_components2=None,
    cca_components=20,
    filter_prop=0.1,
    n_iters=1,
    randomized_svd=False, 
    svd_runs=1,
    verbose=True
)

spm.filter_bad_matches(target='pivot', filter_prop=0.3, verbose=True)

spm.propagate(
    svd_components1=100, 
    svd_components2=None, 
    randomized_svd=False, 
    svd_runs=1, 
    verbose=True
)

spm.filter_bad_matches(
    target='propagated',
    filter_prop=0.3,
    verbose=True
)

#####
################################## MaxFuse process stopped





##### saving out the results:



## 1. matching information
## first we save out the matching index (full dataset version, direction 2-->1). This means all CODEX cells (after filtering) has one match in the snRNA-seq cell dataset; in this matching index dataframe, unique(codex) will be full dataset (after filtering) while unique(snRNA-seq) might be smaller as snRNA-seq cell could be matched repeatedly.

matching = spm.get_matching(order=(2, 1), target='full_data', return_format='dict')
matching_one2one = [[], [], []]
for idx1, indices2_and_scores in matching.items():
    matching_one2one[0].append(idx1)
    matching_one2one[1].append(indices2_and_scores[0][0])
    matching_one2one[2].append(indices2_and_scores[0][1])
    
dd = {"data1":matching_one2one[1], "data2":matching_one2one[0], "score": matching_one2one[2]}
df = pd.DataFrame.from_dict(dd)
df.to_csv("/home/bkzhu/MaxFuse/hubmap_phase2/match/output/CL_0719/full21.csv") # this is dataframe telling for each codex cell, which snrnaseq cell is considered the best match




## 2. integrated embedding
## for the embedding, here we want the cca embedding for both the codex cells (all after filtering) and all snrnaseq that got matched. Since for visualization, it will not be a good idea to put repeated cells (eg for snrnaseq, there are cells got matched  repeatedly), we first get the unique matched snrnaseq cells, and then plot them.

matching = spm.get_matching(order=(2, 1), target='full_data', return_format='dict')
matching_one2one = [[], [], []]

for idx1, indices2_and_scores in matching.items():
    matching_one2one[0].append(idx1)
    matching_one2one[1].append(indices2_and_scores[0][0])
    matching_one2one[2].append(indices2_and_scores[0][1])

dd = {"data1":np.unique(matching_one2one[1])}
df = pd.DataFrame.from_dict(dd)
df.to_csv("/home/bkzhu/MaxFuse/hubmap_phase2/match/output/CL_0719/full21_unique_data1.csv")


## then we can compute a umap based on the cca scores produced by maxfuse and save them out for plotting
arr1, arr2 = spm.get_embedding(target=[spm.active_arr1[np.unique(matching_one2one[1]),:],
                                      spm.active_arr2[matching_one2one[0],:]], refit=False, cca_components=None)
reducer = umap.UMAP()
embedding = reducer.fit_transform(np.vstack([arr1, arr2]))
ebddf = pd.DataFrame(embedding, columns = ['umap1', 'umap2'])
ebddf.to_csv("/home/bkzhu/MaxFuse/hubmap_phase2/match/output/CL_0719/full21_unique_embedding.csv")


