# full data set run for related spatial analysis
# for maxfuse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay

import anndata as ad
import scanpy as sc

import sys
sys.path.append("../../MaxFuse_devo/09302022V/")
import match
import metrics
from scipy.io import mmread
import os

# read in data
root_dir = '/tonsil_v2/'
rna = mmread(root_dir + '/RNA/tonsil_rna_0510.txt')
meta_rna = pd.read_csv(root_dir + '/RNA/tonsil_rna_0510_meta.csv')
rna_names = pd.read_csv(root_dir + '/RNA/tonsil_rna_0510_names.csv')
#
protein = pd.read_csv(root_dir + '/Codex/FCS_output_DeepCell_extOnly/formatch_clusters_x28_y715V2.csv')
protein.rename(columns={'collagen.IV':'collagen IV'}, inplace=True)
protein.rename(columns={'HLA.DR':'HLA-DR'}, inplace=True)
annotation_rna = meta_rna['cluster.info'].to_numpy()
annotation_pro = protein['cluster.term'].to_numpy()

# convert to adata
rna_adata = ad.AnnData(
    rna.tocsr(), dtype=np.float32
)
rna_adata.var_names = rna_names['names']
rna_adata.obs_names = meta_rna['Unnamed: 0']

protein = protein.drop(['X','Unnamed: 0','cellLabelInImage','cellSize','HOECHST1',
                       'PointNum','cluster.term','seurat_res1.0'], axis=1)
protein_adata = ad.AnnData(
    protein.to_numpy(), dtype=np.float32
)
protein_adata.var_names = protein.columns

# get rna--protein correspondence
correspondence = pd.read_csv('/tonsil_v2/match/protein_rna_name_conversionV11.csv')
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

rna_shared = rna_adata[:, rna_protein_correspondence[:, 0]].X.todense()
protein_shared = protein_adata[:, rna_protein_correspondence[:, 1]].X.copy()

# remove static features for shared
mask = ((rna_shared.std(axis=0) > 0.5) & (protein_shared.std(axis=0) > 0.1)).A1
rna_shared = rna_shared[:, mask]
protein_shared = protein_shared[:, mask]
print([rna_shared.shape,protein_shared.shape])

# normalize shared RNA counts
rna_shared = ad.AnnData(rna_shared)
sc.pp.normalize_total(rna_shared)
sc.pp.log1p(rna_shared)
sc.pp.scale(rna_shared)
rna_shared = rna_shared.X
# normalize shared protein counts
protein_shared = ad.AnnData(protein_shared)
sc.pp.normalize_total(protein_shared)
sc.pp.log1p(protein_shared)
sc.pp.scale(protein_shared)
protein_shared = protein_shared.X

##################################

# normalize RNA counts
sc.pp.normalize_total(rna_adata)
sc.pp.log1p(rna_adata)
sc.pp.highly_variable_genes(rna_adata, n_top_genes=5000)
rna_adata = rna_adata[:, rna_adata.var.highly_variable].copy()
sc.pp.scale(rna_adata)
rna_active = rna_adata.X
# normalize protein counts
sc.pp.normalize_total(protein_adata)
sc.pp.log1p(protein_adata)
sc.pp.highly_variable_genes(protein_adata)
sc.pp.scale(protein_adata)
protein_active = protein_adata.X

mf = match.MaxFuse(
    shared_arr1=rna_shared,
    shared_arr2=protein_shared,
    active_arr1=rna_active,
    active_arr2=protein_active,
    method='centroid_shrinkage',
    labels1=None,
    labels2=None
)

mf.split_into_batches(
    max_outward_size=10000,
    matching_ratio=5,
    metacell_size=2,
    method = 'binning',
    seed=None,
    verbose=True
)

mf.construct_graphs(
    n_neighbors1=15,
    n_neighbors2=15,
    svd_components1=50,
    svd_components2=15,
    resolution1=2,
    resolution2=2,
    randomized_svd=False, 
    svd_runs=1,
    resolution_tol=0.1,
    leiden_runs=1,
    leiden_seed=None,
    verbose=True
)

mf.find_initial_pivots(
    wt1=0.3, wt2=0.3,
    svd_components1=25, svd_components2=20,
    randomized_svd=False, svd_runs=1,
    verbose=True
)

mf.refine_pivots(
    wt1=0.3, wt2=0.3,
    svd_components1=200, svd_components2=None,
    cca_components=45,
    filter_prop=0.,
    n_iters=1,
    randomized_svd=False, 
    svd_runs=1,
    verbose=True
)

mf.filter_bad_matches(target='pivot', filter_prop=0.5, verbose=True)

mf.propagate(
    wt1=0.7, wt2=0.7,
    svd_components1=50, 
    svd_components2=None, 
    randomized_svd=False, 
    svd_runs=1, 
    verbose=True
)

mf.filter_bad_matches(
    target='propagated',
    filter_prop=0.3,
    verbose=True
)

full_matching = mf.get_matching(order=(2, 1), target='full_data')
full_matchingV2 = match_utils.address_matching_redundancy(full_matching, order=(2, 1))

full_matchingV2 = match_utils.address_matching_redundancy(full_matching, order=(2, 1))
full = pd.DataFrame(list(zip(full_matchingV2[0],full_matchingV2[1],full_matchingV2[2])), columns = ["idx1","idx2","score"])
full.to_csv("/tonsil_v2/match/match_output/1205_tonsil_fullID.csv", index=False)

reducer = umap.UMAP()
embedding = reducer.fit_transform(np.vstack([arr1[:,0:15], arr2[:,0:15]]))
ebddf = pd.DataFrame(embedding, columns = ['umap1', 'umap2'])
ebddf.to_csv("/tonsil_v2/match/match_output/full/mf/allrna_vanilla_umap_cca15.csv")
