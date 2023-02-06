import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay

import anndata as ad
import scanpy as sc

import sys
sys.path.append("../../MaxFuse_devo/09302022V/")
import match
from scipy.io import mmread
import os

# input chunk
# read in files
out_root = "/tonsil_v2/match/bench_out/"
in_root = "/tonsil_v2/match/bench_input/"
batch = 5
out_idx = 15

for i in range(batch):
    batch_name = "b" + str(i+1) + "/"
    out_dir =out_root + batch_name + "mf/"
    in_dir = in_root + batch_name
    isExist = os.path.exists(out_dir)
    if not isExist:
        os.makedirs(out_dir)
    # read in files
    rna = mmread(in_dir + "rna.txt")
    protein = pd.read_csv(in_dir + "pro.csv")
    meta_rna = pd.read_csv(in_dir + "meta_rna.csv")
    meta_pro = pd.read_csv(in_dir + "meta_pro.csv")
    
    protein = protein.rename(columns={'HLA.DR': 'HLA-DR','collagen.IV': 'collagen IV'}) # for python name auto changed
    rna_names = pd.read_csv("/tonsil_v2/RNA/tonsil_rna_0510_names.csv")
    #
    rna_adata = ad.AnnData(
        rna.tocsr(), dtype=np.float32
    )
    rna_adata.var_names = rna_names['names']
    rna_adata.obs_names = meta_rna['X']
    
    #protein = protein.drop(['Unnamed: 0','cellLabelInImage','cellSize','HOECHST1',
    #                   'PointNum','cluster.term','seurat_res1.0'], axis=1)
    protein_adata = ad.AnnData(
        protein.to_numpy(), dtype=np.float32
    )
    protein_adata.var_names = protein.columns

    correspondence = pd.read_csv('/tonsil_v2/match/protein_rna_name_conversionV8.csv')
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
    protein_shared = protein_adata[:, rna_protein_correspondence[:, 1]].X

    # remove static features for shared
    mask = ((rna_shared.std(axis=0) > 0.5) & (protein_shared.std(axis=0) > 0.1)).A1
    rna_shared = rna_shared[:, mask]
    protein_shared = protein_shared[:, mask]

    # normalize shared RNA counts
    rna_shared = ad.AnnData(rna_shared)
    sc.pp.normalize_total(rna_shared)
    sc.pp.log1p(rna_shared)
    sc.pp.scale(rna_shared)
    rna_shared = rna_shared.X
    ## protein shared
    protein_shared = ad.AnnData(protein_shared) # value of codex is pre-normalized already
    sc.pp.normalize_total(protein_shared)
    sc.pp.log1p(protein_shared)
    sc.pp.scale(protein_shared)
    protein_shared = protein_shared.X

    # only retain highly variable genes
    # normalize RNA counts
    sc.pp.normalize_total(rna_adata)
    sc.pp.log1p(rna_adata)
    sc.pp.highly_variable_genes(rna_adata, n_top_genes = 5000)
    rna_adata = rna_adata[:, rna_adata.var.highly_variable].copy()
    sc.pp.scale(rna_adata)
    #
    sc.pp.normalize_total(protein_adata)
    sc.pp.log1p(protein_adata)
    sc.pp.highly_variable_genes(protein_adata)
    sc.pp.scale(protein_adata)
    # check active all features
    rna_active = rna_adata.X
    protein_active = protein_adata.X
    rna_active = rna_active[:, rna_active.std(axis=0) > 0.01]
    protein_active = protein_active[:, protein_active.std(axis=0) > 0.01] ## temporary values
    #
    spm = match.MaxFuse(
        shared_arr1=rna_shared,
        shared_arr2=protein_shared,
        active_arr1=rna_active,
        active_arr2=protein_active,
        method='centroid_shrinkage',
        labels1=None, # if None, then use scanpy clustering pipeline
        labels2=None
    )
    #
    spm.split_into_batches(
        max_outward_size=5000,
        matching_ratio=5,
        metacell_size=3,
        verbose=True,
        seed = 42
    )
    #
    spm.construct_graphs(
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
    #
    spm.find_initial_pivots(
        wt1=0.3, wt2=0.3,
        svd_components1=25, svd_components2=20,
        randomized_svd=False, svd_runs=1,
        verbose=True
    )
    #
    spm.refine_pivots(
        wt1=0.3, wt2=0.3,
        svd_components1=200, svd_components2=None,
        cca_components=45,
        filter_prop=0.,
        n_iters=1,
        randomized_svd=False, 
        svd_runs=1,
        verbose=True
    )
    #
    spm.filter_bad_matches(target='pivot', filter_prop=0.5, verbose=True)
    #
    spm.propagate(
        svd_components1=50, 
        svd_components2=None, 
        randomized_svd=False, 
        svd_runs=1, 
        verbose=True
    )
    #
    spm.filter_bad_matches(
        target='propagated',
        filter_prop=0.,
        verbose=True
    )
    matching = spm.get_matching(order=(2, 1), target='full_data')
    # save metrics file
    data = {'method': ['sxc']}  
    pd.DataFrame(data).to_csv(out_dir + "/metrics.csv",index=False)
    #
    full = pd.DataFrame(list(zip(matching[0],matching[1],matching[2])), columns = ["idx1","idx2","score"])
    full.to_csv(out_dir + "/full_idx.csv", index=False)
    #
    arr1_cca, arr2_cca = spm.get_embedding(
        active_arr1 = spm.active_arr1,
        active_arr2 = spm.active_arr2,
        refit=False,
        matching=None,
        order=None,
        cca_components=15,
        cca_max_iter=None
    )
    #
    arr1_df = pd.DataFrame(arr1_cca).iloc[:,0:out_idx]
    arr2_df = pd.DataFrame(arr2_cca).iloc[:,0:out_idx]
    arr1_df.to_csv(out_dir + "/full_embed_x0.csv",index=False)
    arr2_df.to_csv(out_dir + "/full_embed_y0.csv", index=False)

    