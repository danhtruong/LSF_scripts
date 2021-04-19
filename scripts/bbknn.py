import numpy as np
import scanpy as sc
import pandas as pd
import anndata
#import bbknn

results_file = 'outputs/mtl_bbknn.h5ad'  # the file that will store the analysis results

pca = pd.read_csv('inputs/mtl_pca.csv', index_col=0)
batch = pd.read_csv('inputs/mtl_batch.csv', index_col=0)
batch = batch.set_index(pca.index)
print(pca.head())
print(batch.head())
adata = anndata.AnnData(X=pca, obs=batch)
sc.tl.pca(adata)
print(adata)
adata.obsm['X_pca'] = pca
adata = sc.external.pp.bbknn(adata, batch_key='x', metric='euclidean')
adata.write(results_file)
