import os
import gc
import ot
import pickle
import anndata
import scanpy as sc
import pandas as pd
import numpy as np
from scipy import sparse
from scipy.stats import spearmanr, pearsonr
from scipy.spatial import distance_matrix
import matplotlib.pyplot as plt

import commot as ct

indir = "/public/workspace/stu21230110/AMP/01data/"

samples = ["C1","C2","C3","P1","P2","P3"]

outdir = "/public/workspace/stu21230110/AMP/11commot/"

for sample_name in samples:
    outdir = "/public/workspace/stu21230110/AMP/11commot/" + sample_name + "/"
    adata = sc.read(indir + sample_name + ".spatial.h5ad")
    adata.var_names_make_unique()
    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata)
    adata_dis500 = adata.copy()
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata = adata[:, adata.var.highly_variable]
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=0.4)
    sc.pl.spatial(adata, color='leiden')
    plt.savefig(outdir + sample_name + ".umap.pdf",bbox_inches = "tight")
    df_cellchat = ct.pp.ligand_receptor_database(species="human", 
    #                                             signaling_type='Secreted Signaling', 
                                                 database='CellChat')
    df_cellchat_filtered = ct.pp.filter_lr_database(df_cellchat, adata_dis500, min_cell_pct=0.05)
    ct.tl.spatial_communication(adata_dis500,
                                database_name='cellchat', df_ligrec=df_cellchat_filtered, dis_thr=500, heteromeric=True, pathway_sum=True)
    adata_dis500.write(outdir + sample_name + '.commot.spatial.h5ad')
    df_cellchat_filtered.to_csv(outdir + sample_name + ".cellchat.filtered.csv")

for sample_name in samples:
    outdir = "/public/workspace/stu21230110/AMP/11commot/" + sample_name + "/"
    adata_dis500 = sc.read(outdir + sample_name + ".commot.spatial.h5ad")
    df_cellchat_filtered = pd.read_csv(outdir + sample_name + ".cellchat.filtered.csv")
    pathways = list(set(df_cellchat_filtered.iloc[:,3]))
    for pathway in pathways:
        ct.tl.communication_direction(adata_dis500, database_name='cellchat', pathway_name=pathway, k=5)
        ct.pl.plot_cell_communication(adata_dis500, database_name='cellchat', pathway_name=pathway, 
                                      arrow_color = "#ff9900", 
                                      plot_method='grid', background_legend=True,scale=0.000006, ndsize=8, grid_density=0.4, summary='sender', background='image', clustering='leiden', cmap='Alphabet',normalize_v = True, normalize_v_quantile=0.995)
        plt.savefig(outdir + sample_name + '.' + pathway + '.signal.arrow.grid.spatial.pdf',bbox_inches = 'tight')



