conda activate squidpy
pip install pytables
conda install anaconda::pytables

#调整过格式
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

os.chdir("/public/workspace/stu21230110/AMP/01data/rds/")

#调整过格式
adata = sc.read("/public/workspace/stu21230110/AMP/01data/spatial.nmf.niche.frequency.h5ad")
def huage_adata_info(adata):
    mat=pd.DataFrame(data=adata.X.todense(),index=adata.obs_names,columns=adata.var_names)
    mat.to_csv("mat.csv")
    meta=pd.DataFrame(data=adata.obs)
    meta.to_csv('metadata.tsv',sep="\t")
    cord=pd.DataFrame(data=adata.obsm['spatial'],index=adata.obs_names,columns=['x','y'])
    cord.to_csv('position_'+'spatial'+'.tsv',sep="\t")
    umap=pd.DataFrame(data=adata.obsm["X_umap"],index=adata.obs_names,columns=['x','y'])
    umap.to_csv('position_'+"X_umap"+'.tsv',sep="\t")

huage_adata_info(adata)