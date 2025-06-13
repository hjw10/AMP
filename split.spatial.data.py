import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import scanpy as sc
import random
from collections import Counter
import matplotlib
matplotlib.use('Agg')

#change single sample spatial

adata_vis = sc.read("/public/workspace/stu21230110/AMP/01data/spatial.nmf.niche.frequency.h5ad")

outdir = "/public/workspace/stu21230110/AMP/01data/"

samples = ['C1', 'C2', 'C3', 'P1', 'P2', 'P3']

sample_name = samples[5]

adata = adata_vis[adata_vis.obs['sampleID'] == sample_name].copy() 

adata.uns['spatial'] = {
    'P3': {
        'coordinates': adata.obsm['spatial'].copy(),  # 坐标
        'images': adata.uns['spatial']['P3']['images'].copy(),  # 图像
        'scalefactors': adata.uns['spatial']['P3']['scalefactors'].copy()  # 缩放因子
    }
}
adata.write(outdir + sample_name + '.spatial.h5ad')