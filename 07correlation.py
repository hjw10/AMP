import scanpy as sc
import squidpy as sq
import numpy as np
import pandas as pd
import os
import sys
import seaborn as sb
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['figure.figsize'] = (6, 5)
sc.set_figure_params(dpi=100, vector_friendly=True)
import scvi
import anndata as ad
import warnings
warnings.filterwarnings("ignore")
from collections import Counter
import ipywidgets as widgets
from ipywidgets import interact, interact_manual
from IPython.core.display import display, HTML
import random

#Define a colour map for gene expression
colors2 = plt.cm.Reds(np.linspace(0, 1, 128))
colors3 = plt.cm.Greys_r(np.linspace(0.7,0.8,20))
#colorsComb = np.vstack([colors3, colors2])
#mymap = colors.LinearSegmentedColormap.from_list('my_colormap', colorsComb)
from matplotlib import colors
colorsComb = np.vstack([plt.cm.Reds(np.linspace(0, 1, 128)), plt.cm.Greys_r(np.linspace(0.7, 0.8, 0))])
mymap = colors.LinearSegmentedColormap.from_list('my_colormap', colorsComb)

# Helper function to split list in chunks
def chunks(lista, n):
    for i in range(0, len(lista), n):
        yield lista[i:i + n]        
        plt.rcParams['figure.figsize'] = (6, 5)

def mysize(w, h, d):
    fig, ax = plt.subplots(figsize = (w, h), dpi = d)
    return(fig.gca())

## frequently used variables
from matplotlib import colors
import matplotlib.pyplot as plt
colorsComb = np.vstack([plt.cm.Reds(np.linspace(0, 1, 128)), plt.cm.Greys_r(np.linspace(0.7, 0.8, 0))])
mymap = colors.LinearSegmentedColormap.from_list('my_colormap', colorsComb)

## Along these Lines, a colourmap diverging from gray to red
gray_red = colors.LinearSegmentedColormap.from_list("grouping", ["lightgray", "red", "darkred"], N = 128)
## Some more Colour Maps
gray_violet = colors.LinearSegmentedColormap.from_list("grouping", ["lightgray", "mediumvioletred", "indigo"], N = 128)
gray_blue = colors.LinearSegmentedColormap.from_list("grouping", ["lightgray", "cornflowerblue", "darkblue"], N = 128)

def mysize(w, h, d):
    fig, ax = plt.subplots(figsize = (w, h), dpi = d)
    return(fig.gca())
#plt.rcParams['figure.figsize'] = (6, 5)
#sc.set_figure_params(dpi=120, vector_friendly=True)

import matplotlib.colors as colors
c_low = colors.colorConverter.to_rgba('orange', alpha = 0)
c_high = colors.colorConverter.to_rgba('red',alpha = 1)
cmap_transparent = colors.LinearSegmentedColormap.from_list('rb_cmap',[c_low, c_high], 512)
import matplotlib.colors as colors
c_low2 = colors.colorConverter.to_rgba('green', alpha = 0)
c_high2 = colors.colorConverter.to_rgba('darkblue',alpha = 1)
cmap_transparent2 = colors.LinearSegmentedColormap.from_list('rb_cmap',[c_low2, c_high2], 512)


outdir = "/public/workspace/stu21230110/AMP/07correlation/"

adata_vis = sc.read("/public/workspace/stu21230110/AMP/01data/spatial.nmf.niche.frequency.h5ad")

cts = adata_vis.obsm['q05_cell_abundance_w_sf'].columns.values.copy()
cts = ['B_cells', 'Endothelial_cells', 'Epithelial_cells', 'Macrophage',
       'Monocyte', 'Myeloid_progenitor_cells', 'NK_cells',
       'Smooth_muscle_cells', 'Stromal_cells', 'T_cells',
       'Tissue_stem_cells']
#array(['B_cells', 'Endothelial_cells', 'Epithelial_cells', 'Macrophage',
#       'Monocyte', 'Myeloid_progenitor_cells', 'NK_cells',
#       'Smooth_muscle_cells', 'Stromal_cells', 'T_cells',
#       'Tissue_stem_cells'], dtype=object)
#adata_vis.obs = pd.concat([adata_vis.obs,adata_vis.obsm['q05_cell_abundance_w_sf']], axis=1)
sc.pl.matrixplot(adata_vis, var_names= cts, groupby="Niche_NMF",cmap=mymap, swap_axes=False, standard_scale="var",dendrogram=True)#, save="matrix_Niche_NMF_celltypes_filtered.pdf")# figsize=(8,20), )
plt.savefig(outdir + "matrix_Niche_NMF_celltypes_filtered.pdf",bbox_inches = "tight")

cts2 = ['Endothelial_cells',
        'Epithelial_cells','Tissue_stem_cells',
        'Myeloid_progenitor_cells',
        'Macrophage','B_cells','Monocyte','T_cells',
        'Stromal_cells',
        'Smooth_muscle_cells']
sc.pl.matrixplot(adata_vis, var_names= cts2, groupby="Niche_NMF",cmap=mymap, swap_axes=False, standard_scale="var",dendrogram=True)#, save="matrix_Niche_NMF_celltypes_filtered.pdf")# figsize=(8,20), )
plt.savefig(outdir + "matrix_Niche_NMF_celltypes_filtered.pdf",bbox_inches = "tight")

#correlation
test_tab = sc.get.obs_df(adata_vis, keys = cts)
#test_tab


import seaborn as sns
from scipy.spatial import distance
from scipy.cluster import hierarchy

correlations = test_tab.corr()
correlations_array = np.asarray(test_tab.corr())

row_linkage = hierarchy.linkage(
    distance.pdist(correlations_array), method='average')

col_linkage = hierarchy.linkage(
    distance.pdist(correlations_array.T), method='average')

sns.set(font_scale=0.975)
p=sns.clustermap(correlations, row_linkage=row_linkage, col_linkage=col_linkage, method="average", figsize=(11, 11), cmap=sns.diverging_palette(250, 10, n=100),
               vmin=-1, vmax=1, center=0, cbar_pos=(-0.1, .2, .03, .4))

p.savefig(outdir + "clustermap_celltype_correlation.pdf",bbox_inches = "tight")
