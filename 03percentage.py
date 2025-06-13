import scanpy as sc
import squidpy as sq
import numpy as np
import pandas as pd
import os
import sys
import seaborn as sb
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42
from matplotlib import colors
import scvi
import anndata as ad
import warnings
warnings.filterwarnings("ignore")
from collections import Counter
import ipywidgets as widgets
from ipywidgets import interact, interact_manual
plt.rcParams['figure.figsize'] = (6, 6)
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

sc.set_figure_params(dpi=100, vector_friendly=True)
def mysize(w, h, d):
    fig, ax = plt.subplots(figsize = (w, h), dpi = d)
    return(fig.gca())

plt.rcParams['figure.figsize'] = (6, 5)
sc.set_figure_params(dpi=100, vector_friendly=True)

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

print(f"numpy=={np.__version__}")
#numpy==1.26.4
print(f"squidpy=={sq.__version__}")
#squidpy==1.6.3

import cell2location as c2l
from cell2location.utils import select_slide


outdir = "/public/workspace/stu21230110/AMP/03percentage/"

adata_vis = sc.read("/public/workspace/stu21230110/AMP/01data/spatial.nmf.plots.h5ad")

adata_vis.obs['treatment'] = adata_vis.obs['sample'].apply(lambda x: 'pain' if x.startswith('P') else 'control' if x.startswith('C') else None)
adata_vis.obs['sampleID'] = adata_vis.obs['sample'].copy()

cts = ['B_cells', 'Endothelial_cells', 'Epithelial_cells', 'Macrophage',
       'Monocyte', 'Myeloid_progenitor_cells', 'NK_cells',
       'Smooth_muscle_cells', 'Stromal_cells', 'T_cells', 'Tissue_stem_cells']

tab = sc.get.obs_df(adata_vis, keys= ["sampleID","Niche_NMF","treatment"]+ cts)

ct_label = "cell_type"
group_by = "sampleID"
xlabel = "Niche_NMF"
ctss = cts
#cm = "sampleID"

tab = tab.loc[:, ctss + [xlabel]].copy()
tab.head()

test_tab2 = tab.copy()
test_tab = test_tab2.groupby([xlabel]).mean().reset_index()
print(test_tab.shape)
test_tab.head()

test_tab.index = test_tab["Niche_NMF"].tolist()
test_tab = test_tab.drop("Niche_NMF", axis=1)
test_tab

test_tab2 = test_tab.div(test_tab.sum(axis=0), axis=1)
test_tab2

cts2 = ['Stromal_cells',
        'Endothelial_cells',
        'Macrophage', 'Monocyte', 'NK_cells', 'T_cells','B_cells',
        'Epithelial_cells','Tissue_stem_cells',
        'Myeloid_progenitor_cells',
        'Smooth_muscle_cells']

test_tab3 = test_tab2.copy()
test_tab3 = test_tab3[cts2]
test_tab3 = test_tab3.iloc[[0,1,3,4,2,5]] 
test_tab4 = test_tab3*100
sb.reset_defaults()
plt.figure(figsize=(15, 2.25))
sb.heatmap(test_tab4, annot=True, cmap="Reds", fmt=".0f",linewidths=0.25, linecolor='grey', )
plt.rcParams['pdf.fonttype'] = 42
plt.savefig(outdir + "NMF_niches_celltype_percentages.pdf",bbox_inches = "tight")

adata_vis.write("/public/workspace/stu21230110/AMP/01data/spatial.nmf.frequency.h5ad")
