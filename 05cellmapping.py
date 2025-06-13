import scanpy as sc
import squidpy as sq
import numpy as np
import pandas as pd
import os
import sys
import seaborn as sb
import matplotlib.pyplot as plt
import matplotlib as mpl
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


def plot_spatial_cm(adata_vis, gene,group, vmins, vcenters, vmaxs, img_keys, cmaps, palettes):
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(14,4), )
    plt.suptitle("                                                           Control                           ", y=1.05)
    with plt.rc_context({'axes.facecolor':  'white','figure.figsize': [4, 4]}):
        slide = select_slide(adata_vis, "C1")
        sc.pl.spatial(slide, cmap=cmaps,color=gene,size=1.3,img_key=img_keys,vmin=vmins, use_raw=False,layer="log1p",groups=group, vmax=vmaxs, ax=ax1, show=False,vcenter= vcenters,palette=palettes, legend_loc=None)
        slide = select_slide(adata_vis, "C2")
        sc.pl.spatial(slide, cmap=cmaps,color=gene,size=1.3,img_key=img_keys,vmin=vmins, use_raw=False,layer="log1p",groups=group, vmax=vmaxs, ax=ax2, show=False,vcenter= vcenters,palette=palettes, legend_loc=None)
        slide = select_slide(adata_vis, "C3")
        sc.pl.spatial(slide, cmap=cmaps,color=gene,size=1.3,img_key=img_keys,vmin=vmins, use_raw=False,layer="log1p",groups=group, vmax=vmaxs, ax=ax3, show=False,vcenter= vcenters,palette=palettes, legend_loc=None)
    plt.savefig(outdir + gene + ".control.spatial.pdf", bbox_inches = "tight")
    fig, (ax4,ax5,ax6) = plt.subplots(1, 3, figsize=(14,4), )
    plt.suptitle("                                                           Pain                           ", y=1.05)
    with plt.rc_context({'axes.facecolor':  'white','figure.figsize': [4, 4]}):
        slide = select_slide(adata_vis, "P1")
        sc.pl.spatial(slide, cmap=cmaps,color=gene,size=1.3,img_key=img_keys,vmin=vmins, use_raw=False,layer="log1p",groups=group, vmax=vmaxs, ax=ax4, show=False,vcenter= vcenters,palette=palettes, legend_loc=None)
        slide = select_slide(adata_vis, "P2")
        sc.pl.spatial(slide, cmap=cmaps,color=gene,size=1.3,img_key=img_keys,vmin=vmins, use_raw=False,layer="log1p",groups=group, vmax=vmaxs, ax=ax5, show=False,vcenter= vcenters,palette=palettes, legend_loc=None)
        slide = select_slide(adata_vis, "P3")
        sc.pl.spatial(slide, cmap=cmaps,color=gene,size=1.3,img_key=img_keys,vmin=vmins, use_raw=False,layer="log1p",groups=group, vmax=vmaxs, ax=ax6, show=False,vcenter= vcenters,palette=palettes)
    plt.savefig(outdir + gene + ".pain.spatial.pdf", bbox_inches = "tight")


outdir = "/public/workspace/stu21230110/AMP/05cellmapping/"

adata_vis = sc.read("/public/workspace/stu21230110/AMP/01data/spatial.nmf.niche.frequency.h5ad")

gene = "Niche_NMF"
vmins = None
vcenters= None
vmaxs = 1
img_keys='lowres'
cmaps= cmap_transparent2
Niche_NMF_palette = dict(zip(sorted(adata_vis.obs.Niche_NMF.cat.categories), adata_vis.uns["Niche_NMF_colors"]))
palettes=Niche_NMF_palette
group= None
from cell2location.utils import select_slide
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(14,4), )
plt.suptitle("                                                   Painfree                           ", y=1.05)
with plt.rc_context({'axes.facecolor':  'white','figure.figsize': [4, 4]}):
    slide = select_slide(adata_vis, "C1")
    sc.pl.spatial(slide, cmap=cmaps,color=gene,size=1.3,img_key=img_keys,vmin=vmins, use_raw=False,layer="log1p",groups=group, vmax=vmaxs, ax=ax1, show=False,vcenter= vcenters,palette=palettes, legend_loc=None)
    slide = select_slide(adata_vis, "C2")
    sc.pl.spatial(slide, cmap=cmaps,color=gene,size=1.3,img_key=img_keys,vmin=vmins, use_raw=False,layer="log1p",groups=group, vmax=vmaxs, ax=ax2, show=False,vcenter= vcenters,palette=palettes, legend_loc=None)
    slide = select_slide(adata_vis, "C3")
    sc.pl.spatial(slide, cmap=cmaps,color=gene,size=1.3,img_key=img_keys,vmin=vmins, use_raw=False,layer="log1p",groups=group, vmax=vmaxs, ax=ax3, show=False,vcenter= vcenters,palette=palettes, legend_loc=None)

plt.savefig(outdir + "spatial_plot_Niche_NMF_tissue_control.pdf",bbox_inches='tight')

fig, (ax4, ax5, ax6) = plt.subplots(1, 3, figsize=(14,4), )
plt.suptitle("                                                   Pain                               ", y=1.05)
with plt.rc_context({'axes.facecolor':  'white','figure.figsize': [4, 4]}):
    slide = select_slide(adata_vis, "P1")
    sc.pl.spatial(slide, cmap=cmaps,color=gene,size=1.3,img_key=img_keys,vmin=vmins, use_raw=False,layer="log1p",groups=group, vmax=vmaxs, ax=ax4, show=False,vcenter= vcenters,palette=palettes, legend_loc=None)
    slide = select_slide(adata_vis, "P2")
    sc.pl.spatial(slide, cmap=cmaps,color=gene,size=1.3,img_key=img_keys,vmin=vmins, use_raw=False,layer="log1p",groups=group, vmax=vmaxs, ax=ax5, show=False,vcenter= vcenters,palette=palettes, legend_loc=None)
    slide = select_slide(adata_vis, "P3")
    sc.pl.spatial(slide, cmap=cmaps,color=gene,size=1.3,img_key=img_keys,vmin=vmins, use_raw=False,layer="log1p",groups=group, vmax=vmaxs, ax=ax6, show=False,vcenter= vcenters,palette=palettes)

plt.savefig(outdir + "spatial_plot_Niche_NMF_tissue_Pain.pdf",bbox_inches='tight')


sampleID = ['C1','C2','C3','P1','P2','P3']
cts = ['B_cells', 'Epithelial_cells', 'Tissue_stem_cells',"T_cells"]

for sample in sampleID:
    slide = select_slide(adata_vis, sample)
    with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
        sc.pl.spatial(slide, cmap= gray_blue,##bwr(蓝红) hsv(变化明显) PuRd（红）magma
        # show first 8 cell types
        color=cts,
        ncols = 2,
        size=1.3,
        img_key='hires',
        # limit color scale at 99.2% quantile of cell abundance
        vmin=0.05, vmax='p99.2')
    plt.savefig(outdir  + sample + '_spatial.pdf',bbox_inches = 'tight')


# Now we use cell2location plotter that allows showing multiple cell types in one panel
from cell2location.plt import plot_spatial

# select up to 6 clusters
clust_labels = ['B_cells', 'Epithelial_cells', 'Tissue_stem_cells']
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

from cell2location.plt import plot_spatial

for sample in sampleID:
    slide = select_slide(adata_vis, sample)
    with mpl.rc_context({'figure.figsize': (15, 15)}):
        fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
        show_img=True,
        # 'fast' (white background) or 'dark_background'
        style='fast',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=7,
        colorbar_position='right'
#        reorder_cmap=[1,0,3,5], #[4, 3, 2,1,0,6,5],
    )
    plt.savefig(outdir + sample + "spatial_co_mapping.pdf",bbox_inches = "tight")
