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



outdir = "/public/workspace/stu21230110/AMP/04frequency/"

adata_vis = sc.read("/public/workspace/stu21230110/AMP/01data/spatial.nmf.frequency.h5ad")

adata_vis.uns['Niche_NMF_colors'] = ['#1f77b4ff', '#ff7f0eff', '#2ca02cff', '#d62728ff', '#9467bdff', '#8c564bff']

Niche_NMF_palette = dict(zip(adata_vis.obs.Niche_NMF.cat.categories, adata_vis.uns["Niche_NMF_colors"]))

Niche_NMF_palette = {'Niche1': '#1f77b4', 
                     'Niche2': '#ff7f0e', 
                     'Niche3': '#2ca02c', 
                     'Niche4': '#d62728', 
                     'Niche5': '#9467bd', 
                     'Niche6': '#8c564b'}

cts = adata_vis.obsm['q05_cell_abundance_w_sf'].columns.values.copy()
cts

sc.pl.umap(adata_vis, color="Niche_NMF", palette=Niche_NMF_palette,)
plt.savefig(outdir + "Niche.umap.pdf", bbox_inches='tight')


def standart_lineplot(data, order, xlabel, ylabel, typ = None, gene = None, smooth = None, column = None, cols = None,
                      palette = None, width = 1, title = None, rotation = None, figsize = (15, 5), tick_size = None,
                      label_size = None,
                      order_smooth = 3, legend = None, conf_int = None, scatter = None, save = None):
    if smooth: 
        ## Possible to set alpha of scatter with scatter_kws={'alpha': 0.1}
        if typ:
            cat = sb.lmplot(data = data, x = xlabel, y = gene, ci = conf_int, order = order_smooth,
                            scatter = scatter, hue = typ, truncate = True, palette = cols)
        else:
            cat = sb.lmplot(data = data, x = xlabel, y = gene, ci = conf_int, order = order_smooth, scatter = scatter, 
                            palette = cols)  
    else:
        ## Removed Parameter order = order, as order should be given numerically anyways.
        if typ:
            cat = sb.catplot(data = data, x = xlabel, y = gene, linestyles = "-", kind = "point", hue = typ,
                             palette = cols)
        else:
            cat = sb.catplot(data = data, x = xlabel, y = gene, linestyles = "-", kind = "point", palette = cols)
        if scatter:
            cat2 = sb.stripplot(data = data, x = xlabel, y = gene, palette = ["black"], hue = typ, size = 7)
            if typ:
                cat2.legend_.remove()
    cat.set(xticks = np.unique(data.loc[:, xlabel]))
    cat.set_xticklabels(order)
    cat.fig.set_size_inches(figsize)
    if rotation:
        cat.ax.set_xticklabels(order, rotation = 'vertical')
    cat.ax.set_title(title, size = label_size)
    cat.ax.set_xlabel(xlabel, size = label_size)
    cat.ax.set_ylabel(ylabel, size = label_size)
    cat.ax.tick_params(labelsize = tick_size)
    if save:
        cat.fig.savefig("%s" %save, bbox_inches = "tight")
        print("Saving figure to %s" %save)
    plt.show()
    plt.close()
    

def split_boxplot(tab, order, xlabel, ylabel, column = None, value = None, cols = None, width = 1, title = None,
                  figsize = (15, 6), legend_loc = None, jitter = None, save = None):
    fig, ax = plt.subplots()
    fig.set_size_inches(figsize)
    if cols is not None:
        fig = sb.boxplot(data = tab, hue = value, x = xlabel, y = column, order = order, width = width, palette = cols)
    else:
        fig = sb.boxplot(data = tab, hue = value, x = xlabel, y = column, order = order, width = width)
    if jitter is not None:
        fig = sb.swarmplot(data = tab, color = "black", x = xlabel, y = column, order = order)    
    if value is not None:
        plt.legend(loc = legend_loc)
    if title:
        fig.set_title(title, size = 15)
    fig.set_xlabel(xlabel, size = 15)
    fig.set_ylabel(ylabel, size = 15)
    fig.tick_params(labelsize = 12) 
    if save:
        fig.get_figure().savefig("%s" %(save))
    plt.show()
    plt.close()


## Relative Frequencies
def plot_relFreq(relFreqs, cluster, cols, order, xlabel = "days", condition = "batch", legend_loc = "upper right",
                 figsize = (15,6), width = .5, jitter = None, save = None):
    ## Subset according to order
    relFreqs = relFreqs.loc[relFreqs[xlabel].isin(order)]
    split_boxplot(relFreqs, order = order, xlabel = xlabel, ylabel = "relative frequency", value = condition,
                  column = cluster, cols = cols, width = width, title = cluster, figsize = figsize,
                  legend_loc = legend_loc, jitter = jitter, save = save)
    

## New adapted Version
def calc_relFreq(a, group_by = "cell_type", xlabel = "days", condition = "batch"):
    freqs = a.obs.groupby(["sampleID", group_by]).size()
    samples = np.unique(a.obs["sampleID"])
    ind = a.obs[group_by].cat.categories
    relFreqs = [freqs[ident] / sum(freqs[ident]) for ident in samples]
    relFreqs = pd.DataFrame(relFreqs, columns = ind, index = samples).fillna(0)
    #relFreqs[xlabel] = grouping.loc[samples, xlabel]  ## when using Grouping Table
    cell_types = {}
    combis = a.obs.groupby(["sampleID", xlabel]).groups.keys()
    for c in combis:
        cell_types[c[0]] = c[1]
    relFreqs[xlabel] = [cell_types[l] for l in relFreqs.index]
    ## Todo, add for condition
    if condition:
        combis = a.obs.groupby(["sampleID", condition]).groups.keys()
        for c in combis:
            cell_types[c[0]] = c[1]
        relFreqs[condition] = [cell_types[l] for l in relFreqs.index]
    return relFreqs

def plot_gene_boxplot(tab, xlabel = "cell_type", condition = None, figsize = (10, 5), legend = True,
                      palette = ["gray", "red", "blue"], score = "Axin2", size = 4, scatter = None,
                      rotate = False, width = 0.7, save = None):
    sb.set_style("ticks")  ## show ticks
    fig, ax = plt.subplots()
    fig.set_size_inches(figsize)
    sf = False if scatter else True
    if condition:
        fig = sb.boxplot(data = tab, x = xlabel, y = score, width = width, hue = condition,
                         showfliers = sf, palette = palette)
        if scatter:
            fig = sb.stripplot(data = tab, x = xlabel, y = score, palette = ["black"],
                               size = size, hue = condition, dodge = True)
    else:
        fig = sb.boxplot(data = tab, x = xlabel, y = score, width = width, showfliers = sf, palette = palette)
        if scatter:
            fig = sb.stripplot(data = tab, x = xlabel, y = score, palette = ["black"], size = size, dodge = True)
    if rotate:
        fig.set_xticklabels(fig.get_xticklabels(), rotation = 90)
    else:
        fig.set_xticklabels(fig.get_xticklabels())
    if condition:
        ## Remove the dots from the legend if scatter is used
        handles, labels = ax.get_legend_handles_labels()
        n = tab.loc[:, condition].nunique()
        plt.legend(handles[0:n], labels[0:n], bbox_to_anchor=(1.05, 1), loc = 2, borderaxespad = 0.)
    #else:
    #    ax.legend_.remove()
    plt.setp(ax.artists, edgecolor = "black")
    plt.setp(ax.lines, color = "black")
    sb.despine()  ## to not show ouline box
    if save:
        print("Saving to %s" %save)
        plt.savefig(save, bbox_to_anchor = "tight")
    plt.show()

def calc_relFreq_per_cluster(a, group_by = "cell_type", xlabel = "days", condition = None):
    freqs = a.obs.groupby([group_by, xlabel]).size()
    celltypes = np.unique(a.obs[group_by])
    ind = a.obs[xlabel].cat.categories
    relFreqs = [freqs[ident] / sum(freqs[ident]) for ident in celltypes]
    relFreqs = pd.DataFrame(relFreqs, columns = ind, index = celltypes).fillna(0)
    cell_types = {}
    combis = a.obs.groupby([group_by, xlabel]).groups.keys()
    for c in combis:
        cell_types[c[0]] = c[1]
    relFreqs[group_by] = relFreqs.index
    ## Todo, add for condition
    if condition:
        combis = a.obs.groupby([group_by, condition]).groups.keys()
        for c in combis:
            cell_types[c[0]] = c[1]
        relFreqs[condition] = [cell_types[l] for l in relFreqs.index]
    return relFreqs

def plot_cluster_composition(relFreqs, xlabel = "name", figsize = (6, 10), 
                             width = 0.8, order = None, errbar = None, labelsize = 15, ticksize = 13,
                             capsize = None,
                             margins = (0.02, 0.04), cols = None, save = None):
    """
    Given a relative Frequency table, plot the values as stacked barplot.
    Parameters
    ----------
    relFreqs
        Pandas DataFrame containing relative Frequencies
    figsize (default (6, 10)
        Set size of Figure in form of (width, height)
    xlabel
        Label to group by on the x-axis
    order (default None)
        Specify as List if manual ordering is desired
    width (default 0.8)
        Specify width of the Bars
    errbar (default None)
        Set to true to plot on top (only possible when grouping the frequencies)
    capsize (default None)
        Size of the horizontal lines of the errorbar
    labelsize (default 15)
        Set size of x and y asix legends
    ticksize (default 13)
        Set size of x axis ticks and legend
    save (default None)
        Set full file path in order to save Figure (e.g. /path/to/file.pdf)
    """
    import matplotlib.patches as mpatches
    patches = []
    fig, ax = plt.subplots()
    fig.set_size_inches(figsize)
    order = np.unique(relFreqs.loc[:, xlabel]) if order is None else order
    ci = 95 if errbar else None
    ax.margins(margins[0], margins[1])
    cell_types = np.flip([col for col in relFreqs.columns if col not in ["identifier", xlabel]])
    #cell_types = np.flip(np.setdiff1d(relFreqs.columns, ["identifier", xlabel]))
    bars = pd.DataFrame(index = order, data = np.zeros(len(order)))
    plot_data = pd.DataFrame(relFreqs.loc[:, xlabel])
    for i, typ in enumerate(cell_types):
        sum_up = [relFreqs.loc[:, typ].values[i] + bars.loc[g].values[0] for i, g in enumerate(relFreqs.loc[:, xlabel])]
        plot_data[typ] = sum_up
        bars.iloc[:, 0] = bars.iloc[:, 0] + relFreqs.loc[:, [typ, xlabel]].groupby(xlabel).mean().loc[order, typ]
    for i, typ in enumerate(reversed(cell_types)):
        fig = sb.barplot(data = plot_data, x = xlabel, y = typ, order = order,
                         ci = ci, errcolor = "black", color = cols[i], capsize = capsize)
        patches.append(mpatches.Patch(color = cols[i], label = typ))
    ax.set_xlabel(xlabel, size = labelsize)
    ax.set_ylabel("relative frequency", size = labelsize)
    ax.tick_params(labelsize = ticksize)
    ax.set_xticklabels(labels = order, rotation = 'vertical')
    ## In order to change width of bars...
    for bar in fig.patches:
        centre = bar.get_x() + bar.get_width()/2.
        bar.set_x(centre - width/2.)
        bar.set_width(width)
    plt.legend(handles = patches, loc = "center left", bbox_to_anchor=(1.02, 0.5),
               prop = {"size": ticksize}, frameon = False)
    if save:
        plt.savefig("%s" %(save), bbox_inches='tight')
        print("Saving Figure to %s" %save)
    plt.show()
    plt.close()


## Get Relative Frequencies per identifier
adata_vis.uns['treatment_colors'] = ['#1f77b4ff', '#ff7f0eff']

xlabel = "treatment"
cell_type_label = "Niche_NMF"
cols = adata_vis.uns["%s_colors" %xlabel]

relFreqs = calc_relFreq_per_cluster(adata_vis, group_by = cell_type_label, xlabel = xlabel)
         control      pain Niche_NMF
Niche1  0.494816  0.505184    Niche1
Niche2  0.823529  0.176471    Niche2
Niche3  0.684211  0.315789    Niche3
Niche4  0.298421  0.701579    Niche4
Niche5  0.629630  0.370370    Niche5
Niche6  0.579817  0.420183    Niche6

plot_cluster_composition(relFreqs, xlabel = cell_type_label, figsize = (5, 4), cols = cols, 
                               margins = (0.02, 0.04), width = 0.7, order = None,save = outdir + 'relfreq_Niches.pdf')

#sc.pl.umap(adata_vis, color=["treatment",],vmax="p99",ncols=2, wspace=0.25,s=5,save = "UMAP_treatment.pdf")

## Get Relative Frequencies per identifier
xlabel = "treatment"
cell_type_label = "Niche_NMF"
cols = adata_vis.uns["%s_colors" %xlabel]

relFreqs = calc_relFreq_per_cluster(adata_vis, group_by = cell_type_label, xlabel = xlabel)
relFreqs

plt.rcParams['pdf.fonttype'] = 42

#test
xlabel = "treatment"
cell_types_label = "Niche_NMF"
cols = adata_vis.uns["%s_colors" %cell_types_label]
celltypes = adata_vis.obs[cell_types_label].cat.categories

relFreqs = calc_relFreq(adata_vis, group_by = cell_types_label, xlabel = xlabel, condition = None)
relFreqs
>>> relFreqs
      Niche1    Niche2    Niche5    Niche3    Niche4    Niche6 treatment
C1  0.776860  0.000918  0.005051  0.005969  0.110193  0.101010   control
C2  0.727769  0.006001  0.003273  0.000000  0.212220  0.050736   control
C3  0.800266  0.001332  0.000000  0.000000  0.194407  0.003995   control
P1  0.735126  0.000000  0.008150  0.000000  0.254279  0.002445      pain
P2  0.630863  0.000716  0.000000  0.000358  0.313283  0.054780      pain
P3  0.592698  0.000570  0.000000  0.002852  0.362236  0.041643      pain

relFreqs2 = relFreqs.iloc[:,[0,1,3,4,2,5,6]]
relFreqs2
      Niche1    Niche2    Niche3    Niche4    Niche5    Niche6 treatment
C1  0.776860  0.000918  0.005969  0.110193  0.005051  0.101010   control
C2  0.727769  0.006001  0.000000  0.212220  0.003273  0.050736   control
C3  0.800266  0.001332  0.000000  0.194407  0.000000  0.003995   control
P1  0.735126  0.000000  0.000000  0.254279  0.008150  0.002445      pain
P2  0.630863  0.000716  0.000358  0.313283  0.000000  0.054780      pain
P3  0.592698  0.000570  0.002852  0.362236  0.000000  0.041643      pain

plot_cluster_composition(relFreqs2, xlabel = xlabel, figsize = (3, 6), order = None,
                               errbar = False, cols = cols, width = 0.8,save = outdir + "relfreq_treatment.pdf")


xlabel = "sampleID"
cell_types_label = "Niche_NMF"
cols = adata_vis.uns["%s_colors" %cell_types_label]
celltypes = adata_vis.obs[cell_types_label].cat.categories

relFreqs = calc_relFreq(adata_vis, group_by = cell_types_label, xlabel = xlabel, condition = None)
relFreqs
      Niche1    Niche2    Niche5    Niche3    Niche4    Niche6 sampleID
C1  0.776860  0.000918  0.005051  0.005969  0.110193  0.101010       C1
C2  0.727769  0.006001  0.003273  0.000000  0.212220  0.050736       C2
C3  0.800266  0.001332  0.000000  0.000000  0.194407  0.003995       C3
P1  0.735126  0.000000  0.008150  0.000000  0.254279  0.002445       P1
P2  0.630863  0.000716  0.000000  0.000358  0.313283  0.054780       P2
P3  0.592698  0.000570  0.000000  0.002852  0.362236  0.041643       P3

relFreqs2 = relFreqs.iloc[:,[0,1,3,4,2,5,6]]
plt.rcParams['figure.figsize'] = (6, 5)
sc.set_figure_params(dpi=100, vector_friendly=True)

plot_cluster_composition(relFreqs2, xlabel = xlabel, figsize = (5, 4), cols = cols, 
                               margins = (0.02, 0.04), width = 0.7, order = None, save = outdir + "relfreq_sampleID.pdf")

adata_vis.write("/public/workspace/stu21230110/AMP/01data/spatial.nmf.niche.frequency.h5ad")