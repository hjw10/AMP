##NMF非负矩阵分解
##画niche分布
##
import os
import sys
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

import pandas as pd
import scvi
import squidpy as sq
import cell2location as c2l
from cell2location.utils import select_slide
from cell2location.plt.plot_heatmap import clustermap
import warnings
warnings.filterwarnings("ignore")

from collections import Counter
import ipywidgets as widgets
from ipywidgets import interact, interact_manual
from IPython.core.display import display, HTML
import random

print(f"numpy=={np.__version__}")
#numpy==1.26.4
print(f"squidpy=={sq.__version__}")
#squidpy==1.6.3

#Define a colour map for gene expression
colors2 = plt.cm.Reds(np.linspace(0, 1, 128))
colors3 = plt.cm.Greys_r(np.linspace(0.7,0.8,20))
#colorsComb = np.vstack([colors3, colors2])
#mymap = colors.LinearSegmentedColormap.from_list('my_colormap', colorsComb)
from matplotlib import colors
colorsComb = np.vstack([plt.cm.Reds(np.linspace(0, 1, 128)), plt.cm.Greys_r(np.linspace(0.7, 0.8, 0))])
mymap = colors.LinearSegmentedColormap.from_list('my_colormap', colorsComb)

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




outdir = "/public/workspace/stu21230110/AMP/02nmf/"

adata_vis = sc.read("/public/workspace/stu21230110/AMP/01data/20250512adata_vis.h5ad")

adata = sc.read("/public/workspace/stu21230110/AMP/01data/sc.h5ad")

cell_anno = pd.read_csv("/public/workspace/stu21230110/AMP/01data/st_cell2location_res.csv")

samples = ['C1','C2','C3','P1','P2','P3']

adata_vis.obs = adata_vis.obs.rename(columns={
    col: col.replace('means_per_cluster_mu_fg_', '') 
    for col in adata_vis.obs.columns 
    if col.startswith('means_per_cluster_mu_fg_')
})

adata_vis.var_names = adata_vis.var["SYMBOL"].tolist().copy()

adata.obs.singleRnew.cat.categories

cts = ['B_cells', 'Endothelial_cells', 'Epithelial_cells', 'Macrophage',
       'Monocyte', 'Myeloid_progenitor_cells', 'NK_cells',
       'Smooth_muscle_cells', 'Stromal_cells', 'T_cells', 'Tissue_stem_cells']

adata_vis.write("/public/workspace/stu21230110/AMP/01data/spatial.rename.h5ad")





outdir = "/public/workspace/stu21230110/AMP/02nmf/"

adata = sc.read("/public/workspace/stu21230110/AMP/01data/sc.h5ad")

samples = ['C1','C2','C3','P1','P2','P3']

adata_vis = sc.read("/public/workspace/stu21230110/AMP/01data/spatial.rename.h5ad")

results_folder = outdir + "/visium_cytassist_integrated/"

run_name = f'{results_folder}/cell2location_map_facts'

from cell2location import run_colocation
res_dict, adata_vis = run_colocation(
    adata_vis,
    model_name='CoLocatedGroupsSklearnNMF',
    train_args={
      'n_fact':  np.arange(3,20), # IMPORTANT: use a wider range of the number of factors (5-30)
      'sample_name_col': 'sample', # columns in adata_vis.obs that identifies sample
      'n_restarts': 3 # number of training restarts
    },
    export_args={'path': f'{run_name}/CoLocatedComb/'}
)
                  
#continue with 6 factors
results_folder = outdir + "visium_cytassist_integrated/"                  
run_name = f'{results_folder}/cell2location_map_6_facts'

from cell2location import run_colocation
res_dict, adata_vis = run_colocation(
    adata_vis,
    model_name='CoLocatedGroupsSklearnNMF',
    train_args={
      'n_fact':  np.arange(6,7), # IMPORTANT: use a wider range of the number of factors (5-30)
      'sample_name_col': 'sample', # columns in adata_vis.obs that identifies sample
      'n_restarts': 3 # number of training restarts
    },
    export_args={'path': f'{run_name}/CoLocatedComb/'}
)

cts_f = ['B_cells', 'Endothelial_cells', 'Epithelial_cells', 'Macrophage',
       'Monocyte', 'Myeloid_progenitor_cells', 'NK_cells',
       'Smooth_muscle_cells', 'Stromal_cells', 'T_cells', 'Tissue_stem_cells']

import pickle
# import models
def unpickle_model(path, mod_name):
    r""" Unpickle model
    """
    file = path + 'model_' + mod_name + ".p"
    mod1_ann = pickle.load(file = open(file, "rb"))
    return mod1_ann['mod']

n_fact = 6
mod_path = "/public/workspace/stu21230110/AMP/02nmf/visium_cytassist_integrated/cell2location_map_6_facts/CoLocatedComb/CoLocatedGroupsSklearnNMF_10535locations_11factors/models/"

mod_sk = unpickle_model(mod_path, f'n_fact{n_fact}')

filterval= 0.05
pretest_sk = mod_sk.cell_type_fractions[(mod_sk.cell_type_fractions['fact_0']>filterval)|(mod_sk.cell_type_fractions['fact_1']>filterval)|(mod_sk.cell_type_fractions['fact_2']>filterval)|(mod_sk.cell_type_fractions['fact_3']>filterval)|(mod_sk.cell_type_fractions['fact_4']>filterval)|(mod_sk.cell_type_fractions['fact_5']>filterval)]
pretest_sk

b_dev_sel =  cts_nmf
cts_nmf


fact_filt = ['fact_3','fact_2','fact_1','fact_5','fact_0','fact_4']
filterval= 0.05

test_cm = mod_sk.cell_type_fractions
test_sk = test_cm[(test_cm['fact_0']>filterval)|(test_cm['fact_1']>filterval)|(test_cm['fact_2']>filterval)|(test_cm['fact_3']>filterval)|(test_cm['fact_4']>filterval)|(test_cm['fact_5']>filterval)]
#test_sk = (test_sk.T / test_sk.max(1)).T

mpl.rc_file_defaults()
mpl.rcParams['pdf.fonttype'] = 42 # enables correct plotting of text
clustermap(test_sk.loc[b_dev_sel, fact_filt[::-1]],
           cluster_rows=True, cluster_cols=False, 
           figure_size=[5.9 + 0.12 * mod_sk.n_fact, 3 + 0.15 * mod_sk.n_var],
           fun_type='dotplot', array_size=None)

plt.savefig(outdir + "NMF_cytassist_visum_human_6facts.pdf", bbox_inches='tight')


def binarize_score_flexible(a, cts, groupby = "cell_type", percentile = 100, none_label = "unmapped"):
    tab = sc.get.obs_df(a, keys = cts + [groupby])
    tab = tab.groupby([groupby]).mean().reset_index()
    tab.set_index(groupby, inplace = True)
    ## Define cell type specific thresholds
    thresholds = dict.fromkeys(tab.columns)
    for ct in tab.columns:
        thresholds[ct] = np.percentile(tab.loc[:, ct], percentile)
    cluster_map = dict.fromkeys(a.obs.loc[:, groupby].cat.categories)
    for cluster in cluster_map.keys():
        assigned = []
        for ct in tab.columns:
            if tab.loc[cluster, ct] > thresholds[ct]:
                assigned.append(ct)
        cluster_map[cluster] = tab.loc[cluster].idxmax() 
    annot = [cluster_map[leid] for leid in a.obs.loc[:, groupby]]
    annot = [none_label if a == "" else a for a in annot]
    return annot

nfacts_adata_vis= ['mean_nUMI_factorsfact_0', 'mean_nUMI_factorsfact_1', 'mean_nUMI_factorsfact_2', 'mean_nUMI_factorsfact_3', 'mean_nUMI_factorsfact_4', 'mean_nUMI_factorsfact_5']

adata_vis.obs["spots"] = adata_vis.obs.index.copy().astype('category')

adata_vis.obs["mapped_spots"] = binarize_score_flexible(adata_vis, nfacts_adata_vis, groupby = "spots", percentile = 100)
cts_nmf = ['B_cells', 'Endothelial_cells', 'Epithelial_cells', 'Macrophage',
       'Monocyte', 'Myeloid_progenitor_cells', 'NK_cells',
       'Smooth_muscle_cells', 'Stromal_cells', 'T_cells', 'Tissue_stem_cells']
sc.pl.matrixplot(adata_vis, var_names= cts_nmf, groupby="mapped_spots",cmap=mymap, swap_axes=False, standard_scale=None,)# figsize=(8,20), )
plt.savefig(outdir + "celltypes.nmf.heatmap.pdf",bbox_inches = "tight")

dict([(x,x) for x in adata_vis.obs.mapped_spots.cat.categories])

cell_mapping_dict ={
 'mean_nUMI_factorsfact_0': 'Niche1',
 'mean_nUMI_factorsfact_1': 'Niche2',
 'mean_nUMI_factorsfact_3': 'Niche3',
 'mean_nUMI_factorsfact_4': 'Niche4',
 'mean_nUMI_factorsfact_2': 'Niche5',
 'mean_nUMI_factorsfact_5': 'Niche6'}

adata_vis.obs['Niche_NMF'] = adata_vis.obs.mapped_spots.map(cell_mapping_dict).astype('category')

# Get highly variable genes and create preprocessing (pp) copy of data for running the model
sc.pp.highly_variable_genes(adata_vis, flavor='seurat_v3', n_top_genes=6000, subset=False)
adata_pp = adata_vis[:, adata_vis.var["highly_variable"] == True].copy()
adata_pp

sc.tl.pca(adata_pp)
sc.pp.neighbors(adata_pp, use_rep='q05_cell_abundance_w_sf',
                n_neighbors = 15)
# Cluster spots into regions using scanpy
sc.tl.leiden(adata_pp, resolution=1.1)
sc.tl.umap(adata_pp, min_dist = 0.3, spread = 1)
adata_vis.obsm['X_pca'] = adata_pp.obsm['X_pca']
adata_vis.obsm['X_umap'] = adata_pp.obsm['X_umap']
adata_vis.obs["leiden"] = adata_pp.obs["leiden"]

sc.pl.umap(adata_vis, color=["Niche_NMF",],vmax="p99",ncols=2, wspace=0.25,s=7.5)
plt.savefig(outdir + "umap.niche.nmf.pdf",bbox_inches = "tight")

adata_vis.uns["Niche_NMF_colors"]
#['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']

Niche_NMF_palette = dict(zip(adata_vis.obs.Niche_NMF.cat.categories, adata_vis.uns["Niche_NMF_colors"]))
Niche_NMF_palette = {'Niche1': '#1f77b4', 
                     'Niche2': '#ff7f0e', 
                     'Niche3': '#2ca02c', 
                     'Niche4': '#d62728', 
                     'Niche5': '#9467bd', 
                     'Niche6': '#8c564b'}

gene = "Niche_NMF"
vmins = None
vcenters= None
vmaxs = 1
img_keys='lowres'
cmaps= cmap_transparent2
palettes=Niche_NMF_palette
group= None

plot_spatial_cm(adata_vis, gene,group, vmins, vcenters, vmaxs, img_keys, cmaps, palettes)

adata_vis.write("/public/workspace/stu21230110/AMP/01data/spatial.nmf.h5ad")

adata_vis = sc.read("/public/workspace/stu21230110/AMP/01data/spatial.nmf.h5ad")

cts2 = ['Endothelial_cells',
        'Epithelial_cells','Tissue_stem_cells',
        'Myeloid_progenitor_cells',
        'B_cells', 'Macrophage', 'NK_cells', 'Monocyte', 'T_cells',
        'Stromal_cells',
        'Smooth_muscle_cells']

sc.pl.matrixplot(adata_vis, var_names= cts2, groupby="Niche_NMF",cmap=mymap, swap_axes=False, 
                 standard_scale="var",dendrogram=True)

plt.savefig(outdir + "matrix_Niche_NMF_celltypes_filtered.pdf",bbox_inches = "tight")

#correlation
test_tab = sc.get.obs_df(adata_vis, keys = cts2)

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

p.savefig(outdir + "clustermap_celltype_correlation.pdf", bbox_inches = "tight")

#all markers calculation
## All Markers
sc.tl.rank_genes_groups(adata_vis, groupby = 'Niche_NMF', groups = "all", use_raw = False, method = "wilcoxon",n_genes=2000, pts= True)
sc.pl.rank_genes_groups(adata_vis)
plt.savefig(outdir + "Niche_NMF.allMarkers.pdf",bbox_inches = "tight")

## Combine into one Data Frame (comparable to Marker Table)
result = adata_vis.uns['rank_genes_groups']
allMarkers = []
for cluster in result['names'].dtype.names:
    current = pd.DataFrame({"gene": result["names"][cluster], "score": result["scores"][cluster],
                            "logfoldchange": result["logfoldchanges"][cluster], "pval": result["pvals"][cluster],
                            "pval_adj": result["pvals_adj"][cluster],
                            "pct_within":result["pts"].loc[result["names"][cluster]][cluster],
                            "pct_outside":result["pts_rest"].loc[result["names"][cluster]][cluster],  "cluster": cluster})
    allMarkers.append(current)

allMarkers = pd.concat(allMarkers)
allMarkers.head()

allMarkers.to_csv(outdir + "AllMarkers_Niche_NMF_2000genes.txt", sep = "\t", index = False)

marker3 = allMarkers
marker3 = marker3[marker3['pct_outside'] < 0.50]
marker3 = marker3[marker3['pct_within'] > 0.5]
marker3 = marker3[marker3['logfoldchange'] > 0.5]
#marker3=list(marker3['gene'])
## add average expression of markers to table
#marker3_expr = avg_expr_ct['Myofibroblasts']
#marker3_expr = marker3_expr[marker3_expr.index.isin(marker3['gene'])]
#marker3['avg_expr'] = list(marker3_expr)
## check distribution of average expression
#b.histplot(marker3_expr)
#len(marker3)
marker3

#Code zum filtern hast ja schon, die folgenden Zeilen wählen dann die top n gene per cluster aus (im gleichen ordering wie das dendrogram)
marker3.sort_values(by = ["cluster", "logfoldchange"], ascending = [True, False], inplace = True)

#Get ordering as in Dendrogram
#sc.tl.dendrogram(adata_vis, groupby = "leiden_25")
order = adata_vis.uns["dendrogram_Niche_NMF"]["dendrogram_info"]["ivl"]

#Get top genes from Markers Table
n_genes = 8
genes = []
for typ in order:
    curgenes = marker3.loc[marker3.cluster == typ, "gene"].values[0:n_genes]
    genes = genes + [g for g in curgenes if g not in genes]

## Create total-normalized layer of counts (each cell will have a total sum of 10,000)
adata_vis.layers["counts"] = adata_vis.X.copy()
adata_vis.layers["total_normalized"] = sc.pp.normalize_total(adata_vis, target_sum=1e4, inplace=False)["X"]
adata_vis.layers["log1p"] = sc.pp.log1p(adata_vis.layers["total_normalized"], copy=True)

sc.pl.matrixplot(adata_vis, var_names = genes, standard_scale = "var", groupby = "Niche_NMF", dendrogram = True,figsize=(16,2.5),layer="log1p")

plt.savefig(outdir + "top%sgenes_per_Niche_NMF.pdf" %n_genes, bbox_inches = "tight")

#plot normalized n_facts and cell types
import pickle
# import models
def unpickle_model(path, mod_name):
    r""" Unpickle model
    """
    file = path + 'model_' + mod_name + ".p"
    mod1_ann = pickle.load(file = open(file, "rb"))
    return mod1_ann['mod']

n_fact = 6
#mod_path = "/CoLocatedComb/CoLocatedGroupsSklearnNMF_57787locations_37factors/models/"
mod_path = "/public/workspace/stu21230110/AMP/02nmf/visium_cytassist_integrated/cell2location_map_6_facts/CoLocatedComb/CoLocatedGroupsSklearnNMF_10535locations_11factors/models/"

mod_sk = unpickle_model(mod_path, f'n_fact{n_fact}')

filterval= 0.05
pretest_sk = mod_sk.cell_type_fractions[(mod_sk.cell_type_fractions['fact_0']>filterval)|(mod_sk.cell_type_fractions['fact_1']>filterval)|(mod_sk.cell_type_fractions['fact_2']>filterval)|(mod_sk.cell_type_fractions['fact_3']>filterval)|(mod_sk.cell_type_fractions['fact_4']>filterval)|(mod_sk.cell_type_fractions['fact_5']>filterval)]
pretest_sk

import matplotlib as mpl
from cell2location.plt.plot_heatmap import clustermap
b_dev_sel =  cts_nmf
fact_filt = ['fact_3','fact_2','fact_1','fact_5','fact_0','fact_4']
    
filterval= 0.05
test_cm = mod_sk.cell_type_fractions
test_sk = test_cm[(test_cm['fact_0']>filterval)|(test_cm['fact_1']>filterval)|(test_cm['fact_2']>filterval)|(test_cm['fact_3']>filterval)|(test_cm['fact_4']>filterval)|(test_cm['fact_5']>filterval)]
#test_sk = (test_sk.T / test_sk.max(1)).T

mpl.rc_file_defaults()
mpl.rcParams['pdf.fonttype'] = 42 # enables correct plotting of text
clustermap(test_sk.loc[b_dev_sel, fact_filt[::-1]],
           cluster_rows=True, cluster_cols=False, 
           figure_size=[5.9 + 0.12 * mod_sk.n_fact, 3 + 0.15 * mod_sk.n_var],
           fun_type='heatmap', array_size=None)

plt.savefig(outdir + "NMF_cytassist_visum_6facts_new_heatmap.pdf", bbox_inches='tight')

adata_vis.write("/public/workspace/stu21230110/AMP/01data/spatial.nmf.plots.h5ad")