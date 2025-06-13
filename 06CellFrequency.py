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


outdir = "/public/workspace/stu21230110/AMP/06CellFrequency/"

adata_vis = sc.read("/public/workspace/stu21230110/AMP/01data/spatial.nmf.niche.frequency.h5ad")

adata = sc.read("/public/workspace/stu21230110/AMP/01data/sc.h5ad")

adata.obs.singleRnew.value_counts()
#singleRnew
#Stromal_cells               35444
#Epithelial_cells             7568
#T_cells                      6846
#Smooth_muscle_cells          4459
#Monocyte                     2096
#Endothelial_cells            2052
#Tissue_stem_cells            1293
#Myeloid_progenitor_cells      874
#NK_cells                      873
#Macrophage                    651
#B_cells                       638
#Name: count, dtype: int64

#within Niche per treatment
cts = ['B_cells', 'Endothelial_cells', 'Epithelial_cells', 'Macrophage',
       'Monocyte', 'Myeloid_progenitor_cells', 'NK_cells',
       'Smooth_muscle_cells', 'Stromal_cells', 'T_cells',
       'Tissue_stem_cells']
tab = sc.get.obs_df(adata_vis, keys= ["sampleID","Niche_NMF","treatment"]+ cts)

#scale cts data
ct_label = "cell_type"
group_by = "sampleID"
xlabel = "treatment"
ctss = cts
#cm = "sampleID"

tab = tab.loc[:, ctss + [group_by, xlabel]].copy()
tab.head()

## Calculate the mean per cell type and compartment
test_tab = tab.groupby([group_by, xlabel]).mean().reset_index()
print(test_tab.shape)
test_tab.head()

test_tab = test_tab.dropna()
test_tab

test_tab2 = test_tab.copy()
test_tab2

test_tab2 = test_tab2.drop(["sampleID","treatment"], axis=1)
test_tab2

test_tab2 = test_tab2.div(test_tab2.sum(axis=1), axis=0)
test_tab2

test_tab

test_tab_final = test_tab2.join(test_tab["sampleID"])
test_tab_final = test_tab_final.join(test_tab["treatment"])
test_tab_final

test_tab_final["identifier"] = test_tab_final["treatment"].astype(str) + "_spatial"
test_tab_final["identifier"] = test_tab_final["identifier"].str.capitalize()
test_tab_final

## relFreqs table from adata atlas
## Get Relative Frequencies per identifier
xlabel = "singleRnew"
cell_type_label = "orig.ident"

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


relFreqs = calc_relFreq_per_cluster(adata, group_by = cell_type_label, xlabel = xlabel)
relFreqs
#     B_cells  Endothelial_cells  Epithelial_cells  Macrophage  Monocyte  ...  Smooth_muscle_cells  Stromal_cells   T_cells  Tissue_stem_cells  orig.ident
#C1  0.006006           0.101154          0.048163    0.001766  0.014013  ...             0.146020       0.477037  0.184056           0.006712          C1
#C2  0.041014           0.019484          0.098043    0.032918  0.075890  ...             0.049733       0.490302  0.141103           0.025445          C2
#C3  0.001130           0.018174          0.269391    0.003565  0.023043  ...             0.054174       0.545304  0.033913           0.031739          C3
#P1  0.002627           0.017240          0.023315    0.005665  0.030047  ...             0.072818       0.798292  0.026188           0.008210          P1
#P2  0.007716           0.029745          0.249035    0.014188  0.041444  ...             0.055881       0.228998  0.274300           0.025140          P2
#P3  0.001675           0.027851          0.059404    0.003702  0.014102  ...             0.061784       0.709237  0.069099           0.024943          P3

relFreqs["sampleID"] = relFreqs.index
def get_sample_group(row_name):
    if row_name.startswith('C'):
        return 'Control'
    else:
        return 'Pain'

# 应用分组规则
relFreqs['treatment'] = relFreqs.index.map(get_sample_group)
relFreqs["identifier"] = relFreqs["treatment"].astype(str) + "_scRNA"
relFreqs

relFreqs2 =  relFreqs.drop(["orig.ident"], axis=1)
relFreqs2.reset_index(drop=True, inplace=True)
relFreqs2

combined = pd.concat([test_tab_final, relFreqs2], axis = 0)
combined

combined = combined.drop(["treatment"], axis=1)
combined

ct_label = "cell_type"
group_by = "sampleID"
xlabel = "identifier"
ctss = cts
cm = "sampleID"

## in order to plot all cell types, melt the data frame
tab_bleo = pd.melt(combined, id_vars = [group_by, xlabel] , var_name = ct_label)
print(tab_bleo.shape)
tab_bleo.head()

fig, ax = plt.subplots()
fig.set_size_inches(9, 5)
#order = ["Alveolar_AT1",]# "Alveolar_AT2", "Alveolar_Macrophage", "Fibrotic",]

fig = sb.boxplot(data = tab_bleo, hue = "identifier", x = "cell_type", y = "value", order = None,showfliers=False) 
fig.legend(bbox_to_anchor= (1.01, 1), loc='upper left', borderaxespad=0)
plt.xlabel("")
plt.ylim(0,0.1)
plt.ylabel("cell type frequency")
plt.title("Comparision of cell type frequencies between spatial samples and scRNA samples")
fig.set_xticklabels(ax.get_xticklabels(),rotation=90)
plt.savefig(outdir + 'Comparision_celltypefrequencies_scRNA_spatial_cut.pdf',bbox_inches = "tight")

#0-0.3
fig, ax = plt.subplots()
fig.set_size_inches(9, 5)
#order = ["Alveolar_AT1",]# "Alveolar_AT2", "Alveolar_Macrophage", "Fibrotic",]

fig = sb.boxplot(data = tab_bleo, hue = "identifier", x = "cell_type", y = "value", order = None,showfliers=False) 
fig.legend(bbox_to_anchor= (1.01, 1), loc='upper left', borderaxespad=0)
plt.xlabel("")
plt.ylim(0,0.3)
plt.ylabel("cell type frequency")
plt.title("Comparision of cell type frequencies between spatial samples and scRNA samples")
fig.set_xticklabels(ax.get_xticklabels(),rotation=90)
plt.savefig(outdir + 'Comparision_celltypefrequencies_scRNA_spatial_cut0-0.3.pdf',bbox_inches = "tight")







fig, ax = plt.subplots()
fig.set_size_inches(10,5)
#order = ["Alveolar_AT1",]# "Alveolar_AT2", "Alveolar_Macrophage", "Fibrotic",]

fig = sb.boxplot(data = tab_bleo, hue = "identifier", x = "cell_type", y = "value", order = None,showfliers=False) #hue_order = ["AT1", "AT2", "Macrophage FABP4+", "Myofibroblast","Macrophage"])
fig.legend(bbox_to_anchor= (1.01, 1), loc='upper left', borderaxespad=0)
plt.xlabel("")
#plt.ylim(0,0.4)
plt.ylabel("cell type frequency")
plt.title("Comparision of cell type frequencies between spatial samples and scRNA samples")
fig.set_xticklabels(ax.get_xticklabels(),rotation=90)
plt.savefig(outdir + 'Comparision_celltypefrequencies_scRNA_spatial.pdf',bbox_inches = "tight")

import itertools
a = cts
b = np.unique(combined['identifier'].to_numpy())

healthy_pairs = [("Control_spatial", "Control_scRNA")]
pain_pairs = [("Pain_spatial", "Pain_scRNA")]

# 生成最终的 pairs2
pairs2 = []
pairs3 = []
for cell in a:
    # 对每个细胞成分，生成 Healthy 和 pain 的组合
    for healthy_pair in healthy_pairs:
        pairs2.append(((cell, healthy_pair[0]), (cell, healthy_pair[1])))

for cell in a:
    for pain_pair in pain_pairs:
        pairs3.append(((cell, pain_pair[0]), (cell, pain_pair[1])))

pairs4 = pairs2 + pairs3


from statannotations.Annotator import Annotator


fig, ax = plt.subplots()
fig.set_size_inches(10, 5)
#order = ["Alveolar_AT1",]# "Alveolar_AT2", "Alveolar_Macrophage", "Fibrotic",]

fig = sb.boxplot(data = tab_bleo, hue = "identifier", x = "cell_type", y = "value", hue_order = ["Control_spatial","Control_scRNA","Pain_spatial","Pain_scRNA"],showfliers=False) #hue_order = ["AT1", "AT2", "Macrophage FABP4+", "Myofibroblast","Macrophage"])
fig.legend(bbox_to_anchor= (1.01, 1), loc='upper left', borderaxespad=0)
plt.xlabel("")
plt.ylabel("frequency")
plt.title("Comparision of cell type frequencies between spatial samples and scRNA samples")
fig.set_xticklabels(ax.get_xticklabels(),rotation=90)
annot = Annotator(fig, pairs4, data = tab_bleo,hue = "identifier", x = "cell_type", y = "value", hide_non_significant=True)
annot.configure(test='Mann-Whitney', verbose=0) #'t-test_ind' 'Mann-Whitney'
#annot.configure(test='Mann-Whitney', text_format='star', loc='outside')
annot.apply_test()
annot.annotate()
plt.savefig(outdir + 'Comparision_celltypefrequencies_scRNA_spatial_sig.pdf',bbox_inches = "tight")
