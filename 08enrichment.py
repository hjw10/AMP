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
#import scvi
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

outdir = "/public/workspace/stu21230110/AMP/08enrichment/"

adata_vis = sc.read("/public/workspace/stu21230110/AMP/01data/spatial.nmf.niche.frequency.h5ad")

#all marker calculation
## All Markers
sc.tl.rank_genes_groups(adata_vis, groupby = 'Niche_NMF', groups = "all", use_raw = False, method = "wilcoxon",n_genes=2000, pts= True)
sc.pl.rank_genes_groups(adata_vis)
plt.savefig(outdir + "all_markers_calculation.pdf",bbox_inches='tight')

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

## Write to file 
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

sc.pl.matrixplot(adata_vis, var_names = genes, standard_scale = "var", groupby = "Niche_NMF", dendrogram = True,figsize=(16,2.5),layer="counts",)

plt.savefig(outdir + "top%sgenes_per_Niche_NMF.pdf" %n_genes,bbox_inches = "tight")


for i in adata_vis.obs.Niche_NMF.cat.categories:
    allMarkers_subset =allMarkers[allMarkers["cluster"]== i ]
    allMarkers_subset.to_csv(outdir + "AllMarkers_human_Niche_NMF_%s.csv" %i, sep = "\t", index = False)

marker3 = allMarkers
marker3 = marker3[marker3['pct_outside'] < 0.25]
marker3 = marker3[marker3['pct_within'] > 0.25]
marker3 = marker3[marker3['logfoldchange'] > 1]
marker3

#Code zum filtern hast ja schon, die folgenden Zeilen wählen dann die top n gene per cluster aus (im gleichen ordering wie das dendrogram)
marker3.sort_values(by = ["cluster", "logfoldchange"], ascending = [True, False], inplace = True)

#Get ordering as in Dendrogram
sc.tl.dendrogram(adata_vis, groupby = "Niche_NMF")
order = adata_vis.uns["dendrogram_Niche_NMF"]["dendrogram_info"]["ivl"]

#Get top genes from Markers Table
n_genes = 10
genes = []
for typ in order:
    curgenes = marker3.loc[marker3.cluster == typ, "gene"].values[0:n_genes]
    genes = genes + [g for g in curgenes if g not in genes]

sc.pl.matrixplot(adata_vis, var_names = genes, standard_scale = "var", groupby = "Niche_NMF", dendrogram = True,figsize=(20,3),layer="log1p")
plt.savefig(outdir + "Allmarkers_matrixplottop%sgenes_per_Niche_NMF_CM_025background_025celltype10logfc.pdf" %n_genes,bbox_inches = "tight")

#gene set enrichment on DGe genes between Niche_NMF_CM
allMarkers =allMarkers[allMarkers["pval_adj"]<0.05]
allMarkers
allMarkers.cluster.value_counts()


import gseapy
#Available databases : ‘Human’, ‘Mouse’, ‘Yeast’, ‘Fly’, ‘Fish’, ‘Worm’ 
gene_set_names = gseapy.get_library_name(organism='Human')
print(gene_set_names)

allMarkers_subset =allMarkers[allMarkers["cluster"]== "Niche1"]
allMarkers_subset =allMarkers_subset[allMarkers_subset["pval_adj"] < 0.05  ]
allMarkers_subset =allMarkers_subset[allMarkers_subset["logfoldchange"] > 2]
print(len(allMarkers_subset))

#gseapy.enrichr
niches = ['Niche1', 'Niche3', 'Niche4','Niche5','Niche6']
for i in niches:
    genelist = sc.get.rank_genes_groups_df(adata_vis, group=i, 
                                     log2fc_min=2, 
                                    pval_cutoff=0.05)['names'].squeeze().str.strip().tolist()
    enr_res = gseapy.enrichr(gene_list=genelist,
                     organism='Human',
                     background=adata_vis.var_names,
                     gene_sets='GO_Biological_Process_2021',
                     cutoff = 0.05,)
    enr_res.results["Niche"] = "%s" %i
    enr_res.results= enr_res.results[enr_res.results["Adjusted P-value"] <0.05]
    #enr_res.results=  enr_res.results.nsmallest(n=6, columns=['Adjusted P-value'],keep="all")
    enr_res.results.to_csv(outdir + "GO_enrichment/GO_enrichment_%s.csv" %i, sep = "\t", index = False)

gsea_res = pd.concat([pd.read_csv(outdir + "GO_enrichment/GO_enrichment_%s.csv" %i,sep = "\t")for i in niches], axis = 0)
gsea_res

# filter for only those top 6 terms per Niche, calculated in the prevíous filtered step
# allows to get also the non significant values from other Niches for shared categories to fill the empty spaces
#gsea_table_new = gsea_res[gsea_res["Term"].isin(gsea_table["Term"])] 
#gsea_table_new
top_8_per_niche = gsea_res.groupby('Niche').apply(lambda x: x.nlargest(8, 'Combined Score')).reset_index(drop=True)
top_8_per_niche = gsea_res.groupby('Niche').apply(lambda x: x.nsmallest(8, 'Adjusted P-value')).reset_index(drop=True)
gsea_table_new = top_8_per_niche

gsea_table_new = gsea_table_new.drop(["Gene_set","P-value","Old P-value","Old adjusted P-value","Odds Ratio","Genes"], axis=1)
gsea_table_new

gsea_table_new["Adjusted P-value"] = np.log10(gsea_table_new["Adjusted P-value"])
gsea_table_new

gsea_table_new["Adjusted P-value"] = gsea_table_new["Adjusted P-value"]*(-1)
gsea_table_new

gsea_table_new['Adjusted P-value'].values[gsea_table_new['Adjusted P-value'] > 30] = 30
gsea_table_new.sort_values(by=["Adjusted P-value"])
gsea_table_new

#gsea_table_new.replace([np.inf, -np.inf], 1200, inplace=True)
#gsea_table_new.sort_values(by=["Combined Score"])
#gsea_table_new


#gsea_table_new_cm= pd.pivot(gsea_table_new,index="Term", columns="Niche", values="Combined Score")
gsea_table_new_cm= pd.pivot(gsea_table_new,index="Term", columns="Niche", values="Adjusted P-value")
gsea_table_new_cm = gsea_table_new_cm.fillna(0)
gsea_table_new_cm

# 设置图形大小
import seaborn as sns
plt.figure(figsize=(10, 8))
# 绘制热图
sns.heatmap(gsea_table_new_cm, cmap=mymap, annot=False, linewidths=.5)

# 添加标题
plt.title('Top 8 Enriched Pathways per Niche')
plt.xlabel('Niche')
plt.ylabel('Term')

#plt.savefig(outdir + "enrich_heatmap_combined_score.pdf",bbox_inches = "tight")
plt.savefig(outdir + "GO_enrichment/enrich_heatmap_adj_Pval.pdf",bbox_inches = "tight")
gsea_table_new_cm.to_csv(outdir + "GO_enrichment/GO_enrichment_all.csv", sep = "\t", index = True)

#again but now for MSig
#'MSigDB_Hallmark_2020'
#?gseapy.enrichr
niches = ['Niche1','Niche3', 'Niche4', 'Niche5', 'Niche6']
for i in niches:
    genelist = sc.get.rank_genes_groups_df(adata_vis, group=i, 
                                     log2fc_min=0.5, 
                                    pval_cutoff=0.05)['names'].squeeze().str.strip().tolist()
    enr_res = gseapy.enrichr(gene_list=genelist,
                     organism='Human',
                     background=adata_vis.var_names,
                     gene_sets='MSigDB_Hallmark_2020',
                     cutoff = 0.05,)
    enr_res.results["Niche"] = "%s" %i
    enr_res.results= enr_res.results[enr_res.results["Adjusted P-value"] <0.05]
    #enr_res.results=  enr_res.results.nsmallest(n=5, columns=['Adjusted P-value'],keep="all")
    enr_res.results.to_csv(outdir + "/MSigDB_enrichment/MSigDB_Hallmark_2020_enrichment_%s.csv" %i, sep = "\t", index = False)


gsea_res = pd.concat([pd.read_csv(outdir + "/MSigDB_enrichment/MSigDB_Hallmark_2020_enrichment_%s.csv" %i,sep = "\t")for i in niches], axis = 0)
gsea_res

top_8_per_niche = gsea_res.groupby('Niche').apply(lambda x: x.nsmallest(8, 'Adjusted P-value')).reset_index(drop=True)
gsea_table_new = top_8_per_niche
gsea_table_new

gsea_table_new = gsea_table_new.drop(["Gene_set","P-value","Old P-value","Old adjusted P-value","Odds Ratio","Genes","Combined Score"], axis=1)
gsea_table_new

gsea_table_new["Adjusted P-value"] = np.log10(gsea_table_new["Adjusted P-value"])
gsea_table_new

gsea_table_new["Adjusted P-value"] = gsea_table_new["Adjusted P-value"]*(-1)
gsea_table_new

gsea_table_new['Adjusted P-value'].values[gsea_table_new['Adjusted P-value'] > 35] = 35
gsea_table_new.sort_values(by=["Adjusted P-value"])
gsea_table_new

gsea_table_new_cm= pd.pivot(gsea_table_new,index="Term", columns="Niche", values="Adjusted P-value")
gsea_table_new_cm = gsea_table_new_cm.fillna(0)
gsea_table_new_cm

# 设置图形大小
import seaborn as sns
plt.figure(figsize=(10, 8))
# 绘制热图
sns.heatmap(gsea_table_new_cm, cmap=mymap, annot=False, linewidths=.5)

# 添加标题
plt.title('Top 8 Enriched Pathways per Niche')
plt.xlabel('Niche')
plt.ylabel('Term')

plt.savefig(outdir + "MSigDB_enrichment/enrich_heatmap_adj_Pval.pdf",bbox_inches = "tight")
gsea_table_new_cm.to_csv(outdir + "MSigDB_enrichment/MSigDB_Hallmark_2020_enrichment_all.csv", sep = "\t", index = True)

