import os
import numpy as np
import stlearn as st
import pandas as pd
from matplotlib import pyplot as plt
import scanpy as sc
import random
from collections import Counter
import matplotlib
matplotlib.use('Agg')

indir = "/public/workspace/stu21230110/AMP/01data/"

samples = ["C1","C2","C3","P1","P2","P3"]

sample_name = samples[0]

outdir = "/public/workspace/stu21230110/AMP/10stLearn/" + sample_name + '/'

adata = sc.read(indir + sample_name + '.spatial.h5ad')
#samples = adata_vis.obs['sampleID'].unique()
#adatas = {sample: adata_vis[adata_vis.obs['sample'] == sample].copy() for sample in samples}

# Basic normalisation #
st.pp.filter_genes(adata, min_cells=3)
st.pp.normalize_total(adata) # NOTE: no log1p

adata = st.convert_scanpy(adata)

# Loading the LR databases available within stlearn (from NATMI)
lrs = st.tl.cci.load_lrs(['connectomeDB2020_lit'], species='human')
print(len(lrs))

# Running the analysis #
st.tl.cci.run(adata, lrs,
                  min_spots = 5, #Filter out any LR pairs with no scores for less than min_spots
                  distance=0, # None defaults to spot+immediate neighbours; distance=0 for within-spot mode
                  n_pairs=2000, # Number of random pairs to generate; low as example, recommend ~10,000
                  n_cpus=20, # Number of CPUs for parallel. If None, detects & use all available.
                  )

lr_info = adata.uns['lr_summary'] # A dataframe detailing the LR pairs ranked by number of significant spots.
lr_features = adata.uns['lrfeatures']
lr_info.to_csv(outdir + sample_name + '.lr_info.spatial.csv')
lr_features.to_csv(outdir + sample_name + '.lr_features.spatial.csv')

#uns_1 = adata.uns['rank_genes_groups_filtered']
#uns_2 = adata.uns['rank_genes_groups']
#adata.uns['rank_genes_groups_filtered']= {key: str(value) for key, value in uns_1.items()}
#adata.uns['rank_genes_groups']= {key: str(value) for key, value in uns_2.items()}
#adata.obs.columns = adata.obs.columns.astype(str)
#adata.var.columns = adata.var.columns.astype(str)
adata.uns['lrfeatures'] = np.array(adata.uns['lrfeatures']).astype(str)
#adata = sc.read(outdir + sample_name + '.spatial.h5ad')
#lr_features = pd.read_csv(outdir + sample_name + '.lr_features.spatial.csv')
#adata.uns['lrfeatures'] = lr_features
######adata.uns['lr_summary'] = np.array(adata.uns['lr_summary']).astype(str)


print('\n', lr_info)
adata.write(outdir + sample_name + '.stlearn.spatial.h5ad')
adata.uns['lrfeatures'] = lr_features

st.tl.cci.adj_pvals(adata, correct_axis='spot',pval_adj_cutoff=0.01, adj_method='fdr_bh')

st.pl.lr_summary(adata, n_top=500)
plt.savefig(outdir + sample_name + '.LR_Rank_500.pdf',dpi = 300)

st.pl.lr_summary(adata, n_top=50, figsize=(10,7))
plt.savefig(outdir + sample_name + '.LR_Rank_50.pdf',dpi = 300)

st.pl.lr_diagnostics(adata, figsize=(10,4.5))
plt.savefig(outdir + sample_name + '.Diagnostic.pdf',dpi = 300)

st.pl.lr_n_spots(adata, n_top=50, figsize=(11, 8),max_text=100)
plt.savefig(outdir + sample_name + '.LR.sig.top50.pdf',dpi = 300)

st.pl.lr_n_spots(adata, n_top=500, figsize=(11, 5),max_text=100)
plt.savefig(outdir + sample_name + '.LR.sig.top500.pdf',dpi = 300)



adata.uns['lr_summary']['interaction_score'] = sum(np.asarray(adata.obsm['lr_scores']))

lr_summary = pd.DataFrame(adata.uns['lr_summary'])
lr_summary.to_csv(outdir + sample_name + '.lr_summary.spatial.csv')

best_lr = pd.DataFrame(adata.uns['lr_summary']).sort_values(['n_spots_sig','interaction_score'],ascending=[False,False]).index.values[:20]

lr_interaction_score = pd.DataFrame(adata.obsm['lr_scores'],index = adata.obs.index,columns = adata.uns['lr_summary'].index)

lr_interaction_score.to_csv(outdir + sample_name + '.LR.interaction.score.csv')

lr_interaction_pval = pd.DataFrame(adata.obsm['p_adjs'],index = adata.obs.index,columns = adata.uns['lr_summary'].index)

lr_interaction_pval.to_csv(outdir + sample_name + '.LR.interaction.padj.csv')

adata_raw = adata.copy()


best_lr = ["PTPRM_PTPRM","MIF_CD74","MIF_CXCR4","APP_CD74"]

for lr in best_lr:
	stats = ['lr_scores','lr_sig_scores', 'p_vals', 'p_adjs', '-log10(p_adjs)']
	fig, axes = plt.subplots(ncols=len(stats), figsize=(20,10))
	for i, stat in enumerate(stats):
		st.pl.lr_result_plot(adata, use_result=stat, use_lr=lr, show_color_bar=True, ax=axes[i])
		axes[i].set_title(f'{lr} {stat}')
		plt.savefig(outdir + sample_name + '.LR.%s.pdf'%(lr),dpi = 300,bbox_inches = "tight")
	try:
		st.pl.lr_plot(adata, lr, inner_size_prop=1, outer_mode='binary', pt_scale=10,use_label=None,show_image=True,sig_spots=True)
		plt.savefig(outdir + sample_name + '.LR.%s.interaction.pdf'%(lr),dpi = 300,bbox_inches = "tight") ###The receptor is in green, the ligand is in red. The inner-point is the receptor, the outter point is the ligand.
		st.pl.lr_plot(adata,lr,inner_size_prop=0.04, middle_size_prop=.07, outer_size_prop=.4,outer_mode='continuous',
pt_scale=60,use_label=None, show_image=True,sig_spots=True)
		plt.savefig(outdir + sample_name + '.LR.%s.coexpression.pdf'%(lr),dpi = 300,bbox_inches = "tight")
	except:
		print("%s is not a lr pair ~~~"%(lr))

st.pl.gene_plot(adata, gene_symbols="CD74", contour=True,cell_alpha=0.5,figsize=(10, 6))
plt.savefig(outdir + sample_name + ".CD74.pdf",dpi = 300,bbox_inches = "tight")

st.pl.gene_plot(adata, gene_symbols="MIF", contour=True,cell_alpha=0.5,figsize=(10, 6))
plt.savefig(outdir + sample_name + ".MIF.pdf",dpi = 300,bbox_inches = "tight")

st.pl.gene_plot(adata, gene_symbols="PTPRM", contour=True,cell_alpha=0.5,figsize=(10, 6))
plt.savefig(outdir + sample_name + ".PTPRM.pdf",dpi = 300,bbox_inches = "tight")

st.pl.gene_plot(adata, gene_symbols="PTGS2", contour=True,cell_alpha=0.5,figsize=(10, 6))
plt.savefig(outdir + sample_name + ".PTGS2.pdf",dpi = 300,bbox_inches = "tight")

st.pl.gene_plot(adata, gene_symbols="CXCR4", contour=True,cell_alpha=0.5,figsize=(10, 6))
plt.savefig(outdir + sample_name + ".CXCR4.pdf",dpi = 300,bbox_inches = "tight")


senescence = ["MIF","CD74","CXCR4","CD44",
              "PTPRM",
              "APP"]
B_inflammation = ["CD74","AFF3","BANK1","HLA-DRB1","HLA-DPB1","HLA-DQA1","HLA-DRA","REL"]
inflammation = ["CD74","AFF3","BANK1","REL"]
B_anti_infla = ["MANF","PDX4"]
Epi_invasion = ["PIK3R3","SNHG12","SOX4"]

#st.pl.gene_plot(adata, gene_symbols=B_inflammation, method="CumSum")
st.pl.gene_plot(adata, gene_symbols=inflammation, method="CumSum")
plt.savefig(outdir + sample_name + ".inflammation.pdf",dpi = 300,bbox_inches = "tight")

st.pl.gene_plot(adata, gene_symbols=B_anti_infla, method="CumSum")
plt.savefig(outdir + sample_name + ".B.anti.inflammation.pdf",dpi = 300,bbox_inches = "tight")

st.pl.gene_plot(adata, gene_symbols=Epi_invasion, method="CumSum")
plt.savefig(outdir + sample_name + ".Epi.invasion.pdf",dpi = 300,bbox_inches = "tight")

st.pl.gene_plot(adata, gene_symbols=senescence, method="CumSum")
plt.savefig(outdir + sample_name + ".senescence.pdf",dpi = 300,bbox_inches = "tight")

#adata = sc.read(outdir + sample_name + '.spatial.correct.h5ad')

####Predicting significant CCIs
# Running the counting of co-occurence of cell types and LR expression hotspots #
st.tl.cci.run_cci(adata, 'Niche_NMF', # Spot cell information either in data.obs or data.uns
                  min_spots=3, # Minimum number of spots for LR to be tested.
                  spot_mixtures=True, # If True will use the label transfer scores,
                                      # so spots can have multiple cell types if score>cell_prop_cutoff
                  cell_prop_cutoff=0.1, # Spot considered to have cell type if score>0.1
                  sig_spots=True, # Only consider neighbourhoods of spots which had significant LR scores.
                  n_perms=1000 # Permutations of cell information to get background, recommend ~1000
                 )

#adata.uns['lrfeatures'] = np.array(adata.uns['lrfeatures']).astype(str)
#adata.write(outdir + sample_name + '.spatial.interaction.h5ad')
#adata = sc.read(outdir + sample_name + '.spatial.interaction.h5ad')

st.pl.cci_check(adata, 'Niche_NMF',figsize=(16, 14))

plt.savefig(outdir + sample_name + '.CCI.LR.interaction.pdf',dpi = 300,bbox_inches = "tight")
