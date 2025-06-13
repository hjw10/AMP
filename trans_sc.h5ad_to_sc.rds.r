.libPaths("~/tools/seurat4/")
library(Seurat)
library(CellChat)
library(patchwork)
library(scater)
library(tidyverse)
library(cowplot)
library(Matrix.utils)
library(dplyr)
library(magrittr)
library(Matrix)
library(purrr)
library(reshape2)
library(S4Vectors)
library(tibble)
library(RColorBrewer)
library(openxlsx)
options(stringsAsFactors = FALSE)
library(NMF)
library(ggalluvial)


CellChatDB <- CellChatDB.human 

str(CellChatDB$geneInfo)

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

gene_key <- CellChatDB$geneInfo

outdir = "/public/workspace/stu21230110/SPA_result/12sc_cellchat"
#load data
library(sceasy)
reticulate::use_python("/public/workspace/stu21230110/anaconda3/bin/python",required = T)
sceasy::convertFormat("/public/workspace/stu21230110/SPA_result/myh5/scRNA.combined.frequence.h5ad", from="anndata", to="seurat", 
                      outFile="/public/workspace/stu21230110/SPA_result/myRDS/scRNA.combined.frequence.rds")