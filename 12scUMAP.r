options(future.globals.maxSize= 80000*1024^2) 
library(data.table)
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)

setwd("/public/workspace/stu21230110/AMP/12scUMAP/")

#adata = readRDS("/public/workspace/stu21230110/AMP/01data/rds/scRNA.combined.frequence.rds")

load("/public/workspace/stu21230110/AMP/01data/rds/06_scRNA_harmony_singler202400810.Rdata")
adata <- scRNA_harmony_singler

###颜色
sample_color <- c("#F39B7FFF","#E377C2FF","#3C5488FF","#7E6148FF","#00A087FF","#8494FF","#4DBBD5FF", "#91D1C2FF","#9467BDFF","#E64B35FF","#B09C85FF")
pdf("sc.umap.celltype1.pdf")
DimPlot(adata, group.by="Annotation", label=T, label.size=3.5, reduction='umap',cols=sample_color)
dev.off()

pdf("sc.umap.celltype2.pdf")
DimPlot(adata, group.by="Annotation", label=F, label.size=3.5, reduction='umap',cols=sample_color)
dev.off()

