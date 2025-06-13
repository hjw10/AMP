#devtools::install_github("sqjin/CellChat")
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

#load data
#library(sceasy)
#reticulate::use_python("/public/workspace/stu21230110/anaconda3/bin/python",required = T)
#sceasy::convertFormat("/public/workspace/stu21230110/SPA_result/myh5/scRNA.combined.frequence.h5ad", from="anndata", to="seurat", 
#                      outFile="/public/workspace/stu21230110/SPA_result/myRDS/scRNA.combined.frequence.rds")

#adata <- readRDS("/public/workspace/stu21230110/AMP/01data/rds/scRNA.combined.frequence.rds")
load("/public/workspace/stu21230110/AMP/01data/rds/06_scRNA_harmony_singler202400810.Rdata")
adata <- scRNA_harmony_singler

gene_metadata = adata[["RNA"]]@meta.features
gene_metadata["ensembl_id"] = rownames(gene_metadata)
rownames(gene_metadata) <- gene_metadata$gene_name
head(gene_metadata)

table(adata@meta.data$group)

counts = as.matrix(GetAssayData(adata, assay = "RNA", slot = "data"))
metadata = adata@meta.data
features = gene_metadata
features_orig = gene_metadata

gene_key <- CellChatDB$geneInfo
gene_key <- gene_key %>%
  select(Symbol, `Ensembl.Gene.ID`) %>%
  filter(grepl("^ENS", `Ensembl.Gene.ID`)) %>%
  mutate(ensembl_id = `Ensembl.Gene.ID`) 
rownames(gene_key) = gene_key$ensembl_id
gene_key <- gene_key %>%
  select(Symbol,ensembl_id) %>%
  distinct()
head(gene_key)

setwd("/public/workspace/stu21230110/AMP/13cellchat/control")

# 创建cellchat对象
seurat=subset(adata,group=="C")#读取seurat对象
data.input = as.matrix(GetAssayData(seurat, assay = "RNA", slot = "data"))#得到标准化表达矩阵
meta = seurat@meta.data # 获得实验信息
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "singleRnew")#创建cellchat对象

levels(cellchat@idents) # 显示细胞类型
groupSize <- as.numeric(table(cellchat@idents)) # 每种细胞类型的细胞数目统计
groupSize

CellChatDB <- CellChatDB.human # 若是小鼠，请用 CellChatDB.mouse 数据库
#showDatabaseCategory(CellChatDB)#显示数据库中的配受体分类

cellchat@idents = droplevels(cellchat@idents, exclude = setdiff(levels(cellchat@idents),unique(cellchat@idents)))#丢弃未使用的因子水平，在我们这个例子数据中需要该步骤


cellchat@DB <- CellChatDB #数据库赋值
cellchat <- subsetData(cellchat) # 必须，即时是整个数据集
cellchat <- identifyOverExpressedGenes(cellchat)#鉴定过表达基因
cellchat <- identifyOverExpressedInteractions(cellchat)#鉴定过表达互作
cellchat <- projectData(cellchat, PPI.human)#将基因表达数据映射PPI互作网络

#################################Inference of cell-cell communication network#####################################
cellchat <- computeCommunProb(cellchat)#计算互作概率
cellchat <- filterCommunication(cellchat, min.cells = 1)#过滤掉低于min.cells的细胞类型
cellchat <- computeCommunProbPathway(cellchat)#计算pathway水平的互作概率，通过总结所有的配受体对
cellchat <- aggregateNet(cellchat)#得到整体细胞互作网络

#cellchat <- readRDS("cellchat_control.rds")
saveRDS(cellchat,file="cellchat_control.rds")
cellchat@netP$pathways#显示显著的pathways
######一定要看结果

[1] "COLLAGEN"   "MK"         "LAMININ"    "MIF"        "APP"
 [6] "MHC-I"      "FN1"        "MHC-II"     "CD99"       "PARs"
[11] "CLEC"       "PTN"        "CADM"       "PTPRM"      "GALECTIN"
[16] "TENASCIN"   "VISFATIN"   "THY1"       "CCL"        "NRXN"
[21] "GAS"        "CD45"       "SEMA3"      "SPP1"       "CD46"
[26] "ITGB2"      "ICAM"       "IGF"        "TGFb"       "JAM"
[31] "NOTCH"      "PECAM1"     "THBS"       "VEGF"       "TNF"
[36] "ncWNT"      "EDN"        "CDH"        "BMP"        "EPHA"
[41] "ADGRE5"     "HH"         "SEMA4"      "ANNEXIN"    "PDGF"
[46] "COMPLEMENT" "ALCAM"      "CD6"        "ANGPT"      "CDH5"
[51] "CDH1"       "CXCL"       "ESAM"       "NECTIN"     "HSPG"
[56] "ANGPTL"     "RESISTIN"   "MPZ"        "CD40"       "CHEMERIN"
[61] "IFN-II"     "EGF"        "BAFF"       "WNT"        "APRIL"
[66] "EPHB"       "PROS"       "NRG"        "IL16"       "SEMA6"
[71] "CNTN"       "OCLN"





pdf("01-1count_weight.pdf")
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")#互作数目
dev.off()

pdf("01-2count_weight.pdf")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")#互作强度
dev.off()







mat <- cellchat@net$weight#互作强度值
pdf("02allcell_count_weight.pdf")
par(mfrow = c(3,4), xpd=TRUE) # 设定绘图区，3行4列
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))#构建1个0值矩阵
  mat2[i, ] <- mat[i, ]#替换掉矩阵对应行
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])#网络图
}
dev.off()
# 层次图
levels(cellchat@idents) 

 [1] "B_cells"                  "Endothelial_cells"
 [3] "Epithelial_cells"         "Macrophage"
 [5] "Monocyte"                 "Myeloid_progenitor_cells"
 [7] "NK_cells"                 "Smooth_muscle_cells"
 [9] "Stromal_cells"            "T_cells"
[11] "Tissue_stem_cells"

pathways.show <- "CCL"#指定信号通路
vertex.receiver = 6 # 指定目标细胞类型
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver,layout = "circle")#绘图
# 所有显著pathways
pathways.show.all <- cellchat@netP$pathways

for (i in 1:length(pathways.show.all)) {
  # 信号pathway以及对应的每个配受体对可视化，网络图，弦图，热图
  netVisual(cellchat, signaling = pathways.show.all[i], layout = "circle",out.format=c("pdf"))#网络图
  netVisual(cellchat, signaling = pathways.show.all[i], layout = "chord",out.format=c("pdf"),height=12)#弦图
  pdf(file=paste0("03-1",pathways.show.all[i], "_heatmap.pdf"))#热图
  print(netVisual_heatmap(cellchat, signaling = pathways.show.all[i], color.heatmap = "Reds"))
  dev.off()
  # 计算每个配受体对对信号通路的贡献大小
  print("LR contribution!")
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])#绘图
  ggsave(filename=paste0("03-2",pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 8, height = 8, units = 'in', dpi = 300)#保存图片
}

#对感兴趣的细胞类型，气泡图显示所有显著LR对
pdf("03-0netVisual_bubble_B_cells.pdf",width=10,height=20)
netVisual_bubble(cellchat, sources.use = c(1:5,7:11), targets.use = 1, remove.isolate = FALSE)
dev.off()
pdf("03-0netVisual_bubble_Epithelial_cells.pdf",width=10,height=20)
netVisual_bubble(cellchat, sources.use = c(1:5,7:11), targets.use = 3, remove.isolate = FALSE)
dev.off()
pdf("03-0netVisual_bubble_Tissue_stem_cells.pdf",width=10,height=20)
netVisual_bubble(cellchat, sources.use = c(1:5,7:11), targets.use = 11, remove.isolate = FALSE)
dev.off()
#绘制感兴趣信号通路相关配受体对的基因表达图

pdf("03-1-plotGeneExpression.pdf",width=8,height=10)
plotGeneExpression(cellchat, signaling = "COLLAGEN")
dev.off()


# 网络分析
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # 计算信号通路的网络中心得分
# 热图可视化细胞所扮演的主要角色，亦可用散点图来表示
for(i in 1:length(pathways.show.all)){
  pdf(file=paste0("04-1",pathways.show.all[i], "_major_signaling_roles_heatmap.pdf"),width=10,height=8)#热图
  print(netAnalysis_signalingRole_network(cellchat, signaling = pathways.show.all[i], width = 8, height = 2.5, font.size = 10))
  dev.off()
  pdf(file=paste0("04-2",pathways.show.all[i], "_major_signaling_roles_scatter.pdf"),width=8,height=8)#散点图
  print(netAnalysis_signalingRole_scatter(cellchat,signaling = pathways.show.all[i]))
  dev.off()
}

# 所有信号通路汇总分析，看哪个信号通路贡献最大
pdf("05netAnalysis_signalingRole_heatmap.pdf",height=32,width=18)
par(mfrow = c(1,2))#
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing",width=10,height=30)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming",width=10,height=30)
ht1 + ht2
dev.off()

#保存结果
df.net <- subsetCommunication(cellchat)#得到细胞互作数据框
write.csv(df.net,"06-net_lr.csv")#输出文件
df.netp <- subsetCommunication(cellchat,slot.name = "netP")#得到信号通路水平的细胞互作结果
write.csv(df.netp,"07-net_pathway.csv")#输出文件
cellchat_control<-cellchat
saveRDS(cellchat_control, file = "control_cellchat_single_dataset.rds")#保存cellchat对象为rds文件




##########################################################################################11 跑重症组
#cellchat  20221026重新跑screen -S pd11
#单数据集
#load("01RDada/06_adata_singler202400810.Rdata")
####
#adata<-adata_singler
library(Seurat)#载入R包
library(ggplot2)
library(CellChat)
library(patchwork)
library(NMF)


setwd("/public/workspace/stu21230110/AMP/13cellchat/pain")
# 创建cellchat对象
seurat=subset(adata,group=="P")#读取seurat对象
data.input = as.matrix(GetAssayData(seurat, assay = "RNA", slot = "data"))#得到标准化表达矩阵
meta = seurat@meta.data # 获得实验信息
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "singleRnew")#创建cellchat对象

levels(cellchat@idents) # 显示细胞类型
groupSize <- as.numeric(table(cellchat@idents)) # 每种细胞类型的细胞数目统计
groupSize

CellChatDB <- CellChatDB.human # 若是小鼠，请用 CellChatDB.mouse 数据库
#showDatabaseCategory(CellChatDB)#显示数据库中的配受体分类

cellchat@idents = droplevels(cellchat@idents, exclude = setdiff(levels(cellchat@idents),unique(cellchat@idents)))#丢弃未使用的因子水平，在我们这个例子数据中需要该步骤


cellchat@DB <- CellChatDB #数据库赋值
cellchat <- subsetData(cellchat) # 必须，即时是整个数据集
cellchat <- identifyOverExpressedGenes(cellchat)#鉴定过表达基因
cellchat <- identifyOverExpressedInteractions(cellchat)#鉴定过表达互作
cellchat <- projectData(cellchat, PPI.human)#将基因表达数据映射PPI互作网络

#################################Inference of cell-cell communication network#####################################
cellchat <- computeCommunProb(cellchat)#计算互作概率
cellchat <- filterCommunication(cellchat, min.cells = 1)#过滤掉低于min.cells的细胞类型
cellchat <- computeCommunProbPathway(cellchat)#计算pathway水平的互作概率，通过总结所有的配受体对
cellchat <- aggregateNet(cellchat)#得到整体细胞互作网络
saveRDS(cellchat,file="cellchat_P.rds")
cellchat@netP$pathways#显示显著的pathways
######一定要看结果

 [1] "COLLAGEN"   "MIF"        "LAMININ"    "MK"         "APP"
 [6] "MHC-I"      "FN1"        "CD99"       "MHC-II"     "PARs"
[11] "VISFATIN"   "CLEC"       "PTPRM"      "ADGRE5"     "CD45"
[16] "PTN"        "TENASCIN"   "GALECTIN"   "ITGB2"      "THY1"
[21] "CCL"        "ICAM"       "PECAM1"     "NRXN"       "IGF"
[26] "VEGF"       "ANNEXIN"    "CXCL"       "THBS"       "TGFb"
[31] "SPP1"       "CD22"       "JAM"        "EPHA"       "CD46"
[36] "CDH"        "NOTCH"      "ncWNT"      "CADM"       "GAS"
[41] "EGF"        "BMP"        "PDGF"       "IFN-II"     "NECTIN"
[46] "SEMA3"      "ANGPTL"     "CDH1"       "EDN"        "ANGPT"
[51] "ALCAM"      "CD6"        "CDH5"       "SEMA4"      "CD226"
[56] "ESAM"       "SELL"       "COMPLEMENT" "TIGIT"      "GRN"
[61] "NCAM"       "PROS"       "MPZ"        "EPHB"       "WNT"
[66] "FGF"        "HGF"        "BAFF"       "CHEMERIN"   "NRG"
[71] "SEMA6"      "SELPLG"     "SN"         "OCLN"







pdf("01-1count_weight.pdf")
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")#互作数目
dev.off()

pdf("01-2count_weight.pdf")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")#互作强度
dev.off()

mat <- cellchat@net$weight#互作强度值
pdf("02allcell_count_weight.pdf")
par(mfrow = c(3,4), xpd=TRUE) # 设定绘图区，3行4列
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))#构建1个0值矩阵
  mat2[i, ] <- mat[i, ]#替换掉矩阵对应行
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])#网络图
}
dev.off()
# 层次图
levels(cellchat@idents) 

 [1] "B_cells"                  "Endothelial_cells"
 [3] "Epithelial_cells"         "Macrophage"
 [5] "Monocyte"                 "Myeloid_progenitor_cells"
 [7] "NK_cells"                 "Smooth_muscle_cells"
 [9] "Stromal_cells"            "T_cells"
[11] "Tissue_stem_cells"

pathways.show <- "CCL"#指定信号通路
vertex.receiver = 6 # 指定目标细胞类型
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver,layout = "circle")#绘图
# 所有显著pathways
pathways.show.all <- cellchat@netP$pathways

for (i in 1:length(pathways.show.all)) {
  # 信号pathway以及对应的每个配受体对可视化，网络图，弦图，热图
  netVisual(cellchat, signaling = pathways.show.all[i], layout = "circle",out.format=c("pdf"))#网络图
  netVisual(cellchat, signaling = pathways.show.all[i], layout = "chord",out.format=c("pdf"),height=12)#弦图
  pdf(file=paste0("03-1",pathways.show.all[i], "_heatmap.pdf"))#热图
  print(netVisual_heatmap(cellchat, signaling = pathways.show.all[i], color.heatmap = "Reds"))
  dev.off()
  # 计算每个配受体对对信号通路的贡献大小
  print("LR contribution!")
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])#绘图
  ggsave(filename=paste0("03-2",pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 8, height = 8, units = 'in', dpi = 300)#保存图片
}

#对感兴趣的细胞类型，气泡图显示所有显著LR对
pdf("03-0netVisual_bubble_B_cells.pdf",width=10,height=20)
netVisual_bubble(cellchat, sources.use = c(1:5,7:11), targets.use = 1, remove.isolate = FALSE)
dev.off()
pdf("03-0netVisual_bubble_Epithelial_cells.pdf",width=10,height=20)
netVisual_bubble(cellchat, sources.use = c(1:5,7:11), targets.use = 3, remove.isolate = FALSE)
dev.off()
pdf("03-0netVisual_bubble_Tissue_stem_cells.pdf",width=10,height=20)
netVisual_bubble(cellchat, sources.use = c(1:5,7:11), targets.use = 11, remove.isolate = FALSE)
dev.off()
#绘制感兴趣信号通路相关配受体对的基因表达图

pdf("03-1-plotGeneExpression.pdf",width=8,height=10)
plotGeneExpression(cellchat, signaling = "SN")
dev.off()


# 网络分析
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # 计算信号通路的网络中心得分
# 热图可视化细胞所扮演的主要角色，亦可用散点图来表示
for(i in 1:length(pathways.show.all)){
  pdf(file=paste0("04-1",pathways.show.all[i], "_major_signaling_roles_heatmap.pdf"),width=10,height=8)#热图
  print(netAnalysis_signalingRole_network(cellchat, signaling = pathways.show.all[i], width = 8, height = 2.5, font.size = 10))
  dev.off()
  pdf(file=paste0("04-2",pathways.show.all[i], "_major_signaling_roles_scatter.pdf"),width=8,height=8)#散点图
  print(netAnalysis_signalingRole_scatter(cellchat,signaling = pathways.show.all[i]))
  dev.off()
}

# 所有信号通路汇总分析，看哪个信号通路贡献最大
pdf("05netAnalysis_signalingRole_heatmap.pdf",height=32,width=18)
par(mfrow = c(1,2))#
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing",width=10,height=30)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming",width=10,height=30)
ht1 + ht2
dev.off()
#保存结果
df.net <- subsetCommunication(cellchat)#得到细胞互作数据框
write.csv(df.net,"06-net_lr.csv")#输出文件
df.netp <- subsetCommunication(cellchat,slot.name = "netP")#得到信号通路水平的细胞互作结果
write.csv(df.netp,"07-net_pathway.csv")#输出文件
cellchat_pain<-cellchat
saveRDS(cellchat_pain, file = "pain_cellchat_single_dataset.rds")#保存cellchat对象为rds文件




#########################################################12cellchat多数据集
#拥有相同细胞构成的多个数据集比较分析
library(CellChat)#载入R包
library(patchwork)



setwd("/public/workspace/stu21230110/AMP/13cellchat/diff")#切换工作目录

cellchat.C <- readRDS("../control/control_cellchat_single_dataset.rds")#读取ctrl的cellchat对象
cellchat.P <- readRDS("../pain/pain_cellchat_single_dataset.rds")#读取treat的cellchat对象
object.list <- list(C= cellchat.C, P= cellchat.P)#创建列表
cellchat <- mergeCellChat(object.list, add.names = names(object.list))#合并cellchat对象
#cellchat <- readRDS("cellchat_comparisonAnalysis_P_vs_C.rds")

pdf(file="01num_weight.pdf",width=16,height=8)
#互作数目和强度比较
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))#互作数目
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")#互作强度
gg1 + gg2
dev.off()
#circle plot，网络图
pdf(file="02circle_plot.pdf",width=16,height=8)
par(mfrow = c(1,2), xpd=TRUE)#设置绘图区域,1行2列
netVisual_diffInteraction(cellchat, weight.scale = T)#互作数目
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")#互作强度
dev.off()
#heatmap，热图
pdf(file="03heatmap.pdf",width=16,height=8)
gg1 <- netVisual_heatmap(cellchat)#互作数目
gg2 <- netVisual_heatmap(cellchat, measure = "weight")#互作强度
gg1 + gg2
dev.off()

#上述差异分析只适合两组比较，如果超过两组，可以采用下面的方式来比较分析。
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))#获得多个数据集间最大的边权重
#par(mfrow = c(1,2), xpd=TRUE)#设置绘图区域，1行2列
for (i in 1:length(object.list)) {#依次绘制每个数据对象的互作网络图
  p<-netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
  #ggsave(filename=paste0("04-",i, "_network.pdf"), plot=p, width = 8, height = 8, units = 'in', dpi = 300)#保存图片
  pdf(file=paste0("04-",i, "_network.pdf"), width = 8, height = 8)
  p
  dev.off()
}

#发送（outgoing）和接收（incoming）信号角色比较
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})#计算互作数目
weight.MinMax <- c(min(num.link), max(num.link)) # 获取互作数目的最小值和最大值，用来标准化图形的点大小
#gg <- list()#新建一个空列表
for (i in 1:length(object.list)) {#绘制气泡图
  p <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
  ggsave(filename=paste0("05-",i, "_dotplot.pdf"), plot=p, width = 8, height = 8, units = 'in', dpi = 300)#保存图片

}
#patchwork::wrap_plots(plots = gg)
#对感兴趣细胞类型的信号通路变化进行分析
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "B_cells")
pdf("B_cells.change.pdf")
gg1
dev.off()

gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Epithelial_cells")
pdf("Epithelial_cells.change.pdf")
gg2
dev.off()

gg3 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Tissue_stem_cells")
pdf("Tissue_stem_cells.change.pdf")
gg3
dev.off()

#两组信号通路差异分析
pdf(file="06deg.pdf",width=16,height=8)
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2
dev.off()
#发送（outgoing）和接收（incoming）信号通路比较
library(ComplexHeatmap)
i = 1
#合并所有数据集鉴定的信号通路 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 30, height = 25)#第一个对象热图
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 30, height = 25)#第二个对象热图
pdf(file="07two.pdf",width=40,height=26)
draw(ht1 + ht2, ht_gap = unit(1, "cm"))#合并两个绘图对象
dev.off()
#配受体对差异分析
#netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:10),  comparison = c(1, 2), angle.x = 45)
#分别显示互作增加或减少
target <- c(1,3,11)
pdf(file="08interaction_deg_B_cells.pdf",width=16,height=12)
gg1 <- netVisual_bubble(cellchat, sources.use = c(2:11), targets.use = 1,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in patient", angle.x = 45, remove.isolate = T)#互作增加
gg2 <- netVisual_bubble(cellchat, sources.use = c(2:11), targets.use = 1,  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in patient", angle.x = 45, remove.isolate = T)#互作减少
gg1 + gg2
dev.off()
pdf(file="08interaction_deg_Epithelial_cells.pdf",width=16,height=12)
gg1 <- netVisual_bubble(cellchat, sources.use = c(1:2,4:11), targets.use = 3,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in patient", angle.x = 45, remove.isolate = T)#互作增加
gg2 <- netVisual_bubble(cellchat, sources.use = c(1:2,4:11), targets.use = 3,  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in patient", angle.x = 45, remove.isolate = T)#互作减少
gg1 + gg2
dev.off()
pdf(file="08interaction_deg_Tissue_stem_cells.pdf",width=16,height=12)
gg1 <- netVisual_bubble(cellchat, sources.use = c(1:10), targets.use = 11,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in patient", angle.x = 45, remove.isolate = T)#互作增加
gg2 <- netVisual_bubble(cellchat, sources.use = c(1:10), targets.use = 11,  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in patient", angle.x = 45, remove.isolate = T)#互作减少
gg1 + gg2
dev.off()

pdf(file="08-2interaction_deg_B_cells.pdf",width=16,height=12)
gg1 <- netVisual_bubble(cellchat, sources.use = 1, targets.use = c(2:11),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in patient", angle.x = 45, remove.isolate = T)#互作增加
gg2 <- netVisual_bubble(cellchat, sources.use = 1, targets.use = c(2:11),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in patient", angle.x = 45, remove.isolate = T)#互作减少
gg1 + gg2
dev.off()
pdf(file="08-2interaction_deg_Epithelial_cells.pdf",width=16,height=12)
gg1 <- netVisual_bubble(cellchat, sources.use = 3, targets.use = c(1:2,4:11),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in patient", angle.x = 45, remove.isolate = T)#互作增加
gg2 <- netVisual_bubble(cellchat, sources.use = 3, targets.use = c(1:2,4:11),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in patient", angle.x = 45, remove.isolate = T)#互作减少
gg1 + gg2
dev.off()
pdf(file="08-2interaction_deg_Tissue_stem_cells.pdf",width=16,height=12)
gg1 <- netVisual_bubble(cellchat, sources.use = 11, targets.use = c(1:10),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in patient", angle.x = 45, remove.isolate = T)#互作增加
gg2 <- netVisual_bubble(cellchat, sources.use = 11, targets.use = c(1:10),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in patient", angle.x = 45, remove.isolate = T)#互作减少
gg1 + gg2
dev.off()

#c("B_cells","Epithelial_cells","Tissue_stem_cells")
pdf(file="08-3interaction_deg_B_Epi_Tissue.pdf",width=16,height=12)
gg1 <- netVisual_bubble(cellchat, 
                        sources.use = c("B_cells","Epithelial_cells","Tissue_stem_cells"), 
                        targets.use = c("B_cells","Epithelial_cells","Tissue_stem_cells"),  
                        comparison = c(1, 2), 
                        max.dataset = 2, 
                        title.name = "Increased signaling in patient", 
                        angle.x = 45, 
                        remove.isolate = T)#互作增加
gg2 <- netVisual_bubble(cellchat, 
                        sources.use = c("B_cells","Epithelial_cells","Tissue_stem_cells"), 
                        targets.use = c("B_cells","Epithelial_cells","Tissue_stem_cells"),  
                        comparison = c(1, 2), 
                        max.dataset = 1, 
                        title.name = "Decreased signaling in patient", 
                        angle.x = 45,
                        remove.isolate = T)#互作减少
gg1 + gg2
dev.off()




gg1$data#得到配受体对数据
gg2$data
#画基因表达
cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("C", "P")) # 设置因子水平
pdf(file="09gene.MIF.pdf",width=16,height=8)
plotGeneExpression(cellchat, signaling = "MIF", split.by = "datasets", colors.ggplot = T)#画基因表达图
dev.off()
#保存rds文件
saveRDS(cellchat, file = "cellchat_comparisonAnalysis_P_vs_C.rds")







setwd("/public/workspace/stu21230110/AMP/13cellchat/niche4")

Control <- readRDS("../control/control_cellchat_single_dataset.rds")#读取ctrl的cellchat对象
Pain <- readRDS("../pain/pain_cellchat_single_dataset.rds")#读取treat的cellchat对象

groupSize_Control <- as.numeric(table(Control@idents))

# set side-by-side plot displays 

pdf("Number.of.interactions.Control.pdf")
p <- netVisual_circle(
    Control@net$count, 
    #vertex.weight = groupSize_Control, 
    weight.scale = TRUE, 
    label.edge= F, 
    title.name = "Number of interactions - Control",
    sources.use = c("B_cells","Epithelial_cells","Tissue_stem_cells"),
    targets.use = c("B_cells","Epithelial_cells","Tissue_stem_cells"),
    remove.isolate=T
)
p
dev.off()

groupSize_Pain <- as.numeric(table(Pain@idents))

pdf("Number.of.interactions.Pain.pdf")
p <- netVisual_circle(
    Pain@net$count, 
    #vertex.weight = groupSize_Pain, 
    weight.scale = TRUE, 
    label.edge= FALSE, 
    title.name = "Number of interactions - Pain",
    sources.use = c("B_cells","Epithelial_cells","Tissue_stem_cells"),
    targets.use = c("B_cells","Epithelial_cells","Tissue_stem_cells"),
    remove.isolate=T
)
p
dev.off()


pdf(file = "CellChat_Pain_Niche4_Circle_numbered.pdf")
groupSize_Control <- as.numeric(table(Control@idents))
# set side-by-side plot displays 
par(mfrow = c(1,2), xpd=TRUE)
p1 <- netVisual_circle(
    Control@net$count, 
    #vertex.weight = groupSize_Control, 
    weight.scale = TRUE, 
    label.edge= T, 
    title.name = "Number of interactions - Control",
    sources.use = c("B_cells","Epithelial_cells","Tissue_stem_cells"),
    targets.use = c("B_cells","Epithelial_cells","Tissue_stem_cells"),
    remove.isolate=T
)

groupSize_Pain <- as.numeric(table(Pain@idents))
p2 <- netVisual_circle(
    Pain@net$count, 
    #vertex.weight = groupSize_Pain, 
    weight.scale = TRUE, 
    label.edge= T, 
    title.name = "Number of interactions - Pain",
    sources.use = c("B_cells","Epithelial_cells","Tissue_stem_cells"),
    targets.use = c("B_cells","Epithelial_cells","Tissue_stem_cells"),
    remove.isolate=T
)
dev.off()

pdf("Interaction.weights.strength.pdf")
groupSize_Control <- as.numeric(table(Control@idents))
# set side-by-side plot displays 
par(mfrow = c(1,2), xpd=TRUE)
p1 <- netVisual_circle(
    Control@net$weight, 
    #vertex.weight = groupSize_Control, 
    weight.scale = TRUE, 
    label.edge= F, 
    title.name = "Interaction weights/strength - Control",
    sources.use = c("B_cells","Epithelial_cells","Tissue_stem_cells"),
    targets.use = c("B_cells","Epithelial_cells","Tissue_stem_cells"),
    remove.isolate=T
)

groupSize_Pain <- as.numeric(table(Pain@idents))
p2 <- netVisual_circle(
    Pain@net$weight, 
    #vertex.weight = groupSize_Pain, 
    weight.scale = TRUE, 
    label.edge= FALSE, 
    title.name = "Interaction weights/strength - Pain",
    sources.use = c("B_cells","Epithelial_cells","Tissue_stem_cells"),
    targets.use = c("B_cells","Epithelial_cells","Tissue_stem_cells"),
    remove.isolate=T
)
dev.off()

pdf("Interaction.weights.strength.Numbered.pdf")
groupSize_Control <- as.numeric(table(Control@idents))
# set side-by-side plot displays 
par(mfrow = c(1,2), xpd=TRUE)
p1 <- netVisual_circle(
    Control@net$weight, 
    #vertex.weight = groupSize_Control, 
    weight.scale = TRUE, 
    label.edge= TRUE, 
    title.name = "Interaction weights/strength - Control",
    sources.use = c("B_cells","Epithelial_cells","Tissue_stem_cells"),
    targets.use = c("B_cells","Epithelial_cells","Tissue_stem_cells"),
    remove.isolate=T
)
groupSize_Pain <- as.numeric(table(Pain@idents))
p2 <- netVisual_circle(
    Pain@net$weight, 
    #vertex.weight = groupSize_Pain, 
    weight.scale = TRUE, 
    label.edge= TRUE, 
    title.name = "Interaction weights/strength - Pain",
    sources.use = c("B_cells","Epithelial_cells","Tissue_stem_cells"),
    targets.use = c("B_cells","Epithelial_cells","Tissue_stem_cells"),
    remove.isolate=T
)
dev.off()



pdf("cell_group_pathways_Pain.pdf")
netVisual_chord_gene(Pain,
    sources.use = c("B_cells","Epithelial_cells","Tissue_stem_cells"),
    targets.use = c("B_cells","Epithelial_cells","Tissue_stem_cells"),
    slot.name = "netP", legend.pos.x = 10,
    thresh = 0.01)
dev.off()


Niche4_net <- subsetCommunication(Pain,         
    sources.use = c("B_cells","Epithelial_cells","Tissue_stem_cells"),
    targets.use = c("B_cells","Epithelial_cells","Tissue_stem_cells"))

str(Niche4_net)
sorted_Niche4_net <- Niche4_net %>% arrange(prob)
head(sorted_Niche4_net)

## filter for only signifcant interactions
sorted_Niche4_net <- filter(sorted_Niche4_net, pval<0.01)
str(sorted_Niche4_net)

write.csv(sorted_Niche4_net, "Pain_Niche4_CellChat.csv")
unique(sorted_Niche4_net$pathway_name)

## count number of interactions per category... maybe a way to split plots?
counts_annotation <- setNames(data.frame(table(sorted_Niche4_net$annotation)), c("annotation", "LR_pairs"))

pie_labels <- paste0(counts_annotation$LR_pairs," = ", round(100 * counts_annotation$LR_pairs/sum(counts_annotation$LR_pairs), 1), "%")
pdf("Niche4_pie.pdf")
pie(counts_annotation$LR_pairs, labels = paste0(counts_annotation$annotation," (", pie_labels,")"))
dev.off()




#cell cell contact
## filter for only signifcant interactions
Niche4_net_ccc <- filter(sorted_Niche4_net, annotation == "Cell-Cell Contact")
str(Niche4_net_ccc)
write.csv(Niche4_net_ccc, "Pain_Niche4_CCC.csv")
counts_pathway_ccc <- setNames(data.frame(table(Niche4_net_ccc$pathway_name)), c("annotation", "LR_pairs"))
   annotation LR_pairs
1      ADGRE5        2
2         APP        6
3        CADM        4
4        CD22        1
5        CD45        1
6        CD46        2
7         CDH        1
8        CDH1        1
9         JAM        1
10      NOTCH        2
11       OCLN        1
12      PTPRM        4
13       SELL        1
14      SEMA4        1
免疫调控：ADGRE5、CD22、CD45、SELL、SEMA4

pdf(file = "CellChat_Pain_Niche4_CCC_pathways_pie.pdf")
#pie_labels <- paste0(" (",counts_pathway_ccc$LR_pairs," = ", round(100 * counts_pathway_ccc$LR_pairs/sum(counts_pathway_ccc$LR_pairs), 2), "%",")")
pie(counts_pathway_ccc$LR_pairs, labels = paste0(counts_pathway_ccc$annotation," (",counts_pathway_ccc$LR_pairs,")"))
dev.off()

Niche4_net_ccc_pairs <- data.frame(interaction_name = Niche4_net_ccc$interaction_name,interaction_name_2 = Niche4_net_ccc$interaction_name_2)
Niche4_net_ccc_pairs

options(stringsAsFactors = FALSE)
bubble <- netVisual_bubble(
    Pain, 
    #sources.use = c(1:4), 
    #targets.use = c(16,17), 
    # signaling = c("MHC-I","MHC-II"),
    pairLR.use = Niche4_net_ccc_pairs,
    remove.isolate = FALSE,
    thresh=0.05,
    sources.use = c("B_cells","Epithelial_cells","Tissue_stem_cells"),
    targets.use = c("B_cells","Epithelial_cells","Tissue_stem_cells"),
    title.name = "Niche4 Niche - Cell-Cell Contact")
ggsave(filename="CellChat_Pain_Niche4_CCC_pairs_heatmap.pdf", plot=bubble, dpi = 300)
str(Niche4_net_ccc)


## filter for only signifcant interactions
Niche4_net_ccc_test <- filter(Niche4_net_ccc, ligand == "PTPRM")
str(Niche4_net_ccc_test)
pairLR.EPHA <- extractEnrichedLR(Pain, signaling = "PTPRM", geneLR.return = FALSE)
pairLR.EPHA

pdf("CellChat_Pain_CCC_Niche4.PTPRM.pdf")
netAnalysis_contribution(Pain, signaling = "PTPRM")
dev.off()
pdf("CellChat_Pain_Niche4_CCC_PTPRM_PTPRM_PTPRM_chord.pdf")
netVisual_individual(Pain, signaling = "PTPRM", pairLR.use = "PTPRM_PTPRM", layout = "chord")
dev.off()


Niche4_net_ccc_test <- filter(Niche4_net_ccc, ligand == "APP")
str(Niche4_net_ccc_test)
pairLR.EPHA <- extractEnrichedLR(Pain, signaling = "APP", geneLR.return = FALSE)
pairLR.EPHA
pdf("CellChat_Pain_CCC_Niche4.APP.pdf")
netAnalysis_contribution(Pain, signaling = "APP")
dev.off()
pdf("CellChat_Pain_Niche4_CCC_APP_APP_CD74_chord.pdf")
netVisual_individual(Pain, signaling = "APP", pairLR.use = "APP_CD74", layout = "chord")
dev.off()




## filter for only signifcant interactions
Niche4_net_ssg <- filter(sorted_Niche4_net, annotation == "Secreted Signaling")
str(Niche4_net_ssg)
write.csv(Niche4_net_ssg, "Pain_Niche4_SSG.csv")

counts_pathway_ssg <- setNames(data.frame(table(Niche4_net_ssg$pathway_name)), c("annotation", "LR_pairs"))
  annotation LR_pairs
1        EGF        2
2        IGF        1
3        MIF       12
4         MK        2
5   VISFATIN        2
MIF	免疫调节、炎症	自身免疫病、肿瘤免疫逃逸
VISFATIN	NAD+代谢、炎症	代谢综合征、衰老相关疾病

pdf(file = "CellChat_Pain_Niche4_SSG_pathways_pie.pdf")
#pie_labels <- paste0(" (",counts_pathway_ssg$LR_pairs," = ", round(100 * counts_pathway_ssg$LR_pairs/sum(counts_pathway_ssg$LR_pairs), 2), "%",")")
pie_labels <- paste0(" (",counts_pathway_ssg$LR_pairs,")")
pie(counts_pathway_ssg$LR_pairs, labels = paste0(counts_pathway_ssg$annotation,pie_labels))
dev.off()

Niche4_net_ssg_pairs <- data.frame(interaction_name = Niche4_net_ssg$interaction_name,interaction_name_2 = Niche4_net_ssg$interaction_name_2)
Niche4_net_ssg_pairs

bubble <- netVisual_bubble(
    Pain, 
    #sources.use = c(1:4), 
    #targets.use = c(16,17), 
    # signaling = c("MHC-I","MHC-II"),
    pairLR.use = Niche4_net_ssg_pairs,
    remove.isolate = FALSE,
    thresh=0.05,
    sources.use = c("B_cells","Epithelial_cells","Tissue_stem_cells"),
    targets.use = c("B_cells","Epithelial_cells","Tissue_stem_cells"),
    title.name = "Niche4 Niche - Secreted Signaling",)
ggsave(filename="CellChat_Pain_Niche4_SSG_pairs_heatmap.pdf", plot=bubble, dpi = 300)

## filter for only signifcant interactions
Niche4_net_ssg_test <- filter(Niche4_net_ssg, ligand == "MIF")
str(Niche4_net_ssg_test)
pairLR.EPHA <- extractEnrichedLR(Pain, signaling = "MIF", geneLR.return = FALSE)
pairLR.EPHA

pdf("CellChat_Pain_SSG_Niche4.MIF.pdf")
netAnalysis_contribution(Pain, signaling = "MIF")
dev.off()

pdf("CellChat_Pain_Niche4_SSG_MIF_CD74+CXCR4_chord.pdf")
netVisual_individual(Pain, signaling = "MIF", pairLR.use = "MIF_CD74_CXCR4", layout = "chord")
dev.off()

pdf("CellChat_Pain_Niche4_SSG_MIF_CD74+CD44_chord.pdf")
netVisual_individual(Pain, signaling = "MIF", pairLR.use = "MIF_CD74_CD44", layout = "chord")
dev.off()

pdf("CellChat_Pain_Niche4_SSG_MIF_MIF_ACKR3_chord.pdf")
netVisual_individual(Pain, signaling = "MIF", pairLR.use = "MIF_ACKR3", layout = "chord")
dev.off()
#(3) MIF (Macrophage Migration Inhibitory Factor)
#关联疼痛：
#MIF通过结合CD74/CXCR4受体，激活脊髓小胶质细胞，释放IL-1β、TNF-α等促炎因子，导致中枢敏化（慢性疼痛的关键机制）。
#在关节炎和神经损伤模型中，MIF抑制剂可减轻痛觉超敏。



#ECM receptor
## filter for only signifcant interactions
Niche4_net_ecm <- filter(sorted_Niche4_net, annotation == "ECM-Receptor")
str(Niche4_net_ecm)
write.csv(Niche4_net_ecm, "Pain_Niche4_ECM.csv")

counts_pathway_ecm <- setNames(data.frame(table(Niche4_net_ecm$pathway_name)), c("annotation", "LR_pairs"))
  annotation LR_pairs
1   COLLAGEN       17
2    LAMININ       15


pdf(file = "CellChat_Pain_Niche4_ECM_pathways_pie.pdf")
#pie_labels <- paste0(" (",counts_pathway_ssg$LR_pairs," = ", round(100 * counts_pathway_ssg$LR_pairs/sum(counts_pathway_ssg$LR_pairs), 2), "%",")")
pie_labels <- paste0(" (",counts_pathway_ecm$LR_pairs,")")
pie(counts_pathway_ecm$LR_pairs, labels = paste0(counts_pathway_ecm$annotation,pie_labels))
dev.off()

Niche4_net_ecm_pairs <- data.frame(interaction_name = Niche4_net_ecm$interaction_name,interaction_name_2 = Niche4_net_ecm$interaction_name_2)
Niche4_net_ecm_pairs

bubble <- netVisual_bubble(
    Pain, 
    #sources.use = c(1:4), 
    #targets.use = c(16,17), 
    # signaling = c("MHC-I","MHC-II"),
    pairLR.use = Niche4_net_ecm_pairs,
    remove.isolate = FALSE,
    thresh=0.05,
    sources.use = c("B_cells","Epithelial_cells","Tissue_stem_cells"),
    targets.use = c("B_cells","Epithelial_cells","Tissue_stem_cells"),
    title.name = "Niche4 Niche - ECM-Receptor")
ggsave(filename="CellChat_Pain_Niche4_ECM_pairs_heatmap.pdf", plot=bubble, dpi = 300)