#rm(list = ls())
#V5_path = "/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/seurat5/"
#.libPaths(V5_path)
#.libPaths()
library(stringr)
library(Seurat)
library(Augur)
library(dplyr)
library(patchwork)
library(viridis)
library(BiocParallel)
library(ggplot2)

load("/public/workspace/stu21230110/AMP/01data/rds/06_scRNA_harmony_singler202400810.Rdata")
outdir = '/public/workspace/stu21230110/AMP/15augur/'
setwd(outdir)

# check
scRNA = scRNA_harmony_singler
table(scRNA$orig.ident)
table(scRNA$Annotation)


# 此函数的主要目的是评估并计算每种细胞类型在特定条件（如实验组别或时间点）下的表现，通常通过计算区域下曲线（AUC）来量化。它是单一条件下的细胞类型表现的度量，可以帮助识别在特定实验设置下哪些细胞类型表现最为显著或重要。
sc_augur <- calculate_auc(scRNA,cell_type_col="Annotation",label_col="group",n_threads=4)

saveRDS(sc_augur,"scRNA.augur.rds")
#sc_augur <- readRDS("augur.rds")
head(sc_augur$AUC,5)

# A tibble: 11 × 2
   cell_type                  auc
   <chr>                    <dbl>
 1 Epithelial_cells         0.766
 2 B_cells                  0.759
 3 Tissue_stem_cells        0.665
 4 T_cells                  0.652
 5 Smooth_muscle_cells      0.634
 6 Monocyte                 0.615
 7 NK_cells                 0.604
 8 Stromal_cells            0.587
 9 Endothelial_cells        0.582
10 Myeloid_progenitor_cells 0.581
11 Macrophage               0.564

# 明确哪个细胞群体的RNA水平变化最大
# plot_lollipop函数是基于ggplot框架的,所以很容易就可以美化图片
# 还可以提取原始函数进行修改
setwd(outdir)
pdf("sc_augur.pdf")
plot_lollipop(sc_augur) +
    geom_segment(aes(xend = cell_type,yend = 0.5), size = 1) +
    geom_point(size = 3,aes(color = cell_type)) + # 增加点的大小和颜色映射
    scale_color_manual(values = c(
        "B_cells"="#F39B7FFF",
        "Endothelial_cells"="#E377C2FF",
        "Epithelial_cells"="#3C5488FF",
        "Macrophage"="#7E6148FF",
        "Monocyte"="#00A087FF",
        "Myeloid_progenitor_cells"="#8494FF",
        "NK_cells"="#4DBBD5FF",
        "Smooth_muscle_cells"="#9467BDFF",
        "Stromal_cells"="#E64B35FF",
        "T_cells"="#91D1C2FF",
        "Tissue_stem_cells"="#B09C85FF"
))
dev.off()

# 结合umap
# rank模式
pdf('rank_umap.pdf')
plot_umap(sc_augur,
            scRNA,
            mode="rank",
            reduction="umap",
            palette="cividis",
#            augur_mode="default",
            cell_type_col="Annotation")
dev.off()

# default模式
pdf('default_umap.pdf')
plot_umap(sc_augur,
            scRNA,
            mode="default",
            reduction="umap",
            palette="viridis", #"viridis", "plasma", "magma", "inferno" "YlGnBu"
#            augur_mode="default",
            cell_type_col="Annotation")
dev.off()



###颜色
sample_color <- c("#BCBD22FF","#4DBBD5FF","#00A087FF","#91D1C2FF","#F39B7FFF","#F7B6D2FF","#7E6148FF","#E64B35FF","#3C5488FF","#E377C2FF","#8491B4FF","#B09C85FF")
sample_color <- c("#F39B7FFF","#E377C2FF","#3C5488FF","#7E6148FF","#00A087FF","#8494FF","#4DBBD5FF", "#91D1C2FF","#9467BDFF","#E64B35FF","#B09C85FF")
