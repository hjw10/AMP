library(pheatmap)

outdir = "/public/workspace/stu21230110/AMP/08enrichment/"

gsea_table_new_cm = read.csv(paste0(outdir, "GO_enrichment/GO_enrichment_all.csv"),header = T, sep = '\t',row.names = 1)

head(gsea_table_new_cm)

pdf(paste0(outdir,"GO_enrichment/GO_heatmap.pdf"))
my_col <- colorRampPalette(c("white","orange", "red","darkred"))(50)
my_col[0] <- "white"
p<-pheatmap(gsea_table_new_cm,color = my_col)
dev.off()



gsea_table_new_cm = read.csv(paste0(outdir, "MSigDB_enrichment/MSigDB_Hallmark_2020_enrichment_all.csv"),header = T, sep = '\t',row.names = 1)
pdf(paste0(outdir,"MSigDB_enrichment/MSigDB_heatmap.pdf"))
my_col <- colorRampPalette(c("white","orange", "red","darkred"))(50)
my_col[0] <- "white"
p<-pheatmap(gsea_table_new_cm,color = my_col)
dev.off()
