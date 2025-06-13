##系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)

library(rhdf5)
library(dplyr)
library(data.table)
library(Matrix)
library(rjson)
library(Seurat)
setwd("/public/workspace/stu21230110/AMP/01data/rds/")

mydata = read.csv("mat.csv")

meta = read.table("metadata.tsv",sep="\t",header=T,row.names=1)

pos = read.table("position_spatial.tsv",sep="\t",header=T,row.names=1)

umap = read.table("position_X_umap.tsv",sep="\t",header=T,row.names=1)

rownames(mydata) = mydata[,1]
mydata = mydata[,-1]
mat <- Matrix(t(mydata), sparse = TRUE)
exp = (mydata)

obj <- CreateSeuratObject(t(exp),project='Spatial',assay='Spatial',meta.data=meta)

tissue_lowres_image <- matrix(1, max(pos$y), max(pos$x))
tissue_positions_list <- data.frame(row.names = colnames(obj),tissue =1,row =pos$y, col =pos$x,imagerow=pos$y, imagecol =pos$x)
scalefactors_json <- toJSON(list(fiducial_diameter_fullres = 1, tissue_hires_scalef = 1, tissue_lowres_scalef = 1))
seurat_spatialObj <- obj
generate_spatialObj <- function(image, scale.factors,tissue.positions, filter.matrix = TRUE){
    if(filter.matrix){ 
        tissue.positions <-tissue.positions[which(tissue.positions$tissue == 1), , drop = FALSE]
        }
    unnormalized.radius <- scale.factors$fiducial_diameter_fullres 
    scale.factors$tissue_lowres_scalef
    spot.radius <- unnormalized.radius / max(dim(image))
    return(new(Class = 'VisiumV1',image = image,scale.factors = scalefactors(spot = scale.factors$tissue_hires_scalef,
                                                                             fiducial = scale.factors$fiducial_diameter_fullres,
                                                                             hires = scale.factors$tissue_hires_scalef,
                                                                             lowres = scale.factors$tissue_lowres_scalef),
                                                                             coordinates = tissue.positions,
                                                                             spot.radius = spot.radius))
 }
spatialObj <- generate_spatialObj(image=tissue_lowres_image,scale.factors=fromJSON(scalefactors_json),tissue.positions = tissue_positions_list)
spatialObj <- spatialObj[Cells(seurat_spatialObj)]
DefaultAssay(spatialObj) <- 'Spatial'
seurat_spatialObj[['slice1']] <- spatialObj

predict_assay <- CreateAssayObject(t(as.matrix(seurat_spatialObj@meta.data[,20:30])))

# 将此 assay 添加到之前创建的 Seurat 对象中
seurat_spatialObj[['predictions']] <- predict_assay

saveRDS(seurat_spatialObj,"spatial.NMF.fact6.new.frequence.rds")