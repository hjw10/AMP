##Misty
outdir = "/public/workspace/stu21230110/AMP/09Misty/"

defined_cols = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080')

library(argparse)
library(tidyverse)
library(Seurat)
library(mistyR)
source("misty_utilities.R") 
setwd(outdir)

slide = readRDS("/public/workspace/stu21230110/AMP/01data/rds/spatial.NMF.fact6.new.frequence.rds")

run_colocalization <- function(slide, 
                               assay, 
                               useful_features, 
                               out_label, 
                               misty_out_alias = outdir) {  ###输出目录大家自己制定
   
  # Define assay of each view ---------------
  view_assays <- list("main" = assay,
                      "juxta" = assay,
                      "para" = assay)
  # Define features of each view ------------
  view_features <- list("main" = useful_features, 
                        "juxta" = useful_features,
                        "para" = useful_features)
  # Define spatial context of each view -----
  view_types <- list("main" = "intra", 
                     "juxta" = "juxta",
                     "para" = "para")
  # Define additional parameters (l in case of paraview,
  # n of neighbors in case of juxta) --------
  view_params <- list("main" = NULL, 
                      "juxta" = 2,
                      "para" = 5)
  
  misty_out <- paste0(misty_out_alias, 
                      out_label, "_", assay)
  
  run_misty_seurat(visium.slide = slide,
                   view.assays = view_assays,
                   view.features = view_features,
                   view.types = view_types,
                   view.params = view_params,
                   spot.ids = NULL,
                   out.alias = misty_out)
  
  return(misty_out)
}


DefaultAssay(slide) <- 'predictions'

useful_features <- rownames(slide)   ####也可以自我设定感兴趣的细胞类型

#useful_features <- useful_features[! useful_features %in% "prolif"]

slide_control = subset(slide,treatment == "control")

slide_pain = subset(slide,treatment == "pain")

slide_id = "control"
slide_id = "pain"

mout <- run_colocalization(slide = slide_pain,
                     useful_features = useful_features,
                     out_label = slide_id,
                     assay = "predictions",
                     misty_out_alias = outdir)


misty_res_slide <- collect_results(mout)
  
  plot_folder <- paste0(mout, "/plots")
  
  system(paste0("mkdir ", plot_folder))
  
  pdf(file = paste0(plot_folder, "/", slide_id, "_", "summary_plots.pdf"))
  
  mistyR::plot_improvement_stats(misty_res_slide)
  mistyR::plot_view_contributions(misty_res_slide)
  
  mistyR::plot_interaction_heatmap(misty_res_slide, "intra", cutoff = 0)
  mistyR::plot_interaction_communities(misty_res_slide, "intra", cutoff = 0.5)
  
  mistyR::plot_interaction_heatmap(misty_res_slide, "juxta_2", cutoff = 0)
  mistyR::plot_interaction_communities(misty_res_slide, "juxta_2", cutoff = 0.5)
  
  mistyR::plot_interaction_heatmap(misty_res_slide, "para_5", cutoff = 0)
  mistyR::plot_interaction_communities(misty_res_slide, "para_5", cutoff = 0.5)
  
  dev.off()