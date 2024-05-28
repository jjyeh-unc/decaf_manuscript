library("Matrix", lib.loc = "/nas/longleaf/home/changfei/R/x86_64-pc-linux-gnu-library/4.2")
library("SeuratObject", lib.loc = "/nas/longleaf/home/changfei/R/x86_64-pc-linux-gnu-library/4.2")
library("Seurat", lib.loc = "/nas/longleaf/home/changfei/R/x86_64-pc-linux-gnu-library/4.2")
library(VAM)           # single cell GSEA)
library(dplyr)         # tidy data 
library(ggplot2)       # pretty plots
library(janitor)       # clean data
library(magrittr)      # pipes
library(mixtools)      # Gaussian mixture model estimation
library(ggsignif)      # significance bars
library(patchwork)     # align plots
library(latex2exp)     # LaTeX
library(paletteer)     # color palettes
library(kableExtra)    # pretty tables
library(wesanderson)   # more color palettesa
library(SingleCellExperiment)
library(pals)
library(gridExtra)


gc() ## clear the memory

# load datasets------------------------------------------------------
Marker_genes <- readRDS("/work/users/c/h/changfei/Other_Markers/Our_markerGenes/subtypeSchema.231026.rds")

data_dir_list <- list(
  Elyada = "/work/users/c/h/changfei/00_SCISSORS/data/Relabeled_datasets_version2th/Elyada_final.Rds",
  Elyada_fib = "/work/users/c/h/changfei/00_SCISSORS/data/Relabeled_datasets_version2th/Elyada_fibro.Rds",
  Yehseq = "/work/users/c/h/changfei/01_CAPER_Paper/04_scRNA_process/data/UNC_sc_all.rds",
  Yehsesq_fib = "/work/users/c/h/changfei/01_CAPER_Paper/04_scRNA_process/data/UNC_sc_CAF.rds"
)
data_list <- list()
for (name in names(data_dir_list)) {
  data_list[[name]] <- readRDS(data_dir_list[[name]])
}




# marker genes----------------------------------------------------------

## DeCAF genes and Elyada CAF genes

gene_list_1 <- readRDS("/work/users/c/h/changfei/Other_Markers/VAM_gene_list/VAM_Gene_list(DeCAF_Elyada_CAF_fromsubtypeSchema_231026).Rds")
names(gene_list_1) <- names(gene_list_1) %>%
  sub("-CAF","",.)


## SCISSORS CAF genes

all_markers <- readRDS("/work/users/c/h/changfei/00_SCISSORS/data/VAM_Markers_version5th/All/combined/VAM_Scissors_Allmarkers_All_combined.Rds")
CAF_markers <- list()
top <- c(5,10,25,50)
caf_name <- c("SCISSORS-CAF-myCAF","SCISSORS-CAF-iCAF","SCISSORS-CAF-apCAF")
for (i in top) {
  for (name in caf_name) {
    name_tmp <- paste0(name,"-top",i)
    CAF_markers[[name_tmp]] <- head(all_markers[[name]],i)
  }
}
names(CAF_markers)
for (name in caf_name) {
  name_tmp <- paste0(name,"-all")
  CAF_markers[[name_tmp]] <- all_markers[[name]]
}

gene_sets <- CAF_markers
gene_list_2 <- c(CAF_markers["SCISSORS-CAF-myCAF-top25"],CAF_markers["SCISSORS-CAF-iCAF-top25"])
names(gene_list_2) <- c("SCISSORS-permCAF","SCISSORS-restCAF")


## combine gene list

gene_sets <- c(gene_list_1,gene_list_2)
gene_sets <- gene_sets[c("DeCAF-permCAF","DeCAF-restCAF","SCISSORS-permCAF","SCISSORS-restCAF","Elyada-myCAF","Elyada-iCAF")]





# VAM ---------------------------------------------------------------------------
VAM_plot <- list()
source('/work/users/c/h/changfei/R_functions/VAM_Functions.R')
color_mapping <- list(
  "DeCAF-permCAF"    = "violetred1",
  "DeCAF-restCAF"    = "turquoise4",
  "SCISSORS-permCAF" = "violetred1",
  "SCISSORS-restCAF" = "turquoise4",
  "Elyada-myCAF"     = "darkgreen",
  "Elyada-iCAF"      = "coral"
)

for (name in names(data_list)) {
  sc <- data_list[[name]]
  gene.set.collection = createGeneSetCollection_mouse_human (gene.ids=sc@assays$SCT@data@Dimnames[[1]],
                                                             gene.set.collection=gene_sets)
  sc@active.assay <- "SCT"
  sc <- vamForSeurat(sc, 
                     gene.set.collection = gene.set.collection, 
                     center = FALSE, 
                     gamma = TRUE, 
                     sample.cov = FALSE, 
                     return.dist = TRUE)
  DefaultAssay(object = sc) = "VAMcdf"
  
  if(name %in% c("Elyada","Elyada_fib")){
    reduction_VAM  <- "tsne"
  }else{
    reduction_VAM  <- "umap"
  }
  
  for (featuresTmp in names(gene.set.collection)) {
    VAM_plot[[name]][[featuresTmp]] <- FeaturePlot(sc, reduction = reduction_VAM, features = featuresTmp, cols = paletteer_c("grDevices::RdYlBu", 30,direction = -1))+
      theme(
            # legend.position = "none",
            axis.text = element_blank(),
            axis.title.x = element_text(size = 20),  # Adjust the font size for the x-axis label
            axis.title.y = element_text(size = 20),# Adjust the font size for the y-axis label
            plot.title = element_text(size = 25), # Adjust the font size for the plot title
            axis.line = element_blank(),  # Remove both x-axis and y-axis lines
            axis.ticks = element_blank()  # Optionally, remove ticks as well
      )+
      labs(
        y = element_blank(),  # Remove y-axis label
        x = element_blank(),  # Remove x-axis label
        fill = NULL)+
      ggtitle(paste0(featuresTmp))
  }
}


# organize the figures for plotting--------------------------------------------
VAM_plot_Elyada <- c(VAM_plot$Elyada[c(1,3,5)],
                        VAM_plot$Elyada_fib[c(1,3,5)],
                        VAM_plot$Elyada[c(2,4,6)],
                        VAM_plot$Elyada_fib[c(2,4,6)]

)

VAM_plot_Yehseq <- c(
                  VAM_plot$Yehseq[c(1,3,5)],
                  VAM_plot$Yehsesq_fib[c(1,3,5)],
                  VAM_plot$Yehseq[c(2,4,6)],
                  VAM_plot$Yehsesq_fib[c(2,4,6)]
)

## add one figure for key
color_key_figures <- list()
color_key_figures[["key"]] <- FeaturePlot(sc, reduction = reduction_VAM, features = featuresTmp, cols = paletteer_c("grDevices::RdYlBu", 30,direction = -1))


pdf("/work/users/c/h/changfei/01_CAPER_Paper/plot/DeCAFpaper_VAM_DeCAFgenes_SICCORSgenes_Elyadagenes_4th.pdf",width = 30,height =10)
# grid.arrange(grobs = plots[1:12], ncol = 6)
grid.arrange(grobs = VAM_plot_Elyada, ncol = 6)
grid.arrange(grobs = VAM_plot_Yehseq, ncol = 6)
grid.arrange(grobs = color_key_figures[1:12], ncol = 6)
dev.off()









