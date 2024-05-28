library(dplyr)       # tidy data 
library(Seurat)      # single cell infrastructure
library(ggplot2)     # plots
library(dplyr)         # tidy data 
library(Seurat)        # single cell infrastructure
library(ggplot2)       # pretty plots
library(janitor)       # clean data
library(magrittr)      # pipes
library(SCISSORS)      # our package
library(mixtools)      # Gaussian mixture model estimation
library(ggsignif)      # significance bars
library(patchwork)     # align plots
library(latex2exp)     # LaTeX
library(paletteer)     # color palettes
library(kableExtra)    # pretty tables
library(wesanderson)   # more color palettes
library(gridExtra)
library(pals)
gc()



plots_Elyada <- list()


# whole data t-sne--------------------------------------------------------------

sc <- readRDS("/work/users/c/h/changfei/01_CAPER_Paper/00_DeCAF_paper_summary/Fig1/Fig1_A/data/Elyada_final.Rds")

##change label for plot
sc@meta.data %<>% mutate(
  label_paper_plot = case_when(seurat_clusters %in% c(0) ~ "Macrophage", 
                               seurat_clusters %in% c(1) ~ "NK", 
                               seurat_clusters %in% c(2) ~ "Monocyte",  
                               seurat_clusters %in% c(3) ~ "T-1",  
                               seurat_clusters %in% c(4) ~ "B",  
                               seurat_clusters %in% c(5) ~ "Ductal", 
                               seurat_clusters %in% c(6) ~ "DC", 
                               seurat_clusters %in% c(7) ~ "T-2", 
                               seurat_clusters %in% c(8) ~ "Fibroblast", 
                               seurat_clusters %in% c(9) ~ "Acinar",
                               seurat_clusters %in% c(10) ~ "Normal Stroma", 
                               seurat_clusters %in% c(11) ~ "Plasma",
                               seurat_clusters %in% c(12) ~ "pDC", 
                               
  ))

##set color maping
color_wholedata <- c(
  "Macrophage"    = "#1F78C8",
  "NK"            = "#ff0000",
  "Monocyte"      = "#33a02c", 
  "T-1"           = "#6A33C2",  
  "B"             = "#ff7f00", 
  "Ductal"        = "#565656", 
  "DC"            = "#FFD700", 
  "T-2"           = "#a6cee3",
  "Fibroblast"    = "#FB6496", 
  "Acinar"        = "#b2df8a",
  "Normal Stroma" = "#CAB2D6",
  "Plasma"        = "#FDBF6F",
  "pDC"           = "#999999"
)

## t-sne plot
plots_Elyada[[1]] <- DimPlot(sc, reduction = "tsne",group.by = "label_paper_plot", label = TRUE,pt.size = 0.01) + 
  scale_color_manual(values = color_wholedata)+
  theme(legend.position = "none",
        axis.text = element_blank(),
        axis.title.x = element_text(size = 20),  # Adjust the font size for the x-axis label
        axis.title.y = element_text(size = 20),   # Adjust the font size for the y-axis label
        axis.line = element_blank(),  # Remove both x-axis and y-axis lines
        axis.ticks = element_blank()  # Optionally, remove ticks as well
  )+
  labs(title="",
       y = element_blank(),  # Remove y-axis label
       x = element_blank(),  # Remove x-axis label
       # y = element_text("t-SNE-2"), 
       # x = element_text("t-SNE-1"), 
       fill = NULL)
plots_Elyada[[1]]




#fib tsne-----------------------------------------------------------------------

sc <- readRDS("/work/users/c/h/changfei/01_CAPER_Paper/00_DeCAF_paper_summary/Fig1/Fig1_A/data/Elyada_fibro.Rds")

## change label for plot
sc@meta.data %<>% mutate(
  label_paper_plot = case_when(
    label %in% c("myCAF") ~ "permCAF", 
    label %in% c("iCAF") ~ "restCAF",  
    label %in% c("apCAF") ~ "apCAF",  
    label %in% c("Perivascular" ) ~ "Pericyte",  
    label %in% c("Endothelial") ~ "Endothelial"))
###### set color maping
color_mapping <- c(
  "permCAF"               = "#FE46A5",
  "restCAF"               = "#008080",
  "apCAF"                 = "#9683EC",
  "Pericyte"              = "#E4D00A",
  "Endothelial"           = "#98988A"
)
###### t-sne plot
plots_Elyada[[7]] <- DimPlot(sc, reduction = "tsne",label = TRUE,group.by = "label_paper_plot",pt.size = 0.01,label.size = 8) + 
  scale_color_manual(values = color_mapping)+
  theme(legend.position = "none",
        axis.text = element_blank(),
        axis.title.x = element_text(size = 20),  # Adjust the font size for the x-axis label
        axis.title.y = element_text(size = 20),   # Adjust the font size for the y-axis label
        axis.line = element_blank(),  # Remove both x-axis and y-axis lines
        axis.ticks = element_blank()  # Optionally, remove ticks as well
  )+
  labs(title="",
       y = element_blank(),  # Remove y-axis label
       x = element_blank(),  # Remove x-axis label
       # y = element_text("t-SNE-2"), 
       # x = element_text("t-SNE-1"), 
       fill = NULL)

plots_Elyada[[7]]


## plot all figures-------------------------------------------------------------
pdf("/work/users/c/h/changfei/01_CAPER_Paper/00_DeCAF_paper_summary/Fig1/Fig1_A/figures/DeCAFpaper_figure1_Elyada_VAM_final.pdf",width = 28,height =10)
# grid.arrange(grobs = plots[1:12], ncol = 6)
grid.arrange(grobs = plots_Elyada[1:12], ncol = 6)
dev.off()

