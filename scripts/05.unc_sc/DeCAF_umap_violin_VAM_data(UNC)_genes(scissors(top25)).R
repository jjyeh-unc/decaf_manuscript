library(dplyr)       # tidy data 
library(Seurat)      # single cell infrastructure
library(ggplot2)       # pretty plots
library(janitor)       # clean data
library(magrittr)      # pipes
library(mixtools)      # Gaussian mixture model estimation
library(ggsignif)      # significance bars
library(patchwork)     # align plots
library(latex2exp)     # LaTeX
library(paletteer)     # color palettes
library(kableExtra)    # pretty tables
library(wesanderson)   # more color palettes
library(pals)
gc()



#UMAP---------------------------------------------------------------------------
plots <- list()

##Yehseq whole data umap

sc <- readRDS("/work/users/c/h/changfei/01_CAPER_Paper/04_scRNA_process/data/UNC_sc_all.rds")
# 
# ### change label for plot
# sc@meta.data %<>% mutate(
#   # more accurate ref. pop. now that I've annotated the immune cells
#   label_paper_plot = case_when(seurat_clusters %in% c(0) ~ "Fibroblast-1", 
#                                seurat_clusters %in% c(1) ~ "T", 
#                                seurat_clusters %in% c(2) ~ "Macrophage",  
#                                seurat_clusters %in% c(3) ~ "Epithelial",  
#                                seurat_clusters %in% c(4) ~ "Fibroblast-2",  
#                                seurat_clusters %in% c(5) ~ "Endothelial", 
#                                seurat_clusters %in% c(6) ~ "Fibroblast-3", 
#                                seurat_clusters %in% c(7) ~ "NK", 
#                                seurat_clusters %in% c(8) ~ "B", 
#                                seurat_clusters %in% c(9) ~ "Schwann"
#   ))
# 
# saveRDS(sc,"/work/users/c/h/changfei/01_CAPER_Paper/00_DeCAF_paper_summary/Fig1/Fig1_FGHI/data/UNC_sc_all.rds")

###### define color
color_names <- unique(sc@meta.data$label_paper_plot)
color_wholedata <- cols25(n=10)
names(color_wholedata) <- color_names
color_wholedata["Fibroblast-1"]="violetred1"
color_wholedata["Fibroblast-2"]="#008080"
color_wholedata["Fibroblast-3"]="#E4D00A"

###### umap plot
plots[[1]] <- DimPlot(sc, reduction = "umap",group.by = "label_paper_plot", label = TRUE,pt.size = 0.01,label.size = 10) + 
  scale_color_manual(values = color_wholedata)+
  theme(legend.position = "none",
        axis.text = element_blank(),
        axis.title.x = element_text(size = 20),  # Adjust the font size for the x-axis label
        axis.title.y = element_text(size = 20),   # Adjust the font size for the y-axis label
        axis.line = element_blank(),  # Remove both x-axis and y-axis lines
        axis.ticks = element_blank()  # Optionally, remove ticks as well
  )+
  labs(
    y = element_blank(),  # Remove y-axis label
    x = element_blank(),  # Remove x-axis label
    # y = element_text("t-SNE-2"), 
    # x = element_text("t-SNE-1"), 
    fill = NULL)+
  ggtitle("")

plots[[1]]


## Yehseq fib ump

sc <- readRDS("/work/users/c/h/changfei/01_CAPER_Paper/04_scRNA_process/data/UNC_sc_CAF.rds")

unique(sc@meta.data$label_seperate)
###define color
color_mapping <- c(
  "proCAF-1"                 = "palevioletred1",
  "proCAF-2"                 = "violetred1",
  "restCAF-1"                 = "turquoise4",
  "restCAF-2"                 = "lightseagreen",
  "restCAF-3"                 = "steelblue1",
  "apCAF"                     = "#9683EC",
  "Pericyte"                  = "#E4D00A",
  "Endothelial"               = "#98988A"
)

###### umap
plots[[2]] <- DimPlot(sc, reduction = "umap",label = TRUE,group.by = "label_seperate",pt.size = 0.5,label.size = 10) + 
  scale_color_manual(values = color_mapping)+
  theme(legend.position = "none",
        axis.text = element_blank(),
        axis.title.x = element_text(size = 20),  # Adjust the font size for the x-axis label
        axis.title.y = element_text(size = 20),   # Adjust the font size for the y-axis label
        axis.line = element_blank(),  # Remove both x-axis and y-axis lines
        axis.ticks = element_blank()  # Optionally, remove ticks as well
  )+
  labs(
    y = element_blank(),  # Remove y-axis label
    x = element_blank(),  # Remove x-axis label
    # y = element_text("t-SNE-2"), 
    # x = element_text("t-SNE-1"), 
    fill = NULL)+
  ggtitle("")

plots[[2]]




# Yehseq fib VAM---------------------------------------------------------------------------

sc <- readRDS("/work/users/c/h/changfei/01_CAPER_Paper/04_scRNA_process/data/UNC_sc_CAF.rds")

# marker gene list
## Elyada_markers
load("/work/users/c/h/changfei/Other_Markers/Our_markerGenes/cmbSubtypes.230508.RData")
Elyada_genes <- readRDS("/work/users/c/h/changfei/Other_Markers/Our_markerGenes/Genes_Elyada.rds")
schemaList
subtypeGeneList[17]
MarkerGene_list_Elyada <- list()
###panCAF
MarkerGene_list_Elyada[["Elyada-panCAF"]] <- Elyada_genes$`Elyada-PanCAF`
###iCAF,myCAF
marker_info <- subtypeGeneList[[14]]
Cell_type <- colnames(marker_info)[2:3]
for (i in 1:length(Cell_type)) {
  Gene_marker_set_name <- paste0("Elyada-CAF-",Cell_type[i])
  col_index <- which(colnames(marker_info)==Cell_type[i])
  MarkerGene_list_Elyada[[Gene_marker_set_name]] <- marker_info[which(marker_info[,col_index]=="TRUE"),]$geneSymbol
}
###apCAF
marker_info <- subtypeGeneList[[15]]
Cell_type <- colnames(marker_info)[2:3]
Gene_marker_set_name <- paste0("Elyada-CAF-mouse-",Cell_type[1])
col_index <- which(colnames(marker_info)==Cell_type[1])
MarkerGene_list_Elyada[[Gene_marker_set_name]] <- marker_info[which(marker_info[,col_index]=="TRUE"),]$geneSymbol
###Perivascular
MarkerGene_list_Elyada[["Elyada-Perivascular"]] <- Elyada_genes$`Elyada-Perivascular`
MarkerGene_list_Elyada <- MarkerGene_list_Elyada[c(1,3,2,4,5)]

###get gene set
gene_sets1 <- readRDS("/work/users/c/h/changfei/00_SCISSORS/data/VAM_Markers_version5th/Top5/combined/VAM_Scissors_Allmarkers_top5_combined.Rds")
gene_sets2 <- readRDS("/work/users/c/h/changfei/00_SCISSORS/data/VAM_Markers_version5th/Top10/combined/VAM_Scissors_Allmarkers_Top10_combined.Rds")
gene_sets3 <- readRDS("/work/users/c/h/changfei/00_SCISSORS/data/VAM_Markers_version5th/Top25/combined/VAM_Scissors_Allmarkers_Top25_combined.Rds")
gene_sets4 <- readRDS("/work/users/c/h/changfei/00_SCISSORS/data/VAM_Markers_version5th/All/combined/VAM_Scissors_Allmarkers_All_combined.Rds")
gene_sets5 <- readRDS("/users/c/h/changfei/old/Yehseq/Gene_sets.rds")

gene_sets <- c(gene_sets1,gene_sets2,gene_sets3,gene_sets4,gene_sets5,MarkerGene_list_Elyada)
gene_sets$`SCISSORS-CAFTop25-apCAF`
gene_sets$`SCISSORSCAFtop25-apCAF`

### find genes
source('/work/users/c/h/changfei/R_functions/VAM_Functions.R')
gene.set.collection = createGeneSetCollection_mouse_human (gene.ids=sc@assays$SCT@data@Dimnames[[1]],
                                                           gene.set.collection=gene_sets)
length(gene.set.collection)
###### VAM
sc@active.assay <- "SCT"
sc <- vamForSeurat(sc, 
                   gene.set.collection = gene.set.collection, 
                   center = FALSE, 
                   gamma = TRUE, 
                   sample.cov = FALSE, 
                   return.dist = TRUE)
DefaultAssay(object = sc) = "VAMcdf"


###prepare for VAM
plots_VAM <- list()
plots_violin <- list()
### SCISSORS
# figure_title_SCISSORS <- c("SCISSORS-panCAF", "SCISSORS-myCAF-like", "SCISSORS-iCAF-like", "SCISSORS-apCAF-like","SCISSORS-PSC")
figure_title_SCISSORS <- c("panCAF", "proCAF", "restCAF", "apCAF","Pericyte")

figure_color_SCISSORS <- c("red","violetred1","turquoise4","#9683EC", "#E4D00A")

gene.set.collection_SCISSORS <- list(
  
  top25 = c( "SCISSORS-panCAFTop25-panCAF",
             "SCISSORS-CAFTop25-myCAF",
             "SCISSORS-CAFTop25-iCAF" ,
             # "SCISSORS-CAFTop25-apCAF",
             "SCISSORSCAFtop25-apCAF",
             "SCISSORS-StromaTop25-Perivascular"
             )

)

### prepare data for violinplots
data_input <- data.frame(t(sc@assays$VAMcdf@data)) %>%
  bind_cols((sc@meta.data %>% dplyr::select(label_seperate)))

color_mapping <- c(
  "proCAF-1"                 = "palevioletred1",
  "proCAF-2"                 = "violetred1",
  "restCAF-1"                 = "turquoise4",
  "restCAF-2"                 = "lightseagreen",
  "restCAF-3"                 = "steelblue1",
  "apCAF"                     = "#9683EC",
  "Pericyte"                  = "#E4D00A",
  "Endothelial"               = "#98988A"
)

### VAM and violin plot
i=1
for(name in names(gene.set.collection_SCISSORS)){
  features_tmp <- gene.set.collection_SCISSORS[[name]]
  for (featuresTmp in 1:length(features_tmp)) {
    plots_VAM[[i]] <- FeaturePlot(sc, reduction = "umap", features = features_tmp[featuresTmp], cols  =c("lightgrey", figure_color_SCISSORS[featuresTmp]))+
      theme(legend.position = "none",
            axis.text = element_blank(),
            axis.title.x = element_text(size = 20),  # Adjust the font size for the x-axis label
            axis.title.y = element_text(size = 20),# Adjust the font size for the y-axis label
            plot.title = element_text(size = 30), # Adjust the font size for the plot title
            axis.line = element_blank(),  # Remove both x-axis and y-axis lines
            axis.ticks = element_blank()  # Optionally, remove ticks as well
      )+
      labs(
        y = element_blank(),  # Remove y-axis label
        x = element_blank(),  # Remove x-axis label
        # y = element_text("t-SNE-2"), 
        # x = element_text("t-SNE-1"), 
        fill = NULL)+
      ggtitle("")

    ## Violin plot
    plots_violin[[i]] <- ggplot(data_input,aes_string(x = "label_seperate", y = gsub("-",".",features_tmp[featuresTmp]), fill = "label_seperate")) +
      geom_violin(scale = "width", draw_quantiles = 0.5, size = 0.5) +
      scale_fill_manual(values = color_mapping)+
      labs(y = "VAM Score", x = NULL, fill = NULL) +
      ggtitle(paste0(figure_title_SCISSORS[featuresTmp]))+
      # ggtitle(figure_title[j])+
      theme(legend.position = "none",
            panel.background = element_blank(),
            axis.line = element_line(color = "black"),
            plot.title = element_text(hjust = 0.5, size =35 ,face = "bold"),
            axis.title.y = element_text(size = 30),
            axis.text.x = element_text(angle = 45, hjust = 1,size = 30),
            axis.text.y = element_text(size = 30)) 
    
    i=i+1
  }
}



# plot all figures-------------------------------------------------------------
pdf("/work/users/c/h/changfei/01_CAPER_Paper/00_DeCAF_paper_summary/Fig1/Fig1_FGHI/figures/DeCAFpaper_Yehseq_QCSeurat_data_umap_VAM_violin_scissors(top25)_markers.pdf",width = 37,height =12)
# grid.arrange(grobs = plots[1:12], ncol = 6)
grid.arrange(grobs = c(plots[1],plots_violin,plots[2],plots_VAM), ncol = 6)

dev.off()
