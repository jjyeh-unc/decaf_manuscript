library(dplyr)       # tidy data 
library(Seurat)      # single cell infrastructure
library(ggplot2)     # plots
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
library(ComplexHeatmap)
library(circlize )
gc()

# marker genes------------------------------------------------------------------
## read marker gene_list
MarkerGene_list_VAM <- readRDS("/work/users/c/h/changfei/01_CAPER_Paper/00_DeCAF_paper_summary/Fig1/Fig1_EJ/data/MarkerGene_list_combine.Rds")

## exclude Helms
names(MarkerGene_list_VAM)
MarkerGene_list_VAM <- MarkerGene_list_VAM[-4]

## transfer all gene markers into one list
gene_list <- list()
for (top_level_name in names(MarkerGene_list_VAM)) {
  sub_list <- MarkerGene_list_VAM[[top_level_name]]
  gene_list <- c(gene_list, sub_list)
}



# Yehseq_fibro heatmap --------------------------------------------------------------

sc <- readRDS("/work/users/c/h/changfei/01_CAPER_Paper/00_DeCAF_paper_summary/Fig1/Fig1_FGHI/data/Yehseq_fibro_7clusters_labeled.Rds")
unique(sc@meta.data$label)

## prepare median VAM matrix-----------------

### VAM
source('/work/users/c/h/changfei/R_functions/VAM_Functions.R')
gene.set.collection = createGeneSetCollection_mouse_human (gene.ids=sc@assays$SCT@data@Dimnames[[1]],
                                                           gene.set.collection=gene_list)
sc@active.assay <- "SCT"
sc <- vamForSeurat(sc, 
                   gene.set.collection = gene.set.collection, 
                   center = FALSE, 
                   gamma = TRUE, 
                   sample.cov = FALSE, 
                   return.dist = TRUE)
DefaultAssay(object = sc) = "VAMcdf"

### prepare data to median VAM vs celltype
VAM_score_info <- sc@assays$VAMcdf@data
VAM_score_info <- as.data.frame(t(sc@assays$VAMcdf@data))
VAM_score_info$label <-sc@meta.data$label_seperate
colnames(VAM_score_info)

VAM_Median <- VAM_score_info %>%
  group_by(label) %>%
  summarize_all(median)%>%
  as.data.frame()

rownames(VAM_Median) <- VAM_Median$label
VAM_Median <- VAM_Median[,-1]
VAM_Median_matrix <- as.matrix(t(VAM_Median))

VAM_Median_matrix_zscore <- t(scale(t(VAM_Median_matrix)))
VAM_Median_matrix_zscore[is.na(VAM_Median_matrix_zscore)] <- 0

VAM_Median_matrix_zscore_df <- as.data.frame(VAM_Median_matrix_zscore)

schema_groups <- names(MarkerGene_list_VAM)
VAM_Median_matrix_zscore_df$schema_groups <- case_when(rownames(VAM_Median_matrix_zscore_df)%in%names(MarkerGene_list_VAM[[schema_groups[1]]]) ~ schema_groups[1], 
                                                       rownames(VAM_Median_matrix_zscore_df)%in%names(MarkerGene_list_VAM[[schema_groups[2]]]) ~ schema_groups[2], 
                                                       rownames(VAM_Median_matrix_zscore_df)%in%names(MarkerGene_list_VAM[[schema_groups[3]]]) ~ schema_groups[3], 
                                                       rownames(VAM_Median_matrix_zscore_df)%in%names(MarkerGene_list_VAM[[schema_groups[4]]]) ~ schema_groups[4], 
                                                       rownames(VAM_Median_matrix_zscore_df)%in%names(MarkerGene_list_VAM[[schema_groups[5]]]) ~ schema_groups[5], 
                                                       rownames(VAM_Median_matrix_zscore_df)%in%names(MarkerGene_list_VAM[[schema_groups[6]]]) ~ schema_groups[6], 
                                                       rownames(VAM_Median_matrix_zscore_df)%in%names(MarkerGene_list_VAM[[schema_groups[7]]]) ~ schema_groups[7], 
                                                       rownames(VAM_Median_matrix_zscore_df)%in%names(MarkerGene_list_VAM[[schema_groups[8]]]) ~ schema_groups[8], 
                                                       rownames(VAM_Median_matrix_zscore_df)%in%names(MarkerGene_list_VAM[[schema_groups[9]]]) ~ schema_groups[9], 
                                                       rownames(VAM_Median_matrix_zscore_df)%in%names(MarkerGene_list_VAM[[schema_groups[10]]]) ~ schema_groups[10]
                                                       
)

VAM_Median_matrix_zscore_df$schema_groups <- factor(VAM_Median_matrix_zscore_df$schema_groups,levels = c("SCISSORS","Elyada","Ogawa","Chen","Grunwald","Hwang","Oh","Wang","MS","Maurer" ))
VAM_Median_matrix_zscore_df <- VAM_Median_matrix_zscore_df %>%
  arrange(schema_groups)

### change row orders as Laura suggested
new_row_order <- c("SCISSORS-CAF-top10-myCAF", "SCISSORS-CAF-top10-iCAF", "SCISSORS-CAF-top10-apCAF", "SCISSORS-Stroma-Top10-Perivascular",
                   "Elyada-CAF-myCAF", "Elyada-CAF-iCAF", "Elyada-CAF-mouse-apCAF", "Elyada-Perivascular",
                   "Ogawa-F-stroma","Ogawa-A-stroma","Ogawa-C-stroma",
                   "Chen-cCAF-Top10","Chen-csCAF-Top10","Chen-PSC-Top10",
                   "Grunwald-Reactive-Inter-Top10","Grunwald-Deserted-Top10",
                   "Hwang-Fibroblast-Myofibroblastic-Top10","Hwang-Fibroblast-Neurotropic-Top10","Hwang-Fibroblast-Adhesive-Top10","Hwang-Fibroblast-Immunomodulatory-Top10",
                   "Oh-Stroma-myCAF-Top10","Oh-Stroma-IL11-CAF-Top10","Oh-Stroma-Myocyte-Top10","Oh-Stroma-iCAF-Top10","Oh-Stroma-csCAF-Top10","Oh-Stroma-Schwann-Top10","Oh-Stroma-qPSC-Top10","Oh-Stroma-smPSC-Top10",
                   "Wang-myCAF-Top10","Wang-nonTypicalCAF-Top10","Wang-iCAF-Top10","Wang-meCAF-Top10","Wang-Stellate-1-Top10","Wang-Stellate-2-Top10",
                   "MS-Activated","MS-Normal",
                   "Maurer-ECM","Maurer-Immune")
index_order <- match(rownames(VAM_Median_matrix_zscore_df),new_row_order)
VAM_Median_matrix_zscore_df <- VAM_Median_matrix_zscore_df[order(index_order),]

### revers the row and col, get prepared matrix
VAM_Median_matrix_zscore_HM <- as.matrix(VAM_Median_matrix_zscore_df[,c(1:7)])%>% 
  t(.)


## heatmap-----------------
### set group color
color_group <- c()
color_group = c(cols25(n=length(levels(VAM_Median_matrix_zscore_df$schema_groups))))
color_group[1] <- "#ff0000"
color_group[2] <- "#1F78C8"
names(color_group)=levels(VAM_Median_matrix_zscore_df$schema_groups)

### change colnames as laura suggested
new_colnames <- colnames(VAM_Median_matrix_zscore_HM) %>% 
  sub("^[^-]+-", "", .)%>% 
  sub("-Top10", "", .)%>% 
  sub("-top10", "", .)%>%
  sub("CAF-", "", .)%>%
  sub("Fibroblast-", "", .)%>%
  sub("Stroma-", "", .)
new_colnames[7] <- "apCAF-mouse"  
new_colnames[1:4] <- c("myCAF-like","iCAF-like","apCAF-like","PSC")

### set bottom annotation
bottom_annotation <- HeatmapAnnotation(
  marker_group = VAM_Median_matrix_zscore_df$schema_groups,
  colname = anno_text(new_colnames, rot = 45, just = "right", offset = unit(1, "npc") - unit(0, "mm")),
    col = list(marker_group = color_group),
  simple_anno_size = unit(1, "mm"),
  show_annotation_name = FALSE,
  show_legend = c(FALSE)
)

p1<-Heatmap(VAM_Median_matrix_zscore_HM,
            name="YehSeq Fib",
            column_split = VAM_Median_matrix_zscore_df$schema_groups,# Specify the group labels
            show_column_names = FALSE, 
            row_names_side = "left",
            bottom_annotation = bottom_annotation,
            cluster_rows = FALSE,
            cluster_columns = FALSE)

p1




#  Elyada_fibro heatmap--------------------------------------------------------------

## read data-----------------
sc <- readRDS("/work/users/c/h/changfei/CAPER_Paper/00_DeCAF_paper_summary/Fig1/Fig1_A/data/Elyada_fibro.Rds")
sc@meta.data$label <- factor(sc@meta.data$label, levels = c("myCAF","iCAF","apCAF","Perivascular","Endothelial"))
sc@meta.data$label_new <- case_when(sc@meta.data$label == "myCAF" ~ "myCAF-like",
                                    sc@meta.data$label == "iCAF" ~ "iCAF-like",
                                    sc@meta.data$label == "apCAF" ~ "apCAF-like",
                                    sc@meta.data$label == "Perivascular" ~ "PSC",
                                    sc@meta.data$label == "Endothelial" ~ "Endothelial")
  
sc@meta.data$label_new <- factor(sc@meta.data$label_new, levels = c("myCAF-like","iCAF-like","apCAF-like","PSC","Endothelial"))


## prepare median VAM matrix-----------------
### VAM
source('/work/users/c/h/changfei/R_functions/VAM_Functions.R')
gene.set.collection = createGeneSetCollection_mouse_human (gene.ids=sc@assays$SCT@data@Dimnames[[1]],
                                                           gene.set.collection=gene_list)
sc@active.assay <- "SCT"
sc <- vamForSeurat(sc, 
                   gene.set.collection = gene.set.collection, 
                   center = FALSE, 
                   gamma = TRUE, 
                   sample.cov = FALSE, 
                   return.dist = TRUE)
DefaultAssay(object = sc) = "VAMcdf"

### prepare data to median VAM vs celltype
VAM_score_info <- as.data.frame(t(sc@assays$VAMcdf@data))
VAM_score_info$label_new <-sc@meta.data$label_new

VAM_Median <- VAM_score_info %>%
  group_by(label_new) %>%
  summarize_all(median)%>%
  as.data.frame()

rownames(VAM_Median) <- VAM_Median$label_new
VAM_Median <- VAM_Median[,-1]
VAM_Median_matrix <- as.matrix(t(VAM_Median))

VAM_Median_matrix_zscore <- t(scale(t(VAM_Median_matrix)))
VAM_Median_matrix_zscore[is.na(VAM_Median_matrix_zscore)] <- 0

VAM_Median_matrix_zscore_df <- as.data.frame(VAM_Median_matrix_zscore)

schema_groups <- names(MarkerGene_list_VAM)
VAM_Median_matrix_zscore_df$schema_groups <- case_when(rownames(VAM_Median_matrix_zscore_df)%in%names(MarkerGene_list_VAM[[schema_groups[1]]]) ~ schema_groups[1], 
                                                       rownames(VAM_Median_matrix_zscore_df)%in%names(MarkerGene_list_VAM[[schema_groups[2]]]) ~ schema_groups[2], 
                                                       rownames(VAM_Median_matrix_zscore_df)%in%names(MarkerGene_list_VAM[[schema_groups[3]]]) ~ schema_groups[3], 
                                                       rownames(VAM_Median_matrix_zscore_df)%in%names(MarkerGene_list_VAM[[schema_groups[4]]]) ~ schema_groups[4], 
                                                       rownames(VAM_Median_matrix_zscore_df)%in%names(MarkerGene_list_VAM[[schema_groups[5]]]) ~ schema_groups[5], 
                                                       rownames(VAM_Median_matrix_zscore_df)%in%names(MarkerGene_list_VAM[[schema_groups[6]]]) ~ schema_groups[6], 
                                                       rownames(VAM_Median_matrix_zscore_df)%in%names(MarkerGene_list_VAM[[schema_groups[7]]]) ~ schema_groups[7], 
                                                       rownames(VAM_Median_matrix_zscore_df)%in%names(MarkerGene_list_VAM[[schema_groups[8]]]) ~ schema_groups[8], 
                                                       rownames(VAM_Median_matrix_zscore_df)%in%names(MarkerGene_list_VAM[[schema_groups[9]]]) ~ schema_groups[9], 
                                                       rownames(VAM_Median_matrix_zscore_df)%in%names(MarkerGene_list_VAM[[schema_groups[10]]]) ~ schema_groups[10]
                                                       
)

VAM_Median_matrix_zscore_df$schema_groups <- factor(VAM_Median_matrix_zscore_df$schema_groups,levels = c("SCISSORS","Elyada","Ogawa","Chen","Grunwald","Hwang","Oh","Wang","MS","Maurer" ))
VAM_Median_matrix_zscore_df <- VAM_Median_matrix_zscore_df %>%
  arrange(schema_groups)

### change row orders as Laura suggested
new_row_order <- c("SCISSORS-CAF-top10-myCAF", "SCISSORS-CAF-top10-iCAF", "SCISSORS-CAF-top10-apCAF", "SCISSORS-Stroma-Top10-Perivascular",
                   "Elyada-CAF-myCAF", "Elyada-CAF-iCAF", "Elyada-CAF-mouse-apCAF", "Elyada-Perivascular",
                   "Ogawa-F-stroma","Ogawa-A-stroma","Ogawa-C-stroma",
                   "Chen-cCAF-Top10","Chen-csCAF-Top10","Chen-PSC-Top10",
                   "Grunwald-Reactive-Inter-Top10","Grunwald-Deserted-Top10",
                   "Hwang-Fibroblast-Myofibroblastic-Top10","Hwang-Fibroblast-Neurotropic-Top10","Hwang-Fibroblast-Adhesive-Top10","Hwang-Fibroblast-Immunomodulatory-Top10",
                   "Oh-Stroma-myCAF-Top10","Oh-Stroma-IL11-CAF-Top10","Oh-Stroma-Myocyte-Top10","Oh-Stroma-iCAF-Top10","Oh-Stroma-csCAF-Top10","Oh-Stroma-Schwann-Top10","Oh-Stroma-qPSC-Top10","Oh-Stroma-smPSC-Top10",
                   "Wang-myCAF-Top10","Wang-nonTypicalCAF-Top10","Wang-iCAF-Top10","Wang-meCAF-Top10","Wang-Stellate-1-Top10","Wang-Stellate-2-Top10",
                   "MS-Activated","MS-Normal",
                   "Maurer-ECM","Maurer-Immune")
index_order <- match(rownames(VAM_Median_matrix_zscore_df),new_row_order)
VAM_Median_matrix_zscore_df <- VAM_Median_matrix_zscore_df[order(index_order),]

### revers the row and col, get prepared matrix
VAM_Median_matrix_zscore_HM <- as.matrix(VAM_Median_matrix_zscore_df[,c(1:4)])%>% 
  t(.)


## heatmap------------
### set group color
color_group <- c()
color_group = c(cols25(n=length(levels(VAM_Median_matrix_zscore_df$schema_groups))))
color_group[1] <- "#ff0000"
color_group[2] <- "#1F78C8"
names(color_group)=levels(VAM_Median_matrix_zscore_df$schema_groups)

### change colnames as laura suggested
new_colnames <- colnames(VAM_Median_matrix_zscore_HM) %>% 
  sub("^[^-]+-", "", .)%>% 
  sub("-Top10", "", .)%>% 
  sub("-top10", "", .)%>%
  sub("CAF-", "", .)%>%
  sub("Fibroblast-", "", .)%>%
  sub("Stroma-", "", .)
new_colnames[7] <- "apCAF-mouse"  
new_colnames[1:4] <- c("myCAF-like","iCAF-like","apCAF-like","PSC")

### set bottom annotation
bottom_annotation <- HeatmapAnnotation(
  marker_group = VAM_Median_matrix_zscore_df$schema_groups,
  colname = anno_text(new_colnames, rot = 45, just = "right", offset = unit(1, "npc") - unit(0, "mm")),
  col = list(marker_group = color_group),
  simple_anno_size = unit(1, "mm"),
  show_annotation_name = FALSE,
  show_legend = c(FALSE)
)

p2<-Heatmap(VAM_Median_matrix_zscore_HM,
            name="Elyada Fib",
            column_split = VAM_Median_matrix_zscore_df$schema_groups,# Specify the group labels
            show_column_names = FALSE, 
            row_names_side = "left",
            bottom_annotation = bottom_annotation,
            cluster_rows = FALSE,
            cluster_columns = FALSE)
p2



pdf("/work/users/c/h/changfei/01_CAPER_Paper/00_DeCAF_paper_summary/Fig1/Fig1_EJ/figures/DeCAFpaper_Yehseq_fib_heatmap.pdf",width = 12,height = 3.32)
p1
dev.off()
pdf("/work/users/c/h/changfei/01_CAPER_Paper/00_DeCAF_paper_summary/Fig1/Fig1_EJ/figures/DeCAFpaper_Elyada_fib_heatmap.pdf",width = 12,height = 2.6)
p2
dev.off()


