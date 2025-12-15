library(tidyverse)  
library(stringr)    
library(dplyr)      
library(latex2exp) 
library(openxlsx)
library(readxl)
library(ComplexHeatmap)
library(circlize)
gc()

# read 6 public paper markers----------------------------------
folder_path <- "/work/users/c/h/changfei/Other_Markers/MarkerGene_list_public"
file_names <- list.files(path = folder_path, full.names = FALSE)
Other_markers <- data.frame(
  "data_name" = c("Ogawa","Chen","Grunwald","Helms","Hwang","Oh","Wang"),
  "file_names" = file_names
)

MarkerGene_list <- list()

for (i in 1:dim(Other_markers)[1]) {
  MarkerGene_list_tmp <- readRDS(paste0(folder_path,"/",Other_markers$file_names[i]))
  MarkerGene_list[[paste0(Other_markers$data_name[i])]] <- MarkerGene_list_tmp

}
names(MarkerGene_list)
# other markers including:-----------------------------------------
## Moffitt stroma: activated, normal;Maurer: ECM-rich, Immune-rich;Scissors: iCAF-like , myCAF-like, apCAF-like, PSC . Top10;
##readElyada: iCAF, myCAF, perivascular, mouse  apCAF

load("/work/users/c/h/changfei/Other_Markers/Our_markerGenes/cmbSubtypes.230508.RData")
schemaList
subtypeGeneList[5]

##MS
Cell_type <- colnames(subtypeGeneList[[5]])[2:3]
marker_info <- subtypeGeneList[[5]]
MarkerGene_list_tmp <- list()
for (i in 1:length(Cell_type)) {
  Gene_marker_set_name <- paste0("MS-",Cell_type[i])
  col_index <- which(colnames(subtypeGeneList[[5]])==Cell_type[i])
  MarkerGene_list_tmp[[Gene_marker_set_name]] <- marker_info[which(marker_info[,col_index]=="TRUE"),]$geneSymbol
}
MarkerGene_list[["MS"]] <- MarkerGene_list_tmp

##Maurer
marker_info <- subtypeGeneList[[11]]
Cell_type <- colnames(marker_info)[2:3]
MarkerGene_list_tmp <- list()
for (i in 1:length(Cell_type)) {
  Gene_marker_set_name <- paste0("Maurer-",Cell_type[i])
  col_index <- which(colnames(marker_info)==Cell_type[i])
  MarkerGene_list_tmp[[Gene_marker_set_name]] <- marker_info[which(marker_info[,col_index]=="TRUE"),]$geneSymbol
}
MarkerGene_list[["Maurer"]] <- MarkerGene_list_tmp

##Elyada-------
schemaList
MarkerGene_list_tmp <- list()
###iCAF,myCAF
marker_info <- subtypeGeneList[[14]]
Cell_type <- colnames(marker_info)[2:3]
for (i in 1:length(Cell_type)) {
  Gene_marker_set_name <- paste0("Elyada-CAF-",Cell_type[i])
  col_index <- which(colnames(marker_info)==Cell_type[i])
  MarkerGene_list_tmp[[Gene_marker_set_name]] <- marker_info[which(marker_info[,col_index]=="TRUE"),]$geneSymbol
}
###apCAF
marker_info <- subtypeGeneList[[15]]
Cell_type <- colnames(marker_info)[2:3]
Gene_marker_set_name <- paste0("Elyada-CAF-mouse-",Cell_type[1])
col_index <- which(colnames(marker_info)==Cell_type[1])
MarkerGene_list_tmp[[Gene_marker_set_name]] <- marker_info[which(marker_info[,col_index]=="TRUE"),]$geneSymbol
###Perivascular
gene_list <- readRDS("/work/users/c/h/changfei/Other_Markers/Our_markerGenes/Genes_Elyada.rds")
MarkerGene_list_tmp[["Elyada-Perivascular"]] <- gene_list$`Elyada-Perivascular`
###combine to the final gene list
MarkerGene_list[["Elyada"]] <- MarkerGene_list_tmp

## Scissors----
MarkerGene_list_tmp <- list()
###PSC
gene_list <- readRDS("/work/users/c/h/changfei/00_SCISSORS/data/VAM_Markers_version5th/All/combined/VAM_Scissors_Allmarkers_All_combined.Rds")
names(gene_list)
MarkerGene_list_tmp[["SCISSORS-proCAF"]] <- gene_list$`SCISSORS-CAF-myCAF`
MarkerGene_list_tmp[["SCISSORS-restCAF"]] <- gene_list$`SCISSORS-CAF-iCAF`
MarkerGene_list_tmp[["SCISSORS-apCAF"]] <- gene_list$`SCISSORS-CAF-apCAF`
MarkerGene_list_tmp[["SCISSORS-Perivascular"]] <- gene_list$`SCISSORS-Stroma-Perivascular`
###combine to the final gene list
MarkerGene_list[["SCISSORS"]] <- MarkerGene_list_tmp

## DeCAF genes----
gene_list_DeCAF <- readRDS("/work/users/c/h/changfei/Other_Markers/VAM_gene_list/VAM_Gene_list(DeCAF_Elyada_CAF_fromsubtypeSchema_231026).Rds")
gene_list_DeCAF <- gene_list_DeCAF[c(1,2)]
MarkerGene_list[["DeCAF"]] <- gene_list_DeCAF

## selection: Ogawa, chen, gruanwald: all;HWang: fibroblast; Oh: stroma; Wang: all--------------------
MarkerGene_list <- MarkerGene_list
MarkerGene_list$Hwang <- MarkerGene_list$Hwang[15:18]
MarkerGene_list$Oh <- MarkerGene_list$Oh[27:34]

## exclude Helms-----------------------------
names(MarkerGene_list)
MarkerGene_list <- MarkerGene_list[-4]
names(MarkerGene_list)
MarkerGene_list$SCISSORS

##transfer all gene markers into one list------------------------------------
MarkerGene_list_tmp <- list()
for (name in names(MarkerGene_list)) {
  sub_list <- MarkerGene_list[[name]]
  MarkerGene_list_tmp <- c(MarkerGene_list_tmp, sub_list)
}
names(MarkerGene_list)
names(MarkerGene_list_tmp)



# calculate number of overlapping------------------------
## Initialize an empty data frame for results
percentage_df <- data.frame(matrix(NA, nrow = length(MarkerGene_list_tmp), ncol = length(MarkerGene_list_tmp)))
colnames(percentage_df) <- names(MarkerGene_list_tmp)
rownames(percentage_df) <- names(MarkerGene_list_tmp)

## Calculate overlapping percentage based on the minimum length of vectors
for (i in seq_along(MarkerGene_list_tmp)) {
  for (j in seq_along(MarkerGene_list_tmp)) {
    overlap_numb <- sum(MarkerGene_list_tmp[[j]] %in% MarkerGene_list_tmp[[i]])
    percentage_df[i, j] <- overlap_numb
  }
}

## add dataname to the results
schema_groups <- names(MarkerGene_list)
percentage_df$schema_groups <- case_when(rownames(percentage_df)%in%names(MarkerGene_list[[schema_groups[1]]]) ~ schema_groups[1], 
                                         rownames(percentage_df)%in%names(MarkerGene_list[[schema_groups[2]]]) ~ schema_groups[2], 
                                         rownames(percentage_df)%in%names(MarkerGene_list[[schema_groups[3]]]) ~ schema_groups[3], 
                                         rownames(percentage_df)%in%names(MarkerGene_list[[schema_groups[4]]]) ~ schema_groups[4], 
                                         rownames(percentage_df)%in%names(MarkerGene_list[[schema_groups[5]]]) ~ schema_groups[5], 
                                         rownames(percentage_df)%in%names(MarkerGene_list[[schema_groups[6]]]) ~ schema_groups[6], 
                                         rownames(percentage_df)%in%names(MarkerGene_list[[schema_groups[7]]]) ~ schema_groups[7], 
                                         rownames(percentage_df)%in%names(MarkerGene_list[[schema_groups[8]]]) ~ schema_groups[8], 
                                         rownames(percentage_df)%in%names(MarkerGene_list[[schema_groups[9]]]) ~ schema_groups[9],
                                         rownames(percentage_df)%in%names(MarkerGene_list[[schema_groups[10]]]) ~ schema_groups[10],
                                         rownames(percentage_df)%in%names(MarkerGene_list[[schema_groups[11]]]) ~ schema_groups[11]
                                                       
)


percentage_df$schema_groups <- factor(percentage_df$schema_groups,levels = c("DeCAF","SCISSORS","Elyada","Ogawa","Chen","Grunwald","Hwang","Oh","Wang","MS","Maurer" ))
percentage_df <- percentage_df %>%
  arrange(schema_groups)

### change row orders as Laura suggested
new_row_order <- c("DeCAF-proCAF","DeCAF-restCAF ",
                   "SCISSORS-proCAF", "SCISSORS-restCAF", "SCISSORS-apCAF", "SCISSORS-Perivascular",
                   "Elyada-CAF-myCAF", "Elyada-CAF-iCAF", "Elyada-CAF-mouse-apCAF", "Elyada-Perivascular",
                   "Ogawa-F-stroma","Ogawa-A-stroma","Ogawa-C-stroma",
                   "Chen-cCAF","Chen-csCAF","Chen-PSC",
                   "Grunwald-Reactive-Inter","Grunwald-Deserted",
                   "Hwang-Fibroblast-Myofibroblastic","Hwang-Fibroblast-Neurotropic","Hwang-Fibroblast-Adhesive","Hwang-Fibroblast-Immunomodulatory",
                   "Oh-Stroma-myCAF","Oh-Stroma-IL11-CAF","Oh-Stroma-Myocyte","Oh-Stroma-iCAF","Oh-Stroma-csCAF","Oh-Stroma-Schwann","Oh-Stroma-qPSC","Oh-Stroma-smPSC",
                   "Wang-myCAF","Wang-nonTypicalCAF","Wang-iCAF","Wang-meCAF","Wang-Stellate-1","Wang-Stellate-2",
                   "MS-Activated","MS-Normal",
                   "Maurer-ECM","Maurer-Immune")
index_order <- match(rownames(percentage_df),new_row_order)
percentage_df <- percentage_df[order(index_order),]

index_order <- match(colnames(percentage_df),new_row_order)
percentage_df <- percentage_df[,order(index_order)]


colnames(percentage_df)
### revers the row and col, get prepared matrix
percentage_HM <- as.matrix(percentage_df[,c(1:40)])%>% 
  t(.)

library("pals")
# heatmap-----------------------------
# set group color
color_group <- c()
color_group = c(cols25(n=length(levels(percentage_df$schema_groups))))
color_group[1] <- "#ff0000"
color_group[2] <- "#1F78C8"
names(color_group)=levels(percentage_df$schema_groups)

### change colnames as laura suggested
new_colnames <- colnames(percentage_HM) %>% 
  sub("^[^-]+-", "", .)%>% 
  sub("-Top10", "", .)%>% 
  sub("-top10", "", .)%>%
  sub("CAF-", "", .)%>%
  sub("Fibroblast-", "", .)%>%
  sub("Stroma-", "", .)
new_colnames[7] <- "apCAF-mouse" 

### change colnames as laura suggested
new_row <- rownames(percentage_HM) %>% 
  sub("^[^-]+-", "", .)%>% 
  sub("-Top10", "", .)%>% 
  sub("-top10", "", .)%>%
  sub("CAF-", "", .)%>%
  sub("Fibroblast-", "", .)%>%
  sub("Stroma-", "", .)
new_row[7] <- "apCAF-mouse"  

### set bottom annotation
bottom_annotation <- HeatmapAnnotation(
  marker_group = percentage_df$schema_groups,
  colname = anno_text(new_colnames, rot = 45, just = "right", offset = unit(1, "npc") - unit(0, "mm")),
  col = list(marker_group = color_group),
  simple_anno_size = unit(1, "mm"),
  show_annotation_name = FALSE,
  show_legend = c(FALSE)
)
left_annotation <- rowAnnotation(
  marker_group = percentage_df$schema_groups,
  col = list(marker_group = color_group),
  simple_anno_size = unit(1, "mm"),
  show_annotation_name = FALSE,
  show_legend = c(FALSE)
)

color_gradient <- colorRamp2(c(0, 100), c("gray90", "red"))
p1<-Heatmap(percentage_HM,
            name=" ",
            col = color_gradient,
            split = percentage_df$schema_groups,
            column_split = percentage_df$schema_groups,# Specify the group labels
            show_column_names = FALSE, 
            row_labels = new_row,
            row_names_side = "left",
            row_title_side = "right" ,
            row_title = " ",
            # row_title_gp = gpar(margin = unit(c(0, 20, 0, 10), "mm")),
            bottom_annotation = bottom_annotation,
            left_annotation  = left_annotation,
            cluster_rows = FALSE,
            cluster_columns = FALSE,
            cell_fun = function(j, i, x, y, width, height, fill) {
              grid.text(percentage_HM[i, j], x, y, gp = gpar(fontsize = 8))
            }

)

p1
pdf("/work/users/c/h/changfei/01_CAPER_Paper/00_DeCAF_paper_summary/FigS1/figures/DeCAFpaper_heatmap_overlaping_genes-2th.pdf",width = 13,height = 12)
p1
dev.off()

## save the df to excel
rownames(percentage_df) <- new_row_order
saveRDS(percentage_df,"/work/users/c/h/changfei/CAPER_Paper/00_DeCAF_paper_summary/FigS1/data/DeCAFpaper_heatmap_overlaping_genes.Rds" )

