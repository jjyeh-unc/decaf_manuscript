# load libraries
library(gtable)
library(gridExtra)
library(grid)
library(gplots)
library(openxlsx)
library(tidyverse)
library(reshape2)
library(ggtext)

gc()



#subset infromation of "T.cells.CD8","T.cells.regulatory..Tregs."----------------- 
result_dir_list <- list.dirs("/work/users/c/h/changfei/01_CAPER_Paper/00_DeCAF_paper_summary/Fig6/Fig6_A/data/public_data_11_cibersort_results",full.names = FALSE)
result_dir_list<- result_dir_list[-1]
data_list <- list()

## 11 pubulic datasets cibersort results
for (name in result_dir_list) {
  setwd(paste0("/work/users/c/h/changfei/01_CAPER_Paper/00_DeCAF_paper_summary/Fig6/Fig6_A/data/public_data_11_cibersort_results/",name))
  result_tmp <- read.table("CIBERSORTx_Adjusted.txt",head=TRUE, sep = "\t")
  rownames(result_tmp) <- result_tmp$Mixture
  result_tmp <- result_tmp[,which(colnames(result_tmp) %in% c("T.cells.CD8","T.cells.regulatory..Tregs."))]
  colname_final <- colnames(result_tmp)
  
  setwd(paste0("/work/users/c/h/changfei/CIBERSORTX/public_data_11/data_classified_anddataname"))
  data_tmp <- readRDS(paste0(name,".Rds"))
  restCAF_samplenames <-data_tmp$sampInfo[which(data_tmp$sampInfo$PurISS.final=="iCAF"),]$Tumor.Sample.ID
  permCAF_samplenames <-data_tmp$sampInfo[which(data_tmp$sampInfo$PurISS.final=="myCAF"),]$Tumor.Sample.ID
  
  result_tmp$PurISS.final <- case_when(rownames(result_tmp) %in% permCAF_samplenames ~ "permCAF",
                                       rownames(result_tmp) %in% restCAF_samplenames ~ "restCAF")
  result_tmp$data_name <- name
  data_list[[name]] <- result_tmp
}

##  yehseq data cibersort results
yehseq_result <- read.table("/work/users/c/h/changfei/01_CAPER_Paper/00_DeCAF_paper_summary/Fig6/Fig6_A/data/yehseq_pdac/CIBERSORTx_Adjusted.txt",head=TRUE, sep = "\t") 
rownames(yehseq_result) <- yehseq_result$Mixture
yehseq_result <- yehseq_result[,which(colnames(yehseq_result) %in% c("T.cells.CD8","T.cells.regulatory..Tregs."))]


yehseq_data <- readRDS("/work/users/x/l/xlpeng/yehseq_decaf_freeze/yehseq_pdac_pi.decaf_freeze.rds")
restCAF_samplenames <-yehseq_data$subtypeCall[which(yehseq_data$subtypeCall$DeCAF=="restCAF"),]$sampID
permCAF_samplenames <-yehseq_data$subtypeCall[which(yehseq_data$subtypeCall$DeCAF=="permCAF"),]$sampID
yehseq_result$PurISS.final <- case_when(rownames(yehseq_result) %in% permCAF_samplenames ~ "permCAF",
                                        rownames(yehseq_result) %in% restCAF_samplenames ~ "restCAF")
yehseq_result$data_name <- "YehSeq_Pi"
data_list[["YehSeq_Pi"]] <- yehseq_result

data_combined_df <- do.call(rbind, data_list)

## level data_combined_df
level_dataset <- c("TCGA_PAAD",
                   "CPTAC",
                   "Moffitt_GEO_array",
                   "PACA_AU_seq",
                   "PACA_AU_array",
                   "Linehan",
                   "Puleo_array",
                   "Olive",
                   "Dijk",
                   "Grunwald",
                   "Hayashi",
                   "YehSeq_Pi") 
level_PurISS.final <- c("permCAF","restCAF")
data_combined_df$data_name <- factor(data_combined_df$data_name,levels = level_dataset)
data_combined_df$PurISS.final <- factor(data_combined_df$PurISS.final,levels = level_PurISS.final)

## Log transformation log(+min,2) get log(ratio)
all_values <- c(data_combined_df$T.cells.CD8,data_combined_df$T.cells.regulatory..Tregs.)
all_values <- all_values[which(all_values>0)]
min_value <- min(all_values)

log2_T.cells.CD8 <- log(data_combined_df$T.cells.CD8 + min_value,2)
log2_T.cells.regulatory..Tregs. <- log(data_combined_df$T.cells.regulatory..Tregs. + min_value,2)

data_combined_df$ratio <- log2_T.cells.CD8-log2_T.cells.regulatory..Tregs.
summary(data_combined_df)

# violin plot----------------------------
p3 <- ggplot(data = data_combined_df, 
             aes(x = PurISS.final, y = ratio, color = PurISS.final)) +
  geom_violin(size = 0.3,draw_quantiles = c(0.5)) +
  scale_color_manual(values =  c("permCAF" = "#FE46A5", "restCAF" = "#008080")) +
  facet_wrap(~data_name, nrow = 1,strip.position = "bottom") + 
  geom_signif(comparisons = list(c("permCAF","restCAF")), size = 0.2, textsize = 8, margin_top = 0.075, 
              test = wilcox.test, color = "black", map_signif_level = FALSE, step_increase = 0.1) +
  theme(
    # legend.position = "none",
    legend.text = element_text(size = 15),
    axis.text = element_blank(),
    axis.text.y = element_text(size=25),
    strip.text = element_markdown(angle = 50, size = 25, hjust = 0.8, vjust = 0.8),
    # strip.text = element_text(angle = 0,size = 25, hjust = 1,vjust = 1),
    strip.background = element_blank(),
    # strip.switch.pad.grid = unit(1, "lines"),
    panel.background = element_rect(fill = "transparent"),# Remove background color
    panel.border = element_rect(color = "black", fill = NA, size = 0.2),
    axis.title.y = element_text(size = 25),
    axis.ticks.x = element_blank())+
  ylim(-17, 17) +
  labs(title="",
       y = "Log2(CD8 T/Treg)",
       # y = "CD8 T/Treg",
       x = "", 
       fill = NULL)
p3

# boxplot------------------------------
p4 <- ggplot(data = data_combined_df, 
             aes(x = PurISS.final, y = ratio, color = PurISS.final)) +
  geom_boxplot() +
  scale_color_manual(values =  c("permCAF" = "#FE46A5", "restCAF" = "#008080")) +
  facet_wrap(~data_name, nrow = 1,strip.position = "bottom") + 
  geom_signif(comparisons = list(c("permCAF","restCAF")), size = 0.2, textsize = 8, margin_top = 0.075, 
              test = wilcox.test, color = "black", map_signif_level = FALSE, step_increase = 0.1) +
  theme(
        # legend.position = "none",
        legend.text = element_text(size = 15),
        axis.text = element_blank(),
        axis.text.y = element_text(size=25),
        strip.text = element_markdown(angle = 50, size = 25, hjust = 0.8, vjust = 0.8),
        # strip.text = element_text(angle = 0,size = 25, hjust = 1,vjust = 1),
        strip.background = element_blank(),
        # strip.switch.pad.grid = unit(1, "lines"),
        panel.background = element_rect(fill = "transparent"),# Remove background color
        panel.border = element_rect(color = "black", fill = NA, size = 0.2),
        axis.title.y = element_text(size = 25),
        axis.ticks.x = element_blank())+
  ylim(-17, 17) +
  labs(title="",
       y = "Log2(CD8 T/Treg)",
       # y = "CD8 T/Treg",
       x = "", 
       fill = NULL)
p4


pdf("/work/users/c/h/changfei/CAPER_Paper/00_DeCAF_paper_summary/Fig6/Fig6_B/figures/DeCAFpaper_CD8T_Treg_violin_boxplot_20231214.pdf",width = 30,height = 7)
p3
p4
dev.off()



