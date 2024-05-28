# load libraries
library(gtable)
library(gridExtra)
library(grid)
library(gplots)
library(openxlsx)
library(tidyverse)
library(reshape2)
gc()

# public data cibersort results-------------------------------
# get folder names of public data cibersort results
result_dir_list <- list.dirs("/work/users/c/h/changfei/01_CAPER_Paper/00_DeCAF_paper_summary/Fig6/Fig6_A/data/public_data_11_cibersort_results",full.names = FALSE)
result_dir_list<- result_dir_list[-1]

# calculate p_value and logFC
p_value_list <- list()
log_fold_list <- list()
colname_final <- c()
iCAF_mean_list <- list()
myCAF_mean_list <- list()
for (name in result_dir_list) {
  
  setwd(paste0("/work/users/c/h/changfei/01_CAPER_Paper/00_DeCAF_paper_summary/Fig6/Fig6_A/data/public_data_11_cibersort_results/",name))
  result_tmp <- read.table("CIBERSORTx_Adjusted.txt",head=TRUE, sep = "\t")
  rownames(result_tmp) <- result_tmp$Mixture
  result_tmp <- result_tmp[,-which(colnames(result_tmp) %in% c("Mixture","P.value","Correlation","RMSE"))]
  result_tmp <- t(result_tmp)
  colname_final <- colnames(result_tmp)
  
  setwd(paste0("/work/users/c/h/changfei/CIBERSORTX/public_data_11/data_classified_anddataname"))
  data_tmp <- readRDS(paste0(name,".Rds"))
  iCAF_samplenames <-data_tmp$sampInfo[which(data_tmp$sampInfo$PurISS.final=="iCAF"),]$Tumor.Sample.ID
  myCAF_samplenames <-data_tmp$sampInfo[which(data_tmp$sampInfo$PurISS.final=="myCAF"),]$Tumor.Sample.ID
  
  iCAF_result <- result_tmp[,which(colnames(result_tmp) %in% iCAF_samplenames)]
  myCAF_result <- result_tmp[,which(colnames(result_tmp) %in% myCAF_samplenames)]
  
  # p-value
  p_value_tmp <- c()
  for (i in 1:nrow(iCAF_result)) {
    p_value_tmp[i] <- wilcox.test(iCAF_result[i, ], myCAF_result[i, ])$p.value
  }
  # replace NA with 1
  p_value_tmp <- ifelse(is.na(p_value_tmp), 1, p_value_tmp)
  p_value_list[[name]] <- p_value_tmp
  
  # mean, logfold
  iCAF_mean <- rowMeans(iCAF_result)
  myCAF_mean <- rowMeans(myCAF_result)
  iCAF_mean_list[[name]] <- iCAF_mean
  myCAF_mean_list[[name]] <- myCAF_mean
}

myCAF_mean_df <- as.data.frame(myCAF_mean_list)
iCAF_mean_df <- as.data.frame(iCAF_mean_list)

p_value_df <- as.data.frame(p_value_list)
summary(p_value_df)


# yehseq data cibersort results-------------------------------
yehseq_result <- read.table("/work/users/c/h/changfei/01_CAPER_Paper/00_DeCAF_paper_summary/Fig6/Fig6_A/data/yehseq_pdac/CIBERSORTx_Adjusted.txt",head=TRUE, sep = "\t") 
rownames(yehseq_result) <- yehseq_result$Mixture
yehseq_result <- yehseq_result[,-which(colnames(yehseq_result) %in% c("Mixture","P.value","Correlation","RMSE"))]
yehseq_result <- t(yehseq_result)
colname_final <- colnames(yehseq_result)

yehseq_data <- readRDS("/work/users/x/l/xlpeng/yehseq_decaf_freeze/yehseq_pdac_pi.decaf_freeze.rds")
yehseq_data$subtypeCall
iCAF_samplenames <-yehseq_data$subtypeCall[which(yehseq_data$subtypeCall$DeCAF=="restCAF"),]$sampID
myCAF_samplenames <-yehseq_data$subtypeCall[which(yehseq_data$subtypeCall$DeCAF=="permCAF"),]$sampID

iCAF_yehseq_result <- yehseq_result[,which(colnames(yehseq_result) %in% iCAF_samplenames)]
myCAF_yehseq_result <- yehseq_result[,which(colnames(yehseq_result) %in% myCAF_samplenames)]

## get yehseq data p-value
p_value_tmp <- c()
for (i in 1:nrow(iCAF_yehseq_result)) {
  p_value_tmp[i] <- wilcox.test(iCAF_yehseq_result[i, ], myCAF_yehseq_result[i, ])$p.value
}

p_value_tmp <- ifelse(is.na(p_value_tmp), 1, p_value_tmp)

## get yehseq data mean
iCAF_mean <- rowMeans(iCAF_yehseq_result)
myCAF_mean <- rowMeans(myCAF_yehseq_result)

#combine mean and get logfold
iCAF_mean_df_combine  <- data.frame("YehSeq_Pi" = iCAF_mean,iCAF_mean_df)
myCAF_mean_df_combine  <- data.frame("YehSeq_Pi" = myCAF_mean,myCAF_mean_df)

all_values <- c(unlist(myCAF_mean_df_combine),unlist(iCAF_mean_df_combine))
all_values <- all_values[which(all_values>0)]
min_value <- min(all_values)
log_fold_df_combine <- log(myCAF_mean_df_combine+min_value,2)-log( iCAF_mean_df_combine+min_value,2)

#combine p-value
p_value_df_combine  <- data.frame("YehSeq_Pi" = p_value_tmp,p_value_df)
rownames(p_value_df_combine) <- rownames(log_fold_df_combine)

# save result to excel
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
ws <- openxlsx::createWorkbook()

### sheet1 permCAF_samples_mean
sheetName = 1
openxlsx::addWorksheet(ws, sheetName = "permCAF_samples_mean")
openxlsx::writeData(ws, sheet = sheetName, x = myCAF_mean_df_combine[,level_dataset],rowNames = TRUE)

### sheet2 restCAF_samples_mean
sheetName = 2
openxlsx::addWorksheet(ws, sheetName = "restCAF_samples_mean")
openxlsx::writeData(ws, sheet = sheetName, x = iCAF_mean_df_combine[,level_dataset],rowNames = TRUE)

### sheet3 logFC
sheetName = 3
openxlsx::addWorksheet(ws, sheetName = "logFC_permCAFmean-restCAFmean")
openxlsx::writeData(ws, sheet = sheetName, x = log_fold_df_combine[,level_dataset],rowNames = TRUE)

### sheet4 pvalue
sheetName = 4
openxlsx::addWorksheet(ws, sheetName = "pvalue_wilcox.test")
openxlsx::writeData(ws, sheet = sheetName, x = p_value_df_combine[,level_dataset],rowNames = TRUE)

openxlsx::saveWorkbook(ws, file = "/work/users/c/h/changfei/01_CAPER_Paper/00_DeCAF_paper_summary/Fig6/Fig6_A/figures/DeCAFpaper_Cibersort_Dotplot_results.xlsx")



#sort rows by numbers of significant------------------------------------------------------------------------ 
##get icaf p_value
log_fold_df_combine_icaf <- apply(log_fold_df_combine, 2, function(x) ifelse(x < 0, 1, 0))
p_value_df_combine_icaf <- apply(p_value_df_combine, 2, function(x) ifelse(x < 0.05, 1, 0))

p_value_df_icaf <- log_fold_df_combine_icaf*p_value_df_combine_icaf
p_value_df_icaf<- as.data.frame(p_value_df_icaf)
p_value_df_icaf$sig_count_icaf <- rowSums(p_value_df_icaf)

##get mycaf p_value
log_fold_df_combine_mycaf <- apply(log_fold_df_combine, 2, function(x) ifelse(x > 0, 1, 0))
p_value_df_combine_mycaf <- apply(p_value_df_combine, 2, function(x) ifelse(x < 0.05, 1, 0))

p_value_df_mycaf <- log_fold_df_combine_mycaf*p_value_df_combine_mycaf
p_value_df_mycaf<- as.data.frame(p_value_df_mycaf)
p_value_df_mycaf$sig_count_mycaf <- rowSums(p_value_df_mycaf)

##combine and sort by sig count
sort_df <- data_frame(
  "celltype" = rownames(p_value_df_mycaf),
  "sig_count_icaf" = p_value_df_icaf$sig_count_icaf,
  "sig_count_mycaf" = p_value_df_mycaf$sig_count_mycaf
)
sort_df_sorted <- sort_df[order(sort_df$sig_count_icaf, -sort_df$sig_count_mycaf), ]
sort_df_sorted <- as.data.frame(sort_df_sorted)


# sort p_value_df_melt log_fold_df_melt by the above sorted celltype order 
new_order <- match(sort_df_sorted$celltype, rownames(p_value_df_combine))
p_value_df_combine <- p_value_df_combine[new_order,]
new_order <- match(sort_df_sorted$celltype, rownames(log_fold_df_combine))
log_fold_df_combine <- log_fold_df_combine[new_order,]


# melt data for ploting
celltype <- rownames(p_value_df_combine)
p_value_df_combine_melt <- data.frame(celltype,p_value_df_combine)
p_value_df_melt <- melt(p_value_df_combine_melt,id.vars = "celltype", variable.name = "data_name", value.name = "p_value")
colnames(p_value_df_melt) <- c("celltype","data_name","p_value")

log_fold_df_combine_melt <- data.frame(celltype,log_fold_df_combine)
log_fold_df_melt <- melt(log_fold_df_combine_melt,id.vars = "celltype", variable.name = "data_name", value.name = "lfc_value")
colnames(log_fold_df_melt) <- c("celltype","data_name","lfc_value")
all(p_value_df_melt$celltype==log_fold_df_melt$celltype)
all(p_value_df_melt$data_name==log_fold_df_melt$data_name)


# combine p-value and lfc
lfc_value <- log_fold_df_melt$lfc_value
p_lfc_df <-data.frame(p_value_df_melt,"lfc_value"=lfc_value) 
summary(p_lfc_df)
##set level
p_lfc_df$celltype <- factor(p_lfc_df$celltype,levels = sort_df_sorted$celltype )
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
p_lfc_df$data_name <- factor(p_lfc_df$data_name,levels = level_dataset)



#adjust data p-value  p---1/p, 
p_lfc_df_adj <- p_lfc_df
p_lfc_df_adj$p_value_inverse <- ifelse(p_lfc_df_adj$p_value>0.05,0,1/(p_lfc_df_adj$p_value ))
#adjust data lfc: lfc>3  ===3  lfc<-3 ===-3
p_lfc_df_adj$LogFC <- ifelse(p_lfc_df_adj$lfc_value > 2, 2, 
                             ifelse(p_lfc_df_adj$lfc_value < -2, -2, p_lfc_df_adj$lfc_value))
summary(p_lfc_df_adj)

#group p
breaks <- c(0, 20, 100, 1000, Inf)
sizes <- c(">0.05", "0.05~0.01", "0.01~0.001", "<0.001") ##according how p_lfc_df_adj$p_value calculated
p_lfc_df_adj$p_value_range <- cut(p_lfc_df_adj$p_value_inverse, breaks = breaks, labels = sizes)

ggplot(p_lfc_df_adj, aes(x = data_name, y = celltype)) +
  geom_point(aes(size = p_value_range, color = LogFC),stroke = 1) +
  scale_color_gradient2(low ="#008080" , high = "#FE46A5", mid = "white", midpoint = 0)+ 
  scale_size_manual(values = c(">0.05"=0,"0.05~0.01" = 4, "0.01~0.001" = 8, "<0.001" = 12)) +
  # scale_size_continuous(range = c(0, 10),breaks = c(20,100,1000))+
  labs(x = "datasets", y = "celltype_LM22") +
  theme(panel.background = element_rect(fill = "grey"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 15))
# theme_minimal() 


pdf("/work/users/c/h/changfei/01_CAPER_Paper/00_DeCAF_paper_summary/Fig6/Fig6_A/figures/cibersort_lm22_dotplot_20231107.pdf",width = 15,height = 10)
ggplot(p_lfc_df_adj, aes(x = data_name, y = celltype)) +
  geom_point(aes(size = p_value_range, color = LogFC),stroke = 1) +
  scale_color_gradient2(low ="#008080" , high = "#FE46A5", mid = "white", midpoint = 0)+ 
  scale_size_manual(values = c(">0.05"=0,"0.05~0.01" = 4, "0.01~0.001" = 8, "<0.001" = 11)) +
    # scale_size_manual(values = c(">0.05"=0,"0.05~0.01" = 5, "0.01~0.001" = 10, "<0.001" = 14)) +
  # scale_size_continuous(range = c(0, 10),breaks = c(20,100,1000))+
  # labs(x = "datasets", y = "celltype_LM22") +
  labs(x = "", y = "") +
  theme(panel.background = element_rect(fill = "grey"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 20),
        axis.text.y = element_text(size = 20))
# theme_minimal() 
dev.off()


