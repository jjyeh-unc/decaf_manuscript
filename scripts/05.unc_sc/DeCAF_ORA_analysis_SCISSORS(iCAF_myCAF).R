library(openxlsx)
library(DESeq2)
library("tximport")
library("tximportData")
library(stringr)
library("plyr")
library(dplyr)
library(enrichplot)  ##load befor clusterprofiler
library(clusterProfiler)
library(ggplot2)
library(ggrepel)
library(EnhancedVolcano)
library(ComplexHeatmap)
library(DOSE)
library(venn)
gc()

# get myCAF iCAF markers
CAF_markers <- readRDS("/work/users/c/h/changfei/00_SCISSORS/Markers/Markers_refind_allmarkers_version5th/rds_files/all_markers/CAF_bed_markers.Rds")
table(CAF_markers$cluster)
myCAF <- CAF_markers[which(CAF_markers$cluster == "myCAF"),]$gene
iCAF <- CAF_markers[which(CAF_markers$cluster == "iCAF"),]$gene
gene_list <- list(
  myCAF = myCAF,
  iCAF = iCAF
)

# myCAF iCAF_regulated analysis------------------------------------------------------------------
dir_data <- "/work/users/c/h/changfei/01_CAPER_Paper/00_DeCAF_paper_summary/Fig1/Fig1_D/data"
setwd(dir_data)

organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)

kegg_organism = "hsa"

for (i in 1:length(gene_list)) {
  gene_list_ORA <- gene_list[[i]]
 ## ORA_GO analysis
  ego <- enrichGO(gene = gene_list_ORA,
                        OrgDb = organism,
                        keyType = 'SYMBOL', 
                        ont = "ALL",
                        pAdjustMethod = "BH", 
                        pvalueCutoff  = 0.05, 
                        readable      = TRUE ) 
  # saveRDS(ego, file = paste0("ORA_GO_", names(gene_list[i]), ".Rds"))
  write.xlsx(ego, "ORA_analysis_CAPER_SCISSORS.xls",
             sheetName = paste0("ORA_GO_", names(gene_list[i])), append = TRUE)
  ## ORA_KEGG analysis
  gene_list_KEGG <- bitr(gene_list_ORA, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)
  enrich_KEGG <- enrichKEGG(gene = gene_list_KEGG$ENTREZID ,
                                  organism = kegg_organism,
                                  pvalueCutoff = 0.05,
                                  keyType = "kegg")
  enrich_KEGG <- setReadable(enrich_KEGG, 'org.Hs.eg.db', 'ENTREZID')%>%
    filter(p.adjust < 0.05,qvalue<0.2)       ##there is a but in enrichKEGG() all the related pathways will be outputed
  # saveRDS(enrich_KEGG, file = paste0("ORA_KEGG_",names(gene_list[i]), ".Rds"))
  write.xlsx(enrich_KEGG, "ORA_analysis_CAPER_SCISSORS.xls",
             sheetName =  paste0("ORA_KEGG_", names(gene_list[i])), append = TRUE)
}


# ORA figure--------------------------------------------------------------------
library(readxl)

data_go <- read_excel("/work/users/c/h/changfei/01_CAPER_Paper/00_DeCAF_paper_summary/Fig1/Fig1_D/data/ORA_dotplot2.xls",sheet=1)
data_kegg <- read_excel("/work/users/c/h/changfei/01_CAPER_Paper/00_DeCAF_paper_summary/Fig1/Fig1_D/data/ORA_dotplot2.xls",sheet=2)

# deal with the data 
data_go$ONTOLOGY <- paste0("GO_",data_go$ONTOLOGY)
data_kegg$type <- "iCAF"
data_kegg$ONTOLOGY <- "KEGG"

# sub set useful information
data_go_sub <- data_go[,c("ONTOLOGY","Description","GeneRatio","p.adjust","type")]
data_kegg_sub <- data_kegg[,c("ONTOLOGY","Description","GeneRatio","p.adjust","type")]



# combine the two data and arrange them
data_combine <- rbind(data_go_sub,data_kegg_sub)
data_combine$ONTOLOGY <- factor(data_combine$ONTOLOGY,levels = rev(c("GO_BP","GO_MF","GO_CC","KEGG")))
data_combine <- data_combine %>%
  arrange(ONTOLOGY,GeneRatio)

data_combine$Description_ONTOLOGY <- paste0(data_combine$ONTOLOGY,":",data_combine$Description)
data_combine$Description_ONTOLOGY <- factor(data_combine$Description_ONTOLOGY,unique(data_combine$Description_ONTOLOGY))
data_combine$GeneRatio_sign <- ifelse(data_combine$type=="iCAF",data_combine$GeneRatio,-(data_combine$GeneRatio))


data_combine$padj_inverse <- 1/data_combine$p.adjust 
breaks <- c(0, 20, 100, 1000, Inf)
sizes <- c(">0.05", "0.05~0.01", "0.01~0.001", "<0.001") ##according how p_lfc_df_adj$p_value calculated
data_combine$padj_range <- cut(data_combine$padj_inverse, breaks = breaks, labels = sizes)

ggplot(data_combine, aes(x = GeneRatio_sign, y = Description_ONTOLOGY)) +
  geom_point(aes(size = padj_range,color = type),stroke = 1.5)+
  scale_color_manual(values = c("iCAF" = "turquoise4", "myCAF" = "violetred1"))+
  geom_vline(xintercept = 0, linetype = "dashed", color = "black",size = 0.5) +
  labs(x = "Â±GeneRatio", y = "") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15))

p <- ggplot(data_combine, aes(x = GeneRatio_sign, y = Description_ONTOLOGY)) +
  geom_point(aes(size = padj_range,color = type),stroke = 1.5)+
  scale_color_manual(values = c("iCAF" = "turquoise4", "myCAF" = "violetred1"))+
  geom_vline(xintercept = 0, linetype = "dashed", color = "black",size = 0.5) +
  labs(x = "Â±GeneRatio", y = "") +
  theme(axis.text.x = element_text(size = 17),
        axis.text.y = element_text(size = 17),
        axis.title.x = element_text(size = 17),
        legend.text = element_text(size = 17),
        legend.title = element_text(size = 17))

pdf("/work/users/c/h/changfei/01_CAPER_Paper/00_DeCAF_paper_summary/Fig1/Fig1_D/figures/DeCAFpaper_SCISSORS_CAF_ORA.pdf",width = 12,height = 8)
p
dev.off()

