# library
library(VAM)           # single cell GSEA
library(dplyr)         # tidy data 
library(Seurat)        # single cell infrastructure
library(ggplot2)       # pretty plots
library(SingleR)       # cell type assignment
library(janitor)       # clean data
library(magrittr)      # pipes
library(SCISSORS)      # our package
library(mixtools)      # Gaussian mixture model estimation
library(ggsignif)      # significance bars
library(patchwork)     # align plots
library(latex2exp)     # LaTeX
library(paletteer)     # color palettes
library(reticulate)    # Python interface
library(kableExtra)    # pretty tables
library(wesanderson)   # more color palettesa
library(SingleCellExperiment)
gc() ## clear the memory

# set dir for saving processed data
setwd("/work/users/c/h/changfei/01_CAPER_Paper/04_scRNA_process/data")


# read raw counts and create seurat object----------------------------------------

raw_counts <- Read10X(data.dir = "/work/users/c/h/changfei/01_CAPER_Paper/04_scRNA_process/Raw_data/filtered_feature_bc_matrix")

Yehseq <- CreateSeuratObject(raw_counts,
                             min.cells = 3,
                             min.features = 200
)

Yehseq@meta.data$sample.name <- case_when(grepl("-1", rownames(Yehseq@meta.data)) ~ "P190429T1", 
                                      grepl("-2", rownames(Yehseq@meta.data)) ~ "P190722T1", 
                                      grepl("-3", rownames(Yehseq@meta.data)) ~ "P191121T2", 
                                      grepl("-4", rownames(Yehseq@meta.data)) ~ "P200220T1",
                                      grepl("-5", rownames(Yehseq@meta.data)) ~ "P220106T1",
                                      grepl("-6", rownames(Yehseq@meta.data)) ~ "P220330T1")
Yehseq@meta.data$sample.type <- case_when(grepl("-1", rownames(Yehseq@meta.data)) ~ "FOLFIRINOX", 
                                      grepl("-2", rownames(Yehseq@meta.data)) ~ "FOLFIRINOX", 
                                      grepl("-3", rownames(Yehseq@meta.data)) ~ "FOLFIRINOX", 
                                      grepl("-4", rownames(Yehseq@meta.data)) ~ "FOLFIRINOX",
                                      grepl("-5", rownames(Yehseq@meta.data)) ~ "gemcitabine-abraxane",
                                      grepl("-6", rownames(Yehseq@meta.data)) ~ "untreated")
saveRDS(Yehseq, file = "Yehseq_sampleinfo.Rds")



# QC----------------------------------------------------------------------------

## calculate percent_MT
# Yehseq <- readRDS("Yehseq_sampleinfo.Rds")
Yehseq[["percent_MT"]] <- PercentageFeatureSet(Yehseq, pattern = "^MT-|^mt-") 
table(Yehseq$sample.name)

##filter normalization and scaling
Yehseq <- subset(Yehseq, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent_MT < 5)
Yehseq <- NormalizeData(Yehseq)
all.genes <- rownames(Yehseq)
Yehseq <- ScaleData(Yehseq, features = all.genes)
saveRDS(Yehseq, file = "Yehseq_QC.Rds")



# first round clustering using SCISSORS-----------------------------------------

# Yehseq <- readRDS("Yehseq_QC.Rds")

## prepare data
Yehseq <- PrepareData(seurat.object = Yehseq, 
                      use.sct = TRUE,
                      n.HVG = 3000,
                      regress.mt = TRUE,
                      regress.cc = FALSE,
                      n.PC = 20,
                      which.dim.reduc = c("tsne", "umap"),
                      initial.resolution = .3,
                      use.parallel = FALSE,
                      random.seed = 629)
DimPlot(Yehseq, reduction = "umap", label=TRUE) + scale_color_paletteer_d("ggthemes::Classic_20")

saveRDS(Yehseq, file = "Yehseq_prepared.Rds")



#singleR annotation-------------------------------------------------------------

# Yehseq <- readRDS("Yehseq0_prepared.Rds")

## Single Cell RNA-seq Reference Data

sc_ref <- readRDS("/work/users/x/l/xlpeng/scRNAseq_reference/singleR_ref_scPDAC.Rds")

Yehseq_preds <- SingleR(
  test = data.frame(Yehseq@assays$RNA@data), 
  ref = sc_ref, 
  labels = sc_ref$label, 
  method = "cluster", 
  clusters = Yehseq$seurat_clusters, 
  de.method = "wilcox")

Yehseq[["SingleR.labels.sc"]] <- Yehseq_preds$labels[match(Yehseq[[]][["seurat_clusters"]], rownames(Yehseq_preds))]

unique(Yehseq@meta.data$SingleR.labels.sc)
Idents(Yehseq) <- Yehseq@meta.data$SingleR.labels.sc
DimPlot(Yehseq, reduction = "umap", group.by = "SingleR.labels.sc", label=TRUE, pt.size = 0.75) 


## Bulk Tissue RNA-seq Reference Data

bulk_ref <- celldex::HumanPrimaryCellAtlasData()
Yehseq_preds <- SingleR(test = data.frame(Yehseq@assays$RNA@data), 
                        ref = bulk_ref, 
                        labels = bulk_ref$label.main, 
                        method = "cluster", 
                        clusters = Yehseq$seurat_clusters, 
                        de.method = "wilcox")
Yehseq[["SingleR.labels.bulk"]] <- Yehseq_preds$labels[match(Yehseq[[]][["seurat_clusters"]], rownames(Yehseq_preds))]
unique(Yehseq@meta.data$SingleR.labels.bulk)
Idents(Yehseq) <- Yehseq@meta.data$SingleR.labels.bulk
DimPlot(Yehseq, reduction = "umap", group.by = "SingleR.labels.bulk", label=TRUE, pt.size = 0.75) 

saveRDS(Yehseq, file = "Yehseq_Annotated.Rds")


## confirmed annotation

Yehseq@meta.data %<>% mutate(
  # more accurate ref. pop. now that I've annotated the immune cells
  label_paper_plot = case_when(seurat_clusters %in% c(0) ~ "Fibroblast-1", 
                               seurat_clusters %in% c(1) ~ "T", 
                               seurat_clusters %in% c(2) ~ "Macrophage",  
                               seurat_clusters %in% c(3) ~ "Epithelial",  
                               seurat_clusters %in% c(4) ~ "Fibroblast-2",  
                               seurat_clusters %in% c(5) ~ "Endothelial", 
                               seurat_clusters %in% c(6) ~ "Fibroblast-3", 
                               seurat_clusters %in% c(7) ~ "NK", 
                               seurat_clusters %in% c(8) ~ "B", 
                               seurat_clusters %in% c(9) ~ "Schwann"
  ))

saveRDS(Yehseq,"UNC_sc_all.rds")



#Fibro reclustering using SCISSORS----------------------------------------------

# Yehseq <- readRDS("UNC_sc_all.rds")

##reclustering
Idents(Yehseq) <- Yehseq@meta.data$seurat_clusters
Yehseq_fibro_recluster <- ReclusterCells(Yehseq,
                        which.clust = c(0,4,6), 
                        merge.clusters = TRUE,  
                        n.HVG = 3000, 
                        n.PC = 20, 
                        resolution.vals = c(.05, .1, .2, .3, .4), 
                        k.vals = c(10, 20, 30, 40), 
                        nn.metric = "cosine",
                        use.parallel = FALSE,
                        use.sct = TRUE,
                        redo.embedding = TRUE)

Yehseq_fibro<- ReduceDimensions (obj = Yehseq_fibro_recluster,
                                 n.PC = 20,
                                 which.algos = c("umap"),
                                 random.seed = 312) 

DimPlot(Yehseq_fibro, reduction = "umap", label=TRUE, pt.size = 0.75)
saveRDS(Yehseq_fibro, file = "Yehseq_fibro.Rds")

##reclustering to get less clusters
k.val <- round(sqrt(ncol(Yehseq_fibro)))
seurat.object <- FindNeighbors(Yehseq_fibro,
                               reduction = "pca",
                               dims = 1:20,
                               k.param = k.val,
                               annoy.metric = "cosine",
                               nn.method = "annoy",
                               verbose = FALSE)
recluster_Yehseq_fibro <- Seurat::FindClusters(seurat.object,
                                                 resolution = 0.2,
                                                 random.seed = 312,
                                                 algorithm = 1,
                                                 verbose = FALSE)
DimPlot(recluster_Yehseq_fibro, reduction = "umap", label=TRUE, pt.size = 0.75)
saveRDS(recluster_Yehseq_fibro, file = "Yehseq_fibro_7clusters.Rds")


## re-lable fib according the VAM plots

# recluster_Yehseq_fibro <- readRDS("Yehseq_fibro_7clusters.Rds")

recluster_Yehseq_fibro@meta.data %<>% mutate(
  label = case_when(seurat_clusters %in% c(1,5) ~ "permCAF", 
                    seurat_clusters %in% c(0,2,6) ~ "restCAF", 
                    seurat_clusters %in% c(4) ~ "apCAF", 
                    seurat_clusters %in% c(3) ~ "Pericyte"))
recluster_Yehseq_fibro@meta.data %<>% mutate(
  label_seperate = case_when(seurat_clusters %in% c(1) ~ "permCAF-1", 
                             seurat_clusters %in% c(5) ~ "permCAF-2", 
                             seurat_clusters %in% c(0) ~ "restCAF-1", 
                             seurat_clusters %in% c(2) ~ "restCAF-2", 
                             seurat_clusters %in% c(6) ~ "restCAF-3", 
                             seurat_clusters %in% c(4) ~ "apCAF", 
                             seurat_clusters %in% c(3) ~ "Pericyte"))

recluster_Yehseq_fibro@meta.data$label <- factor(recluster_Yehseq_fibro@meta.data$label,levels = c("permCAF","restCAF","apCAF","Pericyte"))
recluster_Yehseq_fibro@meta.data$label_seperate <- factor(recluster_Yehseq_fibro@meta.data$label_seperate,
                                                          levels = c("permCAF-1","permCAF-2","restCAF-1","restCAF-2","restCAF-3","apCAF","Pericyte"))
saveRDS(recluster_Yehseq_fibro, file = "Yehseq_fibro_7clusters_labeled.Rds")

saveRDS(recluster_Yehseq_fibro, file = "YUNC_sc_CAF.rds")


