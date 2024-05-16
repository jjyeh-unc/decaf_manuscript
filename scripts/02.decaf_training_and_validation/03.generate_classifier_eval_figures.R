#####################################################################################

##########     Classifier Evaluation Figure Generating Script             ############

####################################################################################




######################################################################################
############## Step 1: Create Functions for Generating Figures ########################
######################################################################################




# load libraries
library(pheatmap)
library(ggplot2)
library(scales)
library(caret)
library(pROC)
library(kableExtra)
library(nsROC)
library(ROCR)

## Colors for DeCAF
DeCAFList = c("permCAF","restCAF")
DeCAFCol = c("violetred1","turquoise4")

############################################################

## Function for Heatmap
Plot_heatmap_function <- function(subtype_data, gene_data, DeCAF_data, classifier, training, CAF_type = c("K2", "K4"), schemaList, include_PurIST = F){
  
  ## List of Samples To Include
  if(CAF_type == "K4"){
    samples_subset = !(subtype_data$dataSet %in% training) & !(is.na(subtype_data$SCISSORS_CAF_ALL))
  }
  
  if(CAF_type == "K2"){
    samples_subset = !(subtype_data$dataSet %in% training) & !(is.na(subtype_data$SCISSORS_CAF))
  }
  
  
  ## Transform TSP's to ranks
  TSPgeneMat = apply(gene_data[,samples_subset], 2, rank)
  
  ## Transfrom Calls to Ordered Factor
  DeCAF_data$DeCAF_graded = factor(DeCAF_data$DeCAF_graded, levels = c("Strong restCAF", "Likely restCAF", "Lean restCAF",
                                                                                           "Lean permCAF", "Likely permCAF", "Strong permCAF"))
  ## Make ColorScapr
  colscale = colorRampPalette(c("turquoise4","white", "violetred1"))(6)
  colscale2 = colorRampPalette(c("turquoise4","white", "violetred1"))(5)
  colscale3 = colorRampPalette(c("turquoise4","white", "violetred1"))(2)
  colscale4 = c(colscale2, "black")
  
  # Create cleaned up study labels
  Study = subtype_data$dataSet[samples_subset]
  
  
  
  ## Change CAF based on if K2 or K4 required
  if(CAF_type == "K4"){
    cc_subtype = subtype_data$SCISSORS_CAF_ALL[samples_subset]
    cc_subtype_names = c("restCAF", "Mixed.restCAF","Mixed","Mixed.permCAF", "permCAF", "Absent")
    cc_colscale = colscale4
  }
  if(CAF_type == "K2"){
    cc_subtype = subtype_data$SCISSORS_CAF[samples_subset]
    cc_subtype_names = c("restCAF", "permCAF")
    cc_colscale = colscale3
  }
  
  
  # create annotation vector for the samples
  annotation_col = data.frame(Study = Study, DeCAF = DeCAF_data$DeCAF_graded[samples_subset],
                              'CC Subtype' = cc_subtype,
                              # "Elyada CAF" = subtype_data$Elyada_CAF[samples_subset],
                              #"Moffitt Stroma" = subtype_data$MS_K2[samples_subset], 
                              check.names = F)
  if(include_PurIST == T){
    annotation_col$PurIST = subtype_data$PurIST[samples_subset]
  }
  rownames(annotation_col) = colnames(gene_data)[samples_subset]
  
  # create annotation color vector

  index_PurIST = which(schemaList == "PurIST")
  annotation_col_colors = list(
    DeCAF = colscale,
    'CC Subtype' = cc_colscale,
    'Training Labels' = DeCAFCol,
    'Up Subtype' =DeCAFCol
  )
  names(annotation_col_colors$'CC Subtype')  = cc_subtype_names
  names(annotation_col_colors$'Training Labels') = DeCAFList
  names(annotation_col_colors$DeCAF) = levels(DeCAF_data$DeCAF_graded)
  names(annotation_col_colors$'Up Subtype') = DeCAFList
  
  if(include_PurIST == T){
    annotation_col_colors$PurIST = subtypeColList[[index_PurIST]]
    names(annotation_col_colors$PurIST) = subtypeList[[index_PurIST]]
  }
  
  
  
  ## Sort columns and rows of heatmap
  sort = order(DeCAF_data$DeCAF_prob[samples_subset])
  tsp = match(as.vector(t(classifier$TSPs[classifier$fit$beta[-1]!=0,])),rownames(TSPgeneMat))
  
  ## Classifier gene annotation
  annotation_row = data.frame('Up Subtype' = factor((rep(DeCAFList  ,length(tsp)/2))),check.names = F)
  rownames(annotation_row) = rownames(TSPgeneMat)[tsp]
  
  
  ## Generate heatmap
  heatmap = pheatmap(TSPgeneMat[tsp,sort], annotation_col = annotation_col[sort,],annotation_row = annotation_row, 
                     annotation_colors = annotation_col_colors,show_colnames = F, cluster_cols = F,cluster_rows = F,
                     annotation_names_row = F,width = 12, height = 12,fontsize = 11)
  
  
  
  ## Output
  return(heatmap)
  
  
}

###################################################################################################

## Function to Generate Pooled ROC
Plot_Pooled_ROC <- function(subtype_data, DeCAF_data,  training){
  
  ## Testing Set
  testing0 =  unique(subtype_data$dataSet)[!unique(subtype_data$dataSet) %in% training]
  
  ## Sample Subset
  samples_subset = !(subtype_data$dataSet %in% training) & !is.na(subtype_data$SCISSORS_CAF)
  DeCAF_data = DeCAF_data[samples_subset,]
  subtype_data = subtype_data[samples_subset,]
  
  ## Calculate TP, FP, TN, FN
  data = NULL
  for(i in 1:(length(testing0))){
    pred = DeCAF_data$DeCAF[DeCAF_data$dataSet == testing0[i]]
    pred_p = DeCAF_data$DeCAF_prob[DeCAF_data$dataSet == testing0[i]]
    label = subtype_data$SCISSORS_CAF[subtype_data$dataSet == testing0[i]]
    pred = prediction(prediction = pred_p, labels = as.numeric(label == "permCAF"))
    TP = as.numeric(attr(pred,"tp")[[1]])
    FP = as.numeric(attr(pred,"fp")[[1]])
    TN = as.numeric(attr(pred,"tn")[[1]])
    FN = as.numeric(attr(pred,"fn")[[1]])
    Author = rep(i, length(TP))
    data = rbind(data, cbind(Author, TP, FP, TN, FN))
  }
  data  = data.frame(data)
  
  
  
  
  
    ## Make ROC curve
    return(
      metaROC(data,cex.main = 1,cex.lab =1,cex.axis = 1, model="random-effects", plot.Author=TRUE,plot.inter.var = T,cex.Author = .5,lwd.Author = 2,alpha.trans = .25,Ni = 100)
      
    )  
  
  
  
}
  


#####################################################################################################


## Function to Generate Table of Classifier Performance
Table_Classif_Perform <- function(subtype_data, DeCAF_data,  training){
  
  ## Testing Set
  testing0 =  unique(subtype_data$dataSet)[!unique(subtype_data$dataSet) %in% training]
  
  ## Sample Subset
  samples_subset = !(subtype_data$dataSet %in% training) & !is.na(subtype_data$SCISSORS_CAF)
  DeCAF_data = DeCAF_data[samples_subset,]
  subtype_data = subtype_data[samples_subset,]
  
  ## Empty Vector for metrics
  N = c()
  permCAF = c()
  Accuracy = c()
  Sensitivity = c()
  Specificity = c()
  AUC = c()
  
  
  samples_subset_train = (SubtypeCmb$dataSet %in% training) & !is.na(SubtypeCmb$SCISSORS_CAF)
  DeCAF_data_t = DeCAFCmb[samples_subset_train,]
  subtype_data_t = SubtypeCmb[samples_subset_train,]
  
  
  ## Obtain Metric for each study 
  for(i in 1:(length(testing0))){
    pred = DeCAF_data$DeCAF[DeCAF_data$dataSet == testing0[i]]
    pred_p = DeCAF_data$DeCAF_prob[DeCAF_data$dataSet == testing0[i]]
    label = subtype_data$SCISSORS_CAF[subtype_data$dataSet == testing0[i]]
    N[i] = length(label)
    permCAF[i] = sum(label == "permCAF")
    Accuracy[i] = round(mean(pred == label),3)
    Sensitivity[i] = round(sensitivity(as.factor(pred), as.factor(label), positive = "permCAF"),3)
    Specificity[i] = round(specificity(as.factor(pred), as.factor(label), negative = "restCAF"),3)
    roc = roc(response = label, predictor = pred_p, levels = c("restCAF", "permCAF"), quiet = T)
    AUC[i] = round(roc$auc[1],3)
    
    
  }
  
  ## Combien Results into one df
  pmat = data.frame(N, permCAF, Accuracy, Sensitivity, Specificity, AUC)
  
  ## Clean up Study Name
  testing0 = gsub("_PrimaryPDAC","",testing0)
  rownames(pmat) = testing0[1:(length(testing0))]
  
  
  
  return(pmat)
  

  
}



  
######################################################################################
##############          Step 2: Load In Data                 ########################
######################################################################################

## Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## Load Subtype, Gene Data, and DeCAF Calls
SubtypeCmb <- read.csv("../../data/public_PDAC_clean/combined_data/combined_subtype_DeCAF.csv")
GeneDataCmb <- read.csv("../../data/public_PDAC_clean/combined_data/combined_classifier_gene_exp_DeCAF.csv",
                        row.names = 1)
DeCAFCmb <- read.csv("../../data/public_PDAC_clean/combined_data/combined_classifier_calls_DeCAF.csv")

## Load Classifier to obtain training groups
load("../R/DeCAF/decaf_classifier.Rdata")
classifier = decaf_classifier$classifier2
training = unique(decaf_classifier$anno$study)

## Load Schema List - Aka Color Codes for other classifiers
load("../../data/cmbSubtypes.RData")



######################################################################################
##############          Step 3: Generate Figures              ########################
######################################################################################

#CAF K4 Figure Generator
K4_heatmap = Plot_heatmap_function(subtype_data = SubtypeCmb, 
                                   gene_data = GeneDataCmb, 
                                   DeCAF_data = DeCAFCmb, 
                                   classifier = classifier, 
                                   training = training, 
                                   CAF_type = "K4",
                                   schemaList = schemaList, include_PurIST = T)

filename =  paste("../../results/figures/heatmap_DeCAF_v_SCISSORS.pdf", sep ="")
pdf(filename, width = 10, height = 7)
K4_heatmap
dev.off()


## Figure 4E Generator
ROC_curve = Plot_Pooled_ROC(subtype_data = SubtypeCmb, DeCAF_data = DeCAFCmb,  training = training)
filename =  paste("../../results/figures/DeCAF_meta_ROC.pdf", sep ="")

pdf(filename, height = 7, width = 7, pointsize = 25)
Plot_Pooled_ROC(subtype_data = SubtypeCmb, DeCAF_data = DeCAFCmb,  training = training)
dev.off()
dev.off()


## Figure 4G Generator
table = Table_Classif_Perform(subtype_data = SubtypeCmb, DeCAF_data = DeCAFCmb,  training = training)
write.csv(table, "../../results/tables/DeCAF_Metrics_ByStudy.csv")

kable(table, "latex", booktabs = T) %>%
  kable_styling(latex_options = c("striped","scale_down"),full_width = F)  %>% 
  save_kable(file = "../../results/tables/DeCAF_Metrics_ByStudy.csv")









