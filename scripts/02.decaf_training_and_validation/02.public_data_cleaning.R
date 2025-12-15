#####################################################################################

##                             Data Cleaning Script                      ############

####################################################################################




######################################################################################
############## Step 1: Create Functions for Loading Data ########################
######################################################################################




#### Load Raw Subtype Calls ####
Load_data_subtype <- function(rDataName) {
  
 dataSet =  readRDS(sprintf("../../data/public_PDAC/%s.caf_subtype.rds", rDataName))
  

 
  ## Change myCAF -> proCAF & iCAF -> restCAF
 
 dataSet$Subtype$SCISSORS_CAF_K2_top25.vK[dataSet$Subtype$SCISSORS_CAF_K2_top25.vK == "myCAF"] = "proCAF"
 dataSet$Subtype$SCISSORS_CAF_K2_top25.vK[dataSet$Subtype$SCISSORS_CAF_K2_top25.vK == "Mixed.myCAF"] = "Mixed.proCAF"
 dataSet$Subtype$SCISSORS_CAF_K2_top25.vK[dataSet$Subtype$SCISSORS_CAF_K2_top25.vK == "Mixed.iCAF"] = "Mixed.restCAF"
  dataSet$Subtype$SCISSORS_CAF_K2_top25.vK[dataSet$Subtype$SCISSORS_CAF_K2_top25.vK == "iCAF"] = "restCAF"
  dataSet$Subtype$CAF = dataSet$Subtype$SCISSORS_CAF_K2_top25.vK
  dataSet$Subtype$SCISSORS_CAF_K2_top25.K2[dataSet$Subtype$SCISSORS_CAF_K2_top25.K2 == "myCAF"] = "proCAF"
  dataSet$Subtype$SCISSORS_CAF_K2_top25.K2[dataSet$Subtype$SCISSORS_CAF_K2_top25.K2 == "iCAF"] = "restCAF"
  
  ## Change Mixed proCAF -> proCAF
  dataSet$Subtype$CAF[dataSet$Subtype$CAF == "Mixed"] = NA
  dataSet$Subtype$CAF[dataSet$Subtype$CAF == "Absent"] = NA
  dataSet$Subtype$CAF[dataSet$Subtype$CAF == "Mixed.proCAF"] = "proCAF"
  dataSet$Subtype$CAF[dataSet$Subtype$CAF == "Mixed.restCAF"] = "restCAF"
  dataSet$Subtype$CAF_K4 = dataSet$Subtype$SCISSORS_CAF_K2_top25.vK
  
  ## remove Absent Call 
  dataSet$Subtype$CAF_ALL = dataSet$Subtype$CAF_K4
  
  dataSet$Subtype$CAF_K4[dataSet$Subtype$CAF_K4 == "Absent"] = NA
  
  dataSet$Subtype$CAF.classifier = dataSet$Subtype$SCISSORS_CAF_K2_top25.vK.classifier
  
  return(dataSet)
}


########################################


#### Load Raw Survival Data ####
Load_data_survival <- function(rDataName) {
  
 dataSet =  readRDS(sprintf("../../data/public_PDAC/%s.survival_data.rds", rDataName))
  
  
  return(dataSet)
}




##################################

#### Load Raw Gene and Clinical Data ####
Load_data_gene <- function(rDataName, type = c("array", "seq")) {
  
  
  dat = try(readRDS(sprintf("../../data/public_PDAC/%s.rds", rDataName)))
  
  # process expression information
  dat2 = dat$ex
  dat2 = dat2[!is.na(dat$featInfo$SYMBOL),]
  genes = dat$featInfo$SYMBOL[!is.na(dat$featInfo$SYMBOL)]
  if(any(table(genes) > 1)){
    # need to collapse expression
    if(type == "array"){
      dat2 = aggregate(. ~ genes, dat = dat2, FUN = mean)
      rownames(dat2) = dat2[,1]
      dat2 = dat2[,-1]
    }else{
      dat2 = aggregate(. ~ genes, dat = data.frame(dat2), FUN = sum)
      rownames(dat2) = dat2[,1]
      dat2 = dat2[,-1]
      #dat2 = dat2[rowMeans(dat2 > 5) > .2 ,]
    }
  }else{
    rownames(dat2) = dat$featInfo$SYMBOL
  }
  
  colnames(dat2) = colnames(dat$ex)     
  
  ## Resave back to dat
  dat$expression = dat2
  
  
  return(dat)
}



#################################


## Generate Cleaned Survival data
Generate_Survival_Data <- function(survSet, subtypeSet, DeCAF_calls, save, rDataName) {
  
  survDat <- data.frame(ID = survSet$sampID,
                        time = as.numeric(survSet$time),
                        status = as.numeric(survSet$event),
                        surv_whitelist = survSet$whitelist,
                        SCISSORS_CAF = subtypeSet$Subtype$CAF,
                        SCISSORS_CAF_K4 = subtypeSet$Subtype$CAF_K4,
                        SCISSORS_CAF_ALL = subtypeSet$Subtype$CAF_ALL,
                        SCISSORS_CAF_K2 = subtypeSet$Subtype$SCISSORS_CAF_K2_top25.K2,
                        MS_K2 = subtypeSet$Subtype$MS,
                        Elyada_CAF = subtypeSet$Subtype$Elyada_CAF,
                        PurIST = subtypeSet$Subtype$PurIST,
                        Maurer =  subtypeSet$Subtype$Maurer,
                        Puleo = subtypeSet$Subtype$Puleo,
                        DeCAF = DeCAF_calls$DeCAF,
                        DeCAF_prob = DeCAF_calls$DeCAF_prob,
                        stringsAsFactors = FALSE)
  survDat <- survDat[survSet$whitelist == TRUE,]
  
  if(save == TRUE){
    filename = paste("../../data/public_PDAC_clean/survival/", rDataName, "_survival.csv", sep = "")
    write.csv(survDat,filename, row.names = FALSE)
  }
  return(survDat)
}

#################################


## Generate Cleaned Survival data
Generate_SuppTable_Data <- function(survSet, subtypeSet, DeCAF_calls, rDataName) {
  
  SuppTable_Data <- data.frame(Dataset = rDataName,
                               ID = subtypeSet$Subtype$sampID,
                               PDAC_Pi_whitelist = !is.na(subtypeSet$Subtype$CAF_ALL),
                               DeCAF = DeCAF_calls$DeCAF,
                               DeCAF_graded = DeCAF_calls$DeCAF_graded,
                               DeCAF_prob = DeCAF_calls$DeCAF_prob,
                               PurIST = subtypeSet$Subtype$PurIST,
                               PurIST_graded = subtypeSet$Subtype$PurIST_graded,
                               PurIST_prob = subtypeSet$Subtype$PurIST.prob,
                               SCISSORS_K2 = subtypeSet$Subtype$SCISSORS_CAF_K2_top25.K2,
                               SCISSORS_vK = subtypeSet$Subtype$SCISSORS_CAF_K2_top25.vK,
                               SCISSORS_combined = subtypeSet$Subtype$CAF,
                               SCISSORS_labels_whitelist = !is.na(subtypeSet$Subtype$CAF) & subtypeSet$Subtype$SCISSORS_CAF_K2_top25.vK.classifier,
                               Elyada = subtypeSet$Subtype$Elyada_CAF,
                               Moffitt_stroma = subtypeSet$Subtype$MS,
                               Maurer =  subtypeSet$Subtype$Maurer,
                               Puleo = subtypeSet$Subtype$Puleo)
  if(class(survSet) != "data.frame"){
    SuppTable_Data$Survival_analysis_whitelist = F
    SuppTable_Data$OS_time = NA
    SuppTable_Data$OS_event = NA
  }
  if(class(survSet) == "data.frame"){
    SuppTable_Data$Survival_analysis_whitelist = survSet$whitelist
    SuppTable_Data$OS_time = as.numeric(survSet$time)
    SuppTable_Data$OS_event = as.numeric(survSet$event)
  }


  return(SuppTable_Data)
}

## Generate Cleaned Subtype data
Generate_Subtype_Data <- function(subtypeSet, DeCAF_calls, save, rDataName) {
  
  subtype <- data.frame(ID = subtypeSet$Subtype$sampID,
                        SCISSORS_CAF = subtypeSet$Subtype$CAF,
                        SCISSORS_CAF_K4 = subtypeSet$Subtype$CAF_K4,
                        SCISSORS_CAF_ALL = subtypeSet$Subtype$CAF_ALL,
                        SCISSORS_CAF_K2 = subtypeSet$Subtype$SCISSORS_CAF_K2_top25.K2,
                        MS_K2 = subtypeSet$Subtype$MS,
                        Elyada_CAF = subtypeSet$Subtype$Elyada_CAF,
                        PurIST = subtypeSet$Subtype$PurIST,
                        Maurer =  subtypeSet$Subtype$Maurer,
                        Puleo = subtypeSet$Subtype$Puleo,
                        DeCAF = DeCAF_calls$DeCAF,
                        CAF_whitelist = subtypeSet$Subtype$CAF.classifier,
                        stringsAsFactors = FALSE)
  subtype <- subtype[subtype$CAF_whitelist == TRUE,]
  
  if(save == TRUE){
    filename = paste("../../data/public_PDAC_clean/subtype/", rDataName, "_subtype.csv", sep = "")
    write.csv(subtype,filename, row.names = FALSE)
  }
  return(subtype)
}


################################



## Extract Classifier Genes
Extract_Classifier_Genes <- function(data, classifier, subtypeSet, save, rDataName) {
 
  rows = which(rownames(data) %in% classifier$TSPs[classifier$fit$beta[-1] != 0,])
  white_list = subtypeSet$Subtype$CAF.classifier
  if(sum(is.na(white_list)) == nrow(subtypeSet$Subtype)){
    white_list = rep(T, nrow(subtypeSet$Subtype))
  }
  data = data[rows,white_list]
  tsp = match(as.vector(t(classifier$TSPs[classifier$fit$beta[-1]!=0,])),rownames(data))
  data = data[tsp,]
  
  
  if(save == TRUE){
    filename = paste("../../data/public_PDAC_clean/classifier_gene_exp/", rDataName, "_classifier_gene_expression.csv", sep = "")
    write.csv(data,filename, row.names = FALSE)
  }
  return(data)
}



######################################################################################
################ Step 2: Load Classifier #############################################
######################################################################################

## Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## Load Classifier
load("../R/DeCAF/decaf_classifier.Rdata")
classifier = decaf_classifier$classifier2

## Load Evaluate Classifier Function 
source("../R/DeCAF/decaf_functions.R")
training = unique(decaf_classifier$anno$study)





######################################################################################
######## Step 3: Clean Survival, Subtype, and Clinical Raw Data #######################
######################################################################################




## Names of Data
DataNames <- c("CPTAC", "Dijk", "Grunwald", "Hayashi",
               "Linehan", "Moffitt_GEO_array", "Olive", "PACA_AU_array", "PACA_AU_seq",
               "Puleo_array", "TCGA_PAAD")
## Data Types
data_type = c("seq", "seq", "seq", "seq",
              "seq", "array", "seq", "array", "seq",
              "array", "seq")
## Is there survival info?
survival_info = c(T, T, T, F,
              T, T, F, T, T,
              T, T)



## Loop Through All Data
for(i in 1:length(DataNames)){
  
  print(paste(DataNames[i], data_type[i], "survival info", survival_info[i], sep = " - "))
  
  ## Load Data
  rDataName <- DataNames[i]
  subtypeSet <- Load_data_subtype(rDataName)
  dataSet <- Load_data_gene(rDataName, type = data_type[i] )
  print("loading done")
  
  ## Generate DeCAF Calls
  DeCAF_calls <- apply_decaf(dataSet$expression, classifier)
  
  ## Clean and Save Survival Data if applicable
  if(survival_info[i]){
    survSet <- Load_data_survival(rDataName)
    survDat <- Generate_Survival_Data(survSet, subtypeSet, DeCAF_calls, save = TRUE, rDataName)
    survDat$dataSet <- rDataName
    if(i == 1){
      survDatCmb <- survDat
    } else {
      survDatCmb <- rbind(survDatCmb,survDat)
    }
    
  }
  if(!survival_info[i]){
    survSet = NA
  }
  
  ## Save for Supplemental Table
  SuppTable_Data_study = Generate_SuppTable_Data(survSet, subtypeSet, DeCAF_calls, rDataName)
  
  ## Clean and Save Data
  subtypeDat <- Generate_Subtype_Data(subtypeSet, DeCAF_calls, save = TRUE, rDataName)
  classif_exp <- Extract_Classifier_Genes(dataSet$expression, classifier, subtypeSet, save = T, rDataName)
  
  
  ## Add Study Name
  subtypeDat$dataSet <- rDataName
  DeCAF_calls$dataSet <- rDataName
  
  ## Only Subtype True DeCAF Calls
  DeCAF_calls = DeCAF_calls[subtypeSet$Subtype$CAF.classifier,]
  
  ## Combine With Other Studies
  if(i == 1){
    subtypeDatCmb <- subtypeDat
    classif_expCmb <- classif_exp
    DeCAF_callsCmb <- DeCAF_calls
    SuppTable_Data <- SuppTable_Data_study
  } else {
    subtypeDatCmb <- rbind(subtypeDatCmb, subtypeDat)
    classif_expCmb <- cbind(classif_expCmb, classif_exp)
    DeCAF_callsCmb <- rbind(DeCAF_callsCmb, DeCAF_calls)
    SuppTable_Data <- rbind(SuppTable_Data, SuppTable_Data_study)
  }
  
  
  
  
}




## Save Combined Data
write.csv(survDatCmb,"../../data/public_PDAC_clean/combined_data/combined_survival_DeCAF.csv", row.names = FALSE)
write.csv(subtypeDatCmb,"../../data/public_PDAC_clean/combined_data/combined_subtype_DeCAF.csv", row.names = FALSE)
write.csv(classif_expCmb,"../../data/public_PDAC_clean/combined_data/combined_classifier_gene_exp_DeCAF.csv")
write.csv(DeCAF_callsCmb,"../../data/public_PDAC_clean/combined_data/combined_classifier_calls_DeCAF.csv")

## Change Survival of Duplicated ID's to Individual Level Only
idDup <- SuppTable_Data$ID[(duplicated(SuppTable_Data$ID))]
idDup_Paca_array <- which((SuppTable_Data$ID %in% idDup) & (SuppTable_Data$dataSet == "PACA_AU_array"))
SuppTable_Data$Survival_analysis_whitelist[idDup_Paca_array] <- "Individual Study Only - Not Combined"
write.xlsx(SuppTable_Data, "../../results/tables/Supplementary Table - clinical and molecular data - public.xlsx", row.names = F)



