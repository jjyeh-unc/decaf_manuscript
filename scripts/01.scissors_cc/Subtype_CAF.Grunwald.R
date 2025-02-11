############################## Functions and libraries ##############################
# set dir
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# load libraries
library(ConsensusClusterPlus) # R3.xx and R4.xx have different versions of ConsensusClusterPlus
library("RColorBrewer")
library(openxlsx)
library(stringr)

# load functions
file.sources <- list.files("../R/",pattern="*.R")
file.sources <- paste("../R/", file.sources, sep="")
sapply(file.sources, source)

# load PurIST
load("../R/PurIST/fitteds_public_2019-02-12.Rdata")
source("../R/PurIST/functions.R")

# load subtype info
### This is a combined subtype object with
### subtypeColList, subtypeGeneList, subtypeList and schemaList
load("../../data/cmbSubtypes.RData")
print("Subtype schemas available for use:")
print(schemaList)

############################## Subtyping ##############################
rDataName <- "Grunwald"
dataSet <- readRDS(paste("../../data/public_PDAC/",rDataName,".rds",sep=""))
sampSub <- 30:nrow(dataSet$sampInfo)
cafSubtype <- list()
cafSubtype$Subtype <- data.frame(sampID = colnames(dataSet$ex), stringsAsFactors = FALSE)
                 
## PurIST
dataSet <- Call_PurIST(dataSet)
cafSubtype$Subtype["PurIST"] <- dataSet$sampInfo$PurIST
cafSubtype$Subtype["PurIST_graded"] <- dataSet$sampInfo$PurIST_graded
cafSubtype$Subtype["PurIST.prob"] <- dataSet$sampInfo$PurIST.prob

if(FALSE) {
## SCISSORS unscaled top10
tmpSchema <- "SCISSORS_CAF_K2_top10"
tmpRowScale <- FALSE
tmpK <- 2
dataSet <- Call_consensus(dataSet, tmpSchema, sampSub, 2, tmpK, tmpRowScale, TRUE, FALSE, "km","euclidean", "R-4.X.X")
tmpSubtypeLab <- c("myCAF","iCAF")
dataSet <- Relabel_subtype(dataSet, tmpSubtypeLab, tmpSchema, tmpK, tmpRowScale)
cafSubtype <- Save_cafSubtype_v2(cafSubtype, dataSet, tmpSchema, tmpK, tmpRowScale, "K2")
}

## SCISSORS unscaled top25
tmpSchema <- "SCISSORS_CAF_K2_top25"
tmpRowScale <- FALSE
tmpK <- 3
dataSet <- Call_consensus(dataSet, tmpSchema, sampSub, 2, tmpK, tmpRowScale, TRUE, TRUE, "km","euclidean", "R-4.X.X")
tmpSubtypeLab <- c("myCAF","Mixed","iCAF")
dataSet <- Relabel_subtype(dataSet, tmpSubtypeLab, tmpSchema, tmpK, tmpRowScale)
cafSubtype <- Save_cafSubtype(cafSubtype, dataSet, tmpSchema, tmpK, tmpRowScale, "vK")

## SCISSORS unscaled top25
tmpSchema <- "SCISSORS_CAF_K2_top25"
tmpRowScale <- FALSE
tmpK <- 2
dataSet <- Call_consensus(dataSet, tmpSchema, sampSub, 2, tmpK, tmpRowScale, TRUE, TRUE, "km","euclidean", "R-4.X.X")
tmpSubtypeLab <- c("myCAF","iCAF")
dataSet <- Relabel_subtype(dataSet, tmpSubtypeLab, tmpSchema, tmpK, tmpRowScale)
cafSubtype <- Save_cafSubtype(cafSubtype, dataSet, tmpSchema, tmpK, tmpRowScale, "K2")

## Elyada unscaled
tmpSchema <- "Elyada_CAF"
tmpRowScale <- FALSE
tmpK <- 2
dataSet <- Call_consensus(dataSet, tmpSchema, sampSub, 2, tmpK, tmpRowScale, TRUE, FALSE, "km","euclidean", "R-4.X.X")
tmpSubtypeLab <- c("myCAF","iCAF")
dataSet <- Relabel_subtype(dataSet, tmpSubtypeLab, tmpSchema, tmpK, tmpRowScale)
cafSubtype <- Save_cafSubtype(cafSubtype, dataSet, tmpSchema, tmpK, tmpRowScale)

## MS unscaled
tmpSchema <- "MS"
tmpRowScale <- FALSE
tmpK <- 2
dataSet <- Call_consensus(dataSet, tmpSchema, sampSub, 2, tmpK, tmpRowScale, TRUE, TRUE, "km","euclidean", "R-4.X.X")
tmpSubtypeLab <- c("Activated","Normal")
dataSet <- Relabel_subtype(dataSet, tmpSubtypeLab, tmpSchema, tmpK, tmpRowScale)
cafSubtype <- Save_cafSubtype(cafSubtype, dataSet, tmpSchema, tmpK, tmpRowScale)

## Maurer
tmpSchema <- "Maurer"
tmpRowScale <- FALSE
tmpK <- 2
dataSet <- Call_consensus(dataSet, tmpSchema, sampSub, 2, tmpK, tmpRowScale, TRUE, FALSE, "km","euclidean", "R-4.X.X")
tmpSubtypeLab <- c("ECM-rich","Immune-rich")
dataSet <- Relabel_subtype(dataSet, tmpSubtypeLab, tmpSchema, tmpK, tmpRowScale)
cafSubtype <- Save_cafSubtype(cafSubtype, dataSet, tmpSchema, tmpK, tmpRowScale)

## Puleo 
tmpSchema <- "Puleo"
dataSet <- Call_centroid(dataSet = dataSet, 
                         tmpSchema = tmpSchema, # this gene set fit the format
                         sampSub = sampSub,
                         lowVarFlt = FALSE, # filter low vairance genes or not; set this to FALSE first. If you get error for low variance of genes, set to TRUE.
                         distance = "euclidean", # see ?ConsensusClusterPlus() for other options
                         Rversion = "R-4.X.X" # or "R-3.X.X" # different R versions use different clusterAlg (km vs kmdist)
)
cafSubtype$Puleo <- dataSet$Puleo
cafSubtype$Subtype$Puleo <- dataSet$sampInfo$Puleo
cafSubtype$Subtype$Puleo.classifier <- FALSE
cafSubtype$Subtype$Puleo.classifier[which(!is.na(cafSubtype$Subtype$Puleo))] <- TRUE

# save
saveRDS(cafSubtype, file = paste("../../data/public_PDAC/",rDataName,".caf_subtype.rds",sep = ""))

