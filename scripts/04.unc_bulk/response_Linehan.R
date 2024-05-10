# set dir
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

############################################################
#################### Functions and libraries

# load libraries
library(survival)
library(survminer)
library(gtable)
library(gridExtra)
library(grid)
library(gplots)
#library(coin)
library(openxlsx)
library(MASS)
#library(splitstackshape)
library(reshape2)
#library(mclust)
#library(riverplot)
#library(Rmisc)
library(stringr)

# load subtype info
### This is a combined subtype object with
### subtypeColList, subtypeGeneList, subtypeList and schemaList
load("../../data/cmbSubtypes.RData")
print("Subtype schemas available for use:")
print(schemaList)

# load functions
file.sources <- list.files("../R/",pattern="*.R")
file.sources <- paste("../R/", file.sources, sep="")
sapply(file.sources, source)

############################## Analyze survival ################################

dataSet <- readRDS("../../data/public_PDAC/Linehan.rds")
dataSet <- Call_DeCAF(dataSet)
subtypeCall <- readRDS("../../data/public_PDAC/Linehan.caf_subtype.rds")

pdf("../../results/figures/Linehan_response.pdf")

### response
responseDat <- data.frame(sampID = dataSet$sampInfo$patient_id,
                          Treatment = dataSet$sampInfo$Treatment,
                          pre.post = dataSet$sampInfo$treatment.1before.2after,
                          Change = dataSet$sampInfo$Percent.change,
                          Resection = dataSet$sampInfo$Resection,
                          DeCAF = dataSet$sampInfo$DeCAF,
                          DeCAF_prob = dataSet$sampInfo$DeCAF_prob,
                          PurIST = subtypeCall$Subtype$PurIST,
                          PurIST_prob = subtypeCall$Subtype$PurIST.prob,
                          stringsAsFactors = FALSE)

responseDat <- responseDat[which(!(is.na(responseDat$Change))),]
responseDat.pre <-  responseDat[which(responseDat$pre.post == 1),]
names(responseDat.pre)[which(names(responseDat.pre) %in% c("DeCAF","DeCAF_prob","PurIST","PurIST_prob"))] <- paste(names(responseDat.pre)[which(names(responseDat.pre) %in% c("DeCAF","DeCAF_prob","PurIST","PurIST_prob"))],"pre",sep=".")
responseDat.post <-  responseDat[which(responseDat$pre.post == 2),]
responseDat.pre$DeCAF.post <- responseDat.post$DeCAF[match(responseDat.pre$sampID,responseDat.post$sampID)]
responseDat.pre$DeCAF_prob.post <- responseDat.post$DeCAF_prob[match(responseDat.pre$sampID,responseDat.post$sampID)]
responseDat.pre$PurIST.post <- responseDat.post$PurIST[match(responseDat.pre$sampID,responseDat.post$sampID)]
responseDat.pre$PurIST_prob.post <- responseDat.post$PurIST_prob[match(responseDat.pre$sampID,responseDat.post$sampID)]

responseDat.pre$ChangeStatus[responseDat.pre$Change <=0] <- "Shrinkage"
responseDat.pre$ChangeStatus[responseDat.pre$Change >0] <- "Growth"

responseDat.pre$RECIST[responseDat.pre$Change>=20] <- "PD"  
responseDat.pre$RECIST[responseDat.pre$Change<=-30] <- "PR"  
responseDat.pre$RECIST[responseDat.pre$Change>-30 & responseDat.pre$Change<20] <- "SD"

responseDat.all <- responseDat.pre[responseDat.pre$Treatment=="FOLFIRINOX+PF04136309",]

# cibersort -----------------------------------------
cibersort <- read.table("../../data/derived/CIBERSORT_Linehan.txt", sep = "\t", header = T)
tmpSplit <- data.frame(str_split_fixed(gsub("Linehan_S1124.","",cibersort$Mixture), "\\.", 2))
tmpSplit$X2 <- gsub("0","",tmpSplit$X2)

# combine
tmpPre <- cibersort[match(paste0(responseDat.all$sampID,".1"),
                          paste0(tmpSplit$X1,".",tmpSplit$X2)), -1]
names(tmpPre) <- paste0(names(tmpPre),"_pre")

tmpPost <- cibersort[match(paste0(responseDat.all$sampID,".2"),
                           paste0(tmpSplit$X1,".",tmpSplit$X2)), -1]
names(tmpPost) <- paste0(names(tmpPost),"_post")

tmpImmuneChange <- tmpPre
for (i in 1:ncol(tmpPre)) {
  tmpImmuneChange[,i] <- tmpPost[,i]-tmpPre[,i]
}
names(tmpImmuneChange) <- gsub("_pre","_diff",names(tmpImmuneChange))

responseDat.all <- cbind(responseDat.all, tmpPre, tmpPost, tmpImmuneChange)
responseDat.all$M1_M2_ratio_pre <- (responseDat.all$Macrophages.M1_pre+0.001)/(responseDat.all$Macrophages.M2_pre+0.001)
responseDat.all$M1_M2_ratio_post <- (responseDat.all$Macrophages.M1_post+0.001)/(responseDat.all$Macrophages.M2_post+0.001)

# plot -----------------------------------------
#responseDat.pre <- responseDat.pre[order(rank(responseDat.pre$Treatment),responseDat.pre$Change,decreasing = TRUE),]
responseDat.pre <- responseDat.pre[order(responseDat.pre$Change,decreasing = TRUE),]

# plot response vs DeCAF_prob
par(mgp=c(2,0.6,0))
par(mfrow=c(3,3),mar=c(3,3,1,1)+0.1)

# all
#responseDat.all <- responseDat.pre[responseDat.pre$Treatment=="FOLFIRINOX+PF04136309",]
corr <- cor.test(as.numeric(responseDat.all$Change),
                 as.numeric(responseDat.all$DeCAF_prob.pre),
                 method="spearman", use= "complete.obs")

# changed to pearson in actual figures

# permCAF response
corr <- cor.test(as.numeric(responseDat.all$Change[responseDat.all$DeCAF.pre=="permCAF"]),
                 as.numeric(responseDat.all$DeCAF_prob.pre[responseDat.all$DeCAF.pre=="permCAF"]),
                 method="pearson", use= "complete.obs")
plot(as.numeric(responseDat.all$Change[responseDat.all$DeCAF.pre=="permCAF"]),
     as.numeric(responseDat.all$DeCAF_prob.pre[responseDat.all$DeCAF.pre=="permCAF"]),
     main = "",
     xlab = "% tumor size change", ylab = "permCAF probability (pre)",
     cex.lab = 1.2, cex.axis = 1.2,
     xlim = c(-50,30), ylim = c(0.5,1),
     col=c("violetred1","turquoise4")[as.factor(responseDat.all$DeCAF.pre[responseDat.all$DeCAF.pre=="permCAF"])], 
     #bg = c("violetred1","turquoise4")[as.factor(responseDat.all$DeCAF.post[responseDat.all$DeCAF.pre=="permCAF"])], 
     pch=c(22,24,21)[as.factor(as.numeric(as.factor(responseDat.all$RECIST[responseDat.all$DeCAF.pre=="permCAF"])))], 
     cex = 1.5, lwd = 1)
abline(v = 20, lwd=1, lty=2)
abline(v = -30, lwd=1, lty=2)
#abline(h = 0.5, lwd=1, lty=2)
text(-5, 0.6, cex = 1.2,
     pos = 3, paste0("r = ",round(corr$estimate,3)))
text(-5, 0.55,  cex = 1.2,
     pos = 3, paste0("p = ",round(corr$p.value,3)))




# response vs immune change -------------------

### neutrophil
corr <- cor.test(responseDat.all$Neutrophils_diff, responseDat.all$Change, 
                 method = "spearman")

# used pearson
corr <- cor.test(responseDat.all$Neutrophils_diff, responseDat.all$Change, 
                 method = "pearson")

plot(responseDat.all$Change, responseDat.all$Neutrophils_diff, 
     main = "",
     xlab = "% tumor size change", ylab = "Neutrophils % change",
     cex.lab = 1.2, cex.axis = 1.2,
     xlim = c(-50,30), #ylim = c(0.5,1),
     col=c("violetred1","turquoise4")[as.factor(responseDat.all$DeCAF.pre)], 
     #bg = c("violetred1","turquoise4")[as.factor(responseDat.all$DeCAF.post)], 
     pch=c(22,24,21)[as.factor(as.numeric(as.factor(responseDat.all$RECIST)))], 
     cex = 1.5, lwd = 1)
text(10, -0.25, cex = 1.2,
     pos = 3, paste("r = ",round(corr$estimate,3),sep=""))
text(10, -0.3, cex = 1.2,
     pos = 3, paste("p = ",round(corr$p.value,3),sep=""))






wilTest.p <- wilcox.test(responseDat.all$Neutrophils_pre[responseDat.all$DeCAF.pre=="permCAF"], 
            responseDat.all$Neutrophils_post[responseDat.all$DeCAF.pre=="permCAF"], paired = T)
wilTest.r <- wilcox.test(responseDat.all$Neutrophils_pre[responseDat.all$DeCAF.pre=="restCAF"], 
            responseDat.all$Neutrophils_post[responseDat.all$DeCAF.pre=="restCAF"], paired = T)
boxplot(responseDat.all$Neutrophils_pre[responseDat.all$DeCAF.pre=="permCAF"], 
        responseDat.all$Neutrophils_post[responseDat.all$DeCAF.pre=="permCAF"],
        responseDat.all$Neutrophils_pre[responseDat.all$DeCAF.pre=="restCAF"],
        responseDat.all$Neutrophils_post[responseDat.all$DeCAF.pre=="restCAF"],
        at = c(1,2,4,5), names = c("Pre","Post","Pre","Post"), 
        ylab = "Neutrophil %", main = "",
        border = c("violetred1","violetred1","turquoise4","turquoise4"),
        col = "white", ylim =c(0,0.6),
        cex.lab = 1.2, cex.axis = 1.2)
text(1.5, 0.49, cex = 1.2,
     pos = 3, paste0("permCAF\np = ", round(wilTest.p$p.value,3)))
text(4.5, 0.49, cex = 1.2,
     pos = 3, paste0("restCAF\np = ", round(wilTest.r$p.value,3)))




# restCAF response
corr <- cor.test(as.numeric(responseDat.all$Change[responseDat.all$DeCAF.pre=="restCAF"]),
                 as.numeric(responseDat.all$DeCAF_prob.pre[responseDat.all$DeCAF.pre=="restCAF"]),
                 method="pearson", use= "complete.obs")

plot(as.numeric(responseDat.all$Change[responseDat.all$DeCAF.pre=="restCAF"]),
     as.numeric(responseDat.all$DeCAF_prob.pre[responseDat.all$DeCAF.pre=="restCAF"]),
     main = "",
     xlab = "% tumor size change", ylab = "permCAF probability (pre)",
     cex.lab = 1.2, cex.axis = 1.2,
     xlim = c(-50,30), ylim = c(0,0.5),
     col=c("turquoise4")[as.factor(responseDat.all$DeCAF.pre[responseDat.all$DeCAF.pre=="restCAF"])], 
     #bg = c("violetred1","turquoise4")[as.factor(responseDat.all$DeCAF.post[responseDat.all$DeCAF.pre=="restCAF"])], 
     pch=c(22,24,21)[as.factor(as.numeric(as.factor(responseDat.all$RECIST[responseDat.all$DeCAF.pre=="restCAF"])))], 
     cex = 1.5, lwd = 1)
abline(v = 20, lwd=1, lty=2)
abline(v = -30, lwd=1, lty=2)
#abline(h = 0.5, lwd=1, lty=2)
text(-5, 0.4, cex = 1.2,
     pos = 3, paste0("r = ",round(corr$estimate,3)))
text(-5, 0.35,  cex = 1.2,
     pos = 3, paste0("p = ",round(corr$p.value,3)))


# 
### macrophage M2
corr <- cor.test(responseDat.all$Macrophages.M2_diff, responseDat.all$Change, 
                 method = "spearman")
corr <- cor.test(responseDat.all$Macrophages.M2_diff, responseDat.all$Change, 
                 method = "pearson")

plot(responseDat.all$Change, responseDat.all$Macrophages.M2_diff, 
     main = "",
     xlab = "% tumor size change", ylab = "M2 macrophage % change",
     cex.lab = 1.2, cex.axis = 1.2,
     xlim = c(-50,30), #ylim = c(0.5,1),
     col=c("violetred1","turquoise4")[as.factor(responseDat.all$DeCAF.pre)], 
     #bg = c("violetred1","turquoise4")[as.factor(responseDat.all$DeCAF.post)], 
     pch=c(22,24,21)[as.factor(as.numeric(as.factor(responseDat.all$RECIST)))], 
     cex = 1.5, lwd = 1)
text(10, 0.2, cex = 1.2,
     pos = 3, paste("r = ",round(corr$estimate,3),sep=""))
text(10, 0.15, cex = 1.2,
     pos = 3, paste("p = ",round(corr$p.value,3),sep=""))


wilTest.p <- wilcox.test(responseDat.all$Macrophages.M2_pre[responseDat.all$DeCAF.pre=="permCAF"], 
            responseDat.all$Macrophages.M2_post[responseDat.all$DeCAF.pre=="permCAF"], paired = T)
wilTest.r <-wilcox.test(responseDat.all$Macrophages.M2_pre[responseDat.all$DeCAF.pre=="restCAF"], 
            responseDat.all$Macrophages.M2_post[responseDat.all$DeCAF.pre=="restCAF"], paired = T)
boxplot(responseDat.all$Macrophages.M2_pre[responseDat.all$DeCAF.pre=="permCAF"], 
        responseDat.all$Macrophages.M2_post[responseDat.all$DeCAF.pre=="permCAF"],
        responseDat.all$Macrophages.M2_pre[responseDat.all$DeCAF.pre=="restCAF"],
        responseDat.all$Macrophages.M2_post[responseDat.all$DeCAF.pre=="restCAF"],
        at = c(1,2,4,5), names = c("Pre","Post","Pre","Post"), 
        ylab = "M2 macrophage %", main = "",
        border = c("violetred1","violetred1","turquoise4","turquoise4"),
        col = "white", ylim =c(0,0.6),
        cex.lab = 1.2, cex.axis = 1.2)
text(1.5, 0.49, cex = 1.2,
     pos = 3, paste0("permCAF\np = ", round(wilTest.p$p.value,3)))
text(4.5, 0.49, cex = 1.2,
     pos = 3, paste0("restCAF\np = ", round(wilTest.r$p.value,3)))


# linehan pre vs post --------------------------------------
corr <- cor.test(as.numeric(responseDat.all$DeCAF_prob.post),
                 as.numeric(responseDat.all$DeCAF_prob.pre),
                 method="spearman", use= "complete.obs")
corr <- cor.test(as.numeric(responseDat.all$DeCAF_prob.post),
                 as.numeric(responseDat.all$DeCAF_prob.pre),
                 method="pearson", use= "complete.obs")

plot(as.numeric(responseDat.all$DeCAF_prob.post),
     as.numeric(responseDat.all$DeCAF_prob.pre),
     main = "",
     ylab = "permCAF probability (pre)", xlab = "permCAF probability (post)",
     cex.lab = 1.2, cex.axis = 1.2,
     xlim = c(0,1), ylim = c(0,1),
     #col=c("darkred","darkgreen","gold2")[as.numeric(as.factor(responseDat.all$RECIST))], 
     col=c("violetred1","turquoise4")[as.numeric(as.factor(responseDat.all$DeCAF.pre))],
     pch=c(0,2,1)[as.numeric(as.factor(responseDat.all$RECIST))], 
     cex=1.5)
abline(v = 0.5, lwd=1, lty=2)
abline(h = 0.5, lwd=1, lty=2)
text(0.2, 0.85, cex = 1.2,
     pos = 3, paste("r = ",round(corr$estimate,3),sep=""))
text(0.2, 0.75, cex = 1.2,
     pos = 3, 
     "p < 0.001") #paste0("p = ", corr$p.value))

# print legend
#plot.new()
#legend(xy.coords(x=0,y=.98),
#       legend=c("DeCAF by treatment(tx)",
#                "permCAF(pre-tx) & permCAF(post-tx) ",
#                "permCAF(pre-tx) & restCAF(post-tx) ",
#                "restCAF(pre-tx) & permCAF(post-tx) ",
#                "restCAF(pre-tx) & restCAF(post-tx) ",
#                "Best respone","PR","SD","PD"),
#       fill = c("white","violetred1","turquoise4","violetred1","turquoise4", 
#              "white","black","black","black"),
#       col = c("white","violetred1","violetred1","turquoise4","turquoise4",
#               "white","black","black","black"),
#       pch = c(21,21,21, 21,21,
#               21,22,24,21),
#       border=FALSE, bty="n",
#       x.intersp = 1 , cex = 1, horiz=FALSE)

dev.off()
