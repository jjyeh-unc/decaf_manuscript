
library(stringr)
library(openxlsx)
library(survival)
library(survminer)
library(survivalAnalysis)
#library(forestplot)
library(magrittr)
library(dplyr)

dataSet <- readRDS("../../data/newly_generated/UNC_bulk.rds")

## parse clincial ---------------------------------------------------------------------
survDat <- data.frame(sampID = colnames(dataSet$ex), 
                      time = dataSet$clinicalDat$`Overall.Survival.(OR.until.death)`,
                      event = dataSet$clinicalDat$Last.Event,
                      neoadj.tx = dataSet$clinicalDat$Neoadj.Tx,
                      neoadj.tx.clean = dataSet$clinicalDat$neo_reg_clean,
                      M = dataSet$clinicalDat$M,
                      DeCAF.prob = dataSet$subtypeCall$DeCAF.prob,
                      DeCAF = dataSet$subtypeCall$DeCAF,
                      DeCAF.graded = dataSet$subtypeCall$DeCAF.graded,
                      PurIST.prob = dataSet$subtypeCall$PurIST.prob,
                      PurIST = dataSet$subtypeCall$PurIST,
                      PurIST.graded = dataSet$subtypeCall$PurIST.graded,
                      `Myxoid(yes.1,.no.0)` = dataSet$pathCall$`Myxoid(yes.1,.no.0)`,
                      `Collagenized(yes.1,.no.0)` = dataSet$pathCall$`Collagenized(yes.1,.no.0)`,
                      `Myx1/mixed2/Coll3` = dataSet$pathCall$`Myx1/mixed2/Coll3`,
                      #whitelist = FALSE,
                      stringsAsFactors = FALSE)
survDat$event[survDat$event == 500 | 
                 survDat$event ==900 |
                 survDat$event ==600 |
                survDat$event ==650 ] <- 0
survDat$event[survDat$event == 2 | 
                 survDat$event ==3] <- 1
survDat$event <- as.numeric(survDat$event)
survDat$time <- as.numeric(survDat$time)
survDat$DeCAF <- factor(survDat$DeCAF , levels = c("permCAF","restCAF"))
survDat$PurIST <- factor(survDat$PurIST, levels = c("Basal-like","Classical"))
survDat0 <- survDat
write.xlsx(survDat0, "../../results/tables/UNC-bulk_clinical_and_molecular_data.xlsx")


pdf("../../results/figures/yehseq_surv.pdf")

# DeCAF --------------------------------------------------------------------
survDat <- survDat0
km <- with(survDat, Surv(time,event))
p = coxph(km ~ DeCAF, data = survDat)
bic = round(BIC(p),3)
hr <- round(1/summary(p)$coefficients[2],3)
km_fit <- survfit(km ~ DeCAF, data = survDat, type = "kaplan-meier")
p <- ggsurvplot(font.legend = 16, 
                font.title = 16,        
                font.x = 16,  
                font.y = 16, 
                font.tickslab = 16,
                risk.table.fontsize = 7,
                risk.table.height = 0.35,
                surv.median.line = "hv", 
                size = 0.3, 
                censor.size=7, 
                pval.size = 8,
                pval.coord = c(96, 0.75),
                xlim = c(0,144),
                km_fit, conf.int = F, pval = T,
                          legend.title="",break.time.by = 12,
                          legend.labs=c("permCAF","restCAF"),
                          palette = c("violetred1","turquoise4"),
                          xlab = "Time (months)", risk.table = T,
                          title = paste("YehSeq\nperm vs rest HR=",hr,
                                        " BIC=",bic,sep=""))
print(p)

# no neo
survDat <- survDat0[which(survDat0$neoadj.tx %in% c("No","no")), ]
km <- with(survDat, Surv(time,event))
p = coxph(km ~ DeCAF, data = survDat)
bic = round(BIC(p),3)
hr <- round(1/summary(p)$coefficients[2],3)
km_fit <- survfit(km ~ DeCAF, data = survDat, type = "kaplan-meier")
p <- ggsurvplot(font.legend = 16, 
                font.title = 16,        
                font.x = 16,  
                font.y = 16, 
                font.tickslab = 16,
                risk.table.fontsize = 7,
                risk.table.height = 0.35,
                surv.median.line = "hv", 
                size = 0.3, 
                censor.size=7, 
                pval.size = 8,
                pval.coord = c(96, 0.75),
                xlim = c(0,144),
                km_fit, conf.int = F, pval = T,
                          legend.title="",break.time.by = 12,
                          legend.labs=c("permCAF","restCAF"),
                          palette = c("violetred1","turquoise4"),
                          xlab = "Time (months)", risk.table = T,
                          title = paste("YehSeq w/o neo\nperm vs rest HR=",hr,
                                        " BIC=",bic,sep=""))
print(p)

# FFX
survDat <- survDat0[which(survDat0$neoadj.tx.clean %in% c("FOLFIRINOX")), ]
km <- with(survDat, Surv(time,event))
p = coxph(km ~ DeCAF, data = survDat)
bic = round(BIC(p),3)
hr <- round(1/summary(p)$coefficients[2],3)
km_fit <- survfit(km ~ DeCAF, data = survDat, type = "kaplan-meier")
p <- ggsurvplot(font.legend = 16, 
                font.title = 16,        
                font.x = 16,  
                font.y = 16, 
                font.tickslab = 16,
                risk.table.fontsize = 7,
                risk.table.height = 0.35,
                surv.median.line = "hv", 
                size = 0.3, 
                censor.size=7, 
                pval.size = 8,
                pval.coord = c(96, 0.75),
                xlim = c(0,144),
                km_fit, conf.int = F, pval = T,
                          legend.title="",break.time.by = 12,
                          legend.labs=c("permCAF","restCAF"),
                          palette = c("violetred1","turquoise4"),
                          xlab = "Time (months)", risk.table = T,
                          title = paste("YehSeq FOLFIRINOX\nperm vs rest HR=",hr,
                                        " BIC=",bic,sep=""))
print(p)

# interaction ------------------------------------------------------------------------------
# PurIST * DeCAF
survDat <- survDat0
survDat$cmb <- factor(paste(survDat$DeCAF, survDat$PurIST,  sep = "*"),
                    levels = c("permCAF*Basal-like","permCAF*Classical","restCAF*Basal-like","restCAF*Classical"))
km <- with(survDat, Surv(time,event))
p = coxph(km ~ cmb, data = survDat)
bic = round(BIC(p),3)
#km_fit <- survfit(km ~ cmb, data = survDat, type = "kaplan-meier")
km_fit <- survfit(km ~ DeCAF+PurIST, data = survDat, type = "kaplan-meier")
p <- ggsurvplot(font.legend = 16, 
                font.title = 16,        
                font.x = 16,  
                font.y = 16, 
                font.tickslab = 16,
                risk.table.fontsize = 7,
                risk.table.height = 0.35,
                size = 0.3, 
                censor.size=7, 
                pval.size = 8,
                pval.coord = c(96, 0.75),
                xlim = c(0,144),
                km_fit, conf.int = F, pval = T,
                legend.title="",
                break.time.by = 12,
                legend.labs=levels(survDat$cmb),
                palette = c("violetred1","#FE46A5","turquoise4","#008080"),
                linetype = c(1,2,1,2),
                xlab = "Time (months)", 
                risk.table = T,
                title = paste("YehSeq DeCAF*PurIST",sep=""))
print(p)

# interaction
table(survDat[,c("PurIST","DeCAF")])
fisher.test(table(survDat[,c("PurIST","DeCAF")]))

# multi-variate
survDat$PurIST <- factor(survDat$PurIST, levels = c("Classical","Basal-like"))
survDat$DeCAF <- factor(survDat$DeCAF, levels = c("restCAF","permCAF"))
fit <- coxph(formula = Surv(time, event) ~ DeCAF+PurIST , data = survDat)
p <- ggforest(fit, 
              main = "Hazard ratio",
              fontsize = 1.2,
              cpositions = c(0.02, 0.13, 0.32)) 
print(p)

if(FALSE) {
p <- survDat %>%
   analyse_multivariate(vars(time, event),
                        vars(PurIST, DeCAF)) %>%
   forest_plot(label_headers = c(endpoint = "Endpoint", factor = "Variables", n = "n"))
}

# PurIST * DeCAF (no neo) -------------------------------------------
survDat <- survDat0[which(survDat0$neoadj.tx %in% c("No","no")), ]
survDat$cmb <- factor(paste(survDat$DeCAF, survDat$PurIST,  sep = "*"),
                    levels = c("permCAF*Basal-like","permCAF*Classical","restCAF*Basal-like","restCAF*Classical"))
km <- with(survDat, Surv(time,event))
p = coxph(km ~ cmb, data = survDat)
bic = round(BIC(p),3)
#km_fit <- survfit(km ~ cmb, data = survDat, type = "kaplan-meier")
km_fit <- survfit(km ~ DeCAF+PurIST, data = survDat, type = "kaplan-meier")
p <- ggsurvplot(font.legend = 16, 
                font.title = 16,        
                font.x = 16,  
                font.y = 16, 
                font.tickslab = 16,
                risk.table.fontsize = 7,
                risk.table.height = 0.35,
                size = 0.3, 
                censor.size=7, 
                pval.size = 8,
                pval.coord = c(96, 0.75),
                xlim = c(0,144),
                km_fit, conf.int = F, pval = T,
                legend.title="",
                break.time.by = 12,
                legend.labs=levels(survDat$cmb),
                palette = c("violetred1","#FE46A5","turquoise4","#008080"),
                linetype = c(1,2,1,2),
                xlab = "Time (months)", 
                risk.table = T,
                title = paste("YehSeq w/o neo DeCAF*PurIST",sep=""))
print(p)

# interaction
table(survDat[,c("PurIST","DeCAF")])
fisher.test(table(survDat[,c("PurIST","DeCAF")]))

# multi-variate
survDat$PurIST <- factor(survDat$PurIST, levels = c("Classical","Basal-like"))
survDat$DeCAF <- factor(survDat$DeCAF, levels = c("restCAF","permCAF"))
fit <- coxph(formula = Surv(time, event) ~ DeCAF+PurIST , data = survDat)
p <- ggforest(fit, 
              main = "YehSeq w/o neo",
              fontsize = 1.2,
              cpositions = c(0.02, 0.13, 0.32)) 
print(p)

# PurIST * DeCAF (FFX) ------------------------------------------------------
survDat <- survDat0[which(survDat0$neoadj.tx.clean %in% c("FOLFIRINOX")), ]
survDat$cmb <- factor(paste(survDat$DeCAF, survDat$PurIST,  sep = "*"),
                    levels = c("permCAF*Basal-like","permCAF*Classical","restCAF*Basal-like","restCAF*Classical"))
km <- with(survDat, Surv(time,event))
p = coxph(km ~ cmb, data = survDat)
bic = round(BIC(p),3)
#km_fit <- survfit(km ~ cmb, data = survDat, type = "kaplan-meier")
km_fit <- survfit(km ~ DeCAF+PurIST, data = survDat, type = "kaplan-meier")
p <- ggsurvplot(font.legend = 16, 
                font.title = 16,        
                font.x = 16,  
                font.y = 16, 
                font.tickslab = 16,
                risk.table.fontsize = 7,
                risk.table.height = 0.35,
                size = 0.3, 
                censor.size=7, 
                pval.size = 8,
                pval.coord = c(96, 0.75),
                xlim = c(0,144),
                km_fit, conf.int = F, pval = T,
                legend.title="",
                break.time.by = 12,
                legend.labs=levels(survDat$cmb),
                palette = c("violetred1","#FE46A5","turquoise4","#008080"),
                linetype = c(1,2,1,2),
                xlab = "Time (months)", 
                risk.table = T,
                title = paste("YehSeq FOLFIRINOX DeCAF*PurIST",sep=""))
print(p)

# interaction
table(survDat[,c("PurIST","DeCAF")])
fisher.test(table(survDat[,c("PurIST","DeCAF")]))

# multi-variate
survDat$PurIST <- factor(survDat$PurIST, levels = c("Classical","Basal-like"))
survDat$DeCAF <- factor(survDat$DeCAF, levels = c("restCAF","permCAF"))
fit <- coxph(formula = Surv(time, event) ~ DeCAF+PurIST , data = survDat)
p <- ggforest(fit, 
              main = "YehSeq FOLFIRINOX",
              fontsize = 1.2,
              cpositions = c(0.02, 0.13, 0.32)) 
print(p)

dev.off()
