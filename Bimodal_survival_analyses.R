####### HOXA9 BIMODALITY AND SURVIVAL ANALYSES ############


library(multimode)

TPM = read.csv("log_modified_LAML_TPM.csv",sep=",",h=T)
TPM$X <- NULL
row.names(TPM) <- TPM$Hybridization.REF
TPM$Hybridization.REF<- NULL

HOXA9 <- as.numeric(TPM[4607,])
hist(HOXA9,cex=2,cex.lab=2,cex.main=2,xlab="log(HOXA9 TPM+1)",ylab="Frequency (patients)",main="",cex.axis=2) 
modetest(HOXA9)
locmodes(HOXA9, 2, display=T) 
# Is HOXA9 bimodal?
modetest(HOXA9, mod0=2, B=100) 



##### SURVIVAL ANALYSES

library(survival)
library(ggplot2)
library(survminer)


data = read.csv("HOXA9_survival.csv",sep=",",h=T)

kmsurvival <- survfit(Surv(OS,vital2)~HOXA9, data = data)
par(mar=c(5,6,4,1)+.1)
ggsurv <- ggsurvplot(kmsurvival, data = data,xlab="Time (Months)", pval=T, risk.table = TRUE, legend.title = "HOXA9", legend.labs = c("High", "Low"), font.x = 18,font.y=18)
ggsurv$plot <- ggsurv$plot + 
  theme(legend.text = element_text(size = 14, color = "black", face = "bold"), legend.title = element_text(size = 14, color = "black", face = "bold"))
ggsurv
summary(kmsurvival)

survdiff(Surv(data2$OS,data2$vital2)~data2$HOXA9, data = data2,rho=0)

####### FOR M2 AND M4 FAB SUBTYPES IN HOXA9 COHORTS

M2_data = read.csv("M2_HOXA9_survival.csv",sep=",",h=T)

kmsurvival_M2 <- survfit(Surv(OS,vital2)~HOXA9, data = M2_data)
par(mar=c(5,6,4,1)+.1)

ggsurv <- ggsurvplot(kmsurvival_M2, data = M2_data,xlab="Time (Months)", pval=T, risk.table = TRUE, legend.title = "FAB M2 HOXA9", legend.labs = c("High", "Low"), font.x = 18,font.y=18)
ggsurv$plot <- ggsurv$plot + 
  theme(legend.text = element_text(size = 14, color = "black", face = "bold"), legend.title = element_text(size = 14, color = "black", face = "bold"))
ggsurv

summary(kmsurvival_M2)

##

M4_data = read.csv("M4_HOXA9_survival.csv",sep=",",h=T)

kmsurvival_M4 <- survfit(Surv(OS,vital2)~HOXA9, data = M4_data)
par(mar=c(5,6,4,1)+.1)

ggsurv <- ggsurvplot(kmsurvival_M4, data = M4_data,xlab="Time (Months)", pval=T, risk.table = TRUE, legend.title = "FAB M4 HOXA9", legend.labs = c("High", "Low"), font.x = 18,font.y=18)
ggsurv$plot <- ggsurv$plot + 
  theme(legend.text = element_text(size = 14, color = "black", face = "bold"), legend.title = element_text(size = 14, color = "black", face = "bold"))
ggsurv

summary(kmsurvival_M4)

### for hoxa9 with age categories

data2 = read.csv("HOXA9_survival_age.csv",sep=",",h=T)

kmsurvival2 <- survfit(Surv(OS,vital2)~HOXA9, data = data2)
par(mar=c(5,6,4,1)+.1)
ggsurv <- ggsurvplot(kmsurvival2, data = data2,xlab="Time (Months)", pval=T, risk.table = TRUE, legend.title = "HOXA9", legend.labs = c("High", "Low"), font.x = 18,font.y=18)
ggsurv$plot <- ggsurv$plot + 
  theme(legend.text = element_text(size = 14, color = "black", face = "bold"), legend.title = element_text(size = 14, color = "black", face = "bold"))
ggsurv

summary(kmsurvival2)
