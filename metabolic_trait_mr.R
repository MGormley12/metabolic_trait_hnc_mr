setwd("")

install.packages("devtools")
devtools::install_github("MRCIEU/TwoSampleMR")
devtools::install_github("MRCIEU/MRInstruments")
install.packages("devtools")
devtools::install_github("MRCIEU/MRInstruments", force = TRUE)
devtools::install_github('MRCIEU/TwoSampleMR')
install.packages("MendelianRandomization", force = TRUE)
install.packages("LDlinkR")
install.packages("plyr")
install.packages("ggplot2")
install.packages("ggpubr")
install.packages("simex")
devtools::install_github("rondolab/MR-PRESSO")

library(devtools)
library(TwoSampleMR)
library(MRInstruments)
library(plyr) 
library(ggplot2)
library(MendelianRandomization)
library(gridExtra)
library(grid)
library(lattice)
library(LDlinkR)
library(ggpubr)
library(simex)
library(MRPRESSO)
ao <- available_outcomes()

# Power calculations

# GAME-ON oral/oropharyngeal cancer GWAS only- 6034 cases, 6585 controls
n <- (6034+6585) 
ratio <- 6034/6585

sig <- 0.05

Betas <- seq(from=0, to=0.5, by=0.0005)
Powers <- as.data.frame(Betas)
Powers$ORs <- exp(Betas)

Powers$Var0.5 <- (pnorm(sqrt(n*0.005*(ratio/(1+ratio))*(1/(1+ratio)))*Powers$Betas-qnorm(1-sig/2)))*100
Powers$Var1 <- (pnorm(sqrt(n*0.01*(ratio/(1+ratio))*(1/(1+ratio)))*Powers$Betas-qnorm(1-sig/2)))*100
Powers$Var2.5 <- (pnorm(sqrt(n*0.025*(ratio/(1+ratio))*(1/(1+ratio)))*Powers$Betas-qnorm(1-sig/2)))*100
Powers$Var5 <- (pnorm(sqrt(n*0.05*(ratio/(1+ratio))*(1/(1+ratio)))*Powers$Betas-qnorm(1-sig/2)))*100
Powers$Var10 <- (pnorm(sqrt(n*0.1*(ratio/(1+ratio))*(1/(1+ratio)))*Powers$Betas-qnorm(1-sig/2)))*100

PowerPlot <- ggplot(Powers, aes(ORs)) +
  geom_line(aes(y = Var0.5, colour = "9")) +  
  geom_line(aes(y = Var1, colour = "7")) +
  geom_line(aes(y = Var2.5, colour = "5")) + 
  geom_line(aes(y = Var5, colour = "3")) +
  geom_line(aes(y = Var10, colour = "1")) +
  xlab("Odds ratio per unit increase in risk factor") +
  scale_y_continuous("Power (%)", limits=c(0, 100), breaks=c(20,40,60,80,100)) +
  theme(axis.title.x = element_text(face="bold", size=20), axis.text.x  = element_text(vjust=0.5, size=16)) +
  theme(axis.title.y = element_text(face="bold", size=20), axis.text.y  = element_text(vjust=0.5, size=16)) +
  theme(legend.text=element_text(size=12), legend.title=element_text(size=12)) +
  scale_colour_discrete("% Variance", labels= c("10.0", "5.0", "2.5", "1.0", "0.5")) +
  geom_hline(yintercept=70) +
  ggtitle("Power for analysing all GAME-ON cases & controls") +
  theme(plot.title=element_text(lineheight=5, size= rel(1.2), face="bold")) +
  theme_classic()

dest.plot <- "gameon_power_2SMR.png"

png(dest.plot, width = 10*500, height = 5*500, res=500)
PowerPlot
dev.off()


# GAME-ON oral cancer GWAS only-  2990 cases,  6585 controls
n <- (2990+6585) 
ratio <- 2990/6585

sig <- 0.05

Betas <- seq(from=0, to=0.5, by=0.0005)
Powers <- as.data.frame(Betas)
Powers$ORs <- exp(Betas)

Powers$Var0.5 <- (pnorm(sqrt(n*0.005*(ratio/(1+ratio))*(1/(1+ratio)))*Powers$Betas-qnorm(1-sig/2)))*100
Powers$Var1 <- (pnorm(sqrt(n*0.01*(ratio/(1+ratio))*(1/(1+ratio)))*Powers$Betas-qnorm(1-sig/2)))*100
Powers$Var2.5 <- (pnorm(sqrt(n*0.025*(ratio/(1+ratio))*(1/(1+ratio)))*Powers$Betas-qnorm(1-sig/2)))*100
Powers$Var5 <- (pnorm(sqrt(n*0.05*(ratio/(1+ratio))*(1/(1+ratio)))*Powers$Betas-qnorm(1-sig/2)))*100
Powers$Var10 <- (pnorm(sqrt(n*0.1*(ratio/(1+ratio))*(1/(1+ratio)))*Powers$Betas-qnorm(1-sig/2)))*100

PowerPlot <- ggplot(Powers, aes(ORs)) +
  geom_line(aes(y = Var0.5, colour = "9")) +  
  geom_line(aes(y = Var1, colour = "7")) +
  geom_line(aes(y = Var2.5, colour = "5")) + 
  geom_line(aes(y = Var5, colour = "3")) +
  geom_line(aes(y = Var10, colour = "1")) +
  xlab("Odds ratio per unit increase in risk factor") +
  scale_y_continuous("Power (%)", limits=c(0, 100), breaks=c(20,40,60,80,100)) +
  theme(axis.title.x = element_text(face="bold", size=20), axis.text.x  = element_text(vjust=0.5, size=16)) +
  theme(axis.title.y = element_text(face="bold", size=20), axis.text.y  = element_text(vjust=0.5, size=16)) +
  theme(legend.text=element_text(size=12), legend.title=element_text(size=12)) +
  scale_colour_discrete("% Variance", labels= c("10.0", "5.0", "2.5", "1.0", "0.5")) +
  geom_hline(yintercept=70) +
  ggtitle("Power for analysing GAME-ON oral cancer cases & controls") +
  theme(plot.title=element_text(lineheight=5, size= rel(1.2), face="bold")) +
  theme_classic()

dest.plot <- "gameon_oc_power_2SMR.png"

png(dest.plot, width = 10*500, height = 5*500, res=500)
PowerPlot
dev.off()

##GAME-ON oropharyngeal cancer GWAS only- 2641 cases, 6585 controls
n <- (2641+6585) 
ratio <- 2641/6585

sig <- 0.05

Betas <- seq(from=0, to=0.5, by=0.0005)
Powers <- as.data.frame(Betas)
Powers$ORs <- exp(Betas)

Powers$Var0.5 <- (pnorm(sqrt(n*0.005*(ratio/(1+ratio))*(1/(1+ratio)))*Powers$Betas-qnorm(1-sig/2)))*100
Powers$Var1 <- (pnorm(sqrt(n*0.01*(ratio/(1+ratio))*(1/(1+ratio)))*Powers$Betas-qnorm(1-sig/2)))*100
Powers$Var2.5 <- (pnorm(sqrt(n*0.025*(ratio/(1+ratio))*(1/(1+ratio)))*Powers$Betas-qnorm(1-sig/2)))*100
Powers$Var5 <- (pnorm(sqrt(n*0.05*(ratio/(1+ratio))*(1/(1+ratio)))*Powers$Betas-qnorm(1-sig/2)))*100
Powers$Var10 <- (pnorm(sqrt(n*0.1*(ratio/(1+ratio))*(1/(1+ratio)))*Powers$Betas-qnorm(1-sig/2)))*100

PowerPlot <- ggplot(Powers, aes(ORs)) +
  geom_line(aes(y = Var0.5, colour = "9")) +  
  geom_line(aes(y = Var1, colour = "7")) +
  geom_line(aes(y = Var2.5, colour = "5")) + 
  geom_line(aes(y = Var5, colour = "3")) +
  geom_line(aes(y = Var10, colour = "1")) +
  xlab("Odds ratio per unit increase in risk factor") +
  scale_y_continuous("Power (%)", limits=c(0, 100), breaks=c(20,40,60,80,100)) +
  theme(axis.title.x = element_text(face="bold", size=20), axis.text.x  = element_text(vjust=0.5, size=16)) +
  theme(axis.title.y = element_text(face="bold", size=20), axis.text.y  = element_text(vjust=0.5, size=16)) +
  theme(legend.text=element_text(size=12), legend.title=element_text(size=12)) +
  scale_colour_discrete("% Variance", labels= c("10.0", "5.0", "2.5", "1.0", "0.5")) +
  geom_hline(yintercept=70) +
  ggtitle("Power for analysing GAME-ON oropharyngeal cancer cases & controls") +
  theme(plot.title=element_text(lineheight=5, size= rel(1.2), face="bold")) +
  theme_classic()

dest.plot <- "gameon_opc_power_2SMR.png"

png(dest.plot, width = 10*500, height = 5*500, res=500)
PowerPlot
dev.off()


# Body mass index (BMI) on oral and oropharngeal cancer

exposure_dat <- extract_instruments(outcomes='ieu-b-40')
write.csv(exposure_dat,"exposure_dat_bmi.csv")

exposure_dat <- read_exposure_data(
  "exposure_dat_bmi.csv",
  clump = TRUE,
  sep = ",",
  phenotype_col = "dexposure",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  eaf_col = "eaf.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  pval_col = "pval.exposure",
  min_pval = 1e-200,
  log_pval = FALSE
)
write.csv(exposure_dat,'clumped_bmi_exposure.csv')

#Extract SNPs for BMI from Yengo et al. from dbGaP or derive from IEU Open GWAS
outcome_dat <- read_outcome_data("bmi_hnc.txt", sep="\t")
dat <- harmonise_data(exposure_dat, outcome_dat)
dat$outcome <- "Oral and oropharyngeal cancer"
dat$exposure <- "Body Mass Index"
write.csv(dat, "bmi_hnc_harmonised.csv")    
mr_results <- mr(dat)
mr_results
or_results <- generate_odds_ratios(mr_results)
or_results
mr_report(dat, output_path = "~/OneDrive - University of Bristol/Documents/Documents/GW4-CAT PhD/Research and publications/Metabolic_Phenotypes_MR/bmi_mr",author="Mark Gormley", study = paste("Body Mass Index", "Oral and oropharyngeal cancer",sep=""))
results<-cbind.data.frame(or_results$outcome,or_results$nsnp,or_results$method,or_results$b,or_results$se,or_results$pval,or_results$or,or_results$or_lci95,or_results$or_uci95)
write.table(results, "bmi_hnc_results.txt")

#Run scatter plot 
mr_scatter <- mr_scatter_plot(mr_results, dat)
mr_scatter
pdf("bmi_hnc_scatter.pdf",  width = 15, height = 20)
ggarrange(mr_scatter$`hRMkq6.zLc7p5`,ncol=2, nrow=2, widths = 2, heights = 1)
dev.off()

#6. PLOT FOREST AND LEAVE-ONE-OUT
res_single <- mr_singlesnp(dat)
mr_forest <- mr_forest_plot(res_single)
mr_forest

pdf("bmi_hnc_forest.pdf",  width = 15, height = 20)
ggarrange(mr_forest$`hRMkq6.zLc7p5`,ncol=2, nrow=2, widths = 2, heights = 1)
dev.off()

#Leave one out analysis
res_loo <- mr_leaveoneout(dat)
mr_loo <- mr_leaveoneout_plot(res_loo)
mr_loo

pdf("bmi_hnc_loo.pdf",  width = 15, height = 20)
ggarrange(mr_loo$`hRMkq6.zLc7p5`, nrow=2, widths = 2, heights = 1)
dev.off()

#Heterogeneity and pleiotropy analysis 
Het<-mr_heterogeneity(dat)
Plt<-mr_pleiotropy_test(dat)
Sin<-mr_singlesnp(dat)
write.table(Het, "bmi_hnc_heterogeneity.txt")
write.table(Plt, "bmi_hnc_pleiotropy.txt")
write.table(Sin, "bmi_hnc_single_snp.txt")
p1 <- mr_scatter_plot(mr_results, dat)
p1[[1]]

res_single <- mr_singlesnp(dat)
p2 <- mr_forest_plot(res_single)
p2[[1]]

res_loo <- mr_leaveoneout(dat)
p3 <- mr_leaveoneout_plot(res_loo)
p3[[1]]

p4 <- mr_funnel_plot(res_single)
p4[[1]]

#Assess outliers
#Radial plots 
devtools::install_github("WSpiller/RadialMR")
library(RadialMR)
library(MendelianRandomization)

dat <- dat[dat$SNP%in%res_single$SNP,]

raddat <- format_radial(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome, dat$SNP)
ivwrad <- ivw_radial(raddat, alpha=0.05/63, weights=3)
dim(ivwrad$outliers)[1] 
ivwrad <- ivw_radial(raddat, alpha=0.05, weights=3)
dim(ivwrad$outliers)[1] 
eggrad <- egger_radial(raddat, alpha=0.05, weights=3)
eggrad$coef 
dim(eggrad$outliers)[1] 

#plot_radial(ivwrad, TRUE, FALSE, TRUE)
plot_radial(c(ivwrad,eggrad), TRUE, FALSE, TRUE)

ivwrad$qstatistic 
ivwrad$sortoutliers <- ivwrad$outliers[order(ivwrad$outliers$p.value),]
ivwrad$sortoutliers$Qsum <- cumsum(ivwrad$sortoutliers$Q_statistic)
ivwrad$sortoutliers$Qdif <- ivwrad$sortoutliers$Qsum - ivwrad$qstatistic
write.csv(ivwrad$sortoutliers, "bmi_outliers.csv", row.names=F, quote=F)

#Remove top outliers
dat2 <- dat[!dat$SNP %in% ivwrad$outliers$SNP,]
mr_results2 <- mr(dat2)
or_results <- generate_odds_ratios(mr_results2)

results<-cbind.data.frame(or_results$outcome,or_results$nsnp,or_results$method,or_results$b,or_results$se,or_results$pval,or_results$or,or_results$or_lci95,or_results$or_uci95)
write.table(results, "bmi_hnc_nooutliers.txt")

#MR presso 
devtools::install_github("rondolab/MR-PRESSO")
library(MRPRESSO)
mr_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat, NbDistribution = 1000,  SignifThreshold = 0.05)
mr_presso
write.table(mr_presso, "bmi_mr_presso.txt")

#Calculate I-squared, r-squared and F-statistics
#I-squared function
Isq <- function(y,s){
  k          = length(y)
  w          = 1/s^2; sum.w  = sum(w)
  mu.hat     = sum(y*w)/sum.w  
  Q          = sum(w*(y-mu.hat)^2)
  Isq        = (Q - (k-1))/Q
  Isq        = max(0,Isq)
  return(Isq)
}

#Calculate Isq weighted and unweighted
I2<-c()
dat <- harmonise_data(exposure_dat, outcome_dat, action = 1)
str(dat)

#F-statistic
dat$samplesize.exposure <- 681275
dat$samplesize.outcome <- 6034
dat <- steiger_filtering(dat) 

N = dat$samplesize.exposure[1] #sample size
K = length(dat$SNP) #number of SNPs
total_r2 <- sum(dat$rsq.exposure) 
Fstat <- (N-K-1)/K * total_r2 / (1-total_r2)
total_r2
Fstat

#Rename required columns
dat $BetaXG<-dat $beta.exposure
dat $seBetaXG<-dat $se.exposure
BetaXG   = dat $BetaXG
seBetaXG = dat $seBetaXG 
seBetaYG<-dat $se.outcome

BXG = abs(BetaXG)         
# Calculate F statistics
# and I-squared statistics
# to measure Instrument 
# strength for MR-Egger

F = BXG^2/seBetaXG^2
mF = mean(F)
Isq_unweighted <- Isq(BXG,seBetaXG) #unweighted
Isq_weighted <- Isq((BXG/seBetaYG),(seBetaXG/seBetaYG)) #weighted

#Save results
output<-cbind (F, mF, Isq_unweighted, Isq_weighted)
I2<-rbind(I2, output)
colnames(I2) <- c("F", "mF", "Isq_unweighted", "Isq_weighted")
write.csv(I2, file="regression_dilution_isq_weighted_bmi.csv", row.names = FALSE)

# OPC and OC subsite analysis
exposure_dat <- read_exposure_data(
  "exposure_dat_bmi.csv",
  clump = TRUE,
  sep = ",",
  phenotype_col = "dexposure",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  eaf_col = "eaf.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  pval_col = "pval.exposure",
  min_pval = 1e-200,
  log_pval = FALSE
)


outcome_dat <- read_outcome_data("bmi_oc.txt", sep="\t")
dat <- harmonise_data(exposure_dat, outcome_dat)
write.csv(dat, "bmi_oc_harmonised.csv")    
mr <- mr(dat)
or_results <- generate_odds_ratios(mr)
or_results
write.table(or_results, "bmi_oc_results.txt")

outcome_dat <- read_outcome_data("bmi_opc.txt", sep="\t")
dat <- harmonise_data(exposure_dat, outcome_dat)
write.csv(dat, "bmi_opc_harmonised.csv")    
mr <- mr(dat)
or_results <- generate_odds_ratios(mr)
or_results
write.table(or_results, "bmi_opc_results.txt")


# Waist Circumference (WC) on oral and oropharyngeal cancer

exposure_dat <- extract_instruments(outcomes='ieu-a-60')
write.csv(exposure_dat,"exposure_dat_wc.csv")

exposure_dat <- read_exposure_data(
  "exposure_dat_wc.csv",
  clump = TRUE,
  sep = ",",
  phenotype_col = "exposure",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  eaf_col = "eaf.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  pval_col = "p-value",
  min_pval = 1e-200,
  log_pval = FALSE
)
write.csv(exposure_dat,'clumped_wc_exposure.csv')


#Extract SNPs for WC from Shungin et al. from dbGaP or derive from IEU Open GWAS
outcome_dat <- read_outcome_data("wc_hnc.txt", sep="\t")
dat <- harmonise_data(exposure_dat, outcome_dat)
dat$outcome <- "Oral and oropharyngeal cancer"
dat$exposure <- "Waist Circumference"
write.csv(dat, "wc_hnc_harmonised.csv")    
mr_results <- mr(dat)
mr_results
or_results <- generate_odds_ratios(mr_results)
or_results
mr_report(dat, output_path = "~OneDrive - University of Bristol/Documents/Documents/GW4-CAT PhD/Research and publications/Metabolic_Phenotypes_MR/wc_mr",author="Mark Gormley", study = paste("Waist Circumference", "Oral and oropharyngeal cancer",sep=""))
results<-cbind.data.frame(or_results$outcome,or_results$nsnp,or_results$method,or_results$b,or_results$se,or_results$pval,or_results$or,or_results$or_lci95,or_results$or_uci95)
write.table(results, "wc_hnc_results.txt")

#Run scatter plot
mr_scatter <- mr_scatter_plot(mr_results, dat)
mr_scatter
pdf("wc_hnc_scatter.pdf",  width = 15, height = 20)
ggarrange(mr_scatter$`cGVCgr.O1cTGs`,ncol=2, nrow=2, widths = 2, heights = 1)
dev.off()

#6. PLOT FOREST AND LEAVE-ONE-OUT
res_single <- mr_singlesnp(dat)
mr_forest <- mr_forest_plot(res_single)
mr_forest

pdf("wc_hnc_forest.pdf",  width = 15, height = 20)
ggarrange(mr_forest$`cGVCgr.O1cTGs`,ncol=2, nrow=2, widths = 2, heights = 1)
dev.off()

#Leave one out analysis
res_loo <- mr_leaveoneout(dat)
mr_loo <- mr_leaveoneout_plot(res_loo)
mr_loo

pdf("wc_hnc_loo.pdf",  width = 15, height = 20)
ggarrange(mr_loo$`cGVCgr.O1cTGs`, nrow=2, widths = 2, heights = 1)
dev.off()

#Heterogeneity and pleiotropy analysis 
Het<-mr_heterogeneity(dat)
Plt<-mr_pleiotropy_test(dat)
Sin<-mr_singlesnp(dat)
write.table(Het, "wc_hnc_heterogeneity.txt")
write.table(Plt, "wc_hnc_pleiotropy.txt")
write.table(Sin, "wc_hnc_single_snp.txt")
p1 <- mr_scatter_plot(mr_results, dat)
p1[[1]]

res_single <- mr_singlesnp(dat)
p2 <- mr_forest_plot(res_single)
p2[[1]]

res_loo <- mr_leaveoneout(dat)
p3 <- mr_leaveoneout_plot(res_loo)
p3[[1]]

p4 <- mr_funnel_plot(res_single)
p4[[1]]

#Assess outliers
#Radial plots 
devtools::install_github("WSpiller/RadialMR")
library(RadialMR)
library(MendelianRandomization)

dat <- dat[dat$SNP%in%res_single$SNP,]

raddat <- format_radial(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome, dat$SNP)
ivwrad <- ivw_radial(raddat, alpha=0.05/63, weights=3)
dim(ivwrad$outliers)[1] 
ivwrad <- ivw_radial(raddat, alpha=0.05, weights=3)
dim(ivwrad$outliers)[1] 
eggrad <- egger_radial(raddat, alpha=0.05, weights=3)
eggrad$coef 
dim(eggrad$outliers)[1]

#plot_radial(ivwrad, TRUE, FALSE, TRUE)
plot_radial(c(ivwrad,eggrad), TRUE, FALSE, TRUE)

ivwrad$qstatistic 
ivwrad$sortoutliers <- ivwrad$outliers[order(ivwrad$outliers$p.value),]
ivwrad$sortoutliers$Qsum <- cumsum(ivwrad$sortoutliers$Q_statistic)
ivwrad$sortoutliers$Qdif <- ivwrad$sortoutliers$Qsum - ivwrad$qstatistic
write.csv(ivwrad$sortoutliers, "wc_outliers.csv", row.names=F, quote=F)

#Remove top outliers
dat2 <- dat[!dat$SNP %in% ivwrad$outliers$SNP,]
mr_results2 <- mr(dat2)
or_results <- generate_odds_ratios(mr_results2)

results<-cbind.data.frame(or_results$outcome,or_results$nsnp,or_results$method,or_results$b,or_results$se,or_results$pval,or_results$or,or_results$or_lci95,or_results$or_uci95)
write.table(results, "wc_hnc_nooutliers.txt")

#MR presso 
devtools::install_github("rondolab/MR-PRESSO")
library(MRPRESSO)
mr_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat, NbDistribution = 1000,  SignifThreshold = 0.05)
mr_presso
write.table(mr_presso, "wc_mr_presso.txt")

#Calculate I-squared, r-squared and F-statistics
#I-squared function
Isq <- function(y,s){
  k          = length(y)
  w          = 1/s^2; sum.w  = sum(w)
  mu.hat     = sum(y*w)/sum.w  
  Q          = sum(w*(y-mu.hat)^2)
  Isq        = (Q - (k-1))/Q
  Isq        = max(0,Isq)
  return(Isq)
}

#Calculate Isq weighted and unweighted
I2<-c()
dat <- harmonise_data(exposure_dat, outcome_dat, action = 1)
str(dat)

#F-statistic
dat$samplesize.exposure <- 224459
dat$samplesize.outcome <- 6034
dat <- steiger_filtering(dat) 

N = dat$samplesize.exposure[1] #sample size
K = length(dat$SNP) #number of SNPs
total_r2 <- sum(dat$rsq.exposure) 
Fstat <- (N-K-1)/K * total_r2 / (1-total_r2)
total_r2
Fstat

#Rename required columns
dat $BetaXG<-dat $beta.exposure
dat $seBetaXG<-dat $se.exposure
BetaXG   = dat $BetaXG
seBetaXG = dat $seBetaXG 
seBetaYG<-dat $se.outcome

BXG = abs(BetaXG)         
# Calculate F statistics
# and I-squared statistics
# to measure Instrument 
# strength for MR-Egger

F = BXG^2/seBetaXG^2
mF = mean(F)
Isq_unweighted <- Isq(BXG,seBetaXG) #unweighted
Isq_weighted <- Isq((BXG/seBetaYG),(seBetaXG/seBetaYG)) #weighted

#Save results
output<-cbind (F, mF, Isq_unweighted, Isq_weighted)
I2<-rbind(I2, output)
colnames(I2) <- c("F", "mF", "Isq_unweighted", "Isq_weighted")
write.csv(I2, file="regression_dilution_isq_weighted_wc.csv", row.names = FALSE)


#SIMEX correction
#Create empty dataframe to store output
simexegger<-c()

#Run SIMEX
#Rename required columns
dat$BetaXG<-dat$beta.exposure
dat$seBetaXG<-dat$se.exposure
dat$BetaYG<-dat$beta.outcome
dat$seBetaYG<-dat$se.outcome
BetaXG <- dat$BetaXG
BetaYG <- dat$BetaYG
seBetaXG <- dat$seBetaXG
seBetaYG <- dat$seBetaYG

BYG <- BetaYG*sign(BetaXG)#Pre-processing steps to ensure all gene--exposure estimates are positive
BXG <- abs(BetaXG)         

#MR-Egger regression (weighted) 
Fit1 <- lm(BYG ~ BXG,weights=1/seBetaYG^2,x=TRUE,y=TRUE)

#MR-Egger regression (unweighted)
Fit2 <- lm(BYG~BXG,x=TRUE,y=TRUE) 

#Simulation extrapolation 
mod.sim1 <- simex(Fit1,B=1000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE") 
mod.sim2 <- simex(Fit2,B=1000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE") 
mod1<-summary(mod.sim1)
mod2<-summary(mod.sim2)
#Print results
mod1
mod2
#Use BXG result and exponentiate the estimates
OR = exp(-0.440179) 
CIL = exp(-0.440179 - 1.96 * 0.472805) 
CIU = exp(-0.440179 + 1.96 * 0.472805) 


# OPC and OC subsite analysis
exposure_dat <- extract_instruments(outcomes='ieu-a-60')
write.csv(exposure_dat,"exposure_dat_wc.csv")

exposure_dat <- read_exposure_data(
  "exposure_dat_wc.csv",
  clump = TRUE,
  sep = ",",
  phenotype_col = "exposure",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  eaf_col = "eaf.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  pval_col = "p-value",
  min_pval = 1e-200,
  log_pval = FALSE
)

outcome_dat <- read_outcome_data("wc_oc.txt", sep="\t")
dat <- harmonise_data(exposure_dat, outcome_dat)
write.csv(dat, "wc_oc_harmonised.csv")    
mr <- mr(dat)
or_results <- generate_odds_ratios(mr)
or_results
write.table(or_results, "wc_oc_results.txt")

outcome_dat <- read_outcome_data("wc_opc.txt", sep="\t")
dat <- harmonise_data(exposure_dat, outcome_dat)
write.csv(dat, "wc_opc_harmonised.csv")    
mr <- mr(dat)
or_results <- generate_odds_ratios(mr)
or_results
write.table(or_results, "wc_opc_results.txt")


# Waist-to-hip ratio (WHR) on oral and oropharyngeal cancer

exposure_dat <- read_exposure_data("exposure_dat_whr.txt", sep="\t")
# already clumped
write.csv(exposure_dat,'clumped_bmi_exposure.csv')

#Extract SNPs for WHR from Shungin et al. from dbGaP or derive from IEU Open GWAS
outcome_dat <- read_outcome_data("whr_hnc.txt", sep="\t")
dat <- harmonise_data(exposure_dat, outcome_dat)
dat$outcome <- "Oral and oropharyngeal cancer"
dat$exposure <- "Waist-to-hip ratio"
write.csv(dat, "whr_hnc_harmonised.csv")    
mr_results <- mr(dat)
mr_results
or_results <- generate_odds_ratios(mr_results)
or_results
mr_report(dat, output_path = "~OneDrive - University of Bristol/Documents/Documents/GW4-CAT PhD/Research and publications/Metabolic_Phenotypes_MR/whr_mr",author="Mark Gormley", study = paste("Waist-to-hip ratio", "Oral and oropharyngeal cancer",sep=""))
results<-cbind.data.frame(or_results$outcome,or_results$nsnp,or_results$method,or_results$b,or_results$se,or_results$pval,or_results$or,or_results$or_lci95,or_results$or_uci95)
write.table(results, "whr_hnc_results.txt")

#Run scatter plot
mr_scatter <- mr_scatter_plot(mr_results, dat)
mr_scatter
pdf("whr_hnc_scatter.pdf",  width = 15, height = 20)
ggarrange(mr_scatter$`AE3WXe.x1ucWE`,ncol=2, nrow=2, widths = 2, heights = 1)
dev.off()

#6. PLOT FOREST AND LEAVE-ONE-OUT
res_single <- mr_singlesnp(dat)
mr_forest <- mr_forest_plot(res_single)
mr_forest

pdf("whr_hnc_forest.pdf",  width = 15, height = 20)
ggarrange(mr_forest$`AE3WXe.x1ucWE`,ncol=2, nrow=2, widths = 2, heights = 1)
dev.off()

#Leave one out analysis
res_loo <- mr_leaveoneout(dat)
mr_loo <- mr_leaveoneout_plot(res_loo)
mr_loo

pdf("whr_hnc_loo.pdf",  width = 15, height = 20)
ggarrange(mr_loo$`AE3WXe.x1ucWE`, nrow=2, widths = 2, heights = 1)
dev.off()

#Heterogeneity and pleiotropy analysis 
Het<-mr_heterogeneity(dat)
Plt<-mr_pleiotropy_test(dat)
Sin<-mr_singlesnp(dat)
write.table(Het, "whr_hnc_heterogeneity.txt")
write.table(Plt, "whr_hnc_pleiotropy.txt")
write.table(Sin, "whr_hnc_single_snp.txt")
p1 <- mr_scatter_plot(mr_results, dat)
p1[[1]]

res_single <- mr_singlesnp(dat)
p2 <- mr_forest_plot(res_single)
p2[[1]]

res_loo <- mr_leaveoneout(dat)
p3 <- mr_leaveoneout_plot(res_loo)
p3[[1]]

p4 <- mr_funnel_plot(res_single)
p4[[1]]

#Assess outliers
#Radial plots 
devtools::install_github("WSpiller/RadialMR")
library(RadialMR)
library(MendelianRandomization)

dat <- dat[dat$SNP%in%res_single$SNP,]

raddat <- format_radial(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome, dat$SNP)
ivwrad <- ivw_radial(raddat, alpha=0.05/63, weights=3)
dim(ivwrad$outliers)[1]
ivwrad <- ivw_radial(raddat, alpha=0.05, weights=3)
dim(ivwrad$outliers)[1] 
eggrad <- egger_radial(raddat, alpha=0.05, weights=3)
eggrad$coef 
dim(eggrad$outliers)[1]

#plot_radial(ivwrad, TRUE, FALSE, TRUE)
plot_radial(c(ivwrad,eggrad), TRUE, FALSE, TRUE)

ivwrad$qstatistic 
ivwrad$sortoutliers <- ivwrad$outliers[order(ivwrad$outliers$p.value),]
ivwrad$sortoutliers$Qsum <- cumsum(ivwrad$sortoutliers$Q_statistic)
ivwrad$sortoutliers$Qdif <- ivwrad$sortoutliers$Qsum - ivwrad$qstatistic
write.csv(ivwrad$sortoutliers, "whr_outliers.csv", row.names=F, quote=F)

#Remove top outliers
dat2 <- dat[!dat$SNP %in% ivwrad$outliers$SNP,]
mr_results2 <- mr(dat2)
or_results <- generate_odds_ratios(mr_results2)

results<-cbind.data.frame(or_results$outcome,or_results$nsnp,or_results$method,or_results$b,or_results$se,or_results$pval,or_results$or,or_results$or_lci95,or_results$or_uci95)
write.table(results, "whr_hnc_nooutliers.txt")

#MR presso 
devtools::install_github("rondolab/MR-PRESSO")
library(MRPRESSO)
mr_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat, NbDistribution = 1000,  SignifThreshold = 0.05)
mr_presso
write.table(mr_presso, "whr_mr_presso.txt")

#Calculate I-squared, r-squared and F-statistics
#I-squared function
Isq <- function(y,s){
  k          = length(y)
  w          = 1/s^2; sum.w  = sum(w)
  mu.hat     = sum(y*w)/sum.w  
  Q          = sum(w*(y-mu.hat)^2)
  Isq        = (Q - (k-1))/Q
  Isq        = max(0,Isq)
  return(Isq)
}

#Calculate Isq weighted and unweighted
I2<-c()
dat <- harmonise_data(exposure_dat, outcome_dat, action = 1)
str(dat)

#F-statistic
dat$samplesize.exposure <- 694649
dat$samplesize.outcome <- 6034
dat <- steiger_filtering(dat) 

N = dat$samplesize.exposure[1] #sample size
K = length(dat$SNP) #number of SNPs
total_r2 <- sum(dat$rsq.exposure) 
Fstat <- (N-K-1)/K * total_r2 / (1-total_r2)
total_r2
Fstat

#Rename required columns
dat $BetaXG<-dat $beta.exposure
dat $seBetaXG<-dat $se.exposure
BetaXG   = dat $BetaXG
seBetaXG = dat $seBetaXG 
seBetaYG<-dat $se.outcome
BXG = abs(BetaXG)        

# Calculate F statistics
# and I-squared statistics
# to measure Instrument 
# strength for MR-Egger

F = BXG^2/seBetaXG^2
mF = mean(F)
Isq_unweighted <- Isq(BXG,seBetaXG) #unweighted
Isq_weighted <- Isq((BXG/seBetaYG),(seBetaXG/seBetaYG)) #weighted

#Save results
output<-cbind (F, mF, Isq_unweighted, Isq_weighted)
I2<-rbind(I2, output)
colnames(I2) <- c("F", "mF", "Isq_unweighted", "Isq_weighted")
write.csv(I2, file="regression_dilution_isq_weighted_whr.csv", row.names = FALSE)


#SIMEX correction
#Create empty dataframe to store output
simexegger<-c()

#Run SIMEX
#Rename required columns
dat$BetaXG<-dat$beta.exposure
dat$seBetaXG<-dat$se.exposure
dat$BetaYG<-dat$beta.outcome
dat$seBetaYG<-dat$se.outcome
BetaXG <- dat$BetaXG
BetaYG <- dat$BetaYG
seBetaXG <- dat$seBetaXG
seBetaYG <- dat$seBetaYG

BYG <- BetaYG*sign(BetaXG)#Pre-processing steps to ensure all gene--exposure estimates are positive
BXG <- abs(BetaXG)         

#MR-Egger regression (weighted) 
Fit1 <- lm(BYG ~ BXG,weights=1/seBetaYG^2,x=TRUE,y=TRUE)

#MR-Egger regression (unweighted)
Fit2 <- lm(BYG~BXG,x=TRUE,y=TRUE) 

#Simulation extrapolation 
mod.sim1 <- simex(Fit1,B=1000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE") 
mod.sim2 <- simex(Fit2,B=1000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE") 
mod1<-summary(mod.sim1)
mod2<-summary(mod.sim2)
#Print results
mod1
mod2
#Use BXG result and exponentiate the estimates
OR = exp(0.79961) 
CIL = exp(0.79961 - 1.96 * 1.62756) 
CIU = exp(0.79961 + 1.96 * 1.62756) 


# OPC and OC subsite analysis
exposure_dat <- read_exposure_data("exposure_dat_whr.txt", sep="\t")
# already clumped

outcome_dat <- read_outcome_data("whr_oc.txt", sep="\t")
dat <- harmonise_data(exposure_dat, outcome_dat)
write.csv(dat, "whr_oc_harmonised.csv")    
mr <- mr(dat)
or_results <- generate_odds_ratios(mr)
or_results
write.table(or_results, "whr_oc_results.txt")

outcome_dat <- read_outcome_data("whr_opc.txt", sep="\t")
dat <- harmonise_data(exposure_dat, outcome_dat)
write.csv(dat, "whr_opc_harmonised.csv")    
mr <- mr(dat)
or_results <- generate_odds_ratios(mr)
or_results
write.table(or_results, "whr_opc_results.txt")


#Type 2 diabetes (T2D) on oral and oropharyngeal cancer

exposure_dat <- read_exposure_data(
  "exposure_dat_t2d.csv",
  clump = TRUE,
  sep = ",",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  min_pval = 1e-200,
  log_pval = FALSE
)
write.csv(exposure_dat,'clumped_t2d_exposure.csv')

#Extract SNPs for T2D from Mahajan et al. from dbGaP or derive from IEU Open GWAS
outcome_dat <- read_outcome_data("t2d_hnc.txt", sep="\t")
dat <- harmonise_data(exposure_dat, outcome_dat)
dat$outcome <- "Oral and oropharyngeal cancer"
dat$exposure <- "Type 2 diabetes"
write.csv(dat, "t2d_hnc_harmonised.csv")    
mr_results <- mr(dat)
mr_results
or_results <- generate_odds_ratios(mr_results)
or_results
mr_report(dat, output_path = "~OneDrive - University of Bristol/Documents/Documents/GW4-CAT PhD/Research and publications/Metabolic_Phenotypes_MR/t2d_mr",author="Mark Gormley", study = paste("Type 2 diabetes", "Oral and oropharyngeal cancer",sep=""))
results<-cbind.data.frame(or_results$outcome,or_results$nsnp,or_results$method,or_results$b,or_results$se,or_results$pval,or_results$or,or_results$or_lci95,or_results$or_uci95)
write.table(results, "t2d_hnc_results.txt")

#Run scatter plot
mr_scatter <- mr_scatter_plot(mr_results, dat)
mr_scatter
pdf("t2d_hnc_scatter.pdf",  width = 15, height = 20)
ggarrange(mr_scatter$`r923KE.up5Vai`,ncol=2, nrow=2, widths = 2, heights = 1)
dev.off()

#6. PLOT FOREST AND LEAVE-ONE-OUT
res_single <- mr_singlesnp(dat)
mr_forest <- mr_forest_plot(res_single)
mr_forest
pdf("t2d_hnc_forest.pdf",  width = 15, height = 20)
ggarrange(mr_forest$`r923KE.up5Vai`,ncol=2, nrow=2, widths = 2, heights = 1)
dev.off()

#Leave one out analysis
res_loo <- mr_leaveoneout(dat)
mr_loo <- mr_leaveoneout_plot(res_loo)
mr_loo

pdf("t2d_hnc_loo.pdf",  width = 15, height = 20)
ggarrange(mr_loo$`r923KE.up5Vai`, nrow=2, widths = 2, heights = 1)
dev.off()

#Heterogeneity and pleiotropy analysis 
Het<-mr_heterogeneity(dat)
Plt<-mr_pleiotropy_test(dat)
Sin<-mr_singlesnp(dat)
write.table(Het, "t2d_hnc_heterogeneity.txt")
write.table(Plt, "t2d_hnc_pleiotropy.txt")
write.table(Sin, "t2d_hnc_single_snp.txt")
p1 <- mr_scatter_plot(mr_results, dat)
p1[[1]]

res_single <- mr_singlesnp(dat)
p2 <- mr_forest_plot(res_single)
p2[[1]]

res_loo <- mr_leaveoneout(dat)
p3 <- mr_leaveoneout_plot(res_loo)
p3[[1]]

p4 <- mr_funnel_plot(res_single)
p4[[1]]

#Assess outliers
#Radial plots 
devtools::install_github("WSpiller/RadialMR")
library(RadialMR)
library(MendelianRandomization)

dat <- dat[dat$SNP%in%res_single$SNP,]

raddat <- format_radial(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome, dat$SNP)
ivwrad <- ivw_radial(raddat, alpha=0.05/63, weights=3)
dim(ivwrad$outliers)[1]
ivwrad <- ivw_radial(raddat, alpha=0.05, weights=3)
dim(ivwrad$outliers)[1]
eggrad <- egger_radial(raddat, alpha=0.05, weights=3)
eggrad$coef 
dim(eggrad$outliers)[1]

#plot_radial(ivwrad, TRUE, FALSE, TRUE)
plot_radial(c(ivwrad,eggrad), TRUE, FALSE, TRUE)

ivwrad$qstatistic 
ivwrad$sortoutliers <- ivwrad$outliers[order(ivwrad$outliers$p.value),]
ivwrad$sortoutliers$Qsum <- cumsum(ivwrad$sortoutliers$Q_statistic)
ivwrad$sortoutliers$Qdif <- ivwrad$sortoutliers$Qsum - ivwrad$qstatistic
write.csv(ivwrad$sortoutliers, "t2d_outliers.csv", row.names=F, quote=F)

#Remove top outliers
dat2 <- dat[!dat$SNP %in% ivwrad$outliers$SNP,]
mr_results2 <- mr(dat2)
or_results <- generate_odds_ratios(mr_results2)

results<-cbind.data.frame(or_results$outcome,or_results$nsnp,or_results$method,or_results$b,or_results$se,or_results$pval,or_results$or,or_results$or_lci95,or_results$or_uci95)
write.table(results, "t2d_hnc_nooutliers.txt")

#MR presso 
devtools::install_github("rondolab/MR-PRESSO")
library(MRPRESSO)
mr_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat, NbDistribution = 1000,  SignifThreshold = 0.05)
mr_presso
write.table(mr_presso, "t2d_mr_presso.txt")

#Calculate I-squared, r-squared and F-statistics
#I-squared function
Isq <- function(y,s){
  k          = length(y)
  w          = 1/s^2; sum.w  = sum(w)
  mu.hat     = sum(y*w)/sum.w  
  Q          = sum(w*(y-mu.hat)^2)
  Isq        = (Q - (k-1))/Q
  Isq        = max(0,Isq)
  return(Isq)
}

#Calculate Isq weighted and unweighted
I2<-c()
dat <- harmonise_data(exposure_dat, outcome_dat, action = 1)
str(dat)

#F-statistic
dat$samplesize.exposure <- 1407282
dat$samplesize.outcome <- 6034
dat <- steiger_filtering(dat) 

N = dat$samplesize.exposure[1] #sample size
K = length(dat$SNP) #number of SNPs
total_r2 <- sum(dat$rsq.exposure) 
Fstat <- (N-K-1)/K * total_r2 / (1-total_r2)
total_r2
Fstat

#Rename required columns
dat $BetaXG<-dat $beta.exposure
dat $seBetaXG<-dat $se.exposure
BetaXG   = dat $BetaXG
seBetaXG = dat $seBetaXG 
seBetaYG<-dat $se.outcome
BXG = abs(BetaXG)        

# Calculate F statistics
# and I-squared statistics
# to measure Instrument 
# strength for MR-Egger

F = BXG^2/seBetaXG^2
mF = mean(F)
Isq_unweighted <- Isq(BXG,seBetaXG) #unweighted
Isq_weighted <- Isq((BXG/seBetaYG),(seBetaXG/seBetaYG)) #weighted

#Save results
output<-cbind (F, mF, Isq_unweighted, Isq_weighted)
I2<-rbind(I2, output)
colnames(I2) <- c("F", "mF", "Isq_unweighted", "Isq_weighted")
write.csv(I2, file="regression_dilution_isq_weighted_t2d.csv", row.names = FALSE)


# OPC and OC subsite analysis
exposure_dat <- read_exposure_data(
  "exposure_dat_t2d.csv",
  clump = TRUE,
  sep = ",",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  min_pval = 1e-200,
  log_pval = FALSE
)

outcome_dat <- read_outcome_data("t2d_oc.txt", sep="\t")
dat <- harmonise_data(exposure_dat, outcome_dat)
write.csv(dat, "t2d_oc_harmonised.csv")    
mr <- mr(dat)
or_results <- generate_odds_ratios(mr)
or_results
write.table(or_results, "t2d_oc_results.txt")

outcome_dat <- read_outcome_data("t2d_opc.txt", sep="\t")
dat <- harmonise_data(exposure_dat, outcome_dat)
write.csv(dat, "t2d_opc_harmonised.csv")    
mr <- mr(dat)
or_results <- generate_odds_ratios(mr)
or_results
write.table(or_results, "t2d_opc_results.txt")


#Glycated haemoglobin (HbA1c) on oral and oropharyngeal cancer

exposure_dat <- read_exposure_data(
  "exposure_dat_hba1c.txt",
  clump = TRUE,
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  min_pval = 1e-200,
  log_pval = FALSE
)
write.csv(exposure_dat,'clumped_hba1c_exposure.csv')

#Extract SNPs for HbA1c from Wheeler et al. from dbGaP or derive from IEU Open GWAS
outcome_dat <- read_outcome_data("hba1c_hnc.txt", sep="\t")
dat <- harmonise_data(exposure_dat, outcome_dat)
dat$outcome <- "Oral and oropharyngeal cancer"
dat$exposure <- "Glycated haemoglobin"
write.csv(dat, "hba1c_hnc_harmonised.csv")    
mr_results <- mr(dat)
mr_results
or_results <- generate_odds_ratios(mr_results)
or_results
mr_report(dat, output_path = "~OneDrive - University of Bristol/Documents/Documents/GW4-CAT PhD/Research and publications/Metabolic_Phenotypes_MR/hba1c_mr",author="Mark Gormley", study = paste("Glycated haemoglobin", "Oral and oropharyngeal cancer",sep=""))
results<-cbind.data.frame(or_results$outcome,or_results$nsnp,or_results$method,or_results$b,or_results$se,or_results$pval,or_results$or,or_results$or_lci95,or_results$or_uci95)
write.table(results, "hba1c_hnc_results.txt")

#Run scatter plot
mr_scatter <- mr_scatter_plot(mr_results, dat)
mr_scatter

pdf("hba1c_hnc_scatter.pdf",  width = 15, height = 20)
ggarrange(mr_scatter$`5BU4SZ.6w4dEU`,ncol=2, nrow=2, widths = 2, heights = 1)
dev.off()

#6. PLOT FOREST AND LEAVE-ONE-OUT
res_single <- mr_singlesnp(dat)
mr_forest <- mr_forest_plot(res_single)
mr_forest

pdf("hba1c_hnc_forest.pdf",  width = 15, height = 20)
ggarrange(mr_forest$`5BU4SZ.6w4dEU`,ncol=2, nrow=2, widths = 2, heights = 1)
dev.off()

#Leave one out analysis
res_loo <- mr_leaveoneout(dat)
mr_loo <- mr_leaveoneout_plot(res_loo)
mr_loo

pdf("hba1c_hnc_loo.pdf",  width = 15, height = 20)
ggarrange(mr_loo$`5BU4SZ.6w4dEU`, nrow=2, widths = 2, heights = 1)
dev.off()

#Heterogeneity and pleiotropy analysis 
Het<-mr_heterogeneity(dat)
Plt<-mr_pleiotropy_test(dat)
Sin<-mr_singlesnp(dat)
write.table(Het, "hba1c_hnc_heterogeneity.txt")
write.table(Plt, "hba1c_hnc_pleiotropy.txt")
write.table(Sin, "hba1c_hnc_single_snp.txt")
p1 <- mr_scatter_plot(mr_results, dat)
p1[[1]]

res_single <- mr_singlesnp(dat)
p2 <- mr_forest_plot(res_single)
p2[[1]]

res_loo <- mr_leaveoneout(dat)
p3 <- mr_leaveoneout_plot(res_loo)
p3[[1]]

p4 <- mr_funnel_plot(res_single)
p4[[1]]

#Assess outliers
#Radial plots 
devtools::install_github("WSpiller/RadialMR")
library(RadialMR)
library(MendelianRandomization)

dat <- dat[dat$SNP%in%res_single$SNP,]

raddat <- format_radial(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome, dat$SNP)
ivwrad <- ivw_radial(raddat, alpha=0.05/63, weights=3)
dim(ivwrad$outliers)[1]
ivwrad <- ivw_radial(raddat, alpha=0.05, weights=3)
dim(ivwrad$outliers)[1]
eggrad <- egger_radial(raddat, alpha=0.05, weights=3)
eggrad$coef 
dim(eggrad$outliers)[1]

#plot_radial(ivwrad, TRUE, FALSE, TRUE)
plot_radial(c(ivwrad,eggrad), TRUE, FALSE, TRUE)

ivwrad$qstatistic 
ivwrad$sortoutliers <- ivwrad$outliers[order(ivwrad$outliers$p.value),]
ivwrad$sortoutliers$Qsum <- cumsum(ivwrad$sortoutliers$Q_statistic)
ivwrad$sortoutliers$Qdif <- ivwrad$sortoutliers$Qsum - ivwrad$qstatistic
write.csv(ivwrad$sortoutliers, "hba1c_outliers.csv", row.names=F, quote=F)

#Remove top outliers
dat2 <- dat[!dat$SNP %in% ivwrad$outliers$SNP,]
mr_results2 <- mr(dat2)
or_results <- generate_odds_ratios(mr_results2)

results<-cbind.data.frame(or_results$outcome,or_results$nsnp,or_results$method,or_results$b,or_results$se,or_results$pval,or_results$or,or_results$or_lci95,or_results$or_uci95)
write.table(results, "hba1c_hnc_nooutliers.txt")

#MR presso 
devtools::install_github("rondolab/MR-PRESSO")
library(MRPRESSO)
mr_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat, NbDistribution = 1000,  SignifThreshold = 0.05)
mr_presso
write.table(mr_presso, "hba1c_mr_presso.txt")

#Calculate I-squared, r-squared and F-statistics
#I-squared function
Isq <- function(y,s){
  k          = length(y)
  w          = 1/s^2; sum.w  = sum(w)
  mu.hat     = sum(y*w)/sum.w  
  Q          = sum(w*(y-mu.hat)^2)
  Isq        = (Q - (k-1))/Q
  Isq        = max(0,Isq)
  return(Isq)
}

#Calculate Isq weighted and unweighted
I2<-c()
dat <- harmonise_data(exposure_dat, outcome_dat, action = 1)
str(dat)

#F-statistic
dat$samplesize.exposure <- 159940
dat$samplesize.outcome <- 6034
dat <- steiger_filtering(dat) 

N = dat$samplesize.exposure[1] #sample size
K = length(dat$SNP) #number of SNPs
total_r2 <- sum(dat$rsq.exposure) 
Fstat <- (N-K-1)/K * total_r2 / (1-total_r2)
total_r2
Fstat

#Rename required columns
dat $BetaXG<-dat $beta.exposure
dat $seBetaXG<-dat $se.exposure
BetaXG   = dat $BetaXG
seBetaXG = dat $seBetaXG 
seBetaYG<-dat $se.outcome
BXG = abs(BetaXG)        

# Calculate F statistics
# and I-squared statistics
# to measure Instrument 
# strength for MR-Egger

F = BXG^2/seBetaXG^2
mF = mean(F)
Isq_unweighted <- Isq(BXG,seBetaXG) #unweighted
Isq_weighted <- Isq((BXG/seBetaYG),(seBetaXG/seBetaYG)) #weighted

#Save results
output<-cbind (F, mF, Isq_unweighted, Isq_weighted)
I2<-rbind(I2, output)
colnames(I2) <- c("F", "mF", "Isq_unweighted", "Isq_weighted")
write.csv(I2, file="regression_dilution_isq_weighted_hba1c.csv", row.names = FALSE)

# OPC and OC subsite analysis

exposure_dat <- read_exposure_data(
  "exposure_dat_hba1c.txt",
  clump = TRUE,
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  min_pval = 1e-200,
  log_pval = FALSE
)

outcome_dat <- read_outcome_data("hba1c_oc.txt", sep="\t")
dat <- harmonise_data(exposure_dat, outcome_dat)
write.csv(dat, "hba1c_oc_harmonised.csv")    
mr <- mr(dat)
or_results <- generate_odds_ratios(mr)
or_results
write.table(or_results, "hba1c_oc_results.txt")

outcome_dat <- read_outcome_data("hba1c_opc.txt", sep="\t")
dat <- harmonise_data(exposure_dat, outcome_dat)
write.csv(dat, "hba1c_opc_harmonised.csv")    
mr <- mr(dat)
or_results <- generate_odds_ratios(mr)
or_results
write.table(or_results, "hba1c_opc_results.txt")



# Fasting glucose (FG) on oral and oropharyngeal cancer

exposure_dat <- read_exposure_data(
  "exposure_dat_fg.txt",
  clump = TRUE,
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  min_pval = 1e-200,
  log_pval = FALSE
)
write.csv(exposure_dat,'clumped_fg_exposure.csv')

#Extract SNPs from Lagou et al. from dbGaP or derive from IEU Open GWAS
outcome_dat <- read_outcome_data("fg_hnc.txt", sep="\t")
dat <- harmonise_data(exposure_dat, outcome_dat)
dat$outcome <- "Oral and oropharyngeal cancer"
dat$exposure <- "Fasting glucose"
write.csv(dat, "fg_hnc_harmonised.csv")    
mr_results <- mr(dat)
mr_results
or_results <- generate_odds_ratios(mr_results)
or_results
mr_report(dat, output_path = "~OneDrive - University of Bristol/Documents/Documents/GW4-CAT PhD/Research and publications/Metabolic_Phenotypes_MR/whr_mr",author="Mark Gormley", study = paste("Fasting glucose", "Oral and oropharyngeal cancer",sep=""))
results<-cbind.data.frame(or_results$outcome,or_results$nsnp,or_results$method,or_results$b,or_results$se,or_results$pval,or_results$or,or_results$or_lci95,or_results$or_uci95)
write.table(results, "fg_hnc_results.txt")

#Run scatter plot 
mr_scatter <- mr_scatter_plot(mr_results, dat)
mr_scatter
pdf("fg_hnc_scatter.pdf",  width = 15, height = 20)
ggarrange(mr_scatter$`GPfZzQ.8LpKoD`,ncol=2, nrow=2, widths = 2, heights = 1)
dev.off()

#6. PLOT FOREST AND LEAVE-ONE-OUT
res_single <- mr_singlesnp(dat)
mr_forest <- mr_forest_plot(res_single)
mr_forest

pdf("fg_hnc_forest.pdf",  width = 15, height = 20)
ggarrange(mr_forest$`GPfZzQ.8LpKoD`,ncol=2, nrow=2, widths = 2, heights = 1)
dev.off()

#Leave one out analysis
res_loo <- mr_leaveoneout(dat)
mr_loo <- mr_leaveoneout_plot(res_loo)
mr_loo

pdf("fg_hnc_loo.pdf",  width = 15, height = 20)
ggarrange(mr_loo$`GPfZzQ.8LpKoD`, nrow=2, widths = 2, heights = 1)
dev.off()

#Heterogeneity and pleiotropy analysis 
Het<-mr_heterogeneity(dat)
Plt<-mr_pleiotropy_test(dat)
Sin<-mr_singlesnp(dat)
write.table(Het, "fg_hnc_heterogeneity.txt")
write.table(Plt, "fg_hnc_pleiotropy.txt")
write.table(Sin, "fg_hnc_single_snp.txt")
p1 <- mr_scatter_plot(mr_results, dat)
p1[[1]]

res_single <- mr_singlesnp(dat)
p2 <- mr_forest_plot(res_single)
p2[[1]]

res_loo <- mr_leaveoneout(dat)
p3 <- mr_leaveoneout_plot(res_loo)
p3[[1]]

p4 <- mr_funnel_plot(res_single)
p4[[1]]

#Assess outliers
#Radial plots 
devtools::install_github("WSpiller/RadialMR")
library(RadialMR)
library(MendelianRandomization)

dat <- dat[dat$SNP%in%res_single$SNP,]

raddat <- format_radial(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome, dat$SNP)
ivwrad <- ivw_radial(raddat, alpha=0.05/63, weights=3)
dim(ivwrad$outliers)[1]
ivwrad <- ivw_radial(raddat, alpha=0.05, weights=3)
dim(ivwrad$outliers)[1]
eggrad <- egger_radial(raddat, alpha=0.05, weights=3)
eggrad$coef 
dim(eggrad$outliers)[1]

#plot_radial(ivwrad, TRUE, FALSE, TRUE)
plot_radial(c(ivwrad,eggrad), TRUE, FALSE, TRUE)

ivwrad$qstatistic 
ivwrad$sortoutliers <- ivwrad$outliers[order(ivwrad$outliers$p.value),]
ivwrad$sortoutliers$Qsum <- cumsum(ivwrad$sortoutliers$Q_statistic)
ivwrad$sortoutliers$Qdif <- ivwrad$sortoutliers$Qsum - ivwrad$qstatistic
write.csv(ivwrad$sortoutliers, "fg_outliers.csv", row.names=F, quote=F)

#Remove top outliers
dat2 <- dat[!dat$SNP %in% ivwrad$outliers$SNP,]
mr_results2 <- mr(dat2)
or_results <- generate_odds_ratios(mr_results2)

results<-cbind.data.frame(or_results$outcome,or_results$nsnp,or_results$method,or_results$b,or_results$se,or_results$pval,or_results$or,or_results$or_lci95,or_results$or_uci95)
write.table(results, "fg_hnc_nooutliers.txt")

#MR presso 
devtools::install_github("rondolab/MR-PRESSO")
library(MRPRESSO)
mr_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat, NbDistribution = 1000,  SignifThreshold = 0.05)
mr_presso
write.table(mr_presso, "fg_mr_presso.txt")

#Calculate I-squared, r-squared and F-statistics
#I-squared function
Isq <- function(y,s){
  k          = length(y)
  w          = 1/s^2; sum.w  = sum(w)
  mu.hat     = sum(y*w)/sum.w  
  Q          = sum(w*(y-mu.hat)^2)
  Isq        = (Q - (k-1))/Q
  Isq        = max(0,Isq)
  return(Isq)
}

#Calculate Isq weighted and unweighted
I2<-c()
dat <- harmonise_data(exposure_dat, outcome_dat, action = 1)
str(dat)

#F-statistic
dat$samplesize.exposure <- 151188
dat$samplesize.outcome <- 6034
dat <- steiger_filtering(dat) 

N = dat$samplesize.exposure[1] #sample size
K = length(dat$SNP) #number of SNPs
total_r2 <- sum(dat$rsq.exposure) 
Fstat <- (N-K-1)/K * total_r2 / (1-total_r2)
total_r2
Fstat

#Rename required columns
dat $BetaXG<-dat $beta.exposure
dat $seBetaXG<-dat $se.exposure
BetaXG   = dat $BetaXG
seBetaXG = dat $seBetaXG 
seBetaYG<-dat $se.outcome
BXG = abs(BetaXG)        

# Calculate F statistics
# and I-squared statistics
# to measure Instrument 
# strength for MR-Egger

F = BXG^2/seBetaXG^2
mF = mean(F)
Isq_unweighted <- Isq(BXG,seBetaXG) #unweighted
Isq_weighted <- Isq((BXG/seBetaYG),(seBetaXG/seBetaYG)) #weighted

#Save results
output<-cbind (F, mF, Isq_unweighted, Isq_weighted)
I2<-rbind(I2, output)
colnames(I2) <- c("F", "mF", "Isq_unweighted", "Isq_weighted")
write.csv(I2, file="regression_dilution_isq_weighted_fg.csv", row.names = FALSE)

# OPC and OC subsite analysis

exposure_dat <- read_exposure_data(
  "exposure_dat_fg.txt",
  clump = TRUE,
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  min_pval = 1e-200,
  log_pval = FALSE
)

outcome_dat <- read_outcome_data("fg_oc.txt", sep="\t")
dat <- harmonise_data(exposure_dat, outcome_dat)
write.csv(dat, "fg_oc_harmonised.csv")    
mr <- mr(dat)
or_results <- generate_odds_ratios(mr)
or_results
write.table(or_results, "fg_oc_results.txt")

outcome_dat <- read_outcome_data("fg_opc.txt", sep="\t")
dat <- harmonise_data(exposure_dat, outcome_dat)
write.csv(dat, "fg_opc_harmonised.csv")    
mr <- mr(dat)
or_results <- generate_odds_ratios(mr)
or_results
write.table(or_results, "fg_opc_results.txt")


#Fasting insulin (FI) on oral and orophrayngeal cancer

exposure_dat <- read_exposure_data(
  "exposure_dat_fi.txt",
  clump = TRUE,
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  min_pval = 1e-200,
  log_pval = FALSE
)
write.csv(exposure_dat,'clumped_fi_exposure.csv')


#Extract SNPs for FI from Lagou et al. from dbGaP or derive from IEU Open GWAS
outcome_dat <- read_outcome_data("fi_hnc.txt", sep="\t")
dat <- harmonise_data(exposure_dat, outcome_dat)
dat$outcome <- "Oral and oropharyngeal cancer"
dat$exposure <- "Fasting insulin"
write.csv(dat, "fi_hnc_harmonised.csv")    
mr_results <- mr(dat)
mr_results
or_results <- generate_odds_ratios(mr_results)
or_results
mr_report(dat, output_path = "~OneDrive - University of Bristol/Documents/Documents/GW4-CAT PhD/Research and publications/Metabolic_Phenotypes_MR/fi_mr",author="Mark Gormley", study = paste("Fasting insulin", "Oral and oropharyngeal cancer",sep=""))
results<-cbind.data.frame(or_results$outcome,or_results$nsnp,or_results$method,or_results$b,or_results$se,or_results$pval,or_results$or,or_results$or_lci95,or_results$or_uci95)
write.table(results, "fi_hnc_results.txt")

#Run scatter plot
mr_scatter <- mr_scatter_plot(mr_results, dat)
mr_scatter

pdf("fi_hnc_scatter.pdf",  width = 15, height = 20)
ggarrange(mr_scatter$`DiU28S.auPXn3`,ncol=2, nrow=2, widths = 2, heights = 1)
dev.off()

#6. PLOT FOREST AND LEAVE-ONE-OUT
res_single <- mr_singlesnp(dat)
mr_forest <- mr_forest_plot(res_single)
mr_forest

pdf("fi_hnc_forest.pdf",  width = 15, height = 20)
ggarrange(mr_forest$`DiU28S.auPXn3`,ncol=2, nrow=2, widths = 2, heights = 1)
dev.off()

#Leave one out analysis
res_loo <- mr_leaveoneout(dat)
mr_loo <- mr_leaveoneout_plot(res_loo)
mr_loo

pdf("fi_hnc_loo.pdf",  width = 15, height = 20)
ggarrange(mr_loo$`DiU28S.auPXn3`, nrow=2, widths = 2, heights = 1)
dev.off()

#Heterogeneity and pleiotropy analysis 
Het<-mr_heterogeneity(dat)
Plt<-mr_pleiotropy_test(dat)
Sin<-mr_singlesnp(dat)
write.table(Het, "fi_hnc_heterogeneity.txt")
write.table(Plt, "fi_hnc_pleiotropy.txt")
write.table(Sin, "fi_hnc_single_snp.txt")
p1 <- mr_scatter_plot(mr_results, dat)
p1[[1]]

res_single <- mr_singlesnp(dat)
p2 <- mr_forest_plot(res_single)
p2[[1]]

res_loo <- mr_leaveoneout(dat)
p3 <- mr_leaveoneout_plot(res_loo)
p3[[1]]

p4 <- mr_funnel_plot(res_single)
p4[[1]]

#Assess outliers
#Radial plots 
devtools::install_github("WSpiller/RadialMR")
library(RadialMR)
library(MendelianRandomization)

dat <- dat[dat$SNP%in%res_single$SNP,]

raddat <- format_radial(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome, dat$SNP)
ivwrad <- ivw_radial(raddat, alpha=0.05/63, weights=3)
dim(ivwrad$outliers)[1]
ivwrad <- ivw_radial(raddat, alpha=0.05, weights=3)
dim(ivwrad$outliers)[1]
eggrad <- egger_radial(raddat, alpha=0.05, weights=3)
eggrad$coef 
dim(eggrad$outliers)[1]

#plot_radial(ivwrad, TRUE, FALSE, TRUE)
plot_radial(c(ivwrad,eggrad), TRUE, FALSE, TRUE)

ivwrad$qstatistic 
ivwrad$sortoutliers <- ivwrad$outliers[order(ivwrad$outliers$p.value),]
ivwrad$sortoutliers$Qsum <- cumsum(ivwrad$sortoutliers$Q_statistic)
ivwrad$sortoutliers$Qdif <- ivwrad$sortoutliers$Qsum - ivwrad$qstatistic
write.csv(ivwrad$sortoutliers, "fi_outliers.csv", row.names=F, quote=F)

#Remove top outliers
dat2 <- dat[!dat$SNP %in% ivwrad$outliers$SNP,]
mr_results2 <- mr(dat2)
or_results <- generate_odds_ratios(mr_results2)

results<-cbind.data.frame(or_results$outcome,or_results$nsnp,or_results$method,or_results$b,or_results$se,or_results$pval,or_results$or,or_results$or_lci95,or_results$or_uci95)
write.table(results, "fi_hnc_nooutliers.txt")

#MR presso 
devtools::install_github("rondolab/MR-PRESSO")
library(MRPRESSO)
mr_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat, NbDistribution = 1000,  SignifThreshold = 0.05)
mr_presso
write.table(mr_presso, "fi_mr_presso.txt")

#Calculate I-squared, r-squared and F-statistics
#I-squared function
Isq <- function(y,s){
  k          = length(y)
  w          = 1/s^2; sum.w  = sum(w)
  mu.hat     = sum(y*w)/sum.w  
  Q          = sum(w*(y-mu.hat)^2)
  Isq        = (Q - (k-1))/Q
  Isq        = max(0,Isq)
  return(Isq)
}

#Calculate Isq weighted and unweighted
I2<-c()
dat <- harmonise_data(exposure_dat, outcome_dat, action = 1)
str(dat)

#F-statistic
dat$samplesize.exposure <- 105056
dat$samplesize.outcome <- 6034
dat <- steiger_filtering(dat) 

N = dat$samplesize.exposure[1] #sample size
K = length(dat$SNP) #number of SNPs
total_r2 <- sum(dat$rsq.exposure) 
Fstat <- (N-K-1)/K * total_r2 / (1-total_r2)
total_r2
Fstat

#Rename required columns
dat $BetaXG<-dat $beta.exposure
dat $seBetaXG<-dat $se.exposure
BetaXG   = dat $BetaXG
seBetaXG = dat $seBetaXG 
seBetaYG<-dat $se.outcome
BXG = abs(BetaXG)        

# Calculate F statistics
# and I-squared statistics
# to measure Instrument 
# strength for MR-Egger

F = BXG^2/seBetaXG^2
mF = mean(F)
Isq_unweighted <- Isq(BXG,seBetaXG) #unweighted
Isq_weighted <- Isq((BXG/seBetaYG),(seBetaXG/seBetaYG)) #weighted

#Save results
output<-cbind (F, mF, Isq_unweighted, Isq_weighted)
I2<-rbind(I2, output)
colnames(I2) <- c("F", "mF", "Isq_unweighted", "Isq_weighted")
write.csv(I2, file="regression_dilution_isq_weighted_fi.csv", row.names = FALSE)

#SIMEX correction
#Create empty dataframe to store output
simexegger<-c()

#Run SIMEX
#Rename required columns
dat$BetaXG<-dat$beta.exposure
dat$seBetaXG<-dat$se.exposure
dat$BetaYG<-dat$beta.outcome
dat$seBetaYG<-dat$se.outcome
BetaXG <- dat$BetaXG
BetaYG <- dat$BetaYG
seBetaXG <- dat$seBetaXG
seBetaYG <- dat$seBetaYG

BYG <- BetaYG*sign(BetaXG)#Pre-processing steps to ensure all gene--exposure estimates are positive
BXG <- abs(BetaXG)         

#MR-Egger regression (weighted) 
Fit1 <- lm(BYG ~ BXG,weights=1/seBetaYG^2,x=TRUE,y=TRUE)

#MR-Egger regression (unweighted)
Fit2 <- lm(BYG~BXG,x=TRUE,y=TRUE) 

#Simulation extrapolation 
mod.sim1 <- simex(Fit1,B=1000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE") 
mod.sim2 <- simex(Fit2,B=1000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE") 
mod1<-summary(mod.sim1)
mod2<-summary(mod.sim2)
#Print results
mod1
mod2
#Use BXG result and exponentiate the estimates
OR = exp(-4.65039) 
CIL = exp(-4.65039 - 1.96 * 3.12237) 
CIU = exp(-4.65039 + 1.96 * 3.12237) 

# OPC and OC subsite analysis

exposure_dat <- read_exposure_data(
  "exposure_dat_fi.txt",
  clump = TRUE,
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  min_pval = 1e-200,
  log_pval = FALSE
)

outcome_dat <- read_outcome_data("fi_oc.txt", sep="\t")
dat <- harmonise_data(exposure_dat, outcome_dat)
write.csv(dat, "fi_oc_harmonised.csv")    
mr <- mr(dat)
or_results <- generate_odds_ratios(mr)
or_results
write.table(or_results, "fi_oc_results.txt")

outcome_dat <- read_outcome_data("fi_opc.txt", sep="\t")
dat <- harmonise_data(exposure_dat, outcome_dat)
write.csv(dat, "fi_opc_harmonised.csv")    
mr <- mr(dat)
or_results <- generate_odds_ratios(mr)
or_results
write.table(or_results, "fi_opc_results.txt")


# Systolic blood pressure (SBP) on oral and oropharyngeal cancer
exposure_dat <- read_exposure_data(
  "exposure_dat_sbp.csv",
  clump = TRUE,
  sep = ",",
  phenotype_col = "Exposure",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "effect allele",
  other_allele_col = "other allele",
  pval_col = "p-value",
  min_pval = 1e-200,
  log_pval = FALSE
)
write.csv(exposure_dat,'clumped_sbp_exposure.csv')

#Extract SNPs for BMI from Evangelou et al. from dbGaP or derive from IEU Open GWAS
outcome_dat <- read_outcome_data("sbp_hnc.txt", sep="\t")
dat <- harmonise_data(exposure_dat, outcome_dat)
dat$outcome <- "Oral and oropharyngeal cancer"
dat$exposure <- "Systolic blood pressure"
write.csv(dat, "sbp_hnc_harmonised.csv")    
mr_results <- mr(dat)
mr_results
or_results <- generate_odds_ratios(mr_results)
or_results
mr_report(dat, output_path = "~/OneDrive - University of Bristol/Documents/Documents/GW4-CAT PhD/Research and publications/Metabolic_Phenotypes_MR/sbp_mr",author="Mark Gormley", study = paste("Systolic blood pressure", "Oral and oropharyngeal cancer",sep=""))
results<-cbind.data.frame(or_results$outcome,or_results$nsnp,or_results$method,or_results$b,or_results$se,or_results$pval,or_results$or,or_results$or_lci95,or_results$or_uci95)
write.table(results, "sbp_hnc_results.csv")

#Run scatter plot
mr_scatter <- mr_scatter_plot(mr_results, dat)
mr_scatter
pdf("sbp_hnc_scatter.pdf",  width = 15, height = 20)
ggarrange(mr_scatter$`ddN0Vx.7oHMil`,ncol=2, nrow=2, widths = 2, heights = 1)
dev.off()

#6. PLOT FOREST AND LEAVE-ONE-OUT
res_single <- mr_singlesnp(dat)
mr_forest <- mr_forest_plot(res_single)
mr_forest

pdf("sbp_hnc_forest.pdf",  width = 15, height = 20)
ggarrange(mr_forest$`ddN0Vx.7oHMil`,ncol=2, nrow=2, widths = 2, heights = 1)
dev.off()

#Leave one out analysis
res_loo <- mr_leaveoneout(dat)
mr_loo <- mr_leaveoneout_plot(res_loo)
mr_loo

pdf("sbp_hnc_loo.pdf",  width = 15, height = 20)
ggarrange(mr_loo$`ddN0Vx.7oHMil`, nrow=2, widths = 2, heights = 1)
dev.off()

#Heterogeneity and pleiotropy analysis 
Het<-mr_heterogeneity(dat)
Plt<-mr_pleiotropy_test(dat)
Sin<-mr_singlesnp(dat)
write.table(Het, "sbp_hnc_heterogeneity.txt")
write.table(Plt, "sbp_hnc_pleiotropy.txt")
write.table(Sin, "sbp_hnc_single_snp.txt")
p1 <- mr_scatter_plot(mr_results, dat)
p1[[1]]

res_single <- mr_singlesnp(dat)
p2 <- mr_forest_plot(res_single)
p2[[1]]

res_loo <- mr_leaveoneout(dat)
p3 <- mr_leaveoneout_plot(res_loo)
p3[[1]]

p4 <- mr_funnel_plot(res_single)
p4[[1]]

#Assess outliers
#Radial plots 
devtools::install_github("WSpiller/RadialMR")
library(RadialMR)
library(MendelianRandomization)

dat <- dat[dat$SNP%in%res_single$SNP,]

raddat <- format_radial(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome, dat$SNP)
ivwrad <- ivw_radial(raddat, alpha=0.05/63, weights=3)
dim(ivwrad$outliers)[1]
ivwrad <- ivw_radial(raddat, alpha=0.05, weights=3)
dim(ivwrad$outliers)[1]
eggrad <- egger_radial(raddat, alpha=0.05, weights=3)
eggrad$coef 
dim(eggrad$outliers)[1]

#plot_radial(ivwrad, TRUE, FALSE, TRUE)
plot_radial(c(ivwrad,eggrad), TRUE, FALSE, TRUE)

ivwrad$qstatistic 
ivwrad$sortoutliers <- ivwrad$outliers[order(ivwrad$outliers$p.value),]
ivwrad$sortoutliers$Qsum <- cumsum(ivwrad$sortoutliers$Q_statistic)
ivwrad$sortoutliers$Qdif <- ivwrad$sortoutliers$Qsum - ivwrad$qstatistic
write.csv(ivwrad$sortoutliers, "sbp_outliers.csv", row.names=F, quote=F)

#Remove top outliers
dat2 <- dat[!dat$SNP %in% ivwrad$outliers$SNP,]
mr_results2 <- mr(dat2)
or_results <- generate_odds_ratios(mr_results2)

results<-cbind.data.frame(or_results$outcome,or_results$nsnp,or_results$method,or_results$b,or_results$se,or_results$pval,or_results$or,or_results$or_lci95,or_results$or_uci95)
write.table(results, "sbp_hnc_nooutliers.txt")

#MR presso 
devtools::install_github("rondolab/MR-PRESSO")
library(MRPRESSO)
mr_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat, NbDistribution = 1000,  SignifThreshold = 0.05)
mr_presso
write.table(mr_presso, "sbp_mr_presso.csv")

#Calculate I-squared, r-squared and F-statistics
#I-squared function
Isq <- function(y,s){
  k          = length(y)
  w          = 1/s^2; sum.w  = sum(w)
  mu.hat     = sum(y*w)/sum.w  
  Q          = sum(w*(y-mu.hat)^2)
  Isq        = (Q - (k-1))/Q
  Isq        = max(0,Isq)
  return(Isq)
}

#Calculate Isq weighted and unweighted
I2<-c()
dat <- harmonise_data(exposure_dat, outcome_dat, action = 1)
str(dat)

#F-statistic
dat$samplesize.exposure <- 757601
dat$samplesize.outcome <- 6034
dat <- steiger_filtering(dat) 

N = dat$samplesize.exposure[1] #sample size
K = length(dat$SNP) #number of SNPs
total_r2 <- sum(dat$rsq.exposure) 
Fstat <- (N-K-1)/K * total_r2 / (1-total_r2)
total_r2
Fstat

#Rename required columns
dat $BetaXG<-dat $beta.exposure
dat $seBetaXG<-dat $se.exposure
BetaXG   = dat $BetaXG
seBetaXG = dat $seBetaXG 
seBetaYG<-dat $se.outcome

BXG = abs(BetaXG)         
# Calculate F statistics
# and I-squared statistics
# to measure Instrument 
# strength for MR-Egger

F = BXG^2/seBetaXG^2
mF = mean(F)
Isq_unweighted <- Isq(BXG,seBetaXG) #unweighted
Isq_weighted <- Isq((BXG/seBetaYG),(seBetaXG/seBetaYG)) #weighted

#Save results
output<-cbind (F, mF, Isq_unweighted, Isq_weighted)
I2<-rbind(I2, output)
colnames(I2) <- c("F", "mF", "Isq_unweighted", "Isq_weighted")
write.csv(I2, file="regression_dilution_isq_weighted_sbp.csv", row.names = FALSE)

#Run SIMEX
#Rename required columns
dat$BetaXG<-dat$beta.exposure
dat$seBetaXG<-dat$se.exposure
dat$BetaYG<-dat$beta.outcome
dat$seBetaYG<-dat$se.outcome
BetaXG <- dat$BetaXG
BetaYG <- dat$BetaYG
seBetaXG <- dat$seBetaXG
seBetaYG <- dat$seBetaYG

BYG <- BetaYG*sign(BetaXG)#Pre-processing steps to ensure all gene--exposure estimates are positive
BXG <- abs(BetaXG)         

#MR-Egger regression (weighted) 
Fit1 <- lm(BYG ~ BXG,weights=1/seBetaYG^2,x=TRUE,y=TRUE)

#MR-Egger regression (unweighted)
Fit2 <- lm(BYG~BXG,x=TRUE,y=TRUE) 

#Simulation extrapolation 
mod.sim1 <- simex(Fit1,B=1000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE") 
mod.sim2 <- simex(Fit2,B=1000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE") 
mod1<-summary(mod.sim1)
mod2<-summary(mod.sim2)
#Print results
mod1
mod2
#Use BXG result and exponentiate the estimates
OR = exp(0.13994) 
CIL = exp(0.13994 - 1.96 * 0.04788) 
CIU = exp(0.13994 + 1.96 * 0.04788) 

# OPC and OC subsite analysis
exposure_dat <- read_exposure_data(
  "exposure_dat_sbp.csv",
  clump = TRUE,
  sep = ",",
  phenotype_col = "Exposure",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "effect allele",
  other_allele_col = "other allele",
  pval_col = "p-value",
  min_pval = 1e-200,
  log_pval = FALSE
)

outcome_dat <- read_outcome_data("sbp_oc.txt", sep="\t")
dat <- harmonise_data(exposure_dat, outcome_dat)
write.csv(dat, "sbp_oc_harmonised.csv")    
mr <- mr(dat)
or_results <- generate_odds_ratios(mr)
or_results
write.table(or_results, "sbp_oc_results.txt")

outcome_dat <- read_outcome_data("sbp_opc.txt", sep="\t")
dat <- harmonise_data(exposure_dat, outcome_dat)
write.csv(dat, "sbp_opc_harmonised.csv")    
mr <- mr(dat)
or_results <- generate_odds_ratios(mr)
or_results
write.table(or_results, "sbp_opc_results.txt")



# Diastolic blood pressure on oral and oropharyngeal cancer

exposure_dat <- read_exposure_data(
  "exposure_dat_dbp.csv",
  clump = TRUE,
  sep = ",",
  phenotype_col = "Exposure",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "effect allele",
  other_allele_col = "other allele",
  pval_col = "p-value",
  min_pval = 1e-200,
  log_pval = FALSE
)
write.csv(exposure_dat,'clumped_dbp_exposure.csv')

#Extract SNPs for BMI from Evangelou et al. from dbGaP or derive from IEU Open GWAS

outcome_dat <- read_outcome_data("dbp_hnc.txt", sep="\t")
dat <- harmonise_data(exposure_dat, outcome_dat)
dat$outcome <- "Oral and oropharyngeal cancer"
dat$exposure <- "Diastolic blood pressure"
write.csv(dat, "dbp_hnc_harmonised.csv")    
mr_results <- mr(dat)
mr_results
or_results <- generate_odds_ratios(mr_results)
or_results
mr_report(dat, output_path = "~/OneDrive - University of Bristol/Documents/Documents/GW4-CAT PhD/Research and publications/Metabolic_Phenotypes_MR/dbp_mr",author="Mark Gormley", study = paste("Diastolic blood pressure", "Oral and oropharyngeal cancer",sep=""))
results<-cbind.data.frame(or_results$outcome,or_results$nsnp,or_results$method,or_results$b,or_results$se,or_results$pval,or_results$or,or_results$or_lci95,or_results$or_uci95)
write.table(results, "dbp_hnc_results.txt")

#Run scatter plot
mr_scatter <- mr_scatter_plot(mr_results, dat)
mr_scatter
pdf("dbp_hnc_scatter.pdf",  width = 15, height = 20)
ggarrange(mr_scatter$`wWJ9EI.NcEcx5`,ncol=2, nrow=2, widths = 2, heights = 1)
dev.off()

#6. PLOT FOREST AND LEAVE-ONE-OUT
res_single <- mr_singlesnp(dat)
mr_forest <- mr_forest_plot(res_single)
mr_forest

pdf("dbp_hnc_forest.pdf",  width = 15, height = 20)
ggarrange(mr_forest$`wWJ9EI.NcEcx5`,ncol=2, nrow=2, widths = 2, heights = 1)
dev.off()

#Leave one out analysis
res_loo <- mr_leaveoneout(dat)
mr_loo <- mr_leaveoneout_plot(res_loo)
mr_loo

pdf("dbp_hnc_loo.pdf",  width = 15, height = 20)
ggarrange(mr_loo$`wWJ9EI.NcEcx5`, nrow=2, widths = 2, heights = 1)
dev.off()

#Heterogeneity and pleiotropy analysis 
Het<-mr_heterogeneity(dat)
Plt<-mr_pleiotropy_test(dat)
Sin<-mr_singlesnp(dat)
write.table(Het, "dbp_hnc_heterogeneity.txt")
write.table(Plt, "dbp_hnc_pleiotropy.txt")
write.table(Sin, "dbp_hnc_single_snp.txt")
p1 <- mr_scatter_plot(mr_results, dat)
p1[[1]]

res_single <- mr_singlesnp(dat)
p2 <- mr_forest_plot(res_single)
p2[[1]]

res_loo <- mr_leaveoneout(dat)
p3 <- mr_leaveoneout_plot(res_loo)
p3[[1]]

p4 <- mr_funnel_plot(res_single)
p4[[1]]

#Assess outliers
#Radial plots 
devtools::install_github("WSpiller/RadialMR")
library(RadialMR)
library(MendelianRandomization)

dat <- dat[dat$SNP%in%res_single$SNP,]

raddat <- format_radial(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome, dat$SNP)
ivwrad <- ivw_radial(raddat, alpha=0.05/63, weights=3)
dim(ivwrad$outliers)[1] 
#0 outliers at bonf 
ivwrad <- ivw_radial(raddat, alpha=0.05, weights=3)
dim(ivwrad$outliers)[1] 
#7 outliers at 0.05

eggrad <- egger_radial(raddat, alpha=0.05, weights=3)
eggrad$coef 
dim(eggrad$outliers)[1] 
#6 outliers at 0.05 

#plot_radial(ivwrad, TRUE, FALSE, TRUE)
plot_radial(c(ivwrad,eggrad), TRUE, FALSE, TRUE)

ivwrad$qstatistic 
ivwrad$sortoutliers <- ivwrad$outliers[order(ivwrad$outliers$p.value),]
ivwrad$sortoutliers$Qsum <- cumsum(ivwrad$sortoutliers$Q_statistic)
ivwrad$sortoutliers$Qdif <- ivwrad$sortoutliers$Qsum - ivwrad$qstatistic
write.csv(ivwrad$sortoutliers, "dbp_outliers.csv", row.names=F, quote=F)

#Remove top outliers
dat2 <- dat[!dat$SNP %in% ivwrad$outliers$SNP,]
mr_results2 <- mr(dat2)
or_results <- generate_odds_ratios(mr_results2)

results<-cbind.data.frame(or_results$outcome,or_results$nsnp,or_results$method,or_results$b,or_results$se,or_results$pval,or_results$or,or_results$or_lci95,or_results$or_uci95)
write.table(results, "dbp_hnc_nooutliers.txt")

#MR presso 
devtools::install_github("rondolab/MR-PRESSO")
library(MRPRESSO)
mr_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat, NbDistribution = 1000,  SignifThreshold = 0.05)
mr_presso
write.table(mr_presso, "dbp_mr_presso.txt")

#Calculate I-squared, r-squared and F-statistics
#I-squared function
Isq <- function(y,s){
  k          = length(y)
  w          = 1/s^2; sum.w  = sum(w)
  mu.hat     = sum(y*w)/sum.w  
  Q          = sum(w*(y-mu.hat)^2)
  Isq        = (Q - (k-1))/Q
  Isq        = max(0,Isq)
  return(Isq)
}

#Calculate Isq weighted and unweighted
I2<-c()
dat <- harmonise_data(exposure_dat, outcome_dat, action = 1)
str(dat)

#F-statistic
dat$samplesize.exposure <- 681275
dat$samplesize.outcome <- 6034
dat <- steiger_filtering(dat) 

N = dat$samplesize.exposure[1] #sample size
K = length(dat$SNP) #number of SNPs
total_r2 <- sum(dat$rsq.exposure) 
Fstat <- (N-K-1)/K * total_r2 / (1-total_r2)
total_r2
Fstat

#Rename required columns
dat $BetaXG<-dat $beta.exposure
dat $seBetaXG<-dat $se.exposure
BetaXG   = dat $BetaXG
seBetaXG = dat $seBetaXG 
seBetaYG<-dat $se.outcome

BXG = abs(BetaXG)         
# Calculate F statistics
# and I-squared statistics
# to measure Instrument 
# strength for MR-Egger

F = BXG^2/seBetaXG^2
mF = mean(F)
Isq_unweighted <- Isq(BXG,seBetaXG) #unweighted
Isq_weighted <- Isq((BXG/seBetaYG),(seBetaXG/seBetaYG)) #weighted

#Save results
output<-cbind (F, mF, Isq_unweighted, Isq_weighted)
I2<-rbind(I2, output)
colnames(I2) <- c("F", "mF", "Isq_unweighted", "Isq_weighted")
write.csv(I2, file="regression_dilution_isq_weighted_dbp.csv", row.names = FALSE)

#Run SIMEX
#Rename required columns
dat$BetaXG<-dat$beta.exposure
dat$seBetaXG<-dat$se.exposure
dat$BetaYG<-dat$beta.outcome
dat$seBetaYG<-dat$se.outcome
BetaXG <- dat$BetaXG
BetaYG <- dat$BetaYG
seBetaXG <- dat$seBetaXG
seBetaYG <- dat$seBetaYG

BYG <- BetaYG*sign(BetaXG)#Pre-processing steps to ensure all gene--exposure estimates are positive
BXG <- abs(BetaXG)         

#MR-Egger regression (weighted) 
Fit1 <- lm(BYG ~ BXG,weights=1/seBetaYG^2,x=TRUE,y=TRUE)

#MR-Egger regression (unweighted)
Fit2 <- lm(BYG~BXG,x=TRUE,y=TRUE) 

#Simulation extrapolation 
mod.sim1 <- simex(Fit1,B=1000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE") 
mod.sim2 <- simex(Fit2,B=1000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE") 
mod1<-summary(mod.sim1)
mod2<-summary(mod.sim2)
#Print results
mod1
mod2
#Use BXG result and exponentiate the estimates
OR = exp(0.06413) 
CIL = exp(0.06413 - 1.96 * 0.08080) 
CIU = exp(0.06413 + 1.96 * 0.08080) 

# OPC and OC subsite analysis
exposure_dat <- read_exposure_data(
  "exposure_dat_dbp.csv",
  clump = TRUE,
  sep = ",",
  phenotype_col = "Exposure",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "effect allele",
  other_allele_col = "other allele",
  pval_col = "p-value",
  min_pval = 1e-200,
  log_pval = FALSE
)

outcome_dat <- read_outcome_data("dbp_oc.txt", sep="\t")
dat <- harmonise_data(exposure_dat, outcome_dat)
write.csv(dat, "dbp_oc_harmonised.csv")    
mr <- mr(dat)
or_results <- generate_odds_ratios(mr)
or_results
write.table(or_results, "dbp_oc_results.txt")

outcome_dat <- read_outcome_data("dbp_opc.txt", sep="\t")
dat <- harmonise_data(exposure_dat, outcome_dat)
write.csv(dat, "dbp_opc_harmonised.csv")    
mr <- mr(dat)
or_results <- generate_odds_ratios(mr)
or_results
write.table(or_results, "dbp_opc_results.txt")


# Instrument-risk factor analyses

# Insert exposure in place of X i.e., bmi, wc, whr, t2d, hba1c or dbp
exposure_dat <- read_exposure_data(
  "exposure_dat_X.csv",
  clump = TRUE,
  sep = ",",
  phenotype_col = "Exposure",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "p-value",
  min_pval = 1e-200,
  log_pval = FALSE
)
write.csv(exposure_dat,'clumped_X_exposure.csv')


# Smoking
# Extract SNPs for DBP from Liu et al. GWAS
outcome_dat <- read_outcome_data("X_smoking.txt", sep="\t")
dat <- harmonise_data(exposure_dat, outcome_dat)
dat$outcome <- "Smoking initiation"
dat$exposure <- "Exposure"
write.csv(dat, "X_smoking_harmonised.csv")    
mr_results <- mr(dat)
mr_results
or_results <- generate_odds_ratios(mr_results)
or_results
results<-cbind.data.frame(or_results$outcome,or_results$nsnp,or_results$method,or_results$b,or_results$se,or_results$pval,or_results$or,or_results$or_lci95,or_results$or_uci95)
write.table(results, "X_smoking_results.txt")

#Run scatter plot code
mr_scatter <- mr_scatter_plot(mr_results, dat)
mr_scatter
pdf("X_smoking_scatter.pdf",  width = 15, height = 20)
ggarrange(mr_scatter$`I8F99P.r7zXbr`,ncol=2, nrow=2, widths = 2, heights = 1)
dev.off()

#6. PLOT FOREST AND LEAVE-ONE-OUT
res_single <- mr_singlesnp(dat)
mr_forest <- mr_forest_plot(res_single)
mr_forest

pdf("X_smoking_forest.pdf",  width = 15, height = 20)
ggarrange(mr_forest$`I8F99P.r7zXbr`,ncol=2, nrow=2, widths = 2, heights = 1)
dev.off()

#Leave one out analysis
res_loo <- mr_leaveoneout(dat)
mr_loo <- mr_leaveoneout_plot(res_loo)
mr_loo

pdf("X_smoking_loo.pdf",  width = 15, height = 20)
ggarrange(mr_loo$`I8F99P.r7zXbr`, nrow=2, widths = 2, heights = 1)
dev.off()

#Heterogeneity and pleiotropy analysis 
Het<-mr_heterogeneity(dat)
Plt<-mr_pleiotropy_test(dat)
Sin<-mr_singlesnp(dat)
write.table(Het, "X_smoking_heterogeneity.txt")
write.table(Plt, "X_smoking_pleiotropy.txt")
write.table(Sin, "X_smoking_single_snp.txt")
p1 <- mr_scatter_plot(mr_results, dat)
p1[[1]]

res_single <- mr_singlesnp(dat)
p2 <- mr_forest_plot(res_single)
p2[[1]]

res_loo <- mr_leaveoneout(dat)
p3 <- mr_leaveoneout_plot(res_loo)
p3[[1]]

p4 <- mr_funnel_plot(res_single)
p4[[1]]

#Assess outliers
#Radial plots 
devtools::install_github("WSpiller/RadialMR")
library(RadialMR)
library(MendelianRandomization)

dat <- dat[dat$SNP%in%res_single$SNP,]

raddat <- format_radial(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome, dat$SNP)
ivwrad <- ivw_radial(raddat, alpha=0.05/63, weights=3)
dim(ivwrad$outliers)[1] 
ivwrad <- ivw_radial(raddat, alpha=0.05, weights=3)
dim(ivwrad$outliers)[1] 
eggrad <- egger_radial(raddat, alpha=0.05, weights=3)
eggrad$coef 
dim(eggrad$outliers)[1] 


#plot_radial(ivwrad, TRUE, FALSE, TRUE)
plot_radial(c(ivwrad,eggrad), TRUE, FALSE, TRUE)

ivwrad$qstatistic 
ivwrad$sortoutliers <- ivwrad$outliers[order(ivwrad$outliers$p.value),]
ivwrad$sortoutliers$Qsum <- cumsum(ivwrad$sortoutliers$Q_statistic)
ivwrad$sortoutliers$Qdif <- ivwrad$sortoutliers$Qsum - ivwrad$qstatistic
write.csv(ivwrad$sortoutliers, "X_smoking_outliers.csv", row.names=F, quote=F)

#Remove top outliers
dat2 <- dat[!dat$SNP %in% ivwrad$outliers$SNP,]
mr_results2 <- mr(dat2)
or_results <- generate_odds_ratios(mr_results2)

results<-cbind.data.frame(or_results$outcome,or_results$nsnp,or_results$method,or_results$b,or_results$se,or_results$pval,or_results$or,or_results$or_lci95,or_results$or_uci95)
write.table(results, "X_smoking_nooutliers.txt")


# Alcohol
# Extract SNPs for DBP from Liu et al. GWAS
outcome_dat <- read_outcome_data("X_alcohol.txt", sep="\t")
dat <- harmonise_data(exposure_dat, outcome_dat)
dat$outcome <- "Alcohol drinks per week"
dat$exposure <- "Exposure"
write.csv(dat, "X_alcohol_harmonised.csv")    
mr_results <- mr(dat)
mr_results
or_results <- generate_odds_ratios(mr_results)
or_results
results<-cbind.data.frame(or_results$outcome,or_results$nsnp,or_results$method,or_results$b,or_results$se,or_results$pval,or_results$or,or_results$or_lci95,or_results$or_uci95)
write.table(results, "X_alcohol_results.txt")

#Run scatter plot code 
mr_scatter <- mr_scatter_plot(mr_results, dat)
mr_scatter
pdf("X_alcohol_scatter.pdf",  width = 15, height = 20)
ggarrange(mr_scatter$`WrFbNW.KwTVNS`,ncol=2, nrow=2, widths = 2, heights = 1)
dev.off()

#6. PLOT FOREST AND LEAVE-ONE-OUT
res_single <- mr_singlesnp(dat)
mr_forest <- mr_forest_plot(res_single)
mr_forest

pdf("X_alcohol_forest.pdf",  width = 15, height = 20)
ggarrange(mr_forest$`WrFbNW.KwTVNS`,ncol=2, nrow=2, widths = 2, heights = 1)
dev.off()

#Leave one out analysis
res_loo <- mr_leaveoneout(dat)
mr_loo <- mr_leaveoneout_plot(res_loo)
mr_loo

pdf("X_alcohol_loo.pdf",  width = 15, height = 20)
ggarrange(mr_loo$`WrFbNW.KwTVNS`, nrow=2, widths = 2, heights = 1)
dev.off()

#Heterogeneity and pleiotropy analysis 
Het<-mr_heterogeneity(dat)
Plt<-mr_pleiotropy_test(dat)
Sin<-mr_singlesnp(dat)
write.table(Het, "X_alcohol_heterogeneity.txt")
write.table(Plt, "X_alcohol_pleiotropy.txt")
write.table(Sin, "X_alcohol_single_snp.txt")
p1 <- mr_scatter_plot(mr_results, dat)
p1[[1]]

res_single <- mr_singlesnp(dat)
p2 <- mr_forest_plot(res_single)
p2[[1]]

res_loo <- mr_leaveoneout(dat)
p3 <- mr_leaveoneout_plot(res_loo)
p3[[1]]

p4 <- mr_funnel_plot(res_single)
p4[[1]]

#Assess outliers
#Radial plots 
devtools::install_github("WSpiller/RadialMR")
library(RadialMR)
library(MendelianRandomization)

dat <- dat[dat$SNP%in%res_single$SNP,]

raddat <- format_radial(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome, dat$SNP)
ivwrad <- ivw_radial(raddat, alpha=0.05/63, weights=3)
dim(ivwrad$outliers)[1]
ivwrad <- ivw_radial(raddat, alpha=0.05, weights=3)
dim(ivwrad$outliers)[1]
eggrad <- egger_radial(raddat, alpha=0.05, weights=3)
eggrad$coef 
dim(eggrad$outliers)[1]

#plot_radial(ivwrad, TRUE, FALSE, TRUE)
plot_radial(c(ivwrad,eggrad), TRUE, FALSE, TRUE)

ivwrad$qstatistic 
ivwrad$sortoutliers <- ivwrad$outliers[order(ivwrad$outliers$p.value),]
ivwrad$sortoutliers$Qsum <- cumsum(ivwrad$sortoutliers$Q_statistic)
ivwrad$sortoutliers$Qdif <- ivwrad$sortoutliers$Qsum - ivwrad$qstatistic
write.csv(ivwrad$sortoutliers, "X_alcohol_outliers.csv", row.names=F, quote=F)

#Remove top outliers
dat2 <- dat[!dat$SNP %in% ivwrad$outliers$SNP,]
mr_results2 <- mr(dat2)
or_results <- generate_odds_ratios(mr_results2)

results<-cbind.data.frame(or_results$outcome,or_results$nsnp,or_results$method,or_results$b,or_results$se,or_results$pval,or_results$or,or_results$or_lci95,or_results$or_uci95)
write.table(results, "X_alcohol_nooutliers.txt")


# Risk tolerance
# Extract SNPs for DBP from Karlsson Linner et al. GWAS

outcome_dat <- read_outcome_data("X_rt.txt", sep="\t")
dat <- harmonise_data(exposure_dat, outcome_dat)
dat$outcome <- "Risk tolerance"
dat$exposure <- "Exposure"
write.csv(dat, "X_rt_harmonised.csv")    
mr_results <- mr(dat)
mr_results
or_results <- generate_odds_ratios(mr_results)
or_results
results<-cbind.data.frame(or_results$outcome,or_results$nsnp,or_results$method,or_results$b,or_results$se,or_results$pval,or_results$or,or_results$or_lci95,or_results$or_uci95)
write.table(results, "X_rt_results.txt")

#Run scatter plot code
mr_scatter <- mr_scatter_plot(mr_results, dat)
mr_scatter
pdf("X_rt_scatter.pdf",  width = 15, height = 20)
ggarrange(mr_scatter$`P51FYM.UQUwFG`,ncol=2, nrow=2, widths = 2, heights = 1)
dev.off()

#6. PLOT FOREST AND LEAVE-ONE-OUT
res_single <- mr_singlesnp(dat)
mr_forest <- mr_forest_plot(res_single)
mr_forest

pdf("X_rt_forest.pdf",  width = 15, height = 20)
ggarrange(mr_forest$`P51FYM.UQUwFG`,ncol=2, nrow=2, widths = 2, heights = 1)
dev.off()

#Leave one out analysis
res_loo <- mr_leaveoneout(dat)
mr_loo <- mr_leaveoneout_plot(res_loo)
mr_loo

pdf("X_ea_loo.pdf",  width = 15, height = 20)
ggarrange(mr_loo$`P51FYM.UQUwFG`, nrow=2, widths = 2, heights = 1)
dev.off()

#Heterogeneity and pleiotropy analysis 
Het<-mr_heterogeneity(dat)
Plt<-mr_pleiotropy_test(dat)
Sin<-mr_singlesnp(dat)
write.table(Het, "X_rt_heterogeneity.txt")
write.table(Plt, "X_rt_pleiotropy.txt")
write.table(Sin, "X_rt_single_snp.txt")
p1 <- mr_scatter_plot(mr_results, dat)
p1[[1]]

res_single <- mr_singlesnp(dat)
p2 <- mr_forest_plot(res_single)
p2[[1]]

res_loo <- mr_leaveoneout(dat)
p3 <- mr_leaveoneout_plot(res_loo)
p3[[1]]

p4 <- mr_funnel_plot(res_single)
p4[[1]]

#Assess outliers
#Radial plots 
devtools::install_github("WSpiller/RadialMR")
library(RadialMR)
library(MendelianRandomization)

dat <- dat[dat$SNP%in%res_single$SNP,]

raddat <- format_radial(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome, dat$SNP)
ivwrad <- ivw_radial(raddat, alpha=0.05/63, weights=3)
dim(ivwrad$outliers)[1] 
ivwrad <- ivw_radial(raddat, alpha=0.05, weights=3)
dim(ivwrad$outliers)[1]
eggrad <- egger_radial(raddat, alpha=0.05, weights=3)
eggrad$coef 
dim(eggrad$outliers)[1]

#plot_radial(ivwrad, TRUE, FALSE, TRUE)
plot_radial(c(ivwrad,eggrad), TRUE, FALSE, TRUE)

ivwrad$qstatistic 
ivwrad$sortoutliers <- ivwrad$outliers[order(ivwrad$outliers$p.value),]
ivwrad$sortoutliers$Qsum <- cumsum(ivwrad$sortoutliers$Q_statistic)
ivwrad$sortoutliers$Qdif <- ivwrad$sortoutliers$Qsum - ivwrad$qstatistic
write.csv(ivwrad$sortoutliers, "X_rt_outliers.csv", row.names=F, quote=F)

#Remove top outliers
dat2 <- dat[!dat$SNP %in% ivwrad$outliers$SNP,]
mr_results2 <- mr(dat2)
or_results <- generate_odds_ratios(mr_results2)

results<-cbind.data.frame(or_results$outcome,or_results$nsnp,or_results$method,or_results$b,or_results$se,or_results$pval,or_results$or,or_results$or_lci95,or_results$or_uci95)
write.table(results, "X_rt_nooutliers.txt")


# Educational attainment
# Extract SNPs for DBP from Lee et al. GWAS

outcome_dat <- read_outcome_data("X_ea.txt", sep="\t")
dat <- harmonise_data(exposure_dat, outcome_dat)
dat$outcome <- "Educational attainment"
dat$exposure <- "Exposure"
write.csv(dat, "X_ea_harmonised.csv")    
mr_results <- mr(dat)
mr_results
or_results <- generate_odds_ratios(mr_results)
or_results
results<-cbind.data.frame(or_results$outcome,or_results$nsnp,or_results$method,or_results$b,or_results$se,or_results$pval,or_results$or,or_results$or_lci95,or_results$or_uci95)
write.table(results, "X_ea_results.txt")

#Run scatter plot
mr_scatter <- mr_scatter_plot(mr_results, dat)
mr_scatter
pdf("X_ea_scatter.pdf",  width = 15, height = 20)
ggarrange(mr_scatter$`YLCgqT.JbgZA7`,ncol=2, nrow=2, widths = 2, heights = 1)
dev.off()

#6. PLOT FOREST AND LEAVE-ONE-OUT
res_single <- mr_singlesnp(dat)
mr_forest <- mr_forest_plot(res_single)
mr_forest

pdf("X_ea_forest.pdf",  width = 15, height = 20)
ggarrange(mr_forest$`YLCgqT.JbgZA7`,ncol=2, nrow=2, widths = 2, heights = 1)
dev.off()

#Leave one out analysis
res_loo <- mr_leaveoneout(dat)
mr_loo <- mr_leaveoneout_plot(res_loo)
mr_loo

pdf("X_ea_loo.pdf",  width = 15, height = 20)
ggarrange(mr_loo$`YLCgqT.JbgZA7`, nrow=2, widths = 2, heights = 1)
dev.off()

#Heterogeneity and pleiotropy analysis 
Het<-mr_heterogeneity(dat)
Plt<-mr_pleiotropy_test(dat)
Sin<-mr_singlesnp(dat)
write.table(Het, "X_ea_heterogeneity.txt")
write.table(Plt, "X_ea_pleiotropy.txt")
write.table(Sin, "X_ea_single_snp.txt")
p1 <- mr_scatter_plot(mr_results, dat)
p1[[1]]

res_single <- mr_singlesnp(dat)
p2 <- mr_forest_plot(res_single)
p2[[1]]

res_loo <- mr_leaveoneout(dat)
p3 <- mr_leaveoneout_plot(res_loo)
p3[[1]]

p4 <- mr_funnel_plot(res_single)
p4[[1]]

#Assess outliers
#Radial plots 
devtools::install_github("WSpiller/RadialMR")
library(RadialMR)
library(MendelianRandomization)

dat <- dat[dat$SNP%in%res_single$SNP,]

raddat <- format_radial(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome, dat$SNP)
ivwrad <- ivw_radial(raddat, alpha=0.05/63, weights=3)
dim(ivwrad$outliers)[1]
ivwrad <- ivw_radial(raddat, alpha=0.05, weights=3)
dim(ivwrad$outliers)[1]
eggrad <- egger_radial(raddat, alpha=0.05, weights=3)
eggrad$coef 
dim(eggrad$outliers)[1]

#plot_radial(ivwrad, TRUE, FALSE, TRUE)
plot_radial(c(ivwrad,eggrad), TRUE, FALSE, TRUE)

ivwrad$qstatistic 
ivwrad$sortoutliers <- ivwrad$outliers[order(ivwrad$outliers$p.value),]
ivwrad$sortoutliers$Qsum <- cumsum(ivwrad$sortoutliers$Q_statistic)
ivwrad$sortoutliers$Qdif <- ivwrad$sortoutliers$Qsum - ivwrad$qstatistic
write.csv(ivwrad$sortoutliers, "X_ea_outliers.csv", row.names=F, quote=F)

#Remove top outliers
dat2 <- dat[!dat$SNP %in% ivwrad$outliers$SNP,]
mr_results2 <- mr(dat2)
or_results <- generate_odds_ratios(mr_results2)

results<-cbind.data.frame(or_results$outcome,or_results$nsnp,or_results$method,or_results$b,or_results$se,or_results$pval,or_results$or,or_results$or_lci95,or_results$or_uci95)
write.table(results, "X_ea_nooutliers.txt")


