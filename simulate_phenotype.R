library(MASS)
library(corpcor)
source("sim_phenotype.R")
source("estVarComp.R")
source("allele_drop_functions.R")
source("assocTestMixedModel_v7.R")


n <- 1600

TEMP=c(2,1,2,1,1,1,2,2,1,2,1,2,1,2,1,2)
SEX <- rep(NA,length(TEMP))
SEX[TEMP==1]="M"
SEX[TEMP==2]="F"

sex <- rep(SEX,(n/16))

# get x chromosome kinship matrix
# get autosomal kinship matrix

# kinship matrix, autosomes
library(kinship2)
tmp <- read.table("pedigree_16individs.txt")
tmp[tmp[,5]==0,5] <- 2
ped <- pedigree(tmp[,2],dadid=tmp[,3],momid=tmp[,4],sex=tmp[,5])
kin <- 2*kinship(ped)

# make it block diagonal
kinFull <- bdiag(kin,kin,kin,kin,kin,
                 kin,kin,kin,kin,kin)
kinFull <- bdiag(kinFull,kinFull,kinFull,kinFull,kinFull,
                 kinFull,kinFull,kinFull,kinFull,kinFull)
dim(kinFull) # 1600 1600, as expected
kinFull <- as.matrix(kinFull)
save(kinFull,file="100Peds_autoKinship.RData")

# kinship matrix, x chromosome
tmp <- read.table("pedigree_16individs_output")
trueKinX <- matrix(NA,nrow=16, ncol=16)
trueKinX[lower.tri(trueKinX,diag=TRUE)] <- tmp[,"V4"]
trueKinXn <- trueKinX
trueKinX <- data.frame(trueKinX)
trueKinX$SEX <- SEX
trueKinXn[trueKinX$SEX=="F",trueKinX$SEX=="F"] <- 2*trueKinXn[trueKinX$SEX=="F",trueKinX$SEX=="F"]

# is this correct???
trueKinXn[trueKinX$SEX=="M",trueKinX$SEX=="F"] <- sqrt(2)*trueKinXn[trueKinX$SEX=="M",trueKinX$SEX=="F"]
trueKinXn[trueKinX$SEX=="F",trueKinX$SEX=="M"] <- sqrt(2)*trueKinXn[trueKinX$SEX=="F",trueKinX$SEX=="M"]

diag(trueKinXn) <- 1
trueKinXn[upper.tri(trueKinXn,diag=FALSE)] <- t(trueKinXn)[upper.tri(trueKinXn,diag=FALSE)]

kinFullX <- bdiag(trueKinXn,trueKinXn,trueKinXn,trueKinXn,trueKinXn,
                  trueKinXn,trueKinXn,trueKinXn,trueKinXn,trueKinXn)
kinFullX <- bdiag(kinFullX,kinFullX,kinFullX,kinFullX,kinFullX,
                  kinFullX,kinFullX,kinFullX,kinFullX,kinFullX)
dim(kinFullX) # 1600 1600
kinFullX <- as.matrix(kinFullX)
save(kinFullX,file="100Peds_xKinship.RData")

rm(list=ls())

######
# add unrelateds to these kinship matrices
kinFull <- get(load("100Peds_autoKinship.RData"))
# add in 500 unrelateds to the end of this matrix
kinFullUnr <- matrix(0,nrow=2100,ncol=2100)
diag(kinFullUnr) <- 1
kinFullUnr[1:1600,1:1600] <- kinFull
save(kinFullUnr,file="100Peds_500unr_autoKinship.RData")

kinFullX <- get(load("100Peds_xKinship.RData"))
kinFullXUnr <- matrix(0,nrow=2100,ncol=2100)
diag(kinFullXUnr) <- 1
kinFullUnr[1:1600,1:1600] <- kinFullX
save(kinFullUnr,file="100Peds_500unr_xKinship.RData")

rm(list=ls())

#####
# these don't depend on the genotypes
sigmaA <- 0.1
sigmaX <- 1.0
sigmaE <- 1.0

totalSigmaScalar <- sigmaE+sigmaA+sigmaX
totalSigma <- sigmaE+sigmaA*kinFull+sigmaX*kinFullX

# choose effect size based on herit, allele freq
herit <- 0.005
p <- 0.2
(beta1 <- getEffSize(p,totalSigmaScalar,herit=herit)) # 0.0128



#####
# type I error rate
# performed 10 simulations where I generated 1000 SNPs, the first of which is causal
# thus, we have 999 null SNPs per simulation

# read in the results
res <- get(load("simResults_20_p2/simRes_1_1.RData"))
simResults <- data.frame(matrix(NA,nrow=11,ncol=19))
names(simResults) <- c(names(res$metricsSet),"false_X0001","false_0001","false_onlyX_0001",
                       "false_0005X","false_0005","false_0005onlyX",
                       "false_pos_001X","false_pos_001","false_pos_001onlyX",
                       "false_pos_01X","false_pos_01","false_pos_01onlyX","nSim")

for(j in 1:11){
  fals001X <- 0; fals001 <- 0; fals001onlyX <- 0
  fals0005X <- 0; fals0005 <- 0; fals0005onlyX <- 0
  Xfals0001 <- 0; fals0001 <- 0; onlyXfals0001 <- 0
  fals01X <- 0; fals01 <- 0; fals01onlyX <- 0
  
  for(i in 1:10){
  fn <- paste("simResults_20_p2/simRes_",i,"_",j,".RData",sep="")
  res <- get(load(fn))
  fn <- paste("simResults_20_p2/simRes_xOnly_",i,"xOnly_",j,".RData",sep="")
  resXonly <- get(load(fn))
  
  # type I error is proportion of false positives
  alpha <- 1e-04
  Xfals0001 <- Xfals0001 + sum(res$mmRes_adjForX[-1,"pval"]<alpha)
  fals0001 <- fals0001 + sum(res$mmRes_notAdjForX[-1,"pval"]<alpha)
  onlyXfals0001 <- onlyXfals0001 + sum(resXonly$mmRes_onlyAdjForX[-1,"pval"]<alpha)

  alpha <- 5e-04
  fals0005X <- fals0005X + sum(res$mmRes_adjForX[-1,"pval"]<alpha)
  fals0005 <- fals0005 + sum(res$mmRes_notAdjForX[-1,"pval"]<alpha)
  fals0005onlyX <- fals0005onlyX + sum(resXonly$mmRes_onlyAdjForX[-1,"pval"]<alpha)
  
  alpha <- 0.001
  fals001X <- fals001X + sum(res$mmRes_adjForX[-1,"pval"]<alpha)
  fals001 <- fals001 + sum(res$mmRes_notAdjForX[-1,"pval"]<alpha)
  fals001onlyX <- fals001onlyX + sum(resXonly$mmRes_onlyAdjForX[-1,"pval"]<alpha)
    
  alpha <- 0.01
  fals01X <- fals01X + sum(res$mmRes_adjForX[-1,"pval"]<alpha)
  fals01 <- fals01 + sum(res$mmRes_notAdjForX[-1,"pval"]<alpha)
  fals01onlyX <- fals01onlyX + sum(resXonly$mmRes_onlyAdjForX[-1,"pval"]<alpha)
  }
  
  simResults[j,1:6] <- res$metricsSet
  simResults$false_X0001[j] <- Xfals0001
  simResults$false_0001[j] <- fals0001
  simResults$false_onlyX_0001[j] <- onlyXfals0001
  
  simResults$false_0005X[j] <- fals0005X
  simResults$false_0005[j] <- fals0005
  simResults$false_0005onlyX[j] <- fals0005onlyX
  
  simResults$false_pos_001X[j] <- fals001X
  simResults$false_pos_001[j] <- fals001
  simResults$false_pos_001onlyX[j] <- fals001onlyX
  
  simResults$false_pos_01X[j] <- fals01X
  simResults$false_pos_01[j] <- fals01
  simResults$false_pos_01onlyX[j] <- fals01onlyX
  
  simResults$nSim[j] <- length(res$mmRes_adjForX[-1,"pval"])*10
}
simResults <- simResults[2:7,]
ty1 <- colSums(simResults)/sum(simResults$nSim)

sum(simResults$nSim) # 29940 simulations, 5K for each of the param values

ty1[7:18]
#false_X0001         false_0001   false_onlyX_0001        false_0005X 
#0.0002368373       0.0007833850       0.0005465476       0.0030242303 
#false_0005    false_0005onlyX     false_pos_001X      false_pos_001 
#0.0036618692       0.0034614684       0.0054290399       0.0059573693 
#false_pos_001onlyX      false_pos_01X       false_pos_01  false_pos_01onlyX 
#0.0062306431       0.0280925487       0.0325013664       0.0295317909 

write.table(simResults,file="typeIerror_10Ksims.txt",sep="\t")

# make a latex table of these results
library(xtable)
# 3 rows: adj for x & autos, adj for just autos, adj for just x
xtable(as.matrix(t(cbind(ty1[c(7,10,13,16)],ty1[c(8,11,14,17)],ty1[c(9,12,15,18)]))),digits=5)
xtable(simResults[,1:6],digits=c(1,5,2,1,1,0,0))

###
# what is the allele frequency for which the false positives are occurring??
# try filtering out SNPs with MAF<0.01 and see if the type I error rate is the same

res <- get(load("simResults_20_p2/simRes_1_1.RData"))
simResults <- data.frame(matrix(NA,nrow=11,ncol=19))
names(simResults) <- c(names(res$metricsSet),"false_X0001","false_0001","false_onlyX_0001",
                       "false_0005X","false_0005","false_0005onlyX",
                       "false_pos_001X","false_pos_001","false_pos_001onlyX",
                       "false_pos_01X","false_pos_01","false_pos_01onlyX","nSim")

for(j in 1:11){
  fals001X <- 0; fals001 <- 0; fals001onlyX <- 0
  fals5e4X <- 0; fals5 <- 0; false4onlyX <- 0
  Xfals0001 <- 0; fals0001 <- 0; onlyXfals0001 <- 0
  fals01X <- 0; fals01 <- 0; fals01onlyX <- 0
  nSim <- 0
    
  for(i in 1:10){
    fn <- paste("simResults_20_p2/simRes_",i,"_",j,".RData",sep="")
    res <- get(load(fn))
    
    res$mmRes_adjForX <- res$mmRes_adjForX[res$mmRes_adjForX[,"MAF"]>0.01,]
    res$mmRes_notAdjForX <- res$mmRes_notAdjForX[res$mmRes_notAdjForX[,"MAF"]>0.01,]
    
    fn <- paste("simResults_20_p2/simRes_xOnly_",i,"xOnly_",j,".RData",sep="")
    resXonly <- get(load(fn))
    resXonly$mmRes_onlyAdjForX <- resXonly$mmRes_onlyAdjForX[resXonly$mmRes_onlyAdjForX[,"MAF"]>0.01,]
    
    # type I error is proportion of false positives
    alpha <- 1e-04
    Xfals0001 <- Xfals0001 + sum(res$mmRes_adjForX[-1,"pval"]<alpha)
    fals0001 <- fals0001 + sum(res$mmRes_notAdjForX[-1,"pval"]<alpha)
    onlyXfals0001 <- onlyXfals0001 + sum(resXonly$mmRes_onlyAdjForX[-1,"pval"]<alpha)
    
    alpha2 <- 5e-04
    fals5e4X <- fals5e4X + sum(res$mmRes_adjForX[-1,"pval"]<alpha2)
    fals5 <- fals5 + sum(res$mmRes_notAdjForX[-1,"pval"]<alpha2)
    false4onlyX <- false4onlyX + sum(resXonly$mmRes_onlyAdjForX[-1,"pval"]<alpha2)
    
    alpha <- 0.001
    fals001X <- fals001X + sum(res$mmRes_adjForX[-1,"pval"]<alpha)
    fals001 <- fals001 + sum(res$mmRes_notAdjForX[-1,"pval"]<alpha)
    fals001onlyX <- fals001onlyX + sum(resXonly$mmRes_onlyAdjForX[-1,"pval"]<alpha)
    
    alpha <- 0.01
    fals01X <- fals01X + sum(res$mmRes_adjForX[-1,"pval"]<alpha)
    fals01 <- fals01 + sum(res$mmRes_notAdjForX[-1,"pval"]<alpha)
    fals01onlyX <- fals01onlyX + sum(resXonly$mmRes_onlyAdjForX[-1,"pval"]<alpha)
    
    nSim <- nSim + length(res$mmRes_adjForX[-1,"pval"])
  }
  
  simResults[j,1:6] <- res$metricsSet
  simResults$false_X0001[j] <- Xfals0001
  simResults$false_0001[j] <- fals0001
  simResults$false_onlyX_0001[j] <- onlyXfals0001
  
  simResults$false_0005X[j] <- fals5e4X
  simResults$false_0005[j] <- fals5
  simResults$false_0005onlyX[j] <- false4onlyX
  
  simResults$false_pos_001X[j] <- fals001X
  simResults$false_pos_001[j] <- fals001
  simResults$false_pos_001onlyX[j] <- fals001onlyX
  
  simResults$false_pos_01X[j] <- fals01X
  simResults$false_pos_01[j] <- fals01
  simResults$false_pos_01onlyX[j] <- fals01onlyX
  
  simResults$nSim[j] <- nSim #instead of: length(res$mmRes_adjForX[-1,"pval"])*10
}

simResults <- simResults[2:7,]
ty1 <- colSums(simResults)/sum(simResults$nSim)
sum(simResults$nSim) # 25770

library(xtable)
# 3 rows: adj for x & autos, adj for just autos, adj for just x
xtable(as.matrix(t(cbind(ty1[c(7,10,13,16)],ty1[c(8,11,14,17)],ty1[c(9,12,15,18)]))),digits=5)
xtable(simResults[,1:6],digits=c(1,5,2,1,1,0,0))




#### 
# get the power results
# all analyses included 2 SNPs
dat <- get(load("simResults_20_p2/simRes_truePos_loop1_1_1.RData"))
powerRes <- data.frame(matrix(NA,nrow=11,ncol=16))
colnames(powerRes) <- c(names(dat$metricsSet),"p_x01","p_01","p_onlyX_01","x_005","nX_005","xOnly_005",
                        "p_x001","p_001","p_001onlyX","nSims")

for(k in 1:11){
powerX <- 0; power <- 0; poweronlyX <- 0
p_X1 <- 0; p_1 <- 0; p_1onlyX <- 0
pX <- 0; p <- 0; ponlyX <- 0

denom <- 0
for(i in 1:1000){
  for(j in 1:10){
    fn <- paste("simResults_20_p2/simRes_truePos_loop",j,"_",i,"_",k,".RData",sep="")
    dat <- get(load(fn))
    
    fn <- paste("simResults_20_p2/simRes_truePos_xOnly_loop",j,"_",i,"_",k,".RData",sep="")
    datonlyX <- get(load(fn))
    
    alpha <- 0.01
    powerX <- powerX + sum(dat$mmRes_adjForX[1,"pval"]<alpha)
    power <- power + sum(dat$mmRes_notAdjForX[1,"pval"]<alpha)
    poweronlyX <- poweronlyX + sum(datonlyX$mmRes_adjForXonly[1,"pval"]<alpha)
    denom <- denom + 1
  
    alpha <- 0.001
    p_X1 <- p_X1 + sum(dat$mmRes_adjForX[1,"pval"]<alpha)
    p_1 <- p_1 + sum(dat$mmRes_notAdjForX[1,"pval"]<alpha)
    p_1onlyX <- p_1onlyX + sum(datonlyX$mmRes_adjForXonly[1,"pval"]<alpha)
    
    alpha <- 0.005
    pX <- pX + sum(dat$mmRes_adjForX[1,"pval"]<alpha)
    p <- p + sum(dat$mmRes_notAdjForX[1,"pval"]<alpha)
    ponlyX <- ponlyX + sum(datonlyX$mmRes_adjForXonly[1,"pval"]<alpha)
    
  }
  powerRes[k,1:6] <- dat$metricsSet

  powerRes$p_x01[k] <- powerX
  powerRes$p_01[k] <- power
  powerRes$p_onlyX_01[k] <- poweronlyX
  
  powerRes$x_005[k] <- pX
  powerRes$nX_005[k] <- p
  powerRes$xOnly_005[k] <- ponlyX
  
  powerRes$p_x001[k] <- p_X1
  powerRes$p_001[k] <- p_1
  powerRes$p_001onlyX[k] <- p_1onlyX

  powerRes$nSims[k] <- denom
}
}

write.table(powerRes,file="powerRes_10000Iters.txt",sep="\t")

# make a plot of these values
powerRes <- powerRes[2:7,]
pdf("powerRes_10000Iters.pdf")
plot(x=powerRes$beta,y=powerRes$p_x01/powerRes$nSims,ylim=c(0,1),type="l",col="blue",
     xlab="Effect Size",ylab="Power",cex.lab=1.2,cex.axis=1.2,lwd=1.2)
points(x=powerRes$beta,y=powerRes$p_01/powerRes$nSims,type="l",lty=3,col="blue",lwd=1.5)
points(x=powerRes$beta,y=powerRes$p_onlyX_01/powerRes$nSims,type="l",lty=5,col="blue",lwd=1.2)

points(x=powerRes$beta,y=powerRes$x_005/powerRes$nSims,type="l",col="gold",lwd=1.2)
points(x=powerRes$beta,y=powerRes$nX_005/powerRes$nSims,type="l",lty=3,col="gold",lwd=1.5)
points(x=powerRes$beta,y=powerRes$xOnly_005/powerRes$nSims,type="l",lty=5,col="gold",lwd=1.2)

points(x=powerRes$beta,y=powerRes$p_x001/powerRes$nSims,type="l",col="magenta",lwd=1.2)
points(x=powerRes$beta,y=powerRes$p_001/powerRes$nSims,type="l",lty=3,col="magenta",lwd=1.5)
points(x=powerRes$beta,y=powerRes$p_001onlyX/powerRes$nSims,type="l",lty=5,col="magenta",lwd=1.2)

abline(h=0.8,col="gray",lwd=0.8)
legend("bottomright",c(expression(paste(alpha," = 0.001")),expression(paste(alpha," = 0.005")),
                       expression(paste(alpha," = 0.01")),"Auto + X adj","Auto adj","X adj"),
       col=c("magenta","gold","blue","black","black","black"),lty=c(1,1,1,1,3,5),lwd=2,cex=1.5)
dev.off()






#################
# process results for different param values

library(BiocParallel)
dat <- get(load("simResults_vars/simRes_1.RData"))
length(dat) # 24; these are the different param combos
# the _X.RData in the object name is the iteration

metrics <- unlist(bplapply(dat,function(x){x$metrics}))
metrics <- data.frame(t(matrix(metrics,ncol=24)))
colnames(metrics) <- c("beta1","herit_snp","p","sigmaA","sigmaX","sigmaE")

# now loop through the list and get the type i error for each of the three models fit, for varying alpha levels
alpha1 <- 0.01
alpha2 <- 0.005
alpha3 <- 0.001
alpha4 <- 5e-04
alpha5 <- 1e-04

metrics$xOnly_1 <- NA
metrics$xOnly_2 <- NA
metrics$xOnly_3 <- NA
metrics$xOnly_4 <- NA
metrics$xOnly_5 <- NA

metrics$autoOnly_1 <- NA
metrics$autoOnly_2 <- NA
metrics$autoOnly_3 <- NA
metrics$autoOnly_4 <- NA
metrics$autoOnly_5 <- NA

metrics$xAuto_1 <- NA
metrics$xAuto_2 <- NA
metrics$xAuto_3 <- NA
metrics$xAuto_4 <- NA
metrics$xAuto_5 <- NA


for(i in 1:length(dat)){

  x1=0; x2=0; x3=0; x4=0; x5=0
  a1=0; a2=0; a3=0; a4=0; a5=0
  xa1=0; xa2=0; xa3=0; xa4=0; xa5=0
  simN=0
  
  for(j in 1:10){
  dat <- get(load(paste("simResults_vars/simRes_",j,".RData",sep="")))
    
  thisRes <- dat[[i]]
  
  x1= x1+sum(thisRes$mmRes_onlyAdjForX[-1,"pval"]<alpha1)
  x2=x2+sum(thisRes$mmRes_onlyAdjForX[-1,"pval"]<alpha2)
  x3=x3+sum(thisRes$mmRes_onlyAdjForX[-1,"pval"]<alpha3)
  x4=x4+sum(thisRes$mmRes_onlyAdjForX[-1,"pval"]<alpha4)
  x5=x5+sum(thisRes$mmRes_onlyAdjForX[-1,"pval"]<alpha5)
  
  a1=a1+sum(thisRes$mmRes_autoOnly[-1,"pval"]<alpha1)
  a2=a2+sum(thisRes$mmRes_autoOnly[-1,"pval"]<alpha2)
  a3=a3+sum(thisRes$mmRes_autoOnly[-1,"pval"]<alpha3)
  a4=a4+sum(thisRes$mmRes_autoOnly[-1,"pval"]<alpha4)
  a5=a5+sum(thisRes$mmRes_autoOnly[-1,"pval"]<alpha5)
  
  xa1=xa1+sum(thisRes$mmRes_both_Xauto[-1,"pval"]<alpha1)
  xa2=xa2+sum(thisRes$mmRes_both_Xauto[-1,"pval"]<alpha2)
  xa3=xa3+sum(thisRes$mmRes_both_Xauto[-1,"pval"]<alpha3)
  xa4=xa4+sum(thisRes$mmRes_both_Xauto[-1,"pval"]<alpha4)
  xa5=xa5+sum(thisRes$mmRes_both_Xauto[-1,"pval"]<alpha5)
  
  simN=simN+nrow(thisRes$mmRes_both_Xauto[-1,])
  
  }  # end the loop through the iterations
  
  metrics$xOnly_1[i] <- x1
  metrics$xOnly_2[i] <- x2
  metrics$xOnly_3[i] <- x3
  metrics$xOnly_4[i] <- x4
  metrics$xOnly_5[i] <- x5
  
  metrics$autoOnly_1[i] <- a1
  metrics$autoOnly_2[i] <- a2
  metrics$autoOnly_3[i] <- a3
  metrics$autoOnly_4[i] <- a4
  metrics$autoOnly_5[i] <- a5
  
  metrics$xAuto_1[i] <- xa1
  metrics$xAuto_2[i] <- xa2
  metrics$xAuto_3[i] <- xa3
  metrics$xAuto_4[i] <- xa4
  metrics$xAuto_5[i] <- xa5
  
  metrics$nSim[i] <- simN
}

write.table(metrics,file="metrics_vars_typeIresults.txt",sep="\t")

rm(list=ls())

metrics <- read.table("metrics_vars_typeIresults.txt",sep="\t")
colSums(metrics)/sum(metrics$nSim)

# make a table of these results for latex
library(xtable)
ty1 <- colSums(metrics)/sum(metrics$nSim)
xtable(t(matrix(c(ty1[7:11],ty1[12:16],ty1[17:21]),nrow=5)),digits=5)

xtable(metrics[,c(1,2,3,4,5,6,22,7,12,17,8,13,18,9,14,19,10,15,20,
                  11,16,21)],digits=c(1,5,3,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))

tabLat <- metrics[,c(1,2,3,4,5,6,22,7,12,17,8,13,18,9,14,19,10,15,20,
                     11,16,21)]
xtable(matrix(colSums(tabLat),nrow=1),digits=0) # totals row
xtable(matrix(colSums(tabLat)/sum(tabLat$nSim),nrow=1),digits=5) # percentage row

# calculate the herit for all x chromosome SNPs
attach(metrics)
herit <- (beta1^2*2*p*(1-p)+sigmaX)/(beta1^2*2*p*(1-p)+sigmaX+sigmaE+sigmaA)
detach(metrics)

metrics$herit_xchr <- herit



