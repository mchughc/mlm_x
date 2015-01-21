library(MASS)
library(corpcor)
source("sim_phenotype.R")
source("estVarComp.R")
source("allele_drop_functions.R")
source("assocTestMixedModel_v7.R")


## Contents:
# 1. Batch results processing: type I error
#  KC matrices
# 2. Make table of type I error results
# 3. Parse estimates for var components using KC matrices
# 3b. Make plots of individual metrics and estimate for 500 iters
# 4. Process results using 2KC for auto and 2,sqrt(2),1 KC for X chr
# 4b. Make plots of individual metrics and estimate for 500 iters
# 5. Create table of metrics for the 24 iterations
# 6. Parse the 2KC for both auto and X chr results
# 7. Look at type I error rate with mf KC matrices
# 8. Get results for simple model with 10K pedigrees
# 9. Get results for more iterations of model A: y=beta1*SNP_x + e
# 10. More simulations
# 11. Type I error 
# 12. Power calcs: check phenotype changes with different beta1 values
# 13. Power simulations
# 14. Power simulation results
# 15. Make power graphs
# 16. Get type I error from the same results
# 17. More power simulations



#####
# 1. Batch results processing: type I error

# cd /projects/geneva/geneva_sata/caitlin/mlm_x
# qsub -N vars batch_simulations_equalvars.sh

# process results for different param values
library(BiocParallel)

dat <- get(load("simResults_vars/simRes_1.RData"))
length(dat) # 24; these are the different param combos
# the _X.RData in the object name is the iteration

metrics <- unlist(bplapply(dat,function(x){x$metrics}))
metrics <- data.frame(t(matrix(metrics,ncol=24)))
colnames(metrics) <- c("beta1","herit2_snp","p","sigmaA","sigmaX","sigmaE")

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


#####
# 2. Make table of type I error results

metrics <- read.table("metrics_vars_typeIresults.txt",sep="\t")
colSums(metrics)/sum(metrics$nSim)

# make a table of these results for latex
library(xtable)
ty1 <- colSums(metrics)/sum(metrics$nSim)
xtable(t(matrix(c(ty1[12:16],ty1[7:11],ty1[17:21]),nrow=5)),digits=5)

# calculate the herit for all x chromosome SNPs
attach(metrics)
herit <- (beta1^2*3*p*(1-p)+sigmaX)/(beta1^2*3*p*(1-p)+sigmaX+sigmaE+sigmaA)
detach(metrics)

metrics$herit_xchr <- herit

xtable(metrics[,c(23,1,2,4,5,7,12,17,8,13,18,9,14,19,10,15,20,
                  11,16,21)],digits=c(1,4,4,3,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))

tabLat <- metrics[,c(1,2,4,5,7,12,17,8,13,18,9,14,19,10,15,20,
                     11,16,21)]
xtable(matrix(colSums(tabLat),nrow=1),digits=0) # totals row
xtable(matrix(colSums(tabLat)/sum(metrics$nSim),nrow=1),digits=5) # percentage row


# a smaller table for the type I error results
xtable(as.matrix(cbind(ty1[c(7,12,17)],ty1[c(8,13,18)],ty1[c(9,14,19)],
                       ty1[c(10,15,20)],ty1[c(11,16,21)])),digits=5)

rm(list=ls())


#####
# 3. Parse estimates for var components using KC matrices

# submitted to cluster:
# qsub -N comps -t 1-50 batch_simulations_equalvars_v2.sh

# for each metric value, calculate the expected var component
# then compare to estimated var comp

meanMatrix <- matrix(0,nrow=3,ncol=13)
colnames(meanMatrix) <- c("true","all3_mean","all3_sd","justAuto_mean","justAuto_sd","justX_mean","justX_sd",
                          "all3_upperMean","all3_lowerMean","justAuto_upperMean","justAuto_lowerMean",
                          "justX_upperMean","justX_lowerMean")
rownames(meanMatrix) <- c("sigmaX","sigmaA","sigmaE")

meanML <- list(meanMatrix,meanMatrix,meanMatrix,meanMatrix,
               meanMatrix,meanMatrix,meanMatrix,meanMatrix,
               meanMatrix,meanMatrix,meanMatrix,meanMatrix,
               meanMatrix,meanMatrix,meanMatrix,meanMatrix,
               meanMatrix,meanMatrix,meanMatrix,meanMatrix,
               meanMatrix,meanMatrix,meanMatrix,meanMatrix)

dat <- get(load("simResults_vars/simRes_varCompsOnly_1_1.RData"))

for(i in 1:length(dat)){
  
  sumX <- NULL; sumE <- NULL
  sumA <- NULL; sumEAuto <- NULL
  sumAbo <- NULL; sumXbo <- NULL; sumEbo <- NULL
  
  sumLabo <- NULL; sumUabo <- NULL
  sumLxbo <- NULL; sumUxbo <- NULL
  sumLebo <- NULL; sumUebo <- NULL
  
  sumLa <- NULL; sumUa <- NULL
  sumLeauto <- NULL; sumUeauto <- NULL
  
  sumLx <- NULL; sumUx <- NULL
  sumLex <- NULL; sumUex <- NULL
  
  for(j in 1:50){
  for(k in 1:10){
    dat <- get(load(paste("simResults_vars/simRes_varCompsOnly_",j,"_",k,".RData",sep="")))
    thisRes <- dat[[i]]
    
    x <- thisRes$ci_inclXonly
    sumX <- c(sumX,x["V_kinshipX","Est"])
    sumE <- c(sumE,x["V_E","Est"])
    sumLx <- c(sumLx,x["V_kinshipX","Lower.95"])
    sumUx <- c(sumUx,x["V_kinshipX","Upper.95"])
    sumLex <- c(sumLx,x["V_E","Lower.95"])
    sumUex <- c(sumUx,x["V_E","Upper.95"])
    
    auto <- thisRes$ci_autoOnly
    sumA <- c(sumA,auto["V_kinshipAuto","Est"])
    sumEAuto <- c(sumEAuto,auto["V_E","Est"])
    sumLa <- c(sumLa,auto["V_kinshipAuto","Lower.95"])
    sumUa <- c(sumUa,auto["V_kinshipAuto","Upper.95"])
    sumLeauto <- c(sumLeauto,auto["V_E","Lower.95"])
    sumUeauto <- c(sumUeauto,auto["V_E","Upper.95"])
    
    bo <- thisRes$ci_autoX
    sumAbo <- c(sumAbo,bo["V_kinshipAuto","Est"])
    sumXbo <- c(sumXbo,bo["V_kinshipX","Est"])
    sumEbo <- c(sumEbo,bo["V_E","Est"])
    
    sumLabo <- c(sumLabo,bo["V_kinshipAuto","Lower.95"])
    sumUabo <- c(sumUabo,bo["V_kinshipAuto","Upper.95"])
    sumLxbo <- c(sumLxbo,bo["V_kinshipX","Lower.95"])
    sumUxbo <- c(sumUxbo,bo["V_kinshipX","Upper.95"])
    sumLebo <- c(sumLebo,bo["V_E","Lower.95"])
    sumUebo <- c(sumUebo,bo["V_E","Upper.95"])
        
  }}  # end the loop through the iterations
  
  m <- thisRes$metricsSet
  hx <- m$beta^2*3*m$p*(1-m$p)+m$sigmaX
  
  meanML[[i]]["sigmaX","true"] <- hx
  meanML[[i]]["sigmaE","true"] <- m$sigmaE
  meanML[[i]]["sigmaA","true"] <- m$sigmaA
  
  meanML[[i]]["sigmaX","all3_mean"] <- mean(sumXbo)
  meanML[[i]]["sigmaX","all3_sd"] <- sd(sumXbo)

  meanML[[i]]["sigmaA","all3_mean"] <- mean(sumAbo)
  meanML[[i]]["sigmaA","all3_sd"] <- sd(sumAbo)
  
  meanML[[i]]["sigmaE","all3_mean"] <- mean(sumEbo)
  meanML[[i]]["sigmaE","all3_sd"] <- sd(sumEbo)
  
  meanML[[i]]["sigmaA","justAuto_mean"] <- mean(sumA)
  meanML[[i]]["sigmaA","justAuto_sd"] <- sd(sumA)
  
  meanML[[i]]["sigmaE","justAuto_mean"] <- mean(sumEAuto)
  meanML[[i]]["sigmaE","justAuto_sd"] <- sd(sumEAuto)
  
  meanML[[i]]["sigmaX","justX_mean"] <- mean(sumX)
  meanML[[i]]["sigmaX","justX_sd"] <- sd(sumX)
  
  meanML[[i]]["sigmaE","justX_mean"] <- mean(sumE)
  meanML[[i]]["sigmaE","justX_sd"] <- sd(sumE)
  
  meanML[[i]]["sigmaX","all3_upperMean"] <- paste("(",format(mean(sumLxbo),digits=5),", ",format(mean(sumUxbo),digits=5),")",sep="")
  meanML[[i]]["sigmaX","all3_lowerMean"] <- mean(sumLxbo)
  meanML[[i]]["sigmaA","all3_upperMean"] <- paste("(",format(mean(sumLabo),digits=5),", ",format(mean(sumUabo),digits=5),")",sep="")
  meanML[[i]]["sigmaA","all3_lowerMean"] <- mean(sumLabo)  
  meanML[[i]]["sigmaE","all3_upperMean"] <- paste("(",format(mean(sumLebo),digits=5),", ",format(mean(sumUebo),digits=5),")",sep="")
  meanML[[i]]["sigmaE","all3_lowerMean"] <- mean(sumLebo)
  
  meanML[[i]]["sigmaA","justAuto_upperMean"] <-  paste("(",format(mean(sumLa),digits=5),", ",format(mean(sumUa),digits=5),")",sep="")
  meanML[[i]]["sigmaA","justAuto_lowerMean"] <- mean(sumLa)
  meanML[[i]]["sigmaE","justAuto_upperMean"] <- paste("(",format(mean(sumLeauto),digits=5),", ",format(mean(sumUeauto),digits=5),")",sep="")
  meanML[[i]]["sigmaE","justAuto_lowerMean"] <- mean(sumLeauto)
  
  meanML[[i]]["sigmaX","justX_upperMean"] <- paste("(",format(mean(sumLx),digits=5),", ",format(mean(sumUx),digits=5),")",sep="")
  meanML[[i]]["sigmaX","justX_lowerMean"] <- mean(sumLx)
  meanML[[i]]["sigmaE","justX_upperMean"] <- paste("(",format(mean(sumLex),digits=5),", ",format(mean(sumUex),digits=5),")",sep="")
  meanML[[i]]["sigmaE","justX_lowerMean"] <- mean(sumLex)
  
}

# calculate 95% ci's for the mean/sd
for(i in 1:length(meanML)){
  iter <- meanML[[i]]
  meanML[[i]] <- cbind(meanML[[i]],rep(NA,3))
  colnames(meanML[[i]])[ncol(meanML[[i]])] <- "all3_ci"
  meanML[[i]] <- cbind(meanML[[i]],rep(NA,3))
  colnames(meanML[[i]])[ncol(meanML[[i]])] <- "justAuto_ci"
  meanML[[i]] <- cbind(meanML[[i]],rep(NA,3))
  colnames(meanML[[i]])[ncol(meanML[[i]])] <- "justX_ci"
  
  meanML[[i]][,"all3_ci"] <- paste("(",format(as.numeric(iter[,"all3_mean"])-1.96*as.numeric(iter[,"all3_sd"])/sqrt(500),
                                              digits=4),", ",
                               format(as.numeric(iter[,"all3_mean"])+1.96*as.numeric(iter[,"all3_sd"])/sqrt(500),
                                      digits=4),")",sep="")
  meanML[[i]][,"justAuto_ci"] <- paste("(",format(as.numeric(iter[,"justAuto_mean"])-1.96*as.numeric(iter[,"justAuto_sd"])/sqrt(500),
                                                  digits=4),", ",
                                   format(as.numeric(iter[,"justAuto_mean"])+1.96*as.numeric(iter[,"justAuto_sd"])/sqrt(500),
                                          digits=4),")",sep="")
  meanML[[i]][,"justX_ci"] <- paste("(",format(as.numeric(iter[,"justX_mean"])-1.96*as.numeric(iter[,"justX_sd"])/sqrt(500),
                                               digits=4),", ",
                                format(as.numeric(iter[,"justX_mean"])+1.96*as.numeric(iter[,"justX_sd"])/sqrt(500),
                                       digits=4),")",sep="")
}

save(meanML,file="var_comp_estimates_500iter.RData")

library(xtable)
for(i in 1:length(meanML)){
  meanML[[i]][,"true"] <- format(as.numeric(meanML[[i]][,"true"]),digits=4)
  meanML[[i]][,"all3_mean"] <- format(as.numeric(meanML[[i]][,"all3_mean"]),digits=4)
  meanML[[i]][,"justX_mean"] <- format(as.numeric(meanML[[i]][,"justX_mean"]),digits=4)
  meanML[[i]][,"justAuto_mean"] <- format(as.numeric(meanML[[i]][,"justAuto_mean"]),digits=4)
  iter <- meanML[[i]]
  cisM1 <- paste(iter[,c("all3_mean")],iter[,c("all3_upperMean")])
  cisM2 <- paste(iter[,"justX_mean"],iter[,"justX_upperMean"])
  cisM3 <- paste(iter[,"justAuto_mean"],iter[,"justAuto_upperMean"])
  print(xtable(matrix(c(iter[,c("true")],cisM1,cisM2,cisM3),nrow=3)))
}


#####
# 3b. Make plots of individual metrics and estimate for 500 iters

## make a plot for one param of the 500 iters vs the estimate (95% ci)
i=17
dat <- get(load("simResults_vars/simRes_varCompsOnly_1_1.RData"))
thisRes <- dat[[i]]
metrics <- thisRes$metricsSet

# want to store thisRes$ci_autoX
v_auto <- data.frame(matrix(NA,nrow=500,ncol=3))
colnames(v_auto) <- c("Est","Lower.95","Upper.95")

v_x <- data.frame(matrix(NA,nrow=500,ncol=3))
colnames(v_x) <- c("Est","Lower.95","Upper.95")

v_e <- data.frame(matrix(NA,nrow=500,ncol=3))
colnames(v_e) <- c("Est","Lower.95","Upper.95")

thisRow <- 0

  for(j in 1:50){
    for(k in 1:10){
      dat <- get(load(paste("simResults_vars/simRes_varCompsOnly_",j,"_",k,".RData",sep="")))
      thisRes <- dat[[i]]
      
      thisRow <- thisRow+1
        
      v_auto[thisRow,"Est"] <- thisRes$ci_autoX["V_kinshipAuto","Est"]
      v_auto[thisRow,"Lower.95"] <- thisRes$ci_autoX["V_kinshipAuto","Lower.95"]
      v_auto[thisRow,"Upper.95"] <- thisRes$ci_autoX["V_kinshipAuto","Upper.95"]
      
      v_x[thisRow,"Est"] <- thisRes$ci_autoX["V_kinshipX","Est"]
      v_x[thisRow,"Lower.95"] <- thisRes$ci_autoX["V_kinshipX","Lower.95"]
      v_x[thisRow,"Upper.95"] <- thisRes$ci_autoX["V_kinshipX","Upper.95"]
      
      v_e[thisRow,"Est"] <- thisRes$ci_autoX["V_E","Est"]
      v_e[thisRow,"Lower.95"] <- thisRes$ci_autoX["V_E","Lower.95"]
      v_e[thisRow,"Upper.95"] <- thisRes$ci_autoX["V_E","Upper.95"]
      
    }}  # end the loop through the iterations
  m <- metrics
  hx <- m$beta^2*3*m$p*(1-m$p)+m$sigmaX
  

pdf("metrics_17_sigmaX_est.pdf",width=14)
areas <- 1:500
minY <- min(v_x$Lower.95)+0.01; maxY <- max(v_x$Upper.95)+0.01
par(mar=c(4,4,1,1)+0.2)
plot(areas, v_x$Est, type='n', xlab="Iteration",
     ylab="Estimate", cex.lab=1.2,ylim=c(minY,maxY))
polygon(x=c(areas, rev(areas)),
        y=c(v_x$Lower.95, rev(v_x$Upper.95)),
        col='grey90',border=NA)
lines(areas, v_x$Est)
abline(h=hx,col="red",lwd=2,lty=2)
dev.off()

pdf("metrics_17_sigmaA_est.pdf",width=14)
areas <- 1:500
minY <- min(v_auto$Lower.95)+0.01; maxY <- max(v_auto$Upper.95)+0.01
par(mar=c(4,4,1,1)+0.2)
plot(areas, v_auto$Est, type='n', xlab="Iteration",
     ylab="Estimate", cex.lab=1.2,ylim=c(minY,maxY))
polygon(x=c(areas, rev(areas)),
        y=c(v_auto$Lower.95, rev(v_auto$Upper.95)),
        col='grey90',border=NA)
lines(areas, v_auto$Est)
abline(h=m$sigmaA,col="red",lwd=2,lty=2)
dev.off()

##### 
# try again with a different set of metrics
i=24
dat <- get(load("simResults_vars/simRes_varCompsOnly_1_1.RData"))
thisRes <- dat[[i]]
metrics <- thisRes$metricsSet

# want to store thisRes$ci_autoX
v_auto <- data.frame(matrix(NA,nrow=500,ncol=3))
colnames(v_auto) <- c("Est","Lower.95","Upper.95")

v_x <- data.frame(matrix(NA,nrow=500,ncol=3))
colnames(v_x) <- c("Est","Lower.95","Upper.95")

v_e <- data.frame(matrix(NA,nrow=500,ncol=3))
colnames(v_e) <- c("Est","Lower.95","Upper.95")

thisRow <- 0

for(j in 1:50){
  for(k in 1:10){
    dat <- get(load(paste("simResults_vars/simRes_varCompsOnly_",j,"_",k,".RData",sep="")))
    thisRes <- dat[[i]]
    
    thisRow <- thisRow+1
    
    v_auto[thisRow,"Est"] <- thisRes$ci_autoX["V_kinshipAuto","Est"]
    v_auto[thisRow,"Lower.95"] <- thisRes$ci_autoX["V_kinshipAuto","Lower.95"]
    v_auto[thisRow,"Upper.95"] <- thisRes$ci_autoX["V_kinshipAuto","Upper.95"]
    
    v_x[thisRow,"Est"] <- thisRes$ci_autoX["V_kinshipX","Est"]
    v_x[thisRow,"Lower.95"] <- thisRes$ci_autoX["V_kinshipX","Lower.95"]
    v_x[thisRow,"Upper.95"] <- thisRes$ci_autoX["V_kinshipX","Upper.95"]
    
    v_e[thisRow,"Est"] <- thisRes$ci_autoX["V_E","Est"]
    v_e[thisRow,"Lower.95"] <- thisRes$ci_autoX["V_E","Lower.95"]
    v_e[thisRow,"Upper.95"] <- thisRes$ci_autoX["V_E","Upper.95"]
    
  }}  # end the loop through the iterations
m <- metrics
hx <- m$beta^2*3*m$p*(1-m$p)+m$sigmaX
newherit <- (m$beta^2*3*m$p*(1-m$p)+m$sigmaX)/(m$beta^2*3*m$p*(1-m$p)+m$sigmaX+m$sigmaE+m$sigmaX)
newbeta <- sqrt(((m$sigmaX+m$sigmaE+m$sigmaA)*newherit-m$sigmaX)/(3*m$p*(1-m$p)*(1-newherit)))

pdf("metrics_24_sigmaX_est.pdf",width=14)
areas <- 1:500
minY <- min(v_x$Lower.95)+0.01; maxY <- max(v_x$Upper.95)+0.01
par(mar=c(4,4,1,1)+0.2)
plot(areas, v_x$Est, type='n', xlab="Iteration",
     ylab="Estimate", cex.lab=1.2,ylim=c(minY,maxY))
polygon(x=c(areas, rev(areas)),
        y=c(v_x$Lower.95, rev(v_x$Upper.95)),
        col='grey90',border=NA)
lines(areas, v_x$Est)
abline(h=hx,col="red",lwd=2,lty=2)
dev.off()

pdf("metrics_24_sigmaA_est.pdf",width=14)
areas <- 1:500
minY <- min(v_auto$Lower.95)+0.01; maxY <- max(v_auto$Upper.95)+0.01
par(mar=c(4,4,1,1)+0.2)
plot(areas, v_auto$Est, type='n', xlab="Iteration",
     ylab="Estimate", cex.lab=1.2,ylim=c(minY,maxY))
polygon(x=c(areas, rev(areas)),
        y=c(v_auto$Lower.95, rev(v_auto$Upper.95)),
        col='grey90',border=NA)
lines(areas, v_auto$Est)
abline(h=m$sigmaA,col="red",lwd=2,lty=2)
dev.off()


##### 
# try again with a different set of metrics
i=1
dat <- get(load("simResults_vars/simRes_varCompsOnly_1_1.RData"))
thisRes <- dat[[i]]
metrics <- thisRes$metricsSet

# want to store thisRes$ci_autoX
v_auto <- data.frame(matrix(NA,nrow=500,ncol=3))
colnames(v_auto) <- c("Est","Lower.95","Upper.95")

v_x <- data.frame(matrix(NA,nrow=500,ncol=3))
colnames(v_x) <- c("Est","Lower.95","Upper.95")

v_e <- data.frame(matrix(NA,nrow=500,ncol=3))
colnames(v_e) <- c("Est","Lower.95","Upper.95")

thisRow <- 0

for(j in 1:50){
  for(k in 1:10){
    dat <- get(load(paste("simResults_vars/simRes_varCompsOnly_",j,"_",k,".RData",sep="")))
    thisRes <- dat[[i]]
    
    thisRow <- thisRow+1
    
    v_auto[thisRow,"Est"] <- thisRes$ci_autoX["V_kinshipAuto","Est"]
    v_auto[thisRow,"Lower.95"] <- thisRes$ci_autoX["V_kinshipAuto","Lower.95"]
    v_auto[thisRow,"Upper.95"] <- thisRes$ci_autoX["V_kinshipAuto","Upper.95"]
    
    v_x[thisRow,"Est"] <- thisRes$ci_autoX["V_kinshipX","Est"]
    v_x[thisRow,"Lower.95"] <- thisRes$ci_autoX["V_kinshipX","Lower.95"]
    v_x[thisRow,"Upper.95"] <- thisRes$ci_autoX["V_kinshipX","Upper.95"]
    
    v_e[thisRow,"Est"] <- thisRes$ci_autoX["V_E","Est"]
    v_e[thisRow,"Lower.95"] <- thisRes$ci_autoX["V_E","Lower.95"]
    v_e[thisRow,"Upper.95"] <- thisRes$ci_autoX["V_E","Upper.95"]
    
  }}  # end the loop through the iterations
m <- metrics
hx <- m$beta^2*3*m$p*(1-m$p)+m$sigmaX
newherit <- (m$beta^2*3*m$p*(1-m$p)+m$sigmaX)/(m$beta^2*3*m$p*(1-m$p)+m$sigmaX+m$sigmaE+m$sigmaX)
newbeta <- sqrt(((m$sigmaX+m$sigmaE+m$sigmaA)*newherit-m$sigmaX)/(3*m$p*(1-m$p)*(1-newherit)))

pdf("metrics_1_sigmaX_est.pdf",width=14)
areas <- 1:500
minY <- min(v_x$Lower.95)+0.01; maxY <- max(v_x$Upper.95)+0.01
par(mar=c(4,4,1,1)+0.2)
plot(areas, v_x$Est, type='n', xlab="Iteration",
     ylab="Estimate", cex.lab=1.2,ylim=c(minY,maxY))
polygon(x=c(areas, rev(areas)),
        y=c(v_x$Lower.95, rev(v_x$Upper.95)),
        col='grey90',border=NA)
lines(areas, v_x$Est)
abline(h=hx,col="red",lwd=2,lty=2)
dev.off()

pdf("metrics_1_sigmaA_est.pdf",width=14)
areas <- 1:500
minY <- min(v_auto$Lower.95)+0.01; maxY <- max(v_auto$Upper.95)+0.01
par(mar=c(4,4,1,1)+0.2)
plot(areas, v_auto$Est, type='n', xlab="Iteration",
     ylab="Estimate", cex.lab=1.2,ylim=c(minY,maxY))
polygon(x=c(areas, rev(areas)),
        y=c(v_auto$Lower.95, rev(v_auto$Upper.95)),
        col='grey90',border=NA)
lines(areas, v_auto$Est)
abline(h=m$sigmaA,col="red",lwd=2,lty=2)
dev.off()

##### 
# try again with a different set of metrics
i=15
dat <- get(load("simResults_vars/simRes_varCompsOnly_1_1.RData"))
thisRes <- dat[[i]]
metrics <- thisRes$metricsSet

# want to store thisRes$ci_autoX
v_auto <- data.frame(matrix(NA,nrow=500,ncol=3))
colnames(v_auto) <- c("Est","Lower.95","Upper.95")

v_x <- data.frame(matrix(NA,nrow=500,ncol=3))
colnames(v_x) <- c("Est","Lower.95","Upper.95")

v_e <- data.frame(matrix(NA,nrow=500,ncol=3))
colnames(v_e) <- c("Est","Lower.95","Upper.95")

thisRow <- 0

for(j in 1:50){
  for(k in 1:10){
    dat <- get(load(paste("simResults_vars/simRes_varCompsOnly_",j,"_",k,".RData",sep="")))
    thisRes <- dat[[i]]
    
    thisRow <- thisRow+1
    
    v_auto[thisRow,"Est"] <- thisRes$ci_autoX["V_kinshipAuto","Est"]
    v_auto[thisRow,"Lower.95"] <- thisRes$ci_autoX["V_kinshipAuto","Lower.95"]
    v_auto[thisRow,"Upper.95"] <- thisRes$ci_autoX["V_kinshipAuto","Upper.95"]
    
    v_x[thisRow,"Est"] <- thisRes$ci_autoX["V_kinshipX","Est"]
    v_x[thisRow,"Lower.95"] <- thisRes$ci_autoX["V_kinshipX","Lower.95"]
    v_x[thisRow,"Upper.95"] <- thisRes$ci_autoX["V_kinshipX","Upper.95"]
    
    v_e[thisRow,"Est"] <- thisRes$ci_autoX["V_E","Est"]
    v_e[thisRow,"Lower.95"] <- thisRes$ci_autoX["V_E","Lower.95"]
    v_e[thisRow,"Upper.95"] <- thisRes$ci_autoX["V_E","Upper.95"]
    
  }}  # end the loop through the iterations
m <- metrics
hx <- m$beta^2*3*m$p*(1-m$p)+m$sigmaX


pdf("metrics_15_sigmaX_est.pdf",width=14)
areas <- 1:500
minY <- min(v_x$Lower.95)+0.01; maxY <- max(v_x$Upper.95)+0.01
par(mar=c(4,4,1,1)+0.2)
plot(areas, v_x$Est, type='n', xlab="Iteration",
     ylab="Estimate", cex.lab=1.2,ylim=c(minY,maxY))
polygon(x=c(areas, rev(areas)),
        y=c(v_x$Lower.95, rev(v_x$Upper.95)),
        col='grey90',border=NA)
lines(areas, v_x$Est)
abline(h=hx,col="red",lwd=2,lty=2)
dev.off()

pdf("metrics_15_sigmaA_est.pdf",width=14)
areas <- 1:500
minY <- min(v_auto$Lower.95)+0.01; maxY <- max(v_auto$Upper.95)+0.01
par(mar=c(4,4,1,1)+0.2)
plot(areas, v_auto$Est, type='n', xlab="Iteration",
     ylab="Estimate", cex.lab=1.2,ylim=c(minY,maxY))
polygon(x=c(areas, rev(areas)),
        y=c(v_auto$Lower.95, rev(v_auto$Upper.95)),
        col='grey90',border=NA)
lines(areas, v_auto$Est)
abline(h=m$sigmaA,col="red",lwd=2,lty=2)
dev.off()

allVals$beta <- NA
for(i in 1:nrow(allVals)){
  allVals$beta[i] <- getEffSize(0.2,allVals$sigmaE[i]+allVals$sigmaA[i]+allVals$sigmaX[i],
                                herit=allVals$herit[i])
}
allVals$true_sigmaX <- allVals$beta^2*3*0.2*(1-0.2)+allVals$sigmaX
allVals$est_sigmaX <- sapply(vals,function(x){as.numeric(x[1,2])})
allVals$est_cbeta12 <- allVals$est_sigmaX-allVals$sigmaX
allVals$est_beta12 <- allVals$est_cbeta12/(3*0.2*0.8)
allVals$beta2 <- allVals$beta^2
allVals$est_beta <- sqrt(allVals$est_beta12)
allVals$beta_totalH <- (allVals$herit*(allVals$sigmaX+allVals$sigmaE+allVals$sigmaA)-allVals$sigmaX)/((1-allVals$herit)*3*0.2*0.8)
allVals$true_sigmaX_totalH <- allVals$sigmaX+3*0.2*0.8*allVals$beta_totalH^2
attach(allVals)

rm(list=ls())


#####
# 4. Process results using 2KC for auto and 2,sqrt(2),1 KC for X chr

# for each metric value, calculate the expected var component
# then compare to estimated var comp

meanMatrix <- matrix(0,nrow=3,ncol=13)
colnames(meanMatrix) <- c("true","all3_mean","all3_sd","justAuto_mean","justAuto_sd","justX_mean","justX_sd",
                          "all3_upperMean","all3_lowerMean","justAuto_upperMean","justAuto_lowerMean",
                          "justX_upperMean","justX_lowerMean")
rownames(meanMatrix) <- c("sigmaX","sigmaA","sigmaE")

meanML <- list(meanMatrix,meanMatrix,meanMatrix,meanMatrix,
               meanMatrix,meanMatrix,meanMatrix,meanMatrix,
               meanMatrix,meanMatrix,meanMatrix,meanMatrix,
               meanMatrix,meanMatrix,meanMatrix,meanMatrix,
               meanMatrix,meanMatrix,meanMatrix,meanMatrix,
               meanMatrix,meanMatrix,meanMatrix,meanMatrix)

dat <- get(load("simResults_vars/simRes_varCompsOnly_2kc_mf_1_1.RData"))

for(i in 1:length(dat)){
  
  sumX <- NULL; sumE <- NULL
  sumA <- NULL; sumEAuto <- NULL
  sumAbo <- NULL; sumXbo <- NULL; sumEbo <- NULL
  
  sumLabo <- NULL; sumUabo <- NULL
  sumLxbo <- NULL; sumUxbo <- NULL
  sumLebo <- NULL; sumUebo <- NULL
  
  sumLa <- NULL; sumUa <- NULL
  sumLeauto <- NULL; sumUeauto <- NULL
  
  sumLx <- NULL; sumUx <- NULL
  sumLex <- NULL; sumUex <- NULL
  
  for(j in 1:50){
    for(k in 1:10){
      dat <- get(load(paste("simResults_vars/simRes_varCompsOnly_2kc_mf_",j,"_",k,".RData",sep="")))
      thisRes <- dat[[i]]
      
      x <- thisRes$ci_inclXonly
      sumX <- c(sumX,x["V_kinshipX","Est"])
      sumE <- c(sumE,x["V_E","Est"])
      sumLx <- c(sumLx,x["V_kinshipX","Lower.95"])
      sumUx <- c(sumUx,x["V_kinshipX","Upper.95"])
      sumLex <- c(sumLx,x["V_E","Lower.95"])
      sumUex <- c(sumUx,x["V_E","Upper.95"])
      
      auto <- thisRes$ci_autoOnly
      sumA <- c(sumA,auto["V_kinshipAuto","Est"])
      sumEAuto <- c(sumEAuto,auto["V_E","Est"])
      sumLa <- c(sumLa,auto["V_kinshipAuto","Lower.95"])
      sumUa <- c(sumUa,auto["V_kinshipAuto","Upper.95"])
      sumLeauto <- c(sumLeauto,auto["V_E","Lower.95"])
      sumUeauto <- c(sumUeauto,auto["V_E","Upper.95"])
      
      bo <- thisRes$ci_autoX
      sumAbo <- c(sumAbo,bo["V_kinshipAuto","Est"])
      sumXbo <- c(sumXbo,bo["V_kinshipX","Est"])
      sumEbo <- c(sumEbo,bo["V_E","Est"])
      
      sumLabo <- c(sumLabo,bo["V_kinshipAuto","Lower.95"])
      sumUabo <- c(sumUabo,bo["V_kinshipAuto","Upper.95"])
      sumLxbo <- c(sumLxbo,bo["V_kinshipX","Lower.95"])
      sumUxbo <- c(sumUxbo,bo["V_kinshipX","Upper.95"])
      sumLebo <- c(sumLebo,bo["V_E","Lower.95"])
      sumUebo <- c(sumUebo,bo["V_E","Upper.95"])
      
    }}  # end the loop through the iterations
  
  m <- thisRes$metricsSet
  hx <- m$beta^2*3*m$p*(1-m$p)+m$sigmaX
  
  meanML[[i]]["sigmaX","true"] <- hx
  meanML[[i]]["sigmaE","true"] <- m$sigmaE
  meanML[[i]]["sigmaA","true"] <- m$sigmaA
  
  meanML[[i]]["sigmaX","all3_mean"] <- mean(sumXbo)
  meanML[[i]]["sigmaX","all3_sd"] <- sd(sumXbo)
  
  meanML[[i]]["sigmaA","all3_mean"] <- mean(sumAbo)
  meanML[[i]]["sigmaA","all3_sd"] <- sd(sumAbo)
  
  meanML[[i]]["sigmaE","all3_mean"] <- mean(sumEbo)
  meanML[[i]]["sigmaE","all3_sd"] <- sd(sumEbo)
  
  meanML[[i]]["sigmaA","justAuto_mean"] <- mean(sumA)
  meanML[[i]]["sigmaA","justAuto_sd"] <- sd(sumA)
  
  meanML[[i]]["sigmaE","justAuto_mean"] <- mean(sumEAuto)
  meanML[[i]]["sigmaE","justAuto_sd"] <- sd(sumEAuto)
  
  meanML[[i]]["sigmaX","justX_mean"] <- mean(sumX)
  meanML[[i]]["sigmaX","justX_sd"] <- sd(sumX)
  
  meanML[[i]]["sigmaE","justX_mean"] <- mean(sumE)
  meanML[[i]]["sigmaE","justX_sd"] <- sd(sumE)
  
  meanML[[i]]["sigmaX","all3_upperMean"] <- paste("(",format(mean(sumLxbo),digits=5),", ",format(mean(sumUxbo),digits=5),")",sep="")
  meanML[[i]]["sigmaX","all3_lowerMean"] <- mean(sumLxbo)
  meanML[[i]]["sigmaA","all3_upperMean"] <- paste("(",format(mean(sumLabo),digits=5),", ",format(mean(sumUabo),digits=5),")",sep="")
  meanML[[i]]["sigmaA","all3_lowerMean"] <- mean(sumLabo)  
  meanML[[i]]["sigmaE","all3_upperMean"] <- paste("(",format(mean(sumLebo),digits=5),", ",format(mean(sumUebo),digits=5),")",sep="")
  meanML[[i]]["sigmaE","all3_lowerMean"] <- mean(sumLebo)
  
  meanML[[i]]["sigmaA","justAuto_upperMean"] <-  paste("(",format(mean(sumLa),digits=5),", ",format(mean(sumUa),digits=5),")",sep="")
  meanML[[i]]["sigmaA","justAuto_lowerMean"] <- mean(sumLa)
  meanML[[i]]["sigmaE","justAuto_upperMean"] <- paste("(",format(mean(sumLeauto),digits=5),", ",format(mean(sumUeauto),digits=5),")",sep="")
  meanML[[i]]["sigmaE","justAuto_lowerMean"] <- mean(sumLeauto)
  
  meanML[[i]]["sigmaX","justX_upperMean"] <- paste("(",format(mean(sumLx),digits=5),", ",format(mean(sumUx),digits=5),")",sep="")
  meanML[[i]]["sigmaX","justX_lowerMean"] <- mean(sumLx)
  meanML[[i]]["sigmaE","justX_upperMean"] <- paste("(",format(mean(sumLex),digits=5),", ",format(mean(sumUex),digits=5),")",sep="")
  meanML[[i]]["sigmaE","justX_lowerMean"] <- mean(sumLex)
  
}

# calculate 95% ci's for the mean/sd
for(i in 1:length(meanML)){
  iter <- meanML[[i]]
  meanML[[i]] <- cbind(meanML[[i]],rep(NA,3))
  colnames(meanML[[i]])[ncol(meanML[[i]])] <- "all3_ci"
  meanML[[i]] <- cbind(meanML[[i]],rep(NA,3))
  colnames(meanML[[i]])[ncol(meanML[[i]])] <- "justAuto_ci"
  meanML[[i]] <- cbind(meanML[[i]],rep(NA,3))
  colnames(meanML[[i]])[ncol(meanML[[i]])] <- "justX_ci"
  
  meanML[[i]][,"all3_ci"] <- paste("(",format(as.numeric(iter[,"all3_mean"])-1.96*as.numeric(iter[,"all3_sd"])/sqrt(500),
                                              digits=4),", ",
                                   format(as.numeric(iter[,"all3_mean"])+1.96*as.numeric(iter[,"all3_sd"])/sqrt(500),
                                          digits=4),")",sep="")
  meanML[[i]][,"justAuto_ci"] <- paste("(",format(as.numeric(iter[,"justAuto_mean"])-1.96*as.numeric(iter[,"justAuto_sd"])/sqrt(500),
                                                  digits=4),", ",
                                       format(as.numeric(iter[,"justAuto_mean"])+1.96*as.numeric(iter[,"justAuto_sd"])/sqrt(500),
                                              digits=4),")",sep="")
  meanML[[i]][,"justX_ci"] <- paste("(",format(as.numeric(iter[,"justX_mean"])-1.96*as.numeric(iter[,"justX_sd"])/sqrt(500),
                                               digits=4),", ",
                                    format(as.numeric(iter[,"justX_mean"])+1.96*as.numeric(iter[,"justX_sd"])/sqrt(500),
                                           digits=4),")",sep="")
}

save(meanML,file="var_comp_estimates_500iter_2kc_mf.RData")

library(xtable)
for(i in 1:length(meanML)){
  meanML[[i]][,"true"] <- format(as.numeric(meanML[[i]][,"true"]),digits=4)
  meanML[[i]][,"all3_mean"] <- format(as.numeric(meanML[[i]][,"all3_mean"]),digits=4)
  meanML[[i]][,"justX_mean"] <- format(as.numeric(meanML[[i]][,"justX_mean"]),digits=4)
  meanML[[i]][,"justAuto_mean"] <- format(as.numeric(meanML[[i]][,"justAuto_mean"]),digits=4)
  iter <- meanML[[i]]
  cisM1 <- paste(iter[,c("all3_mean")],iter[,c("all3_upperMean")])
  cisM2 <- paste(iter[,"justX_mean"],iter[,"justX_upperMean"])
  cisM3 <- paste(iter[,"justAuto_mean"],iter[,"justAuto_upperMean"])
  print(xtable(matrix(c(iter[,c("true")],cisM1,cisM2,cisM3),nrow=3)))
}

rm(list=ls())

# 4b. Make plots of individual metrics and estimate for 500 iters

##### 
# try again with a different set of metrics
i=24
dat <- get(load("simResults_vars/simRes_varCompsOnly_2kc_mf_1_1.RData"))
thisRes <- dat[[i]]
metrics <- thisRes$metricsSet

# want to store thisRes$ci_autoX
v_auto <- data.frame(matrix(NA,nrow=500,ncol=3))
colnames(v_auto) <- c("Est","Lower.95","Upper.95")

v_x <- data.frame(matrix(NA,nrow=500,ncol=3))
colnames(v_x) <- c("Est","Lower.95","Upper.95")

v_e <- data.frame(matrix(NA,nrow=500,ncol=3))
colnames(v_e) <- c("Est","Lower.95","Upper.95")

thisRow <- 0

for(j in 1:50){
  for(k in 1:10){
    dat <- get(load(paste("simResults_vars/simRes_varCompsOnly_2kc_mf_",j,"_",k,".RData",sep="")))
    thisRes <- dat[[i]]
    
    thisRow <- thisRow+1
    
    v_auto[thisRow,"Est"] <- thisRes$ci_autoX["V_kinshipAuto","Est"]
    v_auto[thisRow,"Lower.95"] <- thisRes$ci_autoX["V_kinshipAuto","Lower.95"]
    v_auto[thisRow,"Upper.95"] <- thisRes$ci_autoX["V_kinshipAuto","Upper.95"]
    
    v_x[thisRow,"Est"] <- thisRes$ci_autoX["V_kinshipX","Est"]
    v_x[thisRow,"Lower.95"] <- thisRes$ci_autoX["V_kinshipX","Lower.95"]
    v_x[thisRow,"Upper.95"] <- thisRes$ci_autoX["V_kinshipX","Upper.95"]
    
    v_e[thisRow,"Est"] <- thisRes$ci_autoX["V_E","Est"]
    v_e[thisRow,"Lower.95"] <- thisRes$ci_autoX["V_E","Lower.95"]
    v_e[thisRow,"Upper.95"] <- thisRes$ci_autoX["V_E","Upper.95"]
    
  }}  # end the loop through the iterations
m <- metrics
hx <- m$beta^2*3*m$p*(1-m$p)+m$sigmaX
newherit <- (m$beta^2*3*m$p*(1-m$p)+m$sigmaX)/(m$beta^2*3*m$p*(1-m$p)+m$sigmaX+m$sigmaE+m$sigmaX)
newbeta <- sqrt(((m$sigmaX+m$sigmaE+m$sigmaA)*newherit-m$sigmaX)/(3*m$p*(1-m$p)*(1-newherit)))

pdf("metrics_2kc_mf_24_sigmaX_est.pdf",width=14)
areas <- 1:500
minY <- min(v_x$Lower.95)+0.01; maxY <- max(v_x$Upper.95)+0.01
par(mar=c(4,4,1,1)+0.2)
plot(areas, v_x$Est, type='n', xlab="Iteration",
     ylab="Estimate", cex.lab=1.2,ylim=c(minY,maxY))
polygon(x=c(areas, rev(areas)),
        y=c(v_x$Lower.95, rev(v_x$Upper.95)),
        col='grey90',border=NA)
lines(areas, v_x$Est)
abline(h=hx,col="red",lwd=2,lty=2)
dev.off()

pdf("metrics_2kc_mf_24_sigmaA_est.pdf",width=14)
areas <- 1:500
minY <- min(v_auto$Lower.95)+0.01; maxY <- max(v_auto$Upper.95)+0.01
par(mar=c(4,4,1,1)+0.2)
plot(areas, v_auto$Est, type='n', xlab="Iteration",
     ylab="Estimate", cex.lab=1.2,ylim=c(minY,maxY))
polygon(x=c(areas, rev(areas)),
        y=c(v_auto$Lower.95, rev(v_auto$Upper.95)),
        col='grey90',border=NA)
lines(areas, v_auto$Est)
abline(h=m$sigmaA,col="red",lwd=2,lty=2)
dev.off()

# try again with a different set of metrics
i=17
dat <- get(load("simResults_vars/simRes_varCompsOnly_2kc_mf_1_1.RData"))
thisRes <- dat[[i]]
metrics <- thisRes$metricsSet

# want to store thisRes$ci_autoX
v_auto <- data.frame(matrix(NA,nrow=500,ncol=3))
colnames(v_auto) <- c("Est","Lower.95","Upper.95")

v_x <- data.frame(matrix(NA,nrow=500,ncol=3))
colnames(v_x) <- c("Est","Lower.95","Upper.95")

v_e <- data.frame(matrix(NA,nrow=500,ncol=3))
colnames(v_e) <- c("Est","Lower.95","Upper.95")

thisRow <- 0

for(j in 1:50){
  for(k in 1:10){
    dat <- get(load(paste("simResults_vars/simRes_varCompsOnly_2kc_mf_",j,"_",k,".RData",sep="")))
    thisRes <- dat[[i]]
    
    thisRow <- thisRow+1
    
    v_auto[thisRow,"Est"] <- thisRes$ci_autoX["V_kinshipAuto","Est"]
    v_auto[thisRow,"Lower.95"] <- thisRes$ci_autoX["V_kinshipAuto","Lower.95"]
    v_auto[thisRow,"Upper.95"] <- thisRes$ci_autoX["V_kinshipAuto","Upper.95"]
    
    v_x[thisRow,"Est"] <- thisRes$ci_autoX["V_kinshipX","Est"]
    v_x[thisRow,"Lower.95"] <- thisRes$ci_autoX["V_kinshipX","Lower.95"]
    v_x[thisRow,"Upper.95"] <- thisRes$ci_autoX["V_kinshipX","Upper.95"]
    
    v_e[thisRow,"Est"] <- thisRes$ci_autoX["V_E","Est"]
    v_e[thisRow,"Lower.95"] <- thisRes$ci_autoX["V_E","Lower.95"]
    v_e[thisRow,"Upper.95"] <- thisRes$ci_autoX["V_E","Upper.95"]
    
  }}  # end the loop through the iterations
m <- metrics
hx <- m$beta^2*3*m$p*(1-m$p)+m$sigmaX
newherit <- (m$beta^2*3*m$p*(1-m$p)+m$sigmaX)/(m$beta^2*3*m$p*(1-m$p)+m$sigmaX+m$sigmaE+m$sigmaX)
newbeta <- sqrt(((m$sigmaX+m$sigmaE+m$sigmaA)*newherit-m$sigmaX)/(3*m$p*(1-m$p)*(1-newherit)))

pdf("metrics_2kc_mf_17_sigmaX_est.pdf",width=14)
areas <- 1:500
minY <- min(v_x$Lower.95)+0.01; maxY <- max(v_x$Upper.95)+0.01
par(mar=c(4,4,1,1)+0.2)
plot(areas, v_x$Est, type='n', xlab="Iteration",
     ylab="Estimate", cex.lab=1.2,ylim=c(minY,maxY))
polygon(x=c(areas, rev(areas)),
        y=c(v_x$Lower.95, rev(v_x$Upper.95)),
        col='grey90',border=NA)
lines(areas, v_x$Est)
abline(h=hx,col="red",lwd=2,lty=2)
dev.off()

pdf("metrics_2kc_mf_17_sigmaA_est.pdf",width=14)
areas <- 1:500
minY <- min(v_auto$Lower.95)+0.01; maxY <- max(v_auto$Upper.95)+0.01
par(mar=c(4,4,1,1)+0.2)
plot(areas, v_auto$Est, type='n', xlab="Iteration",
     ylab="Estimate", cex.lab=1.2,ylim=c(minY,maxY))
polygon(x=c(areas, rev(areas)),
        y=c(v_auto$Lower.95, rev(v_auto$Upper.95)),
        col='grey90',border=NA)
lines(areas, v_auto$Est)
abline(h=m$sigmaA,col="red",lwd=2,lty=2)
dev.off()

rm(list=ls())


#####
# 5. Create table of metrics for the 24 iterations

dat <- get(load("simResults_vars/simRes_varCompsOnly_2kc_mf_1_1.RData"))
m <- lapply(dat,function(x){x$metricsSet})     
metrics <- data.frame(matrix(unlist(m),nrow=24,ncol=6,byrow=TRUE))
colnames(metrics) <- colnames(dat[[1]]$metricsSet)
metrics$sigmaXT <- metrics$beta^2*3*metrics$p*(1-metrics$p)+metrics$sigmaX

xtable(metrics[,c(1,2,3,5,7,4,6)],digits=c(0,4,2,1,1,4,1,0))

rm(list=ls())


#####
# 6. Parse the 2KC for both auto and X chr results

# for each metric value, calculate the expected var component
# then compare to estimated var comp

meanMatrix <- matrix(0,nrow=3,ncol=13)
colnames(meanMatrix) <- c("true","all3_mean","all3_sd","justAuto_mean","justAuto_sd","justX_mean","justX_sd",
                          "all3_upperMean","all3_lowerMean","justAuto_upperMean","justAuto_lowerMean",
                          "justX_upperMean","justX_lowerMean")
rownames(meanMatrix) <- c("sigmaX","sigmaA","sigmaE")

meanML <- list(meanMatrix,meanMatrix,meanMatrix,meanMatrix,
               meanMatrix,meanMatrix,meanMatrix,meanMatrix,
               meanMatrix,meanMatrix,meanMatrix,meanMatrix,
               meanMatrix,meanMatrix,meanMatrix,meanMatrix,
               meanMatrix,meanMatrix,meanMatrix,meanMatrix,
               meanMatrix,meanMatrix,meanMatrix,meanMatrix)

dat <- get(load("simResults_vars/simRes_varCompsOnly_2kc_1_1.RData"))

for(i in 1:length(dat)){
  
  sumX <- NULL; sumE <- NULL
  sumA <- NULL; sumEAuto <- NULL
  sumAbo <- NULL; sumXbo <- NULL; sumEbo <- NULL
  
  sumLabo <- NULL; sumUabo <- NULL
  sumLxbo <- NULL; sumUxbo <- NULL
  sumLebo <- NULL; sumUebo <- NULL
  
  sumLa <- NULL; sumUa <- NULL
  sumLeauto <- NULL; sumUeauto <- NULL
  
  sumLx <- NULL; sumUx <- NULL
  sumLex <- NULL; sumUex <- NULL
  
  for(j in 1:50){
    for(k in 1:10){
      dat <- get(load(paste("simResults_vars/simRes_varCompsOnly_2kc_",j,"_",k,".RData",sep="")))
      thisRes <- dat[[i]]
      
      x <- thisRes$ci_inclXonly
      sumX <- c(sumX,x["V_kinshipX","Est"])
      sumE <- c(sumE,x["V_E","Est"])
      sumLx <- c(sumLx,x["V_kinshipX","Lower.95"])
      sumUx <- c(sumUx,x["V_kinshipX","Upper.95"])
      sumLex <- c(sumLx,x["V_E","Lower.95"])
      sumUex <- c(sumUx,x["V_E","Upper.95"])
      
      auto <- thisRes$ci_autoOnly
      sumA <- c(sumA,auto["V_kinshipAuto","Est"])
      sumEAuto <- c(sumEAuto,auto["V_E","Est"])
      sumLa <- c(sumLa,auto["V_kinshipAuto","Lower.95"])
      sumUa <- c(sumUa,auto["V_kinshipAuto","Upper.95"])
      sumLeauto <- c(sumLeauto,auto["V_E","Lower.95"])
      sumUeauto <- c(sumUeauto,auto["V_E","Upper.95"])
      
      bo <- thisRes$ci_autoX
      sumAbo <- c(sumAbo,bo["V_kinshipAuto","Est"])
      sumXbo <- c(sumXbo,bo["V_kinshipX","Est"])
      sumEbo <- c(sumEbo,bo["V_E","Est"])
      
      sumLabo <- c(sumLabo,bo["V_kinshipAuto","Lower.95"])
      sumUabo <- c(sumUabo,bo["V_kinshipAuto","Upper.95"])
      sumLxbo <- c(sumLxbo,bo["V_kinshipX","Lower.95"])
      sumUxbo <- c(sumUxbo,bo["V_kinshipX","Upper.95"])
      sumLebo <- c(sumLebo,bo["V_E","Lower.95"])
      sumUebo <- c(sumUebo,bo["V_E","Upper.95"])
      
    }}  # end the loop through the iterations
  
  m <- thisRes$metricsSet
  hx <- m$beta^2*3*m$p*(1-m$p)+m$sigmaX
  
  meanML[[i]]["sigmaX","true"] <- hx
  meanML[[i]]["sigmaE","true"] <- m$sigmaE
  meanML[[i]]["sigmaA","true"] <- m$sigmaA
  
  meanML[[i]]["sigmaX","all3_mean"] <- mean(sumXbo)
  meanML[[i]]["sigmaX","all3_sd"] <- sd(sumXbo)
  
  meanML[[i]]["sigmaA","all3_mean"] <- mean(sumAbo)
  meanML[[i]]["sigmaA","all3_sd"] <- sd(sumAbo)
  
  meanML[[i]]["sigmaE","all3_mean"] <- mean(sumEbo)
  meanML[[i]]["sigmaE","all3_sd"] <- sd(sumEbo)
  
  meanML[[i]]["sigmaA","justAuto_mean"] <- mean(sumA)
  meanML[[i]]["sigmaA","justAuto_sd"] <- sd(sumA)
  
  meanML[[i]]["sigmaE","justAuto_mean"] <- mean(sumEAuto)
  meanML[[i]]["sigmaE","justAuto_sd"] <- sd(sumEAuto)
  
  meanML[[i]]["sigmaX","justX_mean"] <- mean(sumX)
  meanML[[i]]["sigmaX","justX_sd"] <- sd(sumX)
  
  meanML[[i]]["sigmaE","justX_mean"] <- mean(sumE)
  meanML[[i]]["sigmaE","justX_sd"] <- sd(sumE)
  
  meanML[[i]]["sigmaX","all3_upperMean"] <- paste("(",format(mean(sumLxbo),digits=5),", ",format(mean(sumUxbo),digits=5),")",sep="")
  meanML[[i]]["sigmaX","all3_lowerMean"] <- mean(sumLxbo)
  meanML[[i]]["sigmaA","all3_upperMean"] <- paste("(",format(mean(sumLabo),digits=5),", ",format(mean(sumUabo),digits=5),")",sep="")
  meanML[[i]]["sigmaA","all3_lowerMean"] <- mean(sumLabo)  
  meanML[[i]]["sigmaE","all3_upperMean"] <- paste("(",format(mean(sumLebo),digits=5),", ",format(mean(sumUebo),digits=5),")",sep="")
  meanML[[i]]["sigmaE","all3_lowerMean"] <- mean(sumLebo)
  
  meanML[[i]]["sigmaA","justAuto_upperMean"] <-  paste("(",format(mean(sumLa),digits=5),", ",format(mean(sumUa),digits=5),")",sep="")
  meanML[[i]]["sigmaA","justAuto_lowerMean"] <- mean(sumLa)
  meanML[[i]]["sigmaE","justAuto_upperMean"] <- paste("(",format(mean(sumLeauto),digits=5),", ",format(mean(sumUeauto),digits=5),")",sep="")
  meanML[[i]]["sigmaE","justAuto_lowerMean"] <- mean(sumLeauto)
  
  meanML[[i]]["sigmaX","justX_upperMean"] <- paste("(",format(mean(sumLx),digits=5),", ",format(mean(sumUx),digits=5),")",sep="")
  meanML[[i]]["sigmaX","justX_lowerMean"] <- mean(sumLx)
  meanML[[i]]["sigmaE","justX_upperMean"] <- paste("(",format(mean(sumLex),digits=5),", ",format(mean(sumUex),digits=5),")",sep="")
  meanML[[i]]["sigmaE","justX_lowerMean"] <- mean(sumLex)
  
}

# calculate 95% ci's for the mean/sd
for(i in 1:length(meanML)){
  iter <- meanML[[i]]
  meanML[[i]] <- cbind(meanML[[i]],rep(NA,3))
  colnames(meanML[[i]])[ncol(meanML[[i]])] <- "all3_ci"
  meanML[[i]] <- cbind(meanML[[i]],rep(NA,3))
  colnames(meanML[[i]])[ncol(meanML[[i]])] <- "justAuto_ci"
  meanML[[i]] <- cbind(meanML[[i]],rep(NA,3))
  colnames(meanML[[i]])[ncol(meanML[[i]])] <- "justX_ci"
  
  meanML[[i]][,"all3_ci"] <- paste("(",format(as.numeric(iter[,"all3_mean"])-1.96*as.numeric(iter[,"all3_sd"])/sqrt(500),
                                              digits=4),", ",
                                   format(as.numeric(iter[,"all3_mean"])+1.96*as.numeric(iter[,"all3_sd"])/sqrt(500),
                                          digits=4),")",sep="")
  meanML[[i]][,"justAuto_ci"] <- paste("(",format(as.numeric(iter[,"justAuto_mean"])-1.96*as.numeric(iter[,"justAuto_sd"])/sqrt(500),
                                                  digits=4),", ",
                                       format(as.numeric(iter[,"justAuto_mean"])+1.96*as.numeric(iter[,"justAuto_sd"])/sqrt(500),
                                              digits=4),")",sep="")
  meanML[[i]][,"justX_ci"] <- paste("(",format(as.numeric(iter[,"justX_mean"])-1.96*as.numeric(iter[,"justX_sd"])/sqrt(500),
                                               digits=4),", ",
                                    format(as.numeric(iter[,"justX_mean"])+1.96*as.numeric(iter[,"justX_sd"])/sqrt(500),
                                           digits=4),")",sep="")
}

save(meanML,file="var_comp_estimates_500iter_2kc.RData")

# yep, not good at all...we are needing the 2, sqrt(2) and 1 KC matrix

rm(list=ls())


#####
# 7. Look at type I error rate with mf KC matrices
# process results for different param values

library(BiocParallel)

dat <- get(load("simResults_vars/simRes_mf_1.RData"))
length(dat) # 24; these are the different param combos
# the _X.RData in the object name is the iteration

metrics <- unlist(bplapply(dat,function(x){x$metrics}))
metrics <- data.frame(t(matrix(metrics,ncol=24)))
colnames(metrics) <- c("beta1","herit2_snp","p","sigmaA","sigmaX","sigmaE")

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
    dat <- get(load(paste("simResults_vars/simRes_mf_",j,".RData",sep="")))
    
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

write.table(metrics,file="metrics_vars_mf_typeIresults.txt",sep="\t")

rm(list=ls())

metrics <- read.table("metrics_vars_mf_typeIresults.txt",sep="\t")
colSums(metrics)/sum(metrics$nSim)

# make a table of these results for latex
library(xtable)
ty1 <- colSums(metrics)/sum(metrics$nSim)
xtable(t(matrix(c(ty1[12:16],ty1[7:11],ty1[17:21]),nrow=5)),digits=5)

# calculate the herit for all x chromosome SNPs
attach(metrics)
herit <- (beta1^2*3*p*(1-p)+sigmaX)/(beta1^2*3*p*(1-p)+sigmaX+sigmaE+sigmaA)
detach(metrics)

metrics$herit_xchr <- herit

xtable(metrics[,c(23,1,2,4,5,7,12,17,8,13,18,9,14,19,10,15,20,
                  11,16,21)],digits=c(1,4,4,3,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))

tabLat <- metrics[,c(1,2,4,5,7,12,17,8,13,18,9,14,19,10,15,20,
                     11,16,21)]
xtable(matrix(colSums(tabLat),nrow=1),digits=0) # totals row
xtable(matrix(colSums(tabLat)/sum(metrics$nSim),nrow=1),digits=5) # percentage row


# a smaller table for the type I error results
xtable(as.matrix(cbind(ty1[c(7,12,17)],ty1[c(8,13,18)],ty1[c(9,14,19)],
                       ty1[c(10,15,20)],ty1[c(11,16,21)])),digits=5)

rm(list=ls())


#####
# 8. Get results for simple model with 10K pedigrees

dat <- get(load("500iter_bySex_both_varComps_X.RData"))
dim(dat) # 500 9
paste(format(mean(dat[,1]),digits=4), " (", format(mean(dat[,2]),digits=4), ", ",
      format(mean(dat[,3]),digits=4), ")",sep="")

dat <- get(load("500iter_bySex_both_varComps_E.RData"))
dim(dat) # 500 9
paste(format(mean(dat[,1]),digits=4), " (", format(mean(dat[,2]),digits=4), ", ",
      format(mean(dat[,3]),digits=4), ")",sep="")

## next model
dat <- get(load("500iter_bySex_both_gx_varComps_X.RData"))
dim(dat) # 500 9
paste(format(mean(dat[,1]),digits=4), " (", format(mean(dat[,2]),digits=4), ", ",
      format(mean(dat[,3]),digits=4), ")",sep="")

dat <- get(load("500iter_bySex_both_gx_varComps_E.RData"))
dim(dat) # 500 9
paste(format(mean(dat[,1]),digits=4), " (", format(mean(dat[,2]),digits=4), ", ",
      format(mean(dat[,3]),digits=4), ")",sep="")

## next model
dat <- get(load("500iter_bySex_both_gx_ga_varComps_X.RData"))
dim(dat) # 500 9
paste(format(mean(dat[,1],na.rm=TRUE),digits=4), " (", format(mean(dat[,2],na.rm=TRUE),digits=4), ", ",
      format(mean(dat[,3],na.rm=TRUE),digits=4), ")",sep="")

dat <- get(load("500iter_bySex_both_gx_ga_varComps_A.RData"))
dim(dat) # 500 9
paste(format(mean(dat[,1],na.rm=TRUE),digits=4), " (", format(mean(dat[,2],na.rm=TRUE),digits=4), ", ",
      format(mean(dat[,3],na.rm=TRUE),digits=4), ")",sep="")

dat <- get(load("500iter_bySex_both_gx_ga_varComps_E.RData"))
dim(dat) # 500 9
paste(format(mean(dat[,1],na.rm=TRUE),digits=4), " (", format(mean(dat[,2],na.rm=TRUE),digits=4), ", ",
      format(mean(dat[,3],na.rm=TRUE),digits=4), ")",sep="")

rm(list=ls())


#####
# 9. Get results for more iterations of model A: y=beta1*SNP_x + e

# (do this on the servers)

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x/")

dat <- get(load("/projects/geneva/geneva_sata/caitlin/mlm_x/modelA_results/mmRes_varComps_modelA_SNP1.RData"))
length(dat); names(dat) # 6
# it holds var comps, pvalues for the mixed model
# for 3 models, the "x" version of which is true
# the true var comp for x is 0.4096 
# the true var comp for e is 1

head(dat[[4]])

# look at the true model and see how many true positives are found
# how many false positives are found
# loop through the 500 .RData files in the folder
# the first SNP in the dat[[4]] slot is the truly associated SNP
# the remaining 500 SNPs in the dat[[4]] slot are not associated with the phenotype
true_ps <- rep(NA,500)
for(i in 1:500){
  thisRes <- get(load(paste("modelA_results/mmRes_varComps_modelA_SNP",i,".RData",sep="")))
  pv <- thisRes[[4]]
  true_ps[i] <- pv$pval[1]
}
true_ps # all are 0, well, great, but boring

# check the false positives
false_ps <- rep(NA,1000)
i=1 # this will hold results for SNPs 2-501
thisRes <- get(load(paste("modelA_results/mmRes_varComps_modelA_SNP",i,".RData",sep="")))
pv <- thisRes[[4]]
false_ps[2:501] <- pv$pval[2:501]

i=500 # this will hold results for SNPs 501-1000
thisRes <- get(load(paste("modelA_results/mmRes_varComps_modelA_SNP",i,".RData",sep="")))
pv <- thisRes[[4]]
false_ps[501:1000] <- pv$pval[pv$snpID<=1000&pv$snpID>500]

sum(is.na(false_ps)) # 1; since we never tested the first SNP. it's ok
false_ps <- false_ps[-1]
length(false_ps) # 999

sum(false_ps<0.001) # 1
sum(false_ps<0.01) # 7
sum(false_ps<0.005) # 3
sum(false_ps<0.0001) # 0

#####
# ok, look at these results for when we test the autosomal effects only, but actu it's x chr effects
# these results are stored in the 6th slot of the results files

true_ps <- rep(NA,500)
for(i in 1:500){
  thisRes <- get(load(paste("modelA_results/mmRes_varComps_modelA_SNP",i,".RData",sep="")))
  pv <- thisRes[[6]]
  true_ps[i] <- pv$pval[1]
}
true_ps # all are 0, well, great, but boring

# check the false positives
false_ps <- rep(NA,1000)
i=1 # this will hold results for SNPs 2-501
thisRes <- get(load(paste("modelA_results/mmRes_varComps_modelA_SNP",i,".RData",sep="")))
pv <- thisRes[[6]]
false_ps[2:501] <- pv$pval[2:501]

i=500 # this will hold results for SNPs 501-1000
thisRes <- get(load(paste("modelA_results/mmRes_varComps_modelA_SNP",i,".RData",sep="")))
pv <- thisRes[[6]]
false_ps[501:1000] <- pv$pval[pv$snpID<=1000&pv$snpID>500]

sum(is.na(false_ps)) # 1; since we never tested the first SNP. it's ok
false_ps <- false_ps[-1]
length(false_ps) # 999

sum(false_ps<0.001) # 1
sum(false_ps<0.01) # 11
sum(false_ps<0.005) # 5
sum(false_ps<0.0001) # 0
# so higher type I error than the other method. great!

###
# now look at the model c results, when fit with just autosomal
# this is the realistic case when we think x chr effects are there but we just fit the autosomal ones

true_ps <- rep(NA,500)
for(i in 1:500){
  thisRes <- get(load(paste("modelC_results/mmRes_varComps_modelA_SNP",i,".RData",sep="")))
  pv <- thisRes[[6]]
  true_ps[i] <- pv$pval[1]
}
true_ps # all are 0, well, great, but boring
# not identically zero
summary(true_ps) # basically zero!
#Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
#0.000e+00  0.000e+00  0.000e+00 1.591e-204  0.000e+00 7.954e-202 

# check the false positives
false_ps <- rep(NA,1000)
i=1 # this will hold results for SNPs 2-501
thisRes <- get(load(paste("modelA_results/mmRes_varComps_modelA_SNP",i,".RData",sep="")))
pv <- thisRes[[6]]
false_ps[2:501] <- pv$pval[2:501]

i=500 # this will hold results for SNPs 501-1000
thisRes <- get(load(paste("modelA_results/mmRes_varComps_modelA_SNP",i,".RData",sep="")))
pv <- thisRes[[6]]
false_ps[501:1000] <- pv$pval[pv$snpID<=1000&pv$snpID>500]

sum(is.na(false_ps)) # 1; since we never tested the first SNP. it's ok
false_ps <- false_ps[-1]
length(false_ps) # 999

sum(false_ps<0.001) # 1
sum(false_ps<0.01) # 11
sum(false_ps<0.05)
sum(false_ps<0.005) # 5
sum(false_ps<0.0001) # 0

## what happens if we fit the true model?
true_ps <- rep(NA,500)
for(i in 1:500){
  thisRes <- get(load(paste("modelC_results/mmRes_varComps_modelA_SNP",i,".RData",sep="")))
  pv <- thisRes[[2]]
  true_ps[i] <- pv$pval[1]
}
true_ps # all are 0, well, great, but boring
# not identically zero
summary(true_ps) # basically zero!
#Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
#0.000e+00  0.000e+00  0.000e+00 1.263e-174  0.000e+00 3.785e-172 

# check the false positives
false_ps <- rep(NA,1000)
i=1 # this will hold results for SNPs 2-501
thisRes <- get(load(paste("modelA_results/mmRes_varComps_modelA_SNP",i,".RData",sep="")))
pv <- thisRes[[2]]
false_ps[2:501] <- pv$pval[2:501]

i=500 # this will hold results for SNPs 501-1000
thisRes <- get(load(paste("modelA_results/mmRes_varComps_modelA_SNP",i,".RData",sep="")))
pv <- thisRes[[2]]
false_ps[501:1000] <- pv$pval[pv$snpID<=1000&pv$snpID>500]

sum(is.na(false_ps)) # 1; since we never tested the first SNP. it's ok
false_ps <- false_ps[-1]
length(false_ps) # 999

sum(false_ps<0.001) # 1
sum(false_ps<0.05) # 49
sum(false_ps<0.01) # 6
sum(false_ps<0.005) # 3
sum(false_ps<0.0001) # 0
# ok, so lower type I error rate than the other model fit

rm(list=ls())


#####
# 10. More simulations

## called these functions to simulate
# true phenotype: y=beta1*SNPx + gx + ga + e
# sigmaA=0.3
# sigmaX=0.8
# beta1=0.8

# test true model
# test y=beta1*SNPx + gx + e
# test y=beta1*SNPx + ga + e

# to run, call 
# cd /projects/geneva/geneva_sata/caitlin/mlm_x/
# qsub batch_model_c_sims.sh
# which calls model_batch_sims.R
# which calls model_batch_sims_function.R

# and saves model_c_batch_sims.Rout
# and saves results in modelC_results/mmRes_varComps_modelA_SNP*.RData


## called these functions to simulate
# true phenotype: y=beta1*SNPx + gx + ga + e
# sigmaA=0.3
# sigmaX=0.8
# beta1=0.3

# test true model
# test y=beta1*SNPx + gx + e
# test y=beta1*SNPx + ga + e

# to run, call 
# cd /projects/geneva/geneva_sata/caitlin/mlm_x/
# qsub batch_model_c_sims_smBeta_1.sh
# which calls model_batch_sims_smBeta.R
# which calls model_batch_sims_smBeta_modelC_function.R

# and saves model_c_batch_sims_smBeta_1.Rout
# and saves results in modelC_results_smBeta/mmRes_varComps_modelC_SNP*.RData


## called these functions to simulate
# true phenotype: y=beta1*SNPx + gx + ga + e
# sigmaA=0.3
# sigmaX=0.3
# beta1=0.3

# test true model
# test y=beta1*SNPx + gx + e
# test y=beta1*SNPx + ga + e

# to run, call 
# cd /projects/geneva/geneva_sata/caitlin/mlm_x/
# qsub batch_model_c_sims_smBeta_2.sh
# which calls model_batch_sims_smBeta2.R
# which calls model_batch_sims_smBeta_modelC_function2.R

# and saves model_c_batch_sims_smBeta_2.Rout
# and saves results in modelC_results_smBeta2/mmRes_varComps_modelC_SNP*.RData


#####
# 11. Type I error 

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x/")
source("typeIErr.R")

nSets <- 23
alpha <- 0.005

(res <- typeIErr(nSets,alpha))
res/res[1]
# totalIterations          false_pos_trueModel 
# 1.000000000                  0.005391304 
# false_pos_falseModel  false_pos_trueModel_smBeta2 
# 0.011565217                  0.005217391 
# false_pos_falseModel_smBeta2 
# 0.007739130 

rm(list=ls())


#####
# 12. Power calcs: check phenotype changes with different beta1 values

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")

library(MASS); library(parallel)
library(GWASTools)
library(corpcor); library(doSNOW)
source("sim_phenotype.R")
source("estVarComp.R")
source("allele_drop_functions.R")
source("assocTestMixedModel_v7.R")

##
# load in all the supporting docs

n <- 8000

TEMP=c(2,1,2,1,1,1,2,2,1,2,1,2,1,2,1,2)
SEX <- rep(NA,length(TEMP))
SEX[TEMP==1]="M"
SEX[TEMP==2]="F"

sex <- rep(SEX,(n/16))

kinX <- get(load("500Peds_xKinship.RData"))
kinAuto <- get(load("500Peds_autoKinship.RData"))

kinAuto <- matrix(kinAuto,nrow=n,ncol=n)
kinX <- matrix(kinX,nrow=n,ncol=n)

dat <- get(load("tmp_scanAnnot_pheno_3models_genotypes_p2_8000samples_smBeta.RData"))
geno <- get(load("genotypes_p02_8000samples.RData"))
dim(geno) # 1000 8000

beta1 <- 0.3
pheno1 <- dat[,3]
recov1 <- pheno1-beta1*geno[1,]

# get the 2 new phenotype variables
pheno2 <- recov1+0.1*geno[1,]  
pheno3 <- recov1+0.2*geno[1,]

# now take the new phenotypes, est var comp and get the assoc test results
scan <- data.frame(scanID=1:n,sex=sex,pheno=pheno2)
scan <- ScanAnnotationDataFrame(scan)

genoMt <- MatrixGenotypeReader(genotype=geno,snpID=as.integer(1:nrow(geno)),
                               chromosome=as.integer(rep(23,nrow(geno))),
                               position=as.integer(1:nrow(geno)), scanID=scan$scanID)
genoData <- GenotypeData(genoMt,scanAnnot=scan)

colnames(kinX) <- scan$scanID
colnames(kinAuto) <- scan$scanID

covMatListBoth <- list(kinX,kinAuto) # this is the X chr kinship matrix for all samples
names(covMatListBoth) <- c("kinshipX","kinshipAutos")

varCompBoth <- estVarComp(scan,covMatList=covMatListBoth,"pheno")
varBothci <- estVarCompCI(varCompBoth,prop=FALSE)
#                     Est   Lower 95  Upper 95
#V_kinshipX     0.7525647 0.63094920 0.8741801
#V_kinshipAutos 0.2638003 0.06181687 0.4657838
#V_E            1.0050004 0.94452549 1.0654752

cholSig <- varCompBoth[["cholSigmaInv"]]
mmResBoth <- assocTestMixedModel(genoData,snpStart=1,snpEnd=2,cholSig,outcome="pheno")
#snpID chr    n       MAF minor.allele         Est         SE      Stat
#1     1  23 8000 0.2000000            A  0.07777293 0.02435305 10.198827
#2     2  23 8000 0.2018333            A -0.03059302 0.02428695  1.586714
#pval
#1 0.0014053
#2 0.2077962

# do with the other phenotype
scan <- data.frame(scanID=1:n,sex=sex,pheno=pheno3)
scan <- ScanAnnotationDataFrame(scan)

genoMt <- MatrixGenotypeReader(genotype=geno,snpID=as.integer(1:nrow(geno)),
                               chromosome=as.integer(rep(23,nrow(geno))),
                               position=as.integer(1:nrow(geno)), scanID=scan$scanID)
genoData <- GenotypeData(genoMt,scanAnnot=scan)

varCompBoth3 <- estVarComp(scan,covMatList=covMatListBoth,"pheno")
varBothci3 <- estVarCompCI(varCompBoth3,prop=FALSE)
#                     Est   Lower 95  Upper 95
#V_kinshipX     0.7704196 0.64763614 0.8932031
#V_kinshipAutos 0.2671046 0.06332076 0.4708884
#V_E            1.0030805 0.94236168 1.0637993

cholSig3 <- varCompBoth3[["cholSigmaInv"]]
mmResBoth3 <- assocTestMixedModel(genoData,snpStart=1,snpEnd=2,cholSig3,outcome="pheno")
#snpID chr    n       MAF minor.allele         Est         SE      Stat
#1     1  23 8000 0.2000000            A  0.17764338 0.02440689 52.975286
#2     2  23 8000 0.2018333            A -0.03244745 0.02440208  1.768103
#pval
#1 3.377717e-13
#2 1.836171e-01

rm(list=ls())


#####
# 13. Power simulations

# ok, great, so i can just use the SNP and make the beta values for each one
# call in batch:
# cd /projects/geneva/geneva_sata/caitlin/mlm_x/
# qsub -q thornton.q -t 1-442 -N pow005 batch_model_c_powerCalc_function_005.sh
# model_c_powerCalc_function.R with arguments for the SNP (1-999) and the beta value (0.05 here)
# have other beta values of 0.1, 0.13, 0.15, 0.18, 0.2, 0.25, 0.3

# saves files:
# modelC_results_power/mmRes_SNPXX_beta0.005.RData
# for each SNP and each beta value

###
# only had 442 phenotypes simulated with that scanAnnotation file
# for SNPs 443-999, need to use the other scanAnnotation file
# update the model_c_powerCalc_function.R to _v2.R

# call in batch:
# cd /projects/geneva/geneva_sata/caitlin/mlm_x/
# qsub -q thornton.q -t 443-999 -N pow005 batch_model_c_powerCalc_function_v2_005.sh
# model_c_powerCalc_function_v2.R with arguments for the SNP (1-999) and the beta value (0.05 here)
# have other beta values of 0.1, 0.13, 0.15, 0.18, 0.2, 0.25, 0.3

# saves files:
# modelC_results_power/mmRes_SNPXX_beta0.005.RData
# for each SNP and each beta value


#####
# 14. Power simulation results

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")

# loop through results and make them into one overall file

ct <- 1; ctCI <- 1
files <- list.files("modelC_results_power")
length(files) # first iteration it was 3109; next iteration it's 889
# next iteration it's SNPs 1-442 with a different phenotype 'base' 3536 = 442*8
# next iteration, did more beta values, snps 443-500 with smBeta phenotype and 1-500 with new beta values; 2462
# 2000; this is the 5th iter with smaller beta values

res <- get(load(paste("modelC_results_power/",files[1],sep="")))
totalRes <- data.frame(matrix(NA,nrow=(6*length(files)),ncol=12))
colnames(totalRes) <- c("beta1",colnames(res[[3]]),"model")

varCI <- data.frame(matrix(NA,nrow=(7*length(files)),ncol=5))
colnames(varCI) <- c("beta1",colnames(res[[2]]),"model")

for(i in 1:length(files)){
  
  fn <- paste("modelC_results_power/",files[i],sep="")
  res <- get(load(fn))
  
  totalRes[ct:(ct+1),] <- cbind("beta1"=rep(res[[1]],2),res[[3]],"model"=rep(3,2))
  totalRes[(ct+2):(ct+3),] <- cbind("beta1"=rep(res[[1]],2),res[[5]],"model"=rep(4,2))
  totalRes[(ct+4):(ct+5),] <- cbind("beta1"=rep(res[[1]],2),res[[7]],"model"=rep(2,2))

  varCI[ctCI:(ctCI+2),] <- cbind("beta1"=rep(res[[1]],3),res[[2]],"model"=rep(3,3))
  varCI[(ctCI+3):(ctCI+4),] <- cbind("beta1"=rep(res[[1]],2),res[[4]],"model"=rep(4,2))
  varCI[(ctCI+5):(ctCI+6),] <- cbind("beta1"=rep(res[[1]],2),res[[6]],"model"=rep(2,2))
  
  ct <- ct+6
  ctCI <- ctCI+7
}

dim(totalRes) # 5334 12 | 3rd iter: 21216 12 | 4th iter: 14772 12 | 5th: 12000 12
dim(varCI) # 6223 5 | 3rd iter: 24752 5 | 4th iter: 17234 5 | 5th: 14000 5
# 1778 of each model, 2667 of var ci for model 3, since adj for both x, auto, e
table(totalRes$model) # 3rd iter: 7072 of each model | 4th iter: 4294 | 5th: 4000 of each

save(totalRes,file="pvals_power_results_jan13_v5.RData")
save(varCI,file="varEst_power_results_jan13_v5.RData")

rm(list=ls())

# removed the individual files in the model_c_results_power/ folder

tmp1 <- get(load("pvals_power_results_jan13.RData"))
tmp <- get(load("pvals_power_results_jan13_v2.RData"))
dim(tmp1); dim(tmp) # 18684 12 | 5334 12

totalRes <- rbind(tmp1,tmp)
head(totalRes); dim(totalRes) # 24018 12

save(totalRes,file="pvals_power_results_jan13_all.RData")

rm(list=ls())

tmp1 <- get(load("pvals_power_results_jan13_v3.RData"))
tmp <- get(load("pvals_power_results_jan13_all.RData"))

totalRes <- rbind(tmp1,tmp)
head(totalRes); dim(totalRes) # 45234 12
ftable(totalRes$causal,totalRes$model,totalRes$beta1) # nearly 1000 for each! yesss

# check duplicate rows
totalRes <- totalRes[order(totalRes$pval),]
sum(duplicated(totalRes$pval)) # 72
dups <- totalRes$pval[duplicated(totalRes$pval)]
dups <- totalRes[is.element(totalRes$pval,dups),]
othe <- seq(from=1,to=nrow(dups),by=2)
for(i in othe){
  stopifnot(all(dups[i,]==dups[(i+1),]))
}
# no stops
# remove these dups

totalRes <- totalRes[!duplicated(totalRes$pval),]
dim(totalRes) # 45162 12; 45234-72=45162

save(totalRes,file="pvals_power_results_jan13_all.RData")

# add the 4th iteration
tmp <- get(load("pvals_power_results_jan13_v4.RData"))
tmp1 <- get(load("pvals_power_results_jan13_all.RData"))

totalRes <- rbind(tmp,tmp1)
head(totalRes); dim(totalRes) # 59934 12
ftable(totalRes$causal,totalRes$model,totalRes$beta1) # nearly 1000 for each! yesss

# check duplicate rows
totalRes <- totalRes[order(totalRes$pval),]
sum(duplicated(totalRes$pval)) # 0, good!

save(totalRes,file="pvals_power_results_jan13_all.RData")

##
# add the 5th iteration
tmp <- get(load("pvals_power_results_jan13_v5.RData"))
tmp1 <- get(load("pvals_power_results_jan13_all.RData"))

totalRes <- rbind(tmp,tmp1)
head(totalRes); dim(totalRes) # 71934 12
ftable(totalRes$causal,totalRes$model,totalRes$beta1) # basically 1000 for each!

# check duplicate rows
totalRes <- totalRes[order(totalRes$pval),]
sum(duplicated(totalRes$pval)) # 0, good!

save(totalRes,file="pvals_power_results_jan13_all.RData")

rm(list=ls())


#####
# 15. Make power graphs

source("powerGraph.R")

totalRes <- get(load("pvals_power_results_jan13_all.RData"))
dim(totalRes) # 71934 12

powerGraph(totalRes,b=0.05,fn="power_3models_beta05.pdf")
powerGraph(totalRes,b=0.1,fn="power_3models_beta1.pdf",xlims=c(0,0.5))
powerGraph(totalRes,b=0.13,fn="power_3models_beta13.pdf")
powerGraph(totalRes,b=0.15,fn="power_3models_beta15.pdf")
powerGraph(totalRes,b=0.18,fn="power_3models_beta18.pdf")
powerGraph(totalRes,b=0.3,fn="power_3models_beta3.pdf")

powerGraph(totalRes,b=0.06,fn="power_3models_beta06.pdf")
powerGraph(totalRes,b=0.07,fn="power_3models_beta07.pdf")
powerGraph(totalRes,b=0.08,fn="power_3models_beta08.pdf")
powerGraph(totalRes,b=0.09,fn="power_3models_beta09.pdf")

rm(list=ls())


#####
# 16. Get type I error from the same results

totalRes <- get(load("pvals_power_results_jan13_all.RData"))
dim(totalRes); head(totalRes) # 71934 12
table(totalRes$causal) # 35967 false SNPs
table(totalRes$model[totalRes$causal==FALSE]) # 11989 of each model

nullS_true <- totalRes[totalRes$model==3&totalRes$causal==FALSE,]
nullS_misS <- totalRes[totalRes$model==4&totalRes$causal==FALSE,]
nullS_justX <- totalRes[totalRes$model==2&totalRes$causal==FALSE,]
dim(nullS_true); dim(nullS_misS) # 11989 12 for both
dim(nullS_justX) # 11989 12

alpha <- c(0.05,0.01,0.005,0.001,5e-04)
tyIerr <- data.frame("alpha"=alpha,"trueModel"=NA,"trueModelL"=NA,"trueModelU"=NA,
                     "misSModel"=NA,"misSModelL"=NA,"misSModelU"=NA,
                     "justX"=NA,"justXL"=NA,"justXU"=NA,stringsAsFactors=FALSE)

for(a in alpha){
  
  stderr <- sqrt(a*(1-a)/nrow(nullS_true))
  
  tyIerr$trueModel[tyIerr$alpha==a] <- sum(nullS_true$pval<a)/nrow(nullS_true)
  #stderr <- sqrt(tyIerr$trueModel[tyIerr$alpha==a]*(1-tyIerr$trueModel[tyIerr$alpha==a])/nrow(nullS_true))
  tyIerr$trueModelL[tyIerr$alpha==a] <- tyIerr$trueModel[tyIerr$alpha==a]-1.96*stderr
  tyIerr$trueModelU[tyIerr$alpha==a] <- tyIerr$trueModel[tyIerr$alpha==a]+1.96*stderr
  
  tyIerr$misSModel[tyIerr$alpha==a] <- sum(nullS_misS$pval<a)/nrow(nullS_misS)
  #stderr <- sqrt(tyIerr$misSModel[tyIerr$alpha==a]*(1-tyIerr$misSModel[tyIerr$alpha==a])/nrow(nullS_misS))
  tyIerr$misSModelL[tyIerr$alpha==a] <- tyIerr$misSModel[tyIerr$alpha==a]-1.96*stderr
  tyIerr$misSModelU[tyIerr$alpha==a] <- tyIerr$misSModel[tyIerr$alpha==a]+1.96*stderr
  
  tyIerr$justX[tyIerr$alpha==a] <- sum(nullS_justX$pval<a)/nrow(nullS_justX)
  #stderr <- sqrt(tyIerr$justX[tyIerr$alpha==a]*(1-tyIerr$justX[tyIerr$alpha==a])/nrow(nullS_justX))
  tyIerr$justXL[tyIerr$alpha==a] <- tyIerr$justX[tyIerr$alpha==a]-1.96*stderr
  tyIerr$justXU[tyIerr$alpha==a] <- tyIerr$justX[tyIerr$alpha==a]+1.96*stderr
}

# get CI for true alpha too
tyIerr$alphaL <- tyIerr$alpha-1.96*stderr
tyIerr$alphaU <- tyIerr$alpha+1.96*stderr

tyIerr$cis_true <- paste(format(tyIerr$trueModel,digits=3)," (",
                    format(tyIerr$trueModelL,digits=3),", ",
                    format(tyIerr$trueModelU,digits=3),")",sep="")

tyIerr$cis_misS <- paste(format(tyIerr$misSModel,digits=3)," (",
                         format(tyIerr$misSModelL,digits=3),", ",
                         format(tyIerr$misSModelU,digits=3),")",sep="")

tyIerr$cis_justX <- paste(format(tyIerr$justX,digits=3)," (",
                         format(tyIerr$justXL,digits=3),", ",
                         format(tyIerr$justXU,digits=3),")",sep="")

tyIerr$cis_alpha <- paste(format(tyIerr$alpha,digits=3)," (",
                          format(tyIerr$alphaL,digits=3),", ",
                          format(tyIerr$alphaU,digits=3),")",sep="")

library(xtable)
xtable(t(tyIerr[,c("cis_alpha","cis_true","cis_justX","cis_misS")]),digits=5)

rm(list=ls())


#####
# 17. More power simulations

## make matrix of random effects, ie gx+ga+epsilon
# called in batch:
# cd /projects/geneva/geneva_sata/caitlin/mlm_x/
# qsub -N effRecov batch_randomEffects_recovery.sh
# which calls randomEffects_recovery.R and writes .Rout

# now we have a data frame with 1500 random effects stored:
# randomEffects.RData

# with these, we can generate genotypes and beta values, call var comp and assoc test

# write function that takes a random effect, generates genotype, with given beta calculates phenotype
# then calls var comp and assoc test

# qsub -t 1-1500 -N pow005 batch_model_c_powerCalc_function_005.sh
# which calls model_c_powerCalc_function_v3.R and writes .Rout for beta1=0.05

# qsub -t 1-1500 -N pow007 -q thornton.q batch_model_c_powerCalc_function_007.sh
# which calls model_c_powerCalc_function_v3.R and writes .Rout for beta1=0.07






