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
# 18. Power simulation results
# 19. Compile all varComp estimates
# 20. New type I error results
# 21. Process final round of power simulations
# 22. Power graph for 10K beta1=0.05 simulations
# 23. More simulations!
# 24. Power graph + type I error for 80K beta1=0.05 simulations
# 25. 1M type I error runs
# 26. Power graph for beta1=0.07 simulations
# 27. Type I error rate from 600K runs (for comparison)



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

# removed the individual files in the modelC_results_power/ folder

# now do this on my computer, after moving files from server to my laptop
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

# qsub -t 1-1500 -N pow006 -q thornton.q batch_model_c_powerCalc_function_006.sh
# qsub -t 1-1500 -N pow008 -q thornton.q batch_model_c_powerCalc_function_008.sh
# qsub -t 1-1500 -N pow009 -q thornton.q batch_model_c_powerCalc_function_009.sh
# qsub -t 1-1500 -N pow01 -q thornton.q batch_model_c_powerCalc_function_01.sh
# qsub -t 1-1500 -N pow013 -q thornton.q batch_model_c_powerCalc_function_013.sh

## all results for these are stored in
# modelC_results_power_1mill/

##
## do more iterations with 5 causal SNPs per run
# qsub -t 1-1500 -N pow005 batch_model_c_powerCalc_function_005.sh
# which calls model_c_powerCalc_function_v4.R and writes .Rout for beta1=0.05

# qsub -t 1-1500 -N pow007 -q thornton.q batch_model_c_powerCalc_function_007.sh
# which calls model_c_powerCalc_function_v4.R and writes .Rout for beta1=0.07

# qsub -q olga.q -t 1-100 batch_model_c_powerCalc_function_005.sh # started at 9.11; 
# qsub -t 101-1500 -N pw5_005 -q thornton.q batch_model_c_powerCalc_function_005.sh
# qsub -t 1-100 -N pw5_006 -q olga.q batch_model_c_powerCalc_function_006.sh
# qsub -t 101-1500 -N pw5_006 -q bigmem.q batch_model_c_powerCalc_function_006.sh
# qsub -t 1-1500 -N pw5_007 -q all.q batch_model_c_powerCalc_function_007.sh
# qsub -t 1-1500 -N pw5_008 -q thornton.q batch_model_c_powerCalc_function_008.sh
# qsub -t 1-900 -N pw5_009 -q olga.q batch_model_c_powerCalc_function_009.sh
# qsub -t 901-1500 -N pw5_009 -q thornton.q batch_model_c_powerCalc_function_009.sh
# qsub -t 1-800 -N pw5_01 -q olga.q batch_model_c_powerCalc_function_01.sh
# qsub -t 801-1500 -N pw5_01 -q thornton.q batch_model_c_powerCalc_function_01.sh
# qsub -t 1-1500 -N pw5_013 -q thornton.q batch_model_c_powerCalc_function_013.sh

## all results for these are stored in
# modelC_results_power/


#####
# 18. Power simulation results

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")

# loop through results and make them into one overall file

ct <- 1; ctCI <- 1
files <- list.files("modelC_results_power_1mill/")
length(files) # 10467

res <- get(load(paste("modelC_results_power_1mill/",files[1],sep="")))
totalRes <- data.frame(matrix(NA,nrow=(6*length(files)),ncol=12))
colnames(totalRes) <- c("beta1",colnames(res[[3]]),"model")

varCI <- data.frame(matrix(NA,nrow=(7*length(files)),ncol=5))
colnames(varCI) <- c("beta1",colnames(res[[2]]),"model")

for(i in 1:length(files)){
  
  fn <- paste("modelC_results_power_1mill/",files[i],sep="")
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

dim(totalRes) # 62802 12 
dim(varCI) # 73269 5
table(totalRes$model) # 20934 of each of the three models
table(totalRes$snpID,totalRes$model) # 10467 of each SNP id within each model
table(totalRes$beta1,totalRes$model) # missing the first SNP for beta=0.05, but weird that it's only 2936 for beta=0.13
# that's the least interesting beta value, though
#        2    3    4
#0.05 2998 2998 2998
#0.06 3000 3000 3000
#0.07 3000 3000 3000
#0.08 3000 3000 3000
#0.09 3000 3000 3000
#0.1  3000 3000 3000
#0.13 2936 2936 2936
# 2 rows for each beta, one causal and one null

sum(is.na(totalRes$pval)) # 0
summary(totalRes$pval[totalRes$causal])
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#0.0000000 0.0000129 0.0008365 0.0441100 0.0179600 0.9983000 

summary(totalRes$pval[!totalRes$causal])
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#0.0000188 0.2373000 0.4884000 0.4914000 0.7477000 0.9999000 

head(varCI); tail(varCI)
table(varCI$beta1,varCI$model) # 2 rows for models 2,4, 3 rows for model 3
#        2    3    4
#0.05 2998 4497 2998
#0.06 3000 4500 3000
#0.07 3000 4500 3000
#0.08 3000 4500 3000
#0.09 3000 4500 3000
#0.1  3000 4500 3000
#0.13 2936 4404 2936

# all looks good!
save(totalRes,file="pvals_power_results_jan22.RData")
save(varCI,file="varEst_power_results_jan22.RData")

rm(list=ls())

# removed the individual files in the modelC_results_power_1mill/ folder

## now do this on my laptop, after moving file from server to my laptop
tmp1 <- get(load("pvals_power_results_jan22.RData"))
tmp <- get(load("pvals_power_results_jan13_all.RData")) # these are the results i already have done
dim(tmp1); dim(tmp) # 62802 12 | 71934 12

totalRes <- rbind(tmp1,tmp)
head(totalRes); dim(totalRes) # 134736 12

save(totalRes,file="pvals_power_results_jan22_all.RData")

rm(list=ls())


#####
# 19. Compile all varComp estimates

## now do this on my laptop, after moving file from server to my laptop
tmp1 <- get(load("varEst_power_results_jan13.RData"))
tmp2 <- get(load("varEst_power_results_jan13_v2.RData"))
tmp3 <- get(load("varEst_power_results_jan13_v3.RData"))
tmp4 <- get(load("varEst_power_results_jan13_v4.RData"))
tmp5 <- get(load("varEst_power_results_jan13_v5.RData"))
tmp <- get(load("varEst_power_results_jan22.RData")) # these are the results i already have done
dim(tmp1); dim(tmp) # 21798 5 | 73269 5
dim(tmp2); dim(tmp3); dim(tmp4); dim(tmp5) # 6223 5 | 24752 5 | 17234 5 | 14000 5

totalRes <- rbind(tmp1,tmp2,tmp3,tmp4,tmp5,tmp)
head(totalRes); dim(totalRes) # 157276 5

table(totalRes$beta1,totalRes$model)
# 2 entries for models 2,4, 3 entries for model 3
# the higher beta values have less iterations, that's fine
#         2    3    4
# 0.05 5008 7512 5008
# 0.06 4998 7497 4998
# 0.07 5000 7500 5000
# 0.08 5000 7500 5000
# 0.09 5000 7500 5000
# 0.1  5002 7503 5002
# 0.13 4936 7404 4936
# 0.15 1990 2985 1990
# 0.18 2002 3003 2002
# 0.2  2002 3003 2002
# 0.25 1996 2994 1996
# 0.3  2002 3003 2002

save(totalRes,file="varCI_results_jan22_all.RData")

rm(list=ls())


#####
# 20. New type I error results

totalRes <- get(load("pvals_power_results_jan22_all.RData"))
alpha <- c(0.05,0.01,0.005,0.001)

tyIerr <- data.frame("alpha"=alpha,"both"=NA,"justX"=NA,"justAuto"=NA)

(numBoth <- sum(totalRes$model==3&!totalRes$causal)) # 22456
(numX <- sum(totalRes$model==2&!totalRes$causal)) # 22456
(numAuto <- sum(totalRes$model==4&!totalRes$causal)) # 22456

for(a in 1:length(alpha)){
  tyIerr$both[a] <- sum(totalRes[totalRes$model==3&!totalRes$causal,"pval"]<alpha[a])/numBoth
  tyIerr$justX[a] <- sum(totalRes[totalRes$model==2&!totalRes$causal,"pval"]<alpha[a])/numX
  tyIerr$justAuto[a] <- sum(totalRes[totalRes$model==4&!totalRes$causal,"pval"]<alpha[a])/numAuto
}

library(xtable)
xtable(tyIerr,digits=5)

# get CIs
stdErr <- sqrt(alpha*(1-alpha)/numBoth)
(ci_both <- cbind(tyIerr$both-1.96*stdErr,tyIerr$both+1.96*stdErr)) # all contain the true value
#[1,] 0.0469801776 0.052681383
#[2,] 0.0080056994 0.010608479
#[3,] 0.0040204571 0.005865542
#[4,] 0.0005217601 0.001348564

xtable(ci_both,digits=5)
bothRes <- cbind(alpha,tyIerr$both,ci_both)
xtable(bothRes,digits=5)

(ci_justX <- cbind(tyIerr$justX-1.96*stdErr,tyIerr$justX+1.96*stdErr)) # look good!
# [1,] 0.0458223579 0.051523563
# [2,] 0.0078275733 0.010430353
# [3,] 0.0044212408 0.006266326
# [4,] 0.0006108231 0.001437627

xtable(ci_justX,digits=5)
justXRes <- cbind(alpha,tyIerr$justX,ci_justX)
xtable(justXRes,digits=5)

(ci_justAuto <- cbind(tyIerr$justAuto-1.96*stdErr,tyIerr$justAuto+1.96*stdErr)) # NONE contain true value!!
# [1,] 0.070403761 0.076104967
# [2,] 0.019494834 0.022097614
# [3,] 0.010032213 0.011877299
# [4,] 0.001501454 0.002328258

justAutoRes <- cbind(alpha,tyIerr$justAuto,ci_justAuto)
xtable(justAutoRes,digits=5)

rm(list=ls())


#####
# 21. Process final round of power simulations

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x/")

# loop through results and make them into one overall file

ct <- 1; ctCI <- 1
files <- list.files("modelC_results_power/")
length(files) # 9301 as of 1pm 1/27/14

res <- get(load(paste("modelC_results_power/",files[1],sep=""))) # it's a list of length 3
# first element is fitting x + auto adj on 5 SNPs
# second element is fitting auto adj on 5 SNPs
# third element is fitting x adj on 5 SNPs
justAuto <- res[[2]]
totalRes <- data.frame(matrix(NA,nrow=(length(files)*10*3),ncol=12),stringsAsFactors=FALSE)
colnames(totalRes) <- c("beta1",colnames(justAuto[[3]]$mm_res),"model")
totalRes$model <- as.character(totalRes$model)

varCI <- data.frame(matrix(NA,nrow=(35*length(files)),ncol=5))
colnames(varCI) <- c("beta1",colnames(justAuto[[3]]$var_est_ci),"model")

for(i in 1:length(files)){
  
  fn <- paste("modelC_results_power/",files[i],sep="")
  res <- get(load(fn))
  
  # there are 5 SNPs for each of the three models
  justAuto <- res[[2]]
  justX <- res[[3]]
  both <- res[[1]]
  
  ## look at results for fitting gX and gA
  if(length(both)==5){
  if(!is.null(both[[1]])&length(both[[1]])==4){
    totalRes[ct:(ct+1),1:11] <- cbind("beta1"=rep(both[[1]]$beta1,2),both[[1]]$mm_res)
    totalRes[ct:(ct+1),"model"] <- rep(both[[1]]$model,2)
    
    varCI[ctCI:(ctCI+2),1:4] <- cbind("beta1"=rep(both[[1]]$beta1,3),both[[1]]$var_est_ci)
    varCI[ctCI:(ctCI+2),"model"] <- rep(both[[1]]$model,3)
  }
  if(!is.null(both[[2]])&length(both[[2]])==4){
    totalRes[(ct+2):(ct+3),1:11] <- cbind("beta1"=rep(both[[2]]$beta1,2),both[[2]]$mm_res)
    totalRes[(ct+2):(ct+3),"model"] <- rep(both[[2]]$model,2)
    
    varCI[(ctCI+3):(ctCI+5),1:4] <- cbind("beta1"=rep(both[[2]]$beta1,3),both[[2]]$var_est_ci)
    varCI[(ctCI+3):(ctCI+5),"model"] <- rep(both[[2]]$model,3)
  }
  if(!is.null(both[[3]])&length(both[[3]])==4){
    totalRes[(ct+4):(ct+5),1:11] <- cbind("beta1"=rep(both[[3]]$beta1,2),both[[3]]$mm_res)
    totalRes[(ct+4):(ct+5),"model"] <- rep(both[[3]]$model,2)
    
    varCI[(ctCI+6):(ctCI+8),1:4] <- cbind("beta1"=rep(both[[3]]$beta1,3),both[[3]]$var_est_ci)
    varCI[(ctCI+6):(ctCI+8),"model"] <- rep(both[[3]]$model,3)    
  }
  if(!is.null(both[[4]])&length(both[[4]])==4){
    totalRes[(ct+6):(ct+7),1:11] <- cbind("beta1"=rep(both[[4]]$beta1,2),both[[4]]$mm_res)
    totalRes[(ct+6):(ct+7),"model"] <- rep(both[[4]]$model,2)
    
    varCI[(ctCI+9):(ctCI+11),1:4] <- cbind("beta1"=rep(both[[4]]$beta1,3),both[[4]]$var_est_ci)
    varCI[(ctCI+9):(ctCI+11),"model"] <- rep(both[[4]]$model,3)    
  }
  if(!is.null(both[[5]])&length(both[[5]])==4){
    totalRes[(ct+8):(ct+9),1:11] <- cbind("beta1"=rep(both[[5]]$beta1,2),both[[5]]$mm_res)
    totalRes[(ct+8):(ct+9),"model"] <- rep(both[[5]]$model,2)
    
    varCI[(ctCI+12):(ctCI+14),1:4] <- cbind("beta1"=rep(both[[5]]$beta1,3),both[[5]]$var_est_ci)
    varCI[(ctCI+12):(ctCI+14),"model"] <- rep(both[[5]]$model,3)    
  }}

  ## now results for fitting gA
  if(length(justAuto)==5){
  if(!is.null(justAuto[[1]])&length(justAuto[[1]])==4){
    totalRes[(ct+10):(ct+11),1:11] <- cbind("beta1"=rep(justAuto[[1]]$beta1,2),justAuto[[1]]$mm_res)
    totalRes[(ct+10):(ct+11),"model"] <- rep(justAuto[[1]]$model,2)
    
    varCI[(ctCI+15):(ctCI+16),1:4] <- cbind("beta1"=rep(justAuto[[1]]$beta1,2),justAuto[[1]]$var_est_ci)
    varCI[(ctCI+15):(ctCI+16),"model"] <- rep(justAuto[[1]]$model,2)
  }
  if(!is.null(justAuto[[2]])&length(justAuto[[2]])==4){
    totalRes[(ct+12):(ct+13),1:11] <- cbind("beta1"=rep(justAuto[[2]]$beta1,2),justAuto[[2]]$mm_res)
    totalRes[(ct+12):(ct+13),"model"] <- rep(justAuto[[2]]$model,2)
    
    varCI[(ctCI+17):(ctCI+18),1:4] <- cbind("beta1"=rep(justAuto[[2]]$beta1,2),justAuto[[2]]$var_est_ci)
    varCI[(ctCI+17):(ctCI+18),"model"] <- rep(justAuto[[2]]$model,2)
  }
  if(!is.null(justAuto[[3]])&length(justAuto[[3]])==4){
    totalRes[(ct+14):(ct+15),1:11] <- cbind("beta1"=rep(justAuto[[3]]$beta1,2),justAuto[[3]]$mm_res)
    totalRes[(ct+14):(ct+15),"model"] <- rep(justAuto[[3]]$model,2)
    
    varCI[(ctCI+19):(ctCI+20),1:4] <- cbind("beta1"=rep(justAuto[[3]]$beta1,2),justAuto[[3]]$var_est_ci)
    varCI[(ctCI+19):(ctCI+20),"model"] <- rep(justAuto[[3]]$model,2)
  }
  if(!is.null(justAuto[[4]])&length(justAuto[[4]])==4){
    totalRes[(ct+16):(ct+17),1:11] <- cbind("beta1"=rep(justAuto[[4]]$beta1,2),justAuto[[4]]$mm_res)
    totalRes[(ct+16):(ct+17),"model"] <- rep(justAuto[[4]]$model,2)
    
    varCI[(ctCI+21):(ctCI+22),1:4] <- cbind("beta1"=rep(justAuto[[4]]$beta1,2),justAuto[[4]]$var_est_ci)
    varCI[(ctCI+21):(ctCI+22),"model"] <- rep(justAuto[[4]]$model,2)
  }
  if(!is.null(justAuto[[5]])&length(justAuto[[5]])==4){
    totalRes[(ct+18):(ct+19),1:11] <- cbind("beta1"=rep(justAuto[[5]]$beta1,2),justAuto[[5]]$mm_res)
    totalRes[(ct+18):(ct+19),"model"] <- rep(justAuto[[5]]$model,2)
    
    varCI[(ctCI+23):(ctCI+24),1:4] <- cbind("beta1"=rep(justAuto[[5]]$beta1,2),justAuto[[5]]$var_est_ci)
    varCI[(ctCI+23):(ctCI+24),"model"] <- rep(justAuto[[5]]$model,2)
  }}
  
  ## now results for fitting gX
  if(length(justX)==5){
  if(!is.null(justX[[1]])&length(justX[[1]])==4){
    totalRes[(ct+20):(ct+21),1:11] <- cbind("beta1"=rep(justX[[1]]$beta1,2),justX[[1]]$mm_res)
    totalRes[(ct+20):(ct+21),"model"] <- rep(justX[[1]]$model,2)
    
    varCI[(ctCI+25):(ctCI+26),1:4] <- cbind("beta1"=rep(justX[[1]]$beta1,2),justX[[1]]$var_est_ci)
    varCI[(ctCI+25):(ctCI+26),"model"] <- rep(justX[[1]]$model,2)
  }
  if(!is.null(justX[[2]])&length(justX[[2]])==4){
    totalRes[(ct+22):(ct+23),1:11] <- cbind("beta1"=rep(justX[[2]]$beta1,2),justX[[2]]$mm_res)
    totalRes[(ct+22):(ct+23),"model"] <- rep(justX[[2]]$model,2)
    
    varCI[(ctCI+27):(ctCI+28),1:4] <- cbind("beta1"=rep(justX[[2]]$beta1,2),justX[[2]]$var_est_ci)
    varCI[(ctCI+27):(ctCI+28),"model"] <- rep(justX[[2]]$model,2)
  }
  if(!is.null(justX[[3]])&length(justX[[3]])==4){
    totalRes[(ct+24):(ct+25),1:11] <- cbind("beta1"=rep(justX[[3]]$beta1,2),justX[[3]]$mm_res)
    totalRes[(ct+24):(ct+25),"model"] <- rep(justX[[3]]$model,2)
    
    varCI[(ctCI+29):(ctCI+30),1:4] <- cbind("beta1"=rep(justX[[3]]$beta1,2),justX[[3]]$var_est_ci)
    varCI[(ctCI+29):(ctCI+30),"model"] <- rep(justX[[3]]$model,2)
  }
  if(!is.null(justX[[4]])&length(justX[[4]])==4){
    totalRes[(ct+26):(ct+27),1:11] <- cbind("beta1"=rep(justX[[4]]$beta1,2),justX[[4]]$mm_res)
    totalRes[(ct+26):(ct+27),"model"] <- rep(justX[[4]]$model,2)
    
    varCI[(ctCI+31):(ctCI+32),1:4] <- cbind("beta1"=rep(justX[[4]]$beta1,2),justX[[4]]$var_est_ci)
    varCI[(ctCI+31):(ctCI+32),"model"] <- rep(justX[[4]]$model,2)
  }
  if(!is.null(justX[[5]])&length(justX[[5]])==4){
    totalRes[(ct+28):(ct+29),1:11] <- cbind("beta1"=rep(justX[[5]]$beta1,2),justX[[5]]$mm_res)
    totalRes[(ct+28):(ct+29),"model"] <- rep(justX[[5]]$model,2)
    
    varCI[(ctCI+33):(ctCI+34),1:4] <- cbind("beta1"=rep(justX[[5]]$beta1,2),justX[[5]]$var_est_ci)
    varCI[(ctCI+33):(ctCI+34),"model"] <- rep(justX[[5]]$model,2)
  }}
  
  ct <- ct+30
  ctCI <- ctCI+35
  
  if(i%%100==0){print(i)}
}

dim(totalRes) # 161070 12 
dim(varCI) # 187915 5

table(totalRes$model,exclude=NULL) 
#auto  both     x  <NA> 
#  41240 27952 47614 44264 

table(totalRes$beta1,totalRes$model,exclude=NULL) 
# auto  both     x  <NA>
#   0.05 13192  9980 14134     0
# 0.06  5678  3878  6604     0
# 0.07   138    20   270     0
# 0.08 11786  8742 13206     0
# 0.09  4828  2248  6838     0
# 0.1   5618  3084  6562     0
# <NA>     0     0     0 44264


ftable(totalRes$causal,totalRes$beta1,totalRes$model,exclude=NULL)
# auto  both     x    NA
# FALSE 0.05   6596  4990  7067     0
# 0.06   2839  1939  3302     0
# 0.07     69    10   135     0
# 0.08   5893  4371  6603     0
# 0.09   2414  1124  3419     0
# 0.1    2809  1542  3281     0
# NA        0     0     0     0
# TRUE  0.05   6596  4990  7067     0
# 0.06   2839  1939  3302     0
# 0.07     69    10   135     0
# 0.08   5893  4371  6603     0
# 0.09   2414  1124  3419     0
# 0.1    2809  1542  3281     0
# NA        0     0     0     0
# NA    0.05      0     0     0     0
# 0.06      0     0     0     0
# 0.07      0     0     0     0
# 0.08      0     0     0     0
# 0.09      0     0     0     0
# 0.1       0     0     0     0
# NA        0     0     0 44264


summary(totalRes$pval[totalRes$causal])
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.000   0.000   0.009   0.094   0.086   0.999   14446 

summary(totalRes$pval[!totalRes$causal])
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.000   0.232   0.483   0.487   0.737   1.000   14446 

head(varCI); tail(varCI)
table(varCI$beta1,varCI$model,exclude=NULL) # 2 rows for models 2,4, 3 rows for model 3
# auto  both     x  <NA>
#   0.05 13192 14970 14134     0
# 0.06  5678  5817  6604     0
# 0.07   138    30   270     0
# 0.08 11786 13113 13206     0
# 0.09  4828  3372  6838     0
# 0.1   5618  4626  6562     0
# <NA>     0     0     0 57133

# all looks good!
save(totalRes,file="pvals_power_results_jan27.RData")
save(varCI,file="varEst_power_results_jan27.RData")

rm(list=ls())


#####
# 22. Power graph for 10K beta1=0.05 simulations

source("powerGraph.R")
tmp1 <- get(load("pvals_power_results_jan22_v2.RData"))
tmp2 <- get(load("pvals_power_results_jan22_all.RData"))
dim(tmp1); dim(tmp2) # 277980 12 | 134736 12

totalRes <- rbind(tmp1,tmp2)
dim(totalRes) # 412716 12

# be sure to exclude the NA rows.
sum(is.na(totalRes$pval)) # 73428
totalRes <- totalRes[!is.na(totalRes$pval),]
dim(totalRes) # 339288 12

totalRes$model[totalRes$model==2] <- "x"
totalRes$model[totalRes$model==3] <- "both"
totalRes$model[totalRes$model==4] <- "auto"

ftable(totalRes$causal,totalRes$beta1,totalRes$model)
# auto both    x
# FALSE 0.05  9095 7489 9566
# 0.06  5338 4438 5801
# 0.07  2569 2510 2635
# 0.08  8393 6871 9103
# 0.09  4914 3624 5919
# 0.1   5309 4042 5781
# 0.13  2467 2467 2467
# 0.15   994  994  994
# 0.18  1000 1000 1000
# 0.2   1000 1000 1000
# 0.25   997  997  997
# 0.3   1000 1000 1000
# TRUE  0.05  9095 7489 9566
# 0.06  5338 4438 5801
# 0.07  2569 2510 2635
# 0.08  8393 6871 9103
# 0.09  4914 3624 5919
# 0.1   5309 4042 5781
# 0.13  2467 2467 2467
# 0.15   994  994  994
# 0.18  1000 1000 1000
# 0.2   1000 1000 1000
# 0.25   997  997  997
# 0.3   1000 1000 1000

save(totalRes,file="pvals_power_results_jan26.RData")

powerGraph(totalRes,b=0.05,fn="power_3models_beta05.pdf")
powerGraph(totalRes,b=0.08,fn="power_3models_beta08.pdf")

rm(list=ls())


#####
# 23. More simulations!

library(GWASTools)
setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")

# calls model_c_powerCalc_function_v5.R
# which calls only one SNP at a time
# submits start values for each varCompEst function

# qsub -N pow005 -q thornton.q -t 1501-5000 batch_model_c_powerCalc_function_005.sh
# see how long it takes to do 2500 iterations, for beta1=0.05

# qsub -N res -q thornton.q -t 1-1500 -p -50 batch_model_c_powerCalc_function_005.sh

# have sim results stored in modelC_results_power, modelC_results_power_1mill..._1mill_10
folders <- c("modelC_results_power","modelC_results_power_1mill",paste("modelC_results_power_1mill_",c(2,10),sep=""))
# read in files in all these folders


res <- get(load(paste(folders[1],"/mmRes_SNP1_beta0.05.RData",sep="")))

totalRes <- data.frame(matrix(NA,nrow=(length(folders)*1500*6),ncol=12),stringsAsFactors=FALSE)
colnames(totalRes) <- c("beta1",colnames(res$auto[[4]]),"model")
totalRes$model <- as.character(totalRes$model)

varCI <- data.frame(matrix(NA,nrow=(length(folders)*1500*7),ncol=5))
colnames(varCI) <- c("beta1",colnames(res$auto[[3]]),"model")

ct <- 1; ctCI <- 1

for(i in 1:length(folders)){
  fns <- list.files(folders[i])
  print(paste(folders[i],length(fns)))
  for(j in 1:length(fns)){
    res <- get(load(paste(folders[i],"/",fns[j],sep="")))
    totalRes[ct:(ct+1),1:11] <- cbind("beta1"=rep(res$both[[1]],2),res$both[[4]])
    totalRes$model[ct:(ct+1)] <- rep("both",2)
    
    totalRes[(ct+2):(ct+3),1:11] <- cbind("beta1"=rep(res$x[[1]],2),res$x[[4]])
    totalRes$model[(ct+2):(ct+3)] <- rep("x",2)
    
    totalRes[(ct+4):(ct+5),1:11] <- cbind("beta1"=rep(res$auto[[1]],2),res$auto[[4]])
    totalRes$model[(ct+4):(ct+5)] <- rep("auto",2)
    
    varCI[ctCI:(ctCI+2),1:4] <- cbind("beta1"=rep(res$both[[1]],3),res$both[[3]])
    varCI$model[ctCI:(ctCI+2)] <- rep("both",3)
    
    varCI[(ctCI+3):(ctCI+4),1:4] <- cbind("beta1"=rep(res$x[[1]],2),res$x[[3]])
    varCI$model[(ctCI+3):(ctCI+4)] <- rep("x",2)
    
    varCI[(ctCI+5):(ctCI+6),1:4] <- cbind("beta1"=rep(res$auto[[1]],2),res$auto[[3]])
    varCI$model[(ctCI+5):(ctCI+6)] <- rep("auto",2)  
  
    if(j%%100==0){print(j)}
    
    ct <- ct+6
    ctCI <- ctCI+7
  }
}

ct; ctCI # 23899 | 27882
head(totalRes); head(varCI)

dim(totalRes) # 36000 12
dim(varCI) # 42000 5

totalRes <- totalRes[!is.na(totalRes$pval),]
varCI <- varCI[!is.na(varCI$beta1),]

table(totalRes$model,exclude=NULL) # 7966 for all

table(totalRes$beta1,totalRes$model,exclude=NULL)  # all beta1=0.05
ftable(totalRes$causal,totalRes$beta1,totalRes$model,exclude=NULL)
# 3983 each model, each causal, for beta1=0.05

summary(totalRes$pval[totalRes$causal])
#   Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.000000 0.005847 0.038810 0.141700 0.177300 0.999300 

summary(totalRes$pval[!totalRes$causal])
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#0.0000136 0.2328000 0.4828000 0.4871000 0.7363000 0.9999000 

head(varCI); tail(varCI)
table(varCI$beta1,varCI$model,exclude=NULL) # 2 rows for models 2,4, 3 rows for model 3
# 7966 or 11949 for each model

# all looks good!
save(totalRes,file="pvals_power_results_jan29.RData")
save(varCI,file="varEst_power_results_jan29.RData")

rm(list=ls())

tmp <- get(load("pvals_power_results_jan26.RData"))
tmp2 <- get(load("pvals_power_results_jan29.RData"))

dim(tmp); dim(tmp2) # 339288 12 | 23898 12

totalRes <- rbind(tmp,tmp2)
dim(totalRes) # 363186 12

ftable(totalRes$model,totalRes$causal,totalRes$beta1)
# we have 13078 for each model for beta1=0.05, causal=T
# just need 90K more... :/

save(totalRes,file="pvals_power_results_ALL.RData")

rm(list=ls())


#####
# 24. Process simulations

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x/")

# loop through results and make them into one overall file

ct <- 1; ctCI <- 1
files <- c(list.files(paste("modelC_results_power_1mill_",7:10,sep=""),full.names=TRUE))
length(files) # 13500; 1500 in each

res <- get(load(files[1])) # it's a list of length 3
# first element is fitting x + auto adj on 5 SNPs
# second element is fitting auto adj on 5 SNPs
# third element is fitting x adj on 5 SNPs
justAuto <- res[[2]]
totalRes <- data.frame(matrix(NA,nrow=(length(files)*10*3),ncol=12),stringsAsFactors=FALSE)
colnames(totalRes) <- c("beta1",colnames(justAuto[[3]]$mm_res),"model")
totalRes$model <- as.character(totalRes$model)

varCI <- data.frame(matrix(NA,nrow=(35*length(files)),ncol=5))
colnames(varCI) <- c("beta1",colnames(justAuto[[3]]$var_est_ci),"model")

for(i in 1:length(files)){
  
  res <- get(load(files[i]))
  
  # there are 5 SNPs for each of the three models
  justAuto <- res[[2]]
  justX <- res[[3]]
  both <- res[[1]]
  
  ## look at results for fitting gX and gA
  if(length(both)==5){
    if(!is.null(both[[1]])&length(both[[1]])==4){
      totalRes[ct:(ct+1),1:11] <- cbind("beta1"=rep(both[[1]]$beta1,2),both[[1]]$mm_res)
      totalRes[ct:(ct+1),"model"] <- rep(both[[1]]$model,2)
      
      varCI[ctCI:(ctCI+2),1:4] <- cbind("beta1"=rep(both[[1]]$beta1,3),both[[1]]$var_est_ci)
      varCI[ctCI:(ctCI+2),"model"] <- rep(both[[1]]$model,3)
    }
    if(!is.null(both[[2]])&length(both[[2]])==4){
      totalRes[(ct+2):(ct+3),1:11] <- cbind("beta1"=rep(both[[2]]$beta1,2),both[[2]]$mm_res)
      totalRes[(ct+2):(ct+3),"model"] <- rep(both[[2]]$model,2)
      
      varCI[(ctCI+3):(ctCI+5),1:4] <- cbind("beta1"=rep(both[[2]]$beta1,3),both[[2]]$var_est_ci)
      varCI[(ctCI+3):(ctCI+5),"model"] <- rep(both[[2]]$model,3)
    }
    if(!is.null(both[[3]])&length(both[[3]])==4){
      totalRes[(ct+4):(ct+5),1:11] <- cbind("beta1"=rep(both[[3]]$beta1,2),both[[3]]$mm_res)
      totalRes[(ct+4):(ct+5),"model"] <- rep(both[[3]]$model,2)
      
      varCI[(ctCI+6):(ctCI+8),1:4] <- cbind("beta1"=rep(both[[3]]$beta1,3),both[[3]]$var_est_ci)
      varCI[(ctCI+6):(ctCI+8),"model"] <- rep(both[[3]]$model,3)    
    }
    if(!is.null(both[[4]])&length(both[[4]])==4){
      totalRes[(ct+6):(ct+7),1:11] <- cbind("beta1"=rep(both[[4]]$beta1,2),both[[4]]$mm_res)
      totalRes[(ct+6):(ct+7),"model"] <- rep(both[[4]]$model,2)
      
      varCI[(ctCI+9):(ctCI+11),1:4] <- cbind("beta1"=rep(both[[4]]$beta1,3),both[[4]]$var_est_ci)
      varCI[(ctCI+9):(ctCI+11),"model"] <- rep(both[[4]]$model,3)    
    }
    if(!is.null(both[[5]])&length(both[[5]])==4){
      totalRes[(ct+8):(ct+9),1:11] <- cbind("beta1"=rep(both[[5]]$beta1,2),both[[5]]$mm_res)
      totalRes[(ct+8):(ct+9),"model"] <- rep(both[[5]]$model,2)
      
      varCI[(ctCI+12):(ctCI+14),1:4] <- cbind("beta1"=rep(both[[5]]$beta1,3),both[[5]]$var_est_ci)
      varCI[(ctCI+12):(ctCI+14),"model"] <- rep(both[[5]]$model,3)    
    }}
  
  ## now results for fitting gA
  if(length(justAuto)==5){
    if(!is.null(justAuto[[1]])&length(justAuto[[1]])==4){
      totalRes[(ct+10):(ct+11),1:11] <- cbind("beta1"=rep(justAuto[[1]]$beta1,2),justAuto[[1]]$mm_res)
      totalRes[(ct+10):(ct+11),"model"] <- rep(justAuto[[1]]$model,2)
      
      varCI[(ctCI+15):(ctCI+16),1:4] <- cbind("beta1"=rep(justAuto[[1]]$beta1,2),justAuto[[1]]$var_est_ci)
      varCI[(ctCI+15):(ctCI+16),"model"] <- rep(justAuto[[1]]$model,2)
    }
    if(!is.null(justAuto[[2]])&length(justAuto[[2]])==4){
      totalRes[(ct+12):(ct+13),1:11] <- cbind("beta1"=rep(justAuto[[2]]$beta1,2),justAuto[[2]]$mm_res)
      totalRes[(ct+12):(ct+13),"model"] <- rep(justAuto[[2]]$model,2)
      
      varCI[(ctCI+17):(ctCI+18),1:4] <- cbind("beta1"=rep(justAuto[[2]]$beta1,2),justAuto[[2]]$var_est_ci)
      varCI[(ctCI+17):(ctCI+18),"model"] <- rep(justAuto[[2]]$model,2)
    }
    if(!is.null(justAuto[[3]])&length(justAuto[[3]])==4){
      totalRes[(ct+14):(ct+15),1:11] <- cbind("beta1"=rep(justAuto[[3]]$beta1,2),justAuto[[3]]$mm_res)
      totalRes[(ct+14):(ct+15),"model"] <- rep(justAuto[[3]]$model,2)
      
      varCI[(ctCI+19):(ctCI+20),1:4] <- cbind("beta1"=rep(justAuto[[3]]$beta1,2),justAuto[[3]]$var_est_ci)
      varCI[(ctCI+19):(ctCI+20),"model"] <- rep(justAuto[[3]]$model,2)
    }
    if(!is.null(justAuto[[4]])&length(justAuto[[4]])==4){
      totalRes[(ct+16):(ct+17),1:11] <- cbind("beta1"=rep(justAuto[[4]]$beta1,2),justAuto[[4]]$mm_res)
      totalRes[(ct+16):(ct+17),"model"] <- rep(justAuto[[4]]$model,2)
      
      varCI[(ctCI+21):(ctCI+22),1:4] <- cbind("beta1"=rep(justAuto[[4]]$beta1,2),justAuto[[4]]$var_est_ci)
      varCI[(ctCI+21):(ctCI+22),"model"] <- rep(justAuto[[4]]$model,2)
    }
    if(!is.null(justAuto[[5]])&length(justAuto[[5]])==4){
      totalRes[(ct+18):(ct+19),1:11] <- cbind("beta1"=rep(justAuto[[5]]$beta1,2),justAuto[[5]]$mm_res)
      totalRes[(ct+18):(ct+19),"model"] <- rep(justAuto[[5]]$model,2)
      
      varCI[(ctCI+23):(ctCI+24),1:4] <- cbind("beta1"=rep(justAuto[[5]]$beta1,2),justAuto[[5]]$var_est_ci)
      varCI[(ctCI+23):(ctCI+24),"model"] <- rep(justAuto[[5]]$model,2)
    }}
  
  ## now results for fitting gX
  if(length(justX)==5){
    if(!is.null(justX[[1]])&length(justX[[1]])==4){
      totalRes[(ct+20):(ct+21),1:11] <- cbind("beta1"=rep(justX[[1]]$beta1,2),justX[[1]]$mm_res)
      totalRes[(ct+20):(ct+21),"model"] <- rep(justX[[1]]$model,2)
      
      varCI[(ctCI+25):(ctCI+26),1:4] <- cbind("beta1"=rep(justX[[1]]$beta1,2),justX[[1]]$var_est_ci)
      varCI[(ctCI+25):(ctCI+26),"model"] <- rep(justX[[1]]$model,2)
    }
    if(!is.null(justX[[2]])&length(justX[[2]])==4){
      totalRes[(ct+22):(ct+23),1:11] <- cbind("beta1"=rep(justX[[2]]$beta1,2),justX[[2]]$mm_res)
      totalRes[(ct+22):(ct+23),"model"] <- rep(justX[[2]]$model,2)
      
      varCI[(ctCI+27):(ctCI+28),1:4] <- cbind("beta1"=rep(justX[[2]]$beta1,2),justX[[2]]$var_est_ci)
      varCI[(ctCI+27):(ctCI+28),"model"] <- rep(justX[[2]]$model,2)
    }
    if(!is.null(justX[[3]])&length(justX[[3]])==4){
      totalRes[(ct+24):(ct+25),1:11] <- cbind("beta1"=rep(justX[[3]]$beta1,2),justX[[3]]$mm_res)
      totalRes[(ct+24):(ct+25),"model"] <- rep(justX[[3]]$model,2)
      
      varCI[(ctCI+29):(ctCI+30),1:4] <- cbind("beta1"=rep(justX[[3]]$beta1,2),justX[[3]]$var_est_ci)
      varCI[(ctCI+29):(ctCI+30),"model"] <- rep(justX[[3]]$model,2)
    }
    if(!is.null(justX[[4]])&length(justX[[4]])==4){
      totalRes[(ct+26):(ct+27),1:11] <- cbind("beta1"=rep(justX[[4]]$beta1,2),justX[[4]]$mm_res)
      totalRes[(ct+26):(ct+27),"model"] <- rep(justX[[4]]$model,2)
      
      varCI[(ctCI+31):(ctCI+32),1:4] <- cbind("beta1"=rep(justX[[4]]$beta1,2),justX[[4]]$var_est_ci)
      varCI[(ctCI+31):(ctCI+32),"model"] <- rep(justX[[4]]$model,2)
    }
    if(!is.null(justX[[5]])&length(justX[[5]])==4){
      totalRes[(ct+28):(ct+29),1:11] <- cbind("beta1"=rep(justX[[5]]$beta1,2),justX[[5]]$mm_res)
      totalRes[(ct+28):(ct+29),"model"] <- rep(justX[[5]]$model,2)
      
      varCI[(ctCI+33):(ctCI+34),1:4] <- cbind("beta1"=rep(justX[[5]]$beta1,2),justX[[5]]$var_est_ci)
      varCI[(ctCI+33):(ctCI+34),"model"] <- rep(justX[[5]]$model,2)
    }}
  
  ct <- ct+30
  ctCI <- ctCI+35
  
  if(i%%100==0){print(i)}
}

dim(totalRes) # 135030 12
dim(varCI) # 157535 5

totalRes <- totalRes[!is.na(totalRes$pval),]
varCI <- varCI[!is.na(varCI$model),]

table(totalRes$model,exclude=NULL)# 37450 for each
table(totalRes$beta1,totalRes$model,exclude=NULL) # just beta1=0.05
ftable(totalRes$causal,totalRes$beta1,totalRes$model,exclude=NULL) # 18725 for each causal, each model, beta=0.05

summary(totalRes$pval[totalRes$causal])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.00000 0.00567 0.03996 0.14040 0.17550 0.99980 

summary(totalRes$pval[!totalRes$causal])
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#0.0000248 0.2329000 0.4832000 0.4888000 0.7430000 1.0000000 

head(varCI); tail(varCI)
table(varCI$beta1,varCI$model,exclude=NULL) # 2 rows for models 2,4, 3 rows for model 3

# all looks good!
save(totalRes,file="pvals_power_results_feb6.RData")
save(varCI,file="varEst_power_results_feb6.RData")

rm(list=ls())

tmp <- get(load("pvals_power_results_mar23.RData"))
tmp2 <- get(load("pvals_power_results_ALL.RData"))

dim(tmp); dim(tmp2) # 900000 12 | 3590976 12

totalRes <- rbind(tmp,tmp2)
dim(totalRes) # 4490976 12

sum(duplicated(totalRes$pval)) # 490080; dang.
totalRes <- totalRes[!duplicated(totalRes$pval),]
sum(duplicated(totalRes$pval)) # 0; good
dim(totalRes) # 4040976 12

ftable(totalRes$model,totalRes$causal,totalRes$beta1)
# we have 326043 for each model for beta1=0.05, causal=T
# added some beta=0.07 iterations in there too, 303K for causal=T

save(totalRes,file="pvals_power_results_ALL.RData")

rm(list=ls())


#####
# 24. Power graph + type I error for 80K beta1=0.05 simulations

source("powerGraph.R")

totalRes <- get(load("pvals_power_results_ALL.RData"))
dim(totalRes) # 1835976 12

# be sure to exclude the NA rows.
sum(is.na(totalRes$pval)) # 0

ftable(totalRes$causal,totalRes$beta1,totalRes$model)
# auto  both     x
# FALSE 0.05  80578 78972 81049
# 0.06   7995  6142  8568
# 0.07   3012  2618  3208
# 0.08   8819  7158  9584
# 0.09   7567  5548  8799
# 0.1    8295  6293  9100
# 0.13   8949  7242  9624
# 0.15    994   994   994
# 0.18   1000  1000  1000
# 0.2    1000  1000  1000
# 0.25    997   997   997
# 0.3    1000  1000  1000
# TRUE  0.05  80578 78972 81049
# 0.06   7995  6142  8568
# 0.07   3012  2618  3208
# 0.08   8819  7158  9584
# 0.09   7567  5548  8799
# 0.1    8295  6293  9100
# 0.13   8949  7242  9624
# 0.15    994   994   994
# 0.18   1000  1000  1000
# 0.2    1000  1000  1000
# 0.25    997   997   997
# 0.3    1000  1000  1000

powerGraph(totalRes,b=0.05,fn="power_3models_beta05.pdf")

alpha <- c(0.05,0.01,0.005,0.001,5e-04,1e-04,5e-05,1e-05)

tyIerr <- data.frame("alpha"=alpha,"both"=NA,"justX"=NA,"justAuto"=NA)

(numBoth <- sum(totalRes$model=="both"&!totalRes$causal)) # 296929
(numX <- sum(totalRes$model=="x"&!totalRes$causal)) # 312888
(numAuto <- sum(totalRes$model=="auto"&!totalRes$causal)) # 308171

for(a in 1:length(alpha)){
  tyIerr$both[a] <- sum(totalRes[totalRes$model=="both"&!totalRes$causal,"pval"]<alpha[a])/numBoth
  tyIerr$justX[a] <- sum(totalRes[totalRes$model=="x"&!totalRes$causal,"pval"]<alpha[a])/numX
  tyIerr$justAuto[a] <- sum(totalRes[totalRes$model=="auto"&!totalRes$causal,"pval"]<alpha[a])/numAuto
}

library(xtable)
xtable(tyIerr,digits=6)

# get CIs
stdErr <- sqrt(alpha*(1-alpha)/numBoth)
(ci_both <- cbind(tyIerr$both-1.96*stdErr,tyIerr$both+1.96*stdErr)) # all contain the true value
# [1,]  4.887750e-02 5.135450e-02
# [2,]  9.345148e-03 1.047597e-02
# [3,]  4.583885e-03 5.385517e-03
# [4,]  8.038806e-04 1.163101e-03
# [5,]  3.268835e-04 5.809542e-04
# [6,]  3.564163e-05 1.492883e-04
# [7,]  1.866017e-05 9.902249e-05
# [8,] -9.564015e-06 2.637582e-05

xtable(ci_both,digits=5)
bothRes <- cbind(alpha,tyIerr$both,ci_both)
xtable(bothRes,digits=5)

(ci_justX <- cbind(tyIerr$justX-1.96*stdErr,tyIerr$justX+1.96*stdErr)) # look good!
# [1,]  4.905686e-02 5.153385e-02
# [2,]  9.180887e-03 1.031171e-02
# [3,]  4.779917e-03 5.581549e-03
# [4,]  7.987256e-04 1.157946e-03
# [5,]  3.102511e-04 5.643218e-04
# [6,]  4.693957e-05 1.605862e-04
# [7,]  1.170029e-05 9.206260e-05
# [8,] -1.796992e-05 1.796992e-05

xtable(ci_justX,digits=5)
justXRes <- cbind(alpha,tyIerr$justX,ci_justX)
xtable(justXRes,digits=5)

(ci_justAuto <- cbind(tyIerr$justAuto-1.96*stdErr,tyIerr$justAuto+1.96*stdErr)) # NONE contain true value!!
# [1,] 7.402685e-02 7.650385e-02
# [2,] 1.940294e-02 2.053377e-02
# [3,] 1.058178e-02 1.138341e-02
# [4,] 2.439317e-03 2.798537e-03
# [5,] 1.301470e-03 1.555541e-03
# [6,] 3.195034e-04 4.331501e-04
# [7,] 1.671826e-04 2.475449e-04
# [8,] 4.347118e-05 7.941102e-05

justAutoRes <- cbind(alpha,tyIerr$justAuto,ci_justAuto)
xtable(justAutoRes,digits=5)

rm(list=ls())


#####
# 25. 1M type I error runs

# called 
# qsub -q thornton.q -N testPheno batch_test_phenotypes.sh
# which saves test_phenotypes_250k_2_typeIerr_sigmaA05_sigmaX05.RData

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")
library(GWASTools)

dat1 <- get(load("test_phenotypes_250k_typeIerr_sigmaA05_sigmaX05.RData"))
dat2 <- get(load("test_phenotypes_250k_2_typeIerr_sigmaA05_sigmaX05.RData"))
dat3 <- get(load("test_phenotypes_250k_3_typeIerr_sigmaA05_sigmaX05.RData"))
dat4 <- get(load("test_phenotypes_150k_typeIerr_sigmaA05_sigmaX05.RData"))
dat5 <- get(load("test_phenotypes_100k_typeIerr_sigmaA05_sigmaX05.RData"))

totalRes <- rbind(dat1$mm_res,dat2$mm_res,dat3$mm_res,dat4$mm_res,dat5$mm_res)
dim(totalRes); head(totalRes)
table(totalRes$model) # 1e6 of each

# get the type I error rate from these, for each of the three methods
alpha <- c(0.05,0.01,0.005,0.001,5e-04,1e-04,5e-05,1e-05,5e-06,1e-06)

tyIerr <- data.frame("alpha"=alpha,"both"=NA,"justX"=NA,"justAuto"=NA)

(numBoth <- sum(totalRes$model=="both")) # 1000000
(numX <- sum(totalRes$model=="x")) # 1000000
(numAuto <- sum(totalRes$model=="auto")) # 1000000

for(a in 1:length(alpha)){
  tyIerr$both[a] <- sum(totalRes[totalRes$model=="both","pval"]<alpha[a])/numBoth
  tyIerr$justX[a] <- sum(totalRes[totalRes$model=="x","pval"]<alpha[a])/numX
  tyIerr$justAuto[a] <- sum(totalRes[totalRes$model=="auto","pval"]<alpha[a])/numAuto
}
tyIerr

stdErr <- sqrt(alpha*(1-alpha)/numBoth)
(ci_both <- cbind(tyIerr$both-1.96*stdErr,tyIerr$both+1.96*stdErr)) # all contain the true value

stdErr <- sqrt(alpha*(1-alpha)/numX)
(ci_x <- cbind(tyIerr$justX-1.96*stdErr,tyIerr$justX+1.96*stdErr)) # all contain the true value

stdErr <- sqrt(alpha*(1-alpha)/numAuto)
(ci_auto <- cbind(tyIerr$justAuto-1.96*stdErr,tyIerr$justAuto+1.96*stdErr)) # all contain the true value

library(xtable)
xtable(cbind(tyIerr[,c(1,2)],paste("(",format(ci_both[,1],digits=6),", ",format(ci_both[,2],digits=6),
                                   ")",sep=""),
             tyIerr[,3],paste("(",format(ci_x[,1],digits=6),", ",format(ci_x[,2],digits=6),
                              ")",sep=""),tyIerr[,4],
             paste("(",format(ci_auto[,1],digits=6),", ",format(ci_auto[,2],digits=6),
                   ")",sep="")),digits=6)


####
# called qsub -q thornton.q -N testPheno batch_test_phenotypes.sh
# which saves test_phenotypes_250k_typeIerr_sigmaA03_sigmaX05.RData

dat1 <- get(load("test_phenotypes_250k_typeIerr_sigmaA03_sigmaX05.RData"))
dat2 <- get(load("test_phenotypes_250k_2_typeIerr_sigmaA03_sigmaX05.RData"))
dat3 <- get(load("test_phenotypes_250k_3_typeIerr_sigmaA03_sigmaX05.RData"))
dat4 <- get(load("test_phenotypes_250k_4_typeIerr_sigmaA03_sigmaX05.RData"))

totalRes <- rbind(dat1$mm_res,dat2$mm_res,dat3$mm_res,dat4$mm_res)
dim(totalRes); head(totalRes)
table(totalRes$model) # 1e6 of each

# get the type I error rate from these, for each of the three methods
alpha <- c(0.05,0.01,0.005,0.001,5e-04,1e-04,5e-05,1e-05,5e-06,1e-06)

tyIerr <- data.frame("alpha"=alpha,"both"=NA,"justX"=NA,"justAuto"=NA)

(numBoth <- sum(totalRes$model=="both")) # 1000000
(numX <- sum(totalRes$model=="x")) # 1000000
(numAuto <- sum(totalRes$model=="auto")) # 1000000

for(a in 1:length(alpha)){
  tyIerr$both[a] <- sum(totalRes[totalRes$model=="both","pval"]<alpha[a])/numBoth
  tyIerr$justX[a] <- sum(totalRes[totalRes$model=="x","pval"]<alpha[a])/numX
  tyIerr$justAuto[a] <- sum(totalRes[totalRes$model=="auto","pval"]<alpha[a])/numAuto
}
tyIerr

stdErr <- sqrt(alpha*(1-alpha)/numBoth)
(ci_both <- cbind(tyIerr$both-1.96*stdErr,tyIerr$both+1.96*stdErr)) # all contain the true value

stdErr <- sqrt(alpha*(1-alpha)/numX)
(ci_x <- cbind(tyIerr$justX-1.96*stdErr,tyIerr$justX+1.96*stdErr)) # all contain the true value

stdErr <- sqrt(alpha*(1-alpha)/numAuto)
(ci_auto <- cbind(tyIerr$justAuto-1.96*stdErr,tyIerr$justAuto+1.96*stdErr)) # all contain the true value

library(xtable)
xtable(cbind(tyIerr[,c(1,2)],paste("(",format(ci_both[,1],digits=3),", ",format(ci_both[,2],digits=3),
                                   ")",sep=""),
             tyIerr[,3],paste("(",format(ci_x[,1],digits=3),", ",format(ci_x[,2],digits=3),
                              ")",sep=""),tyIerr[,4],
             paste("(",format(ci_auto[,1],digits=3),", ",format(ci_auto[,2],digits=3),
                   ")",sep="")),digits=6)


# called qsub -q thornton.q -N testPheno batch_test_phenotypes.sh
# which saves test_phenotypes_250k_typeIerr_sigmaA05_sigmaX03.RData

dat1 <- get(load("test_phenotypes_250k_typeIerr_sigmaA05_sigmaX03.RData"))
dat2 <- get(load("test_phenotypes_250k_2_typeIerr_sigmaA05_sigmaX03.RData"))
dat3 <- get(load("test_phenotypes_250k_3_typeIerr_sigmaA05_sigmaX03.RData"))
dat4 <- get(load("test_phenotypes_250k_4_typeIerr_sigmaA05_sigmaX03.RData"))

totalRes <- rbind(dat1$mm_res,dat2$mm_res,dat3$mm_res,dat4$mm_res)
dim(totalRes); head(totalRes)
table(totalRes$model) # 1e6 of each

# get the type I error rate from these, for each of the three methods
alpha <- c(0.05,0.01,0.005,0.001,5e-04,1e-04,5e-05,1e-05,5e-06,1e-06)

tyIerr <- data.frame("alpha"=alpha,"both"=NA,"justX"=NA,"justAuto"=NA)

(numBoth <- sum(totalRes$model=="both")) # 1000000
(numX <- sum(totalRes$model=="x")) # 1000000
(numAuto <- sum(totalRes$model=="auto")) # 1000000

for(a in 1:length(alpha)){
  tyIerr$both[a] <- sum(totalRes[totalRes$model=="both","pval"]<alpha[a])/numBoth
  tyIerr$justX[a] <- sum(totalRes[totalRes$model=="x","pval"]<alpha[a])/numX
  tyIerr$justAuto[a] <- sum(totalRes[totalRes$model=="auto","pval"]<alpha[a])/numAuto
}
tyIerr

stdErr <- sqrt(alpha*(1-alpha)/numBoth)
(ci_both <- cbind(tyIerr$both-1.96*stdErr,tyIerr$both+1.96*stdErr)) # all contain the true value

stdErr <- sqrt(alpha*(1-alpha)/numX)
(ci_x <- cbind(tyIerr$justX-1.96*stdErr,tyIerr$justX+1.96*stdErr)) # all contain the true value

stdErr <- sqrt(alpha*(1-alpha)/numAuto)
(ci_auto <- cbind(tyIerr$justAuto-1.96*stdErr,tyIerr$justAuto+1.96*stdErr)) # all contain the true value

library(xtable)
xtable(cbind(tyIerr[,c(1)],format(tyIerr[,2],digits=3,scientific=TRUE),paste("(",format(ci_both[,1],digits=3),", ",format(ci_both[,2],digits=3),
                                   ")",sep=""),
             format(tyIerr[,3],digits=3,scientific=TRUE),
             paste("(",format(ci_x[,1],digits=3),", ",format(ci_x[,2],digits=3),
                              ")",sep=""),
             format(tyIerr[,4],digits=3,scientific=TRUE),
             paste("(",format(ci_auto[,1],digits=3),", ",format(ci_auto[,2],digits=3),
                   ")",sep="")),digits=6)

## make a plot of these results
library(ggplot2)

######
# sigmaA05_sigmaX03

ci_both <- data.frame(ci_both)
ci_auto <- data.frame(ci_auto)
ci_x <- data.frame(ci_x)
colnames(ci_both) <- colnames(ci_auto) <- colnames(ci_x) <- c("Lower","Upper")
ci_both$model <- "both"
ci_x$model <- "x"
ci_auto$model <- "auto"
ci_both$est <- tyIerr$both; ci_both$alpha <- tyIerr$alpha
ci_auto$est <- tyIerr$justAuto; ci_auto$alpha <- tyIerr$alpha
ci_x$est <- tyIerr$justX; ci_x$alpha <- tyIerr$alpha
dfc <- rbind(ci_both[1:5,],ci_auto[1:5,],ci_x[1:5,])

pdf("tyIerr_1M_cis_a05x03.pdf")
ggplot(dfc,aes(x=model,y=est,color=as.character(alpha),ymin = Lower, ymax=Upper))  +
  geom_abline(intercept=log10(dfc$alpha[1]),slope=0,color="gray") +
  geom_abline(intercept=log10(dfc$alpha[2]), slope=0,color="gray") +
  geom_abline(intercept=log10(dfc$alpha[3]),slope=0,color="gray") +
  geom_abline(intercept=log10(dfc$alpha[4]),slope=0,color="gray") + 
  geom_abline(intercept=log10(dfc$alpha[5]),slope=0,color="gray") + 
  scale_y_log10() +
  #geom_errorbar(aes(ymin=Lower,ymax=Upper),width=0.1,position=position_jitter(width=0.1,height=0)) +
  #geom_point(size=3,position=position_jitter(width=0.1,height=0)) + 
  geom_pointrange(position=position_jitter(width=0.1,height=0)) +
  xlab("Model") + ylab("alpha") + 
  ggtitle("Estimate (95% CI) of Type I Error Rate\nSigma_A=0.5, Sigma_X=0.3") +
  theme_bw() 
  dev.off()
  
######
# sigmaA05_sigmaX05
dat1 <- get(load("test_phenotypes_250k_typeIerr_sigmaA05_sigmaX05.RData"))
dat2 <- get(load("test_phenotypes_250k_2_typeIerr_sigmaA05_sigmaX05.RData"))
dat3 <- get(load("test_phenotypes_250k_3_typeIerr_sigmaA05_sigmaX05.RData"))
dat4 <- get(load("test_phenotypes_150k_typeIerr_sigmaA05_sigmaX05.RData"))
dat5 <- get(load("test_phenotypes_100k_typeIerr_sigmaA05_sigmaX05.RData"))

totalRes <- rbind(dat1$mm_res,dat2$mm_res,dat3$mm_res,dat4$mm_res,dat5$mm_res)
dim(totalRes); head(totalRes)
table(totalRes$model) # 1e6 of each

# get the type I error rate from these, for each of the three methods
alpha <- c(0.05,0.01,0.005,0.001,5e-04,1e-04,5e-05,1e-05,5e-06,1e-06)

tyIerr <- data.frame("alpha"=alpha,"both"=NA,"justX"=NA,"justAuto"=NA)

(numBoth <- sum(totalRes$model=="both")) # 1000000
(numX <- sum(totalRes$model=="x")) # 1000000
(numAuto <- sum(totalRes$model=="auto")) # 1000000

for(a in 1:length(alpha)){
  tyIerr$both[a] <- sum(totalRes[totalRes$model=="both","pval"]<alpha[a])/numBoth
  tyIerr$justX[a] <- sum(totalRes[totalRes$model=="x","pval"]<alpha[a])/numX
  tyIerr$justAuto[a] <- sum(totalRes[totalRes$model=="auto","pval"]<alpha[a])/numAuto
}
tyIerr

stdErr <- sqrt(alpha*(1-alpha)/numBoth)
(ci_both <- cbind(tyIerr$both-1.96*stdErr,tyIerr$both+1.96*stdErr)) # all contain the true value

stdErr <- sqrt(alpha*(1-alpha)/numX)
(ci_x <- cbind(tyIerr$justX-1.96*stdErr,tyIerr$justX+1.96*stdErr)) # all contain the true value

stdErr <- sqrt(alpha*(1-alpha)/numAuto)
(ci_auto <- cbind(tyIerr$justAuto-1.96*stdErr,tyIerr$justAuto+1.96*stdErr)) # all contain the true value


ci_both <- data.frame(ci_both)
ci_auto <- data.frame(ci_auto)
ci_x <- data.frame(ci_x)
colnames(ci_both) <- colnames(ci_auto) <- colnames(ci_x) <- c("Lower","Upper")
ci_both$model <- "both"
ci_x$model <- "x"
ci_auto$model <- "auto"
ci_both$est <- tyIerr$both; ci_both$alpha <- tyIerr$alpha
ci_auto$est <- tyIerr$justAuto; ci_auto$alpha <- tyIerr$alpha
ci_x$est <- tyIerr$justX; ci_x$alpha <- tyIerr$alpha
dfc <- rbind(ci_both[1:5,],ci_auto[1:5,],ci_x[1:5,])

pdf("tyIerr_1M_cis_a05x05.pdf")
ggplot(dfc,aes(x=model,y=est,color=as.character(alpha),ymin = Lower, ymax=Upper))  +
  geom_abline(intercept=log10(dfc$alpha[1]),slope=0,color="gray") +
  geom_abline(intercept=log10(dfc$alpha[2]), slope=0,color="gray") +
  geom_abline(intercept=log10(dfc$alpha[3]),slope=0,color="gray") +
  geom_abline(intercept=log10(dfc$alpha[4]),slope=0,color="gray") + 
  geom_abline(intercept=log10(dfc$alpha[5]),slope=0,color="gray") + 
  scale_y_log10() +
  #geom_errorbar(aes(ymin=Lower,ymax=Upper),width=0.1,position=position_jitter(width=0.1,height=0)) +
  #geom_point(size=3,position=position_jitter(width=0.1,height=0)) + 
  geom_pointrange(position=position_jitter(width=0.1,height=0)) +
  xlab("Model") + ylab("alpha") + 
  ggtitle("Estimate (95% CI) of Type I Error Rate\nSigma_A=0.5, Sigma_X=0.5") +
  theme_bw() 
dev.off()


######
# sigmaA03_sigmaX05

dat1 <- get(load("test_phenotypes_250k_typeIerr_sigmaA03_sigmaX05.RData"))
dat2 <- get(load("test_phenotypes_250k_2_typeIerr_sigmaA03_sigmaX05.RData"))
dat3 <- get(load("test_phenotypes_250k_3_typeIerr_sigmaA03_sigmaX05.RData"))
dat4 <- get(load("test_phenotypes_250k_4_typeIerr_sigmaA03_sigmaX05.RData"))

totalRes <- rbind(dat1$mm_res,dat2$mm_res,dat3$mm_res,dat4$mm_res)
dim(totalRes); head(totalRes)
table(totalRes$model) # 1e6 of each

# get the type I error rate from these, for each of the three methods
alpha <- c(0.05,0.01,0.005,0.001,5e-04,1e-04,5e-05,1e-05,5e-06,1e-06)

tyIerr <- data.frame("alpha"=alpha,"both"=NA,"justX"=NA,"justAuto"=NA)

(numBoth <- sum(totalRes$model=="both")) # 1000000
(numX <- sum(totalRes$model=="x")) # 1000000
(numAuto <- sum(totalRes$model=="auto")) # 1000000

for(a in 1:length(alpha)){
  tyIerr$both[a] <- sum(totalRes[totalRes$model=="both","pval"]<alpha[a])/numBoth
  tyIerr$justX[a] <- sum(totalRes[totalRes$model=="x","pval"]<alpha[a])/numX
  tyIerr$justAuto[a] <- sum(totalRes[totalRes$model=="auto","pval"]<alpha[a])/numAuto
}
tyIerr

stdErr <- sqrt(alpha*(1-alpha)/numBoth)
(ci_both <- cbind(tyIerr$both-1.96*stdErr,tyIerr$both+1.96*stdErr)) # all contain the true value

stdErr <- sqrt(alpha*(1-alpha)/numX)
(ci_x <- cbind(tyIerr$justX-1.96*stdErr,tyIerr$justX+1.96*stdErr)) # all contain the true value

stdErr <- sqrt(alpha*(1-alpha)/numAuto)
(ci_auto <- cbind(tyIerr$justAuto-1.96*stdErr,tyIerr$justAuto+1.96*stdErr)) # all contain the true value


ci_both <- data.frame(ci_both)
ci_auto <- data.frame(ci_auto)
ci_x <- data.frame(ci_x)
colnames(ci_both) <- colnames(ci_auto) <- colnames(ci_x) <- c("Lower","Upper")
ci_both$model <- "both"
ci_x$model <- "x"
ci_auto$model <- "auto"
ci_both$est <- tyIerr$both; ci_both$alpha <- tyIerr$alpha
ci_auto$est <- tyIerr$justAuto; ci_auto$alpha <- tyIerr$alpha
ci_x$est <- tyIerr$justX; ci_x$alpha <- tyIerr$alpha
dfc <- rbind(ci_both[1:5,],ci_auto[1:5,],ci_x[1:5,])

pdf("tyIerr_1M_cis_a03x05.pdf")
ggplot(dfc,aes(x=model,y=est,color=as.character(alpha),ymin = Lower, ymax=Upper))  +
  geom_abline(intercept=log10(dfc$alpha[1]),slope=0,color="gray") +
  geom_abline(intercept=log10(dfc$alpha[2]), slope=0,color="gray") +
  geom_abline(intercept=log10(dfc$alpha[3]),slope=0,color="gray") +
  geom_abline(intercept=log10(dfc$alpha[4]),slope=0,color="gray") + 
  geom_abline(intercept=log10(dfc$alpha[5]),slope=0,color="gray") + 
  scale_y_log10() +
  #geom_errorbar(aes(ymin=Lower,ymax=Upper),width=0.1,position=position_jitter(width=0.1,height=0)) +
  #geom_point(size=3,position=position_jitter(width=0.1,height=0)) + 
  geom_pointrange(position=position_jitter(width=0.1,height=0)) +
  xlab("Model") + ylab("alpha") + 
  ggtitle("Estimate (95% CI) of Type I Error Rate\nSigma_A=0.3, Sigma_X=0.5") +
  theme_bw() 
dev.off()

rm(list=ls())


#####
# 26. Power graph for beta1=0.07 simulations

source("powerGraph.R")

totalRes <- get(load("pvals_power_results_ALL.RData"))
dim(totalRes) # 4040976 12

# be sure to exclude the NA rows.
sum(is.na(totalRes$pval)) # 0

ftable(totalRes$causal,totalRes$beta1,totalRes$model)
#               auto   both      x
# FALSE 0.05  326043 324437 326514
#       0.06    7995   6142   8568
#       0.07  303012 302618 303208
#       0.08    8819   7158   9584
#       0.09    7567   5548   8799
#       0.1     8295   6293   9100
#       0.13    8949   7242   9624
#       0.15     994    994    994
#       0.18    1000   1000   1000
#       0.2     1000   1000   1000
#       0.25     997    997    997
#       0.3     1000   1000   1000
# TRUE  0.05  326043 324437 326514
#       0.06    7995   6142   8568
#       0.07  303012 302618 303208
#       0.08    8819   7158   9584
#       0.09    7567   5548   8799
#       0.1     8295   6293   9100
#       0.13    8949   7242   9624
#       0.15     994    994    994
#       0.18    1000   1000   1000
#       0.2     1000   1000   1000
#       0.25     997    997    997
#       0.3     1000   1000   1000

powerGraph(totalRes,b=0.07,fn="power_3models_beta07.pdf")
powerGraph(totalRes,b=0.05,fn="power_3models_beta05.pdf")

rm(list=ls())


#####
# 27. Type I error rate from 600K runs (for comparison)

# get the type I error rate from these, for each of the three methods
alpha <- c(0.05,0.01,0.005,0.001,5e-04,1e-04,5e-05,1e-05,5e-06,1e-06)

totalRes <- get(load("pvals_power_results_ALL.RData"))
dim(totalRes) # 4040976 12

tyIerr <- data.frame("alpha"=alpha,"both"=NA,"justX"=NA,"justAuto"=NA)

(numBoth <- sum(totalRes$model=="both"&!totalRes$causal)) # 589429
(numX <- sum(totalRes$model=="x"&!totalRes$causal)) # 605388
(numAuto <- sum(totalRes$model=="auto"&!totalRes$causal)) # 600671

for(a in 1:length(alpha)){
  tyIerr$both[a] <- sum(totalRes[!totalRes$causal&totalRes$model=="both","pval"]<alpha[a])/numBoth
  tyIerr$justX[a] <- sum(totalRes[!totalRes$causal&totalRes$model=="x","pval"]<alpha[a])/numX
  tyIerr$justAuto[a] <- sum(totalRes[!totalRes$causal&totalRes$model=="auto","pval"]<alpha[a])/numAuto
}
tyIerr
#    alpha         both        justX     justAuto
# 1  5e-02 5.004899e-02 5.001558e-02 7.486336e-02
# 2  1e-02 9.948392e-03 9.923749e-03 1.941921e-02
# 3  5e-03 4.977206e-03 5.013316e-03 1.085588e-02
# 4  1e-03 9.707584e-04 9.847322e-04 2.764659e-03
# 5  5e-04 4.801115e-04 4.747291e-04 1.521451e-03
# 6  1e-04 8.578795e-05 9.259423e-05 4.025628e-04
# 7  5e-05 4.214145e-05 3.968324e-05 2.294016e-04
# 8  1e-05 9.030310e-06 8.818498e-06 5.772040e-05
# 9  5e-06 6.020207e-06 5.878998e-06 3.552025e-05
# 10 1e-06 3.010103e-06 2.939499e-06 5.920041e-06

stdErr <- sqrt(alpha*(1-alpha)/numBoth)
(ci_both <- cbind(tyIerr$both-1.96*stdErr,tyIerr$both+1.96*stdErr)) 

stdErr <- sqrt(alpha*(1-alpha)/numX)
(ci_x <- cbind(tyIerr$justX-1.96*stdErr,tyIerr$justX+1.96*stdErr))

stdErr <- sqrt(alpha*(1-alpha)/numAuto)
(ci_auto <- cbind(tyIerr$justAuto-1.96*stdErr,tyIerr$justAuto+1.96*stdErr)) 

library(xtable)
xtable(cbind(tyIerr[,c(1,2)],paste("(",format(ci_both[,1],digits=3),", ",format(ci_both[,2],digits=3),
                                   ")",sep=""),
             tyIerr[,3],paste("(",format(ci_x[,1],digits=3),", ",format(ci_x[,2],digits=3),
                              ")",sep=""),tyIerr[,4],
             paste("(",format(ci_auto[,1],digits=3),", ",format(ci_auto[,2],digits=3),
                   ")",sep="")),digits=6)


# make a plot of these results
# have 4 panels, a grid for sigmax, sigmaa = low or medium
# plot log(type I error rate) for each, with cis
library(ggplot2)
library(reshape)

tyIerr$ci_bothL <- ci_both[,1]
tyIerr$ci_bothU <- ci_both[,2]

tyIerr$ci_xL <- ci_x[,1]
tyIerr$ci_xU <- ci_x[,2]

tyIerr$ci_autoL <- ci_auto[,1]
tyIerr$ci_autoU <- ci_auto[,2]

tyIerr <- rbind(tyIerr,tyIerr,tyIerr)
tyIerr$model <- c(rep("both",10),rep("auto",10),rep("x",10))
tyIerr$est <- c(tyIerr$both[1:10],tyIerr$justAuto[1:10],tyIerr$justX[1:10])

tyIerr <- tyIerr[,-c(2:4)]
tyIerr$lower <- c(tyIerr$ci_bothL[1:10],tyIerr$ci_autoL[1:10],tyIerr$ci_xL[1:10])
tyIerr$upper <- c(tyIerr$ci_bothU[1:10],tyIerr$ci_autoU[1:10],tyIerr$ci_xU[1:10])
tyIerr <- tyIerr[,-c(2:7)]

tyIerr$lower[tyIerr$lower==1e-20] <- 0

kp <- c(5e-02,5e-03,5e-04,5e-05,5e-06)
tyIerrSm <- tyIerr[is.element(tyIerr$alpha,kp),]
ggplot(tyIerrSm,aes(x=model,y=log(est),fill=as.factor(alpha),ymin=log(lower), ymax=log(upper)))+geom_pointrange()+
  coord_flip()+theme_bw()+xlab(expression(alpha)) +
  ylab("log(Type I Error Rate)")





