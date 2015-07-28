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
# 28. Check a few iterations of model; what's going on with the 'x' model?
# 29. Rerun some null simulations
# 30. Look at power graph zoomed in on x-axis
# 31. Rerun power calculations
# 32. Simulate phenotypes with different pedigree structure
# 33. Look at new power results
# 34. New null simulations with extreme sigma values
# 35. Boxplots of var comp estimates
# 36. Add new null sims to var comp estimates
# 37. Process simulations excl g_x term
# 38. Null simulations with different 8 person pedigree
# 39. Var comp estimates for 8 person pedigree (fem)
# 40. Var comp estimates for 8 ped excl gX term
# 41. Process null simulations for 16 ped
# 42. Var comp estimates for 16 person pedigree
# 43. Look at new power results, beta=0.07
# 44. New power results, extreme sigma, beta=0.05
# 45. Type I error from causal=F SNPs in power results
# 46. Var comp estiamtes for 16 person ped, extreme sigma
# 47. Null simulations with auto SNP testing, 8 person ped, 'fem'
# 48. Null simulations with auto SNP testing, 8 person ped
# 49. Null simulations with auto SNP testing, 16 person ped
# 50. Look at beta estimates, 8 person ped
# 51. Var of beta estimates
# 52. Look at beta estimates, 8 person ped fem
# 53. Null simulations with auto SNP, 8 person ped fem, extreme sigma
# 54. Null simulations with x chr SNP, 8 person ped fem, extreme sigma
# 55. Look at beta estimates, 8 person ped - different sigma values
# 56. Look at beta estimates, 8 person ped fem - different sigma values

# 57. Power estimates, 8 ped: process results
# 58. Power estimates, 8 ped: make power graph
# 59. Power estimates, 8 ped fem: process results
# 60. Power estimates, 8 ped fem: make power graph
# 61. Power simulations, auto SNP
# 62. Make power graph powerSims_8ped/autoSNP results
# 63. Combine more var comp estimates to make updated boxplots, 8ped_fem
# 64. Power simulations, auto SNP, extreme sigma
# 65. Make power graph powerSims_8ped/autoSNP/extremeSigma results
# 66. Null simulations, 10K, 8 person ped fem
# 67. Make power graph powerSims_8ped_fem/autoSNP/extremeSigma results
# 68. Type I error, 8ped_fem autoSNP
# 69. Type I error, 8ped autoSNP, extreme sigma
# 70. Type I error, 8ped, autoSNP
# 71. Type I error, 8ped fem, X SNP
# 72. Null simulations, 8ped fem, X SNP -- but for var comp estimates
# 73. Make pdf for paper that includes all var comp boxplots
# 74. Power simulations, all configs
# 75. Parse extreme sig type I error results for paper
# 76. Make tyI error plots for JSM poster (x + auto SNPs)
# 77. Make power plots for JSM poster (x + auto SNPs)


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

tmp <- get(load("pvals_power_results_apr7.RData"))
tmp2 <- get(load("pvals_power_results_ALL.RData"))

dim(tmp); dim(tmp2) # 450000 12 | 4490976 12

totalRes <- rbind(tmp,tmp2)
dim(totalRes) # 4940976 12

sum(duplicated(totalRes$pval)) # 0; great!
#totalRes <- totalRes[!duplicated(totalRes$pval),]
#sum(duplicated(totalRes$pval)) # 0; good
#dim(totalRes) # 4040976 12

ftable(totalRes$model,totalRes$causal,totalRes$beta1)
# we have 326043 for each model for beta1=0.05, causal=T
# we have 378012 for each model for beta1=0.07, causal=T

save(totalRes,file="pvals_power_results_ALL.RData")

rm(list=ls())


#####
# 24. Power graph + type I error for 80K beta1=0.05 simulations

source("powerGraph.R")

totalRes <- get(load("pvals_power_results_ALL.RData"))
dim(totalRes) # 4490976 12

# be sure to exclude the NA rows.
sum(is.na(totalRes$pval)) # 0

ftable(totalRes$causal,totalRes$beta1,totalRes$model)
# auto  both     x
# FALSE 0.05  326043 324437 326514
# 0.06    7995   6142   8568
# 0.07  378012 377618 378208
# 0.08    8819   7158   9584
# 0.09    7567   5548   8799
# 0.1     8295   6293   9100
# 0.13    8949   7242   9624
# 0.15     994    994    994
# 0.18    1000   1000   1000
# 0.2     1000   1000   1000
# 0.25     997    997    997
# 0.3     1000   1000   1000
# TRUE  0.05  326043 324437 326514
# 0.06    7995   6142   8568
# 0.07  378012 377618 378208
# 0.08    8819   7158   9584
# 0.09    7567   5548   8799
# 0.1     8295   6293   9100
# 0.13    8949   7242   9624
# 0.15     994    994    994
# 0.18    1000   1000   1000
# 0.2     1000   1000   1000
# 0.25     997    997    997
# 0.3     1000   1000   1000

powerGraph(totalRes,b=0.05,fn="power_3models_beta05.pdf")

alpha <- c(0.05,0.01,0.005,0.001,5e-04,1e-04,5e-05,1e-05)

tyIerr <- data.frame("alpha"=alpha,"both"=NA,"justX"=NA,"justAuto"=NA)

(numBoth <- sum(totalRes$model=="both"&!totalRes$causal)) # 739429
(numX <- sum(totalRes$model=="x"&!totalRes$causal)) # 755388
(numAuto <- sum(totalRes$model=="auto"&!totalRes$causal)) # 750671

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
# [1,] 4.953778e-02 5.053132e-02
# [2,] 9.722778e-03 1.017636e-02
# [3,] 4.793053e-03 5.114592e-03
# [4,] 8.949195e-04 1.039005e-03
# [5,] 4.304979e-04 5.324073e-04
# [6,] 6.240868e-05 1.079931e-04
# [7,] 2.445493e-05 5.668876e-05
# [8,] 9.065207e-07 1.532222e-05

xtable(ci_both,digits=5)
bothRes <- cbind(alpha,tyIerr$both,ci_both)
xtable(bothRes,digits=5)

(ci_justX <- cbind(tyIerr$justX-1.96*stdErr,tyIerr$justX+1.96*stdErr)) # look good!
# [1,] 4.948284e-02 5.047638e-02
# [2,] 9.715120e-03 1.016870e-02
# [3,] 4.836662e-03 5.158202e-03
# [4,] 9.142053e-04 1.058291e-03
# [5,] 4.295930e-04 5.315024e-04
# [6,] 6.590395e-05 1.114883e-04
# [7,] 2.227395e-05 5.450778e-05
# [8,] 7.350893e-07 1.515079e-05

xtable(ci_justX,digits=5)
justXRes <- cbind(alpha,tyIerr$justX,ci_justX)
xtable(justXRes,digits=5)

(ci_justAuto <- cbind(tyIerr$justAuto-1.96*stdErr,tyIerr$justAuto+1.96*stdErr)) # NONE contain true value!!
# [1,] 7.443486e-02 7.542840e-02
# [2,] 1.921581e-02 1.966940e-02
# [3,] 1.069219e-02 1.101373e-02
# [4,] 2.697479e-03 2.841565e-03
# [5,] 1.465022e-03 1.566932e-03
# [6,] 3.755181e-04 4.211025e-04
# [7,] 2.076829e-04 2.399167e-04
# [8,] 4.874209e-05 6.315779e-05

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
dfc <- rbind(ci_both,ci_auto,ci_x)
dfc$sigma_a <- 0.5
dfc$sigma_x <- 0.3

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
tmpdfc <- rbind(ci_both,ci_auto,ci_x)
tmpdfc$sigma_a <- 0.5
tmpdfc$sigma_x <- 0.5
dfc <- rbind(dfc,tmpdfc)

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
tmpdfc <- rbind(ci_both,ci_auto,ci_x)
tmpdfc$sigma_a <- 0.3
tmpdfc$sigma_x <- 0.5
dfc <- rbind(dfc,tmpdfc)


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


## make a plot of all values, for alpha=5e-02

dfc$sigma_a <- paste0("sigma_a=",dfc$sigma_a)
dfc$sigma_x <- paste0("sigma_x=",dfc$sigma_x)

dfcSm <- dfc[dfc$alpha==5e-02,]
ggplot(dfcSm,aes(x=model,y=est)) +geom_point()+ geom_pointrange(aes(ymin=Lower,ymax=Upper))+
  facet_grid(sigma_x~sigma_a) + theme_bw() + geom_hline(aes(yintercept=5e-02,color="gray"))+
  ylab("Type I Error Rate")
dev.off()

dfcSm <- dfc[dfc$alpha==5e-03,]
pdf("typeIerror_alpha5e-3.pdf")
ggplot(dfcSm,aes(x=model,y=est)) +geom_point()+ geom_pointrange(aes(ymin=Lower,ymax=Upper))+
  facet_grid(sigma_x~sigma_a) + theme_bw() + geom_hline(aes(yintercept=5e-03,color="gray"))+
  ylab("Type I Error Rate") + ggtitle(expression(paste("Type I Error Rate, ", alpha, "=0.005")))
dev.off()

pdf("typeIerror_alpha5e-4.pdf")
dfcSm <- dfc[dfc$alpha==5e-04,]
ggplot(dfcSm,aes(x=model,y=est)) +geom_point()+ geom_pointrange(aes(ymin=Lower,ymax=Upper))+
  facet_grid(sigma_x~sigma_a) + theme_bw() + geom_hline(aes(yintercept=5e-04,color="gray"))+
  ylab("Type I Error Rate") + ggtitle(expression(paste("Type I Error Rate, ", alpha, "=0.0005")))
dev.off()

pdf("typeIerror_alpha5e-5.pdf")
dfcSm <- dfc[dfc$alpha==5e-05,]
ggplot(dfcSm,aes(x=model,y=est)) +geom_point()+ geom_pointrange(aes(ymin=Lower,ymax=Upper))+
  facet_grid(sigma_x~sigma_a) + theme_bw() + geom_hline(aes(yintercept=5e-05,color="gray"))+
  ylab("Type I Error Rate") + ggtitle(expression(paste("Type I Error Rate, ", alpha, "=5e-05")))
dev.off()


pdf("typeIerror_alpha5e-6.pdf")
dfcSm <- dfc[dfc$alpha==5e-06,]
ggplot(dfcSm,aes(x=model,y=est)) +geom_point()+ geom_pointrange(aes(ymin=Lower,ymax=Upper))+
  facet_grid(sigma_x~sigma_a) + theme_bw() + geom_hline(aes(yintercept=5e-06,color="gray"))+
  ylab("Type I Error Rate") + ggtitle(expression(paste("Type I Error Rate, ", alpha, "=5e-06")))
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
dim(totalRes) # 4490976 12

tyIerr <- data.frame("alpha"=alpha,"both"=NA,"justX"=NA,"justAuto"=NA)

(numBoth <- sum(totalRes$model=="both"&!totalRes$causal)) # 739429
(numX <- sum(totalRes$model=="x"&!totalRes$causal)) # 755388
(numAuto <- sum(totalRes$model=="auto"&!totalRes$causal)) # 750671

for(a in 1:length(alpha)){
  tyIerr$both[a] <- sum(totalRes[!totalRes$causal&totalRes$model=="both","pval"]<alpha[a])/numBoth
  tyIerr$justX[a] <- sum(totalRes[!totalRes$causal&totalRes$model=="x","pval"]<alpha[a])/numX
  tyIerr$justAuto[a] <- sum(totalRes[!totalRes$causal&totalRes$model=="auto","pval"]<alpha[a])/numAuto
}
tyIerr
# alpha         both        justX     justAuto
# 1 5e-02 5.003455e-02 4.997961e-02 7.493163e-02
# 2 1e-02 9.949569e-03 9.941911e-03 1.944261e-02
# 3 5e-03 4.953822e-03 4.997432e-03 1.085296e-02
# 4 1e-03 9.669623e-04 9.862481e-04 2.769522e-03
# 5 5e-04 4.814526e-04 4.805477e-04 1.515977e-03
# 6 1e-04 8.520088e-05 8.869614e-05 3.983103e-04
# 7 5e-05 4.057185e-05 3.839087e-05 2.237998e-04
# 8 1e-05 8.114369e-06 7.942938e-06 5.594994e-05


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

dfcSm <- tyIerr[tyIerr$alpha==5e-03,]
ggplot(dfcSm,aes(x=model,y=est)) +geom_point()+ geom_pointrange(aes(ymin=lower,ymax=upper))+
      theme_bw() + geom_hline(aes(yintercept=5e-03,color="gray"))+
     ylab("Type I Error Rate") + ggtitle(expression(paste("Type I Error Rate, ", alpha, "=5e-03")))
dev.off()

dfcSm <- tyIerr[tyIerr$alpha==5e-04,]
ggplot(dfcSm,aes(x=model,y=est)) +geom_point()+ geom_pointrange(aes(ymin=lower,ymax=upper))+
  theme_bw() + geom_hline(aes(yintercept=5e-04,color="gray"))+
  ylab("Type I Error Rate") + ggtitle(expression(paste("Type I Error Rate, ", alpha, "=5e-04")))
dev.off()

dfcSm <- tyIerr[tyIerr$alpha==5e-05,]
ggplot(dfcSm,aes(x=model,y=est)) +geom_point()+ geom_pointrange(aes(ymin=lower,ymax=upper))+
  theme_bw() + geom_hline(aes(yintercept=5e-05,color="gray"))+
  ylab("Type I Error Rate") + ggtitle(expression(paste("Type I Error Rate, ", alpha, "=5e-05")))
dev.off()

rm(list=ls())


#####
# 28. Check a few iterations of model; what's going on with the 'x' model?

# simulate 100 genotypes
# simulate 100 phenotypes, all with beta1=0
# perform the mlm_x model
# see how many significant results

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")

library(MASS); library(parallel)
library(GWASTools)
library(corpcor); library(doSNOW)
source("sim_phenotype.R")
source("estVarComp.R")
source("allele_drop_functions.R")
source("assocTestMixedModel_v7.R")

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

genoMat <- matrix(NA,nrow=10000,ncol=n)
colsToGeno <- split(1:n,rep(1:n,each=16,length=n))
geno <- do.call(cbind,lapply(colsToGeno,function(x){genoMat[,x]=Family_alleles_NmarkerX(rep(0.2,10000),10000,SEX)}))
dim(geno) # 100 8000; so we have 100 genotypes simulated

##
# sigmaXMat <- mvrnorm(500,mu=rep(0,16),Sigma=kinX[1:16,1:16])
# sigmaAMat <- mvrnorm(500,mu=rep(0,16),Sigma=kinAuto[1:16,1:16])
#
# noise <- rnorm(n,mean=0,sd=1)
#
# sigmaA <- 0.5
# sigmaX <- 0.5
# pheno <- noise+sigmaA*as.vector(t(sigmaAMat))+sigmaX*as.vector(t(sigmaXMat))

## perform MLM-X on the 100 simulated genotypes
scan1 <- data.frame(scanID=1:n,sex=sex,pheno=pheno)
scan1 <- ScanAnnotationDataFrame(scan1)

genoMt <- MatrixGenotypeReader(genotype=geno,snpID=as.integer(1:nrow(geno)),
                               chromosome=as.integer(rep(23,nrow(geno))),
                               position=as.integer(1:nrow(geno)), scanID=scan1$scanID)

genoData <- GenotypeData(genoMt,scanAnnot=scan1)

snpStart <- 1
snpEnd <- 10000

#colnames(kinX) <- getScanID(genoData)
#colnames(kinAuto) <- getScanID(genoData)

# covMatListBoth <- list(kinX,kinAuto) # this is the X chr kinship matrix for all samples
# names(covMatListBoth) <- c("kinshipX","kinshipAutos")

# varCompBoth <- estVarComp(scan1,covMatList=covMatListBoth,"pheno")
# varBothci <- estVarCompCI(varCompBoth,prop=FALSE)
# varBothci$model <- "both"
#
# cholSigBoth <- varCompBoth[["cholSigmaInv"]]
mmResBoth <- assocTestMixedModel(genoData,snpStart=snpStart,snpEnd=snpEnd,cholSigBoth,outcome="pheno")
mmResBoth$causal <- FALSE
mmResBoth$model <- "both"

# varCompA <- estVarComp(scan1,covMatList=list(kinAuto),"pheno")
# varAutoci <- estVarCompCI(varCompA,prop=FALSE)
# varAutoci$model <- "auto"
#
# cholSigA <- varCompA[["cholSigmaInv"]]
mmResA <- assocTestMixedModel(genoData,snpStart=snpStart,snpEnd=snpEnd,cholSigA,outcome="pheno")
mmResA$causal <- FALSE
mmResA$model <- "auto"

# varCompX <- estVarComp(scan1,covMatList=list(kinX),"pheno")
# varXci <- estVarCompCI(varCompX,prop=FALSE)
# varXci$model <- "x"
#
# cholSigX <- varCompX[["cholSigmaInv"]]
mmResX <- assocTestMixedModel(genoData,snpStart=snpStart,snpEnd=snpEnd,cholSigX,outcome="pheno")
mmResX$causal <- FALSE
mmResX$model <- "x"

mmResF <- rbind(mmResBoth,mmResA,mmResX)
res <- list(sigmaA,sigmaX,mmResF)
names(res) <- c("sigma_auto","sigma_x","mm_res")

#####
sigmaA <- 0.5
sigmaX <- 0.3
pheno <- noise+sigmaA*as.vector(t(sigmaAMat))+sigmaX*as.vector(t(sigmaXMat))

## perform MLM-X on the 100 simulated genotypes
scan1 <- data.frame(scanID=1:n,sex=sex,pheno=pheno)
scan1 <- ScanAnnotationDataFrame(scan1)

genoData <- GenotypeData(genoMt,scanAnnot=scan1)

varCompBoth <- estVarComp(scan1,covMatList=covMatListBoth,"pheno")
varBothci <- estVarCompCI(varCompBoth,prop=FALSE)
varBothci$model <- "both"

cholSigBoth <- varCompBoth[["cholSigmaInv"]]
mmResBoth <- assocTestMixedModel(genoData,snpStart=snpStart,snpEnd=snpEnd,cholSigBoth,outcome="pheno")
mmResBoth$causal <- FALSE
mmResBoth$model <- "both"

varCompA <- estVarComp(scan1,covMatList=list(kinAuto),"pheno")
varAutoci <- estVarCompCI(varCompA,prop=FALSE)
varAutoci$model <- "auto"

cholSigA <- varCompA[["cholSigmaInv"]]
mmResA <- assocTestMixedModel(genoData,snpStart=snpStart,snpEnd=snpEnd,cholSigA,outcome="pheno")
mmResA$causal <- FALSE
mmResA$model <- "auto"

varCompX <- estVarComp(scan1,covMatList=list(kinX),"pheno")
varXci <- estVarCompCI(varCompX,prop=FALSE)
varXci$model <- "x"

cholSigX <- varCompX[["cholSigmaInv"]]
mmResX <- assocTestMixedModel(genoData,snpStart=snpStart,snpEnd=snpEnd,cholSigX,outcome="pheno")
mmResX$causal <- FALSE
mmResX$model <- "x"

mmResF2 <- rbind(mmResBoth,mmResA,mmResX)
mmResF2$sigma_a <- sigmaA
mmResF2$sigma_x <- sigmaX
res2 <- list(sigmaA,sigmaX,mmResF2)
names(res2) <- c("sigma_auto","sigma_x","mm_res")

#####
sigmaA <- 0.3
sigmaX <- 0.3
pheno <- noise+sigmaA*as.vector(t(sigmaAMat))+sigmaX*as.vector(t(sigmaXMat))

## perform MLM-X on the 100 simulated genotypes
scan1 <- data.frame(scanID=1:n,sex=sex,pheno=pheno)
scan1 <- ScanAnnotationDataFrame(scan1)

genoData <- GenotypeData(genoMt,scanAnnot=scan1)

varCompBoth <- estVarComp(scan1,covMatList=covMatListBoth,"pheno")
varBothci <- estVarCompCI(varCompBoth,prop=FALSE)
varBothci$model <- "both"

cholSigBoth <- varCompBoth[["cholSigmaInv"]]
mmResBoth <- assocTestMixedModel(genoData,snpStart=snpStart,snpEnd=snpEnd,cholSigBoth,outcome="pheno")
mmResBoth$causal <- FALSE
mmResBoth$model <- "both"

varCompA <- estVarComp(scan1,covMatList=list(kinAuto),"pheno")
varAutoci <- estVarCompCI(varCompA,prop=FALSE)
varAutoci$model <- "auto"

cholSigA <- varCompA[["cholSigmaInv"]]
mmResA <- assocTestMixedModel(genoData,snpStart=snpStart,snpEnd=snpEnd,cholSigA,outcome="pheno")
mmResA$causal <- FALSE
mmResA$model <- "auto"

varCompX <- estVarComp(scan1,covMatList=list(kinX),"pheno")
varXci <- estVarCompCI(varCompX,prop=FALSE)
varXci$model <- "x"

cholSigX <- varCompX[["cholSigmaInv"]]
mmResX <- assocTestMixedModel(genoData,snpStart=snpStart,snpEnd=snpEnd,cholSigX,outcome="pheno")
mmResX$causal <- FALSE
mmResX$model <- "x"

mmResF3 <- rbind(mmResBoth,mmResA,mmResX)
mmResF3$sigma_a <- sigmaA
mmResF3$sigma_x <- sigmaX
res3 <- list(sigmaA,sigmaX,mmResF3)
names(res3) <- c("sigma_auto","sigma_x","mm_res")

#####
sigmaA <- 0.3
sigmaX <- 0.5
pheno <- noise+sigmaA*as.vector(t(sigmaAMat))+sigmaX*as.vector(t(sigmaXMat))

## perform MLM-X on the 100 simulated genotypes
scan1 <- data.frame(scanID=1:n,sex=sex,pheno=pheno)
scan1 <- ScanAnnotationDataFrame(scan1)

genoData <- GenotypeData(genoMt,scanAnnot=scan1)

varCompBoth <- estVarComp(scan1,covMatList=covMatListBoth,"pheno")
varBothci <- estVarCompCI(varCompBoth,prop=FALSE)
varBothci$model <- "both"

cholSigBoth <- varCompBoth[["cholSigmaInv"]]
mmResBoth <- assocTestMixedModel(genoData,snpStart=snpStart,snpEnd=snpEnd,cholSigBoth,outcome="pheno")
mmResBoth$causal <- FALSE
mmResBoth$model <- "both"

varCompA <- estVarComp(scan1,covMatList=list(kinAuto),"pheno")
varAutoci <- estVarCompCI(varCompA,prop=FALSE)
varAutoci$model <- "auto"

cholSigA <- varCompA[["cholSigmaInv"]]
mmResA <- assocTestMixedModel(genoData,snpStart=snpStart,snpEnd=snpEnd,cholSigA,outcome="pheno")
mmResA$causal <- FALSE
mmResA$model <- "auto"

varCompX <- estVarComp(scan1,covMatList=list(kinX),"pheno")
varXci <- estVarCompCI(varCompX,prop=FALSE)
varXci$model <- "x"

cholSigX <- varCompX[["cholSigmaInv"]]
mmResX <- assocTestMixedModel(genoData,snpStart=snpStart,snpEnd=snpEnd,cholSigX,outcome="pheno")
mmResX$causal <- FALSE
mmResX$model <- "x"

mmResF4 <- rbind(mmResBoth,mmResA,mmResX)
mmResF4$sigma_a <- sigmaA
mmResF4$sigma_x <- sigmaX
res4 <- list(sigmaA,sigmaX,mmResF4)
names(res4) <- c("sigma_auto","sigma_x","mm_res")



## make some plots
library(ggplot2)

mmResF$sigma_a <- res[["sigma_auto"]]
mmResF$sigma_x <- res[["sigma_x"]]

tyIerr <- expand.grid(model=c("auto","both","x"),alpha=c(5e-03,1e-03),sigma_x=c(0.5,0.3),sigma_a=c(0.5,0.3))
tyIerr$est <- NA
tyIerr$lower <- NA; tyIerr$upper <- NA
for(i in 1:6){
  tyIerr$est[i] <- sum(mmResF$pval[mmResF$model==tyIerr$model[i]]<tyIerr$alpha[i])/sum(mmResF$model==tyIerr$model[i])
  stdErr <- sqrt(tyIerr$alpha[i]*(1-tyIerr$alpha[i])/sum(mmResF$model==tyIerr$model[i]))
  tyIerr$lower[i] <- tyIerr$est[i]-1.96*stdErr
  tyIerr$upper[i] <- tyIerr$est[i]+1.96*stdErr
}

for(i in 7:12){
  tyIerr$est[i] <- sum(mmResF2$pval[mmResF2$model==tyIerr$model[i]]<tyIerr$alpha[i])/sum(mmResF2$model==tyIerr$model[i])
  stdErr <- sqrt(tyIerr$alpha[i]*(1-tyIerr$alpha[i])/sum(mmResF2$model==tyIerr$model[i]))
  tyIerr$lower[i] <- tyIerr$est[i]-1.96*stdErr
  tyIerr$upper[i] <- tyIerr$est[i]+1.96*stdErr
}

for(i in 13:18){
  tyIerr$est[i] <- sum(mmResF4$pval[mmResF4$model==tyIerr$model[i]]<tyIerr$alpha[i])/sum(mmResF4$model==tyIerr$model[i])
  stdErr <- sqrt(tyIerr$alpha[i]*(1-tyIerr$alpha[i])/sum(mmResF4$model==tyIerr$model[i]))
  tyIerr$lower[i] <- tyIerr$est[i]-1.96*stdErr
  tyIerr$upper[i] <- tyIerr$est[i]+1.96*stdErr
}

for(i in 19:24){
  tyIerr$est[i] <- sum(mmResF3$pval[mmResF3$model==tyIerr$model[i]]<tyIerr$alpha[i])/sum(mmResF3$model==tyIerr$model[i])
  stdErr <- sqrt(tyIerr$alpha[i]*(1-tyIerr$alpha[i])/sum(mmResF3$model==tyIerr$model[i]))
  tyIerr$lower[i] <- tyIerr$est[i]-1.96*stdErr
  tyIerr$upper[i] <- tyIerr$est[i]+1.96*stdErr
}


totalRes <- tyIerr[tyIerr$alpha==5e-03,]
totalRes$sigma_x <- paste0("sigma_x=",totalRes$sigma_x)
totalRes$sigma_a <- paste0("sigma_a=",totalRes$sigma_a)

pdf("typeIErr_7apr_5e03.pdf")
ggplot(totalRes,aes(x=model,y=est)) +geom_point()+ geom_pointrange(aes(ymin=lower,ymax=upper))+
  theme_bw() + geom_hline(aes(yintercept=alpha,color="gray"))+
  facet_grid(sigma_x~sigma_a) +
  ylab("Type I Error Rate") + ggtitle(expression(paste("Type I Error Rate, ", alpha, "=0.005")))
dev.off()

rm(list=ls())


#####
# 29. Rerun some null simulations

# called in batch:
# cd /projects/geneva/geneva_sata/caitlin/mlm_x/
# qsub -q thornton.q -t 1-10 -N nullSims batch_null_simulations.sh

# which saves res1.RData, res2.RData, ..., res10.RData
# calls gx+ga, ga, gx, linear model for 4 combos of sigma_a, sigma_x \in 0.3, 0.5, sigma_e=1
# 15K null snps each
# saves the mm res and the var components estimates and cis

# parse the results
# resX = which sigma_x, sigma_a config
# res2_X = which iteration

## moved these results .RData files to nullSims_16ped/

library(GWASTools)
setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")

totalRes <- NULL
for(i in seq_len(10)){
  dat <- getobj(paste0("res_",i,".RData"))
  tosv <- dat[["mm_res"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  totalRes <- rbind(totalRes,tosv)

  dat <- getobj(paste0("res2_",i,".RData"))
  tosv <- dat[["mm_res"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  totalRes <- rbind(totalRes,tosv)

  dat <- getobj(paste0("res3_",i,".RData"))
  tosv <- dat[["mm_res"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  totalRes <- rbind(totalRes,tosv)

  dat <- getobj(paste0("res4_",i,".RData"))
  tosv <- dat[["mm_res"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  totalRes <- rbind(totalRes,tosv)
}

dim(totalRes) # 2400000 13
ftable(totalRes$model,totalRes$sigma_a,totalRes$sigma_x)
#                0.3    0.5
# auto   0.3  150000 150000
#        0.5  150000 150000
# both   0.3  150000 150000
#        0.5  150000 150000
# linear 0.3  150000 150000
#        0.5  150000 150000
# x      0.3  150000 150000
#        0.5  150000 150000

# great, so 150K of each type

## make some plots
library(ggplot2)
tyIerr <- expand.grid(model=c("auto","both","x","linear"),alpha=c(5e-03,5e-04),
                      sigma_x=c(0.5,0.3),sigma_a=c(0.5,0.3))
tyIerr$est <- NA
tyIerr$lower <- NA; tyIerr$upper <- NA
for(i in seq_len(nrow(tyIerr))){
  mmResF <- totalRes[totalRes$model==tyIerr$model[i]&totalRes$sigma_x==tyIerr$sigma_x[i]&
                       totalRes$sigma_a==tyIerr$sigma_a[i],]
  tyIerr$est[i] <- sum(mmResF$pval<tyIerr$alpha[i])/nrow(mmResF)
  stdErr <- sqrt(tyIerr$alpha[i]*(1-tyIerr$alpha[i])/nrow(mmResF))
  tyIerr$lower[i] <- tyIerr$est[i]-1.96*stdErr
  tyIerr$upper[i] <- tyIerr$est[i]+1.96*stdErr
}

tyIerr$sigma_a <- paste0("sigma_a=",tyIerr$sigma_a)
tyIerr$sigma_x <- paste0("sigma_x=",tyIerr$sigma_x)

# make the model an ordered factor so it prints both, x, auto, linear
tyIerr$model <- ordered(tyIerr$model,levels=c("both","x","auto","linear"))
tyIerrSm <- tyIerr[tyIerr$alpha==5e-03,]

pdf("typeIErr_5e03.pdf")
ggplot(tyIerrSm,aes(x=model,y=est)) +geom_point()+ geom_pointrange(aes(ymin=lower,ymax=upper))+
  theme_bw() + geom_hline(aes(yintercept=5e-03,color="gray"))+
  facet_grid(sigma_x~sigma_a) +
  ylab("Type I Error Rate") + ggtitle(expression(paste("Type I Error Rate, ", alpha, "=0.005")))
dev.off()

tyIerrSm <- tyIerr[tyIerr$alpha==5e-04,]

pdf("typeIErr_5e04.pdf")
ggplot(tyIerrSm,aes(x=model,y=est)) +geom_point()+ geom_pointrange(aes(ymin=lower,ymax=upper))+
  theme_bw() + geom_hline(aes(yintercept=5e-04,color="gray"))+
  facet_grid(sigma_x~sigma_a) +
  ylab("Type I Error Rate") + ggtitle(expression(paste("Type I Error Rate, ", alpha, "=5e-04")))
dev.off()

## plot the var comp estimates for these
cis_auto <- NULL
cis_x <- NULL
cis_both <- NULL
for(i in seq_len(10)){
  for(j in c("","2","3","4")){
  dat <- getobj(paste0("res",j,"_",i,".RData"))
  tosv <- dat[["ci_x"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  tosv$comp <- c("x","resid")

  cis_x <- rbind(cis_x,tosv)

  tosv <- dat[["ci_auto"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  tosv$comp <- c("auto","resid")

  cis_auto <- rbind(cis_auto,tosv)

  tosv <- dat[["ci_both"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  tosv$comp <- c("x","auto","resid")

  cis_both <- rbind(cis_both,tosv)
}}

# look at kc values comparing x to autos
kinX <- get(load("500Peds_xKinship.RData"))
kinAuto <- get(load("500Peds_autoKinship.RData"))

kinAuto <- matrix(kinAuto,nrow=8000,ncol=8000)
kinX <- matrix(kinX,nrow=8000,ncol=8000)

summary(as.vector(kinAuto[1:16,1:16]))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.0000  0.0000  0.1250  0.1353  0.2500  0.5000
summary(as.vector(kinX[1:16,1:16]))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.0000  0.0000  0.1094  0.1787  0.2500  1.0000
summary(as.vector(kinX[1:16,1:16])-as.vector(kinAuto[1:16,1:16])) # mean is >0...
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#-0.25000  0.00000  0.00000  0.04346  0.12500  0.50000

sum(kinAuto[1:16,1:16]==0) # 70
sum(kinX[1:16,1:16]==0) # 118
sum(kinAuto[1:16,1:16]==0)/length(kinAuto[1:16,1:16]) # 0.2734375
sum(kinX[1:16,1:16]==0)/length(kinAuto[1:16,1:16]) # 0.4609375
# so xKC is 46% 0's and autoKC is only 27%

plot(as.vector(kinX[1:16,1:16]),as.vector(kinAuto[1:16,1:16]),xlab="X KC",ylab="Auto KC")
abline(0,1)
dev.off()
# looks like more smaller kc values, but not true from summary?

cis_autoPl <- cis_auto[cis_auto$comp=="auto",]
mean(cis_autoPl$Est[cis_autoPl$sigma_a==0.3&cis_autoPl$sigma_x==0.3]) # 0.6990352; sum=0.6
mean(cis_autoPl$Est[cis_autoPl$sigma_a==0.3&cis_autoPl$sigma_x==0.5]) # 0.950669; sum=0.8

mean(cis_autoPl$Est[cis_autoPl$sigma_a==0.5&cis_autoPl$sigma_x==0.3]) # 0.8662233; sum=0.8
mean(cis_autoPl$Est[cis_autoPl$sigma_a==0.5&cis_autoPl$sigma_x==0.5]) # 1.130354; sum=1.0

cis_xPl <- cis_x[cis_x$comp=="x",]
mean(cis_xPl$Est[cis_xPl$sigma_a==0.3&cis_xPl$sigma_x==0.3]) # 0.447846; sum=0.6
mean(cis_xPl$Est[cis_xPl$sigma_a==0.3&cis_xPl$sigma_x==0.5]) # 0.6360744; sum=0.8

mean(cis_xPl$Est[cis_xPl$sigma_a==0.5&cis_xPl$sigma_x==0.3]) # 0.5244501; sum=0.8
mean(cis_xPl$Est[cis_xPl$sigma_a==0.5&cis_xPl$sigma_x==0.5]) # 0.712929; sum=1.0

# make histograms of the estimates of each var component, for each model type
allCis <- rbind(cis_x,cis_auto,cis_both)

allCisPl <- allCis[!is.element(allCis$comp,"resid")&allCis$sigma_a==0.3&allCis$sigma_x==0.3,]
allCisPl$model[is.element(allCisPl$model,c("x","auto"))] <- "one comp"
allCisPl$comp[allCisPl$comp=="auto"] <- "sigma^2_a"
allCisPl$comp[allCisPl$comp=="x"] <- "sigma^2_x"
pdf("hist_varComp_s2A03s2X03.pdf")
ggplot(allCisPl,aes(x=Est)) + geom_histogram() + facet_grid(model~comp) +
 ggtitle("sigma^{2}_{x}=0.3, sigma^{2}_{a}=0.3") + theme_bw()
dev.off()

allCisPl <- allCis[!is.element(allCis$comp,"resid")&allCis$sigma_a==0.5&allCis$sigma_x==0.3,]
allCisPl$model[is.element(allCisPl$model,c("x","auto"))] <- "one comp"
allCisPl$comp[allCisPl$comp=="auto"] <- "sigma^2_a"
allCisPl$comp[allCisPl$comp=="x"] <- "sigma^2_x"
pdf("hist_varComp_s2A05s2X03.pdf")
ggplot(allCisPl,aes(x=Est)) + geom_histogram() + facet_grid(model~comp) +
  ggtitle("sigma^{2}_{x}=0.3, sigma^{2}_{a}=0.5") + theme_bw()
dev.off()

allCisPl <- allCis[!is.element(allCis$comp,"resid")&allCis$sigma_a==0.3&allCis$sigma_x==0.5,]
allCisPl$model[is.element(allCisPl$model,c("x","auto"))] <- "one comp"
allCisPl$comp[allCisPl$comp=="auto"] <- "sigma^2_a"
allCisPl$comp[allCisPl$comp=="x"] <- "sigma^2_x"
pdf("hist_varComp_s2A03s2X05.pdf")
ggplot(allCisPl,aes(x=Est)) + geom_histogram() + facet_grid(model~comp) +
  ggtitle("sigma^{2}_{x}=0.5, sigma^{2}_{a}=0.3") + theme_bw()
dev.off()

allCisPl <- allCis[!is.element(allCis$comp,"resid")&allCis$sigma_a==0.5&allCis$sigma_x==0.5,]
allCisPl$model[is.element(allCisPl$model,c("x","auto"))] <- "one comp"
allCisPl$comp[allCisPl$comp=="auto"] <- "sigma^2_a"
allCisPl$comp[allCisPl$comp=="x"] <- "sigma^2_x"
pdf("hist_varComp_s2A05s2X05.pdf")
ggplot(allCisPl,aes(x=Est)) + geom_histogram() + facet_grid(model~comp) +
  ggtitle("sigma^{2}_{x}=0.5, sigma^{2}_{a}=0.5") + theme_bw()
dev.off()

# type I error is worst when ignoring the x chr effects
# var comps are ok when ignorning the x chr effects

# type I error is ok when only fitting x chr effects, but ignoring auto effects
# var comps are underestimated when only fitting x chr effects, but ignoring auto effects

rm(list=ls())


#####
# 30. Look at power graph zoomed in on x-axis

source("powerGraph.R")

totalRes <- get(load("pvals_power_results_ALL.RData"))
dim(totalRes) # 4940976 12

# be sure to exclude the NA rows.
sum(is.na(totalRes$pval)) # 0

ftable(totalRes$causal,totalRes$beta1,totalRes$model)
#               auto   both      x
# FALSE 0.03   75000  75000  75000
# 0.05  326043 324437 326514
# 0.06    7995   6142   8568
# 0.07  378012 377618 378208
# 0.08    8819   7158   9584
# 0.09    7567   5548   8799
# 0.1     8295   6293   9100
# 0.13    8949   7242   9624
# 0.15     994    994    994
# 0.18    1000   1000   1000
# 0.2     1000   1000   1000
# 0.25     997    997    997
# 0.3     1000   1000   1000
# TRUE  0.03   75000  75000  75000
# 0.05  326043 324437 326514
# 0.06    7995   6142   8568
# 0.07  378012 377618 378208
# 0.08    8819   7158   9584
# 0.09    7567   5548   8799
# 0.1     8295   6293   9100
# 0.13    8949   7242   9624
# 0.15     994    994    994
# 0.18    1000   1000   1000
# 0.2     1000   1000   1000
# 0.25     997    997    997
# 0.3     1000   1000   1000

powerGraph(totalRes,b=0.07,fn="power_3models_beta07_zoom.pdf",xlims=c(0,0.05))
powerGraph(totalRes,b=0.05,fn="power_3models_beta05_zoom.pdf",xlims=c(0,0.05))

rm(list=ls())


#####
# 31. Rerun power calculations

# qsub -q thornton.q -t 1-2000 batch_model_c_powerCalc_function_005.sh
# this will produce 10K iterations for each sigma_x, sigma_a value, beta1=0.05 from which to calculate power


#####
# 32. Simulate phenotypes with different pedigree structure

## simulate some with x KC >> auto KC
# mother-son, father-daughter relationships
# full brothers
# uncle-nephew, maternal

# make autosomal and x chr KC matrices, for 8000 people = 1000x8 person pedigree
# see script organize_kinship_matrices.R
# creates 1000Peds_8ped_xKinship.RData and 1000Peds_8ped_autoKinship.RData

# now see what the type I error and var comp estimates look like
# called in batch:
# cd /projects/geneva/geneva_sata/caitlin/mlm_x/
# qsub -q olga.q -t 1-10 -N nullSims batch_null_simulations.sh
# which calls null_simulations_v2.R (pedigree is different)

## moved these results to nullSims_8ped/

library(GWASTools)
setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")

kinX <- getobj("1000Peds_8ped_xKinship.RData")
kinAuto <- getobj("1000Peds_8ped_autoKinship.RData")

kinAuto <- matrix(kinAuto,nrow=8000,ncol=8000)
kinX <- matrix(kinX,nrow=8000,ncol=8000)

summary(as.vector(kinX[1:8,1:8]))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.0000  0.0000  0.2500  0.3438  0.5000  1.0000

summary(as.vector(kinAuto[1:8,1:8]))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.000   0.125   0.250   0.207   0.250   0.500

summary(as.vector(kinX[1:8,1:8])-as.vector(kinAuto[1:8,1:8])) # mean is >0...
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#-0.2500  0.0000  0.1250  0.1367  0.2500  0.5000

kinAuto[1:8,1:8]
kinX[1:8,1:8]

sum(kinAuto[1:8,1:8]==0) # 10
sum(kinX[1:8,1:8]==0) # 18

plot(as.vector(kinX[1:8,1:8]),as.vector(kinAuto[1:8,1:8]),xlab="X KC",ylab="Auto KC")
abline(0,1)
dev.off()
# looks like x kc is larger...

totalRes <- NULL
for(i in seq_len(10)){
  dat <- getobj(paste0("nullSims_8ped/res_",i,".RData"))
  tosv <- dat[["mm_res"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  totalRes <- rbind(totalRes,tosv)

  dat <- getobj(paste0("nullSims_8ped/res2_",i,".RData"))
  tosv <- dat[["mm_res"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  totalRes <- rbind(totalRes,tosv)

  dat <- getobj(paste0("nullSims_8ped/res3_",i,".RData"))
  tosv <- dat[["mm_res"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  totalRes <- rbind(totalRes,tosv)

  dat <- getobj(paste0("nullSims_8ped/res4_",i,".RData"))
  tosv <- dat[["mm_res"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  totalRes <- rbind(totalRes,tosv)
}

dim(totalRes) # 2400000 13
ftable(totalRes$model,totalRes$sigma_a,totalRes$sigma_x)
#                0.3    0.5
# auto   0.3  150000 150000
#        0.5  150000 150000
# both   0.3  150000 150000
#        0.5  150000 150000
# linear 0.3  150000 150000
#        0.5  150000 150000
# x      0.3  150000 150000
#        0.5  150000 150000

# great, so 150K of each type

## make some plots
library(ggplot2)
tyIerr <- expand.grid(model=c("auto","both","x","linear"),alpha=c(5e-03,5e-04,5e-05),
                      sigma_x=c(0.5,0.3),sigma_a=c(0.5,0.3))
tyIerr$est <- NA
tyIerr$lower <- NA; tyIerr$upper <- NA
for(i in seq_len(nrow(tyIerr))){
  mmResF <- totalRes[totalRes$model==tyIerr$model[i]&totalRes$sigma_x==tyIerr$sigma_x[i]&
                       totalRes$sigma_a==tyIerr$sigma_a[i],]
  tyIerr$est[i] <- sum(mmResF$pval<tyIerr$alpha[i])/nrow(mmResF)
  stdErr <- sqrt(tyIerr$alpha[i]*(1-tyIerr$alpha[i])/nrow(mmResF))
  tyIerr$lower[i] <- tyIerr$est[i]-1.96*stdErr
  tyIerr$upper[i] <- tyIerr$est[i]+1.96*stdErr
}

tyIerr$sigma_a <- paste0("sigma_a=",tyIerr$sigma_a)
tyIerr$sigma_x <- paste0("sigma_x=",tyIerr$sigma_x)

# make the model an ordered factor so it prints both, x, auto, linear
tyIerr$model <- ordered(tyIerr$model,levels=c("both","x","auto","linear"))
tyIerrSm <- tyIerr[tyIerr$alpha==5e-03,]

pdf("typeIErr_8ped_5e03.pdf")
ggplot(tyIerrSm,aes(x=model,y=est)) +geom_point()+ geom_pointrange(aes(ymin=lower,ymax=upper))+
  theme_bw() + geom_hline(aes(yintercept=5e-03,color="gray"))+
  facet_grid(sigma_x~sigma_a) +
  ylab("Type I Error Rate") + ggtitle(expression(paste("Type I Error Rate, ", alpha, "=0.005")))
dev.off()

tyIerrSm <- tyIerr[tyIerr$alpha==5e-04,]

pdf("typeIErr_8ped_5e04.pdf")
ggplot(tyIerrSm,aes(x=model,y=est)) +geom_point()+ geom_pointrange(aes(ymin=lower,ymax=upper))+
  theme_bw() + geom_hline(aes(yintercept=5e-04,color="gray"))+
  facet_grid(sigma_x~sigma_a) +
  ylab("Type I Error Rate") + ggtitle(expression(paste("Type I Error Rate, ", alpha, "=5e-04")))
dev.off()

tyIerrSm <- tyIerr[tyIerr$alpha==5e-05,]
pdf("typeIErr_8ped_5e05.pdf")
ggplot(tyIerrSm,aes(x=model,y=est)) +geom_point()+ geom_pointrange(aes(ymin=lower,ymax=upper))+
  theme_bw() + geom_hline(aes(yintercept=5e-05,color="gray"))+
  facet_grid(sigma_x~sigma_a) +
  theme(axis.text=element_text(size=16),axis.title=element_text(size=18),title=element_text(size=20),
                   legend.text=element_text(size=16),strip.text = element_text(size=16)) +
  ylab("Type I Error Rate") + ggtitle(expression(paste("Type I Error Rate, ", alpha, "=5e-05")))
dev.off()

## extreme sigma, now
tyIerr <- expand.grid(model=c("auto","both","x","linear"),alpha=c(5e-03,5e-04,5e-05),
                      sigma_x=c(0.3,3),sigma_a=c(0.3,3))
tyIerr$est <- NA
tyIerr$lower <- NA; tyIerr$upper <- NA
for(i in seq_len(nrow(tyIerr))){
  mmResF <- totalRes[totalRes$model==tyIerr$model[i]&totalRes$sigma_x==tyIerr$sigma_x[i]&
                       totalRes$sigma_a==tyIerr$sigma_a[i],]
  tyIerr$est[i] <- sum(mmResF$pval<tyIerr$alpha[i])/nrow(mmResF)
  stdErr <- sqrt(tyIerr$alpha[i]*(1-tyIerr$alpha[i])/nrow(mmResF))
  tyIerr$lower[i] <- tyIerr$est[i]-1.96*stdErr
  tyIerr$upper[i] <- tyIerr$est[i]+1.96*stdErr
}

tyIerr$sigma_a <- paste0("sigma_a=",tyIerr$sigma_a)
tyIerr$sigma_x <- paste0("sigma_x=",tyIerr$sigma_x)

tyIerr <- tyIerr[!is.nan(tyIerr$est),]
tyIerr <- tyIerr[!(tyIerr$sigma_x=="sigma_x=0.3"&tyIerr$sigma_a=="sigma_a=0.3"),]

# make the model an ordered factor so it prints both, x, auto, linear
tyIerr$model <- ordered(tyIerr$model,levels=c("both","x","auto","linear"))
tyIerrSm <- tyIerr[tyIerr$alpha==5e-05,]

pdf("typeIErr_8ped_5e05_extremeSigma.pdf")
ggplot(tyIerrSm,aes(x=model,y=est)) +geom_point()+ geom_pointrange(aes(ymin=lower,ymax=upper))+
  theme_bw() + geom_hline(aes(yintercept=5e-05,color="gray"))+
  facet_wrap(sigma_x~sigma_a) +
  theme(axis.text=element_text(size=16),axis.title=element_text(size=18),title=element_text(size=20),
          legend.text=element_text(size=16),strip.text = element_text(size=16)) +
  ylab("Type I Error Rate") + ggtitle(expression(paste("Type I Error Rate, ", alpha, "=5e-05")))
dev.off()



## plot the var comp estimates for these
cis_auto <- NULL
cis_x <- NULL
cis_both <- NULL
for(i in seq_len(10)){
  for(j in c("","2","3","4")){
    dat <- getobj(paste0("nullSims_8ped/res",j,"_",i,".RData"))
    tosv <- dat[["ci_x"]]
    tosv$sigma_a <- dat[["sigma_auto"]]
    tosv$sigma_x <- dat[["sigma_x"]]
    tosv$comp <- c("x","resid")

    cis_x <- rbind(cis_x,tosv)

    tosv <- dat[["ci_auto"]]
    tosv$sigma_a <- dat[["sigma_auto"]]
    tosv$sigma_x <- dat[["sigma_x"]]
    tosv$comp <- c("auto","resid")

    cis_auto <- rbind(cis_auto,tosv)

    tosv <- dat[["ci_both"]]
    tosv$sigma_a <- dat[["sigma_auto"]]
    tosv$sigma_x <- dat[["sigma_x"]]
    tosv$comp <- c("x","auto","resid")

    cis_both <- rbind(cis_both,tosv)
  }}

cis_autoPl <- cis_auto[cis_auto$comp=="auto",]
mean(cis_autoPl$Est[cis_autoPl$sigma_a==0.3&cis_autoPl$sigma_x==0.3]) # 0.6817052; sum=0.6
mean(cis_autoPl$Est[cis_autoPl$sigma_a==0.3&cis_autoPl$sigma_x==0.5]) # 0.9328929; sum=0.8

mean(cis_autoPl$Est[cis_autoPl$sigma_a==0.5&cis_autoPl$sigma_x==0.3]) # 0.8571642; sum=0.8
mean(cis_autoPl$Est[cis_autoPl$sigma_a==0.5&cis_autoPl$sigma_x==0.5]) # 1.106899; sum=1.0

cis_autoPl <- cis_auto[cis_auto$comp=="resid",]
mean(cis_autoPl$Est[cis_autoPl$sigma_a==0.3&cis_autoPl$sigma_x==0.3]) # 1.070525; true=1.0
mean(cis_autoPl$Est[cis_autoPl$sigma_a==0.3&cis_autoPl$sigma_x==0.5]) # 1.108168; true=1.0

mean(cis_autoPl$Est[cis_autoPl$sigma_a==0.5&cis_autoPl$sigma_x==0.3]) # 1.073937; true=1.0
mean(cis_autoPl$Est[cis_autoPl$sigma_a==0.5&cis_autoPl$sigma_x==0.5]) # 1.112177; true=1.0



cis_xPl <- cis_x[cis_x$comp=="x",]
mean(cis_xPl$Est[cis_xPl$sigma_a==0.3&cis_xPl$sigma_x==0.3]) # 0.4301809; sum=0.6
mean(cis_xPl$Est[cis_xPl$sigma_a==0.3&cis_xPl$sigma_x==0.5]) # 0.6184052; sum=0.8

mean(cis_xPl$Est[cis_xPl$sigma_a==0.5&cis_xPl$sigma_x==0.3]) # 0.4965284; sum=0.8
mean(cis_xPl$Est[cis_xPl$sigma_a==0.5&cis_xPl$sigma_x==0.5]) # 0.6948854; sum=1.0

cis_xPl <- cis_x[cis_x$comp=="resid",]
mean(cis_xPl$Est[cis_xPl$sigma_a==0.3&cis_xPl$sigma_x==0.3]) # 1.036361; true=1.0
mean(cis_xPl$Est[cis_xPl$sigma_a==0.3&cis_xPl$sigma_x==0.5]) # 1.037935; true=1.0

mean(cis_xPl$Est[cis_xPl$sigma_a==0.5&cis_xPl$sigma_x==0.3]) # 1.069274; true=1.0
mean(cis_xPl$Est[cis_xPl$sigma_a==0.5&cis_xPl$sigma_x==0.5]) # 1.062538; true=1.0



cis_bothPl <- cis_both[cis_both$comp=="x",]
mean(cis_bothPl$Est[cis_bothPl$sigma_a==0.3&cis_bothPl$sigma_x==0.3]) # 0.3358279; true=0.3
mean(cis_bothPl$Est[cis_bothPl$sigma_a==0.3&cis_bothPl$sigma_x==0.5]) # 0.5062398; true=0.5

mean(cis_bothPl$Est[cis_bothPl$sigma_a==0.5&cis_bothPl$sigma_x==0.3]) # 0.3241463; true=0.3
mean(cis_bothPl$Est[cis_bothPl$sigma_a==0.5&cis_bothPl$sigma_x==0.5]) # 0.5173022; true=0.5

cis_bothPl <- cis_both[cis_both$comp=="auto",]
mean(cis_bothPl$Est[cis_bothPl$sigma_a==0.3&cis_bothPl$sigma_x==0.3]) # 0.2394941; true=0.3
mean(cis_bothPl$Est[cis_bothPl$sigma_a==0.3&cis_bothPl$sigma_x==0.5]) # 0.2912639; true=0.3

mean(cis_bothPl$Est[cis_bothPl$sigma_a==0.5&cis_bothPl$sigma_x==0.3]) # 0.4447845; true=0.5
mean(cis_bothPl$Est[cis_bothPl$sigma_a==0.5&cis_bothPl$sigma_x==0.5]) # 0.4630382; true=0.5

cis_bothPl <- cis_both[cis_both$comp=="resid",]
mean(cis_bothPl$Est[cis_bothPl$sigma_a==0.3&cis_bothPl$sigma_x==0.3]) # 1.000223; true=1.0
mean(cis_bothPl$Est[cis_bothPl$sigma_a==0.3&cis_bothPl$sigma_x==0.5]) # 0.9922439; true=1.0

mean(cis_bothPl$Est[cis_bothPl$sigma_a==0.5&cis_bothPl$sigma_x==0.3]) # 1.000157; true=1.0
mean(cis_bothPl$Est[cis_bothPl$sigma_a==0.5&cis_bothPl$sigma_x==0.5]) # 0.9894086; true=1.0


# make histograms of the estimates of each var component, for each model type
allCis <- rbind(cis_x,cis_auto,cis_both)

allCisPl <- allCis[!is.element(allCis$comp,"resid")&allCis$sigma_a==0.3&allCis$sigma_x==0.3,]
allCisPl$model[is.element(allCisPl$model,c("x","auto"))] <- "one comp"
allCisPl$comp[allCisPl$comp=="auto"] <- "sigma^2_a"
allCisPl$comp[allCisPl$comp=="x"] <- "sigma^2_x"
pdf("hist_varComp_8ped_s2A03s2X03.pdf")
ggplot(allCisPl,aes(x=Est)) + geom_histogram() + facet_grid(model~comp) +
  ggtitle("sigma^{2}_{x}=0.3, sigma^{2}_{a}=0.3") + theme_bw()
dev.off()

allCisPl <- allCis[!is.element(allCis$comp,"resid")&allCis$sigma_a==0.5&allCis$sigma_x==0.3,]
allCisPl$model[is.element(allCisPl$model,c("x","auto"))] <- "one comp"
allCisPl$comp[allCisPl$comp=="auto"] <- "sigma^2_a"
allCisPl$comp[allCisPl$comp=="x"] <- "sigma^2_x"
pdf("hist_varComp_8ped_s2A05s2X03.pdf")
ggplot(allCisPl,aes(x=Est)) + geom_histogram() + facet_grid(model~comp) +
  ggtitle("sigma^{2}_{x}=0.3, sigma^{2}_{a}=0.5") + theme_bw()
dev.off()

allCisPl <- allCis[!is.element(allCis$comp,"resid")&allCis$sigma_a==0.3&allCis$sigma_x==0.5,]
allCisPl$model[is.element(allCisPl$model,c("x","auto"))] <- "one comp"
allCisPl$comp[allCisPl$comp=="auto"] <- "sigma^2_a"
allCisPl$comp[allCisPl$comp=="x"] <- "sigma^2_x"
pdf("hist_varComp_8ped_s2A03s2X05.pdf")
ggplot(allCisPl,aes(x=Est)) + geom_histogram() + facet_grid(model~comp) +
  ggtitle("sigma^2[x]=0.5, sigma^2[a]=0.3") + theme_bw()
dev.off()

allCisPl <- allCis[!is.element(allCis$comp,"resid")&allCis$sigma_a==0.5&allCis$sigma_x==0.5,]
allCisPl$model[is.element(allCisPl$model,c("x","auto"))] <- "one comp"
allCisPl$comp[allCisPl$comp=="auto"] <- "sigma^2_a"
allCisPl$comp[allCisPl$comp=="x"] <- "sigma^2_x"
pdf("hist_varComp_8ped_s2A05s2X05.pdf")
ggplot(allCisPl,aes(x=Est)) + geom_histogram() + facet_grid(model~comp) +
  ggtitle("sigma^{2}_{x}=0.5, sigma^{2}_{a}=0.5") + theme_bw()
dev.off()

# var comps are ok when ignorning the x chr effects
# type I error is worst when ignoring the x chr effects

# var comps are underestimated when only fitting x chr effects, but ignoring auto effects
# type I error is ok when only fitting x chr effects, but ignoring auto effects

rm(list=ls())


#####
# 33. Look at new power results

#modelC_results_power_sA03sX03
#modelC_results_power_sA05sX03
#modelC_results_power_sA03sX05
#modelC_results_power_sA05sX05/mmRes_SNP1992_beta0.05.RData

# called from batch:
# cd /projects/geneva/geneva_sata/caitlin/mlm_x/olga_application
# qsub -q thornton.q -N readRes batch_read_results.sh

source("powerGraph.R")

dat <- get(load("pvals_newPower_results_apr10.RData"))
newcs <- get(load("pvals_newPower_results_apr14.RData"))
dat <- rbind(dat,newcs)
dim(dat) # 493350 14
table(dat$beta1) # 480K 0.05, 13350 0.07
dat <- dat[dat$beta1==0.05,]

toPl <- dat[dat$sigma_x=="03"&dat$sigma_a=="03"&!is.na(dat$pval),]
dim(toPl) # 120000 14
table(toPl$causal) # 60000 true
sa3sx3=powerGraph(toPl,b=0.05,fn="powerGraph_sA03sX03.pdf",xlims=c(0,0.05),values=TRUE)
sa3sx3$sigma_a=paste0("sigma_a=",0.3)
sa3sx3$sigma_x=paste0("sigma_x=",0.3)

toPl <- dat[dat$sigma_x=="03"&dat$sigma_a=="05"&!is.na(dat$pval),]
dim(toPl) # 120K 14
table(toPl$causal) # 60K true
sa5sx3=powerGraph(toPl,b=0.05,fn="powerGraph_sA05sX03.pdf",xlims=c(0,0.05),values=TRUE)
sa5sx3$sigma_a=paste0("sigma_a=",0.5)
sa5sx3$sigma_x=paste0("sigma_x=",0.3)

toPl <- dat[dat$sigma_x=="05"&dat$sigma_a=="03"&!is.na(dat$pval),]
dim(toPl) # 120K 14
table(toPl$causal) # 60K
sa3sx5=powerGraph(toPl,b=0.05,fn="powerGraph_sA03sX05.pdf",xlims=c(0,0.05),values=TRUE)
sa3sx5$sigma_a=paste0("sigma_a=",0.3)
sa3sx5$sigma_x=paste0("sigma_x=",0.5)

toPl <- dat[dat$sigma_x=="05"&dat$sigma_a=="05"&!is.na(dat$pval),]
dim(toPl) # 120K 14
table(toPl$causal) # 60K t
sa5sx5=powerGraph(toPl,b=0.05,fn="powerGraph_sA05sX05.pdf",xlims=c(0,0.05),values=TRUE)
sa5sx5$sigma_a=paste0("sigma_a=",0.5)
sa5sx5$sigma_x=paste0("sigma_x=",0.5)

library(ggplot2); library(reshape)
tot <- rbind(sa3sx3,sa3sx5,sa5sx3,sa5sx5)
tot$x_fp_rate <- tot$x_fp/tot$xiter
tot$x_tp_rate <- tot$x_tp/tot$xiter
tot$auto_fp_rate <- tot$auto_fp/tot$autoiter
tot$auto_tp_rate <- tot$auto_tp/tot$autoiter
tot$both_fp_rate <- tot$both_fp/tot$bothiter
tot$both_tp_rate <- tot$both_tp/tot$bothiter

both <- tot[,c("sigma_a","sigma_x","both_fp_rate","both_tp_rate")]
both$model <- "both"
both$fp <- both$both_fp_rate
both$tp <- both$both_tp_rate

auto <- tot[,c("sigma_a","sigma_x","auto_fp_rate","auto_tp_rate")]
auto$model <- "auto"
auto$fp <- auto$auto_fp_rate
auto$tp <- auto$auto_tp_rate

x <- tot[,c("sigma_a","sigma_x","x_fp_rate","x_tp_rate")]
x$model <- "x"
x$fp <- x$x_fp_rate
x$tp <- x$x_tp_rate

toPl <- rbind(x[,c("sigma_x","sigma_a","model","tp","fp")],
              auto[,c("sigma_x","sigma_a","model","tp","fp")],
              both[,c("sigma_x","sigma_a","model","tp","fp")])

toPl <- toPl[toPl$fp<=0.05,]

pdf("power_beta05.pdf")
ggplot(toPl,aes(x=fp,y=tp,color=model)) + geom_line(size=1.5,aes(linetype=model)) + facet_grid(sigma_a~sigma_x) + theme_bw() +
  ylab("True Positive Rate") + xlab("False Positive Rate") +ggtitle(expression(paste(beta,"=0.05")))+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=18),title=element_text(size=20),
       legend.text=element_text(size=16),strip.text = element_text(size=16))
dev.off()

rm(list=ls())


#####
# 34. New null simulations with extreme sigma values

# now see what the type I error and var comp estimates look like
# called in batch:
# cd /projects/geneva/geneva_sata/caitlin/mlm_x/
# qsub -q olga.q -t 1-10 -N nullSims batch_null_simulations.sh
# which calls null_simulations_v3.R (pedigree is different)

## moved these results to nullSims_8ped/

library(GWASTools)
setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")

totalRes <- NULL
for(i in seq_len(10)){
  dat <- getobj(paste0("nullSims_8ped/res_",i,".RData"))
  tosv <- dat[["mm_res"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  totalRes <- rbind(totalRes,tosv)

  dat <- getobj(paste0("nullSims_8ped/res2_",i,".RData"))
  tosv <- dat[["mm_res"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  totalRes <- rbind(totalRes,tosv)

  dat <- getobj(paste0("nullSims_8ped/res3_",i,".RData"))
  tosv <- dat[["mm_res"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  totalRes <- rbind(totalRes,tosv)

  dat <- getobj(paste0("nullSims_8ped/res4_",i,".RData"))
  tosv <- dat[["mm_res"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  totalRes <- rbind(totalRes,tosv)

  dat <- getobj(paste0("nullSims_8ped/extreme_sigma/res_",i,".RData"))
  tosv <- dat[["mm_res"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  totalRes <- rbind(totalRes,tosv)

  dat <- getobj(paste0("nullSims_8ped/extreme_sigma/res2_",i,".RData"))
  tosv <- dat[["mm_res"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  totalRes <- rbind(totalRes,tosv)
}

save(totalRes,file="nullSims_8ped/mmRes_nullSims.RData")

dim(totalRes) # 3600000 13
ftable(totalRes$model,totalRes$sigma_a,totalRes$sigma_x)
#               0.3    0.5      3
# auto   0.3  150000 150000 150000
# 0.5  150000 150000      0
# 3    150000      0      0
# both   0.3  150000 150000 150000
# 0.5  150000 150000      0
# 3    150000      0      0
# linear 0.3  150000 150000 150000
# 0.5  150000 150000      0
# 3    150000      0      0
# x      0.3  150000 150000 150000
# 0.5  150000 150000      0
# 3    150000      0      0

# great, so 150K of each type

## make some plots
library(ggplot2)
tyIerr <- expand.grid(model=c("auto","both","x","linear"),alpha=c(5e-03,5e-04),
                      sigma_x=c(3,0.5,0.3),sigma_a=c(3,0.5,0.3))
tyIerr$est <- NA
tyIerr$lower <- NA; tyIerr$upper <- NA
for(i in seq_len(nrow(tyIerr))){
  mmResF <- totalRes[totalRes$model==tyIerr$model[i]&totalRes$sigma_x==tyIerr$sigma_x[i]&
                       totalRes$sigma_a==tyIerr$sigma_a[i],]
  tyIerr$est[i] <- sum(mmResF$pval<tyIerr$alpha[i])/nrow(mmResF)
  stdErr <- sqrt(tyIerr$alpha[i]*(1-tyIerr$alpha[i])/nrow(mmResF))
  tyIerr$lower[i] <- tyIerr$est[i]-1.96*stdErr
  tyIerr$upper[i] <- tyIerr$est[i]+1.96*stdErr
}

tyIerr$sigma_a <- paste0("sigma_a=",tyIerr$sigma_a)
tyIerr$sigma_x <- paste0("sigma_x=",tyIerr$sigma_x)

# make the model an ordered factor so it prints both, x, auto, linear
tyIerr$model <- ordered(tyIerr$model,levels=c("both","x","auto","linear"))
tyIerrSm <- tyIerr[tyIerr$alpha==5e-03,]

pdf("typeIErr_8ped_5e03_inclExtreme.pdf")
ggplot(tyIerrSm,aes(x=model,y=est)) +geom_point()+ geom_pointrange(aes(ymin=lower,ymax=upper))+
  theme_bw() + geom_hline(aes(yintercept=5e-03,color="gray"))+
  facet_grid(sigma_x~sigma_a) +
  ylab("Type I Error Rate") + ggtitle(expression(paste("Type I Error Rate, ", alpha, "=0.005")))
dev.off()


## plot the var comp estimates for these
cis_auto <- NULL
cis_x <- NULL
cis_both <- NULL
for(i in seq_len(10)){
  for(j in c("res","res2","res3","res4","extreme_sigma/res","extreme_sigma/res2")){
    dat <- getobj(paste0("nullSims_8ped/",j,"_",i,".RData"))
    tosv <- dat[["ci_x"]]
    tosv$sigma_a <- dat[["sigma_auto"]]
    tosv$sigma_x <- dat[["sigma_x"]]
    tosv$comp <- c("x","resid")

    cis_x <- rbind(cis_x,tosv)

    tosv <- dat[["ci_auto"]]
    tosv$sigma_a <- dat[["sigma_auto"]]
    tosv$sigma_x <- dat[["sigma_x"]]
    tosv$comp <- c("auto","resid")

    cis_auto <- rbind(cis_auto,tosv)

    tosv <- dat[["ci_both"]]
    tosv$sigma_a <- dat[["sigma_auto"]]
    tosv$sigma_x <- dat[["sigma_x"]]
    tosv$comp <- c("x","auto","resid")

    cis_both <- rbind(cis_both,tosv)
  }}

vi <- get(load("varEst_newPower_results_apr10.RData"))

vi$comp <- NA
vi$comp[vi$model=="auto"] <- rep(c("auto","resid"),times=sum(vi$model=="auto")/2)
vi$comp[vi$model=="x"] <- rep(c("x","resid"),times=sum(vi$model=="x")/2)
vi$comp[vi$model=="both"] <- rep(c("x","auto","resid"),times=sum(vi$model=="both")/3)

vi$sigma_a[vi$sigma_a=="03"] <- 0.3
vi$sigma_x[vi$sigma_x=="03"] <- 0.3
vi$sigma_a[vi$sigma_a=="05"] <- 0.5
vi$sigma_x[vi$sigma_x=="05"] <- 0.5

cis_both <- rbind(cis_both,vi[vi$model=="both",c("Est","Lower 95","Upper 95","model","sigma_a","sigma_x","comp")])
cis_x <- rbind(cis_x,vi[vi$model=="x",c("Est","Lower 95","Upper 95","model","sigma_a","sigma_x","comp")])
cis_auto <- rbind(cis_auto,vi[vi$model=="auto",c("Est","Lower 95","Upper 95","model","sigma_a","sigma_x","comp")])

cis <- list("cis_x"=cis_x,"cis_auto"=cis_auto,"cis_both"=cis_both)
save(cis,file="nullSims_8ped/varCompCIs_nullSims.RData")

ftable(cis_auto$comp,cis_auto$sigma_a,cis_auto$sigma_x)
#              0.3   0.5     3
# auto  0.3  10010 10010    10
#       0.5  10010 10010     0
#       3       10     0     0
# resid 0.3  10010 10010    10
#       0.5  10010 10010     0
#       3       10     0     0

cis_autoPl <- cis_auto[cis_auto$comp=="auto",]
mean(cis_autoPl$Est[cis_autoPl$sigma_a==0.3&cis_autoPl$sigma_x==0.3]) # 0.6817052; sum=0.6
mean(cis_autoPl$Est[cis_autoPl$sigma_a==0.3&cis_autoPl$sigma_x==0.5]) # 0.9328929; sum=0.8

mean(cis_autoPl$Est[cis_autoPl$sigma_a==0.5&cis_autoPl$sigma_x==0.3]) # 0.8571642; sum=0.8
mean(cis_autoPl$Est[cis_autoPl$sigma_a==0.5&cis_autoPl$sigma_x==0.5]) # 1.106899; sum=1.0

mean(cis_autoPl$Est[cis_autoPl$sigma_a==0.3&cis_autoPl$sigma_x==3]) # 3.705486; sum=3.3
mean(cis_autoPl$Est[cis_autoPl$sigma_a==3&cis_autoPl$sigma_x==0.3]) # 3.262381; sum=3.3


cis_autoPl <- cis_auto[cis_auto$comp=="resid",]
mean(cis_autoPl$Est[cis_autoPl$sigma_a==0.3&cis_autoPl$sigma_x==0.3]) # 1.070525; true=1.0
mean(cis_autoPl$Est[cis_autoPl$sigma_a==0.3&cis_autoPl$sigma_x==0.5]) # 1.108168; true=1.0

mean(cis_autoPl$Est[cis_autoPl$sigma_a==0.5&cis_autoPl$sigma_x==0.3]) # 1.073937; true=1.0
mean(cis_autoPl$Est[cis_autoPl$sigma_a==0.5&cis_autoPl$sigma_x==0.5]) # 1.112177; true=1.0

mean(cis_autoPl$Est[cis_autoPl$sigma_a==0.3&cis_autoPl$sigma_x==3]) # 1.806813; true=1.0
mean(cis_autoPl$Est[cis_autoPl$sigma_a==3&cis_autoPl$sigma_x==0.3]) # 1.106329; true=1.0



cis_xPl <- cis_x[cis_x$comp=="x",]
mean(cis_xPl$Est[cis_xPl$sigma_a==0.3&cis_xPl$sigma_x==0.3]) # 0.4301809; sum=0.6
mean(cis_xPl$Est[cis_xPl$sigma_a==0.3&cis_xPl$sigma_x==0.5]) # 0.6184052; sum=0.8

mean(cis_xPl$Est[cis_xPl$sigma_a==0.5&cis_xPl$sigma_x==0.3]) # 0.4965284; sum=0.8
mean(cis_xPl$Est[cis_xPl$sigma_a==0.5&cis_xPl$sigma_x==0.5]) # 0.6948854; sum=1.0

mean(cis_xPl$Est[cis_xPl$sigma_a==3&cis_xPl$sigma_x==0.3]) # 1.404451; sum=3.3
mean(cis_xPl$Est[cis_xPl$sigma_a==0.3&cis_xPl$sigma_x==3]) # 3.093944; sum=3.3


cis_xPl <- cis_x[cis_x$comp=="resid",]
mean(cis_xPl$Est[cis_xPl$sigma_a==0.3&cis_xPl$sigma_x==0.3]) # 1.036361; true=1.0
mean(cis_xPl$Est[cis_xPl$sigma_a==0.3&cis_xPl$sigma_x==0.5]) # 1.037935; true=1.0

mean(cis_xPl$Est[cis_xPl$sigma_a==0.5&cis_xPl$sigma_x==0.3]) # 1.069274; true=1.0
mean(cis_xPl$Est[cis_xPl$sigma_a==0.5&cis_xPl$sigma_x==0.5]) # 1.062538; true=1.0

mean(cis_xPl$Est[cis_xPl$sigma_a==0.3&cis_xPl$sigma_x==3]) # 1.035329; true=1.0
mean(cis_xPl$Est[cis_xPl$sigma_a==3&cis_xPl$sigma_x==0.3]) # 1.49195; true=1.0


ftable(cis_both$comp,cis_both$sigma_a,cis_both$sigma_x)
#              0.3   0.5     3
# auto  0.3  10010 10010    10
#       0.5  10010 10010     0
#       3       10     0     0
# resid 0.3  10010 10010    10
#       0.5  10010 10010     0
#       3       10     0     0
# x     0.3  10010 10010    10
#       0.5  10010 10010     0
#       3       10     0     0

cis_bothPl <- cis_both[cis_both$comp=="x",]
mean(cis_bothPl$Est[cis_bothPl$sigma_a==0.3&cis_bothPl$sigma_x==0.3]) # 0.300924; true=0.3
mean(cis_bothPl$Est[cis_bothPl$sigma_a==0.3&cis_bothPl$sigma_x==0.5]) # 0.5008091; true=0.5

mean(cis_bothPl$Est[cis_bothPl$sigma_a==0.5&cis_bothPl$sigma_x==0.3]) # 0.3008164; true=0.3
mean(cis_bothPl$Est[cis_bothPl$sigma_a==0.5&cis_bothPl$sigma_x==0.5]) # 0.5022274; true=0.5

mean(cis_bothPl$Est[cis_bothPl$sigma_a==0.3&cis_bothPl$sigma_x==3]) # 2.997379; true=3
mean(cis_bothPl$Est[cis_bothPl$sigma_a==3&cis_bothPl$sigma_x==0.3]) # 0.3256127; true=0.3

cis_bothPl <- cis_both[cis_both$comp=="auto",]
mean(cis_bothPl$Est[cis_bothPl$sigma_a==0.3&cis_bothPl$sigma_x==0.3]) # 0.304379; true=0.3
mean(cis_bothPl$Est[cis_bothPl$sigma_a==0.3&cis_bothPl$sigma_x==0.5]) # 0.2984428; true=0.3

mean(cis_bothPl$Est[cis_bothPl$sigma_a==0.5&cis_bothPl$sigma_x==0.3]) # 0.5038336; true=0.5
mean(cis_bothPl$Est[cis_bothPl$sigma_a==0.5&cis_bothPl$sigma_x==0.5]) # 0.4988112; true=0.5

mean(cis_bothPl$Est[cis_bothPl$sigma_a==0.3&cis_bothPl$sigma_x==3]) # 0.2634433; true=0.3
mean(cis_bothPl$Est[cis_bothPl$sigma_a==3&cis_bothPl$sigma_x==0.3]) # 2.899841; true=3

cis_bothPl <- cis_both[cis_both$comp=="resid",]
mean(cis_bothPl$Est[cis_bothPl$sigma_a==0.3&cis_bothPl$sigma_x==0.3]) # 0.9984585; true=1.0
mean(cis_bothPl$Est[cis_bothPl$sigma_a==0.3&cis_bothPl$sigma_x==0.5]) # 1.000878; true=1.0

mean(cis_bothPl$Est[cis_bothPl$sigma_a==0.5&cis_bothPl$sigma_x==0.3]) # 0.9995569; true=1.0
mean(cis_bothPl$Est[cis_bothPl$sigma_a==0.5&cis_bothPl$sigma_x==0.5]) # 1.000251; true=1.0

mean(cis_bothPl$Est[cis_bothPl$sigma_a==0.3&cis_bothPl$sigma_x==3]) # 0.9909418; true=1.0
mean(cis_bothPl$Est[cis_bothPl$sigma_a==3&cis_bothPl$sigma_x==0.3]) # 1.013769; true=1.0


# make histograms of the estimates of each var component, for each model type
allCis <- rbind(cis_x,cis_auto,cis_both)

# first plot all sigma_resid, since true value is always=1
justE <- allCis[allCis$comp=="resid",]
dim(justE) # 180 7
justE$ind <- 1:nrow(justE)
justE$sigma_a <- paste0("sigma_a=",justE$sigma_a)
justE$sigma_x <- paste0("sigma_x=",justE$sigma_x)
justE$params <- paste(justE$sigma_a,justE$sigma_x,sep=", ")

pdf("sigma_e_estimates.pdf",width=11,height=9)
ggplot(justE,aes(x=ind,y=Est,color=model)) +geom_point(aes(shape=params)) +
  geom_hline(yintercept=1) + xlab("Iteration") + ylab("Estimate") + theme_bw()
dev.off()

# now plot sigma_a; exclude the extremes
justA <- allCis[allCis$comp=="auto"&!is.element(allCis$sigma_a,3)&!is.element(allCis$sigma_x,3),]
dim(justA) # 80 7
justA <- justA[order(justA$model,justA$sigma_a,justA$sigma_x),]
justA$ind <- 1:nrow(justA)
justA$sa <- justA$sigma_a
justA$sigma_a <- paste0("sigma_a=",justA$sigma_a)
justA$sigma_x <- paste0("sigma_x=",justA$sigma_x)
justA$params <- paste(justA$sigma_a,justA$sigma_x,sep=", ")

pdf("sigma_a_estimates.pdf",width=11,height=9)
ggplot(justA,aes(x=ind,y=Est,color=model)) +geom_point() + ggtitle(bquote(sigma[a]^{2}~Estimates)) +
   xlab("Iteration") + ylab("Estimate") + theme_bw() + facet_grid(sigma_x~sigma_a) +
  geom_hline(aes(yintercept=sa))
dev.off()

allCisPl <- allCis[!is.element(allCis$comp,"resid")&allCis$sigma_a==0.3&allCis$sigma_x==0.3,]
allCisPl$model[is.element(allCisPl$model,c("x","auto"))] <- "one comp"
allCisPl$comp[allCisPl$comp=="auto"] <- "sigma^2_a"
allCisPl$comp[allCisPl$comp=="x"] <- "sigma^2_x"
pdf("hist_varComp_8ped_s2A03s2X03.pdf")
ggplot(allCisPl,aes(x=Est)) + geom_histogram() + facet_grid(model~comp) +
  ggtitle("sigma^{2}_{x}=0.3, sigma^{2}_{a}=0.3") + theme_bw()
dev.off()

allCisPl <- allCis[!is.element(allCis$comp,"resid")&allCis$sigma_a==0.5&allCis$sigma_x==0.3,]
allCisPl$model[is.element(allCisPl$model,c("x","auto"))] <- "one comp"
allCisPl$comp[allCisPl$comp=="auto"] <- "sigma^2_a"
allCisPl$comp[allCisPl$comp=="x"] <- "sigma^2_x"
pdf("hist_varComp_8ped_s2A05s2X03.pdf")
ggplot(allCisPl,aes(x=Est)) + geom_histogram() + facet_grid(model~comp) +
  ggtitle("sigma^{2}_{x}=0.3, sigma^{2}_{a}=0.5") + theme_bw()
dev.off()

allCisPl <- allCis[!is.element(allCis$comp,"resid")&allCis$sigma_a==0.3&allCis$sigma_x==0.5,]
allCisPl$model[is.element(allCisPl$model,c("x","auto"))] <- "one comp"
allCisPl$comp[allCisPl$comp=="auto"] <- "sigma^2_a"
allCisPl$comp[allCisPl$comp=="x"] <- "sigma^2_x"
pdf("hist_varComp_8ped_s2A03s2X05.pdf")
ggplot(allCisPl,aes(x=Est)) + geom_histogram() + facet_grid(model~comp) +
  ggtitle("sigma^2[x]=0.5, sigma^2[a]=0.3") + theme_bw()
dev.off()

allCisPl <- allCis[!is.element(allCis$comp,"resid")&allCis$sigma_a==0.5&allCis$sigma_x==0.5,]
allCisPl$model[is.element(allCisPl$model,c("x","auto"))] <- "one comp"
allCisPl$comp[allCisPl$comp=="auto"] <- "sigma^2_a"
allCisPl$comp[allCisPl$comp=="x"] <- "sigma^2_x"
pdf("hist_varComp_8ped_s2A05s2X05.pdf")
ggplot(allCisPl,aes(x=Est)) + geom_histogram() + facet_grid(model~comp) +
  ggtitle("sigma^{2}_{x}=0.5, sigma^{2}_{a}=0.5") + theme_bw()
dev.off()

# var comps are ok when ignorning the x chr effects
# type I error is worst when ignoring the x chr effects

# var comps are underestimated when only fitting x chr effects, but ignoring auto effects
# type I error is ok when only fitting x chr effects, but ignoring auto effects

rm(list=ls())


#####
# 35. Boxplots of var comp estimates

# added new iterations:
# cd /projects/geneva/geneva_sata/caitlin/mlm_x
# qsub -t 1-100 -N nullSims batch_null_simulations.sh
# which calls null_simulations_v4.R
# pheno only has g_a effects, and sigma values are extreme

library(GWASTools)
setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")

cis <- getobj("nullSims_8ped/varCompCIs_nullSims.RData")

newcs <- getobj("varEst_newPower_results_apr14.RData")
dim(newcs); head(newcs) # 295575 7

newcs$sigma_x[newcs$sigma_x=="03"] <- 0.3
newcs$sigma_x[newcs$sigma_x=="05"] <- 0.5
newcs$sigma_a[newcs$sigma_a=="03"] <- 0.3
newcs$sigma_a[newcs$sigma_a=="05"] <- 0.5

cib <- newcs[newcs$model=="both",]
cix <- newcs[newcs$model=="x",]
cia <- newcs[newcs$model=="auto",]

cib$comp <- rep(c("x","auto","resid"),times=nrow(cib)/3)
cix$comp <- rep(c("x","resid"),times=nrow(cix)/2)
cia$comp <- rep(c("auto","resid"),times=nrow(cia)/2)

head(cis[[1]])
cis[[1]] <- rbind(cis[[1]],cix[,c("Est","Lower 95","Upper 95","model","sigma_a","sigma_x","comp")])

head(cis[[2]])
cis[[2]] <- rbind(cis[[2]],cia[,c("Est","Lower 95","Upper 95","model","sigma_a","sigma_x","comp")])

head(cis[[3]])
cis[[3]] <- rbind(cis[[3]],cib[,c("Est","Lower 95","Upper 95","model","sigma_a","sigma_x","comp")])

sapply(cis,dim)
#cis_x cis_auto cis_both
#[1,] 164570   164570   246855
#[2,]      7        7        7

save(cis,file="varCompCIs.RData")

library(ggplot2)
cis <- getobj("nullSims_8ped/varCompCIs_nullSims.RData")
cib <- cis[["cis_both"]]
head(cib)
# exclude the extreme values

cib <- cib[!is.element(cib$sigma_a,3)&!is.element(cib$sigma_x,3),]
cib$sigma_a <- paste0("sigma_a=",cib$sigma_a)
cib$sigma_x <- paste0("sigma_x=",cib$sigma_x)
cib$comp <- ordered(cib$comp,levels=c("auto","x","resid"))
ftable(cib$comp,cib$sigma_a,cib$sigma_x)
#                    sigma_x=0.3 sigma_x=0.5
# auto  sigma_a=0.3        10010       10010
#       sigma_a=0.5        10010       10010
# x     sigma_a=0.3        10010       10010
#       sigma_a=0.5        10010       10010
# resid sigma_a=0.3        10010       10010
#       sigma_a=0.5        10010       10010

# so 10K iterations of each composition
# add in mean for each component
mns <- expand.grid(comp=c("auto","x","resid"),sigma_a=c("sigma_a=0.3","sigma_a=0.5"),
                         sigma_x=c("sigma_x=0.3","sigma_x=0.5"))
mns$comp <- ordered(mns$comp,levels=c("auto","x","resid"))
mns$y <- -0.1
mns$value <- NA
for(i in seq_len(nrow(mns))){
  mns$value[i] <- mean(cib$Est[cib$comp==mns$comp[i]&cib$sigma_a==mns$sigma_a[i]&cib$sigma_x==mns$sigma_x[i]])
}

mns$value <- format(mns$value,digits=3)


pdf("boxplots_varComp_bothEst_allIters.pdf")
ggplot(cib,aes(x=comp,y=Est)) + geom_boxplot() + facet_grid(sigma_a~sigma_x) + theme_bw() +
  xlab("Variance Component") + ylab("Estimate") +
  ggtitle("Variance Component Estimates of 10K Iterations") +
  geom_text(data=mns, aes(x=comp, y=y, label=value), size=4)
dev.off()

# take only 180 iterations of these
pdf("boxplots_varComp_bothEst.pdf")
ggplot(cib,aes(x=comp,y=Est)) + geom_boxplot() + facet_grid(sigma_x~sigma_a) + theme_bw() +
  xlab("Variance Component") + ylab("Estimate") +
  ggtitle("Variance Component Estimates of 180 Iterations") +
  geom_text(data=mns, aes(x=comp, y=y, label=value), size=4)
dev.off()




##
cib <- cis[["cis_auto"]]
head(cib)
# exclude the extreme values

cib <- cib[!is.element(cib$sigma_a,3)&!is.element(cib$sigma_x,3),]
cib$sigma_a <- paste0("sigma_a=",cib$sigma_a)
cib$sigma_x <- paste0("sigma_x=",cib$sigma_x)
cib$comp <- ordered(cib$comp,levels=c("auto","resid"))
ftable(cib$comp,cib$sigma_a,cib$sigma_x)
#                    sigma_x=0.3 sigma_x=0.5
# auto  sigma_a=0.3        10010       10010
#       sigma_a=0.5        10010       10010
# resid sigma_a=0.3        10010       10010
#       sigma_a=0.5        10010       10010
# so 10K iterations of each composition
# add in mean for each component

mns <- expand.grid(comp=c("auto","resid"),sigma_a=c("sigma_a=0.3","sigma_a=0.5"),
                   sigma_x=c("sigma_x=0.3","sigma_x=0.5"))
mns$comp <- ordered(mns$comp,levels=c("auto","resid"))
mns$y <- 0.33
mns$value <- NA
for(i in seq_len(nrow(mns))){
  mns$value[i] <- mean(cib$Est[cib$comp==mns$comp[i]&cib$sigma_a==mns$sigma_a[i]&cib$sigma_x==mns$sigma_x[i]])
}

mns$value <- format(mns$value,digits=3)


pdf("boxplots_varComp_autoEst_allIters.pdf")
ggplot(cib,aes(x=comp,y=Est)) + geom_boxplot() + facet_grid(sigma_a~sigma_x) + theme_bw() +
  xlab("Variance Component") + ylab("Estimate") +
  ggtitle("Variance Component Estimates of 10K Iterations") +
  geom_text(data=mns, aes(x=comp, y=y, label=value), size=4)
dev.off()

pdf("boxplots_varComp_autoEst.pdf")
ggplot(cib,aes(x=comp,y=Est)) + geom_boxplot() + facet_grid(sigma_x~sigma_a) + theme_bw() +
  xlab("Variance Component") + ylab("Estimate") +
  ggtitle("Variance Component Estimates of 180 Iterations") +
  geom_text(data=mns, aes(x=comp, y=y, label=value), size=4)
dev.off()



##
cib <- cis[["cis_x"]]
head(cib)
# exclude the extreme values

cib <- cib[!is.element(cib$sigma_a,3)&!is.element(cib$sigma_x,3),]
cib$sigma_a <- paste0("sigma_a=",cib$sigma_a)
cib$sigma_x <- paste0("sigma_x=",cib$sigma_x)
cib$comp <- ordered(cib$comp,levels=c("x","resid"))
ftable(cib$comp,cib$sigma_a,cib$sigma_x)
#                    sigma_x=0.3 sigma_x=0.5
# x     sigma_a=0.3        10010       10010
#       sigma_a=0.5        10010       10010
# resid sigma_a=0.3        10010       10010
#       sigma_a=0.5        10010       10010
# so 10K iterations of each composition
# add in mean for each component

mns <- expand.grid(comp=c("x","resid"),sigma_a=c("sigma_a=0.3","sigma_a=0.5"),
                   sigma_x=c("sigma_x=0.3","sigma_x=0.5"))
mns$comp <- ordered(mns$comp,levels=c("x","resid"))
mns$y <- 0.2
mns$value <- NA
for(i in seq_len(nrow(mns))){
  mns$value[i] <- mean(cib$Est[cib$comp==mns$comp[i]&cib$sigma_a==mns$sigma_a[i]&cib$sigma_x==mns$sigma_x[i]])
}

mns$value <- format(mns$value,digits=3)

pdf("boxplots_varComp_xEst_allIters.pdf")
ggplot(cib,aes(x=comp,y=Est)) + geom_boxplot() + facet_grid(sigma_a~sigma_x) + theme_bw() +
  xlab("Variance Component") + ylab("Estimate") +
  ggtitle("Variance Component Estimates of 10K Iterations") +
  geom_text(data=mns, aes(x=comp, y=y, label=value), size=4)
dev.off()

pdf("boxplots_varComp_xEst.pdf")
ggplot(cib,aes(x=comp,y=Est)) + geom_boxplot() + facet_grid(sigma_x~sigma_a) + theme_bw() +
  xlab("Variance Component") + ylab("Estimate") +
  ggtitle("Variance Component Estimates of 180 Iterations") +
  geom_text(data=mns, aes(x=comp, y=y, label=value), size=4)
dev.off()


##
# look at the difference between the truth and mean estimate, for the different scenarios
# when only fitting the AUTO effect:
# sX=0.3, difference=
# sX=0.5, difference=
# so about 1/3 increase

# when only fitting the X effect:
# sA=0.3, difference=
# xA=0.5, difference=
# so about 1/2 decrease

rm(list=ls())


#####
# 36. Add new null sims to var comp estimates

# called in batch:
# cd /projects/geneva/geneva_sata/caitlin/mlm_x/
# qsub -q olga.q -t 1-100 -N nullSims batch_null_simulations.sh
# which calls null_simulations_v2.R (pedigree is different)

library(GWASTools)
setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")

totalRes <- NULL
fns <- list.files(pattern=".RData")
resFn <- substr(fns,1,3)
fns <- fns[resFn=="res"]
for(i in fns){
  dat <- getobj(i)
  tosv <- dat[["mm_res"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  totalRes <- rbind(totalRes,tosv)
}

oldR <- getobj("nullSims_8ped/mmRes_nullSims.RData")
tR <- rbind(oldR,totalRes)

dim(tR) # 3600000 13
ftable(tR$model,tR$sigma_a,tR$sigma_x)
#               0.3    0.5      3
# auto   0.3  150000 150000 150000
# 0.5  150000 150000      0
# 3    150000      0      0
# both   0.3  150000 150000 150000
# 0.5  150000 150000      0
# 3    150000      0      0
# linear 0.3  150000 150000 150000
# 0.5  150000 150000      0
# 3    150000      0      0
# x      0.3  150000 150000 150000
# 0.5  150000 150000      0
# 3    150000      0      0

# with these values we can calculate type I error rate
save(tR,file="nullSims_8ped/mmRes_nullSims_v2.RData")

cis_auto <- NULL
cis_x <- NULL
cis_both <- NULL
for(i in fns){
    dat <- getobj(i)
    tosv <- dat[["ci_x"]]
    tosv$sigma_a <- dat[["sigma_auto"]]
    tosv$sigma_x <- dat[["sigma_x"]]
    tosv$comp <- c("x","resid")

    cis_x <- rbind(cis_x,tosv)

    tosv <- dat[["ci_auto"]]
    tosv$sigma_a <- dat[["sigma_auto"]]
    tosv$sigma_x <- dat[["sigma_x"]]
    tosv$comp <- c("auto","resid")

    cis_auto <- rbind(cis_auto,tosv)

    tosv <- dat[["ci_both"]]
    tosv$sigma_a <- dat[["sigma_auto"]]
    tosv$sigma_x <- dat[["sigma_x"]]
    tosv$comp <- c("x","auto","resid")

    cis_both <- rbind(cis_both,tosv)
  }

oldC <- get(load("nullSims_8ped/varCompCIs_nullSims.RData"))

oldC[["cis_x"]] <- rbind(oldC[["cis_x"]],cis_x)
oldC[["cis_auto"]] <- rbind(oldC[["cis_auto"]],cis_auto)
oldC[["cis_both"]] <- rbind(oldC[["cis_both"]],cis_both)

sum(duplicated(oldC[["cis_x"]]$Est)) # 6
sum(duplicated(oldC[["cis_auto"]]$Est)) # 6
sum(duplicated(oldC[["cis_both"]]$Est)) # 9

cisx <- oldC[["cis_x"]]
dim(cisx) # 80804 7
oldC[["cis_x"]] <- cisx[!duplicated(cisx$Est),]
dim(oldC[["cis_x"]]) # 80798 7

cis <- oldC[["cis_auto"]]
dim(cis) # 80804 7
oldC[["cis_auto"]] <- cis[!duplicated(cis$Est),]
dim(oldC[["cis_auto"]]) # 80798 7

cis <- oldC[["cis_both"]]
dim(cis) # 121206 7
oldC[["cis_both"]] <- cis[!duplicated(cis$Est),]
dim(oldC[["cis_both"]]) # 121197 7

save(oldC,file="nullSims_8ped/varCompCIs_nullSims_v2.RData")

rm(list=ls())


#####
# 37. Process simulations excl g_x term

library(GWASTools)
setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")

totalRes <- NULL
fn <- list.files("nullSims_8ped_exclgX/")
for(i in 1:length(fn)){
  dat <- getobj(paste0("nullSims_8ped_exclgX/",fn[i]))
  tosv <- dat[["mm_res"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  totalRes <- rbind(totalRes,tosv)
}

save(totalRes,file="nullSims_8ped_exclgX/mmRes_nullSims.RData")
dim(totalRes) # 11700000 13
table(totalRes$causal) # all F
ftable(totalRes$model,totalRes$sigma_a,totalRes$sigma_x)
#0.3       3
#auto   0.3        0 1500000
#3    1425000       0
#both   0.3        0 1500000
#3    1425000       0
#linear 0.3        0 1500000
#3    1425000       0
#x      0.3        0 1500000
#3    1425000       0

# so only the extreme values
# check type I error quickly
sum(totalRes$pval[totalRes$model=="auto"]<5e-05)/sum(totalRes$model=="auto")
sum(totalRes$pval[totalRes$model=="both"]<5e-05)/sum(totalRes$model=="both")
sum(totalRes$pval[totalRes$model=="x"]<5e-05)/sum(totalRes$model=="x")
sum(totalRes$pval[totalRes$model=="linear"]<5e-05)/sum(totalRes$model=="linear")

tyIerr <- expand.grid(model=c("auto","both","x","linear"),alpha=c(5e-04,5e-05),
                      sigma_x=c(3,0.3),sigma_a=c(3,0.3))
tyIerr$est <- NA
tyIerr$lower <- NA; tyIerr$upper <- NA
for(i in seq_len(nrow(tyIerr))){
  mmResF <- totalRes[totalRes$model==tyIerr$model[i]&totalRes$sigma_x==tyIerr$sigma_x[i]&
                       totalRes$sigma_a==tyIerr$sigma_a[i],]
  tyIerr$est[i] <- sum(mmResF$pval<tyIerr$alpha[i])/nrow(mmResF)
  stdErr <- sqrt(tyIerr$alpha[i]*(1-tyIerr$alpha[i])/nrow(mmResF))
  tyIerr$lower[i] <- tyIerr$est[i]-1.96*stdErr
  tyIerr$upper[i] <- tyIerr$est[i]+1.96*stdErr
}

tyIerr$sigma_a <- paste0("sigma_a=",tyIerr$sigma_a)
tyIerr$sigma_x <- NULL # doesn't include sigma_x in the phenotype

# make the model an ordered factor so it prints both, x, auto, linear
tyIerr$model <- ordered(tyIerr$model,levels=c("both","x","auto","linear"))
tyIerrSm <- tyIerr[tyIerr$alpha==5e-05&!is.nan(tyIerr$est),]

pdf("typeIErr_8ped_5e05_exclgX.pdf",width=11)
ggplot(tyIerrSm,aes(x=model,y=est)) +geom_point()+ geom_pointrange(aes(ymin=lower,ymax=upper))+
  theme_bw() + geom_hline(aes(yintercept=5e-05,color="gray"))+
  facet_wrap(~sigma_a) +
  ylab("Type I Error Rate, 1.5M Iterations, Pheno Excl g_X Term") + ggtitle(expression(paste("Type I Error Rate, ", alpha, "=5e-05")))
dev.off()

tyIerrSm <- tyIerr[tyIerr$alpha==5e-05&!is.nan(tyIerr$est)&!is.element(tyIerr$model,"linear"),]
pdf("typeIErr_8ped_5e05_exclgX_zoom.pdf",width=11)
ggplot(tyIerrSm,aes(x=model,y=est)) +geom_point()+ geom_pointrange(aes(ymin=lower,ymax=upper))+
  theme_bw() + geom_hline(aes(yintercept=5e-05,color="gray"))+
  facet_wrap(~sigma_a) +
  ylab("Type I Error Rate, 1.5M Iterations, Pheno Excl g_X Term") + ggtitle(expression(paste("Type I Error Rate, ", alpha, "=5e-05")))
dev.off()

##
cis_auto <- NULL
cis_x <- NULL
cis_both <- NULL
for(i in 1:length(fn)){
  dat <- getobj(paste0("nullSims_8ped_exclgX/",fn[i]))
  tosv <- dat[["ci_x"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  tosv$comp <- c("x","resid")

  cis_x <- rbind(cis_x,tosv)

  tosv <- dat[["ci_auto"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  tosv$comp <- c("auto","resid")

  cis_auto <- rbind(cis_auto,tosv)

  tosv <- dat[["ci_both"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  tosv$comp <- c("x","auto","resid")

  cis_both <- rbind(cis_both,tosv)
}

cis <- list("cis_x"=cis_x,"cis_auto"=cis_auto,"cix_both"=cis_both)
save(cis,file="nullSims_8ped_exclgX/varCompCIs_nullSims.RData")

rm(list=ls())


#####
# 38. Null simulations with different 8 person pedigree

# cd /projects/geneva/geneva_sata/caitlin/mlm_x
# qsub -t 1-250 batch_null_simulations.sh
# which calls null_simulations_v5.R, calls with new pedigree file xKC<=autoKC

library(GWASTools)
setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")

fn <- list.files("nullSims_8ped_fem/")
length(fn) # 826

cis_auto <- NULL
cis_x <- NULL
cis_both <- NULL
for(i in 1:length(fn)){
  dat <- getobj(paste0("nullSims_8ped_fem/",fn[i]))
  tosv <- dat[["ci_x"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  tosv$comp <- c("x","resid")

  cis_x <- rbind(cis_x,tosv)

  tosv <- dat[["ci_auto"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  tosv$comp <- c("auto","resid")

  cis_auto <- rbind(cis_auto,tosv)

  tosv <- dat[["ci_both"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  tosv$comp <- c("x","auto","resid")

  cis_both <- rbind(cis_both,tosv)
}

cis <- list("cis_x"=cis_x,"cis_auto"=cis_auto,"cis_both"=cis_both)
save(cis,file="nullSims_8ped_fem/varCompCIs_nullSims.RData")

res1 <- fn[substr(fn,1,4)=="res_"]
res2 <- fn[substr(fn,1,5)=="res2_"]
res3 <- fn[substr(fn,1,5)=="res3_"]
res4 <- fn[substr(fn,1,5)=="res4_"]
totalRes <- NULL
# read in first 25 iters of res, res2, res3, res4
for(i in 1:25){
  dat <- getobj(paste0("nullSims_8ped_fem/",res1[i]))
  tosv <- dat[["mm_res"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  totalRes <- rbind(totalRes,tosv)

  dat <- getobj(paste0("nullSims_8ped_fem/",res2[i]))
  tosv <- dat[["mm_res"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  totalRes <- rbind(totalRes,tosv)

  dat <- getobj(paste0("nullSims_8ped_fem/",res3[i]))
  tosv <- dat[["mm_res"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  totalRes <- rbind(totalRes,tosv)

  dat <- getobj(paste0("nullSims_8ped_fem/",res4[i]))
  tosv <- dat[["mm_res"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  totalRes <- rbind(totalRes,tosv)
}

dim(totalRes) # 6000000 13
ftable(totalRes$model,totalRes$sigma_a,totalRes$sigma_x) # 375K for each config

save(totalRes,file="nullSims_8ped_fem/mmRes_nullSims.RData")

tyIerr <- expand.grid(model=c("auto","both","x","linear"),alpha=c(5e-04,5e-05),
                      sigma_x=c(0.5,0.3),sigma_a=c(0.5,0.3))
tyIerr$est <- NA
tyIerr$lower <- NA; tyIerr$upper <- NA
for(i in seq_len(nrow(tyIerr))){
  mmResF <- totalRes[totalRes$model==tyIerr$model[i]&totalRes$sigma_x==tyIerr$sigma_x[i]&
                       totalRes$sigma_a==tyIerr$sigma_a[i],]
  tyIerr$est[i] <- sum(mmResF$pval<tyIerr$alpha[i])/nrow(mmResF)
  stdErr <- sqrt(tyIerr$alpha[i]*(1-tyIerr$alpha[i])/nrow(mmResF))
  tyIerr$lower[i] <- tyIerr$est[i]-1.96*stdErr
  tyIerr$upper[i] <- tyIerr$est[i]+1.96*stdErr
}

tyIerr$sigma_a <- paste0("sigma_a=",tyIerr$sigma_a)
tyIerr$sigma_x <- paste0("sigma_x=",tyIerr$sigma_x)

# make the model an ordered factor so it prints both, x, auto, linear
tyIerr$model <- ordered(tyIerr$model,levels=c("both","x","auto","linear"))
tyIerrSm <- tyIerr[tyIerr$alpha==5e-05&!is.nan(tyIerr$est),]

pdf("typeIErr_8ped_fem_5e05.pdf")
ggplot(tyIerrSm,aes(x=model,y=est)) +geom_point()+ geom_pointrange(aes(ymin=lower,ymax=upper))+
  theme_bw() + geom_hline(aes(yintercept=5e-05,color="gray"))+
  facet_grid(sigma_x~sigma_a) +
  theme(axis.text=element_text(size=16),axis.title=element_text(size=18),title=element_text(size=20),
       legend.text=element_text(size=16),strip.text = element_text(size=16)) +
  ylab("Type I Error Rate") + ggtitle(expression(paste("Type I Error Rate, ", alpha, "=5e-05")))
dev.off()

tyIerrSm <- tyIerr[tyIerr$alpha==5e-05&!is.nan(tyIerr$est)&!is.element(tyIerr$model,"linear"),]
pdf("typeIErr_8ped_fem_5e05_zoom.pdf")
ggplot(tyIerrSm,aes(x=model,y=est)) +geom_point()+ geom_pointrange(aes(ymin=lower,ymax=upper))+
  theme_bw() + geom_hline(aes(yintercept=5e-05,color="gray"))+
  facet_grid(sigma_x~sigma_a) +
  ylab("Type I Error Rate") + ggtitle(expression(paste("Type I Error Rate, ", alpha, "=5e-05")))
dev.off()

rm(list=ls())


#####
# 39. Var comp estimates for 8 person pedigree (fem)

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")
library(GWASTools); library(ggplot2)

cis <- getobj("nullSims_8ped_fem/varCompCIs_nullSims.RData")

cib <- cis[["cis_both"]]
ftable(cib$comp,cib$sigma_a,cib$sigma_x) # between 183-250 of each configuration

cib$sigma_a=paste0("sigma_a=",cib$sigma_a)
cib$sigma_x=paste0("sigma_x=",cib$sigma_x)
cib$comp=ordered(cib$comp,levels=c("auto","x","resid"))

mns <- expand.grid(comp=c("auto","x","resid"),sigma_a=c("sigma_a=0.3","sigma_a=0.5"),
                   sigma_x=c("sigma_x=0.3","sigma_x=0.5"))
mns$comp <- ordered(mns$comp,levels=c("auto","x","resid"))
mns$y <- -0.1
mns$value <- NA
for(i in seq_len(nrow(mns))){
  mns$value[i] <- mean(cib$Est[cib$comp==mns$comp[i]&cib$sigma_a==mns$sigma_a[i]&cib$sigma_x==mns$sigma_x[i]])
}

mns$value <- format(mns$value,digits=3)

pdf("boxplots_varComp_bothEst_8ped_fem.pdf")
ggplot(cib,aes(x=comp,y=Est)) + geom_boxplot() + facet_grid(sigma_x~sigma_a) + theme_bw() +
   xlab("Variance Component") + ylab("Estimate") +
  ggtitle("Variance Component Estimates of 180 Iterations") +
  geom_text(data=mns, aes(x=comp, y=y, label=value), size=4)
dev.off()

##
cib <- cis[["cis_auto"]]
ftable(cib$comp,cib$sigma_a,cib$sigma_x) # 180-200 of each configuration

cib$sigma_a=paste0("sigma_a=",cib$sigma_a)
cib$sigma_x=paste0("sigma_x=",cib$sigma_x)
cib$comp=ordered(cib$comp,levels=c("auto","resid"))

mns <- expand.grid(comp=c("auto","resid"),sigma_a=c("sigma_a=0.3","sigma_a=0.5"),
                   sigma_x=c("sigma_x=0.3","sigma_x=0.5"))
mns$comp <- ordered(mns$comp,levels=c("auto","resid"))
mns$y <- -0.1
mns$value <- NA
for(i in seq_len(nrow(mns))){
  mns$value[i] <- mean(cib$Est[cib$comp==mns$comp[i]&cib$sigma_a==mns$sigma_a[i]&cib$sigma_x==mns$sigma_x[i]])
}

mns$value <- format(mns$value,digits=3)

pdf("boxplots_varComp_autoEst_8ped_fem.pdf")
ggplot(cib,aes(x=comp,y=Est)) + geom_boxplot() + facet_grid(sigma_x~sigma_a) + theme_bw() +
  xlab("Variance Component") + ylab("Estimate") +
  ggtitle("Variance Component Estimates of 180 Iterations") +
  geom_text(data=mns, aes(x=comp, y=y, label=value), size=4)
dev.off()


##
cib <- cis[["cis_x"]]
ftable(cib$comp,cib$sigma_a,cib$sigma_x) # 180-250 of each configuration

cib$sigma_a=paste0("sigma_a=",cib$sigma_a)
cib$sigma_x=paste0("sigma_x=",cib$sigma_x)
cib$comp=ordered(cib$comp,levels=c("x","resid"))

mns <- expand.grid(comp=c("x","resid"),sigma_a=c("sigma_a=0.3","sigma_a=0.5"),
                   sigma_x=c("sigma_x=0.3","sigma_x=0.5"))
mns$comp <- ordered(mns$comp,levels=c("x","resid"))
mns$y <- -0.1
mns$value <- NA
for(i in seq_len(nrow(mns))){
  mns$value[i] <- mean(cib$Est[cib$comp==mns$comp[i]&cib$sigma_a==mns$sigma_a[i]&cib$sigma_x==mns$sigma_x[i]])
}

mns$value <- format(mns$value,digits=3)

pdf("boxplots_varComp_xEst_8ped_fem.pdf")
ggplot(cib,aes(x=comp,y=Est)) + geom_boxplot() + facet_grid(sigma_x~sigma_a) + theme_bw() +
  xlab("Variance Component") + ylab("Estimate") +
  ggtitle("Variance Component Estimates of 180 Iterations") +
  geom_text(data=mns, aes(x=comp, y=y, label=value), size=4)
dev.off()

rm(list=ls())


#####
# 40. Var comp estimates for 8 ped excl gX term

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")
library(GWASTools); library(ggplot2)

cis <- getobj("nullSims_8ped_exclgX/varCompCIs_nullSims.RData")

cib <- cis[["cix_both"]]
ftable(cib$comp,cib$sigma_a,cib$sigma_x) # only 100 or 95 of each configuration

cib$sigma_a=paste0("sigma_a=",cib$sigma_a)
cib$sigma_x=paste0("sigma_x=",cib$sigma_x)
cib$comp=ordered(cib$comp,levels=c("auto","x","resid"))

mns <- expand.grid(comp=c("auto","x","resid"),sigma_a=c("sigma_a=0.3","sigma_a=3"),
                   sigma_x=c("sigma_x=0.3","sigma_x=3"))
mns$comp <- ordered(mns$comp,levels=c("auto","x","resid"))
mns$y <- -0.1
mns$value <- NA
for(i in seq_len(nrow(mns))){
  mns$value[i] <- mean(cib$Est[cib$comp==mns$comp[i]&cib$sigma_a==mns$sigma_a[i]&cib$sigma_x==mns$sigma_x[i]])
}
mns <- mns[!is.nan(mns$value),]
mns$value <- format(mns$value,digits=3)
mns$sigma_x=NULL # no sigma_x variable

pdf("boxplots_varComp_bothEst_8ped_exclgX.pdf")
ggplot(cib,aes(x=comp,y=Est)) + geom_boxplot() + facet_wrap(~sigma_a) + theme_bw() +
  xlab("Variance Component") + ylab("Estimate") +
  ggtitle("Variance Component Estimates of 100 Iterations") +
  geom_text(data=mns, aes(x=comp, y=y, label=value), size=4)
dev.off()

##
cib <- cis[["cis_auto"]]
ftable(cib$comp,cib$sigma_a,cib$sigma_x) # only 250 of each configuration

cib$sigma_a=paste0("sigma_a=",cib$sigma_a)
cib$sigma_x=paste0("sigma_x=",cib$sigma_x)
cib$comp=ordered(cib$comp,levels=c("auto","resid"))

mns <- expand.grid(comp=c("auto","resid"),sigma_a=c("sigma_a=0.3","sigma_a=0.5"),
                   sigma_x=c("sigma_x=0.3","sigma_x=0.5"))
mns$comp <- ordered(mns$comp,levels=c("auto","resid"))
mns$y <- -0.1
mns$value <- NA
for(i in seq_len(nrow(mns))){
  mns$value[i] <- mean(cib$Est[cib$comp==mns$comp[i]&cib$sigma_a==mns$sigma_a[i]&cib$sigma_x==mns$sigma_x[i]])
}

mns$value <- format(mns$value,digits=3)

pdf("boxplots_varComp_autoEst_8ped_fem.pdf")
ggplot(cib,aes(x=comp,y=Est)) + geom_boxplot() + facet_grid(sigma_a~sigma_x) + theme_bw() +
  xlab("Variance Component") + ylab("Estimate") +
  ggtitle("Variance Component Estimates of 250 Iterations") +
  geom_text(data=mns, aes(x=comp, y=y, label=value), size=4)
dev.off()


##
cib <- cis[["cis_x"]]
ftable(cib$comp,cib$sigma_a,cib$sigma_x) # only 250 of each configuration

cib$sigma_a=paste0("sigma_a=",cib$sigma_a)
cib$sigma_x=paste0("sigma_x=",cib$sigma_x)
cib$comp=ordered(cib$comp,levels=c("x","resid"))

mns <- expand.grid(comp=c("x","resid"),sigma_a=c("sigma_a=0.3","sigma_a=0.5"),
                   sigma_x=c("sigma_x=0.3","sigma_x=0.5"))
mns$comp <- ordered(mns$comp,levels=c("x","resid"))
mns$y <- -0.1
mns$value <- NA
for(i in seq_len(nrow(mns))){
  mns$value[i] <- mean(cib$Est[cib$comp==mns$comp[i]&cib$sigma_a==mns$sigma_a[i]&cib$sigma_x==mns$sigma_x[i]])
}

mns$value <- format(mns$value,digits=3)

pdf("boxplots_varComp_xEst_8ped_fem.pdf")
ggplot(cib,aes(x=comp,y=Est)) + geom_boxplot() + facet_grid(sigma_a~sigma_x) + theme_bw() +
  xlab("Variance Component") + ylab("Estimate") +
  ggtitle("Variance Component Estimates of 250 Iterations") +
  geom_text(data=mns, aes(x=comp, y=y, label=value), size=4)
dev.off()

rm(list=ls())


#####
# 41. Process null simulations for 16 ped

library(GWASTools)
setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")

totalRes <- NULL
for(i in 1:25){
  dat <- getobj(paste0("nullSims_16ped/res_",i,".RData"))
  tosv <- dat[["mm_res"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  totalRes <- rbind(totalRes,tosv)

  dat <- getobj(paste0("nullSims_16ped/res2_",i,".RData"))
  tosv <- dat[["mm_res"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  totalRes <- rbind(totalRes,tosv)

  dat <- getobj(paste0("nullSims_16ped/res3_",i,".RData"))
  tosv <- dat[["mm_res"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  totalRes <- rbind(totalRes,tosv)

  dat <- getobj(paste0("nullSims_16ped/res4_",i,".RData"))
  tosv <- dat[["mm_res"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  totalRes <- rbind(totalRes,tosv)

}

save(totalRes,file="nullSims_16ped/mmRes_nullSims.RData")
dim(totalRes) # 11700000 13
table(totalRes$causal) # all F
ftable(totalRes$model,totalRes$sigma_a,totalRes$sigma_x) # 375K iters for each config

# check type I error quickly
sum(totalRes$pval[totalRes$model=="auto"]<5e-05)/sum(totalRes$model=="auto")
sum(totalRes$pval[totalRes$model=="both"]<5e-05)/sum(totalRes$model=="both")
sum(totalRes$pval[totalRes$model=="x"]<5e-05)/sum(totalRes$model=="x")
sum(totalRes$pval[totalRes$model=="linear"]<5e-05)/sum(totalRes$model=="linear")

tyIerr <- expand.grid(model=c("auto","both","x","linear"),alpha=c(5e-04,5e-05),
                      sigma_x=c(0.5,0.3),sigma_a=c(0.5,0.3))
tyIerr$est <- NA
tyIerr$lower <- NA; tyIerr$upper <- NA
for(i in seq_len(nrow(tyIerr))){
  mmResF <- totalRes[totalRes$model==tyIerr$model[i]&totalRes$sigma_x==tyIerr$sigma_x[i]&
                       totalRes$sigma_a==tyIerr$sigma_a[i],]
  tyIerr$est[i] <- sum(mmResF$pval<tyIerr$alpha[i])/nrow(mmResF)
  stdErr <- sqrt(tyIerr$alpha[i]*(1-tyIerr$alpha[i])/nrow(mmResF))
  tyIerr$lower[i] <- tyIerr$est[i]-1.96*stdErr
  tyIerr$upper[i] <- tyIerr$est[i]+1.96*stdErr
}

tyIerr$sigma_a <- paste0("sigma_a=",tyIerr$sigma_a)
tyIerr$sigma_x <- paste0("sigma_x=",tyIerr$sigma_x)

# make the model an ordered factor so it prints both, x, auto, linear
tyIerr$model <- ordered(tyIerr$model,levels=c("both","x","auto","linear"))
tyIerrSm <- tyIerr[tyIerr$alpha==5e-05&!is.nan(tyIerr$est),]

pdf("typeIErr_16ped_5e05.pdf")
ggplot(tyIerrSm,aes(x=model,y=est)) +geom_point()+ geom_pointrange(aes(ymin=lower,ymax=upper))+
  theme_bw() + geom_hline(aes(yintercept=5e-05,color="gray"))+
  facet_grid(sigma_x~sigma_a) +
  ylab("Type I Error Rate, 375K Iterations") + ggtitle(expression(paste("Type I Error Rate, ", alpha, "=5e-05")))
dev.off()

tyIerrSm <- tyIerr[tyIerr$alpha==5e-05&!is.nan(tyIerr$est)&!is.element(tyIerr$model,"linear"),]
pdf("typeIErr_16ped_5e05_zoom.pdf")
ggplot(tyIerrSm,aes(x=model,y=est)) +geom_point()+ geom_pointrange(aes(ymin=lower,ymax=upper))+
  theme_bw() + geom_hline(aes(yintercept=5e-05,color="gray"))+
  facet_grid(sigma_x~sigma_a) +
  ylab("Type I Error Rate, 375K Iterations") + ggtitle(expression(paste("Type I Error Rate, ", alpha, "=5e-05")))
dev.off()

##
cis_auto <- NULL
cis_x <- NULL
cis_both <- NULL
fn <- list.files("nullSims_16ped")
fn <- fn[-1] # this is the "oldRes" folder

for(i in 1:length(fn)){
  dat <- getobj(paste0("nullSims_16ped/",fn[i]))
  tosv <- dat[["ci_x"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  tosv$comp <- c("x","resid")

  cis_x <- rbind(cis_x,tosv)

  tosv <- dat[["ci_auto"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  tosv$comp <- c("auto","resid")

  cis_auto <- rbind(cis_auto,tosv)

  tosv <- dat[["ci_both"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  tosv$comp <- c("x","auto","resid")

  cis_both <- rbind(cis_both,tosv)
}

cis <- list("cis_x"=cis_x,"cis_auto"=cis_auto,"cis_both"=cis_both)
save(cis,file="nullSims_16ped/varCompCIs_nullSims.RData")

rm(list=ls())


#####
# 42. Var comp estimates for 16 person pedigree

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")
library(GWASTools); library(ggplot2)

cis <- getobj("nullSims_16ped/varCompCIs_nullSims.RData")

cib <- cis[["cis_both"]]
ftable(cib$comp,cib$sigma_a,cib$sigma_x) # 250 of each configuration

cib$sigma_a=paste0("sigma_a=",cib$sigma_a)
cib$sigma_x=paste0("sigma_x=",cib$sigma_x)
cib$comp=ordered(cib$comp,levels=c("auto","x","resid"))

mns <- expand.grid(comp=c("auto","x","resid"),sigma_a=c("sigma_a=0.3","sigma_a=0.5"),
                   sigma_x=c("sigma_x=0.3","sigma_x=0.5"))
mns$comp <- ordered(mns$comp,levels=c("auto","x","resid"))
mns$y <- -0.1
mns$value <- NA
for(i in seq_len(nrow(mns))){
  mns$value[i] <- mean(cib$Est[cib$comp==mns$comp[i]&cib$sigma_a==mns$sigma_a[i]&cib$sigma_x==mns$sigma_x[i]])
}
mns <- mns[!is.nan(mns$value),]
mns$value <- format(mns$value,digits=3)

pdf("boxplots_varComp_bothEst_16ped.pdf")
ggplot(cib,aes(x=comp,y=Est)) + geom_boxplot() + facet_wrap(sigma_x~sigma_a) + theme_bw() +
  xlab("Variance Component") + ylab("Estimate") +
  ggtitle("Variance Component Estimates of 250 Iterations") +
  geom_text(data=mns, aes(x=comp, y=y, label=value), size=4)
dev.off()

##
cib <- cis[["cis_auto"]]
ftable(cib$comp,cib$sigma_a,cib$sigma_x) # only 250 of each configuration

cib$sigma_a=paste0("sigma_a=",cib$sigma_a)
cib$sigma_x=paste0("sigma_x=",cib$sigma_x)
cib$comp=ordered(cib$comp,levels=c("auto","resid"))

mns <- expand.grid(comp=c("auto","resid"),sigma_a=c("sigma_a=0.3","sigma_a=0.5"),
                   sigma_x=c("sigma_x=0.3","sigma_x=0.5"))
mns$comp <- ordered(mns$comp,levels=c("auto","resid"))
mns$y <- -0.1
mns$value <- NA
for(i in seq_len(nrow(mns))){
  mns$value[i] <- mean(cib$Est[cib$comp==mns$comp[i]&cib$sigma_a==mns$sigma_a[i]&cib$sigma_x==mns$sigma_x[i]])
}

mns$value <- format(mns$value,digits=3)

pdf("boxplots_varComp_autoEst_16ped.pdf")
ggplot(cib,aes(x=comp,y=Est)) + geom_boxplot() + facet_grid(sigma_a~sigma_x) + theme_bw() +
  xlab("Variance Component") + ylab("Estimate") +
  ggtitle("Variance Component Estimates of 250 Iterations") +
  geom_text(data=mns, aes(x=comp, y=y, label=value), size=4)
dev.off()


##
cib <- cis[["cis_x"]]
ftable(cib$comp,cib$sigma_a,cib$sigma_x) # only 250 of each configuration

cib$sigma_a=paste0("sigma_a=",cib$sigma_a)
cib$sigma_x=paste0("sigma_x=",cib$sigma_x)
cib$comp=ordered(cib$comp,levels=c("x","resid"))

mns <- expand.grid(comp=c("x","resid"),sigma_a=c("sigma_a=0.3","sigma_a=0.5"),
                   sigma_x=c("sigma_x=0.3","sigma_x=0.5"))
mns$comp <- ordered(mns$comp,levels=c("x","resid"))
mns$y <- -0.1
mns$value <- NA
for(i in seq_len(nrow(mns))){
  mns$value[i] <- mean(cib$Est[cib$comp==mns$comp[i]&cib$sigma_a==mns$sigma_a[i]&cib$sigma_x==mns$sigma_x[i]])
}

mns$value <- format(mns$value,digits=3)

pdf("boxplots_varComp_xEst_16ped.pdf")
ggplot(cib,aes(x=comp,y=Est)) + geom_boxplot() + facet_grid(sigma_a~sigma_x) + theme_bw() +
  xlab("Variance Component") + ylab("Estimate") +
  ggtitle("Variance Component Estimates of 250 Iterations") +
  geom_text(data=mns, aes(x=comp, y=y, label=value), size=4)
dev.off()

rm(list=ls())


#####
# 43. Look at new power results, beta=0.07

#modelC_results_power_sA03sX03
#modelC_results_power_sA05sX03
#modelC_results_power_sA03sX05
#modelC_results_power_sA05sX05/mmRes_SNP1992_beta0.05.RData

# called from batch:
# cd /projects/geneva/geneva_sata/caitlin/mlm_x/olga_application
# qsub -q thornton.q -N readRes batch_read_results.sh

source("powerGraph.R")

dat <- get(load("pvals_newPower_results_apr14.RData"))
newcs <- get(load("pvals_newPower_results_apr16.RData"))
dat <- rbind(dat,newcs)
dim(dat) # 493350 14
table(dat$beta1) # 240K 0.05, 253350 0.07
dat <- dat[dat$beta1==0.07,]

# make sure no duplicates
sum(duplicated(dat$pval))
dat <- dat[!duplicated(dat$pval),]
dim(dat) # 240K 14

toPl <- dat[dat$sigma_x=="03"&dat$sigma_a=="03"&!is.na(dat$pval),]
dim(toPl) # 60000 14
table(toPl$causal) # 30000 true
sa3sx3=powerGraph(toPl,b=0.07,fn="powerGraph_sA03sX03.pdf",xlims=c(0,0.05),values=TRUE)
sa3sx3$sigma_a=paste0("sigma_a=",0.3)
sa3sx3$sigma_x=paste0("sigma_x=",0.3)

toPl <- dat[dat$sigma_x=="03"&dat$sigma_a=="05"&!is.na(dat$pval),]
dim(toPl) # 60K 14
table(toPl$causal) # 30K true
sa5sx3=powerGraph(toPl,b=0.07,fn="powerGraph_sA05sX03.pdf",xlims=c(0,0.05),values=TRUE)
sa5sx3$sigma_a=paste0("sigma_a=",0.5)
sa5sx3$sigma_x=paste0("sigma_x=",0.3)

toPl <- dat[dat$sigma_x=="05"&dat$sigma_a=="03"&!is.na(dat$pval),]
dim(toPl) # 60K 14
table(toPl$causal) # 30K
sa3sx5=powerGraph(toPl,b=0.07,fn="powerGraph_sA03sX05.pdf",xlims=c(0,0.05),values=TRUE)
sa3sx5$sigma_a=paste0("sigma_a=",0.3)
sa3sx5$sigma_x=paste0("sigma_x=",0.5)

toPl <- dat[dat$sigma_x=="05"&dat$sigma_a=="05"&!is.na(dat$pval),]
dim(toPl) # 60K 14
table(toPl$causal) # 30K t
sa5sx5=powerGraph(toPl,b=0.07,fn="powerGraph_sA05sX05.pdf",xlims=c(0,0.05),values=TRUE)
sa5sx5$sigma_a=paste0("sigma_a=",0.5)
sa5sx5$sigma_x=paste0("sigma_x=",0.5)

library(ggplot2); library(reshape)
tot <- rbind(sa3sx3,sa3sx5,sa5sx3,sa5sx5)
tot$x_fp_rate <- tot$x_fp/tot$xiter
tot$x_tp_rate <- tot$x_tp/tot$xiter
tot$auto_fp_rate <- tot$auto_fp/tot$autoiter
tot$auto_tp_rate <- tot$auto_tp/tot$autoiter
tot$both_fp_rate <- tot$both_fp/tot$bothiter
tot$both_tp_rate <- tot$both_tp/tot$bothiter

both <- tot[,c("sigma_a","sigma_x","both_fp_rate","both_tp_rate")]
both$model <- "both"
both$fp <- both$both_fp_rate
both$tp <- both$both_tp_rate

auto <- tot[,c("sigma_a","sigma_x","auto_fp_rate","auto_tp_rate")]
auto$model <- "auto"
auto$fp <- auto$auto_fp_rate
auto$tp <- auto$auto_tp_rate

x <- tot[,c("sigma_a","sigma_x","x_fp_rate","x_tp_rate")]
x$model <- "x"
x$fp <- x$x_fp_rate
x$tp <- x$x_tp_rate

toPl <- rbind(x[,c("sigma_x","sigma_a","model","tp","fp")],
              auto[,c("sigma_x","sigma_a","model","tp","fp")],
              both[,c("sigma_x","sigma_a","model","tp","fp")])

toPl <- toPl[toPl$fp<=0.05,]

pdf("power_beta07.pdf")
ggplot(toPl,aes(x=fp,y=tp,color=model)) + geom_line(size=1.5,aes(linetype=model)) + facet_grid(sigma_a~sigma_x) + theme_bw() +
  ylab("True Positive Rate") + xlab("False Positive Rate") +ggtitle(expression(paste(beta,"=0.07")))+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=18),title=element_text(size=20),
        legend.text=element_text(size=16),strip.text = element_text(size=16))
dev.off()

rm(list=ls())


#####
# 44. New power results, extreme sigma, beta=0.05

#modelC_results_power_sA03sX3
#modelC_results_power_sA3sX03

# called from batch:
# cd /projects/geneva/geneva_sata/caitlin/mlm_x/olga_application
# qsub -q thornton.q -N readRes batch_read_results.sh

source("powerGraph.R")

dat <- get(load("pvals_newPower_results_apr17.RData"))
dim(dat) # 120K 14

toPl <- dat[dat$sigma_x=="03"&dat$sigma_a=="3"&!is.na(dat$pval),]
dim(toPl) # 60000 14
table(toPl$causal) # 30000 true
sa3sx3=powerGraph(toPl,b=0.05,fn="powerGraph_sA03sX03.pdf",xlims=c(0,0.05),values=TRUE)
sa3sx3$sigma_a=paste0("sigma_a=",3)
sa3sx3$sigma_x=paste0("sigma_x=",0.3)

toPl <- dat[dat$sigma_x=="3"&dat$sigma_a=="03"&!is.na(dat$pval),]
dim(toPl) # 60K 14
table(toPl$causal) # 30K true
sa5sx3=powerGraph(toPl,b=0.05,fn="powerGraph_sA05sX03.pdf",xlims=c(0,0.05),values=TRUE)
sa5sx3$sigma_a=paste0("sigma_a=",0.3)
sa5sx3$sigma_x=paste0("sigma_x=",3)

library(ggplot2); library(reshape)
tot <- rbind(sa3sx3,sa5sx3)
tot$x_fp_rate <- tot$x_fp/tot$xiter
tot$x_tp_rate <- tot$x_tp/tot$xiter
tot$auto_fp_rate <- tot$auto_fp/tot$autoiter
tot$auto_tp_rate <- tot$auto_tp/tot$autoiter
tot$both_fp_rate <- tot$both_fp/tot$bothiter
tot$both_tp_rate <- tot$both_tp/tot$bothiter

both <- tot[,c("sigma_a","sigma_x","both_fp_rate","both_tp_rate")]
both$model <- "both"
both$fp <- both$both_fp_rate
both$tp <- both$both_tp_rate

auto <- tot[,c("sigma_a","sigma_x","auto_fp_rate","auto_tp_rate")]
auto$model <- "auto"
auto$fp <- auto$auto_fp_rate
auto$tp <- auto$auto_tp_rate

x <- tot[,c("sigma_a","sigma_x","x_fp_rate","x_tp_rate")]
x$model <- "x"
x$fp <- x$x_fp_rate
x$tp <- x$x_tp_rate

toPl <- rbind(x[,c("sigma_x","sigma_a","model","tp","fp")],
              auto[,c("sigma_x","sigma_a","model","tp","fp")],
              both[,c("sigma_x","sigma_a","model","tp","fp")])

toPl <- toPl[toPl$fp<=0.05,]

pdf("power_beta05_extremeSigma.pdf",width=11)
ggplot(toPl,aes(x=fp,y=tp,color=model)) + geom_line(size=1.5,aes(linetype=model)) + facet_wrap(sigma_a~sigma_x) + theme_bw() +
  ylab("True Positive Rate") + xlab("False Positive Rate") +ggtitle(expression(paste(beta,"=0.05")))+
  theme(axis.text=element_text(size=16),axis.title=element_text(size=18),title=element_text(size=20),
        legend.text=element_text(size=16),strip.text = element_text(size=16))
dev.off()

rm(list=ls())


#####
# 45. Type I error from causal=F SNPs in power results

library(GWASTools); library(ggplot2)
dat1 <- getobj("pvals_newPower_results_apr17.RData")
dat2 <- getobj("pvals_newPower_results_apr14.RData")
dat3 <- getobj("pvals_newPower_results_apr10.RData")

dat <- rbind(dat1,dat2,dat3)
sum(duplicated(dat$pval))

dat$sigma_x[is.element(dat$sigma_x,"03")] <- 0.3
dat$sigma_x[is.element(dat$sigma_x,"05")] <- 0.5
dat$sigma_a[is.element(dat$sigma_a,"03")] <- 0.3
dat$sigma_a[is.element(dat$sigma_a,"05")] <- 0.5

### look at type I error rate from these causal==F results
nullS <- dat[!dat$causal,]
ftable(nullS$model,nullS$sigma_x,nullS$sigma_a) # 20K of each configuration
# 10K for the extreme sigma values

# first just find type I error for non-extreme sigma values
toPl <- nullS[!is.element(nullS$sigma_x,3)&!is.element(nullS$sigma_a,3),]
dim(toPl) # 246675 14

alpha <- 5e-04
tyIerr <- expand.grid(model=c("both","auto","x"),sigma_x=c(0.3,0.5),sigma_a=c(0.3,0.5))
tyIerr$fp <- NA
tyIerr$n <- NA

for(i in seq_len(nrow(tyIerr))){
  tyIerr$fp[i] <- sum(toPl$pval[toPl$model==tyIerr$model[i]&toPl$sigma_x==tyIerr$sigma_x[i]&
                                          toPl$sigma_a==tyIerr$sigma_a[i]]<alpha)
  tyIerr$n[i] <- sum(toPl$model==tyIerr$model[i]&toPl$sigma_x==tyIerr$sigma_x[i]&
                                  toPl$sigma_a==tyIerr$sigma_a[i])

}
tyIerr$ty <- tyIerr$fp/tyIerr$n
tyIerr$lower <- tyIerr$ty-1.96*sqrt((alpha*(1-alpha))/tyIerr$n)
tyIerr$upper <- tyIerr$ty+1.96*sqrt((alpha*(1-alpha))/tyIerr$n)
tyIerr$sigma_a <- paste0("sigma_a=",tyIerr$sigma_a)
tyIerr$sigma_x <- paste0("sigma_x=",tyIerr$sigma_x)

ggplot(tyIerr,aes(x=model,y=ty)) + geom_point() + facet_grid(sigma_a~sigma_x) +
  geom_hline(yintercept=alpha) + geom_segment(aes(x=model,y=lower,xend=model,yend=upper))


tyIerr <- expand.grid(model=c("both","auto","x"),sigma_x=c(0.3,3),sigma_a=c(0.3,3))
tyIerr$fp <- NA
tyIerr$n <- NA

for(i in seq_len(nrow(tyIerr))){
  tyIerr$fp[i] <- sum(nullS$pval[nullS$model==tyIerr$model[i]&nullS$sigma_x==tyIerr$sigma_x[i]&
                                  nullS$sigma_a==tyIerr$sigma_a[i]]<alpha)
  tyIerr$n[i] <- sum(nullS$model==tyIerr$model[i]&nullS$sigma_x==tyIerr$sigma_x[i]&
                       nullS$sigma_a==tyIerr$sigma_a[i])

}
tyIerr$ty <- tyIerr$fp/tyIerr$n
tyIerr$lower <- tyIerr$ty-1.96*sqrt((alpha*(1-alpha))/tyIerr$n)
tyIerr$upper <- tyIerr$ty+1.96*sqrt((alpha*(1-alpha))/tyIerr$n)
tyIerr$sigma_a <- paste0("sigma_a=",tyIerr$sigma_a)
tyIerr$sigma_x <- paste0("sigma_x=",tyIerr$sigma_x)

tyIerr <- tyIerr[4:9,]

ggplot(tyIerr,aes(x=model,y=ty)) + geom_point() + facet_wrap(sigma_a~sigma_x) +
  geom_hline(yintercept=alpha) + geom_segment(aes(x=model,y=lower,xend=model,yend=upper))




toPl <- nullS[!is.element(nullS$sigma_x,3)&!is.element(nullS$sigma_a,3),]
dim(toPl) # 246675 14

alpha <- 1e-03
tyIerr <- expand.grid(model=c("both","auto","x"),sigma_x=c(0.3,0.5),sigma_a=c(0.3,0.5))
tyIerr$fp <- NA
tyIerr$n <- NA

for(i in seq_len(nrow(tyIerr))){
  tyIerr$fp[i] <- sum(toPl$pval[toPl$model==tyIerr$model[i]&toPl$sigma_x==tyIerr$sigma_x[i]&
                                  toPl$sigma_a==tyIerr$sigma_a[i]]<alpha)
  tyIerr$n[i] <- sum(toPl$model==tyIerr$model[i]&toPl$sigma_x==tyIerr$sigma_x[i]&
                       toPl$sigma_a==tyIerr$sigma_a[i])

}
tyIerr$ty <- tyIerr$fp/tyIerr$n
tyIerr$lower <- tyIerr$ty-1.96*sqrt((alpha*(1-alpha))/tyIerr$n)
tyIerr$upper <- tyIerr$ty+1.96*sqrt((alpha*(1-alpha))/tyIerr$n)
tyIerr$sigma_a <- paste0("sigma_a=",tyIerr$sigma_a)
tyIerr$sigma_x <- paste0("sigma_x=",tyIerr$sigma_x)

ggplot(tyIerr,aes(x=model,y=ty)) + geom_point() + facet_grid(sigma_a~sigma_x) +
  geom_hline(yintercept=alpha) + geom_segment(aes(x=model,y=lower,xend=model,yend=upper))


tyIerr <- expand.grid(model=c("both","auto","x"),sigma_x=c(0.3,3),sigma_a=c(0.3,3))
tyIerr$fp <- NA
tyIerr$n <- NA

for(i in seq_len(nrow(tyIerr))){
  tyIerr$fp[i] <- sum(nullS$pval[nullS$model==tyIerr$model[i]&nullS$sigma_x==tyIerr$sigma_x[i]&
                                   nullS$sigma_a==tyIerr$sigma_a[i]]<alpha)
  tyIerr$n[i] <- sum(nullS$model==tyIerr$model[i]&nullS$sigma_x==tyIerr$sigma_x[i]&
                       nullS$sigma_a==tyIerr$sigma_a[i])

}
tyIerr$ty <- tyIerr$fp/tyIerr$n
tyIerr$lower <- tyIerr$ty-1.96*sqrt((alpha*(1-alpha))/tyIerr$n)
tyIerr$upper <- tyIerr$ty+1.96*sqrt((alpha*(1-alpha))/tyIerr$n)
tyIerr$sigma_a <- paste0("sigma_a=",tyIerr$sigma_a)
tyIerr$sigma_x <- paste0("sigma_x=",tyIerr$sigma_x)

tyIerr <- tyIerr[4:9,]

pdf("typeIErr_16ped_1e03_extremeSigma.pdf")
ggplot(tyIerr,aes(x=model,y=ty)) + geom_point() + facet_wrap(sigma_a~sigma_x) + theme_bw() +
  geom_hline(yintercept=alpha) + geom_segment(aes(x=model,y=lower,xend=model,yend=upper)) +
  ylab("Type I Error Rate, 10K Iterations") +
  ggtitle(expression(paste("Type I Error Rate, ",alpha, "=0.001",sep="")))
dev.off()

rm(list=ls())


#####
# # 46. Var comp estiamtes for 16 person ped, extreme sigma

library(GWASTools); library(ggplot2)

totalRes <- getobj("varEst_newPower_results_apr17.RData")

table(totalRes$sigma_x,totalRes$sigma_a) # 70K of extreme sigma values
totalRes$sigma_x[totalRes$sigma_x=="03"] <- 0.3
totalRes$sigma_a[totalRes$sigma_a=="03"] <- 0.3

totalRes$comp <- NA
totalRes$comp[totalRes$model=="both"] <- rep(c("x","auto","resid"),times=sum(totalRes$model=="both")/3)
totalRes$comp[totalRes$model=="x"] <- rep(c("x","resid"),times=sum(totalRes$model=="x")/2)
totalRes$comp[totalRes$model=="auto"] <- rep(c("auto","resid"),times=sum(totalRes$model=="auto")/2)

cib <- totalRes[totalRes$model=="both",]
ftable(cib$comp,cib$sigma_a,cib$sigma_x) # 10K of each configuration

cib$sigma_a=paste0("sigma_a=",cib$sigma_a)
cib$sigma_x=paste0("sigma_x=",cib$sigma_x)
cib$comp=ordered(cib$comp,levels=c("auto","x","resid"))

mns <- expand.grid(comp=c("auto","x","resid"),sigma_a=c("sigma_a=0.3","sigma_a=3"),
                   sigma_x=c("sigma_x=0.3","sigma_x=3"))
mns$comp <- ordered(mns$comp,levels=c("auto","x","resid"))
mns$y <- -0.1
mns$value <- NA
for(i in seq_len(nrow(mns))){
  mns$value[i] <- mean(cib$Est[cib$comp==mns$comp[i]&cib$sigma_a==mns$sigma_a[i]&cib$sigma_x==mns$sigma_x[i]])
}
mns <- mns[!is.nan(mns$value),]
mns$value <- format(mns$value,digits=3)

pdf("boxplots_varComp_bothEst_16ped_extremeSigma.pdf",width=11)
ggplot(cib,aes(x=comp,y=Est)) + geom_boxplot() + facet_wrap(sigma_x~sigma_a) + theme_bw() +
  xlab("Variance Component") + ylab("Estimate") +
  ggtitle("Variance Component Estimates of 10K Iterations") +
  geom_text(data=mns, aes(x=comp, y=y, label=value), size=4)
dev.off()

##
cib <- totalRes[totalRes$model=="auto",]
ftable(cib$comp,cib$sigma_a,cib$sigma_x) # 10K of each configuration

cib$sigma_a=paste0("sigma_a=",cib$sigma_a)
cib$sigma_x=paste0("sigma_x=",cib$sigma_x)
cib$comp=ordered(cib$comp,levels=c("auto","resid"))

mns <- expand.grid(comp=c("auto","resid"),sigma_a=c("sigma_a=0.3","sigma_a=3"),
                   sigma_x=c("sigma_x=0.3","sigma_x=3"))
mns$comp <- ordered(mns$comp,levels=c("auto","resid"))
mns$y <- -0.1
mns$value <- NA
for(i in seq_len(nrow(mns))){
  mns$value[i] <- mean(cib$Est[cib$comp==mns$comp[i]&cib$sigma_a==mns$sigma_a[i]&cib$sigma_x==mns$sigma_x[i]])
}
mns <- mns[!is.nan(mns$value),]
mns$value <- format(mns$value,digits=3)

pdf("boxplots_varComp_autoEst_16ped_extremeSigma.pdf",width=11)
ggplot(cib,aes(x=comp,y=Est)) + geom_boxplot() + facet_wrap(sigma_x~sigma_a) + theme_bw() +
  xlab("Variance Component") + ylab("Estimate") +
  ggtitle("Variance Component Estimates of 10K Iterations") +
  geom_text(data=mns, aes(x=comp, y=y, label=value), size=4)
dev.off()


##
cib <- totalRes[totalRes$model=="x",]
ftable(cib$comp,cib$sigma_a,cib$sigma_x) # 10K of each configuration

cib$sigma_a=paste0("sigma_a=",cib$sigma_a)
cib$sigma_x=paste0("sigma_x=",cib$sigma_x)
cib$comp=ordered(cib$comp,levels=c("x","resid"))

mns <- expand.grid(comp=c("x","resid"),sigma_a=c("sigma_a=0.3","sigma_a=3"),
                   sigma_x=c("sigma_x=0.3","sigma_x=3"))
mns$comp <- ordered(mns$comp,levels=c("x","resid"))
mns$y <- -0.1
mns$value <- NA
for(i in seq_len(nrow(mns))){
  mns$value[i] <- mean(cib$Est[cib$comp==mns$comp[i]&cib$sigma_a==mns$sigma_a[i]&cib$sigma_x==mns$sigma_x[i]])
}
mns <- mns[!is.nan(mns$value),]
mns$value <- format(mns$value,digits=3)

pdf("boxplots_varComp_xEst_16ped_extremeSigma.pdf",width=11)
ggplot(cib,aes(x=comp,y=Est)) + geom_boxplot() + facet_wrap(sigma_x~sigma_a) + theme_bw() +
  xlab("Variance Component") + ylab("Estimate") +
  ggtitle("Variance Component Estimates of 10K Iterations") +
  geom_text(data=mns, aes(x=comp, y=y, label=value), size=4)
dev.off()

rm(list=ls())


#####
# 47. Null simulations with auto SNP testing, 8person ped, 'fem'

# cd /projects/geneva/geneva_sata/caitlin/mlm_x
# qsub -t 1-250 batch_null_simulations.sh
# which calls null_simulations_v6.R - this is 8person ped 'fem', nullSims_8ped_fem_autoSNP

library(GWASTools); library(ggplot2)
setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")

fn <- list.files("nullSims_8ped_fem_autoSNP/")
length(fn) # 830

cis_auto <- NULL
cis_x <- NULL
cis_both <- NULL
for(i in 1:length(fn)){
  dat <- getobj(paste0("nullSims_8ped_fem_autoSNP/",fn[i]))
  tosv <- dat[["ci_x"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  tosv$comp <- c("x","resid")

  cis_x <- rbind(cis_x,tosv)

  tosv <- dat[["ci_auto"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  tosv$comp <- c("auto","resid")

  cis_auto <- rbind(cis_auto,tosv)

  tosv <- dat[["ci_both"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  tosv$comp <- c("x","auto","resid")

  cis_both <- rbind(cis_both,tosv)
}

cis <- list("cis_x"=cis_x,"cis_auto"=cis_auto,"cis_both"=cis_both)
save(cis,file="nullSims_8ped_fem_autoSNP/varCompCIs_nullSims.RData")

res1 <- fn[substr(fn,1,4)=="res_"]
res2 <- fn[substr(fn,1,5)=="res2_"]
res3 <- fn[substr(fn,1,5)=="res3_"]
res4 <- fn[substr(fn,1,5)=="res4_"]
totalRes <- NULL
# read in first 25 iters of res, res2, res3, res4
for(i in 1:25){
  dat <- getobj(paste0("nullSims_8ped_fem_autoSNP/",res1[i]))
  tosv <- dat[["mm_res"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  totalRes <- rbind(totalRes,tosv)

  dat <- getobj(paste0("nullSims_8ped_fem_autoSNP/",res2[i]))
  tosv <- dat[["mm_res"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  totalRes <- rbind(totalRes,tosv)

  dat <- getobj(paste0("nullSims_8ped_fem_autoSNP/",res3[i]))
  tosv <- dat[["mm_res"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  totalRes <- rbind(totalRes,tosv)

  dat <- getobj(paste0("nullSims_8ped_fem_autoSNP/",res4[i]))
  tosv <- dat[["mm_res"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  totalRes <- rbind(totalRes,tosv)
}

dim(totalRes) # 6000000 13
ftable(totalRes$model,totalRes$sigma_a,totalRes$sigma_x) # 375K for each config

save(totalRes,file="nullSims_8ped_fem_autoSNP/mmRes_nullSims.RData")

tyIerr <- expand.grid(model=c("auto","both","x","linear"),alpha=c(5e-04,1e-04,5e-05),
                      sigma_x=c(0.5,0.3),sigma_a=c(0.3,0.5))
tyIerr$est <- NA
tyIerr$lower <- NA; tyIerr$upper <- NA
for(i in seq_len(nrow(tyIerr))){
  mmResF <- totalRes[totalRes$model==tyIerr$model[i]&totalRes$sigma_x==tyIerr$sigma_x[i]&
                       totalRes$sigma_a==tyIerr$sigma_a[i],]
  tyIerr$est[i] <- sum(mmResF$pval<tyIerr$alpha[i])/nrow(mmResF)
  stdErr <- sqrt(tyIerr$alpha[i]*(1-tyIerr$alpha[i])/nrow(mmResF))
  tyIerr$lower[i] <- tyIerr$est[i]-1.96*stdErr
  tyIerr$upper[i] <- tyIerr$est[i]+1.96*stdErr
}

tyIerr$sigma_a <- paste0("sigma_a=",tyIerr$sigma_a)
tyIerr$sigma_x <- paste0("sigma_x=",tyIerr$sigma_x)

# make the model an ordered factor so it prints both, x, auto, linear
tyIerr$model <- ordered(tyIerr$model,levels=c("both","x","auto","linear"))
tyIerrSm <- tyIerr[tyIerr$alpha==5e-04&!is.nan(tyIerr$est),]

pdf("typeIErr_8ped_fem_autoSNP_5e04.pdf")
ggplot(tyIerrSm,aes(x=model,y=est)) +geom_point()+ geom_pointrange(aes(ymin=lower,ymax=upper))+
  theme_bw() + geom_hline(aes(yintercept=5e-04,color="gray"))+
  facet_grid(sigma_x~sigma_a) +
  ylab("Type I Error Rate") + ggtitle(expression(paste("Type I Error Rate, ", alpha, "=5e-04")))
dev.off()

tyIerrSm <- tyIerr[tyIerr$alpha==5e-05&!is.nan(tyIerr$est)&!is.element(tyIerr$model,"linear"),]
pdf("typeIErr_8ped_fem_autoSNP_5e05_zoom.pdf")
ggplot(tyIerrSm,aes(x=model,y=est)) +geom_point()+ geom_pointrange(aes(ymin=lower,ymax=upper))+
  theme_bw() + geom_hline(aes(yintercept=5e-05,color="gray"))+
  facet_grid(sigma_x~sigma_a) +
  ylab("Type I Error Rate") + ggtitle(expression(paste("Type I Error Rate, ", alpha, "=5e-05")))
dev.off()

rm(list=ls())


#####
# 48. Null simulations with auto SNP testing, 8 person ped

# cd /projects/geneva/geneva_sata/caitlin/mlm_x
# qsub -t 1-250 batch_null_simulations.sh
# call null_simulations_v8.R - 8 person ped, nullSims_8ped_autoSNP

library(GWASTools); library(ggplot2)
setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")

fn <- list.files("nullSims_8ped_autoSNP/")
length(fn) # 993

cis_auto <- NULL
cis_x <- NULL
cis_both <- NULL
for(i in 1:length(fn)){
  dat <- getobj(paste0("nullSims_8ped_autoSNP/",fn[i]))
  tosv <- dat[["ci_x"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  tosv$comp <- c("x","resid")

  cis_x <- rbind(cis_x,tosv)

  tosv <- dat[["ci_auto"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  tosv$comp <- c("auto","resid")

  cis_auto <- rbind(cis_auto,tosv)

  tosv <- dat[["ci_both"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  tosv$comp <- c("x","auto","resid")

  cis_both <- rbind(cis_both,tosv)
}

cis <- list("cis_x"=cis_x,"cis_auto"=cis_auto,"cis_both"=cis_both)
save(cis,file="nullSims_8ped_autoSNP/varCompCIs_nullSims.RData")

res1 <- fn[substr(fn,1,4)=="res_"]
res2 <- fn[substr(fn,1,5)=="res2_"]
res3 <- fn[substr(fn,1,5)=="res3_"]
res4 <- fn[substr(fn,1,5)=="res4_"]
totalRes <- NULL
# read in first 25 iters of res, res2, res3, res4
for(i in 1:25){
  dat <- getobj(paste0("nullSims_8ped_autoSNP/",res1[i]))
  tosv <- dat[["mm_res"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  totalRes <- rbind(totalRes,tosv)

  dat <- getobj(paste0("nullSims_8ped_autoSNP/",res2[i]))
  tosv <- dat[["mm_res"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  totalRes <- rbind(totalRes,tosv)

  dat <- getobj(paste0("nullSims_8ped_autoSNP/",res3[i]))
  tosv <- dat[["mm_res"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  totalRes <- rbind(totalRes,tosv)

  dat <- getobj(paste0("nullSims_8ped_autoSNP/",res4[i]))
  tosv <- dat[["mm_res"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  totalRes <- rbind(totalRes,tosv)
}

dim(totalRes) # 6000000 13
ftable(totalRes$model,totalRes$sigma_a,totalRes$sigma_x) # 375K for each config

save(totalRes,file="nullSims_8ped_autoSNP/mmRes_nullSims.RData")

tyIerr <- expand.grid(model=c("auto","both","x","linear"),alpha=c(5e-04,5e-05),
                      sigma_x=c(0.5,0.3),sigma_a=c(0.3,0.5))
tyIerr$est <- NA
tyIerr$lower <- NA; tyIerr$upper <- NA
for(i in seq_len(nrow(tyIerr))){
  mmResF <- totalRes[totalRes$model==tyIerr$model[i]&totalRes$sigma_x==tyIerr$sigma_x[i]&
                       totalRes$sigma_a==tyIerr$sigma_a[i],]
  tyIerr$est[i] <- sum(mmResF$pval<tyIerr$alpha[i])/nrow(mmResF)
  stdErr <- sqrt(tyIerr$alpha[i]*(1-tyIerr$alpha[i])/nrow(mmResF))
  tyIerr$lower[i] <- tyIerr$est[i]-1.96*stdErr
  tyIerr$upper[i] <- tyIerr$est[i]+1.96*stdErr
}

tyIerr$sigma_a <- paste0("sigma_a=",tyIerr$sigma_a)
tyIerr$sigma_x <- paste0("sigma_x=",tyIerr$sigma_x)

# make the model an ordered factor so it prints both, x, auto, linear
tyIerr$model <- ordered(tyIerr$model,levels=c("both","x","auto","linear"))
tyIerrSm <- tyIerr[tyIerr$alpha==5e-05&!is.nan(tyIerr$est),]

pdf("typeIErr_8ped_autoSNP_5e05.pdf")
ggplot(tyIerrSm,aes(x=model,y=est)) +geom_point()+ geom_pointrange(aes(ymin=lower,ymax=upper))+
  theme_bw() + geom_hline(aes(yintercept=5e-05,color="gray"))+
  facet_grid(sigma_x~sigma_a) +
  ylab("Type I Error Rate") + ggtitle(expression(paste("Type I Error Rate, ", alpha, "=5e-05")))
dev.off()

tyIerrSm <- tyIerr[tyIerr$alpha==5e-04&!is.nan(tyIerr$est),]
pdf("typeIErr_8ped_autoSNP_5e04.pdf")
ggplot(tyIerrSm,aes(x=model,y=est)) +geom_point()+ geom_pointrange(aes(ymin=lower,ymax=upper))+
  theme_bw() + geom_hline(aes(yintercept=5e-04,color="gray"))+
  facet_grid(sigma_x~sigma_a) +
  ylab("Type I Error Rate") + ggtitle(expression(paste("Type I Error Rate, ", alpha, "=5e-04")))
dev.off()


tyIerrSm <- tyIerr[tyIerr$alpha==5e-05&!is.nan(tyIerr$est)&!is.element(tyIerr$model,"linear"),]
pdf("typeIErr_8ped_autoSNP_5e05_zoom.pdf")
ggplot(tyIerrSm,aes(x=model,y=est)) +geom_point()+ geom_pointrange(aes(ymin=lower,ymax=upper))+
  theme_bw() + geom_hline(aes(yintercept=5e-05,color="gray"))+
  facet_grid(sigma_x~sigma_a) +
  ylab("Type I Error Rate") + ggtitle(expression(paste("Type I Error Rate, ", alpha, "=5e-05")))
dev.off()

rm(list=ls())


#####
# 49. Null simulations with auto SNP testing, 16 person ped

# cd /projects/geneva/geneva_sata/caitlin/mlm_x
# qsub -t 1-250 batch_null_simulations.sh
# call null_simulations_v7.R - 16 person ped, nullSims_16ped_autoSNP

library(GWASTools); library(ggplot2)
setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")

fn <- list.files("nullSims_16ped_autoSNP/")
length(fn) # 994

cis_auto <- NULL
cis_x <- NULL
cis_both <- NULL
for(i in 1:length(fn)){
  dat <- getobj(paste0("nullSims_16ped_autoSNP/",fn[i]))
  tosv <- dat[["ci_x"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  tosv$comp <- c("x","resid")

  cis_x <- rbind(cis_x,tosv)

  tosv <- dat[["ci_auto"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  tosv$comp <- c("auto","resid")

  cis_auto <- rbind(cis_auto,tosv)

  tosv <- dat[["ci_both"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  tosv$comp <- c("x","auto","resid")

  cis_both <- rbind(cis_both,tosv)
}

cis <- list("cis_x"=cis_x,"cis_auto"=cis_auto,"cis_both"=cis_both)
save(cis,file="nullSims_16ped_autoSNP/varCompCIs_nullSims.RData")


res1 <- fn[substr(fn,1,4)=="res_"]
res2 <- fn[substr(fn,1,5)=="res2_"]
res3 <- fn[substr(fn,1,5)=="res3_"]
res4 <- fn[substr(fn,1,5)=="res4_"]
totalRes <- NULL
# read in first 25 iters of res, res2, res3, res4
for(i in 1:25){
  dat <- getobj(paste0("nullSims_16ped_autoSNP/",res1[i]))
  tosv <- dat[["mm_res"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  totalRes <- rbind(totalRes,tosv)

  dat <- getobj(paste0("nullSims_16ped_autoSNP/",res2[i]))
  tosv <- dat[["mm_res"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  totalRes <- rbind(totalRes,tosv)

  dat <- getobj(paste0("nullSims_16ped_autoSNP/",res3[i]))
  tosv <- dat[["mm_res"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  totalRes <- rbind(totalRes,tosv)

  dat <- getobj(paste0("nullSims_16ped_autoSNP/",res4[i]))
  tosv <- dat[["mm_res"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  totalRes <- rbind(totalRes,tosv)
}

dim(totalRes) # 6000000 13
ftable(totalRes$model,totalRes$sigma_a,totalRes$sigma_x) # 375K for each config

save(totalRes,file="nullSims_16ped_autoSNP/mmRes_nullSims.RData")

tyIerr <- expand.grid(model=c("auto","both","x","linear"),alpha=c(5e-04,5e-05),
                      sigma_x=c(0.5,0.3),sigma_a=c(0.3,0.5))
tyIerr$est <- NA
tyIerr$lower <- NA; tyIerr$upper <- NA
for(i in seq_len(nrow(tyIerr))){
  mmResF <- totalRes[totalRes$model==tyIerr$model[i]&totalRes$sigma_x==tyIerr$sigma_x[i]&
                       totalRes$sigma_a==tyIerr$sigma_a[i],]
  tyIerr$est[i] <- sum(mmResF$pval<tyIerr$alpha[i])/nrow(mmResF)
  stdErr <- sqrt(tyIerr$alpha[i]*(1-tyIerr$alpha[i])/nrow(mmResF))
  tyIerr$lower[i] <- tyIerr$est[i]-1.96*stdErr
  tyIerr$upper[i] <- tyIerr$est[i]+1.96*stdErr
}

tyIerr$sigma_a <- paste0("sigma_a=",tyIerr$sigma_a)
tyIerr$sigma_x <- paste0("sigma_x=",tyIerr$sigma_x)

# make the model an ordered factor so it prints both, x, auto, linear
tyIerr$model <- ordered(tyIerr$model,levels=c("both","x","auto","linear"))
tyIerrSm <- tyIerr[tyIerr$alpha==5e-05&!is.nan(tyIerr$est),]

pdf("typeIErr_16ped_autoSNP_5e05.pdf")
ggplot(tyIerrSm,aes(x=model,y=est)) +geom_point()+ geom_pointrange(aes(ymin=lower,ymax=upper))+
  theme_bw() + geom_hline(aes(yintercept=5e-05,color="gray"))+
  facet_grid(sigma_x~sigma_a) +
  ylab("Type I Error Rate, 375K Iterations") + ggtitle(expression(paste("Type I Error Rate, ", alpha, "=5e-05")))
dev.off()

tyIerrSm <- tyIerr[tyIerr$alpha==5e-05&!is.nan(tyIerr$est)&!is.element(tyIerr$model,"linear"),]
pdf("typeIErr_16ped_autoSNP_5e05_zoom.pdf")
ggplot(tyIerrSm,aes(x=model,y=est)) +geom_point()+ geom_pointrange(aes(ymin=lower,ymax=upper))+
  theme_bw() + geom_hline(aes(yintercept=5e-05,color="gray"))+
  facet_grid(sigma_x~sigma_a) +
  ylab("Type I Error Rate, 375K Iterations") + ggtitle(expression(paste("Type I Error Rate, ", alpha, "=5e-05")))
dev.off()

rm(list=ls())


#####
# 50. Look at beta estimates, 8 person ped

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")
library(GWASTools); library(ggplot2)
library(scales)

resAuto <- getobj("nullSims_8ped_autoSNP/mmRes_nullSims.RData")
resX <- getobj("nullSims_8ped/mmRes_nullSims.RData")

### make a 4x4 plot, with both adj on the x axis and all the others on the y
png("beta_ests_xSNP.png",width=720,height=720)
par(mfrow=c(2,2))
plot(resX$Est[resX$sigma_a==0.5&resX$sigma_x==0.3&resX$model=="both"],xlab="Adj for Both Model",
     resX$Est[resX$sigma_a==0.5&resX$sigma_x==0.3&resX$model=="linear"],ylab="Linear Model",
     main="Beta Estimates for Null SNPs\nTesting X Chr SNP",pch=19,col=alpha("gray",0.5),
     cex.lab=1.8)
abline(0,1)
plot(resX$Est[resX$sigma_a==0.5&resX$sigma_x==0.3&resX$model=="both"],xlab="Adj for Both Model",
     resX$Est[resX$sigma_a==0.5&resX$sigma_x==0.3&resX$model=="auto"],ylab="Auto Model",
     pch=19,col=alpha("gray",0.5),
     cex.lab=2)
abline(0,1)
plot(resX$Est[resX$sigma_a==0.5&resX$sigma_x==0.3&resX$model=="both"],xlab="Adj for Both Model",
     resX$Est[resX$sigma_a==0.5&resX$sigma_x==0.3&resX$model=="x"],ylab="X Model",
     pch=19,col=alpha("gray",0.5),
     cex.lab=1.8)
abline(0,1)
dev.off()

png("beta_ses_xSNP.png",width=720,height=720)
par(mfrow=c(2,2))
plot(resX$SE[resX$sigma_a==0.5&resX$sigma_x==0.3&resX$model=="both"],xlab="Adj for Both Model",
     resX$SE[resX$sigma_a==0.5&resX$sigma_x==0.3&resX$model=="linear"],ylab="Linear Model",
     main="Beta SE Estimates for Null SNPs\nTesting X Chr SNP",pch=19,col=alpha("gray",0.5),
     cex.lab=1.8)
abline(0,1)
plot(resX$SE[resX$sigma_a==0.5&resX$sigma_x==0.3&resX$model=="both"],xlab="Adj for Both Model",
     resX$SE[resX$sigma_a==0.5&resX$sigma_x==0.3&resX$model=="auto"],ylab="Auto Model",
     pch=19,col=alpha("gray",0.5),
     cex.lab=2)
abline(0,1)
plot(resX$SE[resX$sigma_a==0.5&resX$sigma_x==0.3&resX$model=="both"],xlab="Adj for Both Model",
     resX$SE[resX$sigma_a==0.5&resX$sigma_x==0.3&resX$model=="x"],ylab="X Model",
     pch=19,col=alpha("gray",0.5),
     cex.lab=1.8)
abline(0,1)
dev.off()

png("beta_ests_autoSNP.png",width=720,height=720)
par(mfrow=c(2,2))
plot(resAuto$Est[resAuto$sigma_a==0.5&resAuto$sigma_x==0.3&resAuto$model=="both"],xlab="Adj for Both Model",
     resAuto$Est[resAuto$sigma_a==0.5&resAuto$sigma_x==0.3&resAuto$model=="linear"],ylab="Linear Model",
     main="Beta Estimates for Null SNPs\nTesting Autosomal SNP",pch=19,col=alpha("gray",0.5),
     cex.lab=1.5)
abline(0,1)
plot(resAuto$Est[resAuto$sigma_a==0.5&resAuto$sigma_x==0.3&resAuto$model=="both"],xlab="Adj for Both Model",
     resAuto$Est[resAuto$sigma_a==0.5&resAuto$sigma_x==0.3&resAuto$model=="auto"],ylab="Auto Model",
     pch=19,col=alpha("gray",0.5),
     cex.lab=1.5)
abline(0,1)
plot(resAuto$Est[resAuto$sigma_a==0.5&resAuto$sigma_x==0.3&resAuto$model=="both"],xlab="Adj for Both Model",
     resAuto$Est[resAuto$sigma_a==0.5&resAuto$sigma_x==0.3&resAuto$model=="x"],ylab="X Model",
     pch=19,col=alpha("gray",0.5),
     cex.lab=1.5)
abline(0,1)
dev.off()

png("beta_ses_autoSNP.png",width=720,height=720)
par(mfrow=c(2,2))
plot(resAuto$SE[resAuto$sigma_a==0.5&resAuto$sigma_x==0.3&resAuto$model=="both"],xlab="Adj for Both Model",
     resAuto$SE[resAuto$sigma_a==0.5&resAuto$sigma_x==0.3&resAuto$model=="linear"],ylab="Linear Model",
     main="Beta SE Estimates for Null SNPs\nTesting Autosomal SNP",pch=19,col=alpha("gray",0.5),
     cex.lab=1.5)
abline(0,1)
plot(resAuto$SE[resAuto$sigma_a==0.5&resAuto$sigma_x==0.3&resAuto$model=="both"],xlab="Adj for Both Model",
     resAuto$SE[resAuto$sigma_a==0.5&resAuto$sigma_x==0.3&resAuto$model=="auto"],ylab="Auto Model",
     pch=19,col=alpha("gray",0.5),
     cex.lab=1.5)
abline(0,1)
plot(resAuto$SE[resAuto$sigma_a==0.5&resAuto$sigma_x==0.3&resAuto$model=="both"],xlab="Adj for Both Model",
     resAuto$SE[resAuto$sigma_a==0.5&resAuto$sigma_x==0.3&resAuto$model=="x"],ylab="X Model",
     pch=19,col=alpha("gray",0.5),
     cex.lab=1.5)
abline(0,1)
dev.off()

rm(list=ls())


#####
# 51. Var of beta estimates

# generate some genotypes, get the estimates of our cov/var matrices, then look at the products of the
# genotypes with our matrices, see how they vary and compare with eachother

library(BiocParallel)
library(MASS); library(parallel)
library(GWASTools)
library(corpcor); library(doSNOW)
source("sim_phenotype.R")
source("estVarComp.R")
source("allele_drop_functions.R")
source("assocTestMixedModel_v7.R")

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x/")
n <- 8000

SEX <- c("F","M","M","M","F","M","M","M")
sex <- rep(SEX,(n/8))

kinX <- get(load("1000Peds_8ped_xKinship.RData"))
kinAuto <- get(load("1000Peds_8ped_autoKinship.RData"))

kinAuto <- matrix(kinAuto,nrow=n,ncol=n)
kinX <- matrix(kinX,nrow=n,ncol=n)
ident <- diag(nrow(kinAuto))

genoMat <- matrix(NA,nrow=150,ncol=n)
colsToGeno <- split(1:n,rep(1:n,each=8,length=n))
genoAuto <- do.call(cbind,lapply(colsToGeno,function(x){genoMat[,x]=Family_alleles_Nmarker_8ped(rep(0.2,150),150,SEX)}))
dim(genoAuto) # 150 8000; so we have 150 genotypes simulated

genoX <- do.call(cbind,lapply(colsToGeno,function(x){genoMat[,x]=Family_alleles_NmarkerX_8ped(rep(0.2,150),150,SEX)}))
dim(genoX) # 150 8000; so we have 150 genotypes simulated

## read in some var comp estimates
varCI <- getobj("nullSims_8ped/varCompCIs_nullSims.RData")
tmpX <- varCI[["cis_x"]]
tmpX <- tmpX[tmpX$sigma_a==0.5&tmpX$sigma_x==0.3,]
dim(tmpX) # 20020 7

tmpA <- varCI[["cis_auto"]]
tmpA <- tmpA[tmpA$sigma_a==0.5&tmpA$sigma_x==0.3,]
dim(tmpA) # 20020 7

tmpB <- varCI[["cis_both"]]
tmpB <- tmpB[tmpB$sigma_a==0.5&tmpB$sigma_x==0.3,]
dim(tmpB) # 30030 7

# get one var comp estimate
vx <- tmpX$Est[5]*kinX+tmpX$Est[6]*ident
va <- tmpA$Est[5]*kinAuto+tmpA$Est[6]*ident
vt <- tmpB$Est[1]*kinX+tmpB$Est[2]*kinAuto+tmpB$Est[3]*ident

sigt <- 0.3*kinX+0.5*kinAuto+ident

summary(apply(genoX,1,var)) # mean 0.5611
summary(apply(genoAuto,1,var)) # mean 0.3201

vxinv <- solve(vx)
vainv <- solve(va)
vtinv <- solve(vt)

# do for all 150 snps
interc <- rep(1, n)
gavx <- apply(genoAuto,1,function(x){desig <- rbind(interc,x)
                                     solve(desig%*% vxinv %*%t(desig))})
gava <- apply(genoAuto,1,function(x){desig <- rbind(interc,x)
                                     solve(desig%*% vainv %*%t(desig))})
gat <- apply(genoAuto,1,function(x){desig <- rbind(interc,x)
                                    solve(desig%*% vtinv %*%t(desig))})

gxvx <- apply(genoX,1,function(x){desig <- rbind(interc,x)
                                  solve(desig%*% vxinv %*%t(desig))})
gxva <- apply(genoX,1,function(x){desig <- rbind(interc,x)
                                  solve(desig%*% vainv %*%t(desig))})
gxt <- apply(genoX,1,function(x){desig <- rbind(interc,x)
                                 solve(desig%*% vtinv %*%t(desig))})

summary(gavx[4,]) # mean 0.0004858
summary(gava[4,]) # mean 0.0005198
summary(gat[4,]) # mean 0.0005183

summary(gxvx[4,]) # mean 0.0003344
summary(gxva[4,]) # mean 0.0003188
summary(gxt[4,]) # mean 0.0003366

# get middle term too
gavxSig <- apply(genoAuto,1,function(x){desig<-rbind(interc,x)
                                        desig%*% vxinv %*% sigt %*% vxinv %*%t(desig)})
gavaSig <- apply(genoAuto,1,function(x){desig <- rbind(interc,x)
                                        desig%*% vainv %*% sigt %*% vainv %*%t(desig)})
gatSig <- apply(genoAuto,1,function(x){desig <- rbind(interc,x)
                                       desig%*% vtinv %*% sigt %*% vtinv %*%t(desig)})

gxvxSig <- apply(genoX,1,function(x){desig <- rbind(interc,x)
                                     desig%*% vxinv %*% sigt %*% vxinv %*%t(desig)})
gxvaSig <- apply(genoX,1,function(x){desig <- rbind(interc,x)
                                     desig%*% vainv %*% sigt %*% vainv %*%t(desig)})
gxtSig <- apply(genoX,1,function(x){desig <- rbind(interc,x)
                                    desig %*% vtinv %*% sigt %*% vtinv %*%t(desig)})

gxvxMult <- gxvx * gxvxSig
gxvaMult <- gxva * gxvaSig
gavxMult <- gavx * gavxSig
gavaMult <- gava * gavaSig

gxvxProd <- gxvx * gxvxSig * gxvx
summary(gxvxProd[4,] / gxvx[4,])

gxvaProd <- gxva * gxvaSig * gxva
summary(gxvaProd[4,] / gxva[4,])

summary(gxvxMult[4,])
summary(gavxMult[4,])
summary(gxvaMult[4,])
summary(gavaMult[4,])

gavaRat <- gava / gat
gavxRat <- gavx / gat
gxvaRat <- gxva / gxt
gxvxRat <- gxvx / gxt

summary(gavaRat[4,])
summary(gxvaRat[4,])
summary(gxvxRat[4,])
summary(gavxRat[4,])

plot(gava[4,],gat[4,],xlab="gava",ylab="ga true")
abline(0,1)
dev.off()

gxvxF <- gxvxMult * gxvxRat
gxvaF <- gxvaMult * gxvaRat
gavxF <- gavxMult * gavxRat
gavaF <- gavaMult * gavaRat

summary(gavaF[4,])
summary(gxvaF[4,])
summary(gavxF[4,])
summary(gxvxF[4,])

gxvxsq <- sqrt(gxvx * gxvxSig * gxvx)/(sqrt(gxt))
gavxsq <- sqrt(gavx * gavxSig * gavx)/(sqrt(gat))
gavasq <- sqrt(gava * gavaSig * gava)/(sqrt(gat))
gxvasq <- sqrt(gxva * gxvaSig * gxva)/(sqrt(gxt))

summary(gxvxsq[4,]) #>1?
summary(gavxsq[4,]) #<1?
summary(gavasq[4,]) #>1?
summary(gxvasq[4,]) #<1?

### wondering why the type of snp, ie auto vs x chr, changes the ratio of our SEs
# how do the SEs change from gxvx to gavx? - this is same model but x chr vs auto SNP
# same for gava to gxva
plot(gava[4,],gxva[4,],xlab="auto model, auto SNP",ylab="auto model, x snp")
abline(0,1)
dev.off()
# no pattern at all!

plot(gxvx[4,],gavx[4,],xlab="x model, x SNP",ylab="x model, auto snp")
abline(0,1)
dev.off() # no pattern

plot(diag(vtinv),diag(vxinv),xlab="sigma inv",ylab="vx inverse")
abline(0,1)
dev.off() # sigma inverse diags are larger

plot(diag(vtinv),diag(vainv),xlab="sigma inv",ylab="va inverse")
abline(0,1)
dev.off() # sigma inverse diags are larger

plot(as.numeric(vtinv[1:8,1:8]),as.numeric(vainv[1:8,1:8]),xlab="sigma inv",ylab="va inverse",xlim=c(-0.15,0.15),
     ylim=c(-0.15,0.15))
abline(0,1); dev.off() # hmm, not really informative

## look at ratio of true/misspecified model to misspecified model
truegava <- gava * gavaSig * gava
truegxvx <- gxvx * gxvxSig * gxvx
truegavx <- gavx * gavxSig * gavx
truegxva <- gxva * gxvaSig * gxva
truegat <- gat * gatSig * gat
truegxt <- gxt * gxtSig * gxt

estRatgava <- truegava/gava
estRatgavx <- truegavx/gavx
estRatgxva <- truegxva/gxva
estRatgxvx <- truegxvx/gxvx
estRatgxt <- truegxt/gxt
estRatgat <- truegat/gat

estRatgxt <- (gxt * gxvt * sigt * gxvt * gxvt / gxt)

summary(estRatgava[4,]) # mean 1.384
summary(estRatgavx[4,]) # mean 1.443
summary(estRatgxvx[4,]) # mean 1.237
summary(estRatgxva[4,]) # mean 1.275

summary(estRatgxt[4,]) # mean 1.197
summary(estRatgat[4,]) # mean 1.334

rm(list=ls())


#####
# 52. Look at beta estimates, 8 person ped fem

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")
library(GWASTools); library(ggplot2)
library(scales)

resAuto <- getobj("nullSims_8ped_fem_autoSNP/mmRes_nullSims.RData")
resX <- getobj("nullSims_8ped_fem/mmRes_nullSims.RData")

## plot est and ses of same sigma_x,sigma_a values for x and auto models
toPl <- resAuto[resAuto$sigma_a==0.5&resAuto$sigma_x==0.3&is.element(resAuto$model,c("auto","x")),]
png("beta_ests_autoSNP_8ped_fem.png")
plot(toPl$Est[toPl$model=="x"],toPl$Est[toPl$model=="auto"],xlab="X Model",ylab="Auto Model",
     main="Beta Estimates for Null SNPs, Auto vs X Model\nTesting Autosomal SNP",pch=19,col=alpha("gray", 0.5),
     cex.lab=1.5,cex.main=1.5)
legend("topleft",c(paste("cor=",format(cor(toPl$Est[toPl$model=="x"],toPl$Est[toPl$model=="auto"]),digits=3))))
abline(0,1)
dev.off()

## do for x chr snp
toPl <- resX[resX$sigma_a==0.5&resX$sigma_x==0.3&is.element(resX$model,c("auto","x")),]
png("beta_ests_xSNP_8ped_fem.png")
plot(toPl$Est[toPl$model=="x"],toPl$Est[toPl$model=="auto"],xlab="X Model",ylab="Auto Model",
     main="Beta Estimates for Null SNPs, Auto vs X Model\nTesting X Chr SNP",pch=19,col = alpha("gray", 0.5))
legend("topleft",c(paste("cor=",format(cor(toPl$Est[toPl$model=="x"],toPl$Est[toPl$model=="auto"]),digits=3))))
abline(0,1)
dev.off()

toPl <- resX[resX$sigma_a==0.5&resX$sigma_x==0.3&is.element(resX$model,c("x","both")),]
png("beta_ests_xBoth_xSNP_8ped_fem.png")
plot(toPl$Est[toPl$model=="x"],toPl$Est[toPl$model=="both"],xlab="X Model",ylab="Adj For Both Model",
     main="Beta Estimates for Null SNPs, Both vs X Model\nTesting X Chr SNP",pch=19,col = alpha("gray", 0.5))
legend("topleft",c(paste("cor=",format(cor(toPl$Est[toPl$model=="x"],toPl$Est[toPl$model=="both"]),digits=3))))
abline(0,1)
dev.off()

# compare the linear model to the x adj model
toPl <- resX[resX$sigma_a==0.5&resX$sigma_x==0.3&is.element(resX$model,c("x","linear")),]
png("beta_ests_xLin_xSNP_8ped_fem.png")
plot(toPl$Est[toPl$model=="x"],toPl$Est[toPl$model=="linear"],xlab="X Model",ylab="Linear Model",
     main="Beta Estimates for Null SNPs, Linear vs X Model\nTesting X Chr SNP",pch=19,col = alpha("gray", 0.5))
legend("topleft",c(paste("cor=",format(cor(toPl$Est[toPl$model=="x"],toPl$Est[toPl$model=="linear"]),digits=3))))
abline(0,1)
dev.off()

### make a 4x4 plot, with both adj on the x axis and all the others on the y
png("beta_ests_xSNP_8ped_fem.png",width=720,height=720)
par(mfrow=c(2,2))
plot(resX$Est[resX$sigma_a==0.5&resX$sigma_x==0.3&resX$model=="both"],xlab="Adj for Both Model",
     resX$Est[resX$sigma_a==0.5&resX$sigma_x==0.3&resX$model=="linear"],ylab="Linear Model",
     main="Beta Estimates for Null SNPs\nTesting X Chr SNP",pch=19,col=alpha("gray",0.5),
     cex.main=1.5, cex.lab=1.8)
abline(0,1)
plot(resX$Est[resX$sigma_a==0.5&resX$sigma_x==0.3&resX$model=="both"],xlab="Adj for Both Model",
     resX$Est[resX$sigma_a==0.5&resX$sigma_x==0.3&resX$model=="auto"],ylab="Auto Model",
     pch=19,col=alpha("gray",0.5),
     cex.lab=2)
abline(0,1)
plot(resX$Est[resX$sigma_a==0.5&resX$sigma_x==0.3&resX$model=="both"],xlab="Adj for Both Model",
     resX$Est[resX$sigma_a==0.5&resX$sigma_x==0.3&resX$model=="x"],ylab="X Model",
     pch=19,col=alpha("gray",0.5),
     cex.lab=1.8)
abline(0,1)
dev.off()

png("beta_ses_xSNP_8ped_fem.png",width=720,height=720)
par(mfrow=c(2,2))
plot(resX$SE[resX$sigma_a==0.5&resX$sigma_x==0.3&resX$model=="both"],xlab="Adj for Both Model",
     resX$SE[resX$sigma_a==0.5&resX$sigma_x==0.3&resX$model=="linear"],ylab="Linear Model",
     main="Beta SE Estimates for Null SNPs\nTesting X Chr SNP",pch=19,col=alpha("gray",0.5),
     cex.lab=1.8)
abline(0,1)
plot(resX$SE[resX$sigma_a==0.5&resX$sigma_x==0.3&resX$model=="both"],xlab="Adj for Both Model",
     resX$SE[resX$sigma_a==0.5&resX$sigma_x==0.3&resX$model=="auto"],ylab="Auto Model",
     pch=19,col=alpha("gray",0.5),
     cex.lab=2)
abline(0,1)
plot(resX$SE[resX$sigma_a==0.5&resX$sigma_x==0.3&resX$model=="both"],xlab="Adj for Both Model",
     resX$SE[resX$sigma_a==0.5&resX$sigma_x==0.3&resX$model=="x"],ylab="X Model",
     pch=19,col=alpha("gray",0.5),
     cex.lab=1.8)
abline(0,1)
dev.off()

png("beta_ests_autoSNP_8ped_fem.png",width=720,height=720)
par(mfrow=c(2,2))
plot(resAuto$Est[resAuto$sigma_a==0.5&resAuto$sigma_x==0.3&resAuto$model=="both"],xlab="Adj for Both Model",
     resAuto$Est[resAuto$sigma_a==0.5&resAuto$sigma_x==0.3&resAuto$model=="linear"],ylab="Linear Model",
     main="Beta Estimates for Null SNPs\nTesting Autosomal SNP",pch=19,col=alpha("gray",0.5),
     cex.lab=1.8)
abline(0,1)
plot(resAuto$Est[resAuto$sigma_a==0.5&resAuto$sigma_x==0.3&resAuto$model=="both"],xlab="Adj for Both Model",
     resAuto$Est[resAuto$sigma_a==0.5&resAuto$sigma_x==0.3&resAuto$model=="auto"],ylab="Auto Model",
     pch=19,col=alpha("gray",0.5),
     cex.lab=2)
abline(0,1)
plot(resAuto$Est[resAuto$sigma_a==0.5&resAuto$sigma_x==0.3&resAuto$model=="both"],xlab="Adj for Both Model",
     resAuto$Est[resAuto$sigma_a==0.5&resAuto$sigma_x==0.3&resAuto$model=="x"],ylab="X Model",
     pch=19,col=alpha("gray",0.5),
     cex.lab=1.8)
abline(0,1)
dev.off()

png("beta_ses_autoSNP_8ped_fem.png",width=720,height=720)
par(mfrow=c(2,2))
plot(resAuto$SE[resAuto$sigma_a==0.5&resAuto$sigma_x==0.3&resAuto$model=="both"],xlab="Adj for Both Model",
     resAuto$SE[resAuto$sigma_a==0.5&resAuto$sigma_x==0.3&resAuto$model=="linear"],ylab="Linear Model",
     main="Beta SE Estimates for Null SNPs\nTesting Autosomal SNP",pch=19,col=alpha("gray",0.5),
     cex.lab=1.5)
abline(0,1)
plot(resAuto$SE[resAuto$sigma_a==0.5&resAuto$sigma_x==0.3&resAuto$model=="both"],xlab="Adj for Both Model",
     resAuto$SE[resAuto$sigma_a==0.5&resAuto$sigma_x==0.3&resAuto$model=="auto"],ylab="Auto Model",
     pch=19,col=alpha("gray",0.5),
     cex.lab=1.5)
abline(0,1)
plot(resAuto$SE[resAuto$sigma_a==0.5&resAuto$sigma_x==0.3&resAuto$model=="both"],xlab="Adj for Both Model",
     resAuto$SE[resAuto$sigma_a==0.5&resAuto$sigma_x==0.3&resAuto$model=="x"],ylab="X Model",
     pch=19,col=alpha("gray",0.5),
     cex.lab=1.5)
abline(0,1)
dev.off()

rm(list=ls())


#####
# 53. Null simulations with auto SNP, 8 person ped fem, extreme sigma

library(GWASTools)
setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")

fn <- list.files("nullSims_8ped_fem/")
length(fn) # 826

fn <- fn[grep("extremeSig",fn)]
length(fn) # 500

cis_auto <- NULL
cis_x <- NULL
cis_both <- NULL
for(i in 1:length(fn)){
  dat <- getobj(paste0("nullSims_8ped_fem/",fn[i]))
  tosv <- dat[["ci_x"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  tosv$comp <- c("x","resid")

  cis_x <- rbind(cis_x,tosv)

  tosv <- dat[["ci_auto"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  tosv$comp <- c("auto","resid")

  cis_auto <- rbind(cis_auto,tosv)

  tosv <- dat[["ci_both"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  tosv$comp <- c("x","auto","resid")

  cis_both <- rbind(cis_both,tosv)
}

cis <- list("cis_x"=cis_x,"cis_auto"=cis_auto,"cis_both"=cis_both)
save(cis,file="nullSims_8ped_fem/varCompCIs_extremeSig_autoSNP_nullSims.RData")

res1 <- fn[substr(fn,1,4)=="res_"]
res2 <- fn[substr(fn,1,5)=="res2_"]

totalRes <- NULL
# read in first 25 iters of res, res2, res3, res4
for(i in 1:50){
  dat <- getobj(paste0("nullSims_8ped_fem/",res1[i]))
  tosv <- dat[["mm_res"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  totalRes <- rbind(totalRes,tosv)

  dat <- getobj(paste0("nullSims_8ped_fem/",res2[i]))
  tosv <- dat[["mm_res"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  totalRes <- rbind(totalRes,tosv)
}

dim(totalRes) # 6000000 13
ftable(totalRes$model,totalRes$sigma_a,totalRes$sigma_x) # 375K for each config

save(totalRes,file="nullSims_8ped_fem/mmRes_extremeSig_autoSNP_nullSims.RData")

## make a type I error plot for the extreme sigma values

## extreme sigma, now
tyIerr <- expand.grid(model=c("auto","both","x","linear"),alpha=c(5e-03,5e-04,5e-05),
                      sigma_x=c(0.3,3),sigma_a=c(0.3,3))
tyIerr$est <- NA
tyIerr$lower <- NA; tyIerr$upper <- NA
for(i in seq_len(nrow(tyIerr))){
  mmResF <- totalRes[totalRes$model==tyIerr$model[i]&totalRes$sigma_x==tyIerr$sigma_x[i]&
                       totalRes$sigma_a==tyIerr$sigma_a[i],]
  tyIerr$est[i] <- sum(mmResF$pval<tyIerr$alpha[i])/nrow(mmResF)
  stdErr <- sqrt(tyIerr$alpha[i]*(1-tyIerr$alpha[i])/nrow(mmResF))
  tyIerr$lower[i] <- tyIerr$est[i]-1.96*stdErr
  tyIerr$upper[i] <- tyIerr$est[i]+1.96*stdErr
}

tyIerr$sigma_a <- paste0("sigma_a=",tyIerr$sigma_a)
tyIerr$sigma_x <- paste0("sigma_x=",tyIerr$sigma_x)

tyIerr <- tyIerr[!is.nan(tyIerr$est),]
tyIerr <- tyIerr[!(tyIerr$sigma_x=="sigma_x=0.3"&tyIerr$sigma_a=="sigma_a=0.3"),]

# make the model an ordered factor so it prints both, x, auto, linear
tyIerr$model <- ordered(tyIerr$model,levels=c("both","x","auto","linear"))
tyIerrSm <- tyIerr[tyIerr$alpha==5e-05,]

pdf("typeIErr_8ped_fem_5e05_autoSNP_extremeSigma.pdf")
ggplot(tyIerrSm,aes(x=model,y=est)) +geom_point()+ geom_pointrange(aes(ymin=lower,ymax=upper))+
  theme_bw() + geom_hline(aes(yintercept=5e-05,color="gray"))+
  facet_wrap(sigma_x~sigma_a) +
  ylab("Type I Error Rate") + ggtitle(expression(paste("Type I Error Rate, ", alpha, "=5e-05")))
dev.off()

rm(list=ls())


#####
# 54. Null simulations with x chr SNP, 8 person ped fem, extreme sigma

library(GWASTools); library(ggplot2)
setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")

fn <- list.files("nullSims_8ped_fem/")
length(fn) # 826

fn <- fn[grep("extremeSig",fn)]
length(fn) # 500

cis_auto <- NULL
cis_x <- NULL
cis_both <- NULL
for(i in 1:length(fn)){
  dat <- getobj(paste0("nullSims_8ped_fem/",fn[i]))
  tosv <- dat[["ci_x"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  tosv$comp <- c("x","resid")

  cis_x <- rbind(cis_x,tosv)

  tosv <- dat[["ci_auto"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  tosv$comp <- c("auto","resid")

  cis_auto <- rbind(cis_auto,tosv)

  tosv <- dat[["ci_both"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  tosv$comp <- c("x","auto","resid")

  cis_both <- rbind(cis_both,tosv)
}

cis <- list("cis_x"=cis_x,"cis_auto"=cis_auto,"cis_both"=cis_both)
save(cis,file="nullSims_8ped_fem/varCompCIs_extremeSig_nullSims.RData")

res1 <- fn[substr(fn,1,4)=="res_"]
res2 <- fn[substr(fn,1,5)=="res2_"]

totalRes <- NULL
# read in first 25 iters of res, res2, res3, res4
for(i in 1:50){
  dat <- getobj(paste0("nullSims_8ped_fem/",res1[i]))
  tosv <- dat[["mm_res"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  totalRes <- rbind(totalRes,tosv)

  dat <- getobj(paste0("nullSims_8ped_fem/",res2[i]))
  tosv <- dat[["mm_res"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  totalRes <- rbind(totalRes,tosv)
}

dim(totalRes) # 6000000 13
ftable(totalRes$model,totalRes$sigma_a,totalRes$sigma_x) # 375K for each config

save(totalRes,file="nullSims_8ped_fem/mmRes_extremeSig_nullSims.RData")

## make a type I error plot for the extreme sigma values

## extreme sigma, now
tyIerr <- expand.grid(model=c("auto","both","x","linear"),alpha=c(5e-03,5e-04,5e-05),
                      sigma_x=c(0.3,3),sigma_a=c(0.3,3))
tyIerr$est <- NA
tyIerr$lower <- NA; tyIerr$upper <- NA
for(i in seq_len(nrow(tyIerr))){
  mmResF <- totalRes[totalRes$model==tyIerr$model[i]&totalRes$sigma_x==tyIerr$sigma_x[i]&
                       totalRes$sigma_a==tyIerr$sigma_a[i],]
  tyIerr$est[i] <- sum(mmResF$pval<tyIerr$alpha[i])/nrow(mmResF)
  stdErr <- sqrt(tyIerr$alpha[i]*(1-tyIerr$alpha[i])/nrow(mmResF))
  tyIerr$lower[i] <- tyIerr$est[i]-1.96*stdErr
  tyIerr$upper[i] <- tyIerr$est[i]+1.96*stdErr
}

tyIerr$sigma_a <- paste0("sigma_a=",tyIerr$sigma_a)
tyIerr$sigma_x <- paste0("sigma_x=",tyIerr$sigma_x)

tyIerr <- tyIerr[!is.nan(tyIerr$est),]
tyIerr <- tyIerr[!(tyIerr$sigma_x=="sigma_x=0.3"&tyIerr$sigma_a=="sigma_a=0.3"),]

# make the model an ordered factor so it prints both, x, auto, linear
tyIerr$model <- ordered(tyIerr$model,levels=c("both","x","auto","linear"))
tyIerrSm <- tyIerr[tyIerr$alpha==5e-05,]

pdf("typeIErr_8ped_fem_5e05_extremeSigma.pdf")
ggplot(tyIerrSm,aes(x=model,y=est)) +geom_point()+ geom_pointrange(aes(ymin=lower,ymax=upper))+
  theme_bw() + geom_hline(aes(yintercept=5e-05,color="gray"))+
  facet_wrap(sigma_x~sigma_a) +
  theme(axis.text=element_text(size=16),axis.title=element_text(size=18),title=element_text(size=20),
        legend.text=element_text(size=16),strip.text = element_text(size=16)) +
  ylab("Type I Error Rate") + ggtitle(expression(paste("Type I Error Rate, ", alpha, "=5e-05")))
dev.off()

rm(list=ls())


#####
# 55. Look at beta estimates, 8 person ped - different sigma values

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")
library(GWASTools); library(ggplot2)
library(scales)

resAuto <- getobj("nullSims_8ped_autoSNP/mmRes_nullSims.RData")
resX <- getobj("nullSims_8ped/mmRes_nullSims.RData")

### make a 4x4 plot, with both adj on the x axis and all the others on the y
png("beta_ests_xSNP_s2A03s2X05.png",width=720,height=720)
par(mfrow=c(2,2))
plot(resX$Est[resX$sigma_x==0.5&resX$sigma_a==0.3&resX$model=="both"],xlab="Adj for Both Model",
     resX$Est[resX$sigma_x==0.5&resX$sigma_a==0.3&resX$model=="linear"],ylab="Linear Model",
     main="Beta Estimates for Null SNPs\nTesting X Chr SNP",pch=19,col=alpha("gray",0.5),
     cex.lab=1.5)
abline(0,1)
plot(resX$Est[resX$sigma_x==0.5&resX$sigma_a==0.3&resX$model=="both"],xlab="Adj for Both Model",
     resX$Est[resX$sigma_x==0.5&resX$sigma_a==0.3&resX$model=="auto"],ylab="Auto Model",
     pch=19,col=alpha("gray",0.5),
     cex.lab=1.5)
abline(0,1)
plot(resX$Est[resX$sigma_x==0.5&resX$sigma_a==0.3&resX$model=="both"],xlab="Adj for Both Model",
     resX$Est[resX$sigma_x==0.5&resX$sigma_a==0.3&resX$model=="x"],ylab="X Model",
     pch=19,col=alpha("gray",0.5),
     cex.lab=1.5)
abline(0,1)
dev.off()

png("beta_ses_xSNP_s2A03s2X05.png",width=720,height=720)
par(mfrow=c(2,2))
plot(resX$SE[resX$sigma_x==0.5&resX$sigma_a==0.3&resX$model=="both"],xlab="Adj for Both Model",
     resX$SE[resX$sigma_x==0.5&resX$sigma_a==0.3&resX$model=="linear"],ylab="Linear Model",
     main="Beta SE Estimates for Null SNPs\nTesting X Chr SNP",pch=19,col=alpha("gray",0.5),
     cex.lab=1.5)
abline(0,1)
plot(resX$SE[resX$sigma_x==0.5&resX$sigma_a==0.3&resX$model=="both"],xlab="Adj for Both Model",
     resX$SE[resX$sigma_x==0.5&resX$sigma_a==0.3&resX$model=="auto"],ylab="Auto Model",
     pch=19,col=alpha("gray",0.5),
     cex.lab=1.5)
abline(0,1)
plot(resX$SE[resX$sigma_x==0.5&resX$sigma_a==0.3&resX$model=="both"],xlab="Adj for Both Model",
     resX$SE[resX$sigma_x==0.5&resX$sigma_a==0.3&resX$model=="x"],ylab="X Model",
     pch=19,col=alpha("gray",0.5),
     cex.lab=1.5)
abline(0,1)
dev.off()

png("beta_ests_autoSNP_s2A03s2X05.png",width=720,height=720)
par(mfrow=c(2,2))
plot(resAuto$Est[resAuto$sigma_x==0.5&resAuto$sigma_a==0.3&resAuto$model=="both"],xlab="Adj for Both Model",
     resAuto$Est[resAuto$sigma_x==0.5&resAuto$sigma_a==0.3&resAuto$model=="linear"],ylab="Linear Model",
     main="Beta Estimates for Null SNPs\nTesting Autosomal SNP",pch=19,col=alpha("gray",0.5),
     cex.lab=1.5)
abline(0,1)
plot(resAuto$Est[resAuto$sigma_x==0.5&resAuto$sigma_a==0.3&resAuto$model=="both"],xlab="Adj for Both Model",
     resAuto$Est[resAuto$sigma_x==0.5&resAuto$sigma_a==0.3&resAuto$model=="auto"],ylab="Auto Model",
     pch=19,col=alpha("gray",0.5),
     cex.lab=1.5)
abline(0,1)
plot(resAuto$Est[resAuto$sigma_x==0.5&resAuto$sigma_a==0.3&resAuto$model=="both"],xlab="Adj for Both Model",
     resAuto$Est[resAuto$sigma_x==0.5&resAuto$sigma_a==0.3&resAuto$model=="x"],ylab="X Model",
     pch=19,col=alpha("gray",0.5),
     cex.lab=1.5)
abline(0,1)
dev.off()

png("beta_ses_autoSNP_s2A03s2X05.png",width=720,height=720)
par(mfrow=c(2,2))
plot(resAuto$SE[resAuto$sigma_x==0.5&resAuto$sigma_a==0.3&resAuto$model=="both"],xlab="Adj for Both Model",
     resAuto$SE[resAuto$sigma_x==0.5&resAuto$sigma_a==0.3&resAuto$model=="linear"],ylab="Linear Model",
     main="Beta SE Estimates for Null SNPs\nTesting Autosomal SNP",pch=19,col=alpha("gray",0.5),
     cex.lab=1.5)
abline(0,1)
plot(resAuto$SE[resAuto$sigma_x==0.5&resAuto$sigma_a==0.3&resAuto$model=="both"],xlab="Adj for Both Model",
     resAuto$SE[resAuto$sigma_x==0.5&resAuto$sigma_a==0.3&resAuto$model=="auto"],ylab="Auto Model",
     pch=19,col=alpha("gray",0.5),
     cex.lab=1.5)
abline(0,1)
plot(resAuto$SE[resAuto$sigma_x==0.5&resAuto$sigma_a==0.3&resAuto$model=="both"],xlab="Adj for Both Model",
     resAuto$SE[resAuto$sigma_x==0.5&resAuto$sigma_a==0.3&resAuto$model=="x"],ylab="X Model",
     pch=19,col=alpha("gray",0.5),
     cex.lab=1.5)
abline(0,1)
dev.off()

rm(list=ls())


#####
# 56. Look at beta estimates, 8 person ped fem - different sigma values

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")
library(GWASTools); library(ggplot2)
library(scales)

resAuto <- getobj("nullSims_8ped_fem_autoSNP/mmRes_nullSims.RData")
resX <- getobj("nullSims_8ped_fem/mmRes_nullSims.RData")

### make a 4x4 plot, with both adj on the x axis and all the others on the y
png("beta_ests_xSNP_8ped_fem_s2A03s2X05.png",width=720,height=720)
par(mfrow=c(2,2))
plot(resX$Est[resX$sigma_x==0.5&resX$sigma_a==0.3&resX$model=="both"],xlab="Adj for Both Model",
     resX$Est[resX$sigma_x==0.5&resX$sigma_a==0.3&resX$model=="linear"],ylab="Linear Model",
     main="Beta Estimates for Null SNPs\nTesting X Chr SNP",pch=19,col=alpha("gray",0.5),
     cex.lab=1.5)
abline(0,1)
plot(resX$Est[resX$sigma_x==0.5&resX$sigma_a==0.3&resX$model=="both"],xlab="Adj for Both Model",
     resX$Est[resX$sigma_x==0.5&resX$sigma_a==0.3&resX$model=="auto"],ylab="Auto Model",
     pch=19,col=alpha("gray",0.5),
     cex.lab=1.5)
abline(0,1)
plot(resX$Est[resX$sigma_x==0.5&resX$sigma_a==0.3&resX$model=="both"],xlab="Adj for Both Model",
     resX$Est[resX$sigma_x==0.5&resX$sigma_a==0.3&resX$model=="x"],ylab="X Model",
     pch=19,col=alpha("gray",0.5),
     cex.lab=1.5)
abline(0,1)
dev.off()

png("beta_ses_xSNP_8ped_fem_s2A03s2X05.png",width=720,height=720)
par(mfrow=c(2,2))
plot(resX$SE[resX$sigma_x==0.5&resX$sigma_a==0.3&resX$model=="both"],xlab="Adj for Both Model",
     resX$SE[resX$sigma_x==0.5&resX$sigma_a==0.3&resX$model=="linear"],ylab="Linear Model",
     main="Beta SE Estimates for Null SNPs\nTesting X Chr SNP",pch=19,col=alpha("gray",0.5),
     cex.lab=1.5)
abline(0,1)
plot(resX$SE[resX$sigma_x==0.5&resX$sigma_a==0.3&resX$model=="both"],xlab="Adj for Both Model",
     resX$SE[resX$sigma_x==0.5&resX$sigma_a==0.3&resX$model=="auto"],ylab="Auto Model",
     pch=19,col=alpha("gray",0.5),
     cex.lab=1.5)
abline(0,1)
plot(resX$SE[resX$sigma_x==0.5&resX$sigma_a==0.3&resX$model=="both"],xlab="Adj for Both Model",
     resX$SE[resX$sigma_x==0.5&resX$sigma_a==0.3&resX$model=="x"],ylab="X Model",
     pch=19,col=alpha("gray",0.5),
     cex.lab=1.5)
abline(0,1)
dev.off()

png("beta_ests_autoSNP_8ped_fem_s2A03s2X05.png",width=720,height=720)
par(mfrow=c(2,2))
plot(resAuto$Est[resAuto$sigma_x==0.5&resAuto$sigma_a==0.3&resAuto$model=="both"],xlab="Adj for Both Model",
     resAuto$Est[resAuto$sigma_x==0.5&resAuto$sigma_a==0.3&resAuto$model=="linear"],ylab="Linear Model",
     main="Beta Estimates for Null SNPs\nTesting Autosomal SNP",pch=19,col=alpha("gray",0.5),
     cex.lab=1.5)
abline(0,1)
plot(resAuto$Est[resAuto$sigma_x==0.5&resAuto$sigma_a==0.3&resAuto$model=="both"],xlab="Adj for Both Model",
     resAuto$Est[resAuto$sigma_x==0.5&resAuto$sigma_a==0.3&resAuto$model=="auto"],ylab="Auto Model",
     pch=19,col=alpha("gray",0.5),
     cex.lab=1.5)
abline(0,1)
plot(resAuto$Est[resAuto$sigma_x==0.5&resAuto$sigma_a==0.3&resAuto$model=="both"],xlab="Adj for Both Model",
     resAuto$Est[resAuto$sigma_x==0.5&resAuto$sigma_a==0.3&resAuto$model=="x"],ylab="X Model",
     pch=19,col=alpha("gray",0.5),
     cex.lab=1.5)
abline(0,1)
dev.off()

png("beta_ses_autoSNP_8ped_fem_s2A03s2X05.png",width=720,height=720)
par(mfrow=c(2,2))
plot(resAuto$SE[resAuto$sigma_x==0.5&resAuto$sigma_a==0.3&resAuto$model=="both"],xlab="Adj for Both Model",
     resAuto$SE[resAuto$sigma_x==0.5&resAuto$sigma_a==0.3&resAuto$model=="linear"],ylab="Linear Model",
     main="Beta SE Estimates for Null SNPs\nTesting Autosomal SNP",pch=19,col=alpha("gray",0.5),
     cex.lab=1.5)
abline(0,1)
plot(resAuto$SE[resAuto$sigma_x==0.5&resAuto$sigma_a==0.3&resAuto$model=="both"],xlab="Adj for Both Model",
     resAuto$SE[resAuto$sigma_x==0.5&resAuto$sigma_a==0.3&resAuto$model=="auto"],ylab="Auto Model",
     pch=19,col=alpha("gray",0.5),
     cex.lab=1.5)
abline(0,1)
plot(resAuto$SE[resAuto$sigma_x==0.5&resAuto$sigma_a==0.3&resAuto$model=="both"],xlab="Adj for Both Model",
     resAuto$SE[resAuto$sigma_x==0.5&resAuto$sigma_a==0.3&resAuto$model=="x"],ylab="X Model",
     pch=19,col=alpha("gray",0.5),
     cex.lab=1.5)
abline(0,1)
dev.off()

rm(list=ls())


#####
# 57. Power estimates, 8 ped: process results

# called power_simulations_v2.R
# stored results in powerSims_8ped/

library(GWASTools); library(ggplot2)
setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")

fn <- list.files("powerSims_8ped/")
length(fn) # 21000

cis_auto <- NULL
cis_x <- NULL
cis_both <- NULL

tmp <- getobj(paste0("powerSims_8ped/",fn[1]))
tmp <- tmp[["mm_res"]]
totalRes <- data.frame(matrix(NA,nrow=8*length(fn),ncol=ncol(tmp)))
colnames(totalRes) <- colnames(tmp)

for(i in 1:length(fn)){
  dat <- getobj(paste0("powerSims_8ped/",fn[i]))
  tosv <- dat[["ci_x"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  tosv$comp <- c("x","resid")

  cis_x <- rbind(cis_x,tosv)

  tosv <- dat[["ci_auto"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  tosv$comp <- c("auto","resid")

  cis_auto <- rbind(cis_auto,tosv)

  tosv <- dat[["ci_both"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  tosv$comp <- c("x","auto","resid")

  cis_both <- rbind(cis_both,tosv)

  ## and get pvalue info
  totalRes[(8*i-7):(8*i),] <- dat[["mm_res"]]
}

cis <- list("cis_x"=cis_x,"cis_auto"=cis_auto,"cis_both"=cis_both)
save(cis,file="powerSims_8ped/varCompCIs_powerSims.RData")

dim(totalRes) # 6000000 13
ftable(totalRes$model,totalRes$sigma_a,totalRes$sigma_x) # 375K for each config

save(totalRes,file="powerSims_8ped/mmRes_powerSims.RData")

rm(list=ls())


#####
# 58. Power estimates, 8 ped: make power graph

library(GWASTools); library(ggplot2)
source("powerGraph.R")

### NOTE: this is the same file name as below
dat <- getobj("mmRes_powerSims.RData")
head(dat)

dat2 <- getobj("mmRes_powerSims_day2.RData")
head(dat2)

datF <- rbind(dat,dat2)
sum(duplicated(datF$pval)) # 64064
datF <- datF[!duplicated(datF$pval),]
dim(datF) # 263936 13

dat <- datF

saveas(dat,"mmRes_8ped_powerSims.RData")

toPl <- dat[dat$sigma_x=="0.3"&dat$sigma_a=="0.3",]
dim(toPl) # 65992 13
table(toPl$causal) # 32996 0, 31996 0.05 and 1000 of 0.07
toPl <- toPl[!is.element(toPl$causal,"0.07"),]
toPl$causal[toPl$causal=="0.05"] <- TRUE
toPl$causal[toPl$causal==0] <- FALSE
toPl$beta1=0.05
sa3sx3=powerGraph(toPl,b=0.05,fn="powerGraph_sA03sX03.pdf",xlims=c(0,0.05),values=TRUE)
sa3sx3$sigma_a=paste0("sigma_a=",0.3)
sa3sx3$sigma_x=paste0("sigma_x=",0.3)

toPl <- dat[dat$sigma_x=="0.5"&dat$sigma_a=="0.3",]
dim(toPl) # 42K 13
table(toPl$causal)
toPl <- toPl[!is.element(toPl$causal,"0.07"),]
toPl$causal[toPl$causal=="0.05"] <- TRUE
toPl$causal[toPl$causal==0] <- FALSE
toPl$beta1=0.05
sx5sa3=powerGraph(toPl,b=0.05,fn="powerGraph_sA05sX03.pdf",xlims=c(0,0.05),values=TRUE)
sx5sa3$sigma_a=paste0("sigma_a=",0.3)
sx5sa3$sigma_x=paste0("sigma_x=",0.5)

toPl <- dat[dat$sigma_x=="0.5"&dat$sigma_a=="0.5",]
dim(toPl) # 42K 13
table(toPl$causal)
toPl <- toPl[!is.element(toPl$causal,"0.07"),]
toPl$causal[toPl$causal=="0.05"] <- TRUE
toPl$causal[toPl$causal==0] <- FALSE
toPl$beta1=0.05
sx5sa5=powerGraph(toPl,b=0.05,fn="powerGraph_sA05sX03.pdf",xlims=c(0,0.05),values=TRUE)
sx5sa5$sigma_a=paste0("sigma_a=",0.5)
sx5sa5$sigma_x=paste0("sigma_x=",0.5)

toPl <- dat[dat$sigma_x=="0.3"&dat$sigma_a=="0.5",]
dim(toPl) # 42K 13
table(toPl$causal)
toPl <- toPl[!is.element(toPl$causal,"0.07"),]
toPl$causal[toPl$causal=="0.05"] <- TRUE
toPl$causal[toPl$causal==0] <- FALSE
toPl$beta1=0.05
sx3sa5=powerGraph(toPl,b=0.05,fn="powerGraph_sA05sX03.pdf",xlims=c(0,0.05),values=TRUE)
sx3sa5$sigma_a=paste0("sigma_a=",0.5)
sx3sa5$sigma_x=paste0("sigma_x=",0.3)

library(ggplot2); library(reshape)
tot <- rbind(sa3sx3,sx5sa3,sx3sa5,sx5sa5)

summary(tot$autoiter); summary(tot$xiter); summary(tot$bothiter)
# all have 8K iterations

tot$x_fp_rate <- tot$x_fp/tot$xiter
tot$x_tp_rate <- tot$x_tp/tot$xiter
tot$auto_fp_rate <- tot$auto_fp/tot$autoiter
tot$auto_tp_rate <- tot$auto_tp/tot$autoiter
tot$both_fp_rate <- tot$both_fp/tot$bothiter
tot$both_tp_rate <- tot$both_tp/tot$bothiter

both <- tot[,c("sigma_a","sigma_x","both_fp_rate","both_tp_rate")]
both$model <- "both"
both$fp <- both$both_fp_rate
both$tp <- both$both_tp_rate

auto <- tot[,c("sigma_a","sigma_x","auto_fp_rate","auto_tp_rate")]
auto$model <- "auto"
auto$fp <- auto$auto_fp_rate
auto$tp <- auto$auto_tp_rate

x <- tot[,c("sigma_a","sigma_x","x_fp_rate","x_tp_rate")]
x$model <- "x"
x$fp <- x$x_fp_rate
x$tp <- x$x_tp_rate

toPl <- rbind(x[,c("sigma_x","sigma_a","model","tp","fp")],
              auto[,c("sigma_x","sigma_a","model","tp","fp")],
              both[,c("sigma_x","sigma_a","model","tp","fp")])

toPl <- toPl[toPl$fp<=0.05,]

pdf("power_beta05_8ped_8kiters.pdf",width=11,height=11)
ggplot(toPl,aes(x=fp,y=tp,color=model)) + geom_line(size=1.5,aes(linetype=model)) + facet_wrap(sigma_a~sigma_x) + theme_bw() +
  ylab("True Positive Rate") + xlab("False Positive Rate") +ggtitle(expression(paste(beta,"=0.05")))+
  theme(axis.text=element_text(size=16),axis.title=element_text(size=18),title=element_text(size=20),
        legend.text=element_text(size=16),strip.text = element_text(size=16))
dev.off()

rm(list=ls())


#####
# 59. Power estimates, 8 ped fem: process results

# called power_simulations_v1.R
# stored results in powerSims_8ped_fem/

library(GWASTools); library(ggplot2)
setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")

fn <- list.files("powerSims_8ped_fem/")
length(fn)

cis_auto <- NULL
cis_x <- NULL
cis_both <- NULL

tmp <- getobj(paste0("powerSims_8ped_fem/",fn[1]))
tmp <- tmp[["mm_res"]]
totalRes <- data.frame(matrix(NA,nrow=8*length(fn),ncol=ncol(tmp)))
colnames(totalRes) <- colnames(tmp)

for(i in 1:length(fn)){
  dat <- getobj(paste0("powerSims_8ped_fem/",fn[i]))
  tosv <- dat[["ci_x"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  tosv$comp <- c("x","resid")

  cis_x <- rbind(cis_x,tosv)

  tosv <- dat[["ci_auto"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  tosv$comp <- c("auto","resid")

  cis_auto <- rbind(cis_auto,tosv)

  tosv <- dat[["ci_both"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  tosv$comp <- c("x","auto","resid")

  cis_both <- rbind(cis_both,tosv)

  ## and get pvalue info
  totalRes[(8*i-7):(8*i),] <- dat[["mm_res"]]

  if(i%%100==0){print(i)}
}

cis <- list("cis_x"=cis_x,"cis_auto"=cis_auto,"cis_both"=cis_both)
save(cis,file="powerSims_8ped_fem/varCompCIs_powerSims.RData")

dim(totalRes) # 6000000 13
ftable(totalRes$model,totalRes$sigma_a,totalRes$sigma_x) # 375K for each config

save(totalRes,file="powerSims_8ped_fem/mmRes_powerSims.RData")

rm(list=ls())


#####
# 60. Power estimates, 8 ped fem: make power graph

library(GWASTools); library(ggplot2)
source("powerGraph.R")

### NOTE: this is the same file name as above
dat <- getobj("mmRes_powerSims.RData")
head(dat)

toPl <- dat[dat$sigma_x=="0.3"&dat$sigma_a=="0.3",]
dim(toPl) # 22608 13
table(toPl$causal) # 11304 0, 10304 0.05 and 1000 of 0.07
toPl <- toPl[!is.element(toPl$causal,"0.07"),]
toPl$causal[toPl$causal=="0.05"] <- TRUE
toPl$causal[toPl$causal==0] <- FALSE
toPl$beta1=0.05
sa3sx3=powerGraph(toPl,b=0.05,fn="powerGraph_sA03sX03.pdf",xlims=c(0,0.05),values=TRUE)
sa3sx3$sigma_a=paste0("sigma_a=",0.3)
sa3sx3$sigma_x=paste0("sigma_x=",0.3)

toPl <- dat[dat$sigma_x=="0.5"&dat$sigma_a=="0.3",]
dim(toPl) # 22288 13
table(toPl$causal)
toPl <- toPl[!is.element(toPl$causal,"0.07"),]
toPl$causal[toPl$causal=="0.05"] <- TRUE
toPl$causal[toPl$causal==0] <- FALSE
toPl$beta1=0.05
sx5sa3=powerGraph(toPl,b=0.05,fn="powerGraph_sA05sX03.pdf",xlims=c(0,0.05),values=TRUE)
sx5sa3$sigma_a=paste0("sigma_a=",0.3)
sx5sa3$sigma_x=paste0("sigma_x=",0.5)

toPl <- dat[dat$sigma_x=="0.5"&dat$sigma_a=="0.5",]
dim(toPl) # 23328 13
table(toPl$causal)
toPl <- toPl[!is.element(toPl$causal,"0.07"),]
toPl$causal[toPl$causal=="0.05"] <- TRUE
toPl$causal[toPl$causal==0] <- FALSE
toPl$beta1=0.05
sx5sa5=powerGraph(toPl,b=0.05,fn="powerGraph_sA05sX03.pdf",xlims=c(0,0.05),values=TRUE)
sx5sa5$sigma_a=paste0("sigma_a=",0.5)
sx5sa5$sigma_x=paste0("sigma_x=",0.5)

toPl <- dat[dat$sigma_x=="0.3"&dat$sigma_a=="0.5",]
dim(toPl) # 42K 13
table(toPl$causal)
toPl <- toPl[!is.element(toPl$causal,"0.07"),]
toPl$causal[toPl$causal=="0.05"] <- TRUE
toPl$causal[toPl$causal==0] <- FALSE
toPl$beta1=0.05
sx3sa5=powerGraph(toPl,b=0.05,fn="powerGraph_sA05sX03.pdf",xlims=c(0,0.05),values=TRUE)
sx3sa5$sigma_a=paste0("sigma_a=",0.5)
sx3sa5$sigma_x=paste0("sigma_x=",0.3)

library(ggplot2); library(reshape)
tot <- rbind(sa3sx3,sx5sa3,sx3sa5,sx5sa5)
tot$x_fp_rate <- tot$x_fp/tot$xiter
tot$x_tp_rate <- tot$x_tp/tot$xiter
tot$auto_fp_rate <- tot$auto_fp/tot$autoiter
tot$auto_tp_rate <- tot$auto_tp/tot$autoiter
tot$both_fp_rate <- tot$both_fp/tot$bothiter
tot$both_tp_rate <- tot$both_tp/tot$bothiter

both <- tot[,c("sigma_a","sigma_x","both_fp_rate","both_tp_rate")]
both$model <- "both"
both$fp <- both$both_fp_rate
both$tp <- both$both_tp_rate

auto <- tot[,c("sigma_a","sigma_x","auto_fp_rate","auto_tp_rate")]
auto$model <- "auto"
auto$fp <- auto$auto_fp_rate
auto$tp <- auto$auto_tp_rate

x <- tot[,c("sigma_a","sigma_x","x_fp_rate","x_tp_rate")]
x$model <- "x"
x$fp <- x$x_fp_rate
x$tp <- x$x_tp_rate

toPl <- rbind(x[,c("sigma_x","sigma_a","model","tp","fp")],
              auto[,c("sigma_x","sigma_a","model","tp","fp")],
              both[,c("sigma_x","sigma_a","model","tp","fp")])

toPl <- toPl[toPl$fp<=0.05,]

pdf("power_beta05_8ped_fem_2500iters.pdf",width=11,height=11)
ggplot(toPl,aes(x=fp,y=tp,color=model)) + geom_line(size=1.5,aes(linetype=model)) + facet_wrap(sigma_a~sigma_x) + theme_bw() +
  ylab("True Positive Rate") + xlab("False Positive Rate") +ggtitle(expression(paste(beta,"=0.05")))+
  theme(axis.text=element_text(size=16),axis.title=element_text(size=18),title=element_text(size=20),
        legend.text=element_text(size=16),strip.text = element_text(size=16))
dev.off()

rm(list=ls())



# call 22500 sims for the non-fem pedigree, ie male centric pedigree
# qsub -q olga.q -p -50 -t 1-22500 -N pow8ped batch_power_simulations.sh


#####
# 61. Power simulations, auto SNP

# call 5000 sims for the fem pedigree, ie balanced pedigree
# qsub -q thornton.q -p -50 -t 1-2500 -N powBalAuto batch_power_simulations.sh
# which calls power_simulations_v3.R: power sims for auto snp with fem/balanced pedigree
# stores results in powerSims_8ped_fem/autoSNP

# call 2500 sims for the male pedigree
# qsub -q olga.q -p -50 -t 1-2000 -N powMalAuto batch_power_simulations.sh
# which calls power_simulations_v4.R: power sims for auto snp with male-centric pedigree
# stores results in powerSims_8ped/autoSNP .txt files for the mmRes, can join these by unix


#####
# 62. Make power graph powerSims_8ped/autoSNP results

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x/")
source("powerGraph.R")

res <- read.table("powerSims_8ped/autoSNP/allSims.txt",header=TRUE,as.is=TRUE)
dim(res) # 268718     13
# remove the header rows from all the files

headRw <- which(res$snpID=="snpID")
res <- res[-headRw,]
dim(res) # 260576     13

table(res$causal,res$model)
#      auto  both linear     x
#0    32572 32572  32572 32572
#0.05 32572 32572  32572 32572

dat <- res

toPl <- dat[dat$sigma_x=="0.3"&dat$sigma_a=="0.3",]
dim(toPl) # 65144 13
table(toPl$causal) # 32572 0, 32572 0.05
toPl$causal[toPl$causal=="0.05"] <- TRUE
toPl$causal[toPl$causal==0] <- FALSE
toPl$beta1=0.05
toPl$causal <- as.logical(toPl$causal)

sa3sx3=powerGraph(toPl,b=0.05,fn="powerGraph_sA03sX03.pdf",xlims=c(0,0.05),values=TRUE)
sa3sx3$sigma_a=paste0("sigma_a=",0.3)
sa3sx3$sigma_x=paste0("sigma_x=",0.3)

toPl <- dat[dat$sigma_x=="0.5"&dat$sigma_a=="0.3",]
dim(toPl) # 65144 13
table(toPl$causal)
toPl$causal[toPl$causal=="0.05"] <- TRUE
toPl$causal[toPl$causal==0] <- FALSE
toPl$beta1=0.05
toPl$causal <- as.logical(toPl$causal)
sx5sa3=powerGraph(toPl,b=0.05,fn="powerGraph_sA05sX03.pdf",xlims=c(0,0.05),values=TRUE)
sx5sa3$sigma_a=paste0("sigma_a=",0.3)
sx5sa3$sigma_x=paste0("sigma_x=",0.5)

toPl <- dat[dat$sigma_x=="0.5"&dat$sigma_a=="0.5",]
dim(toPl) # 65144 13
table(toPl$causal)

toPl$causal[toPl$causal=="0.05"] <- TRUE
toPl$causal[toPl$causal==0] <- FALSE
toPl$beta1=0.05
toPl$causal <- as.logical(toPl$causal)
sx5sa5=powerGraph(toPl,b=0.05,fn="powerGraph_sA05sX03.pdf",xlims=c(0,0.05),values=TRUE)
sx5sa5$sigma_a=paste0("sigma_a=",0.5)
sx5sa5$sigma_x=paste0("sigma_x=",0.5)

toPl <- dat[dat$sigma_x=="0.3"&dat$sigma_a=="0.5",]
dim(toPl) # 65144 13
table(toPl$causal)
toPl <- toPl[!is.element(toPl$causal,"0.07"),]
toPl$causal[toPl$causal=="0.05"] <- TRUE
toPl$causal[toPl$causal==0] <- FALSE
toPl$beta1=0.05
toPl$causal <- as.logical(toPl$causal)
sx3sa5=powerGraph(toPl,b=0.05,fn="powerGraph_sA05sX03.pdf",xlims=c(0,0.05),values=TRUE)
sx3sa5$sigma_a=paste0("sigma_a=",0.5)
sx3sa5$sigma_x=paste0("sigma_x=",0.3)

library(ggplot2); library(reshape)
tot <- rbind(sa3sx3,sx5sa3,sx3sa5,sx5sa5)
tot$x_fp_rate <- tot$x_fp/tot$xiter
tot$x_tp_rate <- tot$x_tp/tot$xiter
tot$auto_fp_rate <- tot$auto_fp/tot$autoiter
tot$auto_tp_rate <- tot$auto_tp/tot$autoiter
tot$both_fp_rate <- tot$both_fp/tot$bothiter
tot$both_tp_rate <- tot$both_tp/tot$bothiter

both <- tot[,c("sigma_a","sigma_x","both_fp_rate","both_tp_rate")]
both$model <- "both"
both$fp <- both$both_fp_rate
both$tp <- both$both_tp_rate

auto <- tot[,c("sigma_a","sigma_x","auto_fp_rate","auto_tp_rate")]
auto$model <- "auto"
auto$fp <- auto$auto_fp_rate
auto$tp <- auto$auto_tp_rate

x <- tot[,c("sigma_a","sigma_x","x_fp_rate","x_tp_rate")]
x$model <- "x"
x$fp <- x$x_fp_rate
x$tp <- x$x_tp_rate

toPl <- rbind(x[,c("sigma_x","sigma_a","model","tp","fp")],
              auto[,c("sigma_x","sigma_a","model","tp","fp")],
              both[,c("sigma_x","sigma_a","model","tp","fp")])

toPl <- toPl[toPl$fp<=0.05,]

pdf("power_beta05_8ped_autoSNP_8kiters.pdf",width=11,height=11)
ggplot(toPl,aes(x=fp,y=tp,color=model)) + geom_line(size=1.5,aes(linetype=model)) + facet_wrap(sigma_a~sigma_x) + theme_bw() +
  ylab("True Positive Rate") + xlab("False Positive Rate") +ggtitle(expression(paste(beta,"=0.05")))+
  theme(axis.text=element_text(size=16),axis.title=element_text(size=18),title=element_text(size=20),
        legend.text=element_text(size=16),strip.text = element_text(size=16))
dev.off()

rm(list=ls())


#####
# 63. Combine more var comp estimates to make updated boxplots, 8ped_fem

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x/nullSims_8ped_fem")

# want for null SNPs
# need var comp ests for mlm-x, x, and auto models

# set of results stored in mlm_x/nullSims_8ped_fem/

library(GWASTools); library(ggplot2)

cis <- getobj("nullSims_8ped_fem/varCompCIs_nullSims.RData")

f <- list.files()
remv <- grep("extremeSig",f)
f <- f[-remv]
f <- f[-1]
f <- f[-length(f)]
length(f) # 826

newResBoth <- data.frame(matrix(NA,nrow=3*826,ncol=7))
newResX <- data.frame(matrix(NA,nrow=2*826,ncol=7))
newResAuto <- data.frame(matrix(NA,nrow=2*826,ncol=7))

ctb=1
ctoth=1

for(i in f){
  dat <- get(load(i))
  toadB <- dat[["ci_both"]]
  toadB$sigma_a <- dat[["sigma_auto"]]
  toadB$sigma_x <- dat[["sigma_x"]]
  toadB$comp <- c("x","auto","resid")

  newResBoth[ctb:(ctb+2),] <- toadB

  toadA <- dat[["ci_auto"]]
  toadA$sigma_a <- dat[["sigma_auto"]]
  toadA$sigma_x <- dat[["sigma_x"]]
  toadA$comp <- c("auto","resid")

  newResAuto[ctoth:(ctoth+1),] <- toadA

  toadX <- dat[["ci_x"]]
  toadX$sigma_a <- dat[["sigma_auto"]]
  toadX$sigma_x <- dat[["sigma_x"]]
  toadX$comp <- c("x","resid")

  newResX[ctoth:(ctoth+1),] <- toadX

  ctb <- ctb+3
  ctoth <- ctoth+2
}

colnames(newResBoth) <- colnames(cis[["cis_both"]])
colnames(newResX) <- colnames(cis[["cis_x"]])
colnames(newResAuto) <- colnames(cis[["cis_auto"]])

cis[["cis_both"]] <- rbind(cis[["cis_both"]],newResBoth)
cis[["cis_x"]] <- rbind(cis[["cis_x"]],newResX)
cis[["cis_auto"]] <- rbind(cis[["cis_auto"]],newResAuto)

cib <- cis[["cis_both"]]
ftable(cib$comp,cib$sigma_a,cib$sigma_x) # between 375-500 of each configuration

cib$sigma_a=paste0("sigma_a=",cib$sigma_a)
cib$sigma_x=paste0("sigma_x=",cib$sigma_x)
cib$comp=ordered(cib$comp,levels=c("auto","x","resid"))

mns <- expand.grid(comp=c("auto","x","resid"),sigma_a=c("sigma_a=0.3","sigma_a=3"),
                   sigma_x=c("sigma_x=0.3","sigma_x=3"))
mns$comp <- ordered(mns$comp,levels=c("auto","x","resid"))
mns$y <- -0.1
mns$value <- NA
for(i in seq_len(nrow(mns))){
  mns$value[i] <- mean(cib$Est[cib$comp==mns$comp[i]&cib$sigma_a==mns$sigma_a[i]&cib$sigma_x==mns$sigma_x[i]])
}
mns <- mns[!is.nan(mns$value),]
mns$value <- format(mns$value,digits=3)
mns$sigma_x=NULL # no sigma_x variable

pdf("boxplots_varComp_autoEst_8ped_fem.pdf")
ggplot(cib,aes(x=comp,y=Est)) + geom_boxplot() + facet_grid(sigma_a~sigma_x) + theme_bw() +
  xlab("Variance Component") + ylab("Estimate") +
  ggtitle("Variance Component Estimates of 500 Iterations") +
  geom_text(data=mns, aes(x=comp, y=y, label=value), size=4)
dev.off()

rm(list=ls())


#####
# 64. Power simulations, auto SNP, extreme sigma

# call 5000 sims for the fem pedigree, ie balanced pedigree
# qsub -q thornton.q -p -50 -t 1-5000 -N balAutoExt batch_power_simulations.sh
# which calls power_simulations_v3_extremeSigma.R: power sims for auto snp with fem/balanced pedigree
# stores results in powerSims_8ped_fem/autoSNP/extremeSigma* .txt files, can join using cat *mmOnly.txt >> allmmRes.txt


# call 5000 sims for the male pedigree
# qsub -q thornton.q -p -50 -t 1-5000 -N femAutoExt batch_power_simulations.sh
# which calls power_simulations_v4_extremeSigma.R: power sims for auto snp with male-centric pedigree
# stores results in powerSims_8ped/autoSNP/extremeSigma* .txt files for the mmRes, can join these by unix


#####
# 65. Make power graph powerSims_8ped/autoSNP/extremeSigma results

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x/")
source("powerGraph.R")

dat <- read.table("powerSims_8ped/autoSNP/allmmRes.txt",header=TRUE,as.is=TRUE)
dim(dat) # 84999 13

headRw <- which(dat$snpID=="snpID")
dat <- dat[-headRw,]
dim(dat) # 80000 13

table(dat$causal, dat$model) # 10K of each

toPl <- dat[dat$sigma_x=="3"&dat$sigma_a=="0.3",]
dim(toPl) # 40000 13
table(toPl$causal) # 20K 0, 20K 0.05
toPl$causal[toPl$causal=="0.05"] <- TRUE
toPl$causal[toPl$causal==0] <- FALSE
toPl$beta1=0.05
toPl$causal <- as.logical(toPl$causal)

sa3sx3=powerGraph(toPl,b=0.05,fn="powerGraph_sA03sX03.pdf",xlims=c(0,0.05),values=TRUE)
sa3sx3$sigma_a=paste0("sigma_a=",0.3)
sa3sx3$sigma_x=paste0("sigma_x=",3)

toPl <- dat[dat$sigma_x=="0.3"&dat$sigma_a=="3",]
dim(toPl) # 65144 13
table(toPl$causal)
toPl$causal[toPl$causal=="0.05"] <- TRUE
toPl$causal[toPl$causal==0] <- FALSE
toPl$beta1=0.05
toPl$causal <- as.logical(toPl$causal)
sx5sa3=powerGraph(toPl,b=0.05,fn="powerGraph_sA05sX03.pdf",xlims=c(0,0.05),values=TRUE)
sx5sa3$sigma_a=paste0("sigma_a=",3)
sx5sa3$sigma_x=paste0("sigma_x=",0.3)


library(ggplot2); library(reshape)
tot <- rbind(sa3sx3,sx5sa3)
tot$x_fp_rate <- tot$x_fp/tot$xiter
tot$x_tp_rate <- tot$x_tp/tot$xiter
tot$auto_fp_rate <- tot$auto_fp/tot$autoiter
tot$auto_tp_rate <- tot$auto_tp/tot$autoiter
tot$both_fp_rate <- tot$both_fp/tot$bothiter
tot$both_tp_rate <- tot$both_tp/tot$bothiter

both <- tot[,c("sigma_a","sigma_x","both_fp_rate","both_tp_rate")]
both$model <- "both"
both$fp <- both$both_fp_rate
both$tp <- both$both_tp_rate

auto <- tot[,c("sigma_a","sigma_x","auto_fp_rate","auto_tp_rate")]
auto$model <- "auto"
auto$fp <- auto$auto_fp_rate
auto$tp <- auto$auto_tp_rate

x <- tot[,c("sigma_a","sigma_x","x_fp_rate","x_tp_rate")]
x$model <- "x"
x$fp <- x$x_fp_rate
x$tp <- x$x_tp_rate

toPl <- rbind(x[,c("sigma_x","sigma_a","model","tp","fp")],
              auto[,c("sigma_x","sigma_a","model","tp","fp")],
              both[,c("sigma_x","sigma_a","model","tp","fp")])

toPl <- toPl[toPl$fp<=0.05,]

pdf("power_beta05_8ped_autoSNP_extrememSigma_10kiters.pdf",width=11,height=11)
ggplot(toPl,aes(x=fp,y=tp,color=model)) + geom_line(size=1.5,aes(linetype=model)) + facet_wrap(sigma_a~sigma_x) + theme_bw() +
  ylab("True Positive Rate") + xlab("False Positive Rate") +ggtitle(expression(paste(beta,"=0.05")))+
  theme(axis.text=element_text(size=16),axis.title=element_text(size=18),title=element_text(size=20),
        legend.text=element_text(size=16),strip.text = element_text(size=16))
dev.off()

rm(list=ls())


#####
# 66. Null simulations, 10K, 8 person ped fem, 8 person ped, auto + x chr SNP

# need to do for both auto and x chr SNP
# null_simulation_v11.R: extreme sigma, 8 ped, auto SNP
# null_simulations_v10.R: extreme sigma, 8 ped fem, X SNP
# null_simulations_v9.R: extreme sigma, 8 ped fem, auto SNP
# null_simulations_v8.R: null sims, 8 ped - auto SNP
# null_simulations_v7.R: 16 person ped
# null_simulations_v6.R: 8 person ped, fem - auto SNP
# null_simulations_v5.R: 8 person ped, fem - X chr SNP

# submitted: qsub -q thornton.q -p -50 -t 1-500 -N nullFemAuto batch_null_simulations.sh
# which calls null_simulations_v6.R and stores results in
# mlm_x/nullSims_8ped_fem/autoSNP/

# and: qsub -q thornton.q -p -50 -t 1-500 -N nullFemX batch_null_simulations.sh
# which calls null_simulations_v5.R and stores results in
# mlm_x/nullSims_8ped_fem/

# qsub -q thornton.q -p -50 -t 1-500 -N nullAuto batch_null_simulations.sh
# which calls null_simulations_v8.R and stores results in:
# mlm_x/nullSims_8ped/autoSNP

# qsub -q thornton.q -p -50 -t 1-500 -N nullAutoFemExt batch_null_simulations.sh
# which calls null_simulations_v11.R and stores results in:
# mlm_x/nullSims_8ped/autoSNP


#####
# 67. Make power graph powerSims_8ped_fem/autoSNP/extremeSigma results

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x/")
source("powerGraph.R")

dat <- read.table("powerSims_8ped_fem/autoSNP/allmmRes_extremeSigma.txt",header=TRUE,as.is=TRUE)
dim(dat) # 84999 13

headRw <- which(dat$snpID=="snpID")
dat <- dat[-headRw,]
dim(dat) # 80000 13

table(dat$causal, dat$model) # 10K of each

toPl <- dat[dat$sigma_x=="3"&dat$sigma_a=="0.3",]
dim(toPl) # 40000 13
table(toPl$causal) # 20K 0, 20K 0.05
toPl$causal[toPl$causal=="0.05"] <- TRUE
toPl$causal[toPl$causal==0] <- FALSE
toPl$beta1=0.05
toPl$causal <- as.logical(toPl$causal)

sa3sx3=powerGraph(toPl,b=0.05,fn="powerGraph_sA03sX03.pdf",xlims=c(0,0.05),values=TRUE)
sa3sx3$sigma_a=paste0("sigma_a=",0.3)
sa3sx3$sigma_x=paste0("sigma_x=",3)

toPl <- dat[dat$sigma_x=="0.3"&dat$sigma_a=="3",]
dim(toPl) # 65144 13
table(toPl$causal)
toPl$causal[toPl$causal=="0.05"] <- TRUE
toPl$causal[toPl$causal==0] <- FALSE
toPl$beta1=0.05
toPl$causal <- as.logical(toPl$causal)
sx5sa3=powerGraph(toPl,b=0.05,fn="powerGraph_sA05sX03.pdf",xlims=c(0,0.05),values=TRUE)
sx5sa3$sigma_a=paste0("sigma_a=",3)
sx5sa3$sigma_x=paste0("sigma_x=",0.3)

library(ggplot2); library(reshape)
tot <- rbind(sa3sx3,sx5sa3)
tot$x_fp_rate <- tot$x_fp/tot$xiter
tot$x_tp_rate <- tot$x_tp/tot$xiter
tot$auto_fp_rate <- tot$auto_fp/tot$autoiter
tot$auto_tp_rate <- tot$auto_tp/tot$autoiter
tot$both_fp_rate <- tot$both_fp/tot$bothiter
tot$both_tp_rate <- tot$both_tp/tot$bothiter

both <- tot[,c("sigma_a","sigma_x","both_fp_rate","both_tp_rate")]
both$model <- "both"
both$fp <- both$both_fp_rate
both$tp <- both$both_tp_rate

auto <- tot[,c("sigma_a","sigma_x","auto_fp_rate","auto_tp_rate")]
auto$model <- "auto"
auto$fp <- auto$auto_fp_rate
auto$tp <- auto$auto_tp_rate

x <- tot[,c("sigma_a","sigma_x","x_fp_rate","x_tp_rate")]
x$model <- "x"
x$fp <- x$x_fp_rate
x$tp <- x$x_tp_rate

toPl <- rbind(x[,c("sigma_x","sigma_a","model","tp","fp")],
              auto[,c("sigma_x","sigma_a","model","tp","fp")],
              both[,c("sigma_x","sigma_a","model","tp","fp")])

toPl <- toPl[toPl$fp<=0.05,]

pdf("power_beta05_8ped_fem_autoSNP_extrememSigma_10kiters.pdf",width=11,height=11)
ggplot(toPl,aes(x=fp,y=tp,color=model)) + geom_line(size=1.5,aes(linetype=model)) + facet_wrap(sigma_a~sigma_x) + theme_bw() +
  ylab("True Positive Rate") + xlab("False Positive Rate") +ggtitle(expression(paste(beta,"=0.05")))+
  theme(axis.text=element_text(size=16),axis.title=element_text(size=18),title=element_text(size=20),
        legend.text=element_text(size=16),strip.text = element_text(size=16))
dev.off()

rm(list=ls())


#####
# 68. Type I error, 8ped_fem autoSNP

library(GWASTools)
setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")

fn <- list.files("nullSims_8ped_fem/autoSNP/")
length(fn) # 1933

cis_auto <- NULL
cis_x <- NULL
cis_both <- NULL
for(i in 1:length(fn)){
  dat <- getobj(paste0("nullSims_8ped_fem/autoSNP/",fn[i]))
  tosv <- dat[["ci_x"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  tosv$comp <- c("x","resid")

  cis_x <- rbind(cis_x,tosv)

  tosv <- dat[["ci_auto"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  tosv$comp <- c("auto","resid")

  cis_auto <- rbind(cis_auto,tosv)

  tosv <- dat[["ci_both"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  tosv$comp <- c("x","auto","resid")

  cis_both <- rbind(cis_both,tosv)
}

cis <- list("cis_x"=cis_x,"cis_auto"=cis_auto,"cis_both"=cis_both)
save(cis,file="nullSims_8ped_fem/varCompCIs_autoSNP_nullSims.RData")

res1 <- fn[substr(fn,1,4)=="res_"]
res2 <- fn[substr(fn,1,5)=="res2_"]
res3 <- fn[substr(fn,1,5)=="res3_"]
res4 <- fn[substr(fn,1,5)=="res4_"]

totalRes <- NULL
# read in first 25 iters of res, res2, res3, res4
for(i in 1:100){
  dat <- getobj(paste0("nullSims_8ped_fem/autoSNP/",res1[i]))
  tosv <- dat[["mm_res"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  totalRes <- rbind(totalRes,tosv)

  dat <- getobj(paste0("nullSims_8ped_fem/autoSNP/",res2[i]))
  tosv <- dat[["mm_res"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  totalRes <- rbind(totalRes,tosv)

  dat <- getobj(paste0("nullSims_8ped_fem/autoSNP/",res3[i]))
  tosv <- dat[["mm_res"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  totalRes <- rbind(totalRes,tosv)

  dat <- getobj(paste0("nullSims_8ped_fem/autoSNP/",res4[i]))
  tosv <- dat[["mm_res"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  totalRes <- rbind(totalRes,tosv)
  }

dim(totalRes) # 120000000 13
ftable(totalRes$model,totalRes$sigma_a,totalRes$sigma_x) # 750K for each config

save(totalRes,file="nullSims_8ped_fem/1M_mmRes_autoSNP_nullSims.RData")

# get type I error for various alpha values
tyIerr <- expand.grid(alpha=c(1e-3,1e-4,1e-5),model=c("both","auto","x"),sigma_a=c(0.3,0.5),sigma_x=c(0.3,0.5))
tyIerr$tyIerrRate <- NA
tyIerr$falsePos <- NA
tyIerr$iters <- NA

for(i in 1:nrow(tyIerr)){
  smallRes <- totalRes[is.element(totalRes$model,tyIerr$model[i])&is.element(totalRes$sigma_a,tyIerr$sigma_a[i])&
                       is.element(totalRes$sigma_x,tyIerr$sigma_x[i]),]
  iters <- nrow(smallRes)
  tyIerr$tyIerrRate[i] <- sum(as.numeric(smallRes$pval)<tyIerr$alpha[i])/iters
  tyIerr$iters[i] <- iters
  tyIerr$falsePos[i] <- sum(as.numeric(smallRes$pval)<tyIerr$alpha[i])
}

library(xtable)
xtable(tyIerr,digits=5)

save(tyIerr,file="nullSims_8ped_fem/1M_tyIerr_autoSNP.RData")


## make a type I error plot

tyIerr <- expand.grid(model=c("auto","both","x","linear"),alpha=c(5e-03,5e-04,5e-05),
                      sigma_x=c(0.3,0.5),sigma_a=c(0.3,0.5))
tyIerr$est <- NA
tyIerr$lower <- NA; tyIerr$upper <- NA
for(i in seq_len(nrow(tyIerr))){
  mmResF <- totalRes[totalRes$model==tyIerr$model[i]&totalRes$sigma_x==tyIerr$sigma_x[i]&
                       totalRes$sigma_a==tyIerr$sigma_a[i],]
  tyIerr$est[i] <- sum(mmResF$pval<tyIerr$alpha[i])/nrow(mmResF)
  stdErr <- sqrt(tyIerr$alpha[i]*(1-tyIerr$alpha[i])/nrow(mmResF))
  tyIerr$lower[i] <- tyIerr$est[i]-1.96*stdErr
  tyIerr$upper[i] <- tyIerr$est[i]+1.96*stdErr
}

tyIerr$sigma_a <- paste0("sigma_a=",tyIerr$sigma_a)
tyIerr$sigma_x <- paste0("sigma_x=",tyIerr$sigma_x)

# make the model an ordered factor so it prints both, x, auto, linear
tyIerr$model <- ordered(tyIerr$model,levels=c("both","x","auto","linear"))
tyIerrSm <- tyIerr[tyIerr$alpha==5e-05,]

pdf("typeIErr_8ped_fem_5e05_autoSNP.pdf")
ggplot(tyIerrSm,aes(x=model,y=est)) +geom_point()+ geom_pointrange(aes(ymin=lower,ymax=upper))+
  theme_bw() + geom_hline(aes(yintercept=5e-05,color="gray"))+
  facet_grid(sigma_x~sigma_a) +
  ylab("Type I Error Rate") + ggtitle(expression(paste("Type I Error Rate, ", alpha, "=5e-04")))
dev.off()

tyIerrSm <- tyIerr[tyIerr$alpha==5e-05&!is.nan(tyIerr$est)&!is.element(tyIerr$model,"linear"),]
pdf("typeIErr_8ped_fem_autoSNP_5e05_zoom.pdf")
ggplot(tyIerrSm,aes(x=model,y=est)) +geom_point()+ geom_pointrange(aes(ymin=lower,ymax=upper))+
  theme_bw() + geom_hline(aes(yintercept=5e-05,color="gray"))+
  facet_grid(sigma_x~sigma_a) +
  ylab("Type I Error Rate") + ggtitle(expression(paste("Type I Error Rate, ", alpha, "=5e-05")))
dev.off()

rm(list=ls())


#####
# 69. Type I error, 8ped autoSNP, extreme sigma

library(GWASTools)
setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")

fn <- list.files("nullSims_8ped/autoSNP/")
fn <- fn[grep("extremeSig",fn)]
length(fn) # 497

cis_auto <- NULL
cis_x <- NULL
cis_both <- NULL
for(i in 1:length(fn)){
  dat <- getobj(paste0("nullSims_8ped/autoSNP/",fn[i]))
  tosv <- dat[["ci_x"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  tosv$comp <- c("x","resid")

  cis_x <- rbind(cis_x,tosv)

  tosv <- dat[["ci_auto"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  tosv$comp <- c("auto","resid")

  cis_auto <- rbind(cis_auto,tosv)

  tosv <- dat[["ci_both"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  tosv$comp <- c("x","auto","resid")

  cis_both <- rbind(cis_both,tosv)
}

cis <- list("cis_x"=cis_x,"cis_auto"=cis_auto,"cis_both"=cis_both)
save(cis,file="nullSims_8ped/varCompCIs_autoSNP_extremeSigma_nullSims.RData")

res1 <- fn[substr(fn,1,4)=="res_"]
res2 <- fn[substr(fn,1,5)=="res2_"]

totalRes <- NULL
# read in first 25 iters of res, res2, res3, res4
for(i in 1:50){
  dat <- getobj(paste0("nullSims_8ped/autoSNP/",res1[i]))
  tosv <- dat[["mm_res"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  totalRes <- rbind(totalRes,tosv)

  dat <- getobj(paste0("nullSims_8ped/autoSNP/",res2[i]))
  tosv <- dat[["mm_res"]]
  tosv$sigma_a <- dat[["sigma_auto"]]
  tosv$sigma_x <- dat[["sigma_x"]]
  totalRes <- rbind(totalRes,tosv)
}

dim(totalRes) # 6000000 13
ftable(totalRes$model,totalRes$sigma_a,totalRes$sigma_x) # 750K for each config

save(totalRes,file="nullSims_8ped/mmRes_autoSNP_extremeSigma_nullSims.RData")

## make a type I error plot

tyIerr <- expand.grid(model=c("auto","both","x","linear"),alpha=c(5e-03,5e-04,5e-05),
                      sigma_x=c(0.3,3),sigma_a=c(0.3,3))
tyIerr$est <- NA
tyIerr$lower <- NA; tyIerr$upper <- NA
for(i in seq_len(nrow(tyIerr))){
  mmResF <- totalRes[totalRes$model==tyIerr$model[i]&totalRes$sigma_x==tyIerr$sigma_x[i]&
                       totalRes$sigma_a==tyIerr$sigma_a[i],]
  tyIerr$est[i] <- sum(mmResF$pval<tyIerr$alpha[i])/nrow(mmResF)
  stdErr <- sqrt(tyIerr$alpha[i]*(1-tyIerr$alpha[i])/nrow(mmResF))
  tyIerr$lower[i] <- tyIerr$est[i]-1.96*stdErr
  tyIerr$upper[i] <- tyIerr$est[i]+1.96*stdErr
}

tyIerr$sigma_a <- paste0("sigma_a=",tyIerr$sigma_a)
tyIerr$sigma_x <- paste0("sigma_x=",tyIerr$sigma_x)

tyIerr <- tyIerr[!is.nan(tyIerr$est),]

# make the model an ordered factor so it prints both, x, auto, linear
tyIerr$model <- ordered(tyIerr$model,levels=c("both","x","auto","linear"))
tyIerrSm <- tyIerr[tyIerr$alpha==5e-05,]

pdf("typeIErr_8ped_extremeSig_5e05_autoSNP.pdf")
ggplot(tyIerrSm,aes(x=model,y=est)) +geom_point()+ geom_pointrange(aes(ymin=lower,ymax=upper))+
  theme_bw() + geom_hline(aes(yintercept=5e-05,color="gray"))+
  facet_wrap(sigma_x~sigma_a) +
  ylab("Type I Error Rate") + ggtitle(expression(paste("Type I Error Rate, ", alpha, "=5e-05")))
dev.off()

tyIerrSm <- tyIerr[tyIerr$alpha==5e-05&!is.nan(tyIerr$est)&!is.element(tyIerr$model,"linear"),]
pdf("typeIErr_8ped_extremeSig_autoSNP_5e05_zoom.pdf")
ggplot(tyIerrSm,aes(x=model,y=est)) +geom_point()+ geom_pointrange(aes(ymin=lower,ymax=upper))+
  theme_bw() + geom_hline(aes(yintercept=5e-05,color="gray"))+
  facet_wrap(sigma_x~sigma_a) +
  ylab("Type I Error Rate") + ggtitle(expression(paste("Type I Error Rate, ", alpha, "=5e-05")))
dev.off()

rm(list=ls())


#####
# 70. Type I error, 8ped, autoSNP

# read in various results tables, combine into 1

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")

mmRes1 <- read.table("nullSims_8ped/autoSNP/mmRes_nullSims_autoSNP.txt",header=T,as.is=T)
mmRes2 <- read.table("nullSims_8ped/autoSNP/allSigmammRes.txt",header=T,as.is=T)

mmRes_autoSNP_nullSims.RData


#####
# 71. Type I error, 8ped fem, X SNP

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")

mmRes <- get(load("nullSims_8ped_fem/mmRes_nullSims.RData"))
tmp <- read.table("nullSims_8ped/allmmRes.txt",header=T,as.is=T,nrows=1e7)
tmp <- tmp[!is.element(tmp$causal,"causal"),] # remove those header rows from combining the results via cat >>
dim(tmp) # 9999959

allRes <- rbind(tmp,mmRes)
save(allRes,file="nullSims_8ped_fem/1M_mmRes_nullSims_8ped_fem.RData")

# get type I error for various alpha values
tyIerr <- expand.grid(alpha=c(1e-3,1e-4,1e-5),model=c("both","auto","x"),sigma_a=c(0.3,0.5),sigma_x=c(0.3,0.5))
tyIerr$tyIerrRate <- NA
tyIerr$falsePos <- NA
tyIerr$iters <- NA

for(i in 1:nrow(tyIerr)){
  smallRes <- allRes[is.element(allRes$model,tyIerr$model[i])&is.element(allRes$sigma_a,tyIerr$sigma_a[i])&
                       is.element(allRes$sigma_x,tyIerr$sigma_x[i]),]
  iters <- nrow(smallRes)
  tyIerr$tyIerrRate[i] <- sum(as.numeric(smallRes$pval)<tyIerr$alpha[i])/iters
  tyIerr$iters[i] <- iters
  tyIerr$falsePos[i] <- sum(as.numeric(smallRes$pval)<tyIerr$alpha[i])
}

library(xtable)
xtable(tyIerr,digits=5)

save(tyIerr,file="nullSims_8ped_fem/1M_tyIerr.RData")

rm(list=ls())


#####
# 72. Null simulations, 8ped fem, X SNP -- but for var comp estimates

# null_simulations_v5_varComps.R
# call 10K of these
# cat varCIsOnly* >> allRes_cisOnly.txt
# rm varCIsOnly*

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x/")
dat <- read.table("nullSims_8ped_fem/allRes_cisOnly.txt",header=FALSE,skip=1,as.is=TRUE,fill=TRUE)
head(dat)
colnames(dat)  <- c("Est","Lower95","Upper95","model","comp","sigma_auto","sigma_x")

table(dat$model,dat$comp)
dat <- dat[!is.element(dat$model,"Upper"),]
dat <- dat[!is.element(dat$model,""),]

table(dat$comp,dat$sigma_auto,dat$sigma_x)

# ok, so want boxplots of the estimate for each var component, stratified by model and component

library(ggplot2)
cib <- dat[dat$model=="both",]
ftable(cib$comp,cib$sigma_a,cib$sigma_x) # 9984 of each configuration

cib$sigma_a=paste0("sigma_a=",cib$sigma_auto)
cib$sigma_x=paste0("sigma_x=",cib$sigma_x)
cib$comp=ordered(cib$comp,levels=c("auto","x","resid"))

cib$Est <- as.numeric(cib$Est)

mns <- expand.grid(comp=c("auto","x","resid"),sigma_a=c("sigma_a=0.3","sigma_a=0.5"),
                   sigma_x=c("sigma_x=0.3","sigma_x=0.5"))
mns$comp <- ordered(mns$comp,levels=c("auto","x","resid"))
mns$y <- -0.1
mns$value <- NA
for(i in seq_len(nrow(mns))){
  mns$value[i] <- mean(cib$Est[cib$comp==mns$comp[i]&cib$sigma_a==mns$sigma_a[i]&cib$sigma_x==mns$sigma_x[i]])
}
mns$value <- format(mns$value,digits=3)

pdf("boxplots_varComp_bothEst_8ped_fem_10Kiters.pdf")
ggplot(cib,aes(x=comp,y=Est)) + geom_boxplot() + facet_grid(sigma_a~sigma_x) + theme_bw() +
  xlab("Variance Component") + ylab("Estimate") +
  ggtitle("Variance Component Estimates of 10K Iterations") +
  geom_text(data=mns, aes(x=comp, y=y, label=value), size=4)
dev.off()
##
cib <- dat[dat$model=="auto",]
ftable(cib$comp,cib$sigma_a,cib$sigma_x) # 9984 of each configuration

cib$sigma_a=paste0("sigma_a=",cib$sigma_auto)
cib$sigma_x=paste0("sigma_x=",cib$sigma_x)
cib$comp=ordered(cib$comp,levels=c("auto","resid"))

cib$Est <- as.numeric(cib$Est)

mns <- expand.grid(comp=c("auto","resid"),sigma_a=c("sigma_a=0.3","sigma_a=0.5"),
                   sigma_x=c("sigma_x=0.3","sigma_x=0.5"))
mns$comp <- ordered(mns$comp,levels=c("auto","resid"))
mns$y <- 0.33
mns$value <- NA
for(i in seq_len(nrow(mns))){
  mns$value[i] <- mean(cib$Est[cib$comp==mns$comp[i]&cib$sigma_a==mns$sigma_a[i]&cib$sigma_x==mns$sigma_x[i]])
}
mns$value <- format(mns$value,digits=3)

pdf("boxplots_varComp_autoEst_8ped_fem_10Kiters.pdf")
ggplot(cib,aes(x=comp,y=Est)) + geom_boxplot() + facet_grid(sigma_a~sigma_x) + theme_bw() +
  xlab("Variance Component") + ylab("Estimate") +
  ggtitle("Variance Component Estimates of 10K Iterations") +
  geom_text(data=mns, aes(x=comp, y=y, label=value), size=4)
dev.off()
##
cib <- dat[dat$model=="x",]
ftable(cib$comp,cib$sigma_a,cib$sigma_x) # 9984 of each configuration

cib$sigma_a=paste0("sigma_a=",cib$sigma_auto)
cib$sigma_x=paste0("sigma_x=",cib$sigma_x)
cib$comp=ordered(cib$comp,levels=c("x","resid"))

cib$Est <- as.numeric(cib$Est)

mns <- expand.grid(comp=c("x","resid"),sigma_a=c("sigma_a=0.3","sigma_a=0.5"),
                   sigma_x=c("sigma_x=0.3","sigma_x=0.5"))
mns$comp <- ordered(mns$comp,levels=c("x","resid"))
mns$y <- 0.2
mns$value <- NA
for(i in seq_len(nrow(mns))){
  mns$value[i] <- mean(cib$Est[cib$comp==mns$comp[i]&cib$sigma_a==mns$sigma_a[i]&cib$sigma_x==mns$sigma_x[i]])
}
mns$value <- format(mns$value,digits=3)

pdf("boxplots_varComp_xEst_8ped_fem_10Kiters.pdf")
ggplot(cib,aes(x=comp,y=Est)) + geom_boxplot() + facet_grid(sigma_a~sigma_x) + theme_bw() +
  xlab("Variance Component") + ylab("Estimate") +
  ggtitle("Variance Component Estimates of 10K Iterations") +
  geom_text(data=mns, aes(x=comp, y=y, label=value), size=4)
dev.off()

rm(list=ls())


#####
# 73. Make pdf for paper that includes all var comp boxplots

library(GWASTools)
library(ggplot2)

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x/")
dat <- read.table("nullSims_8ped_fem/allRes_cisOnly.txt",header=FALSE,skip=1,as.is=TRUE,fill=TRUE)
head(dat)
colnames(dat)  <- c("Est","Lower95","Upper95","model","comp","sigma_auto","sigma_x")

table(dat$model,dat$comp)
dat <- dat[!is.element(dat$model,"Upper"),]
dat <- dat[!is.element(dat$model,""),]

table(dat$comp,dat$sigma_auto,dat$sigma_x)

cis <- getobj("nullSims_8ped/varCompCIs_nullSims.RData")


mf_labeller <- function(var, value){
  value <- as.character(value)
  if (var=="sigma_x") {
    value[value=="sigma_x=0.3"] <- "sigma[X]^{2}==0.3"
    value[value=="sigma_x=0.5"]  <- "sigma[X]^{2}==0.5"
    value <- lapply(value, function(x) parse(text=x))
  }
  if (var=="sigma_a"){
    value[value=="sigma_a=0.3"] <- "sigma[A]^{2}==0.3"
    value[value=="sigma_a=0.5"] <- "sigma[A]^{2}==0.5"
    value <- lapply(value,function(x) parse(text=x))
  }
  return(value)
}


## male-centric, MLM-X results
cib <- cis[["cis_both"]]
head(cib)
# exclude the extreme values

cib <- cib[!is.element(cib$sigma_a,3)&!is.element(cib$sigma_x,3),]
cib$sigma_a <- paste0("sigma_a=",cib$sigma_a)
cib$sigma_x <- paste0("sigma_x=",cib$sigma_x)
cib$comp <- ordered(cib$comp,levels=c("auto","x","resid"))
ftable(cib$comp,cib$sigma_a,cib$sigma_x)
#                    sigma_x=0.3 sigma_x=0.5
# auto  sigma_a=0.3        10010       10010
#       sigma_a=0.5        10010       10010
# x     sigma_a=0.3        10010       10010
#       sigma_a=0.5        10010       10010
# resid sigma_a=0.3        10010       10010
#       sigma_a=0.5        10010       10010

# so 10K iterations of each composition
# add in mean for each component
mns <- expand.grid(comp=c("auto","x","resid"),sigma_a=c("sigma_a=0.3","sigma_a=0.5"),
                   sigma_x=c("sigma_x=0.3","sigma_x=0.5"))
mns$comp <- ordered(mns$comp,levels=c("auto","x","resid"))
mns$y <- -0.1
mns$value <- NA
for(i in seq_len(nrow(mns))){
  mns$value[i] <- mean(cib$Est[cib$comp==mns$comp[i]&cib$sigma_a==mns$sigma_a[i]&cib$sigma_x==mns$sigma_x[i]])
}

mns$value <- format(mns$value,digits=3)

p1 <- ggplot(cib,aes(x=comp,y=Est)) + geom_boxplot() + facet_grid(sigma_a~sigma_x, labeller=mf_labeller) + theme_bw() +
  xlab("Variance Component") + ylab("Estimate") +
  #ggtitle("Variance Component Estimates of 10K Iterations") +
  ggtitle('A') + theme(plot.title=element_text(hjust=0)) +
  geom_text(data=mns, aes(x=comp, y=y, label=value), size=4)

#mtext("A", side=3, line=0.75,adj=0,cex=1.3)

## male-centric, simple MLM results (ie. auto)
cib <- cis[["cis_auto"]]
head(cib)
# exclude the extreme values

cib <- cib[!is.element(cib$sigma_a,3)&!is.element(cib$sigma_x,3),]
cib$sigma_a <- paste0("sigma_a=",cib$sigma_a)
cib$sigma_x <- paste0("sigma_x=",cib$sigma_x)
cib$comp <- ordered(cib$comp,levels=c("auto","resid"))
ftable(cib$comp,cib$sigma_a,cib$sigma_x)
#                    sigma_x=0.3 sigma_x=0.5
# auto  sigma_a=0.3        10010       10010
#       sigma_a=0.5        10010       10010
# resid sigma_a=0.3        10010       10010
#       sigma_a=0.5        10010       10010
# so 10K iterations of each composition
# add in mean for each component

mns <- expand.grid(comp=c("auto","resid"),sigma_a=c("sigma_a=0.3","sigma_a=0.5"),
                   sigma_x=c("sigma_x=0.3","sigma_x=0.5"))
mns$comp <- ordered(mns$comp,levels=c("auto","resid"))
mns$y <- 0.33
mns$value <- NA
for(i in seq_len(nrow(mns))){
  mns$value[i] <- mean(cib$Est[cib$comp==mns$comp[i]&cib$sigma_a==mns$sigma_a[i]&cib$sigma_x==mns$sigma_x[i]])
}

mns$value <- format(mns$value,digits=3)

p2 <- ggplot(cib,aes(x=comp,y=Est)) + geom_boxplot() + facet_grid(sigma_a~sigma_x, labeller=mf_labeller) + theme_bw() +
  xlab("Variance Component") + ylab("Estimate") +
  #ggtitle("Variance Component Estimates of 10K Iterations") +
  ggtitle('B') + theme(plot.title=element_text(hjust=0)) +
  geom_text(data=mns, aes(x=comp, y=y, label=value), size=4)

#mtext("B", side=3, line=0.75,adj=0,cex=1.3)

## male-centric, g_x only results (ie. x)
cib <- cis[["cis_x"]]
head(cib)
# exclude the extreme values

cib <- cib[!is.element(cib$sigma_a,3)&!is.element(cib$sigma_x,3),]
cib$sigma_a <- paste0("sigma_a=",cib$sigma_a)
cib$sigma_x <- paste0("sigma_x=",cib$sigma_x)
cib$comp <- ordered(cib$comp,levels=c("x","resid"))
ftable(cib$comp,cib$sigma_a,cib$sigma_x)
#                    sigma_x=0.3 sigma_x=0.5
# x     sigma_a=0.3        10010       10010
#       sigma_a=0.5        10010       10010
# resid sigma_a=0.3        10010       10010
#       sigma_a=0.5        10010       10010
# so 10K iterations of each composition
# add in mean for each component

mns <- expand.grid(comp=c("x","resid"),sigma_a=c("sigma_a=0.3","sigma_a=0.5"),
                   sigma_x=c("sigma_x=0.3","sigma_x=0.5"))
mns$comp <- ordered(mns$comp,levels=c("x","resid"))
mns$y <- 0.2
mns$value <- NA
for(i in seq_len(nrow(mns))){
  mns$value[i] <- mean(cib$Est[cib$comp==mns$comp[i]&cib$sigma_a==mns$sigma_a[i]&cib$sigma_x==mns$sigma_x[i]])
}

mns$value <- format(mns$value,digits=3)

p3 <- ggplot(cib,aes(x=comp,y=Est)) + geom_boxplot() + facet_grid(sigma_a~sigma_x, labeller=mf_labeller) + theme_bw() +
  xlab("Variance Component") + ylab("Estimate") +
  #ggtitle("Variance Component Estimates of 10K Iterations") +
  ggtitle('C') + theme(plot.title=element_text(hjust=0)) +
  geom_text(data=mns, aes(x=comp, y=y, label=value), size=4)

#mtext("C", side=3, line=0.75,adj=0,cex=1.3)

## balanced (fem), mlm-x results (ie both)
cib <- dat[dat$model=="both",]
ftable(cib$comp,cib$sigma_a,cib$sigma_x) # 9984 of each configuration

cib$sigma_a=paste0("sigma_a=",cib$sigma_auto)
cib$sigma_x=paste0("sigma_x=",cib$sigma_x)
cib$comp=ordered(cib$comp,levels=c("auto","x","resid"))

cib$Est <- as.numeric(cib$Est)

mns <- expand.grid(comp=c("auto","x","resid"),sigma_a=c("sigma_a=0.3","sigma_a=0.5"),
                   sigma_x=c("sigma_x=0.3","sigma_x=0.5"))
mns$comp <- ordered(mns$comp,levels=c("auto","x","resid"))
mns$y <- -0.1
mns$value <- NA
for(i in seq_len(nrow(mns))){
  mns$value[i] <- mean(cib$Est[cib$comp==mns$comp[i]&cib$sigma_a==mns$sigma_a[i]&cib$sigma_x==mns$sigma_x[i]])
}
mns$value <- format(mns$value,digits=3)

p4 <- ggplot(cib,aes(x=comp,y=Est)) + geom_boxplot() + facet_grid(sigma_a~sigma_x, labeller=mf_labeller) + theme_bw() +
  xlab("Variance Component") + ylab("Estimate") +
  #ggtitle("Variance Component Estimates of 10K Iterations") +
  ggtitle('D') + theme(plot.title=element_text(hjust=0)) +
  geom_text(data=mns, aes(x=comp, y=y, label=value), size=4)

#mtext("D", side=3, line=0.75,adj=0,cex=1.3)

## balanced (fem), simple MLM results (ie. auto)

cib <- dat[dat$model=="auto",]
ftable(cib$comp,cib$sigma_a,cib$sigma_x) # 9984 of each configuration

cib$sigma_a=paste0("sigma_a=",cib$sigma_auto)
cib$sigma_x=paste0("sigma_x=",cib$sigma_x)
cib$comp=ordered(cib$comp,levels=c("auto","resid"))

cib$Est <- as.numeric(cib$Est)

mns <- expand.grid(comp=c("auto","resid"),sigma_a=c("sigma_a=0.3","sigma_a=0.5"),
                   sigma_x=c("sigma_x=0.3","sigma_x=0.5"))
mns$comp <- ordered(mns$comp,levels=c("auto","resid"))
mns$y <- 0.33
mns$value <- NA
for(i in seq_len(nrow(mns))){
  mns$value[i] <- mean(cib$Est[cib$comp==mns$comp[i]&cib$sigma_a==mns$sigma_a[i]&cib$sigma_x==mns$sigma_x[i]])
}
mns$value <- format(mns$value,digits=3)

p5 <- ggplot(cib,aes(x=comp,y=Est)) + geom_boxplot() + facet_grid(sigma_a~sigma_x, labeller=mf_labeller) + theme_bw() +
  xlab("Variance Component") + ylab("Estimate") +
  #ggtitle("Variance Component Estimates of 10K Iterations") +
  ggtitle('E') + theme(plot.title=element_text(hjust=0)) +
  geom_text(data=mns, aes(x=comp, y=y, label=value), size=4)

#mtext("E", side=3, line=0.75,adj=0,cex=1.3)

## balanced (fem), g_x only results (ie. x)

cib <- dat[dat$model=="x",]
ftable(cib$comp,cib$sigma_a,cib$sigma_x) # 9984 of each configuration

cib$sigma_a=paste0("sigma_a=",cib$sigma_auto)
cib$sigma_x=paste0("sigma_x=",cib$sigma_x)
cib$comp=ordered(cib$comp,levels=c("x","resid"))

cib$Est <- as.numeric(cib$Est)

mns <- expand.grid(comp=c("x","resid"),sigma_a=c("sigma_a=0.3","sigma_a=0.5"),
                   sigma_x=c("sigma_x=0.3","sigma_x=0.5"))
mns$comp <- ordered(mns$comp,levels=c("x","resid"))
mns$y <- 0.2
mns$value <- NA
for(i in seq_len(nrow(mns))){
  mns$value[i] <- mean(cib$Est[cib$comp==mns$comp[i]&cib$sigma_a==mns$sigma_a[i]&cib$sigma_x==mns$sigma_x[i]])
}
mns$value <- format(mns$value,digits=3)

p6 <- ggplot(cib,aes(x=comp,y=Est)) + geom_boxplot() + facet_grid(sigma_a~sigma_x, labeller=mf_labeller) +
  theme_bw() +
  xlab("Variance Component") + ylab("Estimate") +
  #ggtitle("Variance Component Estimates of 10K Iterations") +
  ggtitle('F') + theme(plot.title=element_text(hjust=0)) +
  geom_text(data=mns, aes(x=comp, y=y, label=value), size=4)

#mtext("F", side=3, line=0.75,adj=0,cex=1.3)

library(gridExtra); library(RGraphics)

pdf("boxplots_varComps_allModels_allPeds_10Kiters.pdf",width=14,height=10.5)
#par( mai=c(0.5, 0.65, 0.4, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 ,mfrow=c(2,3))
grid.arrange(p1,p2,p3,p4,p5,p6,ncol=3)
dev.off()

rm(list=ls())


#####
# 74. Power simulations, all configs

# call 10K sims for the fem pedigree, ie balanced pedigree
# qsub -q olga.q -p -50 -t 1-2500 -N powBalAuto batch_power_simulations.sh
# which calls power_simulations_v3.R: power sims for auto snp with fem/balanced pedigree
# stores results in powerSims_8ped_fem/autoSNP

# call 10K sims for the male-centric pedigree
# qsub -q olga.q -p -50 -t 1-2000 -N powMalAuto batch_power_simulations.sh
# which calls power_simulations_v4.R: power sims for auto snp with male-centric pedigree
# stores results in powerSims_8ped/autoSNP .txt files for the mmRes, can join these by unix


#####
# 75. Parse extreme sig type I error results for paper

library(tidyr)
library(readr)
library(dplyr)
setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")

dat <- get(load("nullSims_8ped_fem/mmRes_extremeSig_autoSNP_nullSims.RData"))
dat <- tbl_df(dat)

tyI_group <- dat %>%
  filter(model!="linear") %>%
  group_by(sigma_x,sigma_a,model)
summarize(tyI_group,n()) # 750K for each variable combo

summarize(tyI_group,tyI_err=sum(pval<1e-3)/n())
#   sigma_x sigma_a model      tyI_err
# 1     0.3     3.0  auto 0.0011013333
# 2     0.3     3.0  both 0.0010920000
# 3     0.3     3.0     x 0.0021640000
# 4     3.0     0.3  auto 0.0010373333
# 5     3.0     0.3  both 0.0009626667
# 6     3.0     0.3     x 0.0010386667

summarize(tyI_group,tyI_err=sum(pval<1e-04)/n())
#   sigma_x sigma_a model      tyI_err
# 1     0.3     3.0  auto 1.173333e-04
# 2     0.3     3.0  both 1.226667e-04
# 3     0.3     3.0     x 3.173333e-04
# 4     3.0     0.3  auto 8.400000e-05
# 5     3.0     0.3  both 9.333333e-05
# 6     3.0     0.3     x 9.733333e-05

summarize(tyI_group,tyI_err=sum(pval<1e-05)/n())
#   sigma_x sigma_a model      tyI_err
# 1     0.3     3.0  auto 1.200000e-05
# 2     0.3     3.0  both 1.200000e-05
# 3     0.3     3.0     x 3.200000e-05
# 4     3.0     0.3  auto 1.066667e-05
# 5     3.0     0.3  both 8.000000e-06
# 6     3.0     0.3     x 9.333333e-06


dat <- get(load("nullSims_8ped_fem/mmRes_extremeSig_nullSims.RData"))
dat <- tbl_df(dat)

tyI_group <- dat %>%
  filter(model!="linear") %>%
  group_by(sigma_x,sigma_a,model)
summarize(tyI_group,n()) # 750K for each variable combo

summarize(tyI_group,tyI_err=sum(pval<1e-3)/n())
#   sigma_x sigma_a model      tyI_err
# 1     0.3     3.0  auto 0.0012933333
# 2     0.3     3.0  both 0.0009826667
# 3     0.3     3.0     x 0.0010280000
# 4     3.0     0.3  auto 0.0046386667
# 5     3.0     0.3  both 0.0009906667
# 6     3.0     0.3     x 0.0009893333

summarize(tyI_group,tyI_err=sum(pval<1e-4)/n())
#   sigma_x sigma_a model      tyI_err
# 1     0.3     3.0  auto 0.0001613333
# 2     0.3     3.0  both 0.0001133333
# 3     0.3     3.0     x 0.0000960000
# 4     3.0     0.3  auto 0.0008746667
# 5     3.0     0.3  both 0.0001200000
# 6     3.0     0.3     x 0.0001213333

summarize(tyI_group,tyI_err=sum(pval<1e-5)/n())
#   sigma_x sigma_a model      tyI_err
# 1     0.3     3.0  auto 1.600000e-05
# 2     0.3     3.0  both 9.333333e-06
# 3     0.3     3.0     x 1.333333e-05
# 4     3.0     0.3  auto 1.600000e-04
# 5     3.0     0.3  both 1.066667e-05
# 6     3.0     0.3     x 1.066667e-05

rm(list=ls())


#####
# 76. Make tyI error plots for JSM poster (x + auto SNPs)

library(ggplot2)
library(dplyr); library(tidyr)
library(readr)

dat <- get(load("/projects/geneva/geneva_sata/caitlin/mlm_x/nullSims_8ped_autoSNP/mmRes_nullSims.RData"))
dim(dat) # 6000000      13
head(dat)

dat$model[dat$model=="both"] <- "MLM-X"
dat$model[dat$model=="auto"] <- "Simple MLM"
dat$model[dat$model=="x"] <- "X only"
dat$model[dat$model=="linear"] <- "Linear"

dat <- tbl_df(dat)
dat <- mutate(dat,causal=NULL)
dat <- group_by(dat,model,sigma_a,sigma_x)

alpha <- 5e-05
# 2x2 of sigma_a, sigma_x values
# 4 entries along the x axis for the 4 models considered
# change the x axis values

summarize(dat,n()) # 375K of each combo
tyI <- summarize(dat,est=sum(pval<alpha)/n(),
                 stdErr=sqrt(alpha*(1-alpha)/n()),
                 lower=est-1.96*stdErr,
                 upper=est+1.96*stdErr)
tyI$sigma_a <- paste0("sigma_a=",tyI$sigma_a)
tyI$sigma_x <- paste0("sigma_x=",tyI$sigma_x)

tyI$model <- ordered(tyI$model,levels=c("MLM-X","X only","Simple MLM","Linear"))

tyIerrSm <- tyI
pdf("/projects/geneva/geneva_sata/caitlin/mlm_x/typeIErr_8ped_autoSNP_5e05.pdf")
ggplot(tyIerrSm,aes(x=model,y=est)) +geom_point()+ geom_pointrange(aes(ymin=lower,ymax=upper))+
  theme_bw() + geom_hline(aes(yintercept=5e-05,color="gray"))+
  facet_grid(sigma_x~sigma_a) + xlab("Model, Testing Autosomal SNP") +
  theme(axis.text=element_text(size=10),axis.title=element_text(size=18),title=element_text(size=20),
        legend.text=element_text(size=16),strip.text = element_text(size=16)) +
  ylab("Type I Error Rate") + ggtitle(expression(paste("Type I Error Rate, ", alpha, "=5e-05")))
dev.off()


dat <- read_delim("/projects/geneva/geneva_sata/caitlin/mlm_x/nullSims_8ped/allmmRes.txt",
                  delim=" ")

dat$model[dat$model=="both"] <- "MLM-X"
dat$model[dat$model=="auto"] <- "Simple MLM"
dat$model[dat$model=="x"] <- "X only"
dat$model[dat$model=="linear"] <- "Linear"

dat <- group_by(dat,model,sigma_a,sigma_x)
dat <- filter(dat,!is.na(sigma_a))

alpha <- 5e-05
# 2x2 of sigma_a, sigma_x values
# 4 entries along the x axis for the 4 models considered
# change the x axis values

summarize(dat,n()) # 375K of each combo
tyI <- summarize(dat,est=sum(pval<alpha)/n(),
                 stdErr=sqrt(alpha*(1-alpha)/n()),
                 lower=est-1.96*stdErr,
                 upper=est+1.96*stdErr)
tyI$sigma_a <- paste0("sigma_a=",tyI$sigma_a)
tyI$sigma_x <- paste0("sigma_x=",tyI$sigma_x)

tyI$model <- ordered(tyI$model,levels=c("MLM-X","X only","Simple MLM","Linear"))

tyIerrSm <- tyI

pdf("/projects/geneva/geneva_sata/caitlin/mlm_x/typeIErr_8ped_5e05.pdf")
ggplot(tyIerrSm,aes(x=model,y=est)) +geom_point()+ geom_pointrange(aes(ymin=lower,ymax=upper))+
  theme_bw() + geom_hline(aes(yintercept=5e-05,color="gray"))+
  facet_grid(sigma_x~sigma_a) + xlab("Model, Testing X Chr SNP") +
  theme(axis.text=element_text(size=10),axis.title=element_text(size=18),title=element_text(size=20),
        legend.text=element_text(size=16),strip.text = element_text(size=16)) +
  ylab("Type I Error Rate") + ggtitle(expression(paste("Type I Error Rate, ", alpha, "=5e-05")))
dev.off()

rm(list=ls())


#####
# 77. Make power plots for JSM poster (x + auto SNPs)

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x/")
source("powerGraph.R")
library(ggplot2); library(reshape)
library(readr); library(tidyr)
library(dplyr)

dat <- read.table("powerSims_8ped/autoSNP/allmmRes.txt",header=TRUE,as.is=TRUE)
dim(dat) # 84999 13

headRw <- which(dat$snpID=="snpID")
dat <- dat[-headRw,]
dim(dat) # 80000 13

dat <- tbl_df(dat)
dat$causal[dat$causal==0.05] <- TRUE
dat$causal[dat$causal==0] <- FALSE
toPl <- dat %>%
  filter(sigma_x==3&sigma_a==0.3) %>%
  mutate(beta1=0.05) %>%
  mutate(causal=as.logical(causal))

sa3sx3=powerGraph(toPl,b=0.05,fn="powerGraph_sA03sX03.pdf",xlims=c(0,0.05),values=TRUE)
sa3sx3$sigma_a=paste0("sigma_a=",0.3)
sa3sx3$sigma_x=paste0("sigma_x=",3)

toPl <- dat %>%
  filter(sigma_x==0.3&sigma_a==3) %>%
  mutate(beta1=0.05) %>%
  mutate(causal=as.logical(causal))
sx5sa3=powerGraph(toPl,b=0.05,fn="powerGraph_sA05sX03.pdf",xlims=c(0,0.05),values=TRUE)
sx5sa3$sigma_a=paste0("sigma_a=",3)
sx5sa3$sigma_x=paste0("sigma_x=",0.3)

tot <- rbind(sa3sx3,sx5sa3)
tot$x_fp_rate <- tot$x_fp/tot$xiter
tot$x_tp_rate <- tot$x_tp/tot$xiter
tot$auto_fp_rate <- tot$auto_fp/tot$autoiter
tot$auto_tp_rate <- tot$auto_tp/tot$autoiter
tot$both_fp_rate <- tot$both_fp/tot$bothiter
tot$both_tp_rate <- tot$both_tp/tot$bothiter

both <- tot[,c("sigma_a","sigma_x","both_fp_rate","both_tp_rate")]
both$model <- "MLM-X"
both$fp <- both$both_fp_rate
both$tp <- both$both_tp_rate

auto <- tot[,c("sigma_a","sigma_x","auto_fp_rate","auto_tp_rate")]
auto$model <- "Simple MLM"
auto$fp <- auto$auto_fp_rate
auto$tp <- auto$auto_tp_rate

x <- tot[,c("sigma_a","sigma_x","x_fp_rate","x_tp_rate")]
x$model <- "X only"
x$fp <- x$x_fp_rate
x$tp <- x$x_tp_rate

toPl <- rbind(x[,c("sigma_x","sigma_a","model","tp","fp")],
              auto[,c("sigma_x","sigma_a","model","tp","fp")],
              both[,c("sigma_x","sigma_a","model","tp","fp")])

toPl <- toPl[toPl$fp<=0.05,]
toPlauto <- toPl

pdf("power_beta05_8ped_autoSNP_extrememSigma_10kiters.pdf",width=11,height=11)
ggplot(toPl,aes(x=fp,y=tp,color=model)) + geom_line(size=1.5,aes(linetype=model)) + facet_wrap(sigma_a~sigma_x) + theme_bw() +
  ylab("True Positive Rate") + xlab("False Positive Rate") +ggtitle(expression(paste(beta,"=0.05, Testing Autosomal SNP")))+
  theme(axis.text=element_text(size=16),axis.title=element_text(size=18),title=element_text(size=20),
        legend.text=element_text(size=16),strip.text = element_text(size=16))
dev.off()

# add in results testing X chr SNP
dat <- read_delim("/projects/geneva/geneva_sata/caitlin/mlm_x/powerSims_8ped/all_mmRes.txt",
                  delim=" ")

dat$causal[dat$causal=="causal"] <- NA
dat$causal[dat$causal==0.05] <- TRUE
dat$causal[dat$causal==0] <- FALSE
toPl <- dat %>%
  filter(sigma_x=="3"&sigma_a=="0.3") %>%
  mutate(beta1=0.05) %>%
  mutate(causal=as.logical(causal))

sa3sx3=powerGraph(toPl,b=0.05,fn="powerGraph_sA03sX03.pdf",xlims=c(0,0.05),values=TRUE)
sa3sx3$sigma_a=paste0("sigma_a=",0.3)
sa3sx3$sigma_x=paste0("sigma_x=",3)

toPl <- dat %>%
  filter(sigma_x==0.3&sigma_a==3) %>%
  mutate(beta1=0.05) %>%
  mutate(causal=as.logical(causal))
sx5sa3=powerGraph(toPl,b=0.05,fn="powerGraph_sA05sX03.pdf",xlims=c(0,0.05),values=TRUE)
sx5sa3$sigma_a=paste0("sigma_a=",3)
sx5sa3$sigma_x=paste0("sigma_x=",0.3)

tot <- rbind(sa3sx3,sx5sa3)
tot$x_fp_rate <- tot$x_fp/tot$xiter
tot$x_tp_rate <- tot$x_tp/tot$xiter
tot$auto_fp_rate <- tot$auto_fp/tot$autoiter
tot$auto_tp_rate <- tot$auto_tp/tot$autoiter
tot$both_fp_rate <- tot$both_fp/tot$bothiter
tot$both_tp_rate <- tot$both_tp/tot$bothiter

both <- tot[,c("sigma_a","sigma_x","both_fp_rate","both_tp_rate")]
both$model <- "MLM-X"
both$fp <- both$both_fp_rate
both$tp <- both$both_tp_rate

auto <- tot[,c("sigma_a","sigma_x","auto_fp_rate","auto_tp_rate")]
auto$model <- "Simple MLM"
auto$fp <- auto$auto_fp_rate
auto$tp <- auto$auto_tp_rate

x <- tot[,c("sigma_a","sigma_x","x_fp_rate","x_tp_rate")]
x$model <- "X only"
x$fp <- x$x_fp_rate
x$tp <- x$x_tp_rate

toPl <- rbind(x[,c("sigma_x","sigma_a","model","tp","fp")],
              auto[,c("sigma_x","sigma_a","model","tp","fp")],
              both[,c("sigma_x","sigma_a","model","tp","fp")])

toPl <- toPl[toPl$fp<=0.05,]

pdf("power_beta05_8ped_extrememSigma_10kiters.pdf",width=11,height=11)
ggplot(toPl,aes(x=fp,y=tp,color=model)) + geom_line(size=1.5,aes(linetype=model)) + facet_wrap(sigma_a~sigma_x) + theme_bw() +
  ylab("True Positive Rate") + xlab("False Positive Rate") +ggtitle(expression(paste(beta,"=0.05, Testing X Chromosome SNP")))+
  theme(axis.text=element_text(size=16),axis.title=element_text(size=18),title=element_text(size=20),
        legend.text=element_text(size=16),strip.text = element_text(size=16))
dev.off()

toPlauto$chr <- "Autosomal SNP"
toPl$chr <- "X Chromosome SNP"
toPlBoth <- rbind(toPl,toPlauto)

pdf("power_beta05_8ped_bothSNP_extrememSigma_10kiters.pdf",width=11,height=11)
ggplot(toPlBoth,aes(x=fp,y=tp,color=model)) + 
  geom_line(size=1.5,aes(linetype=model)) + 
  facet_grid(chr~sigma_a + sigma_x) + theme_bw() + 
  ylab("True Positive Rate") + xlab("False Positive Rate") +
  ggtitle(expression(paste(beta,"=0.05")))+
  theme(axis.text=element_text(size=16),axis.title=element_text(size=18),title=element_text(size=20),
        legend.text=element_text(size=16),strip.text = element_text(size=16))
dev.off()




rm(list=ls())

