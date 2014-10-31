
# script to perform MLM simulations

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")

library(MASS)
library(corpcor)
source("sim_phenotype.R")
source("estVarComp.R")
source("allele_drop_functions.R")
source("assocTestMixedModel_v7.R")


oneSim <- function(kinFull,kinFullX,sigmaA,sigmaX,sigmaE,herit,p,i,fn){
  
  n <- 1600 # number of samples

  TEMP=c(2,1,2,1,1,1,2,2,1,2,1,2,1,2,1,2)
  SEX <- rep(NA,length(TEMP))
  SEX[TEMP==1]="M"
  SEX[TEMP==2]="F"
  
  sex <- rep(SEX,(n/16))
  
  # simulate genotypes now; only depends on the value of p
  otherp <- c(0.01,0.05,0.1,0.2,0.25)
  aFreqs <- rep(otherp,each=700)
  aFreqs[1:5] <- p

  geno <- matrix(NA,nrow=3500,ncol=2100)
  for(i in 1:100){
      geno[,(i*16-15):(i*16)] <- 2*Family_alleles_NmarkerX(aFreqs,3500,SEX)
    }

  for(i in 1:250){
      geno[,(1600+i)] <- 2*INDIV_alleles_Nmarker(aFreqs,3500)
        geno[,(1600+250+i)] <- 2*INDIV_alleles_NmarkerX(aFreqs,3500)
    } # cols are samples, rows are SNPs

  
  
  #genoXtmp <- 2*Family_alleles_NmarkerX(aFreqs,100*500,sex) # 100 peds * 500 SNPs
  #genoX <- t(genoXtmp)
  #genoF <- t(matrix(genoX,nrow=n))

  # add on 500 unrelated samples
#  aFreqs <- rep(otherp,each=100)
#  genoUnr <- matrix(0,nrow=500,ncol=500)
#  for(i in 1:250){
#    genoUnr[,i] <- 2*INDIV_alleles_Nmarker(aFreqs,500) 
#    genoUnr[,(i+250-1)] <- 2*INDIV_alleles_NmarkerX(aFreqs,500)
#  } # cols are samples, rows are SNPs

#  geno <- cbind(genoF,genoUnr)
  
  causalSNP <- as.vector(geno[1,])
  ##

  n <- 2100
  sex <- c(sex,rep("F",250),rep("M",250))

  colnames(kinFull) <- 1:n
  colnames(kinFullX) <- 1:n
  
  # loop through the different values for each of the sigmaE, sigmaA, sigmaX, herit
  allRes <- vector("list",length(sigmaE))
  for(i in 1:length(sigmaE)){
    totalSigmaScalar <- sigmaE[i]+sigmaA[i]+sigmaX[i]
    totalSigma <- sigmaE[i]+sigmaA[i]*kinFull+sigmaX[i]*kinFullX
  
    # calculate effect size based on herit, allele freq
    beta1 <- getEffSize(p,totalSigmaScalar,herit=herit[i])
  
    # simulate the phenotype based on these metrics
    pheno <- simulatePhenotype(kinFull,kinFullX,sigmaA[i],sigmaX[i],sigmaE[i],beta1,causalSNP,seed=i)
  
    df <- data.frame(scanID=1:n,sex=sex,pheno=pheno)
    scan <- ScanAnnotationDataFrame(df)
  
    genoMt <- MatrixGenotypeReader(genotype=geno,snpID=as.integer(1:nrow(geno)),
                                   chromosome=as.integer(rep(23,nrow(geno))),
                                   position=as.integer(1:nrow(geno)), scanID=scan$scanID)
    genoData <- GenotypeData(genoMt,scanAnnot=scan)
  
    ## now estimate the variance components
    #colnames(kinFull) <- scan$scanID
    #colnames(kinFullX) <- scan$scanID
  
    covMatList <- list(kinFull,kinFullX)
    names(covMatList) <- c("kinshipAuto","kinshipX")
    varCompInclX <- estVarComp(scan,covMatList=covMatList,"pheno")
  
    # estimate without adjusting for x chr
    covMatList <- list(kinFull)
    names(covMatList) <- c("kinshipAuto")
    varCompExclX <- estVarComp(scan,covMatList=covMatList,"pheno")
  
    # estimate with only adj for x chr
    covMatList <- list(kinFullX)
    names(covMatList) <- c("kinshipX")
    varComp <- estVarComp(scan,covMatList=covMatList,"pheno")

    
    # call MLM with these results
    cholSig <- varCompInclX[["cholSigmaInv"]]
    mmResXAuto <- assocTestMixedModel(genoData,snpStart=1,snpEnd=nrow(geno),cholSig,outcome="pheno")
  
    # estimate without adjusting for x chr
    cholSig <- varCompExclX[["cholSigmaInv"]]
    mmResAuto <- assocTestMixedModel(genoData,snpStart=1,snpEnd=nrow(geno),cholSig,outcome="pheno")

    # estimate only adj for x chr
    cholSig <- varComp[["cholSigmaInv"]]
    mmResX <- assocTestMixedModel(genoData,snpStart=1,snpEnd=nrow(geno),cholSig,outcome="pheno")
    
    #fname <- paste(fn,"_",i,".RData",sep="")
    metrics <- data.frame(beta=beta1,herit=herit[i],p=p,sigmaA=sigmaA[i],sigmaX=sigmaX[i],sigmaE=sigmaE[i])
    varCompsX <- data.frame(varComp$varComp)
    varCompAutoX <- data.frame(varCompInclX$varComp)
    varCompAuto <- data.frame(varCompExclX$varComp)
    simRes <- list(mmResX,mmResAuto,mmResXAuto,metrics,varCompsX,varCompAutoX,varCompAuto)
    names(simRes) <- c("mmRes_onlyAdjForX","mmRes_autoOnly","mmRes_both_Xauto","metricsSet","varComps_est_inclXonly","varComps_est_autoX","varComps_justAuto")
  
    #save(simRes,file=fname)
    allRes[[i]] <- simRes
    print(paste("********** ended iteration",i," ************"))
  }
  fname <- paste(fn,".RData",sep="")
  save(allRes,file=fname)
  return(NULL)
}


## do a quicker version of the function that does 1000 simulations, saving only the pvalues
