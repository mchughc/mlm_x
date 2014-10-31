
# script to perform MLM simulations

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")

library(MASS)
library(corpcor)
source("sim_phenotype.R")
source("estVarComp.R")
source("allele_drop_functions.R")
source("assocTestMixedModel_v7.R")


oneSim <- function(kinFull,kinFullX,sigmaA,sigmaX,sigmaE,herit,p,fn){
  
  n <- nrow(kinFull) # number of samples

  TEMP=c(2,1,2,1,1,1,2,2,1,2,1,2,1,2,1,2)
  SEX <- rep(NA,length(TEMP))
  SEX[TEMP==1]="M"
  SEX[TEMP==2]="F"
  
  sex <- rep(SEX,(n/16))

  colnames(kinFull) <- 1:n
  colnames(kinFullX) <- 1:n

  
  # simulate genotypes now; only depends on the value of p
  otherp <- c(0.2,0.01)
  aFreqs <- rep(otherp,each=100*500)
  aFreqs[1:1000] <- p
  
  genoXtmp <- 2*Family_alleles_NmarkerX(aFreqs,100*2,sex) # 100 peds * 2 SNPs
  genoX <- matrix(genoXtmp,nrow=n)
  geno <- t(genoX)
  
  totalSigmaScalar <- sigmaE+sigmaA+sigmaX
  
  # loop through the different values for each of the sigmaE, sigmaA, sigmaX, herit
  for(i in 1:length(herit)){
    
    totalSigma <- sigmaE[i]+sigmaA[i]*kinFull+sigmaX[i]*kinFullX
  
    # calculate effect size based on herit, allele freq
    beta1 <- getEffSize(p,totalSigmaScalar[i],herit=herit[i])
  
    # simulate the phenotype based on these metrics
    
      seed <- ceiling(runif(1,min=0,max=10000))
      print(seed)
      causalSNP <- as.vector(genoX[,1])
      pheno <- simulatePhenotype(kinFull,kinFullX,sigmaA[i],sigmaX[i],sigmaE[i],beta1,causalSNP,seed=seed)
    
      df <- data.frame(scanID=1:n,sex=sex,pheno=pheno,family=rep(1:100,each=16))
      scan <- ScanAnnotationDataFrame(df)
  
      genoMt <- MatrixGenotypeReader(genotype=geno,snpID=as.integer(1:nrow(geno)),
                                     chromosome=as.integer(rep(23,nrow(geno))),
                                     position=as.integer(1:nrow(geno)), scanID=scan$scanID)
      genoData <- GenotypeData(genoMt,scanAnnot=scan)
  
      ## now estimate the variance components
#      covMatList <- list(kinFull,kinFullX)
#      names(covMatList) <- c("kinshipAuto","kinshipX")
#      varCompInclX <- estVarComp(scan,covMatList=covMatList,"pheno")
  
      # estimate without adjusting for x chr
#      covMatList <- list(kinFull)
#      names(covMatList) <- c("kinshipAuto")
#      varCompExclX <- estVarComp(scan,covMatList=covMatList,"pheno")

      # estimate with adj for x chr only
    covMatList <- list(kinFullX)
    names(covMatList) <- c("kinshipX")
    varComp <- estVarComp(scan,covMatList=covMatList,"pheno")
  
      # call MLM with these results
#      cholSig <- varCompInclX[["cholSigmaInv"]]
#      mmResX <- assocTestMixedModel(genoData,snpStart=1,snpEnd=nrow(geno),cholSig,outcome="pheno")
  
      # estimate without adjusting for x chr
#      cholSig <- varCompExclX[["cholSigmaInv"]]
#      mmRes <- assocTestMixedModel(genoData,snpStart=1,snpEnd=nrow(geno),cholSig,outcome="pheno")

    cholSig <- varComp[["cholSigmaInv"]]
    mmRes <- assocTestMixedModel(genoData,snpStart=1,snpEnd=nrow(geno),cholSig,outcome="pheno")
    
      fname <- paste(fn,"_",i,".RData",sep="")
      metrics <- data.frame(beta=beta1,herit=herit[i],p=p,sigmaA=sigmaA[i],sigmaX=sigmaX[i],sigmaE=sigmaE[i])
      varComps <- data.frame(varComp$varComp)
#      varCompsExclX <- data.frame(varCompExclX$varComp)
      simRes <- list(mmRes,metrics,varComps)
      names(simRes) <- c("mmRes_adjForXonly","metricsSet","varComps_est_inclXonly")
  
      save(simRes,file=fname)
    
      print(paste("********** ended iteration",i," ************"))
  }
  return(NULL)
}


## do a quicker version of the function that does 1000 simulations, saving only the pvalues
