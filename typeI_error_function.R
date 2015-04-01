
## function to perform var comp estimation and assoc testing
doAnalysis <- function(input){
  
  n <- 8000

  kinX <- get(load("500Peds_xKinship.RData"))
  kinAuto <- get(load("500Peds_autoKinship.RData"))
  
  kinAuto <- matrix(kinAuto,nrow=n,ncol=n)
  kinX <- matrix(kinX,nrow=n,ncol=n)

  genoData <- input[[1]]
  scan <- input[[2]]
  snpStart <- input[[3]]
  snpEnd <- input[[4]]
  newBeta <- input[[5]]
  model <- input[[6]]
  
  colnames(kinX) <- getScanID(genoData)
  colnames(kinAuto) <- getScanID(genoData)

  covMatListBoth <- list(kinX,kinAuto) # this is the X chr kinship matrix for all samples
  names(covMatListBoth) <- c("kinshipX","kinshipAutos")

  if(model=="both"){
    startBoth <- c(newBeta^2*4*0.2*(1-0.2)+0.8,0.3,1)
    varCompBoth <- estVarComp(scan,covMatList=covMatListBoth,"pheno",start=startBoth,AIREML.tol=1e-4)
    varBothci <- estVarCompCI(varCompBoth,prop=FALSE)

    cholSig <- varCompBoth[["cholSigmaInv"]]
    mmResBoth <- assocTestMixedModel(genoData,snpStart=snpStart,snpEnd=snpEnd,cholSig,outcome="pheno")
    mmResBoth$causal <- c(TRUE,FALSE)

    mmRes <- mmResBoth
    varCi <- varBothci
  }
  if(model=="auto"){
    startAuto <- c(newBeta^2*4*0.2*(1-0.2)+0.8 + 0.3, 1)
    varCompA <- estVarComp(scan,covMatList=list(kinAuto),"pheno",start=startAuto,AIREML.tol=1e-4)
    varAutoci <- estVarCompCI(varCompA,prop=FALSE)

    cholSig <- varCompA[["cholSigmaInv"]]
    mmResA <- assocTestMixedModel(genoData,snpStart=snpStart,snpEnd=snpEnd,cholSig,outcome="pheno")
    mmResA$causal <- c(TRUE,FALSE)

    mmRes <- mmResA
    varCi <- varAutoci
  }
  if(model=="x"){
    startX <- c(newBeta^2*4*0.2*(1-0.2)+0.8, 1)
    varCompX <- estVarComp(scan,covMatList=list(kinX),"pheno",start=startX,AIREML.tol=1e-4)
    varXci <- estVarCompCI(varCompX,prop=FALSE)

    cholSig <- varCompX[["cholSigmaInv"]]
    mmResX <- assocTestMixedModel(genoData,snpStart=snpStart,snpEnd=snpEnd,cholSig,outcome="pheno")
    mmResX$causal <- c(TRUE,FALSE)

    mmRes <- mmResX
    varCi <- varXci
  }
  
  res <- list(newBeta,model,varCi,mmRes)
  names(res) <- c("beta1","model","var_est_ci","mm_res")
  
  return(res)
}

# v6: tried to optimize the code for 1M runs quickly!

# v5: added start value for the var comp estimation methods
## NOTE: depends on values of sigma_x and sigma_A: here, sigma_x=0.8 and sigma_a=0.3

# v4:
# does for 5 SNPs per iteration rather than just one

## v3: reads in the random effects
# simulates a new causal SNP
# calculates the new phenotype with beta value and snp
# calls var comp and assoc test

powerCalc <- function(i,newBeta,outputFn){

#  source("/projects/geneva/geneva_sata/caitlin/mlm_x/allele_drop_functions.R")
source("allele_drop_functions.R")  
##
# load in all the supporting docs

  thisCol <- as.integer(i)
  newBeta <- as.numeric(newBeta)
  
  n <- 8000
  
  TEMP=c(2,1,2,1,1,1,2,2,1,2,1,2,1,2,1,2)
  SEX <- rep(NA,length(TEMP))
  SEX[TEMP==1]="M"
  SEX[TEMP==2]="F"

  sex <- rep(SEX,(n/16))

#  kinX <- get(load("500Peds_xKinship.RData"))
#  kinAuto <- get(load("500Peds_autoKinship.RData"))
  
#  kinAuto <- matrix(kinAuto,nrow=n,ncol=n)
#  kinX <- matrix(kinX,nrow=n,ncol=n)

  randomEffects <- get(load("/Users/c8linmch/Dropbox/randomEffects.RData"))

  genoMat <- matrix(NA,nrow=10,ncol=n)
  colsToGeno <- split(1:n,rep(1:n,each=16,length=n))
  geno <- do.call(cbind,lapply(colsToGeno,function(x){genoMat[,x]=Family_alleles_NmarkerX(rep(0.2,10),10,SEX)}))
  
  #simGenos <- vector("list",500)
  #for(i in 1:500){
  #  simGenos[[i]] <- Family_alleles_NmarkerX(rep(0.2,10),10,SEX) # simulate 10 SNPs with MAF 0.2
  #}

  #for(i in 1:length(simGenos)){
  #  geno[,(i*16-15):(i*16)] <- simGenos[[i]]
  #}

  
  ##
  # first make causal SNP geno 1
  causalSNP <- geno[1,]
  recov1 <- randomEffects[,thisCol]
  pheno2 <- recov1+newBeta*causalSNP

  # now take the new phenotypes, est var comp and get the assoc test results
  scan1 <- data.frame(scanID=1:n,sex=sex,pheno=pheno2)
  scan1 <- ScanAnnotationDataFrame(scan1)

  genoMt <- MatrixGenotypeReader(genotype=geno,snpID=as.integer(1:nrow(geno)),
                                 chromosome=as.integer(rep(23,nrow(geno))),
                                 position=as.integer(1:nrow(geno)), scanID=scan1$scanID)
  
  genoData1 <- GenotypeData(genoMt,scanAnnot=scan1)

  ##
  # now make causal SNP geno 3
  causalSNP <- geno[3,]
  pheno2 <- recov1+newBeta*causalSNP

  # now take the new phenotypes, est var comp and get the assoc test results
  scan2 <- data.frame(scanID=1:n,sex=sex,pheno=pheno2)
  scan2 <- ScanAnnotationDataFrame(scan2)
  genoData2 <- GenotypeData(genoMt,scanAnnot=scan2)
    
  ##
  # now make causal SNP geno 5
  causalSNP <- geno[5,]
  pheno2 <- recov1+newBeta*causalSNP

    # now take the new phenotypes, est var comp and get the assoc test results
  scan3 <- data.frame(scanID=1:n,sex=sex,pheno=pheno2)
  scan3 <- ScanAnnotationDataFrame(scan3)
  genoData3 <- GenotypeData(genoMt,scanAnnot=scan3)

  ##
  # now make causal SNP geno 7
  causalSNP <- geno[7,]
  pheno2 <- recov1+newBeta*causalSNP

  # now take the new phenotypes, est var comp and get the assoc test results
  scan4 <- data.frame(scanID=1:n,sex=sex,pheno=pheno2)
  scan4 <- ScanAnnotationDataFrame(scan4)
  genoData4 <- GenotypeData(genoMt,scanAnnot=scan4)
  
  ##
  # now make causal SNP geno 9
  causalSNP <- geno[9,]
  pheno2 <- recov1+newBeta*causalSNP

  # now take the new phenotypes, est var comp and get the assoc test results
  scan5 <- data.frame(scanID=1:n,sex=sex,pheno=pheno2)
  scan5 <- ScanAnnotationDataFrame(scan5)
  genoData5 <- GenotypeData(genoMt,scanAnnot=scan5)


  geno1 <- list(genoData1,scan1,1,2,newBeta,"both")
  geno2 <- list(genoData2,scan2,3,4,newBeta,"both")
  geno3 <- list(genoData3,scan3,5,6,newBeta,"both")
  geno4 <- list(genoData4,scan4,7,8,newBeta,"both")
  geno5 <- list(genoData5,scan5,9,10,newBeta,"both")

  #### do the function calls in parallel
#  numCores <- detectCores()
#  inputs <- list(geno1,geno2,geno3,geno4,geno5)
#  resBoth <- mclapply(inputs,doAnalysis,mc.cores=numCores)

#  resBoth <- bplapply(inputs,doAnalysis)
  resBoth <- vector("list",5)
  resBoth[[1]] <- doAnalysis(geno1)
  resBoth[[2]] <- doAnalysis(geno2)
  resBoth[[3]] <- doAnalysis(geno3)
  resBoth[[4]] <- doAnalysis(geno4)
  resBoth[[5]] <- doAnalysis(geno5)
  
  geno1 <- list(genoData1,scan1,1,2,newBeta,"auto")
  geno2 <- list(genoData2,scan2,3,4,newBeta,"auto")
  geno3 <- list(genoData3,scan3,5,6,newBeta,"auto")
  geno4 <- list(genoData4,scan4,7,8,newBeta,"auto")
  geno5 <- list(genoData5,scan5,9,10,newBeta,"auto")
  
#  numCores <- detectCores()
#  inputs <- list(geno1,geno2,geno3,geno4,geno5)
#  resAuto <- mclapply(inputs,doAnalysis,mc.cores=numCores)

#  resAuto <- bplapply(inputs,doAnalysis)
  resAuto <- vector("list",5)
  resAuto[[1]] <- doAnalysis(geno1)
  resAuto[[2]] <- doAnalysis(geno2)
  resAuto[[3]] <- doAnalysis(geno3)
  resAuto[[4]] <- doAnalysis(geno4)
  resAuto[[5]] <- doAnalysis(geno5)
  
  geno1 <- list(genoData1,scan1,1,2,newBeta,"x")
  geno2 <- list(genoData2,scan2,3,4,newBeta,"x")
  geno3 <- list(genoData3,scan3,5,6,newBeta,"x")
  geno4 <- list(genoData4,scan4,7,8,newBeta,"x")
  geno5 <- list(genoData5,scan5,9,10,newBeta,"x")

#  numCores <- detectCores()
  #inputs <- list(geno1,geno2,geno3,geno4,geno5)
#  resX <- mclapply(inputs,doAnalysis,mc.cores=numCores)

  #resX <- bplapply(inputs,doAnalysis)

  resX <- vector("list",5)
  resX[[1]] <- doAnalysis(geno1)
  resX[[2]] <- doAnalysis(geno2)
  resX[[3]] <- doAnalysis(geno3)
  resX[[4]] <- doAnalysis(geno4)
  resX[[5]] <- doAnalysis(geno5)

  res <- list(resBoth,resAuto,resX)
  names(res) <- c("both","auto","x")
  
  fn <- paste("/projects/geneva/geneva_sata/caitlin/mlm_x/",outputFn,"/mmRes_SNP",thisCol,"_beta",newBeta,".RData",sep="")
  save(res,file=fn)

  return(NULL)
}

#####

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")

#library(BiocParallel)
library(MASS); library(parallel)
library(GWASTools)
library(corpcor); library(doSNOW)
source("sim_phenotype.R")
source("estVarComp.R")
source("allele_drop_functions.R")
source("assocTestMixedModel_v7.R")

args <- commandArgs(TRUE)
#args <- c(1,0.05,"tmp")

##

n <- 8000

TEMP=c(2,1,2,1,1,1,2,2,1,2,1,2,1,2,1,2)
SEX <- rep(NA,length(TEMP))
SEX[TEMP==1]="M"
SEX[TEMP==2]="F"

sex <- rep(SEX,(n/16))

## simulate 1million genotypes first
geno <- Family_alleles_NmarkerX(rep(0.2,(1e6+1)),(1e6+1),SEX)

# make the first SNP the causal genotype
randomEffects <- get(load("/Users/c8linmch/Dropbox/randomEffects.RData"))

randomEffects <- get(load("phenotypes_genotypes_1million_sigmaX05_sigmaA05.RData"))

powerCalc(args[1],args[2],args[3])

q("no")
