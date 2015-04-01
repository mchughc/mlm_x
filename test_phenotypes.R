
# simulate 1million + 1 genotypes
# make phenotype from one SNP, with varying sigma values
# then test all other genotypes against the phenotype

test_phenotypes <- function(numSNPs,sigmaA, sigmaX){
  
  n <- 8000
  
  kinX <- get(load("500Peds_xKinship.RData"))
  kinAuto <- get(load("500Peds_autoKinship.RData"))
  
  kinAuto <- matrix(kinAuto,nrow=n,ncol=n)
  kinX <- matrix(kinX,nrow=n,ncol=n)
  
  sigmaXMat <- mvrnorm(500,mu=rep(0,16),Sigma=kinX[1:16,1:16])
  sigmaAMat <- mvrnorm(500,mu=rep(0,16),Sigma=kinAutos[1:16,1:16])
  
  noise <- rnorm(n,mean=0,sd=1)
  
  TEMP=c(2,1,2,1,1,1,2,2,1,2,1,2,1,2,1,2)
  SEX <- rep(NA,length(TEMP))
  SEX[TEMP==1]="M"
  SEX[TEMP==2]="F"
    
  causalSNP <- Family_alleles_NmarkerX(0.2,1,SEX)
  
  pheno <- noise+0.05*causalSNP+sigmaA*as.vector(t(sigmaAMat))+sigmaX*as.vector(t(sigmaXMat))
  
  aFreqs <- rep(0.2,numSNPs)
  geno <- matrix(NA,nrow=numSNPs,ncol=n)
  
  simGenos <- vector("list",500)
  for(i in 1:500){
    geno[,(i*16-15):(i*16)] <- Family_alleles_NmarkerX(rep(0.2,numSNPs),numSNPs,SEX) # simulate numSNPs with MAF 0.2
  }
  
  scan <- data.frame(scanID=1:n,sex=sex,pheno=pheno)
  scan <- ScanAnnotationDataFrame(scan)
  
  genoMt <- MatrixGenotypeReader(genotype=geno,snpID=as.integer(1:nrow(geno)),
                                 chromosome=as.integer(rep(23,nrow(geno))),
                                 position=as.integer(1:nrow(geno)), scanID=scan$scanID)
  
  genoData <- GenotypeData(genoMt,scanAnnot=scan)
    
  
  snpStart <- 1
  snpEnd <- numSNPs
  
  colnames(kinX) <- getScanID(genoData)
  colnames(kinAuto) <- getScanID(genoData)
  
  covMatListBoth <- list(kinX,kinAuto) # this is the X chr kinship matrix for all samples
  names(covMatListBoth) <- c("kinshipX","kinshipAutos")
  
    startBoth <- c(0.05^2*4*0.2*(1-0.2)+sigmaX,sigmaA,1)
    varCompBoth <- estVarComp(scan,covMatList=covMatListBoth,"pheno",start=startBoth,AIREML.tol=1e-4)
    varBothci <- estVarCompCI(varCompBoth,prop=FALSE)
    
    cholSig <- varCompBoth[["cholSigmaInv"]]
    mmResBoth <- assocTestMixedModel(genoData,snpStart=snpStart,snpEnd=snpEnd,cholSig,outcome="pheno")
    mmResBoth$causal <- FALSE
  mmResBoth$model <- "both"
    
  startAuto <- c(0.05^2*4*0.2*(1-0.2)+sigmaX + sigmaA, 1)
    varCompA <- estVarComp(scan,covMatList=list(kinAuto),"pheno",start=startAuto,AIREML.tol=1e-4)
    varAutoci <- estVarCompCI(varCompA,prop=FALSE)
    
    cholSig <- varCompA[["cholSigmaInv"]]
    mmResA <- assocTestMixedModel(genoData,snpStart=snpStart,snpEnd=snpEnd,cholSig,outcome="pheno")
    mmResA$causal <- FALSE
  mmResA$model <- "auto"
    
    startX <- c(0.05^2*4*0.2*(1-0.2)+sigmaX, 1)
    varCompX <- estVarComp(scan,covMatList=list(kinX),"pheno",start=startX,AIREML.tol=1e-4)
    varXci <- estVarCompCI(varCompX,prop=FALSE)
    
    cholSig <- varCompX[["cholSigmaInv"]]
    mmResX <- assocTestMixedModel(genoData,snpStart=snpStart,snpEnd=snpEnd,cholSig,outcome="pheno")
    mmResX$causal <- FALSE
  mmResX$model <- "x"
      
  res <- list(sigmaA,sigmaX,mmResBoth,varBothci,mmResA,varAutoci,mmResX,varXci)
  names(res) <- c("sigma_auto","sigma_x","mm_res_both","var_est_ci_both","mm_res_auto","var_est_ci_auto",
                  "mm_res_x","var_est_ci_x")
  
  return(res)
}



setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")
library(GWASTools)
library(MASS)
source("allele_drop_functions.R")

res <- test_phenotypes(2,0.5,0.5)
res

q("no")
  