

# simulate phenotypes for 500 SNPs
# for each SNP, simulate 3 phenotypes:
# A: y=beta1*SNPx + noise
# B: y=beta1*SNPx + g_x + noise
# C: y=beta1*SNPx + g_x + g_a + noise

# where
# noise ~ N(0,1)
# g_x ~ MVN(0,kinXchr*sigmaX)
# g_a ~ MVN(0,kinAutos*sigmaA)

# can use the same noise each iteration
# can use the same g_x for models B, C

phenotypes_randEffects <- function(fn,fnX,fnA){

  n <- 8000
  mu <- rep(0,n)

  kinAutos <- get(load("500Peds_autoKinship.RData"))
  kinX <- get(load("500Peds_xKinship.RData"))

  kinAutos <- matrix(kinAutos,nrow=n,ncol=n)
  kinX <- matrix(kinX,nrow=n,ncol=n)

  TEMP=c(2,1,2,1,1,1,2,2,1,2,1,2,1,2,1,2)
  SEX <- rep(NA,length(TEMP))
  SEX[TEMP==1]="M"
  SEX[TEMP==2]="F"
  
  sex <- rep(SEX,(n/16))

  scan <- data.frame(matrix(NA,nrow=8000,ncol=500000+2))
  colnames(scan) <- c("scanID","sex",paste("noise",1:500000,"_modelC",sep=""))
  scan$scanID <- 1:n
  scan$sex <- sex
  
  scanSigmaX <- scan
  scanSigmaA <- scan
  
  save(scan,file=paste("tmp_",fn,sep=""))

  #kinAutosSigma <- kinAutos*as.numeric(sigmaA)
  #kinXSigma <- kinX*sigmaX

# do this for 1M phenotypes
# break into blocks of 5000

  ct <- 1

  for(i in 1:1000){
  # first draw the 500 block diagonal entries from a mvn with covar structure
    sigmaXMat <- mvrnorm(500*500,mu=rep(0,16),Sigma=kinX[1:16,1:16])
    sigmaAMat <- mvrnorm(500*500,mu=rep(0,16),Sigma=kinAutos[1:16,1:16])

    noise <- rnorm(500*n,mean=0,sd=1)
  
    scan[,ct:(ct+500-1)] <- noise
    scanSigmaX[,ct:(ct+500-1)] <- sigmaXMat
    scanSigmaA[,ct:(ct+500-1)] <- sigmaAMat
    ct <- ct+500
  
    print(ct)
    print(i)
  }

  save(scan,file=fn)
save(scanSigmaX,file=fnX)
save(scanSigmaA,file=fnA)
  return(NULL)
  
}




setwd("/projects/geneva/geneva_sata/caitlin/mlm_x/")
library(MASS)
library(GWASTools)

phenotypes_randEffects(fn="phenotypes_genotypes_halfMillion2_noise.RData",
                       fnX="phenotypes_genotypes_halfMillion2_sigmaX05.RData",
                       fnA="phenotypes_genotypes_halfMillion2_sigmaA05.RData")

rm(list=ls())

