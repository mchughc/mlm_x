
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

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x/")

sigmaA <- 0.3
sigmaX <- 0.8
beta1 <- 0.8

genotypes <- get(load("genotypes_p02_8000samples.RData"))

n <- 8000
mu <- rep(0,n)

kinAutos <- get(load("500Peds_autoKinship.RData"))
kinX <- get(load("500Peds_xKinship.RData"))

TEMP=c(2,1,2,1,1,1,2,2,1,2,1,2,1,2,1,2)
SEX <- rep(NA,length(TEMP))
SEX[TEMP==1]="M"
SEX[TEMP==2]="F"

sex <- rep(SEX,(n/16))
scan <- data.frame(scanID=1:n,sex=sex)

phenotypes <- data.frame(matrix(NA,nrow=8000,ncol=500*3))
colnames(phenotypes) <- paste("pheno_SNP",rep(1:500,each=3),"_model",rep(c("A","B","C"),500),sep="")

scan <- cbind(scan,phenotypes)

kinAutosSigma <- kinAutos*sigmaA
kinXSigma <- kinX*sigmaX

for(i in 1:500){

  causalSNP <- genotypes[i,]
  
  noise <- rnorm(n,mean=0,sd=1)
  sigmaXMat <- mvrnorm(1,mu=mu,Sigma=kinXSigma)
  sigmaAMat <- mvrnorm(1,mu=mu,Sigma=kinAutosSigma)
  pheno <- beta1*causalSNP + noise
  
  thisCol <- paste("pheno_SNP",i,"_modelA",sep="")
  scan[,thisCol] <- pheno
  
  thisCol <- paste("pheno_SNP",i,"_modelB",sep="")
  scan[,thisCol] <- pheno + sigmaXMat
  
  thisCol <- paste("pheno_SNP",i,"_modelC",sep="")
  scan[,thisCol] <- pheno + sigmaXMat + sigmaAMat

}

scanAnnot <- ScanAnnotationDataFrame(scan)
save(scanAnnot,file="scanAnnot_pheno_3models_genotypes_p2_8000samples.RData")

rm(list=ls())
