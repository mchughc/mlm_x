
## pedigree simulations using estimated allele frequencies, rather than the simulated values

source("allele_drop_functions.R")
source("assocTestMixedModel_v7.R")
source("estVarComp.R")
library(GWASTools)

############################


FAMNUM1=60;
unrelateds=200
nloci=100
people=FAMNUM1*16+unrelateds

TEMP=c(2,1,2,1,1,1,2,2,1,2,1,2,1,2,1,2)
SEX=rep(NA,16)
SEX[TEMP==1]="M"
SEX[TEMP==2]="F"


f1=.4

Family_alleles_Nmarker(f1,1)
Family_alleles_NmarkerX(f1,1,SEX)
cbind(Family_alleles_NmarkerX(f1,1,SEX),SEX)

##
# what do we expect for these individuals?
# read in true kinship coefficients, multiply by 2 or 1.5 to get relatedness coef
tmp <- read.table("pedigree_16individs_output")
trueKinX <- matrix(NA,nrow=16, ncol=16)
trueKinX[lower.tri(trueKinX,diag=TRUE)] <- tmp[,"V4"]
trueKinXn <- trueKinX
trueKinX <- data.frame(trueKinX)
trueKinX$SEX <- SEX
diag(trueKinX) <- 1
trueKinXn[trueKinX$SEX=="F",trueKinX$SEX=="F"] <- 2*trueKinXn[trueKinX$SEX=="F",trueKinX$SEX=="F"]

# is this correct???
trueKinXn[trueKinX$SEX=="M",trueKinX$SEX=="F"] <- 1.5*trueKinXn[trueKinX$SEX=="M",trueKinX$SEX=="F"]
trueKinXn[trueKinX$SEX=="F",trueKinX$SEX=="M"] <- 1.5*trueKinXn[trueKinX$SEX=="F",trueKinX$SEX=="M"]
###########################


# estimate 1000 xchr loci
set.seed(6422)
genoX <- 2*Family_alleles_NmarkerX(rep(f1,1000),1000,SEX)
aFreq <- apply(genoX,2,function(x){sum(x)/24})
kinshipf <- matrix(NA,nrow=16,ncol=16)
for(i in 1:nrow(genoX)){
  for(j in 1:nrow(genoX)){
    kinshipf[i,j] <- getKinX(genoX[i,],genoX[j,],SEX[i],SEX[j],aFreq)
  }
}
summary(diag(kinshipf))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.5341  0.7530  0.8283  0.8981  1.1060  1.4620 

# do 10K snps
genoX <- 2*Family_alleles_NmarkerX(rep(f1,10000),10000,SEX)
aFreq <- apply(genoX,2,function(x){sum(x)/24})
summary(aFreq)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.0000  0.2500  0.3750  0.3998  0.5417  1.0000 

kinship10Kf <- matrix(NA,nrow=16,ncol=16)
for(i in 1:nrow(genoX)){
  for(j in 1:nrow(genoX)){
    kinship10Kf[i,j] <- getKinX(genoX[i,],genoX[j,],SEX[i],SEX[j],aFreq)
  }
}
summary(diag(kinship10Kf))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.5382  0.7599  0.8189  0.8962  1.1570  1.4000 

# do 100K snps
genoX <- 2*Family_alleles_NmarkerX(rep(f1,100000),100000,SEX)
aFreq <- apply(genoX,2,function(x){sum(x)/24})
summary(aFreq)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.0000  0.2500  0.3750  0.4008  0.5417  1.0000 

kinship100Kf <- matrix(NA,nrow=16,ncol=16)
for(i in 1:nrow(genoX)){
  for(j in 1:nrow(genoX)){
    kinship100Kf[i,j] <- getKinX(genoX[i,],genoX[j,],SEX[i],SEX[j],aFreq)
  }
}
summary(diag(kinship100Kf))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.5362  0.7566  0.8207  0.8992  1.1720  1.3800 


library(grDevices)
pdf("xchr_kc_estimatedVsTrue_usingEstAlleleFreqs.pdf")
x <- trueKinXn[lower.tri(trueKinXn,diag=FALSE)]
plot(jitter(x),
     x-kinshipf[lower.tri(kinshipf,diag=FALSE)],xlab="Expected",
     ylab="Expected - Observed",ylim=c(-0.2,0.6),
     xlim=c(0,1),col=adjustcolor("purple",alpha=0.6),pch=19)
points(jitter(x),
       x-kinship10Kf[lower.tri(kinship10Kf,diag=FALSE)],col=adjustcolor("red",alpha=0.6),pch=19)
points(jitter(x),
       x-kinship100Kf[lower.tri(kinship100Kf,diag=FALSE)],col=adjustcolor("cyan",alpha=0.6),pch=19)
abline(h=0)
legend("topright",c("1k SNPs","10k SNPs","100k SNPs"),col=c("purple","red","cyan"),pch=19)
dev.off()

## what are the pairs/relationships that are the most underestimated?
y <- x-kinship100Kf[lower.tri(kinship100Kf,diag=FALSE)]
which(y>0.4&trueKinXn[lower.tri(trueKinXn,diag=FALSE)]==0.75)
# 76 and 91; maternal uncle-niece and father-son

# what pairs are performing well?
which(y<0.05&trueKinXn[lower.tri(trueKinXn,diag=FALSE)]==0.5)
# 94 -- maternal uncle-niece
which(y<0.05&trueKinXn[lower.tri(trueKinXn,diag=FALSE)]==0.75)
# 37,35,50

tmp <- matrix(1:length(tmp),nrow=16,ncol=16)
# 94 is 6-14, 76 is 12-5 and 91 is 6-11
# 37 is 5-3, 35 is 3-3 and 50 is 2-4


######
# now simulate autosomal SNPs

# kinship matrix
library(kinship2)
tmp <- read.table("pedigree_16individs.txt")
tmp[tmp[,5]==0,5] <- 2
ped <- pedigree(tmp[,2],dadid=tmp[,3],momid=tmp[,4],sex=tmp[,5])
kin <- 2*kinship(ped)

# estimate 1000 SNPs first
set.seed(684271)
geno <- 2*Family_alleles_Nmarker(rep(f1,1000),1000)
aFreq <- apply(geno,2,function(x){sum(x)/32})
kinship <- matrix(NA,nrow=16,ncol=16)
for(i in 1:nrow(geno)){
  for(j in 1:nrow(geno)){
    kinship[i,j] <- getKin(geno[i,],geno[j,],aFreq)
  }
}
summary(diag(kinship))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.7581  0.9886  1.0620  1.1000  1.1340  1.5720 

# now do 10000 SNPs
# estimate 10000 xchr loci now
geno <- 2*Family_alleles_Nmarker(rep(f1,10000),10000)
aFreq <- apply(geno,2,function(x){sum(x)/32})
kinship10K <- matrix(NA,nrow=16,ncol=16)
for(i in 1:nrow(geno)){
  for(j in 1:nrow(geno)){
    kinship10K[i,j] <- getKin(geno[i,],geno[j,],aFreq)
  }
}
summary(diag(kinship10K))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.841   1.050   1.063   1.129   1.131   1.634 

# now 100K SNPs
geno <- 2*Family_alleles_Nmarker(rep(f1,100000),100000)
aFreq <- apply(geno,2,function(x){sum(x)/32})
kinship100K <- matrix(NA,nrow=16,ncol=16)
for(i in 1:nrow(geno)){
  for(j in 1:nrow(geno)){
    kinship100K[i,j] <- getKin(geno[i,],geno[j,],aFreq)
  }
}
summary(diag(kinship100K))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.824   1.061   1.066   1.131   1.145   1.641 

pdf("auto_kc_estimatedVsTrue_usingEstAlleleFreqs.pdf")
plot(jitter(kin[lower.tri(kin,diag=FALSE)]),
     -kinship[lower.tri(kinship,diag=FALSE)]+kin[lower.tri(kin,diag=FALSE)],xlab="Expected",
   ylab="Expected - Observed",ylim=c(-0.5,0.5),
     xlim=c(0,1),col=adjustcolor("purple",alpha=0.6),pch=19)
points(jitter(kin[lower.tri(kin,diag=FALSE)]),
       -kinship10K[lower.tri(kinship10K,diag=FALSE)]+kin[lower.tri(kin,diag=FALSE)],col=adjustcolor("red",alpha=0.6),pch=19)
points(jitter(kin[lower.tri(kin,diag=FALSE)]),
       -kinship100K[lower.tri(kinship100K,diag=FALSE)]+kin[lower.tri(kin,diag=FALSE)],col=adjustcolor("cyan",alpha=0.6),pch=19)
abline(h=0)
legend("topright",c("1k SNPs","10k SNPs","100k SNPs"),col=c("purple","red","cyan"),pch=19)
dev.off()


############
# make a few pedigrees and rbind them together, maybe sample size is too small?
tmp <- rbind(tmp,tmp,tmp,tmp,tmp)
tmp$V1 <- c(rep(1:5,each=16))
ped <- pedigree(tmp[,2],dadid=tmp[,3],momid=tmp[,4],sex=tmp[,5],famid=tmp[,1])
kin <- 2*kinship(ped)
kin <- as.matrix(kin)

geno1 <- 2*Family_alleles_Nmarker(rep(f1,100000),100000)
geno2 <- 2*Family_alleles_Nmarker(rep(f1,100000),100000)
geno3 <- 2*Family_alleles_Nmarker(rep(f1,100000),100000)
geno4 <- 2*Family_alleles_Nmarker(rep(f1,100000),100000)
geno5 <- 2*Family_alleles_Nmarker(rep(f1,100000),100000)
geno <- rbind(geno1,geno2,geno3,geno4,geno5)
dim(geno) # 80 100000
aFreq <- apply(geno,2,function(x){sum(x)/(16*2*5)})
kinship100K <- matrix(NA,nrow=nrow(geno),ncol=nrow(geno))
for(i in 1:nrow(geno)){
  for(j in i:nrow(geno)){
    kinship100K[i,j] <- getKin(geno[i,],geno[j,],aFreq)
  }
}
summary(diag(kinship100K))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.9280  0.9611  0.9664  0.9726  0.9724  1.0400 


geno1 <- 2*Family_alleles_Nmarker(rep(f1,10000),10000)
geno2 <- 2*Family_alleles_Nmarker(rep(f1,10000),10000)
geno3 <- 2*Family_alleles_Nmarker(rep(f1,10000),10000)
geno4 <- 2*Family_alleles_Nmarker(rep(f1,10000),10000)
geno5 <- 2*Family_alleles_Nmarker(rep(f1,10000),10000)
geno <- rbind(geno1,geno2,geno3,geno4,geno5)
dim(geno) # 80 10000
aFreq <- apply(geno,2,function(x){sum(x)/(16*2*5)})
kinship10Kf <- matrix(NA,nrow=nrow(geno),ncol=nrow(geno))
for(i in 1:nrow(geno)){
  for(j in i:nrow(geno)){
    kinship10Kf[i,j] <- getKin(geno[i,],geno[j,],aFreq)
  }
}
summary(diag(kinship10Kf))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.9095  0.9484  0.9635  0.9689  0.9746  1.0530 

geno1 <- 2*Family_alleles_Nmarker(rep(f1,1000),1000)
geno2 <- 2*Family_alleles_Nmarker(rep(f1,1000),1000)
geno3 <- 2*Family_alleles_Nmarker(rep(f1,1000),1000)
geno4 <- 2*Family_alleles_Nmarker(rep(f1,1000),1000)
geno5 <- 2*Family_alleles_Nmarker(rep(f1,1000),1000)
geno <- rbind(geno1,geno2,geno3,geno4,geno5)
dim(geno) # 80 1000
aFreq <- apply(geno,2,function(x){sum(x)/(16*2*5)})
kinship1Kf <- matrix(NA,nrow=nrow(geno),ncol=nrow(geno))
for(i in 1:nrow(geno)){
  for(j in i:nrow(geno)){
    kinship1Kf[i,j] <- getKin(geno[i,],geno[j,],aFreq)
  }
}
summary(diag(kinship1Kf))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.8771  0.9378  0.9709  0.9674  0.9964  1.0560 

pdf("auto_kc_estimatedVsTrue_usingEstAlleleFreqs_80samp.pdf")
plot(jitter(kin[lower.tri(kin,diag=FALSE)]),
     -kinship1Kf[upper.tri(kinship1Kf,diag=FALSE)]+kin[lower.tri(kin,diag=FALSE)],xlab="Expected",
     ylab="Expected - Observed",ylim=c(-0.75,0.75),
     xlim=c(0,1),col=adjustcolor("purple",alpha=0.6),pch=19)
points(jitter(kin[lower.tri(kin,diag=FALSE)]),
       -kinship10Kf[upper.tri(kinship10Kf,diag=FALSE)]+kin[lower.tri(kin,diag=FALSE)],col=adjustcolor("red",alpha=0.6),pch=19)
points(jitter(kin[lower.tri(kin,diag=FALSE)]),
       -kinship100K[upper.tri(kinship100K,diag=FALSE)]+kin[lower.tri(kin,diag=FALSE)],col=adjustcolor("cyan",alpha=0.6),pch=19)
abline(h=0)
legend("topright",c("1k SNPs","10k SNPs","100k SNPs"),col=c("purple","red","cyan"),pch=19)
dev.off()


#### 
# do for x chr snps now

SEX <- rep(SEX,5)

geno1 <- 2*Family_alleles_NmarkerX(rep(f1,100000),100000,SEX)
geno2 <- 2*Family_alleles_NmarkerX(rep(f1,100000),100000,SEX)
geno3 <- 2*Family_alleles_NmarkerX(rep(f1,100000),100000,SEX)
geno4 <- 2*Family_alleles_NmarkerX(rep(f1,100000),100000,SEX)
geno5 <- 2*Family_alleles_NmarkerX(rep(f1,100000),100000,SEX)
geno <- rbind(geno1,geno2,geno3,geno4,geno5)
dim(geno) # 80 100000
aFreq <- apply(geno,2,function(x){sum(x)/(8*2*5+8*5)})
kinship100K <- matrix(NA,nrow=nrow(geno),ncol=nrow(geno))
for(i in 1:nrow(geno)){
  for(j in i:nrow(geno)){
    kinship100K[i,j] <- getKinX(geno[i,],geno[j,],SEX[i],SEX[j],aFreq)
  }
}
summary(diag(kinship100K))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.9140  0.9610  0.9736  0.9817  1.0250  1.0480 


geno1 <- 2*Family_alleles_NmarkerX(rep(f1,10000),10000,SEX)
geno2 <- 2*Family_alleles_NmarkerX(rep(f1,10000),10000,SEX)
geno3 <- 2*Family_alleles_NmarkerX(rep(f1,10000),10000,SEX)
geno4 <- 2*Family_alleles_NmarkerX(rep(f1,10000),10000,SEX)
geno5 <- 2*Family_alleles_NmarkerX(rep(f1,10000),10000,SEX)
geno <- rbind(geno1,geno2,geno3,geno4,geno5)
dim(geno) # 80 10000
aFreq <- apply(geno,2,function(x){sum(x)/(8*2*5+8*5)})
kinship10Kf <- matrix(NA,nrow=nrow(geno),ncol=nrow(geno))
for(i in 1:nrow(geno)){
  for(j in i:nrow(geno)){
    kinship10Kf[i,j] <- getKinX(geno[i,],geno[j,],SEX[i],SEX[j],aFreq)
  }
}
summary(diag(kinship10Kf))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.9114  0.9598  0.9719  0.9811  1.0200  1.0470 

geno1 <- 2*Family_alleles_NmarkerX(rep(f1,1000),1000,SEX)
geno2 <- 2*Family_alleles_NmarkerX(rep(f1,1000),1000,SEX)
geno3 <- 2*Family_alleles_NmarkerX(rep(f1,1000),1000,SEX)
geno4 <- 2*Family_alleles_NmarkerX(rep(f1,1000),1000,SEX)
geno5 <- 2*Family_alleles_NmarkerX(rep(f1,1000),1000,SEX)
geno <- rbind(geno1,geno2,geno3,geno4,geno5)
dim(geno) # 80 1000
aFreq <- apply(geno,2,function(x){sum(x)/(8*5+8*2*5)})
kinship1Kf <- matrix(NA,nrow=nrow(geno),ncol=nrow(geno))
for(i in 1:nrow(geno)){
  for(j in i:nrow(geno)){
    kinship1Kf[i,j] <- getKinX(geno[i,],geno[j,],SEX[i],SEX[j],aFreq)
  }
}
summary(diag(kinship1Kf))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.8866  0.9573  0.9848  0.9842  1.0110  1.0940 

tmp <- read.table("pedigree_80individs_output.txt")
trueKinX <- matrix(0,nrow=80, ncol=80)
for(i in 1:5){
  tmpMat <- trueKinX[(i*16-15):(i*16),(i*16-15):(i*16)]
  tmpMat[lower.tri(tmpMat,diag=TRUE)] <- tmp[tmp[,"V1"]==i,"V4"]
  trueKinX[(i*16-15):(i*16),(i*16-15):(i*16)] <- tmpMat
}
trueKinXn <- trueKinX
trueKinX <- data.frame(trueKinX)
trueKinX$SEX <- SEX
trueKinXn[trueKinX$SEX=="F",trueKinX$SEX=="F"] <- 2*trueKinXn[trueKinX$SEX=="F",trueKinX$SEX=="F"]

# is this correct???
trueKinXn[trueKinX$SEX=="M",trueKinX$SEX=="F"] <- 1.5*trueKinXn[trueKinX$SEX=="M",trueKinX$SEX=="F"]
trueKinXn[trueKinX$SEX=="F",trueKinX$SEX=="M"] <- 1.5*trueKinXn[trueKinX$SEX=="F",trueKinX$SEX=="M"]


pdf("xchr_kc_estimatedVsTrue_usingEstAlleleFreqs_80samp.pdf")
plot(jitter(trueKinXn[lower.tri(trueKinXn,diag=FALSE)]),
     -kinship1Kf[upper.tri(kinship1Kf,diag=FALSE)]+trueKinXn[lower.tri(trueKinXn,diag=FALSE)],xlab="Expected",
     ylab="Expected - Observed",ylim=c(-1,1),
     xlim=c(0,1),col=adjustcolor("purple",alpha=0.6),pch=19)
points(jitter(trueKinXn[lower.tri(trueKinXn,diag=FALSE)]),
       -kinship10Kf[upper.tri(kinship10Kf,diag=FALSE)]+trueKinXn[lower.tri(trueKinXn,diag=FALSE)],col=adjustcolor("red",alpha=0.6),pch=19)
points(jitter(trueKinXn[lower.tri(trueKinXn,diag=FALSE)]),
       -kinship100K[upper.tri(kinship100K,diag=FALSE)]+trueKinXn[lower.tri(trueKinXn,diag=FALSE)],col=adjustcolor("cyan",alpha=0.6),pch=19)
abline(h=0)
legend("topright",c("1k SNPs","10k SNPs","100k SNPs"),col=c("purple","red","cyan"),pch=19)
dev.off()