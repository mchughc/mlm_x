
source("allele_drop_functions.R")
source("assocTestMixedModel_v7.R")
source("estVarComp.R")
library(GWASTools)

############################
## Contents:
# 1. Estimate X chr relatedness
# 2. Estimate for autosomal relatedness 
# 10. Make histograms of X chr KC estimates 




#####
# 1. Estimate X chr relatedness

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

# estimate 1000 xchr loci
set.seed(6422)
genoX <- 2*Family_alleles_NmarkerX(rep(f1,1000),1000,SEX)

kinship <- matrix(NA,nrow=16,ncol=16)
for(i in 1:nrow(genoX)){
  for(j in 1:nrow(genoX)){
    kinship[i,j] <- getKinX(genoX[i,],genoX[j,],SEX[i],SEX[j],f1)
  }
}

summary(diag(kinship))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.9579  0.9862  0.9979  0.9977  1.0100  1.0320 

##
# what do we expect for these individuals?
# read in true kinship coefficients, multiply by 2 to get relatedness coef
tmp <- read.table("pedigree_16individs_output")
trueKinX <- matrix(NA,nrow=16, ncol=16)
trueKinX[lower.tri(trueKinX,diag=TRUE)] <- tmp[,"V4"]
trueKinXn <- trueKinX
trueKinX <- data.frame(trueKinX)
trueKinX$SEX <- SEX
trueKinXn[trueKinX$SEX=="F",trueKinX$SEX=="F"] <- 2*trueKinXn[trueKinX$SEX=="F",trueKinX$SEX=="F"]

# is this correct???
trueKinXn[trueKinX$SEX=="M",trueKinX$SEX=="F"] <- trueKinXn[trueKinX$SEX=="M",trueKinX$SEX=="F"]*sqrt(2)
trueKinXn[trueKinX$SEX=="F",trueKinX$SEX=="M"] <- trueKinXn[trueKinX$SEX=="F",trueKinX$SEX=="M"]*sqrt(2)

diag(trueKinXn) <- 1

# estimate 10000 xchr loci now
set.seed(6422)
genoX <- 2*Family_alleles_NmarkerX(rep(f1,10000),10000,SEX)

kinship10K <- matrix(NA,nrow=16,ncol=16)
for(i in 1:nrow(genoX)){
  for(j in 1:nrow(genoX)){
    kinship10K[i,j] <- getKinX(genoX[i,],genoX[j,],SEX[i],SEX[j],f1)
  }
}

summary(diag(kinship10K))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.9705  0.9957  0.9986  0.9961  1.0010  1.0110 

genoX <- 2*Family_alleles_NmarkerX(rep(f1,100000),100000,SEX)
kinship100K <- matrix(NA,nrow=16,ncol=16)
for(i in 1:nrow(genoX)){
  for(j in 1:nrow(genoX)){
    kinship100K[i,j] <- getKinX(genoX[i,],genoX[j,],SEX[i],SEX[j],f1)
  }
}
summary(diag(kinship100K))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.9968  1.0000  1.0010  1.0010  1.0020  1.0040 
 
genoX <- 2*Family_alleles_NmarkerX(rep(f1,3500),3500,SEX)
kinship35 <- matrix(NA,nrow=16,ncol=16)
for(i in 1:nrow(genoX)){
  for(j in 1:nrow(genoX)){
    kinship35[i,j] <- getKinX(genoX[i,],genoX[j,],SEX[i],SEX[j],f1)
  }
}
summary(diag(kinship35))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.9726  0.9922  0.9999  0.9992  1.0050  1.0340 


library(grDevices)
pdf("xchr_kc_estimatedVsTrue.pdf")
diffMat <- trueKinXn[lower.tri(trueKinXn,diag=FALSE)]-kinship[lower.tri(kinship,diag=FALSE)]
x <- trueKinXn[lower.tri(trueKinXn,diag=FALSE)]
plot(jitter(x), diffMat,
     xlab="Theoretical KC",ylab="Theoretical KC - Estimated KC",
     xlim=c(0,0.6),col=adjustcolor("purple",alpha=0.6),pch=19)
points(jitter(x),
       x-kinship35[lower.tri(kinship100K,diag=FALSE)], col=adjustcolor("orange",alpha=0.6),pch=19)
points(jitter(x),
       x-kinship10K[lower.tri(kinship10K,diag=FALSE)], col=adjustcolor("red",alpha=0.6),pch=19)
points(jitter(x),
       x-kinship100K[lower.tri(kinship100K,diag=FALSE)], col=adjustcolor("cyan",alpha=0.6),pch=19)
abline(h=0)
legend(0.08,-0.043,c("1k SNPs","3500 SNPs","10k SNPs","100k SNPs"),col=c("purple","orange","red","cyan"),pch=19)
dev.off()


# plot stratified on m-m, m-f and f-f pairs
pdf("xchr_kc_estimatedVsTrue_FF.pdf")
plot(jitter(trueKinXn[trueKinX$SEX=="F",trueKinX$SEX=="F"]),
     trueKinXn[trueKinX$SEX=="F",trueKinX$SEX=="F"]-kinship[trueKinX$SEX=="F",trueKinX$SEX=="F"],
     xlab="Expected KC",ylab="Expected KC - Estimated KC",
     xlim=c(0,0.5),col=adjustcolor("purple",alpha=0.6),pch=19)
points(jitter(trueKinXn[trueKinX$SEX=="F",trueKinX$SEX=="F"]),
       trueKinXn[trueKinX$SEX=="F",trueKinX$SEX=="F"]-kinship35[trueKinX$SEX=="F",trueKinX$SEX=="F"],
       col=adjustcolor("orange",alpha=0.6),pch=19)
points(jitter(trueKinXn[trueKinX$SEX=="F",trueKinX$SEX=="F"]),
       trueKinXn[trueKinX$SEX=="F",trueKinX$SEX=="F"]-kinship10K[trueKinX$SEX=="F",trueKinX$SEX=="F"], 
       col=adjustcolor("red",alpha=0.6),pch=19)
points(jitter(trueKinXn[trueKinX$SEX=="F",trueKinX$SEX=="F"]),
       trueKinXn[trueKinX$SEX=="F",trueKinX$SEX=="F"]-kinship100K[trueKinX$SEX=="F",trueKinX$SEX=="F"],
       col=adjustcolor("cyan",alpha=0.6),pch=19)
abline(h=0)
legend("topleft",c("1k SNPs","3500 SNPs","10k SNPs","100k SNPs"),col=c("purple","orange","red","cyan"),pch=19)
dev.off()


pdf("xchr_kc_estimatedVsTrue_FM.pdf")
plot(jitter(trueKinXn[trueKinX$SEX=="F",trueKinX$SEX=="M"]),
     trueKinXn[trueKinX$SEX=="F",trueKinX$SEX=="M"]-kinship[trueKinX$SEX=="F",trueKinX$SEX=="M"],
     xlab="Expected KC",ylab="Expected KC - Estimated KC",
     col=adjustcolor("purple",alpha=0.6),pch=19)
points(jitter(trueKinXn[trueKinX$SEX=="M",trueKinX$SEX=="F"]),
       trueKinXn[trueKinX$SEX=="M",trueKinX$SEX=="F"]-kinship[trueKinX$SEX=="M",trueKinX$SEX=="F"], 
       col=adjustcolor("purple",alpha=0.6),pch=19)
points(jitter(trueKinXn[trueKinX$SEX=="F",trueKinX$SEX=="M"]),
       trueKinXn[trueKinX$SEX=="F",trueKinX$SEX=="M"]-kinship35[trueKinX$SEX=="F",trueKinX$SEX=="M"], 
       col=adjustcolor("orange",alpha=0.6),pch=19)
points(jitter(trueKinXn[trueKinX$SEX=="M",trueKinX$SEX=="F"]),
       trueKinXn[trueKinX$SEX=="M",trueKinX$SEX=="F"]-kinship35[trueKinX$SEX=="M",trueKinX$SEX=="F"], 
       col=adjustcolor("orange",alpha=0.6),pch=19)
points(jitter(trueKinXn[trueKinX$SEX=="F",trueKinX$SEX=="M"]),
       trueKinXn[trueKinX$SEX=="F",trueKinX$SEX=="M"]-kinship10K[trueKinX$SEX=="F",trueKinX$SEX=="M"], 
       col=adjustcolor("red",alpha=0.6),pch=19)
points(jitter(trueKinXn[trueKinX$SEX=="M",trueKinX$SEX=="F"]),
       trueKinXn[trueKinX$SEX=="M",trueKinX$SEX=="F"]-kinship10K[trueKinX$SEX=="M",trueKinX$SEX=="F"], 
       col=adjustcolor("red",alpha=0.6),pch=19)
points(jitter(trueKinXn[trueKinX$SEX=="F",trueKinX$SEX=="M"]),
       trueKinXn[trueKinX$SEX=="F",trueKinX$SEX=="M"]-kinship100K[trueKinX$SEX=="F",trueKinX$SEX=="M"],
       col=adjustcolor("cyan",alpha=0.6),pch=19)
points(jitter(trueKinXn[trueKinX$SEX=="M",trueKinX$SEX=="F"]),
       trueKinXn[trueKinX$SEX=="M",trueKinX$SEX=="F"]-kinship100K[trueKinX$SEX=="M",trueKinX$SEX=="F"],
       col=adjustcolor("cyan",alpha=0.6),pch=19)
abline(h=0)
legend("bottomright",c("1k SNPs","3500 SNPs","10k SNPs","100k SNPs"),col=c("purple","orange","red","cyan"),pch=19)
dev.off()

pdf("xchr_kc_estimatedVsTrue_MM.pdf")
plot(jitter(trueKinXn[trueKinX$SEX=="M",trueKinX$SEX=="M"]),
     trueKinXn[trueKinX$SEX=="M",trueKinX$SEX=="M"]-kinship[trueKinX$SEX=="M",trueKinX$SEX=="M"],
     xlab="Expected KC",ylab="Expected KC - Estimated KC",#ylim=c(-0.1,0.1),
     xlim=c(-0.01,0.51),col=adjustcolor("purple",alpha=0.6),pch=19)
points(jitter(trueKinXn[trueKinX$SEX=="M",trueKinX$SEX=="M"]),
       trueKinXn[trueKinX$SEX=="M",trueKinX$SEX=="M"]-kinship35[trueKinX$SEX=="M",trueKinX$SEX=="M"], 
       col=adjustcolor("orange",alpha=0.6),pch=19)
points(jitter(trueKinXn[trueKinX$SEX=="M",trueKinX$SEX=="M"]),
       trueKinXn[trueKinX$SEX=="M",trueKinX$SEX=="M"]-kinship10K[trueKinX$SEX=="M",trueKinX$SEX=="M"], 
       col=adjustcolor("red",alpha=0.6),pch=19)
points(jitter(trueKinXn[trueKinX$SEX=="M",trueKinX$SEX=="M"]),
       trueKinXn[trueKinX$SEX=="M",trueKinX$SEX=="M"]-kinship100K[trueKinX$SEX=="M",trueKinX$SEX=="M"],
       col=adjustcolor("cyan",alpha=0.6),pch=19)
abline(h=0)
legend("bottomright",c("1k SNPs","3500 SNPs","10k SNPs","100k SNPs"),col=c("purple","orange","red","cyan"),pch=19)
dev.off()


save(kinship,file="xchr_kinshipMatrix/kinship_1kSNPs.RData")
save(kinship35, file="xchr_kinshipMatrix/kinship_3500SNPs.RData")
save(kinship10K,file="xchr_kinshipMatrix/kinship_10kSNPs.RData")
save(kinship100K,file="xchr_kinshipMatrix/kinship_100kSNPs.RData")


####
# 2. Estimate for autosomal relatedness 

# kinship matrix
library(kinship2)
tmp <- read.table("pedigree_16individs.txt")
tmp[tmp[,5]==0,5] <- 2
ped <- pedigree(tmp[,2],dadid=tmp[,3],momid=tmp[,4],sex=tmp[,5])
kin <- 2*kinship(ped)

# estimate 1000 SNPs first
set.seed(684271)
geno <- 2*Family_alleles_Nmarker(rep(f1,1000),1000)
kinship <- matrix(NA,nrow=16,ncol=16)
for(i in 1:nrow(geno)){
  for(j in 1:nrow(geno)){
    kinship[i,j] <- getKin(geno[i,],geno[j,],f1)
  }
}

summary(diag(kinship))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.9196  0.9704  0.9802  0.9867  1.0010  1.0390 

# now do 10000 SNPs
# estimate 10000 xchr loci now
geno <- 2*Family_alleles_Nmarker(rep(f1,10000),10000)
kinship10K <- matrix(NA,nrow=16,ncol=16)
for(i in 1:nrow(geno)){
  for(j in 1:nrow(geno)){
    kinship10K[i,j] <- getKin(geno[i,],geno[j,],f1)
  }
}

summary(diag(kinship10K))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.9865  0.9906  0.9968  0.9970  1.0010  1.0150 

# now 100K SNPs
geno <- 2*Family_alleles_Nmarker(rep(f1,100000),100000)
kinship100K <- matrix(NA,nrow=16,ncol=16)
for(i in 1:nrow(geno)){
  for(j in 1:nrow(geno)){
    kinship100K[i,j] <- getKin(geno[i,],geno[j,],f1)
  }
}
summary(diag(kinship100K))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.9980  0.9996  1.0010  1.0020  1.0030  1.0070 

# finally 3500 SNPs
geno <- 2*Family_alleles_Nmarker(rep(f1,3500),3500)
kinship3500 <- matrix(NA,nrow=16,ncol=16)
for(i in 1:nrow(geno)){
  for(j in 1:nrow(geno)){
    kinship3500[i,j] <- getKin(geno[i,],geno[j,],f1)
  }
}
summary(diag(kinship3500))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.9688  0.9735  0.9833  0.9850  0.9959  1.0070 

pdf("auto_kc_estimatedVsTrue.pdf")
plot(jitter(kin[lower.tri(kin,diag=FALSE)]),
     -kinship[lower.tri(kinship,diag=FALSE)]+kin[lower.tri(kin,diag=FALSE)],xlab="Expected",
     ylab="Expected - Observed",
     xlim=c(0,1),col=adjustcolor("purple",alpha=0.6),pch=19)
points(jitter(kin[lower.tri(kin,diag=FALSE)]),
       -kinship10K[lower.tri(kinship10K,diag=FALSE)]+kin[lower.tri(kin,diag=FALSE)],
       col=adjustcolor("red",alpha=0.6),pch=19)
points(jitter(kin[lower.tri(kin,diag=FALSE)]),
       -kinship100K[lower.tri(kinship100K,diag=FALSE)]+kin[lower.tri(kin,diag=FALSE)],
       col=adjustcolor("cyan",alpha=0.6),pch=19)
abline(h=0)
legend("topright",c("1k SNPs","10k SNPs","100k SNPs"),col=c("purple","red","cyan"),pch=19)
dev.off()

pdf("auto_kc_estimatedVsTrue_v2.pdf")
plot(jitter(kin[lower.tri(kin,diag=FALSE)]),
     -kinship[lower.tri(kinship,diag=FALSE)]+kin[lower.tri(kin,diag=FALSE)],xlab="Expected",
     ylab="Expected - Observed",
     xlim=c(0,1),col=adjustcolor("purple",alpha=0.6),pch=19)
points(jitter(kin[lower.tri(kin,diag=FALSE)]),
       -kinship3500[lower.tri(kinship3500,diag=FALSE)]+kin[lower.tri(kin,diag=FALSE)],
       col=adjustcolor("orange",alpha=0.6),pch=19)
points(jitter(kin[lower.tri(kin,diag=FALSE)]),
       -kinship10K[lower.tri(kinship10K,diag=FALSE)]+kin[lower.tri(kin,diag=FALSE)],
       col=adjustcolor("red",alpha=0.6),pch=19)
points(jitter(kin[lower.tri(kin,diag=FALSE)]),
       -kinship100K[lower.tri(kinship100K,diag=FALSE)]+kin[lower.tri(kin,diag=FALSE)],
       col=adjustcolor("cyan",alpha=0.6),pch=19)
abline(h=0)
legend("topright",c("1k SNPs","3500 SNPs","10k SNPs","100k SNPs"),col=c("purple","orange","red","cyan"),pch=19)
dev.off()

save(kinship,file="auto_kinshipMatrix/kinship_1kSNPs.RData")
save(kinship10K,file="auto_kinshipMatrix/kinship_10kSNPs.RData")
save(kinship100K,file="auto_kinshipMatrix/kinship_100kSNPs.RData")
save(kinship3500,file="auto_kinshipMatrix/kinship_3500SNPs.RData")


#######################
####
# simulate 1000 unrelated individuals with 100K snps
nloci <- 100000
freq1 <- 0.4
geno <- matrix(0,1000,nloci)
for(count in 1:1000){
  geno[count,] <- 2*INDIV_alleles_Nmarker(rep(freq1,nloci),nloci)
  if(count %% 10==0){ print(paste("on count", count))}
}

kinship100Kunrel <- matrix(NA,nrow=nrow(geno),ncol=nrow(geno))
for(i in 1:nrow(geno)){
  for(j in i:nrow(geno)){
    kinship100Kunrel[i,j] <- getKin(geno[i,],geno[j,],f1)
  }
}
summary(diag(kinship100Kunrel))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.9980  0.9996  1.0010  1.0020  1.0030  1.0070 


###
# simulate phenotype
source("sim_phenotype.R")

set.seed(5321)
causalSNPs <- rbinom(ncol(geno),1,0.0005)
sum(causalSNPs) # 56, ok
causalSNPs <- as.logical(causalSNPs)
sum(causalSNPs)

effSNPs <- rnorm(sum(causalSNPs),sd=0.5,mean=4)
hist(effSNPs)
herit <- 0.1

# geno is 100K autosomal genotypes for all 16 individs
pheno <- sim_phenotype(geno, causalSNPs, effSNPs, herit)

####
# need to call assoc test with true kinship matrices and phenotype, get type I error

df <- data.frame(scanID=tmp[,2],sex=tmp[,5],pheno=pheno)
df$sex[df$sex=="2"] <- "F"
df$sex[df$sex=="1"] <- "M"
scan <- ScanAnnotationDataFrame(df)

# use matt's code
kinship_autos <- get(load("auto_kinshipMatrix/kinship_100kSNPs.RData"))
kinship_x <- get(load("xchr_kinshipMatrix/kinship_100kSNPs.RData"))
dim(kinship_autos); dim(kinship_x) # 16 16 | 16 16

colnames(kinship_autos) <- scan$scanID
colnames(kinship_x) <- scan$scanID

covList <- list(kinship_autos,kinship_x)
names(covList) <- c("kinshipAuto","kinshipX")
varComp <- estVarComp(scan,covMatList=covList,"pheno")

# call MLM with these results
tmpG <- geno
tmpGt <- t(geno)
nSNPs <- 1:ncol(geno)
chromosome <- as.integer(rep(1,length(nSNPs)))
position <- as.integer(rep(0,length(nSNPs)))
genoMt <- MatrixGenotypeReader(genotype=tmpGt,snpID=nSNPs,chromosome=chromosome,
                               position=position, scanID=scan$scanID)
genoData <- GenotypeData(genoMt,scanAnnot=scan)
cholSig <- varComp[["cholSigmaInv"]]
mmRes <- assocTestMixedModel(genoData,snpStart=1,snpEnd=ncol(geno),cholSig,outcome="pheno")


#####
library(kinship2)

ped <- read.table("pedigree_16individs.txt")
ped$V5[ped$V5==0] <- 2
pedPlt <- pedigree(ped$V2,ped$V3,ped$V4,ped$V5)
pdf("pedigree_16individs.pdf")
plot(pedPlt,cex=2,lwd=2)
dev.off()



#######
# look at how many x chr SNPs remained after ld-pruning SNPs in OLGA

# want to merge in chromosome info
snpAnnot <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_HCHS_Custom_15041502_B3_all37_v25_AMS.RData"))

snps <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/amstilp/results/pca/v08_study_unrelated_pcair/pca_snps_sel.v08_study_unrel.RData"))
length(snps) # 153411
snps <- data.frame("snpID"=snps)

snps <- merge(snps,pData(snpAnnot)[,c("snpID","chromosome")],by="snpID",all.x=TRUE)
dim(snps) # 153411 2
table(snps$chromosome) # filtered out all the x chr snps!!!


snps <- get(load(
  "/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/amstilp/results/ibd/v02_after_id_issues/ibd_snp_sel.v2.RData"))
length(snps) # 154906
snps <- data.frame("snpID"=snps)

snps <- merge(snps,pData(snpAnnot)[,c("snpID","chromosome")],by="snpID",all.x=TRUE)
dim(snps) # 153411 2
table(snps$chromosome) # wow, only 102 SNPs on the X chromosome
# these automatically remove the non-autosomal SNPs

# need to do some pruning on my own, i guess...
snp.ids <- snpAnnot$snpID[snpAnnot$chromosome==23]
length(snp.ids) # 55905

scanAnnot <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_scanAnnot_v12_AMS.RData"))
scan.ids <- scanAnnot$scanID
gdsobj <- openfn.gds("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/netCDF/samples/OLGA_geno.gds")
library(SNPRelate)
snpset <- snpgdsLDpruning(gdsobj, sample.id=scan.ids, snp.id=snp.ids,
                          autosome.only=FALSE, maf=0.05, missing.rate=0.05,
                          method="corr", slide.max.bp=10*1e6, ld.threshold=0.32, 
                          num.thread=1)
# 3760 SNPs selected

rm(list=ls())


#####
# Make histograms of X chr KC estimates 

TEMP=c(2,1,2,1,1,1,2,2,1,2,1,2,1,2,1,2)
SEX=rep(NA,16)
SEX[TEMP==1]="M"
SEX[TEMP==2]="F"

tmp <- read.table("pedigree_16individs_output")
trueKinX <- matrix(NA,nrow=16, ncol=16)
trueKinX[lower.tri(trueKinX,diag=TRUE)] <- tmp[,"V4"]

diag(trueKinX)[SEX=="F"] <- 0.5
diag(trueKinX)[SEX=="M"] <- 1

est35 <- get(load("xchr_kinshipMatrix/kinship_3500SNPs.RData"))
est35[SEX=="F",SEX=="F"] <- 0.5*est35[SEX=="F",SEX=="F"]
est35[SEX=="F",SEX=="M"] <- (1/sqrt(2))*est35[SEX=="F",SEX=="M"]
est35[SEX=="M",SEX=="F"] <- (1/sqrt(2))*est35[SEX=="M",SEX=="F"]

# ok, now est35 is transformed GRM estimates to X chr KC estimates
# trueKinX is lower.tri matrix of theoretical X chr KC estimates
# make histogram of each relationship type

pdf("hist_xchrKC_byRelType.pdf",width=11)
par(mfrow=c(2,4))
toSel <- trueKinX==0
hist(est35[toSel],xlab="Estimated X Chr KC",main=paste("Histogram of ", sum(toSel,na.rm=T)," Pairs\nwith Theoretical KC 0",sep=""))
abline(v=0,col="red",lty=2,lwd=2)

# toSel <- trueKinX==(1/16) # 0.0625
# hist(est35[toSel],xlab="Estimated X Chr KC",main=paste("Histogram of ", sum(toSel,na.rm=T)," Pairs\nwith Theoretical KC 1/16",sep=""))
# abline(v=0.0625,col="red",lty=2,lwd=2)

# toSel <- trueKinX==(3/32) # 0.09375
# hist(est35[toSel],xlab="Estimated X Chr KC",main=paste("Histogram of ", sum(toSel,na.rm=T)," Pairs\nwith Theoretical KC 3/32",sep=""))
# abline(v=0.09375,col="red",lty=2,lwd=2)

toSel <- trueKinX==(1/8) # 0.125
hist(est35[toSel],xlab="Estimated X Chr KC",main=paste("Histogram of ", sum(toSel,na.rm=T)," Pairs\nwith Theoretical KC 1/8",sep=""))
abline(v=0.125,col="red",lty=2,lwd=2)

toSel <- trueKinX==(3/16) # 0.1875
hist(est35[toSel],xlab="Estimated X Chr KC",main=paste("Histogram of ", sum(toSel,na.rm=T)," Pairs\nwith Theoretical KC 3/16",sep=""))
abline(v=0.1875,col="red",lty=2,lwd=2)

toSel <- trueKinX==0.25
hist(est35[toSel],xlab="Estimated X Chr KC",main=paste("Histogram of ", sum(toSel,na.rm=T)," Pairs\nwith Theoretical KC 1/4",sep=""))
abline(v=0.25,col="red",lty=2,lwd=2)

toSel <- trueKinX==(6/16) # 0.375
hist(est35[toSel],xlab="Estimated X Chr KC",main=paste("Histogram of ", sum(toSel,na.rm=T)," Pairs\nwith Theoretical KC 6/16",sep=""))
abline(v=0.375,col="red",lty=2,lwd=2)

toSel <- trueKinX==0.5
hist(est35[toSel],xlab="Estimated X Chr KC",main=paste("Histogram of ", sum(toSel,na.rm=T)," Pairs\nwith Theoretical KC 1/2",sep=""))
abline(v=0.5,col="red",lty=2,lwd=2)

toSel <- trueKinX==1
hist(est35[toSel],xlab="Estimated X Chr KC",main=paste("Histogram of ", sum(toSel,na.rm=T)," Pairs\nwith Theoretical KC 1",sep=""))
abline(v=1,col="red",lty=2,lwd=2)

dev.off()




