

# estimate the GRM (genetic relatedness matrix) using only X chr SNPs in hapmap MXL samples
# script started 11 sept 2014
# cpm


## Contents:
# 1. Estimate X chr and autosomal GRM for HM3 samples
# 2. Simulate quantitative trait on X chr
# 3. Create gds file from HM .bed, .bim, .fam
# 4. Run PCA on the samples to use for adjustment
# 5. Run the MLM adjusting for X chr relatedness
# 6. Run the MLM adjusting for X chr relatedness with OR 3.0 phenotype
# 7. Run the MLM not adj for X chr with OR 3.0 phenotype on X chr

# 8. Read in results from running on 100 phenotype simulations


#####

# got this from the gcta website:
#http://www.complextraitgenomics.com/software/gcta/estimate_grm.html

# R script to read the GRM binary file
ReadGRMBin=function(prefix, AllN=F, size=4){
  sum_i=function(i){
    return(sum(1:i))
  }
  BinFileName=paste(prefix,".grm.bin",sep="")
  NFileName=paste(prefix,".grm.N.bin",sep="")
  IDFileName=paste(prefix,".grm.id",sep="")
  id = read.table(IDFileName)
  n=dim(id)[1]
  BinFile=file(BinFileName, "rb");
  grm=readBin(BinFile, n=n*(n+1)/2, what=numeric(0), size=size)
  NFile=file(NFileName, "rb");
  if(AllN==T){
    N=readBin(NFile, n=n*(n+1)/2, what=numeric(0), size=size)
  }
  else N=readBin(NFile, n=1, what=numeric(0), size=size)
  i=sapply(1:n, sum_i)
  return(list(diag=grm[i], off=grm[-i], id=id, N=N))
}


#####
# 1. Estimate X chr and autosomal GRM for HM3 samples

# get hapmap plink files from /projects/geneva/geneva_sata/HapMap/HapMap3_r3/hapmap3_r3_b37_fwd.consensus.qc.poly_Ilmn1M 
# .bed, .bim, .fam

# ran this from command line on fisher:
# /projects/geneva/geneva_sata/apps/gcta_1.24.4/gcta64 --bfile /projects/geneva/geneva_sata/HapMap/HapMap3_r3/hapmap3_r3_b37_fwd.consensus.qc.poly_Ilmn1M --maf 0.01 --make-grm-xchr --out /projects/geneva/geneva_sata/caitlin/grm_xchr/hapmap3_r3_b37_fwd.consensus.qc.poly_Ilmn1M --thread-num 10
# /projects/geneva/geneva_sata/apps/gcta_1.24.4/gcta64 --bfile /projects/geneva/geneva_sata/HapMap/HapMap3_r3/hapmap3_r3_b37_fwd.consensus.qc.poly_Ilmn1M --maf 0.01 --autosome --make-grm --out /projects/geneva/geneva_sata/caitlin/grm_xchr/hapmap3_r3_b37_fwd.consensus.qc.poly_Ilmn1M_autos --thread-num 10

# run the .gz to make sure i'm loading the matrices in correctly
# /projects/geneva/geneva_sata/apps/gcta_1.24.4/gcta64 --bfile /projects/geneva/geneva_sata/HapMap/HapMap3_r3/hapmap3_r3_b37_fwd.consensus.qc.poly_Ilmn1M --maf 0.01 --make-grm-xchr-gz --out /projects/geneva/geneva_sata/caitlin/grm_xchr/hapmap3_r3_b37_fwd.consensus.qc.poly_Ilmn1M --thread-num 10
# yep -- looks good! fewf.

ids <- read.table("/projects/geneva/geneva_sata/caitlin/grm_xchr/hapmap3_r3_b37_fwd.consensus.qc.poly_Ilmn1M.grm.id",header=FALSE,as.is=TRUE)
dim(ids); head(ids) # 1397 2

grmMat <- ReadGRMBin("/projects/geneva/geneva_sata/caitlin/grm_xchr/hapmap3_r3_b37_fwd.consensus.qc.poly_Ilmn1M")
class(grmMat); length(grmMat); names(grmMat) # list w 4 elements: diag, off, id, N
# we have ids, check them

all(grmMat$id==ids) # TRUE

summary(grmMat$diag)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.6485  0.8556  1.0260  1.0770  1.2940  2.2990 
summary(grmMat$off)
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#-0.277500 -0.114300 -0.015810 -0.000749  0.104200  1.305000 

hist(grmMat$diag); dev.off()
hist(grmMat$off); dev.off()

table(table(ids$V1))
idTab <- table(ids$V1)
fam <- names(idTab)[idTab>3]
length(fam) # 25

summary(grmMat$diag[is.element(ids$V1,fam)])

rm(list=ls())


#####
# 2. Simulate quantitative trait on X chr

# Simulate a quantitative trait with the heritability of 0.1 for XX SNP ids with XX effect sizes
# make a list of SNPs with their effect sizes

# include SNPs with varying MAF
# calculate MAF for these samples on x chr
# plink --bfile /projects/geneva/geneva_sata/HapMap/HapMap3_r3/hapmap3_r3_b37_fwd.consensus.qc.poly_Ilmn1M --freq --out /projects/geneva/geneva_sata/caitlin/grm_xchr/hapmap3_r3_b37_fwd.consensus.qc.poly_Ilmn1M

afreq <- read.table("/projects/geneva/geneva_sata/caitlin/grm_xchr/hapmap3_r3_b37_fwd.consensus.qc.poly_Ilmn1M.frq",header=TRUE,as.is=TRUE)
head(afreq)
table(afreq$CHR)

afreq <- afreq[afreq$CHR==23,]
dim(afreq)
summary(afreq$MAF)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.00000 0.02501 0.16110 0.18880 0.32880 0.50000 

# get rsIDs for some of the SNPs with MAF values within what we want
caus <- afreq$SNP[afreq$MAF<0.051&afreq$MAF>0.049]
length(caus) # 41

causInpt <- data.frame("SNP"=caus[15:24],"effectSz"=rep(1.3,10))
write.table(causInpt,file="/projects/geneva/geneva_sata/caitlin/grm_xchr/causSNPs_maf05_effSz13.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)

# Simulate 500 cases and 500 controls with the heritability of liability of 0.5 and disease prevalence of 0.1 for 3 times
# /projects/geneva/geneva_sata/apps/gcta_1.24.4/gcta64  --bfile /projects/geneva/geneva_sata/HapMap/HapMap3_r3/hapmap3_r3_b37_fwd.consensus.qc.poly_Ilmn1M  --simu-qt  --simu-causal-loci /projects/geneva/geneva_sata/caitlin/grm_xchr/causSNPs_maf05_effSz13.txt --simu-rep 100  --out /projects/geneva/geneva_sata/caitlin/grm_xchr/hapmap3_r3_b37_fwd.consensus.qc.poly_Ilmn1M_causSNPs_maf05_effSz13 

## do again with OR of 3.0
caus <- afreq$SNP[afreq$MAF<0.202&afreq$MAF>0.198]
length(caus) # 74

# look at their locations
afreq[is.element(afreq$SNP, caus),]

causInpt <- data.frame("SNP"=caus[11:(11+12)],"effectSz"=rep(3.0,13))
write.table(causInpt,file="/projects/geneva/geneva_sata/caitlin/grm_xchr/causSNPs_maf2_effSz3.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)

# Simulate a trait with the heritability of liability of 0.5 and disease prevalence of 0.1 for 3 times
# /projects/geneva/geneva_sata/apps/gcta_1.24.4/gcta64  --bfile /projects/geneva/geneva_sata/HapMap/HapMap3_r3/hapmap3_r3_b37_fwd.consensus.qc.poly_Ilmn1M  --simu-qt  --simu-causal-loci /projects/geneva/geneva_sata/caitlin/grm_xchr/causSNPs_maf2_effSz3.txt --simu-rep 100  --out /projects/geneva/geneva_sata/caitlin/grm_xchr/hapmap3_r3_b37_fwd.consensus.qc.poly_Ilmn1M_causSNPs_maf2_effSz3 

rm(list=ls())


#####
# 3. Create gds file from HM .bed, .bim, .fam

library(SNPRelate)
bed <- "/projects/geneva/geneva_sata/HapMap/HapMap3_r3/hapmap3_r3_b37_fwd.consensus.qc.poly_Ilmn1M.bed"
bim <- "/projects/geneva/geneva_sata/HapMap/HapMap3_r3/hapmap3_r3_b37_fwd.consensus.qc.poly_Ilmn1M.bim"
fam <- "/projects/geneva/geneva_sata/HapMap/HapMap3_r3/hapmap3_r3_b37_fwd.consensus.qc.poly_Ilmn1M.fam"
outfn <- "/projects/geneva/geneva_sata/caitlin/grm_xchr/hapmap3_r3_b37_fwd.consensus.qc.poly_Ilmn1M.gds"
snpgdsBED2GDS(bed, fam, bim, outfn)

rm(list=ls())


#####
# 4. Run PCA on the samples to use for adjustment

gdsfn <- "/projects/geneva/geneva_sata/caitlin/grm_xchr/hapmap3_r3_b37_fwd.consensus.qc.poly_Ilmn1M.gds"
gdsobj <- openfn.gds(gdsfn)

library(SNPRelate)
rv <- snpgdsPCA(gdsobj, maf=0.05, missing.rate=0.05)

save(rv,file="/projects/geneva/geneva_sata/caitlin/grm_xchr/pca_results.RData")
rm(list=ls())


#####
# 5. Run the MLM adjusting for X chr relatedness

# read in the GRM matrices
library(calibrator)
grmMat <- ReadGRMBin("/projects/geneva/geneva_sata/caitlin/grm_xchr/hapmap3_r3_b37_fwd.consensus.qc.poly_Ilmn1M")
class(grmMat); length(grmMat); names(grmMat) # list w 4 elements: diag, off, id, N

length(grmMat$off); length(grmMat$diag) # need to make this into a matrix

kinship_x <- diag(grmMat$diag)
dim(kinship_x)

kinship_x[upper.tri(kinship_x)] <- grmMat$off
kinship_x <- symmetrize(kinship_x)
isSymmetric(kinship_x) # TRUE
kinship_x[1:10,1:10] # appears to have worked

grmMat <- ReadGRMBin("/projects/geneva/geneva_sata/caitlin/grm_xchr/hapmap3_r3_b37_fwd.consensus.qc.poly_Ilmn1M_autos")
class(grmMat); length(grmMat); names(grmMat) # list w 4 elements: diag, off, id, N

length(grmMat$off); length(grmMat$diag) # need to make this into a matrix

kinship_autos <- diag(grmMat$diag)
dim(kinship_autos)

kinship_autos[upper.tri(kinship_autos)] <- grmMat$off
kinship_autos <- symmetrize(kinship_autos)
isSymmetric(kinship_autos) # TRUE
kinship_autos[1:10,1:10] # appears to have worked

# need to make a scanAnnot object
scan <- get(load("/projects/geneva/geneva_sata/HapMap/HapMap3_r3/sample_snp_annot/hapmap3_r3_b36.sample.annot.v6.RData"))
scan

# make a phenotype column, as simulated with GCTA
sims <- read.table("/projects/geneva/geneva_sata/caitlin/grm_xchr/hapmap3_r3_b37_fwd.consensus.qc.poly_Ilmn1M_causSNPs_maf05_effSz13.phen", header=FALSE,as.is=TRUE)
dim(sims) # 1397 102
head(sims) # ok, so 100 simulated phenotypes, first 2 cols are fam and individ ids

scanMerged <- merge(pData(scan),sims[,c("V2","V3")],by.x="coriell.id",by.y="V2",all.x=TRUE)

scan <- scan[order(scan$coriell.id),]
all(scanMerged$coriell.id==scan$coriell.id) # TRUE

scan$pheno_sim <- scanMerged$V3
scan <- scan[order(scan$scanID),]

# need row and col names of the kinship matrices to be the scanIDs
idMap <- grmMat$id
idMap$order <- 1:nrow(idMap)
idMap <- merge(idMap,pData(scan)[,c("coriell.id","scanID")],by.x="V2",by.y="coriell.id")
idMap <- idMap[order(idMap$order),]
rownames(kinship_autos) <- colnames(kinship_autos) <- idMap$scanID
rownames(kinship_x) <- colnames(kinship_x) <- idMap$scanID

# merge in PC1-5
pcRes <- get(load("/projects/geneva/geneva_sata/caitlin/grm_xchr/pca_results.RData"))
names(pcRes)
evs <- pcRes$eigenvect
dim(evs) # 1397 32

idMap <- idMap[order(idMap$scanID),]
scan <- scan[order(idMap$order),]
all(pcRes$sample.id==scan$coriell.id) # TRUE
scan$PC1 <- evs[,1]
scan$PC2 <- evs[,2]
scan$PC3 <- evs[,3]
scan$PC4 <- evs[,4]
scan$PC5 <- evs[,5]

# now they need to be in the same order
idMap <- idMap[order(idMap$scanID),]
scan <- scan[order(scan$scanID),]
scan <- scan[order(idMap$order),]
head(pData(scan))
head(rownames(kinship_autos))

# use matt's code
source("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/mconomos/R/estVarComp.R")
covList <- list(kinship_autos,kinship_x)
names(covList) <- c("kinshipAuto","kinshipX")
varComp <- estVarComp(scan,covMatList=covList,"pheno_sim",covar.vec=c("PC1","PC2","PC3","PC4","PC5"))
varComp$varComp
#V_kinshipAuto    V_kinshipX           V_E 
#4.732034e+01  1.588115e-08  2.174592e+02 

# see what i get when just adjusting for autosomes
covList <- list(kinship_autos)
names(covList) <- c("kinshipAuto")
varComp <- estVarComp(scan,covMatList=covList,"pheno_sim",covar.vec=c("PC1","PC2","PC3","PC4","PC5"))
varComp$varComp
#V_kinshipAuto           V_E 
#42.57832     222.67430 


## 
# now need to call the MM function

# make genotype data object
data <- GdsGenotypeReader("/projects/geneva/geneva_sata/caitlin/grm_xchr/hapmap3_r3_b37_fwd.consensus.qc.poly_Ilmn1M.gds")
idMap <- idMap[order(idMap$order),]
scan <- scan[order(idMap$order),]
scan$scanID <- scan$coriell.id
genoData <- GenotypeData(data,scanAnnot=scan)

chr <- getChromosome(data)
table(chr)
which.max(chr==25) # 944676
source("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/mconomos/R/assocTestMixedModel_v7.R")

# need to change col/rownames back to coriell ids
cholSig <- varComp[["cholSigmaInv"]]
all(idMap$scanID==colnames(cholSig)) # TRUE
colnames(cholSig) <- rownames(cholSig) <- idMap$V2

mmRes <- assocTestMixedModel(genoData,snpStart=1,snpEnd=(which.max(chr==25)-1),cholSig,outcome="pheno_sim",covar.vec=c("PC1","PC2","PC3","PC4","PC5"))

png("manh_withPCA_res.png")
manhattanPlot(mmRes$pval,mmRes$chr)
dev.off()

# look at the pvals for the causal snps
caus <- read.table("/projects/geneva/geneva_sata/caitlin/grm_xchr/causSNPs_maf05_effSz13.txt",as.is=TRUE)
head(caus)

mmRes[is.element(mmRes$snpID,caus$V1),]

rm(list=ls())


#####
# 6. Run the MLM adjusting for X chr relatedness with OR 3.0 phenotype

# read in the GRM matrices
library(calibrator)
grmMat <- ReadGRMBin("/projects/geneva/geneva_sata/caitlin/grm_xchr/hapmap3_r3_b37_fwd.consensus.qc.poly_Ilmn1M")
class(grmMat); length(grmMat); names(grmMat) # list w 4 elements: diag, off, id, N

length(grmMat$off); length(grmMat$diag) # need to make this into a matrix

kinship_x <- diag(grmMat$diag)
dim(kinship_x)

kinship_x[upper.tri(kinship_x)] <- grmMat$off
kinship_x <- symmetrize(kinship_x)
isSymmetric(kinship_x) # TRUE
kinship_x[1:10,1:10] # appears to have worked

grmMat <- ReadGRMBin("/projects/geneva/geneva_sata/caitlin/grm_xchr/hapmap3_r3_b37_fwd.consensus.qc.poly_Ilmn1M_autos")
class(grmMat); length(grmMat); names(grmMat) # list w 4 elements: diag, off, id, N

length(grmMat$off); length(grmMat$diag) # need to make this into a matrix

kinship_autos <- diag(grmMat$diag)
dim(kinship_autos)

kinship_autos[upper.tri(kinship_autos)] <- grmMat$off
kinship_autos <- symmetrize(kinship_autos)
isSymmetric(kinship_autos) # TRUE
kinship_autos[1:10,1:10] # appears to have worked

# need to make a scanAnnot object
scan <- get(load("/projects/geneva/geneva_sata/HapMap/HapMap3_r3/sample_snp_annot/hapmap3_r3_b36.sample.annot.v6.RData"))
scan

# make a phenotype column, as simulated with GCTA
sims <- read.table("/projects/geneva/geneva_sata/caitlin/grm_xchr/hapmap3_r3_b37_fwd.consensus.qc.poly_Ilmn1M_causSNPs_maf2_effSz3.phen", header=FALSE,as.is=TRUE)
dim(sims) # 1397 102
head(sims) # ok, so 100 simulated phenotypes, first 2 cols are fam and individ ids

scanMerged <- merge(pData(scan),sims[,c("V2","V3")],by.x="coriell.id",by.y="V2",all.x=TRUE)

scan <- scan[order(scan$coriell.id),]
all(scanMerged$coriell.id==scan$coriell.id) # TRUE

scan$pheno_sim <- scanMerged$V3
scan <- scan[order(scan$scanID),]

# need row and col names of the kinship matrices to be the scanIDs
idMap <- grmMat$id
idMap$order <- 1:nrow(idMap)
idMap <- merge(idMap,pData(scan)[,c("coriell.id","scanID")],by.x="V2",by.y="coriell.id")
idMap <- idMap[order(idMap$order),]
rownames(kinship_autos) <- colnames(kinship_autos) <- idMap$scanID
rownames(kinship_x) <- colnames(kinship_x) <- idMap$scanID

# merge in PC1-5
pcRes <- get(load("/projects/geneva/geneva_sata/caitlin/grm_xchr/pca_results.RData"))
names(pcRes)
evs <- pcRes$eigenvect
dim(evs) # 1397 32

idMap <- idMap[order(idMap$scanID),]
scan <- scan[order(idMap$order),]
all(pcRes$sample.id==scan$coriell.id) # TRUE
scan$PC1 <- evs[,1]
scan$PC2 <- evs[,2]
scan$PC3 <- evs[,3]
scan$PC4 <- evs[,4]
scan$PC5 <- evs[,5]

# now they need to be in the same order
idMap <- idMap[order(idMap$scanID),]
scan <- scan[order(scan$scanID),]
scan <- scan[order(idMap$order),]
head(pData(scan))
head(rownames(kinship_autos))

# use matt's code
source("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/mconomos/R/estVarComp.R")
covList <- list(kinship_autos,kinship_x)
names(covList) <- c("kinshipAuto","kinshipX")
varComp <- estVarComp(scan,covMatList=covList,"pheno_sim",covar.vec=c("PC1","PC2","PC3","PC4","PC5"))

varComp$varComp
#V_kinshipAuto    V_kinshipX           V_E 
#1.996201e-07  3.143519e+02  2.016340e+03 
# hmm, weird that it's the auto that is zero. wasn't expecting that.

# switch the order of the arguments and see if i still get the same thing
covList <- list(kinship_x,kinship_autos)
names(covList) <- c("kinshipX","kinshipAuto")
varComp2 <- estVarComp(scan,covMatList=covList,"pheno_sim",covar.vec=c("PC1","PC2","PC3","PC4","PC5"))

varComp2$varComp
#V_kinshipX V_kinshipAuto           V_E 
#3.143519e+02  1.996188e-07  2.016340e+03 
# same. weird.



## 
# now need to call the MM function

# make genotype data object
data <- GdsGenotypeReader("/projects/geneva/geneva_sata/caitlin/grm_xchr/hapmap3_r3_b37_fwd.consensus.qc.poly_Ilmn1M.gds")
idMap <- idMap[order(idMap$order),]
scan <- scan[order(idMap$order),]
scan$scanID <- scan$coriell.id
genoData <- GenotypeData(data,scanAnnot=scan)

chr <- getChromosome(data)
table(chr)
which.max(chr==25) # 944676
source("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/mconomos/R/assocTestMixedModel_v7.R")

# need to change col/rownames back to coriell ids
cholSig <- varComp[["cholSigmaInv"]]
all(idMap$scanID==colnames(cholSig)) # TRUE
colnames(cholSig) <- rownames(cholSig) <- idMap$V2

mmRes <- assocTestMixedModel(genoData,snpStart=1,snpEnd=(which.max(chr==25)-1),cholSig,outcome="pheno_sim",covar.vec=c("PC1","PC2","PC3","PC4","PC5"))

png("/projects/geneva/geneva_sata/caitlin/grm_xchr/manh_withPCA_res_OR30.png")
manhattanPlot(mmRes$pval,mmRes$chr)
dev.off()

# look at the pvals for the causal snps
caus <- read.table("/projects/geneva/geneva_sata/caitlin/grm_xchr/causSNPs_maf2_effSz3.txt",as.is=TRUE)
head(caus)

mmRes[is.element(mmRes$snpID,caus$V1),]

# save these results
save(mmRes,file="/projects/geneva/geneva_sata/caitlin/grm_xchr/res_withPCA_OR30.RData")

rm(list=ls())


#####
# 7. Run the MLM not adj for X chr with OR 3.0 phenotype on X chr

# read in the GRM matrices
library(calibrator)

grmMat <- ReadGRMBin("/projects/geneva/geneva_sata/caitlin/grm_xchr/hapmap3_r3_b37_fwd.consensus.qc.poly_Ilmn1M_autos")
class(grmMat); length(grmMat); names(grmMat) # list w 4 elements: diag, off, id, N

length(grmMat$off); length(grmMat$diag) # need to make this into a matrix

kinship_autos <- diag(grmMat$diag)
dim(kinship_autos)

kinship_autos[upper.tri(kinship_autos)] <- grmMat$off
kinship_autos <- symmetrize(kinship_autos)
isSymmetric(kinship_autos) # TRUE
kinship_autos[1:10,1:10] # appears to have worked

# need to make a scanAnnot object
scan <- get(load("/projects/geneva/geneva_sata/HapMap/HapMap3_r3/sample_snp_annot/hapmap3_r3_b36.sample.annot.v6.RData"))
scan

# make a phenotype column, as simulated with GCTA
sims <- read.table("/projects/geneva/geneva_sata/caitlin/grm_xchr/hapmap3_r3_b37_fwd.consensus.qc.poly_Ilmn1M_causSNPs_maf2_effSz3.phen", header=FALSE,as.is=TRUE)
dim(sims) # 1397 102
head(sims) # ok, so 100 simulated phenotypes, first 2 cols are fam and individ ids

scanMerged <- merge(pData(scan),sims[,c("V2","V3")],by.x="coriell.id",by.y="V2",all.x=TRUE)

scan <- scan[order(scan$coriell.id),]
all(scanMerged$coriell.id==scan$coriell.id) # TRUE

scan$pheno_sim <- scanMerged$V3
scan <- scan[order(scan$scanID),]

# need row and col names of the kinship matrices to be the scanIDs
idMap <- grmMat$id
idMap$order <- 1:nrow(idMap)
idMap <- merge(idMap,pData(scan)[,c("coriell.id","scanID")],by.x="V2",by.y="coriell.id")
idMap <- idMap[order(idMap$order),]
rownames(kinship_autos) <- colnames(kinship_autos) <- idMap$scanID

# merge in PC1-5
pcRes <- get(load("/projects/geneva/geneva_sata/caitlin/grm_xchr/pca_results.RData"))
names(pcRes)
evs <- pcRes$eigenvect
dim(evs) # 1397 32

idMap <- idMap[order(idMap$scanID),]
scan <- scan[order(idMap$order),]
all(pcRes$sample.id==scan$coriell.id) # TRUE
scan$PC1 <- evs[,1]
scan$PC2 <- evs[,2]
scan$PC3 <- evs[,3]
scan$PC4 <- evs[,4]
scan$PC5 <- evs[,5]

# now they need to be in the same order
idMap <- idMap[order(idMap$scanID),]
scan <- scan[order(scan$scanID),]
scan <- scan[order(idMap$order),]
head(pData(scan))
head(rownames(kinship_autos))

# use matt's code
source("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/mconomos/R/estVarComp.R")
covList <- list(kinship_autos)
names(covList) <- c("kinshipAuto")
varComp <- estVarComp(scan,covMatList=covList,"pheno_sim",covar.vec=c("PC1","PC2","PC3","PC4","PC5"))

varComp$varComp
#V_kinshipAuto           V_E 
#280.9011     2061.6981 

## 
# now need to call the MM function

# make genotype data object
data <- GdsGenotypeReader("/projects/geneva/geneva_sata/caitlin/grm_xchr/hapmap3_r3_b37_fwd.consensus.qc.poly_Ilmn1M.gds")
idMap <- idMap[order(idMap$order),]
scan <- scan[order(idMap$order),]
scan$scanID <- scan$coriell.id
genoData <- GenotypeData(data,scanAnnot=scan)

chr <- getChromosome(data)
table(chr)
which.max(chr==25) # 944676
source("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/mconomos/R/assocTestMixedModel_v7.R")

# need to change col/rownames back to coriell ids
cholSig <- varComp[["cholSigmaInv"]]
all(idMap$scanID==colnames(cholSig)) # TRUE
colnames(cholSig) <- rownames(cholSig) <- idMap$V2

mmRes <- assocTestMixedModel(genoData,snpStart=1,snpEnd=(which.max(chr==25)-1),cholSig,outcome="pheno_sim",covar.vec=c("PC1","PC2","PC3","PC4","PC5"))

png("/projects/geneva/geneva_sata/caitlin/grm_xchr/manh_withPCA_res_OR30_unadjX.png")
manhattanPlot(mmRes$pval,mmRes$chr)
dev.off()
# umm, looks like i have more power when not adj for X chr relatedness?

# look at the pvals for the causal snps
caus <- read.table("/projects/geneva/geneva_sata/caitlin/grm_xchr/causSNPs_maf2_effSz3.txt",as.is=TRUE)
head(caus)

mmRes[is.element(mmRes$snpID,caus$V1),]
summary(mmRes$pval[is.element(mmRes$snpID,caus$V1)])
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#0.0000000 0.0000677 0.0060340 0.0712600 0.0135100 0.6726000 
# two SNPs (rs5944680, rs5944697) with pval xxe-11, all the others range around xxe-03 or so

# save these results
save(mmRes,file="/projects/geneva/geneva_sata/caitlin/grm_xchr/res_withPCA_OR30_unadjX.RData")

rm(list=ls())


#####
# 8. Read in results from running on 100 phenotype simulations

setwd("/projects/geneva/geneva_sata/caitlin/grm_xchr/")

res_inclX <- get(load("varCompEst_maf2_effSz3.RData"))
res_exclX <- get(load("varCompEst_maf2_effSz3_exclX.RData"))
dim(res_inclX); dim(res_exclX) # both 100 4
# of course, the exclX doesn't have values for the estimated var comp on x chr

apply(res_inclX,2,summary)
apply(res_exclX,2,summary)

library(ggplot2)
resFull <- rbind(res_inclX,res_exclX)
resFull$model <- "inclX"
resFull$model[is.na(resFull$varX)] <- "exclX"

pdf("auto_var_comp_est.pdf")
ggplot(resFull,aes(x=varAuto,fill=model)) +
  geom_histogram(alpha=0.2)
dev.off()

# calculate type I error and power for the two models
pval_inclX <- get(load("pvalues_maf2_effSz3.RData"))
pval_exclX <- get(load("pvalues_maf2_effSz3_exclX.RData"))

dim(pval_inclX); dim(pval_exclX) # 944675 403 for both

# make manhattan plots for one of the simulated phenotype results
library(GWASTools)
png("manh_sim34_inclExclX_maf2_effSz3.png",width=720)
par(mfrow=c(2,1))
manhattanPlot(pval_inclX[,"pval_sim41"],pval_inclX[,"CHR"],main="Results when adjusing for X Chr effects")
manhattanPlot(pval_exclX[,"pval_sim41"],pval_exclX[,"CHR"],main="Results when not adjusting for X Chr effects")
dev.off()

# get the snps that are associated, with effect size 3
caus <- read.table("/projects/geneva/geneva_sata/caitlin/grm_xchr/causSNPs_maf2_effSz3.txt",header=FALSE,as.is=TRUE)
dim(caus); head(caus) # 13 2

pval_inclX$caus <- FALSE
pval_inclX$caus[is.element(pval_inclX$SNP,caus$V1)] <- TRUE

pval_exclX$caus <- FALSE
pval_exclX$caus[is.element(pval_exclX$SNP,caus$V1)] <- TRUE

table(pval_inclX$caus,exclude=NULL); table(pval_exclX$caus,exclude=NULL) # 13 TRUE for both

power <- data.frame(alpha=c(0.05,0.005,5e-04,5e-05,5e-06,5e-07,5e-08))
power$inclX <- NA
power$exclX <- NA
##### POWER, including X adjustment #####
trueHits <- pval_inclX[pval_inclX$caus,c(4:103)]
dim(trueHits) # 13 100; so these are the pvalues for all causal SNPs over all simulations
allPs <- as.vector(as.matrix(trueHits))
for(i in 1:nrow(power)){
  power$inclX[i] <- sum(allPs<power$alpha[i])/length(allPs)
}

##### POWER, excluding X adjustment #####
trueHits <- pval_exclX[pval_exclX$caus,c(4:103)]
dim(trueHits) # 13 100; so these are the pvalues for all causal SNPs over all simulations
allPs <- as.vector(as.matrix(trueHits))
for(i in 1:nrow(power)){
  power$exclX[i] <- sum(allPs<power$alpha[i])/length(allPs)
}

library(xtable)
xtable(power,digits=6)


##### now examine typeI error
# we have many more SNPs to test this

typeIerr <- power
typeIerr$inclX <- NA
typeIerr$exclX <- NA


falseHits <- pval_inclX[!pval_inclX$caus,c(4:103)]
allPs <- as.vector(as.matrix(falseHits))
for(i in 1:nrow(typeIerr)){
  typeIerr$inclX[i] <- sum(allPs<typeIerr$alpha[i])/length(allPs)
}

falseHits <- pval_exclX[!pval_exclX$caus,c(4:103)]
allPs <- as.vector(as.matrix(falseHits))
for(i in 1:nrow(typeIerr)){
  typeIerr$exclX[i] <- sum(allPs<typeIerr$alpha[i])/length(allPs)
}

xtable(typeIerr,digits=6)