

#####
## Contents:
# 1.
# 2.
# 3. Create x kinship with 10K pedigrees
# 4. Create auto kinship with 10K pedigrees
# 5. Create x kinship with 500 pedigrees
# 6. Create auto kinship with 500 pedigrees
# 7. Create x kinship with 1000 pedigrees, 8 person pedigree
# 8. Create auto kinship with 1000 pedigrees, 8 person pedigree
# 9. Create x kinship with 1000 pedigrees, 8 person pedigree w/ X kc <= auto KC
# 10. Create auto kinship with 1000 pedigrees, 8 person pedigree w/ X kc <= auto KC
# 11. Create unrelated kinship w 8000 samples



#####
# 1. 

library(MASS)
library(corpcor)
source("sim_phenotype.R")
source("estVarComp.R")
source("allele_drop_functions.R")
source("assocTestMixedModel_v7.R")

n <- 1600

TEMP=c(2,1,2,1,1,1,2,2,1,2,1,2,1,2,1,2)
SEX <- rep(NA,length(TEMP))
SEX[TEMP==1]="M"
SEX[TEMP==2]="F"

sex <- rep(SEX,(n/16))

# get x chromosome kinship matrix
# get autosomal kinship matrix

# kinship matrix, autosomes
library(kinship2)
tmp <- read.table("pedigree_16individs.txt")
tmp[tmp[,5]==0,"V5"] <- 2
ped <- pedigree(tmp[,2],dadid=tmp[,3],momid=tmp[,4],sex=tmp[,5])
kin <- kinship(ped)

# make it block diagonal
kinFull <- bdiag(kin,kin,kin,kin,kin,
                 kin,kin,kin,kin,kin)
kinFull <- bdiag(kinFull,kinFull,kinFull,kinFull,kinFull,
                 kinFull,kinFull,kinFull,kinFull,kinFull)
dim(kinFull) # 1600 1600, as expected
kinFull <- as.matrix(kinFull)
save(kinFull,file="100Peds_autoKinship.RData")

### save 2*KC and try the analyses with that 
kinFull <- 2*kinFull
kinFull[1:10,1:10]
save(kinFull,file="100Peds_2autoKinship.RData")


# kinship matrix, x chromosome
tmp <- read.table("pedigree_16individs_output")
trueKinX <- matrix(NA,nrow=16, ncol=16)
trueKinX[lower.tri(trueKinX,diag=TRUE)] <- tmp[,"V4"]
trueKinXn <- trueKinX
trueKinX <- data.frame(trueKinX)
trueKinX$SEX <- SEX

#trueKinXn[trueKinX$SEX=="F",trueKinX$SEX=="F"] <- 2*trueKinXn[trueKinX$SEX=="F",trueKinX$SEX=="F"]

# is this correct???
#trueKinXn[trueKinX$SEX=="M",trueKinX$SEX=="F"] <- sqrt(2)*trueKinXn[trueKinX$SEX=="M",trueKinX$SEX=="F"]
#trueKinXn[trueKinX$SEX=="F",trueKinX$SEX=="M"] <- sqrt(2)*trueKinXn[trueKinX$SEX=="F",trueKinX$SEX=="M"]

#diag(trueKinXn) <- 1
diag(trueKinXn)[SEX=="M"] <- 1
diag(trueKinXn)[SEX=="F"] <- 0.5
trueKinXn[upper.tri(trueKinXn,diag=FALSE)] <- t(trueKinXn)[upper.tri(trueKinXn,diag=FALSE)]

kinFullX <- bdiag(trueKinXn,trueKinXn,trueKinXn,trueKinXn,trueKinXn,
                  trueKinXn,trueKinXn,trueKinXn,trueKinXn,trueKinXn)
kinFullX <- bdiag(kinFullX,kinFullX,kinFullX,kinFullX,kinFullX,
                  kinFullX,kinFullX,kinFullX,kinFullX,kinFullX)
dim(kinFullX) # 1600 1600
kinFullX <- as.matrix(kinFullX)
save(kinFullX,file="100Peds_xKinship.RData")

## save 2*KC for some analyses
kinFullX <- 2*kinFullX
kinFullX[1:10,1:10]
save(kinFullX,file="100Peds_2xKinship.RData")

rm(list=ls())


######
# 2.

# add unrelateds to these kinship matrices
kinFull <- get(load("100Peds_autoKinship.RData"))
# add in 500 unrelateds to the end of this matrix
kinFullUnr <- matrix(0,nrow=2100,ncol=2100)
diag(kinFullUnr) <- 0.5
kinFullUnr[1:1600,1:1600] <- kinFull
save(kinFullUnr,file="100Peds_500unr_autoKinship.RData")

kinFullUnr <- 2*kinFullUnr
kinFullUnr[1:10,1:10]
save(kinFullUnr,file="100Peds_500unr_2autoKinship.RData")

kinFullX <- get(load("100Peds_xKinship.RData"))
kinFullXUnr <- matrix(0,nrow=2100,ncol=2100)
# make the unrelated samples to be first 250: females, last 250: males
diag(kinFullXUnr)[1601:1850] <- 0.5 # females
diag(kinFullXUnr)[1851:2100] <- 1 # males
kinFullXUnr[1:1600,1:1600] <- kinFullX
save(kinFullXUnr,file="100Peds_500unr_xKinship.RData")

kinFullXUnr <- 2*kinFullXUnr
kinFullXUnr[1:10,1:10]
save(kinFullXUnr,file="100Peds_500unr_2xKinships.RData")

rm(list=ls())


#####
# 3. Create x kinship with 10K pedigrees

n <- 16000

TEMP=c(2,1,2,1,1,1,2,2,1,2,1,2,1,2,1,2)
SEX <- rep(NA,length(TEMP))
SEX[TEMP==1]="M"
SEX[TEMP==2]="F"

sex <- rep(SEX,(n/16))
length(sex) # 160000
table(sex) # 80K fem, 80K males

tmp <- read.table("pedigree_16individs_output")
trueKinX <- matrix(NA,nrow=16, ncol=16)
trueKinX[lower.tri(trueKinX,diag=TRUE)] <- tmp[,"V4"]
trueKinXn <- trueKinX
trueKinX <- data.frame(trueKinX)
trueKinX$SEX <- SEX

diag(trueKinXn)[SEX=="M"] <- 1
diag(trueKinXn)[SEX=="F"] <- 0.5
trueKinXn[upper.tri(trueKinXn,diag=FALSE)] <- t(trueKinXn)[upper.tri(trueKinXn,diag=FALSE)]

kinFullX <- bdiag(trueKinXn,trueKinXn,trueKinXn,trueKinXn,trueKinXn,
                  trueKinXn,trueKinXn,trueKinXn,trueKinXn,trueKinXn)
kinFullX <- bdiag(kinFullX,kinFullX,kinFullX,kinFullX,kinFullX,
                  kinFullX,kinFullX,kinFullX,kinFullX,kinFullX)
kinFullX <- bdiag(kinFullX,kinFullX,kinFullX,kinFullX,kinFullX,
                  kinFullX,kinFullX,kinFullX,kinFullX,kinFullX)
dim(kinFullX) # 16000 16000; good

save(kinFullX,file="1KPeds_xKinship.RData")
rm(list=ls())


#####
# 4. Create auto kinship with 10K pedigrees

library(kinship2)
tmp <- read.table("pedigree_16individs.txt")
tmp[tmp[,5]==0,"V5"] <- 2
ped <- pedigree(tmp[,2],dadid=tmp[,3],momid=tmp[,4],sex=tmp[,5])
kin <- kinship(ped)

kinFull <- bdiag(kin,kin,kin,kin,kin,
                 kin,kin,kin,kin,kin)
kinFull <- bdiag(kinFull,kinFull,kinFull,kinFull,kinFull,
                 kinFull,kinFull,kinFull,kinFull,kinFull)
kinFull <- bdiag(kinFull,kinFull,kinFull,kinFull,kinFull,
                 kinFull,kinFull,kinFull,kinFull,kinFull)

dim(kinFull) # 16000 16000, as expected
save(kinFull,file="1KPeds_autoKinship.RData")

rm(list=ls())


#####
# 5. Create x kinship with 500 pedigrees

library(kinship2)
n <- 8000

TEMP=c(2,1,2,1,1,1,2,2,1,2,1,2,1,2,1,2)
SEX <- rep(NA,length(TEMP))
SEX[TEMP==1]="M"
SEX[TEMP==2]="F"

sex <- rep(SEX,(n/16))
length(sex) # 8000
table(sex) # 4k fem, 4K males

tmp <- read.table("pedigree_16individs_output")
trueKinX <- matrix(NA,nrow=16, ncol=16)
trueKinX[lower.tri(trueKinX,diag=TRUE)] <- tmp[,"V4"]
trueKinXn <- trueKinX
trueKinX <- data.frame(trueKinX)
trueKinX$SEX <- SEX

diag(trueKinXn)[SEX=="M"] <- 1
diag(trueKinXn)[SEX=="F"] <- 0.5
trueKinXn[upper.tri(trueKinXn,diag=FALSE)] <- t(trueKinXn)[upper.tri(trueKinXn,diag=FALSE)]

kinFullX <- bdiag(trueKinXn,trueKinXn,trueKinXn,trueKinXn,trueKinXn,
                  trueKinXn,trueKinXn,trueKinXn,trueKinXn,trueKinXn)
kinFullX <- bdiag(kinFullX,kinFullX,kinFullX,kinFullX,kinFullX,
                  kinFullX,kinFullX,kinFullX,kinFullX,kinFullX)
kinFullX <- bdiag(kinFullX,kinFullX,kinFullX,kinFullX,kinFullX)
dim(kinFullX) # 8000 8000; good

save(kinFullX,file="500Peds_xKinship.RData")
rm(list=ls())


#####
# 6. Create auto kinship with 500 pedigrees

library(kinship2)
tmp <- read.table("pedigree_16individs.txt")
tmp[tmp[,5]==0,"V5"] <- 2
ped <- pedigree(tmp[,2],dadid=tmp[,3],momid=tmp[,4],sex=tmp[,5])
kin <- kinship(ped)

kinFull <- bdiag(kin,kin,kin,kin,kin,
                 kin,kin,kin,kin,kin)
kinFull <- bdiag(kinFull,kinFull,kinFull,kinFull,kinFull,
                 kinFull,kinFull,kinFull,kinFull,kinFull)
kinFull <- bdiag(kinFull,kinFull,kinFull,kinFull,kinFull)

dim(kinFull) # 8000 8000, as expected
save(kinFull,file="500Peds_autoKinship.RData")

rm(list=ls())


#####
# 7. Create x kinship with 1000 pedigrees, 8 person pedigree

n <- 8000

SEX <- c("F","M","M","M","F","M","M","M")
sex <- rep(SEX,(n/8))
length(sex) # 8000
table(sex) # 2k fem, 6K males

sex_ped <- SEX
sex_ped[sex_ped=="F"] <- 2
sex_ped[sex_ped=="M"] <- 1
ped <- cbind(rep(1,8),1:8,c(0,0,0,2,2,3,3,3),c(0,0,0,1,1,5,5,5),sex_ped)
write.table(ped,file="KinInbcoefX/input_8ped.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

listf <- ped[,c(1,2)]
write.table(listf,file="KinInbcoefX/listfile_8ped.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
# called KinInbcoefX on the servers

###

tmp <- read.table("KinInbcoefX/outfile_8ped")
trueKinX <- matrix(NA,nrow=8, ncol=8)
trueKinX[lower.tri(trueKinX,diag=TRUE)] <- tmp[,"V4"]
trueKinXn <- trueKinX
trueKinX <- data.frame(trueKinX)
trueKinX$SEX <- SEX

diag(trueKinXn)[SEX=="M"] <- 1
diag(trueKinXn)[SEX=="F"] <- 0.5
trueKinXn[upper.tri(trueKinXn,diag=FALSE)] <- t(trueKinXn)[upper.tri(trueKinXn,diag=FALSE)]

kinFullX <- bdiag(trueKinXn,trueKinXn,trueKinXn,trueKinXn,trueKinXn,
                  trueKinXn,trueKinXn,trueKinXn,trueKinXn,trueKinXn)
kinFullX <- bdiag(kinFullX,kinFullX,kinFullX,kinFullX,kinFullX,
                  kinFullX,kinFullX,kinFullX,kinFullX,kinFullX)
kinFullX <- bdiag(kinFullX,kinFullX,kinFullX,kinFullX,kinFullX,
                  kinFullX,kinFullX,kinFullX,kinFullX,kinFullX)
dim(kinFullX) # 8000 8000; good

save(kinFullX,file="1000Peds_8ped_xKinship.RData")

rm(list=ls())


#####
# 8. Create auto kinship with 1000 pedigrees, 8 person pedigree

library(kinship2)
tmp <- read.table("KinInbcoefX/input_8ped.txt")

ped <- pedigree(tmp[,2],dadid=tmp[,3],momid=tmp[,4],sex=tmp[,5])
kin <- kinship(ped)

kinFull <- bdiag(kin,kin,kin,kin,kin,
                 kin,kin,kin,kin,kin)
kinFull <- bdiag(kinFull,kinFull,kinFull,kinFull,kinFull,
                 kinFull,kinFull,kinFull,kinFull,kinFull)
kinFull <- bdiag(kinFull,kinFull,kinFull,kinFull,kinFull,
                 kinFull,kinFull,kinFull,kinFull,kinFull)

dim(kinFull) # 8000 8000, as expected
save(kinFull,file="1000Peds_8ped_autoKinship.RData")

pdf("pedigree_8ped.pdf")
plot.pedigree(ped,cex=2,lwd=2)
dev.off()

rm(list=ls())


#####
# 9. Create x kinship with 1000 pedigrees, 8 person pedigree w/ X kc <= auto KC

n <- 8000

SEX <- c("F","M","F","F","M","M","M","F")
sex <- rep(SEX,(n/8))
length(sex) # 8000
table(sex) # 4k fem, 4K males

sex_ped <- SEX
sex_ped[sex_ped=="F"] <- 2
sex_ped[sex_ped=="M"] <- 1
ped <- cbind(rep(1,8),1:8,c(0,0,0,0,2,2,6,6),c(0,0,0,0,1,1,4,4),sex_ped)
write.table(ped,file="KinInbcoefX/input_8ped_fem.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

listf <- ped[,c(1,2)]
write.table(listf,file="KinInbcoefX/listfile_8ped_fem.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
# called KinInbcoefX on the servers

###

tmp <- read.table("KinInbcoefX/outfile_8ped_fem")
trueKinX <- matrix(NA,nrow=8, ncol=8)
trueKinX[lower.tri(trueKinX,diag=TRUE)] <- tmp[,"V4"]
trueKinXn <- trueKinX
trueKinX <- data.frame(trueKinX)
trueKinX$SEX <- SEX

diag(trueKinXn)[SEX=="M"] <- 1
diag(trueKinXn)[SEX=="F"] <- 0.5
trueKinXn[upper.tri(trueKinXn,diag=FALSE)] <- t(trueKinXn)[upper.tri(trueKinXn,diag=FALSE)]

kinFullX <- bdiag(trueKinXn,trueKinXn,trueKinXn,trueKinXn,trueKinXn,
                  trueKinXn,trueKinXn,trueKinXn,trueKinXn,trueKinXn)
kinFullX <- bdiag(kinFullX,kinFullX,kinFullX,kinFullX,kinFullX,
                  kinFullX,kinFullX,kinFullX,kinFullX,kinFullX)
kinFullX <- bdiag(kinFullX,kinFullX,kinFullX,kinFullX,kinFullX,
                  kinFullX,kinFullX,kinFullX,kinFullX,kinFullX)
dim(kinFullX) # 8000 8000; good

save(kinFullX,file="1000Peds_8ped_fem_xKinship.RData")

rm(list=ls())


#####
# 10. Create auto kinship with 1000 pedigrees, 8 person pedigree w/ X kc <= auto KC

library(kinship2)
tmp <- read.table("input_8ped_fem.txt")

ped <- pedigree(tmp[,2],dadid=tmp[,3],momid=tmp[,4],sex=tmp[,5])
kin <- kinship(ped)

kinFull <- bdiag(kin,kin,kin,kin,kin,
                 kin,kin,kin,kin,kin)
kinFull <- bdiag(kinFull,kinFull,kinFull,kinFull,kinFull,
                 kinFull,kinFull,kinFull,kinFull,kinFull)
kinFull <- bdiag(kinFull,kinFull,kinFull,kinFull,kinFull,
                 kinFull,kinFull,kinFull,kinFull,kinFull)

dim(kinFull) # 8000 8000, as expected
save(kinFull,file="1000Peds_8ped_fem_autoKinship.RData")

pdf("pedigree_8ped_fem.pdf")
plot.pedigree(ped,cex=2,lwd=2) # didn't plot individ 3 since just a founder and nothing else
dev.off()

rm(list=ls())


#####
# 11. Create unrelated kinship w 8000 samples

SEX <- c("F","F","F","F","M","M","M","M")

xDiag <- rep(SEX,1000)
xDiag[xDiag=="F"]	<- 0.5
xDiag[xDiag=="M"] <- 1
xDiag <- as.numeric(xDiag)
kinX <- diag(xDiag)
kinAuto <- 0.5*diag(1000)

save(kinX,file="1000Peds_8000unrel_xKinship.RData")
save(kinAuto,file="1000Peds_8000unrel_autoKinship.RData")

rm(list=ls())