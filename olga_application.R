
## 
# data application of the mlm_x method to the OLGA data

## Contents:
# 1. Get theoretical KC for X chromosome in OLGA samples
 # 1a. Organize .fam file
 # 1b. Call KinInbcoefX
 # 1c. Create matrix from output
# 2. Call PC-Air-X using OLGA samples
# 3. Get pruned set of X chr SNPs
# 4. Run PC-AiR on pruned set of SNPs; smaller sample set
# 5. X chr relatedness estimates
# 6. Prune auto set of SNPs, incl and excl X chr
# 7. PCA on the autosomes+X incl proj to rel
# 8. PCA on the autosomes+X, unrel only
# 9. PCA with auto SNPs
# 10. PC-Relate on X chr SNPs, corr for auto+X PC1-5
# 11. PC-Relate on X chr SNPs, corr for X chr PC1-2
# 12. Parse KC results
# 13. PC-AiR on X chr SNPs adj for X chr KC (which is adj for auto+x ancestry)
# 14. Parse PC-AiR results on X chr after adj for X chr KC (adj for auto+x PC 1-5)
# 15. Run PC-AiR on 3500 autosomal SNPs
# 16. PC-Relate on X chr SNP, corr for auto PC1-5
# 17. Parse PC-AiR results on chr 20 SNPs
# 18. PC-AiR on X chr SNPs adj for X chr KC (which is adj for auto+x ancestry)
# 19. Parse KC results, including #16. results too
# 20. Parse PC-AiR results on X chr SNPs adj for X chr KC (which is adj for auto ancestry)
# 21. PCA-SNP correlation plots
# 22. All PC-Relate runs again, using pruned X chr SNPs
# 23. Look at X-KC for some autosomal relationships, using KC from pruned X chr SNPs
# 24. Compare KC X using all SNPs vs pruned SNPs
# 25. PC-AiR on X chr SNPs adj for X chr KC (which is adj for auto ancestry)
# 26. PC-AiR on X chr SNPs adj for X chr pruned KC (which is adj for auto ancestry)
# 27. Parse results for 25., 26.
# 28. Do eigen() on output from 22. (KC_X matrix adj for auto PC 1-5)
# 29. Do eigen() on KC_X matrix adj for X chr PC 1-2, KC_auto matrix adj for auto PC 1-5, KC_chr20 matrix adj for auto PC 1-5
# 30. Parse chr 20 KC results
# 31. Parse results from #29.
# 32. Rerun X chr PCA
# 33. Rerun pcrelate on xchr, adj for x chr PCs from 32.

# 34. Est VC for LABA2 BCC trait, analysis 316987
# 35. Look at var comp estimates for the 4 models, analysis 316987
# 36. Assoc test for LABA2 BCC trait, analysis 316987
# 37. Parse assoc test results for X chr
# 38. Run assoc test on autosomes now
# 39. Summary statistics for trait LABA2 RBC

# 40. Calculate X chr EVs by gengrp6.strat
# 41. Parse X EVs by gengrp6.strat results
# 42. Est VC for BMI trait, analysis 298888
# 43. Assoc test for BMI trait, analysis 298888
# 44. Look at corr of X PC 1-2 with BMI

# 45. Run ADMIXTURE on the X chromosome
# 46. Plot X chr PC by ADMIXTURE proportions
# 47. Run ADMIXTURE by chromosome
# 48. CAnD on the x vs auto proportions

# 49. Run assoc with LABA2 BCC trait, excluding all PCs

# 50. CAnD on autosomes using local ancestry estimates
# 51. Summaries for CAnD; non param CAnD

# 52. New PO kc estimates, auto vs x chr
# 53. Recalculate lambda on X chr excl Xq28 region
# 54. Redo plots with local ancestry estimates
# 55. Perform CAnD and nonP CAnD and pooled t-test on local ancestry estimates
# 56. Parse CAnD results
# 57. Look at PR chr 2 results in depth
# 58. Look at PR chr 2 results in depth, EUR ancestry

# 59. Re-run ADMIXTURE on the X chromosome
# 60. CAnD on the x vs auto proportions
# 61. Get k0, k1, k2 for X KC F-F pairs
# 62. Make some plots of KC on X chr for MLM-X paper
# 63. Calculate genomic inf across autos only, compare pvalues

# 64. Rerun PC-AiR, PC-Relate iterations w 3rd deg relatives
# 65. Manh of MLM-X, simple MLM X chr results in one png file
# 66. Run pcairPartition() incl & excl the autosomal KING divMat
# 67. Run PC-Relate on set of 3600 SNPs, not on X chr
# 68. Run regression to check corr of auto SNPs with X PCs
# 69. Parse chr 19 PC-Relate results
# 70. Plot distribution of X chr markers across the chr
# 71. QQ plots for simple MLM and MLM-X on RBC trait, with lambda_GC in legend
# 72. Run regression to check corr of X chr SNPs with auto PCs
# 73. Make barplots of var comp estimates w and w/o X chr KC for RBC
# 74. Compare chr 19 results to X chr KC FF pairs
# 75. Rerun CAnD with new p-value combining method



#####
# 1. Get theoretical KC for X chromosome in OLGA samples

## 1a. Organize .fam file
# from command line:
# cp /projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/plink/subjects/OLGA_TOP_subject_level_filtered_update.fam /projects/geneva/geneva_sata/caitlin/mlm_x/olga_application
# cut -f1,2,3,4,5 -d ' ' OLGA_TOP_subject_level_filtered_update.fam > OLGA_TOP_subject_level_filtered_update_sm.fam

# fam id must be consecutive integers from 1-100 (max is 100)
# families must be listed together
setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")
fam <- read.table("olga_application/OLGA_TOP_subject_level_filtered_update_sm.fam",header=FALSE,as.is=TRUE,sep=" ")
dim(fam) # 13204 5
colnames(fam) <- c("family","individ","father","mother","sex")

length(unique(fam$family)) # 11605; eek!
# remove the singletons?
table(table(fam$family))
#    1     2     3     4     5     6 
#10508   693   333    52    11     8 

t <- table(fam$family)
length(names(t)[t==1]) # 10508; these are the ones we want to exclude
singl <- names(t)[t==1]
fam <- fam[!is.element(fam$family,singl),]
dim(fam) # 2696 5; as expected
fam <- fam[order(as.integer(fam$family)),] # need to make the family ids integers 
# got warning for NAs, no worries
length(unique(fam$family)) # 1097

# make family id dictionary
# some of the integers 1-1097 are already used up, so fill in the blanks with the other ones
fam_id_key <- data.frame("orig_famID"=unique(fam$family),stringsAsFactors=FALSE)
fam_id_key$kinIbd_famID <- fam_id_key$orig_famID
toChg <- as.integer(fam_id_key$orig_famID)>1097
toChg[is.na(toChg)] <- TRUE
fam_id_key[toChg,]
# ok, now need to figure out the gaps to fill these with
usedId <- fam_id_key$orig_famID[!toChg]
allIds <- 1:1097
notUsedId <- allIds[!is.element(allIds,as.integer(usedId))]
length(notUsedId); sum(toChg) # both 412! great!

fam_id_key$kinIbd_famID[toChg] <- notUsedId

for(i in 1:nrow(fam_id_key)){
  if(fam_id_key$orig_famID[i]==fam_id_key$kinIbd_famID[i]){ next }
  fam$family[is.element(fam$family,fam_id_key$orig_famID[i])] <- fam_id_key$kinIbd_famID[i]
}

fam <- fam[order(as.integer(fam$family)),]

# individ must also be an integer -- just use integers 1-nrows
individIds <- unique(c(fam$individ,fam$father,fam$mother))
zeroInd <- which(individIds==0)
individIds <- individIds[-zeroInd] # don't include the 0s in there
for(i in 1:length(individIds)){
  fam$individ[fam$individ==individIds[i]] <- i
  fam$father[fam$father==individIds[i]] <- i
  fam$mother[fam$mother==individIds[i]] <- i
}

individ_id_key <- cbind("orig_individID"=individIds,"kinIbd_individID"=1:length(individIds))

## read in scan annot to see which samples are genotyped
library(GWASTools)
scan <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_scanAnnot_v54_AMS.RData"))
head(pData(scan)) # it's col subjectID
sum(is.element(individ_id_key[,"orig_individID"],scan$subjectID)) # 2696

individ_id_key <- data.frame(individ_id_key,stringsAsFactors=FALSE)
individ_id_key$genotyped <- FALSE
individ_id_key$genotyped[is.element(individ_id_key$orig_individID,scan$subjectID)] <- TRUE
table(individ_id_key$genotyped) # 2696 TRUE

# make a list of individs now for which we want kin coef, b/f adding fake samples to the fam file
list_file <- fam[,c("family","individ")]
head(list_file); dim(list_file) # 2696 2
sum(is.element(list_file$individ,individ_id_key$kinIbd_individID[individ_id_key$genotyped])) # 2696
all(is.element(list_file$individ,individ_id_key$kinIbd_individID[individ_id_key$genotyped])) # TRUE

# need to list the father/mothers on their own lines
for(i in 1:nrow(fam)){
  if(fam$mother[i]==fam$father[i] & fam$mother[i]==0){ next }
  if(!is.element(fam$mother[i],fam$individ)){
    newR <- c(fam$family[i],fam$mother[i],0,0,2)
    fam <- rbind(fam,newR)
  }
  if(!is.element(fam$father[i],fam$individ)){
    newR <- c(fam$family[i],fam$father[i],0,0,1)
    fam <- rbind(fam,newR)
  }
}
fam <- fam[order(as.integer(fam$family)),] # great

fam$family <- as.integer(fam$family)
fam$individ <- as.integer(fam$individ)
fam$mother <- as.integer(fam$mother)
fam$father <- as.integer(fam$father)

dim(fam) # 3627 5
length(unique(fam$family)) # 1097
# need to make 5 sub-files with 100 families

keys <- list(individ_id_key,fam_id_key)
names(keys) <- c("individ_id_key","fam_id_key")
save(keys,file="olga_application/famFile_keys.RData")

write.table(fam,file="olga_application/OLGA_TOP_subject_level_filtered_update_allFams.fam",
            sep=" ",quote=FALSE,col.names=FALSE,row.names=FALSE)
write.table(list_file,
            file="olga_application/OLGA_TOP_subject_level_filtered_update_allFams.list",sep=" ",quote=FALSE,
            col.names=FALSE,row.names=FALSE)

rm(list=ls())

## 1b. Call KinInbcoefX
# changed max families in KinInbcoefX.c to be 1200
# then called: g++ KinInbcoefX.c -o KinInbcoefX
# cd /projects/geneva/geneva_sata/caitlin/mlm_x/KinInbcoefX/
# ./KinInbcoefX ../olga_application/OLGA_TOP_subject_level_filtered_update_allFams.fam ../olga_application/OLGA_TOP_subject_level_filtered_update_allFams.list ../olga_application/allFams ../olga_application/error_allFams
# wc -l ../olga_application/allFams || output was 4930

## 1c. Create matrix from output
# kinship matrix, x chromosome
setwd("/projects/geneva/geneva_sata/caitlin/mlm_x/olga_application/")
tmp <- read.table("allFams")
head(tmp); dim(tmp) # 4930 4
# should only have output for the genotyped samples

keys <- get(load("famFile_keys.RData"))
ind_ids <- keys[[1]]

fam <- read.table("OLGA_TOP_subject_level_filtered_update_sm.fam",header=FALSE,as.is=TRUE,sep=" ")
ind_ids <- merge(ind_ids,fam[,c("V2","V5")],by.x="orig_individID",by.y="V2")
dim(ind_ids) # 2696 4; so these are the samples we want to keep

# sort the ind_ids$kinIbd_individID to be in increasing order
ind_ids$kinIbd_individID <- as.integer(ind_ids$kinIbd_individID)
ind_ids <- ind_ids[order(ind_ids$kinIbd_individID),]

trueKinX <- matrix(0, nrow=sum(ind_ids$genotyped), ncol=sum(ind_ids$genotyped))
colnames(trueKinX) <- rownames(trueKinX) <- ind_ids[ind_ids$genotyped,"kinIbd_individID"]
dim(trueKinX) # 2696 2696

diag(trueKinX)[ind_ids$V5==1] <- 1 # males
diag(trueKinX)[ind_ids$V5==2] <- 0.5 # females

table(ind_ids$V5)
#  1    2 
#980 1716 
table(diag(trueKinX))
# 0.5    1 
#1716  980 
# good! looking as expected

# remove the self-rows
# remove the rows that aren't included in the kin matrix
tmp <- tmp[tmp$V2!=tmp$V3,]
dim(tmp) # 2234 4

sum(!is.element(tmp$V3,colnames(trueKinX))) # 0
sum(!is.element(tmp$V2,colnames(trueKinX))) # 0
# great!

tmp <- tmp[order(tmp$V2),]
all(tmp$V2<tmp$V3) # TRUE; so this means it'll be an upper triangular matrix

library(GWASTools)
scan <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_scanAnnot_v54_AMS.RData"))
head(pData(scan)) # it's col subjectID

## make trueKinX large enough to hold all samples, pre-populate with 0, then add in the non-zero values
trueKinX <- matrix(0,nrow=sum(scan$subj.plink),ncol=sum(scan$subj.plink))
rownames(trueKinX) <- colnames(trueKinX) <- scan$subjectID[scan$subj.plink]

# order the ind_ids to be the same order as the scan ids
ind_ids_all <- merge(ind_ids,pData(scan)[,c("subj.plink","subjectID")],by.x="orig_individID",
                     by.y="subjectID",all.y=TRUE,sort=FALSE)
table(ind_ids_all$subj.plink)
ind_ids_all <- ind_ids_all[ind_ids_all$subj.plink,]

ind_ids_Sorted <- ind_ids_all[match(colnames(trueKinX),ind_ids_all$orig_individID),]
allequal(ind_ids_Sorted$orig_individID,colnames(trueKinX)) # TRUE!!!

ind_ids_Sorted$kinIbd_individID[is.na(ind_ids_Sorted$kinIbd_individID)] <- ind_ids_Sorted$orig_individID[is.na(ind_ids_Sorted$kinIbd_individID)]
colnames(trueKinX) <- ind_ids_Sorted$kinIbd_individID
rownames(trueKinX) <- ind_ids_Sorted$kinIbd_individID


for(i in 1:nrow(tmp)){
  rowInd <- which(rownames(trueKinX)==tmp$V2[i])
  colInd <- which(colnames(trueKinX)==tmp$V3[i])
  trueKinX[rowInd,colInd] <- tmp$V4[i]

  if(i %% 100 ==0){ print(i) }
}

# make the col/rownames the original subjectID col/rownames 
colnames(trueKinX) <- ind_ids_Sorted$orig_individID
rownames(trueKinX) <- ind_ids_Sorted$orig_individID

table(trueKinX)
#        0    0.0625   0.09375     0.125    0.1875      0.25     0.375       0.5 
#174343911         3         2        48         6       809        66       771 

all(trueKinX[lower.tri(trueKinX,diag=FALSE)]==0) # TRUE
table(trueKinX[lower.tri(trueKinX,diag=FALSE)]) # all zeroes

trueKinX[lower.tri(trueKinX,diag=FALSE)] <- t(trueKinX)[lower.tri(trueKinX,diag=FALSE)]
isSymmetric(trueKinX) # TRUE

# oh, need to add in the diag elements for the unrel samples
sex_info <- merge(ind_ids_Sorted,pData(scan)[,c("subjectID","subj.plink","sex")],all.x=TRUE,sort=FALSE,
                  by.x="orig_individID",by.y="subjectID")
sex_info <- sex_info[sex_info$subj.plink.y,]

allequal(sex_info$orig_individID,colnames(trueKinX)) # TRUE!

diag(trueKinX)[sex_info$sex=="M"] <- 1
diag(trueKinX)[sex_info$sex=="F"] <- 0.5

table(diag(trueKinX)) # 7749 0.5, 5455 1's
table(scan$sex[scan$subj.plink])
# yay!

save(trueKinX,file="xChr_kc_matrix.RData")

rm(list=ls())


#####
# 2. Call PC-Air-X using OLGA samples

# in email 23 jan 15, matt sent me pcair_X.R function and corresponding manual
setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")
source("pcair_X.R")
library(GWASTools)

scan <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_scanAnnot_v54_AMS.RData"))
head(pData(scan)) # it's col subjectID

snp <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_HCHS_Custom_15041502_B3_all37_v25_AMS.RData"))
head(pData(snp))
snpX <- snp$snpID[snp$chromosome==23]

kcX <- get(load("olga_application/xChr_kc_matrix.RData"))
# check that this is in the same order as the scanAnnot
kcX[1:10,1:10]; dim(kcX) # 13204 13204
allequal(colnames(kcX),rownames(kcX)) # TRUE
allequal(colnames(kcX),scan$subjectID[scan$subj.plink]) # TRUE!

# make genoData object
gdsobj <- GdsGenotypeReader("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/netCDF/subjects/OLGA_subj_geno.gds")
genoData <- GenotypeData(gdsobj,scanAnnot=scan[scan$subj.plink,],snpAnnot=snp)

unrel.set <- scan$scanID[scan$unrelated.pcair.deg4&scan$subj.plink]
scanIncl <- scan$scanID[is.na(scan$gengrp6.outliers)&scan$geno.cntl==0&scan$subj.plink]
pc <- pcair(genoData,unrel.set=unrel.set,Xchr=TRUE,snp.include=snpX,scan.include=scanIncl)
# note this inclues HM samples

pdf("tmp_pc.pdf")
plot(pc$vectors[,1],pc$vectors[,2],xlab="PC1",ylab="PC2")
dev.off()

scan_tmp <- scan
scan <- scan_tmp[scan_tmp$subj.plink,]

pdf("tmp_pc_col_v2.pdf")
plot(-pc$vectors[,1],pc$vectors[,2],xlab="PC1",ylab="PC2",type="n")
points(-pc$vectors[scan$race.cat=="Mexican",1],pc$vectors[scan$race.cat=="Mexican",2],pch=20,col="goldenrod")
points(-pc$vectors[scan$race.cat=="CentralAmerican",1],pc$vectors[scan$race.cat=="CentralAmerican",2],pch=20,col="red")
points(-pc$vectors[scan$race.cat=="SouthAmerican",1],pc$vectors[scan$race.cat=="SouthAmerican",2],pch=20,col="magenta")
points(-pc$vectors[scan$race.cat=="PuertoRican",1],pc$vectors[scan$race.cat=="PuertoRican",2],pch=20,col="green")
points(-pc$vectors[scan$race.cat=="Cuban",1],pc$vectors[scan$race.cat=="Cuban",2],pch=20,col="blue")
points(-pc$vectors[scan$race.cat=="Dominican",1],pc$vectors[scan$race.cat=="Dominican",2],pch=20,col="burlywood")
points(-pc$vectors[scan$race.cat=="Other",1],pc$vectors[scan$race.cat=="Other",2],pch=20,col="black")
points(-pc$vectors[scan$race.cat=="Unknown",1],pc$vectors[scan$race.cat=="Unknown",2],pch=20,col="black")
#legend("bottomright",c("Mexican"))
dev.off()


pdf("tmp_pc_col_ev34.pdf")
plot(-pc$vectors[,3],pc$vectors[,4],xlab="PC1",ylab="PC2",type="n")
points(-pc$vectors[scan$race.cat=="Mexican",3],pc$vectors[scan$race.cat=="Mexican",4],pch=20,col="goldenrod")
points(-pc$vectors[scan$race.cat=="CentralAmerican",3],pc$vectors[scan$race.cat=="CentralAmerican",4],pch=20,col="red")
points(-pc$vectors[scan$race.cat=="SouthAmerican",3],pc$vectors[scan$race.cat=="SouthAmerican",4],pch=20,col="magenta")
points(-pc$vectors[scan$race.cat=="PuertoRican",3],pc$vectors[scan$race.cat=="PuertoRican",4],pch=20,col="green")
points(-pc$vectors[scan$race.cat=="Cuban",3],pc$vectors[scan$race.cat=="Cuban",4],pch=20,col="blue")
points(-pc$vectors[scan$race.cat=="Dominican",3],pc$vectors[scan$race.cat=="Dominican",4],pch=20,col="burlywood")
points(-pc$vectors[scan$race.cat=="Other",3],pc$vectors[scan$race.cat=="Other",4],pch=20,col="black")
points(-pc$vectors[scan$race.cat=="Unknown",3],pc$vectors[scan$race.cat=="Unknown",4],pch=20,col="black")
#legend("bottomright",c("Mexican"))
dev.off()

pc_colors <- rep("black",sum(scan$geno.cntl==0&is.na(scan$gengrp6.outliers)))
pc_colors[scan$race.cat=="Mexican"] <- "goldenrod"
pc_colors[scan$race.cat=="CentralAmerican"] <- "red"
pc_colors[scan$race.cat=="SouthAmerican"] <- "magenta"
pc_colors[scan$race.cat=="PuertoRican"] <- "green"
pc_colors[scan$race.cat=="Cuban"] <- "blue"
pc_colors[scan$race.cat=="Dominican"] <- "burlywood"

library(MASS)
pdf("tmp_parCords.pdf",width=14)
#pc$vectors[,1] <- pc$vectors[,1]*-1
parcoord(pc$vectors[scan$geno.cntl==0&is.na(scan$gengrp6.outliers),],col=pc_colors[scan$geno.cntl==0&is.na(scan$gengrp6.outliers)])
dev.off()

# Mexican - yellow/goldenrod
# CentralAmerican - red
# SouthAmerican - magenta
# PuertoRican - green 
# Cuban - blue
# Dominican - burlywood
# Other - black
# Multi

rm(list=ls())


#####
# 3. Get pruned set of X chr SNPs

library(GWASTools)
library(gdsfmt); library(SNPRelate)
setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")

scanAnnot <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_scanAnnot_v54_AMS.RData"))
scanAnnot # 14160 95
head(pData(scanAnnot)) # it's col subjectID

snp <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_HCHS_Custom_15041502_B3_all37_v25_AMS.RData"))
head(pData(snp))
snp.sel <- snp$snpID[snp$chromosome==23&snp$quality.filter]
length(snp.sel) # 50781

scan.sel <- scanAnnot$scanID[is.na(scanAnnot$gengrp6.outliers)&scanAnnot$geno.cntl==0&scanAnnot$subj.plink&scanAnnot$unrelated.pcair.deg4]
length(scan.sel) # 10272

gdsobj <- snpgdsOpen("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/netCDF/subjects/OLGA_subj_geno.gds")
snpset <- snpgdsLDpruning(gdsobj, sample.id=scan.sel, snp.id=snp.sel,
                          autosome.only=FALSE, maf=0.05, missing.rate=0.05,
                          method="corr", slide.max.bp=10*1e6, ld.threshold=0.32, num.thread=1)
#Chromosome 23: 6.44%, 3600/55905
#3600 SNPs are selected in total.

closefn.gds(gdsobj)

snp.pruned <- unlist(snpset, use.names=FALSE)
length(snp.pruned) # 3600
save(snp.pruned, file="olga_application/snp_sel_xChr_ldPruned.RData")

rm(list=ls())


#####
# 4. Run PC-AiR on pruned set of SNPs; smaller sample set

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")
source("pcair_X_scan.include.R")
library(GWASTools)

scan <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_scanAnnot_v54_AMS.RData"))
head(pData(scan)) 

snp.pruned <- get(load("olga_application/snp_sel_xChr_ldPruned.RData"))
head(snp.pruned); length(snp.pruned) # 3600

# make genoData object
gdsobj <- GdsGenotypeReader("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/netCDF/subjects/OLGA_subj_geno.gds")
genoData <- GenotypeData(gdsobj,scanAnnot=scan[scan$subj.plink,])

unrel.set <- scan$scanID[scan$unrelated.pcair.deg4&is.na(scan$gengrp6.outliers)]
scanIncl <- scan$scanID[is.na(scan$gengrp6.outliers)&scan$geno.cntl==0&scan$subj.plink]

head(unrel.set); head(scanIncl)

pc <- pcair(genoData,unrel.set=unrel.set,Xchr=TRUE,snp.include=snp.pruned,scan.include=scanIncl)
# note this inclues HM samples

save(pc,file="olga_application/pca_prunedXsnps_10272unrel.RData")

pdf("pca_x_ev12.pdf")
plot(pc$vectors[,1],pc$vectors[,2],xlab="PC1",ylab="PC2")
dev.off()
# wow -- looks almost exactly the same as when using non-pruned SNP set

scan <- scan[scan$subj.plink&is.na(scan$gengrp6.outliers)&scan$geno.cntl==0,]

pc$values/pc$sum.values # 0.034817269 0.018669074 0.005356154 0.004933726

pdf("pca_x_ev12_col.pdf")
plot(-pc$vectors[,1],pc$vectors[,2],xlab="EV1 (3.48%)",ylab="EV2 (1.87%)",type="n")
points(-pc$vectors[scan$race.cat=="Mexican",1],pc$vectors[scan$race.cat=="Mexican",2],pch=20,col="goldenrod")
points(-pc$vectors[scan$race.cat=="CentralAmerican",1],pc$vectors[scan$race.cat=="CentralAmerican",2],pch=20,col="red")
points(-pc$vectors[scan$race.cat=="SouthAmerican",1],pc$vectors[scan$race.cat=="SouthAmerican",2],pch=20,col="magenta")
points(-pc$vectors[scan$race.cat=="PuertoRican",1],pc$vectors[scan$race.cat=="PuertoRican",2],pch=20,col="green")
points(-pc$vectors[scan$race.cat=="Cuban",1],pc$vectors[scan$race.cat=="Cuban",2],pch=20,col="blue")
points(-pc$vectors[scan$race.cat=="Dominican",1],pc$vectors[scan$race.cat=="Dominican",2],pch=20,col="burlywood")
points(-pc$vectors[scan$race.cat=="Other",1],pc$vectors[scan$race.cat=="Other",2],pch=20,col="black")
points(-pc$vectors[scan$race.cat=="Unknown",1],pc$vectors[scan$race.cat=="Unknown",2],pch=20,col="black")
legend("topleft",c("CentralAmerican","Cuban","Dominican","Mexican","PuertoRican","SouthAmerican","Other/Unknown"),pch=20,
       col=c("red","blue","burlywood","goldenrod","green","magenta","black"))
dev.off()

png("pca_x_ev12_col.png")
plot(-pc$vectors[,1],pc$vectors[,2],xlab="EV1 (3.48%)",ylab="EV2 (1.87%)",type="n",
     main="PCA with 3,600 LD-pruned X Chromosome SNPs",cex.main=1.5)
points(-pc$vectors[scan$race.cat=="Mexican",1],pc$vectors[scan$race.cat=="Mexican",2],pch=20,col="goldenrod")
points(-pc$vectors[scan$race.cat=="CentralAmerican",1],pc$vectors[scan$race.cat=="CentralAmerican",2],pch=20,col="red")
points(-pc$vectors[scan$race.cat=="SouthAmerican",1],pc$vectors[scan$race.cat=="SouthAmerican",2],pch=20,col="magenta")
points(-pc$vectors[scan$race.cat=="PuertoRican",1],pc$vectors[scan$race.cat=="PuertoRican",2],pch=20,col="green")
points(-pc$vectors[scan$race.cat=="Cuban",1],pc$vectors[scan$race.cat=="Cuban",2],pch=20,col="blue")
points(-pc$vectors[scan$race.cat=="Dominican",1],pc$vectors[scan$race.cat=="Dominican",2],pch=20,col="burlywood")
points(-pc$vectors[scan$race.cat=="Other",1],pc$vectors[scan$race.cat=="Other",2],pch=20,col="black")
points(-pc$vectors[scan$race.cat=="Unknown",1],pc$vectors[scan$race.cat=="Unknown",2],pch=20,col="black")
legend("topleft",c("CentralAmerican","Cuban","Dominican","Mexican","PuertoRican","SouthAmerican","Other/Unknown"),pch=20,
       col=c("red","blue","burlywood","goldenrod","green","magenta","black"))
dev.off()

# do ev1 vs ev2 for all 7 entries, then a legend in the 8th panel
pdf("pca_x_ev12_eachCol.pdf",width=14)
par(mfrow=c(2,4))
plot(-pc$vectors[scan$race.cat=="Mexican",1],pc$vectors[scan$race.cat=="Mexican",2],pch=20,col="goldenrod",xlab="EV1",ylab="EV2",main="Mexican")
plot(-pc$vectors[scan$race.cat=="CentralAmerican",1],pc$vectors[scan$race.cat=="CentralAmerican",2],pch=20,col="red",xlab="EV1",ylab="EV2",main="Central Am")
plot(-pc$vectors[scan$race.cat=="SouthAmerican",1],pc$vectors[scan$race.cat=="SouthAmerican",2],pch=20,col="magenta",xlab="EV1",ylab="EV2",main="South Am")
plot(-pc$vectors[scan$race.cat=="PuertoRican",1],pc$vectors[scan$race.cat=="PuertoRican",2],pch=20,col="green",xlab="EV1",ylab="EV2",main="Puerto Rican")
plot(-pc$vectors[scan$race.cat=="Cuban",1],pc$vectors[scan$race.cat=="Cuban",2],pch=20,col="blue",xlab="EV1",ylab="EV2",main="Cuban")
plot(-pc$vectors[scan$race.cat=="Dominican",1],pc$vectors[scan$race.cat=="Dominican",2],pch=20,col="burlywood",xlab="EV1",ylab="EV2",main="Dominican")
plot(-pc$vectors[is.element(scan$race.cat,c("Unknown","Other")),1],pc$vectors[is.element(scan$race.cat,c("Unknown","Other")),2],pch=20,col="black",xlab="EV1",ylab="EV2",main="Other/Unk")
plot(1:10,1:10,axes=FALSE,xlab="",ylab="",type="n")
legend("topleft",c("CentralAmerican","Cuban","Dominican","Mexican","PuertoRican","SouthAmerican","Other/Unknown"),pch=20,
       col=c("red","blue","burlywood","goldenrod","green","magenta","black"),cex=1.5)
dev.off()

pdf("pca_x_ev34_col.pdf")
plot(-pc$vectors[,3],pc$vectors[,4],xlab="EV3 (0.54%)",ylab="EV4 (0.49%)",type="n")
points(-pc$vectors[scan$race.cat=="Mexican",3],pc$vectors[scan$race.cat=="Mexican",4],pch=20,col="goldenrod")
points(-pc$vectors[scan$race.cat=="CentralAmerican",3],pc$vectors[scan$race.cat=="CentralAmerican",4],pch=20,col="red")
points(-pc$vectors[scan$race.cat=="SouthAmerican",3],pc$vectors[scan$race.cat=="SouthAmerican",4],pch=20,col="magenta")
points(-pc$vectors[scan$race.cat=="PuertoRican",3],pc$vectors[scan$race.cat=="PuertoRican",4],pch=20,col="green")
points(-pc$vectors[scan$race.cat=="Cuban",3],pc$vectors[scan$race.cat=="Cuban",4],pch=20,col="blue")
points(-pc$vectors[scan$race.cat=="Dominican",3],pc$vectors[scan$race.cat=="Dominican",4],pch=20,col="burlywood")
points(-pc$vectors[scan$race.cat=="Other",3],pc$vectors[scan$race.cat=="Other",4],pch=20,col="black")
points(-pc$vectors[scan$race.cat=="Unknown",3],pc$vectors[scan$race.cat=="Unknown",4],pch=20,col="black")
#legend("topleft",c("CentralAmerican","Cuban","Dominican","Mexican","PuertoRican","SouthAmerican","Other/Unknown"),pch=20,
#       col=c("red","blue","burlywood","goldenrod","green","magenta","black"))
dev.off()

png("pca_x_ev34_col.png")
plot(-pc$vectors[,3],pc$vectors[,4],xlab="EV3 (0.54%)",ylab="EV4 (0.49%)",type="n")
points(-pc$vectors[scan$race.cat=="Mexican",3],pc$vectors[scan$race.cat=="Mexican",4],pch=20,col="goldenrod")
points(-pc$vectors[scan$race.cat=="CentralAmerican",3],pc$vectors[scan$race.cat=="CentralAmerican",4],pch=20,col="red")
points(-pc$vectors[scan$race.cat=="SouthAmerican",3],pc$vectors[scan$race.cat=="SouthAmerican",4],pch=20,col="magenta")
points(-pc$vectors[scan$race.cat=="PuertoRican",3],pc$vectors[scan$race.cat=="PuertoRican",4],pch=20,col="green")
points(-pc$vectors[scan$race.cat=="Cuban",3],pc$vectors[scan$race.cat=="Cuban",4],pch=20,col="blue")
points(-pc$vectors[scan$race.cat=="Dominican",3],pc$vectors[scan$race.cat=="Dominican",4],pch=20,col="burlywood")
points(-pc$vectors[scan$race.cat=="Other",3],pc$vectors[scan$race.cat=="Other",4],pch=20,col="black")
points(-pc$vectors[scan$race.cat=="Unknown",3],pc$vectors[scan$race.cat=="Unknown",4],pch=20,col="black")
#legend("topleft",c("CentralAmerican","Cuban","Dominican","Mexican","PuertoRican","SouthAmerican","Other/Unknown"),pch=20,
#       col=c("red","blue","burlywood","goldenrod","green","magenta","black"))
dev.off()



# read in olga results
olga <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/amstilp/results/pca/v08_study_unrelated_pcair/pca.v08_study_unrel.RData"))
dim(olga$eigenvect) # 10642    32

# plot olga ev1 auto vs xchr ev1
# calculate corr

# subset pc$vectors to be samples that are included in OLGA results
sum(is.element(rownames(pc$vectors),olga$sample.id)) # 10594
pc$vectors <- pc$vectors[is.element(rownames(pc$vectors),olga$sample.id),]
dim(pc$vectors) # 10594    10
olga$eigenvect <- olga$eigenvect[is.element(olga$sample.id,rownames(pc$vectors)),]
dim(olga$eigenvect) # 10594 32

cor(olga$eigenvect[,1],-pc$vectors[,1]) # 0.8704225
pdf("pca_x_auto_ev1.pdf")
plot(olga$eigenvect[,1],-pc$vectors[,1],xlab="Autosomal EV1",ylab="X Chr EV1")
legend("topleft",c(expression(paste(rho,"= 0.870"))))
dev.off()

# color the points
scan <- scan[is.element(scan$scanID,rownames(pc$vectors)),]
dim(scan); dim(pc$vectors) # both have 10594 samples

pdf("pca_x_auto_ev1_col.pdf")
plot(-pc$vectors[,1]~olga$eigenvect[,1],xlab="Autosomal EV1",ylab="X Chr EV1",type="n")
points(-pc$vectors[scan$race.cat=="Mexican",1]~olga$eigenvect[scan$race.cat=="Mexican",1],pch=20,col="goldenrod")
points(-pc$vectors[scan$race.cat=="CentralAmerican",1]~olga$eigenvect[scan$race.cat=="CentralAmerican",1],pch=20,col="red")
points(-pc$vectors[scan$race.cat=="SouthAmerican",1]~olga$eigenvect[scan$race.cat=="SouthAmerican",1],pch=20,col="magenta")
points(-pc$vectors[scan$race.cat=="PuertoRican",1]~olga$eigenvect[scan$race.cat=="PuertoRican",1],pch=20,col="green")
points(-pc$vectors[scan$race.cat=="Cuban",1]~olga$eigenvect[scan$race.cat=="Cuban",1],pch=20,col="blue")
points(-pc$vectors[scan$race.cat=="Dominican",1]~olga$eigenvect[scan$race.cat=="Dominican",1],pch=20,col="burlywood")
points(-pc$vectors[scan$race.cat=="Other",1]~olga$eigenvect[scan$race.cat=="Other",1],pch=20,col="black")
points(-pc$vectors[scan$race.cat=="Unknown",1]~olga$eigenvect[scan$race.cat=="Unknown",1],pch=20,col="black")
legend("topleft",c("CentralAmerican","Cuban","Dominican","Mexican","PuertoRican","SouthAmerican","Other/Unknown"),pch=20,
       col=c("red","blue","burlywood","goldenrod","green","magenta","black"),cex=0.9)
legend("bottomright",c(expression(paste(rho,"= 0.870"))),bty="n")
dev.off()

# try colored by sex and see what happens
pdf("pca_x_auto_ev1_sexCol.pdf")
plot(-pc$vectors[,1]~olga$eigenvect[,1],xlab="Autosomal EV1",ylab="X Chr EV1",type="n")
points(-pc$vectors[scan$sex=="M",1]~olga$eigenvect[scan$sex=="M",1],pch=20,col="cyan")
points(-pc$vectors[scan$sex=="F",1]~olga$eigenvect[scan$sex=="F",1],pch=20,col="magenta")
legend("topleft",c("Male","Female"),col=c("cyan","magenta"),pch=19,bty="n")
dev.off()

##
# plot ev2 auto vs xchr ev2
# calculate corr
cor(olga$eigenvect[,2],pc$vectors[,2]) # 0.7441571
pdf("pca_x_auto_ev2.pdf")
plot(olga$eigenvect[,2],pc$vectors[,2],xlab="Autosomal EV2",ylab="X Chr EV2")
legend("topleft",c(expression(paste(rho,"= 0.744"))))
dev.off()

# color the points
scan <- scan[is.element(scan$scanID,rownames(pc$vectors)),]
dim(scan); dim(pc$vectors) # both have 10594 samples

pdf("pca_x_auto_ev2_col.pdf")
plot(pc$vectors[,2]~olga$eigenvect[,2],xlab="Autosomal EV2",ylab="X Chr EV2",type="n")
points(pc$vectors[scan$race.cat=="Mexican",2]~olga$eigenvect[scan$race.cat=="Mexican",2],pch=20,col="goldenrod")
points(pc$vectors[scan$race.cat=="CentralAmerican",2]~olga$eigenvect[scan$race.cat=="CentralAmerican",2],pch=20,col="red")
points(pc$vectors[scan$race.cat=="SouthAmerican",2]~olga$eigenvect[scan$race.cat=="SouthAmerican",2],pch=20,col="magenta")
points(pc$vectors[scan$race.cat=="PuertoRican",2]~olga$eigenvect[scan$race.cat=="PuertoRican",2],pch=20,col="green")
points(pc$vectors[scan$race.cat=="Cuban",2]~olga$eigenvect[scan$race.cat=="Cuban",2],pch=20,col="blue")
points(pc$vectors[scan$race.cat=="Dominican",2]~olga$eigenvect[scan$race.cat=="Dominican",2],pch=20,col="burlywood")
points(pc$vectors[scan$race.cat=="Other",2]~olga$eigenvect[scan$race.cat=="Other",2],pch=20,col="black")
points(pc$vectors[scan$race.cat=="Unknown",2]~olga$eigenvect[scan$race.cat=="Unknown",2],pch=20,col="black")
legend("topleft",c("CentralAmerican","Cuban","Dominican","Mexican","PuertoRican","SouthAmerican","Other/Unknown"),pch=20,
       col=c("red","blue","burlywood","goldenrod","green","magenta","black"),cex=0.9)
legend("bottomright",c(expression(paste(rho,"= 0.744"))),bty="n")
dev.off()

# try colored by sex and see what happens
pdf("pca_x_auto_ev2_sexCol.pdf")
plot(pc$vectors[,2]~olga$eigenvect[,2],xlab="Autosomal EV2",ylab="X Chr EV2",type="n")
points(pc$vectors[scan$sex=="M",2]~olga$eigenvect[scan$sex=="M",2],pch=20,col="cyan")
points(pc$vectors[scan$sex=="F",2]~olga$eigenvect[scan$sex=="F",2],pch=20,col="magenta")
legend("topleft",c("Male","Female"),col=c("cyan","magenta"),pch=19,bty="n")
dev.off()

rm(list=ls())


#####
# 5. X chr relatedness estimates

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")
library(GWASTools)
source("pcrelate_X.R")
# use unrel.set to be the same unrelated set used for PCA on the X chr
# use scan.include to be the same set of individs used for PCA on the X chr

scan <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_scanAnnot_v54_AMS.RData"))
head(pData(scan)) 

# make genoData object
gdsobj <- GdsGenotypeReader("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/netCDF/subjects/OLGA_subj_geno.gds")
genoData <- GenotypeData(gdsobj,scanAnnot=scan[scan$subj.plink,])

unrel.set <- scan$scanID[scan$unrelated.pcair.deg4&is.na(scan$gengrp6.outliers)]
scanIncl <- scan$scanID[is.na(scan$gengrp6.outliers)&scan$geno.cntl==0&scan$subj.plink]

head(unrel.set); head(scanIncl)
length(unrel.set); length(scanIncl) # 10272 | 12747

rel <- pcrelate(genoData, unrel.set=unrel.set, scan.include=scanIncl, Xchr=TRUE)
names(rel)

save(rel,file="olga_application/unadj_pcRelate_Xchr.RData")


## run adj for first 5 PCs
pc <- get(load("olga_application/pca_prunedXsnps_10272unrel.RData"))
pcMat <- pc$vectors[,1:6]
dim(pcMat)

rel <- pcrelate(genoData,unrel.set=unrel.set,scan.include=scanIncl,Xchr=TRUE,pcMat=pcMat)
save(rel,file="olga_application/adj_pcRelate_Xchr.RData")
# results from pcrelate in KC are what i want to run pca on

rm(list=ls())


#####
# 6. Prune auto set of SNPs, incl and excl X chr

library(GWASTools)
library(gdsfmt); library(SNPRelate)
setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")

scanAnnot <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_scanAnnot_v54_AMS.RData"))
scanAnnot # 14160 95
head(pData(scanAnnot)) # it's col subjectID

snp <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_HCHS_Custom_15041502_B3_all37_v25_AMS.RData"))
head(pData(snp))
snpID <- getSnpID(snp)
snp.sel <- snp$snpID[snp$quality.filter]
length(snp.sel) # 2365803

# remove SNPs in spike regions
 filt <- get(data(list=paste("pcaSnpFilters", "hg19", sep=".")))
 chrom <- getChromosome(snp)
 pos <- getPosition(snp)
 snp.filt <- rep(TRUE, length(snpID))
for (f in 1:nrow(filt)) {
     snp.filt[chrom == filt$chrom[f] & filt$start.base[f] < pos & pos < filt$end.base[f]] <- FALSE
  }
snp.sel <- intersect(snp.sel, snpID[snp.filt])
length(snp.sel) # 2325793

scan.sel <- scanAnnot$scanID[is.na(scanAnnot$gengrp6.outliers)&scanAnnot$geno.cntl==0&scanAnnot$subj.plink&scanAnnot$unrelated.pcair.deg4]
length(scan.sel) # 10272

gdsobj <- snpgdsOpen("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/netCDF/subjects/OLGA_subj_geno.gds")
snpset <- snpgdsLDpruning(gdsobj, sample.id=scan.sel, snp.id=snp.sel,
                          autosome.only=TRUE, maf=0.05, missing.rate=0.05,
                          method="corr", slide.max.bp=10*1e6, ld.threshold=0.32, num.thread=1)
#Working space: 10272 samples, 1304735 SNPs
#  SNPs are selected in total.
closefn.gds(gdsobj)

snp.pruned <- unlist(snpset, use.names=FALSE)
length(snp.pruned) # 158778
save(snp.pruned, file="olga_application/snp_sel_auto_ldPruned.RData")

##

gdsobj <- snpgdsOpen("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/netCDF/subjects/OLGA_subj_geno.gds")
snpset <- snpgdsLDpruning(gdsobj, sample.id=scan.sel, snp.id=snp.sel,
                          autosome.only=FALSE, maf=0.05, missing.rate=0.05,
                          method="corr", slide.max.bp=10*1e6, ld.threshold=0.32, num.thread=1)
#Working space: 10272 samples, 1340345 SNPs
# 158778 SNPs are selected in total.

closefn.gds(gdsobj)

snp.pruned <- unlist(snpset, use.names=FALSE)
length(snp.pruned) # 158778
save(snp.pruned, file="olga_application/snp_sel_autoX_ldPruned.RData")

rm(list=ls())


#####
# 7. PCA on the autosomes+X incl proj to rel

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x/")
source("pcair_X_scan.include.R")
library(GWASTools)

scan <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_scanAnnot_v54_AMS.RData"))
head(pData(scan)) 

snp.pruned <- get(load("olga_application/snp_sel_autoX_ldPruned.RData"))
head(snp.pruned); length(snp.pruned) # 158778

# make genoData object
gdsobj <- GdsGenotypeReader("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/netCDF/subjects/OLGA_subj_geno.gds")
genoData <- GenotypeData(gdsobj,scanAnnot=scan[scan$subj.plink,])

unrel.set <- scan$scanID[scan$unrelated.pcair.deg4&is.na(scan$gengrp6.outliers)]
scanIncl <- scan$scanID[is.na(scan$gengrp6.outliers)&scan$geno.cntl==0&scan$subj.plink]

head(unrel.set); head(scanIncl)

pc <- pcair(genoData,unrel.set=unrel.set,Xchr=FALSE,snp.include=snp.pruned,scan.include=scanIncl)
# note this inclues HM samples

save(pc,file="olga_application/pca_prunedAutoXsnps_10272unrelPlusRel.RData")

scan <- scan[scan$subj.plink&is.na(scan$gengrp6.outliers)&scan$geno.cntl==0,]

pc$values/pc$sum.values # 0.0255098027 0.0108610336 0.0009410187 0.0007580892

pdf("olga_application/pca_autoX_ev12_col.pdf")
plot(-pc$vectors[,1],-pc$vectors[,2],xlab="EV1 (2.55%)",ylab="EV2 (1.09%)",type="n")
points(-pc$vectors[scan$race.cat=="Mexican",1],-pc$vectors[scan$race.cat=="Mexican",2],pch=20,col="goldenrod")
points(-pc$vectors[scan$race.cat=="CentralAmerican",1],-pc$vectors[scan$race.cat=="CentralAmerican",2],pch=20,col="red")
points(-pc$vectors[scan$race.cat=="SouthAmerican",1],-pc$vectors[scan$race.cat=="SouthAmerican",2],pch=20,col="magenta")
points(-pc$vectors[scan$race.cat=="PuertoRican",1],-pc$vectors[scan$race.cat=="PuertoRican",2],pch=20,col="green")
points(-pc$vectors[scan$race.cat=="Cuban",1],-pc$vectors[scan$race.cat=="Cuban",2],pch=20,col="blue")
points(-pc$vectors[scan$race.cat=="Dominican",1],-pc$vectors[scan$race.cat=="Dominican",2],pch=20,col="burlywood")
points(-pc$vectors[scan$race.cat=="Other",1],-pc$vectors[scan$race.cat=="Other",2],pch=20,col="black")
points(-pc$vectors[scan$race.cat=="Unknown",1],-pc$vectors[scan$race.cat=="Unknown",2],pch=20,col="black")
legend("topright",c("CentralAmerican","Cuban","Dominican","Mexican","PuertoRican","SouthAmerican","Other/Unknown"),pch=20,
       col=c("red","blue","burlywood","goldenrod","green","magenta","black"))
dev.off()

# do ev1 vs ev2 for all 7 entries, then a legend in the 8th panel
pdf("olga_application/pca_autoX_ev12_eachCol.pdf",width=14)
par(mfrow=c(2,4))
plot(-pc$vectors[scan$race.cat=="Mexican",1],-pc$vectors[scan$race.cat=="Mexican",2],pch=20,col="goldenrod",xlab="EV1",ylab="EV2",main="Mexican")
plot(-pc$vectors[scan$race.cat=="CentralAmerican",1],-pc$vectors[scan$race.cat=="CentralAmerican",2],pch=20,col="red",xlab="EV1",ylab="EV2",main="Central Am")
plot(-pc$vectors[scan$race.cat=="SouthAmerican",1],-pc$vectors[scan$race.cat=="SouthAmerican",2],pch=20,col="magenta",xlab="EV1",ylab="EV2",main="South Am")
plot(-pc$vectors[scan$race.cat=="PuertoRican",1],-pc$vectors[scan$race.cat=="PuertoRican",2],pch=20,col="green",xlab="EV1",ylab="EV2",main="Puerto Rican")
plot(-pc$vectors[scan$race.cat=="Cuban",1],-pc$vectors[scan$race.cat=="Cuban",2],pch=20,col="blue",xlab="EV1",ylab="EV2",main="Cuban")
plot(-pc$vectors[scan$race.cat=="Dominican",1],-pc$vectors[scan$race.cat=="Dominican",2],pch=20,col="burlywood",xlab="EV1",ylab="EV2",main="Dominican")
plot(-pc$vectors[is.element(scan$race.cat,c("Unknown","Other")),1],-pc$vectors[is.element(scan$race.cat,c("Unknown","Other")),2],pch=20,col="black",xlab="EV1",ylab="EV2",main="Other/Unk")
plot(1:10,1:10,axes=FALSE,xlab="",ylab="",type="n")
legend("topleft",c("CentralAmerican","Cuban","Dominican","Mexican","PuertoRican","SouthAmerican","Other/Unknown"),pch=20,
       col=c("red","blue","burlywood","goldenrod","green","magenta","black"),cex=1.5)
dev.off()

pdf("olga_application/pca_autoX_ev34_col.pdf")
plot(-pc$vectors[,3],pc$vectors[,4],xlab="EV3 (0.09%)",ylab="EV4 (0.07%)",type="n")
points(-pc$vectors[scan$race.cat=="Mexican",3],pc$vectors[scan$race.cat=="Mexican",4],pch=20,col="goldenrod")
points(-pc$vectors[scan$race.cat=="CentralAmerican",3],pc$vectors[scan$race.cat=="CentralAmerican",4],pch=20,col="red")
points(-pc$vectors[scan$race.cat=="SouthAmerican",3],pc$vectors[scan$race.cat=="SouthAmerican",4],pch=20,col="magenta")
points(-pc$vectors[scan$race.cat=="PuertoRican",3],pc$vectors[scan$race.cat=="PuertoRican",4],pch=20,col="green")
points(-pc$vectors[scan$race.cat=="Cuban",3],pc$vectors[scan$race.cat=="Cuban",4],pch=20,col="blue")
points(-pc$vectors[scan$race.cat=="Dominican",3],pc$vectors[scan$race.cat=="Dominican",4],pch=20,col="burlywood")
points(-pc$vectors[scan$race.cat=="Other",3],pc$vectors[scan$race.cat=="Other",4],pch=20,col="black")
points(-pc$vectors[scan$race.cat=="Unknown",3],pc$vectors[scan$race.cat=="Unknown",4],pch=20,col="black")
#legend("topleft",c("CentralAmerican","Cuban","Dominican","Mexican","PuertoRican","SouthAmerican","Other/Unknown"),pch=20,
#       col=c("red","blue","burlywood","goldenrod","green","magenta","black"))
dev.off()



# read in olga results
olga <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/amstilp/results/pca/v08_study_unrelated_pcair/pca.v08_study_unrel.RData"))
dim(olga$eigenvect) # 10642    32

olga$eigenvect[,1:5]/olga$eigenval[,1:5] # 0.0255098027 0.0108610336 0.0009410187 0.0007580892


# plot olga ev1 auto vs xchr ev1
# calculate corr

# subset pc$vectors to be samples that are included in OLGA results
sum(is.element(rownames(pc$vectors),olga$sample.id)) # 10594
pc$vectors <- pc$vectors[is.element(rownames(pc$vectors),olga$sample.id),]
dim(pc$vectors) # 10594    10
olga$eigenvect <- olga$eigenvect[is.element(olga$sample.id,rownames(pc$vectors)),]
dim(olga$eigenvect) # 10594 32

cor(olga$eigenvect[,1],-pc$vectors[,1]) # 0.9997652
pdf("olga_application/pca_autoX_auto_ev1.pdf")
plot(olga$eigenvect[,1],-pc$vectors[,1],xlab="Autosomal EV1",ylab="Auto+X Chr EV1")
legend("topleft",c(expression(paste(rho,"= 0.999"))))
dev.off()

# color the points
scan <- scan[is.element(scan$scanID,rownames(pc$vectors)),]
dim(scan); dim(pc$vectors) # both have 10594 samples

pdf("pca_autoX_auto_ev1_col.pdf")
plot(-pc$vectors[,1]~olga$eigenvect[,1],xlab="Autosomal EV1",ylab="Auto+X Chr EV1",type="n")
points(-pc$vectors[scan$race.cat=="Mexican",1]~olga$eigenvect[scan$race.cat=="Mexican",1],pch=20,col="goldenrod")
points(-pc$vectors[scan$race.cat=="CentralAmerican",1]~olga$eigenvect[scan$race.cat=="CentralAmerican",1],pch=20,col="red")
points(-pc$vectors[scan$race.cat=="SouthAmerican",1]~olga$eigenvect[scan$race.cat=="SouthAmerican",1],pch=20,col="magenta")
points(-pc$vectors[scan$race.cat=="PuertoRican",1]~olga$eigenvect[scan$race.cat=="PuertoRican",1],pch=20,col="green")
points(-pc$vectors[scan$race.cat=="Cuban",1]~olga$eigenvect[scan$race.cat=="Cuban",1],pch=20,col="blue")
points(-pc$vectors[scan$race.cat=="Dominican",1]~olga$eigenvect[scan$race.cat=="Dominican",1],pch=20,col="burlywood")
points(-pc$vectors[scan$race.cat=="Other",1]~olga$eigenvect[scan$race.cat=="Other",1],pch=20,col="black")
points(-pc$vectors[scan$race.cat=="Unknown",1]~olga$eigenvect[scan$race.cat=="Unknown",1],pch=20,col="black")
legend("topleft",c("CentralAmerican","Cuban","Dominican","Mexican","PuertoRican","SouthAmerican","Other/Unknown"),pch=20,
       col=c("red","blue","burlywood","goldenrod","green","magenta","black"),cex=0.9)
legend("bottomright",c(expression(paste(rho,"= 0.999"))),bty="n")
dev.off()

# plot auto pc1 vs auto pc2
scan <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_scanAnnot_v54_AMS.RData"))
scan <- scan[is.element(scan$scanID,olga$sample.id),]
dim(scan)

png("pca_auto_ev12_col.png")
plot(olga$eigenvect[,1],olga$eigenvect[,2],xlab="EV1 (2.62%)",ylab="EV2 (1.23%)",main="Autosomal PCA",
     cex.main=1.5,type="n")
points(olga$eigenvect[scan$race.cat=="Mexican"&is.na(scan$gengrp6.outliers),1],
       olga$eigenvect[scan$race.cat=="Mexican"&is.na(scan$gengrp6.outliers),2],
       pch=20,col="goldenrod")
points(olga$eigenvect[scan$race.cat=="CentralAmerican"&is.na(scan$gengrp6.outliers),1],
       olga$eigenvect[scan$race.cat=="CentralAmerican"&is.na(scan$gengrp6.outliers),2],
       pch=20,col="red")
points(olga$eigenvect[scan$race.cat=="SouthAmerican"&is.na(scan$gengrp6.outliers),1],
       olga$eigenvect[scan$race.cat=="SouthAmerican"&is.na(scan$gengrp6.outliers),2],
       pch=20,col="magenta")
points(olga$eigenvect[scan$race.cat=="PuertoRican"&is.na(scan$gengrp6.outliers),1],
       olga$eigenvect[scan$race.cat=="PuertoRican"&is.na(scan$gengrp6.outliers),2],
       pch=20,col="green")
points(olga$eigenvect[scan$race.cat=="Cuban"&is.na(scan$gengrp6.outliers),1],
       olga$eigenvect[scan$race.cat=="Cuban"&is.na(scan$gengrp6.outliers),2],pch=20,col="blue")
points(olga$eigenvect[scan$race.cat=="Dominican"&is.na(scan$gengrp6.outliers),1],
       olga$eigenvect[scan$race.cat=="Dominican"&is.na(scan$gengrp6.outliers),2],
       pch=20,col="burlywood")
points(olga$eigenvect[scan$race.cat=="Other"&is.na(scan$gengrp6.outliers),1],
       olga$eigenvect[scan$race.cat=="Other"&is.na(scan$gengrp6.outliers),2],
       pch=20,col="black")
points(olga$eigenvect[scan$race.cat=="Unknown"&is.na(scan$gengrp6.outliers),1],
       olga$eigenvect[scan$race.cat=="Unknown"&is.na(scan$gengrp6.outliers),2],
       pch=20,col="black")
legend("topleft",c("CentralAmerican","Cuban","Dominican","Mexican","PuertoRican","SouthAmerican","Other/Unknown"),pch=20,
       col=c("red","blue","burlywood","goldenrod","green","magenta","black"),cex=0.9)
dev.off()



##
# plot ev2 auto vs xchr ev2
# calculate corr
cor(olga$eigenvect[,2],-pc$vectors[,2]) # 0.9991056
pdf("pca_autoX_auto_ev2.pdf")
plot(olga$eigenvect[,2],-pc$vectors[,2],xlab="Autosomal EV2",ylab="Auto+X Chr EV2")
legend("topleft",c(expression(paste(rho,"= 0.999"))))
dev.off()

# look at pcs further down
library(MASS)
pdf("olga_application/pca_autoX_parCoords.pdf",width=14)
plotCol <- rep("black",nrow(scan))
plotCol[scan$race.cat=="SouthAmerican"] <- "magenta"
plotCol[scan$race.cat=="Mexican"] <- "goldenrod"
plotCol[scan$race.cat=="CentralAmerican"] <- "red"
plotCol[scan$race.cat=="PuertoRican"] <- "green"
plotCol[scan$race.cat=="Cuban"] <- "blue"
plotCol[scan$race.cat=="Dominican"] <- "burlywood"
parcoord(pc$vectors[,1:10],col=plotCol)
dev.off()

png("olga_application/pca_autoX_parCoords.png",width=960)
plotCol <- rep("black",nrow(scan))
plotCol[scan$race.cat=="SouthAmerican"] <- "magenta"
plotCol[scan$race.cat=="Mexican"] <- "goldenrod"
plotCol[scan$race.cat=="CentralAmerican"] <- "red"
plotCol[scan$race.cat=="PuertoRican"] <- "green"
plotCol[scan$race.cat=="Cuban"] <- "blue"
plotCol[scan$race.cat=="Dominican"] <- "burlywood"
parcoord(pc$vectors[,1:10],col=plotCol)
dev.off()



rm(list=ls())


#####
# 8. PCA on the autosomes+X, unrel only

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")

nSlots <- Sys.getenv("NSLOTS")
nThreads <- ifelse(is.na(strtoi(nSlots) >= 1), 1, strtoi(nSlots))
print(paste("Running with", nThreads,"thread(s)."))

scan <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_scanAnnot_v54_AMS.RData"))
head(pData(scan)) 
unrel.set <- scan$scanID[scan$unrelated.pcair.deg4&is.na(scan$gengrp6.outliers)]

snp.pruned <- get(load("olga_application/snp_sel_autoX_ldPruned.RData"))
head(snp.pruned); length(snp.pruned) # 158778

gdsobj <- snpgdsOpen("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/netCDF/subjects/OLGA_subj_geno.gds")
pca <- snpgdsPCA(gdsobj, sample.id=unrel.set, snp.id=snp.pruned, num.thread=nThreads)
closefn.gds(gdsobj)

save(pca,file="olga_application/pca_prunedAutoXsnps_10272unrel.RData")

rm(list=ls())


#####
# 9. PCA with auto SNPs

library(GWASTools)
library(SNPRelate)

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x/")

scan <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_scanAnnot_v54_AMS.RData"))
head(pData(scan)) 

snp.pruned <- get(load("olga_application/snp_sel_auto_ldPruned.RData"))
length(snp.pruned) # 152391

unrel.set <- scan$scanID[scan$unrelated.pcair.deg4&is.na(scan$gengrp6.outliers)]
scanIncl <- scan$scanID[is.na(scan$gengrp6.outliers)&scan$geno.cntl==0&scan$subj.plink]

head(unrel.set); head(scanIncl) 
length(unrel.set); length(scanIncl) # 10272 | 12747

gdsobj <- snpgdsOpen("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/netCDF/subjects/OLGA_subj_geno.gds")
pca <- snpgdsPCA(gdsobj, sample.id=unrel.set, snp.id=snp.pruned, num.thread=2)
save(pca,file="olga_application/pca_prunedAutosnps_10272unrel.RData")

relStudy <- scanIncl[!is.element(scanIncl,unrel.set)]
length(relStudy) # 2475

SnpLoad <- snpgdsPCASNPLoading(pca,gdsobj,num.thread=nThreads)
eigsamp <- snpgdsPCASampLoading(SnpLoad,gdsobj,sample.id=relStudy,num.thread=nThreads,verbose=TRUE)
closefn.gds(gdsobj)

save(eigsamp,file="olga_application/pca_prunedAutosnps_2475rel.RData")


# put both sets of results into one data.frame
dir <- data.frame("sample.id"=pca$sample.id,"direct"=TRUE)
indir <- data.frame("sample.id"=eigsamp$sample.id,"direct"=FALSE)
sampleInfo <- rbind(dir,indir)
res <- rbind(pca$eigenvect,eigsamp$eigenvect)
dim(res) # 12747 32

totalRes <- list("sampleInfo"=sampleInfo,"eigenvect"=res)

save(totalRes,file="olga_application/pca_prunedAutosnps_10272unrelPlusRel.RData")

rm(list=ls())


#####
# 10. PC-Relate on X chr SNPs, corr for auto PC1-5

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x/")
library(GWASTools)
library(SNPRelate)
source("pcrelate_X.R")

autoPC <- get(load("olga_application/pca_prunedAutoXsnps_10272unrelPlusRel.RData"))

pcMat <- autoPC$vectors[,c(1:5)]
dim(pcMat) # 12747 5

scan <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_scanAnnot_v54_AMS.RData"))
head(pData(scan)) 

# make genoData object
gdsobj <- GdsGenotypeReader("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/netCDF/subjects/OLGA_subj_geno.gds")
genoData <- GenotypeData(gdsobj,scanAnnot=scan[scan$subj.plink,])

unrel.set <- scan$scanID[scan$unrelated.pcair.deg4&is.na(scan$gengrp6.outliers)]
scanIncl <- scan$scanID[is.na(scan$gengrp6.outliers)&scan$geno.cntl==0&scan$subj.plink]

head(unrel.set); head(scanIncl)
length(unrel.set); length(scanIncl) # 10272 | 12747

### need to exclude 13 people with an entirely missing x chr 
# from chromosome.anomalies.csv file, they are study samples with x chr anomaly
# thus, their entire x chr is filtered out. 

chromAnoms <- read.csv("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/dbGaP/To_dbGaP/chromosome_anomalies/chromosome_anomalies.csv",
                       header=TRUE,as.is=T)
ids <- chromAnoms$scanID[chromAnoms$chromosome=="X"&chromAnoms$filter&chromAnoms$whole.chr&
                           !is.element(substr(chromAnoms$SUBJECT_ID,1,2),"NA")]
length(ids) # 13
chromAnoms[is.element(chromAnoms$scanID,ids),] # perfect! exactly what we want

length(scanIncl) # 12747
scanIncl <- scanIncl[!is.element(scanIncl,ids)]
length(scanIncl) # 12734

# take these individs out of pcMat too
pcMat <- pcMat[is.element(rownames(pcMat),scanIncl),]
dim(pcMat) # 12734 5

rel <- pcrelate(genoData, unrel.set=unrel.set, scan.include=scanIncl, Xchr=TRUE, pcMat=pcMat)
names(rel)

save(rel,file="olga_application/pcRelate_Xchr_autoXPC15adj.RData")

## make into matrix
tmp <- matrix(NA,nrow=12734,ncol=12734)
length(tmp[lower.tri(tmp,diag=FALSE)]) # 81071011

kc <- rel$kinship
length(kc[,"kin"]) # 81071011
tmp[lower.tri(tmp,diag=FALSE)] <- kc[,"kin"]

# mirror to upper triangle
tmp[upper.tri(tmp,diag=FALSE)] <- t(tmp)[upper.tri(tmp,diag=FALSE)]

rownames(tmp) <- scanIncl
colnames(tmp) <- scanIncl

diag(tmp) <- 0 # diag elements don't matter

save(tmp,file="olga_application/kcMat_xchr_autoXPC15adj.RData")

rm(list=ls())


#####
# 11. PC-Relate on X chr SNPs, corr for X chr PC1-2

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x/")
library(GWASTools)
library(SNPRelate)
source("pcrelate_X.R")

xPC <- get(load("olga_application/pca_prunedXsnps_10272unrel.RData"))

pcMat <- xPC$vectors[,c(1:2)]
dim(pcMat) # 12747 2

scan <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_scanAnnot_v54_AMS.RData"))
head(pData(scan)) 

# make genoData object
gdsobj <- GdsGenotypeReader("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/netCDF/subjects/OLGA_subj_geno.gds")
genoData <- GenotypeData(gdsobj,scanAnnot=scan[scan$subj.plink,])

unrel.set <- scan$scanID[scan$unrelated.pcair.deg4&is.na(scan$gengrp6.outliers)]
scanIncl <- scan$scanID[is.na(scan$gengrp6.outliers)&scan$geno.cntl==0&scan$subj.plink]

head(unrel.set); head(scanIncl)
length(unrel.set); length(scanIncl) # 10272 | 12747

### need to exclude 13 people with an entirely missing x chr 
# from chromosome.anomalies.csv file, they are study samples with x chr anomaly
# thus, their entire x chr is filtered out. 

chromAnoms <- read.csv("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/dbGaP/To_dbGaP/chromosome_anomalies/chromosome_anomalies.csv",
                       header=TRUE,as.is=T)
ids <- chromAnoms$scanID[chromAnoms$chromosome=="X"&chromAnoms$filter&chromAnoms$whole.chr&
                           !is.element(substr(chromAnoms$SUBJECT_ID,1,2),"NA")]
length(ids) # 13
chromAnoms[is.element(chromAnoms$scanID,ids),] # perfect! exactly what we want

length(scanIncl) # 12747
scanIncl <- scanIncl[!is.element(scanIncl,ids)]
length(scanIncl) # 12734

# take these individs out of pcMat too
pcMat <- pcMat[is.element(rownames(pcMat),scanIncl),]
dim(pcMat) # 12734 2

rel <- pcrelate(genoData, unrel.set=unrel.set, scan.include=scanIncl, Xchr=TRUE, pcMat=pcMat)
names(rel)

save(rel,file="olga_application/pcRelate_Xchr_xPC12adj.RData")

## make into matrix
tmp <- matrix(NA,nrow=12734,ncol=12734)
length(tmp[lower.tri(tmp,diag=FALSE)]) # 81071011

kc <- rel$kinship
length(kc[,"kin"]) # 81071011
tmp[lower.tri(tmp,diag=FALSE)] <- kc[,"kin"]

# mirror to upper triangle
tmp[upper.tri(tmp,diag=FALSE)] <- t(tmp)[upper.tri(tmp,diag=FALSE)]

rownames(tmp) <- scanIncl
colnames(tmp) <- scanIncl

diag(tmp) <- 0 # diag elements don't matter

save(tmp,file="olga_application/kcMat_xchr_xPC12adj.RData")

rm(list=ls())


#####
# 12. Parse KC results

library(GWASTools)
setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")

scan <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_scanAnnot_v54_AMS.RData"))
head(pData(scan))

kcX_autoAdj <- get(load("olga_application/pcRelate_Xchr_autoXPC15adj.RData"))
kcX_autoAdj <- kcX_autoAdj$kinship
kcX_xAdj <- get(load("olga_application/pcRelate_Xchr_xPC12adj.RData"))
kcX_xAdj <- kcX_xAdj$kinship
kcAuto <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/mconomos/results/relatedness/grm/grm_allchr.RData"))
kcAuto <- kcAuto$kinship

kcX_unadj <- get(load("olga_application/unadj_pcRelate_Xchr.RData"))
kcX_unadj <- kcX_unadj$kinship

summary(kcX_unadj$kin) # ranges from -0.18 to 1.03
summary(kcX_autoAdj$kin) # ranges from -0.32 to 0.98
summary(kcX_xAdj$kin) # ranges from -0.16 to 0.97

kcX_xAdj[kcX_xAdj$kin>0.95,]
# ID1    ID2  nsnp       kin
# 646259 669782 31157 0.9656194
# 684581 700299 30186 0.9511280
# 703739 878890 25689 0.9679540
# 757248 821317 33791 0.9621796
# 799718 866519 30746 0.9611269
# 883106 888083 34155 0.9717185

kcX_autoAdj[kcX_autoAdj$kin>0.95,]
# ID1    ID2  nsnp       kin
# 707996 875457 31859 0.9566378
# 757248 821317 33427 0.9710367
# 799718 866519 32527 0.9586097
# 883106 888083 33550 0.9838231
# 896943 901284 28101 0.9513448
# 957205 972362 29310 0.9726953

# what?? only 3 of them are the same. 

pData(scan)[is.element(scan$scanID,c(646259,669782)),]
# woah! full brothers that have the same x chr from their mother!!

pData(scan)[is.element(scan$scanID,c(684581,700299)),] # same household but different family ids.
kcAuto[kcAuto$ID1==684581&kcAuto$ID2==700299,] # kc is 0.2886158, so brothers
kcX_autoAdj[kcX_autoAdj$ID1==684581&kcX_autoAdj$ID2==700299,] # kc is 0.9298
# why not same family id in the scan annot?

##
kcAuto_unr <- kcAuto[kcAuto$kin<0.025,]
dim(kcAuto_unr) # 78527239 7

# make hist of x chr kc estimates for pairs that are unrel on auto
unrel.set <- scan$scanID[scan$unrelated.pcair.deg4&is.na(scan$gengrp6.outliers)]

toHist <- kcX_autoAdj[is.element(kcX_autoAdj$ID1,unrel.set)&is.element(kcX_autoAdj$ID2,unrel.set),]
dim(toHist) # 52649191 4
length(unique(toHist$ID1)) # 10261; great! length of unrel.set

toHistXadj <- kcX_xAdj[is.element(kcX_xAdj$ID1,unrel.set)&is.element(kcX_xAdj$ID2,unrel.set),]
dim(toHistXadj) # 52649191 4
length(unique(toHistXadj$ID1)) # 10261

toHistNoAdj <- kcX_unadj[is.element(kcX_unadj$ID1,unrel.set)&is.element(kcX_unadj$ID2,unrel.set),]
dim(toHistNoAdj) # 52751856 4
length(unique(toHistNoAdj$ID1)) # 10271

pdf("olga_application/hist_xchrKC_autoXPC15adj_autoUnrel.pdf")
hist(toHist$kin,main="X chr KC for auto unrel samples",breaks=50,
     xlab="X Chr KC adj For Auto+X PC 1-5")
# the mean and ci lines were so small/together that they didn't add anything
#abline(v=mean(toHist$kin),col="red")
#abline(v=mean(toHist$kin)-1.96*sd(toHist$kin)/sqrt(nrow(toHist)),col="gray")
#abline(v=mean(toHist$kin)+1.96*sd(toHist$kin)/sqrt(nrow(toHist)),col="gray")
dev.off()

pdf("olga_application/hist_xchrKC_xPC12adj_autoUnrel.pdf")
hist(toHistXadj$kin,main="X chr KC for auto unrel samples",
     xlab="X Chr KC adj For X PC 1-2")

dev.off()

pdf("olga_application/hist_xchrKC_autoXPC15adj_autoUnrel_trunc.pdf")
hist(toHist$kin,main="X chr KC for auto unrel samples",breaks=50,
     xlab="X Chr KC adj For Auto+X PC 1-5",ylim=c(0,5000))
dev.off()

pdf("olga_application/hist_xchrKC_xPC12adj_autoUnrel_trunc.pdf")
hist(toHistXadj$kin,main="X chr KC for auto unrel samples",
     xlab="X Chr KC adj For X PC 1-2",ylim=c(0,5000))
dev.off()

pdf("olga_application/hist_autoKC.pdf")
hist(kcAuto$kin,main="Auto KC",breaks=50,
     xlab="Auto KC")
dev.off()

pdf("olga_application/hist_xchrKC_unadj_autoUnrel.pdf")
hist(toHistNoAdj$kin,main="X chr KC for auto unrel samples",breaks=50,
     xlab="X Chr KC unadj")
dev.off()

pdf("olga_application/hist_xchrKC_unadj_autoUnrel_trunc.pdf")
hist(toHistNoAdj$kin,main="X chr KC for auto unrel samples",breaks=50,
     xlab="X Chr KC unadj",ylim=c(0,5000))
dev.off()

## overlay all these histograms on eachother
pdf("olga_application/hist_allXchrKC.pdf",height=14,width=14)
par(mfrow=c(2,2))
hist(toHist$kin,main="X chr KC for auto unrel samples",breaks=50,
     col=rgb(1,0,0,0.5),xlab="X Chr KC",cex.lab=1.5,cex.main=1.5)
hist(toHistXadj$kin,add=TRUE,col=rgb(0,0,1,0.5),breaks=50,cex.lab=1.5)
hist(toHistNoAdj$kin,add=TRUE,col=rgb(0,1,0,0.5),breaks=50,cex.lab=1.5)
legend("topright",c("Auto+X PC 1-5 Adj","X Chr PC 1-2 Adj","No Adj"),col=c(rgb(1,0,0,0.5),rgb(0,0,1,0.5),rgb(0,1,0,0.5)),
       lty=1,lwd=5,cex=1.5)
hist(toHist$kin,main="",breaks=50,cex.lab=1.5,
     col=rgb(1,0,0,0.5),xlab="X Chr KC")
hist(toHistXadj$kin,col=rgb(0,0,1,0.5),breaks=50,main="",xlab="X Chr KC",cex.lab=1.5)
hist(toHistNoAdj$kin,col=rgb(0,1,0,0.5),breaks=50,main="",xlab="X Chr KC",cex.lab=1.5)
dev.off()


sd(toHist$kin); sd(toHistXadj$kin) # 0.02608632 | 0.02356558
# what is the 95% CI?
c(mean(toHist$kin)-1.96*sd(toHist$kin)/sqrt(nrow(toHist)),mean(toHist$kin)+1.96*sd(toHist$kin)/sqrt(nrow(toHist))) # 0.0138247 0.0138388
c(mean(toHistXadj$kin)-1.96*sd(toHistXadj$kin)/sqrt(nrow(toHistXadj)),mean(toHistXadj$kin)+1.96*sd(toHistXadj$kin)/sqrt(nrow(toHistXadj))) # 0.008863328 0.008876059

# get histograms of all estimates for each of the runs
pdf("olga_application/hist_xkc_autoXadj.pdf")
hist(kcX_autoAdj$kin,cex.lab=1.5,cex.main=1.5,main="X Chr KC for 12,734 Study Samples",
     breaks=50,xlab="X Chr KC")
legend("topright",c(paste("mean =",format(mean(kcX_autoAdj$kin),digits=4)),paste("sd =",format(sd(kcX_autoAdj$kin),digits=4))),
                    bty="n",cex=1.5)
dev.off()

pdf("olga_application/hist_xkc_autoXadj_trunc.pdf")
hist(kcX_autoAdj$kin,cex.lab=1.5,cex.main=1.5,main="X Chr KC for 12,734 Study Samples, Truncated",
     breaks=50,xlab="X Chr KC",ylim=c(0,500000))
dev.off()

pdf("olga_application/hist_xkc_Xadj.pdf")
hist(kcX_xAdj$kin,cex.lab=1.5,cex.main=1.5,main="X Chr KC for 12,734 Study Samples",
     breaks=50,xlab="X Chr KC")
legend("topright",c(paste("mean =",format(mean(kcX_xAdj$kin),digits=4)),paste("sd =",format(sd(kcX_xAdj$kin),digits=4))),
       bty="n",cex=1.5)
dev.off()

pdf("olga_application/hist_xkc_Xadj_trunc.pdf")
hist(kcX_xAdj$kin,cex.lab=1.5,cex.main=1.5,main="X Chr KC for 12,734 Study Samples, Truncated",
     breaks=50,xlab="X Chr KC",ylim=c(0,500000))
dev.off()

pdf("olga_application/hist_xkc_unadj.pdf")
hist(kcX_unadj$kin,cex.lab=1.5,cex.main=1.5,main="X Chr KC for 12,734 Study Samples",
     breaks=50,xlab="X Chr KC")
legend("topright",c(paste("mean =",format(mean(kcX_unadj$kin,na.rm=TRUE),digits=4)),
                    paste("sd =",format(sd(kcX_unadj$kin,na.rm=TRUE),digits=4))),
       bty="n",cex=1.5)
dev.off()

pdf("olga_application/hist_xkc_unadj_trunc.pdf")
hist(kcX_unadj$kin,cex.lab=1.5,cex.main=1.5,main="X Chr KC for 12,734 Study Samples, Truncated",
     breaks=50,xlab="X Chr KC",ylim=c(0,500000))
dev.off()


pdf("olga_application/hist_autokc.pdf")
hist(kcAuto$kin,cex.lab=1.5,cex.main=1.5,main="Auto KC for 12,734 Study Samples",
     breaks=50,xlab="Auto KC")
legend("topright",c(paste("mean =",format(mean(kcAuto$kin,na.rm=TRUE),digits=4)),
                    paste("sd =",format(sd(kcAuto$kin,na.rm=TRUE),digits=4))),
       bty="n",cex=1.5)
dev.off()

pdf("olga_application/hist_autokc_trunc.pdf")
hist(kcAuto$kin,cex.lab=1.5,cex.main=1.5,main="Auto KC for 12,734 Study Samples, Truncated",
     breaks=50,xlab="Auto KC",ylim=c(0,500000))
dev.off()


# plot x_kc vs auto_kc
# first examine those with kcX \in [0.5,1]
sum(kcX_xAdj$kin>0.5) # 495
sum(kcX_autoAdj$kin>0.5) # 511

smXadj <- kcX_xAdj[kcX_xAdj$kin>0.5,]
smXadj$ID12 <- paste(smXadj$ID1,smXadj$ID2)
smAutoadj <- kcX_autoAdj[kcX_autoAdj$kin>0.5,]
smAutoadj$ID12 <- paste(smAutoadj$ID1,smAutoadj$ID2)

toPl <- merge(smXadj,smAutoadj,by="ID12",all=TRUE)
head(toPl); dim(toPl) # 574 9
kcAuto$ID12 <- paste(kcAuto$ID1,kcAuto$ID2)

toPl <- merge(toPl,kcAuto,by="ID12",all.x=TRUE)
# .x is x adj, .y is auto adj, no ending is autoKC
head(toPl); dim(toPl) # 574 16

pdf("olga_application/xKC_autoKC.pdf")
plot(toPl$kin,toPl$kin.x,type="n",xlab="KC, Autosomes",ylab="KC, X Chromosome",ylim=c(0.5,1),xlim=c(0,0.5))
points(toPl$kin,toPl$kin.x,col="cyan",pch=19)
points(toPl$kin,toPl$kin.y,col="magenta",pch=19)
legend("topright",c("Adj for Auto+X PC 1-5","Adj for X chr PC 1-2"),col=c("magenta","cyan"),pch=19)
dev.off()

## look at some of these toPl pairs
toPl <- toPl[order(toPl$kin,decreasing=FALSE),]
# what are the sexes of each individ?
toPlSx <- merge(toPl,pData(scan)[,c("scanID","sex")],by.x="ID1.x",by.y="scanID",all.x=TRUE)
toPlSx2 <- merge(toPlSx,pData(scan)[,c("scanID","sex")],by.x="ID2.x",by.y="scanID",all.x=TRUE)

toPlSx2 <- toPlSx2[order(toPlSx2$kin,decreasing=FALSE),]

# first pair is 2 males, x chr kc 0.5, auto kc=0.006
pData(scan)[is.element(scan$scanID,c(toPlSx2$ID1[1],toPlSx2$ID2[1])),]
# one has family id 708
pData(scan)[is.element(scan$family,708),] # a pair of brothers, plus mom & dad
toPlSx2[is.element(toPlSx2$ID1,c(689417,691920,859470,988420)),]
toPlSx2[is.element(toPlSx2$ID2,c(689417,691920,859470,988420)),]
kcX_xAdj[is.element(kcX_xAdj$ID12,"691920 988420"),]
# what is xkc between random related 723230 and mother 689417?
kcX_xAdj[is.element(kcX_xAdj$ID1,689417)&is.element(kcX_xAdj$ID2,723230),]
kcX_xAdj[is.element(kcX_xAdj$ID1,691920)&is.element(kcX_xAdj$ID2,723230),]

kcX_xAdj$ID12 <- paste(kcX_xAdj$ID1,kcX_xAdj$ID2)
kcX_autoAdj$ID12 <- paste(kcX_autoAdj$ID1,kcX_autoAdj$ID2)

allKC <- merge(kcX_xAdj,kcX_autoAdj,by="ID12",all=TRUE,suffixes=c(".xAdj",".autoAdj"))
dim(allKC) # 81071011 9

allKC2 <- merge(allKC,kcAuto,by="ID12",all=TRUE)
head(allKC2); dim(allKC2) # 81952003 16

save(allKC2,file="olga_application/all_kc_estimates.RData")

### this code took too long to run...
# # now allKC2 holds all pairwise ibd results, using autos, using x chr adj for auto PC, using x chr adj for x chr PC
# allKC2 <- allKC2[,-c(2,3,6,7,10,11)]
# 
# # compare x chr pc adj for auto to x chr pc adj for x chr
# unrelX <- allKC2[allKC2$kin<0.025,] # this is now pairs of individs w auto KC < 0.025, ie unrel on autos
# dim(unrelX) # 78527239 16
# summary(unrelX$kin.autoAdj)
# summary(unrelX$kin.xAdj)
# summary(unrelX$kin)
# 
# # subset by xchrKin large
# unrelX <- unrelX[unrelX$kin.autoAdj>0.05|unrelX$kin.xAdj>0.05,]
# dim(unrelX) # 
# # so these are pairs that have auto KC < 0.025 but x KC > 0.05
# 
# pdf("olga_application/xKC_autoKC_autoKCunrel.pdf")
# plot(unrelX$kin,unrelX$kin.autoAdj,type="n",xlab="Autosomal KC",ylab="X Chr KC",
#      main="All Pairs of Samples with Auto KC < 0.025",ylim=c(0,1),xlim=c(-0.1,0.025))
# points(unrelX$kin,unrelX$kin.xAdj,pch=19,col="cyan")
# points(unrelX$kin,unrelX$kin.autoAdj,pch=19,col="magenta")
# legend("topleft",c("Adj for Auto PC 1-5","Adj for X chr PC 1-2"),col=c("magenta","cyan"),pch=19)
# dev.off()

rm(list=ls())


#####
# 13. PC-AiR on X chr SNPs adj for X chr KC (which is adj for auto ancestry)

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")
source("pcair_X_scan.include.R")
library(GWASTools)

scan <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_scanAnnot_v54_AMS.RData"))
head(pData(scan)) 

snp.pruned <- get(load("olga_application/snp_sel_xChr_ldPruned.RData"))
head(snp.pruned); length(snp.pruned) # 3600

# make genoData object
gdsobj <- GdsGenotypeReader("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/netCDF/subjects/OLGA_subj_geno.gds")
genoData <- GenotypeData(gdsobj,scanAnnot=scan[scan$subj.plink,])

unrel.set <- scan$scanID[scan$unrelated.pcair.deg4&is.na(scan$gengrp6.outliers)]
scanIncl <- scan$scanID[is.na(scan$gengrp6.outliers)&scan$geno.cntl==0&scan$subj.plink]

head(unrel.set); head(scanIncl)
length(unrel.set); length(scanIncl) # 10272 | 12747

kinMat <- get(load("olga_application/kcMat_xchr_autoPC15adj.RData"))
dim(kinMat) # 12734 12734

sum(is.na(kinMat)) # 0

scanIncl <- rownames(kinMat)

source("pcairPartition.R")
pc <- pcair(genoData,Xchr=TRUE,snp.include=snp.pruned,scan.include=scanIncl,kinMat=kinMat,kin.thresh=0.025)
# 24 samples unrelated at this thresh...hmm.
save(pc,file="olga_application/pca_prunedXsnps_unrelRel_autoPC15adj.RData")

rm(list=ls())


#####
# 14. Parse PC-AiR results on X chr after adj for X chr KC (adj for auto+x PC 1-5)

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")
library(GWASTools)

pc <- get(load("olga_application/pca_prunedXsnps_unrelRel_autoPC15adj.RData"))

dim(pc$vectors) # 12734 10
pc$values/pc$sum.values # 0.12556233 0.06968364 0.05613335 0.04533621

scan <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_scanAnnot_v54_AMS.RData"))
head(pData(scan)) 

unrel.set <- scan$scanID[scan$unrelated.pcair.deg4&is.na(scan$gengrp6.outliers)]
scanIncl <- scan$scanID[is.na(scan$gengrp6.outliers)&scan$geno.cntl==0&scan$subj.plink]

head(unrel.set); head(scanIncl)

scan <- scan[is.element(scan$scanID,rownames(pc$vectors)),]

pdf("olga_application/pca_X_adjXkc_adjAutoXPC15_ev12_col.pdf")
plot(-pc$vectors[,1],-pc$vectors[,2],xlab="EV1 (12.56%)",ylab="EV2 (6.97%)",type="n")
points(-pc$vectors[scan$race.cat=="Mexican",1],-pc$vectors[scan$race.cat=="Mexican",2],pch=20,col="goldenrod")
points(-pc$vectors[scan$race.cat=="CentralAmerican",1],-pc$vectors[scan$race.cat=="CentralAmerican",2],pch=20,col="red")
points(-pc$vectors[scan$race.cat=="SouthAmerican",1],-pc$vectors[scan$race.cat=="SouthAmerican",2],pch=20,col="magenta")
points(-pc$vectors[scan$race.cat=="PuertoRican",1],-pc$vectors[scan$race.cat=="PuertoRican",2],pch=20,col="green")
points(-pc$vectors[scan$race.cat=="Cuban",1],-pc$vectors[scan$race.cat=="Cuban",2],pch=20,col="blue")
points(-pc$vectors[scan$race.cat=="Dominican",1],-pc$vectors[scan$race.cat=="Dominican",2],pch=20,col="burlywood")
points(-pc$vectors[scan$race.cat=="Other",1],-pc$vectors[scan$race.cat=="Other",2],pch=20,col="black")
points(-pc$vectors[scan$race.cat=="Unknown",1],-pc$vectors[scan$race.cat=="Unknown",2],pch=20,col="black")
legend("topright",c("CentralAmerican","Cuban","Dominican","Mexican","PuertoRican","SouthAmerican","Other/Unknown"),pch=20,
       col=c("red","blue","burlywood","goldenrod","green","magenta","black"))
dev.off()
# OMG wtf???!!

which(-pc$vectors[,1]>0.6) # 630863
pData(scan)[is.element(scan$scanID,630863),] # male, cuban, unrel at all degrees

# do ev1 vs ev2 for all 7 entries, then a legend in the 8th panel
pdf("olga_application/pca_X_adjXkc_adjAutoXPC15_ev12_eachCol.pdf",width=14)
par(mfrow=c(2,4))
plot(-pc$vectors[scan$race.cat=="Mexican",1],-pc$vectors[scan$race.cat=="Mexican",2],pch=20,col="goldenrod",xlab="EV1",ylab="EV2",main="Mexican")
plot(-pc$vectors[scan$race.cat=="CentralAmerican",1],-pc$vectors[scan$race.cat=="CentralAmerican",2],pch=20,col="red",xlab="EV1",ylab="EV2",main="Central Am")
plot(-pc$vectors[scan$race.cat=="SouthAmerican",1],-pc$vectors[scan$race.cat=="SouthAmerican",2],pch=20,col="magenta",xlab="EV1",ylab="EV2",main="South Am")
plot(-pc$vectors[scan$race.cat=="PuertoRican",1],-pc$vectors[scan$race.cat=="PuertoRican",2],pch=20,col="green",xlab="EV1",ylab="EV2",main="Puerto Rican")
plot(-pc$vectors[scan$race.cat=="Cuban",1],-pc$vectors[scan$race.cat=="Cuban",2],pch=20,col="blue",xlab="EV1",ylab="EV2",main="Cuban")
plot(-pc$vectors[scan$race.cat=="Dominican",1],-pc$vectors[scan$race.cat=="Dominican",2],pch=20,col="burlywood",xlab="EV1",ylab="EV2",main="Dominican")
plot(-pc$vectors[is.element(scan$race.cat,c("Unknown","Other")),1],-pc$vectors[is.element(scan$race.cat,c("Unknown","Other")),2],pch=20,col="black",xlab="EV1",ylab="EV2",main="Other/Unk")
plot(1:10,1:10,axes=FALSE,xlab="",ylab="",type="n")
legend("topleft",c("CentralAmerican","Cuban","Dominican","Mexican","PuertoRican","SouthAmerican","Other/Unknown"),pch=20,
       col=c("red","blue","burlywood","goldenrod","green","magenta","black"),cex=1.5)
dev.off()

pdf("olga_application/pca_X_adjXkc_adjAutoXPC15_ev34_col.pdf")
plot(-pc$vectors[,3],pc$vectors[,4],xlab="EV3 (5.61%)",ylab="EV4 (4.53%)",type="n")
points(-pc$vectors[scan$race.cat=="Mexican",3],pc$vectors[scan$race.cat=="Mexican",4],pch=20,col="goldenrod")
points(-pc$vectors[scan$race.cat=="CentralAmerican",3],pc$vectors[scan$race.cat=="CentralAmerican",4],pch=20,col="red")
points(-pc$vectors[scan$race.cat=="SouthAmerican",3],pc$vectors[scan$race.cat=="SouthAmerican",4],pch=20,col="magenta")
points(-pc$vectors[scan$race.cat=="PuertoRican",3],pc$vectors[scan$race.cat=="PuertoRican",4],pch=20,col="green")
points(-pc$vectors[scan$race.cat=="Cuban",3],pc$vectors[scan$race.cat=="Cuban",4],pch=20,col="blue")
points(-pc$vectors[scan$race.cat=="Dominican",3],pc$vectors[scan$race.cat=="Dominican",4],pch=20,col="burlywood")
points(-pc$vectors[scan$race.cat=="Other",3],pc$vectors[scan$race.cat=="Other",4],pch=20,col="black")
points(-pc$vectors[scan$race.cat=="Unknown",3],pc$vectors[scan$race.cat=="Unknown",4],pch=20,col="black")
#legend("topleft",c("CentralAmerican","Cuban","Dominican","Mexican","PuertoRican","SouthAmerican","Other/Unknown"),pch=20,
#       col=c("red","blue","burlywood","goldenrod","green","magenta","black"))
dev.off()

## hmm. clearly there are families that are being singled out.
# how to change the kin.thresh so it doesn't include relatives? maybe try 0.2 and see what happens?

rm(list=ls())


#####
# 15. Run PC-AiR on 3500 autosomal SNPs

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x/")
source("pcair_X_scan.include.R")
library(GWASTools)

scan <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_scanAnnot_v54_AMS.RData"))
head(pData(scan)) 

snp.pruned <- get(load("olga_application/snp_sel_autoX_ldPruned.RData"))
head(snp.pruned); length(snp.pruned) # 158778

# make genoData object
gdsobj <- GdsGenotypeReader("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/netCDF/subjects/OLGA_subj_geno.gds")
genoData <- GenotypeData(gdsobj,scanAnnot=scan[scan$subj.plink,])

unrel.set <- scan$scanID[scan$unrelated.pcair.deg4&is.na(scan$gengrp6.outliers)]
scanIncl <- scan$scanID[is.na(scan$gengrp6.outliers)&scan$geno.cntl==0&scan$subj.plink]

head(unrel.set); head(scanIncl)
length(unrel.set); length(scanIncl) # 10272 | 12747

chr <- getChromosome(genoData)
snpID <- getSnpID(genoData)

idx <- which(is.element(snpID,snp.pruned))
head(snpID); head(snp.pruned)

table(chr[idx]) # so chr 22 has 2638; chr 20 has 4413
# use chr 20

chr.pruned <- chr[idx]
snp.pruned <- snp.pruned[chr.pruned==20]
length(snp.pruned) # 4413

pc <- pcair(genoData,unrel.set=unrel.set,Xchr=FALSE,snp.include=snp.pruned,scan.include=scanIncl)
save(pc,file="olga_application/pca_chr20snps_10272unrelPlusRel.RData")

rm(list=ls())


#####
# 16. PC-Relate on X chr SNP, corr for auto PC1-5

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x/")
library(GWASTools)
library(SNPRelate)
source("pcrelate_X.R")

autoPC <- get(load("olga_application/pca_prunedAutosnps_10272unrelPlusRel.RData"))

pcMat <- autoPC$eigenvect[,c(1:5)]
dim(pcMat) # 12747 5

scan <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_scanAnnot_v54_AMS.RData"))
head(pData(scan)) 

# make genoData object
gdsobj <- GdsGenotypeReader("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/netCDF/subjects/OLGA_subj_geno.gds")
genoData <- GenotypeData(gdsobj,scanAnnot=scan[scan$subj.plink,])

unrel.set <- scan$scanID[scan$unrelated.pcair.deg4&is.na(scan$gengrp6.outliers)]
scanIncl <- scan$scanID[is.na(scan$gengrp6.outliers)&scan$geno.cntl==0&scan$subj.plink]

head(unrel.set); head(scanIncl)
length(unrel.set); length(scanIncl) # 10272 | 12747

### need to exclude 13 people with an entirely missing x chr 
# from chromosome.anomalies.csv file, they are study samples with x chr anomaly
# thus, their entire x chr is filtered out. 

chromAnoms <- read.csv("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/dbGaP/To_dbGaP/chromosome_anomalies/chromosome_anomalies.csv",
                       header=TRUE,as.is=T)
ids <- chromAnoms$scanID[chromAnoms$chromosome=="X"&chromAnoms$filter&chromAnoms$whole.chr&
                           !is.element(substr(chromAnoms$SUBJECT_ID,1,2),"NA")]
length(ids) # 13
chromAnoms[is.element(chromAnoms$scanID,ids),] # perfect! exactly what we want

length(scanIncl) # 12747
scanIncl <- scanIncl[!is.element(scanIncl,ids)]
length(scanIncl) # 12734

# take these individs out of pcMat too
rownames(pcMat) <- autoPC$sampleInfo[,"sample.id"]
pcMat <- pcMat[is.element(rownames(pcMat),scanIncl),]
dim(pcMat) # 12734 2

rel <- pcrelate(genoData, unrel.set=unrel.set, scan.include=scanIncl, Xchr=TRUE, pcMat=pcMat)
names(rel)

save(rel,file="olga_application/pcRelate_Xchr_autoPC15adj.RData")

## make into matrix
tmp <- matrix(NA,nrow=12734,ncol=12734)
length(tmp[lower.tri(tmp,diag=FALSE)]) # 81071011

kc <- rel$kinship
length(kc[,"kin"]) # 81071011
tmp[lower.tri(tmp,diag=FALSE)] <- kc[,"kin"]

# mirror to upper triangle
tmp[upper.tri(tmp,diag=FALSE)] <- t(tmp)[upper.tri(tmp,diag=FALSE)]

rownames(tmp) <- scanIncl
colnames(tmp) <- scanIncl

diag(tmp) <- 0 # diag elements don't matter

save(tmp,file="olga_application/kcMat_xchr_autoPC15adj.RData")

rm(list=ls())


#####
# 17. Parse PC-AiR results on chr 20 SNPs

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")
library(GWASTools)

pc <- get(load("olga_application/pca_chr20snps_10272unrelPlusRel.RData"))

dim(pc$vectors) # 12747 10
pc$values/pc$sum.values # 0.030891699 0.015642867 0.004148236 0.004050736

scan <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_scanAnnot_v54_AMS.RData"))
head(pData(scan)) 

unrel.set <- scan$scanID[scan$unrelated.pcair.deg4&is.na(scan$gengrp6.outliers)]
scanIncl <- scan$scanID[is.na(scan$gengrp6.outliers)&scan$geno.cntl==0&scan$subj.plink]

head(unrel.set); head(scanIncl)

scan <- scan[is.element(scan$scanID,rownames(pc$vectors)),]

pdf("olga_application/pca_chr20_ev12_col.pdf")
plot(pc$vectors[,1],-pc$vectors[,2],xlab="EV1 (3.09%)",ylab="EV2 (1.56%)",type="n")
points(pc$vectors[scan$race.cat=="Mexican",1],-pc$vectors[scan$race.cat=="Mexican",2],pch=20,col="goldenrod")
points(pc$vectors[scan$race.cat=="CentralAmerican",1],-pc$vectors[scan$race.cat=="CentralAmerican",2],pch=20,col="red")
points(pc$vectors[scan$race.cat=="SouthAmerican",1],-pc$vectors[scan$race.cat=="SouthAmerican",2],pch=20,col="magenta")
points(pc$vectors[scan$race.cat=="PuertoRican",1],-pc$vectors[scan$race.cat=="PuertoRican",2],pch=20,col="green")
points(pc$vectors[scan$race.cat=="Cuban",1],-pc$vectors[scan$race.cat=="Cuban",2],pch=20,col="blue")
points(pc$vectors[scan$race.cat=="Dominican",1],-pc$vectors[scan$race.cat=="Dominican",2],pch=20,col="burlywood")
points(pc$vectors[scan$race.cat=="Other",1],-pc$vectors[scan$race.cat=="Other",2],pch=20,col="black")
points(pc$vectors[scan$race.cat=="Unknown",1],-pc$vectors[scan$race.cat=="Unknown",2],pch=20,col="black")
legend("topright",c("CentralAmerican","Cuban","Dominican","Mexican","PuertoRican","SouthAmerican","Other/Unknown"),pch=20,
       col=c("red","blue","burlywood","goldenrod","green","magenta","black"))
dev.off()

pdf("olga_application/pca_chr20_ev34_col.pdf")
plot(-pc$vectors[,3],pc$vectors[,4],xlab="EV3 (0.415%)",ylab="EV4 (0.405%)",type="n")
points(-pc$vectors[scan$race.cat=="Mexican",3],pc$vectors[scan$race.cat=="Mexican",4],pch=20,col="goldenrod")
points(-pc$vectors[scan$race.cat=="CentralAmerican",3],pc$vectors[scan$race.cat=="CentralAmerican",4],pch=20,col="red")
points(-pc$vectors[scan$race.cat=="SouthAmerican",3],pc$vectors[scan$race.cat=="SouthAmerican",4],pch=20,col="magenta")
points(-pc$vectors[scan$race.cat=="PuertoRican",3],pc$vectors[scan$race.cat=="PuertoRican",4],pch=20,col="green")
points(-pc$vectors[scan$race.cat=="Cuban",3],pc$vectors[scan$race.cat=="Cuban",4],pch=20,col="blue")
points(-pc$vectors[scan$race.cat=="Dominican",3],pc$vectors[scan$race.cat=="Dominican",4],pch=20,col="burlywood")
points(-pc$vectors[scan$race.cat=="Other",3],pc$vectors[scan$race.cat=="Other",4],pch=20,col="black")
points(-pc$vectors[scan$race.cat=="Unknown",3],pc$vectors[scan$race.cat=="Unknown",4],pch=20,col="black")
#legend("topleft",c("CentralAmerican","Cuban","Dominican","Mexican","PuertoRican","SouthAmerican","Other/Unknown"),pch=20,
#       col=c("red","blue","burlywood","goldenrod","green","magenta","black"))
dev.off()

pc_colors <- rep("black",nrow(scan))
pc_colors[scan$race.cat=="Mexican"] <- "goldenrod"
pc_colors[scan$race.cat=="CentralAmerican"] <- "red"
pc_colors[scan$race.cat=="SouthAmerican"] <- "magenta"
pc_colors[scan$race.cat=="PuertoRican"] <- "green"
pc_colors[scan$race.cat=="Cuban"] <- "blue"
pc_colors[scan$race.cat=="Dominican"] <- "burlywood"

library(MASS)
pdf("olga_application/pca_chr20_parCords.pdf",width=14)
#pc$vectors[,1] <- pc$vectors[,1]*-1
parcoord(pc$vectors,col=pc_colors)
dev.off()

# compare with x chr results

xchr <- get(load("olga_application/pca_prunedXsnps_10272unrel.RData"))

# subset pc$vectors to be samples that are included in OLGA results
all(is.element(rownames(pc$vectors),rownames(xchr$vectors))) # TRUE

cor(xchr$vectors[,1],-pc$vectors[,1]) # 0.7968513
pdf("olga_application/pca_chr20_x_ev1.pdf")
plot(-pc$vectors[,1],xchr$vectors[,1],xlab="Chr 20 EV1",ylab="X Chr EV1")
legend("topleft",c(expression(paste(rho,"= 0.797"))))
dev.off()

# color the points
pdf("olga_application/pca_chr20_x_ev1_col.pdf")
plot(-pc$vectors[,1]~xchr$vectors[,1],xlab="X Chr EV1",ylab="Chr 20 EV1",type="n")
points(-pc$vectors[scan$race.cat=="Mexican",1]~xchr$vectors[scan$race.cat=="Mexican",1],pch=20,col="goldenrod")
points(-pc$vectors[scan$race.cat=="CentralAmerican",1]~xchr$vectors[scan$race.cat=="CentralAmerican",1],pch=20,col="red")
points(-pc$vectors[scan$race.cat=="SouthAmerican",1]~xchr$vectors[scan$race.cat=="SouthAmerican",1],pch=20,col="magenta")
points(-pc$vectors[scan$race.cat=="PuertoRican",1]~xchr$vectors[scan$race.cat=="PuertoRican",1],pch=20,col="green")
points(-pc$vectors[scan$race.cat=="Cuban",1]~xchr$vectors[scan$race.cat=="Cuban",1],pch=20,col="blue")
points(-pc$vectors[scan$race.cat=="Dominican",1]~xchr$vectors[scan$race.cat=="Dominican",1],pch=20,col="burlywood")
points(-pc$vectors[scan$race.cat=="Other",1]~xchr$vectors[scan$race.cat=="Other",1],pch=20,col="black")
points(-pc$vectors[scan$race.cat=="Unknown",1]~xchr$vectors[scan$race.cat=="Unknown",1],pch=20,col="black")
legend("topleft",c("CentralAmerican","Cuban","Dominican","Mexican","PuertoRican","SouthAmerican","Other/Unknown"),pch=20,
       col=c("red","blue","burlywood","goldenrod","green","magenta","black"),cex=0.9)
legend("bottomright",c(expression(paste(rho,"= 0.797"))),bty="n")
dev.off()


##
# plot ev2 auto vs xchr ev2
# calculate corr
cor(xchr$vectors[,2],-pc$vectors[,2]) # 0.5950656

pdf("olga_application/pca_chr20_x_ev2_col.pdf")
plot(-pc$vectors[,2]~xchr$vectors[,2],xlab="X Chr EV2",ylab="Chr 20 EV2",type="n")
points(-pc$vectors[scan$race.cat=="Mexican",2]~xchr$vectors[scan$race.cat=="Mexican",2],pch=20,col="goldenrod")
points(-pc$vectors[scan$race.cat=="CentralAmerican",2]~xchr$vectors[scan$race.cat=="CentralAmerican",2],pch=20,col="red")
points(-pc$vectors[scan$race.cat=="SouthAmerican",2]~xchr$vectors[scan$race.cat=="SouthAmerican",2],pch=20,col="magenta")
points(-pc$vectors[scan$race.cat=="PuertoRican",2]~xchr$vectors[scan$race.cat=="PuertoRican",2],pch=20,col="green")
points(-pc$vectors[scan$race.cat=="Cuban",2]~xchr$vectors[scan$race.cat=="Cuban",2],pch=20,col="blue")
points(-pc$vectors[scan$race.cat=="Dominican",2]~xchr$vectors[scan$race.cat=="Dominican",2],pch=20,col="burlywood")
points(-pc$vectors[scan$race.cat=="Other",2]~xchr$vectors[scan$race.cat=="Other",2],pch=20,col="black")
points(-pc$vectors[scan$race.cat=="Unknown",2]~xchr$vectors[scan$race.cat=="Unknown",2],pch=20,col="black")
#legend("topleft",c("CentralAmerican","Cuban","Dominican","Mexican","PuertoRican","SouthAmerican","Other/Unknown"),pch=20,
#       col=c("red","blue","burlywood","goldenrod","green","magenta","black"),cex=0.9)
legend("bottomright",c(expression(paste(rho,"= 0.595"))),bty="n")
dev.off()

rm(list=ls())


#####
# 18. PC-AiR on X chr SNPs adj for X chr KC (which is adj for auto+x ancestry)

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")
source("pcair_X_scan.include.R")
library(GWASTools)

scan <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_scanAnnot_v54_AMS.RData"))
head(pData(scan))

snp.pruned <- get(load("olga_application/snp_sel_xChr_ldPruned.RData"))
head(snp.pruned); length(snp.pruned) # 3600

# make genoData object
gdsobj <- GdsGenotypeReader("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/netCDF/subjects/OLGA_subj_geno.gds")
genoData <- GenotypeData(gdsobj,scanAnnot=scan[scan$subj.plink,])

unrel.set <- scan$scanID[scan$unrelated.pcair.deg4&is.na(scan$gengrp6.outliers)]
scanIncl <- scan$scanID[is.na(scan$gengrp6.outliers)&scan$geno.cntl==0&scan$subj.plink]

head(unrel.set); head(scanIncl)
length(unrel.set); length(scanIncl) # 10272 | 12747

kinMat <- get(load("olga_application/kcMat_xchr_autoXPC15adj.RData"))
dim(kinMat) # 12734 12734

sum(is.na(kinMat)) # 0

scanIncl <- rownames(kinMat)

source("pcairPartition.R")
pc <- pcair(genoData,Xchr=TRUE,snp.include=snp.pruned,scan.include=scanIncl,kinMat=kinMat,kin.thresh=0.2)

save(pc,file="olga_application/pca_prunedXsnps_unrelRel_autoXPC15adj_kc2.RData")

rm(list=ls())


#####
# 19. Parse KC results, including #16. results too

library(GWASTools)
setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")

scan <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_scanAnnot_v54_AMS.RData"))
head(pData(scan))

kcX_autoAdj <- get(load("olga_application/pcRelate_Xchr_autoPC15adj.RData"))
kcX_autoAdj <- kcX_autoAdj$kinship
kcX_autoXAdj <- get(load("olga_application/pcRelate_Xchr_autoXPC15adj.RData"))
kcX_autoXAdj <- kcX_autoXAdj$kinship
kcX_xAdj <- get(load("olga_application/pcRelate_Xchr_xPC12adj.RData"))
kcX_xAdj <- kcX_xAdj$kinship
kcAuto <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/mconomos/results/relatedness/grm/grm_allchr.RData"))
kcAuto <- kcAuto$kinship

kcX_unadj <- get(load("olga_application/unadj_pcRelate_Xchr.RData"))
kcX_unadj <- kcX_unadj$kinship

range(kcX_autoAdj$kin) # ranges from -0.13 to 1.07
summary(kcX_unadj$kin) # ranges from -0.18 to 1.03
summary(kcX_autoXAdj$kin) # ranges from -0.32 to 0.98
summary(kcX_xAdj$kin) # ranges from -0.16 to 0.97

kcX_autoAdj[kcX_autoAdj$kin>0.95,]
# ID1    ID2  nsnp       kin
# 646259 669782 32830 0.9918486
# 660099 845060 34091 1.0132738
# 681836 683348 33857 0.9798943
# 684581 700299 34059 0.9774033
# 707996 875457 33913 0.9574957
# 757248 821317 33710 1.0468203
# 787326 893031 34157 1.0703116
# 799718 866519 34104 0.9586034
# 808658 971603 33236 0.9516704
# 883106 888083 33907 1.0323723
# 896943 901284 33458 0.9557623

kcX_xAdj[kcX_xAdj$kin>0.95,]
# ID1    ID2  nsnp       kin
# 646259 669782 31157 0.9656194
# 684581 700299 30186 0.9511280
# 703739 878890 25689 0.9679540
# 757248 821317 33791 0.9621796
# 799718 866519 30746 0.9611269
# 883106 888083 34155 0.9717185

kcX_autoXAdj[kcX_autoXAdj$kin>0.95,]
# ID1    ID2  nsnp       kin
# 707996 875457 31859 0.9566378
# 757248 821317 33427 0.9710367
# 799718 866519 32527 0.9586097
# 883106 888083 33550 0.9838231
# 896943 901284 28101 0.9513448
# 957205 972362 29310 0.9726953

pData(scan)[is.element(scan$scanID,c(646259,669782)),]
# woah! full brothers that have the same x chr from their mother!!

pData(scan)[is.element(scan$scanID,c(684581,700299)),] # same household but different family ids.
kcAuto[kcAuto$ID1==684581&kcAuto$ID2==700299,] # kc is 0.2886158, so brothers
kcX_autoAdj[kcX_autoAdj$ID1==684581&kcX_autoAdj$ID2==700299,] # kc is 0.9774
kcX_autoXAdj[kcX_autoXAdj$ID1==684581&kcX_autoXAdj$ID2==700299,] # kc is 0.9298

# why not same family id in the scan annot?

##
#kcAuto_unr <- kcAuto[kcAuto$kin<0.025,]
#dim(kcAuto_unr) # 78527239 7

# make hist of x chr kc estimates for pairs that are unrel on auto
unrel.set <- scan$scanID[scan$unrelated.pcair.deg4&is.na(scan$gengrp6.outliers)]

toHist <- kcX_autoXAdj[is.element(kcX_autoXAdj$ID1,unrel.set)&is.element(kcX_autoXAdj$ID2,unrel.set),]
dim(toHist) # 52649191 4
length(unique(toHist$ID1)) # 10261; great! length of unrel.set

toHistXadj <- kcX_xAdj[is.element(kcX_xAdj$ID1,unrel.set)&is.element(kcX_xAdj$ID2,unrel.set),]
dim(toHistXadj) # 52649191 4
length(unique(toHistXadj$ID1)) # 10261

toHistNoAdj <- kcX_unadj[is.element(kcX_unadj$ID1,unrel.set)&is.element(kcX_unadj$ID2,unrel.set),]
dim(toHistNoAdj) # 52751856 4
length(unique(toHistNoAdj$ID1)) # 10271

toHistAutoAdj <- kcX_autoAdj[is.element(kcX_autoAdj$ID1,unrel.set)&is.element(kcX_autoAdj$ID2,unrel.set),]
dim(toHistAutoAdj)
length(unique(toHistAutoAdj$ID1)) # 10261

pdf("olga_application/hist_xchrKC_autoPC15adj_autoUnrel.pdf")
hist(toHistAutoAdj$kin,breaks=50,main="X Chr KC for 12,734 Study Samples",
     xlab="X Chr KC adj For Auto PC 1-5",cex.lab=1.5,cex.main=1.5)
legend("topright",c(paste("mean =",format(mean(kcX_autoAdj$kin),digits=4)),paste("sd =",format(sd(kcX_autoAdj$kin),digits=4))),
       bty="n",cex=1.5)
dev.off()

pdf("olga_application/hist_xchrKC_autoPC15adj_autoUnrel_trunc.pdf")
hist(toHistAutoAdj$kin,main="X chr KC for auto unrel samples",breaks=50,
     xlab="X Chr KC adj For Auto PC 1-5",ylim=c(0,5000))
dev.off()

## overlay all these histograms on eachother
pdf("olga_application/hist_allXchrKC_overlaid.pdf",height=14,width=14)
hist(toHistAutoAdj$kin,main="X chr KC for auto unrel samples",breaks=50,
     col=rgb(1,1,0,0.5),xlab="X Chr KC",cex.lab=1.5,cex.main=1.5)
hist(toHistXadj$kin,add=TRUE,col=rgb(0,0,1,0.5),breaks=50,cex.lab=1.5)
hist(toHistNoAdj$kin,add=TRUE,col=rgb(0,1,0,0.5),breaks=50,cex.lab=1.5)
hist(toHist$kin,add=TRUE,col=rgb(1,0,0,0.5),breaks=50,cex.lab=1.5)
legend("topright",c("Auto PC 1-5 Adj","Auto+X PC 1-5 Adj","X Chr PC 1-2 Adj","No Adj"),
       col=c(rgb(1,1,0,0.5),rgb(1,0,0,0.5),rgb(0,0,1,0.5),rgb(0,1,0,0.5)),
       lty=1,lwd=5,cex=1.5)
dev.off()

pdf("olga_application/hist_allXchrKC.pdf",height=14,width=14)
par(mfrow=c(2,2))
hist(toHistAutoAdj$kin,main="",col=rgb(1,1,0,0.5),breaks=50,xlab="X Chr KC",cex.lab=1.5,xlim=c(-0.2,0.8))
legend("bottomright",c(paste("mean =",format(mean(toHistAutoAdj$kin),digits=4)),paste("sd =",format(sd(toHistAutoAdj$kin),digits=4))),
       bty="n",cex=1.5)
hist(toHist$kin,main="",breaks=50,cex.lab=1.5,col=rgb(1,0,0,0.5),xlab="X Chr KC",xlim=c(-0.2,0.8))
legend("bottomright",c(paste("mean =",format(mean(toHist$kin),digits=4)),paste("sd =",format(sd(toHist$kin),digits=4))),
       bty="n",cex=1.5)
hist(toHistXadj$kin,col=rgb(0,0,1,0.5),breaks=50,main="",xlab="X Chr KC",cex.lab=1.5,xlim=c(-0.2,0.8))
legend("bottomright",c(paste("mean =",format(mean(toHistXadj$kin),digits=4)),paste("sd =",format(sd(toHistXadj$kin),digits=4))),
       bty="n",cex=1.5)
hist(toHistNoAdj$kin,col=rgb(0,1,0,0.5),breaks=50,main="",xlab="X Chr KC",cex.lab=1.5,xlim=c(-0.2,0.8))
legend("bottomright",c(paste("mean =",format(mean(toHistNoAdj$kin,na.rm=T),digits=4)),paste("sd =",format(sd(toHistNoAdj$kin,na.rm=T),digits=4))),
       bty="n",cex=1.5)
legend("topright",c("Auto PC 1-5 Adj","Auto+X PC 1-5 Adj","X Chr PC 1-2 Adj","No Adj"),
       col=c(rgb(1,1,0,0.5),rgb(1,0,0,0.5),rgb(0,0,1,0.5),rgb(0,1,0,0.5)),
       lty=1,lwd=5,cex=1.5)
dev.off()


###
# plot x vs auto KC for parent-offspring pairs
# auto KC should be 1/4
obsRel <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/amstilp/results/ibd/v04_build37_study_control/ibd_obsrel.v4.RData"))
table(obsRel$obs.rel)
#Deg2  Deg3   Dup    FS    PO     U 
# 652   440 20350   773  1807   926 
poIds1 <- obsRel$ID1[obsRel$obs.rel=="PO"]
poIds2 <- obsRel$ID2[obsRel$obs.rel=="PO"]
obsRel$ID12 <- paste(obsRel$ID1,obsRel$ID2)

kcPO_autoAdj <- kcX_autoAdj[is.element(kcX_autoAdj$ID1,poIds1),]
obsPOids <- paste(kcPO_autoAdj$ID1,kcPO_autoAdj$ID2)
kcPO_autoAdj <- kcPO_autoAdj[is.element(obsPOids,paste(poIds1,poIds2)),]
dim(kcPO_autoAdj) # 1442 4
kcPO_autoAdj$ID12 <- paste(kcPO_autoAdj$ID1,kcPO_autoAdj$ID2)

kcPO_autoAdj <- merge(kcPO_autoAdj,obsRel[obsRel$obs.rel=="PO",],by="ID12",all.x=TRUE)
# need to merge in sex, too
kcPO_autoAdj <- merge(kcPO_autoAdj,pData(scan)[,c("scanID","sex")],by.x="ID1.x",by.y="scanID")
kcPO_autoAdj <- merge(kcPO_autoAdj,pData(scan)[,c("scanID","sex")],by.x="ID2.x",by.y="scanID",suffixes=c(".1",".2"))

summary(kcPO_autoAdj$kin) 
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-0.04058  0.29860  0.37140  0.39950  0.54530  0.84310
summary(kcPO_autoAdj$kinship)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.1787  0.2458  0.2483  0.2472  0.2501  0.2655

pdf("olga_application/kc_xvsAuto_poPairs_autoPC15adj.pdf")
plot(kcPO_autoAdj$kin,kcPO_autoAdj$kinship,type="n",xlab="X Chr KC",ylab="Auto KC",
     main="Auto vs X Chr KC for 1,442 PO pairs",cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
points(kcPO_autoAdj$kin[kcPO_autoAdj$sex.1=="M"&kcPO_autoAdj$sex.2=="M"],
       kcPO_autoAdj$kinship[kcPO_autoAdj$sex.1=="M"&kcPO_autoAdj$sex.2=="M"],col=rgb(0,1,0,0.5),pch=19)
points(kcPO_autoAdj$kin[kcPO_autoAdj$sex.1=="F"&kcPO_autoAdj$sex.2=="F"],
       kcPO_autoAdj$kinship[kcPO_autoAdj$sex.1=="F"&kcPO_autoAdj$sex.2=="F"],col=rgb(1,0,0,0.5),pch=19)
points(kcPO_autoAdj$kin[kcPO_autoAdj$sex.1=="M"&kcPO_autoAdj$sex.2=="F"],
       kcPO_autoAdj$kinship[kcPO_autoAdj$sex.1=="M"&kcPO_autoAdj$sex.2=="F"],col=rgb(0,0,1,0.5),pch=19)
points(kcPO_autoAdj$kin[kcPO_autoAdj$sex.1=="F"&kcPO_autoAdj$sex.2=="M"],
       kcPO_autoAdj$kinship[kcPO_autoAdj$sex.1=="F"&kcPO_autoAdj$sex.2=="M"],col=rgb(0,0,1,0.5),pch=19)
legend("bottomright",c("M-M pairs","F-M pairs","F-F pairs"),col=c(rgb(0,1,0,0.5),rgb(0,0,1,0.5),rgb(1,0,0,0.5)),
       pch=19,cex=1.3)
abline(v=0,col="gray"); abline(v=0.25,col="gray"); abline(v=0.5,col="gray")
dev.off()

# do same for unadj
kcPO_unadj <- kcX_unadj[is.element(kcX_unadj$ID1,poIds1),]
obsPOids <- paste(kcPO_unadj$ID1,kcPO_unadj$ID2)
kcPO_unadj <- kcPO_unadj[is.element(obsPOids,paste(poIds1,poIds2)),]
dim(kcPO_unadj) # 1445 4
kcPO_unadj$ID12 <- paste(kcPO_unadj$ID1,kcPO_unadj$ID2)

kcPO_unadj <- merge(kcPO_unadj,obsRel[obsRel$obs.rel=="PO",],by="ID12",all.x=TRUE)
# need to merge in sex, too
kcPO_unadj <- merge(kcPO_unadj,pData(scan)[,c("scanID","sex")],by.x="ID1.x",by.y="scanID")
kcPO_unadj <- merge(kcPO_unadj,pData(scan)[,c("scanID","sex")],by.x="ID2.x",by.y="scanID",suffixes=c(".1",".2"))

summary(kcPO_unadj$kin) 
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
#-0.08891  0.24920  0.32200  0.35040  0.49510  0.80580        3 

pdf("olga_application/kc_xvsAuto_poPairs_unadj.pdf")
plot(kcPO_unadj$kin,kcPO_unadj$kinship,type="n",xlab="X Chr KC",ylab="Auto KC",
     main="Auto vs X Chr KC for 1,442 PO pairs",cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
points(kcPO_unadj$kin[kcPO_unadj$sex.1=="M"&kcPO_unadj$sex.2=="M"],
       kcPO_unadj$kinship[kcPO_unadj$sex.1=="M"&kcPO_unadj$sex.2=="M"],col=rgb(0,1,0,0.5),pch=19)
points(kcPO_unadj$kin[kcPO_unadj$sex.1=="F"&kcPO_unadj$sex.2=="F"],
       kcPO_unadj$kinship[kcPO_unadj$sex.1=="F"&kcPO_unadj$sex.2=="F"],col=rgb(1,0,0,0.5),pch=19)
points(kcPO_unadj$kin[kcPO_unadj$sex.1=="M"&kcPO_unadj$sex.2=="F"],
       kcPO_unadj$kinship[kcPO_unadj$sex.1=="M"&kcPO_unadj$sex.2=="F"],col=rgb(0,0,1,0.5),pch=19)
points(kcPO_unadj$kin[kcPO_unadj$sex.1=="F"&kcPO_unadj$sex.2=="M"],
       kcPO_unadj$kinship[kcPO_unadj$sex.1=="F"&kcPO_unadj$sex.2=="M"],col=rgb(0,0,1,0.5),pch=19)
legend("bottomright",c("M-M pairs","F-M pairs","F-F pairs"),col=c(rgb(0,1,0,0.5),rgb(0,0,1,0.5),rgb(1,0,0,0.5)),
       pch=19,cex=1.3)
abline(v=0,col="gray"); abline(v=0.25,col="gray"); abline(v=0.5,col="gray")
dev.off()

# do same for x chr adj
kcPO_xAdj <- kcX_xAdj[is.element(kcX_xAdj$ID1,poIds1),]
obsPOids <- paste(kcPO_xAdj$ID1,kcPO_xAdj$ID2)
kcPO_xAdj <- kcPO_xAdj[is.element(obsPOids,paste(poIds1,poIds2)),]
dim(kcPO_xAdj) # 1442 4
kcPO_xAdj$ID12 <- paste(kcPO_xAdj$ID1,kcPO_xAdj$ID2)

kcPO_xAdj <- merge(kcPO_xAdj,obsRel[obsRel$obs.rel=="PO",],by="ID12",all.x=TRUE)
# need to merge in sex, too
kcPO_xAdj <- merge(kcPO_xAdj,pData(scan)[,c("scanID","sex")],by.x="ID1.x",by.y="scanID")
kcPO_xAdj <- merge(kcPO_xAdj,pData(scan)[,c("scanID","sex")],by.x="ID2.x",by.y="scanID",suffixes=c(".1",".2"))

summary(kcPO_xAdj$kin) 
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-0.08469  0.24940  0.27890  0.33700  0.50160  0.70920 

pdf("olga_application/kc_xvsAuto_poPairs_xAdj.pdf")
plot(kcPO_xAdj$kin,kcPO_xAdj$kinship,type="n",xlab="X Chr KC",ylab="Auto KC",
     main="Auto vs X Chr KC for 1,442 PO pairs",cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
points(kcPO_xAdj$kin[kcPO_xAdj$sex.1=="M"&kcPO_xAdj$sex.2=="M"],
       kcPO_xAdj$kinship[kcPO_xAdj$sex.1=="M"&kcPO_xAdj$sex.2=="M"],col=rgb(0,1,0,0.5),pch=19)
points(kcPO_xAdj$kin[kcPO_xAdj$sex.1=="F"&kcPO_xAdj$sex.2=="F"],
       kcPO_xAdj$kinship[kcPO_xAdj$sex.1=="F"&kcPO_xAdj$sex.2=="F"],col=rgb(1,0,0,0.5),pch=19)
points(kcPO_xAdj$kin[kcPO_xAdj$sex.1=="M"&kcPO_xAdj$sex.2=="F"],
       kcPO_xAdj$kinship[kcPO_xAdj$sex.1=="M"&kcPO_xAdj$sex.2=="F"],col=rgb(0,0,1,0.5),pch=19)
points(kcPO_xAdj$kin[kcPO_xAdj$sex.1=="F"&kcPO_xAdj$sex.2=="M"],
       kcPO_xAdj$kinship[kcPO_xAdj$sex.1=="F"&kcPO_xAdj$sex.2=="M"],col=rgb(0,0,1,0.5),pch=19)
legend("bottomright",c("M-M pairs","F-M pairs","F-F pairs"),col=c(rgb(0,1,0,0.5),rgb(0,0,1,0.5),rgb(1,0,0,0.5)),
       pch=19,cex=1.3)
abline(v=0,col="gray"); abline(v=0.25,col="gray"); abline(v=0.5,col="gray")
dev.off()


# do same for auto + x chr adj
kcPO_xAdj <- kcX_autoXAdj[is.element(kcX_autoXAdj$ID1,poIds1),]
obsPOids <- paste(kcPO_xAdj$ID1,kcPO_xAdj$ID2)
kcPO_xAdj <- kcPO_xAdj[is.element(obsPOids,paste(poIds1,poIds2)),]
dim(kcPO_xAdj) # 1442 4
kcPO_xAdj$ID12 <- paste(kcPO_xAdj$ID1,kcPO_xAdj$ID2)

kcPO_xAdj <- merge(kcPO_xAdj,obsRel[obsRel$obs.rel=="PO",],by="ID12",all.x=TRUE)
# need to merge in sex, too
kcPO_xAdj <- merge(kcPO_xAdj,pData(scan)[,c("scanID","sex")],by.x="ID1.x",by.y="scanID")
kcPO_xAdj <- merge(kcPO_xAdj,pData(scan)[,c("scanID","sex")],by.x="ID2.x",by.y="scanID",suffixes=c(".1",".2"))

summary(kcPO_xAdj$kin) 
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-0.07931  0.25380  0.29060  0.34300  0.50250  0.68550 

pdf("olga_application/kc_xvsAuto_poPairs_autoXAdj.pdf")
plot(kcPO_xAdj$kin,kcPO_xAdj$kinship,type="n",xlab="X Chr KC",ylab="Auto KC",
     main="Auto vs X Chr KC for 1,442 PO pairs",cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
points(kcPO_xAdj$kin[kcPO_xAdj$sex.1=="M"&kcPO_xAdj$sex.2=="M"],
       kcPO_xAdj$kinship[kcPO_xAdj$sex.1=="M"&kcPO_xAdj$sex.2=="M"],col=rgb(0,1,0,0.5),pch=19)
points(kcPO_xAdj$kin[kcPO_xAdj$sex.1=="F"&kcPO_xAdj$sex.2=="F"],
       kcPO_xAdj$kinship[kcPO_xAdj$sex.1=="F"&kcPO_xAdj$sex.2=="F"],col=rgb(1,0,0,0.5),pch=19)
points(kcPO_xAdj$kin[kcPO_xAdj$sex.1=="M"&kcPO_xAdj$sex.2=="F"],
       kcPO_xAdj$kinship[kcPO_xAdj$sex.1=="M"&kcPO_xAdj$sex.2=="F"],col=rgb(0,0,1,0.5),pch=19)
points(kcPO_xAdj$kin[kcPO_xAdj$sex.1=="F"&kcPO_xAdj$sex.2=="M"],
       kcPO_xAdj$kinship[kcPO_xAdj$sex.1=="F"&kcPO_xAdj$sex.2=="M"],col=rgb(0,0,1,0.5),pch=19)
legend("bottomright",c("M-M pairs","F-M pairs","F-F pairs"),col=c(rgb(0,1,0,0.5),rgb(0,0,1,0.5),rgb(1,0,0,0.5)),
       pch=19,cex=1.3)
abline(v=0,col="gray"); abline(v=0.25,col="gray"); abline(v=0.5,col="gray")
dev.off()

###
# now look at FS relationships
# auto KC should be 1/4
poIds1 <- obsRel$ID1[obsRel$obs.rel=="FS"]
poIds2 <- obsRel$ID2[obsRel$obs.rel=="FS"]

kcPO_autoAdj <- kcX_autoAdj[is.element(kcX_autoAdj$ID1,poIds1),]
obsPOids <- paste(kcPO_autoAdj$ID1,kcPO_autoAdj$ID2)
kcPO_autoAdj <- kcPO_autoAdj[is.element(obsPOids,paste(poIds1,poIds2)),]
dim(kcPO_autoAdj) # 696 4
kcPO_autoAdj$ID12 <- paste(kcPO_autoAdj$ID1,kcPO_autoAdj$ID2)

kcPO_autoAdj <- merge(kcPO_autoAdj,obsRel[obsRel$obs.rel=="FS",],by="ID12",all.x=TRUE)
# need to merge in sex, too
kcPO_autoAdj <- merge(kcPO_autoAdj,pData(scan)[,c("scanID","sex")],by.x="ID1.x",by.y="scanID")
kcPO_autoAdj <- merge(kcPO_autoAdj,pData(scan)[,c("scanID","sex")],by.x="ID2.x",by.y="scanID",suffixes=c(".1",".2"))

summary(kcPO_autoAdj$kin) 
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.03009 0.33820 0.41430 0.42750 0.48500 1.07000 
summary(kcPO_autoAdj$kinship)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.1818  0.2373  0.2521  0.2515  0.2652  0.3073 

pdf("olga_application/kc_xvsAuto_fsPairs_autoPC15adj.pdf")
plot(kcPO_autoAdj$kin,kcPO_autoAdj$kinship,type="n",xlab="X Chr KC",ylab="Auto KC",
     main="Auto vs X Chr KC for 696 FS pairs",cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
points(kcPO_autoAdj$kin[kcPO_autoAdj$sex.1=="M"&kcPO_autoAdj$sex.2=="M"],
       kcPO_autoAdj$kinship[kcPO_autoAdj$sex.1=="M"&kcPO_autoAdj$sex.2=="M"],col=rgb(0,1,0,0.5),pch=19)
points(kcPO_autoAdj$kin[kcPO_autoAdj$sex.1=="F"&kcPO_autoAdj$sex.2=="F"],
       kcPO_autoAdj$kinship[kcPO_autoAdj$sex.1=="F"&kcPO_autoAdj$sex.2=="F"],col=rgb(1,0,0,0.5),pch=19)
points(kcPO_autoAdj$kin[kcPO_autoAdj$sex.1=="M"&kcPO_autoAdj$sex.2=="F"],
       kcPO_autoAdj$kinship[kcPO_autoAdj$sex.1=="M"&kcPO_autoAdj$sex.2=="F"],col=rgb(0,0,1,0.5),pch=19)
points(kcPO_autoAdj$kin[kcPO_autoAdj$sex.1=="F"&kcPO_autoAdj$sex.2=="M"],
       kcPO_autoAdj$kinship[kcPO_autoAdj$sex.1=="F"&kcPO_autoAdj$sex.2=="M"],col=rgb(0,0,1,0.5),pch=19)
legend("bottomright",c("M-M pairs","F-M pairs","F-F pairs"),col=c(rgb(0,1,0,0.5),rgb(0,0,1,0.5),rgb(1,0,0,0.5)),
       pch=19,cex=1.3)
abline(v=(6/16),col="gray"); abline(v=0.25,col="gray"); abline(v=0.5,col="gray")
dev.off()


# do same for unadj
kcPO_unadj <- kcX_unadj[is.element(kcX_unadj$ID1,poIds1),]
obsPOids <- paste(kcPO_unadj$ID1,kcPO_unadj$ID2)
kcPO_unadj <- kcPO_unadj[is.element(obsPOids,paste(poIds1,poIds2)),]
dim(kcPO_unadj) # 697 4
kcPO_unadj$ID12 <- paste(kcPO_unadj$ID1,kcPO_unadj$ID2)

kcPO_unadj <- merge(kcPO_unadj,obsRel[obsRel$obs.rel=="FS",],by="ID12",all.x=TRUE)
# need to merge in sex, too
kcPO_unadj <- merge(kcPO_unadj,pData(scan)[,c("scanID","sex")],by.x="ID1.x",by.y="scanID")
kcPO_unadj <- merge(kcPO_unadj,pData(scan)[,c("scanID","sex")],by.x="ID2.x",by.y="scanID",suffixes=c(".1",".2"))

summary(kcPO_unadj$kin) 
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
#-0.01749  0.28920  0.36390  0.37830  0.43600  1.03100        1 

pdf("olga_application/kc_xvsAuto_fsPairs_unadj.pdf")
plot(kcPO_unadj$kin,kcPO_unadj$kinship,type="n",xlab="X Chr KC",ylab="Auto KC",
     main="Auto vs X Chr KC for 696 FS pairs",cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
points(kcPO_unadj$kin[kcPO_unadj$sex.1=="M"&kcPO_unadj$sex.2=="M"],
       kcPO_unadj$kinship[kcPO_unadj$sex.1=="M"&kcPO_unadj$sex.2=="M"],col=rgb(0,1,0,0.5),pch=19)
points(kcPO_unadj$kin[kcPO_unadj$sex.1=="F"&kcPO_unadj$sex.2=="F"],
       kcPO_unadj$kinship[kcPO_unadj$sex.1=="F"&kcPO_unadj$sex.2=="F"],col=rgb(1,0,0,0.5),pch=19)
points(kcPO_unadj$kin[kcPO_unadj$sex.1=="M"&kcPO_unadj$sex.2=="F"],
       kcPO_unadj$kinship[kcPO_unadj$sex.1=="M"&kcPO_unadj$sex.2=="F"],col=rgb(0,0,1,0.5),pch=19)
points(kcPO_unadj$kin[kcPO_unadj$sex.1=="F"&kcPO_unadj$sex.2=="M"],
       kcPO_unadj$kinship[kcPO_unadj$sex.1=="F"&kcPO_unadj$sex.2=="M"],col=rgb(0,0,1,0.5),pch=19)
legend("bottomright",c("M-M pairs","F-M pairs","F-F pairs"),col=c(rgb(0,1,0,0.5),rgb(0,0,1,0.5),rgb(1,0,0,0.5)),
       pch=19,cex=1.3)
abline(v=(6/16),col="gray"); abline(v=0.25,col="gray"); abline(v=0.5,col="gray")
dev.off()


# do same for x chr adj
kcPO_xAdj <- kcX_xAdj[is.element(kcX_xAdj$ID1,poIds1),]
obsPOids <- paste(kcPO_xAdj$ID1,kcPO_xAdj$ID2)
kcPO_xAdj <- kcPO_xAdj[is.element(obsPOids,paste(poIds1,poIds2)),]
dim(kcPO_xAdj) # 696 4
kcPO_xAdj$ID12 <- paste(kcPO_xAdj$ID1,kcPO_xAdj$ID2)

kcPO_xAdj <- merge(kcPO_xAdj,obsRel[obsRel$obs.rel=="FS",],by="ID12",all.x=TRUE)
# need to merge in sex, too
kcPO_xAdj <- merge(kcPO_xAdj,pData(scan)[,c("scanID","sex")],by.x="ID1.x",by.y="scanID")
kcPO_xAdj <- merge(kcPO_xAdj,pData(scan)[,c("scanID","sex")],by.x="ID2.x",by.y="scanID",suffixes=c(".1",".2"))

summary(kcPO_xAdj$kin) 
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.01809 0.28140 0.36220 0.36370 0.42290 0.97170

pdf("olga_application/kc_xvsAuto_fsPairs_xAdj.pdf")
plot(kcPO_xAdj$kin,kcPO_xAdj$kinship,type="n",xlab="X Chr KC",ylab="Auto KC",
     main="Auto vs X Chr KC for 696 FS pairs",cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
points(kcPO_xAdj$kin[kcPO_xAdj$sex.1=="M"&kcPO_xAdj$sex.2=="M"],
       kcPO_xAdj$kinship[kcPO_xAdj$sex.1=="M"&kcPO_xAdj$sex.2=="M"],col=rgb(0,1,0,0.5),pch=19)
points(kcPO_xAdj$kin[kcPO_xAdj$sex.1=="F"&kcPO_xAdj$sex.2=="F"],
       kcPO_xAdj$kinship[kcPO_xAdj$sex.1=="F"&kcPO_xAdj$sex.2=="F"],col=rgb(1,0,0,0.5),pch=19)
points(kcPO_xAdj$kin[kcPO_xAdj$sex.1=="M"&kcPO_xAdj$sex.2=="F"],
       kcPO_xAdj$kinship[kcPO_xAdj$sex.1=="M"&kcPO_xAdj$sex.2=="F"],col=rgb(0,0,1,0.5),pch=19)
points(kcPO_xAdj$kin[kcPO_xAdj$sex.1=="F"&kcPO_xAdj$sex.2=="M"],
       kcPO_xAdj$kinship[kcPO_xAdj$sex.1=="F"&kcPO_xAdj$sex.2=="M"],col=rgb(0,0,1,0.5),pch=19)
legend("bottomright",c("M-M pairs","F-M pairs","F-F pairs"),col=c(rgb(0,1,0,0.5),rgb(0,0,1,0.5),rgb(1,0,0,0.5)),
       pch=19,cex=1.3)
abline(v=(6/16),col="gray"); abline(v=0.25,col="gray"); abline(v=0.5,col="gray")
dev.off()


# do same for auto + x chr adj
kcPO_xAdj <- kcX_autoXAdj[is.element(kcX_autoXAdj$ID1,poIds1),]
obsPOids <- paste(kcPO_xAdj$ID1,kcPO_xAdj$ID2)
kcPO_xAdj <- kcPO_xAdj[is.element(obsPOids,paste(poIds1,poIds2)),]
dim(kcPO_xAdj) # 1442 4
kcPO_xAdj$ID12 <- paste(kcPO_xAdj$ID1,kcPO_xAdj$ID2)

kcPO_xAdj <- merge(kcPO_xAdj,obsRel[obsRel$obs.rel=="FS",],by="ID12",all.x=TRUE)
# need to merge in sex, too
kcPO_xAdj <- merge(kcPO_xAdj,pData(scan)[,c("scanID","sex")],by.x="ID1.x",by.y="scanID")
kcPO_xAdj <- merge(kcPO_xAdj,pData(scan)[,c("scanID","sex")],by.x="ID2.x",by.y="scanID",suffixes=c(".1",".2"))

summary(kcPO_xAdj$kin) 
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-0.01573  0.28340  0.36840  0.36890  0.43090  0.98380 

pdf("olga_application/kc_xvsAuto_fsPairs_autoXAdj.pdf")
plot(kcPO_xAdj$kin,kcPO_xAdj$kinship,type="n",xlab="X Chr KC",ylab="Auto KC",
     main="Auto vs X Chr KC for 696 FS pairs",cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
points(kcPO_xAdj$kin[kcPO_xAdj$sex.1=="M"&kcPO_xAdj$sex.2=="M"],
       kcPO_xAdj$kinship[kcPO_xAdj$sex.1=="M"&kcPO_xAdj$sex.2=="M"],col=rgb(0,1,0,0.5),pch=19)
points(kcPO_xAdj$kin[kcPO_xAdj$sex.1=="F"&kcPO_xAdj$sex.2=="F"],
       kcPO_xAdj$kinship[kcPO_xAdj$sex.1=="F"&kcPO_xAdj$sex.2=="F"],col=rgb(1,0,0,0.5),pch=19)
points(kcPO_xAdj$kin[kcPO_xAdj$sex.1=="M"&kcPO_xAdj$sex.2=="F"],
       kcPO_xAdj$kinship[kcPO_xAdj$sex.1=="M"&kcPO_xAdj$sex.2=="F"],col=rgb(0,0,1,0.5),pch=19)
points(kcPO_xAdj$kin[kcPO_xAdj$sex.1=="F"&kcPO_xAdj$sex.2=="M"],
       kcPO_xAdj$kinship[kcPO_xAdj$sex.1=="F"&kcPO_xAdj$sex.2=="M"],col=rgb(0,0,1,0.5),pch=19)
legend("bottomright",c("M-M pairs","F-M pairs","F-F pairs"),col=c(rgb(0,1,0,0.5),rgb(0,0,1,0.5),rgb(1,0,0,0.5)),
       pch=19,cex=1.3)
abline(v=0,col="gray"); abline(v=0.25,col="gray"); abline(v=0.5,col="gray")
dev.off()



###
# now look at Deg2 relationships
# auto KC should be 1/8
poIds1 <- obsRel$ID1[obsRel$obs.rel=="Deg2"]
poIds2 <- obsRel$ID2[obsRel$obs.rel=="Deg2"]

kcPO_autoAdj <- kcX_autoAdj[is.element(kcX_autoAdj$ID1,poIds1),]
obsPOids <- paste(kcPO_autoAdj$ID1,kcPO_autoAdj$ID2)
kcPO_autoAdj <- kcPO_autoAdj[is.element(obsPOids,paste(poIds1,poIds2)),]
dim(kcPO_autoAdj) # 564 4
kcPO_autoAdj$ID12 <- paste(kcPO_autoAdj$ID1,kcPO_autoAdj$ID2)

kcPO_autoAdj <- merge(kcPO_autoAdj,obsRel[obsRel$obs.rel=="Deg2",],by="ID12",all.x=TRUE)
# need to merge in sex, too
kcPO_autoAdj <- merge(kcPO_autoAdj,pData(scan)[,c("scanID","sex")],by.x="ID1.x",by.y="scanID")
kcPO_autoAdj <- merge(kcPO_autoAdj,pData(scan)[,c("scanID","sex")],by.x="ID2.x",by.y="scanID",suffixes=c(".1",".2"))

summary(kcPO_autoAdj$kin) 
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-0.03551  0.15740  0.25170  0.27900  0.36860  0.95750 
summary(kcPO_autoAdj$kinship)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.08917 0.11180 0.12210 0.12150 0.13080 0.16810 

pdf("olga_application/kc_xvsAuto_deg2Pairs_autoPC15adj.pdf")
plot(kcPO_autoAdj$kin,kcPO_autoAdj$kinship,type="n",xlab="X Chr KC",ylab="Auto KC",
     main="Auto vs X Chr KC for 564 Deg2 pairs",cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
points(kcPO_autoAdj$kin[kcPO_autoAdj$sex.1=="M"&kcPO_autoAdj$sex.2=="M"],
       kcPO_autoAdj$kinship[kcPO_autoAdj$sex.1=="M"&kcPO_autoAdj$sex.2=="M"],col=rgb(0,1,0,0.5),pch=19)
points(kcPO_autoAdj$kin[kcPO_autoAdj$sex.1=="F"&kcPO_autoAdj$sex.2=="F"],
       kcPO_autoAdj$kinship[kcPO_autoAdj$sex.1=="F"&kcPO_autoAdj$sex.2=="F"],col=rgb(1,0,0,0.5),pch=19)
points(kcPO_autoAdj$kin[kcPO_autoAdj$sex.1=="M"&kcPO_autoAdj$sex.2=="F"],
       kcPO_autoAdj$kinship[kcPO_autoAdj$sex.1=="M"&kcPO_autoAdj$sex.2=="F"],col=rgb(0,0,1,0.5),pch=19)
points(kcPO_autoAdj$kin[kcPO_autoAdj$sex.1=="F"&kcPO_autoAdj$sex.2=="M"],
       kcPO_autoAdj$kinship[kcPO_autoAdj$sex.1=="F"&kcPO_autoAdj$sex.2=="M"],col=rgb(0,0,1,0.5),pch=19)
legend("bottomright",c("M-M pairs","F-M pairs","F-F pairs"),col=c(rgb(0,1,0,0.5),rgb(0,0,1,0.5),rgb(1,0,0,0.5)),
       pch=19,cex=1.3)
abline(v=0,col="gray"); abline(v=0.25,col="gray"); abline(v=0.5,col="gray")
abline(v=(1/8),col="gray"); abline(v=(6/16),col="gray"); abline(v=(3/16),col="gray")
dev.off()


# do same for unadj
kcPO_unadj <- kcX_unadj[is.element(kcX_unadj$ID1,poIds1),]
obsPOids <- paste(kcPO_unadj$ID1,kcPO_unadj$ID2)
kcPO_unadj <- kcPO_unadj[is.element(obsPOids,paste(poIds1,poIds2)),]
dim(kcPO_unadj) # 565 4
kcPO_unadj$ID12 <- paste(kcPO_unadj$ID1,kcPO_unadj$ID2)

kcPO_unadj <- merge(kcPO_unadj,obsRel[obsRel$obs.rel=="Deg2",],by="ID12",all.x=TRUE)
# need to merge in sex, too
kcPO_unadj <- merge(kcPO_unadj,pData(scan)[,c("scanID","sex")],by.x="ID1.x",by.y="scanID")
kcPO_unadj <- merge(kcPO_unadj,pData(scan)[,c("scanID","sex")],by.x="ID2.x",by.y="scanID",suffixes=c(".1",".2"))

summary(kcPO_unadj$kin) 
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
#-0.08319  0.10890  0.20180  0.23010  0.31940  0.90940        1 

pdf("olga_application/kc_xvsAuto_deg2Pairs_unadj.pdf")
plot(kcPO_unadj$kin,kcPO_unadj$kinship,type="n",xlab="X Chr KC",ylab="Auto KC",
     main="Auto vs X Chr KC for 564 Deg2 pairs",cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
points(kcPO_unadj$kin[kcPO_unadj$sex.1=="M"&kcPO_unadj$sex.2=="M"],
       kcPO_unadj$kinship[kcPO_unadj$sex.1=="M"&kcPO_unadj$sex.2=="M"],col=rgb(0,1,0,0.5),pch=19)
points(kcPO_unadj$kin[kcPO_unadj$sex.1=="F"&kcPO_unadj$sex.2=="F"],
       kcPO_unadj$kinship[kcPO_unadj$sex.1=="F"&kcPO_unadj$sex.2=="F"],col=rgb(1,0,0,0.5),pch=19)
points(kcPO_unadj$kin[kcPO_unadj$sex.1=="M"&kcPO_unadj$sex.2=="F"],
       kcPO_unadj$kinship[kcPO_unadj$sex.1=="M"&kcPO_unadj$sex.2=="F"],col=rgb(0,0,1,0.5),pch=19)
points(kcPO_unadj$kin[kcPO_unadj$sex.1=="F"&kcPO_unadj$sex.2=="M"],
       kcPO_unadj$kinship[kcPO_unadj$sex.1=="F"&kcPO_unadj$sex.2=="M"],col=rgb(0,0,1,0.5),pch=19)
legend("bottomright",c("M-M pairs","F-M pairs","F-F pairs"),col=c(rgb(0,1,0,0.5),rgb(0,0,1,0.5),rgb(1,0,0,0.5)),
       pch=19,cex=1.3)
abline(v=0,col="gray"); abline(v=0.25,col="gray"); abline(v=0.5,col="gray")
abline(v=(1/8),col="gray"); abline(v=(6/16),col="gray"); abline(v=(3/16),col="gray")
dev.off()


# do same for x chr adj
kcPO_xAdj <- kcX_xAdj[is.element(kcX_xAdj$ID1,poIds1),]
obsPOids <- paste(kcPO_xAdj$ID1,kcPO_xAdj$ID2)
kcPO_xAdj <- kcPO_xAdj[is.element(obsPOids,paste(poIds1,poIds2)),]
dim(kcPO_xAdj) # 564 4
kcPO_xAdj$ID12 <- paste(kcPO_xAdj$ID1,kcPO_xAdj$ID2)

kcPO_xAdj <- merge(kcPO_xAdj,obsRel[obsRel$obs.rel=="Deg2",],by="ID12",all.x=TRUE)
# need to merge in sex, too
kcPO_xAdj <- merge(kcPO_xAdj,pData(scan)[,c("scanID","sex")],by.x="ID1.x",by.y="scanID")
kcPO_xAdj <- merge(kcPO_xAdj,pData(scan)[,c("scanID","sex")],by.x="ID2.x",by.y="scanID",suffixes=c(".1",".2"))

summary(kcPO_xAdj$kin) 
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-0.07308  0.08913  0.18550  0.21180  0.27860  0.96800 

pdf("olga_application/kc_xvsAuto_deg2Pairs_xAdj.pdf")
plot(kcPO_xAdj$kin,kcPO_xAdj$kinship,type="n",xlab="X Chr KC",ylab="Auto KC",
     main="Auto vs X Chr KC for 564 Deg2 pairs",cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
points(kcPO_xAdj$kin[kcPO_xAdj$sex.1=="M"&kcPO_xAdj$sex.2=="M"],
       kcPO_xAdj$kinship[kcPO_xAdj$sex.1=="M"&kcPO_xAdj$sex.2=="M"],col=rgb(0,1,0,0.5),pch=19)
points(kcPO_xAdj$kin[kcPO_xAdj$sex.1=="F"&kcPO_xAdj$sex.2=="F"],
       kcPO_xAdj$kinship[kcPO_xAdj$sex.1=="F"&kcPO_xAdj$sex.2=="F"],col=rgb(1,0,0,0.5),pch=19)
points(kcPO_xAdj$kin[kcPO_xAdj$sex.1=="M"&kcPO_xAdj$sex.2=="F"],
       kcPO_xAdj$kinship[kcPO_xAdj$sex.1=="M"&kcPO_xAdj$sex.2=="F"],col=rgb(0,0,1,0.5),pch=19)
points(kcPO_xAdj$kin[kcPO_xAdj$sex.1=="F"&kcPO_xAdj$sex.2=="M"],
       kcPO_xAdj$kinship[kcPO_xAdj$sex.1=="F"&kcPO_xAdj$sex.2=="M"],col=rgb(0,0,1,0.5),pch=19)
legend("bottomright",c("M-M pairs","F-M pairs","F-F pairs"),col=c(rgb(0,1,0,0.5),rgb(0,0,1,0.5),rgb(1,0,0,0.5)),
       pch=19,cex=1.3)
abline(v=0,col="gray"); abline(v=0.25,col="gray"); abline(v=0.5,col="gray")
abline(v=(1/8),col="gray"); abline(v=(6/16),col="gray"); abline(v=(3/16),col="gray")
dev.off()

# do same for auto + x chr adj
kcPO_xAdj <- kcX_autoXAdj[is.element(kcX_autoXAdj$ID1,poIds1),]
obsPOids <- paste(kcPO_xAdj$ID1,kcPO_xAdj$ID2)
kcPO_xAdj <- kcPO_xAdj[is.element(obsPOids,paste(poIds1,poIds2)),]
dim(kcPO_xAdj) # 1442 4
kcPO_xAdj$ID12 <- paste(kcPO_xAdj$ID1,kcPO_xAdj$ID2)

kcPO_xAdj <- merge(kcPO_xAdj,obsRel[obsRel$obs.rel=="Deg2",],by="ID12",all.x=TRUE)
# need to merge in sex, too
kcPO_xAdj <- merge(kcPO_xAdj,pData(scan)[,c("scanID","sex")],by.x="ID1.x",by.y="scanID")
kcPO_xAdj <- merge(kcPO_xAdj,pData(scan)[,c("scanID","sex")],by.x="ID2.x",by.y="scanID",suffixes=c(".1",".2"))

summary(kcPO_xAdj$kin) 
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-0.08066  0.09400  0.19030  0.21510  0.28700  0.95660 

pdf("olga_application/kc_xvsAuto_deg2Pairs_autoXAdj.pdf")
plot(kcPO_xAdj$kin,kcPO_xAdj$kinship,type="n",xlab="X Chr KC",ylab="Auto KC",
     main="Auto vs X Chr KC for 594 Deg2 pairs",cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
points(kcPO_xAdj$kin[kcPO_xAdj$sex.1=="M"&kcPO_xAdj$sex.2=="M"],
       kcPO_xAdj$kinship[kcPO_xAdj$sex.1=="M"&kcPO_xAdj$sex.2=="M"],col=rgb(0,1,0,0.5),pch=19)
points(kcPO_xAdj$kin[kcPO_xAdj$sex.1=="F"&kcPO_xAdj$sex.2=="F"],
       kcPO_xAdj$kinship[kcPO_xAdj$sex.1=="F"&kcPO_xAdj$sex.2=="F"],col=rgb(1,0,0,0.5),pch=19)
points(kcPO_xAdj$kin[kcPO_xAdj$sex.1=="M"&kcPO_xAdj$sex.2=="F"],
       kcPO_xAdj$kinship[kcPO_xAdj$sex.1=="M"&kcPO_xAdj$sex.2=="F"],col=rgb(0,0,1,0.5),pch=19)
points(kcPO_xAdj$kin[kcPO_xAdj$sex.1=="F"&kcPO_xAdj$sex.2=="M"],
       kcPO_xAdj$kinship[kcPO_xAdj$sex.1=="F"&kcPO_xAdj$sex.2=="M"],col=rgb(0,0,1,0.5),pch=19)
legend("bottomright",c("M-M pairs","F-M pairs","F-F pairs"),col=c(rgb(0,1,0,0.5),rgb(0,0,1,0.5),rgb(1,0,0,0.5)),
       pch=19,cex=1.3)
abline(v=0,col="gray"); abline(v=0.25,col="gray"); abline(v=0.5,col="gray")
abline(v=(1/8),col="gray"); abline(v=(6/16),col="gray"); abline(v=(3/16),col="gray")
dev.off()

rm(list=ls())


#####
# 20. Parse PC-AiR results on X chr SNPs adj for X chr KC (which is adj for auto ancestry)

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")
library(GWASTools)

pc <- get(load("olga_application/pca_prunedXsnps_unrelRel_autoXPC15adj_kc2.RData"))

dim(pc$vectors) # 12734 10
pc$values/pc$sum.values # 0.034522276 0.018280395 0.005324013 0.004744577

scan <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_scanAnnot_v54_AMS.RData"))
head(pData(scan)) 

unrel.set <- scan$scanID[scan$unrelated.pcair.deg4&is.na(scan$gengrp6.outliers)]
scanIncl <- scan$scanID[is.na(scan$gengrp6.outliers)&scan$geno.cntl==0&scan$subj.plink]

head(unrel.set); head(scanIncl)

scan <- scan[is.element(scan$scanID,rownames(pc$vectors)),]

pdf("olga_application/pca_X_adjXkc_adjAutoXPC15_kc2_ev12_col.pdf")
plot(-pc$vectors[,1],pc$vectors[,2],xlab="EV1 (3.45%)",ylab="EV2 (1.82%)",type="n")
points(-pc$vectors[scan$race.cat=="Mexican",1],pc$vectors[scan$race.cat=="Mexican",2],pch=20,col="goldenrod")
points(-pc$vectors[scan$race.cat=="CentralAmerican",1],pc$vectors[scan$race.cat=="CentralAmerican",2],pch=20,col="red")
points(-pc$vectors[scan$race.cat=="SouthAmerican",1],pc$vectors[scan$race.cat=="SouthAmerican",2],pch=20,col="magenta")
points(-pc$vectors[scan$race.cat=="PuertoRican",1],pc$vectors[scan$race.cat=="PuertoRican",2],pch=20,col="green")
points(-pc$vectors[scan$race.cat=="Cuban",1],pc$vectors[scan$race.cat=="Cuban",2],pch=20,col="blue")
points(-pc$vectors[scan$race.cat=="Dominican",1],pc$vectors[scan$race.cat=="Dominican",2],pch=20,col="burlywood")
points(-pc$vectors[scan$race.cat=="Other",1],pc$vectors[scan$race.cat=="Other",2],pch=20,col="black")
points(-pc$vectors[scan$race.cat=="Unknown",1],pc$vectors[scan$race.cat=="Unknown",2],pch=20,col="black")
legend("topright",c("CentralAmerican","Cuban","Dominican","Mexican","PuertoRican","SouthAmerican","Other/Unknown"),pch=20,
       col=c("red","blue","burlywood","goldenrod","green","magenta","black"))
dev.off()

pdf("olga_application/pca_X_adjXkc_adjAutoXPC15_kc2_ev34_col.pdf")
plot(-pc$vectors[,3],pc$vectors[,4],xlab="EV3 (5.61%)",ylab="EV4 (4.53%)",type="n")
points(-pc$vectors[scan$race.cat=="Mexican",3],pc$vectors[scan$race.cat=="Mexican",4],pch=20,col="goldenrod")
points(-pc$vectors[scan$race.cat=="CentralAmerican",3],pc$vectors[scan$race.cat=="CentralAmerican",4],pch=20,col="red")
points(-pc$vectors[scan$race.cat=="SouthAmerican",3],pc$vectors[scan$race.cat=="SouthAmerican",4],pch=20,col="magenta")
points(-pc$vectors[scan$race.cat=="PuertoRican",3],pc$vectors[scan$race.cat=="PuertoRican",4],pch=20,col="green")
points(-pc$vectors[scan$race.cat=="Cuban",3],pc$vectors[scan$race.cat=="Cuban",4],pch=20,col="blue")
points(-pc$vectors[scan$race.cat=="Dominican",3],pc$vectors[scan$race.cat=="Dominican",4],pch=20,col="burlywood")
points(-pc$vectors[scan$race.cat=="Other",3],pc$vectors[scan$race.cat=="Other",4],pch=20,col="black")
points(-pc$vectors[scan$race.cat=="Unknown",3],pc$vectors[scan$race.cat=="Unknown",4],pch=20,col="black")
#legend("topleft",c("CentralAmerican","Cuban","Dominican","Mexican","PuertoRican","SouthAmerican","Other/Unknown"),pch=20,
#       col=c("red","blue","burlywood","goldenrod","green","magenta","black"))
dev.off()

# looks pretty much the same as the previous results
rm(list=ls())


#####
# 21. PCA-SNP correlation plots

library(GWASTools)
library(SNPRelate)

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x/")

snp.pruned <- get(load("olga_application/snp_sel_xChr_ldPruned.RData"))
head(snp.pruned); length(snp.pruned) # 3600

pca <- get(load("olga_application/pca_prunedXsnps_10272unrel.RData"))
pcaObj <- list("sample.id"=rownames(pca$vectors),"snp.id"=snp.pruned,"eigenval"=pca$values,
               "eigenvect"=pca$vectors,"TraceXTX"=3191314629,"Bayesian"=FALSE, "genmat"=NULL)
class(pcaObj) <- "snpgdsPCAClass"

gdsobj <- snpgdsOpen("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/netCDF/subjects/OLGA_subj_geno.gds")
pca.corr <- snpgdsPCACorr(pcaObj, gdsobj, eig.which=1:6, num.thread=2)

save(pca.corr,file="olga_application/pca_corr_prunedXsnps_10272unrel.RData")

chr <- read.gdsn(index.gdsn(gdsobj, "snp.chromosome")) 

png("olga_application/pca_x_corrManh_ev1.png",width=720)
  plot(abs(pca.corr$snpcorr[1,]), ylim=c(0,1), xlab="", ylab="PC 1",
       col=chr, pch=19,cex=0.8,cex.lab=1.5)
dev.off()

png("olga_application/pca_x_corrManh_ev2.png",width=720)
  plot(abs(pca.corr$snpcorr[2,]), ylim=c(0,1), xlab="", ylab="PC 2",
       col=chr, pch=19,cex=0.8,cex.lab=1.5)
dev.off()

png("olga_application/pca_x_corrManh_ev3.png",width=720)
  plot(abs(pca.corr$snpcorr[3,]), ylim=c(0,1), xlab="", ylab="PC 3",
       col=chr, pch=19,cex=0.8,cex.lab=1.5)
dev.off()


png("olga_application/pca_x_corrManh_ev4.png",width=720)
  plot(abs(pca.corr$snpcorr[4,]), ylim=c(0,1), xlab="", ylab="PC 4",
       col=chr, pch=19,cex=0.8,cex.lab=1.5)
dev.off()

png("olga_application/pca_x_corrManh_ev5.png",width=720)
plot(abs(pca.corr$snpcorr[5,]),ylim=c(0,1),xlab="",ylab="PC 5",col=chr,pch=19,cex=0.8,cex.lab=1.5)
dev.off()

png("olga_application/pca_x_corrManh_ev6.png",width=720)
plot(abs(pca.corr$snpcorr[6,]),ylim=c(0,1),xlab="",ylab="PC 6",col=chr,pch=19,cex=0.8,cex.lab=1.5)
dev.off()

rm(list=ls())


#####
# 22. All PC-Relate runs again, using pruned X chr SNPs

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")
library(GWASTools)
source("pcrelate_X.R")
# use unrel.set to be the same unrelated set used for PCA on the X chr
# use scan.include to be the same set of individs used for PCA on the X chr

scan <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_scanAnnot_v54_AMS.RData"))
head(pData(scan)) 

# make genoData object
gdsobj <- GdsGenotypeReader("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/netCDF/subjects/OLGA_subj_geno.gds")
genoData <- GenotypeData(gdsobj,scanAnnot=scan[scan$subj.plink,])

unrel.set <- scan$scanID[scan$unrelated.pcair.deg4&is.na(scan$gengrp6.outliers)]
scanIncl <- scan$scanID[is.na(scan$gengrp6.outliers)&scan$geno.cntl==0&scan$subj.plink]

head(unrel.set); head(scanIncl)
length(unrel.set); length(scanIncl) # 10272 | 12747

### need to exclude 13 people with an entirely missing x chr 
# from chromosome.anomalies.csv file, they are study samples with x chr anomaly
# thus, their entire x chr is filtered out. 

chromAnoms <- read.csv("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/dbGaP/To_dbGaP/chromosome_anomalies/chromosome_anomalies.csv",
                       header=TRUE,as.is=T)
ids <- chromAnoms$scanID[chromAnoms$chromosome=="X"&chromAnoms$filter&chromAnoms$whole.chr&
                           !is.element(substr(chromAnoms$SUBJECT_ID,1,2),"NA")]
length(ids) # 13
chromAnoms[is.element(chromAnoms$scanID,ids),] # perfect! exactly what we want

length(scanIncl) # 12747
scanIncl <- scanIncl[!is.element(scanIncl,ids)]
length(scanIncl) # 12734

snp.pruned <- get(load("olga_application/snp_sel_xChr_ldPruned.RData"))
head(snp.pruned); length(snp.pruned) # 3600


###
# pcrelate unadj for ancestry

rel <- pcrelate(genoData, unrel.set=unrel.set, scan.include=scanIncl, Xchr=TRUE, snp.include=snp.pruned)
names(rel)

save(rel,file="olga_application/unadj_pcRelate_prunedXchr.RData")


###
# pcrelate adj for auto + x chr pcs 1-5

autoPC <- get(load("olga_application/pca_prunedAutoXsnps_10272unrelPlusRel.RData"))

pcMat <- autoPC$vectors[,c(1:5)]
dim(pcMat) # 12747 5

pcMat <- pcMat[is.element(rownames(pcMat),scanIncl),]
dim(pcMat) # 12734 5

rel <- pcrelate(genoData, unrel.set=unrel.set, scan.include=scanIncl, Xchr=TRUE, pcMat=pcMat, snp.include=snp.pruned)
names(rel)

save(rel,file="olga_application/pcRelate_prunedXchr_autoXPC15adj.RData")

## make into matrix
tmp <- matrix(NA,nrow=12734,ncol=12734)
length(tmp[lower.tri(tmp,diag=FALSE)]) # 81071011

kc <- rel$kinship
length(kc[,"kin"]) # 81071011
tmp[lower.tri(tmp,diag=FALSE)] <- kc[,"kin"]

# mirror to upper triangle
tmp[upper.tri(tmp,diag=FALSE)] <- t(tmp)[upper.tri(tmp,diag=FALSE)]

rownames(tmp) <- scanIncl
colnames(tmp) <- scanIncl

diag(tmp) <- 0 # diag elements don't matter

save(tmp,file="olga_application/kcMat_prunedxchr_autoXPC15adj.RData")


###
# pcrelate adj for x chr pcs 1-2

xPC <- get(load("olga_application/pca_prunedXsnps_10272unrel.RData"))

pcMat <- xPC$vectors[,c(1:2)]
dim(pcMat) # 12747 2

pcMat <- pcMat[is.element(rownames(pcMat),scanIncl),]
dim(pcMat) # 12734 2

rel <- pcrelate(genoData, unrel.set=unrel.set, scan.include=scanIncl, Xchr=TRUE, pcMat=pcMat, snp.include=snp.pruned)
names(rel)

save(rel,file="olga_application/pcRelate_prunedXchr_xPC12adj.RData")

## make into matrix
tmp <- matrix(NA,nrow=12734,ncol=12734)
length(tmp[lower.tri(tmp,diag=FALSE)]) # 81071011

kc <- rel$kinship
length(kc[,"kin"]) # 81071011
tmp[lower.tri(tmp,diag=FALSE)] <- kc[,"kin"]

# mirror to upper triangle
tmp[upper.tri(tmp,diag=FALSE)] <- t(tmp)[upper.tri(tmp,diag=FALSE)]

rownames(tmp) <- scanIncl
colnames(tmp) <- scanIncl

diag(tmp) <- 0 # diag elements don't matter

save(tmp,file="olga_application/kcMat_prunedxchr_xPC12adj.RData")

###
# pcrelate adj for auto pcs 1-5

autoPC <- get(load("olga_application/pca_prunedAutosnps_10272unrelPlusRel.RData"))

pcMat <- autoPC$eigenvect[,c(1:5)]
dim(pcMat) # 12747 5

rownames(pcMat) <- autoPC$sampleInfo[,"sample.id"]
pcMat <- pcMat[is.element(rownames(pcMat),scanIncl),]
dim(pcMat) # 12734 2

rel <- pcrelate(genoData, unrel.set=unrel.set, scan.include=scanIncl, Xchr=TRUE, pcMat=pcMat, snp.include=snp.pruned)
names(rel)

save(rel,file="olga_application/pcRelate_prunedXchr_autoPC15adj.RData")

## make into matrix
tmp <- matrix(NA,nrow=12734,ncol=12734)
length(tmp[lower.tri(tmp,diag=FALSE)]) # 81071011

kc <- rel$kinship
length(kc[,"kin"]) # 81071011
tmp[lower.tri(tmp,diag=FALSE)] <- kc[,"kin"]

# mirror to upper triangle
tmp[upper.tri(tmp,diag=FALSE)] <- t(tmp)[upper.tri(tmp,diag=FALSE)]

rownames(tmp) <- scanIncl
colnames(tmp) <- scanIncl

diag(tmp) <- 0 # diag elements don't matter

save(tmp,file="olga_application/kcMat_prunedxchr_autoPC15adj.RData")

rm(list=ls())


#####
# 23. Look at X-KC for some autosomal relationships, using KC from pruned X chr SNPs

library(GWASTools)
setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")

scan <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_scanAnnot_v54_AMS.RData"))
head(pData(scan))

kcX_autoAdj <- get(load("olga_application/pcRelate_prunedXchr_autoPC15adj.RData"))
kcX_autoAdj <- kcX_autoAdj$kinship
kcX_autoXAdj <- get(load("olga_application/pcRelate_prunedXchr_autoXPC15adj.RData"))
kcX_autoXAdj <- kcX_autoXAdj$kinship
kcX_xAdj <- get(load("olga_application/pcRelate_prunedXchr_xPC12adj.RData"))
kcX_xAdj <- kcX_xAdj$kinship
kcAuto <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/mconomos/results/relatedness/grm/grm_allchr.RData"))
kcAuto <- kcAuto$kinship

kcX_unadj <- get(load("olga_application/unadj_pcRelate_prunedXchr.RData"))
kcX_unadj <- kcX_unadj$kinship

# make hist of x chr kc estimates for pairs that are unrel on auto
unrel.set <- scan$scanID[scan$unrelated.pcair.deg4&is.na(scan$gengrp6.outliers)]

toHist <- kcX_autoXAdj[is.element(kcX_autoXAdj$ID1,unrel.set)&is.element(kcX_autoXAdj$ID2,unrel.set),]
dim(toHist) # 52649191 4
length(unique(toHist$ID1)) # 10261; great! length of unrel.set

toHistXadj <- kcX_xAdj[is.element(kcX_xAdj$ID1,unrel.set)&is.element(kcX_xAdj$ID2,unrel.set),]
dim(toHistXadj) # 52649191 4
length(unique(toHistXadj$ID1)) # 10261

toHistNoAdj <- kcX_unadj[is.element(kcX_unadj$ID1,unrel.set)&is.element(kcX_unadj$ID2,unrel.set),]
dim(toHistNoAdj) # 52751856 4
length(unique(toHistNoAdj$ID1)) # 10271

toHistAutoAdj <- kcX_autoAdj[is.element(kcX_autoAdj$ID1,unrel.set)&is.element(kcX_autoAdj$ID2,unrel.set),]
dim(toHistAutoAdj)
length(unique(toHistAutoAdj$ID1)) # 10261

pdf("olga_application/hist_xPrunedKC_autoPC15adj_autoUnrel.pdf")
hist(toHistAutoAdj$kin,breaks=50,main="X Chr KC for 12,734 Study Samples",
     xlab="X Chr KC adj For Auto PC 1-5",cex.lab=1.5,cex.main=1.5)
legend("topright",c(paste("mean =",format(mean(kcX_autoAdj$kin),digits=4)),paste("sd =",format(sd(kcX_autoAdj$kin),digits=4))),
       bty="n",cex=1.5)
dev.off()

pdf("olga_application/hist_xPrunedKC_autoPC15adj_autoUnrel_trunc.pdf")
hist(toHistAutoAdj$kin,main="X chr KC for auto unrel samples",breaks=50,
     xlab="X Chr KC adj For Auto PC 1-5",ylim=c(0,5000))
dev.off()

## overlay all these histograms on eachother
pdf("olga_application/hist_allXPrunedKC_overlaid.pdf",height=14,width=14)
hist(toHistAutoAdj$kin,main="X chr KC for auto unrel samples",breaks=50,
     col=rgb(1,1,0,0.5),xlab="X Chr KC",cex.lab=1.5,cex.main=1.5)
hist(toHistXadj$kin,add=TRUE,col=rgb(0,0,1,0.5),breaks=50,cex.lab=1.5)
hist(toHistNoAdj$kin,add=TRUE,col=rgb(0,1,0,0.5),breaks=50,cex.lab=1.5)
hist(toHist$kin,add=TRUE,col=rgb(1,0,0,0.5),breaks=50,cex.lab=1.5)
legend("topright",c("Auto PC 1-5 Adj","Auto+X PC 1-5 Adj","X Chr PC 1-2 Adj","No Adj"),
       col=c(rgb(1,1,0,0.5),rgb(1,0,0,0.5),rgb(0,0,1,0.5),rgb(0,1,0,0.5)),
       lty=1,lwd=5,cex=1.5)
dev.off()


pdf("olga_application/hist_allXPrunedKC_justautoXadj.pdf")
hist(toHist$kin,col=rgb(1,0,0,0.5),breaks=50,cex.lab=1.5,
     main="X chr KC for 10,261 Unrelated Samples",xlab="X Chr KC",cex.main=1.5)
dev.off()
pdf("olga_application/hist_allXPrunedKC_justautoXadj_trunc.pdf")
hist(toHist$kin,col=rgb(1,0,0,0.5),breaks=50,cex.lab=1.5,ylim=c(0,1e5),
     main="X chr KC for 10,261 Unrelated Samples",xlab="X Chr KC",cex.main=1.5)
dev.off()



pdf("olga_application/hist_allXPrunedKC.pdf",height=14,width=14)
par(mfrow=c(2,2))
hist(toHistAutoAdj$kin,main="",col=rgb(1,1,0,0.5),breaks=50,xlab="X Chr KC",cex.lab=1.5,xlim=c(-0.2,0.7))
legend("bottomright",c(paste("mean =",format(mean(toHistAutoAdj$kin),digits=4)),paste("sd =",format(sd(toHistAutoAdj$kin),digits=4))),
       bty="n",cex=1.5)
hist(toHist$kin,main="",breaks=50,cex.lab=1.5,col=rgb(1,0,0,0.5),xlab="X Chr KC",xlim=c(-0.2,0.7))
legend("bottomright",c(paste("mean =",format(mean(toHist$kin),digits=4)),paste("sd =",format(sd(toHist$kin),digits=4))),
       bty="n",cex=1.5)
hist(toHistXadj$kin,col=rgb(0,0,1,0.5),breaks=50,main="",xlab="X Chr KC",cex.lab=1.5,xlim=c(-0.2,0.7))
legend("bottomright",c(paste("mean =",format(mean(toHistXadj$kin),digits=4)),paste("sd =",format(sd(toHistXadj$kin),digits=4))),
       bty="n",cex=1.5)
hist(toHistNoAdj$kin,col=rgb(0,1,0,0.5),breaks=50,main="",xlab="X Chr KC",cex.lab=1.5,xlim=c(-0.2,0.7))
legend("bottomright",c(paste("mean =",format(mean(toHistNoAdj$kin,na.rm=T),digits=4)),paste("sd =",format(sd(toHistNoAdj$kin,na.rm=T),digits=4))),
       bty="n",cex=1.5)
legend("topright",c("Auto PC 1-5 Adj","Auto+X PC 1-5 Adj","X Chr PC 1-2 Adj","No Adj"),
       col=c(rgb(1,1,0,0.5),rgb(1,0,0,0.5),rgb(0,0,1,0.5),rgb(0,1,0,0.5)),
       lty=1,lwd=5,cex=1.5)
dev.off()


###
# plot x vs auto KC for parent-offspring pairs
# auto KC should be 1/4
obsRel <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/amstilp/results/ibd/v04_build37_study_control/ibd_obsrel.v4.RData"))
table(obsRel$obs.rel)
#Deg2  Deg3   Dup    FS    PO     U 
# 652   440 20350   773  1807   926 
poIds1 <- obsRel$ID1[obsRel$obs.rel=="PO"]
poIds2 <- obsRel$ID2[obsRel$obs.rel=="PO"]
obsRel$ID12 <- paste(obsRel$ID1,obsRel$ID2)

kcPO_autoAdj <- kcX_autoAdj[is.element(kcX_autoAdj$ID1,poIds1),]
obsPOids <- paste(kcPO_autoAdj$ID1,kcPO_autoAdj$ID2)
kcPO_autoAdj <- kcPO_autoAdj[is.element(obsPOids,paste(poIds1,poIds2)),]
dim(kcPO_autoAdj) # 1442 4
kcPO_autoAdj$ID12 <- paste(kcPO_autoAdj$ID1,kcPO_autoAdj$ID2)

kcPO_autoAdj <- merge(kcPO_autoAdj,obsRel[obsRel$obs.rel=="PO",],by="ID12",all.x=TRUE)
# need to merge in sex, too
kcPO_autoAdj <- merge(kcPO_autoAdj,pData(scan)[,c("scanID","sex")],by.x="ID1.x",by.y="scanID")
kcPO_autoAdj <- merge(kcPO_autoAdj,pData(scan)[,c("scanID","sex")],by.x="ID2.x",by.y="scanID",suffixes=c(".1",".2"))

summary(kcPO_autoAdj$kin) 
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-0.04063  0.27350  0.32960  0.36900  0.51470  0.75530 
summary(kcPO_autoAdj$kinship)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.1787  0.2458  0.2483  0.2472  0.2501  0.2655

pdf("olga_application/kc_xPrunedvsAuto_poPairs_autoPC15adj.pdf")
plot(kcPO_autoAdj$kin,kcPO_autoAdj$kinship,type="n",xlab="X Chr KC",ylab="Auto KC",
     main="Auto vs X Chr KC for 1,442 PO pairs",cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
points(kcPO_autoAdj$kin[kcPO_autoAdj$sex.1=="M"&kcPO_autoAdj$sex.2=="M"],
       kcPO_autoAdj$kinship[kcPO_autoAdj$sex.1=="M"&kcPO_autoAdj$sex.2=="M"],col=rgb(0,1,0,0.5),pch=19)
points(kcPO_autoAdj$kin[kcPO_autoAdj$sex.1=="F"&kcPO_autoAdj$sex.2=="F"],
       kcPO_autoAdj$kinship[kcPO_autoAdj$sex.1=="F"&kcPO_autoAdj$sex.2=="F"],col=rgb(1,0,0,0.5),pch=19)
points(kcPO_autoAdj$kin[kcPO_autoAdj$sex.1=="M"&kcPO_autoAdj$sex.2=="F"],
       kcPO_autoAdj$kinship[kcPO_autoAdj$sex.1=="M"&kcPO_autoAdj$sex.2=="F"],col=rgb(0,0,1,0.5),pch=19)
points(kcPO_autoAdj$kin[kcPO_autoAdj$sex.1=="F"&kcPO_autoAdj$sex.2=="M"],
       kcPO_autoAdj$kinship[kcPO_autoAdj$sex.1=="F"&kcPO_autoAdj$sex.2=="M"],col=rgb(0,0,1,0.5),pch=19)
legend("bottomright",c("M-M pairs","F-M pairs","F-F pairs"),col=c(rgb(0,1,0,0.5),rgb(0,0,1,0.5),rgb(1,0,0,0.5)),
       pch=19,cex=1.3)
abline(v=0,col="gray"); abline(v=0.25,col="gray"); abline(v=0.5,col="gray")
dev.off()

# do same for unadj
kcPO_unadj <- kcX_unadj[is.element(kcX_unadj$ID1,poIds1),]
obsPOids <- paste(kcPO_unadj$ID1,kcPO_unadj$ID2)
kcPO_unadj <- kcPO_unadj[is.element(obsPOids,paste(poIds1,poIds2)),]
dim(kcPO_unadj) # 1442 4
kcPO_unadj$ID12 <- paste(kcPO_unadj$ID1,kcPO_unadj$ID2)

kcPO_unadj <- merge(kcPO_unadj,obsRel[obsRel$obs.rel=="PO",],by="ID12",all.x=TRUE)
# need to merge in sex, too
kcPO_unadj <- merge(kcPO_unadj,pData(scan)[,c("scanID","sex")],by.x="ID1.x",by.y="scanID")
kcPO_unadj <- merge(kcPO_unadj,pData(scan)[,c("scanID","sex")],by.x="ID2.x",by.y="scanID",suffixes=c(".1",".2"))

summary(kcPO_unadj$kin) 
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-0.06461  0.24780  0.30540  0.34390  0.48800  0.74250 

pdf("olga_application/kc_xPrunedvsAuto_poPairs_unadj.pdf")
plot(kcPO_unadj$kin,kcPO_unadj$kinship,type="n",xlab="X Chr KC",ylab="Auto KC",
     main="Auto vs X Chr KC for 1,442 PO pairs",cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
points(kcPO_unadj$kin[kcPO_unadj$sex.1=="M"&kcPO_unadj$sex.2=="M"],
       kcPO_unadj$kinship[kcPO_unadj$sex.1=="M"&kcPO_unadj$sex.2=="M"],col=rgb(0,1,0,0.5),pch=19)
points(kcPO_unadj$kin[kcPO_unadj$sex.1=="F"&kcPO_unadj$sex.2=="F"],
       kcPO_unadj$kinship[kcPO_unadj$sex.1=="F"&kcPO_unadj$sex.2=="F"],col=rgb(1,0,0,0.5),pch=19)
points(kcPO_unadj$kin[kcPO_unadj$sex.1=="M"&kcPO_unadj$sex.2=="F"],
       kcPO_unadj$kinship[kcPO_unadj$sex.1=="M"&kcPO_unadj$sex.2=="F"],col=rgb(0,0,1,0.5),pch=19)
points(kcPO_unadj$kin[kcPO_unadj$sex.1=="F"&kcPO_unadj$sex.2=="M"],
       kcPO_unadj$kinship[kcPO_unadj$sex.1=="F"&kcPO_unadj$sex.2=="M"],col=rgb(0,0,1,0.5),pch=19)
legend("bottomright",c("M-M pairs","F-M pairs","F-F pairs"),col=c(rgb(0,1,0,0.5),rgb(0,0,1,0.5),rgb(1,0,0,0.5)),
       pch=19,cex=1.3)
abline(v=0,col="gray"); abline(v=0.25,col="gray"); abline(v=0.5,col="gray")
dev.off()

# do same for x chr adj
kcPO_xAdj <- kcX_xAdj[is.element(kcX_xAdj$ID1,poIds1),]
obsPOids <- paste(kcPO_xAdj$ID1,kcPO_xAdj$ID2)
kcPO_xAdj <- kcPO_xAdj[is.element(obsPOids,paste(poIds1,poIds2)),]
dim(kcPO_xAdj) # 1442 4
kcPO_xAdj$ID12 <- paste(kcPO_xAdj$ID1,kcPO_xAdj$ID2)

kcPO_xAdj <- merge(kcPO_xAdj,obsRel[obsRel$obs.rel=="PO",],by="ID12",all.x=TRUE)
# need to merge in sex, too
kcPO_xAdj <- merge(kcPO_xAdj,pData(scan)[,c("scanID","sex")],by.x="ID1.x",by.y="scanID")
kcPO_xAdj <- merge(kcPO_xAdj,pData(scan)[,c("scanID","sex")],by.x="ID2.x",by.y="scanID",suffixes=c(".1",".2"))

summary(kcPO_xAdj$kin) 
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-0.05793  0.24960  0.27390  0.33490  0.50080  0.68000 

pdf("olga_application/kc_xPrunedvsAuto_poPairs_xAdj.pdf")
plot(kcPO_xAdj$kin,kcPO_xAdj$kinship,type="n",xlab="X Chr KC",ylab="Auto KC",
     main="Auto vs X Chr KC for 1,442 PO pairs",cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
points(kcPO_xAdj$kin[kcPO_xAdj$sex.1=="M"&kcPO_xAdj$sex.2=="M"],
       kcPO_xAdj$kinship[kcPO_xAdj$sex.1=="M"&kcPO_xAdj$sex.2=="M"],col=rgb(0,1,0,0.5),pch=19)
points(kcPO_xAdj$kin[kcPO_xAdj$sex.1=="F"&kcPO_xAdj$sex.2=="F"],
       kcPO_xAdj$kinship[kcPO_xAdj$sex.1=="F"&kcPO_xAdj$sex.2=="F"],col=rgb(1,0,0,0.5),pch=19)
points(kcPO_xAdj$kin[kcPO_xAdj$sex.1=="M"&kcPO_xAdj$sex.2=="F"],
       kcPO_xAdj$kinship[kcPO_xAdj$sex.1=="M"&kcPO_xAdj$sex.2=="F"],col=rgb(0,0,1,0.5),pch=19)
points(kcPO_xAdj$kin[kcPO_xAdj$sex.1=="F"&kcPO_xAdj$sex.2=="M"],
       kcPO_xAdj$kinship[kcPO_xAdj$sex.1=="F"&kcPO_xAdj$sex.2=="M"],col=rgb(0,0,1,0.5),pch=19)
legend("bottomright",c("M-M pairs","F-F pairs","F-M pairs"),col=c(rgb(0,1,0,1),rgb(1,0,0,1),rgb(0,0,1,1)),
       pch=19,cex=1.3)
abline(v=0,col="gray"); abline(v=0.25,col="gray"); abline(v=0.5,col="gray")
dev.off()


# do same for auto + x chr adj
kcPO_xAdj <- kcX_autoXAdj[is.element(kcX_autoXAdj$ID1,poIds1),]
obsPOids <- paste(kcPO_xAdj$ID1,kcPO_xAdj$ID2)
kcPO_xAdj <- kcPO_xAdj[is.element(obsPOids,paste(poIds1,poIds2)),]
dim(kcPO_xAdj) # 1442 4
kcPO_xAdj$ID12 <- paste(kcPO_xAdj$ID1,kcPO_xAdj$ID2)

kcPO_xAdj <- merge(kcPO_xAdj,obsRel[obsRel$obs.rel=="PO",],by="ID12",all.x=TRUE)
# need to merge in sex, too
kcPO_xAdj <- merge(kcPO_xAdj,pData(scan)[,c("scanID","sex")],by.x="ID1.x",by.y="scanID")
kcPO_xAdj <- merge(kcPO_xAdj,pData(scan)[,c("scanID","sex")],by.x="ID2.x",by.y="scanID",suffixes=c(".1",".2"))

summary(kcPO_xAdj$kin) 
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-0.06666  0.25000  0.27750  0.33620  0.49990  0.67260 

pdf("olga_application/kc_xPrunedvsAuto_poPairs_autoXAdj.pdf")
plot(kcPO_xAdj$kin,kcPO_xAdj$kinship,type="n",xlab="X Chr KC",ylab="Auto KC",
     main="Auto vs X Chr KC for 1,442 PO pairs",cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
points(kcPO_xAdj$kin[kcPO_xAdj$sex.1=="M"&kcPO_xAdj$sex.2=="M"],
       kcPO_xAdj$kinship[kcPO_xAdj$sex.1=="M"&kcPO_xAdj$sex.2=="M"],col=rgb(0,1,0,0.5),pch=19)
points(kcPO_xAdj$kin[kcPO_xAdj$sex.1=="F"&kcPO_xAdj$sex.2=="F"],
       kcPO_xAdj$kinship[kcPO_xAdj$sex.1=="F"&kcPO_xAdj$sex.2=="F"],col=rgb(1,0,0,0.5),pch=19)
points(kcPO_xAdj$kin[kcPO_xAdj$sex.1=="M"&kcPO_xAdj$sex.2=="F"],
       kcPO_xAdj$kinship[kcPO_xAdj$sex.1=="M"&kcPO_xAdj$sex.2=="F"],col=rgb(0,0,1,0.5),pch=19)
points(kcPO_xAdj$kin[kcPO_xAdj$sex.1=="F"&kcPO_xAdj$sex.2=="M"],
       kcPO_xAdj$kinship[kcPO_xAdj$sex.1=="F"&kcPO_xAdj$sex.2=="M"],col=rgb(0,0,1,0.5),pch=19)
legend("bottomright",c("M-M pairs","F-M pairs","F-F pairs"),col=c(rgb(0,1,0,0.5),rgb(0,0,1,0.5),rgb(1,0,0,0.5)),
       pch=19,cex=1.3)
abline(v=0,col="gray"); abline(v=0.25,col="gray"); abline(v=0.5,col="gray")
dev.off()

###
# now look at FS relationships
# auto KC should be 1/4
poIds1 <- obsRel$ID1[obsRel$obs.rel=="FS"]
poIds2 <- obsRel$ID2[obsRel$obs.rel=="FS"]

kcPO_autoAdj <- kcX_autoAdj[is.element(kcX_autoAdj$ID1,poIds1),]
obsPOids <- paste(kcPO_autoAdj$ID1,kcPO_autoAdj$ID2)
kcPO_autoAdj <- kcPO_autoAdj[is.element(obsPOids,paste(poIds1,poIds2)),]
dim(kcPO_autoAdj) # 696 4
kcPO_autoAdj$ID12 <- paste(kcPO_autoAdj$ID1,kcPO_autoAdj$ID2)

kcPO_autoAdj <- merge(kcPO_autoAdj,obsRel[obsRel$obs.rel=="FS",],by="ID12",all.x=TRUE)
# need to merge in sex, too
kcPO_autoAdj <- merge(kcPO_autoAdj,pData(scan)[,c("scanID","sex")],by.x="ID1.x",by.y="scanID")
kcPO_autoAdj <- merge(kcPO_autoAdj,pData(scan)[,c("scanID","sex")],by.x="ID2.x",by.y="scanID",suffixes=c(".1",".2"))

summary(kcPO_autoAdj$kin) 
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.02164 0.31310 0.38850 0.39520 0.44920 1.05800 
summary(kcPO_autoAdj$kinship)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.1818  0.2373  0.2521  0.2515  0.2652  0.3073 

pdf("olga_application/kc_xPrunedvsAuto_fsPairs_autoPC15adj.pdf")
plot(kcPO_autoAdj$kin,kcPO_autoAdj$kinship,type="n",xlab="X Chr KC",ylab="Auto KC",
     main="Auto vs X Chr KC for 696 FS pairs",cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
points(kcPO_autoAdj$kin[kcPO_autoAdj$sex.1=="M"&kcPO_autoAdj$sex.2=="M"],
       kcPO_autoAdj$kinship[kcPO_autoAdj$sex.1=="M"&kcPO_autoAdj$sex.2=="M"],col=rgb(0,1,0,0.5),pch=19)
points(kcPO_autoAdj$kin[kcPO_autoAdj$sex.1=="F"&kcPO_autoAdj$sex.2=="F"],
       kcPO_autoAdj$kinship[kcPO_autoAdj$sex.1=="F"&kcPO_autoAdj$sex.2=="F"],col=rgb(1,0,0,0.5),pch=19)
points(kcPO_autoAdj$kin[kcPO_autoAdj$sex.1=="M"&kcPO_autoAdj$sex.2=="F"],
       kcPO_autoAdj$kinship[kcPO_autoAdj$sex.1=="M"&kcPO_autoAdj$sex.2=="F"],col=rgb(0,0,1,0.5),pch=19)
points(kcPO_autoAdj$kin[kcPO_autoAdj$sex.1=="F"&kcPO_autoAdj$sex.2=="M"],
       kcPO_autoAdj$kinship[kcPO_autoAdj$sex.1=="F"&kcPO_autoAdj$sex.2=="M"],col=rgb(0,0,1,0.5),pch=19)
legend("bottomright",c("M-M pairs","F-M pairs","F-F pairs"),col=c(rgb(0,1,0,0.5),rgb(0,0,1,0.5),rgb(1,0,0,0.5)),
       pch=19,cex=1.3)
abline(v=(6/16),col="gray"); abline(v=0.25,col="gray"); abline(v=0.5,col="gray")
dev.off()


# do same for unadj
kcPO_unadj <- kcX_unadj[is.element(kcX_unadj$ID1,poIds1),]
obsPOids <- paste(kcPO_unadj$ID1,kcPO_unadj$ID2)
kcPO_unadj <- kcPO_unadj[is.element(obsPOids,paste(poIds1,poIds2)),]
dim(kcPO_unadj) # 696 4
kcPO_unadj$ID12 <- paste(kcPO_unadj$ID1,kcPO_unadj$ID2)

kcPO_unadj <- merge(kcPO_unadj,obsRel[obsRel$obs.rel=="FS",],by="ID12",all.x=TRUE)
# need to merge in sex, too
kcPO_unadj <- merge(kcPO_unadj,pData(scan)[,c("scanID","sex")],by.x="ID1.x",by.y="scanID")
kcPO_unadj <- merge(kcPO_unadj,pData(scan)[,c("scanID","sex")],by.x="ID2.x",by.y="scanID",suffixes=c(".1",".2"))

summary(kcPO_unadj$kin) 
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#-0.004463  0.285900  0.362100  0.370100  0.422300  1.033000 

pdf("olga_application/kc_xPrunedvsAuto_fsPairs_unadj.pdf")
plot(kcPO_unadj$kin,kcPO_unadj$kinship,type="n",xlab="X Chr KC",ylab="Auto KC",
     main="Auto vs X Chr KC for 696 FS pairs",cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
points(kcPO_unadj$kin[kcPO_unadj$sex.1=="M"&kcPO_unadj$sex.2=="M"],
       kcPO_unadj$kinship[kcPO_unadj$sex.1=="M"&kcPO_unadj$sex.2=="M"],col=rgb(0,1,0,0.5),pch=19)
points(kcPO_unadj$kin[kcPO_unadj$sex.1=="F"&kcPO_unadj$sex.2=="F"],
       kcPO_unadj$kinship[kcPO_unadj$sex.1=="F"&kcPO_unadj$sex.2=="F"],col=rgb(1,0,0,0.5),pch=19)
points(kcPO_unadj$kin[kcPO_unadj$sex.1=="M"&kcPO_unadj$sex.2=="F"],
       kcPO_unadj$kinship[kcPO_unadj$sex.1=="M"&kcPO_unadj$sex.2=="F"],col=rgb(0,0,1,0.5),pch=19)
points(kcPO_unadj$kin[kcPO_unadj$sex.1=="F"&kcPO_unadj$sex.2=="M"],
       kcPO_unadj$kinship[kcPO_unadj$sex.1=="F"&kcPO_unadj$sex.2=="M"],col=rgb(0,0,1,0.5),pch=19)
legend("bottomright",c("M-M pairs","F-M pairs","F-F pairs"),col=c(rgb(0,1,0,0.5),rgb(0,0,1,0.5),rgb(1,0,0,0.5)),
       pch=19,cex=1.3)
abline(v=(6/16),col="gray"); abline(v=0.25,col="gray"); abline(v=0.5,col="gray")
dev.off()


# do same for x chr adj
kcPO_xAdj <- kcX_xAdj[is.element(kcX_xAdj$ID1,poIds1),]
obsPOids <- paste(kcPO_xAdj$ID1,kcPO_xAdj$ID2)
kcPO_xAdj <- kcPO_xAdj[is.element(obsPOids,paste(poIds1,poIds2)),]
dim(kcPO_xAdj) # 696 4
kcPO_xAdj$ID12 <- paste(kcPO_xAdj$ID1,kcPO_xAdj$ID2)

kcPO_xAdj <- merge(kcPO_xAdj,obsRel[obsRel$obs.rel=="FS",],by="ID12",all.x=TRUE)
# need to merge in sex, too
kcPO_xAdj <- merge(kcPO_xAdj,pData(scan)[,c("scanID","sex")],by.x="ID1.x",by.y="scanID")
kcPO_xAdj <- merge(kcPO_xAdj,pData(scan)[,c("scanID","sex")],by.x="ID2.x",by.y="scanID",suffixes=c(".1",".2"))

summary(kcPO_xAdj$kin) 
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.002151 0.280300 0.363600 0.361000 0.415000 0.971200 

pdf("olga_application/kc_xPrunedvsAuto_fsPairs_xAdj.pdf")
plot(kcPO_xAdj$kin,kcPO_xAdj$kinship,type="n",xlab="X Chr KC",ylab="Auto KC",
     main="Auto vs X Chr KC for 696 FS pairs",cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
points(kcPO_xAdj$kin[kcPO_xAdj$sex.1=="M"&kcPO_xAdj$sex.2=="M"],
       kcPO_xAdj$kinship[kcPO_xAdj$sex.1=="M"&kcPO_xAdj$sex.2=="M"],col=rgb(0,1,0,0.5),pch=19)
points(kcPO_xAdj$kin[kcPO_xAdj$sex.1=="F"&kcPO_xAdj$sex.2=="F"],
       kcPO_xAdj$kinship[kcPO_xAdj$sex.1=="F"&kcPO_xAdj$sex.2=="F"],col=rgb(1,0,0,0.5),pch=19)
points(kcPO_xAdj$kin[kcPO_xAdj$sex.1=="M"&kcPO_xAdj$sex.2=="F"],
       kcPO_xAdj$kinship[kcPO_xAdj$sex.1=="M"&kcPO_xAdj$sex.2=="F"],col=rgb(0,0,1,0.5),pch=19)
points(kcPO_xAdj$kin[kcPO_xAdj$sex.1=="F"&kcPO_xAdj$sex.2=="M"],
       kcPO_xAdj$kinship[kcPO_xAdj$sex.1=="F"&kcPO_xAdj$sex.2=="M"],col=rgb(0,0,1,0.5),pch=19)
legend("bottomright",c("M-M pairs","F-M pairs","F-F pairs"),col=c(rgb(0,1,0,0.5),rgb(0,0,1,0.5),rgb(1,0,0,0.5)),
       pch=19,cex=1.3)
abline(v=(6/16),col="gray"); abline(v=0.25,col="gray"); abline(v=0.5,col="gray")
dev.off()


# do same for auto + x chr adj
kcPO_xAdj <- kcX_autoXAdj[is.element(kcX_autoXAdj$ID1,poIds1),]
obsPOids <- paste(kcPO_xAdj$ID1,kcPO_xAdj$ID2)
kcPO_xAdj <- kcPO_xAdj[is.element(obsPOids,paste(poIds1,poIds2)),]
dim(kcPO_xAdj) # 1442 4
kcPO_xAdj$ID12 <- paste(kcPO_xAdj$ID1,kcPO_xAdj$ID2)

kcPO_xAdj <- merge(kcPO_xAdj,obsRel[obsRel$obs.rel=="FS",],by="ID12",all.x=TRUE)
# need to merge in sex, too
kcPO_xAdj <- merge(kcPO_xAdj,pData(scan)[,c("scanID","sex")],by.x="ID1.x",by.y="scanID")
kcPO_xAdj <- merge(kcPO_xAdj,pData(scan)[,c("scanID","sex")],by.x="ID2.x",by.y="scanID",suffixes=c(".1",".2"))

summary(kcPO_xAdj$kin) 
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-0.01616  0.28210  0.36620  0.36180  0.41490  0.98570 

pdf("olga_application/kc_xPrunedvsAuto_fsPairs_autoXAdj.pdf")
plot(kcPO_xAdj$kin,kcPO_xAdj$kinship,type="n",xlab="X Chr KC",ylab="Auto KC",
     main="Auto vs X Chr KC for 696 FS pairs",cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
points(kcPO_xAdj$kin[kcPO_xAdj$sex.1=="M"&kcPO_xAdj$sex.2=="M"],
       kcPO_xAdj$kinship[kcPO_xAdj$sex.1=="M"&kcPO_xAdj$sex.2=="M"],col=rgb(0,1,0,0.5),pch=19)
points(kcPO_xAdj$kin[kcPO_xAdj$sex.1=="F"&kcPO_xAdj$sex.2=="F"],
       kcPO_xAdj$kinship[kcPO_xAdj$sex.1=="F"&kcPO_xAdj$sex.2=="F"],col=rgb(1,0,0,0.5),pch=19)
points(kcPO_xAdj$kin[kcPO_xAdj$sex.1=="M"&kcPO_xAdj$sex.2=="F"],
       kcPO_xAdj$kinship[kcPO_xAdj$sex.1=="M"&kcPO_xAdj$sex.2=="F"],col=rgb(0,0,1,0.5),pch=19)
points(kcPO_xAdj$kin[kcPO_xAdj$sex.1=="F"&kcPO_xAdj$sex.2=="M"],
       kcPO_xAdj$kinship[kcPO_xAdj$sex.1=="F"&kcPO_xAdj$sex.2=="M"],col=rgb(0,0,1,0.5),pch=19)
legend("bottomright",c("M-M pairs","F-M pairs","F-F pairs"),col=c(rgb(0,1,0,0.5),rgb(0,0,1,0.5),rgb(1,0,0,0.5)),
       pch=19,cex=1.3)
abline(v=0,col="gray"); abline(v=0.25,col="gray"); abline(v=0.5,col="gray")
dev.off()



###
# now look at Deg2 relationships
# auto KC should be 1/8
poIds1 <- obsRel$ID1[obsRel$obs.rel=="Deg2"]
poIds2 <- obsRel$ID2[obsRel$obs.rel=="Deg2"]

kcPO_autoAdj <- kcX_autoAdj[is.element(kcX_autoAdj$ID1,poIds1),]
obsPOids <- paste(kcPO_autoAdj$ID1,kcPO_autoAdj$ID2)
kcPO_autoAdj <- kcPO_autoAdj[is.element(obsPOids,paste(poIds1,poIds2)),]
dim(kcPO_autoAdj) # 564 4
kcPO_autoAdj$ID12 <- paste(kcPO_autoAdj$ID1,kcPO_autoAdj$ID2)

kcPO_autoAdj <- merge(kcPO_autoAdj,obsRel[obsRel$obs.rel=="Deg2",],by="ID12",all.x=TRUE)
# need to merge in sex, too
kcPO_autoAdj <- merge(kcPO_autoAdj,pData(scan)[,c("scanID","sex")],by.x="ID1.x",by.y="scanID")
kcPO_autoAdj <- merge(kcPO_autoAdj,pData(scan)[,c("scanID","sex")],by.x="ID2.x",by.y="scanID",suffixes=c(".1",".2"))

summary(kcPO_autoAdj$kin) 
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-0.02011  0.13220  0.21850  0.24600  0.31870  0.92600 
summary(kcPO_autoAdj$kinship)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.08917 0.11180 0.12210 0.12150 0.13080 0.16810 

pdf("olga_application/kc_xPrunedvsAuto_deg2Pairs_autoPC15adj.pdf")
plot(kcPO_autoAdj$kin,kcPO_autoAdj$kinship,type="n",xlab="X Chr KC",ylab="Auto KC",
     main="Auto vs X Chr KC for 564 Deg2 pairs",cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
points(kcPO_autoAdj$kin[kcPO_autoAdj$sex.1=="M"&kcPO_autoAdj$sex.2=="M"],
       kcPO_autoAdj$kinship[kcPO_autoAdj$sex.1=="M"&kcPO_autoAdj$sex.2=="M"],col=rgb(0,1,0,0.5),pch=19)
points(kcPO_autoAdj$kin[kcPO_autoAdj$sex.1=="F"&kcPO_autoAdj$sex.2=="F"],
       kcPO_autoAdj$kinship[kcPO_autoAdj$sex.1=="F"&kcPO_autoAdj$sex.2=="F"],col=rgb(1,0,0,0.5),pch=19)
points(kcPO_autoAdj$kin[kcPO_autoAdj$sex.1=="M"&kcPO_autoAdj$sex.2=="F"],
       kcPO_autoAdj$kinship[kcPO_autoAdj$sex.1=="M"&kcPO_autoAdj$sex.2=="F"],col=rgb(0,0,1,0.5),pch=19)
points(kcPO_autoAdj$kin[kcPO_autoAdj$sex.1=="F"&kcPO_autoAdj$sex.2=="M"],
       kcPO_autoAdj$kinship[kcPO_autoAdj$sex.1=="F"&kcPO_autoAdj$sex.2=="M"],col=rgb(0,0,1,0.5),pch=19)
legend("bottomright",c("M-M pairs","F-M pairs","F-F pairs"),col=c(rgb(0,1,0,0.5),rgb(0,0,1,0.5),rgb(1,0,0,0.5)),
       pch=19,cex=1.3)
abline(v=0,col="gray"); abline(v=0.25,col="gray"); abline(v=0.5,col="gray")
abline(v=(1/8),col="gray"); abline(v=(6/16),col="gray"); abline(v=(3/16),col="gray")
dev.off()


# do same for unadj
kcPO_unadj <- kcX_unadj[is.element(kcX_unadj$ID1,poIds1),]
obsPOids <- paste(kcPO_unadj$ID1,kcPO_unadj$ID2)
kcPO_unadj <- kcPO_unadj[is.element(obsPOids,paste(poIds1,poIds2)),]
dim(kcPO_unadj) # 565 4
kcPO_unadj$ID12 <- paste(kcPO_unadj$ID1,kcPO_unadj$ID2)

kcPO_unadj <- merge(kcPO_unadj,obsRel[obsRel$obs.rel=="Deg2",],by="ID12",all.x=TRUE)
# need to merge in sex, too
kcPO_unadj <- merge(kcPO_unadj,pData(scan)[,c("scanID","sex")],by.x="ID1.x",by.y="scanID")
kcPO_unadj <- merge(kcPO_unadj,pData(scan)[,c("scanID","sex")],by.x="ID2.x",by.y="scanID",suffixes=c(".1",".2"))

summary(kcPO_unadj$kin) 
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.    
#-0.04523  0.10660  0.19420  0.22100  0.29550  0.89690 

pdf("olga_application/kc_xPrunedvsAuto_deg2Pairs_unadj.pdf")
plot(kcPO_unadj$kin,kcPO_unadj$kinship,type="n",xlab="X Chr KC",ylab="Auto KC",
     main="Auto vs X Chr KC for 564 Deg2 pairs",cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
points(kcPO_unadj$kin[kcPO_unadj$sex.1=="M"&kcPO_unadj$sex.2=="M"],
       kcPO_unadj$kinship[kcPO_unadj$sex.1=="M"&kcPO_unadj$sex.2=="M"],col=rgb(0,1,0,0.5),pch=19)
points(kcPO_unadj$kin[kcPO_unadj$sex.1=="F"&kcPO_unadj$sex.2=="F"],
       kcPO_unadj$kinship[kcPO_unadj$sex.1=="F"&kcPO_unadj$sex.2=="F"],col=rgb(1,0,0,0.5),pch=19)
points(kcPO_unadj$kin[kcPO_unadj$sex.1=="M"&kcPO_unadj$sex.2=="F"],
       kcPO_unadj$kinship[kcPO_unadj$sex.1=="M"&kcPO_unadj$sex.2=="F"],col=rgb(0,0,1,0.5),pch=19)
points(kcPO_unadj$kin[kcPO_unadj$sex.1=="F"&kcPO_unadj$sex.2=="M"],
       kcPO_unadj$kinship[kcPO_unadj$sex.1=="F"&kcPO_unadj$sex.2=="M"],col=rgb(0,0,1,0.5),pch=19)
legend("bottomright",c("M-M pairs","F-M pairs","F-F pairs"),col=c(rgb(0,1,0,0.5),rgb(0,0,1,0.5),rgb(1,0,0,0.5)),
       pch=19,cex=1.3)
abline(v=0,col="gray"); abline(v=0.25,col="gray"); abline(v=0.5,col="gray")
abline(v=(1/8),col="gray"); abline(v=(6/16),col="gray"); abline(v=(3/16),col="gray")
dev.off()


# do same for x chr adj
kcPO_xAdj <- kcX_xAdj[is.element(kcX_xAdj$ID1,poIds1),]
obsPOids <- paste(kcPO_xAdj$ID1,kcPO_xAdj$ID2)
kcPO_xAdj <- kcPO_xAdj[is.element(obsPOids,paste(poIds1,poIds2)),]
dim(kcPO_xAdj) # 564 4
kcPO_xAdj$ID12 <- paste(kcPO_xAdj$ID1,kcPO_xAdj$ID2)

kcPO_xAdj <- merge(kcPO_xAdj,obsRel[obsRel$obs.rel=="Deg2",],by="ID12",all.x=TRUE)
# need to merge in sex, too
kcPO_xAdj <- merge(kcPO_xAdj,pData(scan)[,c("scanID","sex")],by.x="ID1.x",by.y="scanID")
kcPO_xAdj <- merge(kcPO_xAdj,pData(scan)[,c("scanID","sex")],by.x="ID2.x",by.y="scanID",suffixes=c(".1",".2"))

summary(kcPO_xAdj$kin) 
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-0.06192  0.09534  0.18020  0.20820  0.27230  0.98130

pdf("olga_application/kc_xPrunedvsAuto_deg2Pairs_xAdj.pdf")
plot(kcPO_xAdj$kin,kcPO_xAdj$kinship,type="n",xlab="X Chr KC",ylab="Auto KC",
     main="Auto vs X Chr KC for 564 Deg2 pairs",cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
points(kcPO_xAdj$kin[kcPO_xAdj$sex.1=="M"&kcPO_xAdj$sex.2=="M"],
       kcPO_xAdj$kinship[kcPO_xAdj$sex.1=="M"&kcPO_xAdj$sex.2=="M"],col=rgb(0,1,0,0.5),pch=19)
points(kcPO_xAdj$kin[kcPO_xAdj$sex.1=="F"&kcPO_xAdj$sex.2=="F"],
       kcPO_xAdj$kinship[kcPO_xAdj$sex.1=="F"&kcPO_xAdj$sex.2=="F"],col=rgb(1,0,0,0.5),pch=19)
points(kcPO_xAdj$kin[kcPO_xAdj$sex.1=="M"&kcPO_xAdj$sex.2=="F"],
       kcPO_xAdj$kinship[kcPO_xAdj$sex.1=="M"&kcPO_xAdj$sex.2=="F"],col=rgb(0,0,1,0.5),pch=19)
points(kcPO_xAdj$kin[kcPO_xAdj$sex.1=="F"&kcPO_xAdj$sex.2=="M"],
       kcPO_xAdj$kinship[kcPO_xAdj$sex.1=="F"&kcPO_xAdj$sex.2=="M"],col=rgb(0,0,1,0.5),pch=19)
legend("bottomright",c("M-M pairs","F-M pairs","F-F pairs"),col=c(rgb(0,1,0,0.5),rgb(0,0,1,0.5),rgb(1,0,0,0.5)),
       pch=19,cex=1.3)
abline(v=0,col="gray"); abline(v=0.25,col="gray"); abline(v=0.5,col="gray")
abline(v=(1/8),col="gray"); abline(v=(6/16),col="gray"); abline(v=(3/16),col="gray")
dev.off()

# do same for auto + x chr adj
kcPO_xAdj <- kcX_autoXAdj[is.element(kcX_autoXAdj$ID1,poIds1),]
obsPOids <- paste(kcPO_xAdj$ID1,kcPO_xAdj$ID2)
kcPO_xAdj <- kcPO_xAdj[is.element(obsPOids,paste(poIds1,poIds2)),]
dim(kcPO_xAdj) # 564 4
kcPO_xAdj$ID12 <- paste(kcPO_xAdj$ID1,kcPO_xAdj$ID2)

kcPO_xAdj <- merge(kcPO_xAdj,obsRel[obsRel$obs.rel=="Deg2",],by="ID12",all.x=TRUE)
# need to merge in sex, too
kcPO_xAdj <- merge(kcPO_xAdj,pData(scan)[,c("scanID","sex")],by.x="ID1.x",by.y="scanID")
kcPO_xAdj <- merge(kcPO_xAdj,pData(scan)[,c("scanID","sex")],by.x="ID2.x",by.y="scanID",suffixes=c(".1",".2"))

summary(kcPO_xAdj$kin) 
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-0.06276  0.09521  0.18230  0.20790  0.27190  0.92820 

pdf("olga_application/kc_xPrunedvsAuto_deg2Pairs_autoXAdj.pdf")
plot(kcPO_xAdj$kin,kcPO_xAdj$kinship,type="n",xlab="X Chr KC",ylab="Auto KC",
     main="Auto vs X Chr KC for 594 Deg2 pairs",cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
points(kcPO_xAdj$kin[kcPO_xAdj$sex.1=="M"&kcPO_xAdj$sex.2=="M"],
       kcPO_xAdj$kinship[kcPO_xAdj$sex.1=="M"&kcPO_xAdj$sex.2=="M"],col=rgb(0,1,0,0.5),pch=19)
points(kcPO_xAdj$kin[kcPO_xAdj$sex.1=="F"&kcPO_xAdj$sex.2=="F"],
       kcPO_xAdj$kinship[kcPO_xAdj$sex.1=="F"&kcPO_xAdj$sex.2=="F"],col=rgb(1,0,0,0.5),pch=19)
points(kcPO_xAdj$kin[kcPO_xAdj$sex.1=="M"&kcPO_xAdj$sex.2=="F"],
       kcPO_xAdj$kinship[kcPO_xAdj$sex.1=="M"&kcPO_xAdj$sex.2=="F"],col=rgb(0,0,1,0.5),pch=19)
points(kcPO_xAdj$kin[kcPO_xAdj$sex.1=="F"&kcPO_xAdj$sex.2=="M"],
       kcPO_xAdj$kinship[kcPO_xAdj$sex.1=="F"&kcPO_xAdj$sex.2=="M"],col=rgb(0,0,1,0.5),pch=19)
legend("bottomright",c("M-M pairs","F-M pairs","F-F pairs"),col=c(rgb(0,1,0,0.5),rgb(0,0,1,0.5),rgb(1,0,0,0.5)),
       pch=19,cex=1.3)
abline(v=0,col="gray"); abline(v=0.25,col="gray"); abline(v=0.5,col="gray")
abline(v=(1/8),col="gray"); abline(v=(6/16),col="gray"); abline(v=(3/16),col="gray")
dev.off()

rm(list=ls())


#####
# 24. Compare KC X using all SNPs vs pruned SNPs

library(ggplot2)
library(GWASTools)
setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")

## try for auto adj
kcX_autoAdj <- get(load("olga_application/pcRelate_Xchr_autoPC15adj.RData"))
kcX_autoAdj <- kcX_autoAdj$kinship

kcX_P <- get(load("olga_application/pcRelate_prunedXchr_autoPC15adj.RData"))
kcX_P <- kcX_P$kinship

pdf("olga_application/kc_xPrunedvsNot_autoPC15adj.pdf")
df <- data.frame(x=kcX_P$kin,y=kcX_autoAdj$kin)
ggplot(df,aes(x=x,y=y))+stat_binhex()+geom_abline(slope=1,intercept=0)+labs(x="Pruned Set of X Chr SNPs")+
  labs(y="All X Chr SNPs")+labs(title="X KC Estimates using All X Chr SNPs vs Pruned Set of X Chr SNPs")
dev.off()
cor(df$x,df$y) # 0.9022165

## try for adj auto+x
kcX_autoXAdj <- get(load("olga_application/pcRelate_Xchr_autoXPC15adj.RData"))
kcX_autoXAdj <- kcX_autoXAdj$kinship

kcX_P <- get(load("olga_application/pcRelate_prunedXchr_autoXPC15adj.RData"))
kcX_P <- kcX_P$kinship
pdf("olga_application/kc_xPrunedvsNot_autoXPC15adj.pdf")
df <- data.frame(x=kcX_P$kin,y=kcX_autoXAdj$kin)
ggplot(df,aes(x=x,y=y))+stat_binhex()+geom_abline(slope=1,intercept=0)+labs(x="Pruned Set of X Chr SNPs")+
  labs(y="All X Chr SNPs")+labs(title="X KC Estimates using All X Chr SNPs vs Pruned Set of X Chr SNPs")
dev.off()
cor(df$x,df$y) # 0.7116353

## try for x adj  
kcX_xAdj <- get(load("olga_application/pcRelate_Xchr_xPC12adj.RData"))
kcX_xAdj <- kcX_xAdj$kinship

kcX_P <- get(load("olga_application/pcRelate_prunedXchr_xPC12adj.RData"))
kcX_P <- kcX_P$kinship
pdf("olga_application/kc_xPrunedvsNot_XAdj.pdf")
df <- data.frame(x=kcX_P$kin,y=kcX_xAdj$kin)
ggplot(df,aes(x=x,y=y))+stat_binhex()+geom_abline(slope=1,intercept=0)+labs(x="Pruned Set of X Chr SNPs")+
  labs(y="All X Chr SNPs")+labs(title="X KC Estimates using All X Chr SNPs vs Pruned Set of X Chr SNPs")
dev.off()
cor(df$x,df$y) # 0.6698762

## try for unadj
kcX_unadj <- get(load("olga_application/unadj_pcRelate_Xchr.RData"))
kcX_unadj <- kcX_unadj$kinship

kcX_P <- get(load("olga_application/unadj_pcRelate_prunedXchr.RData"))
kcX_P <- kcX_P$kinship

ids <- unique(kcX_P$ID1); length(ids) # 12733
kcX_unadj <- kcX_unadj[is.element(kcX_unadj$ID1,ids),]
kcX_unadj <- kcX_unadj[is.element(kcX_unadj$ID2,ids),]

kcX_P <- kcX_P[is.element(kcX_P$ID1,ids),]
kcX_P <- kcX_P[is.element(kcX_P$ID2,ids),]

df <- data.frame(x=kcX_P$kin,y=kcX_unadj$kin)

pdf("olga_application/kc_xPrunedvsNot_undj.pdf")
ggplot(df,aes(x=x,y=y))+stat_binhex()+geom_abline(slope=1,intercept=0)+labs(x="Pruned Set of X Chr SNPs")+
  labs(y="All X Chr SNPs")+labs(title="X KC Estimates using All X Chr SNPs vs Pruned Set of X Chr SNPs")
dev.off()
cor(df$x,df$y) # 0.9039277

rm(list=ls())


#####
# 25. PC-AiR on X chr SNPs adj for X chr KC (which is adj for auto ancestry)

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")
source("pcair_X_scan.include.R")
library(GWASTools)

scan <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_scanAnnot_v54_AMS.RData"))
head(pData(scan))

snp.pruned <- get(load("olga_application/snp_sel_xChr_ldPruned.RData"))
head(snp.pruned); length(snp.pruned) # 3600

# make genoData object
gdsobj <- GdsGenotypeReader("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/netCDF/subjects/OLGA_subj_geno.gds")
genoData <- GenotypeData(gdsobj,scanAnnot=scan[scan$subj.plink,])

unrel.set <- scan$scanID[scan$unrelated.pcair.deg4&is.na(scan$gengrp6.outliers)]
scanIncl <- scan$scanID[is.na(scan$gengrp6.outliers)&scan$geno.cntl==0&scan$subj.plink]

head(unrel.set); head(scanIncl)
length(unrel.set); length(scanIncl) # 10272 | 12747

kinMat <- get(load("olga_application/kcMat_xchr_autoPC15adj.RData"))
dim(kinMat) # 12734 12734

sum(is.na(kinMat)) # 0

scanIncl <- rownames(kinMat)

source("pcairPartition.R")
pc <- pcair(genoData,Xchr=TRUE,snp.include=snp.pruned,scan.include=scanIncl,kinMat=kinMat,kin.thresh=0.2)

save(pc,file="olga_application/pca_prunedXsnps_unrelRel_autoPC15adj_kc2.RData")

rm(list=ls())


#####
# 26. PC-AiR on X chr SNPs adj for X chr pruned KC (which is adj for auto ancestry)

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")
source("pcair_X_scan.include.R")
library(GWASTools)

scan <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_scanAnnot_v54_AMS.RData"))
head(pData(scan))

snp.pruned <- get(load("olga_application/snp_sel_xChr_ldPruned.RData"))
head(snp.pruned); length(snp.pruned) # 3600

# make genoData object
gdsobj <- GdsGenotypeReader("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/netCDF/subjects/OLGA_subj_geno.gds")
genoData <- GenotypeData(gdsobj,scanAnnot=scan[scan$subj.plink,])

unrel.set <- scan$scanID[scan$unrelated.pcair.deg4&is.na(scan$gengrp6.outliers)]
scanIncl <- scan$scanID[is.na(scan$gengrp6.outliers)&scan$geno.cntl==0&scan$subj.plink]

head(unrel.set); head(scanIncl)
length(unrel.set); length(scanIncl) # 10272 | 12747

kinMat <- get(load("olga_application/kcMat_prunedxchr_autoPC15adj.RData"))
dim(kinMat) # 12734 12734

sum(is.na(kinMat)) # 0

scanIncl <- rownames(kinMat)

source("pcairPartition.R")
pc <- pcair(genoData,Xchr=TRUE,snp.include=snp.pruned,scan.include=scanIncl,kinMat=kinMat,kin.thresh=0.2)

save(pc,file="olga_application/pca_prunedXsnps_unrelRel_autoPC15adj_xPrunedKC_kc2.RData")

rm(list=ls())


#####
# 27. Parse results for 25., 26.

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")
library(GWASTools)

pc <- get(load("olga_application/pca_prunedXsnps_unrelRel_autoPC15adj_xPrunedKC_kc2.RData"))

dim(pc$vectors) # 12734 10
pc$values/pc$sum.values # 0.031704634 0.013357990 0.005618091 0.005037939

scan <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_scanAnnot_v54_AMS.RData"))
head(pData(scan)) 

unrel.set <- scan$scanID[scan$unrelated.pcair.deg4&is.na(scan$gengrp6.outliers)]
scanIncl <- scan$scanID[is.na(scan$gengrp6.outliers)&scan$geno.cntl==0&scan$subj.plink]

head(unrel.set); head(scanIncl)

scan <- scan[is.element(scan$scanID,rownames(pc$vectors)),]

pdf("olga_application/pca_prunedX_adjXkc_adjAutoPC15_xPrunedKC_ev12_col.pdf")
plot(-pc$vectors[,1],-pc$vectors[,2],xlab="EV1 (3.17%)",ylab="EV2 (1.34%)",type="n")
points(-pc$vectors[scan$race.cat=="Mexican",1],-pc$vectors[scan$race.cat=="Mexican",2],pch=20,col="goldenrod")
points(-pc$vectors[scan$race.cat=="CentralAmerican",1],-pc$vectors[scan$race.cat=="CentralAmerican",2],pch=20,col="red")
points(-pc$vectors[scan$race.cat=="SouthAmerican",1],-pc$vectors[scan$race.cat=="SouthAmerican",2],pch=20,col="magenta")
points(-pc$vectors[scan$race.cat=="PuertoRican",1],-pc$vectors[scan$race.cat=="PuertoRican",2],pch=20,col="green")
points(-pc$vectors[scan$race.cat=="Cuban",1],-pc$vectors[scan$race.cat=="Cuban",2],pch=20,col="blue")
points(-pc$vectors[scan$race.cat=="Dominican",1],-pc$vectors[scan$race.cat=="Dominican",2],pch=20,col="burlywood")
points(-pc$vectors[scan$race.cat=="Other",1],-pc$vectors[scan$race.cat=="Other",2],pch=20,col="black")
points(-pc$vectors[scan$race.cat=="Unknown",1],-pc$vectors[scan$race.cat=="Unknown",2],pch=20,col="black")
legend("topright",c("CentralAmerican","Cuban","Dominican","Mexican","PuertoRican","SouthAmerican","Other/Unknown"),pch=20,
       col=c("red","blue","burlywood","goldenrod","green","magenta","black"))
dev.off()

# do ev1 vs ev2 for all 7 entries, then a legend in the 8th panel
pdf("olga_application/pca_prunedX_adjXkc_adjAutoPC15_xPrunedKC_ev12_eachCol.pdf",width=14)
par(mfrow=c(2,4))
plot(-pc$vectors[scan$race.cat=="Mexican",1],-pc$vectors[scan$race.cat=="Mexican",2],pch=20,col="goldenrod",xlab="EV1",ylab="EV2",main="Mexican")
plot(-pc$vectors[scan$race.cat=="CentralAmerican",1],-pc$vectors[scan$race.cat=="CentralAmerican",2],pch=20,col="red",xlab="EV1",ylab="EV2",main="Central Am")
plot(-pc$vectors[scan$race.cat=="SouthAmerican",1],-pc$vectors[scan$race.cat=="SouthAmerican",2],pch=20,col="magenta",xlab="EV1",ylab="EV2",main="South Am")
plot(-pc$vectors[scan$race.cat=="PuertoRican",1],-pc$vectors[scan$race.cat=="PuertoRican",2],pch=20,col="green",xlab="EV1",ylab="EV2",main="Puerto Rican")
plot(-pc$vectors[scan$race.cat=="Cuban",1],-pc$vectors[scan$race.cat=="Cuban",2],pch=20,col="blue",xlab="EV1",ylab="EV2",main="Cuban")
plot(-pc$vectors[scan$race.cat=="Dominican",1],-pc$vectors[scan$race.cat=="Dominican",2],pch=20,col="burlywood",xlab="EV1",ylab="EV2",main="Dominican")
plot(-pc$vectors[is.element(scan$race.cat,c("Unknown","Other")),1],-pc$vectors[is.element(scan$race.cat,c("Unknown","Other")),2],pch=20,col="black",xlab="EV1",ylab="EV2",main="Other/Unk")
plot(1:10,1:10,axes=FALSE,xlab="",ylab="",type="n")
legend("topleft",c("CentralAmerican","Cuban","Dominican","Mexican","PuertoRican","SouthAmerican","Other/Unknown"),pch=20,
       col=c("red","blue","burlywood","goldenrod","green","magenta","black"),cex=1.5)
dev.off()

pdf("olga_application/pca_prunedX_adjXkc_adjAutoPC15_xPrunedKC_ev34_col.pdf")
plot(-pc$vectors[,3],pc$vectors[,4],xlab="EV3 (0.56%)",ylab="EV4 (0.50%)",type="n")
points(-pc$vectors[scan$race.cat=="Mexican",3],pc$vectors[scan$race.cat=="Mexican",4],pch=20,col="goldenrod")
points(-pc$vectors[scan$race.cat=="CentralAmerican",3],pc$vectors[scan$race.cat=="CentralAmerican",4],pch=20,col="red")
points(-pc$vectors[scan$race.cat=="SouthAmerican",3],pc$vectors[scan$race.cat=="SouthAmerican",4],pch=20,col="magenta")
points(-pc$vectors[scan$race.cat=="PuertoRican",3],pc$vectors[scan$race.cat=="PuertoRican",4],pch=20,col="green")
points(-pc$vectors[scan$race.cat=="Cuban",3],pc$vectors[scan$race.cat=="Cuban",4],pch=20,col="blue")
points(-pc$vectors[scan$race.cat=="Dominican",3],pc$vectors[scan$race.cat=="Dominican",4],pch=20,col="burlywood")
points(-pc$vectors[scan$race.cat=="Other",3],pc$vectors[scan$race.cat=="Other",4],pch=20,col="black")
points(-pc$vectors[scan$race.cat=="Unknown",3],pc$vectors[scan$race.cat=="Unknown",4],pch=20,col="black")
#legend("topleft",c("CentralAmerican","Cuban","Dominican","Mexican","PuertoRican","SouthAmerican","Other/Unknown"),pch=20,
#       col=c("red","blue","burlywood","goldenrod","green","magenta","black"))
dev.off()

rm(list=ls())


#####
# 28. Do eigen() on output from 22.

library(GWASTools)
setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")
source("pcairPartition.R")

kinRes <- get(load("olga_application/pcRelate_prunedX_autoPC15adj_noCorrect.RData"))
inbd <- kinRes$inbreed

kc_corr <- get(load("olga_application/pcRelate_prunedXchr_autoPC15adj.RData"))

kinMat <- get(load("olga_application/kcMat_prunedx_autoPC15adj_noCorrect.RData"))
dim(kinMat) # 12734 12734

sum(is.na(kinMat)) # 0

diag(kinMat) <- 1+inbd$f

mns <- rowMeans(kinMat)
summary(mns); head(mns) # good! around zero
#Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
#-5.201e-04 -1.343e-04  5.343e-05  7.192e-05  2.624e-04  7.445e-04 

# remove relatives from this matrix, do eigen() on it
unrel <- pcairPartition(kinMat,kin.thresh=0.2)
length(unrel) # 2
names(unrel) # rels | unrels

save(unrel,file="olga_application/partition_kcMat_prunedx_autoPC15adj_noCorrect_thresh2.RData")

unrel.set <- unrel$unrels
length(unrel.set) # 10509

unrelkcT <- kinMat[is.element(rownames(kinMat),unrel.set),]
unrelkc <- unrelkcT[,is.element(colnames(unrelkcT),unrel.set)]
dim(unrelkc) # 10509 10509
allequal(rownames(unrelkc),colnames(unrelkc)) # TRUE

rSs <- rowSums(unrelkc) 
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#-59.650  -4.792   5.922   2.634  15.190  33.210 

mns <- rowMeans(unrelkc)
#Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
#-0.0056760 -0.0004560  0.0005635  0.0002506  0.0014450  0.0031600 

res <- eigen(unrelkc)
names(res) # values vectors
dim(res$vectors) # 10509 10509

save(res,file="olga_application/eigen_unrel_adjXkc_adjAutoPC15_xPrunedKC_noCorrect_thresh2.RData")

pdf("olga_application/eigen_unrel_adjXkc_adjAutoPC15_xPrunedKC_ev12.pdf")
plot(res$vectors[,1],res$vectors[,2],xlab="EV1",ylab="EV2")
dev.off()

scan <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_scanAnnot_v54_AMS.RData"))
scan <- scan[is.element(scan$scanID,rownames(unrelkc)),]
res$values[1:5]/sum(res$values) # 0.017248231 0.007532178 0.003074856 0.002710641 0.002284822

plotPC <- function(col1,col2,xlab,ylab){
  plot(res$vectors[,col1],res$vectors[,col2],xlab=xlab,ylab=ylab,type="n",cex.lab=1.5,cex.axis=1.5)
  points(res$vectors[scan$race.cat=="Mexican",col1],res$vectors[scan$race.cat=="Mexican",col2],
         pch=20,col="goldenrod")
  points(res$vectors[scan$race.cat=="CentralAmerican",col1],res$vectors[scan$race.cat=="CentralAmerican",col2],pch=20,col="red")
  points(res$vectors[scan$race.cat=="SouthAmerican",col1],res$vectors[scan$race.cat=="SouthAmerican",col2],pch=20,col="magenta")
  points(res$vectors[scan$race.cat=="PuertoRican",col1],res$vectors[scan$race.cat=="PuertoRican",col2],pch=20,col="green")
  points(res$vectors[scan$race.cat=="Cuban",col1],res$vectors[scan$race.cat=="Cuban",col2],pch=20,col="blue")
  points(res$vectors[scan$race.cat=="Dominican",col1],res$vectors[scan$race.cat=="Dominican",col2],pch=20,col="burlywood")
  points(res$vectors[scan$race.cat=="Other",col1],res$vectors[scan$race.cat=="Other",col2],pch=20,col="black")
  points(res$vectors[scan$race.cat=="Unknown",col1],res$vectors[scan$race.cat=="Unknown",col2],pch=20,col="black")
  legend("bottomleft",c("CentralAmerican","Cuban","Dominican","Mexican","PuertoRican","SouthAmerican","Other/Unknown"),pch=20,
         col=c("red","blue","burlywood","goldenrod","green","magenta","black"))
}
  
pdf("olga_application/eigen_unrel_adjXkc_adjAutoPC15_xPrunedKC_ev12_col.pdf")
plotPC(1,2,xlab="EV1 (1.72%)",ylab="EV2 (0.75%)")
dev.off()

png("olga_application/eigen_unrel_adjXkc_adjAutoPC15_xPrunedKC_ev12_col.png")
plotPC(1,2,xlab="EV1 (1.72%)",ylab="EV2 (0.75%)")
dev.off()


pdf("olga_application/eigen_unrel_adjXkc_adjAutoPC15_xPrunedKC_pairs14.pdf")
pairs(res$vectors[,1:4],res$vectors[,1:4])#,text.panel=text(0.5, 0.5, "EV1"))
dev.off()

# do ev1 vs ev2 for all 7 entries, then a legend in the 8th panel
png("olga_application/eigen_unrel_adjXkc_adjAutoPC15_xPrunedKC_ev12_eachCol.png",width=960)
par(mfrow=c(2,4))
plot(res$vectors[scan$race.cat=="Mexican",1],res$vectors[scan$race.cat=="Mexican",2],pch=20,col="goldenrod",
     xlab="EV1",ylab="EV2",main="Mexican",cex.lab=1.5,cex.axis=1.5)
plot(res$vectors[scan$race.cat=="CentralAmerican",1],res$vectors[scan$race.cat=="CentralAmerican",2],pch=20,col="red",xlab="EV1",ylab="EV2",main="Central Am")
plot(res$vectors[scan$race.cat=="SouthAmerican",1],res$vectors[scan$race.cat=="SouthAmerican",2],pch=20,col="magenta",xlab="EV1",ylab="EV2",main="South Am")
plot(res$vectors[scan$race.cat=="PuertoRican",1],res$vectors[scan$race.cat=="PuertoRican",2],pch=20,col="green",xlab="EV1",ylab="EV2",main="Puerto Rican")
plot(res$vectors[scan$race.cat=="Cuban",1],res$vectors[scan$race.cat=="Cuban",2],pch=20,col="blue",xlab="EV1",ylab="EV2",main="Cuban")
plot(res$vectors[scan$race.cat=="Dominican",1],res$vectors[scan$race.cat=="Dominican",2],pch=20,col="burlywood",xlab="EV1",ylab="EV2",main="Dominican")
plot(res$vectors[is.element(scan$race.cat,c("Unknown","Other")),1],res$vectors[is.element(scan$race.cat,c("Unknown","Other")),2],pch=20,col="black",xlab="EV1",ylab="EV2",main="Other/Unk")
plot(1:10,1:10,axes=FALSE,xlab="",ylab="",type="n")
legend("topleft",c("CentralAmerican","Cuban","Dominican","Mexican","PuertoRican","SouthAmerican","Other/Unknown"),pch=20,
       col=c("red","blue","burlywood","goldenrod","green","magenta","black"),cex=1.5)
dev.off()


pdf("olga_application/eigen_unrel_adjXkc_adjAutoPC15_xPrunedKC_ev23_col.pdf")
plotPC(2,3,xlab="EV2 (0.75%)",ylab="EV3 (0.31%)")
dev.off()

png("olga_application/eigen_unrel_adjXkc_adjAutoPC15_xPrunedKC_ev34_col.png")
plotPC(3,4,ylab="EV4 (0.27%)",xlab="EV3 (0.31%)")
dev.off()


# compare to x chr PC 1,2 estimates
pcX <- get(load("olga_application/pca_prunedXsnps_10272unrel.RData"))

xevs <- pcX$vectors[is.element(rownames(pcX$vectors),rownames(unrelkc)),]
dim(xevs) # 10509 10

cor(xevs[,1],res$vectors[,1]) # 0.9454644
cor(xevs[,2],res$vectors[,2]) # 0.8680538
cor(xevs[,3],res$vectors[,3]) # 0.6919725

png("olga_application/eigen_unrel_adjXkc_adjAutoPC15_xPrunedKC_ev1vsPCA_prunedX.png")
plot(xevs[,1],res$vectors[,1],xlab="X EV1",ylab="X auto adj EV1",pch=19,cex.axis=1.5,cex.lab=1.5,type="n")
col2=1; col1=1
points(xevs[scan$race.cat=="Mexican",col1],res$vectors[scan$race.cat=="Mexican",col2],pch=20,col="goldenrod")
points(xevs[scan$race.cat=="CentralAmerican",col1],res$vectors[scan$race.cat=="CentralAmerican",col2],pch=20,col="red")
points(xevs[scan$race.cat=="SouthAmerican",col1],res$vectors[scan$race.cat=="SouthAmerican",col2],pch=20,col="magenta")
points(xevs[scan$race.cat=="PuertoRican",col1],res$vectors[scan$race.cat=="PuertoRican",col2],pch=20,col="green")
points(xevs[scan$race.cat=="Cuban",col1],res$vectors[scan$race.cat=="Cuban",col2],pch=20,col="blue")
points(xevs[scan$race.cat=="Dominican",col1],res$vectors[scan$race.cat=="Dominican",col2],pch=20,col="burlywood")
points(xevs[scan$race.cat=="Other",col1],res$vectors[scan$race.cat=="Other",col2],pch=20,col="black")
points(xevs[scan$race.cat=="Unknown",col1],res$vectors[scan$race.cat=="Unknown",col2],pch=20,col="black")

abline(0,1,col="gray",lwd=2)
legend("topleft",expression(paste(rho,"=0.945")),bty="n",cex=1.3)
legend("bottomright",c("CentralAmerican","Cuban","Dominican","Mexican","PuertoRican","SouthAmerican","Other/Unknown"),pch=20,
       col=c("red","blue","burlywood","goldenrod","green","magenta","black"),cex=1.3)
dev.off()



# compare to auto PC 1,2 estimates
pcAuto <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/amstilp/results/pca/v08_study_unrelated_pcair/pca.v08_study_unrel.RData"))
dim(pcAuto$eigenvect) # 10642    32

autoevs <- pcAuto$eigenvect[is.element(pcAuto$sample.id,rownames(unrelkc)),]
dim(autoevs) # 9417 32

resx <- res$vectors[is.element(rownames(unrelkc),pcAuto$sample.id),]
dim(resx) # 9417 10509

# plot olga ev1 auto vs xchr ev1
# calculate corr

cor(autoevs[,1],-resx[,1]) # 0.8518446 | x chr and autos were 0.8704225

scan <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_scanAnnot_v54_AMS.RData"))
ids <- pcAuto$sample.id[is.element(pcAuto$sample.id,rownames(unrelkc))]
scan <- scan[is.element(scan$scanID,ids),]

png("olga_application/eigen_unrel_adjXkc_adjAutoPC15_xPrunedKC_ev1vsPCA_auto.png")
plot(autoevs[,1],-resx[,1],xlab="Autosomal EV1",ylab="- X auto adj EV1",type="n",cex.lab=1.5,cex.axis=1.5)
col2=1; col1=1
points(autoevs[scan$race.cat=="Mexican",col1],-resx[scan$race.cat=="Mexican",col2],pch=20,col="goldenrod")
points(autoevs[scan$race.cat=="CentralAmerican",col1],-resx[scan$race.cat=="CentralAmerican",col2],pch=20,col="red")
points(autoevs[scan$race.cat=="SouthAmerican",col1],-resx[scan$race.cat=="SouthAmerican",col2],pch=20,col="magenta")
points(autoevs[scan$race.cat=="PuertoRican",col1],-resx[scan$race.cat=="PuertoRican",col2],pch=20,col="green")
points(autoevs[scan$race.cat=="Cuban",col1],-resx[scan$race.cat=="Cuban",col2],pch=20,col="blue")
points(autoevs[scan$race.cat=="Dominican",col1],-resx[scan$race.cat=="Dominican",col2],pch=20,col="burlywood")
points(autoevs[scan$race.cat=="Other",col1],-resx[scan$race.cat=="Other",col2],pch=20,col="black")
points(autoevs[scan$race.cat=="Unknown",col1],-resx[scan$race.cat=="Unknown",col2],pch=20,col="black")

abline(0,1,col="gray",lwd=2)
legend("topleft",expression(paste(rho,"=0.852")),bty="n",cex=1.3)
legend("bottomright",c("CentralAmerican","Cuban","Dominican","Mexican","PuertoRican","SouthAmerican","Other/Unknown"),pch=20,
       col=c("red","blue","burlywood","goldenrod","green","magenta","black"),cex=0.8)
dev.off()

rm(list=ls())


#####
# 29. Do eigen() on KC_X matrix adj for X chr PC 1-2, KC_auto matrix adj for auto PC 1-5, KC_chr20 matrix adj for auto PC 1-5

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x/")
library(GWASTools)
library(SNPRelate)
source("pcrelate_X.R")

autoPC <- get(load("olga_application/pca_prunedAutosnps_10272unrelPlusRel.RData"))

pcMat <- autoPC$eigenvect[,c(1:5)]
dim(pcMat) # 12747 5                                                                                                                                                                                                             
rownames(pcMat) <- autoPC$sampleInfo[,"sample.id"]

scan <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_scanAnnot_v54_AMS.RData"))
head(pData(scan))

# make genoData object                                                                                                                                                                                                           
gdsobj <- GdsGenotypeReader("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/netCDF/subjects/OLGA_subj_geno.gds")
genoData <- GenotypeData(gdsobj,scanAnnot=scan[scan$subj.plink,])

unrel.set <- scan$scanID[scan$unrelated.pcair.deg4&is.na(scan$gengrp6.outliers)]
scanIncl <- scan$scanID[is.na(scan$gengrp6.outliers)&scan$geno.cntl==0&scan$subj.plink]

head(unrel.set); head(scanIncl)
length(unrel.set); length(scanIncl) # 10272 | 12747                                                                                                                                                                              

# read in pruned set of auto SNPs                                                                                                                                                                                                
auto.pruned <- get(load("olga_application/snp_sel_auto_ldPruned.RData"))
length(auto.pruned) # 152391                                                                                                                                                                                                     

###                                                                                                                                                                                                                              
# pc relate call on auto SNPs, adj for auto PC 1-5, no corr                                                                                                                                                                      
rel <- pcrelate(genoData,unrel.set=unrel.set,snp.include=auto.pruned,scan.include=scanIncl,Xchr=FALSE,pcMat=pcMat,correct=FALSE)
names(rel)

save(rel,file="olga_application/pcRelate_Auto_autoPC15adj_noCorrect.RData")

## make into matrix                                                                                                                                                                                                              
tmp <- matrix(NA,nrow=length(scanIncl),ncol=length(scanIncl))
kc <- rel$kinship
tmp[lower.tri(tmp,diag=FALSE)] <- kc[,"kin"]

# mirror to upper triangle                                                                                                                                                                                                       
tmp[upper.tri(tmp,diag=FALSE)] <- t(tmp)[upper.tri(tmp,diag=FALSE)]

rownames(tmp) <- scanIncl
colnames(tmp) <- scanIncl

inbd <- rel$inbreed
diag(tmp) <- 1+inbd$f

unrelPart <- pcairPartition(tmp)
unrelPart.set <- unrelPart$unrels
length(unrelPart.set) # 6974

unrelkcT <- tmp[is.element(rownames(tmp),unrelPart.set),]
unrelkc <- unrelkcT[,is.element(colnames(tmp),unrelPart.set)]
dim(unrelkc) # 6974 6974

eigenD <- eigen(unrelkc)
names(eigenD) # values | vectors
dim(eigenD$vectors) # 6974 6974

save(eigenD,file="olga_application/eigen_unrel_adjAutoKC_adjAutoPC15__noCorrect.RData")

#####                                                                                                                                                                                                                            

# get pruned set of chr 20 SNPs                                                                                                                                                                                                  
chr <- getChromosome(genoData)
snpID <- getSnpID(genoData)

idx <- which(is.element(snpID,auto.pruned))
chr.pruned <- chr[idx]
chr20.pruned <- auto.pruned[chr.pruned==20]
length(chr20.pruned) # 4438                                                                                                                                                                                                      

###                                                                                                                                                                                                                              
# pc relate call on chr 20 SNPs, adj for auto PC 1-5, no corr                                                                                                                                                                    
rel2 <- pcrelate(genoData,unrel.set=unrel.set,snp.include=chr20.pruned,scan.include=scanIncl,Xchr=FALSE,pcMat=pcMat,correct=FALSE)
names(rel2)

save(rel2,file="olga_application/pcRelate_chr20_autoPC15adj_noCorrect.RData")

## make into matrix                                                                                                                                                                                                              
tmp <- matrix(NA,nrow=length(scanIncl),ncol=length(scanIncl))
kc <- rel2$kinship
tmp[lower.tri(tmp,diag=FALSE)] <- kc[,"kin"]

# mirror to upper triangle                                                                                                                                                                                                       
tmp[upper.tri(tmp,diag=FALSE)] <- t(tmp)[upper.tri(tmp,diag=FALSE)]

rownames(tmp) <- scanIncl
colnames(tmp) <- scanIncl

inbd <-rel2$inbreed
diag(tmp) <- 1+inbd$f

unrelPart <- pcairPartition(tmp,kin.thresh=0.1)
unrelPart.set <- unrelPart$unrels
length(unrelPart.set) # 9889

unrelkcT <- tmp[is.element(rownames(tmp),unrelPart.set),]
unrelkc<- unrelkcT[,is.element(colnames(tmp),unrelPart.set)]
dim(unrelkc) # 9889 9889

eigenD <- eigen(unrelkc)
names(eigenD)
dim(eigenD$vectors) # 9889 9889

save(eigenD,file="olga_application/eigen_unrel_adjchr20KC_adjAutoPC15_noCorrect.RData")

#####                                                                                                                                                                                                                            

# read in pruned set of x chr SNPs                                                                                                                                                                                               
snp.pruned <- get(load("olga_application/snp_sel_xChr_ldPruned.RData"))
head(snp.pruned); length(snp.pruned) # 3600                                                                                                                                                                                      

# read in X chr PC mat                                                                                                                                                                                                           
xPC <- get(load("olga_application/pca_prunedXsnps_10272unrel.RData"))
pcMat <- xPC$vectors[,c(1:2)]
dim(pcMat) # 12747 2                                                                                                                                                                                                             

### need to exclude 13 people with an entirely missing x chr                                                                                                                                                                     
# from chromosome.anomalies.csv file, they are study samples with x chr anomaly                                                                                                                                                  
# thus, their entire x chr is filtered out.                                                                                                                                                                                      

chromAnoms <- read.csv("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/dbGaP/To_dbGaP/chromosome_anomalies/chromosome_anomalies.csv",
                       header=TRUE,as.is=T)
ids <- chromAnoms$scanID[chromAnoms$chromosome=="X"&chromAnoms$filter&chromAnoms$whole.chr&
                           !is.element(substr(chromAnoms$SUBJECT_ID,1,2),"NA")]
length(ids) # 13                                                                                                                                                                                                                 
chromAnoms[is.element(chromAnoms$scanID,ids),] # perfect! exactly what we want                                                                                                                                                   

length(scanIncl) # 12747                                                                                                                                                                                                         
scanIncl <- scanIncl[!is.element(scanIncl,ids)]
length(scanIncl) # 12734                                                                                                                                                                                                         

# take these individs out of pcMat too                                                                                                                                                                                           
pcMat <- pcMat[is.element(rownames(pcMat),scanIncl),]
dim(pcMat) # 12734 2                                                                                                                                                                                                             

rel <- pcrelate(genoData, unrel.set=unrel.set,snp.include=snp.pruned, scan.include=scanIncl, Xchr=TRUE, pcMat=pcMat,correct=FALSE)
names(rel)

save(rel,file="olga_application/pcRelate_prunedX_xPC12adj_noCorrect.RData")

## make into matrix                                                                                                                                                                                                              
tmp <- matrix(NA,nrow=length(scanIncl),ncol=length(scanIncl))
kc <- rel$kinship
tmp[lower.tri(tmp,diag=FALSE)] <- kc[,"kin"]

# mirror to upper triangle                                                                                                                                                                                                       
tmp[upper.tri(tmp,diag=FALSE)] <- t(tmp)[upper.tri(tmp,diag=FALSE)]

rownames(tmp) <- scanIncl
colnames(tmp) <- scanIncl

inbd <-rel$inbreed
diag(tmp) <- 1+inbd$f

unrelPart <- pcairPartition(tmp,kin.thresh=0.2)
unrelPart.set <- unrelPart$unrels
length(unrelPart.set) # 10894

unrelkcT <- tmp[is.element(rownames(tmp),unrelPart.set),]
unrelkc<- unrelkcT[,is.element(colnames(tmp),unrelPart.set)]
dim(unrelkc) # 10894

eigenD <- eigen(unrelkc)
names(eigenD)
dim(eigenD$vectors) # 10894 10894

save(eigenD,file="olga_application/eigen_unrel_adjxPrunedKC_adjxPC12_noCorrect.RData")

rm(list=ls())


#####
# 30. Parse chr 20 KC results

library(GWASTools)
setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")

scan <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_scanAnnot_v54_AMS.RData"))
head(pData(scan))

kc20 <- get(load("olga_application/pcRelate_chr20_autoPC15adj_noCorrect.RData"))
kc20 <- kc20$kinship

summary(kc20$kin) # ranges from -0.08 to 0.63

# make hist of x chr kc estimates for pairs that are unrel on auto
unrel.set <- scan$scanID[scan$unrelated.pcair.deg4&is.na(scan$gengrp6.outliers)]

toHist <- kc20[is.element(kc20$ID1,unrel.set)&is.element(kc20$ID2,unrel.set),]
dim(toHist) # 52751856 4
length(unique(toHist$ID1)) # 10271; great! length of unrel.set

pdf("olga_application/hist_chr20KC_autoPC15adj_autoUnrel.pdf")
hist(toHist$kin,breaks=50,main="Chr 20 KC for 12,734 Study Samples\nShowing Autosomal Unrelated Pairs",
     xlab="Chr 20 KC Adj For Auto PC 1-5",cex.lab=1.5,cex.main=1.5)
legend("topright",c(paste("mean =",format(mean(kc20$kin),digits=4)),paste("sd =",format(sd(kc20$kin),digits=4))),
       bty="n",cex=1.5)
dev.off()

pdf("olga_application/hist_chr20KC_autoPC15adj_autoUnrel_trunc.pdf")
hist(toHist$kin,main="Chr 20 KC for auto unrel samples",breaks=50,
     xlab="Chr 20 KC Adj For Auto PC 1-5",ylim=c(0,5000))
dev.off()

rm(list=ls())


#####
# 31. Parse results from #29.

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")
library(GWASTools)

## results from auto pc corr for auto pc
res <- get(load("olga_application/eigen_unrel_adjAutoKC_adjAutoPC15_noCorrect.RData"))
res$values[1:5]/sum(res$values) # 0.0067024027 0.0020627673 0.0005386106 0.0005005866

scan <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_scanAnnot_v54_AMS.RData"))
scan <- scan[is.element(scan$scanID,rownames(res$vectors)),]

png("olga_application/eigen_unrel_adjAutokc_adjAutoPC15_ev12.png")
plot(res$vectors[,1],res$vectors[,2],xlab="EV1 (0.67%)",ylab="EV2 (0.21%)")
dev.off()

plotPC <- function(col1,col2,xlab,ylab){
  plot(res$vectors[,col1],res$vectors[,col2],xlab=xlab,ylab=ylab,type="n",cex.lab=1.5,cex.axis=1.5)
  points(res$vectors[scan$race.cat=="Mexican",col1],res$vectors[scan$race.cat=="Mexican",col2],
         pch=20,col="goldenrod")
  points(res$vectors[scan$race.cat=="CentralAmerican",col1],res$vectors[scan$race.cat=="CentralAmerican",col2],pch=20,col="red")
  points(res$vectors[scan$race.cat=="SouthAmerican",col1],res$vectors[scan$race.cat=="SouthAmerican",col2],pch=20,col="magenta")
  points(res$vectors[scan$race.cat=="PuertoRican",col1],res$vectors[scan$race.cat=="PuertoRican",col2],pch=20,col="green")
  points(res$vectors[scan$race.cat=="Cuban",col1],res$vectors[scan$race.cat=="Cuban",col2],pch=20,col="blue")
  points(res$vectors[scan$race.cat=="Dominican",col1],res$vectors[scan$race.cat=="Dominican",col2],pch=20,col="burlywood")
  points(res$vectors[scan$race.cat=="Other",col1],res$vectors[scan$race.cat=="Other",col2],pch=20,col="black")
  points(res$vectors[scan$race.cat=="Unknown",col1],res$vectors[scan$race.cat=="Unknown",col2],pch=20,col="black")
  legend("bottomleft",c("CentralAmerican","Cuban","Dominican","Mexican","PuertoRican","SouthAmerican","Other/Unknown"),pch=20,
         col=c("red","blue","burlywood","goldenrod","green","magenta","black"))
}

png("olga_application/eigen_unrel_adjAutokc_adjAutoPC15_ev12_col.png")
plotPC(1,2,xlab="EV1 (0.67%)",ylab="EV2 (0.21%)")
dev.off()

png("olga_application/eigen_unrel_adjAutokc_adjAutoPC15_pairs14.png")
pairs(res$vectors[,1:4],res$vectors[,1:4],pch=19)#,text.panel=text(0.5, 0.5, "EV1"))
dev.off()

png("olga_application/eigen_unrel_adjAutokc_adjAutoPC15_ev34_col.png")
plotPC(3,4,ylab="EV4 (0.054%)",xlab="EV3 (0.050%)")
dev.off()


#### eigen() on chr 20 snps, adj for auto pcs
res <- get(load("olga_application/eigen_unrel_adjchr20KC_adjAutoPC15_noCorrect.RData"))
res$values[1:5]/sum(res$values) # 0.013322172 0.005416008 0.002540405 0.002452159 0.001849004

scan <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_scanAnnot_v54_AMS.RData"))
scan <- scan[is.element(scan$scanID,rownames(res$vectors)),]

png("olga_application/eigen_unrel_adjchr20KC_adjAutoPC15_ev12_col.png")
plotPC(1,2,xlab="EV1 (1.33%)",ylab="EV2 (0.54%)")
dev.off()

png("olga_application/eigen_unrel_adjchr20KC_adjAutoPC15_ev12.png")
plot(res$vectors[,1],res$vectors[,2],xlab="EV1 (1.33%)",ylab="EV2 (0.54%)")
dev.off()

png("olga_application/eigen_unrel_adjchr20KC_adjAutoPC15_ev34_col.png")
plotPC(3,4,xlab="EV3 (0.254%)",ylab="EV4 (0.245%)")
dev.off()



#### eigen() on x chr snps, adj for x chr structure
res <- get(load("olga_application/eigen_unrel_adjxPrunedKC_adjxPC12_noCorrect.RData"))
res$values[1:5]/sum(res$values) # 0.003292980 0.002903684 0.002430698 0.002341678 0.002253568

scan <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_scanAnnot_v54_AMS.RData"))
scan <- scan[is.element(scan$scanID,rownames(res$vectors)),]

png("olga_application/eigen_unrel_adjxPrunedKC_adjxPC12_ev12_col.png")
plotPC(1,2,xlab="EV1 (0.33%)",ylab="EV2 (0.29%)")
dev.off()
# great! this is what we are hoping for

png("olga_application/eigen_unrel_adjxPrunedKC_adjxPC12_pairs14.png")
pairs(res$vectors[,1:4],res$vectors[,1:4],pch=19)
dev.off()

png("olga_application/eigen_unrel_adjxPrunedKC_adjxPC12_ev34_col.png")
plotPC(3,4,xlab="EV3 (0.24%)",ylab="EV4 (0.23%)")
dev.off()


#### scree plot of all results on one axis
reschr20 <- get(load("olga_application/eigen_unrel_adjchr20KC_adjAutoPC15_noCorrect.RData"))
reschrXauto <- get(load("olga_application/eigen_unrel_adjXkc_adjAutoPC15_xPrunedKC_noCorrect_thresh2.RData"))
reschrXX <- get(load("olga_application/eigen_unrel_adjxPrunedKC_adjxPC12_noCorrect.RData"))
resauto <- get(load("olga_application/eigen_unrel_adjAutoKC_adjAutoPC15_noCorrect.RData"))

pdf("olga_application/eigen_unrel_screes_ev115.pdf")
plot(x=1:15,y=seq(from=0,to=0.02,length=15),xlab="EV",ylab="Proportion of Variance Explained",type="n")
points(1:15,reschr20$values[1:15]/sum(reschr20$values),type="b",col="purple",pch=1,lwd=1.5)
points(1:15,reschrXauto$values[1:15]/sum(reschrXauto$values),type="b",col="cyan",pch=2,lwd=1.5)
points(1:15,reschrXX$values[1:15]/sum(reschrXX$values),type="b",col="goldenrod",pch=3,lwd=1.5)
points(1:15,resauto$values[1:15]/sum(resauto$values),type="b",col="magenta",pch=5,lwd=1.5)
legend("topright",c("Chr 20 adj for Auto","Chr X adj for Auto","Chr X adj for X","Autos adj for Auto"),
       col=c("purple","cyan","goldenrod","magenta"),lty=1,pch=c(1,2,3,5),lwd=1.5)
dev.off()

rm(list=ls())


#####
# 32. Rerun X chr PCA

# evs on all subjects, excl 19 asian outliers
# use the pruned set of snps calculated earlier

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x/")
source("pcair_X_scan.include.R")
library(QCpipeline)
library(OLGApipeline)
library(OLGAanalysis)
library(GWASTools)

scan <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_scanAnnot_v54_AMS.RData"))
scanFreeze <- getobj("/projects/geneva/gcc-fs2/OLGA/genotype/freeze1/gds/freeze1/SOL_freeze1_scanAnnot.RData")

scanBoth <- merge(pData(scan),pData(scanFreeze),all.y=TRUE,by="scanID")
dim(scanBoth); head(scanBoth) # 12803 97

snp.pruned <- get(load("olga_application/snp_sel_xChr_ldPruned.RData"))
head(snp.pruned); length(snp.pruned) # 3600

# make genoData object
gdsobj <- GdsGenotypeReader("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/netCDF/subjects/OLGA_subj_geno.gds")
genoData <- GenotypeData(gdsobj,scanAnnot=scan[scan$subj.plink,])

unrel.set <- scanBoth$scanID[scanBoth$unrelated.pcair.deg4&!is.element(scanBoth$gengrp6.outliers,"AsianOutliers")]
scanIncl <- scanBoth$scanID[!is.element(scanBoth$gengrp6.outliers,"AsianOutliers")&scanBoth$subj.plink]

length(unrel.set); length(scanIncl) # 10287 | 12784
head(unrel.set); head(scanIncl) 

all(is.element(scanIncl,scan$scanID)) # TRUE

pc <- pcair(genoData,unrel.set=unrel.set,Xchr=TRUE,snp.include=snp.pruned,scan.include=scanIncl)
#Unrelated Set: 10287 Samples 
#Related Set: 2497 Samples
#Running Analysis with 3600 SNPs ...

save(pc,file="olga_application/pca_prunedXsnps_10287unrel12784tot.RData")

scan <- scan[is.element(scan$scanID,scanIncl),]
pc$values/pc$sum.values # 0.034951833 0.018769995 0.005358512 0.004926456

pdf("olga_application/pca_X_12784tot_ev12_col.pdf")
plot(-pc$vectors[,1],-pc$vectors[,2],xlab="EV1 (3.50%)",ylab="EV2 (1.88%)",type="n")
points(-pc$vectors[scan$race.cat=="Mexican",1],-pc$vectors[scan$race.cat=="Mexican",2],pch=20,col="goldenrod")
points(-pc$vectors[scan$race.cat=="CentralAmerican",1],-pc$vectors[scan$race.cat=="CentralAmerican",2],pch=20,col="red")
points(-pc$vectors[scan$race.cat=="SouthAmerican",1],-pc$vectors[scan$race.cat=="SouthAmerican",2],pch=20,col="magenta")
points(-pc$vectors[scan$race.cat=="PuertoRican",1],-pc$vectors[scan$race.cat=="PuertoRican",2],pch=20,col="green")
points(-pc$vectors[scan$race.cat=="Cuban",1],-pc$vectors[scan$race.cat=="Cuban",2],pch=20,col="blue")
points(-pc$vectors[scan$race.cat=="Dominican",1],-pc$vectors[scan$race.cat=="Dominican",2],pch=20,col="burlywood")
points(-pc$vectors[scan$race.cat=="Other",1],-pc$vectors[scan$race.cat=="Other",2],pch=20,col="black")
points(-pc$vectors[scan$race.cat=="Unknown",1],-pc$vectors[scan$race.cat=="Unknown",2],pch=20,col="black")
legend("bottomleft",c("CentralAmerican","Cuban","Dominican","Mexican","PuertoRican","SouthAmerican","Other/Unknown"),pch=20,
       col=c("red","blue","burlywood","goldenrod","green","magenta","black"))
dev.off()

# do ev1 vs ev2 for all 7 entries, then a legend in the 8th panel
pdf("olga_application/pca_X_12784tot_ev12_eachCol.pdf",width=14)
par(mfrow=c(2,4))
plot(-pc$vectors[scan$race.cat=="Mexican",1],-pc$vectors[scan$race.cat=="Mexican",2],pch=20,col="goldenrod",xlab="EV1",ylab="EV2",main="Mexican")
plot(-pc$vectors[scan$race.cat=="CentralAmerican",1],-pc$vectors[scan$race.cat=="CentralAmerican",2],pch=20,col="red",xlab="EV1",ylab="EV2",main="Central Am")
plot(-pc$vectors[scan$race.cat=="SouthAmerican",1],-pc$vectors[scan$race.cat=="SouthAmerican",2],pch=20,col="magenta",xlab="EV1",ylab="EV2",main="South Am")
plot(-pc$vectors[scan$race.cat=="PuertoRican",1],-pc$vectors[scan$race.cat=="PuertoRican",2],pch=20,col="green",xlab="EV1",ylab="EV2",main="Puerto Rican")
plot(-pc$vectors[scan$race.cat=="Cuban",1],-pc$vectors[scan$race.cat=="Cuban",2],pch=20,col="blue",xlab="EV1",ylab="EV2",main="Cuban")
plot(-pc$vectors[scan$race.cat=="Dominican",1],-pc$vectors[scan$race.cat=="Dominican",2],pch=20,col="burlywood",xlab="EV1",ylab="EV2",main="Dominican")
plot(-pc$vectors[is.element(scan$race.cat,c("Unknown","Other")),1],-pc$vectors[is.element(scan$race.cat,c("Unknown","Other")),2],pch=20,col="black",xlab="EV1",ylab="EV2",main="Other/Unk")
plot(1:10,1:10,axes=FALSE,xlab="",ylab="",type="n")
legend("topleft",c("CentralAmerican","Cuban","Dominican","Mexican","PuertoRican","SouthAmerican","Other/Unknown"),pch=20,
       col=c("red","blue","burlywood","goldenrod","green","magenta","black"),cex=1.5)
dev.off()

pdf("olga_application/pca_X_12784tot_ev34_col.pdf")
plot(-pc$vectors[,3],pc$vectors[,4],xlab="EV3 (0.54%)",ylab="EV4 (0.49%)",type="n")
points(-pc$vectors[scan$race.cat=="Mexican",3],pc$vectors[scan$race.cat=="Mexican",4],pch=20,col="goldenrod")
points(-pc$vectors[scan$race.cat=="CentralAmerican",3],pc$vectors[scan$race.cat=="CentralAmerican",4],pch=20,col="red")
points(-pc$vectors[scan$race.cat=="SouthAmerican",3],pc$vectors[scan$race.cat=="SouthAmerican",4],pch=20,col="magenta")
points(-pc$vectors[scan$race.cat=="PuertoRican",3],pc$vectors[scan$race.cat=="PuertoRican",4],pch=20,col="green")
points(-pc$vectors[scan$race.cat=="Cuban",3],pc$vectors[scan$race.cat=="Cuban",4],pch=20,col="blue")
points(-pc$vectors[scan$race.cat=="Dominican",3],pc$vectors[scan$race.cat=="Dominican",4],pch=20,col="burlywood")
points(-pc$vectors[scan$race.cat=="Other",3],pc$vectors[scan$race.cat=="Other",4],pch=20,col="black")
points(-pc$vectors[scan$race.cat=="Unknown",3],pc$vectors[scan$race.cat=="Unknown",4],pch=20,col="black")
#legend("topleft",c("CentralAmerican","Cuban","Dominican","Mexican","PuertoRican","SouthAmerican","Other/Unknown"),pch=20,
#       col=c("red","blue","burlywood","goldenrod","green","magenta","black"))
dev.off()

rm(list=ls())


#####
# 33. Rerun pcrelate on xchr, adj for x chr PCs from 32.

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")
library(GWASTools)
source("pcrelate_X.R")

# use unrel.set to be the same unrelated set used for PCA on the X chr
# use scan.include to be the same set of individs used for PCA on the X chr
scan <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_scanAnnot_v54_AMS.RData"))
scanFreeze <- getobj("/projects/geneva/gcc-fs2/OLGA/genotype/freeze1/gds/freeze1/SOL_freeze1_scanAnnot.RData")
scan <- merge(pData(scan),pData(scanFreeze),by="scanID",all.y=TRUE)

unrel.set <- scan$scanID[scan$unrelated.pcair.deg4&!is.element(scan$gengrp6.outliers,"AsianOutliers")]
scanIncl <- scan$scanID[!is.element(scan$gengrp6.outliers,"AsianOutliers")&scan$subj.plink]

head(unrel.set); head(scanIncl)
length(unrel.set); length(scanIncl) # 10287 | 12784

# read in pruned set of x chr SNPs                                                                                                                                                                                               
snp.pruned <- get(load("olga_application/snp_sel_xChr_ldPruned.RData"))
head(snp.pruned); length(snp.pruned) # 3600                                                                                                                                                                                      

# read in X chr PC mat                                                                                                                                                                                                           
xPC <- get(load("olga_application/pca_prunedXsnps_10287unrel12784tot.RData"))
pcMat <- xPC$vectors[,c(1:2)]
dim(pcMat) # 12784 2

### need to exclude 13 people with an entirely missing x chr                                                                                                                                                                     
# from chromosome.anomalies.csv file, they are study samples with x chr anomaly                                                                                                                                                  
# thus, their entire x chr is filtered out.                                                                                                                                                                                      
chromAnoms <- read.csv("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/dbGaP/To_dbGaP/chromosome_anomalies/chromosome_anomalies.csv",
                       header=TRUE,as.is=T)
ids <- chromAnoms$scanID[chromAnoms$chromosome=="X"&chromAnoms$filter&chromAnoms$whole.chr&
                           !is.element(substr(chromAnoms$SUBJECT_ID,1,2),"NA")]
length(ids) # 13                                                                                                                                                                                                                 
chromAnoms[is.element(chromAnoms$scanID,ids),] # perfect! exactly what we want                                                                                                                                                   

length(scanIncl) # 12784
scanIncl <- scanIncl[!is.element(scanIncl,ids)]
length(scanIncl) # 12771                                                                                                                                                                                                         

# take these individs out of pcMat too                                                                                                                                                                                           
pcMat <- pcMat[is.element(rownames(pcMat),scanIncl),]
dim(pcMat) # 12771 2 

# make genoData object                                                                                                                                                                                                           
gdsobj <- GdsGenotypeReader("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/netCDF/subjects/OLGA_subj_geno.gds")
scan <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_scanAnnot_v54_AMS.RData"))
genoData <- GenotypeData(gdsobj,scanAnnot=scan[scan$subj.plink,])

rel <- pcrelate(genoData, unrel.set=unrel.set,snp.include=snp.pruned, scan.include=scanIncl, Xchr=TRUE, pcMat=pcMat,correct=FALSE)
names(rel)

save(rel,file="olga_application/pcRelate_prunedX_12771tot_xPC12adj_noCorrect.RData")

rm(list=ls())


#####
# 34. Est VC for LABA2 BCC trait, analysis 316987

# see script in
# /projects/geneva/geneva_sata/caitlin/olga_xchr_assoc/run_assoc_estVC.R
# qsub -q thornton.q -N ana316987 batch_run_assoc_estVC.sh
# writes run_assoc_estVC_316987.Rout
## this includes adj for: auto PC 1-5; auto KC, block group and HH as random effects

# qsub -q thornton.q -N ana316987 batch_run_assoc_estVC.sh
# writes run_assoc_estVC_316987.Rout
## this includes adj for: auto PC 1-5, x chr PC 1-2; auto KC, xchr KC, block group and HH as random effects

# qsub -q thornton.q -N ana316987 batch_run_assoc_estVC.sh
# writes run_assoc_estVC_noAutoPC_316987.Rout
## this includes adj for: x chr PC 1-2; auto KC, xchr KC, block group and HH as random effects

# qsub -q thornton.q -N ana316987 batch_run_assoc_estVC.sh
# writes run_assoc_estVC_noAutoPCnoAutoKC_316987.Rout
## this includes adj for: x chr PC 1-2; xchr KC, block group and HH as random effects

# qsub -q thornton.q -N ana316987 batch_run_assoc_estVC.sh
# writes run_assoc_estVC_noAutoKC_316987.Rout
## this includes adj for: auto PC 1-5, x chr PC 1-2; xchr KC, block group and HH as random effects


#####
# 35. Look at var comp estimates for the 4 models, analysis 316987

setwd("/projects/geneva/geneva_sata/caitlin/olga_xchr_assoc")
library(GWASTools)

allIncl <- getobj("Assoc/varCompFrac_316987.RData")
noAutoPC <- getobj("Assoc/varCompFrac_noAutoPC_316987.RData")
rownames(allIncl) <- rownames(noAutoPC) <- c("block","hh","auto","xchr","envr")
allIncl$model <- "auto+xchr"
noAutoPC$model <- "xchr"
allIncl$est <- rownames(allIncl)
noAutoPC$est <- rownames(noAutoPC)

noPCKC <- getobj("Assoc/varCompFrac_noAutoPCnoAutoKC_316987.RData")
noKC <- getobj("Assoc/varCompFrac_noAutoKC_316987.RData")
rownames(noPCKC) <- rownames(noKC) <- c("block","hh","xchr","envr")
noPCKC$model <- "xchr"
noKC$model <- "auto+xchr"
noPCKC$est <- rownames(noPCKC)
noKC$est <- rownames(noKC)

auto <- getobj("Assoc/varCompFrac_auto_316987.RData")
rownames(auto) <- c("block","hh","auto","envr")
auto$model <- "autosomal"
auto$est <- rownames(auto)

library(ggplot2)
pdf("varComp_bothKC_compPC.pdf")
dfc <- rbind(allIncl,noAutoPC)
colnames(dfc) <- c("Proportion","Lower","Upper","model","est")
ggplot(dfc, aes(x=model, y=Proportion, colour=est, group=model)) + 
  geom_errorbar(aes(ymin=Lower, ymax=Upper), width=.1) +
#  geom_line(position=pd) +
  geom_point(size=3) + # 21 is filled circle
  xlab("Var Component") +
  ylab("Proportion") +
  ggtitle("Variance Component Estimates Incl X and Auto KC\nWith and Without Auto PC 1-5 Covar Adjustment") +
#  scale_y_continuous(limits=c(0, max(dfc$len + dfc$se)),    # Set y range
#                     breaks=0:20*4) +                       # Set tick every 4
  theme_bw() +
  theme(legend.justification=c(1,0), legend.position=c(1,0)) # Position legend in bottom right
dev.off()


pdf("varComp_xKC_compPC.pdf")
dfc <- rbind(noPCKC,noKC)
colnames(dfc) <- c("Proportion","Lower","Upper","model","est")
ggplot(dfc, aes(x=model, y=Proportion, colour=est, group=model)) + 
  geom_errorbar(aes(ymin=Lower, ymax=Upper), width=.1) +
  #  geom_line(position=pd) +
  geom_point(size=3) + # 21 is filled circle
  xlab("Var Component") +
  ylab("Proportion") +
  ggtitle("Variance Component Estimates Only Incl X KC\nWith and Without Auto PC 1-5 Covar Adjustment") +
  #  scale_y_continuous(limits=c(0, max(dfc$len + dfc$se)),    # Set y range
  #                     breaks=0:20*4) +                       # Set tick every 4
  theme_bw() +
  theme(legend.justification=c(1,0), legend.position=c(1,0)) # Position legend in bottom right
dev.off()


pdf("varComp_compPC.pdf")
noPCKC$model <- "xPC+xKC"
noKC$model <- "bothPC+xKC"
noAutoPC$model <- "xPC+bothKC"
allIncl$model <- "bothPC+bothKC"
auto$model <- "autosomal"
dfc <- rbind(noPCKC,noKC,allIncl,noAutoPC,auto)
colnames(dfc) <- c("Proportion","Lower","Upper","model","VarComp")
ggplot(dfc, aes(x=model, y=Proportion, colour=VarComp, group=model)) + 
  geom_errorbar(aes(ymin=Lower, ymax=Upper), width=.1) +
  #  geom_line(position=pd) +
  geom_point(size=3) + # 21 is filled circle
  xlab("Model") +
  ylab("Proportion") +
  ggtitle("Variance Component Proportion Estimates\nWith 95% CIs") +
  #  scale_y_continuous(limits=c(0, max(dfc$len + dfc$se)),    # Set y range
  #                     breaks=0:20*4) +                       # Set tick every 4
  theme_bw() +
  theme(legend.justification=c(1,0), legend.position=c(1,0.5)) # Position legend in bottom right
dev.off()

# plot just the xchr ones so we can see better
pdf("varComp_compPC_xchrOnly.pdf")
ggplot(dfc[dfc$VarComp=="xchr",], aes(x=model, y=Proportion, colour=VarComp, group=model)) + 
  geom_errorbar(aes(ymin=Lower, ymax=Upper), width=.1) +
  #  geom_line(position=pd) +
  geom_point(size=3) + # 21 is filled circle
  xlab("Model") +
  ylab("Proportion") + 
  ggtitle("Variance Component Proportion Estimates With 95% CIs\nX Chromosome Estimate Only") +
    scale_y_continuous(limits=c(0, 0.05),    # Set y range
                       breaks=seq(from=0,to=0.05,by=0.01)) +                       # Set tick every 4
  theme_bw() +
  theme(legend.justification=c(1,0), legend.position=c(1,0)) # Position legend in bottom right
dev.off()

library(reshape); library(xtable)
dat <- melt(dfc,id=c("model","VarComp"))
xt <- cast(dat[dat$variable=="Proportion",],model~VarComp)
xtable(xt,digits=5)

rm(list=ls())


#####
# 36. Assoc test for LABA2 BCC trait, analysis 316987

# cd /projects/geneva/geneva_sata/caitlin/olga_xchr_assoc
# ran run_assoc_mixed_segment_auto_316987.R: adj for auto pc 1-5, auto kc
# ran run_assoc_mixed_segment.R: adj for x pc 1-2, auto pc 1-5, auto kc, xchr kc
# ran run_assoc_mixed_segment_noAutoPC_316987.R: adj for x pc 1-2, auto kc, xchr kc
# ran run_assoc_mixed_segment_noAutoPCnoAutoKC_316987.R: adj for x pc 1-2, xchr kc
# ran run_assoc_mixed_segment_noAutoKC_316987.R: adj for x pc 1-2, auto pc 1-5, xchr kc
# and combine_assoc.R for all of those, too

# stored: cd Assoc/
# assoc_auto_316987_chr23.RData
# assoc_316987_chr23.RData
# assoc_noAutoPC_316987_chr23.RData
# assoc_noAutoPCnoAutoKC_316987_chr23.RData
# assoc_noAutoKC_316987_chr23.RData


#####
# 37. Parse assoc test results for X chr

# filter to effN>30
# do manh/qq plot on these SNPs

# cd /projects/geneva/geneva_sata/caitlin/olga_xchr_assoc
# ran plot_qq_manh.R for each of the assoc test results

setwd("/projects/geneva/geneva_sata/caitlin/olga_xchr_assoc/")
library(GWASTools)
library(OLGApipeline)

auto <- getobj("Assoc/assoc_auto_316987_chr23.RData")
allIncl <- getobj("Assoc/assoc_316987_chr23.RData")
nopc <- getobj("Assoc/assoc_noAutoPC_316987_chr23.RData")
nopckc <- getobj("Assoc/assoc_noAutoPCnoAutoKC_316987_chr23.RData")
nokc <- getobj("Assoc/assoc_noAutoKC_316987_chr23.RData")

dim(auto); dim(allIncl); dim(nopc)
dim(nopckc); dim(nokc) # all 936453 17

allequal(auto$effN,allIncl$effN) # TRUE
allequal(auto$snpID,nokc$snpID) # TRUE

filt <- !is.na(auto$pval)&auto$effN>30
table(filt) # 669963 true

png("Plots/auto_vs__pvals_316987.png")
plot(-log10(auto$pval[filt&auto$pval<1e-4]),-log10(allIncl$pval[filt&auto$pval<1e-4]),,xlab="-log10(Pvals) Autosomal Model",
     ylab="-log10(Pvals) Incl X PC1-2, X KC",cex.lab=1.5,cex.axis=1.5,col=rgb(1,0,0,0.5),pch=19)
abline(a=0,b=1)
dev.off()

# do pairs plot of all pvalues, filtered at pval<1e-4
filt2 <- auto$pval<1e-4 | allIncl$pval<1e-4 | nopc$pval<1e-4 | nopckc$pval<1e-4 | nokc$pval<1e-4
table(filt2) # 737 true
table(filt&filt2) # exactly 500 true

autos=-log10(auto$pval[filt&filt2])
AdjForAll=-log10(allIncl$pval[filt&filt2])
NoAutoPC=-log10(nopc$pval[filt&filt2])
NoAutoPCorKC=-log10(nopckc$pval[filt&filt2])
NoAutoKC=-log10(nokc$pval[filt&filt2])
allvals <- data.frame(autos=-log10(auto$pval[filt&filt2]),
                      AdjForAll=-log10(allIncl$pval[filt&filt2]),
                      NoAutoPC=-log10(nopc$pval[filt&filt2]),
                      NoAutoPCorKC=-log10(nopckc$pval[filt&filt2]),
                      NoAutoKC=-log10(nokc$pval[filt&filt2]))
dim(allvals) # 500 5

library(GGally)

png("Plots/pairs_316987.png",height=720,width=720)
custom_plot <- ggpairs(allvals,diag=list(continuous="bar"),
        params=c(col=rgb(1,0,0,0.5)),lower="blank",title="PValue Comparison\nCutoff at -log10(pval)>4, effN>30")
plot21 <- ggplot2::ggplot(allvals,ggplot2::aes(x=autos,y=AdjForAll)) + 
  ggplot2::geom_point(color="red",alpha=0.5) + ggplot2::geom_abline(intercept=0,slope=1)
custom_plot <- putPlot(custom_plot,plot21,2,1)
plot31 <-  ggplot2::ggplot(allvals,ggplot2::aes(x=autos,y=NoAutoPC)) + 
  ggplot2::geom_point(color="red",alpha=0.5) + ggplot2::geom_abline(intercept=0,slope=1)
custom_plot <- putPlot(custom_plot,plot31,3,1)
plot32 <-  ggplot2::ggplot(allvals,ggplot2::aes(x=AdjForAll,y=NoAutoPC)) + 
  ggplot2::geom_point(color="red",alpha=0.5) + ggplot2::geom_abline(intercept=0,slope=1)
custom_plot <- putPlot(custom_plot,plot32,3,2)
plot41 <-  ggplot2::ggplot(allvals,ggplot2::aes(x=autos,y=NoAutoPCorKC)) + 
  ggplot2::geom_point(color="red",alpha=0.5) + ggplot2::geom_abline(intercept=0,slope=1)
custom_plot <- putPlot(custom_plot,plot41,4,1)
plot42 <-  ggplot2::ggplot(allvals,ggplot2::aes(x=AdjForAll,y=NoAutoPCorKC)) + 
  ggplot2::geom_point(color="red",alpha=0.5) + ggplot2::geom_abline(intercept=0,slope=1)
custom_plot <- putPlot(custom_plot,plot42,4,2)
plot43 <-  ggplot2::ggplot(allvals,ggplot2::aes(x=NoAutoPC,y=NoAutoPCorKC)) + 
  ggplot2::geom_point(color="red",alpha=0.5) + ggplot2::geom_abline(intercept=0,slope=1)
custom_plot <- putPlot(custom_plot,plot43,4,3)
plot51 <-  ggplot2::ggplot(allvals,ggplot2::aes(x=autos,y=NoAutoKC)) + 
  ggplot2::geom_point(color="red",alpha=0.5) + ggplot2::geom_abline(intercept=0,slope=1)
custom_plot <- putPlot(custom_plot,plot51,5,1)
plot52 <-  ggplot2::ggplot(allvals,ggplot2::aes(x=AdjForAll,y=NoAutoKC)) + 
  ggplot2::geom_point(color="red",alpha=0.5) + ggplot2::geom_abline(intercept=0,slope=1)
custom_plot <- putPlot(custom_plot,plot52,5,2)
plot53 <-  ggplot2::ggplot(allvals,ggplot2::aes(x=NoAutoPC,y=NoAutoKC)) + 
  ggplot2::geom_point(color="red",alpha=0.5) + ggplot2::geom_abline(intercept=0,slope=1)
custom_plot <- putPlot(custom_plot,plot53,5,3)
plot54 <-  ggplot2::ggplot(allvals,ggplot2::aes(x=NoAutoPCorKC,y=NoAutoKC)) + 
  ggplot2::geom_point(color="red",alpha=0.5) + ggplot2::geom_abline(intercept=0,slope=1)
custom_plot <- putPlot(custom_plot,plot54,5,4)

custom_plot
dev.off()

filt2 <- (!is.na(auto$pval)&auto$pval<1e-3&auto$pval>1e-4) | (!is.na(allIncl$pval)&allIncl$pval<1e-3&allIncl$pval>1e-4) | 
  (!is.na(nopc$pval)&nopc$pval<1e-3&nopc$pval>1e-4) |
  (!is.na(nopckc$pval)&nopckc$pval<1e-3&nopckc$pval>1e-4) | (!is.na(nokc$pval)&nokc$pval<1e-3&nokc$pval>1e-4)
table(filt2) # 2475 true
table(filt&filt2) # 1548 true

allvals <- data.frame(autos=-log10(auto$pval[filt&filt2]),
                      AdjForAll=-log10(allIncl$pval[filt&filt2]),
                      NoAutoPC=-log10(nopc$pval[filt&filt2]),
                      NoAutoPCorKC=-log10(nopckc$pval[filt&filt2]),
                      NoAutoKC=-log10(nokc$pval[filt&filt2]))
dim(allvals) # 1548 5

png("Plots/pairs_log103_4_316987.png",width=720,height=720)
custom_plot <- ggpairs(allvals,diag=list(continuous="bar"),
                       params=c(col=rgb(1,0,0,0.5)),lower="blank",title="PValue Comparison\nCutoff at -log10(pval)>3&<4, effN>30")
plot21 <- ggplot2::ggplot(allvals,ggplot2::aes(x=autos,y=AdjForAll)) + 
  ggplot2::geom_point(color="red",alpha=0.5) + ggplot2::geom_abline(intercept=0,slope=1)
custom_plot <- putPlot(custom_plot,plot21,2,1)
plot31 <-  ggplot2::ggplot(allvals,ggplot2::aes(x=autos,y=NoAutoPC)) + 
  ggplot2::geom_point(color="red",alpha=0.5) + ggplot2::geom_abline(intercept=0,slope=1)
custom_plot <- putPlot(custom_plot,plot31,3,1)
plot32 <-  ggplot2::ggplot(allvals,ggplot2::aes(x=AdjForAll,y=NoAutoPC)) + 
  ggplot2::geom_point(color="red",alpha=0.5) + ggplot2::geom_abline(intercept=0,slope=1)
custom_plot <- putPlot(custom_plot,plot32,3,2)
plot41 <-  ggplot2::ggplot(allvals,ggplot2::aes(x=autos,y=NoAutoPCorKC)) + 
  ggplot2::geom_point(color="red",alpha=0.5) + ggplot2::geom_abline(intercept=0,slope=1)
custom_plot <- putPlot(custom_plot,plot41,4,1)
plot42 <-  ggplot2::ggplot(allvals,ggplot2::aes(x=AdjForAll,y=NoAutoPCorKC)) + 
  ggplot2::geom_point(color="red",alpha=0.5) + ggplot2::geom_abline(intercept=0,slope=1)
custom_plot <- putPlot(custom_plot,plot42,4,2)
plot43 <-  ggplot2::ggplot(allvals,ggplot2::aes(x=NoAutoPC,y=NoAutoPCorKC)) + 
  ggplot2::geom_point(color="red",alpha=0.5) + ggplot2::geom_abline(intercept=0,slope=1)
custom_plot <- putPlot(custom_plot,plot43,4,3)
plot51 <-  ggplot2::ggplot(allvals,ggplot2::aes(x=autos,y=NoAutoKC)) + 
  ggplot2::geom_point(color="red",alpha=0.5) + ggplot2::geom_abline(intercept=0,slope=1)
custom_plot <- putPlot(custom_plot,plot51,5,1)
plot52 <-  ggplot2::ggplot(allvals,ggplot2::aes(x=AdjForAll,y=NoAutoKC)) + 
  ggplot2::geom_point(color="red",alpha=0.5) + ggplot2::geom_abline(intercept=0,slope=1)
custom_plot <- putPlot(custom_plot,plot52,5,2)
plot53 <-  ggplot2::ggplot(allvals,ggplot2::aes(x=NoAutoPC,y=NoAutoKC)) + 
  ggplot2::geom_point(color="red",alpha=0.5) + ggplot2::geom_abline(intercept=0,slope=1)
custom_plot <- putPlot(custom_plot,plot53,5,3)
plot54 <-  ggplot2::ggplot(allvals,ggplot2::aes(x=NoAutoPCorKC,y=NoAutoKC)) + 
  ggplot2::geom_point(color="red",alpha=0.5) + ggplot2::geom_abline(intercept=0,slope=1)
custom_plot <- putPlot(custom_plot,plot54,5,4)

custom_plot
dev.off()

# do a 5 panel manh plot
plotTitles <- c("LABA2 ~ covars + EV1 + EV2 + EV3 + EV4 + EV5 + (1|auto KC)",
                "LABA2 ~ covars + EV1 + EV2 + EV3 + EV4 + EV5 + XEV1 + XEV2 + (1|auto KC) + (1|X KC)",
                "LABA2 ~ covars + XEV1 + XEV2 + (1|auto KC) + (1|X KC)",
                "LABA2 ~ covars + XEV1 + XEV2 + (1|X KC)",
                "LABA2 ~ covars + EV1 + EV2 + EV3 + EV4 + EV5 + XEV1 + XEV2 + (1|X KC)")
pvalues <- data.frame(auto=auto$pval[filt],allIncl=allIncl$pval[filt],nopc=nopc$pval[filt],
                      nopckc=nopckc$pval[filt],nokc=nokc$pval[filt])


png("Plots/all_manhUnfilt_316987.png",width=720,height=1200)
main="LABA2 X Chromosome Manhattan Plots, SNPs with MAC>30"
oma <- c(0, 0, length(main) + 0.1, 0)
addText=""
par(mfrow = c(5, 1), mar = c(5, 5, 4, 2) + 
      0.1, lwd = 1.5, cex.lab = 1.5, cex.main = 1.5, oma = oma)
for (i in 1:5) {
  title <- plotTitles[i]
  manhattanPlot(p = pvalues[,i], chromosome = rep(23,nrow(pvalues)), 
                main = title)
  if (i == 1) {
    mtext(side = 3, line = -1, text = addText, padj = 0.9, 
          adj = 0.02, outer = T, cex = 1.5)
  }
  title(main, outer = T)
  
}
dev.off()

## look at 'index SNP' across all 5 models
which.min(auto$pval)
minP <- rbind(auto[which.min(auto$pval),],allIncl[which.min(allIncl$pval),],nokc[which.min(nokc$pval),],
              nopc[which.min(nopc$pval),],nopckc[which.min(nopckc$pval),])
rownames(minP) <- c("auto","allIncl","noKC","noPC","noPCKC")
library(xtable)
xtable(cbind(minP[,c("Beta","SE","Stat")],format(minP$pval,digits=6,scientific=TRUE)),digits=4)

minP
# snpID chromosome     n        MAF minor.allele      Beta         SE
# auto    27883110          X 12489 0.01979151            B 0.1314149 0.01451786
# allIncl 27883110          X 12489 0.01979151            B 0.1300359 0.01483121
# noKC    27883110          X 12489 0.01979151            B 0.1318145 0.01488613
# noPC    27883110          X 12489 0.01979151            B 0.1300585 0.01480597
# noPCKC  27883110          X 12489 0.01979151            B 0.1314804 0.01485779
# Stat         pval     effN  position alleleA alleleB   oevar
# auto    81.93767 1.404471e-19 386.8517 153764217       C       T 1.03777
# allIncl 76.87290 1.823214e-18 386.8517 153764217       C       T 1.03777
# noKC    78.40841 8.379459e-19 386.8517 153764217       C       T 1.03777
# noPC    77.16209 1.574873e-18 386.8517 153764217       C       T 1.03777
# noPCKC  78.30933 8.810494e-19 386.8517 153764217       C       T 1.03777
# segment info type
# auto         28    1    2
# allIncl      28    1    2
# noKC         28    1    2
# noPC         28    1    2
# noPCKC       28    1    2

# read in snp annot to get rsID and verify position is hg19
olgaData <- OlgaGenotypeData("/projects/geneva/gcc-fs2/OLGA/genotype/freeze1/gds/freeze1")
snpAnnot <- getSnpAnnotation(olgaData, chromosome=23)
pData(snpAnnot)[snpAnnot$snpID==minP$snpID[1],]
# snpID         snp        rsID  position alleleA alleleB chromosome
# 931902 27883110 kgp30626516 kgp30626516 153764217       C       T         23
# segment exp_freq_a1 info certainty type info_type0 concord_type0
# 931902      28        0.02    1         1    2      0.978         0.996
# r2_type0   oevar oevar.male oevar.female
# 931902    0.899 1.03777          1        1.064

snp <- getobj("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_HCHS_Custom_15041502_B3_all37_v25_AMS.RData")
varMetadata(snp) # yep, is hg19 position 

rm(list=ls())


#####
# 38. Run assoc test on autosomes now

# cd /projects/geneva/geneva_sata/caitlin/olga_xchr_assoc

library(GWASTools)
library(OLGApipeline)
olgaData <- OlgaGenotypeData("/projects/geneva/gcc-fs2/OLGA/genotype/freeze1/gds/freeze1")
snpAnnot <- getSnpAnnotation(olgaData, chromosome=chromosome)
range(snpAnnot$segment)
# do this for values of chromosome=1,...,22
rm(list=ls())

# figure out how many segments exist for each chromosome
#chr1: 45
#chr2: 47
#chr3: 37
#chr4: 36
#chr5: 35
#chr6: 32
#chr7: 30
#chr8: 27
#chr9: 22
#chr10: 25
#chr11: 26
#chr12: 24
#chr13: 19
#chr14: 17
#chr15: 16
#chr16: 15 -- combine with 17
#chr17: 15
#chr18: 14
#chr19: 11 -- combine with 20
#chr20: 11
#chr21: 7 -- combine w 22
#chr22: 7

# qsub -q thornton.q -N chr1 -t 1-45 batch_run_assoc_mixed_segment.sh
# qsub -q thornton.q -N chr2 -t 1-47 batch_run_assoc_mixed_segment2.sh
# qsub -q thornton.q -N chr3 -t 1-37 batch_run_assoc_mixed_segment3.sh
# qsub -q thornton.q -N chr4 -t 1-36 batch_run_assoc_mixed_segment4.sh
# qsub -q thornton.q -N chr5 -t 1-35 batch_run_assoc_mixed_segment5.sh
# qsub -q thornton.q -N chr6 -t 1-32 batch_run_assoc_mixed_segment6.sh
# qsub -q thornton.q -N chr7 -t 1-30 batch_run_assoc_mixed_segment7.sh
# qsub -q thornton.q -N chr8 -t 1-27 batch_run_assoc_mixed_segment8.sh
# qsub -q thornton.q -N chr9 -t 1-22 batch_run_assoc_mixed_segment9.sh
# qsub -q thornton.q -N chr10 -t 1-25 batch_run_assoc_mixed_segment10.sh
# qsub -q thornton.q -N chr11 -t 1-26 batch_run_assoc_mixed_segment11.sh
# qsub -q thornton.q -N chr12 -t 1-24 batch_run_assoc_mixed_segment12.sh
# qsub -q thornton.q -N chr13 -t 1-19 batch_run_assoc_mixed_segment13.sh
# qsub -q thornton.q -N chr14 -t 1-17 batch_run_assoc_mixed_segment14.sh
# qsub -q thornton.q -N chr15 -t 1-16 batch_run_assoc_mixed_segment15.sh
# qsub -q thornton.q -N chr1617 -t 1-15 batch_run_assoc_mixed_segment1617.sh
# qsub -q thornton.q -N chr18 -t 1-14 batch_run_assoc_mixed_segment18.sh
# qsub -q thornton.q -N chr1920 -t 1-11 batch_run_assoc_mixed_segment1920.sh
# qsub -q thornton.q -N chr2122 -t 1-7 batch_run_assoc_mixed_segment2122.sh

# ran run_assoc_mixed_segment_auto_316987.R: adj for auto pc 1-5, auto kc

# ran run_assoc_mixed_segment.R: adj for x pc 1-2, auto pc 1-5, auto kc, xchr kc
# ran run_assoc_mixed_segment_noAutoPC_316987.R: adj for x pc 1-2, auto kc, xchr kc
# ran run_assoc_mixed_segment_noAutoPCnoAutoKC_316987.R: adj for x pc 1-2, xchr kc
# ran run_assoc_mixed_segment_noAutoKC_316987.R: adj for x pc 1-2, auto pc 1-5, xchr kc
# and combine_assoc.R for all of those, too

# stored: cd Assoc/
# assoc_auto_316987_chrXX.RData for XX=1,2,...,22,23
# assoc_316987_chrXX.RData for XX=1,2,...,22,23
# assoc_noAutoPC_316987_chrXX.RData for XX=1,2,...,22,23
# assoc_noAutoPCnoAutoKC_316987_chrXX.RData for XX=1,2,...,22,23
# assoc_noAutoKC_316987_chrXX.RData for XX=1,2,...,22,23

# cd /projects/geneva/geneva_sata/caitlin/olga_xchr_assoc
# ran plot_qq_manh.R for each of the assoc test results
chromSep <- "_chr"
files <- list.files(directory, pattern = "assoc.+RData")
base <- unique(matrix(unlist((strsplit(files, chromSep))), 
                      ncol = 2, byrow = TRUE)[, 1])

#####
# 39. Summary statistics for trait LABA2 RBC

library(GWASTools)
library(OLGApipeline)

base.path <- "/projects/geneva/gcc-fs2/OLGA"
base.path.local <- "/projects/geneva/geneva_sata/caitlin/olga_xchr_assoc"
 
 ## connect to DB
analysis_id=316987
db <- getDb("olga_analysis")
attr <- dbGetAnalysisAttributes(db,analysis_id)
config <- dbGetConfig(db, analysis_id,pathPrefix=base.path)
dbDisconnect(db)

scanAnnot <- getobj(file.path(base.path, attr$analysis_path, "Input", paste0("scanAnnot_", analysis_id, ".RData")))

library(ggplot2)
pdf("laba2_hist.pdf")
ggplot(pData(scanAnnot[scanAnnot$assoc.include,]), aes(x=LABA2)) + 
  geom_histogram() +
  #  geom_line(position=pd) +
   xlab("LABA2") +
  ggtitle("Histogram of LABA2 for 12,502 Samples Used in Association Test") +
  #  scale_y_continuous(limits=c(0, max(dfc$len + dfc$se)),    # Set y range
  #                     breaks=0:20*4) +                       # Set tick every 4
  theme_bw() 
dev.off()

# regress laba2 on 5 auto pcs, look at r2
l <- lm(scanAnnot$LABA2[scanAnnot$assoc.include]~scanAnnot$EV1[scanAnnot$assoc.include]+scanAnnot$EV2[scanAnnot$assoc.include]+
          scanAnnot$EV3[scanAnnot$assoc.include]+scanAnnot$EV4[scanAnnot$assoc.include]+scanAnnot$EV5[scanAnnot$assoc.include])
summary(l)
# Multiple R-squared:  0.00365,  Adjusted R-squared:  0.003251; so pretty small!

# merge in x chr pc 1,2
xpc <- get(load("/projects/geneva/geneva_sata/caitlin/mlm_x/olga_application/pca_prunedXsnps_10287unrel12784tot.RData"))
xvec <- data.frame(xpc$vectors[,1:2])
xvec$scanID <- rownames(xpc$vectors)
scan <- merge(pData(scanAnnot),xvec,by="scanID",all.x=TRUE)

l <- lm(scan$LABA2[scan$assoc.include]~scan$X1[scan$assoc.include]+scan$X2[scan$assoc.include])
summary(l)
# Multiple R-squared:  0.00101,  Adjusted R-squared:  0.0008503; quite small for the x chr effect

rm(list=ls())


#####
# 40. Calculate X chr EVs by gengrp6.strat

# evs on each gengrp6.strat set
# use the pruned set of snps calculated earlier

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x/")
source("pcair_X_scan.include.R")
library(QCpipeline)
library(OLGApipeline)
library(OLGAanalysis)
library(GWASTools)

scan <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_scanAnnot_v54_AMS.RData"))
scanFreeze <- getobj("/projects/geneva/gcc-fs2/OLGA/genotype/freeze1/gds/freeze1/SOL_freeze1_scanAnnot.RData")

scanBoth <- merge(pData(scan),pData(scanFreeze),all.y=TRUE,by="scanID")
dim(scanBoth); head(scanBoth) # 12803 97

snp.pruned <- get(load("olga_application/snp_sel_xChr_ldPruned.RData"))
head(snp.pruned); length(snp.pruned) # 3600

# make genoData object
gdsobj <- GdsGenotypeReader("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/netCDF/subjects/OLGA_subj_geno.gds")
genoData <- GenotypeData(gdsobj,scanAnnot=scan[scan$subj.plink,])


### CAm
unrel.set <- scanBoth$scanID[scanBoth$unrelated.pcair.deg4&!is.element(scanBoth$gengrp6.outliers,"AsianOutliers")&
                               is.element(scanBoth$gengrp6.strat,"CentralAmerican")]
scanIncl <- scanBoth$scanID[!is.element(scanBoth$gengrp6.outliers,"AsianOutliers")&scanBoth$subj.plink&
                              is.element(scanBoth$gengrp6.strat,"CentralAmerican")]

length(unrel.set); length(scanIncl) # 1159 | 1373
head(unrel.set); head(scanIncl) 

all(is.element(scanIncl,scan$scanID)) # TRUE

pc <- pcair(genoData,unrel.set=unrel.set,Xchr=TRUE,snp.include=snp.pruned,scan.include=scanIncl)
#Unrelated Set: 1159 Samples 
#Related Set: 214 Samples
#Running Analysis with 3600 SNPs ...

pc$values/pc$sum.values # 0.019437057 0.009851550 0.009140147 0.008454571

save(pc,file="olga_application/pca_prunedXsnps_10287unrel1373CAm.RData")


### Cuban
unrel.set <- scanBoth$scanID[scanBoth$unrelated.pcair.deg4&!is.element(scanBoth$gengrp6.outliers,"AsianOutliers")&
                               is.element(scanBoth$gengrp6.strat,"Cuban")]
scanIncl <- scanBoth$scanID[!is.element(scanBoth$gengrp6.outliers,"AsianOutliers")&scanBoth$subj.plink&
                              is.element(scanBoth$gengrp6.strat,"Cuban")]

length(unrel.set); length(scanIncl) # 1928 | 2238
head(unrel.set); head(scanIncl) 

all(is.element(scanIncl,scan$scanID)) # TRUE

pc <- pcair(genoData,unrel.set=unrel.set,Xchr=TRUE,snp.include=snp.pruned,scan.include=scanIncl)
#Unrelated Set: 1928 Samples 
#Related Set: 310 Samples
#Running Analysis with 3600 SNPs ...

pc$values/pc$sum.values # 0.038423784 0.005891484 0.005596158 0.004937166

save(pc,file="olga_application/pca_prunedXsnps_10287unrel2238Cuban.RData")


### Dominican
unrel.set <- scanBoth$scanID[scanBoth$unrelated.pcair.deg4&!is.element(scanBoth$gengrp6.outliers,"AsianOutliers")&
                               is.element(scanBoth$gengrp6.strat,"Dominican")]
scanIncl <- scanBoth$scanID[!is.element(scanBoth$gengrp6.outliers,"AsianOutliers")&scanBoth$subj.plink&
                              is.element(scanBoth$gengrp6.strat,"Dominican")]

length(unrel.set); length(scanIncl) # 902 | 1182
head(unrel.set); head(scanIncl) 

all(is.element(scanIncl,scan$scanID)) # TRUE

pc <- pcair(genoData,unrel.set=unrel.set,Xchr=TRUE,snp.include=snp.pruned,scan.include=scanIncl)
#Unrelated Set: 902 Samples 
#Related Set: 280 Samples
#Running Analysis with 3600 SNPs ...

pc$values/pc$sum.values # 0.026481741 0.010949196 0.009961251 0.007189963

save(pc,file="olga_application/pca_prunedXsnps_10287unrel1182Domin.RData")


### Mexican
unrel.set <- scanBoth$scanID[scanBoth$unrelated.pcair.deg4&!is.element(scanBoth$gengrp6.outliers,"AsianOutliers")&
                               is.element(scanBoth$gengrp6.strat,"Mexican")]
scanIncl <- scanBoth$scanID[!is.element(scanBoth$gengrp6.outliers,"AsianOutliers")&scanBoth$subj.plink&
                              is.element(scanBoth$gengrp6.strat,"Mexican")]

length(unrel.set); length(scanIncl) # 3725 | 4744
head(unrel.set); head(scanIncl) 

all(is.element(scanIncl,scan$scanID)) # TRUE

pc <- pcair(genoData,unrel.set=unrel.set,Xchr=TRUE,snp.include=snp.pruned,scan.include=scanIncl)
#Unrelated Set: 3725 Samples 
#Related Set: 1019 Samples
#Running Analysis with 3600 SNPs ...

pc$values/pc$sum.values # 0.022370852 0.008395218 0.007664308 0.006012739

save(pc,file="olga_application/pca_prunedXsnps_10287unrel4744Mex.RData")


### PR
unrel.set <- scanBoth$scanID[scanBoth$unrelated.pcair.deg4&!is.element(scanBoth$gengrp6.outliers,"AsianOutliers")&
                               is.element(scanBoth$gengrp6.strat,"PuertoRican")]
scanIncl <- scanBoth$scanID[!is.element(scanBoth$gengrp6.outliers,"AsianOutliers")&scanBoth$subj.plink&
                              is.element(scanBoth$gengrp6.strat,"PuertoRican")]

length(unrel.set); length(scanIncl) # 1780 | 2238
head(unrel.set); head(scanIncl) 

all(is.element(scanIncl,scan$scanID)) # TRUE

pc <- pcair(genoData,unrel.set=unrel.set,Xchr=TRUE,snp.include=snp.pruned,scan.include=scanIncl)
#Unrelated Set: 1780 Samples 
#Related Set: 458 Samples
#Running Analysis with 3600 SNPs ...

pc$values/pc$sum.values # 0.021980928 0.008487731 0.007881785 0.006924683

save(pc,file="olga_application/pca_prunedXsnps_10287unrel2238PR.RData")


### SAm
unrel.set <- scanBoth$scanID[scanBoth$unrelated.pcair.deg4&!is.element(scanBoth$gengrp6.outliers,"AsianOutliers")&
                               is.element(scanBoth$gengrp6.strat,"SouthAmerican")]
scanIncl <- scanBoth$scanID[!is.element(scanBoth$gengrp6.outliers,"AsianOutliers")&scanBoth$subj.plink&
                              is.element(scanBoth$gengrp6.strat,"SouthAmerican")]

length(unrel.set); length(scanIncl) # 758 | 908
head(unrel.set); head(scanIncl) 

all(is.element(scanIncl,scan$scanID)) # TRUE

pc <- pcair(genoData,unrel.set=unrel.set,Xchr=TRUE,snp.include=snp.pruned,scan.include=scanIncl)
#Unrelated Set: 758 Samples 
#Related Set: 150 Samples
#Running Analysis with 3600 SNPs ...

pc$values/pc$sum.values # 0.037533511 0.009575404 0.008235663 0.007638856

save(pc,file="olga_application/pca_prunedXsnps_10287unrel908SAm.RData")

rm(list=ls())


#####
# 41. Parse X EVs by gengrp6.strat results

library(GWASTools)
setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")

scan <- getobj("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_scanAnnot_v54_AMS.RData")

# merge in all the pc values, for each of the gengrp6.strat groups

## CAm first; this is what the sig hits were found in
pc <- getobj("olga_application/pca_prunedXsnps_10287unrel1373CAm.RData")
dim(pc$vectors) # 1373 10

vecs <- pc$vectors
vecs <- data.frame(vecs)
vecs$scanID <- rownames(pc$vectors)

scan <- merge(pData(scan),vecs[,c("X1","X2","scanID")],by="scanID",all.x=TRUE)
table(is.na(scan$X1),scan$gengrp6.strat) # is false for the CentralAmerican, 1373. good!

pc$values[1:2]/pc$sum.values # 0.01943706 0.00985155

scanSm <- scan[is.element(scan$scanID,rownames(pc$vectors)),]

pdf("olga_application/pca_X_1373CAm_ev12_col.pdf")
plot(pc$vectors[,1],pc$vectors[,2],xlab="EV1 (1.94%)",ylab="EV2 (0.99%)",pch=20,col="red",
     main="1373 Central American Samples")
dev.off()


## Cuban
pc <- getobj("olga_application/pca_prunedXsnps_10287unrel2238Cuban.RData")
dim(pc$vectors) # 2238 10

vecs <- pc$vectors
vecs <- data.frame(vecs)
vecs$scanID <- rownames(pc$vectors)

scan <- merge(scan,vecs[,c("X1","X2","scanID")],by="scanID",all.x=TRUE)
scan$X1.x[!is.na(scan$X1.y)] <- scan$X1.y[!is.na(scan$X1.y)]
scan$X2.x[!is.na(scan$X2.y)] <- scan$X2.y[!is.na(scan$X2.y)]
table(is.na(scan$X1.x),is.na(scan$X1.y))
#      FALSE  TRUE
#FALSE  2238  1373
#TRUE      0 10549

colnames(scan) # last 2 cols are .y, .y, remove these
colnames(scan)[is.element(colnames(scan),"X1.x")] <- "X1"
colnames(scan)[is.element(colnames(scan),"X2.x")] <- "X2"

scan <- scan[,-c(98,99)]
table(is.na(scan$X1),scan$gengrp6.strat) # is false for the CentralAmerican, Cuban

pc$values[1:2]/pc$sum.values # 0.038423784 0.005891484

scanSm <- scan[is.element(scan$scanID,rownames(pc$vectors)),]

pdf("olga_application/pca_X_2238Cuban_ev12_col.pdf")
plot(pc$vectors[,1],pc$vectors[,2],xlab="EV1 (3.84%)",ylab="EV2 (0.59%)",pch=20,col="blue",
     main="2238 Cuban Samples")
dev.off()


### Dominican
pc <- getobj("olga_application/pca_prunedXsnps_10287unrel1182Domin.RData")
dim(pc$vectors) # 2238 10

vecs <- pc$vectors
vecs <- data.frame(vecs)
vecs$scanID <- rownames(pc$vectors)

scan <- merge(scan,vecs[,c("X1","X2","scanID")],by="scanID",all.x=TRUE)
scan$X1.x[!is.na(scan$X1.y)] <- scan$X1.y[!is.na(scan$X1.y)]
scan$X2.x[!is.na(scan$X2.y)] <- scan$X2.y[!is.na(scan$X2.y)]
table(is.na(scan$X1.x),is.na(scan$X1.y))
#      FALSE  TRUE
#FALSE  1183  3611
#TRUE      0  9367

colnames(scan) # last 2 cols are .y, .y, remove these
colnames(scan)[is.element(colnames(scan),"X1.x")] <- "X1"
colnames(scan)[is.element(colnames(scan),"X2.x")] <- "X2"

scan <- scan[,-c(98,99)]
table(is.na(scan$X1),scan$gengrp6.strat) # is false for the CentralAmerican, Cuban, Domin

pc$values[1:2]/pc$sum.values # 0.02648174 0.01094920

scanSm <- scan[is.element(scan$scanID,rownames(pc$vectors)),]

pdf("olga_application/pca_X_1183Domin_ev12_col.pdf")
plot(pc$vectors[,1],pc$vectors[,2],xlab="EV1 (2.65%)",ylab="EV2 (1.09%)",pch=20,col="burlywood",
     main="1183 Dominican Samples")
dev.off()


### Mexican
pc <- getobj("olga_application/pca_prunedXsnps_10287unrel4744Mex.RData")
dim(pc$vectors) # 4744 10

vecs <- pc$vectors
vecs <- data.frame(vecs)
vecs$scanID <- rownames(pc$vectors)

scan <- merge(scan,vecs[,c("X1","X2","scanID")],by="scanID",all.x=TRUE)
scan$X1.x[!is.na(scan$X1.y)] <- scan$X1.y[!is.na(scan$X1.y)]
scan$X2.x[!is.na(scan$X2.y)] <- scan$X2.y[!is.na(scan$X2.y)]
table(is.na(scan$X1.x),is.na(scan$X1.y))
#      FALSE  TRUE
#FALSE  4744  4793
#TRUE      0  4623

colnames(scan) # last 2 cols are .y, .y, remove these
colnames(scan)[is.element(colnames(scan),"X1.x")] <- "X1"
colnames(scan)[is.element(colnames(scan),"X2.x")] <- "X2"

scan <- scan[,-c(98,99)]
table(is.na(scan$X1),scan$gengrp6.strat) # is false for the CentralAmerican, Cuban, Domin, Mex

pc$values[1:2]/pc$sum.values # 0.022370852 0.008395218

scanSm <- scan[is.element(scan$scanID,rownames(pc$vectors)),]

pdf("olga_application/pca_X_4744Mex_ev12_col.pdf")
plot(pc$vectors[,1],pc$vectors[,2],xlab="EV1 (2.24%)",ylab="EV2 (0.84%)",pch=20,col="goldenrod",
     main="4744 Mexican Samples")
dev.off()


### Puerto Rican
pc <- getobj("olga_application/pca_prunedXsnps_10287unrel2238PR.RData")
dim(pc$vectors) # 2238 10

vecs <- pc$vectors
vecs <- data.frame(vecs)
vecs$scanID <- rownames(pc$vectors)

scan <- merge(scan,vecs[,c("X1","X2","scanID")],by="scanID",all.x=TRUE)
scan$X1.x[!is.na(scan$X1.y)] <- scan$X1.y[!is.na(scan$X1.y)]
scan$X2.x[!is.na(scan$X2.y)] <- scan$X2.y[!is.na(scan$X2.y)]
table(is.na(scan$X1.x),is.na(scan$X1.y))
#      FALSE  TRUE
#FALSE  2238  9537
#TRUE      0  2385

colnames(scan) # last 2 cols are .y, .y, remove these
colnames(scan)[is.element(colnames(scan),"X1.x")] <- "X1"
colnames(scan)[is.element(colnames(scan),"X2.x")] <- "X2"

scan <- scan[,-c(98,99)]
table(is.na(scan$X1),scan$gengrp6.strat) # is false for the CentralAmerican, Cuban, Domin, Mex, PRs

pc$values[1:2]/pc$sum.values # 0.021980928 0.008487731

scanSm <- scan[is.element(scan$scanID,rownames(pc$vectors)),]

pdf("olga_application/pca_X_2238PR_ev12_col.pdf")
plot(pc$vectors[,1],pc$vectors[,2],xlab="EV1 (2.20%)",ylab="EV2 (0.85%)",pch=20,col="green",
     main="2238 Puerto Rican Samples")
dev.off()


### SAm
pc <- getobj("olga_application/pca_prunedXsnps_10287unrel908SAm.RData")
dim(pc$vectors) # 908 10

vecs <- pc$vectors
vecs <- data.frame(vecs)
vecs$scanID <- rownames(pc$vectors)

scan <- merge(scan,vecs[,c("X1","X2","scanID")],by="scanID",all.x=TRUE)
scan$X1.x[!is.na(scan$X1.y)] <- scan$X1.y[!is.na(scan$X1.y)]
scan$X2.x[!is.na(scan$X2.y)] <- scan$X2.y[!is.na(scan$X2.y)]
table(is.na(scan$X1.x),is.na(scan$X1.y))
#      FALSE  TRUE
#FALSE  908  11775
#TRUE     0   1477

colnames(scan) # last 2 cols are .y, .y, remove these
colnames(scan)[is.element(colnames(scan),"X1.x")] <- "X1"
colnames(scan)[is.element(colnames(scan),"X2.x")] <- "X2"

scan <- scan[,-c(98,99)]
table(is.na(scan$X1),scan$gengrp6.strat) # is false for the CentralAmerican, Cuban, Domin, Mex, PRs, SAm

pc$values[1:2]/pc$sum.values # 0.037533511 0.009575404

scanSm <- scan[is.element(scan$scanID,rownames(pc$vectors)),]

pdf("olga_application/pca_X_908SAm_ev12_col.pdf")
plot(pc$vectors[,1],pc$vectors[,2],xlab="EV1 (3.75%)",ylab="EV2 (0.96%)",pch=20,col="magenta",
     main="908 South American Samples")
dev.off()

# save these X EVs
xpc <- scan[,c("X1","X2","scanID")]
colnames(xpc) <- c("XEV1","XEV2","scanID")
save(xpc,file="olga_application/pca_prunedXsnps_stratAnalys.RData")

rm(list=ls())


#####
# 42. Est VC for BMI trait, analysis 298888

# see script in
# /projects/geneva/geneva_sata/caitlin/olga_xchr_assoc/id298888_bmi/run_assoc_estVC.R
# qsub -q thornton.q -N estVC batch_run_assoc_estVC.sh
# writes run_assoc_estVC_auto_298888.Rout
## this includes adj for: auto PC 1-5; auto KC, block group and HH as random effects

# qsub -q thornton.q -N allIncl batch_run_assoc_estVC.sh
# writes run_assoc_estVC_298888.Rout
## this includes adj for: auto PC 1-5, x chr PC 1-2; auto KC, xchr KC, block group and HH as random effects

# qsub -q thornton.q -N estVC_noAutoKC batch_run_assoc_estVC.sh
# writes run_assoc_estVC_noAutoKC_298888.Rout
## this includes adj for: auto PC 1-5, x chr PC 1-2; xchr KC, block group and HH as random effects

# qsub -q thornton.q -N estVC_noAutoPC batch_run_assoc_estVC.sh
# writes run_assoc_estVC_noAutoPC_298888.Rout
## this includes adj for: x chr PC 1-2; auto KC, xchr KC, block group and HH as random effects

# qsub -q thornton.q -N estVC_noPCnoKC batch_run_assoc_estVC.sh
# writes run_assoc_estVC_noAutoPCnoAutoKC_298888.Rout
## this includes adj for: x chr PC 1-2; xchr KC, block group and HH as random effects

## hmm, x chr VC is 0 for all models. dang!


#####
# 43. Assoc test for BMI trait, analysis 298888

# cd /projects/geneva/geneva_sata/caitlin/olga_xchr_assoc/id298888_bmi/
# ran run_assoc_mixed_segment_noAutoPC_298888.R: adj for x pc 1-2, auto kc, xchr kc
# ran run_assoc_mixed_segment_noAutoKC_298888.R: adj for x pc 1-2, auto pc 1-5, xchr kc
# ran run_assoc_mixed_segment_auto_298888.R: adj for auto pc 1-5, auto kc
# ran run_assoc_mixed_segment.R: adj for x pc 1-2, auto pc 1-5, auto kc, xchr kc
# ran run_assoc_mixed_segment_noAutoPCnoAutoKC_298888.R: adj for x pc 1-2, xchr kc
# and combine_assoc.R for all of those, too

# stored: cd Assoc/
# assoc_auto_298888_chr23.RData
# assoc_298888_chr23.RData
# assoc_noAutoPC_298888_chr23.RData
# assoc_noAutoPCnoAutoKC_298888_chr23.RData
# assoc_noAutoKC_298888_chr23.RData


#####
# 44. Look at corr of X PC 1-2 with BMI

library(GWASTools)
library(OLGApipeline); library(QCpipeline)
library(OLGAanalysis)

analysis_id <- 298888
chromosome <- 23; segment <- 1
dbname <- "olga_analysis"

base.path <- "/projects/geneva/gcc-fs2/OLGA"
base.path.local <- "/projects/geneva/geneva_sata/caitlin/olga_xchr_assoc/id298888_bmi"

## connect to DB
db <- getDb(dbname)
attr <- dbGetAnalysisAttributes(db,analysis_id)
config <- dbGetConfig(db, analysis_id,pathPrefix=base.path)

### CHANGE HERE FOR NEW FIXED EFFECT COVARS ###
config$covars <- c(config$covars,"X1","X2")

dbDisconnect(db)

## check config and set defaults
required <- c("in_gds_file", "outcome", "covars", "analysis_type", "logistic")
optional <- c("block.size", "impute.geno", "genotypeInteractions")
default <- c(5000, FALSE, NA)
config <- setConfigDefaults(config, required, optional, default)
print(config)

scanAnnot <- getobj(file.path(base.path, attr$analysis_path, "Input", paste0("scanAnnot_", analysis_id, ".RData")))

xvec <- getobj("/projects/geneva/geneva_sata/caitlin/mlm_x/olga_application/pca_prunedXsnps_stratAnalys.RData")
colnames(xvec) <- c("X1","X2","scanID")
scan <- merge(pData(scanAnnot),xvec,by="scanID",all.x=TRUE)

cor(scan$log.BMI[scan$assoc.include],scan$X1[scan$assoc.include]) # 0.05442797
cor(scan$log.BMI[scan$assoc.include],scan$X2[scan$assoc.include]) # 0.01256959

rm(list=ls())


#####
# 45. Run ADMIXTURE on the X chromosome

# first run in batch prune and subset gds file to pruned x chr SNPs
# cd /projects/geneva/geneva_sata/caitlin/mlm_x
# qsub -q thornton.q -N xfiltgds batch_gds2ped_subset.sh

## ran admixture with three sets of snps:
# 1. no ld pruning, n=9694
# 2. ld pruning with usual thresh of 0.32, n=2233
# 3. ld pruning with less stringent thresh of 0.48, n=3494

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x/olga_application/admixture")
library(GWASTools)

res1 <- read.table("olga_ref_unrel_xpruned.3.Q",as.is=TRUE)
res2 <- read.table("olga_ref_unrel_x.3.Q",as.is=TRUE)
res3 <- read.table("olga_ref_unrel_xpruned_048.3.Q",as.is=TRUE)

pop <- read.table("olga_ref_unrel_xpruned.pop")
fam <- read.table("olga_ref_unrel_xpruned.fam")

colnames(res1) <- paste("est",1:3,"_pruned",sep="")
colnames(res2) <- paste("est",1:3,sep="")
colnames(res3) <- paste("est",1:3,"_pruned048",sep="")

res <- cbind(fam$V2,pop,res1,res2,res3)
colnames(res)[1:2] <- c("scanID","pop")

res[res$pop=="NAM",] # est3 is NAM
head(res[res$pop=="EUR",]) # est1 is EUR 
# est2 is AFR

colnames(res)[3:5] <- paste(c("EUR","AFR","NAM"),"_pruned",sep="")
colnames(res)[6:8] <- paste(c("EUR","AFR","NAM"),"_all",sep="")
colnames(res)[9:11] <- paste(c("EUR","AFR","NAM"),"_pruned048",sep="")

apply(res[,3:5],2,summary)
#               EUR     AFR     NAM
#Min.       0.00001    0.00001    0.00001
#1st Qu.    0.20900    0.01871    0.13030
#Median     0.40200    0.08579    0.35390
#Mean       0.42490    0.17730    0.39770
#3rd Qu.    0.62840    0.23630    0.63800
#Max.       1.00000    1.00000    1.00000

# write a table to use in svn at tower
towrite <- res[,c(1,2,3:5,12,18)]
colnames(towrite)[3:5] <- c("EUR","AFR","NAM")
write.table(towrite,file="olga_admixture_xchr_results.txt")

library(ggplot2); library(reshape)

set <- res[order(res[,3],res[,4],res[,5]),]
set <- set[set$pop=="-",c("AFR_pruned","NAM_pruned","EUR_pruned")]

set$id <- 1:nrow(set)
set <- melt(set,id="id")

pdf("xchr_barplot_admixture.pdf",width=14)
ggplot(set,aes(id,value,fill=variable)) + geom_bar(stat="identity") +
  ylab("X Chromosome Ancestry Proportion") + xlab("Sample") + 
  ggtitle("ADMIXTURE Estimated X Chr Ancestry Proportions") + theme_bw()
dev.off()

# merge in subpop
scan <- getobj("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_scanAnnot_v54_AMS.RData")

res <- merge(res,pData(scan)[,c("scanID","gengrp6")],by="scanID",all.x=TRUE)
res$gengrp6[is.na(res$gengrp6)] <- "Other"
set <- res[order(res[,3],res[,4],res[,5]),]
set <- set[set$pop=="-",c("AFR_pruned","NAM_pruned","EUR_pruned","gengrp6")]

set$id <- 1:nrow(set)
set <- melt(set,id=c("id","gengrp6"))

pdf("xchr_barplot_admixture_byGengrp6.pdf",width=14)
p =ggplot(set,aes(x=id,y=value,fill=variable)) 
p  + geom_bar(stat="identity",position="stack") + facet_grid(.~gengrp6,scales="free",space="free")+
  ylab("X Chromosome Ancestry Proportion") + xlab("Sample") + 
  ggtitle("ADMIXTURE Estimated X Chr Ancestry Proportions") + theme_bw()
dev.off()

# compare to the autosomal estimates
auto <- read.table("/projects/geneva/gcc-fs2//OLGA/genotype/native_american/reich_data/stephanie/results/olga/1000G/admixture/olga_refpanel_k3/olga_admixture_results.txt",
                   header=TRUE,as.is=TRUE)

colnames(auto)[8:10] <- c("EUR_pruned","AFR_pruned","NAM_pruned")
colnames(auto)[colnames(auto)=="pop"] <- "gengrp6"
auto$genome <- "auto"
res$genome <- "x"
  
res <- merge(res,auto,by="scanID")
res$genome.x[is.na(res$genome.x)] <- res$genome.y[is.na(res$genome.x)]
res$gengrp6.x[is.na(res$gengrp6.x)] <- res$gengrp6.y[is.na(res$gengrp6.x)]
colnames(res)[colnames(res)=="gengrp6.x"] <- "gengrp6"
colnames(res)[colnames(res)=="genome.x"] <- "genome"
set <- res[order(res[,3],res[,4],res[,5]),]
set <- set[set$pop=="-",c("AFR_pruned.x","NAM_pruned.x","EUR_pruned.x",
                          "AFR_pruned.y","NAM_pruned.y","EUR_pruned.y","gengrp6","genome")]


set$id <- 1:nrow(set)
set <- melt(set,id=c("id","gengrp6","genome"))


pdf("xchrVsAuto_barplot_admixture_byGengrp6.pdf",width=14)
p =ggplot(set,aes(x=id,y=value,fill=variable)) 
p  + geom_bar(stat="identity",position="stack") + facet_grid(genome~gengrp6,scales="free",space="free")+
  ylab("Ancestry Proportion") + xlab("Sample") + 
  ggtitle("ADMIXTURE Estimated Autosomal and X Chr Ancestry Proportions") + theme_bw()
dev.off()

# make histograms and density plots of the difference between x chr and autosomal ancestry proportions
res$diffEUR <- res$EUR_pruned.x-res$EUR_pruned.y
res$diffAFR <- res$AFR_pruned.x-res$AFR_pruned.y
res$diffNAM <- res$NAM_pruned.x-res$NAM_pruned.y

apply(res[,c("diffEUR","diffAFR","diffNAM")],2,summary)
#          diffEUR   diffAFR   diffNAM
#Min.    -0.902800 -0.587700 -0.673100
#1st Qu. -0.237100 -0.036720 -0.007242
#Median  -0.116900  0.004255  0.079340
#Mean    -0.122400  0.029230  0.093130
#3rd Qu. -0.004291  0.078270  0.189100
#Max.     0.631700  0.864900  0.740000

set <- data.frame(diff=c(res$diffEUR,res$diffAFR,res$diffNAM),
                  ethn=c(rep("EUR",nrow(res)),rep("AFR",nrow(res)),rep("NAM",nrow(res))))

pdf("xchrVsAuto_hist.pdf")
ggplot(set,aes(x=diff,y=..density..)) + theme_bw() + geom_vline(xintercept=0,color="red") +
geom_histogram(color="black",fill="white") + geom_density() +
facet_grid(.~ethn) + xlab("X Chr - Autosomal Ancestry Proportion")
dev.off()

pdf("xchrVsAuto_boxplot.pdf")
ggplot(set,aes(x=ethn,y=diff)) + theme_bw() + 
  geom_boxplot(color="black",fill="white") +
   xlab("X Chr - Autosomal Ancestry Proportion") + ylab("X Chr - Autosomal Ancestry Proportion")
dev.off()

##### 
# compare the results using different SNP sets
library(GGally)

pdf("pairs_nam.pdf")
ggpairs(res[,c("NAM_pruned.x","NAM_all","NAM_pruned048")],diag=list(continuous="bar"))
dev.off()

pdf("pairs_eur.pdf")
ggpairs(res[,c("EUR_pruned.x","EUR_all","EUR_pruned048")],diag=list(continuous="bar"))
dev.off()

pdf("pairs_afr.pdf")
ggpairs(res[,c("AFR_pruned.x","AFR_all","AFR_pruned048")],diag=list(continuous="bar"))
dev.off()

rm(list=ls())


#####
# 46. Plot X chr PC by ADMIXTURE proportions

library(RColorBrewer)
library(ggplot2)
library(reshape2)
library(grid)
library(OLGApipeline)

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x/olga_application/admixture/genetic_diversity")
ref.pops <- c("Africa", "America", "Europe")
pals <- setNames(c("BuPu", "OrRd", "BuGn"), ref.pops)

admix.xchr <- read.csv("plots/data/admixture_xchr_k3.csv", as.is=TRUE, header=TRUE)
evs <- getobj("../../pca_prunedXsnps_10287unrel12784tot.RData")
evs <- evs$vectors
evs <- data.frame(evs)
evs$scanID <- rownames(evs)
evs <- merge(admix.xchr, evs, by="scanID")
colnames(evs)[6:15] <- paste0("EV",1:10)

for (p in ref.pops) {
  for (y in c("EV2", "EV3")) {
    ggplot(evs, aes_string(x="EV1", y=y, colour=p)) + geom_point() + coord_fixed() +
      scale_colour_gradientn(colours=brewer.pal(5, pals[p])) +
      theme(panel.background=element_rect(fill="gray70"))
    ggsave(paste0("plots/plots/pca_heatmap_xchr_", p, "_", y, ".pdf"), width=6, height=6)
  }
}

rm(list=ls())


#####
# 47. Run ADMIXTURE by chromosome

# first run in batch prune and subset gds file to pruned SNPs by chr
# then call ADMIXTURE
# loop through all chrs
# cd /projects/geneva/geneva_sata/caitlin/mlm_x/olga_application/admixture
# qsub -N admixChr batch_admixture_byChr.sh

# num snps for each chr:
# Chromosome 1: 22.51%, 7557/33576
# Chromosome 2: 20.97%, 7421/35389
# Chromosome 3: 21.94%, 6450/29401
# Chromosome 4: 22.42%, 5795/25846
# Chromosome 5: 21.81%, 5911/27108
# Chromosome 6: 20.69%, 5990/28947
# Chromosome 7: 22.50%, 5263/23392
# Chromosome 8: 19.99%, 4951/24773
# Chromosome 9: 21.15%, 4479/21176
# Chromosome 10: 21.48%, 4996/23257
# Chromosome 11: 21.49%, 4605/21429
# Chromosome 12: 22.15%, 4684/21151
# Chromosome 13: 21.85%, 3620/16570
# Chromosome 14: 22.33%, 3268/14635
# Chromosome 15: 22.90%, 3083/13463
# Chromosome 16: 24.08%, 3249/13495
# Chromosome 17: 25.93%, 3005/11587
# Chromosome 18: 23.15%, 3168/13684
# Chromosome 19: 30.28%, 2265/7480
# Chromosome 20: 24.07%, 2827/11745
# Chromosome 21: 24.23%, 1575/6499
# Chromosome 22: 24.86%, 1626/6540


#####
# 48. CAnD on the x vs auto proportions

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x/olga_application/admixture")
library(GWASTools)

res <- read.table("olga_admixture_xchr_results.txt",header=TRUE,as.is=TRUE)
colnames(res)[3:5] <- c("EUR.x","AFR.x","NAM.x")
res <- res[,c(1,2,6,7,3,4,5)]

for(i in 1:22){
  fn <- paste0("olga_ref_unrel_chr",i,"pruned.3.Q")
  dat <- data.frame(read.table(fn))
  fam <- read.table(paste0("olga_ref_unrel_chr",i,"pruned.fam"),as.is=TRUE)
  pop <- read.table(paste0("olga_ref_unrel_chr",i,"pruned.pop"),as.is=TRUE)
  
  # identify which cols are eur, afr, nam
  eurS <- head(dat[pop$V1=="EUR",1:3])
  colnames(dat)[which(eurS[1,]>0.9)] <- paste0("EUR.",i)
  
  namS <- head(dat[pop$V1=="NAM",1:3])
  colnames(dat)[which(namS[1,]>0.9)] <- paste0("NAM.",i)
  
  afrS <- head(dat[pop$V1=="AFR",1:3])
  colnames(dat)[which(afrS[1,]>0.9)] <- paste0("AFR.",i)
    
  dat$scanID <- fam[,"V2"]  
  res <- merge(res,dat,by="scanID",all.x=TRUE)
}

scan <- getobj("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_scanAnnot_v54_AMS.RData")
res <- merge(res,pData(scan)[,c("scanID","unrelated.pcair.deg4")],by="scanID",all.x=TRUE)

# merge in unrelated set of samples
save(res,file="olga_admixture_results.RData")

rm(list=ls())

## run this locally, now
# be sure to filter out all reference samples

install.packages("~/Documents/CAnD_package/CAnD_0.99.0.tar.gz", repos = NULL, type = "source")
library(CAnD)
library(GWASTools)

# compare with the set of samples from admixture plots on the server
oldRes <-getobj("admix_xchr_aftermerge.RData")
res <- getobj("olga_admixture_results.RData")

sum(is.element(oldRes$scanID,res$scanID)) # 10642
sum(!is.element(res$scanID,oldRes$scanID)) # 84

res <- merge(res,oldRes,by="scanID",all.y=TRUE)

res <- res[res$pop=="-",]
res <- res[res$unrelated.pcair.deg4,]
dim(res) # 10073 74

# first do by all samples together
afrCols <- seq(from=6,to=73,by=3)
colnames(res)[afrCols]
afr <- CAnD(res[,afrCols])
colnames(res)[afrCols-1]
eur <- CAnD(res[,(afrCols-1)])
colnames(res)[afrCols+1]
nam <- CAnD(res[,(afrCols+1)])

res$AFR.auto <- rowMeans(res[,afrCols[2:23]])
res$EUR.auto <- rowMeans(res[,(afrCols-1)[2:23]])
res$NAM.auto <- rowMeans(res[,(afrCols+1)[2:23]])

t.test(res$AFR.x,res$AFR.auto) # p-value = 3.767e-08
t.test(res$AFR.x,res$AFR.auto,paired=TRUE) # p-value < 2.2e-16

t.test(res$EUR.x,res$EUR.auto) # p-value < 2.2e-16
t.test(res$EUR.x,res$EUR.auto,paired=TRUE) # p-value < 2.2e-16

t.test(res$NAM.x,res$NAM.auto) # p-value < 2.2e-16
t.test(res$NAM.x,res$NAM.auto,paired=TRUE) # p-value < 2.2e-16

plotPvals(nam,title="Native American CAnD Results")
# wow, just way off the charts

# do CAnD but excl x chr
afrA <- CAnD(res[,(afrCols)[2:23]])
colnames(res)[(afrCols-1)[2:23]]
eurA <- CAnD(res[,(afrCols-1)[2:23]])
colnames(res)[(afrCols+1)[2:23]]
namA <- CAnD(res[,(afrCols+1)[2:23]])


# now stratify by bkgrd value
table(res$bkgrd)
#CentralAmerican           Cuban       Dominican         Mexican           Other     PuertoRican 
#           1094            1732             897            3635             306            1703 
#SouthAmerican         Unknown 
#          686              20 

eurCols <- afrCols-1
namCols <- afrCols+1

res$bkgrd[res$bkgrd=="Unknown"] <- "Other/Unknown"
res$bkgrd[res$bkgrd=="Other"] <- "Other/Unknown"

subgrpRes <- data.frame("bkgrd"=unique(res$bkgrd),"n"=NA)
newCols <- data.frame(matrix(NA,nrow=nrow(subgrpRes),ncol=23*3))
colnames(newCols) <- paste0(paste0("chr",c(23,1:22)),rep(c(".NAM",".EUR",".AFR"),each=23))
subgrpRes <- cbind(subgrpRes,newCols)
colnames(subgrpRes)

for(i in 1:nrow(subgrpRes)){
  pop <- subgrpRes$bkgrd[i]
  subgrpRes$n[i] <- sum(is.element(res$bkgrd,pop))
  
  subgrpRes[i,3] <- pValues(CAnD(res[res$bkgrd==pop,namCols]))[1]
  subgrpRes[i,4:25] <- pValues(CAnD(res[res$bkgrd==pop,namCols[2:23]]))
  
  subgrpRes[i,26] <- pValues(CAnD(res[res$bkgrd==pop,eurCols]))[1]
  subgrpRes[i,27:48] <- pValues(CAnD(res[res$bkgrd==pop,eurCols[2:23]]))
  
  subgrpRes[i,49] <- pValues(CAnD(res[res$bkgrd==pop,afrCols]))[1]
  subgrpRes[i,50:71] <- pValues(CAnD(res[res$bkgrd==pop,afrCols[2:23]]))  
}
# so the x chr is x vs autosomes, each autosome is that chr vs all other autosomes

# plot these results
# facet wrap for each ancestral subpop
# on the x axis have each chromosome
# on the y axis the -log10(pvalue), where the plotting colors are by bkgrd group, so 6 points per x-axis value

melted <- melt(subgrpRes,id.vars=c("bkgrd","n"))
melted$Ancestry <- substr(melted$variable,start=regexpr("\\.",melted$variable)+1,stop=nchar(as.character(melted$variable)))
melted$Chromosome <- substr(melted$variable,start=regexpr("chr",melted$variable)+3,stop=regexpr("\\.",melted$variable)-1)

#melted$Chromosome[melted$Chromosome==23] <- "X"
melted$Ancestry[melted$Ancestry=="NAM"] <- "America"
melted$Ancestry[melted$Ancestry=="EUR"] <- "Europe"
melted$Ancestry[melted$Ancestry=="AFR"] <- "Africa"

melted$bkgrd <- paste0(melted$bkgrd,"\nn=",melted$n)
cols.all <- setNames(c("red","green","gold","blue","gray","magenta","brown"),unique(melted$bkgrd))

melted$pch <- 1
melted$pch[-log10(melted$value)>30] <- 2
melted$value[-log10(melted$value)>30] <- 1e-29

ggplot(melted,aes(x=as.integer(Chromosome),y=-log10(value),fill=bkgrd)) +geom_point(aes(color=bkgrd,shape=factor(pch)))+
  facet_wrap(~Ancestry,nrow=3)+xlab("Chromosome")+theme_bw()+coord_cartesian(ylim=c(-1,30))+
  geom_hline(yintercept=1.30,size=0.7,color="gray")+scale_color_manual(values=cols.all)
ggsave("CAnD_solResults.pdf", width=11, height=8)

## want the pvalues for x chr vs autosomes
xchr <- subgrpRes[,c("bkgrd","n",paste0("chr23.",c("AFR","NAM","EUR")))]

# get mean autosomal ancestry for pooled t-test
res$EUR.auto <- rowMeans(res[,paste0("EUR.",1:22)])
res$NAM.auto <- rowMeans(res[,paste0("NAM.",1:22)])
res$AFR.auto <- rowMeans(res[,paste0("AFR.",1:22)])

xchr$pooled.AFR <- NA
xchr$pooled.NAM <- NA
xchr$pooled.EUR <- NA

for(i in 1:nrow(xchr)){
  pop <- xchr$bkgrd[i]
  resSm <- res[is.element(res$bkgrd,pop),]
  xchr$chrx.NAM[i] <- t.test(resSm$NAM.x,resSm$NAM.auto,paired=TRUE)$p.value
  xchr$chrx.EUR[i] <- t.test(resSm$EUR.x,resSm$EUR.auto,paired=TRUE)$p.value
  xchr$chrx.AFR[i] <- t.test(resSm$AFR.x,resSm$AFR.auto,paired=TRUE)$p.value
  
  xchr$pooled.NAM[i] <- t.test(resSm$NAM.x,resSm$NAM.auto)$p.value
  xchr$pooled.EUR[i] <- t.test(resSm$EUR.x,resSm$EUR.auto)$p.value
  xchr$pooled.AFR[i] <- t.test(resSm$AFR.x,resSm$AFR.auto)$p.value
}

library(xtable)
xchr[,6:11] <- format(xchr[,6:11],digits=3)
xtable(t(xchr[,c(1:2,9:11,7,8,6)]),digits=3)

toSv <- t(xchr[,c(1:2,9:11,7,8,6)])
write.csv(toSv,file="CAnD_results_xvsAuto.csv",quote=FALSE)

rm(list=ls())


#####
# 49. Run assoc with LABA2 BCC trait, excluding all PCs

## first, est VC
# see script in
# /projects/geneva/geneva_sata/caitlin/olga_xchr_assoc/run_assoc_estVC_noPCs.R
# qsub -q thornton.q -N ana316987 batch_run_assoc_estVC.sh
# writes run_assoc_estVC_316987.Rout
## this includes adj for: auto PC 1-5; auto KC, block group and HH as random effects

# cd /projects/geneva/geneva_sata/caitlin/olga_xchr_assoc
# ran run_assoc_mixed_segment_noPCs_316987.R: adj for auto kc, no pcs (auto nor x chr)

# qsub -q thornton.q -N chr1 -t 1-45 batch_run_assoc_mixed_segment.sh
# qsub -q thornton.q -N chr2 -t 1-47 batch_run_assoc_mixed_segment2.sh
# qsub -q thornton.q -N chr3 -t 1-37 batch_run_assoc_mixed_segment3.sh
# qsub -q thornton.q -N chr4 -t 1-36 batch_run_assoc_mixed_segment4.sh
# qsub -q thornton.q -N chr5 -t 1-35 batch_run_assoc_mixed_segment5.sh
# qsub -q thornton.q -N chr6 -t 1-32 batch_run_assoc_mixed_segment6.sh
# qsub -q thornton.q -N chr7 -t 1-30 batch_run_assoc_mixed_segment7.sh
# qsub -q thornton.q -N chr8 -t 1-27 batch_run_assoc_mixed_segment8.sh
# qsub -q thornton.q -N chr9 -t 1-22 batch_run_assoc_mixed_segment9.sh
# qsub -q thornton.q -N chr10 -t 1-25 batch_run_assoc_mixed_segment10.sh
# qsub -q thornton.q -N chr11 -t 1-26 batch_run_assoc_mixed_segment11.sh
# qsub -q thornton.q -N chr12 -t 1-24 batch_run_assoc_mixed_segment12.sh
# qsub -q thornton.q -N chr13 -t 1-19 batch_run_assoc_mixed_segment13.sh
# qsub -q thornton.q -N chr14 -t 1-17 batch_run_assoc_mixed_segment14.sh
# qsub -q thornton.q -N chr15 -t 1-16 batch_run_assoc_mixed_segment15.sh
# qsub -q thornton.q -N chr1617 -t 1-15 batch_run_assoc_mixed_segment1617.sh
# qsub -q thornton.q -N chr18 -t 1-14 batch_run_assoc_mixed_segment18.sh
# qsub -q thornton.q -N chr1920 -t 1-11 batch_run_assoc_mixed_segment1920.sh
# qsub -q thornton.q -N chr2122 -t 1-7 batch_run_assoc_mixed_segment2122.sh
# qsub -q thornton.q -N chr23 -t 1-28 batch_run_assoc_mixed_segment23.sh

# and combine_assoc.R afterward
# qsub -q thornton.q -N chr1 -t 1-23 batch_combine_assoc.sh

# stored: cd Assoc/
# assoc_noPCs_316987_chr23.RData ... for all chrs

# called
# qsub -q thornton.q -N plot batch_plot_qq_manh.sh


#####
# 50. CAnD on autosomes using local ancestry estimates

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")
library(GWASTools)
library(gdsfmt)

# load in local ancestry estimates
gdsfn <- "/projects/geneva/gcc-fs2/OLGA/genotype/native_american/olga_with_refpanel/reich_1000G/local_ancestry/gds/unique/olga_reich_1000G_lai_new_unique.gds"
gdsobj <- openfn.gds(gdsfn)

# need to get the mean values across chromosomes
# figure out how to loop through and take means for each ancestry, by chr
chrs <- read.gdsn(index.gdsn(gdsobj, "snp.chromosome"))
length(chrs) # 14192
(t <- table(chrs))
#   1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
#1116 1072  939  874  846  811  760  688  661  720  665  709  520  479  484  514 
# 17   18   19   20   21   22 
#506  488  402  438  245  255 
stInd <- cumsum(t)

# ok, so the ancestry estimates are matrices of sample x snp
# store the mean ancestry by chromosome per sample, so need a matrix of results sample x chr*3
ancest <- data.frame(matrix(NA,nrow=12803,ncol=(22*3+1)))
colnames(ancest) <- c("scanID",paste0(rep(c("eur","afr","amer"),each=22),1:22))

ancest$scanID <- read.gdsn(index.gdsn(gdsobj,"sample.id"))

for(i in 1:22){
  if(i==1){st <- 1} else{st <- stInd[i-1]+1}
  stopifnot(all(chrs[st:(t[i]+st-1)]==i))
  tmpEur <- read.gdsn(index.gdsn(gdsobj,"dosage_eur"),start=c(1,st),count=c(-1,t[i]))
  # want row sums, divided by 2*ncol
  colnm <- paste0("eur",i)
  ancest[,colnm] <- rowSums(tmpEur)/(2*ncol(tmpEur))
  
  tmpAfr <- read.gdsn(index.gdsn(gdsobj,"dosage_afr"),start=c(1,st),count=c(-1,t[i]))
  colnm <- paste0("afr",i)
  ancest[,colnm] <- rowSums(tmpAfr)/(2*ncol(tmpAfr))
  
  tmpNam <- read.gdsn(index.gdsn(gdsobj,"dosage_amer"),start=c(1,st),count=c(-1,t[i]))
  colnm <- paste0("amer",i)
  ancest[,colnm] <- rowSums(tmpNam)/(2*ncol(tmpNam)) 
}

# check that each ancestry sums to 1
colseq <- seq(from=2,to=ncol(ancest),by=22)
for(i in 1:22){
  stopifnot(sum(ancest[,(colseq+i-1)])==12803)
} # ok, looks good!

# save this and run the rest locally
save(ancest,file="local_ancestry_byChr.RData")

##
# now we can do CAnD on the local ancestry proportions
# want to stratify by background group

install.packages("~/Documents/CAnD_package/CAnD_0.99.0.tar.gz", repos = NULL, type = "source")
library(CAnD)
library(GWASTools)

# compare with the set of samples from admixture plots on the server
oldRes <- getobj("admix_xchr_aftermerge.RData")
res <- getobj("olga_admixture_results.RData")

sum(is.element(oldRes$scanID,res$scanID)) # 10642
sum(!is.element(res$scanID,oldRes$scanID)) # 84

res <- merge(res,oldRes,by="scanID",all.y=TRUE)

res <- res[res$pop=="-",]
res <- res[res$unrelated.pcair.deg4,]
dim(res) # 10073 78

# first do by all samples together
afrCols <- seq(from=6,to=73,by=3)
colnames(res)[afrCols]
afr <- CAnD(res[,afrCols])
colnames(res)[afrCols-1]
eur <- CAnD(res[,(afrCols-1)])
colnames(res)[afrCols+1]
nam <- CAnD(res[,(afrCols+1)])

res$AFR.auto <- rowMeans(res[,afrCols[2:23]])
res$EUR.auto <- rowMeans(res[,(afrCols-1)[2:23]])
res$NAM.auto <- rowMeans(res[,(afrCols+1)[2:23]])

t.test(res$AFR.x,res$AFR.auto)$p.value # 3.369035e-08
t.test(res$AFR.x,res$AFR.auto,paired=TRUE)$p.value # 7.641825e-33

t.test(res$EUR.x,res$EUR.auto)$p.value # 1.644846e-169
t.test(res$EUR.x,res$EUR.auto,paired=TRUE)$p.value # 0

t.test(res$NAM.x,res$NAM.auto)$p.value # 5.53147e-98
t.test(res$NAM.x,res$NAM.auto,paired=TRUE)$p.value # 0

plotPvals(nam,title="Native American CAnD Results")
# wow, just way off the charts

# do CAnD but excl x chr
afrA <- CAnD(res[,(afrCols)[2:23]])
colnames(res)[(afrCols-1)[2:23]]
eurA <- CAnD(res[,(afrCols-1)[2:23]])
colnames(res)[(afrCols+1)[2:23]]
namA <- CAnD(res[,(afrCols+1)[2:23]])

# now stratify by bkgrd value
table(res$bkgrd)
#CentralAmerican           Cuban       Dominican         Mexican           Other     PuertoRican 
#           1094            1732             897            3635             306            1703 
#SouthAmerican         Unknown 
#          686              20 

eurCols <- afrCols-1
namCols <- afrCols+1

res$bkgrd[res$bkgrd=="Unknown"] <- "Other/Unknown"
res$bkgrd[res$bkgrd=="Other"] <- "Other/Unknown"

subgrpRes <- data.frame("bkgrd"=unique(res$bkgrd),"n"=NA)
newCols <- data.frame(matrix(NA,nrow=nrow(subgrpRes),ncol=23*3))
colnames(newCols) <- paste0(paste0("chr",c(23,1:22)),rep(c(".NAM",".EUR",".AFR"),each=23))
subgrpRes <- cbind(subgrpRes,newCols)
colnames(subgrpRes)

for(i in 1:nrow(subgrpRes)){
  pop <- subgrpRes$bkgrd[i]
  subgrpRes$n[i] <- sum(is.element(res$bkgrd,pop))
  
  subgrpRes[i,3] <- pValues(CAnD(res[res$bkgrd==pop,namCols]))[1]
  subgrpRes[i,4:25] <- pValues(CAnD(res[res$bkgrd==pop,namCols[2:23]]))
  
  subgrpRes[i,26] <- pValues(CAnD(res[res$bkgrd==pop,eurCols]))[1]
  subgrpRes[i,27:48] <- pValues(CAnD(res[res$bkgrd==pop,eurCols[2:23]]))
  
  subgrpRes[i,49] <- pValues(CAnD(res[res$bkgrd==pop,afrCols]))[1]
  subgrpRes[i,50:71] <- pValues(CAnD(res[res$bkgrd==pop,afrCols[2:23]]))  
}
# so the x chr is x vs autosomes, each autosome is that chr vs all other autosomes

# plot these results
# facet wrap for each ancestral subpop
# on the x axis have each chromosome
# on the y axis the -log10(pvalue), where the plotting colors are by bkgrd group, so 6 points per x-axis value

melted <- melt(subgrpRes,id.vars=c("bkgrd","n"))
melted$Ancestry <- substr(melted$variable,start=regexpr("\\.",melted$variable)+1,stop=nchar(as.character(melted$variable)))
melted$Chromosome <- substr(melted$variable,start=regexpr("chr",melted$variable)+3,stop=regexpr("\\.",melted$variable)-1)

#melted$Chromosome[melted$Chromosome==23] <- "X"
melted$Ancestry[melted$Ancestry=="NAM"] <- "America"
melted$Ancestry[melted$Ancestry=="EUR"] <- "Europe"
melted$Ancestry[melted$Ancestry=="AFR"] <- "Africa"

melted$bkgrd <- paste0(melted$bkgrd,"\nn=",melted$n)
cols.all <- setNames(c("red","green","gold","blue","gray","magenta","brown"),unique(melted$bkgrd))

melted$pch <- 1
melted$pch[-log10(melted$value)>30] <- 2
melted$value[-log10(melted$value)>30] <- 1e-29

ggplot(melted,aes(x=as.integer(Chromosome),y=-log10(value),fill=bkgrd)) +geom_point(aes(color=bkgrd,shape=factor(pch)))+
  facet_wrap(~Ancestry,nrow=3)+xlab("Chromosome")+theme_bw()+coord_cartesian(ylim=c(-1,30))+
  geom_hline(yintercept=1.30,size=0.7,color="gray")+scale_color_manual(values=cols.all)
ggsave("CAnD_solResults.pdf", width=11, height=8)

rm(list=ls())


#####
# 51. Summaries for CAnD; non param CAnD


#####
# 52. New PO kc estimates, auto vs x chr

library(GWASTools)
setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")

scan <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_scanAnnot_v54_AMS.RData"))
head(pData(scan))

kcX_xAdj <- get(load("olga_application/pcRelate_prunedXchr_xPC12adj.RData"))
kcX_xAdj <- kcX_xAdj$kinship

obsRel <- get(load("pcrelate_PO_estimates.RData"))
head(obsRel)

poIds1 <- obsRel$ID1
poIds2 <- obsRel$ID2
obsRel$ID12 <- paste(obsRel$ID1,obsRel$ID2)

kcPO_xAdj <- kcX_xAdj[is.element(kcX_xAdj$ID1,poIds1),]
obsPOids <- paste(kcPO_xAdj$ID1,kcPO_xAdj$ID2)
kcPO_xAdj <- kcPO_xAdj[is.element(obsPOids,paste(poIds1,poIds2)),]
dim(kcPO_xAdj) # 1442 4
kcPO_xAdj$ID12 <- paste(kcPO_xAdj$ID1,kcPO_xAdj$ID2)

kcPO_xAdj <- merge(kcPO_xAdj,obsRel,by="ID12",all.x=TRUE,suffixes=c(".x",".auto"))
# need to merge in sex, too
kcPO_xAdj <- merge(kcPO_xAdj,pData(scan)[,c("scanID","sex")],by.x="ID1.x",by.y="scanID")
kcPO_xAdj <- merge(kcPO_xAdj,pData(scan)[,c("scanID","sex")],by.x="ID2.x",by.y="scanID",suffixes=c(".1",".2"))

pdf("olga_application/kc_xPrunedvsAuto_poPairs_xAdj.pdf")
plot(kcPO_xAdj$kin.x,kcPO_xAdj$kin.auto,type="n",xlab="X Chr KC",ylab="Auto KC",
     main="Auto vs X Chr KC for 1,442 PO pairs",cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
points(kcPO_xAdj$kin.x[kcPO_xAdj$sex.1=="M"&kcPO_xAdj$sex.2=="M"],
       kcPO_xAdj$kin.auto[kcPO_xAdj$sex.1=="M"&kcPO_xAdj$sex.2=="M"],col=rgb(0,1,0,0.5),pch=19)
points(kcPO_xAdj$kin.x[kcPO_xAdj$sex.1=="F"&kcPO_xAdj$sex.2=="F"],
       kcPO_xAdj$kin.auto[kcPO_xAdj$sex.1=="F"&kcPO_xAdj$sex.2=="F"],col=rgb(1,0,0,0.5),pch=19)
points(kcPO_xAdj$kin.x[kcPO_xAdj$sex.1=="M"&kcPO_xAdj$sex.2=="F"],
       kcPO_xAdj$kin.auto[kcPO_xAdj$sex.1=="M"&kcPO_xAdj$sex.2=="F"],col=rgb(0,0,1,0.5),pch=19)
points(kcPO_xAdj$kin.x[kcPO_xAdj$sex.1=="F"&kcPO_xAdj$sex.2=="M"],
       kcPO_xAdj$kin.auto[kcPO_xAdj$sex.1=="F"&kcPO_xAdj$sex.2=="M"],col=rgb(0,0,1,0.5),pch=19)
abline(v=0,col="gray"); abline(v=0.25,col="gray"); abline(v=0.5,col="gray")
legend("topleft",c("M-M pairs","F-F pairs","F-M pairs"),col=c(rgb(0,1,0,1),rgb(1,0,0,1),rgb(0,0,1,1)),
       pch=19,cex=1.3,bg="white")
dev.off()

rm(list=ls())


#####
# 53. Recalculate lambda on X chr excl Xq28 region

library(QCpipeline)
library(GWASTools)
source("CAnD_pooled.R")

setwd("/projects/geneva/geneva_sata/caitlin/olga_xchr_assoc/Assoc/")

# read in xchr results, the autosomal model
assoc <- getobj("assoc_auto_316987_chr23.RData")
(mac.thresh <- 30)
filt.imp <- which(assoc$type %in% 0)
if (!("composite.filter" %in% names(assoc))) assoc$composite.filter <- TRUE
assoc$maf.filt <- assoc$effN > mac.thresh
n <- max(assoc$n, na.rm=T)
n # 12489

maf.eff <- format(quadSolveMAF(mac.thresh, n), digits=3)
filt.obs <- which(assoc$type %in% c(2,3))

filt.obs.maf <- intersect(filt.obs, which(assoc$maf.filt))
title.obs.maf <- paste("obs + composite.filter + MAC >", mac.thresh, "\n", prettyNum(length(filt.obs.maf), big.mark=","), "SNPs - MAF >", maf.eff) 
# this is the filter i want

assoc <- assoc[filt.obs.maf,]

stat = assoc[,"Stat"]
stat <- pmax(stat, 0)

(lambda <- calculateLambda(stat, df=1)) # 1.115727, as expected

# figure out what position we should truncate the markers at
# take off the last 1e7 bpairs
maxPos <- max(assoc$position)-1e7
assocSm <- assoc[assoc$position<maxPos,]
stat <- assocSm[["Stat"]]
stat <- pmax(stat,0)
length(stat) # 41718
(lambda <- calculateLambda(stat,df=1)) # 1.104407

##
# do the same exercise for the other assoc test results
assoc <- getobj("assoc_316987_chr23.RData")
(mac.thresh <- 30)
filt.imp <- which(assoc$type %in% 0)
if (!("composite.filter" %in% names(assoc))) assoc$composite.filter <- TRUE
assoc$maf.filt <- assoc$effN > mac.thresh
assoc$maf2.filt <- assoc$effN > 2*mac.thresh
n <- max(assoc$n, na.rm=T)
n # 12489

maf.eff <- format(quadSolveMAF(mac.thresh, n), digits=3)

filt.obs.maf <- intersect(filt.obs, which(assoc$maf.filt))
title.obs.maf <- paste("obs + composite.filter + MAC >", mac.thresh, "\n", prettyNum(length(filt.obs.maf), big.mark=","), "SNPs - MAF >", maf.eff) 
# this is the filter i want

assoc <- assoc[filt.obs.maf,]

stat = assoc[,"Stat"]
stat <- pmax(stat, 0)

(lambda <- calculateLambda(stat, df=1)) # 1.04343

# figure out what position we should truncate the markers at
# take off the last 1e7 bpairs
maxPos <- max(assoc$position)-1e7
assocSm <- assoc[assoc$position<maxPos,]
stat <- assocSm[,"Stat"]
stat <- pmax(stat,0)
length(stat) # 41718
(lambda <- calculateLambda(stat,df=1)) # 1.035676

assoc <-  getobj("assoc_316987_chr23.RData")
length(assoc$pval)

assoc$maf.filt <- assoc$effN > mac.thresh
qual.filt <- which(assoc > 0.8)
mafhi.filt <- intersect(qual.filt, which(assoc$maf.filt))
length(mafhi.filt) # 669969, just as we need!

pos.filt <- which(assoc$position<maxPos)

mafhi.pos.filt <- intersect(mafhi.filt,pos.filt)
length(mafhi.pos.filt) # 626933

# make a plot of these
outfile <- paste("../Plots/pval_manh_single__exclXq28_316987",".png", sep="")
     png(outfile, width=1080, height=380)
mtxt <- "assoc~evx1+evx2"
     par(mar=c(5,5,4,3)+0.1, lwd=1.5, cex.lab=1.5, cex.main=1.5, cex.sub=1.5, 
         oma=c(0,0,length(mtxt)+0.1,0))
     manhattanPlot(assoc[,"pval"][mafhi.pos.filt], assoc$chromosome[mafhi.pos.filt], 
                   main="all + composite.filter + MAC*oevar > 30\n626,933 SNPs - MAF > 0.0012")
     #mtext(side=3, line=-1, text=groupTitle, padj=0.9, adj=0.02, outer=T, cex=1.5)
     #title(mtxt, outer=T)
     dev.off()

rm(list=ls())


#####
# 54. Redo plots with local ancestry estimates

library(GWASTools)

res <- getobj("local_ancestry_byChr.RData")
scan <- getobj("admix_xchr_aftermerge.RData")
colnames(scan)[3:5] <- c("AFR.x","NAM.x","EUR.x")               
# merge in admixture results for x chr since don't have local ancestry estimates
res <- merge(res,scan,by="scanID")

# get avg local ancestry across the autosomes
colnames(res)
afrCols <- 24:45
colnames(res)[afrCols]
res$AFR.auto <- rowMeans(res[,afrCols])

eurCols <- 2:23
colnames(res)[eurCols]
res$EUR.auto <- rowMeans(res[,eurCols])

namCols <- 46:67
colnames(res)[namCols]
res$NAM.auto <- rowMeans(res[,namCols])

# plot the local ancestry estimates averaged across autosomes in barplot form
# want res$AFR.auto,EUR.auto,NAM.auto for each bkgrd
factorWithCount <- function(x,levels){
  n <- sapply(levels, function(l) sum(x == l, na.rm=TRUE))
  factor(x, levels=levels, labels=paste(levels, n, sep="\nn="))
}
pops <- c("Cuban", "Dominican", "PuertoRican", "Mexican", "CentralAmerican", "SouthAmerican", "Other/Unknown")
res$bkgrd[is.element(res$bkgrd,c("Unknown","Other"))] <- "Other/Unknown"
res$bkgrd <- factorWithCount(res$bkgrd, levels=pops)

res$n <- 1:nrow(res)
res$EUR.diff <- res$EUR.x-res$EUR.auto
res$AFR.diff <- res$AFR.x-res$AFR.auto
res$NAM.diff <- res$NAM.x-res$NAM.auto

resSm <- res[,c("scanID","bkgrd","n","AFR.auto","NAM.auto","EUR.auto")]
colnames(resSm)[4:6] <- c("Africa","America","Europe")
resSm$genome <- "Autosomal"
resSm <- resSm[order(resSm$bkgrd, resSm$Europe, resSm$America,resSm$Africa),]
resSm$n <- 1:nrow(resSm)

toAd <- res[,c("scanID","bkgrd","n","AFR.x","NAM.x","EUR.x")]
toAd$genome <- "X"
colnames(toAd)[4:6] <- c("Africa","America","Europe")
toAd <- toAd[order(toAd$bkgrd,toAd$Europe,toAd$America,toAd$Africa),]
toAd$n <- 1:nrow(toAd)
resSm <- rbind(resSm,toAd)

melted <- melt(resSm, id.vars=c("n", "scanID", "bkgrd","genome"))

library(RColorBrewer); library(grid)
DARK2 <- setNames(brewer.pal(8, "Dark2"),
                  c("teal", "orange", "purple", "pink", "green", "gold", "brown", "gray"))
levels <- c("Africa", "America", "Europe")
cols <- setNames(DARK2[c("purple", "orange", "teal")], levels)
colnames(melted)[6] <- "Proportion"

png("localAncest_auto_X_barplot.png",width=960)
ggplot(melted, aes(x=n, y=Proportion, fill=variable, color=variable)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=cols, breaks=rev(names(cols))) +
  scale_color_manual(values=cols, breaks=rev(names(cols))) + 
  facet_grid(genome~bkgrd, scales="free_x") +
  theme_bw() + 
  theme(axis.line=element_blank(),
        axis.ticks.x=element_blank(),axis.text=element_text(size=14),
        axis.text.x=element_blank(),legend.text=element_text(size=14),
        axis.title.x=element_blank(),strip.text=element_text(size=15),
        panel.margin=unit(0, "in"))
dev.off()

## make boxplots now
cols.all <- c(setNames(paste0(cols, "55"), paste0(names(cols),".X")),
              setNames(cols, paste0(names(cols),".Autosomal")))
cols.all <- cols.all[c(4,1,5,2,6,3)]

melted$Ancestry.Type <- paste(melted$variable,melted$genome,sep=".")

pdf("localAncest_boxplot_auto_x_woutPvals.pdf",width=13,height=9)
ggplot(melted, aes(x=bkgrd, y=Proportion, fill=Ancestry.Type)) + geom_boxplot() +
  scale_fill_manual(values=cols.all, breaks=names(cols.all)) +
  facet_wrap(~variable, ncol=1) + theme_bw() +
  ylab("Fraction of continental ancestry") + xlab("Self-identified background") +
  #geom_text(data=pvalPlot, aes(x=bkgrd, y=y, label=value), size=4) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=18),
        title=element_text(size=20),strip.text = element_text(size=16),
        legend.text=element_text(size=16)) 
dev.off()

# make a boxplot of the differences, add a line at y=0
resSm <- res[,c("scanID","bkgrd","n","AFR.diff","NAM.diff","EUR.diff")]
colnames(resSm)[4:6] <- c("Africa","America","Europe")
resSm <- resSm[order(resSm$bkgrd, resSm$Europe, resSm$America,resSm$Africa),]
resSm$n <- 1:nrow(resSm)

melted <- melt(resSm, id.vars=c("n", "scanID", "bkgrd"))

pdf("localAncest_boxplot_auto_x_Diff.pdf",width=11,height=9)
ggplot(melted, aes(x=bkgrd, y=value, fill=variable)) + geom_boxplot() +
  scale_fill_manual(values=cols, breaks=names(cols)) +
  facet_wrap(~variable, ncol=1) + 
  ylab("X Chr - Autosomal Average Local Ancestry") + xlab("Self-identified background") +
  #geom_text(data=pvalPlot, aes(x=bkgrd, y=y, label=value), size=4) +
  theme(axis.text.x=element_text(colour="black",size=12),
        axis.text.y=element_text(size=13)) + theme_bw() + 
  geom_hline(yintercept=0,linetype=2,color="gray")
dev.off()

# make boxplots of the avg local ancestry, by chr
library(ggplot2); library(reshape); library(RColorBrewer)
res <- melt(res,c("scanID","bkgrd"))
res$Ancestry <- NA
res$Ancestry[grep("amer",res$variable)] <- "America"
res$Ancestry[grep("eur",res$variable)] <- "Europe"
res$Ancestry[grep("afr",res$variable)] <- "Africa"
res$Ancestry[grep("AFR",res$variable)] <- "Africa"
res$Ancestry[grep("EUR",res$variable)] <- "Europe"
res$Ancestry[grep("NAM",res$variable)] <- "America"

SET1 <- setNames(brewer.pal(9, "Set1"),
                 c("red", "blue", "green", "purple", "orange", "yellow", "brown", "pink", "gray"))

res$bkgrd <- factor(res$bkgrd)

res$bg <- substr(res$bkgrd,start=1,stop=(regexpr("\n",res$bkgrd)-2))
res$bg[res$bg=="CentralAmerica"] <- "CAm"
res$bg[res$bg=="SouthAmerica"] <- "SAm"
res$bg[res$bg=="PuertoRica"] <- "PR"
res$bg[res$bg=="Other/Unknow"] <- "Unk"
res$bg[res$bg=="Dominica"] <- "Dom"
res$bg[res$bg=="Mexica"] <- "Mex"

res$chr <- NA
res$chr <- gsub("eur","",res$variable)
res$chr <- gsub("amer","",res$chr)
res$chr <- gsub("afr","",res$chr)
res$chr <- gsub("AFR.","",res$chr)
res$chr <- gsub("NAM.","",res$chr)
res$chr <- gsub("EUR.","",res$chr)

res <- res[!is.element(res$chr,c("auto","diff","n")),]

res$chr[res$chr=="x"] <- 23

levels <- c(1:23)
col <- setNames(c(SET1[c("green", "orange", "purple", "red", "yellow", "blue", "gray")],
                  SET1[c("green", "orange", "purple", "red", "yellow", "blue", "gray")],
                  SET1[c("green", "orange", "purple", "red", "yellow", "blue", "gray")],
                  SET1[c("green", "orange")]), levels)

resPl <- res[res$bg=="CAm",]
pdf("localAncest_CAm_autos.pdf",width=11,height=9)
ggplot(resPl,aes(x=as.integer(chr),y=value)) + geom_boxplot(aes(color=chr))+
  facet_wrap(~Ancestry,ncol=1)+theme_bw()+
  xlab("Chromosome") + scale_color_manual(values=col,breaks=rev(names(col)),
                                          guide = guide_legend(title = "")) +
  theme(legend.position = "none") +  ylab("Fraction of Average Local Ancestry") +
  ggtitle("Self-Identified Central Americans, n=1163")
dev.off()

resPl <- res[res$bg=="SAm",]
pdf("localAncest_SAm_autos.pdf",width=11,height=9)
ggplot(resPl,aes(x=as.integer(chr),y=value)) + geom_boxplot(aes(color=chr))+
  facet_wrap(~Ancestry,ncol=1)+theme_bw()+
  xlab("Chromosome") + scale_color_manual(values=col,breaks=rev(names(col)),
                                          guide = guide_legend(title = "")) +
  theme(legend.position = "none") +  ylab("Fraction of Average Local Ancestry") +
  ggtitle("Self-Identified South Americans, n=710")
dev.off()

resPl <- res[res$bg=="Cuba",]
pdf("localAncest_Cuba_autos.pdf",width=11,height=9)
ggplot(resPl,aes(x=as.integer(chr),y=value)) + geom_boxplot(aes(color=chr))+
  facet_wrap(~Ancestry,ncol=1)+theme_bw()+
  xlab("Chromosome") + scale_color_manual(values=col,breaks=rev(names(col)),
                                          guide = guide_legend(title = "")) +
  theme(legend.position = "none") +  ylab("Fraction of Average Local Ancestry") +
  ggtitle("Self-Identified Cubans, n=1779")
dev.off()

resPl <- res[res$bg=="Dom",]
pdf("localAncest_Dom_autos.pdf",width=11,height=9)
ggplot(resPl,aes(x=as.integer(chr),y=value)) + geom_boxplot(aes(color=chr))+
  facet_wrap(~Ancestry,ncol=1)+theme_bw()+
  xlab("Chromosome") + scale_color_manual(values=col,breaks=rev(names(col)),
                                          guide = guide_legend(title = "")) +
  theme(legend.position = "none") +  ylab("Fraction of Average Local Ancestry") +
  ggtitle("Self-Identified Dominicans, n=968")
dev.off()

resPl <- res[res$bg=="Mex",]
pdf("localAncest_Mex_autos.pdf",width=11,height=9)
ggplot(resPl,aes(x=as.integer(chr),y=value)) + geom_boxplot(aes(color=chr))+
  facet_wrap(~Ancestry,ncol=1)+theme_bw()+
  xlab("Chromosome") + scale_color_manual(values=col,breaks=rev(names(col)),
                                          guide = guide_legend(title = "")) +
  theme(legend.position = "none") +  ylab("Fraction of Average Local Ancestry") +
  ggtitle("Self-Identified Mexicans, n=3845")
dev.off()

resPl <- res[res$bg=="PR",]
pdf("localAncest_PR_autos.pdf",width=11,height=9)
ggplot(resPl,aes(x=as.integer(chr),y=value)) + geom_boxplot(aes(color=chr))+
  facet_wrap(~Ancestry,ncol=1)+theme_bw()+
  xlab("Chromosome") + scale_color_manual(values=col,breaks=rev(names(col)),
                                          guide = guide_legend(title = "")) +
  theme(legend.position = "none") +  ylab("Fraction of Average Local Ancestry") +
  ggtitle("Self-Identified Puerto Ricans, n=1832")
dev.off()

resPl <- res[res$bg=="Unk",]
pdf("localAncest_Unk_autos.pdf",width=11,height=9)
ggplot(resPl,aes(x=as.integer(chr),y=value)) + geom_boxplot(aes(color=chr))+
  facet_wrap(~Ancestry,ncol=1)+theme_bw()+
  xlab("Chromosome") + scale_color_manual(values=col,breaks=rev(names(col)),
                                          guide = guide_legend(title = "")) +
  theme(legend.position = "none") +  ylab("Fraction of Average Local Ancestry") +
  ggtitle("Self-Identified Other/Unknown, n=345")
dev.off()

rm(list=ls())


#####
# 55. Perform CAnD and nonP CAnD and pooled t-test on local ancestry estimates

install.packages("~/Documents/CAnD_package/CAnD_0.99.4.tar.gz", repos = NULL, type = "source")
library(CAnD)
library(GWASTools)
source("CAnD_pooled.R")

res <- getobj("local_ancestry_byChr.RData")
scan <- getobj("admix_xchr_aftermerge.RData")
colnames(scan)[3:5] <- c("AFR.x","NAM.x","EUR.x")               
# merge in admixture results for x chr since don't have local ancestry estimates
res <- merge(res,scan,by="scanID")

res$bkgrd[is.element(res$bkgrd,c("Other","Unknown"))] <- "Other/Unknown"

## need to subset so that only unrel samples are included
unrel <- getobj("unrelated_pcair_deg4.RData")

sum(is.element(res$scanID,unrel))
sum(!is.element(res$scanID,unrel))

res <- res[is.element(res$scanID,unrel),]

# get avg local ancestry across the autosomes
colnames(res)
afrCols <- 24:45
colnames(res)[afrCols]
res$AFR.auto <- rowMeans(res[,afrCols])

eurCols <- 2:23
colnames(res)[eurCols]
res$EUR.auto <- rowMeans(res[,eurCols])

namCols <- 46:67
colnames(res)[namCols]
res$NAM.auto <- rowMeans(res[,namCols])


## get CAnD and nonParam CAnD pvalues, also pooled t-test pvalues
# make a table of the results
subgrpRes <- data.frame("bkgrd"=unique(res$bkgrd),"n"=NA)
newCols <- data.frame(matrix(NA,nrow=nrow(subgrpRes),ncol=23*3))
colnames(newCols) <- paste0(paste0("chr",c(23,1:22)),rep(c(".NAM",".EUR",".AFR"),each=23))
subgrpRes <- cbind(subgrpRes,newCols)
colnames(subgrpRes)

namCols <- c(70,46:67)
eurCols <- c(71,2:23)
afrCols <- c(69,24:45)

namAt <- c(70,74)
eurAt <- c(71,73)
afrAt <- c(69,72)

colnames(res)[namCols]; colnames(res)[eurCols]; colnames(res)[afrCols]
nonPres <- subgrpRes
pooledres <- subgrpRes

nonParam_diff <- function(tmp){
  diff_means=tmp[,1]-tmp[,2]
  m=sum(diff_means>0)
  n=nrow(tmp)
  #ifelse(pbinom(m,n,0.5)<0.5,2*pbinom(m,n,0.5),2*(1-pbinom(m,n,0.5)+dbinom(m,n,0.5)))
  binom.test(m,n)$p.value
}

param_diff <- function(tmp){
  t.test(tmp[,1],tmp[,2],paired=TRUE)$p.value
}

pooled_diff <- function(tmp){
  t.test(tmp[,1],tmp[,2],paired=FALSE)$p.value
}

for(i in 1:nrow(subgrpRes)){
  pop <- subgrpRes$bkgrd[i]
  subgrpRes$n[i] <- sum(is.element(res$bkgrd,pop))
  
  #subgrpRes[i,3] <- pValues(CAnD(res[res$bkgrd==pop,namCols]))[1]
  subgrpRes[i,3:25] <- pValues(CAnD(res[res$bkgrd==pop,namCols]))
  subgrpRes[i,3] <- param_diff(res[res$bkgrd==pop,namAt])
  
  #subgrpRes[i,26] <- pValues(CAnD(res[res$bkgrd==pop,eurCols]))[1]
  subgrpRes[i,26:48] <- pValues(CAnD(res[res$bkgrd==pop,eurCols]))
  subgrpRes[i,26] <- param_diff(res[res$bkgrd==pop,eurAt])
  
  #subgrpRes[i,49] <- pValues(CAnD(res[res$bkgrd==pop,afrCols]))[1]
  subgrpRes[i,49:71] <- pValues(CAnD(res[res$bkgrd==pop,afrCols]))  
  subgrpRes[i,49] <- param_diff(res[res$bkgrd==pop,afrAt])
  
  #
  nonPres$n[i] <- sum(is.element(res$bkgrd,pop))
  
  #nonPres[i,3] <- pValues(nonParam_CAnD(res[res$bkgrd==pop,namCols]))[1]
  nonPres[i,3:25] <- pValues(nonParam_CAnD(res[res$bkgrd==pop,namCols]))
  nonPres[i,3] <- nonParam_diff(res[res$bkgrd==pop,namAt])
  
  #nonPres[i,26] <- pValues(nonParam_CAnD(res[res$bkgrd==pop,eurCols]))[1]
  nonPres[i,26:48] <- pValues(nonParam_CAnD(res[res$bkgrd==pop,eurCols]))
  nonPres[i,26] <- nonParam_diff(res[res$bkgrd==pop,eurAt])
  
  #nonPres[i,49] <- pValues(nonParam_CAnD(res[res$bkgrd==pop,afrCols]))[1]
  nonPres[i,49:71] <- pValues(nonParam_CAnD(res[res$bkgrd==pop,afrCols]))  
  nonPres[i,49] <- nonParam_diff(res[res$bkgrd==pop,afrAt])
  
  #
  pooledres$n[i] <- sum(is.element(res$bkgrd,pop))
  
  #pooledres[i,3] <- pValues(CAnD_pooled(res[res$bkgrd==pop,namCols]))[1]
  pooledres[i,3:25] <- pValues(CAnD_pooled(res[res$bkgrd==pop,namCols]))
  pooledres[i,3] <- pooled_diff(res[res$bkgrd==pop,namAt])
  
  #pooledres[i,26] <- pValues(CAnD_pooled(res[res$bkgrd==pop,eurCols]))[1]
  pooledres[i,26:48] <- pValues(CAnD_pooled(res[res$bkgrd==pop,eurCols]))
  pooledres[i,26] <- pooled_diff(res[res$bkgrd==pop,eurAt])
  
  #pooledres[i,49] <- pValues(CAnD_pooled(res[res$bkgrd==pop,afrCols]))[1]
  pooledres[i,49:71] <- pValues(CAnD_pooled(res[res$bkgrd==pop,afrCols]))  
  pooledres[i,49] <- pooled_diff(res[res$bkgrd==pop,afrAt])
}



# x vs autos isn't bonf corrected
# autos is just autos vs other autos, IS bonf corrected

cbind(subgrpRes[,c(3,26,49)],nonPres[,c(3,26,49)],pooledres[,c(3,26,49)])
cand_results <- list("CAnD_pvalue"=subgrpRes,"CAnD_nonP"=nonPres,"Pooled_t"=pooledres)
save(cand_results,file="cand_results_inclX_bonfCorr.RData")

# can run a sens analysis where the autosomes are simply looked at against the autosomes, ie excl the x
subgrpRes <- data.frame("bkgrd"=unique(res$bkgrd),"n"=NA)
newCols <- data.frame(matrix(NA,nrow=nrow(subgrpRes),ncol=22*3))
colnames(newCols) <- paste0(paste0("chr",c(1:22)),rep(c(".NAM",".EUR",".AFR"),each=22))
subgrpRes <- cbind(subgrpRes,newCols)
colnames(subgrpRes)

namCols <- c(46:67)
eurCols <- c(2:23)
afrCols <- c(24:45)

colnames(res)[namCols]; colnames(res)[eurCols]; colnames(res)[afrCols]
nonPresL <- subgrpRes
nonPresG <- subgrpRes
pooledres <- subgrpRes

for(i in 1:nrow(subgrpRes)){
  pop <- subgrpRes$bkgrd[i]
  subgrpRes$n[i] <- sum(is.element(res$bkgrd,pop))
  
  #subgrpRes[i,3] <- pValues(CAnD(res[res$bkgrd==pop,namCols]))[1]
  subgrpRes[i,3:24] <- pValues(CAnD(res[res$bkgrd==pop,namCols]))
  
  #subgrpRes[i,26] <- pValues(CAnD(res[res$bkgrd==pop,eurCols]))[1]
  subgrpRes[i,25:46] <- pValues(CAnD(res[res$bkgrd==pop,eurCols]))
  
  #subgrpRes[i,49] <- pValues(CAnD(res[res$bkgrd==pop,afrCols]))[1]
  subgrpRes[i,47:68] <- pValues(CAnD(res[res$bkgrd==pop,afrCols]))  
  #
  nonPresL$n[i] <- sum(is.element(res$bkgrd,pop))
  
  nonPresL[i,3:24] <- pValues(nonParam_CAnD(res[res$bkgrd==pop,namCols]))
  nonPresL[i,25:46] <- pValues(nonParam_CAnD(res[res$bkgrd==pop,eurCols]))
  if(pop=="Cuban"){nonPresL[i,25:46] <- pValues(nonParam_CAnD(res[res$bkgrd==pop,eurCols],alter="greater"))}
  nonPresL[i,47:68] <- pValues(nonParam_CAnD(res[res$bkgrd==pop,afrCols]))  
  #
  nonPresG$n[i] <- sum(is.element(res$bkgrd,pop))
  
  nonPresG[i,3:24] <- pValues(nonParam_CAnD(res[res$bkgrd==pop,namCols],alter="two.sided"))
  nonPresG[i,25:46] <- pValues(nonParam_CAnD(res[res$bkgrd==pop,eurCols],alter="two.sided"))
  nonPresG[i,47:68] <- pValues(nonParam_CAnD(res[res$bkgrd==pop,afrCols],alter="two.sided"))  
  #
  
  pooledres$n[i] <- sum(is.element(res$bkgrd,pop))
  
  #pooledres[i,3] <- pValues(CAnD_pooled(res[res$bkgrd==pop,namCols]))[1]
  pooledres[i,3:24] <- pValues(CAnD_pooled(res[res$bkgrd==pop,namCols]))
  
  #pooledres[i,26] <- pValues(CAnD_pooled(res[res$bkgrd==pop,eurCols]))[1]
  pooledres[i,25:46] <- pValues(CAnD_pooled(res[res$bkgrd==pop,eurCols]))
  
  #pooledres[i,49] <- pValues(CAnD_pooled(res[res$bkgrd==pop,afrCols]))[1]
  pooledres[i,47:68] <- pValues(CAnD_pooled(res[res$bkgrd==pop,afrCols]))  
}
# all of these are bonf corrected
# these are the autosomal results, excluding the x chr from the pool

auto_results <- list("CAnD_pvalue"=subgrpRes,"CAnD_nonP_less"=nonPresL,
                     "CAnD_nonP_greater"=nonPresG,"Pooled_t"=pooledres)
save(auto_results,file="cand_results_autoOnly_bonfCorr.RData")

rm(list=ls())


#####
# 56. Parse CAnD results

cand_results <- get(load("cand_results_inclX_bonfCorr.RData"))
# make a table of the x chr results
cand <- cand_results[["CAnD_pvalue"]]
nonP <- cand_results[["CAnD_nonP"]]
pooled <- cand_results[["Pooled_t"]]

library(xtable)
r1=cand[,c("bkgrd","n","chr23.NAM","chr23.EUR","chr23.AFR")]
r2=nonP[,c("bkgrd","n","chr23.NAM","chr23.EUR","chr23.AFR")]
r3=pooled[,c("bkgrd","n","chr23.NAM","chr23.EUR","chr23.AFR")]

r1$method="CAnD"
r2$method="nonParam_CAnD"
r3$method="pooled_ttest"

library(ggplot2); library(reshape); library(RColorBrewer)
res <- rbind(melt(r1,c("bkgrd","n","method")),melt(r2,c("bkgrd","n","method")),
             melt(r3,c("bkgrd","n","method")))
res$Ancestry[res$variable=="chr23.NAM"] <- "America"
res$Ancestry[res$variable=="chr23.EUR"] <- "Europe"
res$Ancestry[res$variable=="chr23.AFR"] <- "Africa"
res$bkgrd=paste0(res$bkgrd,"\nn=",res$n)
SET1 <- setNames(brewer.pal(9, "Set1"),
                 c("red", "blue", "green", "purple", "orange", "brown", "pink", "gray"))
  levels <- c("Cuban\nn=1732", "Dominican\nn=897", "PuertoRican\nn=1703", "Mexican\nn=3635", 
              "CentralAmerican\nn=1094", "SouthAmerican\nn=686", "Other/Unknown\nn=326")
  col <- setNames(SET1[c("green", "orange", "purple", "red", "blue","green", "gray")], levels)
  col["reference"] <- "gray80"

res$bkgrd <- ordered(res$bkgrd,levels=levels)

res$bg <- substr(res$bkgrd,start=1,stop=(regexpr("\n",res$bkgrd)-2))
res$bg[res$bg=="CentralAmerica"] <- "CAm"
res$bg[res$bg=="SouthAmerica"] <- "SAm"
res$bg[res$bg=="PuertoRica"] <- "PR"
res$bg[res$bg=="Other/Unknow"] <- "Unk"
res$bg[res$bg=="Dominica"] <- "Dom"
res$bg[res$bg=="Mexica"] <- "Mex"

res$bg <- ordered(res$bg,levels=c("Cuba","Dom","PR","Mex","CAm","SAm","Unk"))

pdf("pvalues_CAnD_etal_olga.pdf",width=11,height=9)
ggplot(res,aes(x=bg,y=-log10(value))) + geom_point(size=3,aes(color=factor(bkgrd)))+
  facet_grid(Ancestry~method)+theme_bw()+geom_hline(yintercept=-log10(0.05))+
  xlab("Self-Identified Background") + scale_color_manual(values=col,breaks=rev(names(col)),
                                                          guide = guide_legend(title = "")) +
  theme(legend.position = "bottom",axis.text=element_text(size=13),axis.title=element_text(size=18),title=element_text(size=20),
        strip.text = element_text(size=16),legend.text=element_text(size=14)) +  labs(y=expression(-log["10"]*"(pvalue)")) 
dev.off()

## now take results, plot pvalues for all chrs for each self-id background in turn
cand_results <- get(load("cand_results_autoOnly_bonfCorr.RData"))
# make a table of the x chr results
cand <- cand_results[["CAnD_pvalue"]]
nonP <- cand_results[["CAnD_nonP_less"]]
pooled <- cand_results[["Pooled_t"]]

r1=cand
r2=nonP
r3=pooled

r1$method="CAnD"
r2$method="nonParam_CAnD"
r3$method="pooled_ttest"

res <- rbind(melt(r1,c("bkgrd","n","method")),melt(r2,c("bkgrd","n","method")),
             melt(r3,c("bkgrd","n","method")))
res$Ancestry <- NA
res$Ancestry[grep("NAM",res$variable)] <- "America"
res$Ancestry[grep("EUR",res$variable)] <- "Europe"
res$Ancestry[grep("AFR",res$variable)] <- "Africa"
res$bkgrd=paste0(res$bkgrd,"\nn=",res$n)
SET1 <- setNames(brewer.pal(9, "Set1"),
                 c("red", "blue", "green", "purple", "orange", "yellow", "brown", "pink", "gray"))
levels <- c("Cuban\nn=1732", "Dominican\nn=897", "PuertoRican\nn=1703", "Mexican\nn=3635", 
            "CentralAmerican\nn=1094", "SouthAmerican\nn=686", "Other/Unknown\nn=326")
col <- setNames(SET1[c("green", "orange", "purple", "red", "yellow", "blue", "gray")], levels)
col["reference"] <- "gray80"

res$bkgrd <- factor(res$bkgrd)

res$bg <- substr(res$bkgrd,start=1,stop=(regexpr("\n",res$bkgrd)-2))
res$bg[res$bg=="CentralAmerica"] <- "CAm"
res$bg[res$bg=="SouthAmerica"] <- "SAm"
res$bg[res$bg=="PuertoRica"] <- "PR"
res$bg[res$bg=="Other/Unknow"] <- "Unk"
res$bg[res$bg=="Dominica"] <- "Dom"
res$bg[res$bg=="Mexica"] <- "Mex"

res$chr <- gsub("chr","",res$variable)
res$chr <- gsub(".NAM","",res$chr)
res$chr <- gsub(".EUR","",res$chr)
res$chr <- gsub(".AFR","",res$chr)

levels <- c(1:22)
col <- setNames(c(SET1[c("green", "orange", "purple", "red", "blue","green", "pink")],
                  SET1[c("green", "orange", "purple", "red",  "blue","green", "pink")],
                  SET1[c("green", "orange", "purple", "red",  "blue","green", "pink")],
                  SET1[c("green")]), levels)

resPl <- res[res$bg=="CAm",]
pdf("pvalues_CAnD_etal_CAm_autos.pdf",width=11,height=9)
ggplot(resPl,aes(x=as.integer(chr),y=-log10(value))) + geom_point(size=3,aes(color=factor(chr)))+
  facet_grid(Ancestry~method)+theme_bw()+geom_hline(yintercept=-log10(0.05))+
  xlab("Chromosome") + scale_color_manual(values=col,breaks=rev(names(col)),
                                          guide = guide_legend(title = "")) +
  theme(legend.position="none",axis.text=element_text(size=14),axis.title=element_text(size=18),title=element_text(size=20),
        strip.text = element_text(size=16)) +  labs(y=expression(-log["10"]*"(pvalue)")) +
  scale_x_continuous(breaks=c(1,5,10,15,20),labels=c(1,5,10,15,20)) +
  ggtitle("Self-Identified Central Americans, n=1094")
dev.off()

resPl <- res[res$bg=="SAm",]
pdf("pvalues_CAnD_etal_SAm_autos.pdf",width=11,height=9)
ggplot(resPl,aes(x=as.integer(chr),y=-log10(value))) + geom_point(size=3,aes(color=factor(chr)))+
  facet_grid(Ancestry~method)+theme_bw()+geom_hline(yintercept=-log10(0.05))+
  xlab("Chromosome") + scale_color_manual(values=col,breaks=rev(names(col)),
                                          guide = guide_legend(title = "")) +
  theme(legend.position="none",axis.text=element_text(size=14),axis.title=element_text(size=18),title=element_text(size=20),
        strip.text = element_text(size=16)) +  labs(y=expression(-log["10"]*"(pvalue)")) +
  scale_x_continuous(breaks=c(1,5,10,15,20),labels=c(1,5,10,15,20)) +
  ggtitle("Self-Identified South Americans, n=686")
dev.off()

resPl <- res[res$bg=="PR",]
pdf("pvalues_CAnD_etal_PR_autos.pdf",width=11,height=9)
ggplot(resPl,aes(x=as.integer(chr),y=-log10(value))) + geom_point(size=3,aes(color=factor(chr)))+
  facet_grid(Ancestry~method)+theme_bw()+geom_hline(yintercept=-log10(0.05))+
  xlab("Chromosome") + scale_color_manual(values=col,breaks=rev(names(col)),
                                          guide = guide_legend(title = "")) +
  theme(legend.position="none",axis.text=element_text(size=14),axis.title=element_text(size=18),title=element_text(size=20),
        strip.text = element_text(size=16)) +  labs(y=expression(-log["10"]*"(pvalue)")) +
  scale_x_continuous(breaks=c(1,5,10,15,20),labels=c(1,5,10,15,20)) +
  ggtitle("Self-Identified Puerto Ricans, n=1703")
dev.off()

resPl <- res[res$bg=="Mex",]
pdf("pvalues_CAnD_etal_Mex_autos.pdf",width=11,height=9)
ggplot(resPl,aes(x=as.integer(chr),y=-log10(value))) + geom_point(size=3,aes(color=factor(chr)))+
  facet_grid(Ancestry~method)+theme_bw()+geom_hline(yintercept=-log10(0.05))+
  xlab("Chromosome") + scale_color_manual(values=col,breaks=rev(names(col)),
                                          guide = guide_legend(title = "")) +
  theme(legend.position="none",axis.text=element_text(size=14),axis.title=element_text(size=18),title=element_text(size=20),
        strip.text = element_text(size=16)) +  labs(y=expression(-log["10"]*"(pvalue)")) +
  scale_x_continuous(breaks=c(1,5,10,15,20),labels=c(1,5,10,15,20)) +
  ggtitle("Self-Identified Mexicans, n=3635")
dev.off()

resPl <- res[res$bg=="Dom",]
pdf("pvalues_CAnD_etal_Dom_autos.pdf",width=11,height=9)
ggplot(resPl,aes(x=as.integer(chr),y=-log10(value))) + geom_point(size=3,aes(color=factor(chr)))+
  facet_grid(Ancestry~method)+theme_bw()+geom_hline(yintercept=-log10(0.05))+
  xlab("Chromosome") + scale_color_manual(values=col,breaks=rev(names(col)),
                                          guide = guide_legend(title = "")) +
  theme(legend.position="none",axis.text=element_text(size=14),axis.title=element_text(size=18),title=element_text(size=20),
        strip.text = element_text(size=16)) +  labs(y=expression(-log["10"]*"(pvalue)")) +
  scale_x_continuous(breaks=c(1,5,10,15,20),labels=c(1,5,10,15,20)) +
  ggtitle("Self-Identified Dominicans, n=897")
dev.off()

resPl <- res[res$bg=="Cuba",]
pdf("pvalues_CAnD_etal_Cuba_autos.pdf",width=11,height=9)
ggplot(resPl,aes(x=as.integer(chr),y=-log10(value))) + geom_point(size=3,aes(color=factor(chr)))+
  facet_grid(Ancestry~method)+theme_bw()+geom_hline(yintercept=-log10(0.05))+
  xlab("Chromosome") + scale_color_manual(values=col,breaks=rev(names(col)),
                                          guide = guide_legend(title = "")) +
  theme(legend.position="none",axis.text=element_text(size=14),axis.title=element_text(size=18),title=element_text(size=20),
        strip.text = element_text(size=16)) +  labs(y=expression(-log["10"]*"(pvalue)")) +
  scale_x_continuous(breaks=c(1,5,10,15,20),labels=c(1,5,10,15,20)) +
  ggtitle("Self-Identified Cubans, n=1732")
dev.off()

resPl <- res[res$bg=="Unk",]
pdf("pvalues_CAnD_etal_Unk_autos.pdf",width=11,height=9)
ggplot(resPl,aes(x=as.integer(chr),y=-log10(value))) + geom_point(size=3,aes(color=factor(chr)))+
  facet_grid(Ancestry~method)+theme_bw()+geom_hline(yintercept=-log10(0.05))+
  xlab("Chromosome") + scale_color_manual(values=col,breaks=rev(names(col)),
                                          guide = guide_legend(title = "")) +
  theme(legend.position="none",axis.text=element_text(size=14),axis.title=element_text(size=18),title=element_text(size=20),
        strip.text = element_text(size=16)) +  labs(y=expression(-log["10"]*"(pvalue)")) +
  scale_x_continuous(breaks=c(1,5,10,15,20),labels=c(1,5,10,15,20)) +
  ggtitle("Self-Identified Other/Unknown, n=326")
dev.off()

rm(list=ls())


#####
# 57. Look at PR chr 2 results in depth

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")
library(GWASTools)
library(gdsfmt)

# load in local ancestry estimates
gdsfn <- "/projects/geneva/gcc-fs2/OLGA/genotype/native_american/olga_with_refpanel/reich_1000G/local_ancestry/gds/unique/olga_reich_1000G_lai_new_unique.gds"
gdsobj <- openfn.gds(gdsfn)

# need to get the mean values across chromosomes
# figure out how to loop through and take means for each ancestry, by chr
chrs <- read.gdsn(index.gdsn(gdsobj, "snp.chromosome"))
length(chrs) # 14192
(t <- table(chrs))
#   1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
#1116 1072  939  874  846  811  760  688  661  720  665  709  520  479  484  514 
# 17   18   19   20   21   22 
#506  488  402  438  245  255 
stInd <- cumsum(t)

pos <- read.gdsn(index.gdsn(gdsobj,"snp.position"))
length(pos); head(pos) # 14192

# ok, so the ancestry estimates are matrices of sample x snp
# just read in chr 2 ancestry estimates, for all SNPs/regions

ancest <- data.frame(matrix(NA,nrow=12803,ncol=(1072*3+1)))
colnames(ancest) <- c("scanID",rep(c("eur","afr","amer"),each=1072))

ancest$scanID <- read.gdsn(index.gdsn(gdsobj,"sample.id"))

tmpEur <- read.gdsn(index.gdsn(gdsobj,"dosage_eur"),start=c(1,stInd[2]),count=c(-1,t[2]))
tmpAfr <- read.gdsn(index.gdsn(gdsobj,"dosage_afr"),start=c(1,stInd[2]),count=c(-1,t[2]))
tmpNam <- read.gdsn(index.gdsn(gdsobj,"dosage_amer"),start=c(1,stInd[2]),count=c(-1,t[2]))
# now have the three ancestry estimates across chr2

# just subset to pr samples
scan <- getobj("admix_xchr_aftermerge.RData")
unrel <- getobj("unrelated_pcair_deg4.RData")

scan <- scan[is.element(scan$scanID,unrel),]
dim(scan) # 10073 in scan; this is correct

tmpEur <- data.frame(tmpEur)
tmpEur$scanID <- ancest$scanID
tmpEur <- merge(tmpEur,scan[,c("scanID","bkgrd")],by="scanID",all.y=TRUE)
dim(tmpEur) # 10073 1074

tmpAfr <- data.frame(tmpAfr)
tmpAfr$scanID <- ancest$scanID
tmpAfr <- merge(tmpAfr,scan[,c("scanID","bkgrd")],by="scanID",all.y=TRUE)
dim(tmpAfr) # 10073 1074

tmpNam <- data.frame(tmpNam)
tmpNam$scanID <- ancest$scanID
tmpNam <- merge(tmpNam,scan[,c("scanID","bkgrd")],by="scanID",all.y=TRUE)
dim(tmpNam) # 10073 1074

# get mean ancestry at each position across all PR individs
meanEur <- data.frame(mn=apply(tmpEur[tmpEur$bkgrd=="PuertoRican",c(2:1073)],2,function(x){sum(x)/(2*length(x))}))
meanAfr <- data.frame(mn=apply(tmpAfr[tmpAfr$bkgrd=="PuertoRican",c(2:1073)],2,function(x){sum(x)/(2*length(x))}))
meanNam <- data.frame(mn=apply(tmpNam[tmpNam$bkgrd=="PuertoRican",c(2:1073)],2,function(x){sum(x)/(2*length(x))}))

meanEur$Ancestry <- "Europe"
meanAfr$Ancestry <- "Africa"
meanNam$Ancestry <- "America"

meanNam$Segment <- 1:nrow(meanNam)
meanEur$Segment <- 1:nrow(meanEur)
meanAfr$Segment <- 1:nrow(meanAfr)

res <- rbind(meanEur,meanAfr,meanNam)
pdf("localAncestry_chr2_PR.pdf",width=11)
ggplot(res,aes(x=Segment,y=mn,color=Ancestry)) + geom_line(size=1.5) + theme_bw() +
  ggtitle("Mean Local Ancestry Proportion on Chromosome 2\n1703 Puerto Ricans") +
  theme(axis.text=element_text(size=16),axis.title=element_text(size=16),
        legend.text=element_text(size=16),title=element_text(size=20))+
  ylab("Mean Proportion Local Ancestry")
dev.off()

# perform CAnD in a sliding window fashion, comparing mean local ancestry in a 50 segment window to 
# a mean of the remaining chromosome

# only do for NAM ancestry, since that was sig in the global test
# calculate the means first; it'll be 2 cols per segment
segMeans <- data.frame(matrix(NA,nrow=sum(tmpNam$bkgrd=="PuertoRican"),ncol=21*2+1))
colnames(segMeans) <- c("scanID",paste0(c("s","wo_s"),rep(1:21,each=2)))

# do the first segment by hand
smNam <- tmpNam[tmpNam$bkgrd=="PuertoRican",2:101]
lgNam <- tmpNam[tmpNam$bkgrd=="PuertoRican",102:(ncol(tmpNam)-1)]
segMeans$s1 <- rowSums(smNam)/(2*ncol(smNam))
segMeans$wo_s1 <- rowSums(lgNam)/(2*ncol(lgNam))

for(i in 2:21){
  cl <- paste0("s",i)
  cl2 <- paste0("wo_s",i)
  
  if(i==21){
    smNam <- tmpNam[tmpNam$bkgrd=="PuertoRican",(1+(i-1)*50):(ncol(tmpNam)-1)]
    lgNam <- tmpNam[tmpNam$bkgrd=="PuertoRican",2:((i-1)*50)]
  }else{
    smNam <- tmpNam[tmpNam$bkgrd=="PuertoRican",(1+(i-1)*50):((i-1)*50+100)]
    lgNam <- tmpNam[tmpNam$bkgrd=="PuertoRican",c(2:((i-1)*50),((i-1)*50+100+1):(ncol(tmpNam)-1))]
  }
  segMeans[,cl] <- rowSums(smNam)/(2*ncol(smNam))
  segMeans[,cl2] <- rowSums(lgNam)/(2*ncol(lgNam))
}

pvals <- rep(NA,21)
for(i in seq(from=2,to=ncol(segMeans),by=2)){
  pvals[(i/2)] <- t.test(segMeans[,i],segMeans[,(i+1)],paired=TRUE)$p.value
}

pvals <- pvals*(length(pvals))
pvals[pvals>=1] <- 1

res <- data.frame(pvalue=pvals,id=1:length(pvals))

pdf("cand_pvalues_chr2_NAM_pr_100segmentBlock.pdf")
ggplot(res,aes(x=id,y=-log10(pvalue))) + theme_bw() + geom_point() +xlab("100 Segment Block") + geom_hline(yintercept=-log10(0.05))+
  ylab(expression(paste(-log[10],"(pvalue)"))) + 
  ggtitle("CAnD Results for Native American Ancestry, Chr 2\n1703 Puerto Ricans")
dev.off()

# segments 2, 6 are significant
# look at these on the local ancestry plot
# seg2: 50:150
# seg6: 250:350

res <- rbind(meanNam)
res$nextMn <- c(res$mn[2:nrow(res)],res$mn[nrow(res)])

pdf("localAncestry_chr2_PR_withSigBlocks.pdf",width=11)
par(mfrow=c(2,1))

res$CAnD_Significant <- "No"
res$CAnD_Significant[(res$Segment>50&res$Segment<150)] <- "Block 2"

ggplot(res) + geom_segment(aes(x=Segment,y=mn,color=CAnD_Significant,
                               xend=Segment+1,yend=nextMn),size=1.5) + theme_bw() +
  ggtitle("Mean Native American Ancestry on Chromosome 2\n1703 Puerto Ricans") +
  ylab("Mean Proportion Local Ancestry") + geom_segment(x=50,y=mean(segMeans[,"s2"]),xend=150,size=1.1,
                                                        yend=mean(segMeans[,"s2"]),color="gray") +
  geom_segment(x=1,y=mean(segMeans[,"wo_s2"]),xend=49,yend=mean(segMeans[,"wo_s2"]),color="gray",size=1.1) + 
  geom_segment(x=151,y=mean(segMeans[,"wo_s2"]),xend=1072,yend=mean(segMeans[,"wo_s2"]),color="gray",size=1.1)

res$CAnD_Significant <- "No"
res$CAnD_Significant[(res$Segment>250&res$Segment<350)] <- "Block 6"

ggplot(res) + geom_segment(aes(x=Segment,y=mn,color=CAnD_Significant,
                               xend=Segment+1,yend=nextMn),size=1.5) + theme_bw() +
  ylab("Mean Proportion Local Ancestry") + 
  geom_segment(x=1,xend=249,y=mean(segMeans[,"wo_s6"]),yend=mean(segMeans[,"wo_s6"]),size=1.1,color="gray") +
  geom_segment(x=351,xend=1072,y=mean(segMeans[,"wo_s6"]),yend=mean(segMeans[,"wo_s6"]),size=1.1,color="gray") +
  geom_segment(x=250,xend=350,y=mean(segMeans[,"s6"]),yend=mean(segMeans[,"s6"]),size=1.1,color="gray")
dev.off()

res$CAnD_Significant <- "No"
res$CAnD_Significant[(res$Segment>50&res$Segment<150)] <- "Block 2"

p1 <- ggplot(res) + geom_segment(aes(x=Segment,y=mn,color=CAnD_Significant,
                               xend=Segment+1,yend=nextMn),size=1.5) + theme_bw() +
  ggtitle("Mean Native American Ancestry on Chromosome 2\n1703 Puerto Ricans") +
  ylab("Mean Proportion Local Ancestry") + geom_segment(x=50,y=mean(segMeans[,"s2"]),xend=150,size=1.1,
                                                        yend=mean(segMeans[,"s2"]),color="gray") +
  geom_segment(x=1,y=mean(segMeans[,"wo_s2"]),xend=49,yend=mean(segMeans[,"wo_s2"]),color="gray",size=1.1) + 
  geom_segment(x=151,y=mean(segMeans[,"wo_s2"]),xend=1072,yend=mean(segMeans[,"wo_s2"]),color="gray",size=1.1)

res$CAnD_Significant <- "No"
res$CAnD_Significant[(res$Segment>250&res$Segment<350)] <- "Block 6"

p2 <- ggplot(res) + geom_segment(aes(x=Segment,y=mn,color=CAnD_Significant,
                               xend=Segment+1,yend=nextMn),size=1.5) + theme_bw() +
  ylab("Mean Proportion Local Ancestry") + 
  geom_segment(x=1,xend=249,y=mean(segMeans[,"wo_s6"]),yend=mean(segMeans[,"wo_s6"]),size=1.1,color="gray") +
  geom_segment(x=351,xend=1072,y=mean(segMeans[,"wo_s6"]),yend=mean(segMeans[,"wo_s6"]),size=1.1,color="gray") +
  geom_segment(x=250,xend=350,y=mean(segMeans[,"s6"]),yend=mean(segMeans[,"s6"]),size=1.1,color="gray")


multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
pdf("localAncestry_chr2_PR_withSigBlocks.pdf",width=21,height=14)
multiplot(p1, p2, cols=1)
dev.off()

# what chromosome position do these 'segment' regions map to?
# pos is the vector of chr positions for each local ancestry reading
# want chr 2, segment 50-150 and 250-350
annot <- data.frame(cbind(chrs,pos))
res$pos <- annot$pos[annot$chrs==2]
range(res$pos[res$Segment>50&res$Segment<150]) # 6538054 19958794
# it's really the first half of the segment
range(res$pos[res$Segment>50&res$Segment<100]) # 6538054 11535738


range(res$pos[res$Segment>250&res$Segment<350]) # 40983374 64690052

res$CAnD_Significant <- "No"
res$CAnD_Significant[(res$Segment>50&res$Segment<100)] <- "Block 2"
ggplot(res) + geom_segment(aes(x=Segment,y=mn,color=CAnD_Significant,
                                     xend=Segment+1,yend=nextMn),size=1.5)

rm(list=ls())


#####
# 58. Look at PR chr 2 results in depth, EUR ancestry

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")
library(GWASTools)
library(gdsfmt)

# load in local ancestry estimates
gdsfn <- "/projects/geneva/gcc-fs2/OLGA/genotype/native_american/olga_with_refpanel/reich_1000G/local_ancestry/gds/unique/olga_reich_1000G_lai_new_unique.gds"
gdsobj <- openfn.gds(gdsfn)

# need to get the mean values across chromosomes
# figure out how to loop through and take means for each ancestry, by chr
chrs <- read.gdsn(index.gdsn(gdsobj, "snp.chromosome"))
length(chrs) # 14192
(t <- table(chrs))
#   1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
#1116 1072  939  874  846  811  760  688  661  720  665  709  520  479  484  514 
# 17   18   19   20   21   22 
#506  488  402  438  245  255 
stInd <- cumsum(t)

pos <- read.gdsn(index.gdsn(gdsobj,"snp.position"))
length(pos); head(pos) # 14192

# ok, so the ancestry estimates are matrices of sample x snp
# just read in chr 2 ancestry estimates, for all SNPs/regions

ancest <- data.frame(matrix(NA,nrow=12803,ncol=(1072*3+1)))
colnames(ancest) <- c("scanID",rep(c("eur","afr","amer"),each=1072))

ancest$scanID <- read.gdsn(index.gdsn(gdsobj,"sample.id"))

tmpEur <- read.gdsn(index.gdsn(gdsobj,"dosage_eur"),start=c(1,1),count=c(-1,-1))
tmpAfr <- read.gdsn(index.gdsn(gdsobj,"dosage_afr"),start=c(1,1),count=c(-1,-1))
tmpNam <- read.gdsn(index.gdsn(gdsobj,"dosage_amer"),start=c(1,1),count=c(-1,-1))
# now have the three ancestry estimates across chr2

# just subset to pr samples
scan <- getobj("admix_xchr_aftermerge.RData")
unrel <- getobj("unrelated_pcair_deg4.RData")

scan <- scan[is.element(scan$scanID,unrel),]
dim(scan) # 10073 in scan; this is correct

tmpEur <- data.frame(tmpEur)
tmpEur$scanID <- ancest$scanID
tmpEur <- merge(tmpEur,scan[,c("scanID","bkgrd")],by="scanID",all.y=TRUE)
dim(tmpEur) # 10073 14194

tmpAfr <- data.frame(tmpAfr)
tmpAfr$scanID <- ancest$scanID
tmpAfr <- merge(tmpAfr,scan[,c("scanID","bkgrd")],by="scanID",all.y=TRUE)
dim(tmpAfr) # 10073 14194

tmpNam <- data.frame(tmpNam)
tmpNam$scanID <- ancest$scanID
tmpNam <- merge(tmpNam,scan[,c("scanID","bkgrd")],by="scanID",all.y=TRUE)
dim(tmpNam) # 10073 14194

# get mean ancestry at each position across all PR individs
meanEur <- data.frame(mn=apply(tmpEur[tmpEur$bkgrd=="PuertoRican",c((1116+2):(1116+2+1073))],2,function(x){sum(x)/(2*length(x))}))
meanAfr <- data.frame(mn=apply(tmpAfr[tmpAfr$bkgrd=="PuertoRican",c((1116+2):(1116+2+1073))],2,function(x){sum(x)/(2*length(x))}))
meanNam <- data.frame(mn=apply(tmpNam[tmpNam$bkgrd=="PuertoRican",c((1116+2):(1116+2+1073))],2,function(x){sum(x)/(2*length(x))}))

meanEur$Ancestry <- "Europe"
meanAfr$Ancestry <- "Africa"
meanNam$Ancestry <- "America"

meanNam$Segment <- 1:nrow(meanNam)
meanEur$Segment <- 1:nrow(meanEur)
meanAfr$Segment <- 1:nrow(meanAfr)

# perform CAnD in a sliding window fashion, comparing mean local ancestry in a 50 segment window to 
# a mean of the remaining chromosome

# do for EUR ancestry
# calculate the means first; it'll be 2 cols per segment
segMeans <- data.frame(matrix(NA,nrow=sum(tmpNam$bkgrd=="PuertoRican"),ncol=21*2+1))
colnames(segMeans) <- c("scanID",paste0(c("s","wo_s"),rep(1:21,each=2)))

# do the first segment by hand
smeur <- tmpEur[tmpEur$bkgrd=="PuertoRican",(1116+2):(1116+101)]
lgeur <- tmpEur[tmpEur$bkgrd=="PuertoRican",c((2:1117),(1218:ncol(tmpEur)-1))]
segMeans$s1 <- rowSums(smeur)/(2*ncol(smeur))
segMeans$wo_s1 <- rowSums(lgeur)/(2*ncol(lgeur))

for(i in 2:21){
  cl <- paste0("s",i)
  cl2 <- paste0("wo_s",i)
  
  smeur <- tmpEur[tmpEur$bkgrd=="PuertoRican",(1116+1+(i-1)*50):(1116+(i-1)*50+100)]
  lgeur <- tmpEur[tmpEur$bkgrd=="PuertoRican",c(2:(1116+(i-1)*50),(1116+(i-1)*50+100+1):(ncol(tmpEur)-1))]

  segMeans[,cl] <- rowSums(smeur)/(2*ncol(smeur))
  segMeans[,cl2] <- rowSums(lgeur)/(2*ncol(lgeur))
}

pvals <- rep(NA,21)
for(i in seq(from=2,to=ncol(segMeans),by=2)){
  pvals[(i/2)] <- t.test(segMeans[,i],segMeans[,(i+1)],paired=TRUE)$p.value
}

pvals <- pvals*(length(pvals))
pvals[pvals>=1] <- 1

res <- data.frame(pvalue=pvals,id=1:length(pvals))

pdf("cand_pvalues_chr2_EUR_pr_100segmentBlock.pdf")
ggplot(res,aes(x=id,y=-log10(pvalue))) + theme_bw() + geom_point() +xlab("100 Segment Block") + geom_hline(yintercept=-log10(0.05))+
  ylab(expression(paste(-log[10],"(pvalue)"))) + theme(axis.text=element_text(size=16),axis.title=element_text(size=16)) +
  ggtitle("CAnD Results for European Ancestry, Chr 2 Window vs Autosomal Mean\n1703 Puerto Ricans")
dev.off()

# get mean across genome
meanEur <- data.frame(mn=apply(tmpEur[tmpEur$bkgrd=="PuertoRican",(2:14193)],2,function(x){sum(x)/(2*length(x))}))
res <- rbind(meanEur)
res$nextMn <- c(res$mn[2:nrow(res)],res$mn[nrow(res)])
res$Segment <- 1:nrow(res)
res$Chr <- chrs

pdf("localAncestry_chr2_EUR_withSigBlocks.pdf",width=11)
par(mfrow=c(2,1))

res$CAnD_Significant <- "No"
res$CAnD_Significant[res$Segment<(100+1116)&res$Segment>1116] <- "Yes"
res$CAnD_Significant[(res$Segment>(350+1116)&res$Segment<(500+1116))] <- "Yes"
res$CAnD_Significant[(res$Segment>(800+1116)&res$Segment<(1000+1116))] <- "Yes"

pdf("localAncestry_chr2_EUR_withSigBlocks.pdf",width=14)
ggplot(res) + geom_segment(aes(x=Segment,y=mn,color=CAnD_Significant,
                               xend=Segment+1,yend=nextMn),size=1.5) + theme_bw() +
  ggtitle("Mean Native American Ancestry Autosomal-Wide\n1703 Puerto Ricans") +
  ylab("Mean Proportion Local Ancestry")
dev.off()


rm(list=ls())


#####
# 59. Re-run ADMIXTURE on the X chromosome

## moved all old results to
# /projects/geneva/geneva_sata/caitlin/mlm_x/olga_application/admixture/old_admixture_results/

# first run in batch prune and subset gds file to pruned x chr SNPs
# cd /projects/geneva/geneva_sata/caitlin/mlm_x/olga_application/admixture
# qsub -q thornton.q -N xfiltgds batch_gds2ped_subset.sh
# which calls gds2ped_subset_v2.R

## ran admixture with three sets of snps:
# 1. no ld pruning, n=9694
# 2. ld pruning with usual thresh of 0.32, n=2233
# 3. ld pruning with less stringent thresh of 0.48, n=3494

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x/olga_application/admixture")
library(GWASTools)

res1 <- read.table("olga_ref_unrel_xpruned.3.Q",as.is=TRUE)
res2 <- read.table("olga_ref_unrel_x.3.Q",as.is=TRUE)
res3 <- read.table("olga_ref_unrel_xpruned_048.3.Q",as.is=TRUE)

pop <- read.table("olga_ref_unrel_xpruned.pop")
fam <- read.table("olga_ref_unrel_xpruned.fam")

colnames(res1) <- paste("est",1:3,"_pruned",sep="")
colnames(res2) <- paste("est",1:3,sep="")
colnames(res3) <- paste("est",1:3,"_pruned048",sep="")

res <- cbind(fam$V2,pop,res1,res2,res3)
colnames(res)[1:2] <- c("scanID","pop")

res[res$pop=="NAM",] # est3 is NAM
head(res[res$pop=="EUR",]) # est1 is EUR 
# est2 is AFR

colnames(res)[3:5] <- paste(c("EUR","AFR","NAM"),"_pruned",sep="")
colnames(res)[6:8] <- paste(c("EUR","AFR","NAM"),"_all",sep="")
colnames(res)[9:11] <- paste(c("EUR","AFR","NAM"),"_pruned048",sep="")

apply(res[,3:5],2,summary)
#        EUR_pruned AFR_pruned NAM_pruned
#Min.       0.00001    0.00001    0.00001
#1st Qu.    0.21910    0.01495    0.11280
#Median     0.41320    0.07324    0.34490
#Mean       0.43850    0.17180    0.38960
#3rd Qu.    0.64660    0.22880    0.63550
#Max.       1.00000    1.00000    1.00000
# write a table to use in svn at tower
towrite <- res[,c(1,2,3:5)]
colnames(towrite)[3:5] <- c("EUR","AFR","NAM")
write.table(towrite,file="olga_admixture_xchr_results.txt")

colnames(towrite)[2] <- "ref.pop"
colnames(towrite)[3:5] <- c("Europe","Africa","America")
write.csv(towrite,file="admixture_xchr_k3.csv",row.names=FALSE,quote=FALSE)

library(ggplot2); library(reshape)

set <- res[order(res[,3],res[,4],res[,5]),]
set <- set[set$pop=="-",c("AFR_pruned","NAM_pruned","EUR_pruned")]

set$id <- 1:nrow(set)
set <- melt(set,id="id")

pdf("xchr_barplot_admixture.pdf",width=14)
ggplot(set,aes(id,value,fill=variable)) + geom_bar(stat="identity") +
  ylab("X Chromosome Ancestry Proportion") + xlab("Sample") + 
  ggtitle("ADMIXTURE Estimated X Chr Ancestry Proportions") + theme_bw()
dev.off()

# merge in subpop
scan <- getobj("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_scanAnnot_v54_AMS.RData")

res <- merge(res,pData(scan)[,c("scanID","gengrp6")],by="scanID",all.x=TRUE)
res$gengrp6[is.na(res$gengrp6)] <- "Other"
set <- res[order(res[,3],res[,4],res[,5]),]
set <- set[set$pop=="-",c("AFR_pruned","NAM_pruned","EUR_pruned","gengrp6")]

set$id <- 1:nrow(set)
set <- melt(set,id=c("id","gengrp6"))

pdf("xchr_barplot_admixture_byGengrp6.pdf",width=14)
p =ggplot(set,aes(x=id,y=value,fill=variable)) 
p  + geom_bar(stat="identity",position="stack") + facet_grid(.~gengrp6,scales="free",space="free")+
  ylab("X Chromosome Ancestry Proportion") + xlab("Sample") + 
  ggtitle("ADMIXTURE Estimated X Chr Ancestry Proportions") + theme_bw()
dev.off()

rm(list=ls())


#####
# 60. CAnD on the x vs auto proportions

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x/olga_application/admixture")
library(GWASTools)

res <- read.table("olga_admixture_xchr_results.txt",header=TRUE,as.is=TRUE)
colnames(res)[3:5] <- c("EUR.x","AFR.x","NAM.x")

autos <- read.csv("genetic_diversity/plots/data/admixture_k3.csv",header=TRUE,as.is=TRUE)

scan <- getobj("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_scanAnnot_v54_AMS.RData")
res <- merge(res,pData(scan)[,c("scanID","gengrp6","unrelated.pcair.deg4")],by="scanID",all.x=TRUE)
res$gengrp6[is.na(res$gengrp6)] <- "Other"


# compare mean autosomal proportions to x chr
# paired t-test, stratified self-id bkgrd group

cand_res <- data.frame(matrix(NA,nrow=7,ncol=7))

bkgs <- unique(res$gengrp6)
colnames(cand_res) <- bkgs
rownames(cand_res) <- c("n","chrx.NAM","chrx.EUR","chrx.AFR","pooled.NAM","pooled.EUR",
                        "pooled.AFR")

for(i in bkgs){
  thisbkg <- res$scanID[is.element(res$gengrp6,i)&res$unrelated.pcair.deg4]
  cand_res["n",colnames(cand_res)==i] <- length(thisbkg)
  cand_res["chrx.EUR",colnames(cand_res)==i] <- t.test(res[is.element(res$scanID,thisbkg),"EUR.x"],
       autos[is.element(autos$scanID,thisbkg),"Europe"],paired=TRUE)$p.value
  cand_res["pooled.EUR",colnames(cand_res)==i] <- t.test(res[is.element(res$scanID,thisbkg),"EUR.x"],
                                                       autos[is.element(autos$scanID,thisbkg),"Europe"])$p.value
  
  cand_res["chrx.NAM",colnames(cand_res)==i] <- t.test(res[is.element(res$scanID,thisbkg),"NAM.x"],
                                                     autos[is.element(autos$scanID,thisbkg),"America"],paired=TRUE)$p.value
  cand_res["pooled.NAM",colnames(cand_res)==i] <- t.test(res[is.element(res$scanID,thisbkg),"NAM.x"],
                                                       autos[is.element(autos$scanID,thisbkg),"America"])$p.value
  
  cand_res["chrx.AFR",colnames(cand_res)==i] <- t.test(res[is.element(res$scanID,thisbkg),"AFR.x"],
                                                     autos[is.element(autos$scanID,thisbkg),"Africa"],paired=TRUE)$p.value
  cand_res["pooled.AFR",colnames(cand_res)==i] <- t.test(res[is.element(res$scanID,thisbkg),"AFR.x"],
                                                       autos[is.element(autos$scanID,thisbkg),"Africa"])$p.value
}

write.csv(cand_res,file="CAnD_results_xvsAuto.csv",quote=FALSE)

rm(list=ls())


#####
# 61. Get k0, k1, k2 for X KC F-F pairs

library(GWASTools)
library(SNPRelate)
setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")

source("pcrelate.R")
# call pcrelate with no x chr option to get k0, k1, k2 calculation
# use unrel.set to be the same unrelated set used for PCA on the X chr
# use scan.include to be the same set of individs used for PCA on the X chr

scan <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_scanAnnot_v54_AMS.RData"))
head(pData(scan)) 

# make genoData object
gdsobj <- GdsGenotypeReader("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/netCDF/subjects/OLGA_subj_geno.gds")
genoData <- GenotypeData(gdsobj,scanAnnot=scan[scan$subj.plink,])

unrel.set <- scan$scanID[scan$unrelated.pcair.deg4&is.na(scan$gengrp6.outliers)]
scanIncl <- scan$scanID[is.na(scan$gengrp6.outliers)&scan$geno.cntl==0&scan$subj.plink & scan$sex=="F"]

head(unrel.set); head(scanIncl)
length(unrel.set); length(scanIncl) # 10272 | 7515

### need to exclude 13 people with an entirely missing x chr 
# from chromosome.anomalies.csv file, they are study samples with x chr anomaly
# thus, their entire x chr is filtered out. 

chromAnoms <- read.csv("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/dbGaP/To_dbGaP/chromosome_anomalies/chromosome_anomalies.csv",
                       header=TRUE,as.is=T)
ids <- chromAnoms$scanID[chromAnoms$chromosome=="X"&chromAnoms$filter&chromAnoms$whole.chr&
                           !is.element(substr(chromAnoms$SUBJECT_ID,1,2),"NA")]
length(ids) # 13
chromAnoms[is.element(chromAnoms$scanID,ids),] # perfect! exactly what we want

length(scanIncl) # 7515
scanIncl <- scanIncl[!is.element(scanIncl,ids)]
length(scanIncl) # 7513

snp.pruned <- get(load("olga_application/snp_sel_xChr_ldPruned.RData"))
head(snp.pruned); length(snp.pruned) # 3600

xPC <- get(load("olga_application/pca_prunedXsnps_10272unrel.RData"))

pcMat <- xPC$vectors[,c(1:2)]
dim(pcMat) # 12747 2

pcMat <- pcMat[is.element(rownames(pcMat),scanIncl),]
dim(pcMat) # 7513 2

rel <- pcrelate(genoData, unrel.set=unrel.set, scan.include=scanIncl, Xchr=TRUE, pcMat=pcMat, snp.include=snp.pruned)
names(rel)

save(rel,file="olga_application/pcRelate_notpcRelateX_prunedXchr_xPC12adj.RData")

rm(list=ls())

##
pcRes <- getobj("olga_application/pcRelate_notpcRelateX_prunedXchr_xPC12adj.RData")
head(pcRes$kinship) 

# plot k0 vs kin for kin>2^(-13/2)
png("k2_kc_FFpairs.png")
plot(pcRes$kinship$k2[pcRes$kinship$kin>2^(-13/2)],pcRes$kinship$kin[pcRes$kinship$kin>2^(-13/2)],
     xlab="k2",ylab="kinship coef",main="k2 vs KC on X Chr for F-F pairs")
dev.off()
# hmm.

png("k0_k1_FFpairs.png")
plot(pcRes$kinship$k0[pcRes$kinship$kin>2^(-13/2)],pcRes$kinship$k1[pcRes$kinship$kin>2^(-13/2)],
     xlab="k0",ylab="k1",main="k0 vs k1 on X Chr for F-F pairs")
dev.off()

# plot the expected values for each of the relationships
# mom-daughter: kc=0.25, k2=0, k1=0.25, k0=0.75
# full sisters: kc=3/8, k2=0.5 -- but could range anywhere from 0-1
# maternal aunt-niece: kc=3/16, k2=0
# maternal gma-gdaughter: kc=1/8, k2=0
# paternal aunt-niece: kc=1/8, k2=0
# paternal gma-gdaughter: kc=0.25, k2=0

rm(list=ls())


#####
# 62. Make some plots of KC on X chr for MLM-X paper

library(GWASTools)
library(SNPRelate)
setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")

xpc <- getobj("olga_application/pcRelate_prunedXchr_xPC12adj.RData")
xkc <- xpc$kinship
xkc$ID12 <- paste(xkc$ID1,xkc$ID2)

# load no Asians PC-Relate results
pcres <- getobj("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/mconomos/results/relatedness/v13_without_asian2/pcrelate_allChr.RData")
dim(pcres$kinship) # 81708936        7

pcres <- pcres$kinship

xsamp <- unique(xkc$ID1)
length(xsamp) # 12733

pcresSm <- pcres[is.element(pcres$ID1,xsamp),]
pcresSm <- pcresSm[is.element(pcresSm$ID2,xsamp),]

length(unique(pcresSm$ID1)) # 12732
all(is.element(unique(pcresSm$ID1),xsamp)) # TRUE
all(is.element(xsamp,unique(pcresSm$ID1))) # FALSE
# dang, so we're off by one sample...

xsamp[!is.element(xsamp,pcresSm$ID1)] # 99933

# add in relationship
# build up a new data frame with relationships included

# number of relatives by each degree
pcPO <- pcresSm$kin > 2^(-5/2) & pcresSm$k0 < 0.025; table(pcPO) # 1442 T
pcFS <- pcresSm$kin > 2^(-5/2) & pcresSm$k0 >= 0.025 & pcresSm$k2 > 2^(-7/2); table(pcFS) # 695 T
pcdeg2 <- pcresSm$kin < 2^(-5/2) & pcresSm$kin > 2^(-7/2) & pcresSm$k2 < 2^(-9/2); table(pcdeg2) # 580 T
pcdeg3 <- pcresSm$kin < 2^(-7/2) & pcresSm$kin > 2^(-9/2); table(pcdeg3) # 484 T
pcdeg4 <- pcresSm$kin < 2^(-9/2) & pcresSm$kin > 2^(-11/2); table(pcdeg4) # 940 T
pcUnrel <- pcresSm$kin < 2^(-11/2); table(pcUnrel)

po <- pcresSm[pcPO,]
po$ID12 <- paste(po$ID1,po$ID2)
po <- merge(po,xkc,by="ID12",suffixes=c(".auto",".x"),all.x=TRUE)
po$relationship <- "PO"

fs <- pcresSm[pcFS,]
fs$ID12 <- paste(fs$ID1,fs$ID2)
fs <- merge(fs,xkc,by="ID12",suffixes=c(".auto",".x"),all.x=TRUE)
fs$relationship <- "FS"

d2 <- pcresSm[pcdeg2,]
d2$ID12 <- paste(d2$ID1,d2$ID2)
d2 <- merge(d2,xkc,by="ID12",suffixes=c(".auto",".x"),all.x=TRUE)
d2$relationship <- "Deg2"

d3 <- pcresSm[pcdeg3,]
d3$ID12 <- paste(d3$ID1,d3$ID2)
d3 <- merge(d3,xkc,by="ID12",suffixes=c(".auto",".x"),all.x=TRUE)
d3$relationship <- "Deg3"

d4 <- pcresSm[pcdeg4,]
d4$ID12 <- paste(d4$ID1,d4$ID2)
d4 <- merge(d4,xkc,by="ID12",suffixes=c(".auto",".x"),all.x=TRUE)
d4$relationship <- "Deg4"

# merge in sex info
scan <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_scanAnnot_v54_AMS.RData"))
head(pData(scan)) 

toPl <- rbind(po,fs,d2,d3,d4)

toPl <- merge(toPl,pData(scan)[,c("scanID","sex")],all.x=TRUE,by.x="ID1.auto",by.y="scanID")
toPl <- merge(toPl,pData(scan)[,c("scanID","sex")],all.x=TRUE,by.x="ID2.auto",by.y="scanID",suffixes=c(".1",".2"))

toPl$sex <- paste(toPl$sex.1,toPl$sex.2,sep="-")
toPl$sex[toPl$sex=="M-F"] <- "F-M"

save(toPl,file="kc_res_pairsAuto.RData")

library(ggplot2)

pdf("kc_allRelationships_xkc_violin.pdf")
toPl$relationship <- factor(toPl$relationship,levels=c("PO","FS","Deg2","Deg3","Deg4"))
ggplot(toPl,aes(x=relationship,y=kin.x,col=sex)) + geom_violin() + theme_bw() + xlab("Relationship") + ylab("X Chr KC")
dev.off()

pdf("kc_allRelationships_xkc.pdf")
toPl$relationship <- factor(toPl$relationship,levels=c("PO","FS","Deg2","Deg3","Deg4"))
ggplot(toPl,aes(x=relationship,y=kin.x,col=sex)) + geom_boxplot() + theme_bw() + xlab("Relationship") + ylab("X Chr KC")
dev.off()

# look at po, fs in detail
pdf("kc_pos_xkc.pdf")
toPlSm <- toPl[toPl$relationship=="PO",]
ggplot(toPlSm,aes(y=kin.auto,x=kin.x,col=sex)) + geom_point() + theme_bw() + xlab("X Chr KC") + ylab("Auto KC")
dev.off()

pdf("kc_fss_xkc.pdf")
toPlSm <- toPl[toPl$relationship=="FS",]
ggplot(toPlSm,aes(y=kin.auto,x=kin.x,col=sex)) + geom_point() + theme_bw() + xlab("X Chr KC") + ylab("Auto KC")
dev.off()

rm(list=ls())


# pairs identified as relatives
idx <- which(pcres$kinship$kin > 2^(-11/2)); length(idx) #4240

summary(pcres$kinship$kin[-idx])
#      Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
#-1.690e-02 -1.385e-03 -8.017e-06 -7.000e-09  1.375e-03  2.208e-02 

brks <- seq(from = min(pcres$kinship$kin[-idx]), to = max(pcres$kinship$kin[-idx]), length=11)
unrel.idx <- NULL
for(i in 1:10){
  print(i)
  tmp <- which(pcres$kinship$kin >= brks[i] & pcres$kinship$kin <= brks[(i+1)]); print(length(tmp))
  if(length(tmp) > 100){
    unrel.idx <- c(unrel.idx, sample(tmp, size=100))
  }else{
    unrel.idx <- c(unrel.idx, tmp)
  }
}
length(unrel.idx) # 908
summary(pcres$kinship$kin[unrel.idx])
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#-0.016900 -0.003177  0.003296  0.004419  0.012670  0.022080 


# assign colors for plotting based on my PC relatedness estimates
clr <- rep("#65171775",length(idx))
clr[pcFS[idx]] <- "#0000ff50"
clr[pcPO[idx]] <- "#ff000050"
clr[pcdeg2[idx]] <- "#ff00ff50"
clr[pcdeg3[idx]] <- "#00ff0050"
clr[pcdeg4[idx]] <- "#ff650050"
# add in unrelateds
clr <- c(clr, rep("#00000030", length(unrel.idx)))

pdf("k0_v_kin_new.pdf")
plot(pcres$kinship$k0[c(idx, unrel.idx)], pcres$kinship$kin[c(idx, unrel.idx)], pch=20, col=clr, xlab=expression(paste("PC-Relate ",k[0]," Estimate")), ylab="PC-Relate Kinship Estimate")
segments(0,-1,0,0.25, lty=2, col=2); segments(-1,0.25,0,0.25, lty=2, col=2)
segments(0.25,-1,0.25,0.25, lty=2, col=4); segments(0,0.25,0.25,0.25, lty=2, col=4)
segments(0.5,-1,0.5,0.125, lty=2, col=6); segments(-1,0.125,0.5,0.125, lty=2, col=6)
segments(0.75,-1,0.75,0.0625, lty=2, col=3); segments(-1,0.0625,0.75,0.0625, lty=2, col=3)
segments(0.875,-1,0.875,0.03125, lty=2, col="orange"); segments(-1,0.03125,0.875,0.03125, lty=2, col="orange")
segments(1,-1,1,0, lty=2); segments(-1,0,1,0, lty=2)
abline(0.25,-0.25,col="gray");
dev.off()

plot(pcres$kinship$k0[c(idx,unrel.idx)], pcres$kinship$k1[c(idx,unrel.idx)], pch=20, col=clr, xlab="PC k0", ylab="PC Kinship")
plot(pcres$kinship$k0[c(idx,unrel.idx)], pcres$kinship$k2[c(idx,unrel.idx)], pch=20, col=clr, xlab="PC k0", ylab="PC Kinship")

rm(list=ls())


#####
# 63. Calculate genomic inf across autos only, compare pvalues

library(QCpipeline)
library(GWASTools)

setwd("/projects/geneva/geneva_sata/caitlin/olga_xchr_assoc/Assoc/")

# read in xchr results, the autosomal model
(mac.thresh <- 30)
totalRes <- data.frame(matrix(NA,nrow=2128491,ncol=4))
ct <- 1
n <- 12489
maf.eff <- format(quadSolveMAF(mac.thresh, n), digits=3)

for(i in 1:22){
  assc <- getobj(paste("assoc_auto_316987_chr",i,".RData",sep=""))
  
  filt.imp <- which(assc$type %in% 0)
  assc$composite.filter <- TRUE
  assc$maf.filt <- assc$effN > mac.thresh
  filt.obs <- which(assc$type %in% c(2,3))
  
  filt.obs.maf <- intersect(filt.obs, which(assc$maf.filt))
  # this is the filter i want
  
  assc <- assc[filt.obs.maf,]
  
  totalRes[ct:(ct+nrow(assc)-1),] <- assc[,c("snpID","chromosome","Stat","pval")]
  ct <- ct+nrow(assc)
}

totalRes <- totalRes[!is.na(totalRes$X1),]
stat = totalRes[,"X3"]
stat <- pmax(stat, 0)

(lambda_simpleMLM <- calculateLambda(stat, df=1)) # 1.04951 

##
# do the same exercise for the other assoc test results
mlmX <- data.frame(matrix(NA,nrow=2128491,ncol=4))
ct <- 1

maxN <- rep(NA,22)
for(i in 1:22){
  assc <- getobj(paste("assoc_316987_chr",i,".RData",sep=""))
  filt.imp <- which(assc$type %in% 0)
  assc$composite.filter <- TRUE
  assc$maf.filt <- assc$effN > mac.thresh
  filt.obs <- which(assc$type %in% c(2,3))
  
  filt.obs.maf <- intersect(filt.obs, which(assc$maf.filt))
  # this is the filter i want
  
  assc <- assc[filt.obs.maf,]
  
  mlmX[ct:(ct+nrow(assc)-1),] <- assc[,c("snpID","chromosome","Stat","pval")]
  ct <- ct+nrow(assc)
}

mlmX <- mlmX[!is.na(mlmX$X1),]
stat = mlmX[,"X3"]
stat <- pmax(stat, 0)

(lambda_MLMx <- calculateLambda(stat, df=1)) # 1.04839

mlmX$pval_simple <- totalRes[,"X4"]
colnames(mlmX)[4] <- "pval"

# make a plot of each pval comparing mlm-x to simple mlm-x, across the autosomes
library(ggplot2)
png("../Plots/autoPvals_mlmX_vsAutoModel.png")
ggplot(mlmX,aes(-log10(pval),-log10(pval_simple))) + geom_point(aes(alpha=0.5)) + 
  labs(y=expression(-log[10]*"(Simple MLM p-value)"), x=expression(-log[10]*"(MLM-X p-value)")) + 
  geom_abline(intercept=0,slope=1) + theme_bw() + theme(legend.position="none")
dev.off()

cor(-log10(mlmX$pval),-log10(mlmX$pval_simple)) # 0.9972682

rm(list=ls())


#####
# 64. Rerun PC-AiR, PC-Relate iterations w 3rd deg relatives

library(GWASTools)
library(OLGApipeline)
library(QCpipeline)
library(MASS)
library(GENESIS)

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")

scan <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_scanAnnot_v54_AMS.RData"))
head(pData(scan)) 

snp.pruned <- get(load("olga_application/snp_sel_xChr_ldPruned.RData"))
head(snp.pruned); length(snp.pruned) # 3600

# make genoData object
gdsobj <- GdsGenotypeReader("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/netCDF/subjects/OLGA_subj_geno.gds")
genoData <- GenotypeData(gdsobj,scanAnnot=scan[scan$subj.plink,])

unrel.set <- scan$scanID[scan$unrelated.pcair&!is.element(scan$gengrp6.outliers, "AsianOutliers")&
                           scan$geno.cntl==0&scan$subj.plink]
scanIncl <- scan$scanID[!is.element(scan$gengrp6.outliers, "AsianOutliers")&scan$subj.plink&scan$geno.cntl==0]

head(unrel.set); head(scanIncl)
length(unrel.set); length(scanIncl) # 10624 | 12784

pc <- pcair(genoData,unrel.set=unrel.set,Xchr=TRUE,snp.include=snp.pruned,scan.include=scanIncl)
save(pc,file="olga_application/pca_prunedXsnps_10272unrel12784rel_unrelDeg3.RData")

# make a par coords plot to check still 1-2 to include

pc_colors <- rep("black",length(scanIncl))
scan <- scan[is.element(scan$scanID,scanIncl),]
pc_colors[scan$race.cat=="Mexican"] <- "goldenrod"
pc_colors[scan$race.cat=="CentralAmerican"] <- "red"
pc_colors[scan$race.cat=="SouthAmerican"] <- "magenta"
pc_colors[scan$race.cat=="PuertoRican"] <- "green"
pc_colors[scan$race.cat=="Cuban"] <- "blue"
pc_colors[scan$race.cat=="Dominican"] <- "burlywood"

pdf("tmp_parCords.pdf",width=14)
#pc$vectors[,1] <- pc$vectors[,1]*-1
parcoord(pc$vectors,col=pc_colors)
dev.off()

## run pc-relate adj for the pc's from above
xPC <- get(load("olga_application/pca_prunedXsnps_10272unrel12784rel_unrelDeg3.RData"))

pcMat <- xPC$vectors[,c(1:2)]
dim(pcMat) # 12784 2

png("xpc_12_3rdDeg.png")
plot(-pcMat[,1],pcMat[,2],xlab="EV 1", ylab="EV 2",col=pc_colors,pch=19)
legend("topleft",c("CentralAmerican","Cuban","Dominican","Mexican","PuertoRican","SouthAmerican","Other/Unknown"),pch=20,
       col=c("red","blue","burlywood","goldenrod","green","magenta","black"),cex=0.8)
dev.off()

png("xpc_12.png")
plot(xPC)
dev.off()

png("xpc_34.png")
plot(xPC,3,4)
dev.off()


oldPCs <- get(load("olga_application/pca_prunedXsnps_10272unrel.RData"))
old12 <- oldPCs$vectors[,c(1:2)]

png("xpc_12_3rdDegvs4th_ev1.png")
plot(pcMat[is.element(rownames(pcMat),rownames(old12)),1],old12[,1],xlab="4th Deg",ylab="3rd Deg")
abline(0,1)
dev.off()

png("xpc_12_3rdDegvs4th_ev2.png")
plot(pcMat[is.element(rownames(pcMat),rownames(old12)),2],old12[,2],xlab="4th Deg",ylab="3rd Deg")
abline(0,1)
dev.off()

scan <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_scanAnnot_v54_AMS.RData"))
head(pData(scan)) 

unrel.set <- scan$scanID[scan$unrelated.pcair&!is.element(scan$gengrp6.outliers, "AsianOutliers")&
                           scan$geno.cntl==0&scan$subj.plink]
scanIncl <- scan$scanID[!is.element(scan$gengrp6.outliers, "AsianOutliers")&scan$subj.plink&scan$geno.cntl==0]

head(unrel.set); head(scanIncl)
length(unrel.set); length(scanIncl) # 10624 | 12784

### need to exclude 13 people with an entirely missing x chr 
# from chromosome.anomalies.csv file, they are study samples with x chr anomaly
# thus, their entire x chr is filtered out. 

chromAnoms <- read.csv("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/dbGaP/To_dbGaP/chromosome_anomalies/chromosome_anomalies.csv",
                       header=TRUE,as.is=T)
ids <- chromAnoms$scanID[chromAnoms$chromosome=="X"&chromAnoms$filter&chromAnoms$whole.chr&
                           !is.element(substr(chromAnoms$SUBJECT_ID,1,2),"NA")]
length(ids) # 13
chromAnoms[is.element(chromAnoms$scanID,ids),] # perfect! exactly what we want

length(scanIncl) # 12784
scanIncl <- scanIncl[!is.element(scanIncl,ids)]
length(scanIncl) # 12771

unrel.set <- unrel.set[!is.element(unrel.set,ids)]
length(unrel.set) # 10614

# take these individs out of pcMat too
pcMat <- pcMat[is.element(rownames(pcMat),scanIncl),]
dim(pcMat) # 12771 2

rel <- pcrelate(genoData, unrel.set=unrel.set, scan.include=scanIncl, Xchr=TRUE, pcMat=pcMat,
                snp.include=snp.pruned)
names(rel)

save(rel,file="olga_application/pcRelate_Xchr_xPC12adj_unrel3Deg.RData")

## make into matrix
tmp <- matrix(NA,nrow=length(scanIncl),ncol=length(scanIncl))
length(tmp[lower.tri(tmp,diag=FALSE)]) # 81071011

kc <- rel$kinship
length(kc[,"kin"]) # 81071011
tmp[lower.tri(tmp,diag=FALSE)] <- kc[,"kin"]

# mirror to upper triangle
tmp[upper.tri(tmp,diag=FALSE)] <- t(tmp)[upper.tri(tmp,diag=FALSE)]

rownames(tmp) <- scanIncl
colnames(tmp) <- scanIncl

diag(tmp) <- 1 + rel$inbreed$f

save(tmp,file="olga_application/kcMat_xchr_xPC12adj_unrel3Deg.RData")

rm(list=ls())


#####
# 65. Manh of MLM-X, simple MLM X chr results in one png file

library(QCpipeline); library(OLGApipeline)
library(GWASTools)

setwd("/projects/geneva/geneva_sata/caitlin/olga_xchr_assoc/Assoc/")

# read in xchr results, the autosomal model
(mac.thresh <- 30)
n <- 12490
maf.eff <- format(quadSolveMAF(mac.thresh, n), digits=3)

assc <- getobj("assoc_auto_316987_chr23.RData")
  
  filt.imp <- which(assc$type %in% 0)
  assc$composite.filter <- TRUE
  assc$maf.filt <- assc$effN > mac.thresh
  filt.obs <- which(assc$type %in% c(2,3))
  
  filt.obs.maf <- intersect(filt.obs, which(assc$maf.filt))
  # this is the filter i want
  
  assc <- assc[filt.obs.maf,]
totalRes <- assc
  
##
# do the same exercise for the other assoc test results
assc <- getobj("assoc_316987_chr23.RData")

  filt.imp <- which(assc$type %in% 0)
  assc$composite.filter <- TRUE
  assc$maf.filt <- assc$effN > mac.thresh
  filt.obs <- which(assc$type %in% c(2,3))
  
  filt.obs.maf <- intersect(filt.obs, which(assc$maf.filt))
  # this is the filter i want
  
  assc <- assc[filt.obs.maf,]
  
mlmX <- assc

# make a plot of each pval comparing mlm-x to simple mlm-x, across the autosomes
library(ggplot2)
library(gridExtra); library(RGraphics)

png("../Plots/autoPvals_mlmX_vsAutoModel.png",width=1080, height=760)
mtxt <- ""
par(mfrow=c(2,1),mar=c(5,5,4,3)+0.1, lwd=1.5, cex.lab=1.5, cex.main=1.5, cex.sub=1.5, oma=c(0,0,length(mtxt)+0.1,0))
manhattanPlot(totalRes$pval, totalRes$chromosome)
mtext("A", side=3, line=0.75,adj=0,cex=2)
manhattanPlot(mlmX$pval, mlmX$chromosome)
mtext("B", side=3, line=0.75,adj=0,cex=2)
dev.off()

rm(list=ls())


#####
# 66. Run pcairPartition() incl & excl the autosomal KING divMat

# decide if we should include the divMat in the next PC-AiR run

library(GWASTools)
library(OLGApipeline)
library(QCpipeline)
library(MASS)
library(GENESIS)
library(ggplot2)

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x")

scan <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_scanAnnot_v54_AMS.RData"))
head(pData(scan)) 

unrel.set <- scan$scanID[scan$unrelated.pcair&!is.element(scan$gengrp6.outliers, "AsianOutliers")&
                           scan$geno.cntl==0&scan$subj.plink]
scanIncl <- scan$scanID[!is.element(scan$gengrp6.outliers, "AsianOutliers")&scan$subj.plink&scan$geno.cntl==0]

head(unrel.set); head(scanIncl)
length(unrel.set); length(scanIncl) # 10624 | 12784

# load king estimates
kingibd <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/amstilp/results/ibd/v04_build37_study_control/ibd.v4.RData")); class(kingibd)
dim(kingibd$kinship) # 14160 14160 - different dimensions
# contains a bunch of extras (dups, controls, no post)

# assign sample.id to names of matrix
length(kingibd$sample.id) # 14160
rownames(kingibd$kinship) <- kingibd$sample.id
colnames(kingibd$kinship) <- kingibd$sample.id
kingibd$kinship[1:5,1:5]

thresh <- 2^(-9/2)

# load kinship matrix, adj for x chr PC 1-2
kinMat <- getobj("olga_application/kcMat_xchr_xPC12adj_unrel3Deg.RData")

pdf("xkc_hist.pdf")
hist(as.vector(kinMat[lower.tri(kinMat,diag=FALSE)]),breaks=100,
     main="X KC for 12,771 Pairs",xlab="X KC")
dev.off()

pdf("xkc_hist_trunc.pdf")
hist(as.vector(kinMat[lower.tri(kinMat,diag=FALSE)]),ylim=c(0,1000),breaks=100,
     main="X KC for 12,771 Pairs",xlab="X KC")
dev.off()

rels <- apply(kinMat,1,function(x){sum(x>thresh)})

kinMat.unrelAuto <- kinMat[rownames(kinMat)%in%unrel.set,colnames(kinMat)%in%unrel.set]
dim(kinMat.unrelAuto) # 10614 10614

summary(as.vector(kinMat.unrelAuto[lower.tri(kinMat.unrelAuto,diag=FALSE)]))

pdf("xkc_unrel_hist.pdf")
hist(as.vector(kinMat.unrelAuto[lower.tri(kinMat.unrelAuto,diag=FALSE)]),breaks=100,
     main="X KC for 10,614 Unrel Pairs", xlab="X KC")
dev.off()

pdf("xkc_unrel_hist_trunc.pdf")
hist(as.vector(kinMat.unrelAuto[lower.tri(kinMat.unrelAuto,diag=FALSE)]),ylim=c(0,1000),breaks=100,
     main="X KC for 10,614 Unrel Pairs", xlab="X KC")
dev.off()

# stratify by ff, fm, mm pairs
kinRes <- getobj("olga_application/pcRelate_Xchr_xPC12adj_unrel3Deg.RData")
kc <- kinRes$kinship
kc.unrelAuto <- kc[kc$ID1 %in% unrel.set & kc$ID2 %in% unrel.set,]

kc.unrelAuto <- merge(kc.unrelAuto,pData(scan)[,c("scanID","sex")],by.x="ID1",by.y="scanID",all.x=TRUE)
kc.unrelAuto <- merge(kc.unrelAuto,pData(scan)[,c("scanID","sex")],by.x="ID2",by.y="scanID",all.x=TRUE,suffixes=c(".1",".2"))
kc.unrelAuto$sexPair <- paste(kc.unrelAuto$sex.1,kc.unrelAuto$sex.2,sep="-")
kc.unrelAuto$sexPair[kc.unrelAuto$sexPair=="M-F"] <- "F-M"

pdf("xkc_unrel_hist_bySexPair.pdf",width=14,height=9)
ggplot(kc.unrelAuto,aes(x=kin)) + geom_histogram(binwidth=0.003) + 
  facet_grid(~sexPair) + theme_bw() + xlab("X KC") 
dev.off()

pdf("xkc_unrel_hist_bySexPair_trunc.pdf",width=14,height=9)
ggplot(kc.unrelAuto,aes(x=kin)) + geom_histogram(binwidth=0.003) + 
  facet_grid(~sexPair) + theme_bw() + xlab("X KC") + coord_cartesian(ylim=c(0, 1000))
dev.off()

# only plot the FF pairs
pdf("xkc_unrel_hist_FFpairs_trunc.pdf")
ggplot(kc.unrelAuto[sexPair=="F-F",],aes(x=kin)) + geom_histogram(binwidth=0.003) + 
  theme_bw() + xlab("X KC") + coord_cartesian(ylim=c(0, 1000)) + ggtitle("X KC for F-F Unrelated Pairs")
dev.off()

kinMat4th <- getobj("olga_application/kcMat_prunedxchr_xPC12adj.RData")
kinMatCAm <- kinMat[rownames(kinMat)%in%rownames(kinMat4th),colnames(kinMat)%in%colnames(kinMat4th)]

png("tmp_comparison.png")
plot(kinMat4th[lower.tri(kinMat4th,diag=F)],kinMatCAm[lower.tri(kinMatCAm,diag=F)],
     xlab="4th Deg X KC",ylab="3rd Deg X KC",col="#00000030")
abline(0,1)
dev.off()

divMat <- kingibd$kinship
inclIdx <- is.element(rownames(divMat),rownames(kinMat))
divMat <- divMat[inclIdx,inclIdx]

ids <- pcairPartition(kinMat=kinMat,divMat=divMat,kin.thresh=thresh,div.thresh=-thresh) # took about 2.5 hrs
save(ids,file="olga_application/unrelIDs_kinMat_divMat_unrelDeg3.RData")
names(ids)
length(ids[[1]]); length(ids[[2]]) # 12074 | 697
divMat <- ids

ids <- pcairPartition(kinMat=kinMat,kin.thresh=thresh,div.thresh=-thresh) # took about 2.5 hrs
save(ids,file="olga_application/unrelIDs_kinMatOnly_unrelDeg3.RData")
names(ids) # rels | unrels
length(ids[[1]]); length(ids[[2]]) # 12072 | 699
kinMat <- ids

sum(is.element(divMat$unrel,kinMat$unrel)) # 386
sum(is.element(kinMat$unrel,divMat$unrel)) # 386

library(VennDiagram)
png("tmp_venn.png")
venn.diagram(x=list("kinMatOnly"=as.integer(kinMat$unrel),"inclDivMat"=as.integer(divMat$unrel)),filename="tmp_venn.png")
dev.off()

rm(list=ls())


#####
# 67. Run PC-Relate on set of 3600 SNPs, not on X chr

library(GWASTools)
library(OLGApipeline)
library(QCpipeline)
library(MASS)
library(GENESIS)

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x/")
source("pcrelate.R")

scan <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_scanAnnot_v54_AMS.RData"))
head(pData(scan)) 

snp.pruned <- get(load("olga_application/snp_sel_auto_ldPruned.RData"))
length(snp.pruned) # 152391

snp <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_HCHS_Custom_15041502_B3_all37_v25_AMS.RData"))
snpChrs <- snp$chromosome[is.element(snp$snpID,snp.pruned)]
table(snpChrs) # chr 19 has 4125, that'll work

snp.pruned <- snp.pruned[snpChrs==19] 
length(snp.pruned) # 4125

# make genoData object
gdsobj <- GdsGenotypeReader("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/netCDF/subjects/OLGA_subj_geno.gds")
genoData <- GenotypeData(gdsobj,scanAnnot=scan[scan$subj.plink,])

unrel.set <- scan$scanID[scan$unrelated.pcair&is.na(scan$gengrp6.outliers)&
                           scan$geno.cntl==0&scan$subj.plink]
scanIncl <- scan$scanID[is.na(scan$gengrp6.outliers)&scan$subj.plink&scan$geno.cntl==0]

head(unrel.set); head(scanIncl)
length(unrel.set); length(scanIncl) # 10624 | 12784

## run pc-relate adj for auto pcs 1-5
autoPC <- get(load("olga_application/pca_prunedAutosnps_10272unrelPlusRel.RData"))

pcMat <- autoPC$eigenvect[,c(1:5)]
dim(pcMat) # 12747 2
rownames(pcMat) <- autoPC$sampleInfo$sample.id

### need to exclude 13 people with an entirely missing x chr 
# from chromosome.anomalies.csv file, they are study samples with x chr anomaly
# thus, their entire x chr is filtered out. 

chromAnoms <- read.csv("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/dbGaP/To_dbGaP/chromosome_anomalies/chromosome_anomalies.csv",
                       header=TRUE,as.is=T)
ids <- chromAnoms$scanID[chromAnoms$chromosome=="X"&chromAnoms$filter&chromAnoms$whole.chr&
                           !is.element(substr(chromAnoms$SUBJECT_ID,1,2),"NA")]
length(ids) # 13
chromAnoms[is.element(chromAnoms$scanID,ids),] # perfect! exactly what we want

length(scanIncl) # 12784
scanIncl <- scanIncl[!is.element(scanIncl,ids)]
length(scanIncl) # 12771

unrel.set <- unrel.set[!is.element(unrel.set,ids)]
length(unrel.set) # 10614

# take these individs out of pcMat too
pcMat <- pcMat[is.element(rownames(pcMat),scanIncl),]
dim(pcMat) # 12734 5

# need to order pcMat the same as scanIncl
scanIncl <- scanIncl[is.element(scanIncl,rownames(pcMat))]

pcMat <- cbind(pcMat,as.integer(rownames(pcMat)))
pcMatTmp <- pcMat[order(pcMat[,6]),]

pcMat <- pcMatTmp[,c(1:5)]

unrel.set <- unrel.set[is.element(unrel.set,scanIncl)]

rel <- pcrelate(genoData, unrel.set=unrel.set, scan.include=scanIncl, Xchr=FALSE, pcMat=pcMat,
                snp.include=snp.pruned)
names(rel)

save(rel,file="olga_application/pcRelate_chr19_4125snps_autoPC15adj_unrel3Deg.RData")

## make into matrix
tmp <- matrix(NA,nrow=length(scanIncl),ncol=length(scanIncl))
length(tmp[lower.tri(tmp,diag=FALSE)]) # 81071011

kc <- rel$kinship
length(kc[,"kin"]) # 81071011
tmp[lower.tri(tmp,diag=FALSE)] <- kc[,"kin"]

# mirror to upper triangle
tmp[upper.tri(tmp,diag=FALSE)] <- t(tmp)[upper.tri(tmp,diag=FALSE)]

rownames(tmp) <- scanIncl
colnames(tmp) <- scanIncl

diag(tmp) <- 1 + rel$inbreed$f

save(tmp,file="olga_application/kcMat_chr19_4125snps_autoPC15adj_unrel3Deg.RData")

rm(list=ls())


#####
# 68. Run regression to check corr of auto SNPs with X PCs

library(GWASTools)
library(OLGApipeline)
library(QCpipeline)
library(MASS)
library(GENESIS)
library(lmtest)

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x/")

# run 2 models: 
# auto SNP ~ auto PC 1-5
# auto SNP ~ auto PC 1-5 + X chr PC 1-2

scan <- getobj("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_scanAnnot_v54_AMS.RData")
head(pData(scan)) 

unrel.set <- scan$scanID[scan$unrelated.pcair&!is.element(scan$gengrp6.outliers,"AsianOutliers")&
                           scan$geno.cntl==0&scan$subj.plink]
scanIncl <- scan$scanID[!is.element(scan$gengrp6.outliers,"AsianOutliers")&scan$subj.plink&scan$geno.cntl==0]

head(unrel.set); head(scanIncl)
length(unrel.set); length(scanIncl) # 10624 | 12784

# adj for auto PCs 1-5

db <- getDb('olga_analysis')
evs <- dbGetEvData(db, evSetID=2, evNums=1:5)
dbDisconnect(db)

head(evs); dim(evs) # 12784 5

xpc <- getobj("olga_application/pca_prunedXsnps_10272unrel12784rel_unrelDeg3.RData")
dim(xpc$vectors) # 12784 10

allequal(rownames(xpc$vectors),evs$scan_id) # TRUE; great!

allequal(scan$scanID[is.element(scan$scanID,scanIncl)],rownames(xpc$vectors)) # TRUE
scan$autoEV1[is.element(scan$scanID,scanIncl)] <- evs$EV1
scan$autoEV2[is.element(scan$scanID,scanIncl)] <- evs$EV2
scan$autoEV3[is.element(scan$scanID,scanIncl)] <- evs$EV3
scan$autoEV4[is.element(scan$scanID,scanIncl)] <- evs$EV4
scan$autoEV5[is.element(scan$scanID,scanIncl)] <- evs$EV5

scan$xEV1[is.element(scan$scanID,scanIncl)] <- xpc$vectors[,1]
scan$xEV2[is.element(scan$scanID,scanIncl)] <- xpc$vectors[,2]
  
saveas(scan,"olga_application/scanAnnot_lrtest.RData")

rm(list=ls())

##
# called in batch:
# lmtest_snps.R and .Rout for each chr 1-22 in parallel
# saved lrTest_autoSNPs/lrRes_chrXX.txt for each chr
# need to cat those in unix, make some plots

## 
# called in batch:
# lmtest_plotRes.R and .Rout


#####
# 69. Parse chr 19 PC-Relate results

library(GWASTools)
library(OLGApipeline)
library(QCpipeline)
library(MASS)
library(GENESIS)
library(ggplot2)

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x/")

kc <- get(load("olga_application/pcRelate_chr19_4125snps_autoPC15adj_unrel3Deg.RData"))
names(kc)
kin <- kc$kinship

scan <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_scanAnnot_v54_AMS.RData"))

head(pData(scan)) 
unrel.set <- scan$scanID[scan$unrelated.pcair&is.na(scan$gengrp6.outliers)&
                           scan$geno.cntl==0&scan$subj.plink]

kc.unrelAuto <- kin[is.element(kin$ID1,unrel.set)&is.element(kin$ID2,unrel.set),]
dim(kc.unrelAuto)
head(kc.unrelAuto)
summary(kc.unrelAuto$kin)
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#-0.098210 -0.006858  0.001194  0.001345  0.009370  0.311000 

pdf("kc_chr19_autounrel_hist.pdf")
ggplot(kc.unrelAuto,aes(x=kin)) + geom_histogram(binwidth=0.003) + xlab("KC") + 
  ggtitle("KC calculated from 4,125 Chr 19 SNPs\nAuto Unrel Pairs") + theme_bw() 
dev.off()

pdf("kc_chr19_autounrel_hist_trunc.pdf")
ggplot(kc.unrelAuto,aes(x=kin)) + geom_histogram(binwidth=0.003) + xlab("KC") + coord_cartesian(ylim=c(0, 1000)) +
  ggtitle("KC calculated from 4,125 Chr 19 SNPs\nAuto Unrel Pairs") + theme_bw() 
dev.off()

rm(list=ls())


#####
# 70. Plot distribution of X chr markers across the chr

library(GWASTools)
library(OLGApipeline)
library(QCpipeline)
library(MASS)
library(GENESIS)
library(ggplot2)

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x/")
source("pcrelate.R")

scan <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_scanAnnot_v54_AMS.RData"))
head(pData(scan)) 

snp.pruned <- get(load("olga_application/snp_sel_xChr_ldPruned.RData"))
length(snp.pruned) # 3600

snp <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_HCHS_Custom_15041502_B3_all37_v25_AMS.RData"))

xchrA <- pData(snp[snp$chromosome==23,])
dim(xchrA)
range(xchrA$position) # 2699676 154921994
# it's 155,270,560 bp long, accd to ucsc genome browser

pseudoautosomal.hg19
#       chrom region start.base  end.base
#X.PAR1     X   PAR1      60001   2699520
#X.PAR2     X   PAR2  154931044 155260560
#X.XTR      X    XTR   88395830  92583067


# Plot the densities of snps in the bed file for each chr seperately
snpDensity <- ggplot(xchrA) + 
  geom_rect(aes(xmin=60001, xmax=2699520, ymin=0, ymax=Inf),fill='gray80', alpha=0.8) +
  geom_rect(aes(xmin=154931044, xmax=155260560, ymin=0,ymax=Inf),fill='gray80', alpha=0.8) +
  geom_rect(aes(xmin=88395830, xmax=92583067, ymin=0,ymax=Inf),fill='gray80', alpha=0.8) +
  geom_histogram(aes(x=position),binwidth=1e4) + # pick a binwidth that is not too small 
  ggtitle("Density of 55,905 X Chr SNPs Across hg19\nBinned by 1e4 bases") +
  xlab("Position") + 
  ylab("SNP density") + 
  theme_bw() + 
  annotate("text", x = 1319760, y = 30, label = "PAR1") + 
  annotate("text", x = (164758+154931044), y = 30, label = "PAR2") + 
  annotate("text", x = (2093618+88395830), y = 30, label = "XTR") + 
  annotate("text", x=60.6*1e6, y=10, label="centromere",angle=90)

# save the plot to .png file
png("chrX_density.png",width=960)
print(snpDensity)
dev.off()

# do for chr 19 for comparison
chr19 <- pData(snp[snp$chromosome==19,])
png("chr19_density.png",width=960)
ggplot(chr19) + geom_histogram(aes(x=position),binwidth=1e4) +
  ggtitle("Density of 52,890 Chr 19 SNPs Across hg19\nBinned by 1e4 bases") +
  xlab("Position") + ylab("SNP density") + theme_bw() +
  annotate("text",x=26.5*1e6,y=10,label="centromere",angle=90)
dev.off()

rm(list=ls())


#####
# 71. QQ plots for simple MLM and MLM-X on RBC trait, with lambda_GC in legend

library(QCpipeline); library(OLGApipeline)
library(GWASTools)
library(dplyr); library(readr)
library(tidyr)

setwd("/projects/geneva/geneva_sata/caitlin/olga_xchr_assoc/Assoc/")

# read in xchr results, the autosomal model
(mac.thresh <- 30)
n <- 12490
maf.eff <- format(quadSolveMAF(mac.thresh, n), digits=3)

assc <- getobj("assoc_auto_316987_chr23.RData")

filt.imp <- which(assc$type %in% 0)
assc$composite.filter <- TRUE
assc$maf.filt <- assc$effN > mac.thresh
filt.obs <- which(assc$type %in% c(2,3))

filt.obs.maf <- intersect(filt.obs, which(assc$maf.filt))
# this is the filter i want

assc <- assc[filt.obs.maf,]
totalRes <- assc

##
# do the same exercise for the other assoc test results
assc <- getobj("assoc_316987_chr23.RData")

filt.imp <- which(assc$type %in% 0)
assc$composite.filter <- TRUE
assc$maf.filt <- assc$effN > mac.thresh
filt.obs <- which(assc$type %in% c(2,3))

filt.obs.maf <- intersect(filt.obs, which(assc$maf.filt))
# this is the filter i want

assc <- assc[filt.obs.maf,]

mlmX <- assc

allequal(mlmX$snpID,totalRes$snpID) # TRUE
mlmX$method <- "MLM-X"
totalRes$method <- "simple MLM"

res <- rbind(mlmX[,c("chromosome","snpID","method","pval")],
             totalRes[,c("chromosome","snpID","method","pval")])
res <- tbl_df(res)

res <- mutate(res,pval_10=-log10(pval))
res_method <- group_by(res,method)
res_method <- mutate(res_method,n=n())
n <- 43868
a <- 1:n
b <- a/n
x <- -log10(b)
res_method <- mutate(res_method,spval=-log10(sort(pval)))
res_method <- mutate(res_method,x=x)

# make a plot of each pval comparing mlm-x to simple mlm-x, across the autosomes
library(ggplot2)
library(gridExtra); library(RGraphics)

txt <- paste("lambda[GC] == 1.116~simple~MLM")
txt2 <- paste("lambda[GC] == 1.043~MLM-X")

png("../Plots/xChrQQ_mlmX_vsAutoModel.png",width=1080,height=760)
ggplot(res_method,aes(x,spval,color=method)) + geom_point(aes(color=method,name=" ")) +
  xlab(expression(paste(-log[10],"(expected P)"))) +
  ylab(expression(paste(-log[10], "(observed P)"))) +
  theme_bw() + geom_abline(intercept=0) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=18,face="bold"),
        legend.text = element_text(size = 18),legend.key = element_rect(colour = "white")) +
  annotate("text", x = 0.5, y = 16.8, label = txt,parse=TRUE,size=6) +
  annotate("text", x = 0.42, y = 16, label = txt2, parse=TRUE,size=6) 
dev.off()

rm(list=ls())


#####
# 72. Run regression to check corr of X chr SNPs with auto PCs

# run 2 models: 
# x SNP ~ x PC 1-2
# x SNP ~ x PC 1-2 + auto PC 1-5

##
# called in batch:
# lmtest_snps_xchr.R and .Rout for chr 23 in parallel
# saved lrTest_autoSNPs/lrRes_chrXX.txt for each chr
# need to cat those in unix, make some plots

## 
# called in batch:
# plot_qq_manh_xchr.R


#####
# 73. Make barplots of var comp estimates w and w/o X chr KC for RBC

library(ggplot2)
library(tidyr); library(dplyr)

mlmX <- data.frame("block group"=0.00396,"household"=0.04950,
                   "autosomal kinship"=0.28453,"x kinship"=0.02935,
                   "residual"=0.63266,model="MLM-X")
simp <- data.frame("block group"=0.00363,"household"=0.04945,
                   "autosomal kinship"=0.28473,"x kinship"=0.00,
                   "residual"=0.66218,model="Simple MLM")

mlmX <- gather(mlmX,varComp,varEst,-model)
simp <- gather(simp,varComp,varEst,-model)
toPl <- rbind(mlmX,simp)

pdf("/projects/geneva/geneva_sata/caitlin/olga_xchr_assoc/prop_var_barchart_RBC.pdf")
ggplot(toPl,aes(x=model,y=varEst,fill=varComp)) + geom_bar(position = "fill",stat = "identity") + theme_bw()+
  ylab("Proportion Variance") + scale_fill_brewer()
dev.off()

rm(list=ls())


#####
# 74. Compare chr 19 results to X chr KC FF pairs

library(GWASTools)
library(OLGApipeline)
library(QCpipeline)
library(MASS)
library(GENESIS)
library(ggplot2)

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x/")

kc19 <- get(load("olga_application/pcRelate_chr19_4125snps_autoPC15adj_unrel3Deg.RData"))
names(kc19)
kin19 <- kc19$kinship

# stratify by ff, fm, mm pairs
kinRes <- getobj("olga_application/pcRelate_Xchr_xPC12adj_unrel3Deg.RData")
kinX <- kinRes$kinship

## subset to KC > 0.1
kin19 <- kin19[kin19$kin>0.1,]
kinX <- kinX[kinX$kin>0.1,]
dim(kin19) # 3047 7
dim(kinX) # 21388 7

kin19$id12 <- paste(kin19$ID1,kin19$ID2)
kinX$id12 <- paste(kinX$ID1,kinX$ID2)
kinBoth <- merge(kin19,kinX,all.x=TRUE,by="id12")
kinBoth <- kinBoth[!is.na(kinBoth$kin.y),]
dim(kinBoth) # 2230 15

# merge in sex info
scan <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_scanAnnot_v54_AMS.RData"))

kin19 <- merge(kin19,pData(scan)[,c("scanID","sex")],by.x="ID1",by.y="scanID",all.x=TRUE)
kin19 <- merge(kin19,pData(scan)[,c("scanID","sex")],by.x="ID2",by.y="scanID",all.x=TRUE,suffixes=c(".1",".2"))

kin19 <- kin19[kin19$sex.1=="F"&kin19$sex.2=="F",]
dim(kin19) # 1309 10

kinBoth <- merge(kin19,kinX,by="id12",all.x=TRUE,suffixes=c(".19",".x"))
kinBoth <- kinBoth[!is.na(kinBoth$kin.x),]
dim(kinBoth) # 1090 17

cor(kinBoth$kin.x,kinBoth$kin.19) # 0.38703



ggplot(kinBoth,aes(x=kin.x,y=kin.19)) + geom_point(color=relationship) + geom_abline() +
  theme_bw()
dev.off()


head(pData(scan)) 
unrel.set <- scan$scanID[scan$unrelated.pcair&is.na(scan$gengrp6.outliers)&
                           scan$geno.cntl==0&scan$subj.plink]

kc.unrelAuto <- kin[is.element(kin$ID1,unrel.set)&is.element(kin$ID2,unrel.set),]
dim(kc.unrelAuto)
head(kc.unrelAuto)
summary(kc.unrelAuto$kin)
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#-0.098210 -0.006858  0.001194  0.001345  0.009370  0.311000 

pdf("kc_chr19_autounrel_hist.pdf")
ggplot(kc.unrelAuto,aes(x=kin)) + geom_histogram(binwidth=0.003) + xlab("KC") + 
  ggtitle("KC calculated from 4,125 Chr 19 SNPs\nAuto Unrel Pairs") + theme_bw() 
dev.off()


kinMat.unrelAuto <- kinMat[rownames(kinMat)%in%unrel.set,colnames(kinMat)%in%unrel.set]
dim(kinMat.unrelAuto) # 10614 10614

summary(as.vector(kinMat.unrelAuto[lower.tri(kinMat.unrelAuto,diag=FALSE)]))

pdf("xkc_unrel_hist.pdf")
hist(as.vector(kinMat.unrelAuto[lower.tri(kinMat.unrelAuto,diag=FALSE)]),breaks=100,
     main="X KC for 10,614 Unrel Pairs", xlab="X KC")
dev.off()

pdf("xkc_unrel_hist_trunc.pdf")
hist(as.vector(kinMat.unrelAuto[lower.tri(kinMat.unrelAuto,diag=FALSE)]),ylim=c(0,1000),breaks=100,
     main="X KC for 10,614 Unrel Pairs", xlab="X KC")
dev.off()

# stratify by ff, fm, mm pairs
kinRes <- getobj("olga_application/pcRelate_Xchr_xPC12adj_unrel3Deg.RData")
kc <- kinRes$kinship
kc.unrelAuto <- kc[kc$ID1 %in% unrel.set & kc$ID2 %in% unrel.set,]

kc.unrelAuto <- merge(kc.unrelAuto,pData(scan)[,c("scanID","sex")],by.x="ID1",by.y="scanID",all.x=TRUE)
kc.unrelAuto <- merge(kc.unrelAuto,pData(scan)[,c("scanID","sex")],by.x="ID2",by.y="scanID",all.x=TRUE,suffixes=c(".1",".2"))
kc.unrelAuto$sexPair <- paste(kc.unrelAuto$sex.1,kc.unrelAuto$sex.2,sep="-")
kc.unrelAuto$sexPair[kc.unrelAuto$sexPair=="M-F"] <- "F-M"

rm(list=ls())


#####
# 75. Rerun CAnD with new p-value combining method

calc_combP <- function(dat){
  sig2.i <- apply(dat,1,var) # this is individual level variance
  
  # get the w_cc',i for each individual i
  m <- ncol(dat)
  n <- nrow(dat)
  tmp1 <- expand.grid(1:m,1:m) # all pairwise combos of chrs
  # remove rows where chrs are =
  cc <- tmp1[,1]==tmp1[,2]
  tmp1 <- tmp1[!cc,]
  w <- rep(NA,nrow(dat))
  for(i in 1:nrow(dat)){
    abar <- mean(as.numeric(dat[i,]))
    tmp2 <- cbind(as.numeric(dat[i,tmp1[,1]])-abar,as.numeric(dat[i,tmp1[,2]])-abar)
    w[i] <- mean(tmp2[,1]*tmp2[,2])
  }
  # w_{cc',i} should be zero under the null, since we are adj for the mean ancestry w/in an individ
  # this is relative to an individ, so there will be no correlation
  # if this is relative to the population itself, there will be correlation
  
  # need a sig2_c for each pair of chromosomes
  sig2.c <- rep(NA,ncol(dat))
  tstat <- rep(NA,ncol(dat))
  for(i in 1:m){
    tmp <- dat[,-i]
    pool <- apply(tmp,1,mean)
    d <- dat[,i]-pool
    dbar <- mean(d)
    #  sdd <- sd(d)/sqrt(n)

    sig2.c[i] <- sum((d-dbar)^2)/(n*(n-1))
    tstat[i] <- dbar/sqrt(sig2.c[i])
    #  tstat[i] <- t.test(nam[,i],pool,paired=TRUE)$statistic
  }
  
  sig.matrix <- diag(nrow=m,ncol=m)
  for(i in 1:m){
    for(j in 1:m){
      denom <- n^2*sqrt(sig2.c[i])*sqrt(sig2.c[j])
      sig.matrix[i,j] <- (1/denom)*sum((m/(m-1)^2)*(w-sig2.i))
    }
  }
  diag(sig.matrix) <- 1
  mean(sig.matrix[lower.tri(sig.matrix)]) # -0.054286
  
  # calculate new stat
  (newstat <- tstat%*%solve(sig.matrix)%*%tstat) # 49.37731
  return(pchisq(newstat,df=m,lower.tail=FALSE))
}

library(GWASTools)

res <- getobj("local_ancestry_byChr.RData")
scan <- getobj("admix_xchr_aftermerge.RData")
colnames(scan)[3:5] <- c("AFR.x","NAM.x","EUR.x")               
# merge in admixture results for x chr since don't have local ancestry estimates
res <- merge(res,scan,by="scanID")

res$bkgrd[is.element(res$bkgrd,c("Other","Unknown"))] <- "Other/Unknown"

## need to subset so that only unrel samples are included
unrel <- getobj("unrelated_pcair_deg4.RData")

sum(is.element(res$scanID,unrel))
sum(!is.element(res$scanID,unrel))

res <- res[is.element(res$scanID,unrel),]

# get avg local ancestry across the autosomes
colnames(res)
afrCols <- 24:45
colnames(res)[afrCols]
res$AFR.auto <- rowMeans(res[,afrCols])

eurCols <- 2:23
colnames(res)[eurCols]
res$EUR.auto <- rowMeans(res[,eurCols])

namCols <- 46:67
colnames(res)[namCols]
res$NAM.auto <- rowMeans(res[,namCols])


## get CAnD pvalues across the autosomes for PR samples, combine using new method

subgrpRes <- data.frame("bkgrd"=unique(res$bkgrd),"n"=NA)
newCols <- data.frame(matrix(NA,nrow=nrow(subgrpRes),ncol=23*3))
colnames(newCols) <- paste0(paste0("chr",c(23,1:22)),rep(c(".NAM",".EUR",".AFR"),each=23))
subgrpRes <- cbind(subgrpRes,newCols)
colnames(subgrpRes)

namCols <- c(70,46:67)
eurCols <- c(71,2:23)
afrCols <- c(69,24:45)

namAt <- c(70,74)
eurAt <- c(71,73)
afrAt <- c(69,72)

colnames(res)[namCols]; colnames(res)[eurCols]; colnames(res)[afrCols]
nonPres <- subgrpRes
pooledres <- subgrpRes

i <- which(subgrpRes$bkgrd=="PuertoRican")
pop <- subgrpRes$bkgrd[i]
subgrpRes$n[i] <- sum(is.element(res$bkgrd,pop))

# genome-wide pvalue
calc_combP(res[res$bkgrd==pop,namCols]) # 4.11335e-173
calc_combP(res[res$bkgrd==pop,eurCols]) # 5.288241e-185
calc_combP(res[res$bkgrd==pop,afrCols]) # 1.797041e-26

# autosomal-wide pvalue 
calc_combP(res[res$bkgrd==pop,namCols[2:23]]) # 1.253574e-49
calc_combP(res[res$bkgrd==pop,eurCols[2:23]]) # 2.4385e-36
calc_combP(res[res$bkgrd==pop,afrCols[2:23]]) # 7.859006e-15

rm(list=ls())
