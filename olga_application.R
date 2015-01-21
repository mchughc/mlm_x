
## 
# data application of the mlm_x method to the OLGA data

## Contents:
# 1. Get theoretical KC for X chromosome in OLGA samples
 # 1a. Organize .fam file
 # 1b. Call KinInbcoefX
 # 1c. Create matrix from output
# 2. Estimate relatedness on the x chromosome in the OLGA samples




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
dim(ind_ids) # 3627 4; so these are the samples we want to keep

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

for(i in 1:nrow(tmp)){
  rowInd <- which(rownames(trueKinX)==tmp$V2[i])
  colInd <- which(colnames(trueKinX)==tmp$V3[i])
  trueKinX[rowInd,colInd] <- tmp$V4[i]
}

all(trueKinX[lower.tri(trueKinX,diag=FALSE)]==0) # TRUE
table(trueKinX[upper.tri(trueKinX,diag=FALSE)]) 
#      0  0.0625 0.09375   0.125  0.1875    0.25   0.375     0.5 
#3631155       3       2      48       6     809      66     771 

table(trueKinX[lower.tri(trueKinX,diag=FALSE)]) # all zeroes

trueKinX[lower.tri(trueKinX,diag=FALSE)] <- t(trueKinX)[lower.tri(trueKinX,diag=FALSE)]
isSymmetric(trueKinX) # TRUE

# add in other genotyped samples that are unrel; reorder the columns and rows.
scan <- get(load("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/sample_snp_annot/SoL_scanAnnot_v54_AMS.RData"))
dim(scan) # 14160 95
table(scan$subj.plink) # 13204 true; same num of individs in the fam file

# just add in individs that are in the fam file that aren't listed in trueKinX file
all(ind_ids$kinIbd_individID==rownames(trueKinX)) # TRUE
all(ind_ids$kinIbd_individID==colnames(trueKinX)) # TRUE

rownames(trueKinX) <- colnames(trueKinX) <- ind_ids$orig_individID
# ok, good. now add in unrel samples

singl <- fam$V2[!is.element(fam$V2,colnames(trueKinX))]
length(singl) # 10508
table(table(fam$V1)) # 10508 singletons, so what we expect
unrelKin <- matrix(0,nrow=nrow(fam),ncol=nrow(fam))
dim(unrelKin) # 13204 13204

unrelKin[1:nrow(trueKinX),1:nrow(trueKinX)] <- trueKinX
# ok, should be ready now!

colnames(unrelKin) <- c(colnames(trueKinX),singl)
rownames(unrelKin) <- colnames(unrelKin)

# oh, need to add in the diag elements for the unrel samples
sex_info <- data.frame("ids"=colnames(unrelKin),stringsAsFactors=FALSE)
dim(sex_info) # 13204 1
sex_info2 <- merge(sex_info,pData(scan)[scan$subj.plink,c("subjectID","sex")],by.x="ids",by.y="subjectID",all.x=TRUE,sort=FALSE)
dim(sex_info2) # 13204 2

all(sex_info$ids==sex_info2$ids) # TRUE

diag(unrelKin)[sex_info2$sex=="F"] <- 0.5
diag(unrelKin)[sex_info2$sex=="M"] <- 1

table(diag(unrelKin)) # 5455 1, 7749 0.5
table(scan$sex[scan$subj.plink]) # 5455 m, 7749 f

save(unrelKin,file="xChr_kc_matrix.RData")

rm(list=ls())


#####




