
# script to run var component estimation and MLM assoc test on phenotype values simulated 100x

source("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/mconomos/R/estVarComp.R")
source("/projects/geneva/gcc-fs2/OLGA/genotype/phaseIa/mconomos/R/assocTestMixedModel_v7.R")

###
# R fxn to read the GRM binary file
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


# files used:
sims <- read.table("/projects/geneva/geneva_sata/caitlin/grm_xchr/hapmap3_r3_b37_fwd.consensus.qc.poly_Ilmn1M_causSNPs_maf2_effSz3.phen", header=FALSE,as.is=TRUE)
grmMatX <- ReadGRMBin("/projects/geneva/geneva_sata/caitlin/grm_xchr/hapmap3_r3_b37_fwd.consensus.qc.poly_Ilmn1M")
grmMatA <- ReadGRMBin("/projects/geneva/geneva_sata/caitlin/grm_xchr/hapmap3_r3_b37_fwd.consensus.qc.poly_Ilmn1M_autos")
scan <- get(load("/projects/geneva/geneva_sata/HapMap/HapMap3_r3/sample_snp_annot/hapmap3_r3_b36.sample.annot.v6.RData"))
pcRes <- get(load("/projects/geneva/geneva_sata/caitlin/grm_xchr/pca_results.RData"))

### start section of things i only have to do once
# read in the GRM matrices
library(calibrator)

kinship_x <- diag(grmMatX$diag)
kinship_x[upper.tri(kinship_x)] <- grmMatX$off
kinship_x <- symmetrize(kinship_x)
stopifnot(isSymmetric(kinship_x))

kinship_autos <- diag(grmMatA$diag)
kinship_autos[upper.tri(kinship_autos)] <- grmMatA$off
kinship_autos <- symmetrize(kinship_autos)
stopifnot(isSymmetric(kinship_autos))

# need to make a scanAnnot object
scan

# need row and col names of the kinship matrices to be the scanIDs
idMap <- grmMatA$id
idMap$order <- 1:nrow(idMap)
idMap <- merge(idMap,pData(scan)[,c("coriell.id","scanID")],by.x="V2",by.y="coriell.id")
idMap <- idMap[order(idMap$order),]
rownames(kinship_autos) <- colnames(kinship_autos) <- idMap$scanID
rownames(kinship_x) <- colnames(kinship_x) <- idMap$scanID

# merge in PC1-5
evs <- pcRes$eigenvect

idMap <- idMap[order(idMap$scanID),]
scan <- scan[order(idMap$order),]
stopifnot(all(pcRes$sample.id==scan$coriell.id))
scan$PC1 <- evs[,1]
scan$PC2 <- evs[,2]
scan$PC3 <- evs[,3]
scan$PC4 <- evs[,4]
scan$PC5 <- evs[,5]

# make genotype data object
data <- GdsGenotypeReader("/projects/geneva/geneva_sata/caitlin/grm_xchr/hapmap3_r3_b37_fwd.consensus.qc.poly_Ilmn1M.gds")
idMap <- idMap[order(idMap$order),]
scan <- scan[order(idMap$order),]
scan$oldScanID <- scan$scanID
scan$scanID <- scan$coriell.id
genoData <- GenotypeData(data,scanAnnot=scan)

scan$scanID <- scan$oldScanID

chr <- getChromosome(data)
endSNP <- which.max(chr==25)-1 # 944675

varResults <- data.frame(phenoSim=1:100,varX=NA,varAuto=NA,varErr=NA)
varResults_exclX <- varResults

# merge in all simulated phenotypes to scan annot
scanMerged <- merge(pData(scan),sims,by.x="coriell.id",by.y="V2",all.x=TRUE)
scan <- scan[order(scan$coriell.id),]

### end section of things i only have to do once 

## loop through simulated phenotype values and perform for each one of the 100 simulations
for(i in 1:100){
  thisPheno <- paste("V",(i+2),sep="")
  scan <- scan[order(scan$coriell.id),]
  stopifnot(all(scanMerged$coriell.id==scan$coriell.id))
  scan$pheno_sim <- scanMerged[,thisPheno]

  # now they need to be in the same order
  scan <- scan[order(scan$scanID),]
  idMap <- idMap[order(idMap$scanID),]
  scan <- scan[order(idMap$order),]

  # use matt's code
  covList <- list(kinship_autos,kinship_x)
  names(covList) <- c("kinshipAuto","kinshipX")
  varComp <- estVarComp(scan,covMatList=covList,"pheno_sim",covar.vec=c("PC1","PC2","PC3","PC4","PC5"))
  
  # store varComp$varComp estimates
  varResults$varX[i] <- varComp$varComp["V_kinshipX"]
  varResults$varAuto[i] <- varComp$varComp["V_kinshipAuto"]
  varResults$varErr[i] <- varComp$varComp["V_E"]
  
  # call MLM with these results
  # need to change col/rownames back to coriell ids
  idMap <- idMap[order(idMap$order),]
  cholSig <- varComp[["cholSigmaInv"]]
  stopifnot(all(idMap$scanID==colnames(cholSig)))
  colnames(cholSig) <- rownames(cholSig) <- idMap$V2
  scan$scanID <- scan$coriell.id
  genoData <- GenotypeData(data,scanAnnot=scan)
  mmRes <- assocTestMixedModel(genoData,snpStart=1,snpEnd=endSNP,cholSig,outcome="pheno_sim",covar.vec=c("PC1","PC2","PC3","PC4","PC5"))

  scan$scanID <- scan$oldScanID
  
  if(i==1){ # create results matrix on first loop through
    pvalueResults <- data.frame(matrix(NA,nrow=nrow(mmRes),ncol=(400+3)))
    colnames(pvalueResults) <- c("SNP", "CHR", "MAF", paste("pval_sim",1:100,sep=""),paste("Est_sim",1:100,sep=""),paste("SE_sim",1:100,sep=""),
                                 paste("Stat_sim",1:100,sep=""))
    pvalueResults$SNP <- mmRes$snpID
    pvalueResults$CHR <- mmRes$chr
    pvalueResults$MAF <- mmRes$MAF
    
    pvalueResults_exclX <- pvalueResults
  }
    
  stopifnot(all(pvalueResults$SNP==mmRes$SNP))
  thisCol <- c(paste("pval_sim",i,sep=""),paste("Est_sim",i,sep=""),paste("SE_sim",i,sep=""),paste("Stat_sim",i,sep=""))
  pvalueResults[,thisCol] <- mmRes[,c("pval","Est","SE","Stat")]
  
  
  
  ##### do the same thing excluding the X chr adjustment
  # now they need to be in the same order
  scan <- scan[order(scan$scanID),]
  idMap <- idMap[order(idMap$scanID),]
  scan <- scan[order(idMap$order),]
  
  # use matt's code
  covList <- list(kinship_autos)
  names(covList) <- c("kinshipAuto")
  varComp <- estVarComp(scan,covMatList=covList,"pheno_sim",covar.vec=c("PC1","PC2","PC3","PC4","PC5"))
  
  # store varComp$varComp estimates
  varResults_exclX$varAuto[i] <- varComp$varComp["V_kinshipAuto"]
  varResults_exclX$varErr[i] <- varComp$varComp["V_E"]
  
  # call MLM with these results
  # need to change col/rownames back to coriell ids
  idMap <- idMap[order(idMap$order),]
  cholSig <- varComp[["cholSigmaInv"]]
  stopifnot(all(idMap$scanID==colnames(cholSig)))
  colnames(cholSig) <- rownames(cholSig) <- idMap$V2
  scan$scanID <- scan$coriell.id
  genoData <- GenotypeData(data,scanAnnot=scan)
  mmRes <- assocTestMixedModel(genoData,snpStart=1,snpEnd=endSNP,cholSig,outcome="pheno_sim",covar.vec=c("PC1","PC2","PC3","PC4","PC5"))
  
  scan$scanID <- scan$oldScanID
  
  stopifnot(all(pvalueResults$SNP==mmRes$SNP))
  thisCol <- c(paste("pval_sim",i,sep=""),paste("Est_sim",i,sep=""),paste("SE_sim",i,sep=""),paste("Stat_sim",i,sep=""))
  pvalueResults_exclX[,thisCol] <- mmRes[,c("pval","Est","SE","Stat")]
    
}

save(pvalueResults,file="/projects/geneva/geneva_sata/caitlin/grm_xchr/pvalues_maf2_effSz3.RData")
save(varResults,file="/projects/geneva/geneva_sata/caitlin/grm_xchr/varCompEst_maf2_effSz3.RData")
  
save(pvalueResults_exclX,file="/projects/geneva/geneva_sata/caitlin/grm_xchr/pvalues_maf2_effSz3_exclX.RData")
save(varResults_exclX,file="/projects/geneva/geneva_sata/caitlin/grm_xchr/varCompEst_maf2_effSz3_exclX.RData")

q("no")

