
# script to run mlm on x in batch
# first est var comp then run assoc test

# for each phenotype, fit 3 models:
# 1. incl adj for x only
# 2. incl adj for auto only
# 3. incl adj for x and auto both

model_batch_sims <- function(model){

library(GWASTools)
setwd("/projects/geneva/geneva_sata/caitlin/mlm_x/")
source("assocTestMixedModel_v7.R")
source("estVarComp.R")

# read in simulated phenotypes and geno data
scan <- get(load("scanAnnot_pheno_3models_genotypes_p2_8000samples.RData"))
scan <- ScanAnnotationDataFrame(scan)

geno <- get(load("genotypes_p02_8000samples.RData"))
genoMt <- MatrixGenotypeReader(genotype=geno,snpID=as.integer(1:nrow(geno)),
                               chromosome=as.integer(rep(23,nrow(geno))),
                               position=as.integer(1:nrow(geno)), scanID=scan$scanID)
genoData <- GenotypeData(genoMt,scanAnnot=scan)

kinX <- get(load("500Peds_xKinship.RData"))
kinAuto <- get(load("500Peds_autoKinship.RData"))

n <- 8000

kinAuto <- matrix(kinAuto,nrow=n,ncol=n)
kinX <- matrix(kinX,nrow=n,ncol=n)

colnames(kinX) <- scan$scanID
colnames(kinAuto) <- scan$scanID

# this is var comp including both x and autos
covMatListBoth <- list(kinX,kinAuto) # this is the X chr kinship matrix for all samples
names(covMatListBoth) <- c("kinshipX","kinshipAutos")

covMatListX <- list(kinX)
names(covMatListX) <- c("kinshipX")

covMatListAutos <- list(kinAuto)
names(covMatListAutos) <- c("kinshipAutos")

# loop through the phenotypes for model A: y=beta1*SNPx + noise
if(model=="A"){colSeq <- seq(from=3,to=ncol(scan),by=3)}
if(model=="B"){colSeq <- seq(from=4,to=ncol(scan),by=3)}
if(model=="C"){colSeq <- seq(from=5,to=ncol(scan),by=3)}

ct <- 1

for(i in colSeq){
  # make phenotype col called "pheno" for the causal SNP in question
  scan$pheno <- pData(scan)[,i]

  genoData <- GenotypeData(genoMt,scanAnnot=scan) # after added pheno variable
  
  varCompBoth <- estVarComp(scan,covMatList=covMatListBoth,"pheno")
  varBothci <- estVarCompCI(varCompBoth,prop=FALSE)

  cholSig <- varCompBoth[["cholSigmaInv"]]
  mmResBoth <- assocTestMixedModel(genoData,snpStart=ct,snpEnd=(ct+500),cholSig,outcome="pheno")
  
  varCompX <- estVarComp(scan,covMatList=covMatListX,"pheno")
  varXci <- estVarCompCI(varCompX,prop=FALSE)
  
  cholSig <- varCompX[["cholSigmaInv"]]
  mmResX <- assocTestMixedModel(genoData,snpStart=ct,snpEnd=(ct+500),cholSig,outcome="pheno")
  
  varCompAuto <- estVarComp(scan,covMatList=covMatListAutos,"pheno")
  varAutoci <- estVarCompCI(varCompAuto,prop=FALSE)
  
  cholSig <- varCompAuto[["cholSigmaInv"]]
  mmResAuto <- assocTestMixedModel(genoData,snpStart=ct,snpEnd=(ct+500),cholSig,outcome="pheno")
  
  resF <- list(varBothci,mmResBoth,varXci,mmResX,varAutoci,mmResAuto)
  names(resF) <- c("var_both_ci","mm_res_both","var_x_ci","mm_res_x","var_auto_ci","mm_res_auto")
  fname <- paste("model",model,"_results/mmRes_varComps_modelA_SNP",ct,".RData",sep="")
  save(resF,file=fname)
  
  ct <- ct+1  
  if(ct%%10==0){print(paste("**** ended iteration",ct,"****"))}

}
  return(NULL)
}





