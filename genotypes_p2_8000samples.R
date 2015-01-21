
# simulate genotypes for 8000 individs = 500 16-person pedigrees

setwd("/projects/geneva/geneva_sata/caitlin/mlm_x/")

library(MASS); library(parallel)
library(BiocParallel); library(GWASTools)
library(corpcor); library(doSNOW)

source("allele_drop_functions.R")

# set parameters                                                                                                                                                                    
p <- 0.2
numSNPs <- 1000 # number of times to perform the analysis                                                                                                                           

n <- 8000

TEMP=c(2,1,2,1,1,1,2,2,1,2,1,2,1,2,1,2)
SEX <- rep(NA,length(TEMP))
SEX[TEMP==1]="M"
SEX[TEMP==2]="F"

# simulate some genotypes                                                                                                                                                           
aFreqs <- rep(p,numSNPs)

#geno <- matrix(NA,nrow=numSNPs,ncol=n)

threads <- detectCores()
cl <- makeCluster(threads)
registerDoSNOW(cl)
simGenos <- foreach(i=1:500) %dopar% Family_alleles_NmarkerX(aFreqs,numSNPs,SEX)

stopCluster(cl)

# need to make the males 4*genotype; females 2*genotype to get coding for M:0,2 and for F:0,1,2                                                                                     
males <- SEX=="M"
simGenos <- bplapply(simGenos,function(x){x[,males] <- 4*x[,males];x[,!males] <- 2*x[,!males]; return(x)})

dim(simGenos)

save(simGenos,file="genotypes_p02_8000samples.RData")

rm(list=ls())