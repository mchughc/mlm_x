

## check the simulated genotypes to see that they have the expected relatedness structure
source("allele_drop_functions.R")

kinFull <- get(load("100Peds_500unr_autoKinship.RData"))
kinFullX <- get(load("100Peds_500unr_xKinship.RData"))

n <- 1600 # number of samples

TEMP=c(2,1,2,1,1,1,2,2,1,2,1,2,1,2,1,2)
SEX <- rep(NA,length(TEMP))
SEX[TEMP==1]="M"
SEX[TEMP==2]="F"

sex <- rep(SEX,(n/16))

# simulate genotypes now; only depends on the value of p
otherp <- c(0.01,0.05,0.1,0.2,0.25)
aFreqs <- rep(otherp,each=700)
aFreqs[1:5] <- 0.2

set.seed(6422)
geno <- matrix(NA,nrow=3500,ncol=2100)
for(i in 1:100){
  geno[,(i*16-15):(i*16)] <- 2*Family_alleles_NmarkerX(aFreqs,3500,SEX)
}

for(i in 1:250){
  geno[,(1600+i)] <- 2*INDIV_alleles_Nmarker(aFreqs,3500) 
  geno[,(1600+250+i)] <- 2*INDIV_alleles_NmarkerX(aFreqs,3500)
} # cols are samples, rows are SNPs

apply(geno[,1600:2100],2,table) # make sure the males at least have no 2's; good
dim(geno) # 3500 2100, so snp x sample -- 2000 SNPs and 2100 samples
sex <- c(sex,rep("F",250),rep("M",250)) 
table(sex) # good - 1050 of each

##
# run ibd on these genotypes -- females only
# doesn't work on the entire sample since we need x chr specific code...

# sid <- 1:length(sex)
# fems <- sid[sex=="F"]
# library(SNPRelate)
# snpgdsCreateGeno("test.gds",genmat=geno,snpfirstdim=TRUE,
#                  sample.id=1:ncol(geno),snp.id=1:nrow(geno),
#                  snp.chromosome=rep(1,nrow(geno)),snp.position=rep(0,nrow(geno)))
# (genofile <- openfn.gds("test.gds"))
# closefn.gds(genofile)
# 
# gdsobj <- snpgdsOpen("test.gds")
# ibd <- snpgdsIBDMoM(gdsobj,sample.id=fems)
# ibd.coeff <- snpgdsIBDSelection(ibd)
# closefn.gds(gdsobj)

#pdf("tmp_ibd_results.pdf")
#plot(ibd.coeff$k0,ibd.coeff$k1)
#dev.off()

#pdf("tmp_ibd_resultsKC.pdf")
#hist(ibd.coeff$kinship)
#dev.off()

## use my code and see what happens
# check to be sure that the simulated genotypes are yielding the patterns of relatedness we expect
kinshipCalc <- matrix(NA,nrow=2100,ncol=2100)
for(i in 1:ncol(geno)){
  for(j in i:ncol(geno)){
    kinshipCalc[i,j] <- getKinX(geno[,i],geno[,j],sex[i],sex[j],aFreqs)
  }
}
summary(diag(kinshipCalc))

# ok, looks much better now!!!
# try the simulations again with this method of simulating the genotypes


library(grDevices)
pdf("test_kc.pdf")
diffMat <- trueKinXn[lower.tri(trueKinXn,diag=FALSE)]-kinshipCalc[lower.tri(kinshipCalc,diag=FALSE)]
x <- trueKinXn[lower.tri(trueKinXn,diag=FALSE)]
plot(jitter(x), diffMat,main="Difference between Simulated Geno KC observed and Expected KC\nFor X Chromosome SNPs",
     xlab="Expected KC",ylab="Expected KC - Estimated KC",
     col=adjustcolor("purple",alpha=0.6),pch=19)
abline(h=0)
dev.off()







## check that KC matrices are correct

library(Matrix)
tmp <- as.matrix(bdiag(kinFull))
diag(tmp) <- 1
allequal(kinFull,tmp) # TRUE
tmp <- as.matrix(bdiag(kinFullX))
diag(tmp) <- 1
allequal(kinFullX,tmp) # TRUE


