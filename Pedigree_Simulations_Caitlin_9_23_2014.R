

alleledrop<-function(x1,x2)
{
  x3=0
  if(x1==0) {x3=x3+0;}
  if(x1==1) {x3=x3+.5;}
  if(x1==.5) { 
    myval=runif(1,0,1)
    if(myval<.5) {x3=x3+.5
    }else {x3=x3+0}
  }

  if(x2==0) {x3=x3+0;}
  if(x2==1) {x3=x3+.5;}
  if(x2==.5){ 
    myval=runif(1,0,1)
    if(myval<.5){x3=x3+.5
    } else {x3=x3+0}
  }
  return(x3)
}


alleledropX<-function(x1,sex1,x2,sex2,sex3)
{
  x3=0
  if(sex3=="M")
  {
  	if(sex1=="F"){genoF=x1}
    if(sex2=="F"){genoF=x2}
    
    if(genoF==0){x3=x3+0}
    if(genoF==1){x3=x3+0.5}
    if(genoF==0.5){
      myval=runif(1,0,1)
      if(myval<0.5){
        x3=x3+0.5
      }else{x3=x3+0}
    }
  } # end if sex3==M loop

  if(sex3=="F")
  {
	  if(sex1=="F" & sex2=="M")
    { genoF=x1
      genoM=x2
	  }
    if(sex1=="M" & sex2=="F")
    { genoF=x2
      genoM=x1
    }
    
    if(genoF==0){x3=x3+0}
    if(genoF==1){x3=x3+0.5}
    if(genoF==0.5){
      myval=runif(1,0,1)
      if(myval<0.5){ x3=x3+0.5
      }else{ x3=x3+0 }
    }
    
    x3=x3+genoM  
  } # end if sex3==F loop

return(x3)
}

#####

INDIV_alleles_Nmarker<-function(FREQ1,nloci)
{
  GENO= rep(0,nloci)
  myval=runif(nloci,0,1)
  for(i in 1:nloci)
  {
    val1=(FREQ1[i])^(2)
    val2=val1+2*FREQ1[i]*(1-FREQ1[i])
    val3=val2+(1-FREQ1[i])^(2)

    if(myval[i]<val1)
      {GENO[i]=1
    } else if(myval[i]<val2) {
      GENO[i]=.5
    } else {GENO[i]=0}
  }

  return(GENO)
}


INDIV_alleles_NmarkerX<-function(FREQ1,nloci)
{
  GENO= rep(0,nloci);

  for(i in 1:nloci)
  {
    val1=(FREQ1[i])
    myval=runif(1,0,1)

    if(myval<val1)
      {GENO[i]=.5} else {GENO[i]=0}
  }

  return(GENO)
}


#####

Family_alleles_Nmarker<-function(FREQ1,nloci)
{
  GENO= matrix(0,16,nloci);

  founder1=1
  founder2=2
  founder3=3
  founder4=4
  founder5=5

  for(count in 1:5)
    {GENO[count,]=INDIV_alleles_Nmarker(FREQ1,nloci)}
  
  count=5

  count=count+1
  sib1=count
  for(i in 1:nloci)
    {GENO[count,i]=alleledrop(GENO[founder1,i],GENO[founder2,i])}

  count=count+1
  sib2=count
  for(i in 1:nloci)
    {GENO[count,i]=alleledrop(GENO[founder1,i],GENO[founder2,i])}

  count=count+1
  sib3=count
  for(i in 1:nloci)
    {GENO[count,i]=alleledrop(GENO[founder1,i],GENO[founder2,i])}

#SECOND GENERATION FAM1

count=count+1
for(i in 1:nloci)
{GENO[count,i]=alleledrop(GENO[sib1,i],GENO[founder3,i])}

count=count+1
for(i in 1:nloci)
{GENO[count,i]=alleledrop(GENO[sib1,i],GENO[founder3,i])}

count=count+1
for(i in 1:nloci)
{GENO[count,i]=alleledrop(GENO[sib1,i],GENO[founder3,i])}


#SECOND GENERATION FAM2

count=count+1
for(i in 1:nloci)
{GENO[count,i]=alleledrop(GENO[sib2,i],GENO[founder4,i])}

count=count+1
for(i in 1:nloci)
{GENO[count,i]=alleledrop(GENO[sib2,i],GENO[founder4,i])}


#SECOND GENERATION FAM3

count=count+1
for(i in 1:nloci)
{GENO[count,i]=alleledrop(GENO[sib3,i],GENO[founder5,i])}

count=count+1
for(i in 1:nloci)
{GENO[count,i]=alleledrop(GENO[sib3,i],GENO[founder5,i])}



count=count+1
for(i in 1:nloci)
{GENO[count,i]=alleledrop(GENO[sib3,i],GENO[founder5,i])}

return(GENO)

}




Family_alleles_NmarkerX<-function(FREQ1,nloci,SEX)
{
	GENO= matrix(0,16,nloci);



for(i in 1:5)
{
	if(SEX[i]=="F")
	{GENO[i,]=INDIV_alleles_Nmarker(FREQ1,nloci)}
if(SEX[i]=="M")
{GENO[i,]=INDIV_alleles_NmarkerX(FREQ1,nloci)}
	}



founder1=1
founder2=2
founder3=3
founder4=4
founder5=5


count=5



count=count+1
 sib1=count
for(i in 1:nloci)
{GENO[count,i]=alleledropX(GENO[founder1,i],SEX[founder1],GENO[founder2,i],SEX[founder2],SEX[sib1])}

count=count+1
 sib2=count
for(i in 1:nloci)
{GENO[count,i]=alleledropX(GENO[founder1,i],SEX[founder1],GENO[founder2,i],SEX[founder2],SEX[sib2])}

count=count+1
 sib3=count
for(i in 1:nloci)
{GENO[count,i]=alleledropX(GENO[founder1,i],SEX[founder1],GENO[founder2,i],SEX[founder2],SEX[sib3])}



#SECOND GENERATION FAM1

count=count+1
for(i in 1:nloci)
{GENO[count,i]=alleledropX(GENO[sib1,i],SEX[sib1],GENO[founder3,i],SEX[founder3],SEX[count])}
	


count=count+1
for(i in 1:nloci)
{GENO[count,i]=alleledropX(GENO[sib1,i],SEX[sib1],GENO[founder3,i],SEX[founder3],SEX[count])}
	


count=count+1
for(i in 1:nloci)
{GENO[count,i]=alleledropX(GENO[sib1,i],SEX[sib1],GENO[founder3,i],SEX[founder3],SEX[count])}
	


#SECOND GENERATION FAM2

count=count+1
for(i in 1:nloci)
{GENO[count,i]=alleledropX(GENO[sib2,i],SEX[sib2],GENO[founder4,i],SEX[founder4],SEX[count])}

count=count+1
for(i in 1:nloci)
{GENO[count,i]=alleledropX(GENO[sib2,i],SEX[sib2],GENO[founder4,i],SEX[founder4],SEX[count])}
	

#SECOND GENERATION FAM3

count=count+1
for(i in 1:nloci)
{GENO[count,i]=alleledropX(GENO[sib3,i],SEX[sib3],GENO[founder5,i],SEX[founder5],SEX[count])}
	
count=count+1
for(i in 1:nloci)
{GENO[count,i]=alleledropX(GENO[sib3,i],SEX[sib3],GENO[founder5,i],SEX[founder5],SEX[count])}

count=count+1
for(i in 1:nloci)
{GENO[count,i]=alleledropX(GENO[sib3,i],SEX[sib3],GENO[founder5,i],SEX[founder5],SEX[count])}
	

return(GENO)

}

#####
# function to estimate x chromosome kinship between two individs
getKin <- function(geno1,geno2,sex1,sex2,afreq){
  # takes as arguments the vector of genotypes for individ1, individ2
  # takes the sex for individ1, individ2 as "M" and "F"
  # takes a vector of allele frequencies for all genotypes 
  
  denom = afreq*(1-afreq)
  if(sex1=="F" & sex2=="F"){
    # basic autosomal case
    a = sum((geno1-2*afreq)*(geno2-2*afreq)/(2*denom))
  }
  if(sex1=="M" & sex2=="F"){
    a = sum((geno1-afreq)*(geno2-2*afreq)/(sqrt(2)*denom))
  }
  if(sex1=="F" & sex2=="M"){
    a = sum((geno1-2*afreq)*(geno2-afreq)/(sqrt(2)*denom))    
  }
  if(sex1=="M" & sex2=="M"){
    a = sum((geno1-afreq)*(geno2-afreq)/(denom))
  }
  return(a/length(geno1))
}


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

# estimate 1000 xchr loci
set.seed(6422)
genoX <- 2*Family_alleles_NmarkerX(rep(f1,1000),1000,SEX)

kinship <- matrix(NA,nrow=16,ncol=16)
for(i in 1:nrow(genoX)){
  for(j in 1:nrow(genoX)){
    kinship[i,j] <- getKin(genoX[i,],genoX[j,],SEX[i],SEX[j],f1)
  }
}
isSymmetric(kinship) # TRUE

summary(diag(kinship))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.9775  0.9883  0.9977  1.0020  1.0100  1.0380 

##
# what do we expect for these individuals?
# read in true kinship coefficients, multiply by 2 to get relatedness coef
tmp <- read.table("pedigree_16individs_output")
trueKinX <- matrix(NA,nrow=16, ncol=16)
trueKinX[lower.tri(trueKinX,diag=TRUE)] <- tmp[,"V4"]
trueKinXn <- trueKinX
trueKinX <- data.frame(trueKinX)
trueKinX$SEX <- SEX
diag(trueKinX) <- 1
trueKinXn[trueKinX$SEX=="F",trueKinX$SEX=="F"] <- 2*trueKinXn[trueKinX$SEX=="F",trueKinX$SEX=="F"]

summary(kinship-trueKinX)

# estimate 10000 xchr loci now
set.seed(6422)
genoX <- 2*Family_alleles_NmarkerX(rep(f1,10000),10000,SEX)

kinship10K <- matrix(NA,nrow=16,ncol=16)
for(i in 1:nrow(genoX)){
  for(j in 1:nrow(genoX)){
    kinship10K[i,j] <- getKin(genoX[i,],genoX[j,],SEX[i],SEX[j],f1)
  }
}

summary(diag(kinship10K))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.9857  0.9918  0.9976  0.9969  0.9996  1.0100 

genoX <- 2*Family_alleles_NmarkerX(rep(f1,100000),100000,SEX)

kinship100K <- matrix(NA,nrow=16,ncol=16)
for(i in 1:nrow(genoX)){
  for(j in 1:nrow(genoX)){
    kinship100K[i,j] <- getKin(genoX[i,],genoX[j,],SEX[i],SEX[j],f1)
  }
}

summary(kinship-trueKinX)
# make qq plot of expected, observed
plot(trueKinXn[lower.tri(trueKinXn,diag=FALSE)],
     -kinship[lower.tri(kinship,diag=FALSE)]+trueKinXn[lower.tri(trueKinXn,diag=FALSE)],xlab="Expected",
     ylab="Expected - Observed",
     xlim=c(0,1),ylim=c(-1,1),col="purple",pch=19)
points(trueKinXn[lower.tri(trueKinXn,diag=FALSE)],
       -kinship10K[lower.tri(kinship10K,diag=FALSE)]+trueKinXn[lower.tri(trueKinXn,diag=FALSE)],col="red",pch=19)
points(trueKinXn[lower.tri(trueKinXn,diag=FALSE)],
       -kinship100K[lower.tri(kinship100K,diag=FALSE)]+trueKinXn[lower.tri(trueKinXn,diag=FALSE)],col="cyan",pch=19)
abline(h=0)

####
# need to get kinship matrices
# need to simulate quantitative phenotype
# need to call assoc test with kinship matrices and phenotype, get type I error

# need to estimate kinship matrices
# need to call assoc test with estimated kinship matrices and phenotype, get type I error

#####
library(kinship2)

ped <- read.table("pedigree_16individs.txt")
ped$V5[ped$V5==0] <- 2
pedPlt <- pedigree(ped$V2,ped$V3,ped$V4,ped$V5)
pdf("pedigree_16individs.pdf")
plot(pedPlt,cex=2,lwd=2)
dev.off()
