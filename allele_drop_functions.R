
alleledrop_half <-function(x1,x2)
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

alleledrop <-function(x1,x2)
  {
    x3=rep(0,length(x1))

    hetInds <- x1==1
    myval <- runif(sum(hetInds),0,1)
    if(sum(hetInds)>0){
      x3[hetInds] <- x3[hetInds]+as.integer(myval<0.5)
    }
    x3[!hetInds] <- x3[!hetInds]+x1[!hetInds]/2

    hetInds <- x2==1
    myval <- runif(sum(hetInds),0,1)
    if(sum(hetInds)>0){
      x3[hetInds] <- x3[hetInds]+as.integer(myval<0.5)
    }
    x3[!hetInds] <- x3[!hetInds]+x2[!hetInds]/2

    return(x3)
  }


alleledrop_single <- function(x1,x2)
  {
    x3=0
        if(x1==0) {x3=0}
        if(x1==2) {x3=1}
        if(x1==1) {
              myval=runif(1,0,1)
              x3=as.integer(myval<0.5)
            }

        if(x2==0) {x3=x3+0;}
        if(x2==2) {x3=x3+1;}
        if(x2==1){
              myval=runif(1,0,1)
              x3 <- x3 + as.integer(myval<0.5)
            }
        return(x3)
    }

alleledropX<-function(x1,sex1,x2,sex2,sex3)
  {
    x3=rep(0,length(x1))
    
    if(sex1=="F"){ # setting genoF to the vector of genotypes corresponding to mother
      genoF=x1; genoM=x2
    }else{genoF=x2; genoM=x1}
    
    hetInds <- genoF==1
    myval <- runif(sum(hetInds),0,1)
    
    if(sex3=="M") # if kid is male
      {
        # for SNPs with genoF=1, draw runif and choose
        if(sum(hetInds)>0){
          x3[hetInds]=x3[hetInds]+2*as.integer(myval<0.5)
        }
        # for SNPs with genoF=0,2, transmit that genotype to kid
        x3[!hetInds] <- x3[!hetInds]+genoF[!hetInds]
      }else{ # kid is female
        if(sum(hetInds)>0){
          x3[hetInds] <- x3[hetInds]+as.integer(myval<0.5)
        }
        x3[!hetInds] <- x3[!hetInds]+genoF[!hetInds]/2

        x3=x3+genoM/2
      }
  return(x3)
  }

alleledropX_single<-function(x1,sex1,x2,sex2,sex3)
{
  x3=0
  if(sex1=="F"){
    genoF=x1; genoM=x2
  }else{genoF=x2; genoM=x1}
  
  if(sex3=="M")
  {
    if(genoF==1){
      myval=runif(1,0,1)
      x3=x3+2*as.integer(myval<0.5)
    }else{ x3=x3+genoF }
  } # end if sex3==M loop
  
  if(sex3=="F")
  {
    if(genoF==1){
      myval=runif(1,0,1)
      x3=x3+as.integer(myval<0.5)
    }else{ x3=x3+genoF/2 }
    
    x3=x3+genoM/2  
  } # end if sex3==F loop
  
  return(x3)
}


alleledropX_half <- function(x1,sex1,x2,sex2,sex3)
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

  val1 <- FREQ1^2
#  val2 <- val1+2*FREQ1*(1-FREQ1)

  GENO[myval<val1] <- 2
  GENO[myval<(val1+2*FREQ1*(1-FREQ1))&myval>=val1] <- 1
  
  return(GENO)
}

INDIV_alleles_Nmarker_half <- function(FREQ1,nloci)
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


INDIV_alleles_NmarkerX_half <- function(FREQ1,nloci)
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



INDIV_alleles_NmarkerX<-function(FREQ1,nloci)
{
  GENO= rep(0,nloci);

  randVals <- runif(nloci,0,1)
  GENO[randVals<FREQ1] <- 2
  
  return(GENO)
}


#####

Family_alleles_Nmarker<-function(FREQ1,nloci)
{
  GENO= matrix(0,nloci,16);
  
  founder1=1
  founder2=2
  founder3=3
  founder4=4
  founder5=5
  
  for(count in 1:5){GENO[,count]=INDIV_alleles_Nmarker(FREQ1,nloci)}
  
  count=6
  
  sib1=count
  for(i in 1:nloci)
  {GENO[i,count]=alleledrop(GENO[i,founder1],GENO[i,founder2])}
  
  count=count+1
  sib2=count
  for(i in 1:nloci)
  {GENO[i,count]=alleledrop(GENO[i,founder1],GENO[i,founder2])}
  
  count=count+1
  sib3=count
  for(i in 1:nloci)
  {GENO[i,count]=alleledrop(GENO[i,founder1],GENO[i,founder2])}
  
  #SECOND GENERATION FAM1
  
  count=count+1
  for(i in 1:nloci)
  {GENO[i,count]=alleledrop(GENO[i,sib1],GENO[i,founder3])}
  
  count=count+1
  for(i in 1:nloci)
  {GENO[i,count]=alleledrop(GENO[i,sib1],GENO[i,founder3])}
  
  count=count+1
  for(i in 1:nloci)
  {GENO[i,count]=alleledrop(GENO[i,sib1],GENO[i,founder3])}
  
  
  #SECOND GENERATION FAM2
  
  count=count+1
  for(i in 1:nloci)
  {GENO[i,count]=alleledrop(GENO[i,sib2],GENO[i,founder4])}
  
  count=count+1
  for(i in 1:nloci)
  {GENO[i,count]=alleledrop(GENO[i,sib2],GENO[i,founder4])}
  
  
  #SECOND GENERATION FAM3
  
  count=count+1
  for(i in 1:nloci)
  {GENO[i,count]=alleledrop(GENO[i,sib3],GENO[i,founder5])}
  
  count=count+1
  for(i in 1:nloci)
  {GENO[i,count]=alleledrop(GENO[i,sib3],GENO[i,founder5])}
  
  count=count+1
  for(i in 1:nloci)
  {GENO[i,count]=alleledrop(GENO[i,sib3],GENO[i,founder5])}
  
  return(GENO)
  
}




Family_alleles_NmarkerX<-function(FREQ1,nloci,SEX)
{
  GENO= matrix(0,nloci,16);
  
  GENO[,1]=INDIV_alleles_Nmarker(FREQ1,nloci)
  GENO[,3]=INDIV_alleles_Nmarker(FREQ1,nloci)
  GENO[,2]=INDIV_alleles_NmarkerX(FREQ1,nloci)
  GENO[,4]=INDIV_alleles_NmarkerX(FREQ1,nloci)
  GENO[,5]=INDIV_alleles_NmarkerX(FREQ1,nloci)

#  for(i in 1:5)
#  {
#    if(SEX[i]=="F")
#    {GENO[,i]=INDIV_alleles_Nmarker(FREQ1,nloci)}
#    if(SEX[i]=="M")
#    {GENO[,i]=INDIV_alleles_NmarkerX(FREQ1,nloci)}
#  }
  
  founder1=1
  founder2=2
  founder3=3
  founder4=4
  founder5=5

  sib1=6
  sib2=7
  sib3=8
  
  count=6
  GENO[,count]=alleledropX(GENO[,founder1],SEX[founder1],GENO[,founder2],SEX[founder2],SEX[sib1])
  
#  for(i in 1:nloci)
#  {GENO[i,count]=alleledropX(GENO[i,founder1],SEX[founder1],GENO[i,founder2],SEX[founder2],SEX[sib1])}
  
  count=count+1
  GENO[,count]=alleledropX(GENO[,founder1],SEX[founder1],GENO[,founder2],SEX[founder2],SEX[sib2])
#  for(i in 1:nloci)
#  {GENO[i,count]=alleledropX(GENO[i,founder1],SEX[founder1],GENO[i,founder2],SEX[founder2],SEX[sib2])}
  
  count=count+1
  GENO[,count]=alleledropX(GENO[,founder1],SEX[founder1],GENO[,founder2],SEX[founder2],SEX[sib3])
#  for(i in 1:nloci)
#  {GENO[i,count]=alleledropX(GENO[i,founder1],SEX[founder1],GENO[i,founder2],SEX[founder2],SEX[sib3])}

  
  #SECOND GENERATION FAM1
  count=count+1
  GENO[,count]=alleledropX(GENO[,sib1],SEX[sib1],GENO[,founder3],SEX[founder3],SEX[count])
#  for(i in 1:nloci)
#  {GENO[i,count]=alleledropX(GENO[i,sib1],SEX[sib1],GENO[i,founder3],SEX[founder3],SEX[count])}
  
  count=count+1
  GENO[,count]=alleledropX(GENO[,sib1],SEX[sib1],GENO[,founder3],SEX[founder3],SEX[count])
#  for(i in 1:nloci)
#  {GENO[i,count]=alleledropX(GENO[i,sib1],SEX[sib1],GENO[i,founder3],SEX[founder3],SEX[count])}
  
  count=count+1
  GENO[,count]=alleledropX(GENO[,sib1],SEX[sib1],GENO[,founder3],SEX[founder3],SEX[count])
#  for(i in 1:nloci)
#  {GENO[i,count]=alleledropX(GENO[i,sib1],SEX[sib1],GENO[i,founder3],SEX[founder3],SEX[count])}
  
  #SECOND GENERATION FAM2
  count=count+1
  GENO[,count]=alleledropX(GENO[,sib2],SEX[sib2],GENO[,founder4],SEX[founder4],SEX[count])
#  for(i in 1:nloci)
#  {GENO[i,count]=alleledropX(GENO[i,sib2],SEX[sib2],GENO[i,founder4],SEX[founder4],SEX[count])}
  
  count=count+1
  GENO[,count]=alleledropX(GENO[,sib2],SEX[sib2],GENO[,founder4],SEX[founder4],SEX[count])
#  for(i in 1:nloci)
#  {GENO[i,count]=alleledropX(GENO[i,sib2],SEX[sib2],GENO[i,founder4],SEX[founder4],SEX[count])}
  
  
  #SECOND GENERATION FAM3
  count=count+1
  GENO[,count]=alleledropX(GENO[,sib3],SEX[sib3],GENO[,founder5],SEX[founder5],SEX[count])
#  for(i in 1:nloci)
#  {GENO[i,count]=alleledropX(GENO[i,sib3],SEX[sib3],GENO[i,founder5],SEX[founder5],SEX[count])}
  
  count=count+1
  GENO[,count]=alleledropX(GENO[,sib3],SEX[sib3],GENO[,founder5],SEX[founder5],SEX[count])
#  for(i in 1:nloci)
#  {GENO[i,count]=alleledropX(GENO[i,sib3],SEX[sib3],GENO[i,founder5],SEX[founder5],SEX[count])}
  
  count=count+1
  GENO[,count]=alleledropX(GENO[,sib3],SEX[sib3],GENO[,founder5],SEX[founder5],SEX[count])
#  for(i in 1:nloci)
#  {GENO[i,count]=alleledropX(GENO[i,sib3],SEX[sib3],GENO[i,founder5],SEX[founder5],SEX[count])}
    
  return(GENO)
  
}

#####
# function to estimate x chromosome kinship between two individs
getKinX_half<- function(geno1,geno2,sex1,sex2,afreq){
  # takes as arguments the vector of genotypes for individ1, individ2
  # takes the sex for individ1, individ2 as "M" and "F"
  # takes a vector of allele frequencies for all genotypes 
  
  denom = afreq*(1-afreq)
  incl = afreq>0 & afreq<1
  
  # subset everything to exclude monomorphic SNPs
  geno1 = geno1[incl]
  geno2 = geno2[incl]
  afreq = afreq[incl]
  denom = denom[incl]
  
  if(sex1=="F" & sex2=="F"){
    # basic autosomal case
    a = sum((geno1-2*afreq)*(geno2-2*afreq)/(2*denom))
  }
  else if(sex1=="M" & sex2=="F"){
    a = sum((geno1-afreq)*(geno2-2*afreq)/(sqrt(2)*denom))
  }
  else if(sex1=="F" & sex2=="M"){
    a = sum((geno1-2*afreq)*(geno2-afreq)/(sqrt(2)*denom))    
  }else{
    a = sum((geno1-afreq)*(geno2-afreq))/sum(denom)
  }
  return(a/length(geno1))
}

getKinX <- function(geno1,geno2,sex1,sex2,afreq){
    # takes as arguments the vector of genotypes for individ1, individ2
    # takes the sex for individ1, individ2 as "M" and "F"
    # takes a vector of allele frequencies for all genotypes

  return(getKin(geno1,geno2,afreq))
}



#####
# function to estimate autosomal kinship between two individs
getKin <- function(geno1,geno2,afreq){
  # takes as arguments the vector of genotypes for individ1, individ2
  # takes a vector of allele frequencies for all genotypes 
  
  denom = afreq*(1-afreq)
  incl = afreq>0.05 & afreq<0.95
  
  # basic autosomal case
  a = sum((geno1[incl]-2*afreq[incl])*(geno2[incl]-2*afreq[incl]))/sum((2*denom[incl]))
  return(a/length(geno1[incl]))
}

