
# simulate a quantitative phenotype
# simple additive model
# y_j = sum_i w_ij u_i + e_j
# i causal variants
# w_ij = (x_ij-2p_i)/(sqrt(2p_i(1-p_i)))
# u_i = allelic effect of ith causal variant
# e_j = residual effect ~ N(0, var(sum_i w_ij u_i)/(1/h^2)-1)


#sim_phenotype <- function(geno, causalSNPs, effSNPs, herit){
  # geno = matrix of genotypes of size nSamp x SNPs
  # causalSNPs = vector of length num SNPs, true for those which are causal SNPs
  # effSNPs = vector of length num causal SNPs with effect size for each causal SNP
  # herit = a numeric value of user-specified heritability

#  snpFreqFull <- apply(geno,2,function(x){sum(x)/(2*length(x))})
#  tmpG <- geno[,causalSNPs&snpFreqFull!=0&snpFreqFull!=1] # exclude monomorphic SNPs
#  snpFreq <- snpFreqFull[causalSNPs&snpFreqFull!=0&snpFreqFull!=1]

#  w <- (tmpG-2*snpFreq)/(sqrt(2*snpFreq*(1-snpFreq)))
  # now w is matrix of size nSamp x causalSNPs, excluding monomorphic SNPs

#  snpFreqEff <- snpFreqFull[causalSNPs]
#  sumTerm <- w %*% effSNPs[snpFreqEff!=0&snpFreqEff!=1]

#  varEps <- var(sumTerm)/((1/herit^2)-1)
#  eps <- rnorm(length(sumTerm),mean=0,sd=sqrt(varEps))

#  pheno <- sumTerm+eps
#  return(pheno)
#}



simulatePhenotype <- function(kinAutos,kinX,sigmaA,sigmaX,sigmaE,beta1,SNP,seed=NULL){

  if(!is.null(seed)){set.seed(seed)}
  
  n <- nrow(kinAutos)
  ident <- diag(n)
  sigmaMat <- as.matrix(kinAutos*sigmaA+kinX*sigmaX+ident*sigmaE)
  sigmaMat <- make.positive.definite(sigmaMat)
  noise <- mvrnorm(1, mu=rep(0,n), Sigma=sigmaMat)
  
  beta0 <- rep(1,n)
  
  y <- beta0 + SNP*beta1 + noise
  
  return(as.vector(y))
}


getEffSize <- function(p,totalSigma,herit=NULL,beta1=NULL){
  # c <- 2*p*(1-p)
  # c = 2 * allele freq * (1-allele freq) = 2p(1-p)

  # for an x chromosome SNP with half female, half male
  # c = 3 * allele freq * (1-allele freq) = 3p(1-p)
  c <- 3*p*(1-p)
  
  if(is.null(beta1)){
    # calculate beta1 based on herit and c (allele freq)
    res <- sqrt((herit*totalSigma)/((1-herit)*c))
  }
  if(is.null(herit)){
    # calculate herit based on c and effect size
    res <- (beta1^2*c)/(beta1^2*c+totalSigma)
  }
  return(res)
}


