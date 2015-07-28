assocTestMixedModel <- function(genoData,
                                snpStart,
                                snpEnd,
                                cholSigmaInv,
                                Phi = NULL,
                                outcome,
                                covar.vec = NULL,
                                ivar.vec = NULL,
                                ivar.SE = "model",
                                ivar.return.betas = FALSE,
                                scan.exclude = NULL,
                                method = "EMMAX",
                                impute.geno = TRUE,
                                block.size = 5000,
                                verbose = TRUE){
  
  
  # set which samples to keep
  scanID <- getScanID(genoData)
  keep <- rep(TRUE, nscan(genoData))
    keep <- keep & !(scanID %in% scan.exclude)
  
  # get chromosome information
  chr <- getChromosome(genoData, index=snpStart:snpEnd)
         
    # remove samples with any missing data
    dat <- as.data.frame(dat[keep,]) ## dat holds outcome, covariates for all samples
    # outcome
    Y <- dat[,outcome] ## Y is vector outcome
    
    # create design matrix
    model.formula <- as.formula(paste(paste(outcome,"~"), paste(covar.vec,collapse="+")))
    W <- model.matrix(model.formula, data=dat)    ## W is a matrix to use for regression, includes intercept & all covars
    
k <- ncol(W)
  n <- length(scanID)
 
  # number of SNPs in the segment
  nsnp.seg <- snpEnd - snpStart + 1
  # determine number of SNP blocks
  nblocks <- ceiling(nsnp.seg/block.size)
  
  # loop through blocks
  for(b in 1:nblocks){
       
    # impute missing genotype values
    # impute to frequency value
    
    keep.geno <- check == 0
    
    # sample size
    n <- sum(keep.geno)
    
      # W: covariate matrix
      # C: cholesky decomposition of sigma inverse
      # sigma inverse: inverse phenotype covariance matrix
        CW <- crossprod(C.block, W.block)
        Mt <- C - CC''W (C'' W)-1''W''C 
        Mt <- C.block - tcrossprod(tcrossprod(C.block,tcrossprod(chol2inv(chol(crossprod(CW))),CW)),CW) # this is a matrix used to adjust phenotype and genotype for fixed effect covariates AND "decorrelating" the phenotype and genotype -- ie adjusting out covariance structure given by sigma matrix.
        Ytilde <- Mt''Y # so ytilde is the phenotype adjusted for the covariates/correlation structure
        sY2 <- sum(Ytilde^2)    
        
    # perform regressions
        Xtilde <- crossprod(Mt,geno) # adjust genotypes for correlation structure and fixed effects
        XtX <- colSums(Xtilde^2) # vector of X^T SigmaInv X (for each SNP)
        # filter monomorphic SNPs
        XtX[which(maf==0)] <- NA
        beta <- as.vector(crossprod(Xtilde,Ytilde)/XtX)
        Vbeta <- (sY2/XtX - beta^2)/(n - k - 1) # RSS/XtX
        Stat <- beta^2/Vbeta 
    
    # collect results
      res[bidx,"Est"] <- beta
        res[bidx,"SE"] <- sqrt(Vbeta)
        res[bidx,"Stat"] <- Stat
        res[bidx,"pval"] <- pchisq(Stat, df=1, lower.tail=FALSE)        
        
    endTime <- Sys.time()
    rate <- format(endTime - startTime, digits=4)
    
    keep.previous <- keep.geno
    
    if(verbose) message(paste("Block", b, "of", nblocks, "Completed -", rate))
  } # end block loop
  
  # results data frame
  res <- as.data.frame(res)
  
  # add in snpID
  res$snpID <- getSnpID(genoData, snpStart:snpEnd)
  
  # convert minor.allele coding back to A/B
  res[,"minor.allele"][res[,"minor.allele"] == 1] <- "A"
  res[,"minor.allele"][res[,"minor.allele"] == 0] <- "B"
  
  if(ivar.return.betas){
    res.betas <- as.data.frame(res.betas)
    names(res.Vbetas) <- res$snpID
    res <- list(cbind(res, res.betas), res.Vbetas)
  }  
  
  return(res)
}


.subsetCholSigmaInv <- function(cholSigmaInv, chol.idx) {
  if(length(chol.idx) > 0){
    # subset cholSigmaInv
    SigmaInv <- tcrossprod(cholSigmaInv)
    for(i in sort(chol.idx, decreasing=TRUE)){
      SigmaInv <- SigmaInv[-i,-i] - tcrossprod(SigmaInv[-i,i])/SigmaInv[i,i]
    }
    cholSigmaInv <- t(chol(SigmaInv))
  }
  
  cholSigmaInv
}
