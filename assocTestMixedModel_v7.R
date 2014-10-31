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
  
  # check that method is valid
  if(!is.element(method,c("EMMAX","MASTOR"))){
    stop("method must be one of EMMAX or MASTOR")
  }
  
  if(method == "MASTOR" & is.null(Phi)){
    stop("MASTOR method requires specification of Phi")
  }
  
  if(!is.null(ivar.vec)){
    if(method != "EMMAX"){
      stop("EMMAX method must be used for Tests with Genotype Interaction")
    }
    if(!is.element(ivar.SE, c("model","robust"))){
      stop("ivar.SE must specify model or robust standard errors")
    }
    if(ivar.SE == "robust"){
      stop("robust SE not operational yet!")
    }
  }
  
  # set which samples to keep
  scanID <- getScanID(genoData)
  keep <- rep(TRUE, nscan(genoData))
  
  # samples excluded from entire analysis
  if(!is.null(scan.exclude)){
    keep <- keep & !(scanID %in% scan.exclude)
  }
  
  # get chromosome information
  chr <- getChromosome(genoData, index=snpStart:snpEnd)
  
  # X chromosome check for sex variable
  if(XchromCode(genoData) %in% chr & !hasSex(genoData)){
    stop("Sex values for the samples are required to compute MAF for X chromosome SNPs")
  }
  
  # Y chromosome 
  if(YchromCode(genoData) %in% chr){
    # check for sex variable
    if(!hasSex(genoData)){
      stop("Sex values for the samples are required for Y chromosome SNPs")
    }
    if(!all(chr == YchromCode(genoData))){
      stop("Y chromosome must be analyzed separately")
    }
    # only keep males
    keep <- keep & (getSex(genoData) == "M")
  }
  
#   # filter individuals with anomalies if not imputing genotypes
#   if(!impute.geno){
#     # read in first SNP to get missingness for segment
#     geno <- getGenotype(genoData, snp=c(snpStart, 1), scan=c(1,-1))
#     keep <- keep & !(is.na(geno))
#   }
  
  # read in outcome and covariate data
  if(verbose) message("Reading in Phenotype and Covariate Data...")
  if(!is.null(covar.vec)){
    cvnames <- unique(unlist(strsplit(covar.vec,"[*:]")))
    # read in data
    dat <- as.data.frame(getScanVariable(genoData, c(outcome,cvnames)))
    # identify samples with any missing data
    keep <- keep & apply(dat,1,function(x){ all(!is.na(x)) })
    
    # read in interaction variable data
    if(!is.null(ivar.vec)){
      ivnames <- unique(unlist(strsplit(ivar.vec,"[*:]")))
      # read in data
      idat <- as.data.frame(getScanVariable(genoData, c(outcome,ivnames)))
      # identify samples with any missing data
      keep <- keep & apply(idat,1,function(x){ all(!is.na(x)) })
      # remove samples with any missing data
      idat <- as.data.frame(idat[keep,])
      # create matrix to store interaction variables
      iformula <- as.formula(paste(paste(outcome,"~"), paste(ivar.vec,collapse="+")))
      V <- model.matrix(iformula, data=idat)
      v <- ncol(V)
    }
    
    # remove samples with any missing data
    dat <- as.data.frame(dat[keep,])
    # outcome
    Y <- dat[,outcome]
    # create design matrix
    model.formula <- as.formula(paste(paste(outcome,"~"), paste(covar.vec,collapse="+")))
    W <- model.matrix(model.formula, data=dat)    
    
  }else{
    # read in data
    dat <- getScanVariable(genoData,outcome)
    # identify samples with any missing data
    keep <- keep & !is.na(dat)
    
    # read in interaction variable data
    if(!is.null(ivar.vec)){
      ivnames <- unique(unlist(strsplit(ivar.vec,"[*:]")))
      # read in data
      idat <- as.data.frame(getScanVariable(genoData, c(outcome,ivnames)))
      # identify samples with any missing data
      keep <- keep & apply(idat,1,function(x){ all(!is.na(x)) })
      # remove samples with any missing data
      idat <- as.data.frame(idat[keep,])
      # create matrix to store interaction variables
      iformula <- as.formula(paste(paste(outcome,"~"), paste(ivar.vec,collapse="+")))
      V <- model.matrix(iformula, data=idat)
      v <- ncol(V)
    }
    
    # outcome
    Y <- dat[keep]
    # design matrix
    W <- matrix(1,nrow=length(Y),ncol=1)
  }
  k <- ncol(W)
  scanID <- scanID[keep]

  # which samples to remove from cholSigmaInv
  if(!all(scanID %in% colnames(cholSigmaInv))){
    stop("All of the included Samples must be in the cholSigmaInv matrix")
  }
  chol.idx <- which(!(colnames(cholSigmaInv) %in% scanID))
  cholSigmaInv <- .subsetCholSigmaInv(cholSigmaInv, chol.idx)

  
  if(method == "MASTOR"){
    # subset Phi matrix
    if(!all(scanID %in% colnames(Phi))){
      stop("All of the included Samples must be in the Phi matrix")
    }
    keepPhi <- colnames(Phi) %in% scanID
    Phi <- Phi[keepPhi,keepPhi]
    
    # check that Phi and SigmaInv match
    if(!all(colnames(Phi) == colnames(cholSigmaInv))){
      stop("The Phi and cholSigmaInv matrices do not match")
    }
  }
  ### THIS IS WHERE THE BIG CHANGES START

  # sample size, assuming no missing genotypes
  n <- length(scanID)
  if(verbose) message("Running analysis with ", n, " Samples")

  
  # number of SNPs in the segment
  nsnp.seg <- snpEnd - snpStart + 1
  # determine number of SNP blocks
  nblocks <- ceiling(nsnp.seg/block.size)
  
  # set up results matrix
  if(method == "MASTOR"){
    nv <- c("snpID","chr","n","MAF","minor.allele","Stat","pval")
    
  }else if(method == "EMMAX"){
    nv <- c("snpID","chr","n","MAF","minor.allele")
    if(is.null(ivar.vec)){
      nv <- append(nv, c("Est","SE","Stat","pval"))
    }else{
      nv <- append(nv, c("GxE.Stat","GxE.pval","Joint.Stat","Joint.pval"))      
    }
  }
  res <- matrix(NA, nrow=nsnp.seg, ncol=length(nv), dimnames=list(NULL, nv))
  
  # hold interaction betas
  if(ivar.return.betas){
    res.betas <- matrix(NA, nrow=nsnp.seg, ncol=v, dimnames=list(NULL, c("Est.G", paste("Est.G",colnames(V)[-1],sep=":"))))
    res.Vbetas <- vector("list",nsnp.seg)
  }
  
  
  # chromosome
  res[,"chr"] <- chr
  
  # since we haven't done any loops here yet:
  keep.previous <- rep(TRUE, n)
  
  if(verbose) message("Beginning Calculations...")
  # loop through blocks
  for(b in 1:nblocks){
    
    # keep track of time for rate reporting
    startTime <- Sys.time()
    
    snp.start.pos <- snpStart + (b-1)*block.size
    nsnp.block <- block.size
    if(snp.start.pos + nsnp.block > snpEnd){
      nsnp.block <- snpEnd - snp.start.pos + 1
    }
    snp.end.pos <- snp.start.pos + nsnp.block - 1
    
    bidx <- ((b-1)*block.size+1):((b-1)*block.size+nsnp.block)
    
    # get genotypes for the block
    geno <- getGenotype(genoData, snp=c(snp.start.pos, nsnp.block), scan=c(1,-1), transpose=TRUE)
    
    ####### ADDING THIS LINE #########
    geno <- t(geno)
    ####### ADDING THIS LINE #########
    
    # subset
    geno <- as.matrix(geno)[keep, , drop=F]
    
    # allele frequency
    freq <- 0.5*colMeans(geno, na.rm=T)
    # for X chr
    if(XchromCode(genoData) %in% chr[bidx]){
      # which are on X chr
      Xidx <- chr[bidx] == XchromCode(genoData)
      # males
      m <- (getSex(genoData) == "M")[keep]
      f <- (getSex(genoData) == "F")[keep]
      # calculate allele freq for X
      freq[Xidx] <- (0.5 * colSums(geno[m,Xidx], na.rm=T) + colSums(geno[f,Xidx], na.rm=T)) / 
        (colSums(!is.na(geno[m, Xidx])) + 2*colSums(!is.na(geno[f, Xidx])))
        
    }
    
    # MAF
    maf <- ifelse(freq < 0.5, freq, 1-freq)
    res[bidx,"MAF"] <- maf
    # minor allele coding:  A = 1, B = 0
    res[bidx,"minor.allele"] <- ifelse(freq < 0.5, 1, 0)
    
    
    # impute missing genotype values
    if(impute.geno){      
      miss.idx <- which(is.na(geno))
      if(length(miss.idx) > 0){
        snp.idx <- ceiling(miss.idx/n)
        geno[miss.idx] <- 2*freq[snp.idx]
      }
    }
    
    # check for missingness
    check <- rowSums(is.na(geno))    
    if (!all(check %in% c(0, nsnp.block))) {
      stop("sporadic missing in block size > 1")
    }
    
    keep.geno <- check == 0
    
    # get rid of missing for this block
    # get rid of missing for this block
    # can probably put this in an if statement
    geno <- as.matrix(geno)[keep.geno, , drop=F]
    
    # sample size
    n <- sum(keep.geno)
    res[bidx, "n"] <- n
    
    if (b == 1 | !all(keep.previous == keep.geno)){
      # calculate matrices
      
      if (b > 1) warning("recalculating matrices!")
      #if(verbose) message("Pre-Computing some Matrices for Analysis...")
      
      
      # subsetting
      W.block <- W[keep.geno, , drop=F]
      Y.block <- Y[keep.geno]
      # check if interactions
      if (!is.null(ivar.vec)){
        V.block <- V[keep.geno, , drop=F]
      }
      
      # here we have to subset the matrix, not the inverse
      # this is a fancy way of getting the inverse of the subset without having to get the original matrix
      chol.idx <- which(!(colnames(cholSigmaInv) %in% scanID[keep.geno]))
      C.block <- .subsetCholSigmaInv(cholSigmaInv, chol.idx)
      if (method == "MASTOR"){
        Phi.block <- Phi[keep.geno, keep.geno]
      }
      
      # W: covariate matrix
      # C: cholesky decomposition of sigma inverse
      # sigma inverse: inverse phenotype covariance matrix
      if(method == "EMMAX"){
        #k <- ncol(W.block) # number of covariates
        CW <- crossprod(C.block, W.block)
        Mt <- C.block - tcrossprod(tcrossprod(C.block,tcrossprod(chol2inv(chol(crossprod(CW))),CW)),CW) # this is a matrix used to adjust phenotype and genotype for fixed effect covariates AND "decorrelating" the phenotype and genotype -- ie adjusting out covariance structure given by sigma matrix.
        Ytilde <- crossprod(Mt,Y.block) # so ytilde is the phenotype adjusted for the covariates/correlation structure
        sY2 <- sum(Ytilde^2)    
        
      }else if(method == "MASTOR"){
        CW <- crossprod(C.block, W.block)
        Mt <- C.block - tcrossprod(tcrossprod(C.block,tcrossprod(chol2inv(chol(crossprod(CW))),CW)),CW)
        PY <- tcrossprod(Mt, crossprod(Y.block, C.block)) # V = PY    
        PYtPhiPY <- crossprod(PY, crossprod(Phi.block, PY))

        PhiInv <- tryCatch(chol2inv(chol(Phi.block)), error=function(e) TRUE)
        if(is.logical(PhiInv)){
          # switch to 2p(1-p)
          message("Phi matrix not invertable, using 2p(1-p) to estimate variance scalar")
        }else{
          PhiInvW <- crossprod(PhiInv,W.block)
          U <- PhiInv - tcrossprod(tcrossprod(PhiInvW,chol2inv(chol(crossprod(W.block,PhiInvW)))),PhiInvW)
        }
      }
    }
    
        
    # perform regressions
    if(method == "EMMAX"){
      # no interaction
      if(is.null(ivar.vec)){
        Xtilde <- crossprod(Mt,geno) # adjust genotypes for correlation structure and fixed effects
        XtX <- colSums(Xtilde^2) # vector of X^T SigmaInv X (for each SNP)
        # filter monomorphic SNPs
        XtX[which(maf==0)] <- NA
        beta <- as.vector(crossprod(Xtilde,Ytilde)/XtX)
        Vbeta <- (sY2/XtX - beta^2)/(n - k - 1) # RSS/XtX
        Stat <- beta^2/Vbeta
        
        # interaction
      }else{
        GxE.Stat <- rep(NA, ncol(geno))
        Joint.Stat <- rep(NA, ncol(geno))
        for(g in 1:ncol(geno)){
          # filter monomorphic or missing SNPs
          if(maf[g] == 0 || is.na(maf[g])){
            next
          }else{
            Xtilde <- crossprod(Mt,geno[,g]*V.block)
            XtX <- crossprod(Xtilde)
            XtXinv <- tryCatch( chol2inv(chol(XtX)), error=function(e){TRUE})  # this is inverse A matrix of sandwich
            # check that the error function above hasn't been called (which returns TRUE instead of the inverse matrix)
            if(is.logical(XtXinv)){
              next
            }
            XtY <- crossprod(Xtilde, Ytilde)
            betas <- crossprod(XtXinv, XtY)
            
            if(ivar.SE == "model"){
              # model based
              RSS <- as.numeric((sY2 - crossprod(XtY,betas))/(n - k - v))
              Vbetas <- XtXinv*RSS
              GxE.Stat[g] <- tryCatch( crossprod(betas[-1],crossprod(chol2inv(chol(Vbetas[-1,-1])),betas[-1])), error=function(e){NA})
              Joint.Stat[g] <- tryCatch( crossprod(betas,crossprod(XtX,betas))/RSS, error=function(e){NA})
            }else if(ivar.SE == "robust"){
              # sandwich
              R <- as.vector(Ytilde - Xtilde %*% betas)
              B <- crossprod(Xtilde*R)   # B matrix of sandwich
              Vbetas <- crossprod(XtXinv,tcrossprod(B,XtXinv))
              GxE.Stat[g] <- tryCatch( crossprod(betas[-1],crossprod(chol2inv(chol(Vbetas[-1,-1])),betas[-1])), error=function(e){NA})
              Joint.Stat[g] <- tryCatch( crossprod(betas,crossprod(chol2inv(chol(Vbetas)),betas)), error=function(e){NA})
            }            
            
            # save interaction betas?
            if(ivar.return.betas){
              res.betas[bidx[g],] <- betas
              res.Vbetas[[bidx[g]]] <- Vbetas
            }
            
          }
        } # end SNP loop
      } # end interaction
      
    }else if(method == "MASTOR"){
      if(is.logical(PhiInv)){
        sigma2 <- 2*freq*(1-freq)
        #sigma2 <- apply(geno,2,function(x){ sum(x == 1) })/nrow(geno)
        #sigma2 <- (colSums(geno^2) - n*4*freq^2)/(n-1)
      }else{
        sigma2 <- colSums(crossprod(U,geno)*geno)/(n-k)
      }        
      Stat <- crossprod(geno,PY)^2/(sigma2*PYtPhiPY)
      # filter monomorphic SNPs
      Stat[which(maf==0 | sigma2 == 0)] <- NA
    }
    
    # collect results
    if(method == "MASTOR"){
      res[bidx,"Stat"] <- Stat
      # compute p-values
      res[bidx,"pval"] <- pchisq(Stat, df=1, lower.tail=FALSE)
      
    }else if(method == "EMMAX"){
      if(is.null(ivar.vec)){
        res[bidx,"Est"] <- beta
        res[bidx,"SE"] <- sqrt(Vbeta)
        res[bidx,"Stat"] <- Stat
        res[bidx,"pval"] <- pchisq(Stat, df=1, lower.tail=FALSE)        
      }else{
        res[bidx,"GxE.Stat"] <- GxE.Stat
        res[bidx,"GxE.pval"] <- pchisq(GxE.Stat, df=(v-1), lower.tail=FALSE)
        res[bidx,"Joint.Stat"] <- Joint.Stat
        res[bidx,"Joint.pval"] <- pchisq(Joint.Stat, df=v, lower.tail=FALSE)        
      }
    }
    
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
