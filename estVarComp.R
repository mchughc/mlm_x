estVarComp <- function(scanAnnot,
                       covMatList,
                       outcome,
                       covar.vec = NULL,
                       scan.exclude = NULL,
                       start = NULL,
                       method = "AIREML",
                       AIREML.tol = 1e-6,
                       ngrid = 100,
                       llim = -10,
                       ulim = 10,
                       grid.tol = 1e-10,
                       verbose = TRUE){

  # set which samples to keep
  scanID <- getVariable(scanAnnot, "scanID")
  keep <- rep(TRUE, length(scanID))

  # samples excluded from entire analysis
  if(!is.null(scan.exclude)){
    keep <- keep & !(scanID %in% scan.exclude)
  }

  # read in outcome and covariate data
  if(verbose) message("Reading in Phenotype and Covariate Data...")
  if(!is.null(covar.vec)){
    cvnames <- unique(unlist(strsplit(covar.vec,"[*:]")))
    # read in data
    dat <- as.data.frame(getVariable(scanAnnot, c(outcome,cvnames)))
    # remove samples with any missing data
    keep <- keep & apply(dat,1,function(x){ all(!is.na(x)) })
    dat <- as.data.frame(dat[keep,])
    # outcome
    Y <- dat[,outcome]
    # create design matrix
    model.formula <- as.formula(paste(paste(outcome,"~"), paste(covar.vec,collapse="+")))
    W <- model.matrix(model.formula, data=dat)
  }else{
    # read in data
    dat <- getVariable(scanAnnot,outcome)
    # remove samples with any missing data
    keep <- keep & !is.na(dat)
    # outcome
    Y <- dat[keep]
    # design matrix
    W <- matrix(1,nrow=length(Y),ncol=1)
  }

  # sample size
  n <- sum(keep)

  # if covMatList is a matrix, convert to a list
  if(class(covMatList) == "matrix"){
    covMatList <- list(A = covMatList)
  }

  # number of covariance structure matrices
  m <- length(covMatList)
  if(m > 1 & method != "AIREML"){
    stop("AIREML method must be used with Delta matrix")
  }

  # if covMatList doesn't have names, assign them
  if(is.null(names(covMatList))){
    names(covMatList) <- paste("A",1:m,sep="")
  }

  # check for starting values
  if(!is.null(start)){
    if(length(start) != (m+1)){
      stop("length of start must equal the length of covMatList + 1")
    }
  }
    
  # subest covariance structure matrices
  keepID <- scanID[keep]  
  for(i in 1:m){
    if(!all(keepID %in% colnames(covMatList[[i]]))){
      stop(paste("All of the included Samples must be in matrix", i, "of covMatList"))
    }
    # subset matrix
    keepMat <- colnames(covMatList[[i]]) %in% keepID
    covMatList[[i]] <- covMatList[[i]][keepMat,keepMat]

    # check that names match
    if(!all(colnames(covMatList[[i]]) == scanID[keep])){
      stop("Column and Row names of matrix", i, "of covMatList must match the scanIDs of scanAnnot")
    }
  }
  

  if(method == "AIREML"){
    if(verbose) message("Using AIREML Procedure...")
    if(verbose) message(paste("Sample Size: ", n))
    if(verbose) message("Computing Variance Component Estimates...")
    if(verbose) message(paste("Sigma^2_",c(names(covMatList),"E"),sep="", collapse="     "))

    # initial values
    sigma2.p <- var(Y)
    AIREML.tol <- AIREML.tol*sigma2.p  # set convergence tolerance dependent on trait
    
    if(is.null(start)){      
      sigma2.k <- rep((1/(m+1))*sigma2.p, (m+1))
    }else{
      sigma2.k <- as.vector(start)
    }
    zeroFLAG <- rep(FALSE, m+1) # indicator of elements that have converged to "0"
    

    reps <- 0
    repeat({
      reps <- reps+1
      
      # variance matrix
      V <- diag(rep(sigma2.k[m+1],n))
      for(i in 1:m){
        V <- V + covMatList[[i]]*sigma2.k[i]
      }      
      #V <- Reduce('+', mapply('*', covMatList, sigma2.k[1:m], SIMPLIFY=FALSE), init = sigma2.k[m+1]*diag(n))

      # inverse
      Vinv <- chol2inv(chol(V))

      # projection matrix
      VinvW <- crossprod(Vinv,W)
      P <- Vinv - tcrossprod(tcrossprod(VinvW,chol2inv(chol(crossprod(W,VinvW)))),VinvW)

      # matrices for later use
      PY <- crossprod(P,Y)
      PPY <- crossprod(P,PY)      

      if(reps > 1  || !is.null(start)){
        # Average Information and Scores
        AI <- matrix(NA, nrow=(m+1), ncol=(m+1))
        score <- rep(NA,(m+1))        
        for(i in 1:m){
          PAPY <- crossprod(P,crossprod(covMatList[[i]],PY))
          score[i] <- -0.5*(sum(P*covMatList[[i]]) - crossprod(Y, PAPY)) 
          AI[i,i] <- 0.5*crossprod(PY, crossprod(covMatList[[i]],PAPY)) # YPAPAPY
          if((i+1) <= m){
            for(j in (i+1):m){
              AI[i,j] <- 0.5*crossprod(PY, crossprod(covMatList[[j]],PAPY)) # YPDPAPY, YPSPAPY
              AI[j,i] <- AI[i,j]
            }
          }
          AI[i,(m+1)] <- 0.5*crossprod(PY, crossprod(covMatList[[i]],PPY)) # YPAPPY
          AI[(m+1),i] <- AI[i,(m+1)]
        }
        score[m+1] <- -0.5*(sum(diag(P)) - crossprod(Y, PPY))
        AI[(m+1),(m+1)] <- 0.5*crossprod(PY,PPY) # YPPPY        

        # update
        AIinvScore <- solve(AI, score)
        sigma2.kplus1 <- sigma2.k + AIinvScore
        
        sigma2.kplus1[zeroFLAG & sigma2.kplus1 < AIREML.tol] <- 0 # set elements that were previously "0" and are still < 0 back to 0 (prevents step-halving due to this component)
        
        # step-halving if step too far
        k <- 1
        while(!all(sigma2.kplus1 >= 0)){
          k <- 0.5*k
          sigma2.kplus1 <- sigma2.k + k*AIinvScore
          sigma2.kplus1[zeroFLAG & sigma2.kplus1 < AIREML.tol] <- 0
        }

      }else{
        # EM step
        sigma2.kplus1 <- rep(NA,(m+1))

        for(i in 1:m){
          PAPY <- crossprod(P,crossprod(covMatList[[i]],PY))
          sigma2.kplus1[i] <- (1/n)*((sigma2.k[i])^2*crossprod(Y,PAPY) + (n*sigma2.k[i] - (sigma2.k[i])^2*sum(P*covMatList[[i]])))
        }
        sigma2.kplus1[m+1] <- (1/n)*((sigma2.k[m+1])^2*crossprod(Y,PPY) + (n*sigma2.k[m+1] - (sigma2.k[m+1])^2*sum(diag(P))))
      }

      # print current estimates
      if(verbose) print(sigma2.kplus1)

      # test for convergence
      #stat <- max(abs(sigma2.kplus1 - sigma2.k))
      stat <- sqrt(sum((sigma2.kplus1 - sigma2.k)^2))
      sigma2.k <- sigma2.kplus1
      if(stat < AIREML.tol) break()
      
      zeroFLAG <- sigma2.k < AIREML.tol # which elements have converged to "0"
      sigma2.k[zeroFLAG] <- 0 # set these to 0
    })

    # estimates
    varComp <- sigma2.k
    names(varComp) <- paste("V_",c(names(covMatList),"E"),sep="")
    # get covariance of estimates
    varCompCov <- solve(AI)
    colnames(varCompCov) <- paste("V_",c(names(covMatList),"E"),sep="")
    rownames(varCompCov) <- paste("V_",c(names(covMatList),"E"),sep="")

    
  }else if(method == "GridSearch"){
    if(verbose) message("Computing Eigen Decompositon of Phi Matrix...")
    q <- ncol(W)
    n.q <- n-q
    M <- diag(n) - tcrossprod(tcrossprod(W,chol2inv(chol(crossprod(W)))),W)
    EIG <- eigen(tcrossprod(tcrossprod(M,Phi),M), symmetric = TRUE)
    lambda <- EIG$values[1:n.q]
    eta <- crossprod(EIG$vectors[,1:n.q],Y); rm(EIG)
    etasq <- eta^2
    
    # EMMA grid search
    if(verbose) message("Performing Grid Search...")
    logdelta <- seq(from=llim, to=ulim, length=ngrid)
    delta <- exp(logdelta)
    Lambdas <- matrix(lambda, nrow = n.q, ncol = ngrid) + matrix(delta, nrow = n.q, ncol = ngrid, byrow=TRUE)
    Etasq <- matrix(etasq, nrow = n.q, ncol = ngrid)
    dLL <- 0.5*delta*(n.q*colSums(Etasq/(Lambdas*Lambdas))/colSums(Etasq/Lambdas) - colSums(1/Lambdas))

    LL.fn <- function(val,lambda,etasq){
      newval <- lambda + exp(val)
      0.5*(n.q*log(n.q/(2*pi)) - 1 - log(sum(etasq/newval)) - sum(log(newval)))
    }
    dLL.fn <- function(val,lambda,etasq){
      newval <- lambda + exp(val)
      0.5*(n.q*sum(etasq/newval^2)/sum(etasq/newval) - sum(1/newval))
    }

    optlogdelta <- vector(length=0)
    optLL <- vector(length=0)
    if(dLL[1] < grid.tol){
      optlogdelta <- append(optlogdelta, llim)
      optLL <- append(optLL, LL.fn(llim,lambda,etasq))
    }
    if(dLL[ngrid] > 0-grid.tol){
      optlogdelta <- append(optlogdelta, ulim)
      optLL <- append(optLL, LL.fn(ulim,lambda,etasq))
    }
    grid.tol2 <- grid.tol^2
    for(i in 1:(ngrid-1)){
      if(dLL[i]*dLL[i+1] < 0-grid.tol2 && dLL[i] > 0 && dLL[i+1] < 0){
        r <- uniroot(dLL.fn, lower=logdelta[i], upper=logdelta[i+1], lambda=lambda, etasq=etasq)
        optlogdelta <- append(optlogdelta, r$root)
        optLL <- append(optLL, LL.fn(r$root,lambda,etasq))
      }
    }

    maxdelta <- exp(optlogdelta[which.max(optLL)])
    maxLL <- max(optLL)
    Va <- sum(etasq/(lambda+maxdelta))/n.q
    Ve <- Va*maxdelta
    varComp <- c(Va,Ve)
    varCompCov <- NULL
  }
    
  if(verbose) message("Computing Cholesky Decomposition of Inverse Covariance Matrix of Phenotype...")
  # Covariance Matrix
  Sigma <- diag(n)*varComp[m+1]
  for(i in 1:m){
    Sigma <- Sigma + covMatList[[i]]*varComp[i]
  }
  # Inverse
  SigmaInv <- chol2inv(chol(Sigma))
  # Cholesky Decomposition
  cholSigmaInv <- t(chol(SigmaInv))
  colnames(cholSigmaInv) <- colnames(covMatList[[1]])
  rownames(cholSigmaInv) <- rownames(covMatList[[1]])

  return(list(varComp = varComp, varCompCov = varCompCov, cholSigmaInv = cholSigmaInv))
  
}


# x is the output from estVarComp
estVarCompCI <- function(x, prop=TRUE){
  if(prop){
    ci <- matrix(NA, nrow=length(x$varComp), ncol=2)
    est <- x$varComp/sum(x$varComp)
    for(i in 1:length(x$varComp)){
      deltaH <- rep(-x$varComp[i]/(sum(x$varComp)^2),length(x$varComp))
      deltaH[i] <- deltaH[i] + sum(x$varComp)/(sum(x$varComp)^2)
      varH <- crossprod(deltaH, crossprod(x$varCompCov, deltaH)) 
      ci[i,] <- est[i] + sqrt(varH)*qnorm(c(0.025,0.975))
    }
    res <- as.data.frame(cbind(est, ci))
    names(res) <- c("Proportion", "Lower 95", "Upper 95")
    
  }else{
    ci <- x$varComp + sqrt(diag(x$varCompCov)) %o% qnorm(c(0.025,0.975))
    res <- as.data.frame(cbind(x$varComp, ci))
    names(res) <- c("Est", "Lower 95", "Upper 95")
  }

  print(res)
}
