
if(verbose) message("Using AIREML Procedure...")
if(verbose) message(paste("Sample Size: ", n))
if(verbose) message("Computing Variance Component Estimates...")
if(verbose) message(paste("Sigma^2_",c(names(covMatList),"E"),sep="", collapse="     "))

# initial values
if(is.null(start)){
  sigma2.p <- var(Y)
  sigma2.k <- rep((1/(m+1))*sigma2.p, (m+1))
}else{
  sigma2.k <- as.vector(start)
}

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
    
    # step-halving if step too far
    k <- 1
    while(!all(sigma2.kplus1 > 0)){
      k <- 0.5*k
      sigma2.kplus1 <- sigma2.k + k*AIinvScore
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
  stat <- max(abs(sigma2.kplus1 - sigma2.k))
  sigma2.k <- sigma2.kplus1
  if(stat < AIREML.tol) break()
})

# estimates
varComp <- sigma2.k
names(varComp) <- paste("V_",c(names(covMatList),"E"),sep="")
# get covariance of estimates
varCompCov <- solve(AI)
colnames(varCompCov) <- paste("V_",c(names(covMatList),"E"),sep="")
rownames(varCompCov) <- paste("V_",c(names(covMatList),"E"),sep="")




#####
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

