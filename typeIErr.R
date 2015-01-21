typeIErr <- function(nSets,alpha){
  totalFP_true <- 0
  totalFP_false <- 0
  
  totalFP_T_2 <- 0
  totalFP_F_2 <- 0
  
  totalIter <- 0
  
  for(i in 1:nSets){
    fn <- paste("modelC_results_smBeta/mmRes_varComps_modelC_SNP",i,".RData",sep="")
    fn2 <- paste("modelC_results_smBeta2/mmRes_varComps_modelC_SNP",i,".RData",sep="")
    
    dat <- get(load(fn))
    dat2 <- get(load(fn2))
    
    # true model fit
    pvT <- dat[[2]]
    
    # false model fit
    pvF <- dat[[6]]
    
    stopifnot(nrow(pvT)==nrow(pvF))
    
    totalFP_true <- totalFP_true + sum(pvT$pval[-1]<alpha)
    totalFP_false <- totalFP_false + sum(pvF$pval[-1]<alpha)
    
    ## look at the second set of beta values
    pvT <- dat2[[2]]
    pvF <- dat2[[6]]
    
    stopifnot(nrow(pvT)==nrow(pvF))
    
    totalFP_T_2 <- totalFP_T_2 + sum(pvT$pval[-1]<alpha)
    totalFP_F_2 <- totalFP_F_2 + sum(pvF$pval[-1]<alpha)
    
    totalIter <- totalIter + length(pvT$pval[-1])
  }
  
  res <- c(totalIter,totalFP_true,totalFP_false,totalFP_T_2, totalFP_F_2)
  names(res) <- c("totalIterations","false_pos_trueModel","false_pos_falseModel",
                  "false_pos_trueModel_smBeta2","false_pos_falseModel_smBeta2")
  
  return(res)
}

