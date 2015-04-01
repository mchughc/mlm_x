
powerGraph <- function(totalRes,b,fn,xlims=c(0,1)){
  # want to make an ROC curve of sorts
  # for various alpha levels calculate both false positive and true positive rate
  # then, plot those false pos rates vs true pos rates for various alpha levels
  
  alpha <- seq(from=1e-04,to=1,by=0.01)
  
  null_auto <- totalRes[totalRes$beta1==b&totalRes$model=="auto"&!totalRes$causal,]
  null_xauto <- totalRes[totalRes$beta1==b&totalRes$model=="both"&!totalRes$causal,]
  null_x <- totalRes[totalRes$beta1==b&totalRes$model=="x"&!totalRes$causal,]
  
  true_auto <- totalRes[totalRes$beta1==b&totalRes$model=="auto"&totalRes$causal,]
  true_xauto <- totalRes[totalRes$beta1==b&totalRes$model=="both"&totalRes$causal,]
  true_x <- totalRes[totalRes$beta1==b&totalRes$model=="x"&totalRes$causal,]
  
  iterAuto <- nrow(true_auto)
  iterBoth <- nrow(true_xauto)
  iterX <- nrow(true_x)
  
  fp_x <- rep(NA,length(alpha))
  fp_xauto <- rep(NA,length(alpha))
  fp_auto <- rep(NA,length(alpha))
  
  tp_x <- rep(NA,length(alpha))
  tp_xauto <- rep(NA,length(alpha))
  tp_auto <- rep(NA,length(alpha))
  
  # for each of the model results
  for(i in 1:length(alpha)){
    fp_x[i] <- sum(null_x$pval<alpha[i])
    fp_xauto[i] <- sum(null_xauto$pval<alpha[i])
    fp_auto[i] <- sum(null_auto$pval<alpha[i])
    
    tp_x[i] <- sum(true_x$pval<alpha[i])
    tp_xauto[i] <- sum(true_xauto$pval<alpha[i])
    tp_auto[i] <- sum(true_auto$pval<alpha[i])
  }
  
  
  pdf(fn)
  plot(fp_xauto/iterBoth,tp_xauto/iterBoth,type="l",xlab="False Positive Rate",ylab="True Positive Rate",lwd=2,
       main=bquote(beta * "=" * .(b)),xlim=xlims,cex.lab=1.5,cex.main=1.5,cex.axis=1.5)
  points(fp_x/iterX,tp_x/iterX,type="l",col="cyan",lwd=2)
  points(fp_auto/iterAuto,tp_auto/iterAuto,type="l",col="red",lwd=2)
  legend("bottomright",c("Adj for X and Auto","Adj for Auto","Adj for X"),col=c("black","red","cyan"),
         cex=1.2,lty=1,lwd=2)
  dev.off()
  
  return(NULL)
}
