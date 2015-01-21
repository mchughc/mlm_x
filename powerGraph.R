
powerGraph <- function(totalRes,b,fn,xlims=c(0,1)){
  # want to make an ROC curve of sorts
  # for various alpha levels, or false positive rates, find the number of true positives
  
  alpha <- seq(from=1e-04,to=1,by=0.01)
  
  totalResSm <- totalRes[totalRes$beta1==b&totalRes$causal&totalRes$model==3,]
  totalResF <- totalRes[totalRes$beta1==b&totalRes$causal&totalRes$model==4,]
  totalResX <- totalRes[totalRes$beta1==b&totalRes$causal&totalRes$model==2,]
  
  iter <- sum(totalRes$beta1==b&totalRes$model==3)/2
  stopifnot(iter==sum(totalRes$beta1==b&totalRes$model==4)/2)
  stopifnot(iter==sum(totalRes$beta1==b&totalRes$model==2)/2)
  y1 <- rep(NA,length(alpha)); y2 <- y1
  y3 <- y1
  for(i in 1:length(alpha)){
    y1[i] <- sum(totalResSm$pval<alpha[i])
    y2[i] <- sum(totalResF$pval<alpha[i])
    y3[i] <- sum(totalResX$pval<alpha[i])
  }
  
  pdf(fn)
  plot(alpha,y1/iter,type="l",xlab="False Positive Rate",ylab="True Positive Rate",lwd=2,
       main=bquote(beta * "=" * .(b)),xlim=xlims)
  points(alpha,y2/iter,type="l",col="red",lwd=2)
  points(alpha,y3/iter,type="l",col="cyan",lwd=2)
  legend("bottomright",c("Adj for X and Auto","Adj for Auto","Adj for X"),col=c("black","red","cyan"),lty=1,lwd=2)
  dev.off()
  
  return(NULL)
}
