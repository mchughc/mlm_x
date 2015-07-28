CAnD_pooled <- function (chrAncest, bonfCorr = TRUE) 
{
  numChrs <- ncol(chrAncest)
  if (ncol(chrAncest) <= 1) {
    stop("chrAncest must have at least two columns of \n         ancestry proportions for testing.")
  }
  if (!all(apply(chrAncest, 2, class) == "numeric")) {
    stop("Ancestry proportions must be numeric values.")
  }
  if (!sum(is.na(chrAncest)) == 0) {
    warning("NA values will be excluded from the analysis.")
  }
  if (nrow(chrAncest) <= 20) {
    warning("The number of samples may be too small for assumptions to hold.\n            Consider running 'nonParam_CAnD' instead.")
  }
  diff_means <- getDiffMatrices(chrAncest, diff = FALSE)
  pval <- rep(NA, numChrs)
  for (i in 1:numChrs) {
    pval[i] <- t.test(chrAncest[, i], diff_means[, i], paired = FALSE)$p.value
  }
  cand_stat <- -2 * sum(log(pval))
  cand <- 1 - pchisq(cand_stat, df = 2 * numChrs)
  if (bonfCorr) {
    pval <- pval * numChrs
  }
  pval <- ifelse(pval > 1, 1, pval)
  names(pval) <- colnames(chrAncest)
  return(new("CAnDResult", test = "parametric", pValues = pval, 
             overallStatistic = cand_stat, overallpValue = cand, BonfCorr = bonfCorr))
}