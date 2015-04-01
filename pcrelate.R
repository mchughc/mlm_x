# define generic
pcrelate <- function(genoData, MAF = 0.05, pcMat = NULL, unrel.set = NULL, ibd.prop = TRUE, snp.include = NULL, scan.include = NULL, block.size = 10000, snpfirstdim = TRUE, correct = TRUE) UseMethod("pcrelate")

# defined method for matrix
pcrelate.matrix <- function(genoData, MAF = 0.05, pcMat = NULL, unrel.set = NULL, ibd.prop = TRUE, scan.include = NULL, block.size = 10000, correct = TRUE){
	
	# checks
	if(MAF < 0 | MAF >= 0.5){
		stop("MAF must be > 0 and  < 0.5")
	}
	
	# samples to include
	if(is.null(scan.include)){
		scan.include <- 1:dim(genoData)[2]
	}
	# sample size
	nsamp <- length(scan.include)
	
	# PC Checks
	if(is.null(pcMat)){
		message("pcMat not specified, Calculating Unadjusted Relatedness Estimates")
		method <- "Unadjusted"
	}else{
		pcMat <- as.matrix(pcMat)		
		if(nsamp != nrow(pcMat)){
			stop("The number of Samples in genoData and pcMat do not match.")
		}
		message(paste("Adjusting for",ncol(pcMat),"PC(s)"))
		method <- "PC-REAP"
	}	
	
	# compute Hat matrix
	message("Computing Hat Matrix...")
	# PC matrix
	X <- cbind(rep(1,nsamp),pcMat)
	if(is.null(unrel.set)){
		# hat matrix: X(X'X)^{-1}X'
		H <- tcrossprod( tcrossprod( X, chol2inv(chol(crossprod(X))) ), X )
	}else{
		Xu <- X[unrel.set,]
		# hat matrix: X(X_u'X_u)^{-1}X_u'
		H <- tcrossprod( tcrossprod( X, chol2inv(chol(crossprod(Xu))) ), Xu )
	}
	
	# determine SNP blocks
	# number of snps
	nloci <- dim(genoData)[1]
	message(paste("Running Analysis with",nloci,"SNPs ..."))
	# number of blocks of snps
	nblocks <- ceiling(nloci/block.size)
	# start and end positions for blocks
	if(nblocks==1){
		snp.start <- 1
		snp.end <- nloci
	}else{
		snp.start <- (0:(nblocks-1))*block.size+1
		snp.end <- c((1:(nblocks-1))*block.size, nloci)
	}
	
	
	# results matrices
	message("Creating data.frames to Store Results...")
    if(ibd.prop){
		res <- data.frame(ID1 = rep(scan.include[-nsamp], times=((nsamp-1):1)),
                          ID2 = unlist(lapply(2:nsamp,function(x){ scan.include[x:nsamp] })),
                          nsnp = rep(0, nsamp*(nsamp-1)/2),
                          kin = rep(0, nsamp*(nsamp-1)/2),
                          k2 = rep(0, nsamp*(nsamp-1)/2),
                          k1 = rep(0, nsamp*(nsamp-1)/2),
                          k0 = rep(0, nsamp*(nsamp-1)/2))
        kindenom = rep(0, nsamp*(nsamp-1)/2)
        k2denom = rep(0, nsamp*(nsamp-1)/2)
        k0denom = rep(0, nsamp*(nsamp-1)/2)
                
	}else{
		res <- data.frame(ID1 = rep(scan.include[-nsamp], times=((nsamp-1):1)),
                          ID2 = unlist(lapply(2:nsamp,function(x){ scan.include[x:nsamp] })),
                          nsnp = rep(0, nsamp*(nsamp-1)/2),
                          kin = rep(0, nsamp*(nsamp-1)/2))                
                kindenom = rep(0, nsamp*(nsamp-1)/2)
	}
	inbreed <- data.frame(ID = scan.include,
                          nsnp = rep(0,nsamp),
                          f = rep(0,nsamp))
    fdenom = rep(0,nsamp)
	
	
	for(bn in 1:nblocks){
		message(paste("Computing Individual Specific Allele Frequencies: Block",bn,"of",nblocks,"..."))
		# load genotype data
		geno <- genoData[snp.start[bn]:snp.end[bn],]
                # set any negative or "3" genotype values to NA
		geno[geno < 0 | geno==3] <- NA
		 
		# minor allele frequency filter
		if(is.null(unrel.set)){
			# allele freq from entire sample
			pA <- 0.5*rowMeans(geno, na.rm = TRUE)
		}else{
			# allele freq from unrelated set
			pA <- 0.5*rowMeans(geno[,unrel.set], na.rm = TRUE)
		}
		# remove SNPs with low MAF
		snp.excl <- which(is.na(pA) | pA < MAF | pA > (1-MAF))
		if(length(snp.excl > 0)){
			geno <- geno[-snp.excl,]
			pA <- pA[-snp.excl]
		}
		# number of snps in the block
		nsnp.block <- dim(geno)[1]
		
		# rather than compute new H for each SNP, impute missing genotype values to sample mean
		# which genotype values are missing
		miss.idx <- which(is.na(geno))
		# if there are missing genotype values
		if(length(miss.idx) > 0){
			# index of which snps have the missing values
			snp.idx <- miss.idx %% nsnp.block; snp.idx[snp.idx==0] <- nsnp.block
			# replace missing genotypes with twice the sample allele frequency
			geno[miss.idx] <- 2*pA[snp.idx]
		}
		
		# matrix of individual specific allele frequencies
		if(is.null(unrel.set)){
			phat <- 0.5*tcrossprod(geno,H)
		}else{
			phat <- 0.5*tcrossprod(geno[,unrel.set],H)
		}
		# opposite allele freq
		qhat <- 1-phat
		
		# plug in a temporary value where genotypes are missing (for filtering)
		phat[miss.idx] <- -1
		# which snps to filter (fitted value too big/small)
		filt.idx <- which(phat < MAF | phat > (1-MAF))
			
		# set missing & filtered values to 0 (so no contribution)
		geno[filt.idx] <- 0
		phat[filt.idx] <- 0
		qhat[filt.idx] <- 0
		
		# product (filtered values already 0)
		phatqhat <- phat*qhat
							
		if(ibd.prop){
			# update inbreeding denominator 
			fdenom <- fdenom + colSums(phatqhat)
			
			message(paste("Computing Dominance Covariance Matrix: Block",bn,"of",nblocks,"..."))			
			# create dominance coded matrix			
			Xd <- matrix(0, nrow = nsnp.block, ncol = nsamp)
			idx0 <- which(geno == 0)
			idx2 <- which(geno == 2)
			Xd[idx0] <- phat[idx0]
			Xd[idx2] <- qhat[idx2]
			Xd <- Xd - phatqhat
			# set missing & filtered SNPs to 0 (so no contribution)
			Xd[filt.idx] <- 0
			
			# updated inbreeding numerator
			inbreed$f <- inbreed$f + colSums(Xd)
			
			# update numerator
			k2Mat <- crossprod(Xd); rm(Xd)
			res$k2 <- res$k2 + k2Mat[lower.tri(k2Mat)]; rm(k2Mat)

            # update denominator
			sumpqpqMat <- crossprod(phatqhat)
			k2denom <- k2denom + sumpqpqMat[lower.tri(sumpqpqMat)]; rm(sumpqpqMat)
			
			
			message(paste("Computing Count of Opposite Homozygotes: Block",bn,"of",nblocks,"..."))
			# update denominator
			sump2q2Mat <- crossprod(phat^2,qhat^2)
			k0denom <- k0denom + sump2q2Mat[lower.tri(sump2q2Mat)] + t(sump2q2Mat)[lower.tri(sump2q2Mat)]; rm(sump2q2Mat)
			
			# homozygotes
			IAA <- geno == 2
			Iaa <- geno == 0
			# set missing & filtered SNPs to FALSE (so no contribution)
			IAA[filt.idx] <- FALSE; Iaa[filt.idx] <- FALSE
						
			# update numerator
			NAAaaMat <- crossprod(IAA,Iaa); rm(IAA); rm(Iaa)
			res$k0 <- res$k0 + NAAaaMat[lower.tri(NAAaaMat)] + t(NAAaaMat)[lower.tri(NAAaaMat)]; rm(NAAaaMat)
			
		}else{
			# update inbreeding numerator
			inbreed$f <- inbreed$f + colSums(geno^2) - colSums(geno) - 2*colSums(phat*geno) + 2*colSums(phat^2)
			
			# update inbreeding denominator (filtered values already 0)			
			fdenom <- fdenom + 2*colSums(phatqhat)
		}
		
		message(paste("Computing Additive Covariance Matrix: Block",bn,"of",nblocks,"..."))
		# matrix of regression residuals (filtered values already 0)
		R <- geno-2*phat; rm(geno)
		
		# update numerator
		kinMat <- crossprod(R); rm(R)
		res$kin <- res$kin + kinMat[lower.tri(kinMat)]; rm(kinMat)	
		
		# update denominator
		sumsqrtpqpqMat <- crossprod(sqrt(phatqhat)); rm(phatqhat)
		kindenom <- kindenom + sumsqrtpqpqMat[lower.tri(sumsqrtpqpqMat)]; rm(sumsqrtpqpqMat)
		
		# update denominator
		#sumpiqjMat <- crossprod(phat,qhat)
		#kindenom <- kindenom + sumpiqjMat[lower.tri(sumpiqjMat)] + t(sumpiqjMat)[lower.tri(sumpiqjMat)]; rm(sumpiqjMat)	
		
		# update denominator (piqi + pjqj)
		#sumpqMat <- matrix(colSums(phatqhat), nrow = nsamp, ncol = nsamp)		
		#kindenom <- kindenom + sumpqMat[lower.tri(sumpqMat)] + t(sumpqMat)[lower.tri(sumpqMat)]; rm(sumpqMat)
				
		# determine the number of snps used for each pair
		snpcount <- matrix(1, nrow=nsnp.block, ncol=nsamp); snpcount[filt.idx] <- 0
		nsnpMat <- crossprod(snpcount); rm(snpcount)
		res$nsnp <- res$nsnp + nsnpMat[lower.tri(nsnpMat)]; 
		inbreed$nsnp <- inbreed$nsnp + diag(nsnpMat); rm(nsnpMat)
	}
	
	message(paste("Computing Final Estimates Across All Blocks ..."))
	# calculate inbreeding estimates
	inbreed$f <- inbreed$f/fdenom; rm(fdenom)
	
	# calculate kinship estimates
	res$kin <- res$kin/(4*kindenom); rm(kindenom)
	
	kincorrect <- NULL
	fcorrect <- NULL	
	if(!is.null(pcMat) & correct){
		# correct for possible overadjustment in each PC
		idxf <- which(inbreed$f < 2^(-11/2))
		for(i in 1:ncol(pcMat)){
			Avec <- pcMat[,i]			
			fvals <- lm(inbreed$f[idxf] ~ I(Avec[idxf]) + I(Avec[idxf]^2))$coef
			inbreed$f <- inbreed$f - fvals[1] - fvals[2]*Avec - fvals[3]*Avec^2	
			fcorrect <- append(fcorrect, fvals)
			
			idx <- which(res$kin < 2^(-11/2))
			Acov <- tcrossprod(Avec); Acov <- Acov[lower.tri(Acov)]
			kvals <- lm(res$kin[idx] ~ Acov[idx])$coef
			res$kin <- res$kin - kvals[1] - kvals[2]*Acov; rm(Acov); rm(Avec)
			kincorrect <- append(kincorrect, kvals)
		}
	}
	
	k2correct <- NULL
	if(ibd.prop){
		# pairwise inbreeding product
		fprod <- tcrossprod(inbreed$f)
		fprod <- fprod[lower.tri(fprod)]
		
		# calculate k2 estimates
		res$k2 <- res$k2/k2denom - fprod; rm(k2denom)
		
		# correction
		if(correct){
			if(!(is.null(pcMat))){
				for(i in 1:ncol(pcMat)){
					Acov <- tcrossprod(pcMat[,i]); Acov <- Acov[lower.tri(Acov)]
					idx <- which(res$kin < 2^(-11/2))
					k2vals <- lm(res$k2[idx] ~ I(Acov[idx]) + I(Acov[idx]^2))$coef
					res$k2 <- res$k2 - k2vals[1] - k2vals[2]*Acov - k2vals[3]*Acov^2; rm(Acov)
					k2correct <- append(k2correct, k2vals)
				}
			}			
			idx <- which(res$k2 < 2^(-9/2))
			k2vals <- lm(res$k2[idx] ~ res$kin[idx])$coef
			res$k2 <- res$k2 - k2vals[1] - k2vals[2]*res$kin
			k2correct <- append(k2correct, k2vals)
		}
			
		# calculate k0 estimates
		res$k0 <- res$k0/k0denom; rm(k0denom)
                
		# index for not PO/FS
		rel2idx <- which(res$kin < 2^(-5/2))
		res$k0[rel2idx] <- (1 - 4*res$kin[rel2idx] + res$k2[rel2idx])
		
		# calculate k1 estimates
		res$k1 <- 1 - res$k2 - res$k0
	}
	
		
	# return results
	out <- list(kinship = res,
                    inbreed = inbreed,
                    kincorrect = kincorrect,
                    fcorrect = fcorrect,
                    k2correct = k2correct,
                    call = match.call(),
                    method = method)				
	class(out) <- "pcrelate"
	return(out)
}


# define method for GWASTools GenotypeData
pcrelate.GenotypeData <- function(genoData, MAF = 0.05, pcMat = NULL, unrel.set = NULL, ibd.prop = TRUE, snp.include = NULL, scan.include = NULL, block.size = 10000, correct = TRUE){
	
	# checks
	if(MAF < 0 | MAF >= 0.5){
		stop("MAF must be > 0 and  < 0.5")
	}
	
	# snps to include
	if(!is.null(snp.include)){
		if(!all(is.element(snp.include,getSnpID(genoData)))){
			stop("Not all of the SNP IDs in snp.include are in genoData")
		}
	}else{
		# all autosomal SNPs
		snp.include <- getSnpID(genoData)[getChromosome(genoData) <= 22]
	}
	
    # load scanID
    scanID <- getScanID(genoData)
    
    # samples to include - get the index
    if(!is.null(scan.include)){
        scan.include <- which(is.element(scanID,scan.include))
    }else{
        scan.include <- scanID
    }
    # sample size
    nsamp <- length(scan.include)
    
    # check that scanID matches unrel.set
    if(!is.null(unrel.set)){
        if(!all(unrel.set %in% scanID)){
            stop("All of the samples in unrel.set must be in the scanIDs of genoData")
        }
        # create index for unrealted set
        unrel.idx <- which(scanID %in% unrel.set)
    }
	
	# PC Checks
	if(is.null(pcMat)){
		message("pcMat not specified, Calculating Unadjusted Relatedness Estimates")
		method <- "Unadjusted"
	}else{
		pcMat <- as.matrix(pcMat)		
		if(nsamp != nrow(pcMat)){
			stop("The number of Samples in genoData and pcMat do not match.")
		}
		message(paste("Adjusting for",ncol(pcMat),"PC(s)"))
		method <- "PC-REAP"
	}	
		
	# compute Hat matrix
	message("Computing Hat Matrix...")
	# PC matrix
	X <- cbind(rep(1,nsamp),pcMat)
	if(is.null(unrel.set)){
		# hat matrix: X(X'X)^{-1}X'
		H <- tcrossprod( tcrossprod( X, chol2inv(chol(crossprod(X))) ), X )
	}else{
		Xu <- X[unrel.idx,]
		# hat matrix: X(X_u'X_u)^{-1}X_u'
		H <- tcrossprod( tcrossprod( X, chol2inv(chol(crossprod(Xu))) ), Xu )
	}
	
	# determine SNP blocks	
	# number of snps to be used
	nloci <- length(snp.include)
	message(paste("Running Analysis with",nloci,"SNPs ..."))
	# number of blocks of snps
	nblocks <- ceiling(nloci/block.size)
	# get index of SNP number in netCDF to include
	snp.include.idx <- which(is.element(getSnpID(genoData),snp.include))
	# start and end positions for blocks
	if(nblocks==1){
		snp.start <- snp.include.idx[1]
		snp.end <- snp.include.idx[nloci]
	}else{
		snp.start <- snp.include.idx[(0:(nblocks-1))*block.size+1]
		snp.end <- snp.include.idx[c((1:(nblocks-1))*block.size, nloci)]
	}
	
	# results matrix
	message("Creating data.frames to Store Results...")
        if(ibd.prop){
		res <- data.frame(ID1 = rep(scan.include[-nsamp], times=((nsamp-1):1)),
                                  ID2 = unlist(lapply(2:nsamp,function(x){ scan.include[x:nsamp] })),
                                  nsnp = rep(0, nsamp*(nsamp-1)/2),
                                  kin = rep(0, nsamp*(nsamp-1)/2),
                                  k2 = rep(0, nsamp*(nsamp-1)/2),
                                  k1 = rep(0, nsamp*(nsamp-1)/2),
                                  k0 = rep(0, nsamp*(nsamp-1)/2))
                kindenom = rep(0, nsamp*(nsamp-1)/2)
                k2denom = rep(0, nsamp*(nsamp-1)/2)
                k0denom = rep(0, nsamp*(nsamp-1)/2)
                
	}else{
		res <- data.frame(ID1 = rep(scan.include[-nsamp], times=((nsamp-1):1)),
                                  ID2 = unlist(lapply(2:nsamp,function(x){ scan.include[x:nsamp] })),
                                  nsnp = rep(0, nsamp*(nsamp-1)/2),
                                  kin = rep(0, nsamp*(nsamp-1)/2))                
                kindenom = rep(0, nsamp*(nsamp-1)/2)
	}
	inbreed <- data.frame(ID = scan.include,
                              nsnp = rep(0,nsamp),
                              f = rep(0,nsamp))
        fdenom = rep(0,nsamp)
	
	
	for(bn in 1:nblocks){
		message(paste("Computing Individual Specific Allele Frequencies: Block",bn,"of",nblocks,"..."))
		# load genotype data
		geno <- getGenotype(genoData, scan=c(1,-1), snp=c(snp.start[bn], snp.end[bn]-snp.start[bn]+1) )
		# subset included snps
		geno <- geno[is.element(getSnpID(genoData, snp.start[bn]:snp.end[bn]), snp.include),]
		# remove exclude samples
		if(!is.null(scan.include)){
			geno <- geno[,scan.include]
		}
		# set any negative genotype values to NA
		geno[geno < 0] <- NA
		
		# minor allele frequency filter
		if(is.null(unrel.set)){
			# allele freq from entire sample
			pA <- 0.5*rowMeans(geno, na.rm = TRUE)
		}else{
			# allele freq from unrelated set
			pA <- 0.5*rowMeans(geno[,unrel.idx], na.rm = TRUE)
		}
		# remove SNPs with low MAF
		snp.excl <- which(is.na(pA) | pA < MAF | pA > (1-MAF))
		if(length(snp.excl > 0)){
			geno <- geno[-snp.excl,]
			pA <- pA[-snp.excl]
		}
		# number of snps in the block
		nsnp.block <- dim(geno)[1]
		
		# rather than compute new H for each SNP, impute missing genotype values to sample mean
		# which genotype values are missing
		miss.idx <- which(is.na(geno))
		# if there are missing genotype values
		if(length(miss.idx) > 0){
			# index of which snps have the missing values
			snp.idx <- miss.idx %% nsnp.block; snp.idx[snp.idx==0] <- nsnp.block
			# replace missing genotypes with twice the sample allele frequency
			geno[miss.idx] <- 2*pA[snp.idx]
		}
		
		# matrix of individual specific allele frequencies
		if(is.null(unrel.set)){
			phat <- 0.5*tcrossprod(geno,H)
		}else{
			phat <- 0.5*tcrossprod(geno[,unrel.idx],H)
		}
		# opposite allele freq
		qhat <- 1-phat
		
		# plug in a temporary value where genotypes are missing (for filtering)
		phat[miss.idx] <- -1
		# which snps to filter (fitted value too big/small)
		filt.idx <- which(phat < MAF | phat > (1-MAF))
			
		# set missing & filtered values to 0 (so no contribution)
		geno[filt.idx] <- 0
		phat[filt.idx] <- 0
		qhat[filt.idx] <- 0
		
		# product (filtered values already 0)
		phatqhat <- phat*qhat
		
		if(ibd.prop){
			# update inbreeding denominator 			
			fdenom <- fdenom + colSums(phatqhat)
			
			message(paste("Computing Dominance Covariance Matrix: Block",bn,"of",nblocks,"..."))			
			# create dominance coded matrix			
			Xd <- matrix(0, nrow = nsnp.block, ncol = nsamp)
			idx0 <- which(geno == 0)
			idx2 <- which(geno == 2)
			Xd[idx0] <- phat[idx0]
			Xd[idx2] <- qhat[idx2]
			Xd <- Xd - phatqhat
			# set missing & filtered SNPs to 0 (so no contribution)
			Xd[filt.idx] <- 0
			
			# updated inbreeding numerator
			inbreed$f <- inbreed$f + colSums(Xd)
			
			# update numerator
			k2Mat <- crossprod(Xd); rm(Xd)
			res$k2 <- res$k2 + k2Mat[lower.tri(k2Mat)]; rm(k2Mat)

            # update denominator
			sumpqpqMat <- crossprod(phatqhat)
			k2denom <- k2denom + sumpqpqMat[lower.tri(sumpqpqMat)]; rm(sumpqpqMat)
			
			
			message(paste("Computing Count of Opposite Homozygotes: Block",bn,"of",nblocks,"..."))
			# update denominator
			sump2q2Mat <- crossprod(phat^2,qhat^2)
			k0denom <- k0denom + sump2q2Mat[lower.tri(sump2q2Mat)] + t(sump2q2Mat)[lower.tri(sump2q2Mat)]; rm(sump2q2Mat)
			
			# homozygotes
			IAA <- geno == 2
			Iaa <- geno == 0
			# set missing & filtered SNPs to FALSE (so no contribution)
			IAA[filt.idx] <- FALSE; Iaa[filt.idx] <- FALSE
						
			# update numerator
			NAAaaMat <- crossprod(IAA,Iaa); rm(IAA); rm(Iaa)
			res$k0 <- res$k0 + NAAaaMat[lower.tri(NAAaaMat)] + t(NAAaaMat)[lower.tri(NAAaaMat)]; rm(NAAaaMat)
			
		}else{
			# update inbreeding numerator
			inbreed$f <- inbreed$f + colSums(geno^2) - colSums(geno) - 2*colSums(phat*geno) + 2*colSums(phat^2)
			
			# update inbreeding denominator (filtered values already 0)			
			fdenom <- fdenom + 2*colSums(phatqhat)
		}
		
		message(paste("Computing Additive Covariance Matrix: Block",bn,"of",nblocks,"..."))
		# matrix of regression residuals (filtered values already 0)
		R <- geno-2*phat; rm(geno)
		
		# update numerator
		kinMat <- crossprod(R); rm(R)
		res$kin <- res$kin + kinMat[lower.tri(kinMat)]; rm(kinMat)	
		
		# update denominator
		sumsqrtpqpqMat <- crossprod(sqrt(phatqhat)); rm(phatqhat)
		kindenom <- kindenom + sumsqrtpqpqMat[lower.tri(sumsqrtpqpqMat)]; rm(sumsqrtpqpqMat)
				
		# determine the number of snps used for each pair
		snpcount <- matrix(1, nrow=nsnp.block, ncol=nsamp); snpcount[filt.idx] <- 0
		nsnpMat <- crossprod(snpcount); rm(snpcount)
		res$nsnp <- res$nsnp + nsnpMat[lower.tri(nsnpMat)]; 
		inbreed$nsnp <- inbreed$nsnp + diag(nsnpMat); rm(nsnpMat)
	}
	
	message(paste("Computing Final Estimates Across All Blocks ..."))
	# calculate inbreeding estimates
	inbreed$f <- inbreed$f/fdenom; rm(fdenom)
	
	# calculate kinship estimates
	res$kin <- res$kin/(4*kindenom); rm(kindenom)
			
	kincorrect <- NULL
	fcorrect <- NULL	
	if(!is.null(pcMat) & correct){
		# correct for possible overadjustment in each PC
		idxf <- which(inbreed$f < 2^(-11/2))
		for(i in 1:ncol(pcMat)){
			Avec <- pcMat[,i]			
			fvals <- lm(inbreed$f[idxf] ~ I(Avec[idxf]) + I(Avec[idxf]^2))$coef
			inbreed$f <- inbreed$f - fvals[1] - fvals[2]*Avec - fvals[3]*Avec^2	
			fcorrect <- append(fcorrect, fvals)
			
			idx <- which(res$kin < 2^(-11/2))
			Acov <- tcrossprod(Avec); Acov <- Acov[lower.tri(Acov)]
			kvals <- lm(res$kin[idx] ~ Acov[idx])$coef
			res$kin <- res$kin - kvals[1] - kvals[2]*Acov; rm(Acov); rm(Avec)
			kincorrect <- append(kincorrect, kvals)
		}
	}
	
	k2correct <- NULL
	if(ibd.prop){
		# pairwise inbreeding product
		fprod <- tcrossprod(inbreed$f)
		fprod <- fprod[lower.tri(fprod)]
		
		# calculate k2 estimates
		res$k2 <- res$k2/k2denom - fprod; rm(k2denom)
		
		# correction
		if(correct){
			if(!(is.null(pcMat))){
				for(i in 1:ncol(pcMat)){
					Acov <- tcrossprod(pcMat[,i]); Acov <- Acov[lower.tri(Acov)]
					idx <- which(res$kin < 2^(-11/2))
					k2vals <- lm(res$k2[idx] ~ I(Acov[idx]) + I(Acov[idx]^2))$coef
					res$k2 <- res$k2 - k2vals[1] - k2vals[2]*Acov - k2vals[3]*Acov^2; rm(Acov)
					k2correct <- append(k2correct, k2vals)
				}
			}			
			idx <- which(res$k2 < 2^(-9/2))
			k2vals <- lm(res$k2[idx] ~ res$kin[idx])$coef
			res$k2 <- res$k2 - k2vals[1] - k2vals[2]*res$kin
			k2correct <- append(k2correct, k2vals)
		}
			
		# calculate k0 estimates
		res$k0 <- res$k0/k0denom; rm(k0denom)
                
		# index for not PO/FS
		rel2idx <- which(res$kin < 2^(-5/2))
		res$k0[rel2idx] <- (1 - 4*res$kin[rel2idx] + res$k2[rel2idx])
		
		# calculate k1 estimates
		res$k1 <- 1 - res$k2 - res$k0
	}
		

    # return results
	out <- list(kinship= res,
                    inbreed = inbreed,
                    kincorrect = kincorrect,
                    fcorrect = fcorrect,
                    k2correct = k2correct,
                    call = match.call(),
                    method = method)	
	class(out) <- "pcrelate"
	return(out)
}



# define method for gds.class data
pcrelate.gds.class <- function(genoData, MAF = 0.05, pcMat = NULL, unrel.set = NULL, ibd.prop = TRUE, snp.include = NULL, scan.include = NULL, block.size = 10000, snpfirstdim = TRUE, correct = TRUE){
	
	# checks
	if(MAF < 0 | MAF >= 0.5){
		stop("MAF must be > 0 and  < 0.5")
	}
	
	# snps to include
	if(!is.null(snp.include)){
		if(!all(is.element(snp.include,read.gdsn(index.gdsn(genoData,"snp.id"))))){
			stop("Not all of the SNP IDs in snp.include are in genoData")
		}
	}else{
		# all autosomal SNPs
		snp.include <- read.gdsn(index.gdsn(genoData,"snp.id"))[read.gdsn(index.gdsn(genoData,c("snp.chromosome"))) <= 22]
	}
	
	# samples to include
	if(is.null(scan.include)){
		scan.include <- read.gdsn(index.gdsn(genoData,"sample.id"))
	}
	# sample size
	nsamp <- length(scan.include)

	# PC checks
	# PC Checks
	if(is.null(pcMat)){
		message("pcMat not specified, Calculating Unadjusted Relatedness Estimates")
		method <- "Unadjusted"
	}else{
		pcMat <- as.matrix(pcMat)		
		if(nsamp != nrow(pcMat)){
			stop("The number of Samples in genoData and pcMat do not match.")
		}
		message(paste("Adjusting for",ncol(pcMat),"PC(s)"))
		method <- "PC-REAP"
	}
	
	# compute Hat matrix
	message("Computing Hat Matrix...")
	# PC matrix
	X <- cbind(rep(1,nsamp),pcMat)
	if(is.null(unrel.set)){
		# hat matrix: X(X'X)^{-1}X'
		H <- tcrossprod( tcrossprod( X, chol2inv(chol(crossprod(X))) ), X )
	}else{
		Xu <- X[unrel.set,]
		# hat matrix: X(X_u'X_u)^{-1}X_u'
		H <- tcrossprod( tcrossprod( X, chol2inv(chol(crossprod(Xu))) ), Xu )
	}
	
	# determine SNP blocks	
	# number of snps to be used
	nloci <- length(snp.include)
	message(paste("Running Analysis with",nloci,"SNPs ..."))
	# number of blocks of snps
	nblocks <- ceiling(nloci/block.size)
	# start and end positions for blocks
	if(nblocks==1){
		snp.start <- 1
		snp.end <- nloci
	}else{
		snp.start <- (0:(nblocks-1))*block.size+1
		snp.end <- c((1:(nblocks-1))*block.size, nloci)
	}
	
	
	# results matrix
	message("Creating data.frames to Store Results...")
	if(ibd.prop){
		res <- data.frame(ID1 = rep(scan.include[-nsamp], times=((nsamp-1):1)),
                                  ID2 = unlist(lapply(2:nsamp,function(x){ scan.include[x:nsamp] })),
                                  nsnp = rep(0, nsamp*(nsamp-1)/2),
                                  kin = rep(0, nsamp*(nsamp-1)/2),
                                  k2 = rep(0, nsamp*(nsamp-1)/2),
                                  k1 = rep(0, nsamp*(nsamp-1)/2),
                                  k0 = rep(0, nsamp*(nsamp-1)/2))
                kindenom = rep(0, nsamp*(nsamp-1)/2)
                k2denom = rep(0, nsamp*(nsamp-1)/2)
                k0denom = rep(0, nsamp*(nsamp-1)/2)
                
	}else{
		res <- data.frame(ID1 = rep(scan.include[-nsamp], times=((nsamp-1):1)),
                                  ID2 = unlist(lapply(2:nsamp,function(x){ scan.include[x:nsamp] })),
                                  nsnp = rep(0, nsamp*(nsamp-1)/2),
                                  kin = rep(0, nsamp*(nsamp-1)/2))                
                kindenom = rep(0, nsamp*(nsamp-1)/2)
	}
	inbreed <- data.frame(ID = scan.include,
                              nsnp = rep(0,nsamp),
                              f = rep(0,nsamp))
        fdenom = rep(0,nsamp)
	
	
	for(bn in 1:nblocks){
		message(paste("Computing Individual Specific Allele Frequencies: Block",bn,"of",nblocks,"..."))
		# load genotype data
		geno <- snpgdsGetGeno(genoData, sample.id = scan.include, snp.id = snp.include[snp.start[bn]:snp.end[bn]], verbose = FALSE)
		# transpose if samples stored as first dimension
		if(!snpfirstdim){
			geno <- t(geno)
		}
		# set any negative or "3" genotype values to NA
		geno[geno < 0 | geno==3] <- NA
		
		# minor allele frequency filter
		if(is.null(unrel.set)){
			# allele freq from entire sample
			pA <- 0.5*rowMeans(geno, na.rm = TRUE)
		}else{
			# allele freq from unrelated set
			pA <- 0.5*rowMeans(geno[,unrel.set], na.rm = TRUE)
		}
		# remove SNPs with low MAF
		snp.excl <- which(is.na(pA) | pA < MAF | pA > (1-MAF))
		if(length(snp.excl > 0)){
			geno <- geno[-snp.excl,]
			pA <- pA[-snp.excl]
		}
		# number of snps in the block
		nsnp.block <- dim(geno)[1]
		
		# rather than compute new H for each SNP, impute missing genotype values to sample mean
		# which genotype values are missing
		miss.idx <- which(is.na(geno))
		# if there are missing genotype values
		if(length(miss.idx) > 0){
			# index of which snps have the missing values
			snp.idx <- miss.idx %% nsnp.block; snp.idx[snp.idx==0] <- nsnp.block
			# replace missing genotypes with twice the sample allele frequency
			geno[miss.idx] <- 2*pA[snp.idx]
		}
		
		# matrix of individual specific allele frequencies
		if(is.null(unrel.set)){
			phat <- 0.5*tcrossprod(geno,H)
		}else{
			phat <- 0.5*tcrossprod(geno[,unrel.set],H)
		}
		# opposite allele freq
		qhat <- 1-phat
		
		# plug in a temporary value where genotypes are missing (for filtering)
		phat[miss.idx] <- -1
		# which snps to filter (fitted value too big/small)
		filt.idx <- which(phat < MAF | phat > (1-MAF))

        # set missing & filtered values to 0 (so no contribution)
		geno[filt.idx] <- 0
		phat[filt.idx] <- 0
		qhat[filt.idx] <- 0
		
		# product (filtered values already 0)
		phatqhat <- phat*qhat
		
		if(ibd.prop){
			# update inbreeding denominator			
			fdenom <- fdenom + colSums(phatqhat)
			
			message(paste("Computing Dominance Covariance Matrix: Block",bn,"of",nblocks,"..."))			
			# create dominance coded matrix			
			Xd <- matrix(0, nrow = nsnp.block, ncol = nsamp)
			idx0 <- which(geno == 0)
			idx2 <- which(geno == 2)
			Xd[idx0] <- phat[idx0]
			Xd[idx2] <- qhat[idx2]
			Xd <- Xd - phatqhat
			# set missing & filtered SNPs to 0 (so no contribution)
			Xd[filt.idx] <- 0
			
			# updated inbreeding numerator
			inbreed$f <- inbreed$f + colSums(Xd)
			
			# update numerator
			k2Mat <- crossprod(Xd); rm(Xd)
			res$k2 <- res$k2 + k2Mat[lower.tri(k2Mat)]; rm(k2Mat)

            # update denominator
			sumpqpqMat <- crossprod(phatqhat)
			k2denom <- k2denom + sumpqpqMat[lower.tri(sumpqpqMat)]; rm(sumpqpqMat)
			
			
			message(paste("Computing Count of Opposite Homozygotes: Block",bn,"of",nblocks,"..."))
			# update denominator
			sump2q2Mat <- crossprod(phat^2,qhat^2)
			k0denom <- k0denom + sump2q2Mat[lower.tri(sump2q2Mat)] + t(sump2q2Mat)[lower.tri(sump2q2Mat)]; rm(sump2q2Mat)
			
			# homozygotes
			IAA <- geno == 2
			Iaa <- geno == 0
			# set missing & filtered SNPs to FALSE (so no contribution)
			IAA[filt.idx] <- FALSE; Iaa[filt.idx] <- FALSE
						
			# update numerator
			NAAaaMat <- crossprod(IAA,Iaa); rm(IAA); rm(Iaa)
			res$k0 <- res$k0 + NAAaaMat[lower.tri(NAAaaMat)] + t(NAAaaMat)[lower.tri(NAAaaMat)]; rm(NAAaaMat)
			
		}else{
			# update inbreeding numerator
			inbreed$f <- inbreed$f + colSums(geno^2) - colSums(geno) - 2*colSums(phat*geno) + 2*colSums(phat^2)
			
			# update inbreeding denominator (filtered values already 0)			
			fdenom <- fdenom + 2*colSums(phatqhat)
		}
		
		message(paste("Computing Additive Covariance Matrix: Block",bn,"of",nblocks,"..."))
		# matrix of regression residuals (filtered values already 0)
		R <- geno-2*phat; rm(geno)
		
		# update numerator
		kinMat <- crossprod(R); rm(R)
		res$kin <- res$kin + kinMat[lower.tri(kinMat)]; rm(kinMat)	
		
		# update denominator
		sumsqrtpqpqMat <- crossprod(sqrt(phatqhat)); rm(phatqhat)
		kindenom <- kindenom + sumsqrtpqpqMat[lower.tri(sumsqrtpqpqMat)]; rm(sumsqrtpqpqMat)
				
		# determine the number of snps used for each pair
		snpcount <- matrix(1, nrow=nsnp.block, ncol=nsamp); snpcount[filt.idx] <- 0
		nsnpMat <- crossprod(snpcount); rm(snpcount)
		res$nsnp <- res$nsnp + nsnpMat[lower.tri(nsnpMat)]; 
		inbreed$nsnp <- inbreed$nsnp + diag(nsnpMat); rm(nsnpMat)
	}
	
	message(paste("Computing Final Estimates Across All Blocks ..."))
	# calculate inbreeding estimates
	inbreed$f <- inbreed$f/fdenom; rm(fdenom)
	
	# calculate kinship estimates
	res$kin <- res$kin/(4*kindenom); rm(kindenom)
		
	kincorrect <- NULL
	fcorrect <- NULL	
	if(!is.null(pcMat) & correct){
		# correct for possible overadjustment in each PC
		idxf <- which(inbreed$f < 2^(-11/2))
		for(i in 1:ncol(pcMat)){
			Avec <- pcMat[,i]			
			fvals <- lm(inbreed$f[idxf] ~ I(Avec[idxf]) + I(Avec[idxf]^2))$coef
			inbreed$f <- inbreed$f - fvals[1] - fvals[2]*Avec - fvals[3]*Avec^2	
			fcorrect <- append(fcorrect, fvals)
			
			idx <- which(res$kin < 2^(-11/2))
			Acov <- tcrossprod(Avec); Acov <- Acov[lower.tri(Acov)]
			kvals <- lm(res$kin[idx] ~ Acov[idx])$coef
			res$kin <- res$kin - kvals[1] - kvals[2]*Acov; rm(Acov); rm(Avec)
			kincorrect <- append(kincorrect, kvals)
		}
	}
	
	k2correct <- NULL
	if(ibd.prop){
		# pairwise inbreeding product
		fprod <- tcrossprod(inbreed$f)
		fprod <- fprod[lower.tri(fprod)]
		
		# calculate k2 estimates
		res$k2 <- res$k2/k2denom - fprod; rm(k2denom)
		
		# correction
		if(correct){
			if(!(is.null(pcMat))){
				for(i in 1:ncol(pcMat)){
					Acov <- tcrossprod(pcMat[,i]); Acov <- Acov[lower.tri(Acov)]
					idx <- which(res$kin < 2^(-11/2))
					k2vals <- lm(res$k2[idx] ~ I(Acov[idx]) + I(Acov[idx]^2))$coef
					res$k2 <- res$k2 - k2vals[1] - k2vals[2]*Acov - k2vals[3]*Acov^2; rm(Acov)
					k2correct <- append(k2correct, k2vals)
				}
			}			
			idx <- which(res$k2 < 2^(-9/2))
			k2vals <- lm(res$k2[idx] ~ res$kin[idx])$coef
			res$k2 <- res$k2 - k2vals[1] - k2vals[2]*res$kin
			k2correct <- append(k2correct, k2vals)
		}
			
		# calculate k0 estimates
		res$k0 <- res$k0/k0denom; rm(k0denom)
                
		# index for not PO/FS
		rel2idx <- which(res$kin < 2^(-5/2))
		res$k0[rel2idx] <- (1 - 4*res$kin[rel2idx] + res$k2[rel2idx])
		
		# calculate k1 estimates
		res$k1 <- 1 - res$k2 - res$k0
	}
			
	
	# return results
	out <- list(kinship= res,
                    inbreed = inbreed,
                    kincorrect = kincorrect,
                    fcorrect = fcorrect,
                    k2correct = k2correct,
                    call = match.call(),
                    method = method)	
	class(out) <- "pcrelate"
	return(out)
}



###########################################################################################
# define generic
pcrelateChr <- function(genoData, chromosome, MAF = 0.05, pcMat = NULL, unrel.set = NULL, ibd.prop = TRUE, snp.include = NULL, scan.include = NULL, block.size = 10000, snpfirstdim = TRUE) UseMethod("pcrelateChr")


# define method for gds.class data
pcrelateChr.gds.class <- function(genoData, chromosome, MAF = 0.05, pcMat = NULL, unrel.set = NULL, ibd.prop = TRUE, snp.include = NULL, scan.include = NULL, block.size = 10000, snpfirstdim = TRUE){
	
	# checks
	if(MAF < 0 | MAF >= 0.5){
		stop("MAF must be > 0 and  < 0.5")
	}
	
	# snps to include
	if(!is.null(snp.include)){
		if(!all(is.element(snp.include,read.gdsn(index.gdsn(genoData,"snp.id"))))){
			stop("Not all of the SNP IDs in snp.include are in genoData")
		}
		# subset SNPs on that chromosome
		snp.include <- snp.include[is.element(snp.include,read.gdsn(index.gdsn(genoData,"snp.id"))[read.gdsn(index.gdsn(genoData,"snp.chromosome")) == chromosome])]
	}else{
		# all SNPs on that chromosome
		snp.include <- read.gdsn(index.gdsn(genoData,"snp.id"))[read.gdsn(index.gdsn(genoData,"snp.chromosome")) == chromosome]
	}
	
	# samples to include
	if(is.null(scan.include)){
		scan.include <- read.gdsn(index.gdsn(genoData,"sample.id"))
	}
	# sample size
	nsamp <- length(scan.include)

	# PC checks
	if(is.null(pcMat)){
		message("pcMat not specified, Calculating Unadjusted Relatedness Estimates")
		method <- "Unadjusted"
	}else{
		pcMat <- as.matrix(pcMat)
		if(nsamp != nrow(pcMat)){
			stop("The number of Samples used in genoData and pcMat do not match.")
		}
		message(paste("Adjusting for",ncol(pcMat),"PC(s)"))
		method <- "PC-REAP"
	}
	
	# compute Hat matrix
	message("Computing Hat Matrix...")
	# PC matrix
	X <- cbind(rep(1,nsamp),pcMat)
	if(is.null(unrel.set)){
		# hat matrix: X(X'X)^{-1}X'
		H <- tcrossprod( tcrossprod( X, chol2inv(chol(crossprod(X))) ), X )
	}else{
		Xu <- X[unrel.set,]
		# hat matrix: X(X_u'X_u)^{-1}X_u'
		H <- tcrossprod( tcrossprod( X, chol2inv(chol(crossprod(Xu))) ), Xu )
	}
	
	# determine SNP blocks	
	# number of snps to be used
	nloci <- length(snp.include)
	message(paste("Running Analysis with",nloci,"SNPs ..."))
	# number of blocks of snps
	nblocks <- ceiling(nloci/block.size)
	# start and end positions for blocks
	if(nblocks==1){
		snp.start <- 1
		snp.end <- nloci
	}else{
		snp.start <- (0:(nblocks-1))*block.size+1
		snp.end <- c((1:(nblocks-1))*block.size, nloci)
	}
	
	
	# results matrix
	message("Creating data.frames to Store Results...")
	if(ibd.prop){
		res <- data.frame(ID1 = rep(scan.include[-nsamp], times=((nsamp-1):1)),
                                  ID2 = unlist(lapply(2:nsamp,function(x){ scan.include[x:nsamp] })),
                                  nsnp = rep(0, nsamp*(nsamp-1)/2),
                                  kin = rep(0, nsamp*(nsamp-1)/2),
                                  kindenom = rep(0, nsamp*(nsamp-1)/2),
                                  k2 = rep(0, nsamp*(nsamp-1)/2),
                                  k2denom = rep(0, nsamp*(nsamp-1)/2),
                                  k0 = rep(0, nsamp*(nsamp-1)/2),
                                  k0denom = rep(0, nsamp*(nsamp-1)/2))		
	}else{
		res <- data.frame(ID1 = rep(scan.include[-nsamp], times=((nsamp-1):1)),
                                  ID2 = unlist(lapply(2:nsamp,function(x){ scan.include[x:nsamp] })),
                                  nsnp = rep(0, nsamp*(nsamp-1)/2),
                                  kin = rep(0, nsamp*(nsamp-1)/2),
                                  kindenom = rep(0, nsamp*(nsamp-1)/2))
	}
	inbreed <- data.frame(ID = scan.include,
                              nsnp = rep(0,nsamp),
                              f = rep(0,nsamp),
                              fdenom = rep(0,nsamp))
							
	
	for(bn in 1:nblocks){
		message(paste("Computing Individual Specific Allele Frequencies: Chromosome",chromosome,"- Block",bn,"of",nblocks,"..."))
		# load genotype data
		geno <- snpgdsGetGeno(genoData, sample.id = scan.include, snp.id = snp.include[snp.start[bn]:snp.end[bn]], verbose = FALSE)
        # transpose if samples stored as first dimension
		if(!snpfirstdim){
                  geno <- t(geno)
		}
		# set any negative or "3" genotype values to NA
                geno[geno < 0 | geno==3] <- NA
		
		# minor allele frequency filter
		if(is.null(unrel.set)){
			# allele freq from entire sample
			pA <- 0.5*rowMeans(geno, na.rm = TRUE)
		}else{
			# allele freq from unrelated set
			pA <- 0.5*rowMeans(geno[,unrel.set], na.rm = TRUE)
		}
		# remove SNPs with low MAF
		snp.excl <- which(is.na(pA) | pA < MAF | pA > (1-MAF))
		if(length(snp.excl > 0)){
			geno <- geno[-snp.excl,]
			pA <- pA[-snp.excl]
		}
		# number of snps in the block
		nsnp.block <- dim(geno)[1]
		
		# rather than compute new H for each SNP, impute missing genotype values to sample mean
		# which genotype values are missing
		miss.idx <- which(is.na(geno))
		# if there are missing genotype values
		if(length(miss.idx) > 0){
			# index of which snps have the missing values
			snp.idx <- miss.idx %% nsnp.block; snp.idx[snp.idx==0] <- nsnp.block
			# replace missing genotypes with twice the sample allele frequency
			geno[miss.idx] <- 2*pA[snp.idx]
		}
		
		# matrix of individual specific allele frequencies
		if(is.null(unrel.set)){
			phat <- 0.5*tcrossprod(geno,H)
		}else{
			phat <- 0.5*tcrossprod(geno[,unrel.set],H)
		}
                # opposite allele freq
                qhat <- 1-phat
                
		# plug in a temporary value where genotypes are missing (for filtering)
		phat[miss.idx] <- -1
		# which snps to filter (fitted value too big/small)
		filt.idx <- which(phat < MAF | phat > (1-MAF))
		
		# set missing & filtered values to 0 (so no contribution)
		geno[filt.idx] <- 0
		phat[filt.idx] <- 0
        qhat[filt.idx] <- 0
        
        # product (filtered values already 0)
        phatqhat <- phat*qhat
		
		if(ibd.prop){			
			# update inbreeding denominator 			
			inbreed$fdenom <- inbreed$fdenom + colSums(phatqhat)
			
			message(paste("Computing Dominance Covariance Matrix: Block",bn,"of",nblocks,"..."))			
			# create dominance coded matrix			
			Xd <- matrix(0, nrow = nsnp.block, ncol = nsamp)
			idx0 <- which(geno == 0)
			idx2 <- which(geno == 2)
			Xd[idx0] <- phat[idx0]
			Xd[idx2] <- qhat[idx2]
			Xd <- Xd - phatqhat
			# set missing & filtered SNPs to 0 (so no contribution)
			Xd[filt.idx] <- 0
			
			# updated inbreeding numerator
			inbreed$f <- inbreed$f + colSums(Xd)
			
			# update numerator
			k2Mat <- crossprod(Xd); rm(Xd)
			res$k2 <- res$k2 + k2Mat[lower.tri(k2Mat)]; rm(k2Mat)
			
			# update denominator
			sumpqpqMat <- crossprod(phatqhat)
			res$k2denom <- res$k2denom + sumpqpqMat[lower.tri(sumpqpqMat)]; rm(sumpqpqMat)	

			
			message(paste("Computing Count of Opposite Homozygotes: Block",bn,"of",nblocks,"..."))
			# update denominator
			sump2q2Mat <- crossprod(phat^2,qhat^2)
			res$k0denom <- res$k0denom + sump2q2Mat[lower.tri(sump2q2Mat)] + t(sump2q2Mat)[lower.tri(sump2q2Mat)] ; rm(sump2q2Mat)
			
			# homozygotes
			IAA <- geno == 2
			Iaa <- geno == 0
			# set missing & filtered SNPs to FALSE (so no contribution)
			IAA[filt.idx] <- FALSE; Iaa[filt.idx] <- FALSE
						
			# update numerator
			NAAaaMat <- crossprod(IAA,Iaa); rm(IAA); rm(Iaa)
			res$k0 <- res$k0 + NAAaaMat[lower.tri(NAAaaMat)] + t(NAAaaMat)[lower.tri(NAAaaMat)]; rm(NAAaaMat)
			
		}else{
			# update inbreeding numerator
			inbreed$f <- inbreed$f + colSums(geno^2) - colSums(geno) - 2*colSums(phat*geno) + 2*colSums(phat^2)
			
			# update inbreeding denominator (filtered values already 0)
			inbreed$fdenom <- inbreed$fdenom + 2*colSums(phatqhat)
		}
		
		message(paste("Computing Additive Covariance Matrix: Block",bn,"of",nblocks,"..."))
		# matrix of regression residuals (filtered values already 0)
		R <- geno-2*phat; rm(geno)
		
		# update numerator
		kinMat <- crossprod(R); rm(R)
		res$kin <- res$kin + kinMat[lower.tri(kinMat)]; rm(kinMat)	
		
		# update denominator
		sumsqrtpqpqMat <- crossprod(sqrt(phatqhat)); rm(phatqhat)
		res$kindenom <- res$kindenom + sumsqrtpqpqMat[lower.tri(sumsqrtpqpqMat)]; rm(sumsqrtpqpqMat)

			
		# determine the number of snps used for each pair
		snpcount <- matrix(1, nrow=nsnp.block, ncol=nsamp); snpcount[filt.idx] <- 0
		nsnpMat <- crossprod(snpcount); rm(snpcount)
		res$nsnp <- res$nsnp + nsnpMat[lower.tri(nsnpMat)]; 
		inbreed$nsnp <- inbreed$nsnp + diag(nsnpMat); rm(nsnpMat)
	}
	
	# return results
	out <- list(kinship = res,
                    inbreed = inbreed,
                    call = match.call(),
                    method = method)
	class(out) <- "pcrelateChr"
	return(out)
}

	
pcrelateChrCombine <- function(filelist, pcMat = NULL, ibd.prop = TRUE, correct = TRUE){
	
	# PC checks
	if(is.null(pcMat)){
		message("pcMat not specified, Calculating Unadjusted Relatedness Estimates")
		method <- "Unadjusted"
		
		# load the first file in filelist
		message("Loading Results from Chromosome 1...")
		tmp <- get(load(filelist[1]))
		# check the class
		if(class(tmp) != "pcrelateChr"){
			stop("The objects in filelist should be class 'pcrelateChr'.")
		}
		
		# sample size
		nsamp <- dim(tmp$inbreed)[1]
				
	}else{
		pcMat <- as.matrix(pcMat)
		message(paste("Adjusted for",ncol(pcMat),"PC(s)"))
		method <- "PC-Relate"
		
		# load the first file in filelist
		message("Loading Results from Chromosome 1...")
		tmp <- get(load(filelist[1]))
		# check the class
		if(class(tmp) != "pcrelateChr"){
			stop("The objects in filelist should be class 'pcrelateChr'.")
		}
		
		# sample size
		nsamp <- dim(tmp$inbreed)[1]
		
		# check
		if(nsamp != nrow(pcMat)){
			stop("The number of Samples used in genoData and pcMat do not match.")
		}
	}
	
	# pull out results
	combined.kinship <- tmp$kinship
	combined.inbreed <- tmp$inbreed

	# load the rest of the files
	if(length(filelist) > 1){
		for(i in 2:length(filelist)){
			# load
			message(paste("Loading Results from Chromosome",i,"..."))
			tmp <- get(load(filelist[i]))
			# check the class
			if(class(tmp) != "pcrelateChr"){
				stop("The objects in filelist should be class 'pcrelateChr'.")
			}
			
			# check that IDs match
			if(!all(combined.kinship[,1:2] == tmp$kinship[,1:2]) | !all(combined.inbreed[,1] == tmp$inbreed[,1])){
				stop("The Individual IDs in these files do not match!")
			}
			
			# update results
			combined.kinship[,-(1:2)] <- combined.kinship[,-(1:2)] + tmp$kinship[,-(1:2)]
			combined.inbreed[,-1] <- combined.inbreed[,-1] + tmp$inbreed[,-1]
		}
	}
	
	message("Computing Final Estimates...")
	# calculate inbreeding estimates
	combined.inbreed$f <- combined.inbreed$f/combined.inbreed$fdenom
	# remove denominator
	combined.inbreed <- combined.inbreed[,-4]
	
	# calculate kinship estimates
	combined.kinship$kin <- combined.kinship$kin/(4*combined.kinship$kindenom)
	# remove denominator
	combined.kinship <- combined.kinship[,-5]

	kincorrect <- NULL
	fcorrect <- NULL
	if(!is.null(pcMat) & correct){
		# correct for possible overadjustment in each PC
		idxf <- which(combined.inbreed$f < 2^(-11/2))
		for(i in 1:ncol(pcMat)){
			Avec <- pcMat[,i]			
			fvals <- lm(combined.inbreed$f[idxf] ~ I(Avec[idxf]) + I(Avec[idxf]^2))$coef
			combined.inbreed$f <- combined.inbreed$f - fvals[1] - fvals[2]*Avec - fvals[3]*Avec^2	
			fcorrect <- append(fcorrect, fvals)

			idx <- which(combined.kinship$kin < 2^(-11/2))
			Acov <- tcrossprod(Avec); Acov <- Acov[lower.tri(Acov)]
			kvals <- lm(combined.kinship$kin[idx] ~ Acov[idx])$coef
			combined.kinship$kin <- combined.kinship$kin - kvals[1] - kvals[2]*Acov; rm(Acov); rm(Avec)
			kincorrect <- append(kincorrect, kvals)
		}
	}
	
	k2correct <- NULL
	if(ibd.prop){
		# pairwise inbreeding product
		fprod <- tcrossprod(combined.inbreed$f)
		fprod <- fprod[lower.tri(fprod)]
		
		# calculate k2 estimates
		combined.kinship$k2 <- combined.kinship$k2/combined.kinship$k2denom - fprod
		
		# correction
		if(correct){
			if(!(is.null(pcMat))){
				for(i in 1:ncol(pcMat)){
					Acov <- tcrossprod(pcMat[,i]); Acov <- Acov[lower.tri(Acov)]
					idx <- which(combined.kinship$kin < 2^(-11/2))
					k2vals <- lm(combined.kinship$k2[idx] ~ I(Acov[idx]) + I(Acov[idx]^2))$coef
					combined.kinship$k2 <- combined.kinship$k2 - k2vals[1] - k2vals[2]*Acov - k2vals[3]*Acov^2; rm(Acov)
					k2correct <- append(k2correct, k2vals)
				}
			}
			idx <- which(combined.kinship$k2 < 2^(-9/2))
			k2vals <- lm(combined.kinship$k2[idx] ~ combined.kinship$kin[idx])$coef
			combined.kinship$k2 <- combined.kinship$k2 - k2vals[1] - k2vals[2]*combined.kinship$kin
			k2correct <- append(k2correct, k2vals)
		}
		
		# calculate k0 estimates
		combined.kinship$k0 <- combined.kinship$k0/combined.kinship$k0denom
		# remove denominator
		combined.kinship <- combined.kinship[,-8]

		# index for not PO/FS
		rel2idx <- which(combined.kinship$kin < 2^(-5/2))
		combined.kinship$k0[rel2idx] <- (1 - 4*combined.kinship$kin[rel2idx] + combined.kinship$k2[rel2idx])
		
		# calculate k1 estimates
		names(combined.kinship)[6] <- "k1"
		combined.kinship$k1 <- 1 - combined.kinship$k0 - combined.kinship$k2
	}	
	
	# return results
	out <- list(kinship = combined.kinship,
                    inbreed = combined.inbreed,
                    kincorrect = kincorrect,
                    fcorrect = fcorrect,
                    k2correct = k2correct,
                    call = match.call(),
                    method = method)	
	class(out) <- "pcrelate"
	return(out)
}
