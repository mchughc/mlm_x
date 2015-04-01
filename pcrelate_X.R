# define method for GWASTools GenotypeData

pcrelate <- function(genoData, MAF = 0.05, pcMat = NULL, unrel.set = NULL, snp.include = NULL, scan.include = NULL, Xchr = FALSE, block.size = 10000, correct = TRUE){
	
	# checks
	if(MAF < 0 | MAF >= 0.5){
		stop("MAF must be > 0 and  < 0.5")
	}
    
    # if X chromosome, check that sex is available
    #if(Xchr){
    #    if(!hasSex(genoData)){
    #        stop("sex must be provided in genoData to analyze X chromosome SNPs")
    #    }
    #}
	
	# snps to include
	if(!is.null(snp.include)){
		if(!all(is.element(snp.include,getSnpID(genoData)))){
			stop("Not all of the SNP IDs in snp.include are in genoData")
		}
    }else{
        if(!Xchr){
            # all autosomal SNPs
            snp.include <- getSnpID(genoData)[getChromosome(genoData) <= 22]
        }else{
            # all X chromosome SNPs
            snp.include <- getSnpID(genoData)[getChromosome(genoData) == XchromCode(genoData)]
        }
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
    res <- data.frame(ID1 = rep(scan.include[-nsamp], times=((nsamp-1):1)),
                        ID2 = unlist(lapply(2:nsamp,function(x){ scan.include[x:nsamp] })),
                        nsnp = rep(0, nsamp*(nsamp-1)/2),
                        kin = rep(0, nsamp*(nsamp-1)/2))
    kindenom = rep(0, nsamp*(nsamp-1)/2)
    
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
	
        # update inbreeding denominator (filtered values already 0)
        fdenom <- fdenom + 2*colSums(phatqhat)
				
		message(paste("Computing Additive Covariance Matrix: Block",bn,"of",nblocks,"..."))
		# matrix of regression residuals (filtered values already 0)
		R <- geno-2*phat; rm(geno)
		
		# update numerator
		kinMat <- crossprod(R); rm(R)
		res$kin <- res$kin + kinMat[lower.tri(kinMat)]
        
        # update inbreeding numerator
        inbreed$f <- inbreed$f + diag(kinMat); rm(kinMat)
		
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
	inbreed$f <- inbreed$f/fdenom - 1; rm(fdenom)
	
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

    # return results
	out <- list(kinship= res,
                    inbreed = inbreed,
                    kincorrect = kincorrect,
                    fcorrect = fcorrect,
                    call = match.call(),
                    method = method)	
	class(out) <- "pcrelate"
	return(out)
}


