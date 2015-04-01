pcair <-
function(genoData, v = 10, MAF = 0.05, kinMat = NULL, kin.thresh = 0.025, divMat = NULL, div.thresh = -0.025, unrel.set = NULL, snp.include = NULL, Xchr = FALSE, block.size = 10000){
	
	# checks
	if(MAF < 0 | MAF > 0.5){
		stop("MAF must be between 0 and 0.5")
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
    
    # check that scanID matches kinMat
    if(!is.null(kinMat)){
        if(is.null(colnames(kinMat))){
            stop("colnames and rownames of kinMat must be individual IDs")
        }
        if(!all(scan.include == colnames(kinMat))){
            stop("colnames and rownames of kinMat must match the scanIDs of genoData")
        }
    }
    # check that scanID matches unrel.set
    if(!is.null(unrel.set)){
        if(!all(unrel.set %in% scan.include)){
            stop("All of the samples in unrel.set must be in the scanIDs of genoData")
        }
    }
	
	# get related and unrelated sets
	if(is.null(kinMat)){
		if(is.null(unrel.set)){
			message("kinMat and unrel.set both unspecified, Running Standard Principal Components Analysis")
			rels <- NULL
			unrels <- scan.include
		}else{
			message("kinMat not specified, using unrel.set as the Unrelated Set")
			rels <- scan.include[!(scan.include %in% unrel.set)]
			unrels <- unrel.set
		}
	}else{
		if(dim(kinMat)[1] != nsamp | dim(kinMat)[2] != nsamp){
			stop("The dimension of kinMat must match the number of samples in genoData")
		}
		if(!is.null(divMat)){
			if(dim(divMat)[1] != nsamp | dim(divMat)[2] != nsamp){
				stop("The dimension of divMat must match the number of samples in genoData")
			}
		}
		if(is.null(unrel.set)){
			message("Partitioning Samples into Related and Unrelated Sets...")
		}else{
			message("Partitioning Samples into Related and Unrelated Sets, unrel.set forced into the Unrelated Set")
		}
		part <- pcairPartition(kinMat = kinMat, kin.thresh = kin.thresh, divMat = divMat, div.thresh = div.thresh, unrel.set = unrel.set)
		rels <- part$rels
		unrels <- part$unrels
		if(is.null(rels)){
			message("No relatives identified, Running Standard Principal Components Analysis")
		}
	}
	
	# create index for related and unrealted sets
	rel.idx <- which(scan.include %in% rels)
	unrel.idx <- which(scan.include %in% unrels)
	nr <- length(rel.idx)
	nu <- length(unrel.idx)
	if(nr > 0){
		method <- "PC-AiR"
	}else{
		method <- "Standard PCA"
	}
	message(paste("Unrelated Set:",nu,"Samples \nRelated Set:",nr,"Samples"))
	
	# determine blocks
	# number of autosomal snps
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
	
	# if relatives
	if(nr > 0){
		# correlation matrix for sPCA
		Psiu <- matrix(0, nrow=nu, ncol=nu)
		# number of snps used
		nsnps <- 0
		
		for(bn in 1:nblocks){
			message(paste("Computing Genetic Correlation Matrix for the Unrelated Set: Block",bn,"of",nblocks,"..."))
			# load genotype data for unrelated set
			geno <- getGenotype(genoData, scan=c(1,-1), snp=c(snp.start[bn], snp.end[bn]-snp.start[bn]+1) )
			# subset included snps
			geno <- geno[is.element(getSnpID(genoData, snp.start[bn]:snp.end[bn]), snp.include),]
            # remove exclude samples
            if(!is.null(scan.include)){
                geno <- geno[,scan.include]
            }
			# subset unrelated set
			geno <- geno[,unrel.idx]
			
			# allele freq est from unrelated set
			pA <- 0.5*rowMeans(geno, na.rm = TRUE)
			# remove monomorphic SNPs
			snp.excl <- which(is.na(pA) | pA <= MAF | pA >= (1-MAF))
			if(length(snp.excl > 0)){
				geno <- geno[-snp.excl,]
				pA <- pA[-snp.excl]
			}
			nsnps <- nsnps + length(pA)
			
			# Standardized genotype values
			# estimated variance at each SNP
			sigma.hat <- sqrt(2*pA*(1-pA))
			# z_i = (x_i - 2\hat{p})/sqrt{2*\hat{p}(1-\hat{p})}
			Zu <- (geno-2*pA)/sigma.hat
			Zu[which(is.na(Zu))] <- 0
			
			# unrelated empirical correlation matrix
			Psiu <- Psiu + crossprod(Zu)
		}
		Psiu <- (1/nsnps)*Psiu
        
		
		# sPCA analysis
		message("Performing PCA on the Unrelated Set...")
		eigu <- eigen(Psiu, symmetric=TRUE)
		
		# subset desired number of eigenvectors
		if(is.null(v)){
			v <- nu
		}
		L <- eigu$values[1:v]
		V <- eigu$vectors[,1:v]
		# sum of eigenvalues
		sum.values <- sum(eigu$values)
		
		# matrix of pseudo-eigenvectors
		Q <- matrix(0, nrow=nr, ncol=v)
		
		# project for related set
		message("Predicting PC Values for the Related Set...")
		for(bn in 1:nblocks){
			# load genotype data
			geno <- getGenotype(genoData, scan=c(1,-1), snp=c(snp.start[bn], snp.end[bn]-snp.start[bn]+1) )
			# subset included snps
			geno <- geno[is.element(getSnpID(genoData, snp.start[bn]:snp.end[bn]), snp.include),]
            # remove exclude samples
            if(!is.null(scan.include)){
                geno <- geno[,scan.include]
            }
			# allele freq est from unrelated set
			pA <- 0.5*rowMeans(geno[,unrel.idx], na.rm = TRUE)
			# remove monomorphic SNPs
			snp.excl <- which(is.na(pA) | pA <= MAF | pA >= (1-MAF))
			if(length(snp.excl > 0)){
				geno <- geno[-snp.excl,]
				pA <- pA[-snp.excl]
			}
			
			# Standardized genotype values
			# estimated variance at each SNP
			sigma.hat <- sqrt(2*pA*(1-pA))
			# z_i = (x_i - 2\hat{p})/sqrt{2*\hat{p}(1-\hat{p})}
			Z <- (geno-2*pA)/sigma.hat
			Z[which(is.na(Z))] <- 0
			
			# subset unrelated and related sets
			Zu <- Z[,unrel.idx]
			Zr <- Z[,rel.idx]
			
			# SNP weights
			WT <- tcrossprod(t(V),Zu)
			WL <- (1/L)*WT
			
			# pseudo eigenvectors for relateds
			Q <- Q + crossprod(Zr,t(WL))
		}
		Q <- (1/nsnps)*Q
		
		# concatenate
		message("Concatenating Results...")
		EIG <- matrix(NA, nrow=nsamp, ncol=v)
		EIG[unrel.idx,] <- V
		EIG[rel.idx,] <- Q
		
	# if no relatives
	}else{
		# correlation matrix for sPCA
		Psi <- matrix(0, nrow=nsamp, ncol=nsamp)
		# number of snps used
		nsnps <- 0
		
		for(bn in 1:nblocks){
			message(paste("Computing Genetic Correlation Matrix: Block",bn,"of",nblocks,"..."))
			# load genotype data
			geno <- getGenotype(genoData, scan=c(1,-1), snp=c(snp.start[bn], snp.end[bn]-snp.start[bn]+1) )
			# subset included snps
			geno <- geno[is.element(getSnpID(genoData, snp.start[bn]:snp.end[bn]), snp.include),]
            # remove exclude samples
            if(!is.null(scan.include)){
                geno <- geno[,scan.include]
            }
			# allele freq est from entire sample
			pA <- 0.5*rowMeans(geno, na.rm = TRUE)
			# remove monomorphic SNPs
			snp.excl <- which(is.na(pA) | pA <= MAF | pA >= (1-MAF))
			if(length(snp.excl > 0)){
				geno <- geno[-snp.excl,]
				pA <- pA[-snp.excl]
			}
			nsnps <- nsnps + length(pA)
			
			# Standardized genotype values
			# estimated variance at each SNP
			sigma.hat <- sqrt(2*pA*(1-pA))
			# z_i = (x_i - 2\hat{p})/sqrt{2*\hat{p}(1-\hat{p})}
			Z <- (geno-2*pA)/sigma.hat
			Z[which(is.na(Z))] <- 0
			
			# empirical correlation matrix
			Psi <- Psi + crossprod(Z)
		}
		Psi <- (1/nsnps)*Psi
		
		# sPCA analysis
		message("Performing Standard PCA...")
		eig <- eigen(Psi, symmetric=TRUE)
		
		# subset desired number of eigenvectors
		if(is.null(v)){
			v <- nsamp
		}
		# output
		EIG <- eig$vectors[,1:v]
		L <- eig$values[1:v]
		sum.values <- sum(eig$values)
	}
    
    # add scanIDs as rownames of EIG
    rownames(EIG) <- scan.include
	
	# return results
	out <- list(vectors = EIG, 
				values = L, 
				sum.values = sum.values, 
				rels = rels, 
				unrels = unrels,
				kin.thresh = kin.thresh,
				div.thresh = -abs(div.thresh),
				nsamp = nsamp,
				nsnps = nsnps,
				MAF = MAF,
				call = match.call(),
				method = method)
	class(out) <- "pcair"
	return(out)
}
