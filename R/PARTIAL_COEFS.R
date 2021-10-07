

PARTIAL_COEFS <- function(cormat, modelRsq=NULL, verbose=TRUE) {

	# for cormat, the DV must be in column 1

	NVARs <- nrow(cormat) #- 1
	
	# standardized coefficients
	betas <- solve(cormat[2:NVARs,2:NVARs]) %*% cormat[2:NVARs,1]  
		
	# partial correlations
	partials <- matrix(-9999,(NVARs-1),1)
	Rinv <- solve(cormat)
	for (luper in 2:NVARs) {
		partials[(luper-1),1] <- Rinv[luper,1] / ( -1 * sqrt(Rinv[1,1] * Rinv[luper,luper]) )
	}
		
	# semi-partial (part) correlations
	# Cohen et al 2003 - Applied Multiple Regression/Correlation Analysis for the Behavioral Sciences p. 102
	
	if (is.null(modelRsq))
		modelRsq <- t(solve(cormat[2:NVARs,2:NVARs]) %*% cormat[2:NVARs,1]) %*% cormat[2:NVARs,1]
	
	semipartials <- sqrt( ( partials^2 / (1 - partials^2) ) %*% (1 - modelRsq) )
	semipartials <- semipartials * sign(partials)

	
	Rcoefs <- data.frame(beta = betas, r = cormat[2:nrow(cormat),1],
                         partial.r = partials, semipartial.r = semipartials )

	if (verbose) print(Rcoefs)


return(invisible(Rcoefs))

}


