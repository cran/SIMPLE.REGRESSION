


# rounds numeric columns in a matrix
# numeric columns named 'p' or 'plevel' are rounded to round_p places
# numeric columns not named 'p' are rounded to round_non_p places

round_boc <- function(donnes, round_non_p = 2, round_p = 5) {
	
	# identify the numeric columns
	#	numers <- apply(donnes, 2, is.numeric)  # does not work consistently 
	for (lupec in 1:ncol(donnes)) {

		if (is.numeric(donnes[,lupec]) == TRUE) 
		
			if (colnames(donnes)[lupec] == 'p' | colnames(donnes)[lupec] == 'plevel' | colnames(donnes)[lupec] == 'Pr(>|t|)')  {
				donnes[,lupec] <- round(donnes[,lupec],round_p)
			} else {
				donnes[,lupec] <- round(donnes[,lupec],round_non_p)				
			}		
		# if (is.numeric(donnes[,lupec]) == FALSE) numers[lupec] = 'FALSE'		
		# if (colnames(donnes)[lupec] == 'p') numers[lupec] = 'FALSE'		
	}
	
	# # set the p column to FALSE
	# numers_not_p <- !names(numers) %in% "p"
	
#	donnes[,numers_not_p] = round(donnes[,numers_not_p],round_non_p) 
	
#	if (any(colnames(donnes) == 'p'))  donnes[,'p'] = round(donnes[,'p'],round_p) 

	return(invisible(donnes))
}










# Collinearity Diagnostics

Collinearity <- function(modelmat, verbose=FALSE) {
	

	# Variance Inflation Factors 
	
	cormat <- cor(modelmat[,-1,drop=FALSE])  # ignoring the intercept column
	
	smc <- 1 - (1 / diag(solve(cormat)))
	
	VIF <- 1 / (1 - smc)
	
	Tolerance <- 1 / VIF
	
	VIFtol <- cbind(Tolerance, VIF)
	
	# Eigenvalues, Condition Indices, Condition Number and Condition Index, & Variance Decomposition Proportions (as in SPSS)
	
	# adapted from https://stat.ethz.ch/pipermail/r-help/2003-July/036756.html
		
	modelmat <- scale(modelmat, center=FALSE) / sqrt(nrow(modelmat) - 1)
	
	svd.modelmat <- svd(modelmat)
	
	singular.values <- svd.modelmat$d
	
	evals <- singular.values **2
	
	condition.indices <- max(svd.modelmat$d) / svd.modelmat$d
	
	phi <- sweep(svd.modelmat$v**2, 2, svd.modelmat$d**2, "/")
	
	VP <- t(sweep(phi, 1, rowSums(phi), "/"))
	
	colnames(VP) <- colnames((modelmat))
	rownames(VP) <- 1:nrow(VP)
	
		
	CondInd <- data.frame(1:length(evals), cbind(evals, condition.indices, VP))	
	colnames(CondInd)[1:4] <- c('Dimension','Eigenvalue','Condition Index', 'Intercept')
	
	if (verbose) {
		message('\n\nCollinearity Diagnostics:\n')

#		message('\n\nTolerance and Variance Inflation Factors:\n')
		print(round(VIFtol, 3), print.gap=4)
		message('\nThe mean Variance Inflation Factor = ', round(mean(VIFtol[,2]),3))
		message('\nMulticollinearity is said to exist when VIF values are > 5, or when Tolerance values')
		message('are < 0.1. Multicollinearity may be biasing a regression model, and there is reason')
		message('for serious concern, when VIF values are greater than 10.\n\n')

		print(round(CondInd,3), print.gap=4)
		message('\nThe coefficients for the Intercept and for the predictors are the Variance Proportions.')

		message('\nEigenvalues that differ greatly in size suggest multicollinearity. Multicollinearity is')
		message('said to exist when the largest condition index is between 10 and 30, and evidence for')
		message('the problem is strong when the largest condition index is > 30. Multicollinearity is')
		message('also said to exist when the Variance Proportions (in the same row) for two variables are')
		message('above .80 and when the corresponding condition index for a dimension is higher than 10 to 30.\n\n')
	}	

simple.Model <- list(VIFtol=VIFtol, CondInd=CondInd)

return(invisible(simple.Model))


}

