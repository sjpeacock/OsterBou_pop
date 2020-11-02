###############################################################################
## Convergence diagnostics
###############################################################################
Isigma <- TRUE

includeChains <- list();length(includeChains) <- dim(models2test)[1]*6
dim(includeChains) <- c(dim(models2test)[1], 6)

for(m in 1:dim(models2test)[1]){
	# Model 1
includeChains[[m,1]] <- c(1:6)
includeChains[[m,2]] <- c(1:6)
includeChains[[m,3]] <- c(1:6)
includeChains[[m,4]] <- c(1:6)
includeChains[[m,5]] <- c(1:6)
includeChains[[m,6]] <- c(1:6)
}

if(Isigma == TRUE){
	models2test[, 'nParams'] <- c(28, 13, 15, 12, 14, 15, 17, 14, 16)
	
	includeChains[[3,1]] <- c(1:3, 5)
	includeChains[[3,4]] <- c(1:2, 5)
	includeChains[[3,5]] <- c(1, 4:6)
	includeChains[[3,6]] <- c(2:4, 6)
	
	includeChains[[5,6]] <- c(2:5)
	includeChains[[5,5]] <- c(1, 3:5)
	includeChains[[5,4]] <- c(1:4, 6)
	includeChains[[5,2]] <- c(2:6)
	
	includeChains[[6,1]] <- c(1, 5, 6)
	includeChains[[6,2]] <- c(2:4)
	includeChains[[6,3]] <- c(1:2, 5:6)
	includeChains[[6,4]] <- c(1, 3, 5, 6)
	includeChains[[6,5]] <- c(1:3, 5, 6)
	includeChains[[6,6]] <- c(1, 3:6)
	
	includeChains[[7,1]] <- c(2, 4:6)
	includeChains[[7,2]] <- c(1,2,5,6)
	includeChains[[7,3]] <- c(2:4)
	includeChains[[7,4]] <- c(2:3, 6)
	includeChains[[7,5]] <- c(2, 5, 6)
	includeChains[[7,6]] <- c(1, 3, 6)
	
	includeChains[[8,1]] <- c(1,3,4)
	includeChains[[8,2]] <- c(3,4,5)
	includeChains[[8,3]] <- c(1,2,5)
	includeChains[[8,4]] <- c(2,4,5,6)
	includeChains[[8,5]] <- c(3,4,6)
	includeChains[[8,6]] <- c(2:4)
	
	includeChains[[9,1]] <- c(2:3)
	includeChains[[9,2]] <- c(2:3)
	includeChains[[9,3]] <- c(2:6)
	includeChains[[9,4]] <- c(1,2,6)
	includeChains[[9,5]] <- c(1:2) #* didn't converge at all
	includeChains[[9,6]] <- c(3:5)
	
} else {
includeChains[[3,1]] <- c(1:3, 5)
includeChains[[3,2]] <- c(1:5)
includeChains[[3,3]] <- c(1:5)
includeChains[[3,5]] <- c(2:5)

includeChains[[6,1]] <- c(1,3,4,6)
includeChains[[6,2]] <- c(2,3,5)
includeChains[[6,3]] <- c(1, 3:6)
includeChains[[6,4]] <- c(2:6)
includeChains[[6,5]] <- c(1, 4:6)
includeChains[[6,6]] <- c(1,2, 4:6)

includeChains[[7,1]] <- c(2,4,6)
includeChains[[7,2]] <- c(1,3,5,6)
includeChains[[7,3]] <- c(2,4,6)
includeChains[[7,4]] <- c(3,4,6)
includeChains[[7,5]] <- c(1:4, 6)
includeChains[[7,6]] <- c(1, 4:5)

includeChains[[8,1]] <- c(1,3,6)
includeChains[[8,2]] <- c(1,2,4,6)
includeChains[[8,3]] <- c(2:5)
includeChains[[8,4]] <- c(1,2,5,6)
includeChains[[8,6]] <- c(1,3,4)

includeChains[[9,1]] <- c(3,4,6)
includeChains[[9,2]] <- c(1:5)
includeChains[[9,3]] <- c(1,2,4:6)
includeChains[[9,4]] <- c(2,4:6)
includeChains[[9,5]] <- c(1,2,5)
}

# List with chains that converged to include
in.all <- list();length(in.all) <- dim(models2test)[1]*6
dim(in.all) <- c(dim(models2test)[1], 6)

summary.in <- in.all

meanEst.all <- list(); length(meanEst.all) <- dim(models2test)[1]
Rhat.all <- list(); length(Rhat.all) <- dim(models2test)[1]

for(m in 1:dim(models2test)[1]){
	
	Rhat.all[[m]] <- matrix(nrow = as.numeric(models2test[m, 'nParams']), ncol = 6)
	meanEst.all[[m]] <- matrix(nrow = as.numeric(models2test[m, 'nParams']), ncol = 6)
	
	for(K in 1:6){
		
		includePar <- rep(NA, dim(out.all[[m, K]][[1]])[2])
		
		if(m == 1){
			if(Isigma == TRUE){
				includePar[c(29:49, 162:168)] <- 1
			} else {
					includePar[c(29:49, 162)] <- 1
			}
			
		}else {
			for(i in 1:dim(out.all[[m, K]][[1]])[2]){
				dummy <- strsplit(colnames(out.all[[m, K]][[1]])[i], split = NULL)
				
				if(dummy[[1]][1] == "p" | dummy[[1]][1] == "N") includePar[i] <- 0 else includePar[i] <- 1
			}}
		
		parIndex <- which(includePar*as.numeric(apply(out.all[[m, K]][[1]], 2, mean) != 0) == 1)
		nP <- length(parIndex)
		
		# Summary of model output, including parameters
		in.all[[m, K]] <- list(); length(in.all[[m, K]]) <- length(includeChains[[m, K]])
		for(i in 1:length(includeChains[[m, K]])){
			in.all[[m, K]][[i]] <- as.mcmc(matrix(
				out.all[[m, K]][[includeChains[[m, K]][i]]][, parIndex], 
				nrow = numIter, 
				ncol = nP, 
				dimnames = list(NULL, colnames(out.all[[m, K]][[1]])[parIndex])))
		}
		
		in.all[[m, K]] <- mcmc.list(in.all[[m, K]])
		
		X <- summary(in.all[[m, K]])
			summary.in[[m, K]] <- X
		
		
		meanEst.all[[m]][, K] <- X[[1]][, 1]
		
		Rhat.all[[m]][, K] <- gelman.diag(in.all[[m, K]])[[1]][, 1]
		
		ind <- round(seq(1, dim(out.all[[m, K]][[1]])[1], length.out = 300))
		
		pdf(width = 8.5, height = 11, pointsize = 10, file = paste("dataAnalysis/plots/Isigma/TracePlots_model", m, "_repOut", K, ".pdf", sep=""))
		par(mfrow = c(4, 3), mar = c(3,2,1,1), oma = c(4,4,2,2))
		I <- 0
		for(i in 1:nP){
			plot(ind, in.all[[m, K]][[1]][ind, i], "l", col = c(1:6)[includeChains[[m, K]][1]], lty = 3, xlab = "", ylab = "", main = rownames(X[[1]])[i], ylim = X[[2]][i, c(1, 5)])
			
			for(j in 2:length(includeChains[[m, K]])){
				lines(ind, in.all[[m, K]][[j]][ind, i], col = c(1:6)[includeChains[[m, K]][j]], lty = 3)
			}
			mtext(side = 3, adj = 0, line = -2, paste("  Rhat = ", round(Rhat.all[[m]][i, K], 2)))
			}
		
		plot(1,1, "n", xaxt="n", yaxt="n", bty = "n", xlab="", ylab="")
		legend(1, 1, legend = 1:6, lwd = 1, col = 1:6, title = "chains")
		dev.off()
		
	} # end K
	} # end m


###############################################################################
## How do parameter values differ among models and among leave-out reps?
###############################################################################

tempInd <- matrix(exp(summary.out[[1]][parIndex[1:21], 1]), nrow = 7, ncol = 3, byrow = FALSE)
