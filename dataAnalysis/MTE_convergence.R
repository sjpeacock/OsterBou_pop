###############################################################################
## Convergence diagnostics
###############################################################################

for(m in 1:dim(models2test)[1]){
	for(K in 1:5){
		#	m <- 2
		# K <- 1
		
		# Summary of model output, including parameters
		summary.out <- summary(out.all[[m, K]])
		
		includePar <- rep(NA, dim(summary.out[[1]])[1])
		for(i in 1:dim(summary.out[[1]])[1]){
			dummy <- strsplit(rownames(summary.out[[1]])[i], split = NULL)
			if(dummy[[1]][1] == "p" | dummy[[1]][1] == "N") includePar[i] <- 0 else includePar[i] <- 1
		}
		
		parIndex <- which(includePar == 1 & summary.out[[1]][, 2] != 0)
		nP <- length(parIndex)
		summary.out[[1]][parIndex, ]
		
		Rhat <- gelman.diag(as.mcmc(list(out.all[[m, K]][[1]][, parIndex], out.all[[m, K]][[2]][, parIndex], out.all[[m, K]][[3]][, parIndex], out.all[[m, K]][[4]][, parIndex])))
		
		pdf(width = 8.5, height = 11, pointsize = 10, file = paste("plots/TracePlots_model", m, "_repOut", K, ".pdf", sep=""))
		par(mfrow = c(4, 3), mar = c(1,2,1,1), oma = c(4,4,2,2))
		I <- 0
		for(i in parIndex){
			I <- I + 1
			plot(seq(1, numIter, 10), out.all[[m, K]][[1]][seq(1, numIter, 10), i], "l", col = 1, lty = 3, xlab = "", ylab = "", main = rownames(summary.out[[1]])[i], xaxt="n")
			for(j in 2:4) lines(seq(1, numIter, 10), out.all[[m, K]][[j]][seq(1, numIter, 10), i], col = j, lty = 3)
			mtext(side = 3, adj = 0, line = -2, paste("  Rhat = ", round(Rhat[[1]][I, 1], 2)))
		}
		dev.off()
		
	}}


###############################################################################
## How do parameter values differ among models and among leave-out reps?
###############################################################################

tempInd <- matrix(exp(summary.out[[1]][parIndex[1:21], 1]), nrow = 7, ncol = 3, byrow = FALSE)
