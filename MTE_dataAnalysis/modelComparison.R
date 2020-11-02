# From MTE_convergence.R, we have the following objects:
#
#     in.all[[m, K]] - list for each model/rep with only converged 
#                      chains and relevant parameters
#     summary.in[[m,K]] - list of summaries
#     meanEst.all[[m]][, K] - list for each model with mean parameter estimates
#                             for each rep
#     Rhat.all[[m]][, K] - Gelman and Rubin convergence diagnostic for each par

###############################################################################
## Likelihood
###############################################################################

#------------------------------------------------------------------------------
# Matrix of parameter values for each temperature (rows) and parameter (columns)
paramsPred.all <- list(); length(paramsPred.all) <- dim(models2test)[1]*6; dim(paramsPred.all) <- c(dim(models2test)[1], 6)

for(m in 1:dim(models2test)[1]){ 
	for(K in 1:6){
		
		paramsPred <- matrix(NA, nrow = dat$nT, ncol = 4)
		
		for(j in 1:3){ # for each parameter u0, u1, rho
			
			if(models2test[m, j] == "I"){
				paramsPred[, j] <- as.numeric(exp(summary.in[[m,K]][[1]][paste('p.mean.log[', c(1:dat$nT),",", j, "]", sep=""), 1]))
			
				} else if (models2test[m, j] == "C"){
				paramsPred[, j] <- rep(summary.in[[m, K]][[1]][paste("a[", j, "]", sep=""), 1], dat$nT)
			
				} else if (models2test[m, j] == "A"){
				params.j <- c(
					a = summary.in[[m, K]][[1]][paste("a[", j, "]", sep=""), 1], 
					E = summary.in[[m, K]][[1]][paste("E[", j, "]", sep=""), 1])
				paramsPred[, j] <- MTE.predict(model = "A", params = params.j, T.all = dat$T.obs,	exponent = c(1, 1, -1, NA)[j])
				
			} else if (models2test[m, j] == "SSU"){
				params.j <- c(
					a = summary.in[[m, K]][[1]][paste("a[", j, "]", sep=""), 1], 
					E = summary.in[[m, K]][[1]][paste("E[", j, "]", sep=""), 1], 
					Eh = summary.in[[m, K]][[1]][paste("Eh[", j, "]", sep=""), 1], 
					Th = summary.in[[m, K]][[1]][paste("Th[", j, "]", sep=""), 1])
				paramsPred[, j] <- MTE.predict(model = "SSU", params = params.j, T.all = dat$T.obs,	exponent = c(1, 1, -1, NA)[j])
			}
		} # end j
		
		if(Isigma == TRUE){
			paramsPred[, 4] <- as.numeric(exp(summary.in[[m,K]][[1]][paste('sigma.log[', c(1:dat$nT), "]", sep=""), 1]))
			} else {
				# Constant sigma
		paramsPred[, 4] <- rep(exp(summary.in[[m, K]][[1]]['sigma.log', 1]), dat$nT)
			}
		
		paramsPred.all[[m, K]] <- paramsPred
	}}

#------------------------------------------------------------------------------
# Calculate likelihood
negLogLik.all <- list(); length(negLogLik.all) <- dim(models2test)[1]*6; dim(negLogLik.all) <- c(dim(models2test)[1], 6)

for(m in 1:dim(models2test)[1]){ 
	for(K in 1:6){
		if(K ==6){
			negLogLik.all[[m, K]] <- calcLik(params.temp = paramsPred.all[[m, K]], dat = dat, holdOut = NA)
		} else {
			negLogLik.all[[m, K]] <- calcLik(params.temp = paramsPred.all[[m, K]], dat = dat, holdOut = K)
		}
	}}

#------------------------------------------------------------------------------
# Table of likelihood for each model 
negLogLikTable <- data.frame(
	model = rep(1:dim(models2test)[1], each = 6), 
	u0 = rep(models2test[, 'u0'], each = 6), 
	u1 = rep(models2test[, 'u1'], each = 6), 
	rho = rep(models2test[, 'rho'], each = 6), 
	validationReplicate = rep(c(1:5, "none"), dim(models2test)[1]), 
	training_nLL = rep(NA, dim(models2test)[1]*6),
	validation_nLL = rep(NA, dim(models2test)[1]*6),
	overall_nLL = rep(NA, dim(models2test)[1]*6)) 

for(m in 1:dim(models2test)[1]){
	for(K in 1:5){
		negLogLikTable[negLogLikTable$model == m & negLogLikTable$validationReplicate == K, 'training_nLL'] <- negLogLik.all[[m, K]][[1]]['training']
		negLogLikTable[negLogLikTable$model == m & negLogLikTable$validationReplicate == K, 'validation_nLL'] <- negLogLik.all[[m, K]][[1]]['validation']
		negLogLikTable[negLogLikTable$model == m & negLogLikTable$validationReplicate == K, 'overall_nLL'] <- NA
	}
	
	negLogLikTable[negLogLikTable$model == m & negLogLikTable$validationReplicate == "none", 'overall_nLL'] <- negLogLik.all[[m, 6]][[1]]['overall']
}

negLogLikTable$mKcode <- paste(negLogLikTable$model, negLogLikTable$validationReplicate, sep="-")

sumNLL <- data.frame(
	models2test, 
	training = round(tapply(negLogLikTable$training_nLL, negLogLikTable$model, mean, na.rm = TRUE), 1),
	validation = round(tapply(negLogLikTable$validation_nLL, negLogLikTable$model, mean, na.rm = TRUE), 1), 
	overall = round(tapply(negLogLikTable$overall_nLL, negLogLikTable$model, mean, na.rm = TRUE), 1))

sumNLL <- sumNLL[order(sumNLL$validation, decreasing = FALSE), ]

sumNLL$AICoverall <- round(calcAIC(nLL = sumNLL$overall, nParams = as.numeric(models2test[match(rownames(sumNLL), rownames(models2test)), 'nParams'])))


###############################################################################
# Tables with model comparison
###############################################################################
write.csv(sumNLL, file = "dataAnalysis/output/modelComparison_Isigma.csv")

###############################################################################
# Tables with parameter output
###############################################################################
if(Isigma == TRUE){
	tab <- data.frame(
		temp = c(dat$T.obs, "y[0]", "E"),
		mu0 = c(paste(
			sprintf("%.3f", exp(summary.in[[1, 6]][[1]][c(1:7), 1])), " (", 
			sprintf("%.3f", exp(summary.in[[1, 6]][[2]][c(1:7), 1])), ", ",
			sprintf("%.3f", exp(summary.in[[1, 6]][[2]][c(1:7), 5])), ")", sep=""),
			paste(
				sprintf("%.3f", summary.in[[bestMod, 6]][[1]][c('a[1]', "E[1]"), 1]), " (", 
				sprintf("%.3f", summary.in[[bestMod, 6]][[2]][c('a[1]', "E[1]"), 1]), ", ", 
				sprintf("%.3f", summary.in[[bestMod, 6]][[2]][c('a[1]', "E[1]"), 5]), ")", sep="")),
		mu1 = c(paste(
			sprintf("%.3f", exp(summary.in[[1, 6]][[1]][c(8:14), 1])), " (", 
			sprintf("%.3f", exp(summary.in[[1, 6]][[2]][c(8:14), 1])), ", ",
			sprintf("%.3f", exp(summary.in[[1, 6]][[2]][c(8:14), 5])), ")", sep=""),
			paste(
				sprintf("%.3f", summary.in[[bestMod, 6]][[1]]['a[2]', 1]), " (", 
				sprintf("%.3f", summary.in[[bestMod, 6]][[2]]['a[2]', 1]), ", ", 
				sprintf("%.3f", summary.in[[bestMod, 6]][[2]]['a[2]', 5]), ")", sep=""), NA),
		rho = c(paste(
			sprintf("%.3f", exp(summary.in[[1, 6]][[1]][15:21, 1])), " (", 
			sprintf("%.3f", exp(summary.in[[1, 6]][[2]][15:21, 1])), ", ",
			sprintf("%.3f", exp(summary.in[[1, 6]][[2]][15:21, 5])), ")", sep=""),
			paste(
				sprintf("%.3f", summary.in[[bestMod, 6]][[1]][c('a[3]', 'E[3]'), 1]), " (", 
				sprintf("%.3f", summary.in[[bestMod, 6]][[2]][c('a[3]', 'E[3]'), 1]), ", ", 
				sprintf("%.3f", summary.in[[bestMod, 6]][[2]][c('a[3]', 'E[3]'), 5]), ")", sep="")),
		sigma1 = c(paste(
			sprintf("%.3f", exp(summary.in[[1, 6]][[1]][22:28, 1])), " (", 
			sprintf("%.3f", exp(summary.in[[1, 6]][[2]][22:28, 1])), ", ",
			sprintf("%.3f", exp(summary.in[[1, 6]][[2]][22:28, 5])), ")", sep=""), rep(NA, 2)),
		sigmaMTE = c(paste(
			sprintf("%.3f", exp(summary.in[[bestMod, 6]][[1]][6:12, 1])), " (", 
			sprintf("%.3f", exp(summary.in[[bestMod, 6]][[2]][6:12, 1])), ", ",
			sprintf("%.3f", exp(summary.in[[bestMod, 6]][[2]][6:12, 5])), ")", sep=""), rep(NA, 2))
	)
	
	
}else{

	tab <- data.frame(
	temp = c(dat$T.obs, "sigma", "y[0]", "E", "Eh", "T[h]", "sigma"),
	mu0 = c(paste(
		sprintf("%.3f", exp(summary.in[[1, 6]][[1]][c(1:7,22), 1])), " (", 
		sprintf("%.3f", exp(summary.in[[1, 6]][[2]][c(1:7,22), 1])), ", ",
	  sprintf("%.3f", exp(summary.in[[1, 6]][[2]][c(1:7,22), 5])), ")", sep=""),
		paste(
			sprintf("%.3f", summary.in[[bestMod, 6]][[1]][c('a[1]', "E[1]"), 1]), " (", 
			sprintf("%.3f", summary.in[[bestMod, 6]][[2]][c('a[1]', "E[1]"), 1]), ", ", 
			sprintf("%.3f", summary.in[[bestMod, 6]][[2]][c('a[1]', "E[1]"), 5]), ")", sep=""),
		c(NA, NA),
		paste(
			sprintf("%.3f", exp(summary.in[[bestMod, 6]][[1]]['sigma.log', 1])), " (", 
			sprintf("%.3f", exp(summary.in[[bestMod, 6]][[2]]['sigma.log', 1])), ", ", 
			sprintf("%.3f", exp(summary.in[[bestMod, 6]][[2]]['sigma.log', 5])), ")", sep="")),
	mu1 = c(paste(
		sprintf("%.3f", exp(summary.in[[1, 6]][[1]][c(8:14), 1])), " (", 
		sprintf("%.3f", exp(summary.in[[1, 6]][[2]][c(8:14), 1])), ", ",
		sprintf("%.3f", exp(summary.in[[1, 6]][[2]][c(8:14), 5])), ")", sep=""),
	  NA,
		paste(
			sprintf("%.3f", summary.in[[bestMod, 6]][[1]]['a[2]', 1]), " (", 
			sprintf("%.3f", summary.in[[bestMod, 6]][[2]]['a[2]', 1]), ", ", 
			sprintf("%.3f", summary.in[[bestMod, 6]][[2]]['a[2]', 5]), ")", sep=""),
		rep(NA,4)),
	rho = c(paste(
		sprintf("%.3f", exp(summary.in[[1, 6]][[1]][15:21, 1])), " (", 
		sprintf("%.3f", exp(summary.in[[1, 6]][[2]][15:21, 1])), ", ",
		sprintf("%.3f", exp(summary.in[[1, 6]][[2]][15:21, 5])), ")", sep=""),
		NA, 
		paste(
			sprintf("%.3f", summary.in[[bestMod, 6]][[1]][c('a[3]', 'E[3]', 'Eh[3]', 'Th[3]'), 1]), " (", 
			sprintf("%.3f", summary.in[[bestMod, 6]][[2]][c('a[3]', 'E[3]', 'Eh[3]', 'Th[3]'), 1]), ", ", 
			sprintf("%.3f", summary.in[[bestMod, 6]][[2]][c('a[3]', 'E[3]', 'Eh[3]', 'Th[3]'), 5]), ")", sep=""),
		NA)
)
}
write.csv(tab, file = "dataAnalysis/output/bestModelParameters_Isigma.csv")

###############################################################################
## Plot I estimates and MTE overall
###############################################################################

K <- 6

parI <- list(
	mean = paramsPred.all[[1, K]],
	lower  = exp(matrix(summary.in[[1, K]][[2]][paste('p.mean.log[', rep(c(1:dat$nT), each = 3),",", rep(1:3, dat$nT), "]", sep=""), 1], byrow = TRUE, nrow = dat$nT)),
	upper = exp(matrix(summary.in[[1, K]][[2]][paste('p.mean.log[', rep(c(1:dat$nT), each = 3),",", rep(1:3, dat$nT), "]", sep=""), 5], byrow = TRUE, nrow = dat$nT))
)

bestMod <- 5
T.all <- seq(0, 40, 0.1)


parMod <- list(
	mean = cbind(
		u0 = rep(NA, length(T.all)),
		u1 = rep(NA, length(T.all)),
		rho = rep(NA, length(T.all))),
	lower = cbind(
		u0 = rep(NA, length(T.all)),
		u1 = rep(NA, length(T.all)),
		rho = rep(NA, length(T.all))),
	upper = cbind(
		u0 = rep(NA, length(T.all)),
		u1 = rep(NA, length(T.all)),
		rho = rep(NA, length(T.all)))
)

for(j in 1:3){ # for each parameter (u0, u1, rho)
	
	#-----------------------------------------------------------------------------
	# Constant
	if(models2test[bestMod, j] == "C"){
		parMod$mean[, j] <- rep(summary.in[[bestMod, 6]][[1]][paste('a[', j, ']', sep=""), 1], length(T.all))
		Cdummy <- numeric(1000)
		chains <- sample(c(1:length(in.all[[bestMod, 6]])), 1000, replace = TRUE)
		iter <- sample(1:numIter, 1000, replace = TRUE)
		for(i in 1:1000){
			Cdummy[i] <- in.all[[bestMod, 6]][[chains[i]]][iter[i], paste('a[', j, ']', sep="")]
		}
		parMod$lower[, j] <- rep(quantile(Cdummy, 0.025), length(T.all))
		parMod$upper[, j] <- rep(quantile(Cdummy, 0.975), length(T.all))
		
		#-----------------------------------------------------------------------------
		# Arrhenius
		
	}else if(models2test[bestMod, j] == "A"){
		
		params.j <- summary.in[[bestMod, 6]]
		parMod$mean[, j] <- params.j[[1]][paste("a[", j, "]", sep = ""), 1] * exp(-params.j[[1]][paste("E[", j, "]", sep = ""), 1]/(8.62*10^-5) * (1/(T.all + 273.15) - 1/(15 + 273.15)))
		
		Adummy <- matrix(NA, nrow = length(T.all), ncol = 1000)
		chains <- sample(c(1:length(in.all[[bestMod, 6]])), 1000, replace = TRUE)
		iter <- sample(1:numIter, 1000, replace = TRUE)
		for(i in 1:1000){
			Adummy[, i] <- in.all[[bestMod, 6]][[chains[i]]][iter[i], paste("a[", j, "]", sep = "")] * exp(-in.all[[bestMod, 6]][[chains[i]]][iter[i], paste("E[", j, "]", sep = "")]/(8.62*10^-5) * (1/(T.all + 273.15) - 1/(15 + 273.15)))
		}
		parMod$lower[, j] <- apply(Adummy, 1, quantile, 0.025)
		parMod$upper[, j] <- apply(Adummy, 1, quantile, 0.975)
		
		#-----------------------------------------------------------------------------
		# SSU
	}else if(models2test[bestMod, j] == "SSU"){
		
		params.j <- summary.in[[bestMod, 6]]
		parMod$mean[, j] <- params.j[[1]][paste("a[", j, "]", sep = ""), 1] * exp(-params.j[[1]][paste("E[", j, "]", sep = ""), 1]/(8.62*10^-5) * (1/(T.all + 273.15) - 1/(15 + 273.15))) * (1 + exp(params.j[[1]][paste("Eh[", j, "]", sep = ""), 1] / (8.62*10^-5)*(-1/(T.all+273.15)+1/(params.j[[1]][paste("Th[", j, "]", sep = ""), 1]+273.15))))^c(1, 1, -1)[j]
		
		SSUdummy <- matrix(NA, nrow = length(T.all), ncol = 1000)
		chains <- sample(c(1:length(in.all[[bestMod,6]])), 1000, replace = TRUE) 
		iter <- sample(1:numIter, 1000, replace = TRUE)
		for(i in 1:1000){
			SSUdummy[, i] <- in.all[[bestMod,6]][[chains[i]]][iter[i], paste("a[", j, "]", sep = "")] * exp(-in.all[[bestMod,6]][[chains[i]]][iter[i], paste("E[", j, "]", sep = "")]/(8.62*10^-5) * (1/(T.all + 273.15) - 1/(15 + 273.15))) * (1 + exp(in.all[[bestMod,6]][[chains[i]]][iter[i], paste("Eh[", j, "]", sep = "")] / (8.62*10^-5)*(-1/(T.all+273.15)+1/(in.all[[bestMod,6]][[chains[i]]][iter[i], paste("Th[", j, "]", sep = "")]+273.15))))^c(1, 1, -1)[j]
		}
		
		parMod$lower[, j] <- apply(SSUdummy, 1, quantile, 0.025)
		parMod$upper[, j] <- apply(SSUdummy, 1, quantile, 0.975)
		
	}
} # end j


yMax <- c(0.6, 0.1, 0.1)
quartz(width = 6.3, height = 2.3, pointsize = 10)
par(mfrow = c(1, 3), mar = c(4, 3, 2, 1), oma = c(1, 3, 0, 0))

for(j in 1:3){
	
	# Temperature independent estimates
	plotCI(dat$T.obs, parI[[1]][, j], li = parI[[2]][, j], ui = parI[[3]][, j], gap = 0, ylim = c(0, yMax[j]), xlab = "", bty = "l", ylab = "", las = 1, xlim = c(0, 40), pch = 21, pt.bg = "white")
	mtext(side = 3, adj =0, line = 0.5, c(expression(paste("a) ", mu[0])), expression(paste("b) ", mu[1])), expression(paste("c) ", rho)))[j])
	
	polygon(x = c(T.all, rev(T.all)), y = c(parMod$lower[, j], rev(parMod$upper[, j])), border = NA, col = c("#FF000060", "#00000060", "#00FF0060")[j])
	lines(T.all, parMod$mean[, j], col= c(2, 1, 3)[j])
}

mtext(side = 1, line = 0, expression(paste("Temperature (", degree, "C)", sep = "")), outer = TRUE)
mtext(side = 2, line = 0, expression(paste("Rate (", italic(y(T)), ")", sep = "")), outer = TRUE)

###############################################################################
## Plot against data
###############################################################################

t.all <- c(1:120)

larvPred <- array(NA, dim = c(dat$nT, length(t.all), 2, 2))
# dim = temp (7), time (1200), stage (2), model (2)

# Simulate model predictions for I and bestMod
for(i in 1:dat$nT){
	larvPred[i, , , 1] <- sim.larv(
			times = t.all, 
			params = list(
				N0 = 50, 
				u0 = paramsPred.all[[1, 6]][i,1], 
				u1 = paramsPred.all[[1, 6]][i,2], 
				rho = paramsPred.all[[1, 6]][i,3], 
				sigma = paramsPred.all[[1, 6]][i,4]))

	larvPred[i, , , 2] <- sim.larv(
		times = t.all, 
		params = list(
			N0 = 50, 
			u0 = paramsPred.all[[bestMod, 6]][i,1], 
			u1 = paramsPred.all[[bestMod, 6]][i,2], 
			rho = paramsPred.all[[bestMod, 6]][i,3], 
			sigma = paramsPred.all[[bestMod, 6]][i,4]))
}

# Add confidence intervals??
larvPredCI <- list(
	lower = array(NA, dim = c(dat$nT, length(t.all), 2, 2)),
	upper = array(NA, dim = c(dat$nT, length(t.all), 2, 2)))

for(i in 1:dat$nT){
	
	# Independent estimates
	larvDummy <- array(NA, dim = c(1000, length(t.all), 2))
	chains <- sample(c(1:length(in.all[[bestMod,6]])), 1000, replace = TRUE) 
	iter <- sample(1:numIter, 1000, replace = TRUE)
	for(j in 1:1000){
  	params.j <- exp(in.all[[1, 6]][[chains[i]]][iter[i], c(paste("p.mean.log[", i, ",", c(1:3), "]", sep=""), "sigma.log")])
  	larvDummy[j, , ] <- sim.larv(
  		times = t.all, 
  		params = list(
  			N0 = 50, 
  			u0 = params.j[1], 
  			u1 = params.j[2], 
  			rho = params.j[3], 
  			sigma = params.j[4]))
	} # end j
	
	CIdummy_preInf <- apply(larvDummy[,,1], 2, quantile, c(0.025, 0.975))
	CIdummy_inf <- apply(larvDummy[,,2], 2, quantile, c(0.025, 0.975))
	larvPredCI$lower[i, , 1, 1] <- CIdummy_preInf['2.5%', ]
	larvPredCI$upper[i, , 1, 1] <- CIdummy_preInf['97.5%', ]
	larvPredCI$lower[i, , 2, 1] <- CIdummy_inf['2.5%', ]
	larvPredCI$upper[i, , 2, 1] <- CIdummy_inf['97.5%', ]

	# MTE estimates
	larvDummy <- array(NA, dim = c(1000, length(t.all), 2))
	chains <- sample(c(1:length(in.all[[bestMod,6]])), 1000, replace = TRUE) 
	iter <- sample(1:numIter, 1000, replace = TRUE)
	for(j in 1:1000){
		params.j <- in.all[[bestMod, 6]][[chains[i]]][iter[i], ]
		larvDummy[j, , ] <- sim.larv(
			times = t.all, 
			params = list(
				N0 = 50, 
				u0 = as.numeric(params.j['a[1]'] * exp(-params.j['E[1]']/(8.62*10^-5) * (1/(dat$T.obs[i] + 273.15) - 1/(15 + 273.15)))), 
				u1 = as.numeric(params.j['a[2]']), 
				rho = as.numeric(params.j['a[3]'] * exp(-params.j['E[3]']/(8.62*10^-5) * (1/(dat$T.obs[i] + 273.15) - 1/(15 + 273.15))) * (1 + exp(params.j['Eh[3]'] / (8.62*10^-5)*(-1/(dat$T.obs[i]+273.15)+1/(params.j['Th[3]']+273.15))))^(-1)), 
				sigma = exp(params.j['sigma.log'])))
	} # end j
	
	CIdummy_preInf <- apply(larvDummy[,,1], 2, quantile, c(0.025, 0.975))
	CIdummy_inf <- apply(larvDummy[,,2], 2, quantile, c(0.025, 0.975))
	larvPredCI$lower[i, , 1, 2] <- CIdummy_preInf['2.5%', ]
	larvPredCI$upper[i, , 1, 2] <- CIdummy_preInf['97.5%', ]
	larvPredCI$lower[i, , 2, 2] <- CIdummy_inf['2.5%', ]
	larvPredCI$upper[i, , 2, 2] <- CIdummy_inf['97.5%', ]
  
	} # end temperature i

#----------------------------------------------------------------------------
# Plot data

T2plot <- c(1,2,3,6)
datCol <- grey(0.7) # "#00000050"
quartz(width = 7, height = 3, pointsize = 10)
par(mfcol = c(2,4), mar = c(1, 1, 0, 0), oma = c(4,4,2,1))

for(i in 1:length(T2plot)){
	
	#-----------------------------------------------------------------------------
	# Pre-infective
	plot(1:120, 1:120, "n", ylim = c(0, 50), bty = "l", xlab = "", ylab = "", las = 1, yaxt="n", xaxt="n")
	axis(side = 1, labels = FALSE)
	if(i == 1) axis(side = 2, las = 1) else axis(side = 2, labels = FALSE)
	
	# polygon(x = c(t.all, rev(t.all)), y = c(larvPredCI$lower[T2plot[i], , 1, 1], rev(larvPredCI$upper[T2plot[i], , 1, 1])), col = "#0000FF30", border = NA)
	# polygon(x = c(t.all, rev(t.all)), y = c(larvPredCI$lower[T2plot[i], , 1, 2], rev(larvPredCI$upper[T2plot[i], , 1, 2])), col = "#FF000030", border = NA)
	
	for(j in 1:5) points(dat$t[dat$ind[T2plot[i], ,j]], dat$n_obs[T2plot[i], , j, 1, 1], pch = j, col = datCol, xpd = NA)
	
	lines(t.all, larvPred[T2plot[i], ,1,1], col = 4, lwd = 1.5)
	lines(t.all, larvPred[T2plot[i], ,1,2], col = 2, lwd = 1.5)
	mtext(side = 3, substitute(paste(temp, degree, "C"), list(temp = dat$T.obs[T2plot[i]])))
	if(i == 1) mtext(side = 2, "Pre-infective larvae", line = 3, cex =0.8)
	if(i == length(T2plot)) legend("topright", col = c(datCol, 4,2), pch = c(1, NA, NA), lwd = c(1, 1.5, 1.5), lty = c(NA, 1, 1), legend = c("data", "indep. est.", "MTE est."), bty = "n")
	
	#-----------------------------------------------------------------------------
	# Infective
	plot(1:120, 1:120, "n", ylim = c(0, 20), bty = "l", xlab = "", ylab = "", las = 1, yaxt="n")
	if(i == 1) axis(side = 2, las = 1) else axis(side = 2, labels = FALSE)
	
	# polygon(x = c(t.all, rev(t.all)), y = c(larvPredCI$lower[T2plot[i], , 2, 1], rev(larvPredCI$upper[T2plot[i], , 2, 1])), col = "#0000FF30", border = NA)
	# polygon(x = c(t.all, rev(t.all)), y = c(larvPredCI$lower[T2plot[i], , 2, 2], rev(larvPredCI$upper[T2plot[i], , 2, 2])), col = "#FF000030", border = NA)
	
	for(j in 1:5) points(dat$t[dat$ind[T2plot[i], ,j]], dat$n_obs[T2plot[i], , j, 1, 2], pch = j, col = datCol, xpd = NA)
	
	lines(t.all, larvPred[T2plot[i], ,2,1], col = 4, lwd = 1.5)
	lines(t.all, larvPred[T2plot[i], ,2,2], col = 2, lwd = 1.5)
	if(i == 1) mtext(side = 2, "Infective larvae", line = 3, cex = 0.8)
}

mtext(side = 1, line = 2,outer = TRUE, "Time (days)")

