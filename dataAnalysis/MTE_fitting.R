###############################################################################
# This code fits temperature relationships from the metabolic theory of 
# ecology to development and mortality parameters for free-living Ostertagia
# larvae. The data are from temperature manipulation experiments conducted
# in the lab.
###############################################################################

library(here)
library(dclone)
library(rjags)
library(scales) #for plotting
library(gplots)

# 
source("dataAnalysis/functions.R")

source("dataAnalysis/MTE_models.R")

# Within models, parameter order is always u0, u1, rho
models2test <- rbind(
	model1 = c("I", "I", "I"), 
	model2 = c("A", "A", "A"), 
	model3 = c("A", "A", "SSU"),
	model4 = c("A", "C", "A"),
	model5 = c("A", "C", "SSU"),
	model6 = c("SSU", "A", "A"),
	model7 = c("SSU", "A", "SSU"),
	model8 = c("SSU", "C", "A"),
	model9 = c("SSU", "C", "SSU"))
dimnames(models2test)[[2]] <- c("u0", "u1", "rho")


# Load data
dat<-load.data(excl.temp=40)

###############################################################################
## Fitting
###############################################################################
numChains <- 4 # Number of cores to use = number of chains
numIter <- 10000 # Number of MCMC steps
	
# Store results
procTime <- array(NA, dim = c(dim(models2test)[1], 6), dimnames = list(rownames(models2test), c(paste("rep", 1:5, sep=""), "none")))

out.all <- list(); length(out.all) <- dim(models2test)[1]*6; dim(out.all) <- c(dim(models2test)[1], 6)

#------------------------------------------------------------------------------
for(m in 1:dim(models2test)[1]){ # for each model being tested
	
	# Pick model 
	model <- pickModel(m)
	
	for(K in 1:6){ # for each hold-one-out validation run
		
		# Which reps to include?
		if(K == 6) incl.reps <- c(1:5) else incl.reps <- c(1:5)[which(c(1:5) != K)]
		
		# Create dat list excluding one rep
		dat.mK <- list(
			nT = dat$nT,
			T.obs = dat$T.obs,
			dt = dat$dt, 
			t = dat$t, 	# Times that counts were performed
			nt = dat$nt, 	# Number of times that counts were
			nt_obs = dat$nt_obs[, incl.reps],
			ind = dat$ind[, , incl.reps],
			n_obs = dat$n_obs[, , incl.reps, , ],
			nk = 4,
			sigN0 = 0.3,
			sigp = 0.05
		)
		
		if(m == 1){ # for the I model (independent temperature fits) the parameters are different
			parNames <- c("p.mean.log", "sigma.log", "N0", "p.rep")
		} else {
			parNames <- c("a", "E", "Eh", "Th", "sigma.log", "N0", "p.rep")
		}
		#----------------------------------------------------------------------------
		# Fit model to subsetted data
		cl <- makeCluster(numChains)
		
		t.start<-proc.time()[3]
		
		fit.mK <- jags.parfit(
			cl, 
			data = dat.mK, 
			params = parNames, 
			model = model, 
			n.iter = numIter, 
			n.chains = numChains) 
		
		stopCluster(cl)
		procTime[m, K] <- round((proc.time()[3] - t.start)/(60), 1) # track process time
		cat(paste("Process time for model", m, "and hold-out", K, "=", procTime[m, K], "minutes"))
		
		#----------------------------------------------------------------------------
		# Store output
		out.all[[m, K]] <- fit.mK
		
	}}
	
saveRDS(out.all, file = "MTE_output_kFold.rds")
# readRDS("MTE_output_kFold.rds")


###############################################################################
## Likelihood
###############################################################################
negLogLik.all <- list(); length(negLogLik.all) <- dim(models2test)[1]*6; dim(negLogLik.all) <- c(dim(models2test)[1], 6)
summary.out <- list(); length(summary.out) <- dim(models2test)[1]*6; dim(summary.out) <- c(dim(models2test)[1], 6)
paramsPred.all <- list(); length(paramsPred.all) <- dim(models2test)[1]*6; dim(paramsPred.all) <- c(dim(models2test)[1], 6)

for(m in 1:dim(models2test)[1]){ 
	for(K in 1:5){
		
		# Summary of model output, including parameters
		out <- summary(out.all[[m, K]])
		summary.out[[m, K]] <- out
		
		# matrix of parameter values for each temperature (rows) and parameter (columns)
		paramsPred <- matrix(NA, nrow = dat$nT, ncol = 4)
		for(j in 1:3){ # for each parameter u0, u1, rho
			if(models2test[m, j] == "I"){
				paramsPred[, j] <- as.numeric(exp(out[[1]][paste('p.mean.log[', c(1:dat$nT),",", j, "]", sep=""), 1]))
			} else if(models2test[m, j] == "C"){
				params.j <- c(a = out[[1]][paste("a[", j, "]", sep=""), 1])
				paramsPred[, j] <- MTE.predict(model = "C", params = params.j, T.all = dat$T.obs,	exponent = c(1, 1, -1, NA)[j])
			
			} else if (models2test[m, j] == "A"){
				params.j <- c(
					a = out[[1]][paste("a[", j, "]", sep=""), 1], 
					E = out[[1]][paste("E[", j, "]", sep=""), 1])
				paramsPred[, j] <- MTE.predict(model = "A", params = params.j, T.all = dat$T.obs,	exponent = c(1, 1, -1, NA)[j])
			
			} else if (models2test[m, j] == "SSU"){
				params.j <- c(
					a = out[[1]][paste("a[", j, "]", sep=""), 1], 
					E = out[[1]][paste("E[", j, "]", sep=""), 1], 
					Eh = out[[1]][paste("Eh[", j, "]", sep=""), 1], 
					Th = out[[1]][paste("Th[", j, "]", sep=""), 1])
				paramsPred[, j] <- MTE.predict(model = "SSU", params = params.j, T.all = dat$T.obs,	exponent = c(1, 1, -1, NA)[j])
			}
		} # end j
		
		# Constant sigma
		paramsPred[, 4] <- rep(exp(out[[1]]['sigma.log', 1]), dat$nT)
		
		paramsPred.all[[m, K]] <- paramsPred
		# Calculate liklihood
		# negLogLik.all[[m, K]] <- calcLik(params.temp = paramsPred, dat = dat, holdOut = K)
	}
}

# Table of likelihood for each model 
negLogLikTable <- data.frame(
	model = rep(1:dim(models2test)[1], each = 5), 
	u0 = rep(models2test[, 'u0'], each = 5), 
	u1 = rep(models2test[, 'u1'], each = 5), 
	rho = rep(models2test[, 'rho'], each = 5), 
	sigma = rep(models2test[, 'sigma'], each = 5),
	validationReplicate = rep(c(1:5), dim(models2test)[1]), 
	training_nLL = rep(NA, dim(models2test)[1]*5),
	validation_nLL = rep(NA, dim(models2test)[1]*5),
	overall_nLL = rep(NA, dim(models2test)[1]*5)) 

for(m in 1:dim(models2test)[1]){
	for(K in 1:5){
		negLogLikTable[negLogLikTable$model == m & negLogLikTable$validationReplicate == K, 'training_nLL'] <- negLogLik.all[[m, K]][[1]]['training']
		negLogLikTable[negLogLikTable$model == m & negLogLikTable$validationReplicate == K, 'validation_nLL'] <- negLogLik.all[[m, K]][[1]]['validation']
		negLogLikTable[negLogLikTable$model == m & negLogLikTable$validationReplicate == K, 'overall_nLL'] <- negLogLik.all[[m, K]][[1]]['overall']
	}
}

negLogLikTable$mKcode <- paste(negLogLikTable$model, negLogLikTable$validationReplicate, sep="-")

sumNLL <- data.frame(models2test, training = tapply(negLogLikTable$training_nLL, negLogLikTable$model, mean),validation = tapply(negLogLikTable$validation_nLL, negLogLikTable$model, mean), overall = tapply(negLogLikTable$overall_nLL, negLogLikTable$model, mean))

order(sumNLL$overall, decreasing = FALSE)

sumNLL <- sumNLL[order(sumNLL$validation, decreasing = FALSE), ]

###############################################################################
## Plot I estimates and MTE for each holdout
###############################################################################

K <- 6

parI <- list(
	mean = paramsPred.all[[1, K]],
	lower  = exp(matrix(summary.out[[1, K]][[2]][paste('p.mean.log[', rep(c(1:dat$nT), each = 3),",", rep(1:3, dat$nT), "]", sep=""), 1], byrow = TRUE, nrow = dat$nT)),
	upper = exp(matrix(summary.out[[1, K]][[2]][paste('p.mean.log[', rep(c(1:dat$nT), each = 3),",", rep(1:3, dat$nT), "]", sep=""), 5], byrow = TRUE, nrow = dat$nT))
)

bestMod <- rbind(A = c(7, 3, 2), SSU = c(5, 4, 7))
T.all <- seq(0, 40, 0.1)


parA <- list(
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
parSSU <- parA
parC <- parA

for(j in 1:3){ # for each parameter (u0, u1, rho)
	
	#-----------------------------------------------------------------------------
	# Constant
	if(j == 2){
		parC$mean[, j] <- rep(summary.out[[7, K]][[1]]['a[2]', 1], length(T.all))
		Cdummy <- numeric(1000)
		chains <- sample(c(1:4), 1000, replace = TRUE)
		iter <- sample(1:numIter, 1000, replace = TRUE)
		for(i in 1:1000){
			Cdummy[i] <- out.all[[7, K]][[chains[i]]][iter[i], 'a[2]']
		}
		parC$lower[, j] <- rep(quantile(Cdummy, 0.025), length(T.all))
		parC$upper[, j] <- rep(quantile(Cdummy, 0.975), length(T.all))
	}

	
	#-----------------------------------------------------------------------------
	# Arrhenius
	
	params.j <- summary.out[[bestMod['A', j], K]]
	parA$mean[, j] <- params.j[[1]][paste("a[", j, "]", sep = ""), 1] * exp(-params.j[[1]][paste("E[", j, "]", sep = ""), 1]/(8.62*10^-5) * (1/(T.all + 273.15) - 1/(15 + 273.15)))
	
	Adummy <- matrix(NA, nrow = length(T.all), ncol = 1000)
	chains <- sample(c(1:4), 1000, replace = TRUE)
	iter <- sample(1:numIter, 1000, replace = TRUE)
	for(i in 1:1000){
		Adummy[, i] <- out.all[[bestMod['A', j], K]][[chains[i]]][iter[i], paste("a[", j, "]", sep = "")] * exp(-out.all[[bestMod['A', j], K]][[chains[i]]][iter[i], paste("E[", j, "]", sep = "")]/(8.62*10^-5) * (1/(T.all + 273.15) - 1/(15 + 273.15)))
	}
	parA$lower[, j] <- apply(Adummy, 1, quantile, 0.025)
	parA$upper[, j] <- apply(Adummy, 1, quantile, 0.975)
	
	#-----------------------------------------------------------------------------
	# SSU
	
	params.j <- summary.out[[bestMod['SSU', j], K]]
	parSSU$mean[, j] <- params.j[[1]][paste("a[", j, "]", sep = ""), 1] * exp(-params.j[[1]][paste("E[", j, "]", sep = ""), 1]/(8.62*10^-5) * (1/(T.all + 273.15) - 1/(15 + 273.15))) * (1 + exp(params.j[[1]][paste("Eh[", j, "]", sep = ""), 1] / (8.62*10^-5)*(-1/(T.all+273.15)+1/(params.j[[1]][paste("Th[", j, "]", sep = ""), 1]+273.15))))^c(1, 1, -1)[j]
	
	SSUdummy <- matrix(NA, nrow = length(T.all), ncol = 1000)
	if(j == 1) chains <- sample(c(1:3), 1000, replace = TRUE) else chains <- sample(c(1:4), 1000, replace = TRUE)
	iter <- sample(1:numIter, 1000, replace = TRUE)
	for(i in 1:1000){
		SSUdummy[, i] <- out.all[[bestMod['SSU', j], K]][[chains[i]]][iter[i], paste("a[", j, "]", sep = "")] * exp(-out.all[[bestMod['SSU', j], K]][[chains[i]]][iter[i], paste("E[", j, "]", sep = "")]/(8.62*10^-5) * (1/(T.all + 273.15) - 1/(15 + 273.15))) * (1 + exp(out.all[[bestMod['SSU', j], K]][[chains[i]]][iter[i], paste("Eh[", j, "]", sep = "")] / (8.62*10^-5)*(-1/(T.all+273.15)+1/(out.all[[bestMod['SSU', j], K]][[chains[i]]][iter[i], paste("Th[", j, "]", sep = "")]+273.15))))^c(1, 1, -1)[j]
	}
	
	#plot(T.all, parSSU$mean[, j], "l", ylim = c(0, 0.6));for(i in sample(1:1000, size = 30)) lines(T.all, SSUdummy[, i], col="#00000030")
	# lines(T.all, SSUdummy[, i], col=2); abline(v = out.all[[bestMod['SSU', j], K]][[chains[i]]][iter[i], paste("Th[", j, "]", sep = "")], lty = 2, col = 2)
	parSSU$lower[, j] <- apply(SSUdummy, 1, quantile, 0.025)
	parSSU$upper[, j] <- apply(SSUdummy, 1, quantile, 0.975)
	
}


yMax <- c(0.6, 0.1, 0.1)
quartz(width = 6.3, height = 2.3, pointsize = 10)
par(mfrow = c(1, 3), mar = c(4, 3, 2, 1), oma = c(1, 3, 0, 0))

for(j in 1:3){
	
	# Temperature independent estimates
	plotCI(dat$T.obs, parI[[1]][, j], li = parI[[2]][, j], ui = parI[[3]][, j], gap = 0, ylim = c(0, yMax[j]), xlab = "", bty = "l", ylab = "", las = 1, xlim = c(0, 40), pch = 21, pt.bg = "white")
	mtext(side = 3, adj =0, line = 0.5, c(expression(paste("a) ", mu[0])), expression(paste("b) ", mu[1])), expression(paste("c) ", rho)))[j])
	
	# MTE estimated curves
	if(j == 2){
		polygon(x = c(T.all, rev(T.all)), y = c(parC$lower[, j], rev(parC$upper[, j])), border = NA, col = "#00000060")
		lines(T.all, parC$mean[, j], col = 1)
	}
	
	polygon(x = c(T.all, rev(T.all)), y = c(parA$lower[, j], rev(parA$upper[, j])), border = NA, col = "#FF000060")
	lines(T.all, parA$mean[, j], col = 2)
	
	polygon(x = c(T.all, rev(T.all)), y = c(parSSU$lower[, j], rev(parSSU$upper[, j])), border = NA, col = "#00FF0060")
	lines(T.all, parSSU$mean[, j], col = 3)
	
	
}
mtext(side = 1, line = 0, expression(paste("Temperature (", degree, "C)", sep = "")), outer = TRUE)
mtext(side = 2, line = 0, expression(paste("Rate (", italic(y(T)), ")", sep = "")), outer = TRUE)

###############################################################################
###############################################################################
j <- 1
quartz(width = 5, height = 5, pointsize = 10)
par(mfrow= c(2,2), mar= c(2,2,2,2))

plot(out.all[[bestMod['SSU', j], K]][[4]][seq(1, 10000, 10), paste("a[", j, "]", sep = "")], out.all[[bestMod['SSU', j], K]][[4]][seq(1, 10000, 10), paste("E[", j, "]", sep = "")], xlab = "a", ylab = "E")
points(params.j[[1]][paste("a[", j, "]", sep = ""), 1], params.j[[1]][paste("E[", j, "]", sep = ""), 1], col =2, cex = 2)

plot(out.all[[bestMod['SSU', j], K]][[4]][seq(1, 10000, 10), paste("Eh[", j, "]", sep = "")], out.all[[bestMod['SSU', j], K]][[4]][seq(1, 10000, 10), paste("E[", j, "]", sep = "")], xlab = "a", ylab = "E")
points(params.j[[1]][paste("Eh[", j, "]", sep = ""), 1], params.j[[1]][paste("E[", j, "]", sep = ""), 1], col =2, cex = 2)

plot(out.all[[bestMod['SSU', j], K]][[4]][seq(1, 10000, 10), paste("a[", j, "]", sep = "")], out.all[[bestMod['SSU', j], K]][[4]][seq(1, 10000, 10), paste("Th[", j, "]", sep = "")], xlab = "a", ylab = "E")
points(params.j[[1]][paste("a[", j, "]", sep = ""), 1], params.j[[1]][paste("Th[", j, "]", sep = ""), 1], col =2, cex = 2)

plot(out.all[[bestMod['SSU', j], K]][[2]][seq(1, 10000, 10), paste("Eh[", j, "]", sep = "")], out.all[[bestMod['SSU', j], K]][[2]][seq(1, 10000, 10), paste("Th[", j, "]", sep = "")], xlab = "a", ylab = "E")
points(params.j[[1]][paste("Eh[", j, "]", sep = ""), 1], params.j[[1]][paste("Th[", j, "]", sep = ""), 1], col =2, cex = 2)

plot(seq(1, 10000, 10), out.all[[bestMod['SSU', j], K]][[1]][seq(1, 10000, 10), paste("Eh[", j, "]", sep = "")], "l", lty = 3)
lines(seq(1, 10000, 10), out.all[[bestMod['SSU', j], K]][[2]][seq(1, 10000, 10), paste("Eh[", j, "]", sep = "")], "l", lty = 3, col = 2)
lines(seq(1, 10000, 10), out.all[[bestMod['SSU', j], K]][[3]][seq(1, 10000, 10), paste("Eh[", j, "]", sep = "")], "l", lty = 3, col = 3)
lines(seq(1, 10000, 10), out.all[[bestMod['SSU', j], K]][[4]][seq(1, 10000, 10), paste("Eh[", j, "]", sep = "")], "l", lty = 3, col = 4)

plot(seq(1, 10000, 10), out.all[[bestMod['SSU', j], K]][[1]][seq(1, 10000, 10), paste("Th[", j, "]", sep = "")], "l", lty = 3)
lines(seq(1, 10000, 10), out.all[[bestMod['SSU', j], K]][[2]][seq(1, 10000, 10), paste("Th[", j, "]", sep = "")], "l", lty = 3, col = 2)
lines(seq(1, 10000, 10), out.all[[bestMod['SSU', j], K]][[3]][seq(1, 10000, 10), paste("Th[", j, "]", sep = "")], "l", lty = 3, col = 3)
lines(seq(1, 10000, 10), out.all[[bestMod['SSU', j], K]][[4]][seq(1, 10000, 10), paste("Th[", j, "]", sep = "")], "l", lty = 3, col = 4)

hist(out.all[[bestMod['SSU', j], K]][[4]][, paste("a[", j, "]", sep = "")])
abline(v = params.j[[1]][paste("a[", j, "]", sep = ""), 1], col = 2)

demean <- function(x){ return((x - mean(x))/sd(x))}
xmax <- 1000
plot(c(1:xmax), demean(out.all[[bestMod['SSU', j], K]][[4]][1:xmax, paste("a[", j, "]", sep = "")]), "l")
lines(c(1:xmax), demean(out.all[[bestMod['SSU', j], K]][[4]][1:xmax, paste("E[", j, "]", sep = "")]), col = 2)
lines(c(1:xmax), demean(out.all[[bestMod['SSU', j], K]][[4]][1:xmax, paste("Eh[", j, "]", sep = "")]), col = 3)
lines(c(1:xmax), demean(out.all[[bestMod['SSU', j], K]][[4]][1:xmax, paste("Th[", j, "]", sep = "")]), col = 4)
