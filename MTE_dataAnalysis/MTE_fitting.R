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

# Choose: Fit the models with contant sigma across temperatures
#         or different (Independent) sigma for each temperature
# source("dataAnalysis/MTE_models.R")
# source("dataAnalysis/MTE_models_Isigma.R")

# Within models, parameter order is always u0, u1, rho
models2test <- rbind(
	model1 = c("I", "I", "I", 22), 
	model2 = c("A", "A", "A", 7), 
	model3 = c("A", "A", "SSU", 9),
	model4 = c("A", "C", "A", 6),
	model5 = c("A", "C", "SSU", 8),
	model6 = c("SSU", "A", "A", 9),
	model7 = c("SSU", "A", "SSU", 11),
	model8 = c("SSU", "C", "A", 8),
	model9 = c("SSU", "C", "SSU", 10))
dimnames(models2test)[[2]] <- c("u0", "u1", "rho", 'nParams')

nF <- 6 # number of "folds" (k-fold) = each of five reps plus a fit with all data

# Load data
dat<-load.data(excl.temp=40)

###############################################################################
## Fitting
###############################################################################
numChains <- 6 # Number of cores to use = number of chains
# Increase number of chains so that we can pick and choose convergence.
numIter <- 10000 # Number of MCMC steps
	
# Store results
procTime <- array(NA, dim = c(dim(models2test)[1], nF), dimnames = list(rownames(models2test), c(paste("rep", 1:5, sep=""), "none")))

out.all <- list(); length(out.all) <- dim(models2test)[1]*6; dim(out.all) <- c(dim(models2test)[1], nF)

#------------------------------------------------------------------------------
for(m in 7:dim(models2test)[1]){ # for each model being tested
	
	# Pick model 
	model <- pickModel(m)
	
	for(K in c(1:6)){ # for each hold-one-out validation run
		
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
			n.adapt = 4000,
			n.chains = numChains) 
		
		stopCluster(cl)
		procTime[m, K] <- round((proc.time()[3] - t.start)/(60), 1) # track process time
		cat(paste("Process time for model", m, "and hold-out", K, "=", procTime[m, K], "minutes"))
		
		#----------------------------------------------------------------------------
		# Store output
		out.all[[m, K]] <- fit.mK
		saveRDS(fit.mK, file = paste("dataAnalysis/output/model", m, "_holdOut", K, ".rds", sep=""))
		
	}}

#  saveRDS(out.all, file = "dataAnalysis/outAll_Isigma.rds")
# out.all <- readRDS("dataAnalysis/output/outAll_Isigma.rds")

