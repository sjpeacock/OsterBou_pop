###############################################################################
# Function to load and clean data
###############################################################################

#------------------------------------------------------------------------------
# Inputs:
# excl.temp: the temperature trials to be subsetted out when loading the data
#            Defaults to 40*C, where few larvae developed 
#------------------------------------------------------------------------------
# Outputs:
# dat: list of data objects that can be fed into JAGS. The list includes:
#      nT: number of unique temperature treatments
#      T.obs: the observed temperatures in each unique temperature treatment
#      dt: time step for (defaults to 1 day)
#      t: times to simulate model in fitting (day)
#      nt: length of t
#      nt_obs: matrix with number of timestpes f
#------------------------------------------------------------------------------

load.data <- function(excl.temp=40){
  data <- read.csv("dataAnalysis/Ostertagia_data.csv")
  data <- subset(data, data$temperature!=excl.temp)
	
	#--------------------------------------------
	# Data
	dt <- 1 #small time step for simulating, larger time step for fitting.
	t <- seq(1, 120, dt)
	nt <- length(t)
 
	nT <- length(unique(data$temperature))
                    
	# Number of timesteps for each temperature and rep (not the same for all reps!)
	nt_obs <- matrix(nrow=nT, ncol=5)
	for(i in 1:nT){
		for(k in 1:5){
			dik <- unique(data$day[data$temperature==unique(data$temperature)[i]&data$rep==k])
			nt_obs[i,k] <- length(dik)
			#print(paste(dik))
		}}
                    
    # Index of t vector corresponding to each t_obs for temp i and rep k
    ind <- matrix(rep(NA, nT*max(nt_obs)*5)); dim(ind) <- c(nT, max(nt_obs), 5)
    for(i in 1:nT){
      for(k in 1:5){
        for(j in 2:nt_obs[i,k]){
          ind[i,j,k] <- which(t==unique(data$day[data$temperature==unique(data$temperature)[i]&data$rep==k])[j])
        }}}	
    
    # Number of individuals counted
    #n_obs[i,j,k,l,s]
    n_obs <- matrix(rep(NA, nT*max(nt_obs)*5*3*2)); dim(n_obs) <- c(nT, max(nt_obs), 5, 3, 2)
    for(i in 1:nT){
      for(k in 1:5){
        for(j in 1:nt_obs[i,k]){
          for(l in 1:3){
            ikl <- which(data$temperature==unique(data$temperature)[i]&data$rep==k&data$count==l)
            d.ikl <- unique(data$day[ikl])
            iklj <- which(data$temperature==unique(data$temperature)[i]&data$rep==k&data$count==l&data$day==d.ikl[j])
            if(length(iklj)!=4) stop("Length iklj != 4")
            n_obs[i,j,k,l,1] <- sum(data$n_obs[iklj[1:3]]) # Preinfectives are sum of eggs, L1, L2
            n_obs[i,j,k,l,2] <- data$n_obs[iklj[4]] # Infectives are L3
            
          }}}}
    
    dat <- list(
    	nT=length(unique(data$temperature)),
    	T.obs = unique(data$temperature),
    	dt=dt, 
    	t=t, 	# Times that counts were performed
    	nt=nt, 	# Number of times that counts were
    	nt_obs=nt_obs,
    	ind=ind,
    	n_obs=n_obs,
    	nk=5,
    	sigN0=0.01,
    	sigp=0.05
    )
    
    return(dat)
   }

###############################################################################
# Function that returns four parameters over temperature T.all
###############################################################################

#------------------------------------------------------------------------------
# Inputs:

#------------------------------------------------------------------------------
# Outputs:

#------------------------------------------------------------------------------

MTE.predict <- function(model="SSU", params, T.all, exponent = 1){
	T0 <- 15
	k <- 8.62*10^-5
	K <- 273.15
	
	if(model=="A"){
		y <- params['a']*exp(-params['E']/k*(1/(T.all+K)-1/(T0+K)))
	}else if(model=="SSU"){
			y <- params['a']*exp(-params['E']/k*(1/(T.all+K)-1/(T0+K)))*(1+exp(params['Eh']/k*(-1/(T.all+K)+1/(params['Th']+K))))^exponent
	}else if(model=="C"){
		y <- rep(params['a'], length(T.all))
	}
	return(as.numeric(y))
}



###############################################################################
# Simulate model
###############################################################################

#------------------------------------------------------------------------------
# Inputs:
#     times: vector of time points to simulate the model at. 
#            Defaults to c(1:1200)/10 (max 120 days)
#     params: named list of parameters for a single temperature
#             E.g., N0, u0, u1, rho, sigma
#------------------------------------------------------------------------------
# Outputs:
# N_pred: matrix with rows = timesteps (length times) and two columns for 
#         pre-infectives and infectives.
#------------------------------------------------------------------------------

sim.larv <- function(times = c(1:1200)/10, params){
	
		#par <- c(N0, mu0, mu1, rho, sigma)
		
		nt <- length(times)
		dt <- times[2] - times[1]
		
		N_pred <- matrix(NA, nrow=nt, ncol=2)
		# Pre-infectives
		for(j in 1:nt){
			N_pred[j,1] <- params$N0 * exp(-params$u0 * times[j]) * (1 - pnorm(q = (log(times[j] * params$rho) + params$sigma^2 / 2) / (params$sigma)))
			}
			
		# Infectives
		# - this is tricky because you have to keep track of when parasites develop (N_1start) and whether parasites that develop at time i survive to time j (N_1surv)
		
		N1_start <- numeric(nt)
		# Infectives that develop at time "jj"
		for(j in 1:nt){
			N1_start[j] <- params$N0 / (times[j] * sqrt(2 * 3.141593 * params$sigma^2)) * exp(-1 / (2* params$sigma^2) * (log(times[j] * params$rho) + params$sigma^2 / 2)^2 - params$u0 * times[j])
			}
		
		# ...and survive to time j
		N1_surv <- matrix(NA, nrow = nt, ncol = nt)
		for(j in 1:nt){
			for(jj in 1:j){
				N1_surv[j,jj] <- N1_start[jj] * exp(-params$u1 * (times[j] - times[jj]))
				}
			
			# Summed from 1:j	
			N_pred[j,2] <- sum(N1_surv[j,1:j] * dt)
			
			} #end j
	
	return(N_pred)

}

###############################################################################
# Calculate likelihood
###############################################################################

#------------------------------------------------------------------------------
# Inputs:
#     params.temp: matrix containing predicted parameter values (columns; u1, 
#                  u2, rho, sigma) over temperature (rows; one for each 
#                  temperature in dat) 
#    dat: data list (same as fitting) to use for calculating likelihood
#    holdOut: the rep number that was held out for validation. If all data were
#             used then set to NA
#------------------------------------------------------------------------------
# Outputs:
# A list containing two elements:
#     logLik: the negative log liklihood for the training, validation (if 
#             applicable) and overall data
#     N_pred: a list of the predicted number of larvae, each element 
#             corresponding to a temperature and within each element a matrix
#             with nrow = times and ncol = 2 (pre-infectives and infectives)
#------------------------------------------------------------------------------

calcLik <- function(params.temp, dat, holdOut = NA){
	N_pred <- list(); length(N_pred) <- dat$nT
	N_predk <- N_pred

	for(i in 1:dat$nT){
		par.i <- list(50, params.temp[i,1], params.temp[i,2], params.temp[i,3], params.temp[i,4])
		names(par.i) <- c('N0', 'u0', 'u1', 'rho', 'sigma')
		N_pred[[i]] <- sim.larv(times = dat$t, params = par.i)
		}
	

	ll <- array(NA, dim = dim(dat$n_obs))
	
	for(i in 1:dat$nT){ #for each temperature
		for(k in 1:dat$nk){ # for each replicate
			for(l in 1:3){ # for each of three counts
				# Likelihood od initial pre-infectives
				#ll <- sum(ll, dpois(n_obs[i,1,k,l,1], lambda=N0_pred[i,k], log=TRUE))
				for(s in 1:2){ #for each stage
					for(j in 2:dat$nt_obs[i,k]){
						ll[i,j,k,l,s] <- dpois(dat$n_obs[i,j,k,l,s], lambda=N_pred[[i]][dat$ind[i,j,k], s], log=TRUE)
					}}}}}
	
	if(is.na(holdOut) == TRUE){
		logLik <- c(overall = -sum(ll, na.rm = TRUE))
	} else {
		logLik <- c(
			training = -sum(ll[, , which(c(1:dat$nk) != holdOut), , ], na.rm = TRUE), 
			validation = -sum(ll[, , holdOut, , ], na.rm = TRUE),
			overall = -sum(ll, na.rm = TRUE))
	}
	return(list(logLik, N_pred))
	}