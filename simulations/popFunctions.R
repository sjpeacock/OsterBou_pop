###############################################################################################
# Time derivative of the PDE system
###############################################################################################

# Input variables:
# 	y: matrix with columns equal to the spatial grid and rows equal to the different variables: 
#				calf
#				yearling 
#				adult
#				developingL4 
#				arrestedL4 
#				P (adult worm)
#       L0 (free-living pre-infective larvae)
#       L3 (free-living infective larvae)
#       V_developingL4
#       V_arrestedL4
#       V_P

# params: vector of parameters, including:
#     muC - mortality rate of calves per day
#     muY - mortality rate of yearlings per day
#     muA - mortality rate of adults per day
#     alpha - per-parasite rate of parasite-induced mortality of adults
#     beta - intake rate of parasites 
#            (*constant for now but may consider seasonal variation as with 
#            fecal output)
#     ppnInhibit - the proportion of larvae that go into arrested development (0-1)
#     rhoL4 - development rate of L4 larvae to adults (per day)
#     muP - mortality of adult parasites per day
#     *lambda* - time-varying rate of egg output per adult parasite
#     gamma - density dependence of parasite fecundity (-0.49)
#     *muL0* - mortality rate (per day) of pre-infective larvae
#     *rho0* - development rate (per day) of pre-infectives to infectives
#     *muL3* - mortality
# * parameters that will vary over time (temperature, etc.)

#------------------------------------------------------------------------------
partial_t.Bou <- function(y, p){ 
	
	#*********
	# 1) Number of stationary and moving calves
	d.calf_stat <- (- p$muC * y['calf_stat', ] - p$start * y['calf_stat', ]  + (p$stop + p$parasitStop * y['P_mov', ]) * y['calf_mov', ]) 
	
	d.calf_mov <- (- p$muC * y['calf_mov', ] - (p$stop + p$parasitStop * y['P_mov', ]) * y['calf_mov', ] + p$start * y['calf_stat', ]) 
	
	#*********
	# 2) Number of yearlings
	d.yearling_stat <- (- p$muY * y['yearling_stat', ] - p$start * y['yearling_stat', ] + (p$stop + p$parasitStop * y['P_mov', ]) * y['yearling_mov', ])
	
	d.yearling_mov <- (- p$muY * y['yearling_mov', ] - (p$stop + p$parasitStop * y['P_mov', ]) * y['yearling_mov', ] + p$start * y['yearling_stat', ]) 
	
	#*********
	# 3) Number of adults
	d.adult_stat <- (- p$muA * y['adult_stat', ] - p$alpha * y['P_stat', ] * y['adult_stat', ] - p$start * y['adult_stat', ] + (p$stop + p$parasitStop * y['P_mov', ]) * y['adult_mov', ]) 
	
	d.adult_mov <- (- p$muA * y['adult_mov', ] - p$alpha * y['P_mov', ] * y['adult_mov', ] - (p$stop + p$parasitStop * y['P_mov', ]) * y['adult_mov', ] + p$start * y['adult_stat', ]) 
	
	#*********
	if((p$stop + p$parasitStop) == 0){
		d.L4_stat <- (p$beta * (1 - p$ppnInhibit) * y['L3', ] - (p$rho4 + p$mu4) * y['L4_stat', ])
		
		d.L4A_stat <- (p$beta * p$ppnInhibit * y['L3', ] -  p$mu4 * y['L4A_stat', ])
		
		d.P_stat <- (p$rho4 * y['L4_stat', ] - p$muP * y['P_stat', ] - p$nuP * y['P_stat', ] * (y['P_stat', ] * (p$k + 1)/ p$k + 1) - p$alpha * y['P_stat',] * (y['P_stat', ] + p$k) / p$k)
		
	} else {
		d.L4_stat <- (p$beta * (1 - p$ppnInhibit) * y['L3', ] - (p$rho4 + p$mu4) * y['L4_stat', ] + (p$stop + p$parasitStop * y['P_mov', ]) * y['adult_mov',]/y['adult_stat',] * (y['L4_mov', ] - y['L4_stat', ]))
		
		d.L4A_stat <- (p$beta * p$ppnInhibit * y['L3', ] -  p$mu4 * y['L4A_stat', ]  + (p$stop + p$parasitStop * y['P_mov', ]) * y['adult_mov',]/y['adult_stat',] * (y['L4A_mov', ] - y['L4A_stat', ]))
		
		d.P_stat <- (p$rho4 * y['L4_stat', ] - p$muP * y['P_stat', ] - p$nuP * y['P_stat', ] * (y['P_stat', ] * (p$k + 1)/ p$k + 1) - p$alpha * y['P_stat',] * (y['P_stat', ] + p$k) / p$k + y['adult_mov', ]/y['adult_stat', ] * (p$stop * (y['P_mov', ] - y['P_stat', ]) + p$parasitStop * y['P_mov', ] * (p$k * (y['P_mov', ] - y['P_stat', ] + 1) + y['P_mov', ])/ p$k))
		
	}
	
	#*********
	if(p$start == 0){
		d.L4_mov <- (p$beta * (1 - p$ppnInhibit) * y['L3', ] - (p$rho4 + p$mu4) * y['L4_mov', ])
		
		d.L4A_mov <- (p$beta * p$ppnInhibit * y['L3', ] - p$mu4 * y['L4A_mov', ])
		
		d.P_mov <- (p$rho4 * y['L4_mov',] - p$muP * y['P_mov',]	- p$nuP * y['P_mov', ] * (y['P_mov', ] * (p$k + 1)/ p$k + 1) 	- (p$alpha + p$parasitStop) * y['P_mov',] * (y['P_mov', ] + p$k) / p$k)
		
	} else {
		d.L4_mov <- (p$beta * (1 - p$ppnInhibit) * y['L3', ] - (p$rho4 + p$mu4) * y['L4_mov', ] + p$start *y['adult_stat',] /  y['adult_mov',] * (y['L4_stat', ] - y['L4_mov', ]))
		
		d.L4A_mov <- (p$beta * p$ppnInhibit * y['L3', ] - p$mu4 * y['L4A_mov', ] + p$start * y['adult_stat',] /  y['adult_mov',] * (y['L4A_stat', ] - y['L4A_mov', ]))
		
		d.P_mov <- (p$rho4 * y['L4_mov',] - p$muP * y['P_mov',]	- p$nuP * y['P_mov', ] * (y['P_mov', ] * (p$k + 1)/ p$k + 1) 	- (p$alpha + p$parasitStop) * y['P_mov',] * (y['P_mov', ] + p$k) / p$k + p$start * y['adult_stat', ]/y['adult_mov', ] * (y['P_stat', ] - y['P_mov', ]))
	}

	# 7) Pre-infective free-living larvae
	d.L0 <- (p$lambda * (y['adult_stat',] * y['P_stat',] ^ (1 + p$phi) + y['adult_mov',] * y['P_mov',] ^ (1 + p$phi)) - p$mu0 * y['L0', ] - p$rho0 *  y['L0', ]) # development
	
	# 8) Infective free-living L3 larvae
	d.L3 <- (p$rho0 *  y['L0', ] - p$mu3 * y['L3', ] - p$beta * y['L3', ] * (y['adult_stat', ] + y['adult_mov', ]))
	
	return(rbind(
		d.calf_mov,
		d.yearling_mov,
		d.adult_mov,
		d.L4_mov,
		d.L4A_mov,
		d.P_mov,
		d.calf_stat,
		d.yearling_stat,
		d.adult_stat,
		d.L4_stat,
		d.L4A_stat,
		d.P_stat,
		d.L0,
		d.L3
	))
	
} #end function


###############################################################################################
# Parasite egg output - lambda
###############################################################################################
# From Stien et al. 2002 Int J Parasit

predict.lambda <- function(DOY){
	
	# Eggs per gram feces per worm, not accounting for density dependence
	# Lambda
	alpha1 <- 0.01
	alpha2 <- 0.345
	mu <- 0.52
	sigma <- 0.087
	
	lambda <- alpha1 + alpha2/(sigma*sqrt(2*pi)) * exp(-(DOY/365 - mu)^2/(2*sigma^2))
	
	# Faeces production rate
	muF <- 0.58 # peak plant biomass in august, 58% thorugh the year
	alphaF1 <- 1300 # faecal production rate in winter (min) based on 1 kg dry matter per day
	maxF <- 5400 # g faeces per day in summer
	sigmaF <- 0.085
	alphaF2 <- (maxF - alphaF1) * sigmaF *sqrt(2*pi) 
	
	fpr <- alphaF1 + alphaF2/(sigmaF*sqrt(2*pi)) * exp(-(DOY/365 - muF)^2/(2*sigmaF^2))
	
	return(lambda * fpr)
	
}


###############################################################################################
# MTE predictions for free-living larvae params
###############################################################################################

# predict.mu0 <- function(
# 	params = c(a = 0.068, E = 0.884),
# 	temp = 15){
# 		return(params['a'] * exp(-params['E']/(8.62*10^-5) * (1/(temp + 273.15) - 1/(15 + 273.15))))
# }

# SSUL
# mu_0*exp(-E_mu/k*(1/(temp[i]+273.15)-1/(15+273.15)))*(1+exp(EL_mu/k*(1/(temp[i]+273.15)-1/(TL_mu+273.15)))+exp(EH_mu/k*(-1/(temp[i]+273.15)+1/(TH_mu+273.15))))

# predict.mu0 <- function(
# 	params = c(a = 0.068, E = 0.884, El = 3.25, Tl = 0),
# 	temp = 15){
#  	return(params['a'] * exp(-params['E']/(8.62*10^-5) * (1/(temp + 273.15) - 1/(15 + 273.15))) * (1  + exp(params['El']/(8.62*10^-5)*(1/(temp+273.15) - 1 / (params['Tl'] + 273.15)))))
# }
# 
# # plot(seq(-30, 30, 0.1), predict.mu0(temp = seq(-30, 30, 0.1)), "l", ylim = c(0, 0.5))
# 
# # muL3 constant across temperature according to best-fit MTE relationship
# mu3 = 0.022
# 	
# predict.rho0 <- function(
# 	params = c(a = 0.032, E = 0.686, Eh = 7.957, Th = 30.568),
# 	temp = 15){
# 	return(params['a'] * exp(-params['E']/(8.62*10^-5) * (1/(temp + 273.15) - 1/(15 + 273.15))) * (1 + exp(params['Eh'] / (8.62*10^-5)*(-1/(temp+273.15)+1/(params['Th']+273.15))))^(-1))
# }

# From previous fits:

predict.mu0<-function(temp){
	0.0548*exp(-0.971/(8.62*10^-5)*(1/(temp+273.15)-1/(15+273.15)))*(1+exp(30.244/(8.62*10^-5)*(-1/(temp+273.15)+1/(39.2+273.15))))
}

predict.mu3<-function(temp){
	0.0211*exp(-(-0.208)/(8.62*10^-5)*(1/(temp +273.15)-1/(15+273.15)))*(1+exp(3.5543/(8.62*10^-5)*(-1/(temp +273.15)+1/(27.6+273.15))))
}

predict.rho0<-function(temp){
	0.0287*exp(-0.537/(8.62*10^-5)*(1/(temp+273.15)-1/(15+273.15)))*(1+exp(3.297/(8.62*10^(-5))*(-1/(temp+273.15)+1/(31.6+273.15))))^(-1)
}


###############################################################################################
# Parasite and host dependent fecundity
###############################################################################################
# Inputs:
#     P_mean: mean parasite burden at conception (October, 222 daysprior to calving)
#     P_VMR: variance-to-mean ratio for P_mean (VMR = mean/k + 1, k = mean / (VMR - 1))
#     numAdults: number of adult female caribou that survive to calving season
#     pCalf0: fecundity (i.e., probability that female has a calf) in the absence of parasites

numCalves <- function(P_mean, numFemales, pCalf0){ #, stoch = FALSE)
 
# 	# Stochastic approach
#  if(stoch == TRUE){
#  	dispersionParam <- P_mean/(P_VMR - 1)
# 	numParasitesPerFemale <- rnbinom(n = round(numFemales), size = dispersionParam, mu = P_mean)
# 	calf <- rbinom(
# 		n = round(numFemales),
# 		size = 1,
# 		prob = pCalf0 - 1/(1 + exp(7.025 - 0.000328*numParasitesPerFemale)))
# 	return(sum(calf))
#  
# 	} else{

	# Deterministic approach
	return(numFemales * (pCalf0 - 1/(1 + exp(7.025 - 0.000328 * P_mean))))

}

# # Plot to compare stochastic and deterministic approaches
# P_mean.all <- seq(0, 10000, 100)
# calves.all <- matrix(NA, ncol = 2, nrow = length(P_mean.all))
# for(i in 1:length(P_mean.all)){
# 	calves.all[i, 1] <- numCalves(P_mean = P_mean.all[i], P_VMR = 8, numFemales = 1000, pCalf0 = 0.9, stoch = FALSE)
# 	calves.all[i, 2] <- numCalves(P_mean = P_mean.all[i], P_VMR = 8, numFemales = 1000, pCalf0 = 0.9, stoch = TRUE)
# }
# 
# plot(P_mean.all, calves.all[, 1], "l", ylim = range(calves.all, na.rm = TRUE), xlab = "Mean parasite burden", ylab = "No. calves per 1000 females")
# lines(P_mean.all, calves.all[,2], lty = 3)
# 
	
###############################################################################################
# Annual temperature cycle 
###############################################################################################
# Input variables:
# 	DOY: vector of julien days 1:365 for which predicted temperature is desired
#   range: one of "Summer",  "Winter",  "Spring",  "Fall", "Calving", which determines the 
#          parameters used to predict temperature based on a cosine curve fitted to data from
#          the respective Bathurst caribou ranges (Russell et al. 2013 Rangifer)
# Output variables:
#   vector of predicted temperatures in degrees Celsius for each DOY input

predict.temp <- function(DOY, range){
	rangeNum <- as.numeric(factor(range, levels = c("Summer",  "Winter",  "Spring",  "Fall", "Calving")))
	# mean temperature
	ck <- c(-10.974553,  -8.029248,  -9.808518,  -9.750996, -11.466570)[rangeNum] 
	# half annual temperature range (amplitude)
	dk <- c(22.546752,  23.653142,  22.808474,  22.990450, 21.577856)[rangeNum] 
	# DOY corresponding to max temperature
	t0 <- c(199.922313, 198.146132, 199.304408, 199.167954, 201.066183)[rangeNum] 
	
	return(ck + dk * cos((DOY - t0)* 2 * pi / 365))
}

###############################################################################################
# Calculate parameters based on DOY
###############################################################################################

calcParams <- function(DOY, temp, ppnInhibit = 0){
	# Need to have parameters as a list because the parameters for stationary larvae will vary in space and time
	params <- list(
	# muC - mortality rate of calves per day
	muC = (1 - 0.45)/365, #approx. initial parameter from Boulanger
	
	# muY - mortality rate of yearlings per day
	muY = (1 - 0.86)/365, # annual Sy = 0.86 from Boulanger et al. 2011
	
	# muA - mortality rate of adults per day
	muA = (1 - 0.86)/365, #approx. initial parameter from Boulanger
	
	# alpha - per-parasite rate of parasite-induced mortality of adults
	alpha = 0,
	
	# rate of starting
	start = startMat[DOY, ],
	
	# rate of stopping
	stop = stopMat[DOY, ],
	
	# per-parasite rate of parasite-induced stopping
	parasitStop = 0,
	
	# beta - intake rate of parasites 
	# *********  Need to better resolve this. ****************************
	# Also depends on dry matter intake and will vary thorughout the year?
	# Assume constant for now. Grenfell et al. 1987 assumed three levels (0.0001, 0.001, 0.01)
	# Seems common to do that so we may just need to look at sensitivity to this parameter.
	beta = 10^-6, #0.0000025
	
	# ppnInhibit - the proportion of larvae that go into arrested development (0-1)
	ppnInhibit = ppnInhibit,
	
	# rho4 - development rate of L4 larvae to adults (per day)
	# From Grenfell et al. 1987, development to adults can take 17 days.
	rho4 = 0.06, 
	
	# mu4 - mortality of L4 larvae per day
	mu4 = 0.02,
	
	# muP - mortality of adult parasites per day
	# Likely density-dependent, see Smith and Grenfell 1985 Parasit. Today.
	# From Grenfell et al. 1987: mu5 = a + b * P where 
	#                            a = 0.1713 per day and b = 0.3082 * 10^-6 per worm per day
	# plot(seq(1, 10^18, length.out = 100), 0.1713 + 0.3082 * 10^-6 *  seq(1, 10^18, length.out = 100), "l")
	# Likely insignificant over the ranges of parasites that we see, use mean
	muP = 0.1713,
	nuP = 0.3082e-6,
		
	# *lambda* - time-varying rate of egg output per adult parasite
	lambda = predict.lambda(DOY)*10^-2,
	
	# phi - density dependence of parasite fecundity (-0.49 Stien et al. 2002 Int J Parasit) 
	phi = -0.49,
	
	# *mu0* - mortality rate (per day) of pre-infective larvae
	mu0 = predict.mu0(temp),
	
	# *rho0* - development rate (per day) of pre-infectives to infectives
	rho0 = predict.rho0(temp),
	
	# *mu3* - mortality
	#          Estimated as constant, but apply in matrix to allow for changes.
	mu3 = predict.mu3(temp),#rep(0.022, length(x)),
	
	# Aggregation parameter
	# ** Need to look at what this is in Bathurst data
	k = 0.8
	)
	
	return(params)
}

###############################################################################################
# Function to set up initial distribution
###############################################################################################
# Initial spatial distribution
initDist <- function(totPop, x, x.start.sd = 50){
	# Note the shift so that population always starts at x = 0
	shift.x <- round(length(x)/2)
	return(c(totPop / sqrt(2 * pi * (x.start.sd^2)) * exp(- (x - x[shift.x])^2 / (2 * x.start.sd^2)))[c(shift.x:length(x), 1:(shift.x - 1))])
}

###############################################################################################
# Function to set up space-time grid
###############################################################################################
gridSetup <- function(
	y, # Number of years to run simulation
	dt # Timestep in days
	){
	
	#------------------------------------------------------------------------------
	# Time grid
	#------------------------------------------------------------------------------
	timeDat <- data.frame(
		year = rep(1:y, each = length(DOY)/dt),
		DOY = rep(rep(DOY, each = 1/dt), y),
		time = rep(seq(dt, 365, dt), y),
		step = c(1:(length(DOY) / dt * y))
	)
	
	#------------------------------------------------------------------------------
	# Spatial grid
	#------------------------------------------------------------------------------
	
	xmax <- 2268 #sum(bouMove) # 2268
	dx <- migSpeed * dt
	x <- seq(0, xmax, dx)
	
	# Grid spaces that fall on each "range" based on caribou movement...I did this manually.
	xRange <- character(length(x))
	for(i in 1:length(x)){
		j <- findInterval(x[i], cumsum(bouMove))
		xRange[i] <- bouRange[j]
	}
	length(x) == length(xRange)
	
	#------------------------------------------------------------------------------
	# Temperature grid
	#------------------------------------------------------------------------------
	tempGrid <- matrix(NA, nrow = 365, ncol = length(x))
	for(i in 1:length(x)){
		tempGrid[, i] <- predict.temp(DOY, range = xRange[i])
	}
	
	# Smooth each day
	tempGrid.smooth <- tempGrid
	midX <- round(length(x)/2)
	for(n in 1:length(DOY)){
		smooth1 <- lowess(tempGrid[n, c((midX + 1):length(x), 1:midX)], f = 0.3)$y[order(c((midX+1):length(x), 1:midX))]
		smooth2 <- lowess(tempGrid[n, ], f = 0.3)$y
		dsmooth <- smooth1-smooth2
		ind1 <- which(abs(dsmooth[1:round(length(x)/2)]) < 10^-2)
		ind3 <- which(abs(dsmooth[(length(x) - round(length(x)/2.5)):(length(x) - round(length(x)/4))]) < 10^-2) + (length(x) - round(length(x)/2.5)) - 1
		
		tempGrid.smooth[n, ] <- c(smooth1[1:ind1[1]], smooth2[(ind1[1] + 1): ind3[1]], smooth1[(ind3[1] + 1):length(x)]) 
	}
	# 14 km grid
	
	return(list(
		timeDat = timeDat,
		x = x,
		dx = dx,
		xRange = xRange, 
		tempGrid = tempGrid.smooth
		))
} 
###############################################################################################
# Function to simulate caribou dynamics within a season
###############################################################################################

simBou <- function(
	timeDat,
	x,
	dx, 
	dt,
	migSpeed,
	initBou, 
	tempGrid,
	# Breeding date, when animals move up a class, is June 7 (DOY = 158)
	breedDOY = as.numeric(strftime(as.Date("1985-06-07"), "%j")),
	# L4 resume development at the start of spring migration
	L4startDOY = as.numeric(strftime(as.Date("1985-04-20"), "%j"))
  ){
	
	nt <- dim(timeDat)[1]
	
	# 1) Set up matrices to store solutions every 7 days
	V <- array(NA, dim = c(dim(initBou)[1], length(x), nt), dimnames = list(rownames(initBou), x, paste(timeDat$year, timeDat$time, sep="-")))
	V[, , 1] <- initBou
	V0 <- array(0, dim = c(dim(initBou)[1], length(x)), dimnames = list(rownames(initBou), x))
	
	# 2) Run through each timestep
	for(n in 1:(nt - 1)){
		# Set parameters based on DOY
		params.n <- calcParams(DOY = timeDat$DOY[n], temp = tempGrid[timeDat$DOY[n], ])
		
		# Try manually setting these to 1 just for a single day??
		if(is.element(round(timeDat$time[n], 1), round(c(99, 168, 250), 1)) == TRUE){
			params.n$start <- 1/dt 
			} else {
				params.n$start <- 0
			}
		
		if(is.element(round(timeDat$time[n], 1), round(c(152, 181, 346), 1)) == TRUE){
			params.n$stop <- 1/dt
		} else {
			params.n$stop <- 0
		}
		
		# Advection speed for each variable: number of grid spaces moved
		u <- migSpeed * dt / dx
		
		# Calculate boundary conditions: torus for circular migration
		Vn <- V[, , n]
		Vnp1 <- Vn
		
		# Spatial advection (upstream differencing) for moving stages
		if(u > 0){
			for(j in match(c("calf_mov", "yearling_mov", "adult_mov", "L4_mov", "L4A_mov", "P_mov"), rownames(initBou))){
				Vnp1[j, ] <- Vn[j, c((dim(V)[2] - u + 1) : (dim(V)[2]), 1:(dim(V)[2] - u))]
			}
		}
		
		#---------------------------------------------------------------------------
		# If breeding date, move caribou up and add calves
		# *** we're going to get weird things happening if the population doesn't mix...
		if(timeDat$time[n] == breedDOY){
			cat("breeding")
			Vn.breed <- Vnp1
			
			newCalves_stat <- numCalves(
				P_mean = V['P_stat', , n],
				numFemales = V['adult_stat', , n] * propFemale,
				pCalf0 = 0.85)
			
			newCalves_mov <- numCalves(
				P_mean = V['P_mov', , n],
				numFemales = V['adult_mov', , n] * propFemale,
				pCalf0 = 0.85)
			
			# # All calves start out in stat category
			Vn.breed['calf_stat', ] <- newCalves_stat + newCalves_mov
			Vn.breed['calf_mov', ] <- rep(0, length(x))
			
			# Yearlings and adults stay in respective categories
			Vn.breed['yearling_stat', ] <- Vnp1['calf_stat', ]
			Vn.breed['yearling_mov', ] <- Vnp1['calf_mov', ]
			Vn.breed['adult_stat', ] <- Vnp1['adult_stat', ] + Vnp1['yearling_stat', ]
			Vn.breed['adult_mov', ] <- Vnp1['adult_mov', ] + Vnp1['yearling_mov', ]
			
			# Mean parasite burden decreases as incoming yearlings are parasite-free
			Vn.breed['P_mov', ] <- (Vnp1['adult_mov', ] * Vnp1['P_mov', ]) / Vn.breed['adult_mov', ]
			Vn.breed['P_stat', ] <- (Vnp1['adult_stat', ] * Vnp1['P_stat', ]) / Vn.breed['adult_stat', ]
			
			# Replace Vnp1 with updated matrix
			Vnp1 <- Vn.breed
		} # end if breed
		#---------------------------------------------------------------------------
		
		#---------------------------------------------------------------------------
		# If start of spring migration, L4 resume development
		if(round(timeDat$time[n], 3) == round(L4startDOY, 3)){
			cat("L4 development resuming")
			Vn.start <- Vnp1
			
			Vn.start['L4_stat', ] <- Vnp1['L4_stat', ] + Vnp1['L4A_stat', ]
			Vn.start['L4A_stat', ] <- rep(0, length(x))
			Vn.start['L4_mov', ] <- Vnp1['L4_mov', ] + Vnp1['L4A_mov', ]
			Vn.start['L4A_mov', ] <- rep(0, length(x))
			
			Vnp1 <- Vn.start
		}
		#---------------------------------------------------------------------------
		
		#---------------------------------------------------------------------------
		# Temporal dynamics (4th order Runge Kutta)
		# Was having issues here so just using Euler's formula right now
		k1 <- partial_t.Bou(y = Vnp1, p = params.n)
		# which(is.na(k1) == TRUE, arr.ind = TRUE)
		# k2 <- partial_t.Bou(y = Vnp1 + k1 * (dt / 2), p = params.n)
		# k3 <- partial_t.Bou(y = Vnp1 + k2 * (dt / 2) , p = params.n)
		# k4 <- partial_t.Bou(y = Vnp1 + k3 * dt, p = params.n)
		# 
		# par(mfrow = c(3,1))
		# plot(x, k1['d.L4_stat', ], "l", lty = 2, ylim = range(k1[c('d.L4_stat', 'd.L4_mov', 'd.L4A_stat', 'd.L4A_mov'),], na.rm = TRUE))
		# lines(x, k1['d.L4_mov', ])
		# lines(x, k1['d.L4A_stat', ], col = 2, lty = 2)
		# lines(x, k1['d.L4A_mov', ], col = 2)
		# 
		# plot(x, k1['d.P_stat', ], "l", col = 4, lty = 2)
		# lines(x, k1['d.P_mov', ], col = 4)
		# 
		# plot(x, k1['d.L0', ], "l")
		# lines(x, k1['d.L3', ], "l", col = 2)
		# #
		V[, , n + 1] <- pmax(V0, Vnp1 + dt * k1)#pmax(V0, Vnp1 + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4))#  #
		
		if(timeDat$DOY[n] == 365 & timeDat$time[n] == timeDat$DOY[n]) cat(paste("Year", timeDat$year[n], "complete\n"))
		
		
	} #end all timesteps	
	
	return(V)
}
###############################################################################
# Plot output
###############################################################################

plot.timestep <- function(V, n, Nrange = NA, Prange = NA, Lrange = NA){
	# Nrange <- range(V[c('adult_mov', "adult_stat"), , ]
	if(is.na(Nrange[1]) == TRUE) Nrange <- range(V[c('adult_mov', "adult_stat"), , n])								
	if(is.na(Prange[1]) == TRUE) Prange <- range(V[c('P_mov', "P_stat"), , n])								
	if(is.na(Lrange[1]) == TRUE) Lrange <- range(V[c('L0', "L3"), , n])
	
	par(mfrow = c(3,1), mar = c(2,5,1,4), oma = c(2, 0, 2, 0))
	
	plot(x, V['adult_mov', , n], "l", ylim = Nrange, bty = "l", xaxt="n", yaxt = "n", ylab = "", lwd = 1.5)
	axis(side = 1, labels = FALSE)
	yTick <- pretty(Nrange)
	axis(side = 2, las = 1, at = yTick, labels = yTick/1000)
	lines(x, V['adult_stat', , n], lty = 3, lwd = 1.5)
	lines(x, V['calf_mov', , n], col = seasonCols['Calving'])
	lines(x, V['calf_stat', , n], lty = 3, col = seasonCols['Calving'])
	lines(x, V['yearling_mov', , n], col = seasonCols['Fall'])
	lines(x, V['yearling_stat', , n], lty = 3, col = seasonCols['Fall'])
	mtext(side = 3, adj = 0, "a) Host population density (1000s per km)")
	
	plot(x, V['P_mov', , n], "l", ylim = Prange, bty = "l", ylab = "", lwd = 1.5, las = 1, xaxt="n")
	axis(side = 1, labels = FALSE)
	lines(x, V['P_stat', , n], lty = 3, lwd = 1.5)
	mtext(side = 3, adj = 0, "b) Mean parasite burden per host")
	
	plot(x, V['L3', , n], "l", ylim = Lrange, bty = "l", ylab = "", lwd = 1.5, las = 1)#, yaxt="n")
	# yTick <- pretty(range(V[c('L0', "L3"), , ], na.rm = TRUE))
	# axis(side = 2, las = 1, at = yTick, labels = yTick*10^-10)
	lines(x, V['L0', , n], lwd = 1.5, col = seasonCols['Winter'])
	mtext(side = 3, adj = 0, "c) Density of free-living parasite larvae (per km)")
	
	D <- as.Date(paste((1984 + timeDat$year[n]), timeDat$DOY[n], sep="-"), format = "%Y-%j")
	mtext(side = 3, adj = 1, line = -1, outer = TRUE, paste("Year ", timeDat$year[n], "\n", strftime(D, format = "%b %d"), "\n timestep ", n, sep =""))
}
