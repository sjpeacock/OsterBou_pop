###############################################################################################
# Scenario for proportion inhibiting
###############################################################################################

pArrestCalc <- function(scenario){
	if(scenario == 4){
		pIn <- c(rep(1, 109), 1/(1 + exp(-0.08*(c(110:365) - 172))))
	} else {
		pIn <- rep(c(1, 0.5, 0)[p], 365)
	}
	return(pIn)
}

###############################################################################################
# Time derivative of the PDE system
###############################################################################################

# Input variables:
# 	y: matrix with columns equal to the spatial grid and rows equal to the different variables: 
#				calf_stat
#				calf_mov
#				yearling_stat
#				yearling_mov
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

# #------------------------------------------------------------------------------

partial_t.Bou <- function(y, p){ 
	
	nx <- dim(y)[2]
	dy <- array(0, dim = dim(y), dimnames = dimnames(y))
	
	#-----------------------------------------------------------------------------
	# If there are moving adults, then stopping can be non-zero
	adult_mov.nonzero <- which(y['adult_mov', ] > 0)
	adult_stat.nonzero <- which(y['adult_stat', ] > 0)
	
	#*****************************************************************************
	# 1) Number of stationary and moving calves
	dy['calf_stat', ] <- - (p$muC + p$start) * y['calf_stat', ] + p$stop * y['calf_mov', ]
	dy['calf_mov', ]  <- - (p$muC + p$stop) * y['calf_mov', ]  + p$start * y['calf_stat', ]
	
	#*****************************************************************************
	# 2) Number of yearlings
	dy['yearling_stat', ] <- - (p$muY + p$start) * y['yearling_stat', ] + p$stop * y['yearling_mov', ]
	dy['yearling_mov', ] <-  - (p$muY + p$stop) * y['yearling_mov', ]  + p$start * y['yearling_stat', ] 
	
	#*****************************************************************************
	# 3) Number of adults
	dy['adult_stat', ] <- - (p$muA + p$start) * y['adult_stat', ] + p$stop * y['adult_mov', ]
	dy['adult_mov', ] <- - (p$muA + p$stop) * y['adult_mov', ] + p$start * y['adult_stat', ] 
	
	#*****************************************************************************
	# 4) Arrested larvae
	dy['L4A_stat', ] <- p$beta * p$ppnInhibit * y['L3', ] * y['adult_stat', ] - (p$muA + p$mu4 + p$start) * y['L4A_stat', ] + p$stop * y['L4A_mov', ]
	dy['L4A_mov', ] <- p$beta * p$ppnInhibit * y['L3', ] * y['adult_mov', ] - (p$muA + p$mu4 + p$stop) * y['L4A_mov', ] + p$start * y['L4A_stat', ]
	
	#*****************************************************************************
	# 5) Developing larvae
	dy['L4_stat', ] <- p$beta * (1 - p$ppnInhibit) * y['L3', ] * y['adult_stat', ] - (p$muA + p$mu4 + p$rho4 + p$start) * y['L4_stat', ] + p$stop * y['L4_mov', ]
	dy['L4_mov', ] <- p$beta * (1 - p$ppnInhibit) * y['L3', ] * y['adult_mov', ] -  (p$muA + p$mu4 + p$rho4 + p$stop) * y['L4_mov', ] + p$start * y['L4_stat', ]
	
	#*****************************************************************************
	# 6) Adult parasites
	PMort_stat <- rep(0, nx)
	PMort_stat[adult_stat.nonzero] <- y['P_stat', adult_stat.nonzero] / y['adult_stat', adult_stat.nonzero] * (p$k + 1)/p$k + 1
	
	dy['P_stat', ] <- p$rho4 * y['L4_stat', ] - (p$muA + p$muP + (p$nuP * PMort_stat) + p$start) * y['P_stat', ] + p$stop * y['P_mov', ]
	
	PMort_mov <- rep(0, nx)
	PMort_mov[adult_mov.nonzero] <- y['P_mov', adult_mov.nonzero] / y['adult_mov', adult_mov.nonzero] * (p$k + 1)/p$k + 1
	
	dy['P_mov', ] <- p$rho4 * y['L4_mov', ] - (p$muA + p$muP + (p$nuP * PMort_mov) + p$stop) * y['P_mov', ] + p$start * y['P_stat', ]
	
	#*****************************************************************************
	# 7) Larvae
	dy['L0', ] <- p$lambda * (y['adult_stat', ] * y['P_stat', ]^(1 + p$phi) + y['adult_mov', ] * y['P_mov', ]^(1 + p$phi)) - (p$mu0 + p$rho0) * y['L0', ]
	
	dy['L3', ] <- p$rho0 * y['L0', ] - p$mu3 * y['L3', ] - p$beta * y['L3', ]* (y['adult_stat', ] + y['adult_mov', ])
	
	#*****************************************************************************
	if(sum(is.na(dy)) > 0) stop("\n\nNAs in derivative function.\n\n")
	unique(which(is.na(dy) == TRUE, arr.ind = TRUE)[, 1])
	return(dy)
	
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

predict.mu0 <- function(temp){
 	return(0.068 * exp(-0.884/(8.62*10^-5) * (1/(temp + 273.15) - 1/(15 + 273.15))) * (1  + exp(2.928/(8.62*10^-5)*(1/(temp+273.15) - 1 / (-3.377 + 273.15)))))
}

predict.rho0 <- function(temp){
	return(0.032 * exp(-0.686/(8.62*10^-5) * (1/(temp + 273.15) - 1/(15 + 273.15))) * (1 + exp(7.957/ (8.62*10^-5)*(-1/(temp+273.15)+1/(30.568+273.15))))^(-1))
}

predict.mu3 <- function(temp){
	0.0211*exp(-0.208/(8.62*10^-5)*(1/(temp +273.15)-1/(15+273.15)))*(1 + exp(3.409 / (8.62*10^-5)*(1/(temp + 273.15)-1/(-19.318 + 273.15))) + exp(3.5543/(8.62*10^-5)*(-1/(temp +273.15)+1/(27.6+273.15))))
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
	return(numFemales * pCalf0 * (1 - 1/(1 + exp(7.025 - 0.000328 * P_mean))))

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

predict.temp <- function(timeDat = NULL, climateScenario = NULL, ground = TRUE, stoch = FALSE){ 
	
	# If no timeDat is supplied, calculate temperatures over one year with timestep of 1 day
	if(length(timeDat) == 0){
		DOY <- c(1:365)
		y <- 1
		dt <- 1
		timeDat <- data.frame(
			year = rep(1:y, each = length(DOY)/dt),
			DOY = rep(rep(DOY, each = 1/dt), y),
			time = rep(seq(dt, 365, dt), y),
			step = c(1:(length(DOY) / dt * y))
		)
		}
	
	if(stoch == FALSE){
		# rangeNum <- as.numeric(factor(range, levels = c("Summer",  "Winter",  "Spring",  "Fall", "Calving")))
		# mean temperature
		ck <- c(-10.974553,  -8.029248,  -9.808518,  -9.750996, -11.466570) 
		# half annual temperature range (amplitude)
		dk <- c(22.546752,  23.653142,  22.808474,  22.990450, 21.577856) 
		# DOY corresponding to max temperature
		t0 <- c(199.922313, 198.146132, 199.304408, 199.167954, 201.066183)#[rangeNum] 
		
		tempDOY <- matrix(NA, nrow = length(DOY), ncol = 5)
		for(i in DOY){
			tempDOY[i, ] <- ck + dk * cos((DOY[i] - t0)* 2 * pi / 365)
		}
		
		# If introducing a climate change scenario
		if(length(climateScenario) > 0){
			if(climateScenario == "rcp26"){
				increase <- 2.2176 - 0.6511 * cos((DOY - 168.5002)* 2 * pi / 365)
				tempDOY <- tempDOY + matrix(rep(increase, 5), nrow = length(DOY), ncol = 5)
			} else if(climateScenario == "rcp85"){
				increase <- 7.7867 - 3.1112 * cos((DOY - 179.0068)* 2 * pi / 365)
				tempDOY <- tempDOY + matrix(rep(increase, 5), nrow = length(DOY), ncol = 5)
			}
		}
		
		# If adjusting for ground temperatures
		if(ground == TRUE){
			# Ground temp
			tempDOY_ground <- tempDOY
			
			for(i in 1:5){ # for each season
				# Minimum ground temperatures of -15*C
				tempDOY_ground[which(tempDOY[, i] < -15), i] <- -15
			
				# Delay DD accumulation by 2 weeks
				d0 <- which(tempDOY_ground[, i] > 0)[1] # what is the DOY when temps first > 0
				tempDOY_ground[d0:(d0 + 14), i] <- 0 # Make ground temps =0 for next two weeks
				
				# 5 *C warmer in the summer
				summerDOY <- c((d0 + 15):max(which(tempDOY[, i] > 0)))
				if(length(climateScenario) == 0){
					tempDOY_summer <- ck[i] + (dk[i] + 5) * cos((summerDOY - t0[i])* 2 * pi / (365 - 14))
					} else if(climateScenario == "rcp26"){
					tempDOY_summer <- ck[i] + (dk[i] + 5) * cos((summerDOY - t0[i])* 2 * pi / (365 - 14)) + 2.2176 - 0.6511 * cos((summerDOY - 168.5002)* 2 * pi / 365)
				} else if(climateScenario == "rcp85"){
					tempDOY_summer <- ck[i] + (dk[i] + 5) * cos((summerDOY - t0[i])* 2 * pi / (365 - 14)) +  7.7867 - 3.1112 * cos((summerDOY - 179.0068)* 2 * pi / 365)
				} 
				
				# Smooth out start of summer over first 30 days
				adjustDown <- seq(- tempDOY_summer[1], 0, length.out = 30)
				tempDOY_summer[1:30] <- tempDOY_summer[1:30] + adjustDown
				
				# Smooth out end of summer over last 30 days
				adjustDown2 <- seq(- tempDOY_summer[length(tempDOY_summer)], 0, length.out = 30)
				tempDOY_summer[(length(tempDOY_summer) - 30 + 1): length(tempDOY_summer)] <- tail(tempDOY_summer, 30) + rev(adjustDown2)
				
				tempDOY_ground[summerDOY, i] <- tempDOY_summer
			}
			
				tempDOY <- tempDOY_ground
			}
		
		tempGrid <- tempDOY
		
	
	} else {
		tempDat <- read.csv("parameterization/BAH_temps.csv")

		tempDat$Date <- as.Date(tempDat$Date)
		tempDat$DOY <- strftime(tempDat$Date, format = "%j")
		tempDat$year <- strftime(tempDat$Date, format = "%Y")

		# sample randomly for each timestep, year, DOY, range combination
			# so for each DOY x range combo (365*5 ranges), we need 1/dt * ny samples
			tempGrid <- array(NA, dim = c(dim(timeDat)[1], 5), dimnames = list(timeDat$step, c("Summer",  "Winter",  "Spring",  "Fall", "Calving")))
			for(i in DOY){
				dum <- tapply(tempDat$T2Mean[as.numeric(tempDat$DOY) == DOY[i]], tempDat$Range[as.numeric(tempDat$DOY) == DOY[i]], sample, size = y/dt)
				tempGrid[which(timeDat$DOY == DOY[i]), ] <- cbind(dum$Summer, dum$Winter, dum$Spring, dum$Fall, dum$Calving)
			} # end DOY

		}
	
	return(tempGrid)
}

###############################################################################################
# Calculate parameters based on DOY
###############################################################################################

calcParams <- function(DOY, temp, ppnInhibit = 0, transmission = "base"){
	
	if(transmission == "base") beta <- 10^-6 else if(transmission == "high") beta <- 10^-5 else if(transmission == "low") beta <- 10^-7
	
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
	beta = beta,
	
	# ppnInhibit - the proportion of larvae that go into arrested development (0-1)
	ppnInhibit = ppnInhibit,
	
	# rho4 - development rate of L4 larvae to adults (per day)
	# From Grenfell et al. 1987, development to adults can take 17 days.
	rho4 = 0.06, 
	
	# mu4 - mortality of L4 larvae per day
	mu4 = 0.002,
	
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
	# Estimated from Bathurst data
	k = 0.9940684
	)
	
	return(params)
}

###############################################################################################
# Function to set up initial distribution
###############################################################################################
# Initial spatial distribution
initDist <- function(totPop, x, x.start.sd = 80){
	# Note the shift so that population always starts at x = 0
	shift.x <- round(length(x)/2)
	return(c(totPop / sqrt(2 * pi * (x.start.sd^2)) * exp(- (x - x[shift.x])^2 / (2 * x.start.sd^2)))[c(shift.x:length(x), 1:(shift.x - 1))])
}


###############################################################################################
# Function to simulate caribou dynamics within a season
###############################################################################################

simBou <- function(
	initBou, # Initial conditions for all x for 14 variables
	tempGrid, # space-time grid of temperatures over an annual cycle
	ppnInhibit, # the proportion of larvae undergoing arrested development, between 0 and 1
	transmission = "base", # the transmission coefficient (beta) can be three levels: base, low, or high
	migSpeed = 14, # migration speed of caribou in km/day (can be zero for simualtions without migration) 
	OctP = NULL # the adult parasite burden in October; if supplied then this affects the pregnancy rate of cows in year 1 of the simulation (for use when carrying on from previous sims)
	){
	
	# If only one value of ppnInhibit is supplied, then use that
	# Otherwise use daily estimate in calcParams below
	if(length(ppnInhibit) == 1) ppnInhibit <- rep(ppnInhibit, 365)
	
	# Breeding date, when animals move up a class, is June 7 (DOY = 158)
	breedDOY <- as.numeric(strftime(as.Date("1985-06-07"), "%j"))
	
	# L4 resume development at the start of spring migration
	# Hoar et al. 2012 show resumption in late March
	L4startDOY <- as.numeric(strftime(as.Date("1985-04-20"), "%j"))
	
	# Advection speed for each variable: number of grid spaces moved
	u <- migSpeed * dt / dx
	
	nt <- dim(timeDat)[1]
	
	# # 1) Set up matrices to store solutions every d days
	# d <- 1
	# ntKeep <- floor(dim(timeDat)[1]/d)
	# nKeep <- seq(1, nt, d)
	
	V <- array(NA, dim = c(dim(initBou)[1], length(x), nt), dimnames = list(rownames(initBou), x, paste(timeDat$year, timeDat$time, sep="-")))
	V[, , 1] <- initBou
	V0 <- array(0, dim = c(dim(initBou)[1], length(x)), dimnames = list(rownames(initBou), x))
	
	# 2) Run through each timestep
	# Match x grid spaces to range numbers for indexing gridded temperature data
	rangeMatch <- match(xRange, c("Summer",  "Winter",  "Spring",  "Fall", "Calving"))
	
	for(n in 1:(nt-1)){
		# Set parameters based on DOY
		params.n <- calcParams(
			DOY = timeDat$DOY[n], 
			temp = tempGrid[timeDat$DOY[n], rangeMatch],
			ppnInhibit = ppnInhibit[timeDat$DOY[n]],
			transmission = transmission)

		# Calculate boundary conditions: torus for circular migration
		Vn <- V[, , n]
		Vnp1 <- Vn
		
		# Spatial advection (upstream differencing) for moving stages
		if(u > 0){
			for(j in match(c("calf_mov", "yearling_mov", "adult_mov", "L4_mov", "L4A_mov", "P_mov"), rownames(initBou))){
				Vnp1[j, ] <- Vn[j, c(c((length(x) - u + 1) : length(x)), 1:(length(x) - u))]
			}
		}
		
		#---------------------------------------------------------------------------
		# If breeding date, move caribou up and add calves
		# *** we're going to get weird things happening if the population doesn't mix...
		# ONLY stationary cows have calves, based on average parasite abundance previous October among all
		if(round(timeDat$time[n], 2) == round(breedDOY, 2)){
			# cat("breeding")
			Vn.breed <- Vnp1
			
			if((n - 240/dt) < 0){ # for the first year of the simulation
				if(length(OctP) == 0){
					P_mean <- 0 	# If not supplied, assume zero parasite burden
				} else {
					P_mean <- OctP
			  }
			} else { # for next years
				adult_stat.nonzero <- which(V['adult_stat', , n - 240/dt] > 0)
				adult_mov.nonzero <- which(V['adult_mov', , n - 240/dt] > 0)
				
				# Mean parasite burden across all hosts = #parasites/#hosts
				P_mean <- sum(c(V['P_stat', adult_stat.nonzero, n - 240/dt], V['P_mov', adult_mov.nonzero, n - 240/dt]))/sum(c(V['adult_stat', adult_stat.nonzero, n - 240/dt],  V['adult_mov', adult_mov.nonzero, n - 240/dt]))
			}
			
			newCalves <- numCalves(
				P_mean = P_mean,
				numFemales = (V['adult_stat', , n] + V['adult_mov', , n]) * propFemale,
				pCalf0 = 0.8)
			
			# # All calves start out in stat category
			Vn.breed['calf_stat', ] <- newCalves
			Vn.breed['calf_mov', ] <- rep(0, length(x))
			
			# Yearlings and adults stay in respective categories
			Vn.breed['yearling_stat', ] <- Vnp1['calf_stat', ]
			Vn.breed['yearling_mov', ] <- Vnp1['calf_mov', ]
			Vn.breed['adult_stat', ] <- Vnp1['adult_stat', ] + Vnp1['yearling_stat', ]
			Vn.breed['adult_mov', ] <- Vnp1['adult_mov', ] + Vnp1['yearling_mov', ]
			
			# Replace Vnp1 with updated matrix
			Vnp1 <- Vn.breed
		} # end if breed
		#---------------------------------------------------------------------------
		
		#---------------------------------------------------------------------------
		# If start of spring migration, L4 resume development
		if(round(timeDat$time[n], 2) == round(L4startDOY, 2)){
			# cat("L4 development resuming")
			Vn.start <- Vnp1
			Vn.start['L4_stat', ] <- Vnp1['L4_stat', ] + Vnp1['L4A_stat', ]
			Vn.start['L4A_stat', ] <- rep(0, length(x))
			Vn.start['L4_mov', ] <- Vnp1['L4_mov', ] + Vnp1['L4A_mov', ]
			Vn.start['L4A_mov', ] <- rep(0, length(x))
			
			Vnp1 <- Vn.start
		}
		#---------------------------------------------------------------------------
		
		#---------------------------------------------------------------------------
		# Temporal dynamics (Euler's formula)
		k1 <- partial_t.Bou(y = Vnp1, p = params.n)
		V[, , n + 1]  <- pmax(V0, Vnp1 + dt * k1)
		
		# # Temporal dynamics (4th order Runge Kutta)
		# k1 <- partial_t(Vnp1, params)
		# k2 <- partial_t(Vnp1 + grid$dt / 2 * k1, params)
		# k3 <- partial_t(Vnp1 + grid$dt / 2 * k2, params)
		# k4 <- partial_t(Vnp1 + grid$dt * k3, params)
		# 
		# V[, , n + 1] <- Vnp1 + grid$dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)#, 0) # Added in max to avoid negative values (Mar 14, 2019)
		# if(is.element(n, nKeep) == TRUE) V[, , n] <- Vn
		
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
