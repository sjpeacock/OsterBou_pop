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
	stop.host <- rep(0, nx)
	stop.host[adult_mov.nonzero] <- p$stop[adult_mov.nonzero] + p$parasitStop * y['P_mov', adult_mov.nonzero] / y['adult_mov', adult_mov.nonzero]
	
	stop.parasit <- rep(0, nx)
	stop.parasit[adult_mov.nonzero] <- p$stop[adult_mov.nonzero] + p$parasitStop * (y['P_mov', adult_mov.nonzero] / y['adult_mov', adult_mov.nonzero] * (p$k + 1) / p$k + 1)
	
	#-----------------------------------------------------------------------------
	# If there are stationary adults, then starting can be non-zero
	adult_stat.nonzero <- which(y['adult_stat', ] > 0)
	start.host <- rep(0, nx)
	start.host[adult_stat.nonzero] <- p$start[adult_stat.nonzero]
	
	start.parasit <- rep(0, nx)
	start.parasit[adult_stat.nonzero] <- p$start[adult_stat.nonzero] * y['P_stat', adult_stat.nonzero]
	#-----------------------------------------------------------------------------
	
	#*****************************************************************************
	# 1) Number of stationary and moving calves
	dy['calf_stat', ] <- - p$muC * y['calf_stat', ] - p$start * y['calf_stat', ]  + stop.host * y['calf_mov', ]
	
	dy['calf_mov', ] <- - p$muC * y['calf_mov', ] - stop.host * y['calf_mov', ] + p$start * y['calf_stat', ]
	
	#*****************************************************************************
	# 2) Number of yearlings
	dy['yearling_stat', ] <- - p$muY * y['yearling_stat', ] - p$start * y['yearling_stat', ] + stop.host * y['yearling_mov', ]
	
	dy['yearling_mov', ] <- - p$muY * y['yearling_mov', ] - stop.host * y['yearling_mov', ] + p$start * y['yearling_stat', ]
	
	#*****************************************************************************
	# 3) Number of adults
	adultMort_stat <- rep(0, nx) # Need to account for zero parasite-induced mort when N_stat = 0
	adultMort_stat[adult_stat.nonzero] <- y['adult_stat', adult_stat.nonzero] * (p$muA + p$alpha * y['P_stat', adult_stat.nonzero] / y['adult_stat', adult_stat.nonzero])
	
	dy['adult_stat', ] <- - adultMort_stat - p$start * y['adult_stat', ] + stop.host * y['adult_mov', ]
	
	adultMort_mov <- rep(0, nx) # Need to account for zero parasite-induced mort when N_stat = 0
	adultMort_mov[adult_mov.nonzero] <- y['adult_mov', adult_mov.nonzero] * (p$muA + p$alpha * y['P_mov', adult_mov.nonzero] / y['adult_mov', adult_mov.nonzero])
	dy['adult_mov', ] <- - adultMort_mov - stop.host * y['adult_mov', ] + p$start * y['adult_stat', ] 
	
	
	#*****************************************************************************
	# 4) Arrested larvae
	L4Mort_stat <- rep(0, nx) # Need to account for zero parasite-induced mort when N_stat = 0
	L4Mort_stat[adult_stat.nonzero] <- p$muA + p$mu4 + p$alpha * y['P_stat', adult_stat.nonzero] / y['adult_stat', adult_stat.nonzero]
	
	dy['L4A_stat', ] <- p$beta * p$ppnInhibit * y['L3', ] * y['adult_stat', ] - L4Mort_stat * y['L4A_stat', ] - p$start * y['L4A_stat', ] + stop.host * y['L4A_mov', ]
	
	L4Mort_mov <- rep(0, nx) # Need to account for zero parasite-induced mort when N_stat = 0
	L4Mort_mov[adult_mov.nonzero] <- p$muA + p$mu4 + p$alpha * y['P_mov', adult_mov.nonzero] / y['adult_mov', adult_mov.nonzero]
	
	dy['L4A_mov', ] <- p$beta * p$ppnInhibit * y['L3', ] * y['adult_mov', ] - y['L4A_mov', ] * L4Mort_mov + p$start * y['L4A_stat', ] - stop.host * y['L4A_mov', ]
	
	#*****************************************************************************
	# 5) Developing larvae
	dy['L4_stat', ] <- p$beta * (1 - p$ppnInhibit) * y['L3', ] * y['adult_stat', ] - y['L4_stat', ] * L4Mort_stat - p$rho4 * y['L4_stat', ] - p$start * y['L4_stat', ] + stop.host * y['L4_mov', ]
	
	dy['L4_mov', ] <- p$beta * (1 - p$ppnInhibit) * y['L3', ] * y['adult_mov', ] -  y['L4_mov', ] * L4Mort_mov - p$rho4 * y['L4_mov', ] + p$start * y['L4_stat', ] - stop.host * y['L4_mov', ]
	
	#*****************************************************************************
	# 6) Adult parasites
	PMort_stat <- rep(0, nx)
	PMort_stat[adult_stat.nonzero] <- y['P_stat', adult_stat.nonzero] / y['adult_stat', adult_stat.nonzero] * (p$k + 1)/p$k + 1
	
	dy['P_stat', ] <- p$rho4 * y['L4_stat', ] - (p$muA + p$muP) * y['P_stat', ] - (p$alpha + p$nuP) * y['P_stat', ] * PMort_stat - start.parasit + y['P_mov', ] * stop.parasit
	
	PMort_mov <- rep(0, nx)
	PMort_mov[adult_mov.nonzero] <- y['P_mov', adult_mov.nonzero] / y['adult_mov', adult_mov.nonzero] * (p$k + 1)/p$k + 1
	
	dy['P_mov', ] <- p$rho4 * y['L4_mov', ] - (p$muA + p$muP) * y['P_mov', ] - (p$alpha + p$nuP) * y['P_mov', ] * PMort_mov + start.parasit - y['P_mov', ] * stop.parasit
	
	#*****************************************************************************
	# 7) Larvae
	dy['L0', ] <- p$lambda * (y['adult_stat', ] * y['P_stat', ]^(1 + p$phi) + y['adult_mov', ] * y['P_mov', ]^(1 + p$phi)) - (p$mu0 + p$rho0) * y['L0', ]
	
	dy['L3', ] <- p$rho0 * y['L0', ] - p$mu3 * y['L3', ] - p$beta * y['L3', ]* (y['adult_stat', ] + y['adult_mov', ])
	
	#*****************************************************************************
	if(sum(is.na(dy)) > 0) stop("\n\nNAs in derivative function.\n\n")
	unique(which(is.na(dy) == TRUE, arr.ind = TRUE)[, 1])
	return(dy)
	
} #end function