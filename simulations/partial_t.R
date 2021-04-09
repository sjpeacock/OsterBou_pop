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
