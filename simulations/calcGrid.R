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
if(dx == 0) dx <- 2
if(dx == 14) dx <- 2
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
# tempGrid <- matrix(NA, nrow = 365, ncol = length(x))
# for(i in 1:length(x)){
# 	tempGrid[, i] <- predict.temp(DOY, range = xRange[i])
# }
# 
# # Smooth each day
# tempGrid.smooth <- tempGrid
# midX <- round(length(x)/2)
# for(n in 1:length(DOY)){
# 	smooth1 <- lowess(tempGrid[n, c((midX + 1):length(x), 1:midX)], f = 0.3)$y[order(c((midX+1):length(x), 1:midX))]
# 	smooth2 <- lowess(tempGrid[n, ], f = 0.3)$y
# 	dsmooth <- smooth1-smooth2
# 	ind1 <- which(abs(dsmooth[1:round(length(x)/2)]) < 10^-2)
# 	ind3 <- which(abs(dsmooth[(length(x) - round(length(x)/2.5)):(length(x) - round(length(x)/4))]) < 10^-2) + (length(x) - round(length(x)/2.5)) - 1
# 	
# 	tempGrid.smooth[n, ] <- c(smooth1[1:ind1[1]], smooth2[(ind1[1] + 1): ind3[1]], smooth1[(ind3[1] + 1):length(x)]) 
# }
# # 14 km grid
# 
# tempGrid <- tempGrid.smooth
# 
# tempGrid <- matrix(NA, nrow = dim(timeDat)[1], ncol = length(x))
# for(i in 1:length(x)){
# 	tempGrid[, i] <- predict.temp(DOY = timeDat$DOY, range = xRange[i], stoch = TRUE)
# }

# image(tempDOY[timeDat$year == 1], x = timeDat$time[timeDat$year == 1], y = x)

tempGrid <- predict.temp(timeDat = timeDat, stoch = stochTemp)

#------------------------------------------------------------------------------
# Stopping and starting matrices
#------------------------------------------------------------------------------

if(startstop == "simple"){
	# Starting:
	startMat <- matrix(0, nrow = length(DOY), ncol = length(x), dimnames = list(DOY, x))
	startMat[99:151, ] <- 1 # Spring migration
	startMat[168:180, ] <- 1 # Post-calving migration
	startMat[250:345, ] <- 1 # Fall migration
	
	# Stopping:
	stopMat <- matrix(1, nrow = length(DOY), ncol = length(x), dimnames = list(DOY, x))
	stopMat[99:151, ] <- 0 # Spring migration
	stopMat[168:180, ] <- 0 # Post-calving migration
	stopMat[250:345, ] <- 0 # Fall migration
} else if(startstop == "agg"){
	# Starting:
	startMat <- matrix(0, nrow = length(DOY), ncol = length(x), dimnames = list(DOY, x))
	stopMat <- matrix(1, nrow = length(DOY), ncol = length(x), dimnames = list(DOY, x))
	
	#*********************
	#*# 1) Spring migration
	#    a) Make everyone everywhere start moving on April 20
	startMat[111:153, ] <- 1
	stopMat[111:169, ] <- 0 # (Will change stopping for calving below)
	
	#*********************
	# 2) Calving
	#   a) Make everyone stop by 652 km (stop = 1) ramping up from stop = 0 at 553 km for days 1 through 169
	# Need to start at DOY = 138 to catch leading edge at x = 602 km
	xInd <- min(which(x > 602)):min(which(x > 682))
	sdCalving <- 30
	yCurve <- 1/sqrt(2 * pi * sdCalving^2) * exp(-(x[xInd] - 682)^2/(2*sdCalving^2))
	yCurve <- (yCurve - min(yCurve))
	yCurve <- yCurve/max(yCurve)
	for(i in 138:169){
		stopMat[i, xInd] <- yCurve
		stopMat[i, max(xInd):length(x)] <- 1
	}
	
	# plot(x, stopMat[139, ], "l")
	
	#*********************
	# 3) Post-calving migration
	#    a) Change everyone to moving
	startMat[169:180, ] <- 1
	stopMat[169:180, ] <- 0
	
	#*********************
	# 4) Summer grazing dispersal
	#    a) Stop starting at the back. 
	#       Stopping started at 553 km above, so the back is 553 + 11 days * 14 km/day = 707 km
	#       The front is 652 km + 11 days * 14 km/day = 806 km
	# xInd <- which(x > 707 & x < 806)
	startMat[181:230, ] <- 0.1
	
	sdGrazing <- 80
	xInd <- which(x > 650 & x < 965)#233:345
	
	for(i in 181:230){
		minY <- seq(0, 1, length.out = length(181:230))[i - 180]
		yCurve <- 1/sqrt(2 * pi * sdGrazing^2) * exp(-(x[xInd] - x[232])^2/(2*sdGrazing^2))
		yCurve <- (yCurve - min(yCurve)) 
		yCurve <- yCurve/max(yCurve)
		yCurve <- yCurve*(1- minY) + minY
		stopMat[i, xInd] <- yCurve
		stopMat[i, max(xInd):length(x)] <- minY
	}
	
	# plot(x, stopMat[230, ], "l")
	# for(i in 182:230) lines(x, stopMat[i, ], col = grey((i - 181)/40))
	
	#*********************
	# 5) Late summer grazing (all stopped)
	startMat[231:250, ] <- 0
	stopMat[231:250, ] <- 1
	
	#*********************
	# 6) Fall migration (pre-breeding)
	startMat[251:290, ] <- 1
	stopMat[251:305, ] <- 0
	
	#*********************
	# 7) Rut/Breeding
	# Need to bunch up a bit; use same appraoch as for calving
	startMat[291:305, ] <- 0
	
	xInd <- which(x > 1550 & x <= 1700)
	sdRut <- 50
	yCurve <- 1/sqrt(2 * pi * sdRut^2) * exp(-(x[xInd] - 1700)^2/(2*sdRut^2))
	yCurve <- (yCurve - min(yCurve))
	yCurve <- yCurve/max(yCurve)
	for(i in 260:305){
		stopMat[i, xInd] <- yCurve
		stopMat[i, max(xInd):length(x)] <- 1
	}
	
	# plot(x, stopMat[270, ], "l")
	
	
	#*********************
	# 8) Fall migration post-breeding
	# Start migrating at the front to spread out to the original sd = 50
	startMat[306:355, ] <- 1
	stopMat[306:355, ] <- 0
	
	
	#*********************
	# 9) Winter
	# Spread out with stopping rate across winter grounds
	# startMat[336:365, ] <- 0
	# stopMat[336:365, ] <- 1
	
	startMat[336:340, ] <- 1
	stopMat[336:340, ] <- 0
	
	startMat[341:365, ] <- 0.1
	startMat[1:25, ] <- 0.1
	sdWinter <- 80
	
	# FOr summer grazing, the stopping is 1 at x[232] = 646.8 km
	# and declines to zero (to start) by x[345] = 963.2
	# and the host population is evetaully centered at x[361] = 1008 
	# SO if we want the host population centered at x[1] = 0
	# need to adjust xInd by -360
	# length(x) + 232 - 360
	# length(x) + 345 - 360
	# xInd <- 683:796 # Stopped them too soon.  Shift a bit?
	xInd <- which(x > 1990)# 712:811
	counter <- 0
	for(i in c(341:365, 1:25)){
		counter <- counter + 1
		minY <- seq(0, 1, length.out = 50)[counter]
		yCurve <- 1/sqrt(2 * pi * sdWinter^2) * exp(-(x[xInd] - x[min(xInd)])^2/(2*sdWinter^2))
		yCurve <- (yCurve - min(yCurve)) 
		yCurve <- yCurve/max(yCurve)
		yCurve <- yCurve*(1- minY) + minY
		stopMat[i, ] <- minY
		stopMat[i, min(which(x > 1500)) : length(x)] <- 1
		# stopMat[i, c(((length(x)-50):length(x)),c(1:49))] <- yCurve
		stopMat[i, xInd] <- yCurve
		
	}
	
	# plot(x, stopMat[341, ], "l")
	# for(i in 337:365) lines(x, stopMat[i, ], col = grey((i - 335)/40))
	
}
