predict.MTE <- function(
	params = c(a = 0.068, E = 0.884, Eh = NA, Th = NA, El = NA, Tl = NA, z = 1),
	temp = 15){
	# Arrhenius	
	if(is.na(params['Eh']) == TRUE & is.na(params['El']) == TRUE){
		return(params['a'] * exp(-params['E']/(8.62*10^-5) * (1/(temp + 273.15) - 1/(15 + 273.15))))
		
	}else{
		# SS upper bound
		if(is.na(params['El']) == TRUE){
			return(params['a'] * exp(-params['E']/(8.62*10^-5) * (1/(temp + 273.15) - 1/(15 + 273.15))) * (1 + exp(params['Eh'] / (8.62*10^-5)*(-1/(temp+273.15)+1/(params['Th']+273.15))))^params['z'])
		}
		# SS lower bound
		if(is.na(params['Eh']) == TRUE){
			return(params['a'] * exp(-params['E']/(8.62*10^-5) * (1/(temp + 273.15) - 1/(15 + 273.15))) * (1 + exp(params['El'] / (8.62*10^-5)*(1/(temp+273.15)-1/(params['Tl']+273.15))))^params['z'])
		}
		# SS upper AND lower bound
		if(sum(is.na(params)) == 0){
			return(params['a'] * exp(-params['E']/(8.62*10^-5) * (1/(temp + 273.15) - 1/(15 + 273.15))) * (1 + exp(params['El'] / (8.62*10^-5)*(1/(temp+273.15)-1/(params['Tl']+273.15))) + exp(params['Eh'] / (8.62*10^-5)*(-1/(temp+273.15)+1/(params['Th']+273.15))))^params['z'])
		}
	}
} # end function

paramsMTE <- list(
	mu0 = c(a = 0.068,E = 0.884, El = 3.25, Tl = 0, Eh = NA, Th = NA, z = 1),
	rho0 = c(a = 0.032, E = 0.686, Eh = 7.957, Th = 30.568, El =NA, Tl = NA, z = -1),
	mu1 = c(a = 0.0211, E = 0.208, Eh = 3.5543, Th = 27.6, El = 2, Tl = 0, z = 1)
)


T.all <- seq(-40, 40, 0.1)

plot(T.all, predict.MTE(params = paramsMTE$mu1, temp = T.all), "l", ylim= c(0,1))
lines(T.all, predict.MTE(params = paramsMTE$mu0, temp = T.all), lty = 2)
lines(T.all, predict.MTE(params = paramsMTE$rho0, temp = T.all), col = 2)


ppnSurv <- function(temp, paramsMTE, infectiousTime = 7){
	mu0 <- predict.MTE(params = paramsMTE$mu0, temp = temp)
	rho0 <- predict.MTE(params = paramsMTE$rho0, temp = temp)
	mu1 <- predict.MTE(params = paramsMTE$mu1, temp = temp)
	
	return(exp(-mu0 * (1/rho0) - mu1 * infectiousTime))
}


plot(T.all, ppnSurv(temp = T.all, paramsMTE, infectiousTime = 7), "l")


###############################################################################
# Bias from using mean temp
###############################################################################

# Read in CARMA MERRA dataset
tempDat <- read.csv("parameterization/BAH_temps.csv")
tempDat$Date <- as.Date(tempDat$Date)
tempDat$DOY <- strftime(tempDat$Date, format = "%j")
tempDat$year <- strftime(tempDat$Date, format = "%Y")

# Calculate bias over two weeks in July (peak transmission?)
DOYrange <- as.numeric(strftime(c(as.Date("2000-07-01"), as.Date("2000-07-14")), format = "%j"))

tempDat2 <- subset(tempDat, (tempDat$DOY >= DOYrange[1] & tempDat$DOY <= DOYrange[2]) & tempDat$Range == "Summer")

# Try 2 different scales: (1) hourly, (2) daily

makeHourly <- function(dailyMean, dailyMax, dailyMin){
	
	# dailyMean <- tempDat2$T2Mean[tempDat2$year == 1980]
	# dailyMax <- tempDat2$T2M_MAX[tempDat2$year == 1980]
	# dailyMin <- tempDat2$T2M_MIN[tempDat2$year == 1980]

	n <- length(dailyMean)
	out <- data.frame(day = rep(1:n, each = 24), hour = rep(1:24,n), temp = rep(NA, n*24))
	
	for(i in 1:n){
		
		upper <- dailyMean[i] + (dailyMax[i] - dailyMean[i]) * cos((c(1:24) - 12)* 2 * 3.141593 / 24)
		lower <- dailyMean[i] + (dailyMean[i] - dailyMin[i]) * cos((c(1:24) - 12)* 2 * 3.141593 / 24)
		
		out$temp[out$day == i] <- c(lower[1:6], upper[7:18], lower[19:24])
	}
	
	return(out)
}

# Sample year 1980
y <- 1980
tempDat.y <- subset(tempDat2, tempDat2$year == y)
out <- makeHourly(dailyMean = tempDat.y$T2Mean, dailyMax = tempDat.y$T2M_MAX, dailyMin = tempDat.y$T2M_MIN)

plot(as.numeric(tempDat.y$DOY) + 0.5, tempDat.y$T2Mean, "o", pch = 19, cex = 0.5, ylim = range(c(tempDat.y$T2M_MIN, tempDat.y$T2M_MAX)), xlim = c(DOYrange[1], DOYrange[2] + 1), xlab = "DOY", ylab = expression(paste("Temperature (", degree, "C)")), las = 1)

lines(as.numeric(tempDat.y$DOY) + 0.5, tempDat.y$T2M_MIN, lty = 2)
lines(as.numeric(tempDat.y$DOY) + 0.5, tempDat.y$T2M_MAX, lty = 2)
points(out$day + out$hour/24 + DOYrange[1] - 1, out$temp, "o", col = 2, pch = 19, cex = 0.5)

#------------------------------------------------------------------------------

parOverTime <- list(
	mu0 = list(
		daily = predict.MTE(params = paramsMTE$mu0, temp = tempDat.y$T2Mean),
	  hourly = predict.MTE(params = paramsMTE$mu0, temp = out$temp)
		),
	mu1 = list(
		daily = predict.MTE(params = paramsMTE$mu1, temp = tempDat.y$T2Mean),
		hourly = predict.MTE(params = paramsMTE$mu1, temp = out$temp)
	),
	rho0 = list(
	  daily = predict.MTE(params = paramsMTE$rho0, temp = tempDat.y$T2Mean),
	  hourly = predict.MTE(params = paramsMTE$rho0, temp = out$temp)
	)
)

parLab <- c(expression(paste("Pre-infective mortality (", mu[0], ")")), expression(paste("Infective mortality (", mu[1], ")")), expression(paste("Pre-infective development (", rho[0], ")")))
						
for(i in 1:3){
	
	plot(as.numeric(tempDat.y$DOY) + 0.5, parOverTime[[i]]$daily, "o", pch = 19, cex = 0.5, ylim = range(parOverTime[[i]]), xlim = c(DOYrange[1], DOYrange[2] + 1), xlab = "DOY", ylab = "Rate per day", las = 1)
	
	mtext(side = 3, line = 0, parLab[i])
	
	lines(out$day + out$hour/24 + DOYrange[1] - 1, parOverTime[[i]]$hourly, lty = 3, col = 2)
	
	points(as.numeric(tempDat.y$DOY) + 0.5, tapply(parOverTime[[i]]$hourly, out$day, sum)/24, "o", col = 2, pch = 19, cex = 0.5)
}

#------------------------------------------------------------------------------
# Percent developed or survived

dLarvae <- function(y, par){
	dy1 <- - (par[1] + par[3]) * y[1]
	dy2 <- par[3] * y[1] - par[2] * y[2]
	return(c(dy1, dy2))
}

larv <- list(
	daily = matrix(NA, nrow = dim(tempDat.y)[1] + 1, ncol = 2),
  hourly = matrix(NA, nrow = dim(out)[1] + 1, ncol = 2)	
)

# Daily
larv$daily[1, ] <- c(1, 0)

for(i in 2:(dim(tempDat.y)[1]+1)){
	larv$daily[i, ] <- larv$daily[i - 1, ] + dLarvae(y = larv$daily[i - 1, ], par = c(parOverTime$mu0$daily[i-1], parOverTime$mu1$daily[i-1], parOverTime$rho0$daily[i-1])) * 1
}

# Hourly
larv$hourly[1, ] <- c(1, 0)

for(i in 2:(dim(out)[1] + 1)){
	larv$hourly[i, ] <- larv$hourly[i - 1, ] + dLarvae(y = larv$hourly[i - 1, ], par = c(parOverTime$mu0$hourly[i-1], parOverTime$mu1$hourly[i-1], parOverTime$rho0$hourly[i-1])) * (1/24)
}

plot(c(0:(dim(tempDat.y)[1]))+183, larv$daily[, 1], "l", ylim = c(0,1), xlab="Day", ylab = "% survival", las= 1)
lines(c(0:(dim(tempDat.y)[1]))+183, larv$daily[, 2], lty = 2)
lines(c(0:(dim(out)[1]))/24+183, larv$hourly[, 1], col = 2)
lines(c(0:(dim(out)[1]))/24+183, larv$hourly[, 2], lty = 2, col = 2)
legend("topright", lty = c(1,2), c("pre-infective", "infective"), bty = "n")


###############################################################################
# Throhgout the year

# For each DOY, calculate daily mean change in pre-infective and infective larvae
# using average temp and hourly temps.  Calculate percent bias in the change in larvae.

# Repeat for each year we have data

# Plot range for each DOY.

n.y <- length(unique(tempDat$year))

dL_relativeBias <- array(NA, dim = c(365, n.y, 2))

for(i in 1:365){
	for(j  in 1:n.y){
		
		dat.ij <- subset(tempDat, as.numeric(tempDat$DOY) == i & as.numeric(tempDat$year) == unique(as.numeric(tempDat$year))[j] & tempDat$Range == "Fall")
		
		hour.ij <- makeHourly(dailyMean = dat.ij$T2Mean, dailyMax = dat.ij$T2M_MAX,dailyMin = dat.ij$T2M_MIN) 
		
		par.ij <- list(
			daily = c(predict.MTE(params = paramsMTE$mu0, temp = dat.ij$T2Mean),predict.MTE(params = paramsMTE$mu1, temp = dat.ij$T2Mean),predict.MTE(params = paramsMTE$rho0, temp = dat.ij$T2Mean)),
			hourly = cbind(predict.MTE(params = paramsMTE$mu0, temp = hour.ij$temp), predict.MTE(params = paramsMTE$mu1, temp = hour.ij$temp),	hourly = predict.MTE(params = paramsMTE$rho0, temp = hour.ij$temp))
		)
		
		dailyChange <- 1 + pmax(-1, dLarvae(c(1, 1), par = par.ij$daily))
		
		hourlyChange <- matrix(NA, nrow = 25, ncol = 2)
		hourlyChange[1, ] <- 1
		for(t in 1:24){
			hourlyChange[t + 1, ] <- hourlyChange[t, ] + pmax(- hourlyChange[t, ], dLarvae(hourlyChange[t, ], par.ij$hourly[t, ])*(1/24))
		}
		
		if(hourlyChange[25, 1] == 0){
			dL_relativeBias[i, j, 1] <- (dailyChange[1] - hourlyChange[25, 1])
		} else {
			dL_relativeBias[i, j, 1] <- (dailyChange[1] - hourlyChange[25, 1]) / hourlyChange[25, 1]
		}
		
		if(hourlyChange[25, 2] == 0){
			dL_relativeBias[i, j, 2] <- (dailyChange[2] - hourlyChange[25, 2])
		} else {
			dL_relativeBias[i, j, 2] <- (dailyChange[2] - hourlyChange[25, 2]) / hourlyChange[25, 2]
		}
	
		
	}
}

plot(1:365, apply(dL_relativeBias[,,1], 1, quantile, 0.5, na.rm = TRUE), "l", ylim = c(-1, 3), xlab = "DOY", ylab = "% relative bias", las= 1, main = "1980-2019")
lines(1:365, apply(dL_relativeBias[,,1], 1, quantile, 0.025, na.rm = TRUE), lty = 2)
lines(1:365, apply(dL_relativeBias[,,1], 1, quantile, 0.975, na.rm = TRUE), lty = 2)

##

plot(1:366, tempDat$T2Mean[tempDat$year == unique(tempDat$year)[1] & tempDat$Range == "Fall"], "l", xlim = c(80, 120))
polygon(x = c(c(1:366), rev(c(1:366))), y = c(tempDat$T2M_MIN[tempDat$year == unique(tempDat$year)[1] & tempDat$Range == "Fall"], rev(tempDat$T2M_MAX[tempDat$year == unique(tempDat$year)[1] & tempDat$Range == "Fall"])), col = "#00000060", border = NA)
abline(v = which(dL_relativeBias[, 1, 1] > 1), col = 2)

plot(1:365, dL_relativeBias[, 1, 1], "l", ylim = c(-1.1, 3), xlim = c(80, 120))
# for(j in 2:3) lines(1:365, dL_relativeBias[, j, 1], col = j)
abline(v = which(dL_relativeBias[, 1, 1] > 1), col = 2)
lines(1:365, dL_relativeBias[, 1, 2], lty = 2)


plot(1:365, dL_relativeBias[, 1, 1], "l", ylim = c(-1.1, 3))

# Example of bias
i <- 87


par(mfrow = c(2,3))

plot(hour.ij$hour, hour.ij$temp, "l", xlab = "Hour", ylab = "Temp", col = 2,ylim = c(dat.ij$T2M_MIN, dat.ij$T2M_MAX))
abline(h = c(dat.ij$T2Mean, dat.ij$T2M_MAX, dat.ij$T2M_MIN), lty = c(1,2,2))
mtext(side = 3, line = 0, "Temperature", adj = 0)

for(i in 1:3){
	plot(hour.ij$hour, par.ij$hourly[, i], "l", col = 2, ylab = "Rate per day", xlab = "Hour", las = 1)
	mtext(side = 3, line = 0, parLab[i], adj = 0) 
	abline(h = par.ij$daily[i])
}

for(i in 1:2){
	plot(c(0, hour.ij$hour), hourlyChange[, i], "l", col = 2, las =1, ylan = "Survival", xlab = "Hour")
	lines(c(0, 24), c(1, dailyChange[i]))
	mtext(side = 3, line = 0, c("Proportion of pre-infective", "Proportion of infective")[i], adj = 0)
	
	abline(h = dailyChange[i], lty = 3)
	text(3,  dailyChange[i], round(dailyChange[i], 2), pos = 3)
	
	abline(h = hourlyChange[25, i], lty = 3, col = 2)
	text(3,  hourlyChange[25, i], round(hourlyChange[25, i], 2), col = 2, pos = 3)
	
	text(20, 1, paste(round((dailyChange[i] - hourlyChange[25, i])/hourlyChange[25, i], 2)*100, "% bias", sep = ""), pos = 1)
	
}


upper <- dat.ij$T2Mean + (dat.ij$T2M_MAX - dat.ij$T2Mean) * cos((c(1:24) - 12)* 2 * 3.141593 / 24)
lower <- dat.ij$T2Mean + (dat.ij$T2Mean - dat.ij$T2M_MIN) * cos((c(1:24) - 12)* 2 * 3.141593 / 24)

plot(hour.ij$hour, upper, "l", ylim = range(c(upper, lower)))
abline(h = c(dat.ij$T2Mean, dat.ij$T2M_MAX, dat.ij$T2M_MIN), lty = c(1,2,2))
lines(hour.ij$hour, c(lower[1:6], upper[7:18], lower[19:24]), col = 2)

plot(hour.ij$hour, lower, "l", ylim = range(c(upper, lower)))
abline(h = c(dat.ij$T2Mean, dat.ij$T2M_MAX, dat.ij$T2M_MIN), lty = c(1,2,2))

