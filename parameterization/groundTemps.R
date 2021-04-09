###############################################################################
# Compare air temperatures from CARMA MERRA:
#
# Russell, D.E., Whitfield, P.H., Cai, J., Gunn, A., White, R.G., Poole, K., 
# 2013. CARMA’s MERRA-based caribou range climate database. Rangifer 33, 145. 
# https://doi.org/10.7557/2.33.2.2535
#
# to profiles from 
#
# Convey, P., S. J. Coulson, M. R. Worland, and A. Sjöblom. “The Importance of Understanding Annual and Shorter-Term Temperature Patterns and Variation in the Surface Levels of Polar Soils for Terrestrial Biota.” Polar Biology 41, no. 8 (2018): 1587–1605. https://doi.org/10.1007/s00300-018-2299-0.
#
# and see if we can estimate ground temps for BAH range
###############################################################################

tempDat <- read.csv("parameterization/BAH_temps.csv")

tempDat$Date <- as.Date(tempDat$Date)
tempDat$DOY <- as.numeric(strftime(tempDat$Date, format = "%j"))
tempDat$year <- strftime(tempDat$Date, format = "%Y")

# Consider summer range only since that's where the majority of transmission
# occurs

tempDat <- subset(tempDat, tempDat$Range == "Summer")

min(tempDat$T2M_MIN) # - 54.5
max(tempDat$T2M_MAX) # 27.5
min(tempDat$T2Mean) # - 51.6

fun2 <- function(x){
	length(which(x > 0))
}

posCount <- function(dat){
mean(tapply(dat$T2Mean, dat$year, fun2)) #119 days with mean positive temp
}

y <- 1980
z <- as.numeric(tempDat$T2Mean[tempDat$year == y & tempDat$DOY > 100 & tempDat$DOY < 280] < 0)
x <- tempDat$Date[tempDat$year == y & tempDat$DOY > 100 & tempDat$DOY < 280]
plot(x, z, "l")

# Number of freeze-thaw events per year
fun3 <- function(x){
	z <- as.numeric(x < 0)
	return(sum(diff(z) == -1))
}

FTcount <- function(dat){
	mean(tapply(dat$T2Mean, dat$year, fun3)) 
}


plot(tempDat$Date, tempDat$T2Mean, "l", ylab = "Mean air temperature (*C)", xlab = "")
abline(h = 0)

# Summary statistics
tempSumm <- data.frame(
	season = c("MAM", "JJA", "SON", "DJF", "Annual"), 
	TMean = rep(NA, 5),
	TMean_Min = rep(NA, 5),
	TMean_Max = rep(NA, 5),
	TMin_Min = rep(NA, 5),
	TMax_Max = rep(NA, 5),
	FTevents = rep(NA, 5),
	DaysG0 = rep(NA, 5))

for(i in 1:5){
	if(i == 1){
		tempDat.i <- subset(tempDat, tempDat$DOY >= 61 & tempDat$DOY <= 152)
	} else if(i == 2){
		tempDat.i <- subset(tempDat, tempDat$DOY >= 153 & tempDat$DOY <= 244)
	}else if(i == 3){
		tempDat.i <- subset(tempDat, tempDat$DOY >= 245 & tempDat$DOY <= 335)
	} else if(i == 4){
		tempDat.i <- subset(tempDat, tempDat$DOY >= 336 | tempDat$DOY <= 60)
	} else{
		tempDat.i <- tempDat
	}
	
	tempSumm$TMean[i] <- round(mean(tempDat.i$T2Mean), 2)
	tempSumm$TMean_Min[i] <- round(min(tempDat.i$T2Mean), 2)
	tempSumm$TMean_Max[i] <- round(max(tempDat.i$T2Mean), 2)
	tempSumm$TMin_Min[i] <- round(min(tempDat.i$T2M_MIN), 2)
	tempSumm$TMax_Max[i] <- round(max(tempDat.i$T2M_MAX), 2)
	
	tempSumm$FTevents[i] <- round(FTcount(tempDat.i), 1)
	tempSumm$DaysG0[i] <- round(posCount(tempDat.i), 1)

}

write.csv(tempSumm, file = "parameterization/summerRangeTempSummary.csv")

###############################################################################
# Cumulative degree days
###############################################################################
ny <- length(unique(tempDat$year))

dd <- matrix(NA, nrow = 365, ncol = ny) 

for(y in 1:ny){
	tempDat.y <- subset(tempDat, tempDat$year == unique(tempDat$year)[y])
	
	# Set negative days to zero
	tempDat.y$T2Mean[tempDat.y$T2Mean < 0] <- 0
	
	# Remove leap year day
	if(length(tempDat.y$T2Mean) > 365) tempDat.y <- tempDat.y[c(1:59, 61:366), ]
	
	dd[, y] <- cumsum(tempDat.y$T2Mean)
	
}

# thaw dates
thaw <- c(round(mean(apply(dd, 2, function(x) which(x > max(x))[1]))), range(apply(dd, 2, function(x) which(x > 0)[1])))

freeze <- c(round(mean(apply(dd, 2, function(x) which(x == max(x))[1]))), range(apply(dd, 2, function(x) which(x == max(x))[1])))

plot(1:365, dd[, 1], "l", col = "#00000030", ylim = c(0, max(dd)), las = 1, ylab = "Cumulative degree days", xlab = "DOY")
for(i in 2:ny) lines(1:365, dd[, i], col = "#00000030")
lines(1:365, apply(dd, 1, mean), lwd = 2)
abline(v = thaw, col = 2, lwd = c(2,0.8, 0.8))
abline(v = freeze, col = 4, lwd = c(2,0.8, 0.8))
u <- par('usr')

text(thaw, rep(1300, 3), srt = 90, xpd = NA, col = 2, strftime(as.Date(thaw, origin = "1970-01-01"), format = "%b-%d"), pos = 2)

text(freeze, rep(200, 3), srt = 90, xpd = NA, col = 4, strftime(as.Date(freeze, origin = "1970-01-01"), format = "%b-%d"), pos = 2)

###############################################################################
# High Arctic Shrub data from Convey et al. (2018)
###############################################################################

grDat <- read.csv("parameterization/Convey2018_ArcticB.csv")

# Calculate mean daily temperatures
TGMean <- tapply(grDat$Temp, grDat$Date, mean)
TGMin <- tapply(grDat$Temp, grDat$Date, min)
TGMax <- tapply(grDat$Temp, grDat$Date, max)

GDate <- as.Date(names(TGMean))

plot(GDate, TGMean, "l", ylim = range(TGMean, TGMin, TGMax), las = 1, ylab = "Mean daily ground temperature (*C)", xlab = "")
abline(h = 0)
lines(GDate, TGMin, col = 4)
lines(GDate, TGMax, col = 2)

qDate <- function(n){
	x <- locator(n)
	print(as.Date(x$x, origin = "1970-01-01"))
}

# Freeze up at the end of Sept