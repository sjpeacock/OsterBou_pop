###############################################################################
# Looking at annual temperature cycles over the range of the Bathurst caribou
# As compiled by Don Russell for CARMA from MERRA's climate data:
#
# Russell, D.E., Whitfield, P.H., Cai, J., Gunn, A., White, R.G., Poole, K., 
# 2013. CARMA’s MERRA-based caribou range climate database. Rangifer 33, 145. 
# https://doi.org/10.7557/2.33.2.2535
#
library(gplots)
library(dclone)
###############################################################################
seasonCols <- c(
	"Spring" = "#9FCE62", 
	"Calving" = "#FFE548", 
	"postCalving" = "#FFA251",
	"Summer" = "#DC4A45",
	"Fall" = "#8E4D1D",
	"Winter" = "#618CCF")
seasonCols2 <- paste(seasonCols, "60", sep = "")
names(seasonCols2) <- names(seasonCols)

tempDat <- read.csv("parameterization/BAH_temps.csv")

tempDat$Date <- as.Date(tempDat$Date)
tempDat$DOY <- strftime(tempDat$Date, format = "%j")
tempDat$year <- strftime(tempDat$Date, format = "%Y")

season <- unique(tempDat$Range)

par(mfrow = c(3,1))

for(y in c(1980, 1995, 2010)){
ind <- which(tempDat$year == y)

plot(tempDat$Date[ind][tempDat$Range[ind] == season[1]], tempDat$T2Mean[ind][tempDat$Range[ind] == season[1]], "l", ylim = range(tempDat$T2Mean, an.rm=TRUE), col = seasonCols[season[1]], las = 1)
for(j in 2:5){
	lines(tempDat$Date[ind][tempDat$Range[ind] == season[j]], tempDat$T2Mean[ind][tempDat$Range[ind] == season[j]], col = seasonCols[season[j]], las = 1)
}
mtext(side = 3, adj = 0, line= 0.5, y)

}
mtext(side = 2, outer=TRUE, expression(paste("Temperature (", degree, "C)")), line = 2)
# All the seasons do pretty much look the same!!


#------------------------------------------------------------------------------
# How much interannual variability?
nY <- length(unique(tempDat$year))

for(j in 1:5){
	plot(c(1, 366), c(-50, 30), "n", bty="l", ylab = "Temperature (*C)")
	mtext(side = 3, season[j], col = seasonCols[season[j]])
	for(i in 1:nY){
		Y <- tempDat$T2Mean[tempDat$year == unique(tempDat$year)[i] & tempDat$Range == season[j]]
		lines(1:length(Y), Y, col= seasonCols2[season[j]])
		
	}
}

	
#------------------------------------------------------------------------------
# Fit sine curve to basic data

# 1) Assuming different parameters for summer, winter, spring, etc. ranges
tempDat.range <- list(
	n = dim(tempDat)[1],
	DOY = tempDat$DOY,
	range = as.numeric(factor(tempDat$Range, levels = season)),
  temp = tempDat$T2Mean
)

tempMod.range <- function(){
	# Define parameters
	for(i in 1:5){
		ck[i] ~ dnorm(10, pow(10, -2))
	  dk[i] ~ dnorm(10, pow(10, -2))
	  t0[i] ~ dunif(0, 365)
	}
	
	sig ~ dlnorm(-3, pow(2, -2))
	
	for(z in 1:n){
		tempPred[z] <- ck[range[z]] + dk[range[z]] * cos((DOY[z] - t0[range[z]])* 2 * 3.141593 / 365)
		temp[z] ~ dnorm(tempPred[z], pow(sig, -2))
	}
}

cl <- makeCluster(3)

t.start<-proc.time()[3]

fit.range <- jags.parfit(
	cl, 
	data = tempDat.range, 
	params = c("ck", "dk", "t0", "sig"), 
	model = tempMod.range, 
	n.iter = 1000,
	n.adapt = 1000,
	n.chains = 3) 

stopCluster(cl)
procTime[m, K] <-  # track process time
cat(paste("Process time ", round((proc.time()[3] - t.start)/(60), 1), "minutes"))

out1 <- summary(fit.range)
p <- out1[[1]][,1]

DOY.date <- as.Date(paste("2000", c(1:365), sep="-"), format = "%Y-%j")
par(mfrow=c(1,1), mar=c(4,4,2,1), oma = rep(0,4))
plot(DOY.date, rep(1, length(DOY.date)), ylim = c(-40, 20), "n", bty="l", las= 1, xlab = "", ylab = "Temp")
for(i in 1:5){
	lines(DOY.date, p[paste('ck[', i, ']', sep="")] + p[paste('dk[', i, ']', sep="")] * cos((DOY - p[paste('t0[', i, ']', sep="")]) * 2 * pi / 365), col = seasonCols[i])
	maxT <- p[paste('ck[', i, ']', sep="")] + p[paste('dk[', i, ']', sep="")]
	t0 <- as.Date(paste("2000", p[paste('t0[', i, ']', sep="")], sep="-"), format = "%Y-%j") 

segments(x0 = as.Date("1999-12-15"), x1 = t0, y0 = maxT, y1 = maxT, col = seasonCols[i], lty = 3)
segments(x0 = t0, x1 = t0, y0 = maxT, y1 = -50, col = seasonCols[i], lty = 3)
}

legend("topright", lwd = 1, col = seasonCols, season, bty = "n")

