moveDat <- read.csv("parameterization/Gunn2001_bouMovement.csv")



range(moveDat$DOY)

seasonCols <- c(
	"Spring" = "#9FCE62", 
	"Calving" = "#FFE548", 
	"postCalving" = "#FFA251",
	"Summer" = "#DC4A45",
	"Fall" = "#8E4D1D",
	"Winter" = "#618CCF")

seasonStart <- c(
	"Spring" = 111, # April 20
	"Calving" = 154, # June 2
	"postCalving" = 169,
	"Summer" = 181,
	"Fall" = 251,
	"Winter" = 336
)

# Average movement rates for each period
movRate <- numeric(6)
for(i in 1:6){
	if(i == 6){ # if winter
		ind <- which(moveDat$DOY < seasonStart[i] | moveDat$DOY > seasonStart[6])
	} else {
		ind <- which(moveDat$DOY >= seasonStart[i] & moveDat$DOY < seasonStart[i+1])
	}
	movRate[i] <- mean(moveDat$movementRate[ind])
}

plot(moveDat$DOY, moveDat$movementRate, "n", xlab = "DOY", ylab = "Movement rate (km/day", las = 1, xlim = c(1, 365), xaxs="i")
# abline(v = seasonStart, col = seasonCols, lwd = 4)
for(i in 1:5){
	polygon(x = seasonStart[c(i, i+1, i+1, i)], y = c(0, 0, 25, 25), col= seasonCols[i], border = NA)
	segments(x0 = seasonStart[i], x1 = seasonStart[i+1] - 1, y0 = movRate[i], y1 = movRate[i], col = "white", lwd = 4)
}
polygon(x = c(seasonStart[6], 365, 365, seasonStart[6]), y = c(0, 0, 25, 25), col= seasonCols[6], border = NA)
polygon(x = c(1, seasonStart[1], seasonStart[1], 1), y = c(0, 0, 25, 25), col= seasonCols[6], border = NA)
segments(x0 = c(seasonStart[6], 1), x1 = c(365, seasonStart[1] - 1), y0 = rep(movRate[6], 2), y1 = rep(movRate[6], 2), col = "white", lwd = 4)

points(moveDat$DOY, moveDat$movementRate, pch = 19, cex = 0.5)

L <- loess(movementRate ~ DOY, data = moveDat, span = 60/365)
Lpred <- predict(L1, newdata = data.frame(DOY = c(1:365)))
lines(1:365, Lpred, col = grey(0.8), lwd = 3)

write.csv(data.frame(DOY = c(1:365), movementRate = Lpred), file = "parameterization/dailyMovementRates.csv")


# Breakpoints?
library(BreakPoints)
N_break_point(moveDat$movementRate[order(moveDat$DOY)], n_max = 8, n_period=14)
