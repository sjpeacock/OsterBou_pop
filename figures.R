library(gplots)
library(ggsci)
library(pals)
colPal <- pal_locuszoom()#scico(n = 7, palette = "berlin")
cols <- colPal(7)
pal.safe(cols)

library(pals)
colPal <- pal_locuszoom()#scico(n = 7, palette = "berlin")
cols <- colPal(7)
colPal2 <- pal_locuszoom(alpha = 0.5)#scico(n = 7, palette = "berlin")
cols2 <- colPal2(7)

names(cols) <- c("Calving", "Summer", "Spring", "Fall", "Winter", "purple", "grey")
names(cols2) <- c("Calving", "Summer", "Spring", "Fall", "Winter", "purple", "grey")
###############################################################################
# Conceptual illustration of MTE relationships
###############################################################################

Kelvin <- function(temp){ return(temp + 273.15)}

E <- 1
Eh <- 3.25
El <- 3.25
Tl <- Kelvin(0)
Th <- Kelvin(30)
T0 <- Kelvin(15)
a <- 0.1
k <- 8.62 * 10^-5

Temp <- seq(-10, 40, 0.1)
nT <- length(Temp)

param <- cbind(
	
	constant = rep(a, nT),
	
	A = a * exp(-E / k * (1 / Kelvin(Temp) - 1 / T0)),
	
	SSUL = a * exp(-E/k*(1/Kelvin(Temp) - 1/T0))  * (1 + exp(El/k * (1/Kelvin(Temp) - 1/Tl)) + exp(Eh/k*(-1/Kelvin(Temp) + 1 / Th))),
	
	SSUL.rho = a * exp(-E/k*(1/Kelvin(Temp) - 1/T0))  * (1 + exp(El/k * (1/Kelvin(Temp) - 1/Tl)) + exp(Eh/k*(-1/Kelvin(Temp) + 1 / Th)))^(-1)

	)

quartz(width = 6.3, height = 3, pointsize = 10)
layout(matrix(c(1,1,2), nrow = 1))
 par(mar = c(4,3,1,1))
plot(Temp, param[, 1], "l", ylim = c(0, 1), col= cols[1], las = 1, bty="l", yaxt="n", xaxt="n", xlab = "", ylab = "", lwd = 1.5)
u <- par('usr')
arrows(x0 = u[1], x1 = u[2], y0 = u[3], y1 = u[3], length = 0.1, xpd = NA)
arrows(x0 = u[1], x1 = u[1], y0 = u[3], y1 = u[4], length = 0.1, xpd = NA)
lines(Temp, param[, 2], col = cols[2], lwd = 1.5)
lines(Temp, param[, 3], col = cols[3], lty = 2, lwd = 1.5)
lines(Temp, param[, 4], col = cols[4], lty = 2, lwd = 1.5)
abline(v = c(Tl, Th)-273.15, lty = 2)
axis(side = 1, at = c(Tl, Th)-273.15, labels = c(expression(paste(italic(T[L]))), expression(paste(italic(T[H])))))
mtext(side = 1, line = 3, "Temperature")
mtext(side = 2, line = 1, "Rate (y)")

segments(x0=45, x1=50, y0=u[3] + 0.8*(u[4]-u[3]), y1=u[3] + 0.8*(u[4]-u[3]), xpd=NA, col = cols[1], lwd = 1.5)
text(50, u[3] + 0.8*(u[4]-u[3]), pos = 4,"Constant rate", xpd = NA, cex = 1.2, col = cols[1])

segments(x0=45, x1=50, y0=u[3] + 0.7*(u[4]-u[3]), y1=u[3] + 0.7*(u[4]-u[3]), col = cols[2], xpd=NA, lwd = 1.5)
text(50, u[3] + 0.7*(u[4]-u[3]), pos = 4,"Arrhenius", xpd = NA, cex = 1.2, col = cols[2])
# text(45, u[3] + 0.6*(u[4]-u[3]), expression(y = y[0]*exp(-E/k*(1/T - 1/T[0]))), xpd = NA)

segments(x0=45, x1=50, y0=u[3] + 0.45*(u[4]-u[3]), y1=u[3] + 0.45*(u[4]-u[3]), col = cols[3], lty = 2, xpd=NA, lwd = 1.5)
text(50, u[3] + 0.45*(u[4]-u[3]), pos = 4,"Sharpe-Schoolfield", xpd = NA, cex = 1.2, col = cols[3])

segments(x0=45, x1=50, y0=u[3] + 0.2*(u[4]-u[3]), y1=u[3] + 0.2*(u[4]-u[3]), col = cols[4], lty = 2, xpd=NA, lwd = 1.5)
text(50, u[3] + 0.2*(u[4]-u[3]), pos = 4,"Inverse Sharpe-Schoolfield", xpd = NA, cex = 1.2, col = cols[4])

###############################################################################
# Inhibition scenarios
###############################################################################
z <- as.Date(paste("1985", c(1:365), sep = "-"), format = "%Y-%j")
DOY <- c(1:365)

pI4 <- c(rep(1, 109), 1/(1 + exp(-0.08*(c(110:365) - 172))))

quartz(width = 6.3, height = 3, pointsize = 10)
par(mfrow = c(1,1), mar = c(3,4,5,10))
plot(z, DOY, "n", ylim = c(0,1), ylab  = "Proportion arresting", las = 1, xlab = "", xaxs = "i")
axis(side = 3, at = z[c(1, seq(90, 365, 90))], labels = c(1, seq(90, 365, 90)))
lines(z, pI4, lwd = 3, col = cols[4])
abline(h = 1, col = cols[1], lwd = 1.5, lty = 1)
abline(h = 0.5, col = cols[2], lwd = 1.5, lty = 2)
abline(h = 0, col = cols[3], lwd = 1.5, lty = 3)
u <- par('usr')
arrows(x0 = z[110], x1 = z[110], y1 = u[3], y0 = 1.2, xpd = NA, length = 0.08)
text(z[110], 1.4, "Resumption of development \n Apr 20", xpd = NA)

# arrows(x0 = z[172], x1 = z[172], y1 = 0.5, y0 = 1.3, xpd = NA, length = 0.08)
abline(v = z[172], lty = 2)# text(z[172], 1.45, "Solstice\nJun 21", xpd = NA)


# arrows(x0 = z[232], x1 = z[232], y1 = 1, y0 = 1.3, xpd = NA, length = 0.08)
# text(z[232], 1.4, "Aug 20", xpd = NA)

legend(5857, 0.6, lwd = c(rep(1.5, 3),3), lty = c(1:3,1), col = cols[1:4], c("Scenario 1", "Scenario 2", "Scenario 3", "Scenario 4"), bty = "n", xpd = NA)

###############################################################################
# Temperature data vs model
###############################################################################

tempDat <- read.csv("parameterization/BAH_temps.csv")

tempDat$Date <- as.Date(tempDat$Date)
tempDat$DOY <- strftime(tempDat$Date, format = "%j")
tempDat$year <- strftime(tempDat$Date, format = "%Y")

season <- unique(tempDat$Range)
# mean temperature c("Summer",  "Winter",  "Spring",  "Fall", "Calving")
ck <- c(-10.974553,  -8.029248,  -9.808518,  -9.750996, -11.466570) 
# half annual temperature range (amplitude)
dk <- c(22.546752,  23.653142,  22.808474,  22.990450, 21.577856) 
# DOY corresponding to max temperature
t0 <- c(199.922313, 198.146132, 199.304408, 199.167954, 201.066183)

# tempChanges <- cbind(mean = c(0, 2, 5, 10, rep(0, 4)), range = c(rep(0, 4), c(0, 2, 5, 10)))
# 
# par(mfrow = c(2,2), mar = c(4,4,2,1))
# 
# # Plot two different seasons, different years
# 
# ind1 <- which(tempDat$year == 1980)
# ind2 <- which(tempDat$year == 2000)
# ind3 <- which(tempDat$year == 2019)
# 
# quartz(width = 6.3, height = 4, pointsize = 10)
# par(mfrow = c(2,2), mar = c(3, 1, 0, 0), oma = c(0, 4, 1, 1))
# 
# # for(i in 1:2){
# i <- 2
# 	season.i <- c("Calving", "Winter")[i]
# 	I <- match(season.i, c("Summer",  "Winter",  "Spring",  "Fall", "Calving"))
# 	
# 	plot(tempDat$Date[ind1][tempDat$Range[ind1] == season.i], tempDat$T2Mean[ind1][tempDat$Range[ind1] == season.i], "l", ylim = c(-45, 20), col = cols2[season.i], las = 1, xlab = "", ylab = "", yaxt = "n")
# 	axis(side = 2, at = seq(-40, 20, 20), las = 1)
# 	lines(tempDat$Date[ind1][tempDat$Range[ind1] == season.i], tempDat$T2Mean[ind2][tempDat$Range[ind2] == "Calving"], col= cols2[season.i])
# 	lines(tempDat$Date[ind1][tempDat$Range[ind1] == season.i], tempDat$T2Mean[ind3][tempDat$Range[ind2] == season.i], col= cols2[season.i])
# 	lines(tempDat$Date[ind1][tempDat$Range[ind1] == season.i], ck[I] + dk[I] * cos((1:366 - t0[I])* 2 * pi / 365), col = cols[season.i], lwd = 2)
# 	mtext(side = 3, line = -1.5, paste(" ", LETTERS[1]), adj = 0)
# # }
# # legend("topleft", lwd= c(1, 2), col = c(cols2[season.i], cols[season.i]), c("Daily mean", "Fitted sine curve"), bty = "n")
# 
# plot(tempDat$Date[ind1][tempDat$Range[ind1] == season.i], tempDat$T2Mean[ind1][tempDat$Range[ind1] == season.i], "n", ylim = c(-45, 20), col = cols2[season.i], las = 1, xlab = "", ylab = "", yaxt = "n")
# axis(side = 2, at = seq(-40, 20, 20), labels = FALSE)
# for(i in 1:5){
# 	season.i <- c("Summer",  "Winter",  "Spring",  "Fall", "Calving")[i]
# 	I <- match(season.i, c("Summer",  "Winter",  "Spring",  "Fall", "Calving"))
# 	lines(tempDat$Date[ind1][tempDat$Range[ind1] == season.i], ck[I] + dk[I] * cos((1:366 - t0[I])* 2 * pi / 365), col = cols[season.i])
# }
# mtext(side = 3, line = -1.5, paste(" ", LETTERS[2]), adj = 0)
# legend("bottom", fill = cols[c(5,3,1,2,4)], names(cols)[c(5,3,1,2,4)], cex = 0.8, border = NA, bty = "n", title = "Range")
# 
# season.i <- "Winter"
# I <- match(season.i, c("Summer",  "Winter",  "Spring",  "Fall", "Calving"))
# 
# plot(tempDat$Date[ind1][tempDat$Range[ind1] == season.i], ck[I] + dk[I] * cos((1:366 - t0[I])* 2 * pi / 365), "l", , las = 1, xlab = "", ylab = "", ylim = c(-45, 28), lwd = 2, col = cols[season.i])
# # axis(side = 2, labels = FALSE)
# for(i in 1:3){
# 	lines(tempDat$Date[ind1][tempDat$Range[ind1] == season.i], (ck[I]+ c(2, 5, 10)[i]) + dk[I] * cos((1:366 - t0[I])* 2 * pi / 365), lty = i+1, col = cols[season.i]) 
# }
# mtext(side = 3, line = -1.5, paste(" ", LETTERS[3]), adj = 0)
# legend("bottom", lty = c(1:4), lwd = c(2, rep(1, 3)), c("Current", expression(paste("+2", degree, "C")), expression(paste("+5", degree, "C")), expression(paste("+10", degree, "C"))),cex = 0.8, bty = "n", col = cols[season.i], title = "Change in mean")
# 
# plot(tempDat$Date[ind1][tempDat$Range[ind1] == season.i], ck[I] + dk[I] * cos((1:366 - t0[I])* 2 * pi / 365), "l", , las = 1, xlab = "", ylab = "", ylim = c(-45, 28), lwd = 2, col = cols[season.i], yaxt = "n")
# axis(side = 2, labels = FALSE)
# for(i in 1:3){
# 	lines(tempDat$Date[ind1][tempDat$Range[ind1] == season.i], ck[I] + (dk[I]+ c(2, 5, 10)[i]) * cos((1:366 - t0[I])* 2 * pi / 365), lty = i+1, col = cols[season.i]) 
# }
# mtext(side = 3, line = -1.5, paste(" ", LETTERS[4]), adj = 0)
# legend("bottom", lty = c(1:4), lwd = c(2, rep(1, 3)), c("Current", expression(paste("+2", degree, "C")), expression(paste("+5", degree, "C")), expression(paste("+10", degree, "C"))),cex = 0.8, bty = "n", col = cols[season.i], title = "Change in range")
# 
# 
# mtext(side = 2, line = 2, outer = TRUE, expression(paste("Temperature (", degree, "C)")))
# 

#------------------------------------------------------------------------------
# Including RCP temperature changes

tChange <- read.csv("ClimateChangeData/tasSummaryCalving.csv")
parChange <- read.csv("ClimateChangeData/parChange.csv")

DOY <- c(1:365)
xDate <- as.Date(paste("2100", DOY, sep = "-"), format = "%Y-%j")
xSeason <- as.Date(paste("2100", c(3, 6, 9, 12), 01, sep = "-"))

quartz(width = 3.2, height = 6, pointsize = 10)
par(mfrow = c(3,1), mar = c(3,5,1,1))

# A) Real temperature data from MERRA
ind1 <- which(tempDat$year == 2001) #2000 was leap year!
season1 <- "Calving"
I <- match(season1, c("Summer",  "Winter",  "Spring",  "Fall", "Calving"))
# 	
plot(xDate, tempDat$T2Mean[ind1][tempDat$Range[ind1] == season1], "l", ylim = c(-45, 20), col = grey(0.6), las = 1, xlab = "", ylab = expression(paste("Temperature (", degree, "C)", sep = "")), yaxt = "n")
axis(side = 2, at = seq(-40, 20, 20), las = 1)
lines(xDate, ck[I] + dk[I] * cos((DOY - t0[I])* 2 * pi / 365), col = 1, lwd = 2)
mtext(side = 3, line = -1.5, paste(" ", LETTERS[1]), adj = 0)


# B) Climate change scenarios
plot(xDate, rep(1, 365), "n", col = 4, xlab = "", ylab = expression(paste("Temperature change (", degree, "C)", sep = "")), las = 1, ylim = range(tChange[, c('tas25', 'tas75')]), xaxs = "i")

# Add seasonal data from RCP
for(j in 1:2){ # for each model
	for(i in 2:4){ # for each season except winter (which spans year)
		polygon(x = xSeason[c(i-1, i-1, i, i)], y = c(tChange[tChange$model == c("rcp26", "rcp85")[j] & tChange$season == c("DJF", "MAM", "JJA", "SON")[i], c('tas25', 'tas75', 'tas75', 'tas25')]), col = cols2[c(4, 1)[j]], border = NA)
		segments(x0 = xSeason[i-1], x1 = xSeason[i], y0 = tChange[tChange$model == c("rcp26", "rcp85")[j] & tChange$season == c("DJF", "MAM", "JJA", "SON")[i], 'tas50'], col = cols[c(4, 1)[j]])
	}
	# FOr winter
	polygon(x = c(rep(as.Date("2100-01-01"), 2), xSeason[c(1,1)]), y = c(tChange[tChange$model == c("rcp26", "rcp85")[j] & tChange$season == "DJF", c('tas25', 'tas75', 'tas75', 'tas25')]), col = cols2[c(4, 1)[j]], border = NA)
	segments(x0 = as.Date("2100-01-01"), x1 = xSeason[1], y0 = tChange[tChange$model == c("rcp26", "rcp85")[j] & tChange$season == "DJF", 'tas50'], col = cols[c(4, 1)[j]])
	polygon(x = c(xSeason[c(4, 4)] , rep(as.Date("2100-12-31"), 2)), y = c(tChange[tChange$model == c("rcp26", "rcp85")[j] & tChange$season == "DJF", c('tas25', 'tas75', 'tas75', 'tas25')]), col = cols2[c(4, 1)[j]], border = NA)
	segments(x0 = xSeason[4], x1 = as.Date("2100-12-31"), y0 = tChange[tChange$model == c("rcp26", "rcp85")[j] & tChange$season == "DJF", 'tas50'], col = cols[c(4, 1)[j]])				
	
	lines(xDate, parChange[parChange$X == 'ck', j+1] + parChange[parChange$X == 'dk', j+1] * cos((DOY - parChange[parChange$X == 't0', j+1])* 2 * pi / 365), col = cols[c(4, 1)[j]], lty = 2)
	
}
mtext(side = 3, line = -1.5, paste(" ", LETTERS[2]), adj = 0)

# Bringing it together: temps under two scenarios
tempNow <- ck[I] + dk[I] * cos((DOY - t0[I])* 2 * pi / 365)
temp26 <- tempNow +	parChange[parChange$X == 'ck', 2] + parChange[parChange$X == 'dk', 2] * cos((DOY - parChange[parChange$X == 't0', 2])* 2 * pi / 365)
temp85 <- tempNow +	parChange[parChange$X == 'ck', 3] + parChange[parChange$X == 'dk', 3] * cos((DOY - parChange[parChange$X == 't0', 3])* 2 * pi / 365)
temps <- cbind(tempNow, temp26, temp85)

plot(xDate, tempNow, "l", lwd = 2, col = 1, xlab = "", ylab = expression(paste("Temperature (", degree, "C)", sep = "")), las = 1, ylim = c(-45, 20), yaxt = "n")
axis(side = 2, at = seq(-40, 20, 20), las = 1)
for(j in 1:2) lines(xDate, temps[, j+1], col = cols[c(4, 1)[j]], lwd = 2)
mtext(side = 3, line = -1.5, paste(" ", LETTERS[3]), adj = 0)
 legend("topright", lwd = 2, col = c(1,cols[c(4, 1)]), legend  = c("past", "RCP 2.6", "RCP 8.5"), bty = "n" )
 
#------------------------------------------------------------------------------
# Adjusting for ground temperatures
#------------------------------------------------------------------------------
tempPred <- list(
	air = list(
		current = predict.temp(timeDat = NA, ground = FALSE),
		low = predict.temp(timeDat = NA, climateScenario = "rcp26", ground = FALSE),
		high = predict.temp(timeDat = NA, climateScenario = "rcp85", ground = FALSE)),
	ground = list(
		current = predict.temp(timeDat = NA),
		low = predict.temp(timeDat = NA, climateScenario = "rcp26"),
		high = predict.temp(timeDat = NA, climateScenario = "rcp85"))
)
 	
 
quartz(width = 6.4, height = 4, pointsize = 10)
# layout(matrix(c(1,2,3,3), nrow = 2, byrow = TRUE))
par(mfrow= c(2,2), mar = c(2, 4, 1, 2), oma = c(1,1,0,0))
 
# A) Real temperature data from MERRA
ind1 <- which(tempDat$year == 2001) #2000 was leap year!
season1 <- "Calving"
I <- match(season1, c("Summer",  "Winter",  "Spring",  "Fall", "Calving"))
# 	
plot(xDate, tempDat$T2Mean[ind1][tempDat$Range[ind1] == season1], "l", ylim = c(-40, 20), col = grey(0.6), las = 1, xlab = "", ylab = expression(paste("Air temperature (", degree, "C)", sep = "")), yaxt = "n", xaxs = "i")
abline(h = 0, lty = 2)
axis(side = 2, at = seq(-40, 20, 20), las = 1)
lines(xDate, tempPred$air$current[, I], col = 1)
mtext(side = 3, line = -1.5, paste(" ", LETTERS[1]), adj = 0)


# B) Climate change scenarios
plot(xDate, rep(1, 365), "n", col = 4, xlab = "", ylab = expression(paste("Temperature change (", degree, "C)", sep = "")), las = 1, ylim = range(tChange[, c('tas25', 'tas75')]), xaxs = "i")

# Add seasonal data from RCP
for(j in 1:2){ # for each model
	for(i in 2:4){ # for each season except winter (which spans year)
		polygon(x = xSeason[c(i-1, i-1, i, i)], y = c(tChange[tChange$model == c("rcp26", "rcp85")[j] & tChange$season == c("DJF", "MAM", "JJA", "SON")[i], c('tas25', 'tas75', 'tas75', 'tas25')]), col = cols2[c(4, 1)[j]], border = NA)
		segments(x0 = xSeason[i-1], x1 = xSeason[i], y0 = tChange[tChange$model == c("rcp26", "rcp85")[j] & tChange$season == c("DJF", "MAM", "JJA", "SON")[i], 'tas50'], col = cols[c(4, 1)[j]])
	}
	# FOr winter
	polygon(x = c(rep(as.Date("2100-01-01"), 2), xSeason[c(1,1)]), y = c(tChange[tChange$model == c("rcp26", "rcp85")[j] & tChange$season == "DJF", c('tas25', 'tas75', 'tas75', 'tas25')]), col = cols2[c(4, 1)[j]], border = NA)
	segments(x0 = as.Date("2100-01-01"), x1 = xSeason[1], y0 = tChange[tChange$model == c("rcp26", "rcp85")[j] & tChange$season == "DJF", 'tas50'], col = cols[c(4, 1)[j]])
	polygon(x = c(xSeason[c(4, 4)] , rep(as.Date("2100-12-31"), 2)), y = c(tChange[tChange$model == c("rcp26", "rcp85")[j] & tChange$season == "DJF", c('tas25', 'tas75', 'tas75', 'tas25')]), col = cols2[c(4, 1)[j]], border = NA)
	segments(x0 = xSeason[4], x1 = as.Date("2100-12-31"), y0 = tChange[tChange$model == c("rcp26", "rcp85")[j] & tChange$season == "DJF", 'tas50'], col = cols[c(4, 1)[j]])				
	
	lines(xDate, parChange[parChange$X == 'ck', j+1] + parChange[parChange$X == 'dk', j+1] * cos((DOY - parChange[parChange$X == 't0', j+1])* 2 * pi / 365), col = cols[c(4, 1)[j]], lty = 2)
	
}
mtext(side = 3, line = -1.5, paste(" ", LETTERS[2]), adj = 0)
  
# C) Annual air temperature projected
for(i in 1:2){
	plot(xDate, tempPred$air$current[, I], "n", xlab = "", ylab = c(expression(paste("Air temperature (", degree, "C)", sep = "")), expression(paste("Ground temperature (", degree, "C)", sep = "")))[i], las = 1, ylim = c(-40, 20), yaxt = "n", xaxs = "i")
	if(i==2) abline(h = -15, lty = 3)
	axis(side = 2, at = seq(-40, 20, 20), las = 1)
	for(j in 1:3){
			lines(xDate, tempPred[[i]][[j]][, I], col = c(1, cols[c(4, 1)])[j])
		}
	mtext(side = 3, line = -1.5, paste(" ", LETTERS[2+i]), adj = 0)
	abline(h = 0, lty = 2)
}
legend("bottomright", lwd = 1, col = c(1,cols[c(4, 1)]), legend  = c("past", "low emissions", "high emissions"), bty = "n" )
	
###############################################################################
# Larvae parameters and MTE relationships
###############################################################################

par1 <- read.csv("MTE_dataAnalysis/output/bestModelParameters.csv")

# Function to convert character numbers from .csv above to numeric for plotting
convNum <- function(x, type = "m"){
	x.out <- length(x)
	for(i in 1:length(x)){
		dum <- strsplit(x[i], split = "()")[[1]]
		if(type=="m"){
			x.out[i] <- as.numeric(paste(as.numeric(dum[1]), ".", as.numeric(dum[3]), as.numeric(dum[4]), as.numeric(dum[5]), sep = ""))
		} else if(type == "l"){
			x.out[i] <- as.numeric(paste(as.numeric(dum[8]), ".", as.numeric(dum[10]), as.numeric(dum[11]), as.numeric(dum[12]), sep = ""))
		} else if(type == "u"){
			x.out[i] <- as.numeric(paste(as.numeric(dum[15]), ".", as.numeric(dum[17]), as.numeric(dum[18]), as.numeric(dum[19]), sep = ""))
		}
	} # end i
	
	return(x.out)
}

# Function to predict param over temp from MTE hyperparameters
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

# Function to plot true temperture range and 10th and 90th quantiles
plotTpolys <- function(){
	u <- par('usr')
	polygon(x = c(range(tempDat$T2Mean), rev(range(tempDat$T2Mean))), y = u[c(3,3,4,4)], border = NA, col = cols2[4])#"#00000030")
	polygon(x = c(quantile(tempDat$T2Mean, c(0.1, 0.9)), rev(quantile(tempDat$T2Mean, c(0.1, 0.9)))), y = u[c(3,3,4,4)], border = NA, col = cols2[5])#
}


temp <- seq(5, 35, 5)
temp.all <- seq(-50, 40, 0.2)
quartz(width = 8.3, height = 2.8, pointsize = 10)
par(mfrow = c(1,3), mar = c(3,3,2,1), oma = c(1, 2, 0, 8))

#------------------------------------------------------------------------------
#mu0
#------------------------------------------------------------------------------
plot(temp.all, rep(0, length(temp.all)), "n", las = 1, ylab = "", xlab = "", main = expression(mu[0]), xlim = c(-50, 40), ylim = c(0, 1))
plotTpolys()

lines(temp.all, predict.MTE(params = c(a = convNum(par1[9,3]), E = convNum(par1[10,3]), Eh = NA, Th = NA, z = 1), temp = temp.all), lty = 3, lwd = 1.5)

lines(temp.all, predict.mu0(temp.all))

plotCI(temp, convNum(par1[1:7, 'mu0'], type = "m"), li = convNum(par1[1:7, 'mu0'], type = "l"), ui = convNum(par1[1:7, 'mu0'], type = "u"), gap = 0,  sfrac= 0.008, add = TRUE, pch = 21, pt.bg = "white")

# Points from Aleuy et al. (2020) Marshallagia study
# freezing survival of L1
points(-9, 0.0855, pch = 8, xpd = NA)
points(-20, 0.9857, pch = 8, xpd = NA)
mtext(side = 3, line = -1.5, paste(" ", LETTERS[1]), adj = 0)

#------------------------------------------------------------------------------
# mu1
#------------------------------------------------------------------------------

plot(temp.all, rep(0, length(temp.all)), "n", las = 1, ylab = "", xlab = "", main = expression(mu[1]), xlim = c(-50, 40), ylim = c(0, 0.1))
plotTpolys()

lines(temp.all, rep(convNum(par1[9, 'mu1']), length(temp.all)), lty = 3, lwd = 1.5)

lines(temp.all, predict.mu3(temp.all))

plotCI(temp[1:6], convNum(par1[1:6, 'mu1'], type = "m"), li = convNum(par1[1:6, 'mu1'], type = "l"), ui = convNum(par1[1:6, 'mu1'], type = "u"), gap = 0,  sfrac= 0.008, add = TRUE, pch = 21, pt.bg = "white")

# Points from Aleuy et al. (2020) Marshallagia study
# freezing survival of L3
points(-9, 0.002344, pch = 8, xpd = NA)
points(-20, 0.01762, pch = 8, xpd = NA)
# points(-35, 3, pch = 8, xpd = NA)
mtext(side = 3, line = -1.5, paste(" ", LETTERS[2]), adj = 0)

#------------------------------------------------------------------------------
# rho
#------------------------------------------------------------------------------

plot(temp.all, rep(0, length(temp.all)), "n", las = 1, ylab = "", xlab = "", main = expression(rho), xlim = c(-50, 40), ylim = c(0, 0.1))
plotTpolys()

lines(temp.all, predict.MTE(params = c(a = convNum(par1[9, 'rho']), E = convNum(par1[10, 'rho']), Eh = convNum(par1[11, 'rho']), Th = 30.568, z = -1), temp = temp.all), lty = 3, lwd = 1.5)

lines(temp.all, predict.rho0(temp.all))
plotCI(temp[1:7], convNum(par1[1:7, 'rho'], type = "m"), li = convNum(par1[1:7, 'rho'], type = "l"), ui = convNum(par1[1:7, 'rho'], type = "u"), gap = 0,  sfrac= 0.008, add = TRUE, pch = 21, pt.bg = "white")

legend(48, 0.1, pch = c(21, 8, NA, NA), lty = c(NA, NA, 3, 1), lwd = c(NA, NA, 1.5, 1), legend = c("Exp't estimates", "Marshallagia", "Fitted MTE", "Assumed MTE"), xpd = NA, bty = "n")

mtext(side = 3, line = -1.5, paste(" ", LETTERS[3]), adj = 0)
mtext(side = 1, outer = TRUE, expression(paste("Temperature (", degree, "C)")))
mtext(side = 2, outer = TRUE, "Parameter value")
