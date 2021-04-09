
source("simulations/bouSetup.R")
source("simulations/popFunctions.R")

###############################################################################################
# One-year simulation to look at host population changes
# No parasites
###############################################################################################
y <- 1
dt <- 1
stochTemp <- FALSE
startstop <- "agg"
source("simulations/calcGrid.R")

initBou1 <- rbind(
	calf_mov = rep(0, length(x)),
	yearling_mov = rep(0, length(x)),
	adult_mov = rep(0, length(x)),
	L4_mov = rep(0, length(x)),
	L4A_mov = rep(0, length(x)),
	P_mov = rep(0, length(x)),
	calf_stat = initDist(totPop = bouDist1985['calf'], x),
	yearling_stat = initDist(totPop = bouDist1985['yearling'], x),
	adult_stat = initDist(totPop = c(bouDist1985['cow'] + bouDist1985['bull']), x),
	L4_stat = rep(0, length(x)),
	L4A_stat = rep(0, length(x)),
	P_stat = rep(0, length(x)),
	L0 = rep(0 * dx, length(x)),
	L3 = rep(0 * dx, length(x))
)


V1 <- simBou(
	initBou = initBou1, 
	tempGrid = tempGrid,
	ppnInhibit = 1)

plot.timestep(V1, n = which(timeDat$DOY == as.numeric(strftime(as.Date("1985-06-15"), format = "%j")))[1])

plot(timeDat$time, apply(V1['adult_mov', , ], 2, sum)*dx + apply(V1['adult_stat', , ], 2, sum)*dx, "l", ylim = c(0, 450000))
lines(timeDat$time, apply(V1['adult_mov', , ], 2, sum)*dx, col = 3)
lines(timeDat$time, apply(V1['adult_stat', , ], 2, sum)*dx, col = 2)


#------------------------------------------------------------------------------
# Produce GIF for presentation
#------------------------------------------------------------------------------

Nrange <- range(V1[c('adult_mov', 'adult_stat'), , ])
V <- V1

days2plot <- DOY
nInd <- match(round(days2plot, 2), round(timeDat$time, 2))
intervals <- rep(0.05, length(nInd))
intervals[which(timeDat$DOY[nInd] == 99) + c(-4:4)] <- 0.2 # spring migration
intervals[which(timeDat$DOY[nInd] == 152) + c(-4:4)] <- 0.2 # calving
intervals[which(timeDat$DOY[nInd] == 168) + c(-4:4)] <- 0.2 # post-calving mig
intervals[which(timeDat$DOY[nInd] == 180) + c(-4:4)] <- 0.2 # summer
intervals[which(timeDat$DOY[nInd] == 250) + c(-4:4)] <- 0.2 # summer

Sys.setenv(PATH="/opt/ImageMagick/bin")
ani.options(interval = intervals, ani.width = 600, ani.height = 400, ani.dev="png")

saveGIF(expr = {
	for(i in nInd){
		par(mfrow = c(1,1), mar = c(4,4,4,1), family = "Helvetica")
		# Moving
		plot(x, V['adult_mov', , 1], "n", ylim = Nrange, bty = "l", yaxt = "n", ylab = "", lwd = 1.5, xlab="", yaxs = "i", xaxs = "i")
		
		# u <- par('usr')
		# points(x, rep(u[3] - (u[4] - u[3])*0.01, length(x)), pos = 1, col = seasonCols[xRange], pch = 15, xpd = NA)
		
		axis(side = 1, labels = FALSE, lwd = 1.5)
		yTick <- pretty(Nrange)
		axis(side = 2, las = 1, at = yTick, labels = sprintf("%.1f", yTick/1000))
		
		lines(x, V['adult_mov', , i], col = 1, lwd = 1.5)#DOYcol[timeDat$DOY[i]])
		lines(x, V['adult_stat', , i], col = 1, lty = 2, lwd = 2)# DOYcol[timeDat$DOY[i]],lty = 2, lwd = 2)
		
		mtext(side = 3, line = 1, strftime(as.Date(paste("1985", timeDat$DOY[i], sep = "-"), format = "%Y-%j"), format = "%b %d"))
		
		lines(x, stopMat[i, ]*Nrange[2], col = 2, xpd = NA)
		lines(x, startMat[i, ]*Nrange[2], col = 3, xpd =NA)
		
		legend("topleft", lwd = c(1.5, 2, 1, 1), col = c(1,1,2,3), c("Moving hosts", "Stationary hosts", "Stopping rate", "Starting rate"), lty = c(1, 2, 1, 1), bty = "n")
		
		if(is.element(timeDat$DOY[i], c(111:153))) mtext(side = 3, line = 2, "Spring migration", font =2, xpd = NA)
		
		if(timeDat$DOY[i] == breedDOY) mtext(side = 3, line = 2, "Peak calving", font =2, xpd = NA) else if(is.element(timeDat$DOY[i], c(154:168))) mtext(side = 3, line = 2, "Calving season", font =2, xpd = NA)
		
		if(is.element(timeDat$DOY[i], c(169:180))) mtext(side = 3, line = 2, "Post-calving migration", font =2, xpd = NA)
		if(is.element(timeDat$DOY[i], c(181:250))) mtext(side = 3, line = 2, "Summer range", font =2, xpd = NA)
		if(is.element(timeDat$DOY[i], c(251:290))) mtext(side = 3, line = 2, "Fall migration", font =2, xpd = NA)
		if(is.element(timeDat$DOY[i], c(291:305))) mtext(side = 3, line = 2, "Breeding season", font =2, xpd = NA)
		if(is.element(timeDat$DOY[i], c(306:335))) mtext(side = 3, line = 2, "Fall migration", font =2, xpd = NA)
		if(is.element(timeDat$DOY[i], c(336:365))) mtext(side = 3, line = 2, "Winter range", font =2, xpd = NA)
	}},
	, movie.name = "hostMovememt.gif", 
	img.name = "hostMovement")


###############################################################################################
# Base simulation starting with parasites
###############################################################################################

# Starting abundance of parasites in adults
# From Bathurst surveys, appprox. 
initP <- 914

y <- 30
dt <- 1
stochTemp <- FALSE
startstop <- "agg"

source("simulations/calcGrid.R")

initBou <- rbind(
	calf_mov = rep(0, length(x)),
	yearling_mov = rep(0, length(x)),
	adult_mov = rep(0, length(x)),
	L4_mov = rep(0, length(x)),
	L4A_mov = rep(0, length(x)),
	P_mov = rep(0, length(x)),
	calf_stat = initDist(totPop = bouDist1985['calf'], x),
	yearling_stat = initDist(totPop = bouDist1985['yearling'], x),
	adult_stat = initDist(totPop = c(bouDist1985['cow'] + bouDist1985['bull']), x),
	L4_stat = rep(0, length(x)),
	L4A_stat = rep(0, length(x)),
	P_stat = initP * initDist(totPop = 10^5 * c(bouDist1985['cow'] + bouDist1985['bull'])/sum(bouDist1985), x),
	L0 = rep(0, length(x)),
	L3 = rep(0, length(x))
)

V <- simBou(
	initBou = initBou, 
	tempGrid = tempGrid,
	ppnInhibit = 1)

plotOverTime(years2plot = 1:30, V1 = V)

# Hosts
plot(timeDat$time, apply(V['adult_mov', , ], 2, sum)*dx + apply(V['adult_stat', , ], 2, sum)*dx, "l", ylim = c(0, 450000))
lines(timeDat$time, apply(V['adult_mov', , ], 2, sum)*dx, col = 3)
lines(timeDat$time, apply(V['adult_stat', , ], 2, sum)*dx, col = 2)

# Parasites
plot(timeDat$time, apply(V['P_mov', , ], 2, sum)*dx + apply(V['P_stat', , ], 2, sum)*dx, "l", ylim = c(0, 450000))
lines(timeDat$time, apply(V['P_mov', , ], 2, sum)*dx, col = 3)
lines(timeDat$time, apply(V['P_stat', , ], 2, sum)*dx, col = 2)
lines(timeDat$time, apply(V['L4A_mov', , ], 2, sum)*dx + apply(V['L4A_stat', , ], 2, sum)*dx, lty = 2)
lines(timeDat$time, apply(V['L4A_mov', , ], 2, sum)*dx, col = 3, lty = 2)
lines(timeDat$time, apply(V['L4A_stat', , ], 2, sum)*dx, col = 2, lty = 2)


#------------------------------------------------------------------------------
V_base <- simBou(initBou = initBou, tempGrid = tempGrid, ppnInhibit = 1)
V_noInhibit <- simBou(initBou = initBou, tempGrid = tempGrid, ppnInhibit = 0)
V_plus5 <- simBou(initBou = initBou, tempGrid = tempGrid + 5, ppnInhibit = 0)
V_plus5inhibit <- simBou(initBou = initBou, tempGrid = tempGrid + 5, ppnInhibit = 1)

plotOverTime(years2plot = 1:30, V1 = V)

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------




###############################################################################
# Plot output
###############################################################################
ranges <- list(N = range(V[c('adult_mov', "adult_stat"), , ]), P = range(V[c('P_mov', "P_stat"), , ]), L = range(V[c('L0', "L3"), , ]))
plot.timestep(V, n = which(timeDat$DOY == 300)[1], Nrange = ranges$N, Prange = ranges$P, Lrange = ranges$L)
#------------------------------------------------------------------------------
# Animation
#------------------------------------------------------------------------------

saveHTML(expr = {
	for(n in round(seq(timeDat$step[1], which(timeDat$year == 10)[1], length.out = 120))){
		plot.timestep(V, n)
	}
},
   title = "OsterBou pop dynamics", 
   description = "Animation loops through time", 
   verbose = FALSE)

#------------------------------------------------------------------------------
# Plot change in one location over time
#------------------------------------------------------------------------------

plot(as.Date(paste(1984+timeDat$year, timeDat$DOY, sep="-"), format = "%Y-%j"), V['L0', 181, ], "l")
plot(as.Date(paste(1984+timeDat$year, timeDat$DOY, sep="-"), format = "%Y-%j"), V['P_mov', 181, ], "l")
plot(as.Date(paste(1984+timeDat$year, timeDat$DOY, sep="-"), format = "%Y-%j"), V['P_stat', 181, ], "l")


plot(as.Date(paste(1984+timeDat$year, timeDat$DOY, sep="-"), format = "%Y-%j"), V['yearling_mov', 182, ], "l")

ind1 <- which(paste(1984+timeDat$year, timeDat$DOY, sep="-") == "1986-150")
ind2 <- which(paste(1984+timeDat$year, timeDat$DOY, sep="-") == "1986-165")
ind3 <- which(paste(1984+timeDat$year, timeDat$DOY, sep="-") == "1986-180")

lines(as.Date(paste(1984+timeDat$year, timeDat$DOY, sep="-"), format = "%Y-%j"), V['yearling_mov', 192, ], col = grey(0.8))
abline(v = as.Date(paste(1984+timeDat$year, timeDat$DOY, sep="-"), format = "%Y-%j")[c(ind1, ind2, ind3)], col = c(2:4))

plot(x, V['yearling_mov', , ind1[1]], "l", col = 2)
lines(x, V['yearling_mov', , ind2[1]], col = 3)
lines(x, V['yearling_mov', , ind3[1]], col = 4)


# CHange in total parasite and total host population over time
par(mfrow=c(3,1))
plot(as.Date(paste(1984+timeDat$year, timeDat$DOY, sep="-"), format = "%Y-%j"), apply(V['adult_mov', , ], 2, sum), "l", main = "Total hosts over time")
plot(as.Date(paste(1984+timeDat$year, timeDat$DOY, sep="-"), format = "%Y-%j"), apply(V['P_mov', , ], 2, sum), "l", main = "Total P over time")
plot(as.Date(paste(1984+timeDat$year, timeDat$DOY, sep="-"), format = "%Y-%j"), apply(V['L3', , ], 2, sum), "l", main = "Total L3 over time")
lines(as.Date(paste(1984+timeDat$year, timeDat$DOY, sep="-"), format = "%Y-%j"), apply(V['L0', , ], 2, sum), col = 4)

plot(as.Date(paste(1984+timeDat$year, timeDat$DOY, sep="-"), format = "%Y-%j"), apply(V['L0', , ], 2, sum), , ylim = c(0, 10^3))

#--------------
DOYcol <- colorRampPalette(colors=seasonCols[c("Winter", "Spring", "Calving", "Summer", "Fall", "Winter")])(n = 365)

par(mfrow = c(2,1), mar = c(4,4,2,1))
# Moving
plot(x, V['adult_mov', , 1], "l", col = DOYcol[timeDat$DOY[1]],  ylim = range(V['adult_mov', , ]), bty = "l", xaxt="n", yaxt = "n", ylab = "", lwd = 1.5)
axis(side = 1, labels = FALSE)
yTick <- pretty(Nrange)
axis(side = 2, las = 1, at = yTick, labels = yTick/1000)
for(i in round(seq(1, dim(timeDat)[1], length.out=12))){
	lines(x, V['adult_mov', , i], col = DOYcol[timeDat$DOY[i]])
}

# Stationary
plot(x, V['adult_stat', , 1], "l", col = DOYcol[timeDat$DOY[1]], ylim = range(V['adult_stat', , ]), bty = "l", xaxt="n", yaxt = "n", ylab = "", lwd = 1.5)
axis(side = 1, labels = FALSE)
yTick <- pretty(Nrange)
axis(side = 2, las = 1, at = yTick, labels = yTick/1000)
for(i in round(seq(1, dim(timeDat)[1], length.out=12))){
	lines(x, V['adult_stat', , i], col = DOYcol[timeDat$DOY[i]])
}

# nInd <- which((timeDat$time - timeDat$DOY) == 0)[c(1,31,59,90,99:180,200, 230, 250:345, 365)]

#------------------------------------------------------------------------------
# Compare basic simulation with agg and simple host movement
#------------------------------------------------------------------------------
V_agg; V_simple

Ntot <- list(agg = apply(V_agg['adult_mov', , ] + V_agg['adult_stat', , ], 2, sum), simple = apply(V_simple['adult_mov', , ] + V_simple['adult_stat', , ], 2, sum))
Pmean <- list(agg = apply((V_agg['P_stat', , ] + V_agg['P_mov', , ]), 2, sum)/Ntot$agg, simple = apply((V_simple['P_stat', , ] + V_simple['P_mov', , ]), 2, sum)/Ntot$simple)
L4Amean <- apply((V['L4A_stat', , ] + V['L4A_mov', , ]), 2, sum)/Ntot
L4mean <- apply((V['L4_stat', , ] + V['L4_mov', , ]), 2, sum)/Ntot
L0tot <- list(agg = apply(V_agg['L0', , ], 2, sum), simple = apply(V_simple['L0', , ], 2, sum))
L3tot <- list(agg = apply(V_agg['L3', , ], 2, sum), simple = apply(V_simple['L3', , ], 2, sum))

# varOverTime <- list()
# varOverTime$base <- cbind(Ntot, Pmean, L4Amean, L4mean, L0tot, L3tot)
# varOverTime$noInhibit <- cbind(Ntot, Pmean, L4Amean, L4mean, L0tot, L3tot)
# varOverTime$plus5 <- cbind(Ntot, Pmean, L4Amean, L4mean, L0tot, L3tot)
# varOverTime$plus5inhibit <- cbind(Ntot, Pmean, L4Amean, L4mean, L0tot, L3tot)


y2plot <- 5:10
tInd <- which(is.element(timeDat$year, y2plot))
dates2plot <- as.Date(paste(timeDat$year, timeDat$DOY, sep="-"), format = "%Y-%j")

par(mfrow=c(2,1), mar = c(3,4,1,5), oma = c(1,0,2,0))


# plot(dates2plot, Ntot$agg, type = "l", ylim = range(Ntot))
# lines(dates2plot, Ntot$simple, col = 2)
# 
# par(new= TRUE)
# plot(dates2plot, Pmean$agg, "l", lty = 3, bty = "n", xaxt="n", yaxt="n", xlab = "", ylab = "", ylim = c(0, 3))

plot(dates2plot, Pmean$agg, type = "l", ylim = c(0,3), col = 2, las = 1)
lines(dates2plot, Pmean$simple)

plot(dates2plot, L0tot$agg, type = "l", ylim = range(L0tot), col = 2, las = 1)
lines(dates2plot, L0tot$simple)

###############################################################################
# Stochastic temperature
###############################################################################

y <- 30
dt <- 1
startstop <- "agg"

V_stochTemp <- list(); length(V_stochTemp) <- 2
for(s in 1:2){
	stochTemp <- c(FALSE, TRUE)[s]
	source("simulations/calcGrid.R")
	
	V_stochTemp[[s]] <- simBou(
		initBou = initBou, 
		tempGrid = tempGrid,
		ppnInhibit = 1)

}


#----------------------------------------------------------------------------
j <- 10015
plot(x, apply(V_stochTemp[[2]][c('P_mov', 'P_stat'), , j], 2, sum), "l")
lines(x, apply(V_stochTemp[[1]][c('P_mov', 'P_stat'), , j], 2, sum), col = 2)

plot(x, V_stochTemp[[2]]['L0', , j], "l")
lines(x, V_stochTemp[[1]]['L0', , j], col = 2)

# Differenes are VERY small when using stochastic temps

plotOverTime(years2plot = 1:30, V1 = V_stochTemp[[1]], V2 = V_stochTemp[[2]])

#----------------------------------------------------------------------------
stochTemp <- FALSE
source("simulations/calcGrid.R")
tempGrid_det <- tempGrid

stochTemp <- TRUE
source("simulations/calcGrid.R")
tempGrid_stoch <- tempGrid

# Summaries
Ntot <- list(
	det = apply(V_stochTemp[[1]]['adult_mov', , ] + V_stochTemp[[1]]['adult_stat', , ], 2, sum),
	stoch = apply(V_stochTemp[[2]]['adult_mov', , ] + V_stochTemp[[2]]['adult_stat', , ], 2, sum)
)

Pmean = list(
	det = apply(V_stochTemp[[1]]['P_mov', , ] + V_stochTemp[[1]]['P_stat', , ], 2, sum)/Ntot$det,
	stoch = apply(V_stochTemp[[2]]['P_mov', , ] + V_stochTemp[[2]]['P_stat', , ], 2, sum)/Ntot$stoch
)

L3tot <- list(
	det = apply(V_stochTemp[[1]]['L3', , ], 2, sum) * dx,
	stoch = apply(V_stochTemp[[2]]['L3', , ], 2, sum) * dx
)
	
# Plot
years2plot <- 3:4
tInd <- which(is.element(timeDat$year, years2plot))
dates2plot <- as.Date(paste(timeDat$year[tInd], timeDat$DOY[tInd], sep="-"), format = "%Y-%j")

pdf(file = "simulations/figures/StochTemp.pdf", width = 6.6, height = 7, pointsize = 10)
par(mfrow = c(3,1), mar = c(4,6,2,1), oma = c(3,0,0,0))
plot(as.Date(paste(rep(years2plot, each = 365), rep(c(1:365), length(years2plot)), sep = "-"), format  = "%Y-%j"), tempGrid_stoch[(years2plot[1]*365 + 1):((max(years2plot)+1)*365), 1], "l", col = "#ff6361", las = 1, bty = "l", ylab = "", xlab = "", xaxt = "n")
axis(side = 1, at = as.Date(paste(c(years2plot, tail(years2plot)+1), 1, sep="-"), format = "%Y-%j"), labels = c(years2plot, tail(years2plot)+1))
mtext(side = 2, line = 4.5, expression(paste("Temperature (", degree, "C)", sep = "")))
lines(as.Date(paste(rep(years2plot, each = 365), rep(c(1:365), length(years2plot)), sep = "-"), format  = "%Y-%j"), tempGrid_det[(years2plot[1]*365 + 1):((max(years2plot)+1)*365), 1], col = "#003f5c", lwd = 1.5)
legend("topright", lwd = c(1.5,  1), col = c("#003f5c", "#ff6361"), legend = c("Smoothed", "Stochastic"), bty = "n")
mtext(side= 3, adj = 0, line = -0.5, "   a")

#-----
par(mar = c(2,6,1,1))

plot(dates2plot, Pmean$det[tInd], "n", bty = "l", las = 1, ylab = "", xaxt = "n", xlab = "")
axis(side = 1, at = as.Date(paste(c(years2plot, tail(years2plot)+1), 1, sep="-"), format = "%Y-%j"), labels = c(years2plot, tail(years2plot)+1))
lines(dates2plot, Pmean$stoch[tInd], col = "#ff6361", xpd = NA)
lines(dates2plot, Pmean$det[tInd], lwd = 1.5, col = "#003f5c")
mtext(side = 2, line = 4.5, "Adult parasite burden")
mtext(side= 3, adj = 0, line = -0.5, "   b")

plot(dates2plot, L3tot$det[tInd], "n", bty = "l", las = 1, ylab = "", xaxt = "n", xlab = "")
axis(side = 1, at = as.Date(paste(c(years2plot, tail(years2plot)+1), 1, sep="-"), format = "%Y-%j"), labels = c(years2plot, tail(years2plot)+1))
lines(dates2plot, L3tot$stoch[tInd], col = "#ff6361", xpd = NA)
lines(dates2plot, L3tot$det[tInd], lwd = 1.5, col = "#003f5c")
mtext(side = 2, line = 4.5, "Free-living infectious larvae")
mtext(side= 3, adj = 0, line = -0.5, "   c")

mtext(side = 1, line = 3, "Year in simulation")
dev.off()

###############################################################################
# Temperature changes
###############################################################################
tempChanges <- cbind(mean = c(0, 2, 5, 10, rep(0, 4)), range = c(rep(0, 4), c(0, 2, 5, 10)))

y <- 30
dt <- 1
startstop <- "agg"
stochTemp <- FALSE
source("simulations/calcGrid.R")

# Parallelize over multiple cores
registerDoParallel(cores = 4)

ptime <- system.time({ 
	outTempChanges2 <- foreach (j = 1:dim(tempChanges)[1]) %dopar% {
		
		tempGrid <- predict.temp(timeDat = timeDat, stoch = FALSE, meanChange = tempChanges[j, 1], rangeChange = tempChanges[j,2])
		
		V <- simBou(
			initBou = initBou, 
			tempGrid = tempGrid,
			ppnInhibit = 1)
		
		maxP <- numeric(y)
		for(i in 1:y){
			ind <- which(timeDat$year == i)
			non0 <- list(stat = which(V['adult_stat',,ind] > 0, arr.ind = TRUE), mov = which(V['adult_mov',,ind] > 0, arr.ind = TRUE))
			maxP[i] <- max(c(V['P_stat',,ind][non0$stat]/V['adult_stat',,ind][non0$stat], V['P_mov',,ind][non0$mov]/V['adult_mov',,ind][non0$mov]), na.rm = TRUE)
		}
		
		mean(tail(log(maxP[2:y]/maxP[1:(y-1)]), 5))
		
	} # end over all parList
	
})[3]

print(paste("Run time:", ptime/60))

pdf(file = "simulations/figures/tempChanges.pdf", width = 5, height = 3, pointsize = 8)
par(mfrow = c(1,1), mar = c(5,4,2,1))

plot(c(1:4, 6:9), unlist(outTempChanges2), xaxt="n", xlim= c(0.5, 9.5), col = rep(c("#003f5c", "#7a5195",	"#ef5675", "#ffa600"), 2), cex = 2, pch = 1, ylim = c(-0.3, 0.2), lwd = 2, las = 1, ylab = "Annual change in peak parasite burden", xlab = "")
abline(h = 0)
abline(v = 5)

points(c(1:4, 6:9), unlist(outTempChanges), col = rep(c("#003f5c", "#7a5195",	"#ef5675", "#ffa600"), 2), cex = 2, pch = 19)

axis(side = 1, at = 1:4, labels = c(0, "+2", "+5", "+10"))
axis(side = 1, at = 6:9, labels = c(0, "+2", "+5", "+10"))

text(2.5, -0.45, "Change in\nmean temperature", xpd = NA)
text(7.5, -0.45, "Change in\ntemperature range", xpd = NA)

legend("bottomleft", pch = c(19, 1), pt.cex = 2, col = grey(0.8), c("100% arrested", "0% arrested"), bty = "n")

dev.off()
###############################################################################
# Proportion that inhibit
###############################################################################


# What is the annual rate of change in peak abundance?
maxP <- numeric(y)
for(i in 1:y){
	ind <- which(timeDat$year == i)
	non0 <- list(stat = which(V['adult_stat',,ind] > 0, arr.ind = TRUE), mov = which(V['adult_mov',,ind] > 0, arr.ind = TRUE))
	
	maxP[i] <- max(c(V['P_stat',,ind][non0$stat]/V['adult_stat',,ind][non0$stat], V['P_mov',,ind][non0$mov]/V['adult_mov',,ind][non0$mov]), na.rm = TRUE)
}

plot(1:y, maxP)
plot(2:y, log(maxP[2:y]/maxP[1:(y-1)]), ylim = c(0, 0.5))
abline(h = 0)
abline(h = mean(tail(log(maxP[2:y]/maxP[1:(y-1)]), 5)), col = 2)

#-------------
ppnInhibit.all <- seq(0, 1, 0.02)

# Parallelize over multiple cores
registerDoParallel(cores = 5)

ptime <- system.time({ 
	out5 <- foreach (j = 1:length(ppnInhibit.all)) %dopar% {
		
	V <- simBou(
		initBou = initBou, 
		tempGrid = tempGrid + 5,
		ppnInhibit = ppnInhibit.all[j])
	
	maxP <- numeric(y)
	for(i in 1:y){
		ind <- which(timeDat$year == i)
		non0 <- list(stat = which(V['adult_stat',,ind] > 0, arr.ind = TRUE), mov = which(V['adult_mov',,ind] > 0, arr.ind = TRUE))
		maxP[i] <- max(c(V['P_stat',,ind][non0$stat]/V['adult_stat',,ind][non0$stat], V['P_mov',,ind][non0$mov]/V['adult_mov',,ind][non0$mov]), na.rm = TRUE)
	}
		
	mean(tail(log(maxP[2:y]/maxP[1:(y-1)]), 5))
	
	} # end over all parList
	
})[3]

print(paste("Run time:", ptime/60))

ppnIn_trial1 <- data.frame(ppnInhibit = seq(0,1,0.05), annualChange = unlist(out_ppnInhinit))
ppnIn_trial2 <- data.frame(ppnInhibit = seq(0,1,0.02), annualChange = unlist(out))
ppnIn_trial5 <- data.frame(ppnInhibit = seq(0,1,0.02), annualChange = unlist(out5))

quartz(width = 4, height = 2.8, pointsize = 10)   
par(mfrow = c(1,1), mar = c(4,4,2,1))
plot(ppnIn_trial2$ppnInhibit, ppnIn_trial2$annualChange, "n", xlab = "Proportion of ingested larvae entering arrested state", ylab = "Annual change in peak parasite burden", bty = "l", xaxs = "i", las = 1, ylim = c(-0.25, 0.25))
points(ppnIn_trial2$ppnInhibit, ppnIn_trial2$annualChange, "o", pch = 21, bg = "#003f5c", xpd = NA)
abline(h = 0)
points(ppnIn_trial5$ppnInhibit, ppnIn_trial5$annualChange, "o", pch = 21, bg = "#ffa600", xpd = NA)
legend("topleft", pch = 21, pt.bg = c("#003f5c", "#ffa600"), legend = c("Current", expression(paste("+5", degree, "C"))), bty = "n", lwd = 1)

# Plot of metric

plot(x, V['adult_mov', , 160] + V['adult_stat', , 160], "l", lty = 2)
lines(x, V['adult_mov', , 140] + V['adult_stat', , 140])

par(new=TRUE)
plot(x, V['P_mov', , 140] + V['P_stat', , 140], "l", col = 2, yaxt = "n")
axis(side = 4, col = 2)
lines(x, V['P_mov', , 160] + V['P_stat', , 160], col = 2, lty = 2)

par(new = TRUE)
plot(x[which(V['adult_mov', , 140] > 2)], rep(sum(V[c('P_mov', 'P_stat'), , 140])/sum(V[c('adult_mov', 'adult_stat'), , 140]), length(which(V['adult_mov', , 140] > 2))), "l", col = 3, ylim = c(0, 1.5*10^-5), xlim = range(x))

lines(x[which((V['adult_stat', , 160] + V['adult_mov', , 160]) > 1)], rep(sum(V[c('P_mov', 'P_stat'), , 160])/sum(V[c('adult_mov', 'adult_stat'), , 160]), length(which((V['adult_stat', , 160] + V['adult_mov', , 160]) > 1))), col = 3, lty = 2)
