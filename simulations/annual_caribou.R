
library(animation)
library(gtools)
library(fBasics) # for Heaviside function
library(doParallel)
library(gplots)

source("simulations/popFunctions.R")

seasonCols <- c(
	"Spring" = "#9FCE62", 
	"Calving" = "#FFE548", 
	"postCalving" = "#FFA251",
	"Summer" = "#DC4A45",
	"Fall" = "#8E4D1D",
	"Winter" = "#618CCF")

###############################################################################################
# Parameter definition
###############################################################################################

# Define annual cycle
DOY <- c(1:365)

# What habitat (range) do caribou occupy throughout the year?
# See Fig. 11 of https://www.enr.gov.nt.ca/sites/enr/files/resources/draft_-_caribou_range_assessment_and_technical_information.pdf

# Used for temperature data from MERRA for the different ranges
bouRange <- c(
	rep("Winter", 109), # January - April 19
	rep("Spring", 43), # 20 April - June 1
	rep("Calving", 27), # 2 June - 28 Jun (includes post-calving migration)
	rep("Summer", 70), # 29 Jun - 6 September
	rep("Fall", 85), # 7 September - 30 November
	rep("Winter", 31) #  1 December - 31 December
)

if(length(bouRange) != length(DOY)) warning("\n\n\n\n\n***********\n\n\n\nbouRange not the right length. Seasons don't add up.\n\n\n\n\n***********\n\n\n\n")

#------------------------------------------------------------------------------
# Movement speed
#------------------------------------------------------------------------------
# What is the migration speed of caribou throughout the year

# 1) Based on data
# bouMove <- read.csv("parameterization/dailyMovementRates.csv")
# # How long is the migration based on the Gunn daily movement rates?
# sum(round(bouMove$movementRate)) # 3089 km

# Interpolated for seasons.
migSpeed <- 14

bouMove <- c(
	rep(0, (109 - 11)),# winter
	rep(migSpeed, 11 + 43), #spring
	rep(0, 15), # calving
	rep(migSpeed, 27-15), # post-calving migration
	rep(0, 70), # summer
	rep(migSpeed, 96),# Fall migration
	rep(0, 20)) # winter

###############################################################################################
# Space-time grid
###############################################################################################

    

#------------------------------------------------------------------------------
# Starting and stopping rates
#------------------------------------------------------------------------------

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


# # More complex
# # Starting:
# startMat <- matrix(0, nrow = length(DOY), ncol = length(x), dimnames = list(DOY, x))
# 
# # Build up to spring migration
# for(i in 1:5){
#  startMat[99 - c(5:1)[i], c(201:length(x), 1:200)] <- 1 - pnorm(q = x, mean = (2268 + 50 * seq(-1.96, 1.96, length.out = 5)[i]) - x[200], sd = 50)
#  I <- which(startMat[99 - c(5:1)[i], ] == 0)[1]
#  startMat[99 - c(5:1)[i], I:(which(x > 1000)[1])] <- pnorm(q = x[I:(which(x > 1000)[1])], mean = 800, sd = 50)
# }
# # Spring migration
# startMat[99:151, ] <- 1 
# 
# # plot(x, startMat[94, ], "n"); for(i in 1:5) lines(x, startMat[94+i, ], col = grey(i/7))
# 
# # Slow down at calving grounds
# for(i in 1:5){
# 	startMat[151 - i, ] <- 1 - pnorm(q = x, mean = (756 + 20 * seq(-1.96, 1.96, length.out = 5)[i]), sd = 50)
# }
# plot(x, startMat[151, ], "l", ylim = c(0,1)); for(i in 1:5) lines(x, startMat[151 - i, ], col = grey(i/7))
# abline(v = 756)
# 
# # 
# frontLine <- 2268+ 1.96*50
# plot(x, 1 - pnorm(q = x, mean = frontLine, sd = 50), "l", ylim = c(0,1))
# abline(v = frontLine)
# lines(x, startMat[99, ], col = 2)


###############################################################################################
# Finish setting up grid
###############################################################################################

# Distribution of herd among stages/sex from Boulanger et al. 2011 J Wild Man
bouDist1985 <- c(yearling = 88000, bull = 138000, calf = 176000, cow = 240000)


# Proportion of adults that are female for breeding purposes
propFemale <- bouDist1985['cow']/(bouDist1985['bull'] + bouDist1985['cow'])



###############################################################################################
# One-year simulation to look at host population changes
# No parasites
###############################################################################################
# Testing stop/starting
initBou1 <- rbind(
	calf_mov = rep(0, length(x)),
	yearling_mov = rep(0, length(x)),
	adult_mov = rep(0.1, length(x)),
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


grid1 <- gridSetup(y = 1, dt = 0.2)
V1 <- simBou(
	timeDat = grid1$timeDat, 
	x = grid1$x, 
	grid1$dx, 
	dt = 0.2,
	initBou = initBou1, 
	tempGrid = grid1$tempGrid,
	migSpeed = 14)

plot.timestep(V1, n = which(grid1$timeDat$DOY == as.numeric(strftime(as.Date("1985-06-15"), format = "%j")))[1])

#------------------------------------------------------------------------------
# Produce GIF for presentation
#------------------------------------------------------------------------------

Nrange <- range(V1[c('adult_mov', 'adult_stat'), , ])
x <- grid1$x
timeDat <- grid1$timeDat
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
ani.options(interval = intervals, ani.width = 600, ani.height = 300, ani.dev="png")

saveGIF(expr = {
	for(i in nInd){
		par(mfrow = c(1,1), mar = c(4,4,2,1), family = "Helvetica")
		# Moving
		plot(x, V['adult_mov', , 1], "n", ylim = Nrange, bty = "l", yaxt = "n", ylab = "", lwd = 1.5, xlab="")
		axis(side = 1, labels = FALSE)
		yTick <- pretty(Nrange)
		axis(side = 2, las = 1, at = yTick, labels = sprintf("%.1f", yTick/1000))
		points(x, rep(-0.1, length(x)), col = seasonCols[grid1$xRange], pch = 15, xpd = NA)
		
		lines(x, V['adult_mov', , i], col = 1, lwd = 1.5)#DOYcol[timeDat$DOY[i]])
		lines(x, V['adult_stat', , i], col = 1, lty = 2, lwd = 2)# DOYcol[timeDat$DOY[i]],lty = 2, lwd = 2)
		
		mtext(side = 3, line = 1, strftime(as.Date(paste("1985", timeDat$DOY[i], sep = "-"), format = "%Y-%j"), format = "%b %d"))
		
	}},
	, movie.name = "hostMovememt.gif", 
	img.name = "hostMovement")


###############################################################################################
# Base simulation starting with parasites
###############################################################################################

# Starting abundance of parasites in adults
# From Bathurst surveys, appprox. 
initP <- 100 

grid2 <- gridSetup(y = 3, dt = 0.2)

x <- grid2$x
initBou2 <- rbind(
	calf_mov = rep(0, length(x)),
	yearling_mov = rep(0, length(x)),
	adult_mov = rep(0.1, length(x)),
	L4_mov = rep(0, length(x)),
	L4A_mov = rep(0, length(x)),
	P_mov = rep(initP, length(x)),
	calf_stat = initDist(totPop = 10^3*bouDist1985['calf']/sum(bouDist1985), x),
	yearling_stat = initDist(totPop = 10^3*bouDist1985['yearling']/sum(bouDist1985), x),
	adult_stat = initDist(totPop = 10^3*c(bouDist1985['cow'] + bouDist1985['bull'])/sum(bouDist1985), x),
	L4_stat = rep(0, length(x)),
	L4A_stat = rep(0, length(x)),
	P_stat = rep(initP, length(x)),
	L0 = rep(0 * dx, length(x)),
	L3 = rep(0 * dx, length(x))
)

V2 <- simBou(
	timeDat = grid2$timeDat, 
	x = grid2$x, 
	grid2$dx, 
	dt = 0.2,
	initBou = initBou2, 
	tempGrid = grid2$tempGrid,
	migSpeed = 14)


plot.timestep(V2, n = 158/dt)
plot.timestep(V2, n = (365+158)/dt)
plot.timestep(V2, n = (365+249)/dt)
plot.timestep(V2, n = 4000)


# CHange in total population size over time
# CHange in total parasite and total host population over time
timeDat <- grid2$timeDat
V <- V2
Ntot <- apply(V['adult_mov', , ] + V['adult_stat', , ], 2, sum)

par(mfrow=c(3,1))
plot(as.Date(paste(1984+timeDat$year, timeDat$DOY, sep="-"), format = "%Y-%j"), Ntot, "l", main = "Total hosts over time")

plot(as.Date(paste(1984+timeDat$year, timeDat$DOY, sep="-"), format = "%Y-%j"), apply(V['P_stat', , ], 2, sum), "l", main = "Total P over time", ylim = c(0, 10^3))


plot(as.Date(paste(1984+timeDat$year, timeDat$DOY, sep="-"), format = "%Y-%j"), apply(V['L3', , ], 2, sum), "l", main = "Total L3 over time")
lines(as.Date(paste(1984+timeDat$year, timeDat$DOY, sep="-"), format = "%Y-%j"), apply(V['L0', , ], 2, sum), col = 4)


# CHanges over time in population at breeding grounds
timeDat <- grid2$timeDat

plot(timeDat)
















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


