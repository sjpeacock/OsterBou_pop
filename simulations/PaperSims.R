source("simulations/bouSetup.R")
source("simulations/popFunctions.R")
library(gplots)
library(ggsci)

library(pals)
colPal <- pal_locuszoom()#scico(n = 7, palette = "berlin")
cols <- colPal(7)[c(1,2,3,5,4,6,7)]
colPal2 <- pal_locuszoom(alpha = 0.5)#scico(n = 7, palette = "berlin")
cols2 <- colPal2(7) [c(1,2,3,5,4,6,7)]
###############################################################################
# Set up grid and initial conditions
###############################################################################

# Starting abundance of parasites in adults
# From Bathurst surveys, appprox. 
initP <- 914

y <- 100
dt <- 1
stochTemp <- FALSE
startstop <- "agg"

source("simulations/calcGrid.R")

Bou0 <- rbind(
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

# Simulate a 100-year "burnin" for each level of transmission (3) and parasite inhibition (4)
V.init <- list(); length(V.init) <- 12; dim(V.init) <- c(3, 4)
initBou <- list(); length(initBou) <- 12; dim(initBou) <- c(3, 4)

# Want to store final variables but also mean parasite burdens in October 10
# (breedDOY - 240) which will influence year 1 pregnancy rates in subsequent simulations
n <- which(timeDat$DOY == 283 & timeDat$year == 100)

# Only keep the last 20 years, otherwise too much memory!
for(i in 1:3){ # for 3 levels of transmission
	for(p in 1:4){ # for 3 levels of parasite inhibition
		
		V.dum <- simBou(
			initBou = Bou0, 
			tempGrid = tempGrid,
			ppnInhibit = pArrestCalc(p),
			transmission = c("low", "base", "high")[i])
		
		V.init[[i, p]] <- V.dum[, , which(timeDat$year > 80)]
		
		initBou[[i, p]] <- list(
			var = V.dum[, , dim(V.dum)[3]],
			OctP = sum(V.dum[c('P_stat', 'P_mov'), , n]) / sum(V.dum[c('adult_stat', 'adult_mov'), , n])
		)

		}}

# Check
i <- 3; p <- 1
par(mfrow = c(3,1))
for(j in 1:3){
	J <- grep(c("adult_", "P_", "L3")[j], dimnames(initBou[[i,p]]$var)[[1]])
	plot(x, initBou[[i,p]]$var[J[1], ], "l", ylim = range(initBou[[i,p]]$var[J, ]), main = c("adult_", "P_", "L3")[j])
	if(length(J) > 1) lines(x, initBou[[i,p]]$var[J[2], ], col = 4)
}

# saveRDS(initBou, file = "simulations/output/initBou.rds")
# initBou <- readRDS("simulations/output/initBou.rds")


###############################################################################
# Look at annual dynamics in year 30 of each simulation
###############################################################################
# Note that the broad patterns don't change as we run the simulations for longer,
# just the scale of the 100% inhibition parasite burdens is much higher and
# thus comparisons are more difficult.
# Similarly the annual patterns are the same in the peak and trough of the longer
# term cycles, it's just the scale that shifts.
y <- 30
dt <- 1
stochTemp <- FALSE
startstop <- "agg"

source("simulations/calcGrid.R")

V.annual <- list(); length(V.annual) <- 24; dim(V.annual) <- c(2, 3, 4)
for(m in 1:2){ # for migratory and non-migratory simulations
	for(i in 1:3){ # for 3 levels of transmission
		for(p in 1:4){ # for 4 levels of parasite inhibition
			V.dum <- simBou(
				initBou = Bou0, 
				tempGrid = tempGrid,
				ppnInhibit = pArrestCalc(p),
				migSpeed = c(14, 0)[m],
				transmission = c("low", "base", "high")[i])
			
			V.annual[[m, i, p]] <- V.dum[, , which(timeDat$year == 30)]
		}}
}

#------------------------------------------------------------------------------
# Summarize annual dynamics
#------------------------------------------------------------------------------

oneYear <- array(NA, dim = c(2, 3, 4, 5, 365))
# dimensions = beta (3) * inhibition (4) * variable * day

for(m in 1:2){ # for migratory and resident
	for(i in 1:3){
		for(p in 1:4){ # for each level of ppnInhibit (1, 0.5, 0)
			for(j in 1:365){
				oneYear[m, i, p, 1, j] <- sum(V.annual[[m, i, p]][c('P_mov', 'P_stat'), , j]) / sum(V.annual[[m, i, p]][c('adult_mov', 'adult_stat'), , j])
				oneYear[m, i, p, 2, j] <- sum(V.annual[[m, i, p]][c('L4A_mov', 'L4A_stat'), , j]) / sum(V.annual[[m, i, p]][c('adult_mov', 'adult_stat'), , j])
				oneYear[m, i, p, 3, j] <- sum(V.annual[[m, i, p]][c('L4_mov', 'L4_stat'), , j]) / sum(V.annual[[m, i, p]][c('adult_mov', 'adult_stat'), , j])
				oneYear[m, i, p, 4, j] <- sum(V.annual[[m, i, p]]['L0', , j])
				oneYear[m, i, p, 5, j] <- sum(V.annual[[m, i, p]]['L3', , j])
			} # end day j
			
		}}
} # end m

# saveRDS(oneYear, file = "simulations/output/oneYear30.rds")
# oneYear <- readRDS("simulations/output/oneYear30.rds")

#------------------------------------------------------------------------------
# Plot
#------------------------------------------------------------------------------

# Plot all levels of transmission or just the base level?
allBeta <- FALSE

# Plot resident populations in dashed
plotRes <- TRUE

if(allBeta == TRUE){
	I <- c(1:3) 
	quartz(width = 6.3, height = 6, pointsize = 10)
	layout(mat = cbind(
		matrix(c(rep(1:3, each = 3), 16, rep(4:5, each = 3)), ncol = 1),
		5 + matrix(c(rep(1:3, each = 3), 11, rep(4:5, each = 3)), ncol = 1),
		10 + matrix(c(rep(1:3, each = 3), 6, rep(4:5, each = 3)), ncol = 1)))
	
	} else{
		I <- 2
		quartz(width = 4, height = 6, pointsize = 10)
		layout(mat =	matrix(c(rep(1:3, each = 3), 6, rep(4:5, each = 3)), ncol = 1))
	} 

xDate <- as.Date(paste(30, c(1:365), sep = "-"), format = "%y-%j") 
yPerc <- 0.2

par(oma = c(4, 5, 4, 1))
par(mar = c(0,4,0,1))
# [m, i, p, 5, j]

for(i in I){
	for(s in 1:5){ # For each stage
		if(allBeta == FALSE){
			ylim.s <- extendrange(oneYear[1, i, , s, ], f = c(0, yPerc))
		} else {
			ylim.s <- extendrange(oneYear[1, , , s, ], f = c(0, yPerc))
		}
		
		if(s == 3|s == 5){
			plot(xDate, oneYear[1, i, 1, 1, ], "n", ylim = ylim.s, yaxt = "n", ylab = "", xaxs = "i")
		}else{
			plot(xDate, oneYear[1, i, 1, 1, ], "n", ylim = ylim.s, yaxt = "n", ylab = "", xaxt = "n", xaxs = "i")
		}
		
		for(p in 1:4){
			lines(xDate, oneYear[1, i, p, s, ], col = cols[p], lwd = 1.5)
			if(plotRes == TRUE) lines(xDate, oneYear[2, i, p, s, ], col = cols[p], lty = 2)
		}
		
		if(i == 1 | allBeta == FALSE) mtext(side = 2, line = 1, c("Adult\nparasites", "Arrested\nlarvae", "Developing\nlarvae", "Pre-infective\nlarvae", "Infective\nlarvae")[s])
		if((i == 1 | allBeta == FALSE) & s == 1){
		legend("topleft", col = cols[c(1:4)], title = "Inhibition", lwd = 1.5, legend = c("100%", "50%", "0%", "variable"), bty = "n")
		}
	
	} #  end s
} #end i


###############################################################################
# Run simulations for subsequent 100 years
###############################################################################

y <- 100
dt <- 1
stochTemp <- FALSE
startstop <- "agg"

source("simulations/calcGrid.R")

# Store AnnualP and Annual N
# Store full dynamics for the last 20 years, so peak and trough can be further 
# explored

# Simulate a 100-year simulation for each 
# r: climate scenario (NULL, RCP2.6, RCP8.5), 
# i: level of transmission (low, base, high), and 
# p: parasite inhibition (1, 0.5, 0, variable).

V <- list(); length(V) <- 36; dim(V) <- c(3, 4, 3)
annualSumm <- list(); length(annualSumm) <- 36; dim(annualSumm) <- c(3, 4, 3)

# Only keep the last 20 years, otherwise too much memory!

for(r in 1:3){ # for three climate scenarios: past, RCP 2.6, and RCP 8.5
	
	if(r == 1){
		tempGrid <- predict.temp(timeDat = timeDat, stoch = FALSE)
	} else if(r == 2){
		tempGrid <- predict.temp(timeDat = timeDat, stoch = FALSE, climateScenario = "rcp26")
	} else if(r == 3){
		tempGrid <- predict.temp(timeDat = timeDat, stoch = FALSE, climateScenario = "rcp85")
	}
	
	for(i in 1:3){ # for 3 levels of transmission
		for(p in 1:4){ # for 3 levels of parasite inhibition
			V.dum <- simBou(
				initBou = initBou[[i, p]]$var, 
				tempGrid = tempGrid,
				ppnInhibit = pArrestCalc(p),
				transmission = c("low", "base", "high")[i], 
				OctP = initBou[[i, p]]$OctP)
			
			# Store last 20 years of each simulation
			V[[i, p, r]] <- V.dum[, , which(timeDat$year > 80)]
			
			# Calculate AnnualP and AnnualN
			annualSumm[[i, p, r]] <- matrix(NA, nrow = y, ncol = 2, dimnames = list(c(1:y), c("annualP", "annualN")))
			
			for(n in 1:y){
				V.n <- V.dum[, , which(timeDat$year == n)]
				Pavg <- numeric(365)
				Ntot <- numeric(365)
				for(j in 1:365){
					Pavg[j] <- sum(V.n[c('P_mov', 'P_stat'), , j]) / sum(V.n[c('adult_mov', 'adult_stat'), , j])
					Ntot[j] <- sum(V.n[c('adult_mov', 'adult_stat'), , j])
				} # end day j
				
				annualSumm[[i, p, r]][n, 1] <- sum(Pavg, na.rm = TRUE)
				annualSumm[[i, p, r]][n, 2] <- sum(V.n[c('adult_mov', 'yearling_mov', "calf_mov", "adult_stat", "yearling_stat", "calf_stat"), , 159]) * dx
			} # end n
			
		} # end r climate scenario
			} # end P
} # end i
  

# saveRDS(annualSumm, file = "annualSumm_GT.rds")
# saveRDS(V, file = "V_GT.rds")

###############################################################################
# Summarize cycles in average parasite burden and host abundance over 100 years
###############################################################################

# all climate scenarios have the same y-axis scale for comparison
ylims <- list(
	P = list(
		lowTrans <- range(annualSumm[[1, 3, 1]][, 1]),
		baseTrans <- range(annualSumm[[2, 3, 1]][, 1]),
		highTrans <- range(annualSumm[[3, 3, 1]][, 1])
	),
	N = list(
		lowTrans <- range(annualSumm[[1, 3, 1]][, 2]),
		baseTrans <- range(annualSumm[[2, 3, 1]][, 2]),
		highTrans <- range(annualSumm[[3, 3, 1]][, 2])
	))
for(i in 1:3){
	for(r in 1:3){
		for(p in 2:4){
			ylims[[1]][[i]] <- range(c(ylims[[1]][[i]], annualSumm[[i, p, r]][, 1]))
			ylims[[2]][[i]] <- range(c(ylims[[2]][[i]], annualSumm[[i, p, r]][, 2]))
		}
	}
}



pdf(file = "figures/supplement/annualP_allScenarios_zoomed_GT.pdf", width = 8, height = 8, pointsize = 10)
par(mfrow = c(3,3), mar = c(1,1,1,1), oma = c(4,5,2,2))
for(i in 1:3){ # for each of three levels of transmission
	for(r in 1:3){ # for each climate scenario
		
		# Plot all three inhibition on the same figure using cols[1-3] (red = 100% inhibition)
		plot(1:y, annualSumm[[i, 1, r]][, 1], "n", ylim = ylims[[1]][[i]], yaxt = "n", xaxt = "n")
		for(p in 1:4)	lines(1:y, annualSumm[[i, p, r]][, 1], "o", pch = 19, cex = 0.3, col = cols[p], lwd = 1.2)
		if(r == 1) axis(side = 2, las = 1) else axis(side = 2, labels = FALSE)
		if(i == 3) axis(side = 1) else axis(side = 1, labels = FALSE)
		
		# # Plot host population over top in dashed line
		# par(new = TRUE)
		# plot(1:y, annualSumm[[i, 1, r]][, 2], "n", ylim = ylims[[2]][[i]], yaxt = "n", xaxt = "n")
		# for(p in 1:3)	lines(1:y, annualSumm[[i, p, r]][, 2], col = cols[p], lwd = 0.5)
		# 
		if(i == 1) mtext(side = 3, line = 1, c("Current", "RCP 2.6", "RCP 8.5")[r])
		
		if(r == 3) mtext(side = 4, c(expression(paste("Low ", beta)), expression(paste("Base ", beta)), expression(paste("High ", beta)))[i], line = 1)
		
		# if(r == 3 & i == 1) legend("right", col = cols[1:3], title = "Inhibition", legend = c("100%", "50%", "0%"), bty = "n", lwd = 1.2)
	}
	}
mtext(side = 1, "Year in simulation", outer = TRUE, line = 2)
mtext(side = 2, "Average annual parasite pressure", outer = TRUE, line = 3.5)
dev.off()

#------------------------------------------------------------------------------
# Plot of base transmission under current conditions for main text
#------------------------------------------------------------------------------
i <- 2 # Transmission rate scenario (2 = base)
r <- 1 # Climate scenario (1 = current)

ylimCycles <- array(NA, dim = c(4,2,2))
for(p in 1:4){
	ylimCycles[p, 1, ] <- extendrange(annualSumm[[i, p, r]][60:y, 1], f = 0.15)
	ylimCycles[p, 2, ] <- extendrange(annualSumm[[i, p, r]][60:y, 2], f = 0.15)
}
# ylimCycles[4,1,] <- c(82, 83) 

# quartz(width = 6.8, height = 5, pointsize = 10)
pdf(file = "figures/popCycles_GT.pdf", width = 6.8, height = 5, pointsize = 10)
par(mfrow = c(2, 2), mar = c(1,4,1,4), oma = c(3,1.5,1,4))

for(p in 1:4){
	plot(60:y, annualSumm[[i, p, r]][60:y, 1], "o", pch = c(21,22,24,25)[p], col = cols[p], bg = cols[p], cex = 0.8, lwd = 1.2, xaxt = "n", ylab = "", xlab = "", yaxt = "n", ylim = ylimCycles[p,1,])
	if(p > 2) axis(side = 1, at = seq(60, 100, 20), labels = seq(160, 200, 20)) else axis(side = 1, at = seq(60, 100, 20), labels = FALSE)
	axis(side = 2, at = pretty(annualSumm[[i, p, r]][60:y, 1]), labels = pretty(annualSumm[[i, p, r]][60:y, 1]*10^-6), las = 1)
	if(p == 1){
		arrows(x0 = 76, x1 = 94, y0 = max(annualSumm[[i, p, r]][60:y, 1]), y1 = max(annualSumm[[i, p, r]][60:y, 1]), length = 0.06, cod = 3, lwd = 1.2)
		text(85, max(annualSumm[[i, p, r]][60:y, 1]), pos = 1, "18 years")
	} else if(p == 2){
		arrows(x0 = 71, x1 = 81, y0 = max(annualSumm[[i, p, r]][60:y, 1]), y1 = max(annualSumm[[i, p, r]][60:y, 1]), length = 0.06, cod = 3, lwd = 1.2)
		text(76, max(annualSumm[[i, p, r]][60:y, 1]), pos = 1, "10 years")
	} else if(p == 3){
		arrows(x0 = 86, x1 = 95, y0 = max(annualSumm[[i, p, r]][60:y, 1]), y1 = max(annualSumm[[i, p, r]][60:y, 1]), length = 0.06, cod = 3, lwd = 1.2)
		text(90.5, max(annualSumm[[i, p, r]][60:y, 1]), pos = 1, "9 years")
	}
	
	par(new = TRUE)
	plot(60:y, annualSumm[[i, p, r]][60:y, 2], "o", pch = c(21,22,24,25)[p], col = cols[p], bg = "white", cex = 0.8, lty = 3, yaxt = "n", xaxt = "n", ylim = ylimCycles[p, 2 ,], ylab = "", xlab = "")
	axis(side = 4, at = pretty(annualSumm[[i, p, r]][60:y, 2]), labels = pretty(annualSumm[[i, p, r]][60:y, 2]*10^-6), las = 1)
	if(p == 1){
		arrows(x0 = 80, x1 = 84, y0 = min(annualSumm[[i, p, r]][60:y, 2]), y1 = min(annualSumm[[i, p, r]][60:y, 2]), code = 3, length = 0.06, lwd = 1.2)
		text(82, min(annualSumm[[i, p, r]][60:y, 2]),pos =1, "4 years")
	} else if(p == 2){
		arrows(x0 = 83, x1 = 86, y0 = min(annualSumm[[i, p, r]][60:y, 2]), y1 = min(annualSumm[[i, p, r]][60:y, 2]), code = 3, length = 0.06, lwd = 1.2)
		text(84.5, min(annualSumm[[i, p, r]][60:y, 2]),pos =1, "3 years")
	} else if(p == 3){
		arrows(x0 = 88, x1 = 90, y0 = min(annualSumm[[i, p, r]][60:y, 2]), y1 = min(annualSumm[[i, p, r]][60:y, 2]), code = 3, length = 0.06, lwd = 1.2)
		text(89, min(annualSumm[[i, p, r]][60:y, 2]), pos =1, "2 years")
	}
	
if(p == 4) legend("bottomright", pch = c(25, 25), pt.bg = c(cols[4], "white"), lty = c(1, 3), lwd = c(1.2, 1), legend = c("Parasite pressure", "Total hosts"), col = cols[4], bty = "n")

	mtext(side =3, adj = 0, line = -1.5, paste(" ", LETTERS[p]))
	
	
} # end p

mtext(side = 2, outer = TRUE, expression(paste("Annual parasite pressure (parasite", {}%*%{}, "days (host)", {}^-1 %*% 10^-6, ")")), line = )
mtext(side =4, outer = TRUE, expression(paste("Total host population size (", {}%*%10^-6, ")")), line = 2)
mtext(side = 1, outer = TRUE, "Year in simulation", line = 2)

dev.off()
###############################################################################
# For climate scenarios, summarize 
# (1) average parasite burdens in last 80 years
# (2) average host population in the last 80 years.
###############################################################################
boot <- function(x){
	bootX <- apply(matrix(sample(x, length(x) * 1000, replace = TRUE), nrow = 1000), 1, mean)
	return(c(quantile(bootX, 0.025), mean(bootX), quantile(bootX, 0.975)))
}

avg80 <- array(NA, dim = c(3, 4, 3, 2, 3))
for(i in 1:3){
	for(r in 1:3){
		for(p in 1:4){
			for(m in 1:2){
				avg80[i, p, r, m, ] <- boot(annualSumm[[i, p, r]][21:100, m])
			}
		}
	}
}

# Plot as a percentage of the mean in current scenarios
# quartz(width = 4, height = 4, pointsize = 10)
ylims2 <- matrix(c(0.9, 1.15, 0.85, 1.1), 2, 2)

for(i in 1:3){
	if(i == 2) pdf(file = "figures/climateChangeEffects.pdf", width = 3.2, height = 4, pointsize = 10)	
	if(i == 1) pdf(file = "figures/supplement/climateChangeEffects_lowBeta.pdf", width = 3.2, height = 4, pointsize = 10)	
	if(i == 3) pdf(file = "figures/supplement/climateChangeEffects_highBeta.pdf", width = 3.2, height = 4, pointsize = 10)	

	par(mfrow = c(2,1), mar = c(2,5,2,1), oma = c(2, 0, 0, 6))
	for(m in 1:2){
		plot(1:3, avg80[i, 1, , m, 2], "n", ylim = ylims2[, m], xlim = c(0.5, 3.5), xaxs = "i", xaxt = "n", las = 1, xlab = "", ylab = "")
		axis(side = 1, at = c(0.5, 1.5, 2.5, 3.5), labels = FALSE)
		axis(side = 1, at = c(1, 2, 3), labels = c("Current", "RCP 2.6", "RCP 8.5"), tck = 0)
		abline(h = 1)
		abline(v = c(1.5, 2.5))
		for(p in 1:4){
			segments(x0 = 1:3 + c(-0.3, -0.1, 0.1, 0.3)[p], x1 = 1:3 + c(-0.3, -0.1, 0.1, 0.3)[p], y0 = avg80[i, p, , m, 1]/avg80[i, p, 1, m, 2], y1 = avg80[i, p, , m, 3]/avg80[i, p, 1, m, 2], col = cols2[p], lwd = 2)
			points(1:3 + c(-0.3, -0.1, 0.1, 0.3)[p], avg80[i, p, , m, 2]/avg80[i, p, 1, m, 2], pch = c(21,22,24,25)[p], col = cols[p], bg = cols2[p])
		}
		mtext(side = 3, adj = 0, line = 0, LETTERS[m])
		mtext(side = 2, c("% change in parasites", "% change in host pop'n")[m], line = 3)
	}
	mtext(side = 1, "Climate scenario", line = 3)
	legend(3.8, 1.5, pch = c(21,22,24,25), col = cols[1:4], pt.bg = cols2[1:4], title = "Arrested\ndevelopment", legend = c("100%", "50%", "0%", "variable"), bty = "n", cex = 0.8, pt.cex = 1.0, xpd= NA)

	dev.off()
}
#------------------------------------------------------------------------------
# Does the timing of infection differ and is that why we see declines in host
# abundance with no significant increase in parasite populations?
#------------------------------------------------------------------------------
# Look at relative summary in last year of simulation

year2plot <- 100 #last 80 years in simulation

annual100 <- array(NA, dim = c(3, 3, 3, 6, 365))
for(i in 1:3){
	for(p in 1:3){ # for each level of ppnInhibit (1, 0.5, 0)
		for(r in 1:3){ # for migratory and resident
			V.ipr <- V[[i,p,r]][, , which(timeDat$year[timeDat$year > 80] == year2plot)]
			for(j in 1:365){
				annual100[i, p, r, 1, j] <- sum(V.ipr[c('P_mov', 'P_stat'), , j]) / sum(V.ipr[c('adult_mov', 'adult_stat'), , j])
				annual100[i, p, r, 2, j] <- sum(V.ipr[c('L4A_mov', 'L4A_stat'), , j]) / sum(V.ipr[c('adult_mov', 'adult_stat'), , j])
				annual100[i, p, r, 3, j] <- sum(V.ipr[c('L4_mov', 'L4_stat'), , j]) / sum(V.ipr[c('adult_mov', 'adult_stat'), , j])
				annual100[i, p, r, 4, j] <- sum(V.ipr['L0', , j])
				annual100[i, p, r, 5, j] <- sum(V.ipr['L3', , j])
				annual100[i, p, r, 6, j] <- sum(V.ipr[c('adult_mov', 'adult_stat', 'yearling_mov', 'yearling_stat', 'calf_mov', "calf_stat"), , j])
			} # end day j
			
			
		}}
} # end m


#-----
# Plot 
#-----
xDate <- as.Date(paste(2100, c(1:365), sep = "-"), format = "%Y-%j") 

i <- 2
p <- 1
par(mfrow = c(2,1))
for(s in c(1,6)){
	plot(xDate, annual100[i, p, 1, s, ], "l", col= 1, ylim = extendrange(annual100[i, p, 1, s, ], f = 0.05))
	par(new = TRUE)
	plot(xDate, annual100[i, p, 2, s, ], "l", col= cols[4], ylim = extendrange(annual100[i, p, 2, s, ], f = 0.05))
	par(new = TRUE)
	plot(xDate, annual100[i, p, 3, s, ], "l", col= cols[1], ylim = extendrange(annual100[i, p, 3, s, ], f = 0.05))
	abline(v = as.Date(paste("2100", breedDOY - 240 + 365, sep = "-"), format = "%Y-%j"))
	# points(rep(as.Date(paste("2100", breedDOY - 240 + 365, sep = "-"), format = "%Y-%j"), 2))
}


par(mfrow = c(2,1))
s <- 1
plot(xDate, annual100[i, p, 1, s, ], "l", col= 1, ylim = extendrange(annual100[i, p, 1, s, ], f = 0.05))
par(new = TRUE)
plot(xDate, annual100[i, p, 2, s, ], "l", col= cols[4], ylim = extendrange(annual100[i, p, 2, s, ], f = 0.05))
par(new = TRUE)
plot(xDate, annual100[i, p, 3, s, ], "l", col= cols[1], ylim = extendrange(annual100[i, p, 3, s, ], f = 0.05))
abline(v = as.Date(paste("2100", breedDOY - 240 + 365, sep = "-"), format = "%Y-%j"))


preg <- cbind(seq(0, 60000, 100), 0.8 - 1/(1 + exp(7.025 - 0.000328*seq(0, 60000, 100))))

p <- 2


for(s in c(1,6)){
	plot(xDate, annual100[i, p, 1, s, ], "l", col= 1, ylim = extendrange(annual100[i, p, 1:3, s, ], f = 0.05))
	for(r in 2:3) lines(xDate, annual100[i, p, r, s, ], "l", col= cols[c(4,1)[r-1]])
	abline(v = as.Date(paste("2100", breedDOY - 240 + 365, sep = "-"), format = "%Y-%j"), lty = 3)
	abline(v = as.Date(paste("2100", breedDOY, sep = "-"), format = "%Y-%j"), lty = 3)
if(s == 1) lines(preg[, 2]/max(preg[,2]) * 180))
}










