y <- c(1.9560242163689843, 5.7644653721707515, 3.5834363643880156, 0.968812593714496, 0, 0)
x <- c(0, 1, 7, 30, 3*30, 6*30)

plot(x, y, bty = "l", xlab = "Days at -20*C", ylab = "epg")

out <- nls(y ~ a * exp(-u * x), data= data.frame(x, y), start = list(a = 10, u = 0.1))
summary(out)


xDummy <- 0:180
yDummy <- summary(out)$parameters[1,1] * exp(- summary(out)$parameters[2,1] * xDummy)

lines(xDummy, yDummy)
text(60, 3, )

###############################################################################
# Pandey 1972 studies
###############################################################################

eggSurv <- data.frame(
	day = rep(c(1/24, 2/24, 4/14, 6/24, 12/24, 1, 2, 4, 7*c(1,2,4,6,8,10,12,16,20,24,28,32,36,40,44,46,50,52)), 6),
	temp = rep(c(-10, 1, 4, 40, 45, 50), each = 26),
	surv = c(rep(100, 7), 74, 58, 36, 20, 3, 0, 0, rep(NA, 12),
					 rep(100, 9), c(95, 90, 66, 25, 22, 17, 12, 11, 11, 4, 5, 3, 2, 2, 2, 0, 0),
					 rep(100, 12 )
			
	)
)

par(mfrow = c(1,1))
#--------------------
# EGGS
y <- c(rep(100, 7), 74, 58, 36, 20, 3, 0, 0, rep(NA, 12))
x <- c(1/24, 2/24, 4/14, 6/24, 12/24, 1, 2, 4, 7*c(1,2,4,6,8,10,12,16,20,24,28,32,36,40,44,46,50,52))


plot(x, y, bty = "l", xlab = "Days at -10*C", ylab = "% survival", xlim = c(0, 100))

out <- nls(y ~ 100*exp(-u * x), data= data.frame(x,y), start = list(u = 0.1))
summary(out)


xDummy <- 0:180
yDummy <- 100 * exp(- summary(out)$parameters[1,1] * xDummy)

lines(xDummy, yDummy)
text(60, 60, expression(mu[egg] == 0.06759))
mtext(side = 3, line = 1, "-10*C from Pandey (1972)")

#-------------------
#L1
yL1 <- c(rep(100, 4), 90, 95, 80, 50, 20, 9, 3, 1, rep(0, 3))
xL1 <- c(15/60/24, 30/60/24, 1/24, 2/24, 4/24, 6/24, 1:7, 14, 21)

points(xL1, yL1, col = 2)

outL1 <- nls(yL1 ~ 100*exp(-u * xL1), data= data.frame(xL1,yL1), start = list(u = 0.1))
summary(outL1)
yDummyL1 <- 100 * exp(- summary(outL1)$parameters[1,1] * xDummy)

lines(xDummy, yDummyL1, col = 2)
text(60, 50, expression(mu[L1] == 0.44391), col = 2)

#-------------------
# L3

yL3 <- c(rep(100, 7), 80, 50, 25, 13, 6, 1, 0, 0)
xL3 <- c(c(0.5, 1, 2, 4, 8, 12)/24, 1, 2, 4, c(1, 2, 4, 6, 8, 10)*7)

points(xL3, yL3, col = 4)

outL3 <- nls(yL3 ~ 100*exp(-u * xL3), data= data.frame(xL3, yL3), start = list(u = 0.1))
summary(outL3)
yDummyL3 <- 100 * exp(- summary(outL3)$parameters[1,1] * xDummy)

lines(xDummy, yDummyL3, col = 4)
text(60, 40, expression(mu[L3] == 0.15324), col = 4)

###############################################################################
# Kafle et al. 2019 UP/VE freezing study
###############################################################################

datKafle <- read.csv("/Users/stephaniepeacock/Google\ Drive/Kutz\ Group/Pratap/freezing.csv")

head(datKafle)

# Data with initial = 0 don't add any information so we can remove:
datKafle <- datKafle[which(datKafle$initial > 0), ]


#---------------------------------------------------------------------------------------------------
# 1) DATA: Look at the proportion (+/- 95% SE) for UP & VE
#---------------------------------------------------------------------------------------------------

datKafle$combo <- paste(datKafle$species, datKafle$temp, datKafle$day, sep="-")

datKafle2 <- data.frame(
	initial = tapply(datKafle$initial, datKafle$combo, sum), 
	final = tapply(datKafle$final, datKafle$combo, sum), 
	temp = tapply(datKafle$temp, datKafle$combo, mean), 
	day = tapply(datKafle$day, datKafle$combo, mean), 
	species = as.factor(tapply(datKafle$species, datKafle$combo, unique))
	)

datKafle2$proportion <- datKafle2$final / datKafle2$initial

xDummy <- 1:180

outVE <- nls(proportion ~ a * exp(-u * day), data = subset(datKafle2, datKafle2$species == "VE"), start = list(a = 0.9, u = 0.001))
yDummyVE <- summary(outVE)$parameters[1,1] * exp(- summary(outVE)$parameters[2,1] * xDummy)

outUP <- nls(proportion ~ a * exp(-u * day), data = subset(datKafle2, datKafle2$species == "UP"), start = list(a = 0.9, u = 0.001))
yDummyUP <- summary(outUP)$parameters[1,1] * exp(- summary(outUP)$parameters[2,1] * xDummy)

plot(jitter(datKafle2$day), datKafle2$proportion, col = c(2,4)[as.numeric(datKafle2$species)], pch = as.numeric(as.factor(datKafle2$temp)), ylim = c(0.8, 1), xlab = "Day")
lines(xDummy, yDummyVE, col = 4)
text(100, 0.92, expression(mu[VE]==0.000398), col = 4)
lines(xDummy, yDummyUP, col = 2)
text(70, 0.85, expression(mu[UP]==0.000430), col = 2)

# Linear analysis with temperature

datKafle2$y <- log(datKafle2$proportion)

out.all <- lm(y ~ species + day:as.factor(temp), dat = datKafle2)
summary(out.all)

dumDat <- data.frame(
	day = rep(seq(0, 180, 2), 4), 
	temp = rep(unique(datKafle2$temp), each = 91), 
	VE = predict(out.all, newdata = data.frame(day = rep(seq(0, 180, 2), 4), temp = rep(unique(datKafle2$temp), each = 91), species = rep("VE", 91*4))), 
	UP = predict(out.all, newdata = data.frame(day = rep(seq(0, 180, 2), 4), temp = rep(unique(datKafle2$temp), each = 91), species = rep("UP", 91*4))))
				
lines(dumDat$day[dumDat$temp == -10], exp(dumDat$VE[dumDat$temp == -10]), col = 4)
lines(dumDat$day[dumDat$temp == -10], exp(dumDat$UP[dumDat$temp == -10]), col = 2)

lines(dumDat$day[dumDat$temp == -25], exp(dumDat$VE[dumDat$temp == -25]), col = 4, lty = 2)
lines(dumDat$day[dumDat$temp == -25], exp(dumDat$UP[dumDat$temp == -25]), col = 2, lty = 2)

lines(dumDat$day[dumDat$temp == -40], exp(dumDat$VE[dumDat$temp == -40]), col = 4, lty = 3)
lines(dumDat$day[dumDat$temp == -40], exp(dumDat$UP[dumDat$temp == -40]), col = 2, lty = 3)

lines(dumDat$day[dumDat$temp == -80], exp(dumDat$UP[dumDat$temp == -80]), col = 2, lty = 4)

# Across both species and all temps, freezing mortality of L1 seems to be 0.0004ish

###############################################################################
# deBruyn 2010 MSc thesis on freezing EPG recovery for mixed parasite eggs in 
# caribou feces frozen at -20*C
###############################################################################

datDeB <- data.frame(
	exp = c(rep("short", 3), rep("long", 6)),
	day = c(0, 1, 7, 7, c(1, 3, 6, 9, 12)*30),
	epg = c(24, 14, 8.1, 23, 20.4, 8.6, 2.7, 0.67, 0.67)
)

fitDat <- data.frame(y = datDeB$epg[datDeB$exp == "long"], x = datDeB$day[datDeB$exp == "long"])
outDeB <- nls(y ~ a * exp(u * x), data = fitDat, start = list(a = 25, u = -0.01))
summary(outDeB)

xDummy <- c(0:365)
yDummyDeB <- summary(outDeB)$parameters[1,1] * exp(summary(outDeB)$parameters[2,1] * xDummy)

plot(datDeB$day, datDeB$epg, col = as.numeric(as.factor(datDeB$exp)))
lines(xDummy, yDummyDeB)

summary(outDeB)$parameters[2,1]
# Egg moratlity at -20*C = 0.0117

###############################################################################
# Aleuy et al. 2020 Freeze survival of Marshallagia
###############################################################################

datAleuy <- data.frame(
	stage = rep(c(rep("Egg", 5), rep("L1", 5), rep("L3", 5)), 3),
	temp = rep(c(-9, -20, -35), each = 15),
	day = rep(c(3, 7, 14, 30, 90), 9),
	ppnSurv = c(
		c(0.974, 0.818, 0.978, 1, 0.848), c(0.10, 0.576, 0.513, 0.323, 0.0520), c(0.974, 1, 0.922, 0.970, 0.803),
		c(0.896, 0.918, 0.848, 0.890, 0.253), c(0.052, rep(0, 4)), c(0.896, 0.796, 0.327, 0.550, 0.550),
		c(0.342, rep(0, 4)), rep(0, 5), rep(0, 5)))


xDummy <- c(0:90)
mortRates <- matrix(NA, nrow = 3, ncol = 3, dimnames = list(c(-9, -20, -35), c("Egg", "L1", "L3")))

quartz(width = 6.3, height = 3.2, pointsize = 10)
par(mfrow = c(1,3), mar = c(4,4,2.5,1))
for(i in 1:3){ # temperatures
	Ti <- c(-9, -20, -35)[i]
	plot(datAleuy$day[datAleuy$temp == Ti], datAleuy$ppnSurv[datAleuy$temp == Ti], col = c(2, 1, 4)[as.numeric(as.factor(datAleuy$stage[datAleuy$temp == Ti]))], ylim = c(0, 1), lwd = 1.5, main = substitute(paste(Ti2, degree, "C", sep = ""), list(Ti2 = Ti)), xlab = "Expsoure (days)", ylab = "Survival (%)", las = 1)
	
	for(j in 1:3){# stages
		out.ij <- nls(y ~ exp(u * x), data = data.frame(y = datAleuy$ppnSurv[datAleuy$temp == Ti & datAleuy$stage == c("Egg", "L1", "L3")[j]], x = datAleuy$day[datAleuy$temp == Ti & datAleuy$stage == c("Egg", "L1", "L3")[j]]), start = list(u = -0.01))
		mortRates[i,j] <- summary(out.ij)$parameters[1,1]
		lines(xDummy, exp(mortRates[i,j] * xDummy), col = c(2, 1, 4)[j], lwd = 1.5)
		text(c(45, 55, 65)[j], exp(mortRates[i,j] * c(45, 55, 65)[j]), round(mortRates[i,j], 4), col = c(2, 1, 4)[j], pos = 3, xpd = NA)
	} # end j
	mtext(side = 3, line = 1, LETTERS[i], adj = 0)
} # end i
legend("topright", lwd = 1.5, col = c(2, 1, 4), c("Egg", "L1", "L3"))
