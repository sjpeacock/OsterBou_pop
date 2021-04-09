library(gplots)
par1 <- read.csv("MTE_dataAnalysis/output/bestModelParameters.csv")
par2 <- read.csv("MTE_dataAnalysis/output/bestModelParameters_Isigma.csv")

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

plotTpolys <- function(){
	u <- par('usr')
	polygon(x = c(-52.46723,  24.38447, 24.38447, -52.46723 ), y = u[c(3,3,4,4)], border = NA, col = "#00000030")
	polygon(x = c(-33.52071,  15.62382, 15.62382, -33.52071), y = u[c(3,3,4,4)], border = NA, col = "#00000030")
}
	


temp <- seq(5, 35, 5)
temp.all <- seq(-50, 40, 0.2)
quartz(width = 6.3, height = 8, pointsize = 10)
par(mfrow = c(2, 2), mar = c(4,4,2,1))

#mu0
plotCI(temp, convNum(par1[1:7, 'mu0'], type = "m"), li = convNum(par1[1:7, 'mu0'], type = "l"), ui = convNum(par1[1:7, 'mu0'], type = "u"), gap = 0.3,  sfrac= 0.008, bty = "l", las = 1, ylab = "", xlab = "", main = expression(mu[0]), xlim = c(-50, 40), ylim = c(0, 1))

plotTpolys()
plotCI(temp, convNum(par2[1:7, 'mu0'], type = "m"), li = convNum(par2[1:7, 'mu0'], type = "l"), ui = convNum(par2[1:7, 'mu0'], type = "u"), col = 2, gap = 0, cex=0.7, add=TRUE, pch = 21, pt.bg = "white")
lines(temp.all, predict.MTE(params = c(a = convNum(par1[9,3]), E = convNum(par1[10,3]), Eh = NA, Th = NA, z = 1), temp = temp.all))
lines(temp.all, predict.MTE(params = c(a = convNum(par2[8,3]), E = convNum(par2[9,3]), Eh = NA, Th = NA, z = 1), temp = temp.all), col = 2)

lines(temp.all, predict.MTE(params = c(a = 0.068, E = 0.884, Eh =NA, Th = NA, El = 3.227, Tl = -6.282,  z = 1), temp = temp.all), lty = 3, lwd = 2)

lines(temp.all, predict.MTE(params = c(a = 0.056, E = 0.65, Eh =3.25, Th = 37.5, El = 3.25, Tl = -2.5,  z = 1), temp = temp.all), col = 4, lwd = 1.5)
lines(temp.all, predict.MTE(params = c(a = 0.056, E = 0.65, Eh =NA, Th = NA, El = NA, Tl = NA,  z = 1), temp = temp.all), col = 4, lwd = 1.5, lty = 2)

# Bryanne's freezing of eggs (BUT none developed...)
points(-20, 0.03630041, pch = 1, lwd = 2, col = 7)

# debruyn 2010 epg recovery
points(-20, 0.0117, pch = 19, lwd = 2, col = 7)

# Pandey 1972 egg and L1
points(-10, 0.06759, pch = 1, col = 6, lwd = 2)
points(-10, 0.44391, pch = 3, col = 6, lwd = 2)

# Ale's marshallagia study
points(-9, 0.001783, pch = 1, col = 5, xpd = NA, lwd = 2) #egg
points(-20, 0.01132, pch = 1, col = 5, xpd = NA, lwd = 2) #egg
points(-35, 0.3904, pch = 1, col = 5, xpd = NA, lwd = 2) #egg

points(-9, 0.0855, pch = 3, col = 5, xpd = NA, lwd = 2) #L1
points(-20, 0.9857, pch = 3, col = 5, xpd = NA, lwd = 2) #L1
# points(-35, 3, pch = 3, col = 5, xpd = NA, lwd = 2) #L1

plot(1,1,"n", bty="n", xaxt="n", yaxt="n", xlab = "", ylab = "")
legend("topleft", col = c(1,2,4,4,1,7,7,6,5,1,1), lty = c(1,1,1,2,3,rep(NA, 6)), lwd = c(1,1,1.5,2,2, rep(2, 6)), c(expression(paste("Constant ", sigma)), expression(paste("Different ", sigma)), "Molnar 2013 SSUL", "Molnar 2013 Arr.", "Proposed", "Hoar (2012), O. grue eggs", "deBruyn (2010), O. grue+ eggs", "Pandey (1972), O. ostertagi", "Aleuy et al. (2020), Marshallagia", "Egg", "L1"), pch = c(21, 21, NA, NA, NA, 1, 19, rep(15, 2), 1, 3),  ncol = 1, xpd = NA, bty = "n", pt.bg = "white")

#------------------------------------------------------------------------------
# mu 3
#------------------------------------------------------------------------------

plotCI(temp, convNum(par1[1:7, 'mu1'], type = "m"), li = convNum(par1[1:7, 'mu1'], type = "l"), ui = convNum(par1[1:7, 'mu1'], type = "u"), gap = 0.3,  sfrac= 0.008, bty = "l", las = 1, ylab = "", xlab = "", main = expression(mu[3]), xlim = range(temp.all), ylim = c(0, 0.2))
plotTpolys()
plotCI(temp, convNum(par2[1:7, 'mu1'], type = "m"), li = convNum(par2[1:7, 'mu1'], type = "l"), ui = convNum(par2[1:7, 'mu1'], type = "u"), col = 2, gap = 0, cex=0.7, add=TRUE, pch = 21, pt.bg = "white")
lines(temp.all, rep(convNum(par1[9, 'mu1']), length(temp.all)))
lines(temp.all, rep(convNum(par2[8, 'mu1']), length(temp.all)), col = 2)

lines(temp.all, predict.MTE(params = c(a = 0.0211, E = 0.208, Eh = 3.5543, Th = 27.6, El = 2.16, Tl = -18.8,  z = 1), temp = temp.all), lty = 3, lwd = 2)
# lines(temp.all, predict.MTE(params = c(a = 0.0211, E = 0.208, Eh = 3.5543, Th = 27.6, El = 3.249, Tl = -28.267,  z = 1), temp = temp.all), lty = 3)

lines(temp.all, predict.MTE(params = c(a = 0.056, E = 0.65, Eh =3.25, Th = 37.5, El = 3.25, Tl = -2.5,  z = 1), temp = temp.all), col = 4, lwd = 1.5)
lines(temp.all, predict.MTE(params = c(a = 0.056, E = 0.65, Eh =NA, Th = NA, El = NA, Tl = NA,  z = 1), temp = temp.all), col = 4, lwd = 1.5, lty = 2)

# Pandey 1972 L3
points(-10, 0.15324, pch = 3, col = 6, lwd = 2)

# Ale's marshallagia study
points(-9, 0.002344, pch = 3, col = 5, xpd = NA, lwd = 2) #L3
points(-20, 0.01762, pch = 3, col = 5, xpd = NA, lwd = 2) #L3
points(-35, 3, pch = 3, col = 5, xpd = NA, lwd = 2) #L3

# Anja's marshallagia study (Carlsson et al. 2013)
points(-48, 4, pch = 3, col = 7, lwd = 2)

lines(temp.all, predict.MTE(params = c(a = 0.068, E = 0.884, Eh =NA, Th = NA, El = 3.227, Tl = -6.282,  z = 1), temp = temp.all), lty = 3, lwd = 2, col = "white")
#------------------------------------------------------------------------------
# rho0
#------------------------------------------------------------------------------
plotCI(temp, convNum(par1[1:7, 'rho'], type = "m"), li = convNum(par1[1:7, 'rho'], type = "l"), ui = convNum(par1[1:7, 'rho'], type = "u"), gap = 0.3,  sfrac= 0.008, bty = "l", las = 1, ylab = "", xlab = "", main = expression(rho[0]), xlim = range(temp.all), ylim = c(0, 0.12))
plotTpolys()
plotCI(temp, convNum(par2[1:7, 'rho'], type = "m"), li = convNum(par2[1:7, 'rho'], type = "l"), ui = convNum(par2[1:7, 'rho'], type = "u"), col = 2, gap = 0, cex=0.7, add=TRUE, pch = 21, pt.bg = "white")
lines(temp.all, predict.MTE(params = c(a = convNum(par1[9, 'rho']), E = convNum(par1[10, 'rho']), Eh = convNum(par1[11, 'rho']), Th = 30.568, z = -1), temp = temp.all))
lines(temp.all, predict.MTE(params = c(a = convNum(par2[8,'rho']), E = convNum(par2[9,'rho']), Eh = NA, Th = NA, z = -1), temp = temp.all), col = 2)

lines(temp.all, predict.MTE(params = c(a = 0.032, E = 0.686, Eh = 7.957, Th = 30.568, El = NA, Tl = NA,  z = -1), temp = temp.all), lty = 3, lwd = 2)

lines(temp.all, predict.MTE(params = c(a = (1/29.6), E = 0.65, Eh =3.25, Th = 32.5, El = 3.25, Tl = 2.5,  z = -1), temp = temp.all), col = 4, lwd = 1.5)
lines(temp.all, predict.MTE(params = c(a = (1/29.6), E = 0.65, Eh =NA, Th = NA, El = NA, Tl = NA,  z = -1), temp = temp.all), col = 4, lwd = 1.5, lty = 2)

mtext(side = 1, outer=TRUE, "Temperature (*C)")
# # sigma
# plotCI(temp, convNum(par2[1:7, 'sigma1'], type = "m"), li = convNum(par2[1:7, 'sigma1'], type = "l"), ui = convNum(par2[1:7, 'sigma1'], type = "u"), gap = 0.3,  sfrac= 0.008, bty = "l", las = 1, ylab = "", xlab = "", main = expression(sigma), xlim = range(temp.all), col = NA)
# plotTpolys()
# plotCI(temp, convNum(par2[1:7, 'sigma1'], type = "m"), li = convNum(par2[1:7, 'sigma1'], type = "l"), ui = convNum(par2[1:7, 'sigma1'], type = "u"), col = 2, gap = 0, cex=0.7, add=TRUE, pch = 21, pt.bg = "white")
# lines(temp.all, rep(convNum(par1[13, 3]), length(temp.all)))
# 
##################################################################
# Estimate lower bounds from Marshallagia data
##################################################################
library(dclone)
# Set upper bounds and other parameters based on Bryanne's experiments

mod_u0Lower <- function(){
	
	#Priors on two parameters
	El ~ dnorm(3.25, pow(0.1, -2))
	Tl ~ dunif(-40, 0)
	
	# Model prediction
	for(i in 1:length(rate)){
		rate_pred[i] <- a * exp(- E /(8.62*10^-5) * (1/(temp[i] + 273.15) - 1/(15 + 273.15))) * (1 + exp(El / (8.62*10^-5)*(1 / (temp[i] + 273.15) - 1/(Tl + 273.15))))
	}
	
	# Likelihood
	for(i in 1:length(rate)){
		rate[i] ~ dlnorm(log(rate_pred[i]), pow(0.2, -2))
	}
}

# u0
dat_u0 <- list(
	a = 0.068, 
	E = 0.884,
	temp = c(-9, -20),#, -35),
	rate = c(0.0855312, 0.98579604) #3)
)

fit_u0 <- jags.fit(data = dat_u0, params = c("El", "Tl"), model = mod_u0Lower)
summary(fit_u0)
plot(fit_u0, las = 1)

# # By nonlinear least squares
# fit_u0_nls <- nls(formula = y ~ 0.068 * exp(- 0.884 /(8.62*10^-5) * (1/(temp + 273.15) - 1/(15 + 273.15))) * (1 + exp(El / (8.62*10^-5)*(1 / (temp + 273.15) - 1/(Tl + 273.15)))), data = data.frame(y = c(0.0855312, 0.98579604, 3), temp = c(-9, -20, -35)), start = list(El = 3.25, Tl = 0), lower = c(0, -40), upper = c(10, 0), algorithm = "port")
# summary(fit_u0_nls)
# 

# Plots to Ale's data
temp.all <- seq(-50, 40, 0.2)
plot(temp.all, predict.MTE(params = c(a = 0.068, E = 0.884, Eh =NA, Th = NA, El = summary(fit_u0)[[1]]['El', 1], Tl = summary(fit_u0)[[1]]['Tl', 1],  z = 1), temp.all), "l", ylim = c(0,1))
# lines(temp.all, predict.MTE(params = c(a = 0.068, E = 0.884, Eh =NA, Th = NA, El = summary(fit_u0_nls)$parameters['El', 1], Tl = summary(fit_u0_nls)$parameters['Tl', 1],  z = 1), temp.all), lty =2, col = 2)								 

 # Ale's marshallagia study
points(-9, 0.0855, pch = 3, col = 5, xpd = NA, lwd = 2) #L1
points(-20, 0.9857, pch = 3, col = 5, xpd = NA, lwd = 2) #L1
points(-35, 3, pch = 3, col = 5, xpd = NA, lwd = 2) #L1

#------------------------------------------------------------------------------
# u3 fit
mod_u3Lower <- function(){
	
	#Priors on two parameters
	El ~ dnorm(3.25, pow(1, -2))
	Tl ~ dunif(-40, 0)
	
	# Model prediction
	for(i in 1:length(rate)){
		rate_pred[i] <- 0.0211 * exp(- 0.208 /(8.62*10^-5) * (1/(temp[i] + 273.15) - 1/(15 + 273.15))) * (1 + exp(El / (8.62*10^-5)*(1/(temp[i] + 273.15)-1/(Tl+273.15))) + exp(3.5543 / (8.62*10^-5)*(-1/(temp[i] + 273.15)+1/(27.6+273.15))))
	}
	
	# Likelihood
	for(i in 1:length(rate)){
		rate[i] ~ dlnorm(log(rate_pred[i]), pow(0.2, -2))
	}
}

# u0
dat_u3 <- list(
	temp = c(-9, -20),#, -35),
	rate = c(0.002344, 0.01762) #3)
)

fit_u3 <- jags.fit(data = dat_u3, params = c("El", "Tl"), model = mod_u3Lower)
summary(fit_u3)
plot(fit_u3)

# fit_u3_nls <- nls(formula = y ~ 0.0211 * exp(- 0.208 /(8.62*10^-5) * (1/(temp + 273.15) - 1/(15 + 273.15))) * (1 + exp(El / (8.62*10^-5)*(1/(temp+273.15)-1/(Tl+273.15))) + exp(3.5543 / (8.62*10^-5)*(-1/(temp+273.15)+1/(27.6+273.15)))), data = data.frame(y = c(0.002344, 0.01762), temp = c(-9, -20)), start = list(El = 3.25, Tl = 0))
# summary(fit_u3_nls)


plot(temp.all, predict.MTE(params = c(a = 0.0211, E = 0.208, Eh = 3.5543, Th = 27.6, El = summary(fit_u3)[[1]]['El', 1], Tl = summary(fit_u3)[[1]]['Tl', 1],  z = 1), temp = temp.all), "l", ylim = c(0,0.1))
			
# Ale's marshallagia study
points(-9, 0.002344, pch = 3, col = 5, xpd = NA, lwd = 2) #L3
points(-20, 0.01762, pch = 3, col = 5, xpd = NA, lwd = 2) #L3
points(-35, 3, pch = 3, col = 5, xpd = NA, lwd = 2) #L3

