#
# Look at starting and stopping grids
# 

y <- 10
dt <- 0.2
stochTemp <- FALSE
startstop <- "simple"

source("simulations/calcGrid.R")

initBou <- rbind(
	calf_mov = rep(0, length(x)),
	yearling_mov = rep(0, length(x)),
	adult_mov = rep(0, length(x)),
	L4_mov = rep(0, length(x)),
	L4A_mov = rep(0, length(x)),
	P_mov = rep(0, length(x)),
	calf_stat = initDist(totPop = 10^5 * bouDist1985['calf']/sum(bouDist1985), x),
	yearling_stat = initDist(totPop = 10^5 * bouDist1985['yearling']/sum(bouDist1985), x),
	adult_stat = initDist(totPop = 10^5 * c(bouDist1985['cow'] + bouDist1985['bull'])/sum(bouDist1985), x),
	L4_stat = rep(0, length(x)),
	L4A_stat = rep(0, length(x)),
	P_stat = initP * initDist(totPop = 10^5 * c(bouDist1985['cow'] + bouDist1985['bull'])/sum(bouDist1985), x),
	L0 = rep(0, length(x)),
	L3 = rep(0, length(x))
)

V1 <- simBou(initBou = initBou, tempGrid = tempGrid, ppnInhibit = 1)

y <- 10
dt <- 0.5
stochTemp <- FALSE
startstop <- "simple"

source("simulations/calcGrid.R")

initBou <- rbind(
	calf_mov = rep(0, length(x)),
	yearling_mov = rep(0, length(x)),
	adult_mov = rep(0, length(x)),
	L4_mov = rep(0, length(x)),
	L4A_mov = rep(0, length(x)),
	P_mov = rep(0, length(x)),
	calf_stat = initDist(totPop = 10^5 * bouDist1985['calf']/sum(bouDist1985), x),
	yearling_stat = initDist(totPop = 10^5 * bouDist1985['yearling']/sum(bouDist1985), x),
	adult_stat = initDist(totPop = 10^5 * c(bouDist1985['cow'] + bouDist1985['bull'])/sum(bouDist1985), x),
	L4_stat = rep(0, length(x)),
	L4A_stat = rep(0, length(x)),
	P_stat = initP * initDist(totPop = 10^5 * c(bouDist1985['cow'] + bouDist1985['bull'])/sum(bouDist1985), x),
	L0 = rep(0, length(x)),
	L3 = rep(0, length(x))
)

V2 <- simBou(initBou = initBou, tempGrid = tempGrid, ppnInhibit = 1)

################################################################################
# Plotting
dayo <- "10-100"
	
plot(as.numeric(rownames(V1[1,,])), apply(V1[c('adult_stat', 'adult_mov'), , which(colnames(V1[1,,]) == dayo)], 2, sum), "l", col = 2, ylim = range(c(apply(V1[c('adult_stat', 'adult_mov'), , which(colnames(V1[1,,]) == dayo)], 2, sum), apply(V2[c('adult_stat', 'adult_mov'), , which(colnames(V2[1,,]) == dayo)], 2, sum))))
lines(as.numeric(rownames(V2[1,,])), V2['adult_stat', , which(colnames(V2[1,,]) == dayo)] + V2['adult_mov', , which(colnames(V2[1,,]) == dayo)], col = 4, lty = 3)


plot(as.numeric(rownames(V1[1,,])), apply(V1[c('P_mov', 'P_stat'), , which(colnames(V1[1,,]) == dayo)], 2, sum), "l", col = 2, ylim = range(c(apply(V1[c('P_mov', 'P_stat'), , which(colnames(V1[1,,]) == dayo)], 2, sum), apply(V2[c('P_mov', 'P_stat'), , which(colnames(V2[1,,]) == dayo)], 2, sum))))
lines(as.numeric(rownames(V2[1,,])), V2['P_stat', , which(colnames(V2[1,,]) == dayo)] + V2['P_mov', , which(colnames(V2[1,,]) == dayo)], col = 4, lty = 3)


saveHTML(expr = {
	for(y in 1:10){
		for(d in seq(1, 365, 10)){
			dayo <- paste(y, d, sep = "-")
			
			plot(as.numeric(rownames(V1[1,,])), V1['adult_stat', , which(colnames(V1[1,,]) == dayo)], "l", ylim = range(c(V1['adult_stat', , which(colnames(V1[1,,]) == dayo)], V2['adult_stat', , which(colnames(V2[1,,]) == dayo)])))
			
			lines(as.numeric(rownames(V2[1,,])), V2['adult_stat', , which(colnames(V2[1,,]) == dayo)], lty = 2, lwd = 1.5)
	}
})