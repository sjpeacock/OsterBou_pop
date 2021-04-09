###############################################################################
# Annual dynamics at peak and trough
###############################################################################

# 1) What years to look at?
# year2plot <- array(NA, dim = c(3,3,3,2))
# par(mfrow = c(1,1), mar = c(4,4,2,1), oma = rep(0, 4))
# for(i in 1:3){ #for each transmission
# 	for(r in 1:2){# for each climate scenario
# 		for(p in 1:3){# for each inhibition
# 			plot(80:y, annualSumm[[i, p, r]][80:y, 1], "o", col = cols[p])
# 			year2plot[i, p, r, ] <- round(locator(2)$x)
# 			}}}
# saveRDS(year2plot, file = "year2plot_PeakTrough.rds")
year2plot <- readRDS("year2plot_PeakTrough.rds")

oneYear <- array(NA, dim = c(3, 3, 3, 2, 5, 365))
# dimensions = beta (3) * inhibition (3) * climate (3) * peak/trough (2) * variable * day

for(i in 1:3){
	for(p in 1:3){ # for each level of ppnInhibit (1, 0.5, 0)
		for(r in 1:3){ # for each climate scenario
			for(m in 1:2){ # for peak and trough
				V.i <- V[[i, p, r]][, , which(timeDat$year[which(timeDat$year > 80)] == year2plot[i, p, r, m])]
				for(j in 1:365){
					oneYear[i, p, r, m, 1, j] <- sum(V.i[c('P_mov', 'P_stat'), , j]) / sum(V.i[c('adult_mov', 'adult_stat'), , j])
					oneYear[i, p, r, m, 2, j] <- sum(V.i[c('L4_mov', 'L4_stat'), , j]) / sum(V.i[c('adult_mov', 'adult_stat'), , j])
					oneYear[i, p, r, m, 3, j] <- sum(V.i[c('L4A_mov', 'L4A_stat'), , j]) / sum(V.i[c('adult_mov', 'adult_stat'), , j])
					oneYear[i, p, r, m, 4, j] <- sum(V.i['L0', , j])
					oneYear[i, p, r, m, 5, j] <- sum(V.i['L3', , j])
				} # end day j
				
			}}}}

#------------------------------------------------------------------------------
# Plot for given transmission and climate scenario

i <- 2
r <- 1

xDate <- as.Date(paste(30, c(1:365), sep = "-"), format = "%y-%j") 
yPerc <- 0.1
quartz(width = 6.3, height = 8, pointsize = 10)
layout(mat = cbind(
	matrix(c(rep(1:3, each = 3), 11, rep(4:5, each = 3)), ncol = 1),
	5 + matrix(c(rep(1:3, each = 3), 7, rep(4:5, each = 3)), ncol = 1)))
par(oma = c(4, 5, 4, 1))

par(mar = c(0,4,0,1))

for(m in 1:2){ # For peak and trough
	for(s in 1:5){ # For each stage
		if(s == 3|s == 5){
			plot(xDate, oneYear[i, 1, r, m, s, ], "n", ylim = extendrange(oneYear[i, , r, , s, ], f = c(0, yPerc)), yaxt = "n", ylab = "", xaxs = "i")
		}else{
			plot(xDate, oneYear[i, 1, r, m, s, ], "n", ylim = extendrange(oneYear[i, , r, , s, ], f = c(0, yPerc)), yaxt = "n", ylab = "", xaxt = "n", xaxs = "i")
		}
		
		for(p in 1:3){
			lines(xDate, oneYear[i, p, r, m, s, ], col = cols[p], lwd = 1.5)
		}
		
		if(m == 1) mtext(side = 2, line = 4.5, c("Adult\nparasites", "Arrested\nlarvae", "Developing\nlarvae", "Pre-infective\nlarvae", "Infective\nlarvae")[s])
		if(m == 1 & s == 1){
			legend("topright", col = cols[c(1:3)], title = "Inhibition", lwd = 1.5, legend = c("100%", "50%", "0%"))
		}
		
		if(s == 1) mtext(side = 3, line = 1.3, c("Peak", "Trough")[m])
	}
}






