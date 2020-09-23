cols <- 1:4

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
