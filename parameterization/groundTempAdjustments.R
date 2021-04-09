c("Summer",  "Winter",  "Spring",  "Fall", "Calving")

# mean temperature
ck <- c(-10.974553,  -8.029248,  -9.808518,  -9.750996, -11.466570) [1]
# half annual temperature range (amplitude)
dk <- c(22.546752,  23.653142,  22.808474,  22.990450, 21.577856)[1]
# DOY corresponding to max temperature
t0 <- c(199.922313, 198.146132, 199.304408, 199.167954, 201.066183)[1]

DOY <- c(1:365)

tempDOY <- ck + dk * cos((DOY - t0)* 2 * pi / 365)

# Ground temp
tempDOY2 <- tempDOY

# Minimum ground temperatures of -15*C
tempDOY2[which(tempDOY < -15)] <- -15

# Delay DD accumulation by 2 weeks
d0 <- which(tempDOY > 0)[1]
tempDOY2[d0:(d0 + 14)] <- 0

# 3 *C warmer in the summer
summerDOY <- c((d0 + 15):max(which(tempDOY > 0)))
tempDOY3 <- ck + (dk + 5) * cos((summerDOY - t0)* 2 * pi / (365 - 14))

adjustDown <- seq(-tempDOY3[1], 0, length.out = 30)
tempDOY3[1:30] <- tempDOY3[1:30] + adjustDown
# Smooth out end of summer
adjustDown2 <- seq(- tempDOY3[length(tempDOY3)], 0, length.out = 30)
tempDOY3[(length(tempDOY3) - 30 + 1): length(tempDOY3)] <- tail(tempDOY3, 30) + rev(adjustDown2
)

tempDOY2[summerDOY] <- tempDOY3


# Plot adjustments
quartz(width = 4.8, height = 3, pointsize =10)
par(mar = c(4,4,4,1))
plot(DOY, tempDOY, "l", las = 1, ylab = expression(paste("Temperature (", degree, "C)", sep = "")), lwd = 2, ylim = c(-35, 20), col = grey(0.6), xaxs = "i", lty = 3)
abline(v = t0, lty = 2)
abline(h = 0)
abline(h = -15, lty = 2, col = 4)
abline(h = c(max(tempDOY), max(tempDOY2)), lty = 2, col = 2)
# text(50, max(tempDOY) + 1.5, col = 2, expression(paste(3, degree, "C")))
# arrows(30, max(tempDOY), 30, max(tempDOY2), col = 2, length = 0.08)
lines(DOY, tempDOY2, lwd = 2)
axis(side = 3, at = as.numeric(strftime(as.Date(c(paste(2100, seq(1,12,2), 01, sep = "-"))), format = "%j")), labels = c("Jan", "Mar", "May", "Jul", "Sep", "Nov"))
points(c(145, 200, 343), c(-3.4, 8.4, -18.2), cex = 2, pch = 21, bg = "white")
text(c(145, 200, 343), c(-3.4, 8.4, -18.2), c(1:3), cex =0.8)
