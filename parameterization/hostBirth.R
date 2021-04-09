# WHat is the relationship between calf production and parasite burdens?
# From Albon et al. 2002: based on treated and untreated pregnancy and calving rates
#m Note Hughes et al. 2009 found no effect of lumen parasites on probability
# that DU caribou were pregnancy, though there was an effect on body condition, and
# it is well established that body condition affects pregnancy rates.

AlbonData <- read.csv("parameterization/Albon2002_1a.csv",header = TRUE)

plot(AlbonData$intensity, AlbonData$calfDiff, xlim=c(0, 16000), ylim = c(0, 0.2))
lines(seq(0,16000,10), 1/(1 + exp(7.025 - 0.000328*seq(0,16000,10))), lwd=1.5)

# num. calves per female host per year, accounting for parasites
# pCalf : probability of having a calf in the absence of parasites
pCalf <- 0.85
P <- seq(0, 20000, 100)

pCalfP <- pCalf - 1/(1 + exp(7.025 - 0.000328*P))

par(mar = c(4,4,2,4))
plot(P, pCalfP, "l", bty="l", ylim = c(0.3, 1), ylab="Probability of calf", xlab = "Parasite intensity",  yaxt="n")
abline(h = pCalf, lty = 2)
axis(side = 2, at = seq(0.4, 1, 0.2), las=1)
text(8000, 0.9, paste("Difference = ", round(pCalf - pCalfP[findInterval(7275, P)], 5)), col = 2)

par(new = TRUE)
plot(P, dnbinom(P, size = 0.9940684, mu = 914), 'l', col = 2, xaxt="n", yaxt="n", xlab = "", ylab = "", bty = "u")
axis(side = 4, col = 2, las = 1)
mtext(side = 4, col = 2, "Probability density based on Bathurst survey", line = 3)
abline(v = 7275, lty = 3, col = 2)



# for(i in 1:3){
# 	pCalf <- c(0.8, 0.9, 1)[i]
# 	lines(P, pCalf - 1/(1 + exp(7.025 - 0.000328*P)), col = i)
# 	abline(h = pCalf, lty = 2, col = i)
# }
text(5500, 0.5, "Pr(calf) in absence of parasites (f) \ntemperature dependent??")
