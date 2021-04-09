###############################################################################
# From Siten et al. (2002 Int J Parasit)
# Larval shedding depends on paraiste infection intensity and time of year
###############################################################################

DOY <- c(1:365)/365

beta <- 0.49
a1 <- 0.010
a2 <- 0.345

# Peak feces production time 
mu <- 0.52*365

sig <- 0.087
k1 <- 0.55
k2 <- 0.02

FPR <- a1 + a2/(sig*sqrt(2*pi)) * exp(-(DOY - mu)^2/(2*sig^2))

plot(as.Date(paste("2020", DOY*365, sep="-"), format = "%Y-%j"), FPR, "l")
