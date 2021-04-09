# Map of the Bathurst caribou migration and ranges

library(PBSmapping)
gshhg <- "~/Google Drive/Mapping/gshhg-bin-2.3.7/"
xlim <- c(-128, -124) + 360
ylim <- c(48.6, 51.4)
land <- importGSHHS(paste0(gshhg,"gshhs_i.b"), xlim = xlim, ylim = ylim, maxLevel = 2, useWest = TRUE)
rivers <- importGSHHS(paste0(gshhg,"wdb_rivers_i.b"), xlim = xlim, ylim = ylim, useWest = TRUE)
borders <- importGSHHS(paste0(gshhg,"wdb_borders_i.b"), xlim = xlim, ylim = ylim, useWest = TRUE, maxLevel = 1)

#------------------------------------------------------------------------------
# Sampling site locations
siteDat <- data.frame(uniqueSiteID = unique(paste(dat$region, dat$location, dat$date, sep = "--")))
siteDat$EID <- 1:length(siteDat$uniqueSiteID)
siteDat$nFish <- rep(NA, dim(siteDat)[1])
siteDat$X <- rep(NA, dim(siteDat)[1])
siteDat$Y <- rep(NA, dim(siteDat)[1])
siteDat$region <- rep(NA, dim(siteDat)[1])
siteDat$location <- rep(NA, dim(siteDat)[1])

for(i in 1:dim(siteDat)[1]){
	dum <- which(paste(dat$region, dat$location, dat$date, sep = "--") == siteDat$uniqueSiteID[i])
	siteDat$nFish[i] <- length(dum)
	siteDat$X[i] <- - dat$lon[dum[1]]
	siteDat$Y[i] <- dat$lat[dum[1]]
	siteDat$region[i] <- dat$region[dum[1]]
	siteDat$location[i] <- dat$location[dum[1]]
}

siteDat1 <- as.EventData(siteDat[which(is.na(siteDat$X) == FALSE), ], projection = "LL")


#------------------------------------------------------------------------------
# Farm locations
datFarm <- read.csv('data/lice-audit-verif-pou-2011-ongoing-rpt-pac-dfo-mpo-aquaculture-eng.csv')
datFarm$Audit.Date <- as.Date(datFarm$Audit.Date, format = "%d-%b-%y")
datFarm <- subset(datFarm, datFarm$Audit.Date >= "2019-01-01")	

farmLoc <- data.frame(farmName = unique(datFarm$Site.Common.Name))
farmLoc$EID <- 1:length(farmLoc$farmName)
farmLoc$X <- datFarm$Longitude[match(farmLoc$farmName, datFarm$Site.Common.Name)]
farmLoc$Y <- datFarm$Latitude[match(farmLoc$farmName, datFarm$Site.Common.Name)]
farmLoc$FishHealthZone <- datFarm$Fish.Health.Zone[match(farmLoc$farmName, datFarm$Site.Common.Name)]
farmLoc$Company <- datFarm$Licence.Holder[match(farmLoc$farmName, datFarm$Site.Common.Name)]

farmLoc <- as.EventData(farmLoc, projection = "LL")
# Plot map
quartz(width = 5, height = 6)
plotMap(land, xlim = xlim - 360, ylim = ylim,	col = grey(0.8), bg = "aliceblue", las = 1, lwd = 0.5, border = grey(0.6))
addPoints(farmLoc, pch = 14 + as.numeric(as.factor(c(farmLoc$Company))), cex = 0.6)
addPoints(siteDat1, pch = pchRegion[siteDat1$region], col = colRegion[siteDat1$region])#, col = c(3, 2, 4)[as.numeric(as.factor(siteDat1$region))])
legend("bottomleft", pch = c(15:17, pchRegion), pt.lwd = 1.5, col = c(rep(1, 3), colRegion), legend = c("Farm (Cermaq)", "Farm (Greig)", "Farm (MOWI)", regions), bg = "white", pt.cex = c(rep(0.6, 3), rep(1, length(regions))), cex = 0.8)


#######################################