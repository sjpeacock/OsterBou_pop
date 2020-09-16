
library(animation)
library(gtools)
library(fBasics) # for Heaviside function
library(doParallel)

# Acronyms used:
#	VMR: variance-to-mean ratio. Assuming a negative binomial distribution of 
#			parasites among hosts, the VMR is related to the overdispersion parameter 
#			of the neg binom distribution

###############################################################################
# Functions
###############################################################################

# Here, rho is not within-host development (assumed 0) but development of 
# non-infectious larvae

#------------------------------------------------------------------------------
# Parasite and host dependent fecundity
#------------------------------------------------------------------------------
pCalf <- function(params){
	pCalf - 1/(1 + exp(7.025 - 0.000328*P))

#------------------------------------------------------------------------------
# Time derivative of the PDE system
#------------------------------------------------------------------------------
# Input variables:
# 	Y: matrix with columns equal to the spatial grid and rows equal to the seven
#			different variables: 
#				(1) stationary hosts (H_stat),
#				(2) stationary parasites (P_stat), 
#				(3) stationary VMR (VMR_stat),
#				(4) migratory hosts (H_mig), 
#				(5) migratory parasites (P_mig), 
#				(6) migratory VMR (VMR_mig), 
#				(7) free-living larvae (larv).
# params: vector of parameters, including:
			beta
			
partial_t.caribou<-function(Y, params, migrate="BOTH"){ 
	
	if(dim(V)[1]==8){ #Both stationary and transient sub populations included
		
		d.N<- params['beta']*(V[1,] + V[4,]) - (params['mu'] + params['omega'] + params['alpha']*V[2,])*V[1,]  + (params['gamma'] + params['theta']*V[5,])*V[4,]
		d.m<- params['lambda']*V[8,] - V[2,]*(params['sigma'] + params['alpha']*V[3,] + params['beta']*(V[1,] + V[4,])/V[1,]) + (V[4,]/V[1,])*(params['gamma']*(V[5,] - V[2,]) + params['theta']*V[5,]*(V[6,] + V[5,] - V[2,]))
		d.A <- params['beta']*V[2,]*(V[1,]+V[4,])/V[1,] + (1 - V[3,])*(params['lambda']*V[8,]/V[2,] - params['rho'] + params['sigma'] + V[3,]*params['alpha']) + (V[4,]*V[5,])/(V[1,]*V[2,])*(params['theta']*(V[6,]*(3*V[5,] + 2*V[6,] - 1 - V[3,] - 2*V[2,]) + (V[5,] - V[2,])^2 - V[3,]*V[5,]) + params['gamma']*(V[5,] + V[6,] - V[3,] - 2*V[2,] + V[2,]^2/V[5,]))	
		d.N_hat<- -(params['mu'] + params['gamma'] + (params['alpha'] + params['theta'])*V[5,])*V[4,] + params['omega']*V[1,]		
		d.m_hat<-  params['lambda']*V[8,] - V[5,]*(params['sigma'] + (params['alpha'] + params['theta'])*V[6,]) + params['omega']*V[1,]/V[4,]*(V[2,] - V[5,]) 		
		d.A_hat <- (1 - V[6,])*(params['lambda']*V[8,]/V[5,] + params['sigma'] + V[6,]*(params['alpha'] + params['theta'])) + params['omega']*(V[1,]*V[2,])/(V[4,]*V[5,])*(V[2,] + V[3,] - V[6,] - 2*V[5,] + V[5,]^2/V[2,])		
		d.L0 <- params['kappa']*(V[1,]*V[2,] + V[4,]*V[5,]) - (params['mu0'] + params['rho'])*V[7,] - params['lambda']*V[7,]*(V[1,] + V[4,])		
		d.L1 <- params['rho']*V[7,] - params['mu1']*V[8,] - params['lambda']*V[8,]*(V[1,] + V[4,])		
		return(rbind(d.N, d.m, d.A, d.N_hat, d.m_hat, d.A_hat, d.L0, d.L1))
	
	}else if(dim(V)[1]==5 & migrate=="YES"){ #Only transient hosts; no switching
	
		d.N_hat<- -(params['mu'] + params['gamma'] + (params['alpha'] + params['theta'])*V[2,])*V[1,]		
		d.m_hat<-  params['lambda']*V[5,] - V[2,]*(params['sigma'] + params['alpha']*V[3,]) 	
		d.A_hat <- (1 - V[3,])*(params['lambda']*V[5,]/V[2,] + params['sigma'] + V[3,]*params['alpha'])				
		d.L0 <- params['kappa']*(V[1,]*V[2,]) - (params['mu0'] + params['rho'])*V[4,] - params['lambda']*V[4,]*V[1,]				
		d.L1 <- params['rho']*V[4,] - params['mu1']*V[5,] - params['lambda']*V[5,]*V[1,]		
		return(rbind(d.N_hat, d.m_hat, d.A_hat, d.L0, d.L1))

	}else if(dim(V)[1]==5 & migrate=="NO"){ #Only stationary hosts; no switching
	
		d.N<- (params['beta'] - params['mu'] - params['alpha']*V[2,])*V[1,]  			
		d.m<- params['lambda']*V[5,] - V[2,]*(params['sigma'] + params['alpha']*V[3,] + params['beta']) 			
		d.A <- params['beta']*V[2,] - (V[3,] - 1)*(params['lambda']*V[5,]/V[2,] + params['sigma'] + V[3,]*params['alpha']) 			
		d.L0 <- params['kappa']*(V[1,]*V[2,]) - (params['mu0'] + params['rho'])*V[4,] - params['lambda']*V[4,]*V[1,]					
		d.L1 <- params['rho']*V[4,] - params['mu1']*V[5,] - params['lambda']*V[5,]*V[1,]		
		return(rbind(d.N, d.m, d.A, d.L0, d.L1))

	}else{cat("Something wrong...didn't match a category!")}
	
} #end function

sim.caribou<-function(params, Nx, Nt, dt, dx, init, migrate="BOTH", include.full.end=FALSE){
	x.samp<-round(seq(1,Nx,8))
	t.samp<-round(seq(1,Nt,8))
		
	# 1) Set up matrices to store solutions:
	V<-matrix(NA); length(V)<-dim(init)[1]*(Nx+1)*(Nt+1); dim(V)<-c(dim(init)[1], (Nx+1), (Nt+1))
	V[,,1]<-init
				
	# Advection speed for each variable: number of grid spaces moved
	if(migrate=="BOTH"){ u<-c(rep(0,3), rep(params['c.'], 3), 0, 0)*dt/dx
	}else if(migrate=="YES"){u<-c(rep(params['c.'], 3), 0, 0)*dt/dx
	}else if(migrate=="NO"){u<-rep(0, 5)*dt/dx}
			
	# 2) Loop through timesteps	
	for(n in 1:Nt){ # For each timestep, n
			
		# Calculate boundary conditions: torus for circular migration
		Vn<-V[,,n]
		
		# Spatial advection (upstream differencing)
		Vnp1<-Vn
		for(j in 1:(dim(Vn)[1])){
			if(u[j]>0) Vnp1[j,] <- Vn[j,c((Nx+2-u[j]):(Nx+1), 1:(Nx+1-u[j]))]
		}
		
		# Temporal dynamics (4th order Runge Kutta)
		k1<-partial_t.caribou(Vnp1, params, migrate)
		k2<-partial_t.caribou(Vnp1 + dt/2*k1, params, migrate)
		k3<-partial_t.caribou(Vnp1 + dt/2*k2, params, migrate)
		k4<-partial_t.caribou(Vnp1 + dt*k3, params, migrate)
		
		V[,,n+1]<- Vnp1 + dt/6*(k1 + 2*k2 + 2*k3 + k4)
		
		} #end timestep n
		
		V.out<-V[,x.samp,t.samp]
		
		if(include.full.end == FALSE) return(list(V.out, t=t.samp, x=x.samp)) else return(list(V.out, t=t.samp, x=x.samp, end=V[,,Nt+1]))
	}


###############################################################################################
# Set up
###############################################################################################

# SIX seasons:
# S: Dec - Mar - winter
# M: Apr-May migration to breeding grounds
# S: June - calving
# M: July - post-calving migration
# S: Aug - Sep- summer grounds
# M: Oct - Nov - migration to overwinter grounds
# **May change among seasons??

mu0<-function(temp){
	0.0548*exp(-0.971/(8.62*10^-5)*(1/(temp+273.15)-1/(15+273.15)))*(1+exp(30.244/(8.62*10^-5)*(-1/(temp+273.15)+1/(39.2+273.15))))
}

mu1<-function(temp){
	0.0211*exp(-(-0.208)/(8.62*10^-5)*(1/(temp +273.15)-1/(15+273.15)))*(1+exp(3.5543/(8.62*10^-5)*(-1/(temp +273.15)+1/(27.6+273.15))))
}

rho<-function(temp){
	0.0287*exp(-0.537/(8.62*10^-5)*(1/(temp+273.15)-1/(15+273.15)))*(1+exp(3.297/(8.62*10^(-5))*(-1/(temp+273.15)+1/(31.6+273.15))))^(-1)
}



# Units: km, days
temp.seasons<-c(-20, 5, 20, 25, 15, -5)+5

params.bou<-rbind(
	beta=c(0, 0, 0.35, 0, 0, 0)/365, # birth only in calving season
	mu=c(0.1, 0, 0, 0, 0, 0)/365, #mortality only in winter
	alpha=rep(0, 6), #assume no parasite-induced mortality for now
	sigma=rep(5, 6)/365, # within-host parasite mort is 5
	kappa=c(0, 10, 80, 160, 80, 0), #parasite production per parasite per day** (Stien & Irvine 2002 Int J Parasit Fig. 3 BAM!) ** For final model will have to consider that egg production per host does not increase linearly with number of parasites!! See paper.
	mu0=mu0(temp=temp.seasons), #Mortality of uninfectives**
	mu1=mu1(temp=temp.seasons), #Moratlity of infectives**
	rho=rho(temp=temp.seasons), #Development rate of uninfectives**
	lambda=rep(0.0000025, 6),#c(0.001, 0.0001, 0.001, 0.0001, 0.001, 0.0001)/365,
	gamma=rep(0, 6), # switch to stationary
	omega=rep(0,6),# switch to transient
	theta=rep(0,6),
	c.=c(0, 16, 0, 24, 0, 16)
	)
colnames(params.bou)<-c("winter", "spring.mig", "calving", "post-calving.mig", "summer", "fall.mig")

#------------------------------------------------------------
# Plot of metabolically determined parameters:
#------------------------------------------------------------
theme.cols<-c(dg="#1D9A78", lg="#8BC145", db="#36AFCE", lb="#1D6FA9", r="#B74919", o="#F19D19")
T.all<-seq(-20,35, 0.1)

pdf("MTE_plus5.pdf", width=8, height=3)
# quartz(width=8, height=3)
par(mfrow=c(1,3), mar=c(4,4,5,1))
plot(T.all, mu0(T.all), "l",  bty="l", lwd=1.5, las=1, ylab=expression(paste("Mortality of L0 (", d^-1, ")", sep="")), ylim=c(0, 0.2), xlab=""); 
segments(x0=temp.seasons, y0=rep(0, 6), x1=temp.seasons, y1=mu0(temp.seasons), col=theme.cols[c('db','lg','lg','o', "o", "db")], lty=c(2,1))
arrows(x0=temp.seasons, y0=mu0(temp.seasons), x1=par('usr')[1], y1=mu0(temp.seasons), length=0.08, col=theme.cols[c('db','lg','lg','o', "o", "db")], lty=c(2,1))
mtext(side=3, adj=0, line=0.5 ,"a)")

plot(T.all, rho(T.all), "l",  bty="l", lwd=1.5, las=1, ylab=expression(paste("Development rate of L0 (", d^-1, ")", sep="")), xlab="") 
segments(x0=temp.seasons, y0=rep(0, 6), x1=temp.seasons, y1=rho(temp.seasons), col=theme.cols[c('db','lg','lg','o', "o", "db")], lty=c(2,1))
arrows(x0=temp.seasons, y0=rho(temp.seasons), x1=par('usr')[1], y1=rho(temp.seasons), length=0.08, col=theme.cols[c('db','lg','lg','o', "o", "db")], lty=c(2,1))
mtext(side=3, adj=0, line=0.5 ,"b)")

legend(-40, 0.085, xpd=NA, lty=rep(c(1,2), 3), col=theme.cols[c('db','db','lg','lg','o', "o")], ncol=3, legend=colnames(params.bou)[c(6,1:5)], bty="n")

plot(T.all, mu1(T.all), "l",  bty="l", lwd=1.5, las=1, ylab=expression(paste("Mortality of L1 (", d^-1, ")", sep="")), ylim=c(0, 0.10), xlab="") 
segments(x0=temp.seasons, y0=rep(0, 6), x1=temp.seasons, y1=mu1(temp.seasons), col=theme.cols[c('db','lg','lg','o', "o", "db")], lty=c(2,1))
arrows(x0=temp.seasons, y0=mu1(temp.seasons), x1=par('usr')[1], y1=mu1(temp.seasons), length=0.08, col=theme.cols[c('db','lg','lg','o', "o", "db")], lty=c(2,1))
mtext(side=3, adj=0, line=0.5 ,"c)")

mtext(side=1, outer=TRUE, expression(paste("Temperature (", degree, "C)", sep="")), line=-1)
dev.off()
#------------------------------------------------------------
#------------------------------------------------------------

# Temporal grid, different for the different seasons?
tmax<-c(121, 61, 30, 31, 61, 61)

# How far will they migrate?
sum(params.bou['c.',]*tmax) #2696 km

# Spatial grid in kms
xmin <- 0
xmax<-sum(params['c.',]*tmax)
n.x<-sum(params['c.',]*tmax)
dx <- (xmax-xmin)/n.x
x <- seq(xmin, xmax, dx)

# Temporal grid - want it so that the number of grid spaces moved at speed c in one timestep is a whole number
dt<-1/8 # timestep every 4 hours (1/8 days*24 hours/day = 4 hours)
params.bou['c.',]*dt/dx # Move 2 or 3 grid spaces in a timestep. OK
n.t<-round(tmax/dt)
sum(n.t)*dt #one-year, 365 days

# Initial condition
# k = m/(A-1), A = (m+k)/k
sd.start<-30
mean.start<-130
m.start<-10000 #Albon 2002: Fig. 2 Ostertagia abundance, within range of Stien and Irvine 2002
k.start<-0.043 + 0.69*m.start*10^(-3) #from Irvine et al. 2000 in Stien and Irvine...2002 Int J Parasit


init.bou<-matrix(c(
	c(10^4/sqrt(2*pi*(sd.start^2))*exp(-(x-mean.start)^2/(2*sd.start^2))+1), # N_hat
	rep(m.start, n.x+1), # m_hat
	rep((m.start+k.start)/k.start, n.x+1), # A_hat
	rep(1, n.x+1),
	rep(1, n.x+1)
	), nrow=5, ncol=n.x+1, byrow=TRUE)



###############################################################################################
# Simulation
###############################################################################################


y<-20

n.seasons<-length(tmax)*y
out<-list(); length(out)<-length(n.seasons)

for(i in 1:y){
	
	# 1) WINTER
	if(i==1) out[[1]]<-sim.caribou(params= params.bou[,1], Nx=n.x, Nt=n.t[1], dt=dt, dx=dx, init=init.bou, migrate="NO", include.full.end=TRUE)
	if(i>1) out[[(i-1)*6 + 1]]<-sim.caribou(params= params.bou[,1], Nx=n.x, Nt=n.t[1], dt=dt, dx=dx, init=out[[(i-1)*6]]$end, migrate="NO",include.full.end=TRUE)
	
	# 2) spring migration
	out[[(i-1)*6 + 2]]<-sim.caribou(params= params.bou[,2], Nx=n.x, Nt=n.t[2], dt=dt, dx=dx, init=out[[(i-1)*6 + 1]]$end, migrate="YES", include.full.end=TRUE)
	
	# 3) CALVING
	out[[(i-1)*6 + 3]]<-sim.caribou(params= params.bou[,3], Nx=n.x, Nt=n.t[3], dt=dt, dx=dx, init=out[[(i-1)*6 + 2]]$end, migrate="NO",include.full.end=TRUE)
	
	# 4) post-calving migration
	out[[(i-1)*6 + 4]]<-sim.caribou(params= params.bou[,4], Nx=n.x, Nt=n.t[4], dt=dt, dx=dx, init=out[[(i-1)*6 + 3]]$end, migrate="YES", include.full.end=TRUE)
	
	# 5) SUMMER
	out[[(i-1)*6 + 5]]<-sim.caribou(params= params.bou[,5], Nx=n.x, Nt=n.t[5], dt=dt, dx=dx, init=out[[(i-1)*6 + 4]]$end, migrate="NO", include.full.end=TRUE)
	
	# 6) fall migration
	out[[i*6]]<-sim.caribou(params= params.bou[,6], Nx=n.x, Nt=n.t[6], dt=dt, dx=dx, init=out[[(i-1)*6 + 5]]$end, migrate="YES", include.full.end=TRUE)
	
	
}	

tmaxminus<-c(0, tmax)

date.x<-matrix(NA, nrow=1, ncol=3)
for(i in 1:y){
	for(s in 1:6){
		w<-seq(10, tmax[s], 10)
		date.x<-rbind(date.x, cbind(rep(i, length(w)), rep(s, length(w)), sum(tmaxminus[1:s]) +w)) 
		if(i==1&s==1){
			V.all<-out[[(i-1)*6 + s]][[1]][,,w]
		}else{
			nis<-dim(V.all)[3]
			V.all1<-matrix(rep(NA, 5*dim(V.all)[2]*(nis+length(w)))); dim(V.all1)<-c(5, dim(V.all)[2], nis+length(w))
			V.all1[,,1:nis]<-V.all
			V.all1[,,(nis+1):(nis+length(w))]<-out[[(i-1)*6 + s]][[1]][,,w]
			V.all<-V.all1
			}
		}}
	
dim(V.all)		


 # #Plot 5 years
# for(i in 1:y){
	# for(s in c(1,3,5)){
		# w1<-seq(1, tmax[s], 10)
		# w2<-seq(1, tmax[s+1], 10)
		# par(mfrow=c(4,1), mar=c(4,4,3,1))
		# for(j in c(1,2,4,5)){
			# plot(out[[(i-1)*6 + s]]$x, out[[(i-1)*6 + s]][[1]][j,,1], "n", ylim=range(V.all[j,,]), ylab=c("Hosts", "Parasites", NA, "Uninfectious Larvae", "Infectious Larvae")[j], bty="l", xlab="Distance (km)")
			# if(j==1) mtext(side=3, paste("Year ", i, ", ", colnames(params)[s], sep=""))
			# for(k in 1:length(w1)){
				# lines(out[[(i-1)*6 + s]]$x, out[[(i-1)*6 + s]][[1]][j,,w1[k]], col=rainbow(n=length(w1))[k], lty=2)
			# }
			# for(k in 1:length(w2)){
				# lines(out[[(i-1)*6 + s + 1]]$x, out[[(i-1)*6 + s+1]][[1]][j,,w2[k]], col=rainbow(n=length(w2))[k])			
			# }	
# }	}}	
		

#####################################################################################################
#####################################################################################################
# Track: 
#	- infectious parasite larvae at each location (winter, breeding grounds, summer)
#	- mean parasite burden across all hosts

date.x2<-as.Date(paste(2000+date.x[,1], date.x[,3], sep="-"), format="%Y-%j")
# date.x3<-strftime(date.x2, format="%m-%d")
loc.seasons<-unique(cumsum(params.bou['c.',]*tmax))[1:3]
L0<-matrix(rep(NA, 3*dim(V.all)[3]), nrow=dim(V.all)[3], ncol=3)
L1<-L0
m.avg<-numeric(dim(V.all)[3])

for(i in 1:dim(V.all)[3]){ #every 10 days
	for(j in 1:3){
		L0[i,j]<-V.all[4,findInterval(loc.seasons[j], out[[1]]$x)+1,i]
		L1[i,j]<-V.all[5,findInterval(loc.seasons[j], out[[1]]$x)+1,i]
		}
	m.avg[i]<-sum(V.all[1,,i]/sum(V.all[1,,i])*V.all[2,,i])
	}
	
	
# Plot
plot.seasons<-function(s){
	for(i in 1:y){
	if(s==1)x1<-as.Date(paste(1999+c(1:y)[i], c(12,6,8)[s], 01, sep="-")) else x1<-as.Date(paste(2000+c(1:y)[i], c(12,6,8)[s], 01, sep="-"))
	x2<-as.Date(paste(2000+c(1:y)[i], c(3,6,9)[s], c(31, 30, 30)[s], sep="-"))
	u<-par('usr')
	polygon(x=c(x1,x2,x2,x1), y=c(u[c(3,3,4,4)]), border=NA, col=paste(theme.cols[c('db', 'lg', 'o')[s]], "20", sep=""))
	}
}

theme.cols<-c(dg="#1D9A78", lg="#8BC145", db="#36AFCE", lb="#1D6FA9", r="#B74919", o="#F19D19")
pdf(file="Annual_plus5.pdf", width=10, height=8.5)
par(mfrow=c(4,1), mar=c(3,5,2,1))
for(s in 1:3){
	plot(date.x2[2:length(date.x2)], L1[,s], "l", col=theme.cols[c('db', 'lg', 'o')[s]], las=1, bty="l", ylab="Free-living larvae", xaxt="n", xlab="", main=paste("Parasite larvae on", c("wintering", "calving", "summer")[s], "grounds"))
	lines(date.x2[2:length(date.x2)], L0[,s], lty=2, col=theme.cols[c('db', 'lg', 'o')[s]])
	axis(side=1, at=as.Date(paste(2000+c(1:y), 01, 01, sep="-")), labels=c(1:y))
	plot.seasons(s)
	
	if(s==1) legend('topleft', lty=c(2,1,NA), c("Non-infectious", "Infectious", "Caribou present"), pch=c(NA, NA, 15), col=c(rep(theme.cols['db'], 2), paste(theme.cols['db'], 20, sep="")), bg="white", pt.cex=3)
	}

plot(date.x2[2:length(date.x2)], m.avg, "l", col=theme.cols['r'], las=1, bty="l", ylab="Mean parasites per host", xaxt="n", xlab="Time (years)", main=c("Parasites per host"))
# abline(h=m.start, lty=2, col=theme.cols['r'])
axis(side=1, at=as.Date(paste(2000+c(1:y), 01, 01, sep="-")), labels=c(1:y))
plot.seasons(1); plot.seasons(2); plot.seasons(3)
mtext(side=1, outer=TRUE)

dev.off()
