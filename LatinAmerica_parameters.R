# Loosely following code in: 
# https://osf.io/7dshm/

library(ggplot2)
library(grid)
library(dplyr)
library(pracma)
library(covid19up)
library(rstudioapi)

populations <- c()
source('Reimplement_SEIRfunc2.R')

parameter_space <- read.table('parameter_combo.txt',header=T)
countryID <- parameter_space['country']
j <- as.numeric(parameter_space['btest_index'])
t <- as.numeric(parameter_space['time'])

d <- read.csv(paste('../../LatinAmerica_Data/',countryID,'.csv', sep=''))

btest = seq(0.01,3.5,by=0.01)
a <- 1/8.29 # incubation period; see Qin, J. et al. Estimation of incubation period distribution of COVID-19 using disease onset forward time: a novel cross-sectional and forward follow-up study. Sci. Adv. 6, eabc1202 (2020). 8.29 days.
g <- 1/5 # Recovery rate (days infectious); see Epidemiology and transmission of COVID-19 in 391 cases and 1286 of their close contacts in Shenzhen, China: a retrospective cohort study. Lancet Infect Dis. 2020. says: five days.
ief <- g/a
rho <- 10
dt <- 1
Nruns <- 100

if (countryID=='Tejas_revised')
	{
	I0 <- sum(diff(d$Cumulative.Probable.Cases[1:14] + d$Cumulative.Confirmed.Cases[1:14])) # Confirm logic with Luis
	d <- d[14:nrow(d),]
	d$Date <- as.Date(d$date,'%m/%d/%Y')
	R0 <- d$total_cases[1] - d$Cumulative.Fatalities[1] + d$people_fully_vaccinated[1]
	}
if (countryID!='Tejas_revised')	
	{
	I0 <- sum(diff(d$total_cases[(min(which(!is.na(d[,'new_vaccinations_smoothed'])))-14):(min(which(!is.na(d[,'new_vaccinations_smoothed'])))) ])) # Prevalence based on all infectious from 14 days ago through present.
	R0 <- d$total_cases[min(which(!is.na(d[,'new_vaccinations_smoothed'])))-14]
	d <- d[min(which(!is.na(d[,'new_vaccinations_smoothed']))):nrow(d),]
	d$Date <- as.Date(d[min(which(!is.na(d[,'new_vaccinations_smoothed']))):nrow(d),'date'],'%m/%d/%Y')
	}
d$dt = c(diff(d$Date),0)
d$Time <- 1:length(d$Date)
tdiff <- diff(d$Date)
Time <- length(d$Time)

date1 <- min(d$Date)
date2 <- max(d$Date)

if (t <= length(d$Date))
	{
	tv = bv = ev = dv = c()
	# initial conditions; use only post-vaccination situation
	N <- unique(d$population)
	populations <- c(populations, N)


	E0 = round(I0*(ief+rnorm(1,mean=0,sd=0.25)))
	V0 = ifelse(countryID=='Tejas_revised', d$people_fully_vaccinated[1], 0)
	S0 = N - E0 - I0 - R0 - V0

	# Determine the time-lagged vaccination rate; assume it takes 30 days to become fully vaccinated
	if(countryID=='Tejas_revised')
		{
		vaxRate <-  c(diff(d$people_fully_vaccinated[1:31])/(N-d$total_cases[1:30] - d$Cumulative.Fatalities[1:30] - d$people_fully_vaccinated[1:30]), d$daily_vaccinations[31:nrow(d)]/(N-d$total_cases[31:nrow(d)] - d$Cumulative.Fatalities[31:nrow(d)] - d$people_fully_vaccinated[31:nrow(d)]))
		}
	else
		{
		# impute fully vaccinated data if missing
		for (i in 31:nrow(d))
			{
			if (is.na(d$people_fully_vaccinated[i]))
				{
				d$people_fully_vaccinated[i] <- d$people_fully_vaccinated[i-1]
				}
			}
		vaxRate <-  c(rep(0,30), d[31:nrow(d),'new_vaccinations_smoothed']/(N-d$total_cases[31:nrow(d)]-d$total_deaths[31:nrow(d)]-d$people_fully_vaccinated[31:nrow(d)]))
		}

	# Implement the EnKF:
	numb = length(btest)
	E = matrix(data=rep(0,(Time-1)*numb),ncol=numb,nrow=Time-1)
	best = rep(0,Time-1)
	b = btest[j]
	Z = matrix(data=rep(c(S0,E0,I0,R0),Nruns),ncol=4,nrow=Nruns,byrow=T)
	Z[,2] = round(Z[,2]*(1+rnorm(Nruns,mean=0,sd=0.5)))
	dt = d$dt[t]

	# Simulate ensemble
	sim = seir_custom(N=N,I0=Z[,3],E0=Z[,2],R0=Z[,4], V0=V0, b=b,a=a,g=g,Time=dt,runs=Nruns, vaxRate=vaxRate)
	Z[,] <- as.numeric(as.matrix(sim[,paste('Var',1:4,sep='')]))

	# Ensemble Kalman filter
	H = c(0,0,1,1)
	P = cov(Z)
	K = P%*%H/as.numeric(H%*%P%*%H+rho)

	Zprime = Z
	for ( r in 1:Nruns ) {
		Zprime[r,] = round(Z[r,]- 0.5*K*as.numeric(t(H)%*%Z[r,] + t(H)%*%colMeans(Z) - 2*d$total_cases[t+1])) # Note they use cumulative cases in xinfer.R
		Zprime[r,1] = N - sum(Zprime[r,2:4]) - sim[r,'Var5']    # enforce N = S+E+I+R+V
		}

	E[t,j] = E[t,j] + 1/(2*(H%*%P%*%H+rho))*norm(H%*%colMeans(Z)-d$total_cases[t+1])^2 + 0.5*log(H%*%P%*%H+rho)


	Z = Zprime
	# if ( d$Date[t+1]==date1 )  ZS = Z
	# if ( d$Date[t+1]==date2 )  ZL = Z

	E.df = data.frame(time=t, beta=b, E=E[t,j])
	outname = sprintf('../../LatinAmerica_Data/LKbeta/E_%s_%s_%s.dat',countryID,j,t)
	write.table(file=outname,E.df,quote=F,row.names=F)
	}


num <- 0

#for (countryID in locations) {
#	num <- num + 1
#	inname = sprintf('LatinAmerica_Data/LKbeta/E_%s.dat,countryID')
#	est = read.table(inname,header=T)

#	Time = max(est$time)
#	best = vector(length=Time)
#	for ( i in 1:Time ) {
#		idx = which(est$time==i)
#		E = est$E[idx]
#		idx1 = which.min(E)
#		best[i] = est$beta[idx[idx1]]
#		}
#	date = est$date[1:Time]
#	b.df = data.frame(Time=1:Time, Date=date, b_est=best,
#	LK_ID=rep(countryID, Time))
#	if ( num==1 )  beta.df = b.df
#	else  beta.df = rbind(beta.df,b.df)
#	}

#write.table(file='LatinAmerica_Data/LKbeta/betaVals.dat',beta.df,quote=F,row.names=F)
