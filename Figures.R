library(RColorBrewer)
# constants:
sigma <- 1/5.2 # See LatinAmerica_parameters.R
gamma <- 1/3 # See LatinAmerica_parameters.R
# Quarantine efficacy:
q <- 99.2342 # quarantine efficacy; Assessing the Effectiveness of Mass Testing and Quarantine in the Spread of COVID-19 in Beijing and Xinjiang, 2020; Complexity 2021

# Better measure:
# For CR, Pan., Ur.: Hasell, J., Mathieu, E., Beltekian, D. et al. A cross-country database of COVID-19 testing. Sci Data 7, 345 (2020). https://doi.org/10.1038/s41597-020-00688-8
# Tests per million per day
tests_by_country <- c('CostaRica'=1459/1e6, 'Panama'=2287/1e6, 'Uruguay'= 2477/1e6, 'Tejas_revised'=1909.887/1e6)
# For TX: https://coronavirus.jhu.edu/testing/states-comparison (Accessed Aug. 24): 101415 tests per 100k so daily test total per million is total tests per million/#days in pandemic: (101415*10)/as.numeric(difftime(as.Date('2021-08-24'),as.Date('2020-03-11')))

Days <- 170 # days post first vaccination for which we have records on all regions except TX

recovery_lag <-  5 # See LatinAmerica_parameters.R


Rm <- function(b, B1, sigmaV, q, su, sv)
	{
	# Eurozone way of doing it
	(B1*sigma*(q + sigmaV)*su + b*B1*(q + sigma)*sigmaV*sv)/sqrt(B1*(gamma + q)*(q + sigma)*(q + sigmaV)*(n)*(sigma*(q + sigmaV)*su + b*(q + sigma)*sigmaV*sv))	
	}

bRange <- c(0.05,0.5,1,5) # relative infectivity of the mutant strain

get_Rms <- function(sigmaV, q, dats_by_day, countryID)
	{
	R_ms <- matrix(nrow=nrow(dats_by_day), ncol=length(bRange))
	for (i in 1:length(dats_by_day[,'Date']))
		{
		for (j in 1:length(bRange))
			{
			R_ms[i,j] <- Rm(bRange[j], B1[1], sigmaV, q, dats_by_day[i,'Su'], dats_by_day[i,'Sv'])
			}
		}
	return(R_ms)
	}
	
# Since the sensitivity analysis suggested the delay in incubation had limited effect:

sigma_v <- 1

# Have a directory with subdirectories where input data are organized by region
locations <- c('Panama','CostaRica','Tejas_revised','Uruguay')
B1s <- read.table('LatinAmerica_Data/LKbeta/betaVals_revised.dat',header=T)
heatMaps <- list()
incidence <- list()
susceptibles <- list()
vaxes <- list()
map_to_heatMap <- c()
beta_ests <- list()

counter <- 1

for (countryID in locations)
	{
	d <- read.csv(paste('/home/kokamoto/Documents/Oshigoto/CoronaChan/VaccineResistance/LatinAmerica_Data/',countryID,'.csv', sep=''))
			
	# remove all records from before people were fully vaccinated:
	if (countryID != 'Tejas_revised')
		{
		non_vax_days <- 1:min(which(!is.na(d[,'new_vaccinations_smoothed'])))
		d <- d[-non_vax_days,]
		}

# Determine the time-lagged vaccination rate; assume it takes 30 days to become fully vaccinated
	if(countryID=='Tejas_revised')
		{
		for (i in 1:nrow(d))
			{
			if (is.na(d$people_fully_vaccinated[i]))
				{
				d$people_fully_vaccinated[i] <-  d$people_fully_vaccinated[i-1]
				}
			# As there is a date in Texas with NA deaths
			if (is.na(d$total_deaths[i]))
				{
				d$total_deaths[i] <-  d$total_deaths[i-1]
				}
			}
		Sv <- d$people_fully_vaccinated + d$total_cases
		}
	else
		{
		# impute fully vaccinated data if missing; note the non-Texas data begin when vaccinations began.
		for (i in 31:nrow(d))
			{
			if (is.na(d$people_fully_vaccinated[i]))
				{
				d$people_fully_vaccinated[i] <- d$people_fully_vaccinated[i-1]
				}
			}
		Sv <-  c(rep(0,30), d[31:nrow(d),'new_vaccinations_smoothed']) +  d$total_cases
		}
		
	Su <- d$population - Sv - d$total_deaths
	B1 <- B1s[which(B1s[,'LK_ID']==countryID),'b_est']
	beta_ests[[counter]] <- B1
	dats_by_day <- data.frame(Date=1:Days, Su=Su[1:Days], Sv=Sv[1:Days])
	n <- unique(d$population[1])
	my_q <- tests_by_country[countryID] / (n/1e6) # Convert tests per million into tests per capita
	if (countryID == 'Tejas_revised')
		{
		incidence[[counter]] <- (d$New.Probable.Cases[1:Days] / n) * 1e6
		}
	else
		{
		incidence[[counter]] <- d$new_cases_smoothed_per_million[1:Days]
		} 
	susceptibles[[counter]] <- Su[1:Days] / n * 1e6 # Convert susceptibles to per-million
	vaxes[[counter]] <- Sv[1:Days]/n # Convert vaccinated to per-million

	heatMaps[[counter]] <- get_Rms(sigma_v, my_q, dats_by_day, countryID)
	map_to_heatMap <- c(map_to_heatMap, paste(countryID, '_', sigma_v,sep=''))
	counter <- counter + 1
	}

labs <- c('A','B','C','D')

par(las=1,mfrow=c(2,2),mar=c(6,7,5,5))	

cols <- brewer.pal(8,"Set1")

for (counter in 1:4)
	{
	plot(vaxes[[counter]],t='l',xlab='Days since vaccinations began',ylab='Fraction of population recovered\n or fully vaccinated',main=labs[counter],ylim=range(vaxes), col=cols[1],cex.lab=1.8,lwd=2,cex.main=3,cex.axis=1.5)

	if (counter != 3)
		{
	legend(25,0.5, c(expression(beta[mu]~'='~'0.05'), expression(beta[mu]~'='~'0.5'), expression(beta[mu]~'='~'1'), expression(beta[mu]~'='~'5')), lwd=c(3,3,3,3), lty=rep(2,4),cex=1.5,col=cols[2:5],box.col='WHITE')
		}
	else
		{
		legend(100,0.23, c(expression(beta[mu]~'='~'0.05'), expression(beta[mu]~'='~'0.5'), expression(beta[mu]~'='~'1'), expression(beta[mu]~'='~'5')), lwd=c(3,3,3,3), lty=rep(2,4),cex=1.5,col=cols[2:5],box.col='WHITE')
		}
	
	cols <- brewer.pal(8,"Set1")
	par(new=T,las=1)
	plot(log(heatMaps[[counter]][,3]), axes=F, xlab='', ylab='',t='l',lwd=3,ylim=range(lapply(heatMaps, log)),lty=2,col=cols[2])
	lines(log(heatMaps[[counter]][,1]),lty=2,col=cols[3],lwd=3)
	lines(log(heatMaps[[counter]][,2]),lty=2,col=cols[4],lwd=3)
	lines(log(heatMaps[[counter]][,4]),lty=2,col=cols[5],lwd=3)
	abline(h=0)
	axis(side=4, at=pretty(c(0.25,2)),cex.axis=1.5)
	mtext(expression(atop(Log,R[mu])), side=4, line=1, padj=-2.5,cex=1.8)
	}
	
# Supplementary Time Series for Beta, Su:

par(las=1,mfrow=c(2,2),mar=c(6,7,5,5))	
for (counter in 1:4)
	{
	plot(susceptibles[[counter]],t='l',xlab='Days since vaccinations began',ylab='',main=labs[counter],ylim=range(susceptibles), col=cols[1],cex.lab=1.8,lwd=2,cex.main=3,cex.axis=1.5)
	par(las=0)
	mtext('Susceptible hosts per million', side=2, line=3.4, padj=-1.5,cex=1.4)

	if (counter==1)
		{
		legend(75,7.5e5, c(expression(S[u](t)), expression(beta[A](t))), lwd=c(2,2),col=cols[1:2],cex=1.5)
		}
	par(new=T,las=1)
	plot(beta_ests[[counter]], axes=F, xlab='', ylab='',t='l',lwd=2,ylim=range(beta_ests),col=cols[2])
	axis(side=4, at=pretty(c(0.25,4)),cex.axis=1.4)
	mtext(expression(beta[A]), side=4, line=3, padj=-1.5,cex=1.4)
	}
