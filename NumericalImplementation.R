library(fields)

# constants:
sigma <- 1/8.29 # See LatinAmerica_parameters.R
gamma <- 1/5 # See LatinAmerica_parameters.R
# Quarantine efficacy:
q <- 99.2342 # quarantine efficacy; Assessing the Effectiveness of Mass Testing and Quarantine in the Spread of COVID-19 in Beijing and Xinjiang, 2020; Complexity 2021

# Better measure:
# For CR, Pan., Ur.: Hasell, J., Mathieu, E., Beltekian, D. et al. A cross-country database of COVID-19 testing. Sci Data 7, 345 (2020). https://doi.org/10.1038/s41597-020-00688-8
# Tests per million per day
tests_by_country <- c('CostaRica'=1459/1e6, 'Panama'=2287/1e6, 'Uruguay'= 2477/1e6, 'Tejas_revised'=1909.887/1e6)
# For TX: https://coronavirus.jhu.edu/testing/states-comparison (Accessed Aug. 24): 101415 tests per 100k so daily test total per million is total tests per million/#days in pandemic: (101415*10)/as.numeric(difftime(as.Date('2021-08-24'),as.Date('2020-03-11')))

recovery_lag <-  8.29 + 5 # See LatinAmerica_parameters.R

Rm <- function(b, B1, sigmaV, q, su, sv)
	{
	# Eurozone way of doing it
	(B1*sigma*(q + sigmaV)*su + b*B1*(q + sigma)*sigmaV*sv)/sqrt(B1*(gamma + q)*(q + sigma)*(q + sigmaV)*(n)*(sigma*(q + sigmaV)*su + b*(q + sigma)*sigmaV*sv))	
	}

get_Rms <- function(sigmaV, q, dats_by_day, countryID)
	{
	bRange <- seq(0,1,length=30)
	R_ms <- matrix(nrow=nrow(dats_by_day), ncol=30)
	for (i in 1:length(dats_by_day[,'Date']))
		{
		for (j in 1:30)
			{
			R_ms[i,j] <- Rm(bRange[j], B1[1], sigmaV, q, dats_by_day[i,'Su'], dats_by_day[i,'Sv'])
			}
		}
	return(R_ms)
	}

# Sensitivity analysis
# scenario 1 (drag; much longer incubation period), scenario 2 (identical incubation period), scenario 3 (much shorter incubation period)
sigma_vs <- c(0.1, 1,10)
# increase in 
B_factor 

# Have a directory with subdirectories where input data are organized by region
locations <- c('Panama','CostaRica','Tejas_revised','Uruguay')
B1s <- read.table('LatinAmerica_Data/LKbeta/betaVals.dat',header=T)
heatMaps <- list()
map_to_heatMap <- c()

counter <- 1

# To plot the time series:

par(mfrow=c(2,2))
plot(B1s[which(B1s$LK_ID=='CostaRica'),3],t='l',xlab='Day', ylab=expression('Estimated'~beta), main='Costa Rica')
plot(B1s[which(B1s$LK_ID=='Panama'),3],t='l',xlab='Day', ylab=expression('Estimated'~beta), main='Panama')
plot(B1s[which(B1s$LK_ID=='Tejas_revised'),3],t='l',xlab='Day', ylab=expression('Estimated'~beta), main='Texas')
plot(B1s[which(B1s$LK_ID=='Uruguay'),3],t='l',xlab='Day', ylab=expression('Estimated'~beta), main='Uruguay')

for (sig_v in sigma_vs)
	{
	for (countryID in locations)
		{
		d <- read.csv(paste('/home/kokamoto/Documents/Oshigoto/CoronaChan/VaccineResistance/LatinAmerica_Data/',countryID,'.csv', sep=''))
		
		recovereds <- c(rep(0, as.integer(recovery_lag)), d$total_cases)[1:nrow(d)]
		
		# remove all records from before people were fully vaccinated:
		non_vax_days <- 1:min(which(d$people_fully_vaccinated > 0))
		d <- d[-non_vax_days,]
		recovereds <- recovereds[-non_vax_days]
		
		excluded_dates <- unique(c(which(is.na(d$people_fully_vaccinated)), which(is.na(d$total_deaths)), which(is.na(d$total_cases)), which(is.na(recovereds))))
			
		Su <- d$population - d$total_cases - ifelse(is.na(d$people_fully_vaccinated), 0, d$people_fully_vaccinated) - ifelse(is.na(d$total_deaths), d$total_deaths,0)
		Sv <- d$people_fully_vaccinated + recovereds
		Date <- d$date
		B1 <- B1s[which(B1s[,'LK_ID']==countryID),'b_est']

		if (length(excluded_dates) > 0)
			{
			Su <- Su[-excluded_dates]
			Sv <- Sv[-excluded_dates]
			Date <- d$date[-excluded_dates]
			B1 <- B1s[which(B1s[,'LK_ID']==countryID),'b_est'][-excluded_dates]
			}

		dats_by_day <- data.frame(Date, Su, Sv)
		n <- unique(d$population[1])
		my_q <- tests_by_country[countryID] / (n/1e6) # Convert tests per million into tests per capita

		heatMaps[[counter]] <- get_Rms(sig_v, my_q, dats_by_day, countryID)
		map_to_heatMap <- c(map_to_heatMap, paste(countryID, '_', sig_v,sep=''))
		counter <- counter + 1
		}
	}

draw_heatmap <- function(R_ms, countryID)
	{
	bRange <- seq(0,1,length=30)
	# pretty-fication
	pdf(file=paste('R_m_for_', countryID, '.pdf',sep=''), width=6,height=2.4)
	par(mfrow=c(1,3))
	for (i in 1:length(sigma_vs))
		{
		sigmaV <- sigma_vs[i]
		target <- which(map_to_heatMap==paste(countryID,'_', sigmaV,sep=''))
		# image.plot(x=1:nrow(R_ms[[target]]), y=bRange, ifelse(R_ms[[target]] < 1, NA, R_ms[[target]]), xlab='Days since vaccination begins', ylab='Mutant relative infectivity', main=bquote(paste('Evolutionary Instability; ', sigma[v]==.(sigmaV))), zlim=c(0,max(unlist(R_ms))))
		image.plot(x=1:nrow(R_ms[[target]]), y=bRange, ifelse(R_ms[[target]] < 1, NA, R_ms[[target]]), xlab='Days since vaccination begins', ylab='Mutant relative infectivity', main='Evolutionary Instability', zlim=c(0,max(unlist(R_ms))))
		mtext(bquote(sigma[v]==.(sigmaV)))
		}
	dev.off()
	}

for (countryID in locations)
	{
	draw_heatmap(heatMaps, countryID)
	}

