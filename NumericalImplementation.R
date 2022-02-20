library(fields)

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
	
bRange <- seq(0,5,length=30) # relative infectivity of the mutant strain

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

# Sensitivity analysis
# scenario 1 (drag; much longer incubation period), scenario 2 (identical incubation period), scenario 3 (much shorter incubation period)
sigma_vs <- c(0.1, 1,10)

# Have a directory with subdirectories where input data are organized by region
locations <- c('Panama','CostaRica','Tejas_revised','Uruguay')
B1s <- read.table('LatinAmerica_Data/LKbeta/betaVals.dat',header=T)
heatMaps <- list()
map_to_heatMap <- c()

counter <- 1

for (sig_v in sigma_vs)
	{
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

		dats_by_day <- data.frame(Date=1:Days, Su=Su[1:Days], Sv=Sv[1:Days])
		n <- unique(d$population[1])
		my_q <- tests_by_country[countryID] / (n/1e6) # Convert tests per million into tests per capita

		heatMaps[[counter]] <- get_Rms(sig_v, my_q, dats_by_day, countryID)
		map_to_heatMap <- c(map_to_heatMap, paste(countryID, '_', sig_v,sep=''))
		counter <- counter + 1
		}
	}

draw_heatmap <- function(R_ms, countryID)
	{
	# pretty-fication
	pdf(file=paste('R_m_for_', countryID, '.pdf',sep=''), width=6,height=2.4)
	par(mfrow=c(1,3))
	for (i in 1:length(sigma_vs))
		{
		sigmaV <- sigma_vs[i]
		target <- which(map_to_heatMap==paste(countryID,'_', sigmaV,sep=''))
		#image.plot(x=1:nrow(R_ms[[target]]), y=bRange, ifelse(R_ms[[target]] < 1, NA, R_ms[[target]]), xlab='Days since vaccination begins', ylab='Mutant relative infectivity', main='Evolutionary Instability', zlim=c(0,max(unlist(R_ms))))
		image.plot(x=1:nrow(R_ms[[target]]), y=bRange, R_ms[[target]], xlab='Days since vaccination begins', ylab='Mutant relative infectivity', main='Evolutionary Instability')
		mtext(bquote(sigma[v]==.(sigmaV)))
		}
	dev.off()
	}

for (countryID in locations)
	{
	draw_heatmap(heatMaps, countryID)
	}

