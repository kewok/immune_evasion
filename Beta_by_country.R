locations <- c('Panama','CostaRica','Tejas_revised','Uruguay')

btest_index <- 1:length(seq(0.01,3.5,by=0.01))

Days <- 170 # days for which we have records on all regions
E_vals <- list(locations)


for (countryID in locations)
	{
	E_vals[[countryID]] <- matrix(nrow=Days, ncol=length(btest_index))
	
	for (Date in 1:Days) 
		{
		for (j in 1:length(btest_index))
			{
			inname = sprintf('LatinAmerica_Data/LKbeta/E_%s_%s_%s.dat',countryID, j, Date)
			E_vals[[countryID]][Date,j] <- read.table(inname, header=T)$E
			}
		}
	}
	
beta.df <- data.frame(nrow=Days*length(locations),ncol=3)

current_row <- 1
bvals <- seq(0.01,3.5,by=0.01)
for (countryID in locations)
	{
	for (Date in 1:Days)
		{
		beta.df[current_row, 1] <- Date
		beta.df[current_row, 2] <- bvals[which(E_vals[[countryID]][Date,]==min(E_vals[[countryID]][Date,]))[1]]
		beta.df[current_row, 3] <- countryID
		current_row <- current_row + 1
		}
	}
	
colnames(beta.df) <- c('Date','b_est','LK_ID')
