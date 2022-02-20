locations <- c('Panama','CostaRica','Tejas_revised','Uruguay')
btest_index <- 1:length(seq(0.01,3.5,by=0.01))
parameter_space_global <- expand.grid(locations, btest_index, 1:236)

colnames(parameter_space_global) <- c('country', 'btest_index', 'time')

args <- (commandArgs(trailingOnly=TRUE))
n = as.numeric(args[1])

# Get the value to the nearest 10000th and create a subdirectory

for (rep in 1:300)
	{
	i <- 300*(n-1) + rep
	subdir <- paste('subdir', n, sep='')
	if (sum(subdir!=list.files()))
  		{
		system(paste('mkdir', subdir))
		}

	setwd(subdir)
	directory <- paste('Scenario',i,sep='')
	system(paste('mkdir', directory))
	setwd(directory)
	write.table(parameter_space_global[i,], file='parameter_combo.txt', col.names=T, row.names=F)
	system('cp ../../LatinAmerica_parameters.R .')
	system('cp ../../Reimplement_SEIRfunc2.R .')
	source('LatinAmerica_parameters.R')
	setwd('../..')
	}

