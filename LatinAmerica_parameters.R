# Loosely following code in: 
# https://osf.io/7dshm/

# install.packages('/home/kokamoto/Documents/Oshigoto/CoronaChan/VaccineResistance/RawData_Inferring_beta/R-Code/covid19up_0.1.0.tar.gz', repos=NULL)
library(ggplot2)
library(grid)
library(dplyr)
library(pracma)
library(covid19up)
library(rstudioapi)

locations <- c('Panama','CostaRica','Tejas_revised','Uruguay')
populations <- c()

source('Beta_by_country.R')

write.table(file='LatinAmerica_Data/LKbeta/betaVals_revised.dat',beta.df,quote=F,row.names=F)
