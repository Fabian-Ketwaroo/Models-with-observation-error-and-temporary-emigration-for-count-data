setwd("~/Objective priors for N-mixture models/Parrots")

library(tidyverse)
library(lubridate)


orange_parrots = read.csv("orange_wing_data.csv")

head(orange_parrots)

orange_parrots$Year = as.numeric(format(as.Date(orange_parrots$Date, format="%d/%m/%Y"),"%Y"))
orange_parrots$months = as.numeric(format(as.Date(orange_parrots$Date, format="%d/%m/%Y"),"%m"))
head(orange_parrots)

nb_months = c(9,10,11,12)
nb = length(nb_months)
J= 12

C_y1 = matrix(NA,J , nb) 
Bred_year1 = orange_parrots

a = Bred_year1[ which(Bred_year1$Year==2004 & Bred_year1$months == 9 ), ]
a$Total.Ind
dim(a)

for (i in 1:length(nb_months)) {
  a = Bred_year1[ which(Bred_year1$Year==2004 & Bred_year1$months == nb_months[i] ), ]
  l = length(a$Total.counts)
  C_y1[1:l,i] = a$Total.counts
}

C_y1 # rows are sampling occasions and columns are months

nb_months = 1:9
nb = length(nb_months)

C_y2 = matrix(NA,J , nb) 
Bred_year1 = orange_parrots

a = Bred_year1[ which(Bred_year1$Year==2005 & Bred_year1$months == 1 ), ]
a$Total.Ind
dim(a)

for (i in 1:length(nb_months)) {
  a = Bred_year1[ which(Bred_year1$Year==2005 & Bred_year1$months == nb_months[i] ), ]
  l = length(a$Total.counts)
  C_y2[1:l,i] = a$Total.counts
}

C_y2 # rows are sampling occasions and columns are months

Oparrots = cbind(C_y1, C_y2)
Oparrots # here the rows are the sampling occasions and the columns are the months 

save(Oparrots,file =  "Oparrots.Rdata")

# Covariates
head(orange_parrots)
#Need to convert time to 0,1
S = length(orange_parrots$Total.counts)
orange_parrots$ctime = numeric(length(S))
head(orange_parrots)

for (s in 1:S) {
  
  if(orange_parrots$Time[s]=="PM") orange_parrots$ctime[s] = 1
  else  orange_parrots$ctime[s] = 0
  
}

orange_parrots
sum(orange_parrots$ctime)


orange_parrots$crain = numeric(length(S))
head(orange_parrots)

for (s in 1:S) {
  
  if(orange_parrots$Rain[s]==1) orange_parrots$crain[s] = 1
  else  orange_parrots$crain[s] = 2
  
}

head(orange_parrots)

orange_parrots$pcloud = orange_parrots$cloudy = numeric(length(S))
head(orange_parrots)

for (s in 1:S) {
  
  if(orange_parrots$Cloudiness[s]==2) orange_parrots$pcloud[s] = 1
  else  orange_parrots$pcloud[s] = 0
  
  if(orange_parrots$Cloudiness[s]==3) orange_parrots$cloudy[s] = 1
  else  orange_parrots$cloudy[s] = 0

    
}

head(orange_parrots)

orange_parrots$low_wind = orange_parrots$strong_wind = numeric(length(S))
head(orange_parrots)

for (s in 1:S) {
  
  if(orange_parrots$Wind[s]==2) orange_parrots$low_wind[s] = 1
  else  orange_parrots$low_wind[s] = 0
  
  if(orange_parrots$Wind[s]==3) orange_parrots$strong_wind[s] = 1
  else  orange_parrots$strong_wind[s] = 0
  
  
}

head(orange_parrots)

save(orange_parrots, file =  "orange_parrots.Rdata" )

