
# Cleaning data 
library(readxl)
parrots = read_excel("2017-2019PC.xlsx")
head(parrots)
dim(parrots)

library(tidyverse)
Bred_year1  = filter(parrots, Breeding.year == 1)
head(Bred_year1)
dim(Bred_year1)

Bred_year2  = filter(parrots, Breeding.year == 2)
head(Bred_year2)
dim(Bred_year2)

Bred_year3  = filter(parrots, Breeding.year == 3)
head(Bred_year3)
dim(Bred_year3)


# Breeding year 1
nb_months = 12; J = 10
C_y1 = matrix(NA,J , nb_months ) 
head(Bred_year1)

a = Bred_year1[ which(Bred_year1$Breeding.month == 7), ]
a$Total.Ind
dim(a)

for (i in 1: nb_months) {
 a = Bred_year1[ which(Bred_year1$Breeding.month == i), ]$Total.Ind 
 l = length(a)
 C_y1[1:l,i] = a
}

C_y1 # rows are sampling occasions and columns are months


# Breeding year 2
C_y2 = matrix(NA, J, nb_months ) 

a = Bred_year2[ which(Bred_year1$Breeding.month == 7), ]
a$Total.Ind
dim(a)

for (i in 1: nb_months) {
  a = Bred_year2[ which(Bred_year2$Breeding.month == i), ]$Total.Ind 
  l = length(a)
  C_y2[1:l,i] = a
}

C_y1; C_y2


# Breeding year 3
C_y3 = matrix(NA, J, nb_months ) 

for (i in 1:nb_months) {
  a = Bred_year3[ which(Bred_year3$Breeding.month == i), ]$Total.Ind
  print(length(a))

}
# Looks like there is 10 sampling occasions in some months now

for (i in 1: nb_months) {
  a = Bred_year3[ which(Bred_year3$Breeding.month == i), ]$Total.Ind 
  l = length(a)
  C_y3[1:l,i] = a
}

C_y3

# Different number of sampling (J)
# To account for the different number of years maybe I can create some indicator in the observation process. 
# Somewhere in the loop to say stop at the jth occasion for month t
# I can also add a counter which says that if there is NA in the data then probability of detection is zero.

cparrots = cbind(C_y1, C_y2, C_y3)
cparrots

#save(cparrots, file = "cparrots.Rdata")

# Cleaning data to obtain am/pm information
# Let 1 be am and 2 be pm
T =  nb_months*3
index = numeric(T)

for (t in 1:T) {
  index[t] = sum(1- is.na(cparrots[,t])) # create indicator where to end loop
  
}

index

time <- matrix(NA, J, T)

for (t in 1:T) {
  time[1:index[t],t] <- rep(c(1,2), index[t]/2)
  
}

time

#save(time, file = "parrots_time.Rdata")

#---------------------------------------------------------------------------------------------------------------------------------------------------------------

# Line plots for each year to get an idea of the counts

C_y1;C_y2;C_y3 # counts for each year 

df <- as.data.frame(matrix(0, nb_months*3, 3 ))
dim(df)
colnames(df) <- c("Counts", "Month", "Year")
head(df)

c1  <- colSums(C_y1, na.rm = T)
c2  <- colSums(C_y2, na.rm = T)
c3 <- colSums(C_y3, na.rm = T)

df$Counts = c(c1,c2,c3)
head(df)

m1 <- m2<- m3 <- 1:12 #c("Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct")
df$Month = c(m1,m2,m3)

y1 <- rep("Year 1", 12)
y2 <- rep("Year 2", 12)
y3 <- rep("Year 3", 12)

df$Year <- c(y1,y2,y3)
head(df)

df$Rate =  df$Counts/index
#yr1 <- df[1:12,]

# x axis treated as continuous variable
#df2$dose <- as.numeric(as.vector(df2$dose))
df$Month <- as.factor(df$Month)
ggplot(data=df, aes(x=Month, y=Rate, group=Year, color = Year)) +
  geom_line() + geom_point() +
  theme_minimal()

# Can change Months to Nov-Oct 
m1 <- m2<- m3 <- c("Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct")
df$Month = c(m1,m2,m3)

y1 <- rep("Year 1", 12)
y2 <- rep("Year 2", 12)
y3 <- rep("Year 3", 12)

df$Year <- c(y1,y2,y3)
head(df)
#yr1 <- df[1:12,]

df$Month <- factor(df$Month, levels = c("Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct") )
ggplot(data=df, aes(x=Month, y=Rate, group=Year, color = Year)) +
  geom_line() + geom_point() +
  theme_minimal() 

#ggplot(data=df, aes(x=Month, y=Rate, group=Year, color = Year)) + geom_point() 

C_y1; time[, 1:nb_months]

dim(Bred_year1)
head(Bred_year1)
yr1 <- as.data.frame(matrix(0, 90, 3))
colnames(yr1) <- c("Counts", "Month", "AM.PM")
head(yr1)
yr1$Counts <- Bred_year1$Total.Ind
yr1$Month <- Bred_year1$Month
yr1$AM.PM <- Bred_year1$AM.PM
head(yr1)

yr1$Month <- factor(yr1$Month, levels = c("Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct") )
ggplot(data=yr1, aes(x=Month, y=Counts, group= AM.PM, color = AM.PM, shape= AM.PM)) +  geom_point(size = 3)+ theme_minimal() + ylab("AM/PM Counts") +ggtitle("Year 1")

l = dim(Bred_year2)[1]
yr2 <- as.data.frame(matrix(0, l, 3))
colnames(yr2) <- c("Counts", "Month", "AM.PM")
head(yr2)
yr2$Counts <- Bred_year2$Total.Ind
yr2$Month <- Bred_year2$Month
yr2$AM.PM <- Bred_year2$AM.PM
head(yr2)

yr2$Month <- factor(yr2$Month, levels = c("Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct") )
ggplot(data=yr2, aes(x=Month, y=Counts, group= AM.PM, color = AM.PM, shape= AM.PM)) +  geom_point(size = 3)+ theme_minimal() + ylab("AM/PM Counts") +ggtitle("Year 2")


l = dim(Bred_year3)[1]
yr2 <- as.data.frame(matrix(0, l, 3))
colnames(yr2) <- c("Counts", "Month", "AM.PM")
head(yr2)
yr2$Counts <- Bred_year3$Total.Ind
yr2$Month <- Bred_year3$Month
yr2$AM.PM <- Bred_year3$AM.PM
head(yr2)

yr2$Month <- factor(yr2$Month, levels = c("Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct") )
ggplot(data=yr2, aes(x=Month, y=Counts, group= AM.PM, color = AM.PM, shape= AM.PM)) +  geom_point(size = 3)+ theme_minimal() + ylab("AM/PM Counts") +ggtitle("Year 3")


# Using each point as AM/PM hard to intrepret. I can add all occasions of AM/PM to get sum of counts

C_y1

sumy1 <- as.data.frame(matrix(0, 2*nb_months,3))
dim(sumy1)
colnames(sumy1) <- c("Counts", "Month", "AM.PM")
#M <- seq(1,12,by=2)
for (t in 1:12 ) {
  print(sum(C_y1[c(1,3,5,7,9),t], na.rm = T))
  print(sum(C_y1[c(2,4,6,8,10),t], na.rm = T))
  
}

C_y1
head(yr1)
id2 <- index[1:12]

a <- b<-  numeric(nb_months)
for (t in 1:12 ) {
 a[t] <- sum(C_y1[c(1,3,5,7,9),t], na.rm = T)/(id2[t]/2)
 b[t] <-  sum(C_y1[c(2,4,6,8,10),t], na.rm = T)/(id2[t]/2)
  
}
a;b
sumy1[seq(1,24,2), 1] <- a
sumy1[seq(1,24,2)+1, 1] <- b
sumy1

m1 <- c("Nov", "Nov", "Dec","Dec", "Jan", "Jan", "Feb", "Feb", "Mar", "Mar", "Apr", "Apr", "May", "May", "Jun", "Jun", "Jul", "Jul", "Aug", "Aug", "Sep", "Sep", "Oct", "Oct")
sumy1$Month <- m1
sumy1

sumy1$AM.PM <- Bred_year1$AM.PM[1:24]
sumy1

sumy1$Month <- factor(sumy1$Month, levels = c("Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct") )
ggplot(data=sumy1, aes(x=Month, y= Counts, group= AM.PM, color = AM.PM)) + geom_line() + geom_point() + theme_minimal() +ylab("Rate")+ ggtitle("Year 1")


C_y2

sumy1 <- as.data.frame(matrix(0, 2*nb_months,3))
dim(sumy1)
colnames(sumy1) <- c("Counts", "Month", "AM.PM")
#M <- seq(1,12,by=2)
for (t in 1:12 ) {
  print(sum(C_y2[c(1,3,5,7,9),t], na.rm = T))
  print(sum(C_y2[c(2,4,6,8,10),t], na.rm = T))
  
}

C_y2
id2 <- index[12:24]
a <- b<-  numeric(nb_months)
for (t in 1:12 ) {
  a[t] <- sum(C_y2[c(1,3,5,7,9),t], na.rm = T)/(id2[t]/2)
  b[t] <-  sum(C_y2[c(2,4,6,8,10),t], na.rm = T)/(id2[t]/2)
  
}
a;b
sumy1[seq(1,24,2), 1] <- a
sumy1[seq(1,24,2)+1, 1] <- b
sumy1

m1 <- c("Nov", "Nov", "Dec","Dec", "Jan", "Jan", "Feb", "Feb", "Mar", "Mar", "Apr", "Apr", "May", "May", "Jun", "Jun", "Jul", "Jul", "Aug", "Aug", "Sep", "Sep", "Oct", "Oct")
sumy1$Month <- m1
sumy1

sumy1$AM.PM <- Bred_year1$AM.PM[1:24]
sumy1


sumy1$Month <- factor(sumy1$Month, levels = c("Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct") )
ggplot(data=sumy1, aes(x=Month, y= Counts, group= AM.PM, color = AM.PM)) + geom_line() + geom_point() + theme_minimal() +ylab("Rate")+ ggtitle("Year 2")

C_y3

sumy1 <- as.data.frame(matrix(0, 2*nb_months,3))
dim(sumy1)
colnames(sumy1) <- c("Counts", "Month", "AM.PM")

id2 <- index[24:36]
a <- b<-  numeric(nb_months)
for (t in 1:12 ) {
  a[t] <- sum(C_y3[c(1,3,5,7,9),t], na.rm = T)/(id2[t]/2)
  b[t] <-  sum(C_y3[c(2,4,6,8,10),t], na.rm = T)/(id2[t]/2)
  
}
a;b
sumy1[seq(1,24,2), 1] <- a
sumy1[seq(1,24,2)+1, 1] <- b
sumy1

m1 <- c("Nov", "Nov", "Dec","Dec", "Jan", "Jan", "Feb", "Feb", "Mar", "Mar", "Apr", "Apr", "May", "May", "Jun", "Jun", "Jul", "Jul", "Aug", "Aug", "Sep", "Sep", "Oct", "Oct")
sumy1$Month <- m1
sumy1

sumy1$AM.PM <- Bred_year1$AM.PM[1:24]
sumy1


sumy1$Month <- factor(sumy1$Month, levels = c("Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct") )
ggplot(data=sumy1, aes(x=Month, y= Counts, group= AM.PM, color = AM.PM)) + geom_line() + geom_point() + theme_minimal() +ylab("Rate")+ ggtitle("Year 3")


#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Cleaning covariates for p 
# Becca and Simon think that precipitation, wind speed, and visibility to be the most relevant.

library(readxl)
parrots = read_excel("2017-2019PC.xlsx")
head(parrots)

parrots = read_excel("2017-2019PC.xlsx")
head(parrots)
covariates <- read_excel("Roost.Covariates.xlsx")
head(covariates[,15:21])
dim(covariates)
dim(parrots) # matches 
covariates$Breeding.year = parrots$Breeding.year
covariates$Total.precipitaton = as.numeric(covariates$Total.precipitaton)

library(tidyverse)
Bred_year1  = filter(covariates, Breeding.year == 1)
head(Bred_year1)
dim(Bred_year1)

Bred_year2  = filter(covariates, Breeding.year == 2)
head(Bred_year2)
dim(Bred_year2)

Bred_year3  = filter(covariates, Breeding.year == 3)
head(Bred_year3)
dim(Bred_year3)


# Precipitation
# Breeding year 1
nb_months = 12; J = 10
C_y1 = matrix(NA,J , nb_months ) 
head(Bred_year1)

a = Bred_year1[ which(Bred_year1$Breeding.month == 7), ]
a$Total.precipitaton
dim(a)

for (i in 1: nb_months) {
  a = Bred_year1[ which(Bred_year1$Breeding.month == i), ]$Total.precipitaton
  l = length(a)
  C_y1[1:l,i] = a
}

C_y1 # rows are sampling occasions and columns are months
C_y1[1:2,1] = C_y1[1:2,7] =0 # set missing values to zero so I dont have to specify a prior on it
index[1:12]

# Breeding year 2
C_y2 = matrix(NA, J, nb_months ) 

for (i in 1: nb_months) {
  a = Bred_year2[ which(Bred_year2$Breeding.month == i), ]$Total.precipitaton
  l = length(a)
  C_y2[1:l,i] = a
}

C_y2
C_y2[1:8,4] = 0
index[13:24]

# Breeding year 3
C_y3 = matrix(NA, J, nb_months ) 

for (i in 1: nb_months) {
  a = Bred_year3[ which(Bred_year3$Breeding.month == i), ]$Total.precipitaton
  l = length(a)
  C_y3[1:l,i] = a
}

C_y3
index[25:36]
C_y3[1:2,5]  = C_y3[3:4,11] = 0


precipitation = cbind(C_y1, C_y2, C_y3)
precipitation
dim(precipitation)
save(precipitation, file = "precipitation.Rdata")

# Wind seeed
C_y1 = matrix(NA,J , nb_months ) 
head(Bred_year1)

a = Bred_year1[ which(Bred_year1$Breeding.month == 7), ]
a$Average.Wind.Speed
dim(a)

for (i in 1: nb_months) {
  a = Bred_year1[ which(Bred_year1$Breeding.month == i), ]$Average.Wind.Speed
  l = length(a)
  C_y1[1:l,i] = a
}

C_y1 # rows are sampling occasions and columns are months


# Breeding year 2
C_y2 = matrix(NA, J, nb_months ) 

for (i in 1: nb_months) {
  a = Bred_year2[ which(Bred_year2$Breeding.month == i), ]$Average.Wind.Speed
  l = length(a)
  C_y2[1:l,i] = a
}

C_y1; C_y2


# Breeding year 3
C_y3 = matrix(NA, J, nb_months ) 

for (i in 1: nb_months) {
  a = Bred_year3[ which(Bred_year3$Breeding.month == i), ]$Average.Wind.Speed
  l = length(a)
  C_y3[1:l,i] = a
}

C_y3

wind_speed = cbind(C_y1, C_y2, C_y3)
wind_speed
save(wind_speed, file = "wind_speed.Rdata")

# Visibility
C_y1 = matrix(NA,J , nb_months ) 
head(Bred_year1)

a = Bred_year1[ which(Bred_year1$Breeding.month == 7), ]
a$Visability
dim(a)

for (i in 1: nb_months) {
  a = Bred_year1[ which(Bred_year1$Breeding.month == i), ]$Visability
  l = length(a)
  C_y1[1:l,i] = a
}

C_y1 # rows are sampling occasions and columns are months


# Breeding year 2
C_y2 = matrix(NA, J, nb_months ) 

for (i in 1: nb_months) {
  a = Bred_year2[ which(Bred_year2$Breeding.month == i), ]$Visability
  l = length(a)
  C_y2[1:l,i] = a
}

C_y1; C_y2


# Breeding year 3
C_y3 = matrix(NA, J, nb_months ) 

for (i in 1: nb_months) {
  a = Bred_year3[ which(Bred_year3$Breeding.month == i), ]$Visability
  l = length(a)
  C_y3[1:l,i] = a
}

C_y3

Visability = cbind(C_y1, C_y2, C_y3)
Visability
save(Visability, file = "Visability.Rdata")



#===============================================================================================================================================

#Cleaning covariates for p 

library(readxl)
parrots = read_excel("2017-2019PC.xlsx")
head(parrots)

parrots = read_excel("2017-2019PC.xlsx")
head(parrots)
covariates <- read_excel("Roost.Covariates.xlsx")
head(covariates[,15:21])
dim(covariates)
dim(parrots) # matches 
covariates$Breeding.year = parrots$Breeding.year

library(tidyverse)
Bred_year1  = filter(covariates, Breeding.year == 1)
head(Bred_year1)
dim(Bred_year1)

Bred_year2  = filter(covariates, Breeding.year == 2)
head(Bred_year2)
dim(Bred_year2)

Bred_year3  = filter(covariates, Breeding.year == 3)
head(Bred_year3)
dim(Bred_year3)

#covariates Temp = Temp, Humidity = Humidity, Visability = Visability, AWS=AWS, Rain = Rain, Storm = Storm, Time = Time, weather = weather
nb_months = 12; J = 10
T =  nb_months*3
index = numeric(T)

for (t in 1:T) {
  index[t] = sum(1- is.na(cparrots[,t])) # create indicator where to end loop
  
}

index

time <- matrix(NA, J, T)

for (t in 1:T) {
  time[1:index[t],t] <- rep(c(1,2), index[t]/2)
  
}

time

# Wind seeed
C_y1 = matrix(NA,J , nb_months ) 
head(Bred_year1)

a = Bred_year1[ which(Bred_year1$Breeding.month == 7), ]
a$Average.Wind.Speed
dim(a)

for (i in 1: nb_months) {
  a = Bred_year1[ which(Bred_year1$Breeding.month == i), ]$Average.Wind.Speed
  l = length(a)
  C_y1[1:l,i] = a
}

C_y1 # rows are sampling occasions and columns are months


# Breeding year 2
C_y2 = matrix(NA, J, nb_months ) 

for (i in 1: nb_months) {
  a = Bred_year2[ which(Bred_year2$Breeding.month == i), ]$Average.Wind.Speed
  l = length(a)
  C_y2[1:l,i] = a
}

C_y1; C_y2


# Breeding year 3
C_y3 = matrix(NA, J, nb_months ) 

for (i in 1: nb_months) {
  a = Bred_year3[ which(Bred_year3$Breeding.month == i), ]$Average.Wind.Speed
  l = length(a)
  C_y3[1:l,i] = a
}

C_y3

wind_speed = cbind(C_y1, C_y2, C_y3)
wind_speed
#save(wind_speed, file = "wind_speed.Rdata")

# Visibility
C_y1 = matrix(NA,J , nb_months ) 
head(Bred_year1)

a = Bred_year1[ which(Bred_year1$Breeding.month == 7), ]
a$Visability
dim(a)

for (i in 1: nb_months) {
  a = Bred_year1[ which(Bred_year1$Breeding.month == i), ]$Visability
  l = length(a)
  C_y1[1:l,i] = a
}

C_y1 # rows are sampling occasions and columns are months


# Breeding year 2
C_y2 = matrix(NA, J, nb_months ) 

for (i in 1: nb_months) {
  a = Bred_year2[ which(Bred_year2$Breeding.month == i), ]$Visability
  l = length(a)
  C_y2[1:l,i] = a
}

C_y1; C_y2


# Breeding year 3
C_y3 = matrix(NA, J, nb_months ) 

for (i in 1: nb_months) {
  a = Bred_year3[ which(Bred_year3$Breeding.month == i), ]$Visability
  l = length(a)
  C_y3[1:l,i] = a
}

C_y3

Visability = cbind(C_y1, C_y2, C_y3)
Visability

# Temp
C_y1 = matrix(NA,J , nb_months ) 
head(Bred_year1)

a = Bred_year1[ which(Bred_year1$Breeding.month == 7), ]
a$Median.Temp
dim(a)

for (i in 1: nb_months) {
  a = Bred_year1[ which(Bred_year1$Breeding.month == i), ]$Median.Temp
  l = length(a)
  C_y1[1:l,i] = a
}

C_y1 # rows are sampling occasions and columns are months


# Breeding year 2
C_y2 = matrix(NA, J, nb_months ) 

for (i in 1: nb_months) {
  a = Bred_year2[ which(Bred_year2$Breeding.month == i), ]$Median.Temp
  l = length(a)
  C_y2[1:l,i] = a
}

C_y1; C_y2


# Breeding year 3
C_y3 = matrix(NA, J, nb_months ) 

for (i in 1: nb_months) {
  a = Bred_year3[ which(Bred_year3$Breeding.month == i), ]$Median.Temp
  l = length(a)
  C_y3[1:l,i] = a
}

C_y3

Temp = cbind(C_y1, C_y2, C_y3)
Temp

## Humidity
C_y1 = matrix(NA,J , nb_months ) 
head(Bred_year1)

a = Bred_year1[ which(Bred_year1$Breeding.month == 7), ]
a$Average.relative.humidity
dim(a)

for (i in 1: nb_months) {
  a = Bred_year1[ which(Bred_year1$Breeding.month == i), ]$Average.relative.humidity
  l = length(a)
  C_y1[1:l,i] = a
}

C_y1 # rows are sampling occasions and columns are months


# Breeding year 2
C_y2 = matrix(NA, J, nb_months ) 

for (i in 1: nb_months) {
  a = Bred_year2[ which(Bred_year2$Breeding.month == i), ]$Average.relative.humidity
  l = length(a)
  C_y2[1:l,i] = a
}

C_y1; C_y2


# Breeding year 3
C_y3 = matrix(NA, J, nb_months ) 

for (i in 1: nb_months) {
  a = Bred_year3[ which(Bred_year3$Breeding.month == i), ]$Average.relative.humidity
  l = length(a)
  C_y3[1:l,i] = a
}

C_y3

Humidity = cbind(C_y1, C_y2, C_y3)
Humidity


## Rain
C_y1 = matrix(NA,J , nb_months ) 
head(Bred_year1)

a = Bred_year1[ which(Bred_year1$Breeding.month == 7), ]
a$Rain.Drizzle
dim(a)

for (i in 1: nb_months) {
  a = Bred_year1[ which(Bred_year1$Breeding.month == i), ]$Rain.Drizzle
  l = length(a)
  C_y1[1:l,i] = a
}

C_y1 # rows are sampling occasions and columns are months


# Breeding year 2
C_y2 = matrix(NA, J, nb_months ) 

for (i in 1: nb_months) {
  a = Bred_year2[ which(Bred_year2$Breeding.month == i), ]$Rain.Drizzle
  l = length(a)
  C_y2[1:l,i] = a
}

C_y1; C_y2


# Breeding year 3
C_y3 = matrix(NA, J, nb_months ) 

for (i in 1: nb_months) {
  a = Bred_year3[ which(Bred_year3$Breeding.month == i), ]$Rain.Drizzle
  l = length(a)
  C_y3[1:l,i] = a
}

C_y3

Rain = cbind(C_y1, C_y2, C_y3)
Rain

## Storm
C_y1 = matrix(NA,J , nb_months ) 
head(Bred_year1)

a = Bred_year1[ which(Bred_year1$Breeding.month == 7), ]
a$Storm.Thunder
dim(a)

for (i in 1: nb_months) {
  a = Bred_year1[ which(Bred_year1$Breeding.month == i), ]$Storm.Thunder
  l = length(a)
  C_y1[1:l,i] = a
}

C_y1 # rows are sampling occasions and columns are months


# Breeding year 2
C_y2 = matrix(NA, J, nb_months ) 

for (i in 1: nb_months) {
  a = Bred_year2[ which(Bred_year2$Breeding.month == i), ]$Storm.Thunder
  l = length(a)
  C_y2[1:l,i] = a
}

C_y1; C_y2


# Breeding year 3
C_y3 = matrix(NA, J, nb_months ) 

for (i in 1: nb_months) {
  a = Bred_year3[ which(Bred_year3$Breeding.month == i), ]$Storm.Thunder
  l = length(a)
  C_y3[1:l,i] = a
}

C_y3

Storm = cbind(C_y1, C_y2, C_y3)
Storm


## wather
C_y1 = matrix(NA,J , nb_months ) 
head(Bred_year1)

a = Bred_year1[ which(Bred_year1$Breeding.month == 7), ]
a$Obs.Weather
dim(a)

for (i in 1: nb_months) {
  a = Bred_year1[ which(Bred_year1$Breeding.month == i), ]$Obs.Weather
  l = length(a)
  C_y1[1:l,i] = a
}

C_y1 # rows are sampling occasions and columns are months


# Breeding year 2
C_y2 = matrix(NA, J, nb_months ) 

for (i in 1: nb_months) {
  a = Bred_year2[ which(Bred_year2$Breeding.month == i), ]$Obs.Weather
  l = length(a)
  C_y2[1:l,i] = a
}

C_y1; C_y2


# Breeding year 3
C_y3 = matrix(NA, J, nb_months ) 

for (i in 1: nb_months) {
  a = Bred_year3[ which(Bred_year3$Breeding.month == i), ]$Obs.Weather
  l = length(a)
  C_y3[1:l,i] = a
}

C_y3

weather = cbind(C_y1, C_y2, C_y3)
weather


X = array(c(Temp, Humidity, Visability, wind_speed, Rain, Storm, time), dim =c(J, T,7))
dim(X)
head(X[,,1])
