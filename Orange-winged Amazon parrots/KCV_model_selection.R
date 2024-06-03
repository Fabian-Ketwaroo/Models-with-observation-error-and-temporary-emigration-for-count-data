

# Rscript to analyse Orange parrots datasets to class of TE models developed
{
  library(nimble)
  nimbleOptions(MCMCusePredictiveDependenciesInCalculations = TRUE)
  nimbleOptions(MCMCorderPosteriorPredictiveSamplersLast = FALSE)
}

rm(list=ls())
#=============================================================================================================================================
# DP

#=============================================================================================================================================

# Cleaned data 
#load("~/Objective priors for N-mixture models/Parrots/Oparrots.Rdata")
load("Oparrots.Rdata")

{
  Oparrots
  J = dim(Oparrots)[1]
  nmonths = dim(Oparrots)[2]
  
  nweeks = J/2
  week1 = matrix(0, 2, nweeks )
  head(week1)
  
  for (t in 1:6) {
    
    week1[,t] = Oparrots[( (t*2-1):(t*2)) ,1]
    
  }
  
  week1
  Oparrots
  
  week2 = matrix(0, 2, nweeks )
  
  for (t in 1:6) {
    
    week2[,t] = Oparrots[( (t*2-1):(t*2)) ,2]
    
  }
  
  week2
  Oparrots[,2]
  
  w3 = matrix(0, 2, nweeks)
  i=3
  for (t in 1:6) {
    
    w3[,t] = Oparrots[( (t*2-1):(t*2)) ,i]
    
  }
  
  w3
  Oparrots[,i]
  
  w4 = matrix(0, 2, nweeks)
  i=4
  for (t in 1:6) {
    
    w4[,t] = Oparrots[( (t*2-1):(t*2)) ,i]
    
  }
  
  w4
  Oparrots[,i]
  
  w5 = matrix(0, 2, nweeks)
  i=5
  for (t in 1:6) {
    
    w5[,t] = Oparrots[( (t*2-1):(t*2)) ,i]
    
  }
  
  w5
  Oparrots[,i]
  
  
  w6 = matrix(0, 2, nweeks)
  i=6
  for (t in 1:6) {
    
    w6[,t] = Oparrots[( (t*2-1):(t*2)) ,i]
    
  }
  
  w6
  Oparrots[,i]
  
  
  w7 = matrix(0, 2, nweeks)
  i=7
  for (t in 1:6) {
    
    w7[,t] = Oparrots[( (t*2-1):(t*2)) ,i]
    
  }
  
  w7
  Oparrots[,i]
  
  w8 = matrix(0, 2, nweeks)
  i=8
  for (t in 1:6) {
    
    w8[,t] = Oparrots[( (t*2-1):(t*2)) ,i]
    
  }
  
  w8
  Oparrots[,i]
  
  
  w9 = matrix(0, 2, nweeks)
  i=9
  for (t in 1:6) {
    
    w9[,t] = Oparrots[( (t*2-1):(t*2)) ,i]
    
  }
  
  w9
  Oparrots[,i]
  
  w10 = matrix(0, 2, nweeks)
  i=10
  for (t in 1:6) {
    
    w10[,t] = Oparrots[( (t*2-1):(t*2)) ,i]
    
  }
  
  w10
  Oparrots[,i]
  
  w11 = matrix(0, 2, nweeks)
  i=11
  for (t in 1:6) {
    
    w11[,t] = Oparrots[( (t*2-1):(t*2)) ,i]
    
  }
  
  w11
  Oparrots[,i]
  
  w12 = matrix(0, 2, nweeks)
  i=12
  for (t in 1:6) {
    
    w12[,t] = Oparrots[( (t*2-1):(t*2)) ,i]
    
  }
  
  w12
  Oparrots[,i]
  
  
  w13 = matrix(0, 2, nweeks)
  i=13
  for (t in 1:6) {
    
    w13[,t] = Oparrots[( (t*2-1):(t*2)) ,i]
    
  }
  
  w13
  Oparrots[,i]
  
  
}

weeks_Opar = cbind(week1, week2,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12,w13)
#head(weeks_Opar)

nweeks =dim(weeks_Opar)[2]
index = numeric(nweeks)

for (t in 1:nweeks) {
  index[t] = sum(1- is.na(weeks_Opar[,t])) # create indicator where to end loop
  
}

#index
nweeks_index = index[index!=0] # number of sampling occasions in each month
nweeks_index
nweeks =length(nweeks_index)

#weeks_Opar # need to remvove na columns
dim(weeks_Opar)
Y = weeks_Opar[, index!=0]
#Y

# Covariates 
load("Oparrots_Xmatweekly.Rdata")
#head(X)
#X[,,1]
#Y # matches 

pos <- c(1,1,2,2,3,4) # The ‘pos’ vector gives the position of the factor corresponding to each component of beta
m <- c(2,2,2,2,1,1) # The ‘m’ vector captures the same size group for all levels belonging to the same categorical covariate. 


numVars = dim(X)[3]#ncov+1 #number of covariates (including the intercept)


win.data <- list(X =X,index = nweeks_index, pos= pos, m = m, nb_levels = length(pos), nweeks = nweeks,  numVars = numVars, K= nweeks)
str(win.data)

nimTE <- nimbleCode({
  
  # Priors
  lambda ~ dgamma(0.01,0.01)
  sigma.re ~ dunif(0,100)
  
  lambda2 ~ dgamma(0.001,0.001)
  # sigma.tu ~ dunif(0,10)
  tau.tu <- 1/sigma.re^2
  
  for (l in 1:nb_levels){
    # using the "pos" vector allows us to affect the same  prior
    # inclusion probability for all levels within the same categorical
    pi_g[l] ~ dbern(p_g)
    pi_g_pos[l] <- pi_g[pos[l]]
    
    tau_g[l] ~ dgamma((m[pos[l]] + 1) / 2,pow(lambda2,2) / 2)
    # pos is used to ensure the same tau_coef value for all levels
    # within the same categorical covariate
    tau_g_pos[l] <- tau_g[pos[l]]
    
    b[l] ~ dnorm(0, var =  tau.tu*(1/pow(tau_g_pos[l],2)))
    beta[l+1] <- (1 - pi_g_pos[l]) * 0 + pi_g_pos[l] * b[l]
    
  }
  
  #Intercept
  beta[1] ~ dnorm(0, sd = 2)
  p_g <- 0.5
  
  # DP set up
  for(i in 1:K) {
    gammatilde[i] ~ dgamma(2,0.1)
    betatilde[i] ~ dgamma(2,0.1)
  }
  
  z[1:nweeks] ~ dCRP(alpha, size = nweeks)
  alpha ~ dgamma(1, 1)
  
  
  M ~ dpois(lambda)                    # Super-population size M
  
  # Likelihood
  # Ecological model for true abundance
  for(t in 1:nweeks){                     # Loop over months
    
    theta[t] ~ dbeta(gammatilde[z[t]], betatilde[z[t]])
    N[t] ~ dbin(theta[t], M)            # 'Current' population size N
    
   linpred[1:index[t],t] <- (X[1:index[t],t, 1:numVars] %*% beta[1:numVars])[1:index[t],1]
    
    for(j in 1:index[t]){
      y[j,t] ~ dbin(p[j,t], N[t])
      re[j,t] ~ dnorm(0, sigma.re)
      logit(p[j,t]) <- linpred[j,t] + re[j,t]
      
    }
    
    
  }
  
  
  
})

Nst <- colSums(weeks_Opar, na.rm = TRUE)  #apply(Oparrots, 2, max, na.rm = TRUE)+100 # Inits for latent N
Nst = Nst[Nst !=0]
Mst <- max(Nst) + 500
TEmodel = nimbleModel(nimTE, constants = win.data,data = list(y =Y), inits = list(M = Mst,  N= Nst, lambda =  1000,lambda2 =1, tau_g = rgamma(length(pos), 1,1),b = rep(0, length(pos)), sigma.re = 1, pi_g = rep(1, length(pos)), re =  matrix(rnorm(2*nweeks), 2, nweeks) , w = rep(1, numVars),beta = rep(0, numVars),theta = rbeta(nweeks,0.5, 0.5), z = sample(1:win.data$K, size = nweeks, replace = TRUE), alpha = 1, gammatilde = rgamma(win.data$K, 0.5, 1) ,betatilde = rgamma(win.data$K, 0.5, 1) ))
TEmodel$calculate()

cTEmodel <- compileNimble(TEmodel)
TEconf = configureMCMC(TEmodel, enableWAIC = TRUE)



TEFoldFunction <- function(i){
  foldNodes_i <- paste0('y[ ,', i ,' ]')  # will return 'y[,i]' for i = 1 e.g.
  return(foldNodes_i)
}

# We define our own loss function as well.
# The function below will compute the root mean squared error.

RMSElossFunction <- function(simulatedDataValues, actualDataValues){
  dataLength <- length(simulatedDataValues) # simulatedDataValues is a vector
  SSE <- 0
  for(i in 1:dataLength){
    SSE <- SSE + (simulatedDataValues[i] - actualDataValues[i])^2
  }
  MSE <- SSE / dataLength
  RMSE <- sqrt(MSE)
  return(RMSE)
}


# Need to leave out co-variates 
DPcrossValOutput <- runCrossValidate(MCMCconfiguration = TEconf,
                                      k = nweeks,
                                      foldFunction = TEFoldFunction,
                                      lossFunction = RMSElossFunction,
                                      MCMCcontrol = list(niter = 105000, nburnin = 50000, thin = 20),nBootReps = NA,
                                      silent = F)



save(DPcrossValOutput, file ="OparrotsDP_KCV.Rdata")

####################################################################################################################

# RW1
rm(list = ls())
load("Oparrots.Rdata")

{
  Oparrots
  J = dim(Oparrots)[1]
  nmonths = dim(Oparrots)[2]
  
  nweeks = J/2
  week1 = matrix(0, 2, nweeks )
  head(week1)
  
  for (t in 1:6) {
    
    week1[,t] = Oparrots[( (t*2-1):(t*2)) ,1]
    
  }
  
  week1
  Oparrots
  
  week2 = matrix(0, 2, nweeks )
  
  for (t in 1:6) {
    
    week2[,t] = Oparrots[( (t*2-1):(t*2)) ,2]
    
  }
  
  week2
  Oparrots[,2]
  
  w3 = matrix(0, 2, nweeks)
  i=3
  for (t in 1:6) {
    
    w3[,t] = Oparrots[( (t*2-1):(t*2)) ,i]
    
  }
  
  w3
  Oparrots[,i]
  
  w4 = matrix(0, 2, nweeks)
  i=4
  for (t in 1:6) {
    
    w4[,t] = Oparrots[( (t*2-1):(t*2)) ,i]
    
  }
  
  w4
  Oparrots[,i]
  
  w5 = matrix(0, 2, nweeks)
  i=5
  for (t in 1:6) {
    
    w5[,t] = Oparrots[( (t*2-1):(t*2)) ,i]
    
  }
  
  w5
  Oparrots[,i]
  
  
  w6 = matrix(0, 2, nweeks)
  i=6
  for (t in 1:6) {
    
    w6[,t] = Oparrots[( (t*2-1):(t*2)) ,i]
    
  }
  
  w6
  Oparrots[,i]
  
  
  w7 = matrix(0, 2, nweeks)
  i=7
  for (t in 1:6) {
    
    w7[,t] = Oparrots[( (t*2-1):(t*2)) ,i]
    
  }
  
  w7
  Oparrots[,i]
  
  w8 = matrix(0, 2, nweeks)
  i=8
  for (t in 1:6) {
    
    w8[,t] = Oparrots[( (t*2-1):(t*2)) ,i]
    
  }
  
  w8
  Oparrots[,i]
  
  
  w9 = matrix(0, 2, nweeks)
  i=9
  for (t in 1:6) {
    
    w9[,t] = Oparrots[( (t*2-1):(t*2)) ,i]
    
  }
  
  w9
  Oparrots[,i]
  
  w10 = matrix(0, 2, nweeks)
  i=10
  for (t in 1:6) {
    
    w10[,t] = Oparrots[( (t*2-1):(t*2)) ,i]
    
  }
  
  w10
  Oparrots[,i]
  
  w11 = matrix(0, 2, nweeks)
  i=11
  for (t in 1:6) {
    
    w11[,t] = Oparrots[( (t*2-1):(t*2)) ,i]
    
  }
  
  w11
  Oparrots[,i]
  
  w12 = matrix(0, 2, nweeks)
  i=12
  for (t in 1:6) {
    
    w12[,t] = Oparrots[( (t*2-1):(t*2)) ,i]
    
  }
  
  w12
  Oparrots[,i]
  
  
  w13 = matrix(0, 2, nweeks)
  i=13
  for (t in 1:6) {
    
    w13[,t] = Oparrots[( (t*2-1):(t*2)) ,i]
    
  }
  
  w13
  Oparrots[,i]
  
  
}

weeks_Opar = cbind(week1, week2,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12,w13)
#head(weeks_Opar)

nweeks =dim(weeks_Opar)[2]
index = numeric(nweeks)

for (t in 1:nweeks) {
  index[t] = sum(1- is.na(weeks_Opar[,t])) # create indicator where to end loop
  
}

#index
nweeks_index = index[index!=0] # number of sampling occasions in each month
nweeks_index
nweeks =length(nweeks_index)

#weeks_Opar # need to remvove na columns
dim(weeks_Opar)
Y = weeks_Opar[, index!=0]
#Y

# Covariates 
load("Oparrots_Xmatweekly.Rdata")

# Cleaned data 
pos <- c(1,1,2,2,3,4) # The ‘pos’ vector gives the position of the factor corresponding to each component of beta
m <- c(2,2,2,2,1,1) # The ‘m’ vector captures the same size group for all levels belonging to the same categorical covariate. 


numVars = dim(X)[3]#ncov+1 #number of covariates (including the intercept)

T = nweeks
w= matrix(0,T,T)

for (i in 1:T) {
  if(i == 1) w[i,i+1] = 1
  if(i == T) w[i, i-1] =1
  else{
    w[i,c(i-1,i+1)] = 1
  }
}

#w
CarArg =  as.carAdjacency(w) # arguments for car distriubtion.

win.data <- list(X =X,index = nweeks_index, pos= pos, m = m, nb_levels = length(pos), nweeks = nweeks,  numVars = numVars, K= nweeks, L = length(CarArg$adj), adj = CarArg$adj, weights= CarArg$weights, num = CarArg$num)
str(win.data)

nimTE <- nimbleCode({
  
  # Priors
  lambda ~ dgamma(0.01,0.01)
  sigma.re ~ dunif(0,100)
  sigma ~ dunif(0,15)
  tau <- 1/sigma^2
  
  lambda2 ~ dgamma(0.001,0.001)
  # sigma.tu ~ dunif(0,10)
  tau.tu <- 1/sigma.re^2
  
  for (l in 1:nb_levels){
    # using the "pos" vector allows us to affect the same  prior
    # inclusion probability for all levels within the same categorical
    pi_g[l] ~ dbern(p_g)
    pi_g_pos[l] <- pi_g[pos[l]]
    
    tau_g[l] ~ dgamma((m[pos[l]] + 1) / 2,pow(lambda2,2) / 2)
    # pos is used to ensure the same tau_coef value for all levels
    # within the same categorical covariate
    tau_g_pos[l] <- tau_g[pos[l]]
    
    b[l] ~ dnorm(0, var =  tau.tu*(1/pow(tau_g_pos[l],2)))
    beta[l+1] <- (1 - pi_g_pos[l]) * 0 + pi_g_pos[l] * b[l]
    
  }
  
  #Intercept
  beta[1] ~ dnorm(0, sd = 2)
  p_g <- 0.5
  
  
  M ~ dpois(lambda)                    # Super-population size M
  
  # RW1 approximated by ICAR
  s[1:nweeks] ~ dcar_normal(adj[1:L], weights[1:L], num[1:nweeks], tau, zero_mean = 0)
  
  
  # Likelihood
  # Ecological model for true abundance
  for(t in 1:nweeks){                     # Loop over months
    
    logit(theta[t]) <- s[t]
    N[t] ~ dbin(theta[t], M)            # 'Current' population size N
    
    linpred[1:index[t],t] <- (X[1:index[t],t, 1:numVars] %*% beta[1:numVars])[1:index[t],1]
    
    for(j in 1:index[t]){
      y[j,t] ~ dbin(p[j,t], N[t])
      re[j,t] ~ dnorm(0, sigma.re)
      logit(p[j,t]) <- linpred[j,t] + re[j,t]
      
    }
    
    
  }
  
  
  
})

Nst <- colSums(weeks_Opar, na.rm = TRUE)  #apply(Oparrots, 2, max, na.rm = TRUE)+100 # Inits for latent N
Nst = Nst[Nst !=0]
Mst <- max(Nst) + 500
TEmodel = nimbleModel(nimTE, constants = win.data, data = list(y =Y), inits = list(sigma.re=1,pi_g = rep(1, length(pos)), lambda2 =1, tau_g = rgamma(length(pos), 1,1),b = rep(0, length(pos)), re =  matrix(rnorm(2*nweeks), 2, nweeks), M = Mst,  N= Nst, lambda =  rgamma(1,100,0.1), beta = rep(0, numVars), s = rnorm(nweeks, 0,1), sigma = runif(1) )  )
TEmodel$calculate()

cTEmodel <- compileNimble(TEmodel)
TEconf = configureMCMC(TEmodel, enableWAIC = TRUE)


TEFoldFunction <- function(i){
  foldNodes_i <- paste0('y[ ,', i ,' ]')  # will return 'y[,i]' for i = 1 e.g.
  return(foldNodes_i)
}

# We define our own loss function as well.
# The function below will compute the root mean squared error.

RMSElossFunction <- function(simulatedDataValues, actualDataValues){
  dataLength <- length(simulatedDataValues) # simulatedDataValues is a vector
  SSE <- 0
  for(i in 1:dataLength){
    SSE <- SSE + (simulatedDataValues[i] - actualDataValues[i])^2
  }
  MSE <- SSE / dataLength
  RMSE <- sqrt(MSE)
  return(RMSE)
}


# Need to leave out co-variates 
RW1crossValOutput <- runCrossValidate(MCMCconfiguration = TEconf,
                                     k = nweeks,
                                     foldFunction = TEFoldFunction,
                                     lossFunction = RMSElossFunction,
                                     MCMCcontrol = list(niter = 105000, nburnin = 50000, thin = 20),nBootReps = NA,
                                     silent = TRUE)



save(RW1crossValOutput, file ="OparrotsRW1_KCV.Rdata")



####################################################################################################################

# RW2
rm(list = ls())
# Cleaned data 
load("Oparrots.Rdata")

{
  Oparrots
  J = dim(Oparrots)[1]
  nmonths = dim(Oparrots)[2]
  
  nweeks = J/2
  week1 = matrix(0, 2, nweeks )
  head(week1)
  
  for (t in 1:6) {
    
    week1[,t] = Oparrots[( (t*2-1):(t*2)) ,1]
    
  }
  
  week1
  Oparrots
  
  week2 = matrix(0, 2, nweeks )
  
  for (t in 1:6) {
    
    week2[,t] = Oparrots[( (t*2-1):(t*2)) ,2]
    
  }
  
  week2
  Oparrots[,2]
  
  w3 = matrix(0, 2, nweeks)
  i=3
  for (t in 1:6) {
    
    w3[,t] = Oparrots[( (t*2-1):(t*2)) ,i]
    
  }
  
  w3
  Oparrots[,i]
  
  w4 = matrix(0, 2, nweeks)
  i=4
  for (t in 1:6) {
    
    w4[,t] = Oparrots[( (t*2-1):(t*2)) ,i]
    
  }
  
  w4
  Oparrots[,i]
  
  w5 = matrix(0, 2, nweeks)
  i=5
  for (t in 1:6) {
    
    w5[,t] = Oparrots[( (t*2-1):(t*2)) ,i]
    
  }
  
  w5
  Oparrots[,i]
  
  
  w6 = matrix(0, 2, nweeks)
  i=6
  for (t in 1:6) {
    
    w6[,t] = Oparrots[( (t*2-1):(t*2)) ,i]
    
  }
  
  w6
  Oparrots[,i]
  
  
  w7 = matrix(0, 2, nweeks)
  i=7
  for (t in 1:6) {
    
    w7[,t] = Oparrots[( (t*2-1):(t*2)) ,i]
    
  }
  
  w7
  Oparrots[,i]
  
  w8 = matrix(0, 2, nweeks)
  i=8
  for (t in 1:6) {
    
    w8[,t] = Oparrots[( (t*2-1):(t*2)) ,i]
    
  }
  
  w8
  Oparrots[,i]
  
  
  w9 = matrix(0, 2, nweeks)
  i=9
  for (t in 1:6) {
    
    w9[,t] = Oparrots[( (t*2-1):(t*2)) ,i]
    
  }
  
  w9
  Oparrots[,i]
  
  w10 = matrix(0, 2, nweeks)
  i=10
  for (t in 1:6) {
    
    w10[,t] = Oparrots[( (t*2-1):(t*2)) ,i]
    
  }
  
  w10
  Oparrots[,i]
  
  w11 = matrix(0, 2, nweeks)
  i=11
  for (t in 1:6) {
    
    w11[,t] = Oparrots[( (t*2-1):(t*2)) ,i]
    
  }
  
  w11
  Oparrots[,i]
  
  w12 = matrix(0, 2, nweeks)
  i=12
  for (t in 1:6) {
    
    w12[,t] = Oparrots[( (t*2-1):(t*2)) ,i]
    
  }
  
  w12
  Oparrots[,i]
  
  
  w13 = matrix(0, 2, nweeks)
  i=13
  for (t in 1:6) {
    
    w13[,t] = Oparrots[( (t*2-1):(t*2)) ,i]
    
  }
  
  w13
  Oparrots[,i]
  
  
}

weeks_Opar = cbind(week1, week2,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12,w13)
#head(weeks_Opar)

nweeks =dim(weeks_Opar)[2]
index = numeric(nweeks)

for (t in 1:nweeks) {
  index[t] = sum(1- is.na(weeks_Opar[,t])) # create indicator where to end loop
  
}

#index
nweeks_index = index[index!=0] # number of sampling occasions in each month
nweeks_index
nweeks =length(nweeks_index)

#weeks_Opar # need to remvove na columns
dim(weeks_Opar)
Y = weeks_Opar[, index!=0]
#Y

# Covariates 
load("Oparrots_Xmatweekly.Rdata")

pos <- c(1,1,2,2,3,4) # The ‘pos’ vector gives the position of the factor corresponding to each component of beta
m <- c(2,2,2,2,1,1) # The ‘m’ vector captures the same size group for all levels belonging to the same categorical covariate. 


numVars = dim(X)[3]#ncov+1 #number of covariates (including the intercept)

T = nweeks
w= matrix(0,T,T)

w[1, c(1+1,1+2)] = c(2, -1)
w[T, c(T-1, T-2)] = c(2,-1)
w[2, c(2-1,2+1, 2+2)] = c(2, 4,-1)

for (i in 3:(T-2)) {
  w[i,c(i-2,i-1,i+1,i+2)] = c(-1,4,4,-1)
} 

w[T-1, c(T-1-2, T-1-1, T)] = c(-1,4,2)

#w
CarArg =  as.carAdjacency(w) # arguments for car distriubtion.


win.data <- list(X =X,index = nweeks_index, pos= pos, m = m, nb_levels = length(pos), nweeks = nweeks,  numVars = numVars, K= nweeks, L = length(CarArg$adj), adj = CarArg$adj, weights= CarArg$weights, num = CarArg$num)
str(win.data)

nimTE <- nimbleCode({
  
  # Priors
  lambda ~ dgamma(0.01,0.01)
  sigma.re ~ dunif(0,100)
  sigma ~ dunif(0,15)
  tau <- 1/sigma^2
  
  lambda2 ~ dgamma(0.001,0.001)
  # sigma.tu ~ dunif(0,10)
  tau.tu <- 1/sigma.re^2
  
  for (l in 1:nb_levels){
    # using the "pos" vector allows us to affect the same  prior
    # inclusion probability for all levels within the same categorical
    pi_g[l] ~ dbern(p_g)
    pi_g_pos[l] <- pi_g[pos[l]]
    
    tau_g[l] ~ dgamma((m[pos[l]] + 1) / 2,pow(lambda2,2) / 2)
    # pos is used to ensure the same tau_coef value for all levels
    # within the same categorical covariate
    tau_g_pos[l] <- tau_g[pos[l]]
    
    b[l] ~ dnorm(0, var =  tau.tu*(1/pow(tau_g_pos[l],2)))
    beta[l+1] <- (1 - pi_g_pos[l]) * 0 + pi_g_pos[l] * b[l]
    
  }
  
  #Intercept
  beta[1] ~ dnorm(0, sd = 2)
  p_g <- 0.5
  
  
  M ~ dpois(lambda)                    # Super-population size M
  
  # RW1 approximated by ICAR
  s[1:nweeks] ~ dcar_normal(adj[1:L], weights[1:L], num[1:nweeks], tau, zero_mean = 0)
  
  
  # Likelihood
  # Ecological model for true abundance
  for(t in 1:nweeks){                     # Loop over months
    
    logit(theta[t]) <- s[t]
    N[t] ~ dbin(theta[t], M)            # 'Current' population size N
    
    linpred[1:index[t],t] <- (X[1:index[t],t, 1:numVars] %*% beta[1:numVars])[1:index[t],1]
    
    for(j in 1:index[t]){
      y[j,t] ~ dbin(p[j,t], N[t])
      re[j,t] ~ dnorm(0, sigma.re)
      logit(p[j,t]) <- linpred[j,t] + re[j,t]
      
    }
    
    
  }
  
  
  
})

Nst <- colSums(weeks_Opar, na.rm = TRUE)  #apply(Oparrots, 2, max, na.rm = TRUE)+100 # Inits for latent N
Nst = Nst[Nst !=0]
Mst <- max(Nst) + 500
TEmodel = nimbleModel(nimTE, constants = win.data, data = list(y =Y), inits = list(sigma.re=1,pi_g = rep(1, length(pos)), lambda2 =1, tau_g = rgamma(length(pos), 1,1),b = rep(0, length(pos)), re =  matrix(rnorm(2*nweeks), 2, nweeks), M = Mst,  N= Nst, lambda =  rgamma(1,100,0.1), beta = rep(0, numVars), s = rnorm(nweeks, 0,1), sigma = runif(1) )  )
TEmodel$calculate()

cTEmodel <- compileNimble(TEmodel)
TEconf = configureMCMC(TEmodel, enableWAIC = TRUE)



TEFoldFunction <- function(i){
  foldNodes_i <- paste0('y[ ,', i ,' ]')  # will return 'y[,i]' for i = 1 e.g.
  return(foldNodes_i)
}

# We define our own loss function as well.
# The function below will compute the root mean squared error.

RMSElossFunction <- function(simulatedDataValues, actualDataValues){
  dataLength <- length(simulatedDataValues) # simulatedDataValues is a vector
  SSE <- 0
  for(i in 1:dataLength){
    SSE <- SSE + (simulatedDataValues[i] - actualDataValues[i])^2
  }
  MSE <- SSE / dataLength
  RMSE <- sqrt(MSE)
  return(RMSE)
}


# Need to leave out co-variates 
RW2crossValOutput <- runCrossValidate(MCMCconfiguration = TEconf,
                                      k = nweeks,
                                      foldFunction = TEFoldFunction,
                                      lossFunction = RMSElossFunction,
                                      MCMCcontrol = list(niter = 105000, nburnin = 50000, thin = 20),nBootReps = NA,
                                      silent = TRUE)



save(RW2crossValOutput, file ="OparrotsRW2_KCV.Rdata")





####################################################################################################################

# AR1
rm(list = ls())
# Cleaned data 
load("Oparrots.Rdata")

{
  Oparrots
  J = dim(Oparrots)[1]
  nmonths = dim(Oparrots)[2]
  
  nweeks = J/2
  week1 = matrix(0, 2, nweeks )
  head(week1)
  
  for (t in 1:6) {
    
    week1[,t] = Oparrots[( (t*2-1):(t*2)) ,1]
    
  }
  
  week1
  Oparrots
  
  week2 = matrix(0, 2, nweeks )
  
  for (t in 1:6) {
    
    week2[,t] = Oparrots[( (t*2-1):(t*2)) ,2]
    
  }
  
  week2
  Oparrots[,2]
  
  w3 = matrix(0, 2, nweeks)
  i=3
  for (t in 1:6) {
    
    w3[,t] = Oparrots[( (t*2-1):(t*2)) ,i]
    
  }
  
  w3
  Oparrots[,i]
  
  w4 = matrix(0, 2, nweeks)
  i=4
  for (t in 1:6) {
    
    w4[,t] = Oparrots[( (t*2-1):(t*2)) ,i]
    
  }
  
  w4
  Oparrots[,i]
  
  w5 = matrix(0, 2, nweeks)
  i=5
  for (t in 1:6) {
    
    w5[,t] = Oparrots[( (t*2-1):(t*2)) ,i]
    
  }
  
  w5
  Oparrots[,i]
  
  
  w6 = matrix(0, 2, nweeks)
  i=6
  for (t in 1:6) {
    
    w6[,t] = Oparrots[( (t*2-1):(t*2)) ,i]
    
  }
  
  w6
  Oparrots[,i]
  
  
  w7 = matrix(0, 2, nweeks)
  i=7
  for (t in 1:6) {
    
    w7[,t] = Oparrots[( (t*2-1):(t*2)) ,i]
    
  }
  
  w7
  Oparrots[,i]
  
  w8 = matrix(0, 2, nweeks)
  i=8
  for (t in 1:6) {
    
    w8[,t] = Oparrots[( (t*2-1):(t*2)) ,i]
    
  }
  
  w8
  Oparrots[,i]
  
  
  w9 = matrix(0, 2, nweeks)
  i=9
  for (t in 1:6) {
    
    w9[,t] = Oparrots[( (t*2-1):(t*2)) ,i]
    
  }
  
  w9
  Oparrots[,i]
  
  w10 = matrix(0, 2, nweeks)
  i=10
  for (t in 1:6) {
    
    w10[,t] = Oparrots[( (t*2-1):(t*2)) ,i]
    
  }
  
  w10
  Oparrots[,i]
  
  w11 = matrix(0, 2, nweeks)
  i=11
  for (t in 1:6) {
    
    w11[,t] = Oparrots[( (t*2-1):(t*2)) ,i]
    
  }
  
  w11
  Oparrots[,i]
  
  w12 = matrix(0, 2, nweeks)
  i=12
  for (t in 1:6) {
    
    w12[,t] = Oparrots[( (t*2-1):(t*2)) ,i]
    
  }
  
  w12
  Oparrots[,i]
  
  
  w13 = matrix(0, 2, nweeks)
  i=13
  for (t in 1:6) {
    
    w13[,t] = Oparrots[( (t*2-1):(t*2)) ,i]
    
  }
  
  w13
  Oparrots[,i]
  
  
}

weeks_Opar = cbind(week1, week2,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12,w13)
#head(weeks_Opar)

nweeks =dim(weeks_Opar)[2]
index = numeric(nweeks)

for (t in 1:nweeks) {
  index[t] = sum(1- is.na(weeks_Opar[,t])) # create indicator where to end loop
  
}

#index
nweeks_index = index[index!=0] # number of sampling occasions in each month
nweeks_index
nweeks =length(nweeks_index)

#weeks_Opar # need to remvove na columns
dim(weeks_Opar)
Y = weeks_Opar[, index!=0]
#Y

# Covariates 
#load("C:/Users/Fabian Ketwaroo/Documents/Objective priors for N-mixture models/Parrots/Oparrots_Xmatweekly.Rdata")
load("Oparrots_Xmatweekly.Rdata")

pos <- c(1,1,2,2,3,4) # The ‘pos’ vector gives the position of the factor corresponding to each component of beta
m <- c(2,2,2,2,1,1) # The ‘m’ vector captures the same size group for all levels belonging to the same categorical covariate. 


numVars = dim(X)[3]#ncov+1 #number of covariates (including the intercept)


win.data <- list(X =X,index = nweeks_index, pos= pos, m = m, nb_levels = length(pos), nweeks = nweeks,  numVars = numVars)
str(win.data)

nimTE <- nimbleCode({
  
  # Priors
  lambda ~ dgamma(0.01,0.01)
  sigma.re ~ dunif(0,100)
  sigma ~ dunif(0,15)
  tau <- 1/sigma^2
  
  sigma.e ~ dunif(0,10)
  alpha ~ dnorm(0,1)#dgamma(0.01,0.01)
  rho ~ dunif(-1,1)
  
  
  lambda2 ~ dgamma(0.001,0.001)
  # sigma.tu ~ dunif(0,10)
  tau.tu <- 1/sigma.re^2
  
  for (l in 1:nb_levels){
    # using the "pos" vector allows us to affect the same  prior
    # inclusion probability for all levels within the same categorical
    pi_g[l] ~ dbern(p_g)
    pi_g_pos[l] <- pi_g[pos[l]]
    
    tau_g[l] ~ dgamma((m[pos[l]] + 1) / 2,pow(lambda2,2) / 2)
    # pos is used to ensure the same tau_coef value for all levels
    # within the same categorical covariate
    tau_g_pos[l] <- tau_g[pos[l]]
    
    b[l] ~ dnorm(0, var =  tau.tu*(1/pow(tau_g_pos[l],2)))
    beta[l+1] <- (1 - pi_g_pos[l]) * 0 + pi_g_pos[l] * b[l]
    
  }
  
  #Intercept
  beta[1] ~ dnorm(0, sd = 2)
  p_g <- 0.5
  
  
  M ~ dpois(lambda)                    # Super-population size M
  
  # AR1
  s[1] ~ dnorm(0, sd =(tau*(1-rho^2))^(-0.5) )
  
  for (t in 2:nmonths) {
    s[t] ~ dnorm(rho*s[t-1], sd = sigma.e)
  }
  
  # Likelihood
  # Ecological model for true abundance
  for(t in 1:nweeks){                     # Loop over months
    
    logit(theta[t]) <- alpha + s[t] 
    N[t] ~ dbin(theta[t], M)            # 'Current' population size N
    
    linpred[1:index[t],t] <- (X[1:index[t],t, 1:numVars] %*% beta[1:numVars])[1:index[t],1]
    
    for(j in 1:index[t]){
      y[j,t] ~ dbin(p[j,t], N[t])
      re[j,t] ~ dnorm(0, sigma.re)
      logit(p[j,t]) <- linpred[j,t] + re[j,t]
      
    }
    
    
  }
  
  
  
})

Nst <- colSums(weeks_Opar, na.rm = TRUE)  #apply(Oparrots, 2, max, na.rm = TRUE)+100 # Inits for latent N
Nst = Nst[Nst !=0]
Mst <- max(Nst) + 500
TEmodel = nimbleModel(nimTE, constants = win.data, data = list(y =Y), inits = list(sigma.e = 1, rho = runif(1, -1,1), alpha = 1,sigma.re=1,pi_g = rep(1, length(pos)), lambda2 =1, tau_g = rgamma(length(pos), 1,1),b = rep(0, length(pos)), re =  matrix(rnorm(2*nweeks), 2, nweeks), M = Mst,  N= Nst, lambda =  rgamma(1,100,0.1), beta = rep(0, numVars), s = rnorm(nweeks, 0,1), sigma = runif(1) )  )
TEmodel$calculate()

cTEmodel <- compileNimble(TEmodel)
TEconf = configureMCMC(TEmodel, enableWAIC = TRUE)

TEFoldFunction <- function(i){
  foldNodes_i <- paste0('y[ ,', i ,' ]')  # will return 'y[,i]' for i = 1 e.g.
  return(foldNodes_i)
}

# We define our own loss function as well.
# The function below will compute the root mean squared error.

RMSElossFunction <- function(simulatedDataValues, actualDataValues){
  dataLength <- length(simulatedDataValues) # simulatedDataValues is a vector
  SSE <- 0
  for(i in 1:dataLength){
    SSE <- SSE + (simulatedDataValues[i] - actualDataValues[i])^2
  }
  MSE <- SSE / dataLength
  RMSE <- sqrt(MSE)
  return(RMSE)
}


# Need to leave out co-variates 
AR1crossValOutput <- runCrossValidate(MCMCconfiguration = TEconf,
                                      k = nweeks,
                                      foldFunction = TEFoldFunction,
                                      lossFunction = RMSElossFunction,
                                      MCMCcontrol = list(niter = 105000, nburnin = 50000, thin = 20),nBootReps = NA,
                                      silent = TRUE)



save(AR1crossValOutput, file ="OparrotsAR1_KCV.Rdata")
