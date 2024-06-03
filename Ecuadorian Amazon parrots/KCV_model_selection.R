
{
  library(nimble)
  nimbleOptions(MCMCusePredictiveDependenciesInCalculations = TRUE)
  nimbleOptions(MCMCorderPosteriorPredictiveSamplersLast = FALSE)

}



# Cleaned data 
load("cparrots.Rdata")

cparrots
J = dim(cparrots)[1]
nmonths = dim(cparrots)[2]

index = numeric(nmonths)

for (t in 1:nmonths) {
  index[t] = sum(1- is.na(cparrots[,t])) # create indicator where to end loop
  
}

index # number of sampling occasions in each month

Y =cparrots # need to remvove
Y

# Covariates 
load("Bparrots_Xmat.Rdata")
head(X)
dim(X)
pos <- c(1,2,3,4,5,6,7,8,8,8) # The ‘pos’ vector gives the position of the factor corresponding to each component of beta
m <- c(1,1,1,1,1,1,1,3,3,3) # The ‘m’ vector captures the same size group for all levels belonging to the same categorical covariate. 

numVars = dim(X)[3]

yN_index = rep(1:3, each = 12)

win.data <- list( index = index, pos = pos, m = m, nb_levels = length(pos), nmonths = nmonths,  numVars = numVars, K=36, Y=3, yN_index = yN_index)
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
    gammatilde[i] ~ dgamma(2,0.1)#dgamma(100,4)
    betatilde[i] ~ dgamma(2,0.1)#dgamma(100,4)
  }
  
  z[1:nmonths] ~ dCRP(alpha, size = nmonths)
  alpha ~ dgamma(1, 1)
  
  
  M ~ dpois(lambda)                    # Super-population size M
  
  # Yearly abundance
  for (y1 in 1:Y) {
    
    delta[y1] ~ dunif(0,1)
    yN[y1] ~ dbin(delta[y1], M)
    
  }
  
  # Likelihood
  # Ecological model for true abundance
  for(t in 1:nmonths){                     # Loop over months
    
    theta[t] ~ dbeta(gammatilde[z[t]], betatilde[z[t]])
    N[t] ~ dbin(theta[t], yN[yN_index[t]])            # 'Current' population size N
    
    linpred[1:index[t],t] <- (X[1:index[t],t, 1:numVars] %*% beta[1:numVars])[1:index[t],1]
    
    for (j in 1:index[t]) {
      y[j,t] ~ dbin(p[j,t], N[t])
      re[j,t] ~ dnorm(0, sigma.re)
      logit(p[j,t]) <- linpred[j,t] + re[j,t]
    }
    
    
  }
  
  
})

Nst <- apply(cparrots, 2, max, na.rm = TRUE)+100 # Inits for latent N
Mst <- max(Nst) + 500
Kst <-  rbinom(3, Mst, 0.5)
TEmodel = nimbleModel(nimTE, constants = win.data,data = list(y =Y, X=X), inits = list(pi_g = rep(1, length(pos)),lambda2 =1, tau_g = rgamma(length(pos), 0.01,0.01), b = rep(0, length(pos)), M = Mst,  N= Nst, yN = Kst, sigma.re = 1, re = matrix(rnorm(J*nmonths, 0,1), J, nmonths), lambda =  runif(1, 100, 500),theta = runif(1,0,1), beta = rep(0, numVars), delta = runif(3,0,1),theta = rbeta(nmonths,0.5, 0.5), z = sample(1:win.data$K, size = nmonths, replace = TRUE), alpha = 1, gammatilde = rgamma(win.data$K, 0.5, 1) ,betatilde = rgamma(win.data$K, 0.5, 1) ))
TEmodel$calculate()

cTEmodel <- compileNimble(TEmodel)
TEconf = configureMCMC(TEmodel)



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
                                     k = nmonths,
                                     foldFunction = TEFoldFunction,
                                     lossFunction = RMSElossFunction,
                                     MCMCcontrol = list(niter = 135000, nburnin = 35000, thin = 20), nBootReps = NA,
                                     silent = TRUE)



save(DPcrossValOutput, file ="BparrotsDP_KCV.Rdata")




#==================================================================================================================================================
#==================================================================================================================================================
rm(list = ls())
# Rw1
library(nimble)
nimbleOptions(MCMCusePredictiveDependenciesInCalculations = TRUE)
nimbleOptions(MCMCorderPosteriorPredictiveSamplersLast = FALSE)


# Cleaned data 
load("cparrots.Rdata")
#load("~/Objective priors for N-mixture models/Parrots/parrots_time.Rdata")

cparrots
J = dim(cparrots)[1]
nmonths = dim(cparrots)[2]

index = numeric(nmonths)

for (t in 1:nmonths) {
  index[t] = sum(1- is.na(cparrots[,t])) # create indicator where to end loop
  
}

index # number of sampling occasions in each month

Y =cparrots # need to remvove
Y

# Covariates 
load("Bparrots_Xmat.Rdata")
head(X)
dim(X)
pos <- c(1,2,3,4,5,6,7,8,8,8) # The ‘pos’ vector gives the position of the factor corresponding to each component of beta
m <- c(1,1,1,1,1,1,1,3,3,3) # The ‘m’ vector captures the same size group for all levels belonging to the same categorical covariate. 

numVars = dim(X)[3]

yN_index = rep(1:3, each = 12)

T = nmonths
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

win.data <- list(Y=3,yN_index = yN_index, index = index, pos = pos, m = m, nb_levels = length(pos),  nmonths = T,  numVars = numVars, L = length(CarArg$adj), adj = CarArg$adj, weights= CarArg$weights, num = CarArg$num)
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
  
 
  
  # RW1 approximated by ICAR
  s[1:nmonths] ~ dcar_normal(adj[1:L], weights[1:L], num[1:nmonths], tau, zero_mean = 0)
  
  
  
  M ~ dpois(lambda)                    # Super-population size M
  
  # Yearly abundance
  for (y1 in 1:Y) {
    
    delta[y1] ~ dunif(0,1)
    yN[y1] ~ dbin(delta[y1], M)
    
  }
  
  for(t in 1:nmonths){                     # Loop over months
    
    logit(theta[t]) <- s[t] #~ dbeta(gammatilde[z[t]], betatilde[z[t]])
    N[t] ~ dbin(theta[t], yN[yN_index[t]])   
    
    linpred[1:index[t],t] <- (X[1:index[t],t, 1:numVars] %*% beta[1:numVars])[1:index[t],1]
    
    for (j in 1:index[t]) {
      y[j,t] ~ dbin(p[j,t], N[t])
      re[j,t] ~ dnorm(0, sigma.re)
      logit(p[j,t]) <- linpred[j,t] + re[j,t]
    }
    
  }
  
 
})

Nst <- apply(cparrots, 2, max, na.rm = TRUE)+100 # Inits for latent N
Mst <- max(Nst) + 500
Kst <-  rbinom(3, Mst, 0.5)
TEmodel = nimbleModel(nimTE, constants = win.data, data = list(y =Y, X=X), inits = list(delta = runif(3,0,1), yN =Kst, pi_g = rep(1, length(pos)), lambda2 =1, tau_g = rgamma(length(pos), 1,1), b = rep(0, length(pos)), sigma.re=1, re = matrix(rnorm(J*nmonths, 0,1), J, nmonths), M = Mst,  N= Nst,  lambda =  rgamma(1,10,0.1), beta = rep(0, numVars), s = rnorm(T, 0,1), sigma = runif(1) )  )
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
                                     k = nmonths,
                                     foldFunction = TEFoldFunction,
                                     lossFunction = RMSElossFunction,
                                     MCMCcontrol = list(niter = 135000, nburnin = 35000, thin = 20), nBootReps = NA, silent = F)



save(RW1crossValOutput, file ="BparrotsRW1_KCV.Rdata")



#==================================================================================================================================================

rm(list = ls())
# Rw2
library(nimble)


# Becca Parrots 
# Cleaned data 
load("cparrots.Rdata")
#load("~/Objective priors for N-mixture models/Parrots/parrots_time.Rdata")

cparrots
J = dim(cparrots)[1]
nmonths = dim(cparrots)[2]

index = numeric(nmonths)

for (t in 1:nmonths) {
  index[t] = sum(1- is.na(cparrots[,t])) # create indicator where to end loop
  
}

index # number of sampling occasions in each month

Y =cparrots # need to remvove
Y

# Covariates 
load("Bparrots_Xmat.Rdata")
head(X)
dim(X)
pos <- c(1,2,3,4,5,6,7,8,8,8) # The ‘pos’ vector gives the position of the factor corresponding to each component of beta
m <- c(1,1,1,1,1,1,1,3,3,3) # The ‘m’ vector captures the same size group for all levels belonging to the same categorical covariate. 

numVars = dim(X)[3]

yN_index = rep(1:3, each = 12)

T = nmonths
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

win.data <- list(Y=3,yN_index = yN_index, index = index, pos = pos, m = m, nb_levels = length(pos),  nmonths = T,  numVars = numVars, L = length(CarArg$adj), adj = CarArg$adj, weights= CarArg$weights, num = CarArg$num)
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
  
  
  
  # RW1 approximated by ICAR
  s[1:nmonths] ~ dcar_normal(adj[1:L], weights[1:L], num[1:nmonths], tau, zero_mean = 0)
  
  
  
  M ~ dpois(lambda)                    # Super-population size M
  
  # Yearly abundance
  for (y1 in 1:Y) {
    
    delta[y1] ~ dunif(0,1)
    yN[y1] ~ dbin(delta[y1], M)
    
  }
  
  for(t in 1:nmonths){                     # Loop over months
    
    logit(theta[t]) <- s[t] #~ dbeta(gammatilde[z[t]], betatilde[z[t]])
    N[t] ~ dbin(theta[t], yN[yN_index[t]])   
    
    linpred[1:index[t],t] <- (X[1:index[t],t, 1:numVars] %*% beta[1:numVars])[1:index[t],1]
    
    for (j in 1:index[t]) {
      y[j,t] ~ dbin(p[j,t], N[t])
      re[j,t] ~ dnorm(0, sigma.re)
      logit(p[j,t]) <- linpred[j,t] + re[j,t]
    }
    
  }
  
  
})

Nst <- apply(cparrots, 2, max, na.rm = TRUE)+100 # Inits for latent N
Mst <- max(Nst) + 500
Kst <-  rbinom(3, Mst, 0.5)
TEmodel = nimbleModel(nimTE, constants = win.data, data = list(y =Y, X=X), inits = list(delta = runif(3,0,1), yN =Kst, pi_g = rep(1, length(pos)), lambda2 =1, tau_g = rgamma(length(pos), 1,1), b = rep(0, length(pos)), sigma.re=1, re = matrix(rnorm(J*nmonths, 0,1), J, nmonths), M = Mst,  N= Nst,  lambda =  rgamma(1,10,0.1), beta = rep(0, numVars), s = rnorm(T, 0,1), sigma = runif(1) )  )
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
                                      k = nmonths,
                                      foldFunction = TEFoldFunction,
                                      lossFunction = RMSElossFunction,
                                      MCMCcontrol = list(niter = 135000, nburnin = 35000, thin = 20), nBootReps = NA,
                                      silent = TRUE)



save(RW2crossValOutput, file ="BparrotsRW2_KCV.Rdata")
#=================================================================================================================================================

rm(list = ls())
# Cor
library(nimble)


# Cleaned data 
load("cparrots.Rdata")


cparrots
J = dim(cparrots)[1]
nmonths = dim(cparrots)[2]

index = numeric(nmonths)

for (t in 1:nmonths) {
  index[t] = sum(1- is.na(cparrots[,t])) # create indicator where to end loop
  
}

index # number of sampling occasions in each month

Y =cparrots # need to remvove
Y

# Covariates 
load("Bparrots_Xmat.Rdata")
head(X)
dim(X)
pos <- c(1,2,3,4,5,6,7,8,8,8) # The ‘pos’ vector gives the position of the factor corresponding to each component of beta
m <- c(1,1,1,1,1,1,1,3,3,3) # The ‘m’ vector captures the same size group for all levels belonging to the same categorical covariate. 

numVars = dim(X)[3]

yN_index = rep(1:3, each = 12)


T = nmonths
w= matrix(0,T,T)

for (i in 1:T) {
  if(i == 1) w[i,c(i+1, T/3*1+1, T/3*2+1 )] = c(1,1,1)
  if(i == T) w[i, c(i-1, 24, 12)] = c(1,1,1)
  else{
    w[i,c(i-1,i+1)] = 1
  }
}


for (i in 2:12) {
  w[i, c(T/3+i ,T/3*2+i)] = 1
}


#  w[1:12,]

for (i in 13:24) {
  w[i,c( abs(T/3-i) ,T/3+i)] = 1
}

#  w[13:24,]


for (i in 25:35) {
  w[i,c(abs(T/3*2-i), abs(T/3-i))] = 1
}
CarArg =  as.carAdjacency(w) # arguments for car distriubtion.


win.data <- list(Y=3,yN_index = yN_index, index = index, pos = pos, m = m, nb_levels = length(pos),  nmonths = T,  numVars = numVars, L = length(CarArg$adj), adj = CarArg$adj, weights= CarArg$weights, num = CarArg$num)
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
  
  
  
  # RW1 approximated by ICAR
  s[1:nmonths] ~ dcar_normal(adj[1:L], weights[1:L], num[1:nmonths], tau, zero_mean = 0)
  
  
  
  M ~ dpois(lambda)                    # Super-population size M
  
  # Yearly abundance
  for (y1 in 1:Y) {
    
    delta[y1] ~ dunif(0,1)
    yN[y1] ~ dbin(delta[y1], M)
    
  }
  
  for(t in 1:nmonths){                     # Loop over months
    
    logit(theta[t]) <- s[t] #~ dbeta(gammatilde[z[t]], betatilde[z[t]])
    N[t] ~ dbin(theta[t], yN[yN_index[t]])   
    
    linpred[1:index[t],t] <- (X[1:index[t],t, 1:numVars] %*% beta[1:numVars])[1:index[t],1]
    
    for (j in 1:index[t]) {
      y[j,t] ~ dbin(p[j,t], N[t])
      re[j,t] ~ dnorm(0, sigma.re)
      logit(p[j,t]) <- linpred[j,t] + re[j,t]
    }
    
  }
  
  
})

Nst <- apply(cparrots, 2, max, na.rm = TRUE)+100 # Inits for latent N
Mst <- max(Nst) + 500
Kst <-  rbinom(3, Mst, 0.5)
TEmodel = nimbleModel(nimTE, constants = win.data, data = list(y =Y, X=X), inits = list(delta = runif(3,0,1), yN =Kst, pi_g = rep(1, length(pos)), lambda2 =1, tau_g = rgamma(length(pos), 1,1), b = rep(0, length(pos)), sigma.re=1, re = matrix(rnorm(J*nmonths, 0,1), J, nmonths), M = Mst,  N= Nst,  lambda =  rgamma(1,10,0.1), beta = rep(0, numVars), s = rnorm(T, 0,1), sigma = runif(1) )  )
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
CorcrossValOutput <- runCrossValidate(MCMCconfiguration = TEconf,
                                      k = nmonths,
                                      foldFunction = TEFoldFunction,
                                      lossFunction = RMSElossFunction,
                                      MCMCcontrol = list(niter = 135000, nburnin = 35000, thin = 20),nBootReps = NA,
                                      silent = TRUE)



save(CorcrossValOutput, file ="BparrotsCor_KCV.Rdata")
#==================================================================================================================================================
rm(list = ls())
# Ar1 model

# Cleaned data 
load("cparrots.Rdata")
cparrots
J = dim(cparrots)[1]
nmonths = dim(cparrots)[2]

index = numeric(nmonths)

for (t in 1:nmonths) {
  index[t] = sum(1- is.na(cparrots[,t])) # create indicator where to end loop
  
}

index # number of sampling occasions in each month

Y =cparrots # need to remvove
Y

# Covariates 
load("Bparrots_Xmat.Rdata")
head(X)
dim(X)
pos <- c(1,2,3,4,5,6,7,8,8,8) # The ‘pos’ vector gives the position of the factor corresponding to each component of beta
m <- c(1,1,1,1,1,1,1,3,3,3) # The ‘m’ vector captures the same size group for all levels belonging to the same categorical covariate. 

numVars = dim(X)[3]

yN_index = rep(1:3, each = 12)

T = nmonths

win.data <- list(Y=3,yN_index = yN_index, index = index, pos = pos, m = m, nb_levels = length(pos),  nmonths = T,  numVars = numVars)
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
  
  
  # AR1
  s[1] ~ dnorm(0, sd =(tau*(1-rho^2))^(-0.5) )
  
  for (t in 2:nmonths) {
    s[t] ~ dnorm(rho*s[t-1], sd = sigma.e)
  }
  
  
  M ~ dpois(lambda)                    # Super-population size M
  
  # Yearly abundance
  for (y1 in 1:Y) {
    
    delta[y1] ~ dunif(0,1)
    yN[y1] ~ dbin(delta[y1], M)
    
  }
  
  
  for(t in 1:nmonths){                     # Loop over months
    
    logit(theta[t]) <- alpha + s[t] #~ dbeta(gammatilde[z[t]], betatilde[z[t]])
    N[t] ~ dbin(theta[t], yN[yN_index[t]])   
    
    linpred[1:index[t],t] <- (X[1:index[t],t, 1:numVars] %*% beta[1:numVars])[1:index[t],1]
    
    for (j in 1:index[t]) {
      y[j,t] ~ dbin(p[j,t], N[t])
      re[j,t] ~ dnorm(0, sigma.re)
      logit(p[j,t]) <- linpred[j,t] + re[j,t]
    }
    
    
  }
  
 
  
})

Nst <- apply(cparrots, 2, max, na.rm = TRUE)+100 # Inits for latent N
Mst <- max(Nst) + 500
Kst <-  rbinom(3, Mst, 0.5)
TEmodel = nimbleModel(nimTE, constants = win.data, data = list(y =Y, X=X), inits = list(delta = runif(3,0,1), yN =Kst,pi_g = rep(1, length(pos)), lambda2 =1, tau_g = rgamma(length(pos), 1,1), b = rep(0, length(pos)),sigma.re=1, sigma.e = 1, rho = runif(1, -1,1), alpha = 1, re = matrix(rnorm(J*nmonths, 0,1), J, nmonths), M = Mst,  N= Nst, lambda =  rgamma(1,10,0.1), beta = rep(0, numVars), s = rnorm(T, 0,1), sigma = runif(1) )  )
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
                                      k = nmonths,
                                      foldFunction = TEFoldFunction,
                                      lossFunction = RMSElossFunction,
                                      MCMCcontrol = list(niter = 135000, nburnin = 35000, thin = 20), nBootReps = NA,
                                      silent = TRUE)



save(AR1crossValOutput, file ="BparrotsAR1_KCV.Rdata")


