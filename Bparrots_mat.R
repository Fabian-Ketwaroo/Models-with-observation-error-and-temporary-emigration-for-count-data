{
  library(nimble)
  library(coda)
  #library(mcmcplots)
  #library(ggmcmc)
  #library(mcclust)
  #library(mcclust.ext)
  #library(readxl)
  #library(pgdraw)
  #library(data.table)
  #library(MCMCvis)
  
  
}


# Becca Parrots 

# Cleaned data 
load("~/Objective priors for N-mixture models/Parrots/cparrots.Rdata")
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/Parrots/Model selection/Bparrots_Xmat.Rdata")
#load("~/Misc/cparrots.Rdata")
#load("~/Misc/Bparrots_Xmat.Rdata")
#X = array(c( Intercept, Temp, Humidity, Visability, wind_speed,Rain, Storm, time, weatherCloud, weatherRain, weatherSunshine), dim =c(J, T,11))
X[,,6]
X[,,8]

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
#load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/Parrots/Model selection/Bparrots_Xmat.Rdata")
#head(X)
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

TEconf$monitors
TEconf$resetMonitors()
TEconf$addMonitors(c("lambda", "theta", "beta", "pi_g", "sigma.re", "M", "N", "z", "yN",  "gammatilde", "betatilde" ))
#TEconf$addMonitors2("p")

TEmcmc = buildMCMC(TEconf)
TEmod = compileNimble(TEmcmc, project = TEmodel, resetFunctions = TRUE)
#TEmcmc.out <- runMCMC(TEmod, niter = 55000, nburnin = 10000, thin = 20, samplesAsCodaMCMC = TRUE, nchains = 2, summary = TRUE) 
#TEmcmc.out <- runMCMC(TEmod, niter = 350, samplesAsCodaMCMC = TRUE, nchains = 1, summary = TRUE) 
#TEmcmc.out <- runMCMC(TEmod, niter = 155000, nburnin = 25000, thin = 20, samplesAsCodaMCMC = TRUE, nchains = 1, summary = TRUE) 
TEmcmc.out <- runMCMC(TEmod, niter = 135000, nburnin = 35000, thin = 20, samplesAsCodaMCMC = TRUE, nchains = 2, summary = TRUE) 
TEmcmc.out$summary$all.chains
#head(X)
#TE_WAIC  = TEmcmc.out$WAIC

samples <- TEmcmc.out$samples
save(samples, file = "Bparrots_DP_BGLSS_mat.Rdata")


#==================================================================================================================================================
#==================================================================================================================================================

# Rw1
library(nimble)
#source("parrots_BVS_samplers.R")


# Becca Parrots 
# Cleaned data 
load("~/Objective priors for N-mixture models/Parrots/cparrots.Rdata")
load("~/Objective priors for N-mixture models/Parrots/parrots_time.Rdata")
load("~/Misc/cparrots.Rdata")
load("~/Misc/Bparrots_Xmat.Rdata")


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
#head(X)
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
  for (y in 1:Y) {
    
    delta[y] ~ dunif(0,1)
    yN[y] ~ dbin(delta[y], M)
    
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
TEconf$monitors
TEconf$resetMonitors()
TEconf$addMonitors(c("lambda",  "theta", "beta", "pi_g", "M", "yN", "N", "sigma" ))
#TEconf$addMonitors(c("lambda", "sdBeta", "sigma.re", "theta", "beta", "w", "M", "N", "sigma" ))

TEmcmc = buildMCMC(TEconf)
TEmod = compileNimble(TEmcmc, project = TEmodel, resetFunctions = TRUE)
#TEmcmc.out <- runMCMC(TEmod, niter = 5000, nburnin = 1000, samplesAsCodaMCMC = TRUE, nchains = 1, summary = TRUE) 
#TEmcmc.out <- runMCMC(TEmod, niter = 150000, nburnin = 50000, thin = 20, samplesAsCodaMCMC = TRUE, nchains = 2, WAIC = TRUE, summary = TRUE) 
#TEmcmc.out <- runMCMC(TEmod, niter = 55000, nburnin = 10000, thin = 20, samplesAsCodaMCMC = TRUE, nchains = 2, WAIC = TRUE, summary = TRUE) 
TEmcmc.out <- runMCMC(TEmod, niter = 135000, nburnin = 35000, thin = 20, thin2 = 20, samplesAsCodaMCMC = TRUE, nchains = 2, WAIC = TRUE, summary = TRUE) 
samples <- TEmcmc.out$samples

save(samples, file = "Bparrots_RW1_BGLSS_mat.Rdata")


#==================================================================================================================================================

# Rw2
library(nimble)


# Becca Parrots 
# Cleaned data 
#load("~/Objective priors for N-mixture models/Parrots/cparrots.Rdata")
#load("~/Objective priors for N-mixture models/Parrots/parrots_time.Rdata")
load("~/Misc/cparrots.Rdata")
load("~/Misc/Bparrots_Xmat.Rdata")


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
#load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/Parrots/Model selection/Bparrots_Xmat.Rdata")
#head(X)
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
  for (y in 1:Y) {
    
    delta[y] ~ dunif(0,1)
    yN[y] ~ dbin(delta[y], M)
    
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
TEconf$monitors
TEconf$resetMonitors()
TEconf$addMonitors(c("lambda",  "theta", "beta", "pi_g", "M", "yN", "N", "sigma" ))
#TEconf$addMonitors(c("lambda", "sdBeta", "sigma.re", "theta", "beta", "w", "M", "N", "sigma" ))

TEmcmc = buildMCMC(TEconf)
TEmod = compileNimble(TEmcmc, project = TEmodel, resetFunctions = TRUE)
#TEmcmc.out <- runMCMC(TEmod, niter = 5000, nburnin = 1000, samplesAsCodaMCMC = TRUE, nchains = 1, summary = TRUE) 
#TEmcmc.out <- runMCMC(TEmod, niter = 150000, nburnin = 50000, thin = 20, samplesAsCodaMCMC = TRUE, nchains = 2, WAIC = TRUE, summary = TRUE) 
#TEmcmc.out <- runMCMC(TEmod, niter = 55000, nburnin = 10000, thin = 20, samplesAsCodaMCMC = TRUE, nchains = 2, WAIC = TRUE, summary = TRUE) 
TEmcmc.out <- runMCMC(TEmod, niter = 135000, nburnin = 35000, thin = 20, thin2 = 20, samplesAsCodaMCMC = TRUE, nchains = 2, WAIC = TRUE, summary = TRUE) 
samples <- TEmcmc.out$samples

save(samples, file = "Bparrots_RW2_BGLSS_mat.Rdata")


#=================================================================================================================================================
# Cor
library(nimble)
#source("parrots_BVS_samplers.R")


# Becca Parrots 
# Cleaned data 
load("~/Objective priors for N-mixture models/Parrots/cparrots.Rdata")
load("~/Objective priors for N-mixture models/Parrots/parrots_time.Rdata")
#load("~/Misc/cparrots.Rdata")
#load("~/Misc/Bparrots_Xmat.Rdata")


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
#load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/Parrots/Model selection/Bparrots_Xmat.Rdata")
#head(X)
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
  for (y in 1:Y) {
    
    delta[y] ~ dunif(0,1)
    yN[y] ~ dbin(delta[y], M)
    
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
TEconf$monitors
TEconf$resetMonitors()
TEconf$addMonitors(c("lambda",  "theta", "beta", "pi_g", "M", "yN", "N", "sigma" ))
#TEconf$addMonitors(c("lambda", "sdBeta", "sigma.re", "theta", "beta", "w", "M", "N", "sigma" ))

TEmcmc = buildMCMC(TEconf)
TEmod = compileNimble(TEmcmc, project = TEmodel, resetFunctions = TRUE)
#TEmcmc.out <- runMCMC(TEmod, niter = 5000, nburnin = 1000, samplesAsCodaMCMC = TRUE, nchains = 1, summary = TRUE) 
#TEmcmc.out <- runMCMC(TEmod, niter = 150000, nburnin = 50000, thin = 20, samplesAsCodaMCMC = TRUE, nchains = 2, WAIC = TRUE, summary = TRUE) 
#TEmcmc.out <- runMCMC(TEmod, niter = 55000, nburnin = 10000, thin = 20, samplesAsCodaMCMC = TRUE, nchains = 2, WAIC = TRUE, summary = TRUE) 
TEmcmc.out <- runMCMC(TEmod, niter = 135000, nburnin = 35000, thin = 20, thin2 = 20, samplesAsCodaMCMC = TRUE, nchains = 2, WAIC = TRUE, summary = TRUE) 
samples <- TEmcmc.out$samples

save(samples, file = "Bparrots_Cor_BGLSS_mat.Rdata")

#==================================================================================================================================================
# Ar1 model

library(nimble)
#source("parrots_BVS_samplers.R")


# Becca Parrots 
# Cleaned data 
#load("~/Objective priors for N-mixture models/Parrots/cparrots.Rdata")
#load("~/Objective priors for N-mixture models/Parrots/parrots_time.Rdata")
load("~/Misc/cparrots.Rdata")
load("~/Misc/Bparrots_Xmat.Rdata")


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
#load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/Parrots/Model selection/Bparrots_Xmat.Rdata")
#head(X)
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
  for (y in 1:Y) {
    
    delta[y] ~ dunif(0,1)
    yN[y] ~ dbin(delta[y], M)
    
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
TEconf$monitors
TEconf$resetMonitors()
TEconf$addMonitors(c("lambda",  "theta", "beta", "pi_g", "M", "yN", "N", "sigma", "rho" ))
#TEconf$addMonitors(c("lambda", "sdBeta", "sigma.re", "theta", "beta", "w", "M", "N", "sigma" ))

TEmcmc = buildMCMC(TEconf)
TEmod = compileNimble(TEmcmc, project = TEmodel, resetFunctions = TRUE)
#TEmcmc.out <- runMCMC(TEmod, niter = 5000, nburnin = 1000, samplesAsCodaMCMC = TRUE, nchains = 1, summary = TRUE) 
#TEmcmc.out <- runMCMC(TEmod, niter = 150000, nburnin = 50000, thin = 20, samplesAsCodaMCMC = TRUE, nchains = 2, WAIC = TRUE, summary = TRUE) 
#TEmcmc.out <- runMCMC(TEmod, niter = 55000, nburnin = 10000, thin = 20, samplesAsCodaMCMC = TRUE, nchains = 2, WAIC = TRUE, summary = TRUE) 
TEmcmc.out <- runMCMC(TEmod, niter = 135000, nburnin = 35000, thin = 20, thin2 = 20, samplesAsCodaMCMC = TRUE, nchains = 2, WAIC = TRUE, summary = TRUE) 
samples <- TEmcmc.out$samples

save(samples, file = "Bparrots_AR1_BGLSS_mat.Rdata")


#==================================================================================================================================================
# Results 

# Estimates of N, avaiabiliy and DP clusters
#==================================================================================================================================================

{
  library(nimble)
  library(coda)
  library(mcmcplots)
  library(ggmcmc)
  library(mcclust)
  library(mcclust.ext)
  library(readxl)
  library(pgdraw)
  library(data.table)
  library(MCMCvis)
}

# DP results 
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/Parrots/Model selection/Bparrots_DP_BGLSS_mat.Rdata")
#load("~/Objective priors for N-mixture models/Parrots/Bparrots_BGLSS_psamples.Rdata")
head(samples)

#MCMCsummary(samples[, c(1:45, 110:120) ])

samples = mcmc(do.call(rbind, samples))
head(samples)
dim(samples)
zs <-samples[,172:207]
dim(zs)
head(zs)
thetas <- MCMCsummary(as.matrix(samples[ ,133:168]),Rhat = F, n.eff = F)
tail(thetas)
head(thetas)
dim(thetas)
d = data.frame(thetas[,1])

psm=comp.psm(zs)
z.VI=minVI(psm,zs,method=("all"),include.greedy=T)
summary(z.VI)
#plot(z.VI,data=d)

z.ci <- as.data.frame(z.VI$cl) #matrix(z.VI$cl,5,12)
colnames(z.ci) <- rep(c("Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct"),3)
z.ci
length(z.ci[1,1:12])

z.summary <- data.frame(matrix(0, 3,13))
colnames(z.summary) <- c("Year","Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct")
z.summary[,1] <- c(1,2,3)
z.summary[1,-1] <- z.ci[5,1:12]
z.summary[2,-1] <- z.ci[5,13:24]
z.summary[3,-1] <- z.ci[5,25:36]
print(z.summary, row.names = F)


# yN, N and theta summary 
#samples = mcmc(do.call(rbind, samples))
head(samples)
N <- samples[,2:37]
head(N)
colnames(N) = c("Nov1", "Dec1", "Jan1", "Feb1", "Mar1", "Apr1", "May1", "Jun1", "Jul1", "Aug1", "Sep1", "Oct1", "Nov2", "Dec2", "Jan2", "Feb2", "Mar2", "Apr2", "May2", "Jun2", "Jul2", "Aug2", "Sep2", "Oct2", "Nov3", "Dec3", "Jan3", "Feb3", "Mar3", "Apr3", "May3", "Jun3", "Jul3", "Aug3", "Sep3", "Oct3"  )
Ns = ggs(N, keep_original_order = T)
#head(Ns)
#Ns$Parameter= sort(Ns$Parameter, decreasing = TRUE)
dim(Ns)
#ggs_caterpillar(Ns, horizontal = F, sort = F, thin_ci = c(0.5,0.5)) + xlab("Population size") + ylab("Month") + ggtitle(" ") + theme_classic() + theme( axis.title.x = element_text(size=14), axis.title.y = element_text(size=14) ,axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust=1))
#ggs_caterpillar(Ns, horizontal = F, sort = F, thin_ci = c(0.5,0.5)) + xlab("Population size") + ylab("Month") + ggtitle(" ")  +theme_classic() + theme( axis.title.x = element_text(size=20), axis.title.y = element_text(size=20), axis.text.y = element_text(size = 12) ,axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust=1))  
ggs_caterpillar(Ns, horizontal = F, sort = F, thin_ci = c(0.5,0.5)) + xlab("Population size") + ylab("Month") + ggtitle(" ")  +theme_classic() + theme( axis.title.x = element_text(size=20), axis.title.y = element_text(size=20), axis.text.y = element_text(size = 12) ,axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust=1))  + scale_x_continuous(limits = c(0,275), breaks = seq(0, 275, by=55), expand =c(0,0) )


head(samples)

ynn <- samples[, 169:171]
head(ynn)
colnames(ynn) = c("Year 1", "Year 2", "Year 3")
ggs_caterpillar(ggs(ynn), horizontal = F, sort = F) + xlab("Abundance") + ylab("Year") + ggtitle("Posterior Summary")


the = samples[,133:168]
colnames(the) = c("Nov1", "Dec1", "Jan1", "Feb1", "Mar1", "Apr1", "May1", "Jun1", "Jul1", "Aug1", "Sep1", "Oct1", "Nov2", "Dec2", "Jan2", "Feb2", "Mar2", "Apr2", "May2", "Jun2", "Jul2", "Aug2", "Sep2", "Oct2", "Nov3", "Dec3", "Jan3", "Feb3", "Mar3", "Apr3", "May3", "Jun3", "Jul3", "Aug3", "Sep3", "Oct3"  )
Ns = ggs(the, keep_original_order = T)
head(Ns)
ggs_caterpillar(Ns, horizontal = F, sort = F) + xlab("Availability") + ylab("Month") + ggtitle("Posterior Summary")

# BVS
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/Parrots/Model selection/Bparrots_DP_BGLSS_mat.Rdata")
head(samples)
(betas = MCMCsummary(samples)[38:48,])

X = data.frame(intercept, Temp = Temp, Humidity = Humidity, Visability = Visability, AWS=AWS, Rain = Rain, Storm = Storm, Time = Time, weather = weather)# need to add weather 
(pi_g = MCMCsummary(samples)[122:129,])

exp(betas[c(6,7,8),][,1])

#==================================================================================================================================================
# Rw1
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/Parrots/Model selection/Bparrots_RW1_BGLSS_mat.Rdata")
head(samples)

#MCMCsummary(samples[, c(1:45, 110:120) ])

samples = mcmc(do.call(rbind, samples))

head(samples)
N <- samples[,2:37]
head(N)
colnames(N) = c("Nov1", "Dec1", "Jan1", "Feb1", "Mar1", "Apr1", "May1", "Jun1", "Jul1", "Aug1", "Sep1", "Oct1", "Nov2", "Dec2", "Jan2", "Feb2", "Mar2", "Apr2", "May2", "Jun2", "Jul2", "Aug2", "Sep2", "Oct2", "Nov3", "Dec3", "Jan3", "Feb3", "Mar3", "Apr3", "May3", "Jun3", "Jul3", "Aug3", "Sep3", "Oct3"  )
Ns = ggs(N, keep_original_order = T)
#head(Ns)
#Ns$Parameter= sort(Ns$Parameter, decreasing = TRUE)
dim(Ns)
#ggs_caterpillar(Ns, horizontal = F, sort = F, thin_ci = c(0.5,0.5)) + xlab("Abundance") + ylab("Month") + ggtitle(" ") + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#ggs_caterpillar(Ns, horizontal = F, sort = F, thin_ci = c(0.5,0.5)) + xlab("Population size") + ylab("Month") + ggtitle(" ")  +theme_classic() + theme( axis.title.x = element_text(size=20), axis.title.y = element_text(size=20), axis.text.y = element_text(size = 12) ,axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust=1))  
ggs_caterpillar(Ns, horizontal = F, sort = F, thin_ci = c(0.5,0.5)) + xlab("Population size") + ylab("Month") + ggtitle(" ")  +theme_classic() + theme( axis.title.x = element_text(size=20), axis.title.y = element_text(size=20), axis.text.y = element_text(size = 12) ,axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust=1))  + scale_x_continuous(limits = c(0,250), breaks = seq(0, 250, by=50), expand =c(0,0) )


the = samples[,61:96]
colnames(the) = c("Nov1", "Dec1", "Jan1", "Feb1", "Mar1", "Apr1", "May1", "Jun1", "Jul1", "Aug1", "Sep1", "Oct1", "Nov2", "Dec2", "Jan2", "Feb2", "Mar2", "Apr2", "May2", "Jun2", "Jul2", "Aug2", "Sep2", "Oct2", "Nov3", "Dec3", "Jan3", "Feb3", "Mar3", "Apr3", "May3", "Jun3", "Jul3", "Aug3", "Sep3", "Oct3"  )
Ns = ggs(the, keep_original_order = T)
#head(Ns)
ggs_caterpillar(Ns, horizontal = F, sort = F) + xlab("Availability") + ylab("Month") + ggtitle("Posterior Summary")


# BVs
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/Parrots/Model selection/Bparrots_RW1_BGLSS_mat.Rdata")

(betas = MCMCsummary(samples)[38:48,])

(betas = MCMCsummary(samples, digits = 3)[38:48,c(1,2,3,5)])

#X = data.frame(itercept, Temp = Temp, Humidity = Humidity, Visability = Visability, AWS=AWS, Rain = Rain, Storm = Storm, Time = Time, weather = weather)# need to add weather 

(pi_g = MCMCsummary(samples)[50:57,])

betas[c(6,7,8),]
exp(betas[c(6,7,8),][,1])
#1-exp(betas[c(5),][,1])
(exp(betas[c(6,7,8),][,1])-1)*100

#==================================================================================================================================================
# Rw2
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/Parrots/Model selection/Bparrots_RW2_BGLSS_mat.Rdata")
head(samples)

#MCMCsummary(samples[, c(1:45, 110:120) ])

samples = mcmc(do.call(rbind, samples))

head(samples)
N <- samples[,2:37]
head(N)
colnames(N) = c("Nov1", "Dec1", "Jan1", "Feb1", "Mar1", "Apr1", "May1", "Jun1", "Jul1", "Aug1", "Sep1", "Oct1", "Nov2", "Dec2", "Jan2", "Feb2", "Mar2", "Apr2", "May2", "Jun2", "Jul2", "Aug2", "Sep2", "Oct2", "Nov3", "Dec3", "Jan3", "Feb3", "Mar3", "Apr3", "May3", "Jun3", "Jul3", "Aug3", "Sep3", "Oct3"  )
Ns = ggs(N, keep_original_order = T)
#head(Ns)
#Ns$Parameter= sort(Ns$Parameter, decreasing = TRUE)
dim(Ns)
#ggs_caterpillar(Ns, horizontal = F, sort = F, thin_ci = c(0.5,0.5)) + xlab("Abundance") + ylab("Month") + ggtitle(" ")+ theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#ggs_caterpillar(Ns, horizontal = F, sort = F, thin_ci = c(0.5,0.5)) + xlab("Population size") + ylab("Month") + ggtitle(" ")  +theme_classic() + theme( axis.title.x = element_text(size=20), axis.title.y = element_text(size=20), axis.text.y = element_text(size = 12) ,axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust=1))  
ggs_caterpillar(Ns, horizontal = F, sort = F, thin_ci = c(0.5,0.5)) + xlab("Population size") + ylab("Month") + ggtitle(" ")  +theme_classic() + theme( axis.title.x = element_text(size=20), axis.title.y = element_text(size=20), axis.text.y = element_text(size = 12) ,axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust=1))  + scale_x_continuous(limits = c(0,250), breaks = seq(0, 250, by=50), expand =c(0,0) )



the = samples[,61:96]
colnames(the) = c("Nov1", "Dec1", "Jan1", "Feb1", "Mar1", "Apr1", "May1", "Jun1", "Jul1", "Aug1", "Sep1", "Oct1", "Nov2", "Dec2", "Jan2", "Feb2", "Mar2", "Apr2", "May2", "Jun2", "Jul2", "Aug2", "Sep2", "Oct2", "Nov3", "Dec3", "Jan3", "Feb3", "Mar3", "Apr3", "May3", "Jun3", "Jul3", "Aug3", "Sep3", "Oct3"  )
Ns = ggs(the, keep_original_order = T)
#head(Ns)
ggs_caterpillar(Ns, horizontal = F, sort = F) + xlab("Availability") + ylab("Month") + ggtitle("Posterior Summary")


# BVs
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/Parrots/Model selection/Bparrots_RW2_BGLSS_mat.Rdata")

(betas = MCMCsummary(samples)[38:48,])

X = data.frame(Temp = Temp, Humidity = Humidity, Visability = Visability, AWS=AWS, Rain = Rain, Storm = Storm, Time = Time, weather = weather)# need to add weather 

(pi_g = MCMCsummary(samples)[50:57,])

betas[c(6,7,8),]
exp(betas[c(6,7,8),1])
#==================================================================================================================================================
# cor
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/Parrots/Model selection/Bparrots_Cor_BGLSS_mat.Rdata")
head(samples)

#MCMCsummary(samples[, c(1:45, 110:120) ])

samples = mcmc(do.call(rbind, samples))

head(samples)
N <- samples[,2:37]
head(N)
colnames(N) = c("Nov1", "Dec1", "Jan1", "Feb1", "Mar1", "Apr1", "May1", "Jun1", "Jul1", "Aug1", "Sep1", "Oct1", "Nov2", "Dec2", "Jan2", "Feb2", "Mar2", "Apr2", "May2", "Jun2", "Jul2", "Aug2", "Sep2", "Oct2", "Nov3", "Dec3", "Jan3", "Feb3", "Mar3", "Apr3", "May3", "Jun3", "Jul3", "Aug3", "Sep3", "Oct3"  )
Ns = ggs(N, keep_original_order = T)
#head(Ns)
#Ns$Parameter= sort(Ns$Parameter, decreasing = TRUE)
dim(Ns)
#ggs_caterpillar(Ns, horizontal = F, sort = F, thin_ci = c(0.5,0.5)) + xlab("Abundance") + ylab("Month") + ggtitle(" ")+ theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggs_caterpillar(Ns, horizontal = F, sort = F, thin_ci = c(0.5,0.5)) + xlab("Population size") + ylab("Month") + ggtitle(" ")  +theme_classic() + theme( axis.title.x = element_text(size=20), axis.title.y = element_text(size=20), axis.text.y = element_text(size = 12) ,axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust=1))  

ggs_caterpillar(Ns, horizontal = F, sort = F, thin_ci = c(0.5,0.5)) + xlab("Population size") + ylab("Month") + ggtitle(" ")  +theme_classic() + theme( axis.title.x = element_text(size=20), axis.title.y = element_text(size=20), axis.text.y = element_text(size = 12) ,axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust=1))  + scale_x_continuous(limits = c(0,250), breaks = seq(0, 250, by=50), expand =c(0,0) )
#+ scale_x_continuous(limits = c(0,12000), expand = c(0,0))

the = samples[,61:96]
colnames(the) = c("Nov1", "Dec1", "Jan1", "Feb1", "Mar1", "Apr1", "May1", "Jun1", "Jul1", "Aug1", "Sep1", "Oct1", "Nov2", "Dec2", "Jan2", "Feb2", "Mar2", "Apr2", "May2", "Jun2", "Jul2", "Aug2", "Sep2", "Oct2", "Nov3", "Dec3", "Jan3", "Feb3", "Mar3", "Apr3", "May3", "Jun3", "Jul3", "Aug3", "Sep3", "Oct3"  )
Ns = ggs(the, keep_original_order = T)
#head(Ns)
ggs_caterpillar(Ns, horizontal = F, sort = F) + xlab("Availability") + ylab("Month") + ggtitle("Posterior Summary")


# BVs
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/Parrots/Model selection/Bparrots_Cor_BGLSS_mat.Rdata")
(betas = MCMCsummary(samples)[38:48,])

(betas = MCMCsummary(samples, digits = 3)[38:48,c(1,2,3,5)])

#X = data.frame(Temp = Temp, Humidity = Humidity, Visability = Visability, AWS=AWS, Rain = Rain, Storm = Storm, Time = Time, weather = weather)# need to add weather 

(pi_g = MCMCsummary(samples)[50:57,])

betas[c(6,7,8),]
exp(betas[c(6,7,8),][,1])

#==================================================================================================================================================
# AR1
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/Parrots/Model selection/Bparrots_AR1_BGLSS_mat.Rdata")
head(samples)

#MCMCsummary(samples[, c(1:45, 110:120) ])

samples = mcmc(do.call(rbind, samples))

head(samples)
N <- samples[,2:37]
head(N)
colnames(N) = c("Nov1", "Dec1", "Jan1", "Feb1", "Mar1", "Apr1", "May1", "Jun1", "Jul1", "Aug1", "Sep1", "Oct1", "Nov2", "Dec2", "Jan2", "Feb2", "Mar2", "Apr2", "May2", "Jun2", "Jul2", "Aug2", "Sep2", "Oct2", "Nov3", "Dec3", "Jan3", "Feb3", "Mar3", "Apr3", "May3", "Jun3", "Jul3", "Aug3", "Sep3", "Oct3"  )
Ns = ggs(N, keep_original_order = T)
#head(Ns)
#Ns$Parameter= sort(Ns$Parameter, decreasing = TRUE)
dim(Ns)
#ggs_caterpillar(Ns, horizontal = F, sort = F) + xlab("Abundance") + ylab("Month") + ggtitle("Posterior Summary")
#ggs_caterpillar(Ns, horizontal = F, sort = F, thin_ci = c(0.5,0.5)) + xlab("Abundance") + ylab("Month") + ggtitle(" ") + theme_classic()  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#ggs_caterpillar(Ns, horizontal = F, sort = F, thin_ci = c(0.5,0.5)) + xlab("Population size") + ylab("Month") + ggtitle(" ")  +theme_classic() + theme( axis.title.x = element_text(size=20), axis.title.y = element_text(size=20), axis.text.y = element_text(size = 12) ,axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust=1))  
ggs_caterpillar(Ns, horizontal = F, sort = F, thin_ci = c(0.5,0.5)) + xlab("Population size") + ylab("Month") + ggtitle(" ")  +theme_classic() + theme( axis.title.x = element_text(size=20), axis.title.y = element_text(size=20), axis.text.y = element_text(size = 12) ,axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust=1))  + scale_x_continuous(limits = c(0,250), breaks = seq(0, 250, by=50), expand =c(0,0) )


the = samples[,62:97]
colnames(the) = c("Nov1", "Dec1", "Jan1", "Feb1", "Mar1", "Apr1", "May1", "Jun1", "Jul1", "Aug1", "Sep1", "Oct1", "Nov2", "Dec2", "Jan2", "Feb2", "Mar2", "Apr2", "May2", "Jun2", "Jul2", "Aug2", "Sep2", "Oct2", "Nov3", "Dec3", "Jan3", "Feb3", "Mar3", "Apr3", "May3", "Jun3", "Jul3", "Aug3", "Sep3", "Oct3"  )
Ns = ggs(the, keep_original_order = T)
#head(Ns)
ggs_caterpillar(Ns, horizontal = F, sort = F) + xlab("Availability") + ylab("Month") + ggtitle("Posterior Summary")


# BVs
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/Parrots/Model selection/Bparrots_AR1_BGLSS_mat.Rdata")
(betas = MCMCsummary(samples)[38:48,])

#X = data.frame(Temp = Temp, Humidity = Humidity, Visability = Visability, AWS=AWS, Rain = Rain, Storm = Storm, Time = Time, weather = weather)# need to add weather 

(pi_g = MCMCsummary(samples)[50:57,])


betas[c(6,7,8),]
exp(betas[c(6,7,8),][,1])
#==================================================================================================================================================

#========================================================================================================================================
#PIP plot

rm(list=ls())
# Cor results 
library(MCMCvis)
#library(ggplot2)
load("C:/Users/Fabian Ketwaroo/OneDrive - University of Kent/Parrots/Model selection/Bparrots_Cor_BGLSS_mat.Rdata")

head(samples)
#samples = mcmc(do.call(rbind, samples))
MCMCsummary(samples[, c(50:57) ])
dim(samples)
PIP =  MCMCsummary(samples[, c(50:57) ])[,1]
#PIP = as.matrix(PIP)
Variables = c("Temperature", "Humidity", "Visability", "Wind Speed", "Rain", "Storm", "Time", "Weather")

plot(PIP, pch = 15, xlab = "Variable", xaxt="n", ylab = "Inclusion Pobability")
abline(h=0.5)
axis(1, at=1:8, labels=Variables, las=2, cex.axis=0.75)




load("~/Objective priors for N-mixture models/Parrots/Oparrots_TE_RW1_BGLSS.Rdata")
MCMCsummary(samples[, c(23:26) ])
dim(samples)
PIP = MCMCsummary(samples[, c(23:26) ])[,1]
PIP
PIP = c(0.567, 0.412,0.345,0.44)
x = factor(1:4)
plot(PIP, pch = 15, xlab = "Variable")
abline(h=0.5)
axis(1:4, lwd =1)


x <- seq(10,200,10)
y <- runif(x)

plot(x,y)



