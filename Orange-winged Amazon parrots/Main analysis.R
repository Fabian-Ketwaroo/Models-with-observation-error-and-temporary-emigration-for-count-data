
# Rscript to analyse Orange parrots data sets to class of TE models developed
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

#=============================================================================================================================================
# Clean data
#=============================================================================================================================================

# Cleaned data 
load("Oparrots.Rdata")

Oparrots
J = dim(Oparrots)[1]
nmonths = dim(Oparrots)[2]

{
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
  
  weeks_Opar = cbind(week1, week2,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12,w13)
  nweeks =dim(weeks_Opar)[2]
  
  
}


index = numeric(nweeks)

for (t in 1:nweeks) {
  index[t] = sum(1- is.na(weeks_Opar[,t])) # create indicator where to end loop
  
}

nweks_index = index[index!=0] # number of sampling occasions in each month
nweks_index
nweeks =length(nweks_index)

P =c(weeks_Opar) # need to remvove
P
y = P[!is.na(P)]
y

S = length(y)
Nindex = rep(1:nweeks, nweks_index)
length(Nindex)

# Covariates 
load("orange_parrots.Rdata")
head(orange_parrots)

# Center and standardisze continouous variables 
partly_cloud = as.factor(orange_parrots$pcloud)
cloudy = as.factor(orange_parrots$cloudy)
Rain = as.factor(orange_parrots$crain)
low_wind = as.factor(orange_parrots$low_wind)
strong_wind = as.factor(orange_parrots$strong_wind)
Time = as.factor(orange_parrots$ctime) # 1 represents PM
Wind = as.factor(orange_parrots$Wind)
Cloud = as.factor(orange_parrots$Cloudiness)

ncov = 4 # number of explanatory variables in the model excluding the intercept
X = data.frame(Cloud = Cloud, Wind = Wind, Rain = Rain,Time = Time)

head(X)
originalX <- X # backup the original  matrix of covariates
X <- model.matrix(~ ., X) # creates a design (or model) matrix
# Here the design matrix has a column of ones for the intercept
head(X)

pos <- c(1,1,2,2,3,4) # The ‘pos’ vector gives the position of the factor corresponding to each component of beta
m <- c(2,2,2,2,1,1) # The ‘m’ vector captures the same size group for all levels belonging to the same categorical covariate. 


numVars = dim(X)[2]#ncov+1 #number of covariates (including the intercept)

#===================================================================================================================================================
# DP
#===================================================================================================================================================


win.data <- list(S= S,X =X, Nindex = Nindex, pos= pos, m = m, nb_levels = length(pos), nweeks = nweeks,  numVars = numVars, K= nweeks)
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
  
  linpred[1:S] <- (X[1:S, 1:numVars] %*% beta[1:numVars])[1:S,1]   
  
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
    
  }
  
  for (s in 1:S) {
    y[s] ~ dbin(p[s], N[Nindex[s]])
    re[s] ~ dnorm(0, sigma.re)
    logit(p[s]) <- linpred[s] + re[s]
  }
  
  
})

Nst <- colSums(weeks_Opar, na.rm = TRUE)  #apply(Oparrots, 2, max, na.rm = TRUE)+100 # Inits for latent N
Nst = Nst[Nst !=0]
Mst <- max(Nst) + 500
TEmodel = nimbleModel(nimTE, constants = win.data,data = list(y =y), inits = list(M = Mst,  N= Nst, lambda =  1000,lambda2 =1, tau_g = rgamma(length(pos), 1,1),b = rep(0, length(pos)), sigma.re = 1, pi_g = rep(1, length(pos)), re = rnorm(S, 0,1), w = rep(1, numVars),beta = rep(0, numVars),theta = rbeta(nweeks,0.5, 0.5), z = sample(1:win.data$K, size = nweeks, replace = TRUE), alpha = 1, gammatilde = rgamma(win.data$K, 0.5, 1) ,betatilde = rgamma(win.data$K, 0.5, 1) ))
TEmodel$calculate()

cTEmodel <- compileNimble(TEmodel)
TEconf = configureMCMC(TEmodel, enableWAIC = TRUE)

TEconf$monitors
TEconf$resetMonitors()
TEconf$addMonitors(c("lambda", "pi_g", "theta", "beta", "sigma.re", "M", "N", "z",  "gammatilde", "betatilde" ))

TEmcmc = buildMCMC(TEconf)
TEmod = compileNimble(TEmcmc, project = TEmodel, resetFunctions = TRUE)
TEmcmc.out <- runMCMC(TEmod, niter = 205000, nburnin = 105000, thin = 20, thin2= 20, samplesAsCodaMCMC = TRUE, nchains = 2, WAIC = TRUE, summary = TRUE) 
TEmcmc.out$summary$all.chains
Opar_WAIC = TEmcmc.out$WAIC
samples <- TEmcmc.out$samples
head(samples)

save(samples, file = "weeklyOparrots_TE_DP_BGLSS1.Rdata")

#===================================================================================================================================================
# RW1
#===================================================================================================================================================

T = nweeks
w= matrix(0,T,T)

for (i in 1:T) {
  if(i == 1) w[i,i+1] = 1
  if(i == T) w[i, i-1] =1
  else{
    w[i,c(i-1,i+1)] = 1
  }
}

w
CarArg =  as.carAdjacency(w) # arguments for car distriubtion.


pos <- c(1,1,2,2,3,4) # The ‘pos’ vector gives the position of the factor corresponding to each component of beta
m <- c(2,2,2,2,1,1) # The ‘m’ vector captures the same size group for all levels belonging to the same categorical covariate. 


win.data <- list(S= S, X= X, Nindex = Nindex,pos= pos, m = m, nb_levels = length(pos),  nmonths = nweeks,  numVars = numVars, L = length(CarArg$adj), adj = CarArg$adj, weights= CarArg$weights, num = CarArg$num)
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
  
  linpred[1:S] <- (X[1:S, 1:numVars] %*% beta[1:numVars])[1:S,1]   
  
  # RW1 approximated by ICAR
  s[1:nmonths] ~ dcar_normal(adj[1:L], weights[1:L], num[1:nmonths], tau, zero_mean = 0)
  
  
  
  M ~ dpois(lambda)                    # Super-population size M
  
  
  for(t in 1:nmonths){                     # Loop over months
    
    logit(theta[t]) <- s[t] #~ dbeta(gammatilde[z[t]], betatilde[z[t]])
    N[t] ~ dbin(theta[t], M)            # 'Current' population size N
    
  }
  
  for (s1 in 1:S) {
    y[s1] ~ dbin(p[s1], N[Nindex[s1]])
    re[s1] ~ dnorm(0, sigma.re)
    logit(p[s1]) <- linpred[s1] + re[s1]
  }
  
})


Nst <- colSums(weeks_Opar, na.rm = TRUE)  #apply(Oparrots, 2, max, na.rm = TRUE)+100 # Inits for latent N
Nst = Nst[Nst !=0]
Mst <- max(Nst) + 500
TEmodel = nimbleModel(nimTE, constants = win.data, data = list(y =y), inits = list(sigma.re=1,pi_g = rep(1, length(pos)), lambda2 =1, tau_g = rgamma(length(pos), 1,1),b = rep(0, length(pos)), re = rnorm(S,0,1), M = Mst,  N= Nst, sdBeta = runif(1, 0,4), lambda =  rgamma(1,10,0.1), beta = rep(0, numVars), s = rnorm(nweeks, 0,1), sigma = runif(1) )  )
TEmodel$calculate()

cTEmodel <- compileNimble(TEmodel)
TEconf = configureMCMC(TEmodel, enableWAIC = TRUE)

TEconf$monitors
TEconf$resetMonitors()
TEconf$addMonitors(c("lambda", "pi_g", "theta", "beta",  "M", "N", "sigma" ))


TEmcmc = buildMCMC(TEconf)
TEmod = compileNimble(TEmcmc, project = TEmodel, resetFunctions = TRUE)
TEmcmc.out <- runMCMC(TEmod, niter = 205000, nburnin = 105000, thin = 20, thin2= 20, samplesAsCodaMCMC = TRUE, nchains = 2, WAIC = TRUE, summary = TRUE) 
TEmcmc.out$summary$all.chains
Opar_RW1_WAIC = TEmcmc.out$WAIC
samples <- TEmcmc.out$samples
#MCMCsummary(samples)
head(samples)
MCMCsummary(samples)

save(samples, file = "weeklyOparrots_TE_RW1_BGLSS1.Rdata")

#===================================================================================================================================================
# RW2
#===================================================================================================================================================

T = nweeks
w= matrix(0,T,T)

w[1, c(1+1,1+2)] = c(2, -1)
w[T, c(T-1, T-2)] = c(2,-1)
w[2, c(2-1,2+1, 2+2)] = c(2, 4,-1)

for (i in 3:(T-2)) {
  w[i,c(i-2,i-1,i+1,i+2)] = c(-1,4,4,-1)
} 

w[T-1, c(T-1-2, T-1-1, T)] = c(-1,4,2)

w
CarArg =  as.carAdjacency(w) # arguments for car distriubtion.


pos <- c(1,1,2,2,3,4) # The ‘pos’ vector gives the position of the factor corresponding to each component of beta
m <- c(2,2,2,2,1,1) # The ‘m’ vector captures the same size group for all levels belonging to the same categorical covariate. 


win.data <- list(S= S, X=X,  Nindex = Nindex,pos= pos, m = m, nb_levels = length(pos),  nmonths = nweeks,  numVars = numVars, L = length(CarArg$adj), adj = CarArg$adj, weights= CarArg$weights, num = CarArg$num)
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
  
  linpred[1:S] <- (X[1:S, 1:numVars] %*% beta[1:numVars])[1:S,1]   
  
  
  # RW1 approximated by ICAR
  s[1:nmonths] ~ dcar_normal(adj[1:L], weights[1:L], num[1:nmonths], tau, zero_mean = 0)
  
  
  
  M ~ dpois(lambda)                    # Super-population size M
  
  
  for(t in 1:nmonths){                     # Loop over months
    
    logit(theta[t]) <- s[t] #~ dbeta(gammatilde[z[t]], betatilde[z[t]])
    N[t] ~ dbin(theta[t], M)            # 'Current' population size N
    
  }
  
  
  for (s1 in 1:S) {
    y[s1] ~ dbin(p[s1], N[Nindex[s1]])
    re[s1] ~ dnorm(0, sigma.re)
    logit(p[s1]) <- linpred[s1] + re[s1]
  }
  
})

Nst <- colSums(weeks_Opar, na.rm = TRUE)  #apply(Oparrots, 2, max, na.rm = TRUE)+100 # Inits for latent N
Nst = Nst[Nst !=0]
Mst <- max(Nst) + 500
TEmodel = nimbleModel(nimTE, constants = win.data, data = list(y =y), inits = list(sigma.re=1,pi_g = rep(1, length(pos)), lambda2 =1, tau_g = rgamma(length(pos), 1,1),b = rep(0, length(pos)), re = rnorm(S,0,1), M = Mst,  N= Nst,  lambda =  rgamma(1,10,0.1), beta = rep(0, numVars), s = rnorm(nweeks, 0,1), sigma = runif(1) )  )
TEmodel$calculate()

cTEmodel <- compileNimble(TEmodel)
TEconf = configureMCMC(TEmodel, enableWAIC = TRUE)

TEconf$monitors
TEconf$resetMonitors()
TEconf$addMonitors(c("lambda", "pi_g", "theta", "beta",  "M", "N", "sigma" ))

TEmcmc = buildMCMC(TEconf)
TEmod = compileNimble(TEmcmc, project = TEmodel, resetFunctions = TRUE)
TEmcmc.out <- runMCMC(TEmod, niter = 405000, nburnin = 155000, thin = 30, thin2= 20, samplesAsCodaMCMC = TRUE, nchains = 2, WAIC = TRUE, summary = TRUE) 
TEmcmc.out$summary$all.chains
samples <- TEmcmc.out$samples
save(samples, file = "weeklyOparrots_TE_RW2_RE_BGLSS1.Rdata")

#===================================================================================================================================================
# AR1
#===================================================================================================================================================


win.data <- list(S= S, Nindex = Nindex,pos= pos, m = m, nb_levels = length(pos),  nmonths = nweeks,  numVars = numVars)
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
  
  linpred[1:S] <- (X[1:S, 1:numVars] %*% beta[1:numVars])[1:S,1]   
  
  # AR1
  s[1] ~ dnorm(0, sd =(tau*(1-rho^2))^(-0.5) )
  
  for (t in 2:nmonths) {
    s[t] ~ dnorm(rho*s[t-1], sd = sigma.e)
  }
  
  
  M ~ dpois(lambda)                    # Super-population size M
  
  
  for(t in 1:nmonths){                     # Loop over months
    
    logit(theta[t]) <- alpha + s[t] #~ dbeta(gammatilde[z[t]], betatilde[z[t]])
    N[t] ~ dbin(theta[t], M)            # 'Current' population size N
    
  }
  
  
  for (s1 in 1:S) {
    y[s1] ~ dbin(p[s1], N[Nindex[s1]])
    re[s1] ~ dnorm(0, sigma.re)
    logit(p[s1]) <- linpred[s1] + re[s1]
  }
  
})


Nst <- colSums(weeks_Opar, na.rm = TRUE)  #apply(Oparrots, 2, max, na.rm = TRUE)+100 # Inits for latent N
Nst = Nst[Nst !=0]
Mst <- max(Nst) + 500
TEmodel = nimbleModel(nimTE, constants = win.data, data = list(y =y, X=X), inits = list(sigma.re=1,pi_g = rep(1, length(pos)), lambda2 =1, tau_g = rgamma(length(pos), 1,1),b = rep(0, length(pos)),sigma.e = 1, rho = runif(1, -1,1), alpha = 1, re = rnorm(S,0,1), M = Mst,  N= Nst, lambda =  rgamma(1,10,0.1), beta = rep(0, numVars), s = rnorm(nweeks, 0,1), sigma = runif(1) )  )
TEmodel$calculate()

cTEmodel <- compileNimble(TEmodel)
TEconf = configureMCMC(TEmodel, enableWAIC = TRUE)

TEconf$monitors
TEconf$resetMonitors()
TEconf$addMonitors(c("lambda",  "pi_g", "theta", "beta", "M", "N", "sigma", "rho" ))

TEmcmc = buildMCMC(TEconf)
TEmod = compileNimble(TEmcmc, project = TEmodel, resetFunctions = TRUE)
TEmcmc.out <- runMCMC(TEmod, niter = 175000, nburnin = 75000, thin = 20, thin2 = 20,  samplesAsCodaMCMC = TRUE, nchains = 2, WAIC = TRUE, summary = TRUE) 
TEmcmc.out$summary$all.chains
Opar_AR1_WAIC = TEmcmc.out$WAIC
samples <- TEmcmc.out$samples

save(samples, file = "weeklyOparrots_TE_AR1_RE_BGLSS.Rdata")
#==================================================================================================================================================
