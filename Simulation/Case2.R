
library(nimble)
library(coda)
library(mcmcplots)
library(ggmcmc)
library(mcclust)
library(mcclust.ext)
library(pgdraw)
library(data.table)
library(MCMCvis)

# p approx 0.6

simRW1 <- function(nsurveys, nmonths, sigma, beta, ncov){
  
  S = nsurveys*nmonths
  M = rpois(1, lambda)
  
  s = numeric(nmonths)
  s[1] = 0.5#runif(1, -1,0.5)#runif(1, -10,10)
  
  for (t in 2:T) {
    s[t] <- rnorm(1,s[t-1], sigma)
  }
  
  theta <- ilogit(s)
  
  N <- rbinom(nmonths, M, theta)
  
  # Design matrix
  X = matrix(0, nrow = S, ncol = (ncov-2))  # first 4 continouous and last one categorical with 2 levels
  X[, 1:(ncov-2)] <- rnorm(S*(ncov-2))
  #X[, ncov] <- as.factor(sample(x=c(1,2, 3), size=S, replace=TRUE, prob=rep(1/3, 3))) 
  # cat_var = rbinom(S,1, 0.5)
  # X = cbind(X, cat_var)
  # X[, 5] = as.factor(X[,5])
  #head(X) 
  originalX <- X # backup the original  matrix of covariates
  X = as.data.frame(X)
  X <- model.matrix(~ ., X) # creates a design (or model) matrix
  A = data.frame(x5 = as.factor(sample(x=c(1,2, 3), size=S, replace=TRUE, prob=rep(1/3, 3))) )
  A <- model.matrix(~ ., A) # creates a design (or model) matrix
  B = data.frame(x6 = as.factor(sample(x=c(1,2), size=S, replace=TRUE, prob=rep(1/2, 2))) )
  B <- model.matrix(~ ., B) # creates a design (or model) matri
  X = cbind(X, B, A[,-1] )
  X= X[, -7]
  head(X)
  # Here the design matrix has a column of ones for the intercept
  
  # here the 3rd and 5th beta are non-significant 
  
  #true_betas <- c(c(0.1, 1, 0.4, 0.8, 1.5, 1))
  #true_gamma <- c(1,1,1,1,1,1)
  #true_gammabeta <- gamma*betas
  #re <- rnorm(S, 0, sigma.re)
  p <-  matrix(ilogit((X %*% beta)) ,nsurveys,nmonths)
  #head(p)
  #mean(apply(p, 2, mean)); plot(apply(p, 2, mean))
  
  Y <- matrix(0, nsurveys,nmonths)
  
  for (t in 1:nmonths) {
    #index[t] = sum(1- is.na(cparrots[,t])) # create indicator where to end loop
    Y[, t] <- rbinom(nsurveys,N[t], p[,t] )
  }
  
  y <- c(Y)
  
  return(list(y=y,Y=Y,X=X, N=N, p =p ))
  
}


J = nsurveys = 8#dim(cparrots)[1]
T = nmonths= 36#dim(cparrots)[2]

lambda = 100
sigma = 1
ncov1 = 7 # 5 continous and 2 categorical variables with 2 and 3 levels respectively
betas <- c(0.75, 1.25, 0.2, 2,0,-0.6, 0.5,-1, 0)# p approx 0.6
nsim = 50
#load("~/Objective priors for N-mixture models/Parrots/RW1_seed.Rdata")
load("~/Misc/RW1_seed.Rdata")
numbers = sdata[1:nsim]

# with RE 
betas_sum = array(0, dim = c(length(betas)-2, 5, nsim))
Ns_sum = array(0, dim = c(T, 5, nsim))
sigma_sum = sigmare_sum= matrix(0, nsim, 5)
lambda_sum = pi_g =  matrix(0, nsim, 5)
pi_g = array(0, dim = c(length(betas)-4, 5, nsim))
# effective sample size
ess = matrix(NA, nsim, 53)

for (n in 1:nsim) {
  
  print(n)
  set.seed(numbers[n])
  data <- simRW1(nsurveys, nmonths, sigma, betas, ncov1)
  mean(data$p)
  #plot(density(data$theta))
  # data$theta
  # plot(data$theta)
  plot(data$N, type = "b")
  # 
  
  numVars = length(betas)-2#ncov+1 #number of covariates (including the intercept)
  X = data$X[, -c(2,7)]
  #head(X)
  
  pos <- c(1,2,3,4,5,5) # The ‘pos’ vector gives the position of the factor corresponding to each component of beta
  m <- c(1,1,1,1,2,2) # The ‘m’ vector captures the same size group for all levels belonging to the same categorical covariate. 
  
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
  
  Nindex = rep(1:T, each = J ) 
  S = T*J
  
  win.data <- list(pos = pos, m =m, nb_levels = length(pos), S= S, tau_b0 = 3600, Nindex = Nindex,  nmonths = T,  numVars = numVars, L = length(CarArg$adj), adj = CarArg$adj, weights= CarArg$weights, num = CarArg$num)
  str(win.data)
  
  if(n==1){
    
    
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
      
      for (s in 1:S) {
        y[s] ~ dbin(p[s], N[Nindex[s]])
        re[s] ~ dnorm(0, sigma.re)
        logit(p[s]) <- linpred[s] + re[s]
      }
      
    })
    
    
    Nst <- apply(data$Y, 2, max, na.rm = TRUE)+100 # Inits for latent N
    Mst <- max(Nst)
    TEmodel = nimbleModel(nimTE, constants = win.data, data = list(y =data$y, X=X), inits = list(pi_g = rep(1, length(pos)), lambda2 =1, tau_g = rgamma(length(pos), 1,1), sigma.tu =1, b = rep(0, length(pos)), sigma.re=1, re = rnorm(S,0,1), M = Mst,  N= Nst, lambda =  rgamma(1,10,0.1), beta = rep(0, numVars), s = rnorm(T, 0,1), sigma = runif(1) )  )
    TEmodel$calculate()
    
    cTEmodel <- compileNimble(TEmodel)
    TEconf = configureMCMC(TEmodel, enableWAIC = TRUE)
    
    TEconf$monitors
    TEconf$resetMonitors()
    TEconf$addMonitors(c("lambda","sigma.re", "beta", "pi_g", "M", "N", "sigma"  ))
    #TEconf$addMonitors(c("lambda", "sdBeta", "sigma.re", "theta", "beta", "w", "M", "N", "sigma" ))
    
    TEmcmc = buildMCMC(TEconf)
    TEmod = compileNimble(TEmcmc, project = TEmodel, resetFunctions = TRUE)
    #TEmcmc.out <- runMCMC(TEmod, niter = 5000, nburnin = 1000, samplesAsCodaMCMC = TRUE, nchains = 1, summary = TRUE) 
    #TEmcmc.out <- runMCMC(TEmod, niter = 150000, nburnin = 50000, thin = 20, samplesAsCodaMCMC = TRUE, nchains = 2, WAIC = TRUE, summary = TRUE) 
    TEmcmc.out <- runMCMC(TEmod, niter = 55000, nburnin = 10000, thin = 20, samplesAsCodaMCMC = TRUE, nchains = 2, WAIC = TRUE, summary = TRUE) 
    #TEmcmc.out$summary$all.chains#[1:52,]
    #betas[-c(2,7)]
    betas_sum[,,n] = TEmcmc.out$summary$all.chains[(T+2): (T+(length(betas))-1 ), ]
    pi_g[,,n] = TEmcmc.out$summary$all.chains[46:50,]
    lambda_sum[n,] = TEmcmc.out$summary$all.chains["lambda",]
    Ns_sum[,,n] =  TEmcmc.out$summary$all.chains[2:(T+1),]
    sigma_sum[n,] = TEmcmc.out$summary$all.chains["sigma",]
    sigmare_sum[n, ] = TEmcmc.out$summary$all.chains["sigma.re",]
    samples = TEmcmc.out$samples
    ess[n,] = effectiveSize(samples)
    
  } else{
    
    Nst <- apply(data$Y, 2, max, na.rm = TRUE)+100 # Inits for latent N
    Mst <- max(Nst)
    
    cTEmodel$y= data$y
    cTEmodel$X = X
    cTEmodel$N = Nst
    cTEmodel$M = Mst
    cTEmodel$lambda = rgamma(1,10,0.1)
    cTEmodel$s = rnorm(T,0,1)
    cTEmodel$sigma = runif(1)
    cTEmodel$pi_g = rep(1, length(pos))
    cTEmodel$tau_g = rgamma(length(pos), 1,1)
    cTEmodel$lambda2 = 1#rgamma(1,10,0.1)
    
    #TEmcmc.out <- runMCMC(TEmod, niter = 150000, nburnin = 50000, thin = 20, samplesAsCodaMCMC = TRUE, nchains = 2, summary = TRUE) 
    TEmcmc.out <- runMCMC(TEmod, niter = 55000, nburnin = 10000, thin = 20, samplesAsCodaMCMC = TRUE, nchains = 2, WAIC = TRUE, summary = TRUE) 
    #TEmcmc.out$summary$all.chains#[1:52,]
    betas_sum[,,n] = TEmcmc.out$summary$all.chains[(T+2): (T+(length(betas))-1 ), ]
    pi_g[,,n] = TEmcmc.out$summary$all.chains[46:50,]
    lambda_sum[n,] = TEmcmc.out$summary$all.chains["lambda",]
    Ns_sum[,,n] =  TEmcmc.out$summary$all.chains[2:(T+1),]
    sigma_sum[n,] = TEmcmc.out$summary$all.chains["sigma",]
    sigmare_sum[n, ] = TEmcmc.out$summary$all.chains["sigma.re",]
    samples = TEmcmc.out$samples
    ess[n,] = effectiveSize(samples)
    
    
  }
  
  save(betas_sum, file= "wbetas_BGLSS_re_RW1_p06s.Rdata")
  save(Ns_sum, file= "wNs_BGLSS_re_RW1_p06s.Rdata")
  save(sigma_sum, file= "wsigma_BGLSS_re_RW1_p06s.Rdata")
  save(ess, file = "wess_BGLSS_re_RW1_p06s.Rdata")
  save(lambda_sum, file= "wlambda_BGLSS_re_RW1_p06s.Rdata")
  save(pi_g, file = "wpi_g_BGLSS_reRW1_p06s.Rdata")
  save(sigmare_sum, file = "wsigmare_BGLSS_reRW1_p06s.Rdata" )
  
}

# without RE
betas_sum = array(0, dim = c(length(betas), 5, nsim))
Ns_sum = array(0, dim = c(T, 5, nsim))
sigma_sum = sigmare_sum= matrix(0, nsim, 5)
lambda_sum = pi_g =  matrix(0, nsim, 5)
pi_g = array(0, dim = c(length(betas)-2, 5, nsim))
# effective sample size
ess = matrix(NA, nsim, 56)

for (n in 1:nsim) {
  
  print(n)
  set.seed(numbers[n])
  data <- simRW1(nsurveys, nmonths, sigma, betas, ncov1)
  mean(data$p)
  #plot(density(data$theta))
  # data$theta
  # plot(data$theta)
  plot(data$N, type = "b")
  # 
  
  numVars = length(betas)#-2#ncov+1 #number of covariates (including the intercept)
  X = data$X#[, -c(6,7)]
  head(X)
  pos <- c(1,2,3,4,5,6,7,7) # The ‘pos’ vector gives the position of the factor corresponding to each component of beta
  m <- c(1,1,1,1,1,1,2,2) # The ‘m’ vector captures the same size group for all levels belonging to the same categorical covariate. 
  
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
  
  Nindex = rep(1:T, each = J ) 
  S = T*J
  
  win.data <- list(pos = pos, m =m, nb_levels = length(pos), S= S, tau_b0 = 3600, Nindex = Nindex,  nmonths = T,  numVars = numVars, L = length(CarArg$adj), adj = CarArg$adj, weights= CarArg$weights, num = CarArg$num)
  str(win.data)
  
  if(n==1){
    
    
    nimTE <- nimbleCode({
      
      # Priors
      lambda ~ dgamma(0.01,0.01)
      #sigma.re ~ dunif(0,100)
      sigma ~ dunif(0,15)
      tau <- 1/sigma^2
      
      lambda2 ~ dgamma(0.001,0.001)
      # sigma.tu ~ dunif(0,10)
      #tau.tu <- 1/sigma.re^2
      
      for (l in 1:nb_levels){
        # using the "pos" vector allows us to affect the same  prior
        # inclusion probability for all levels within the same categorical
        pi_g[l] ~ dbern(p_g)
        pi_g_pos[l] <- pi_g[pos[l]]
        
        tau_g[l] ~ dgamma((m[pos[l]] + 1) / 2,pow(lambda2,2) / 2)
        # pos is used to ensure the same tau_coef value for all levels
        # within the same categorical covariate
        tau_g_pos[l] <- tau_g[pos[l]]
        
        b[l] ~ dnorm(0, var =  (1/pow(tau_g_pos[l],2)))
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
      
      for (s in 1:S) {
        y[s] ~ dbin(p[s], N[Nindex[s]])
        # re[s] ~ dnorm(0, sigma.re)
        logit(p[s]) <- linpred[s] #+ re[s]
      }
      
    })
    
    
    Nst <- apply(data$Y, 2, max, na.rm = TRUE)+100 # Inits for latent N
    Mst <- max(Nst)
    TEmodel = nimbleModel(nimTE, constants = win.data, data = list(y =data$y, X=X), inits = list(pi_g = rep(1, length(pos)), lambda2 =1, tau_g = rgamma(length(pos), 1,1), sigma.tu =1, b = rep(0, length(pos)), sigma.re=1, re = rnorm(S,0,1), M = Mst,  N= Nst, lambda =  rgamma(1,10,0.1), beta = rep(0, numVars), s = rnorm(T, 0,1), sigma = runif(1) )  )
    TEmodel$calculate()
    
    cTEmodel <- compileNimble(TEmodel)
    TEconf = configureMCMC(TEmodel, enableWAIC = TRUE)
    
    TEconf$monitors
    TEconf$resetMonitors()
    TEconf$addMonitors(c("lambda","beta", "pi_g", "M", "N", "sigma"  ))
    #TEconf$addMonitors(c("lambda", "sdBeta", "sigma.re", "theta", "beta", "w", "M", "N", "sigma" ))
    
    TEmcmc = buildMCMC(TEconf)
    TEmod = compileNimble(TEmcmc, project = TEmodel, resetFunctions = TRUE)
    #TEmcmc.out <- runMCMC(TEmod, niter = 5000, nburnin = 1000, samplesAsCodaMCMC = TRUE, nchains = 1, summary = TRUE) 
    #TEmcmc.out <- runMCMC(TEmod, niter = 150000, nburnin = 50000, thin = 20, samplesAsCodaMCMC = TRUE, nchains = 2, WAIC = TRUE, summary = TRUE) 
    TEmcmc.out <- runMCMC(TEmod, niter = 55000, nburnin = 10000, thin = 20, samplesAsCodaMCMC = TRUE, nchains = 2, WAIC = TRUE, summary = TRUE) 
    #TEmcmc.out$summary$all.chains#[1:52,]
    #betas[-c(6,7)]
    betas_sum[,,n] = TEmcmc.out$summary$all.chains[(T+2): (T+(length(betas))+1 ), ]
    pi_g[,,n] = TEmcmc.out$summary$all.chains[48:54,]
    lambda_sum[n,] = TEmcmc.out$summary$all.chains["lambda",]
    Ns_sum[,,n] =  TEmcmc.out$summary$all.chains[2:(T+1),]
    sigma_sum[n,] = TEmcmc.out$summary$all.chains["sigma",]
    samples = TEmcmc.out$samples
    ess[n,] = effectiveSize(samples)
    
  } else{
    
    Nst <- apply(data$Y, 2, max, na.rm = TRUE)+100 # Inits for latent N
    Mst <- max(Nst)
    
    cTEmodel$y= data$y
    cTEmodel$X = X
    cTEmodel$N = Nst
    cTEmodel$M = Mst
    cTEmodel$lambda = rgamma(1,10,0.1)
    cTEmodel$s = rnorm(T,0,1)
    cTEmodel$sigma = runif(1)
    cTEmodel$pi_g = rep(1, length(pos))
    cTEmodel$tau_g = rgamma(length(pos), 1,1)
    cTEmodel$lambda2 = 1#rgamma(1,10,0.1)
    
    #TEmcmc.out <- runMCMC(TEmod, niter = 150000, nburnin = 50000, thin = 20, samplesAsCodaMCMC = TRUE, nchains = 2, summary = TRUE) 
    TEmcmc.out <- runMCMC(TEmod, niter = 55000, nburnin = 10000, thin = 20, samplesAsCodaMCMC = TRUE, nchains = 2, WAIC = TRUE, summary = TRUE) 
    #TEmcmc.out$summary$all.chains#[1:52,]
    betas_sum[,,n] = TEmcmc.out$summary$all.chains[(T+2): (T+(length(betas))+1 ), ]
    pi_g[,,n] = TEmcmc.out$summary$all.chains[48:54,]
    lambda_sum[n,] = TEmcmc.out$summary$all.chains["lambda",]
    Ns_sum[,,n] =  TEmcmc.out$summary$all.chains[2:(T+1),]
    sigma_sum[n,] = TEmcmc.out$summary$all.chains["sigma",]
    samples = TEmcmc.out$samples
    ess[n,] = effectiveSize(samples)
    
    
  }
  
  save(betas_sum, file= "wbetas_BGLSS_RW1_p06s.Rdata")
  save(Ns_sum, file= "wNs_BGLSS_RW1_p06s.Rdata")
  save(sigma_sum, file= "wsigma_BGLSS_RW1_p06s.Rdata")
  save(ess, file = "wess_BGLSS_RW1_p06s.Rdata")
  save(lambda_sum, file= "wlambda_BGLSS_RW1_p06s.Rdata")
  save(pi_g, file = "wpi_g_BGLSS_RW1_p06s.Rdata")
  save(sigmare_sum, file = "wsigmare_BGLSS_RW1_p06s.Rdata" )
  
}


#===============================================================================================================================================

library(nimble)
library(coda)
library(mcmcplots)
library(ggmcmc)
library(mcclust)
library(mcclust.ext)
library(pgdraw)
library(data.table)
library(MCMCvis)

# p approx 0.3

# Generate data
simRW1 <- function(nsurveys, nmonths, sigma, beta, ncov){
  
  S = nsurveys*nmonths
  M = rpois(1, lambda)
  
  s = numeric(nmonths)
  s[1] = 0.5#runif(1, -1,0.5)#runif(1, -10,10)
  
  for (t in 2:T) {
    s[t] <- rnorm(1,s[t-1], sigma)
  }
  
  theta <- ilogit(s)
  
  N <- rbinom(nmonths, M, theta)
  
  # Design matrix
  X = matrix(0, nrow = S, ncol = (ncov-2))  # first 4 continouous and last one categorical with 2 levels
  X[, 1:(ncov-2)] <- rnorm(S*(ncov-2))
  #X[, ncov] <- as.factor(sample(x=c(1,2, 3), size=S, replace=TRUE, prob=rep(1/3, 3))) 
  # cat_var = rbinom(S,1, 0.5)
  # X = cbind(X, cat_var)
  # X[, 5] = as.factor(X[,5])
  #head(X) 
  originalX <- X # backup the original  matrix of covariates
  X = as.data.frame(X)
  X <- model.matrix(~ ., X) # creates a design (or model) matrix
  A = data.frame(x5 = as.factor(sample(x=c(1,2, 3), size=S, replace=TRUE, prob=rep(1/3, 3))) )
  A <- model.matrix(~ ., A) # creates a design (or model) matrix
  B = data.frame(x6 = as.factor(sample(x=c(1,2), size=S, replace=TRUE, prob=rep(1/2, 2))) )
  B <- model.matrix(~ ., B) # creates a design (or model) matri
  X = cbind(X, B, A[,-1] )
  X= X[, -7]
  head(X)
  # Here the design matrix has a column of ones for the intercept
  
  # here the 3rd and 5th beta are non-significant 
  
  #true_betas <- c(c(0.1, 1, 0.4, 0.8, 1.5, 1))
  #true_gamma <- c(1,1,1,1,1,1)
  #true_gammabeta <- gamma*betas
  #re <- rnorm(S, 0, sigma.re)
  p <-  matrix(ilogit((X %*% beta)) ,nsurveys,nmonths)
  #head(p)
  #mean(apply(p, 2, mean)); plot(apply(p, 2, mean))
  
  Y <- matrix(0, nsurveys,nmonths)
  
  for (t in 1:nmonths) {
    #index[t] = sum(1- is.na(cparrots[,t])) # create indicator where to end loop
    Y[, t] <- rbinom(nsurveys,N[t], p[,t] )
  }
  
  y <- c(Y)
  
  return(list(y=y,Y=Y,X=X, N=N, p =p ))
  
}


J = nsurveys = 8#dim(cparrots)[1]
T = nmonths= 36#dim(cparrots)[2]

lambda = 100
sigma = 1
ncov1 = 7 # 5 continous and 2 categorical variables with 2 and 3 levels respectively
betas <- c(-1.5, 1.25, 0.2, 2,0,-0.6, 0.5,-1, 0)# p approx 0.3
nsim = 50
#load("~/Objective priors for N-mixture models/Parrots/RW1_seed.Rdata")
load("~/Misc/RW1_seed.Rdata")
numbers = sdata[1:nsim]


# with RE 
betas_sum = array(0, dim = c(length(betas)-2, 5, nsim))
Ns_sum = array(0, dim = c(T, 5, nsim))
sigma_sum = sigmare_sum= matrix(0, nsim, 5)
lambda_sum = pi_g =  matrix(0, nsim, 5)
pi_g = array(0, dim = c(length(betas)-4, 5, nsim))
# effective sample size
ess = matrix(NA, nsim, 53)
for (n in 1:nsim) {
  
  print(n)
  set.seed(numbers[n])
  data <- simRW1(nsurveys, nmonths, sigma, betas, ncov1)
  mean(data$p)
  #plot(density(data$theta))
  # data$theta
  # plot(data$theta)
  plot(data$N, type = "b")
  # 
  
  numVars = length(betas)-2#ncov+1 #number of covariates (including the intercept)
  X = data$X[, -c(2,7)]
  
  pos <- c(1,2,3,4,5,5) # The ‘pos’ vector gives the position of the factor corresponding to each component of beta
  m <- c(1,1,1,1,2,2) # The ‘m’ vector captures the same size group for all levels belonging to the same categorical covariate. 
  
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
  
  Nindex = rep(1:T, each = J ) 
  S = T*J
  
  win.data <- list(pos = pos, m =m, nb_levels = length(pos), S= S, tau_b0 = 3600, Nindex = Nindex,  nmonths = T,  numVars = numVars, L = length(CarArg$adj), adj = CarArg$adj, weights= CarArg$weights, num = CarArg$num)
  str(win.data)
  
  if(n==1){
    
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
      
      for (s in 1:S) {
        y[s] ~ dbin(p[s], N[Nindex[s]])
        re[s] ~ dnorm(0, sigma.re)
        logit(p[s]) <- linpred[s] + re[s]
      }
      
    })
    
    
    Nst <- apply(data$Y, 2, max, na.rm = TRUE)+100 # Inits for latent N
    Mst <- max(Nst)
    TEmodel = nimbleModel(nimTE, constants = win.data, data = list(y =data$y, X=X), inits = list(pi_g = rep(1, length(pos)), lambda2 =1, tau_g = rgamma(length(pos), 1,1), sigma.tu =1, b = rep(0, length(pos)), sigma.re=1, re = rnorm(S,0,1), M = Mst,  N= Nst, lambda =  rgamma(1,10,0.1), beta = rep(0, numVars), s = rnorm(T, 0,1), sigma = runif(1) )  )
    TEmodel$calculate()
    
    cTEmodel <- compileNimble(TEmodel)
    TEconf = configureMCMC(TEmodel, enableWAIC = TRUE)
    
    TEconf$monitors
    TEconf$resetMonitors()
    TEconf$addMonitors(c("lambda","sigma.re", "beta", "pi_g", "M", "N", "sigma"  ))
    #TEconf$addMonitors(c("lambda", "sdBeta", "sigma.re", "theta", "beta", "w", "M", "N", "sigma" ))
    
    TEmcmc = buildMCMC(TEconf)
    TEmod = compileNimble(TEmcmc, project = TEmodel, resetFunctions = TRUE)
    #TEmcmc.out <- runMCMC(TEmod, niter = 5000, nburnin = 1000, samplesAsCodaMCMC = TRUE, nchains = 1, summary = TRUE) 
    TEmcmc.out <- runMCMC(TEmod, niter = 150000, nburnin = 50000, thin = 20, samplesAsCodaMCMC = TRUE, nchains = 2, WAIC = TRUE, summary = TRUE) 
    #TEmcmc.out <- runMCMC(TEmod, niter = 55000, nburnin = 10000, thin = 20, samplesAsCodaMCMC = TRUE, nchains = 2, WAIC = TRUE, summary = TRUE) 
    # TEmcmc.out$summary$all.chains#[1:52,]
    betas_sum[,,n] = TEmcmc.out$summary$all.chains[(T+2): (T+(length(betas))-1 ), ]
    pi_g[,,n] = TEmcmc.out$summary$all.chains[46:50,]
    lambda_sum[n,] = TEmcmc.out$summary$all.chains["lambda",]
    Ns_sum[,,n] =  TEmcmc.out$summary$all.chains[2:(T+1),]
    sigma_sum[n,] = TEmcmc.out$summary$all.chains["sigma",]
    sigmare_sum[n, ] = TEmcmc.out$summary$all.chains["sigma.re",]
    samples = TEmcmc.out$samples
    ess[n,] = effectiveSize(samples)
    
  } else{
    
    Nst <- apply(data$Y, 2, max, na.rm = TRUE)+100 # Inits for latent N
    Mst <- max(Nst)
    
    cTEmodel$y= data$y
    cTEmodel$X = X
    cTEmodel$N = Nst
    cTEmodel$M = Mst
    cTEmodel$lambda = rgamma(1,10,0.1)
    cTEmodel$s = rnorm(T,0,1)
    cTEmodel$sigma = runif(1)
    cTEmodel$pi_g = rep(1, length(pos))
    cTEmodel$tau_g = rgamma(length(pos), 1,1)
    cTEmodel$lambda2 = 1#rgamma(1,10,0.1)
    
    TEmcmc.out <- runMCMC(TEmod, niter = 150000, nburnin = 50000, thin = 20, samplesAsCodaMCMC = TRUE, nchains = 2, summary = TRUE) 
    #TEmcmc.out <- runMCMC(TEmod, niter = 55000, nburnin = 10000, thin = 20, samplesAsCodaMCMC = TRUE, nchains = 2, WAIC = TRUE, summary = TRUE) 
    #TEmcmc.out$summary$all.chains#[1:52,]
    betas_sum[,,n] = TEmcmc.out$summary$all.chains[(T+2): (T+(length(betas))-1 ), ]
    pi_g[,,n] = TEmcmc.out$summary$all.chains[46:50,]
    lambda_sum[n,] = TEmcmc.out$summary$all.chains["lambda",]
    Ns_sum[,,n] =  TEmcmc.out$summary$all.chains[2:(T+1),]
    sigma_sum[n,] = TEmcmc.out$summary$all.chains["sigma",]
    sigmare_sum[n, ] = TEmcmc.out$summary$all.chains["sigma.re",]
    samples = TEmcmc.out$samples
    ess[n,] = effectiveSize(samples)
    
    
  }
  
  save(betas_sum, file= "wbetas_BGLSS_re_RW1_p03s.Rdata")
  save(Ns_sum, file= "wNs_BGLSS_re_RW1_p03s.Rdata")
  save(sigma_sum, file= "wsigma_BGLSS_re_RW1_p03s.Rdata")
  save(ess, file = "wess_BGLSS_re_RW1_p03s.Rdata")
  save(lambda_sum, file= "wlambda_BGLSS_re_RW1_p03s.Rdata")
  save(pi_g, file = "wpi_g_BGLSS_reRW1_p03s.Rdata")
  save(sigmare_sum, file = "wsigmare_BGLSS_reRW1_p03s.Rdata" )
  
}


# without RE 
betas_sum = array(0, dim = c(length(betas), 5, nsim))
Ns_sum = array(0, dim = c(T, 5, nsim))
sigma_sum = sigmare_sum= matrix(0, nsim, 5)
lambda_sum = pi_g =  matrix(0, nsim, 5)
pi_g = array(0, dim = c(length(betas)-2, 5, nsim))
# effective sample size
ess = matrix(NA, nsim, 56)

for (n in 1:nsim) {
  
  print(n)
  set.seed(numbers[n])
  data <- simRW1(nsurveys, nmonths, sigma, betas, ncov1)
  mean(data$p)
  #plot(density(data$theta))
  # data$theta
  # plot(data$theta)
  plot(data$N, type = "b")
  # 
  
  numVars = length(betas)#-2#ncov+1 #number of covariates (including the intercept)
  X = data$X#[, -c(6,7)]
  head(X)
  pos <- c(1,2,3,4,5,6,7,7) # The ‘pos’ vector gives the position of the factor corresponding to each component of beta
  m <- c(1,1,1,1,1,1,2,2) # The ‘m’ vector captures the same size group for all levels belonging to the same categorical covariate. 
  
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
  
  Nindex = rep(1:T, each = J ) 
  S = T*J
  
  win.data <- list(pos = pos, m =m, nb_levels = length(pos), S= S, tau_b0 = 3600, Nindex = Nindex,  nmonths = T,  numVars = numVars, L = length(CarArg$adj), adj = CarArg$adj, weights= CarArg$weights, num = CarArg$num)
  str(win.data)
  
  if(n==1){
    
    nimTE <- nimbleCode({
      
      # Priors
      lambda ~ dgamma(0.01,0.01)
      # sigma.re ~ dunif(0,100)
      sigma ~ dunif(0,15)
      tau <- 1/sigma^2
      
      lambda2 ~ dgamma(0.001,0.001)
      # sigma.tu ~ dunif(0,10)
      # tau.tu <- 1/sigma.re^2
      
      for (l in 1:nb_levels){
        # using the "pos" vector allows us to affect the same  prior
        # inclusion probability for all levels within the same categorical
        pi_g[l] ~ dbern(p_g)
        pi_g_pos[l] <- pi_g[pos[l]]
        
        tau_g[l] ~ dgamma((m[pos[l]] + 1) / 2,pow(lambda2,2) / 2)
        # pos is used to ensure the same tau_coef value for all levels
        # within the same categorical covariate
        tau_g_pos[l] <- tau_g[pos[l]]
        
        b[l] ~ dnorm(0, var =  (1/pow(tau_g_pos[l],2)))
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
      
      for (s in 1:S) {
        y[s] ~ dbin(p[s], N[Nindex[s]])
        #   re[s] ~ dnorm(0, sigma.re)
        logit(p[s]) <- linpred[s] #+ re[s]
      }
      
    })
    
    
    Nst <- apply(data$Y, 2, max, na.rm = TRUE)+100 # Inits for latent N
    Mst <- max(Nst)
    TEmodel = nimbleModel(nimTE, constants = win.data, data = list(y =data$y, X=X), inits = list(pi_g = rep(1, length(pos)), lambda2 =1, tau_g = rgamma(length(pos), 1,1), sigma.tu =1, b = rep(0, length(pos)), sigma.re=1, re = rnorm(S,0,1), M = Mst,  N= Nst, lambda =  rgamma(1,10,0.1), beta = rep(0, numVars), s = rnorm(T, 0,1), sigma = runif(1) )  )
    TEmodel$calculate()
    
    cTEmodel <- compileNimble(TEmodel)
    TEconf = configureMCMC(TEmodel, enableWAIC = TRUE)
    
    TEconf$monitors
    TEconf$resetMonitors()
    TEconf$addMonitors(c("lambda", "beta", "pi_g", "M", "N", "sigma"  ))
    #TEconf$addMonitors(c("lambda", "sdBeta", "sigma.re", "theta", "beta", "w", "M", "N", "sigma" ))
    
    TEmcmc = buildMCMC(TEconf)
    TEmod = compileNimble(TEmcmc, project = TEmodel, resetFunctions = TRUE)
    #TEmcmc.out <- runMCMC(TEmod, niter = 5000, nburnin = 1000, samplesAsCodaMCMC = TRUE, nchains = 2, summary = TRUE) 
    TEmcmc.out <- runMCMC(TEmod, niter = 150000, nburnin = 50000, thin = 20, samplesAsCodaMCMC = TRUE, nchains = 2, WAIC = TRUE, summary = TRUE) 
    #TEmcmc.out <- runMCMC(TEmod, niter = 55000, nburnin = 10000, thin = 20, samplesAsCodaMCMC = TRUE, nchains = 2, WAIC = TRUE, summary = TRUE) 
    # TEmcmc.out$summary$all.chains#[1:52,]
    betas_sum[,,n] = TEmcmc.out$summary$all.chains[(T+2): (T+(length(betas))+1 ), ]
    pi_g[,,n] = TEmcmc.out$summary$all.chains[48:54,]
    lambda_sum[n,] = TEmcmc.out$summary$all.chains["lambda",]
    Ns_sum[,,n] =  TEmcmc.out$summary$all.chains[2:(T+1),]
    sigma_sum[n,] = TEmcmc.out$summary$all.chains["sigma",]
    samples = TEmcmc.out$samples
    ess[n,] = effectiveSize(samples)
    
  } else{
    
    Nst <- apply(data$Y, 2, max, na.rm = TRUE)+100 # Inits for latent N
    Mst <- max(Nst)
    
    cTEmodel$y= data$y
    cTEmodel$X = X
    cTEmodel$N = Nst
    cTEmodel$M = Mst
    cTEmodel$lambda = rgamma(1,10,0.1)
    cTEmodel$s = rnorm(T,0,1)
    cTEmodel$sigma = runif(1)
    cTEmodel$pi_g = rep(1, length(pos))
    cTEmodel$tau_g = rgamma(length(pos), 1,1)
    cTEmodel$lambda2 = 1#rgamma(1,10,0.1)
    
    TEmcmc.out <- runMCMC(TEmod, niter = 150000, nburnin = 50000, thin = 20, samplesAsCodaMCMC = TRUE, nchains = 2, summary = TRUE) 
    #TEmcmc.out <- runMCMC(TEmod, niter = 55000, nburnin = 10000, thin = 20, samplesAsCodaMCMC = TRUE, nchains = 2, WAIC = TRUE, summary = TRUE) 
    #TEmcmc.out$summary$all.chains#[1:52,]
    betas_sum[,,n] = TEmcmc.out$summary$all.chains[(T+2): (T+(length(betas))+1 ), ]
    pi_g[,,n] =TEmcmc.out$summary$all.chains[48:54,]
    lambda_sum[n,] = TEmcmc.out$summary$all.chains["lambda",]
    Ns_sum[,,n] =  TEmcmc.out$summary$all.chains[2:(T+1),]
    sigma_sum[n,] = TEmcmc.out$summary$all.chains["sigma",]
    samples = TEmcmc.out$samples
    ess[n,] = effectiveSize(samples)
    
    
  }
  
  save(betas_sum, file= "wbetas_BGLSS_RW1_p03s.Rdata")
  save(Ns_sum, file= "wNs_BGLSS_RW1_p03s.Rdata")
  save(sigma_sum, file= "wsigma_BGLSS_RW1_p03s.Rdata")
  save(ess, file = "wess_BGLSS_RW1_p03s.Rdata")
  save(lambda_sum, file= "wlambda_BGLSS_RW1_p03s.Rdata")
  save(pi_g, file = "wpi_g_BGLSS_RW1_p03s.Rdata")
  
}

