# Rscript where I show all models work when p apprx 0.6
library(nimble)
library(coda)
library(mcclust)
library(mcclust.ext)

#====================================================================================================================================================
# DP 

# Generate data

# Generate data
simDP <- function(nsurveys, nmonths, beta, ncov){
  
  
  M = rpois(1, lambda)
  
  n = nmonths/2
  S = nsurveys*nmonths
  x = rbeta(n, 10, 10);mean(x); var(x)
  z = rbeta(n, 10, 1); mean(z); var(z)
  #theta <- list(y = c(x, z)) # to clusters
  theta <- list(y =  sample(c(x, z)) ) # to clusters
  
  N <- rbinom(nmonths, M, theta$y)
  
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


J = nsurveys=  8#dim(cparrots)[1]
T = nmonths = 36#dim(cparrots)[2]
lambda = 100
#betas <- c(0.5, 0.3, 0.1, 0.5, 0.1, 1, 0, 0) # p approx 0.6
#true_betas <- c(-1.5, -0.8, -0.6, -0.1, -0.8, -0.1, 0, 0)# p approx 0.3
ncov1 = 7 # 5 continous and 2 categorical variables with 2 and 3 levels respectively
betas <- c(0.75, 1.25, 0.2, 2,0,-0.6, 0.5,-1, 0)# p approx 0.6
nsim = 50
numbers = 1:nsim


betas_sum = array(0, dim = c(length(betas), 5, nsim))
Ns_sum = array(0, dim = c(T, 5, nsim))
lambda_sum = matrix(0, nsim, 5)
zs_sum = matrix(0, nsim, 2)
# effective sample size
ess = matrix(NA, nsim, 46)
Nindex = rep(1:T, each = J ) 

for (n in 1:nsim) {
  
  print(n)  
  set.seed(numbers[n])
  #set.seed(1993)
  data = simDP(J, T,  betas, ncov1)
  mean(data$p)
  plot(data$N, type = "b")
  
  numVars = length(betas)#-2#ncov+1 #number of covariates (including the intercept)
  X =  data$X#[,c(1:4,7,8)]
  dim(X)
  
  ncov <- numVars-1 # number of explanatory variables in the model excluding the intercept
  column_covariate = 1:ncov
  
  # true or false if the covaraites is categorical or numerical - I only have numerical covariates to begin with
  classCovariates <- rep(T, ncov)
  #classCovariates[match(fac_covariates, column_covariate)] <- F
  
  
  indexes_covariates <- c()
  indexes_covariates[1] <- 1
  
  if(any(column_covariate != 0)){
    
    indexes_covariates <- c()
    indexes_covariates[1] <- 1
    k <- 2
    for (i in 1:ncov) {
      if(classCovariates[i]){
        indexes_covariates[k] <- i + 1
        k <- k + 1
      } else {
        num_levels <- length(unique(X[,i]))
        indexes_covariates[k + 0:(num_levels-2)] <- i + 1
        k <- k + num_levels - 1
      }
    }
    
  }
  
  indexes_covariates # if only numerical covariates this is just 1:(ncov+1)
  # I also think tha the 1 is for the intercept
  
  # compute C matrix given only numerical covariates 
  C <- matrix(0, nrow = ncol(X) - 1, ncol = ncol(X) - 1)
  l <- 0 # l loop through the covariates 
  for (i in 1:ncov) {
    if(classCovariates[i]){ # if it's a numerical covariates
      C[l + 1, l + 1] <- 1
      l <- l + 1
    } else { # if it's a categorical covariates
      num_levels <- length(unique(originalX[,i]))
      C[l + 1:(num_levels-1), l + 1:(num_levels-1)] <- .5 * (diag(1, nrow = (num_levels-1)) + 
                                                               matrix(1, nrow = (num_levels-1), ncol = (num_levels-1)))
      l <- l + num_levels - 1
    }
  }
  
  C
  Sigma <- cov(X[,-1, drop = F])
  sumCSigma <- sum(C * Sigma)
  
  d_bar <- ncov      
  d_bar <- ifelse(d_bar <= ncov, d_bar, 2)
  
  mu_psi = 0#0.5
  phi_mu =4
  phi_beta = 1/4
  
  # Prior parameters for beta, beta ~ MVN(b, B)   
  b_psi <- c(mu_psi, rep(0, ncol(X) - 1))
  sigma_beta_psi <- C * phi_beta  
  b_psi #b
  
  B_psi <- matrix(0, nrow = ncol(X), ncol = ncol(X))
  B_psi[1, 1] <- phi_mu
  B_psi[2:(ncol(X)), 2:(ncol(X))] <- sigma_beta_psi
  B_psi # B
  
  S = T*J
  win.data <- list(S= S, tau_b0 = 3600, Nindex = Nindex, b = b_psi, B = B_psi,  nmonths = T,  numVars = numVars, K=36)
  str(win.data)
  
  if(n==1){
    
    nimTE <- nimbleCode({
      
      # Priors
      lambda ~ dgamma(0.01,0.01)
      
      beta[1:numVars] ~ dmnorm(b[1:numVars], B[1:numVars, 1:numVars]) 
      linpred[1:S] <- (X[1:S, 1:numVars] %*% beta[1:numVars])[1:S,1]
      
      # DP set up
      for(i in 1:K) {
        gammatilde[i] ~  dgamma(2,0.1)#dgamma(0.1,0.1)
        betatilde[i] ~  dgamma(2,0.1)#dgamma(0.1,0.1)
      }
      
      z[1:nmonths] ~ dCRP(alpha, size = nmonths)
      alpha ~ dgamma(1, 1)
      
      
      M ~ dpois(lambda)                    # Super-population size M
      
      
      for(t in 1:nmonths){                     # Loop over months
        
        theta[t] ~ dbeta(gammatilde[z[t]], betatilde[z[t]])
        N[t] ~ dbin(theta[t], M)            # 'Current' population size N
        
      }
      
      for (s in 1:S) {
        y[s] ~ dbin(p[s], N[Nindex[s]])
        #  re[s] ~ dnorm(0, sigma.re)
        logit(p[s]) <- linpred[s] #+ re[s]
      }
      
    })
    
    
    Nst <- apply(data$Y, 2, max, na.rm = TRUE)+100 # Inits for latent N
    Mst <- max(Nst)
    TEmodel = nimbleModel(nimTE, constants = win.data, data = list(y =data$y, X=X), inits = list(M = Mst,N =Nst, lambda =  rgamma(1,10,0.1), beta = rep(0, numVars),  theta = rbeta(T,0.5, 0.5), z = sample(1:win.data$K, size = T, replace = TRUE), alpha = 1, gammatilde = rgamma(win.data$K, 0.5, 1) ,betatilde = rgamma(win.data$K, 0.5, 1)  ))
    TEmodel$calculate()
    
    cTEmodel <- compileNimble(TEmodel)
    TEconf = configureMCMC(TEmodel, enableWAIC = TRUE)
    
    TEconf$monitors
    TEconf$resetMonitors()
    TEconf$addMonitors(c("lambda", "theta", "beta", "M", "N", "z" ))
    
    TEmcmc = buildMCMC(TEconf)
    TEmod = compileNimble(TEmcmc, project = TEmodel, resetFunctions = TRUE)
    #TEmcmc.out <- runMCMC(TEmod, niter = 5000, nburnin = 1000, samplesAsCodaMCMC = TRUE, nchains = 1, WAIC = TRUE, summary = TRUE) 
    #TEmcmc.out <- runMCMC(TEmod, niter = 105000, nburnin = 15000, thin = 20, samplesAsCodaMCMC = TRUE, nchains = 1, summary = TRUE) 
    TEmcmc.out <- runMCMC(TEmod, niter = 55000, nburnin = 10000, thin = 20, samplesAsCodaMCMC = TRUE, nchains = 2, WAIC = TRUE, summary = TRUE) 
    betas_sum[,,n] = TEmcmc.out$summary$all.chains[(T+2): (T+(length(betas))+1 ), ]
    lambda_sum[n,] = TEmcmc.out$summary$all.chains["lambda",]
    Ns_sum[,,n] =  TEmcmc.out$summary$all.chains[2:(T+1),]
    samples = TEmcmc.out$samples
    
    ess[n,] = effectiveSize(samples[, c(1:46)])
    
    samples = mcmc(do.call(rbind, samples))
    zs <- samples[,84:119]
    # thetas <- MCMCsummary(as.matrix(samples[ ,47:82]),Rhat = F, n.eff = F)
    # d = data.frame(thetas[,1])
    
    psm=comp.psm(zs)
    z.VI=minVI(psm,zs,method=("all"),include.greedy=F)
    zs_sum[n,] = summary(z.VI)$t[[1]][1:2]
    
    
  } else {
    
    Nst <- apply(data$Y, 2, max, na.rm = TRUE)+100 # Inits for latent N
    Mst <- max(Nst)
    
    cTEmodel$y= data$y
    cTEmodel$X = X
    cTEmodel$N = Nst
    cTEmodel$M = Mst
    cTEmodel$lambda = rgamma(1,10,0.1)
    
    TEmcmc.out <- runMCMC(TEmod, niter = 55000, nburnin = 10000, thin = 20, samplesAsCodaMCMC = TRUE, nchains = 2, WAIC = TRUE, summary = TRUE) 
    betas_sum[,,n] = TEmcmc.out$summary$all.chains[(T+2): (T+(length(betas))+1 ), ]
    lambda_sum[n,] = TEmcmc.out$summary$all.chains["lambda",]
    Ns_sum[,,n] =  TEmcmc.out$summary$all.chains[2:(T+1),]
    samples = TEmcmc.out$samples
    
    ess[n,] = effectiveSize(samples[, c(1:46)])
    
    samples = mcmc(do.call(rbind, samples))
    zs <- samples[,84:119]
    # thetas <- MCMCsummary(as.matrix(samples[ ,47:82]),Rhat = F, n.eff = F)
    # d = data.frame(thetas[,1])
    
    psm=comp.psm(zs)
    z.VI=minVI(psm,zs,method=("all"),include.greedy=F)
    zs_sum[n,] = summary(z.VI)$t[[1]][1:2]
    
    
    
  }
  
  save(betas_sum, file= "betas_TE_DP_p06s.Rdata")
  save(Ns_sum, file= "Ns_TE_DP_p06s.Rdata")
  save(ess, file = "ess_TE_DP_p06s.Rdata")
  save(lambda_sum, file= "lambda_TE_DP_p06s.Rdata")
  save(zs_sum, file= "zs_TE_DP_BVS_p06s.Rdata")
  
}


# p approx 0.3

lambda = 100
ncov1 = 7
betas <- c(-1.5, 1.25, 0.2, 2,0,-0.6, 0.5,-1, 0)# p approx 0.3
nsim = 50
numbers = 1:nsim


betas_sum = array(0, dim = c(length(betas), 5, nsim))
Ns_sum = array(0, dim = c(T, 5, nsim))
lambda_sum = matrix(0, nsim, 5)
zs_sum = matrix(0, nsim, 2)
# effective sample size
ess = matrix(NA, nsim, 46)
Nindex = rep(1:T, each = J ) 

for (n in 1:nsim) {
  
  print(n)  
  set.seed(numbers[n])
  #set.seed(1993)
  data = simDP(J, T,  betas, ncov1)
  mean(data$p)
  plot(data$N, type = "b")
  
  numVars = length(betas)#-2#ncov+1 #number of covariates (including the intercept)
  X =  data$X#[,c(1:4,7,8)]
  dim(X)
  
  ncov <- numVars-1 # number of explanatory variables in the model excluding the intercept
  column_covariate = 1:ncov
  
  # true or false if the covaraites is categorical or numerical - I only have numerical covariates to begin with
  classCovariates <- rep(T, ncov)
  #classCovariates[match(fac_covariates, column_covariate)] <- F
  
  
  indexes_covariates <- c()
  indexes_covariates[1] <- 1
  
  if(any(column_covariate != 0)){
    
    indexes_covariates <- c()
    indexes_covariates[1] <- 1
    k <- 2
    for (i in 1:ncov) {
      if(classCovariates[i]){
        indexes_covariates[k] <- i + 1
        k <- k + 1
      } else {
        num_levels <- length(unique(X[,i]))
        indexes_covariates[k + 0:(num_levels-2)] <- i + 1
        k <- k + num_levels - 1
      }
    }
    
  }
  
  indexes_covariates # if only numerical covariates this is just 1:(ncov+1)
  # I also think tha the 1 is for the intercept
  
  # compute C matrix given only numerical covariates 
  C <- matrix(0, nrow = ncol(X) - 1, ncol = ncol(X) - 1)
  l <- 0 # l loop through the covariates 
  for (i in 1:ncov) {
    if(classCovariates[i]){ # if it's a numerical covariates
      C[l + 1, l + 1] <- 1
      l <- l + 1
    } else { # if it's a categorical covariates
      num_levels <- length(unique(originalX[,i]))
      C[l + 1:(num_levels-1), l + 1:(num_levels-1)] <- .5 * (diag(1, nrow = (num_levels-1)) + 
                                                               matrix(1, nrow = (num_levels-1), ncol = (num_levels-1)))
      l <- l + num_levels - 1
    }
  }
  
  C
  Sigma <- cov(X[,-1, drop = F])
  sumCSigma <- sum(C * Sigma)
  
  d_bar <- ncov      
  d_bar <- ifelse(d_bar <= ncov, d_bar, 2)
  
  mu_psi = 0#0.5
  phi_mu =4
  phi_beta = 1/4
  
  # Prior parameters for beta, beta ~ MVN(b, B)   
  b_psi <- c(mu_psi, rep(0, ncol(X) - 1))
  sigma_beta_psi <- C * phi_beta  
  b_psi #b
  
  B_psi <- matrix(0, nrow = ncol(X), ncol = ncol(X))
  B_psi[1, 1] <- phi_mu
  B_psi[2:(ncol(X)), 2:(ncol(X))] <- sigma_beta_psi
  B_psi # B
  
  S = T*J
  win.data <- list(S= S, tau_b0 = 3600, Nindex = Nindex, b = b_psi, B = B_psi,  nmonths = T,  numVars = numVars, K=36)
  str(win.data)
  
  if(n==1){
    
    nimTE <- nimbleCode({
      
      # Priors
      lambda ~ dgamma(0.01,0.01)
      
      beta[1:numVars] ~ dmnorm(b[1:numVars], B[1:numVars, 1:numVars]) 
      linpred[1:S] <- (X[1:S, 1:numVars] %*% beta[1:numVars])[1:S,1]
      
      # DP set up
      for(i in 1:K) {
        gammatilde[i] ~  dgamma(2,0.1)#dgamma(0.1,0.1)
        betatilde[i] ~  dgamma(2,0.1)#dgamma(0.1,0.1)
      }
      
      z[1:nmonths] ~ dCRP(alpha, size = nmonths)
      alpha ~ dgamma(1, 1)
      
      
      M ~ dpois(lambda)                    # Super-population size M
      
      
      for(t in 1:nmonths){                     # Loop over months
        
        theta[t] ~ dbeta(gammatilde[z[t]], betatilde[z[t]])
        N[t] ~ dbin(theta[t], M)            # 'Current' population size N
        
      }
      
      for (s in 1:S) {
        y[s] ~ dbin(p[s], N[Nindex[s]])
        #  re[s] ~ dnorm(0, sigma.re)
        logit(p[s]) <- linpred[s] #+ re[s]
      }
      
    })
    
    
    Nst <- apply(data$Y, 2, max, na.rm = TRUE)+100 # Inits for latent N
    Mst <- max(Nst)
    TEmodel = nimbleModel(nimTE, constants = win.data, data = list(y =data$y, X=X), inits = list(M = Mst,N =Nst, lambda =  rgamma(1,10,0.1), beta = rep(0, numVars),  theta = rbeta(T,0.5, 0.5), z = sample(1:win.data$K, size = T, replace = TRUE), alpha = 1, gammatilde = rgamma(win.data$K, 0.5, 1) ,betatilde = rgamma(win.data$K, 0.5, 1)  ))
    TEmodel$calculate()
    
    cTEmodel <- compileNimble(TEmodel)
    TEconf = configureMCMC(TEmodel, enableWAIC = TRUE)
    
    TEconf$monitors
    TEconf$resetMonitors()
    TEconf$addMonitors(c("lambda", "theta", "beta", "M", "N", "z" ))
    
    TEmcmc = buildMCMC(TEconf)
    TEmod = compileNimble(TEmcmc, project = TEmodel, resetFunctions = TRUE)
    #TEmcmc.out <- runMCMC(TEmod, niter = 5000, nburnin = 1000, samplesAsCodaMCMC = TRUE, nchains = 1, WAIC = TRUE, summary = TRUE) 
    #TEmcmc.out <- runMCMC(TEmod, niter = 150000, nburnin = 35000, thin = 20, samplesAsCodaMCMC = TRUE, nchains = 2, WAIC = TRUE, summary = TRUE) 
     TEmcmc.out <- runMCMC(TEmod, niter = 55000, nburnin = 10000, thin = 20, samplesAsCodaMCMC = TRUE, nchains = 2, WAIC = TRUE, summary = TRUE) 
    betas_sum[,,n] = TEmcmc.out$summary$all.chains[(T+2): (T+(length(betas))+1 ), ]
    lambda_sum[n,] = TEmcmc.out$summary$all.chains["lambda",]
    Ns_sum[,,n] =  TEmcmc.out$summary$all.chains[2:(T+1),]
    samples = TEmcmc.out$samples
    
    ess[n,] = effectiveSize(samples[, c(1:46)])
    
    samples = mcmc(do.call(rbind, samples))
    zs <- samples[,84:119]
    # thetas <- MCMCsummary(as.matrix(samples[ ,47:82]),Rhat = F, n.eff = F)
    # d = data.frame(thetas[,1])
    
    psm=comp.psm(zs)
    z.VI=minVI(psm,zs,method=("all"),include.greedy=F)
    zs_sum[n,] = summary(z.VI)$t[[1]][1:2]
    
    
  } else {
    
    Nst <- apply(data$Y, 2, max, na.rm = TRUE)+100 # Inits for latent N
    Mst <- max(Nst)
    
    cTEmodel$y= data$y
    cTEmodel$X = X
    cTEmodel$N = Nst
    cTEmodel$M = Mst
    cTEmodel$lambda = rgamma(1,10,0.1)
    
    TEmcmc.out <- runMCMC(TEmod, niter = 55000, nburnin = 10000, thin = 20, samplesAsCodaMCMC = TRUE, nchains = 2, WAIC = TRUE, summary = TRUE) 
    #TEmcmc.out <- runMCMC(TEmod, niter = 150000, nburnin = 35000, thin = 20, samplesAsCodaMCMC = TRUE, nchains = 2, WAIC = TRUE, summary = TRUE) 
    betas_sum[,,n] = TEmcmc.out$summary$all.chains[(T+2): (T+(length(betas))+1 ), ]
    lambda_sum[n,] = TEmcmc.out$summary$all.chains["lambda",]
    Ns_sum[,,n] =  TEmcmc.out$summary$all.chains[2:(T+1),]
    samples = TEmcmc.out$samples
    
    ess[n,] = effectiveSize(samples[, c(1:46)])
    
    samples = mcmc(do.call(rbind, samples))
    zs <- samples[,84:119]
    # thetas <- MCMCsummary(as.matrix(samples[ ,47:82]),Rhat = F, n.eff = F)
    # d = data.frame(thetas[,1])
    
    psm=comp.psm(zs)
    z.VI=minVI(psm,zs,method=("all"),include.greedy=F)
    zs_sum[n,] = summary(z.VI)$t[[1]][1:2]
    
    
    
  }
  
  save(betas_sum, file= "betas_TE_DP_p03s.Rdata")
  save(Ns_sum, file= "Ns_TE_DP_p03s.Rdata")
  save(ess, file = "ess_TE_DP_p03s.Rdata")
  save(lambda_sum, file= "lambda_TE_DP_p03s.Rdata")
  save(zs_sum, file= "zs_TE_DP_BVS_p03s.Rdata")
  
}


#====================================================================================================================================================
# RW1 
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


J = 8#dim(cparrots)[1]
T = 36#dim(cparrots)[2]

lambda = 100
sigma = 1
ncov1 = 7 # 5 continous and 2 categorical variables with 2 and 3 levels respectively
betas <- c(0.75, 1.25, 0.2, 2,0,-0.6, 0.5,-1, 0)# p approx 0.6

nsim = 50
#load("~/Objective priors for N-mixture models/Parrots/RW1_seed.Rdata")
load("~/Nmix/RW1_seed.Rdata")
numbers = sdata[1:nsim]


betas_sum = array(0, dim = c(length(betas), 5, nsim))
Ns_sum = array(0, dim = c(T, 5, nsim))
sigma_sum = matrix(0, nsim, 5)
lambda_sum = matrix(0, nsim, 5)
# effective sample size
ess = matrix(NA, nsim, 84)

for (n in 1:nsim) {
  
  print(n)
  set.seed(numbers[n])
  data =  simRW1(J,T, sigma, betas, ncov1)
  mean(data$p)
  #plot(density(data$theta))
  # data$theta
  # plot(data$theta)
  plot(data$N, type = "b")
  # 
  
  numVars = length(betas)#ncov+1 #number of covariates (including the intercept)
  X = data$X#[,c(1:4,7,8)]
  
  ncov <- numVars-1 # number of explanatory variables in the model excluding the intercept
  column_covariate = 1:ncov
  
  # true or false if the covaraites is categorical or numerical - I only have numerical covariates to begin with
  classCovariates <- rep(T, ncov)
  #classCovariates[match(fac_covariates, column_covariate)] <- F
  
  
  indexes_covariates <- c()
  indexes_covariates[1] <- 1
  
  if(any(column_covariate != 0)){
    
    indexes_covariates <- c()
    indexes_covariates[1] <- 1
    k <- 2
    for (i in 1:ncov) {
      if(classCovariates[i]){
        indexes_covariates[k] <- i + 1
        k <- k + 1
      } else {
        num_levels <- length(unique(X[,i]))
        indexes_covariates[k + 0:(num_levels-2)] <- i + 1
        k <- k + num_levels - 1
      }
    }
    
  }
  
  indexes_covariates # if only numerical covariates this is just 1:(ncov+1)
  # I also think tha the 1 is for the intercept
  
  # compute C matrix given only numerical covariates 
  C <- matrix(0, nrow = ncol(X) - 1, ncol = ncol(X) - 1)
  l <- 0 # l loop through the covariates 
  for (i in 1:ncov) {
    if(classCovariates[i]){ # if it's a numerical covariates
      C[l + 1, l + 1] <- 1
      l <- l + 1
    } else { # if it's a categorical covariates
      num_levels <- length(unique(originalX[,i]))
      C[l + 1:(num_levels-1), l + 1:(num_levels-1)] <- .5 * (diag(1, nrow = (num_levels-1)) + 
                                                               matrix(1, nrow = (num_levels-1), ncol = (num_levels-1)))
      l <- l + num_levels - 1
    }
  }
  
  C
  Sigma <- cov(X[,-1, drop = F])
  sumCSigma <- sum(C * Sigma)
  
  d_bar <- ncov      
  d_bar <- ifelse(d_bar <= ncov, d_bar, 2)
  
  mu_psi = 0#0.5
  phi_mu =4
  phi_beta = 1/4
  
  # Prior parameters for beta, beta ~ MVN(b, B)   
  b_psi <- c(mu_psi, rep(0, ncol(X) - 1))
  sigma_beta_psi <- C * phi_beta  
  b_psi #b
  
  B_psi <- matrix(0, nrow = ncol(X), ncol = ncol(X))
  B_psi[1, 1] <- phi_mu
  B_psi[2:(ncol(X)), 2:(ncol(X))] <- sigma_beta_psi
  B_psi # B
  
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
  
  win.data <- list(S= S, tau_b0 = 3600, Nindex = Nindex,  nmonths = T, b = b_psi, B = B_psi, numVars = numVars, L = length(CarArg$adj), adj = CarArg$adj, weights= CarArg$weights, num = CarArg$num)
  str(win.data)
  
  if(n==1){
    
    
    nimTE <- nimbleCode({
      
      # Priors
      lambda ~ dgamma(0.01,0.01)
      sigma ~ dunif(0,15)
      tau <- 1/sigma^2
      
      beta[1:numVars] ~ dmnorm(b[1:numVars], B[1:numVars, 1:numVars]) 
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
        #  re[s] ~ dnorm(0, sigma.re)
        logit(p[s]) <- linpred[s] #+ re[s]
      }
      
    })
    
    
    Nst <- apply(data$Y, 2, max, na.rm = TRUE)+100 # Inits for latent N
    Mst <- max(Nst)
    TEmodel = nimbleModel(nimTE, constants = win.data, data = list(y =data$y, X=X), inits = list(sigma.re=1, M = Mst,  N= Nst, lambda =  rgamma(1,10,0.1), beta = rep(0, numVars), s = rnorm(T, 0,1), sigma = runif(1) )  )
    TEmodel$calculate()
    
    cTEmodel <- compileNimble(TEmodel)
    TEconf = configureMCMC(TEmodel, enableWAIC = TRUE)
    
    TEconf$monitors
    TEconf$resetMonitors()
    TEconf$addMonitors(c("lambda", "theta", "beta", "M", "N", "sigma" ))
    #TEconf$addMonitors(c("lambda", "sdBeta", "sigma.re", "theta", "beta", "w", "M", "N", "sigma" ))
    
    TEmcmc = buildMCMC(TEconf)
    TEmod = compileNimble(TEmcmc, project = TEmodel, resetFunctions = TRUE)
    #TEmcmc.out <- runMCMC(TEmod, niter = 5000, nburnin = 1000, samplesAsCodaMCMC = TRUE, nchains = 1, summary = TRUE) 
    #TEmcmc.out <- runMCMC(TEmod, niter = 150000, nburnin = 50000, thin = 20, samplesAsCodaMCMC = TRUE, nchains = 2, WAIC = TRUE, summary = TRUE) 
    TEmcmc.out <- runMCMC(TEmod, niter = 55000, nburnin = 10000, thin = 20, samplesAsCodaMCMC = TRUE, nchains = 2, WAIC = TRUE, summary = TRUE) 
    # TEmcmc.out$summary$all.chains#[1:52,]
    betas_sum[,,n] = TEmcmc.out$summary$all.chains[(T+2): (T+(length(betas))+1 ), ]
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
    
    #TEmcmc.out <- runMCMC(TEmod, niter = 150000, nburnin = 50000, thin = 20, samplesAsCodaMCMC = TRUE, nchains = 2, summary = TRUE) 
    TEmcmc.out <- runMCMC(TEmod, niter = 55000, nburnin = 10000, thin = 20, samplesAsCodaMCMC = TRUE, nchains = 2, WAIC = TRUE, summary = TRUE) 
    #TEmcmc.out$summary$all.chains#[1:52,]
    betas_sum[,,n] = TEmcmc.out$summary$all.chains[(T+2): (T+(length(betas))+1 ), ]
    lambda_sum[n,] = TEmcmc.out$summary$all.chains["lambda",]
    Ns_sum[,,n] =  TEmcmc.out$summary$all.chains[2:(T+1),]
    sigma_sum[n,] = TEmcmc.out$summary$all.chains["sigma",]
    samples = TEmcmc.out$samples
    ess[n,] = effectiveSize(samples)
    
    
  }
  
  save(betas_sum, file= "betas_TE_RW1_p06s.Rdata")
  save(Ns_sum, file= "Ns_TE_RW1_p06s.Rdata")
  save(sigma_sum, file= "sigma_TE_RW1_p06s.Rdata")
  save(ess, file = "ess_TE_RW1_p06s.Rdata")
  save(lambda_sum, file= "lambda_TE_RW1_p06s.Rdata")
  
  
}



J = 8#dim(cparrots)[1]
T = 36#dim(cparrots)[2]

lambda = 100
ncov1 = 7
betas <- c(-1.5, 1.25, 0.2, 2,0,-0.6, 0.5,-1, 0)# p approx 0.3
sigma = 1
nsim = 50
#load("~/Objective priors for N-mixture models/Parrots/RW1_seed.Rdata")
load("~/Nmix/RW1_seed.Rdata")
numbers = sdata[1:nsim]


betas_sum = array(0, dim = c(length(betas), 5, nsim))
Ns_sum = array(0, dim = c(T, 5, nsim))
sigma_sum = matrix(0, nsim, 5)
lambda_sum = matrix(0, nsim, 5)
# effective sample size
ess = matrix(NA, nsim, 84)

for (n in 1:nsim) {
  
  print(n)
  set.seed(numbers[n])
  data = simRW1(J,T, sigma, betas, ncov1)
  mean(data$p)
  #plot(density(data$theta))
  # data$theta
  # plot(data$theta)
  plot(data$N, type = "b")
  # 
  
  numVars = length(betas)#ncov+1 #number of covariates (including the intercept)
  X = data$X#[,c(1:4,7,8)]
  
  ncov <- numVars-1 # number of explanatory variables in the model excluding the intercept
  column_covariate = 1:ncov
  
  # true or false if the covaraites is categorical or numerical - I only have numerical covariates to begin with
  classCovariates <- rep(T, ncov)
  #classCovariates[match(fac_covariates, column_covariate)] <- F
  
  
  indexes_covariates <- c()
  indexes_covariates[1] <- 1
  
  if(any(column_covariate != 0)){
    
    indexes_covariates <- c()
    indexes_covariates[1] <- 1
    k <- 2
    for (i in 1:ncov) {
      if(classCovariates[i]){
        indexes_covariates[k] <- i + 1
        k <- k + 1
      } else {
        num_levels <- length(unique(X[,i]))
        indexes_covariates[k + 0:(num_levels-2)] <- i + 1
        k <- k + num_levels - 1
      }
    }
    
  }
  
  indexes_covariates # if only numerical covariates this is just 1:(ncov+1)
  # I also think tha the 1 is for the intercept
  
  # compute C matrix given only numerical covariates 
  C <- matrix(0, nrow = ncol(X) - 1, ncol = ncol(X) - 1)
  l <- 0 # l loop through the covariates 
  for (i in 1:ncov) {
    if(classCovariates[i]){ # if it's a numerical covariates
      C[l + 1, l + 1] <- 1
      l <- l + 1
    } else { # if it's a categorical covariates
      num_levels <- length(unique(originalX[,i]))
      C[l + 1:(num_levels-1), l + 1:(num_levels-1)] <- .5 * (diag(1, nrow = (num_levels-1)) + 
                                                               matrix(1, nrow = (num_levels-1), ncol = (num_levels-1)))
      l <- l + num_levels - 1
    }
  }
  
  C
  Sigma <- cov(X[,-1, drop = F])
  sumCSigma <- sum(C * Sigma)
  
  d_bar <- ncov      
  d_bar <- ifelse(d_bar <= ncov, d_bar, 2)
  
  mu_psi = 0#0.5
  phi_mu =4
  phi_beta = 1/4
  
  # Prior parameters for beta, beta ~ MVN(b, B)   
  b_psi <- c(mu_psi, rep(0, ncol(X) - 1))
  sigma_beta_psi <- C * phi_beta  
  b_psi #b
  
  B_psi <- matrix(0, nrow = ncol(X), ncol = ncol(X))
  B_psi[1, 1] <- phi_mu
  B_psi[2:(ncol(X)), 2:(ncol(X))] <- sigma_beta_psi
  B_psi # B
  
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
  
  win.data <- list(S= S, tau_b0 = 3600, Nindex = Nindex,  nmonths = T, b = b_psi, B = B_psi, numVars = numVars, L = length(CarArg$adj), adj = CarArg$adj, weights= CarArg$weights, num = CarArg$num)
  str(win.data)
  
  if(n==1){
    
    
    nimTE <- nimbleCode({
      
      # Priors
      lambda ~ dgamma(0.01,0.01)
      sigma ~ dunif(0,15)
      tau <- 1/sigma^2
      
      beta[1:numVars] ~ dmnorm(b[1:numVars], B[1:numVars, 1:numVars]) 
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
        #  re[s] ~ dnorm(0, sigma.re)
        logit(p[s]) <- linpred[s] #+ re[s]
      }
      
    })
    
    
    Nst <- apply(data$Y, 2, max, na.rm = TRUE)+100 # Inits for latent N
    Mst <- max(Nst)
    TEmodel = nimbleModel(nimTE, constants = win.data, data = list(y =data$y, X=X), inits = list(sigma.re=1, M = Mst,  N= Nst, lambda =  rgamma(1,10,0.1), beta = rep(0, numVars), s = rnorm(T, 0,1), sigma = runif(1) )  )
    TEmodel$calculate()
    
    cTEmodel <- compileNimble(TEmodel)
    TEconf = configureMCMC(TEmodel, enableWAIC = TRUE)
    
    TEconf$monitors
    TEconf$resetMonitors()
    TEconf$addMonitors(c("lambda", "theta", "beta", "M", "N", "sigma" ))
    #TEconf$addMonitors(c("lambda", "sdBeta", "sigma.re", "theta", "beta", "w", "M", "N", "sigma" ))
    
    TEmcmc = buildMCMC(TEconf)
    TEmod = compileNimble(TEmcmc, project = TEmodel, resetFunctions = TRUE)
    #TEmcmc.out <- runMCMC(TEmod, niter = 5000, nburnin = 1000, samplesAsCodaMCMC = TRUE, nchains = 1, summary = TRUE) 
    #TEmcmc.out <- runMCMC(TEmod, niter = 150000, nburnin = 50000, thin = 20, samplesAsCodaMCMC = TRUE, nchains = 2, WAIC = TRUE, summary = TRUE) 
    TEmcmc.out <- runMCMC(TEmod, niter = 55000, nburnin = 10000, thin = 20, samplesAsCodaMCMC = TRUE, nchains = 2, WAIC = TRUE, summary = TRUE) 
    # TEmcmc.out$summary$all.chains#[1:52,]
    betas_sum[,,n] = TEmcmc.out$summary$all.chains[(T+2): (T+(length(betas))+1 ), ]
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
    
    #TEmcmc.out <- runMCMC(TEmod, niter = 150000, nburnin = 50000, thin = 20, samplesAsCodaMCMC = TRUE, nchains = 2, summary = TRUE) 
    TEmcmc.out <- runMCMC(TEmod, niter = 55000, nburnin = 10000, thin = 20, samplesAsCodaMCMC = TRUE, nchains = 2, WAIC = TRUE, summary = TRUE) 
    #TEmcmc.out$summary$all.chains#[1:52,]
    betas_sum[,,n] = TEmcmc.out$summary$all.chains[(T+2): (T+(length(betas))+1 ), ]
    lambda_sum[n,] = TEmcmc.out$summary$all.chains["lambda",]
    Ns_sum[,,n] =  TEmcmc.out$summary$all.chains[2:(T+1),]
    sigma_sum[n,] = TEmcmc.out$summary$all.chains["sigma",]
    samples = TEmcmc.out$samples
    ess[n,] = effectiveSize(samples)
    
    
  }
  
  save(betas_sum, file= "betas_TE_RW1_p03s.Rdata")
  save(Ns_sum, file= "Ns_TE_RW1_p03s.Rdata")
  save(sigma_sum, file= "sigma_TE_RW1_p03s.Rdata")
  save(ess, file = "ess_TE_RW1_p03s.Rdata")
  save(lambda_sum, file= "lambda_TE_RW1_p03s.Rdata")
  
  
}

#===============================================================================================================================================


#=======================================================================================================================================================
# Results 
#=====================================================================================================================================================

{
  library(nimble)
  library(coda)
  library(pgdraw)
  library(data.table)
  #source("TE_BVS_samplers.R")
  library(mcmcplots)
  library(MCMCvis)
  library(ggmcmc)
  library(truncnorm)
  
}

#=====================================================================================================================================================
# DP
#=====================================================================================================================================================

# p approx 0.6

simDP <- function(nsurveys, nmonths, beta, ncov){
  
  
  M = rpois(1, lambda)
  
  n = nmonths/2
  S = nsurveys*nmonths
  x = rbeta(n, 10, 10);mean(x); var(x)
  z = rbeta(n, 10, 1); mean(z); var(z)
  #theta <- list(y = c(x, z)) # to clusters
  theta <- list(y =  sample(c(x, z)) ) # to clusters
  
  N <- rbinom(nmonths, M, theta$y)
  
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


J = nsurveys=  8#dim(cparrots)[1]
T = nmonths = 36#dim(cparrots)[2]
lambda = 100
#betas <- c(0.5, 0.3, 0.1, 0.5, 0.1, 1, 0, 0) # p approx 0.6
#true_betas <- c(-1.5, -0.8, -0.6, -0.1, -0.8, -0.1, 0, 0)# p approx 0.3
ncov1 = 7 # 5 continous and 2 categorical variables with 2 and 3 levels respectively
betas <- c(0.75, 1.25, 0.2, 2,0,-0.6, 0.5,-1, 0)# p approx 0.6
nsim = 50
numbers = 1:nsim

sumSim <- function(samples, true_value, nsim, gamma = F ){
  lam_sum = samples
  lambda = true_value
  if (gamma == TRUE){
    mean1 = mean(lam_sum[,1]) # mean
    return(mean1)
  } else{
    med = median(lam_sum[,2]) # mean
    cov = mean(lam_sum[,4] < lambda & lam_sum[,5]>lambda) # coverage
    rmse= sqrt(mean((lam_sum[,2]-lambda)^2))/lambda
    
    lam = lam_sum[,2]
    N= nsim
    lam_bias = numeric(N)
    
    for (i in 1:N) {
      lam_bias[i] = (lam[i] - lambda)/lambda
    }
    
    bias = median(lam_bias)
    out = list(med = med, cov = cov, rmse = rmse, bias = bias)
    return(out)
    
  }
  
}

#lam_sum=t(betas_sum[1,,])
#lambda = betas[1]

load("~/Objective priors for N-mixture models/Parrots/sc1/betas_TE_DP_p06s.Rdata")
sumSim(t(betas_sum[1,,]), betas[1], nsim)
sumSim(t(betas_sum[2,,]), betas[2], nsim)
sumSim(t(betas_sum[3,,]), betas[3], nsim)

sumSim(t(betas_sum[4,,]), betas[4], nsim)
sumSim(t(betas_sum[5,,]), betas[5], nsim)
sumSim(t(betas_sum[6,,]), betas[6], nsim)

sumSim(t(betas_sum[7,,]), betas[7], nsim)
sumSim(t(betas_sum[8,,]), betas[8], nsim)
sumSim(t(betas_sum[9,,]), betas[9], nsim)
betas#[c(1:4, 7:8)]


b1 = sumSim(t(betas_sum[1,,]), betas[1], nsim)$bias
b2 = sumSim(t(betas_sum[2,,]), betas[2], nsim)$bias
b3 = sumSim(t(betas_sum[3,,]), betas[3], nsim)$bias
b4 = sumSim(t(betas_sum[4,,]), betas[4], nsim)$bias
b5 = sumSim(t(betas_sum[5,,]), betas[5], nsim)$bias
b6 = sumSim(t(betas_sum[6,,]), betas[6], nsim)$bias
b7 = sumSim(t(betas_sum[7,,]), betas[7], nsim)$bias
b8 = sumSim(t(betas_sum[8,,]), betas[8], nsim)$bias

mean(b2, b3, b4,b5,b6, b7, b8) # mean bias for all coefficients
median(b2, b3, b4,b5,b6, b7, b8) # mean bias for all coefficients

b1 = sumSim(t(betas_sum[1,,]), betas[1], nsim)$cov
b2 = sumSim(t(betas_sum[2,,]), betas[2], nsim)$cov
b3 = sumSim(t(betas_sum[3,,]), betas[3], nsim)$cov
b4 = sumSim(t(betas_sum[4,,]), betas[4], nsim)$cov
b5 = sumSim(t(betas_sum[5,,]), betas[5], nsim)$cov
b6 = sumSim(t(betas_sum[6,,]), betas[6], nsim)$cov
b7 = sumSim(t(betas_sum[7,,]), betas[7], nsim)$cov
b8 = sumSim(t(betas_sum[8,,]), betas[8], nsim)$cov

mean(b2, b3, b4,b5,b6, b7, b8) # mean bias for all coefficients
median(b2, b3, b4,b5,b6, b7, b8) # mean bias for all coefficients


b1 = sumSim(t(betas_sum[1,,]), betas[1], nsim)$rmse
b2 = sumSim(t(betas_sum[2,,]), betas[2], nsim)$rmse
b3 = sumSim(t(betas_sum[3,,]), betas[3], nsim)$rmse
b4 = sumSim(t(betas_sum[4,,]), betas[4], nsim)$rmse
b5 = sumSim(t(betas_sum[5,,]), betas[5], nsim)$rmse
b6 = sumSim(t(betas_sum[6,,]), betas[6], nsim)$rmse
b7 = sumSim(t(betas_sum[7,,]), betas[7], nsim)$rmse
b8 = sumSim(t(betas_sum[8,,]), betas[8], nsim)$rmse

mean(b2, b3, b4,b5,b6, b7, b8) # mean bias for all coefficients
median(b2, b3, b4,b5,b6, b7, b8) # mean bias for all coefficients



# lambda
load("~/Objective priors for N-mixture models/Parrots/sc1/lambda_TE_DP_p06s.Rdata")
sumSim(lambda_sum, lambda, nsim)


load("~/Objective priors for N-mixture models/Parrots/sc1/zs_TE_DP_BVS_p06s.Rdata")
zs_sum
length(which( rowSums(zs_sum) == 36))
length(which(zs_sum[,1] < 36))
mean(zs_sum[,1])
mean(zs_sum[,2])

misclass = numeric(nsim)
for (i in 1:nsim) {
  misclass[i] =  (abs(zs_sum[i,1]- 18) + abs(zs_sum[i,2]-18))/36
}

misclass
mean(misclass)
median(misclass)
plot(density(misclass))
abline(v=median(misclass), col= "red" )
abline(v=mean(misclass), col= "blue" )


# best to use median not mean

# Ns
load("~/Objective priors for N-mixture models/Parrots/sc1/Ns_TE_DP_p06s.Rdata")
RW1_N = Ns_sum

# Bias 
RW1N_bias = matrix(0, nsim, T)

#RMSE
RW1N_rmse = matrix(0, nsim, T)

# Coverage 
RW1_cov = matrix(0, nsim,T)

for (i in 1:nsim) {
  
  set.seed(numbers[i])
  data = simDP(J, T,  betas, ncov1)
  N = data$N
  RW1N_rmse[i,] = sqrt(mean((RW1_N[, 2,i]-N)^2))/N
  for (t in 1:T) {
    RW1N_bias[i,t] = (RW1_N[t,2,i] - N[t])/N[t]
    if( RW1_N[t,4, i] <= N[t] & RW1_N[t,5, i]>=N[t]) RW1_cov[i,t] = 1
    else RW1_cov[i,t] = 0
    
  }
  
  
}


# Absolute bias
y = abs( colMeans(RW1N_bias))
plot(y, col = "blue", pch = 20, ylim = c(min(y),max(y) ) )

#  Bias 
y =  colMeans(RW1N_bias)
plot(y, col = "black", pch = 20, ylim = c(min(y),max(y)), ylab= "Relative bias", xlab = "Sampling occasions" )
abline(h = 0)
mean(y) # Mean bias for all N
median(y)


#  RMSE
y =  colMeans(RW1N_rmse)
plot(y, col = "black", pch = 20, ylim = c(min(y),max(y)), ylab= "Relative bias", xlab = "Sampling occasions" )
abline(h = 0)
mean(y) # Mean bias for all N
median(y)


# Coverage
b=  colMeans(RW1_cov)
plot(b, col = "black", pch = 20, ylim = c(min(b),max(b)), ylab= "Coverage", xlab = "Months" )
mean(b) # Mean coverage for all N
median(b)
#=====================================================================================================================================================

# p approx 0.3
simDP <- function(nsurveys, nmonths, beta, ncov){
  
  
  M = rpois(1, lambda)
  
  n = nmonths/2
  S = nsurveys*nmonths
  x = rbeta(n, 10, 10);mean(x); var(x)
  z = rbeta(n, 10, 1); mean(z); var(z)
  #theta <- list(y = c(x, z)) # to clusters
  theta <- list(y =  sample(c(x, z)) ) # to clusters
  
  N <- rbinom(nmonths, M, theta$y)
  
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


J = nsurveys=  8#dim(cparrots)[1]
T = nmonths = 36#dim(cparrots)[2]
lambda = 100


ncov1 = 7
betas <- c(-1.5, 1.25, 0.2, 2,0,-0.6, 0.5,-1, 0)# p approx 0.3
nsim = 50
numbers = 1:nsim

sumSim <- function(samples, true_value, nsim, gamma = F ){
  lam_sum = samples
  lambda = true_value
  if (gamma == TRUE){
    mean1 = mean(lam_sum[,1]) # mean
    return(mean1)
  } else{
    med = median(lam_sum[,2]) # mean
    cov = mean(lam_sum[,4] < lambda & lam_sum[,5]>lambda) # coverage
    rmse= sqrt(mean((lam_sum[,2]-lambda)^2))/lambda
    
    lam = lam_sum[,2]
    N= nsim
    lam_bias = numeric(N)
    
    for (i in 1:N) {
      lam_bias[i] = (lam[i] - lambda)/lambda
    }
    
    bias = median(lam_bias)
    out = list(med = med, cov = cov, rmse = rmse, bias = bias)
    return(out)
    
  }
  
}

load("~/Objective priors for N-mixture models/Parrots/sc1/betas_TE_DP_p03s.Rdata")
sumSim(t(betas_sum[1,,]), betas[1], nsim)
sumSim(t(betas_sum[2,,]), betas[2], nsim)
sumSim(t(betas_sum[3,,]), betas[3], nsim)

sumSim(t(betas_sum[4,,]), betas[4], nsim)
sumSim(t(betas_sum[5,,]), betas[5], nsim)
sumSim(t(betas_sum[6,,]), betas[6], nsim)

sumSim(t(betas_sum[7,,]), betas[7], nsim)
sumSim(t(betas_sum[8,,]), betas[8], nsim)
sumSim(t(betas_sum[9,,]), betas[9], nsim)
betas#[c(1:4, 7:8)]

b1 = sumSim(t(betas_sum[1,,]), betas[1], nsim)$bias
b2 = sumSim(t(betas_sum[2,,]), betas[2], nsim)$bias
b3 = sumSim(t(betas_sum[3,,]), betas[3], nsim)$bias
b4 = sumSim(t(betas_sum[4,,]), betas[4], nsim)$bias
b5 = sumSim(t(betas_sum[5,,]), betas[5], nsim)$bias
b6 = sumSim(t(betas_sum[6,,]), betas[6], nsim)$bias
b7 = sumSim(t(betas_sum[7,,]), betas[7], nsim)$bias
b8 = sumSim(t(betas_sum[8,,]), betas[8], nsim)$bias


mean(b2, b3, b4,b5,b6, b7, b8) # mean bias for all coefficients
median(b2, b3, b4,b5,b6, b7, b8) # mean bias for all coefficients

b1 = sumSim(t(betas_sum[1,,]), betas[1], nsim)$cov
b2 = sumSim(t(betas_sum[2,,]), betas[2], nsim)$cov
b3 = sumSim(t(betas_sum[3,,]), betas[3], nsim)$cov
b4 = sumSim(t(betas_sum[4,,]), betas[4], nsim)$cov
b5 = sumSim(t(betas_sum[5,,]), betas[5], nsim)$cov
b6 = sumSim(t(betas_sum[6,,]), betas[6], nsim)$cov
b7 = sumSim(t(betas_sum[7,,]), betas[7], nsim)$cov
b8 = sumSim(t(betas_sum[8,,]), betas[8], nsim)$cov

mean(b2, b3, b4,b5,b6, b7, b8) # mean bias for all coefficients
median(b2, b3, b4,b5,b6, b7, b8) # mean bias for all coefficients


b1 = sumSim(t(betas_sum[1,,]), betas[1], nsim)$rmse
b2 = sumSim(t(betas_sum[2,,]), betas[2], nsim)$rmse
b3 = sumSim(t(betas_sum[3,,]), betas[3], nsim)$rmse
b4 = sumSim(t(betas_sum[4,,]), betas[4], nsim)$rmse
b5 = sumSim(t(betas_sum[5,,]), betas[5], nsim)$rmse
b6 = sumSim(t(betas_sum[6,,]), betas[6], nsim)$rmse
b7 = sumSim(t(betas_sum[7,,]), betas[7], nsim)$rmse
b8 = sumSim(t(betas_sum[8,,]), betas[8], nsim)$rmse

mean(b2, b3, b4,b5,b6, b7, b8) # mean bias for all coefficients
median(b2, b3, b4,b5,b6, b7, b8) # mean bias for all coefficients



# lambda
load("~/Objective priors for N-mixture models/Parrots/sc1/lambda_TE_DP_p03s.Rdata")
sumSim(lambda_sum, lambda, nsim)


load("~/Objective priors for N-mixture models/Parrots/sc1/zs_TE_DP_BVS_p03s.Rdata")
zs_sum
length(which( rowSums(zs_sum) == 36))
length(which(zs_sum[,1] < 36))
mean(zs_sum[,1])
mean(zs_sum[,2], na.rm = T)

misclass = numeric(nsim)
for (i in 1:nsim) {
  misclass[i] =  (abs(zs_sum[i,1]- 18) + abs(zs_sum[i,2]-18))/36
}

misclass
mean(misclass)
median(misclass)
plot(density(misclass))
abline(v=median(misclass), col= "red" )
abline(v=mean(misclass), col= "blue" )
# best to use median not mean


# Ns
load("~/Objective priors for N-mixture models/Parrots/sc1/Ns_TE_DP_p03s.Rdata")
RW1_N = Ns_sum

# Bias 
RW1N_bias = matrix(0, nsim, T)

#RMSE
RW1N_rmse = matrix(0, nsim, T)

# Coverage 
RW1_cov = matrix(0, nsim,T)

for (i in 1:nsim) {
  
  set.seed(numbers[i])
  data = simDP(J, T,  betas, ncov1)
  N = data$N
  RW1N_rmse[i,] = sqrt(mean((RW1_N[, 2,i]-N)^2))/N
  for (t in 1:T) {
    RW1N_bias[i,t] = (RW1_N[t,2,i] - N[t])/N[t]
    if( RW1_N[t,4, i] <= N[t] & RW1_N[t,5, i]>=N[t]) RW1_cov[i,t] = 1
    else RW1_cov[i,t] = 0
    
  }
  
  
}


# Absolute bias
y = abs( colMeans(RW1N_bias))
plot(y, col = "blue", pch = 20, ylim = c(min(y),max(y) ) )

#  Bias 
y =  colMeans(RW1N_bias)
plot(y, col = "black", pch = 20, ylim = c(min(y),max(y)), ylab= "Relative bias", xlab = "Sampling occasions" )
abline(h = 0)
mean(y) # Mean bias for all N
median(y)


#  RMSE
y =  colMeans(RW1N_rmse)
plot(y, col = "black", pch = 20, ylim = c(min(y),max(y)), ylab= "Relative bias", xlab = "Sampling occasions" )
abline(h = 0)
mean(y) # Mean bias for all N
median(y)


# Coverage
b=  colMeans(RW1_cov)
plot(b, col = "black", pch = 20, ylim = c(min(b),max(b)), ylab= "Coverage", xlab = "Months" )
mean(b) # Mean coverage for all N
median(b)
#=====================================================================================================================================================

#=====================================================================================================================================================
# RW1
#=====================================================================================================================================================

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


J = 8#dim(cparrots)[1]
T = 36#dim(cparrots)[2]

lambda = 100
sigma = 1
ncov1 = 7 # 5 continous and 2 categorical variables with 2 and 3 levels respectively
betas <- c(0.75, 1.25, 0.2, 2,0,-0.6, 0.5,-1, 0)# p approx 0.6
nsim = 50
load("~/Objective priors for N-mixture models/Parrots/RW1_seed.Rdata")
#load("~/Nmix/RW1_seed.Rdata")
numbers = sdata[1:nsim]

sumSim <- function(samples, true_value, nsim, gamma = F ){
  lam_sum = samples
  lambda = true_value
  if (gamma == TRUE){
    mean1 = mean(lam_sum[,1]) # mean
    return(mean1)
  } else{
    med = median(lam_sum[,2]) # mean
    cov = mean(lam_sum[,4] < lambda & lam_sum[,5]>lambda) # coverage
    rmse= sqrt(mean((lam_sum[,2]-lambda)^2))/lambda
    
    lam = lam_sum[,2]
    N= nsim
    lam_bias = numeric(N)
    
    for (i in 1:N) {
      lam_bias[i] = (lam[i] - lambda)/lambda
    }
    
    bias = median(lam_bias)
    out = list(med = med, cov = cov, rmse = rmse, bias = bias)
    return(out)
    
  }
  
}

load("~/Objective priors for N-mixture models/Parrots/sc1/betas_TE_RW1_p06s.Rdata")
sumSim(t(betas_sum[1,,]), betas[1], nsim)
sumSim(t(betas_sum[2,,]), betas[2], nsim)
sumSim(t(betas_sum[3,,]), betas[3], nsim)

sumSim(t(betas_sum[4,,]), betas[4], nsim)
sumSim(t(betas_sum[5,,]), betas[5], nsim)
sumSim(t(betas_sum[6,,]), betas[6], nsim)

sumSim(t(betas_sum[7,,]), betas[7], nsim)
sumSim(t(betas_sum[8,,]), betas[8], nsim)
sumSim(t(betas_sum[9,,]), betas[9], nsim)

betas#[c(1:4, 7:8)]

b1 = sumSim(t(betas_sum[1,,]), betas[1], nsim)$bias
b2 = sumSim(t(betas_sum[2,,]), betas[2], nsim)$bias
b3 = sumSim(t(betas_sum[3,,]), betas[3], nsim)$bias
b4 = sumSim(t(betas_sum[4,,]), betas[4], nsim)$bias
b5 = sumSim(t(betas_sum[5,,]), betas[5], nsim)$bias
b6 = sumSim(t(betas_sum[6,,]), betas[6], nsim)$bias
b7 = sumSim(t(betas_sum[7,,]), betas[7], nsim)$bias
b8 = sumSim(t(betas_sum[8,,]), betas[8], nsim)$bias

mean(b2, b3, b4,b5,b6, b7, b8) # mean bias for all coefficients
median(b2, b3, b4,b5,b6, b7, b8) # mean bias for all coefficients

b1 = sumSim(t(betas_sum[1,,]), betas[1], nsim)$cov
b2 = sumSim(t(betas_sum[2,,]), betas[2], nsim)$cov
b3 = sumSim(t(betas_sum[3,,]), betas[3], nsim)$cov
b4 = sumSim(t(betas_sum[4,,]), betas[4], nsim)$cov
b5 = sumSim(t(betas_sum[5,,]), betas[5], nsim)$cov
b6 = sumSim(t(betas_sum[6,,]), betas[6], nsim)$cov
b7 = sumSim(t(betas_sum[7,,]), betas[7], nsim)$cov
b8 = sumSim(t(betas_sum[8,,]), betas[8], nsim)$cov

mean(b2, b3, b4,b5,b6, b7, b8) # mean bias for all coefficients
median(b2, b3, b4,b5,b6, b7, b8) # mean bias for all coefficients

b1 = sumSim(t(betas_sum[1,,]), betas[1], nsim)$rmse
b2 = sumSim(t(betas_sum[2,,]), betas[2], nsim)$rmse
b3 = sumSim(t(betas_sum[3,,]), betas[3], nsim)$rmse
b4 = sumSim(t(betas_sum[4,,]), betas[4], nsim)$rmse
b5 = sumSim(t(betas_sum[5,,]), betas[5], nsim)$rmse
b6 = sumSim(t(betas_sum[6,,]), betas[6], nsim)$rmse
b7 = sumSim(t(betas_sum[7,,]), betas[7], nsim)$rmse
b8 = sumSim(t(betas_sum[8,,]), betas[8], nsim)$rmse

mean(b2, b3, b4,b5,b6, b7, b8) # mean bias for all coefficients
median(b2, b3, b4,b5,b6, b7, b8) # mean bias for all coefficients



# lambda
load("~/Objective priors for N-mixture models/Parrots/sc1/lambda_TE_RW1_p06s.Rdata")
sumSim(lambda_sum, lambda, nsim)

load("~/Objective priors for N-mixture models/Parrots/sc1/sigma_TE_RW1_p06s.Rdata")
sumSim(sigma_sum, sigma, nsim)

# Ns
load("~/Objective priors for N-mixture models/Parrots/sc1/Ns_TE_RW1_p06s.Rdata")
RW1_N = Ns_sum

# Bias 
RW1N_bias = matrix(0, nsim, T)

# RMSE
RW1N_rmse = matrix(0, nsim, T)


# Coverage 
RW1_cov = matrix(0, nsim,T)

for (i in 1:nsim) {
  
  set.seed(numbers[i])
  data = simRW1(J,T, sigma, betas, ncov1)
  N = data$N
  RW1N_rmse[i,] = sqrt(mean((RW1_N[, 2,i]-N)^2))/N
  for (t in 1:T) {
    RW1N_bias[i,t] = (RW1_N[t,2,i] - N[t])/N[t]
    if( RW1_N[t,4, i] <= N[t] & RW1_N[t,5, i]>=N[t]) RW1_cov[i,t] = 1
    else RW1_cov[i,t] = 0
    
  }
  
  
}

# Absolute bias
y = abs( colMeans(RW1N_bias))
plot(y, col = "blue", pch = 20, ylim = c(min(y),max(y) ) )

#  Bias 
y =  colMeans(RW1N_bias)
plot(y, col = "black", pch = 20, ylim = c(min(y),max(y)), ylab= "Relative bias", xlab = "Sampling occasions" )
abline(h = 0)
median(y)

#  RMSE 
y =  colMeans(RW1N_rmse)
plot(y, col = "black", pch = 20, ylim = c(min(y),max(y)), ylab= "Relative bias", xlab = "Sampling occasions" )
abline(h = 0)
median(y)


# Coverage
b=  colMeans(RW1_cov)
plot(b, col = "black", pch = 20, ylim = c(min(b),max(b)), ylab= "Coverage", xlab = "Months" )
median(b)
#=====================================================================================================================================================

# p approx 0.3


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


J = 8#dim(cparrots)[1]
T = 36#dim(cparrots)[2]

lambda = 100
sigma = 1
ncov1 = 7 # 5 continous and 2 categorical variables with 2 and 3 levels respectively
betas <- c(-1.5, 1.25, 0.2, 2,0,-0.6, 0.5,-1, 0)# p approx 0.3

nsim = 50
load("~/Objective priors for N-mixture models/Parrots/RW1_seed.Rdata")
#load("~/Nmix/RW1_seed.Rdata")
numbers = sdata[1:nsim]

sumSim <- function(samples, true_value, nsim, gamma = F ){
  lam_sum = samples
  lambda = true_value
  if (gamma == TRUE){
    mean1 = mean(lam_sum[,1]) # mean
    return(mean1)
  } else{
    med = median(lam_sum[,2]) # mean
    cov = mean(lam_sum[,4] < lambda & lam_sum[,5]>lambda) # coverage
    rmse= sqrt(mean((lam_sum[,2]-lambda)^2))/lambda
    
    lam = lam_sum[,2]
    N= nsim
    lam_bias = numeric(N)
    
    for (i in 1:N) {
      lam_bias[i] = (lam[i] - lambda)/lambda
    }
    
    bias = median(lam_bias)
    out = list(med = med, cov = cov, rmse = rmse, bias = bias)
    return(out)
    
  }
  
}

load("~/Objective priors for N-mixture models/Parrots/sc1/betas_TE_RW1_p03s.Rdata")
sumSim(t(betas_sum[1,,]), betas[1], nsim)
sumSim(t(betas_sum[2,,]), betas[2], nsim)
sumSim(t(betas_sum[3,,]), betas[3], nsim)

sumSim(t(betas_sum[4,,]), betas[4], nsim)
sumSim(t(betas_sum[5,,]), betas[5], nsim)
sumSim(t(betas_sum[6,,]), betas[6], nsim)

sumSim(t(betas_sum[7,,]), betas[7], nsim)
sumSim(t(betas_sum[8,,]), betas[8], nsim)
sumSim(t(betas_sum[9,,]), betas[9], nsim)
betas#[c(1:4, 7:8)]

b1 = sumSim(t(betas_sum[1,,]), betas[1], nsim)$bias
b2 = sumSim(t(betas_sum[2,,]), betas[2], nsim)$bias
b3 = sumSim(t(betas_sum[3,,]), betas[3], nsim)$bias
b4 = sumSim(t(betas_sum[4,,]), betas[4], nsim)$bias
b5 = sumSim(t(betas_sum[5,,]), betas[5], nsim)$bias
b6 = sumSim(t(betas_sum[6,,]), betas[6], nsim)$bias
b7 = sumSim(t(betas_sum[7,,]), betas[7], nsim)$bias
b8 = sumSim(t(betas_sum[8,,]), betas[8], nsim)$bias

mean(b2, b3, b4,b5,b6, b7, b8) # mean bias for all coefficients
median(b2, b3, b4,b5,b6, b7, b8) # mean bias for all coefficients

b1 = sumSim(t(betas_sum[1,,]), betas[1], nsim)$cov
b2 = sumSim(t(betas_sum[2,,]), betas[2], nsim)$cov
b3 = sumSim(t(betas_sum[3,,]), betas[3], nsim)$cov
b4 = sumSim(t(betas_sum[4,,]), betas[4], nsim)$cov
b5 = sumSim(t(betas_sum[5,,]), betas[5], nsim)$cov
b6 = sumSim(t(betas_sum[6,,]), betas[6], nsim)$cov
b7 = sumSim(t(betas_sum[7,,]), betas[7], nsim)$cov
b8 = sumSim(t(betas_sum[8,,]), betas[8], nsim)$cov

mean(b2, b3, b4,b5,b6, b7, b8) # mean bias for all coefficients
median(b2, b3, b4,b5,b6, b7, b8) # mean bias for all coefficients


b1 = sumSim(t(betas_sum[1,,]), betas[1], nsim)$rmse
b2 = sumSim(t(betas_sum[2,,]), betas[2], nsim)$rmse
b3 = sumSim(t(betas_sum[3,,]), betas[3], nsim)$rmse
b4 = sumSim(t(betas_sum[4,,]), betas[4], nsim)$rmse
b5 = sumSim(t(betas_sum[5,,]), betas[5], nsim)$rmse
b6 = sumSim(t(betas_sum[6,,]), betas[6], nsim)$rmse
b7 = sumSim(t(betas_sum[7,,]), betas[7], nsim)$rmse
b8 = sumSim(t(betas_sum[8,,]), betas[8], nsim)$rmse

mean(b2, b3, b4,b5,b6, b7, b8) # mean bias for all coefficients
median(b2, b3, b4,b5,b6, b7, b8) # mean bias for all coefficients

# lambda
load("~/Objective priors for N-mixture models/Parrots/sc1/lambda_TE_RW1_p03s.Rdata")
sumSim(lambda_sum, lambda, nsim)

load("~/Objective priors for N-mixture models/Parrots/sc1/sigma_TE_RW1_p03s.Rdata")
sumSim(sigma_sum, sigma, nsim)

# Ns
load("~/Objective priors for N-mixture models/Parrots/sc1/Ns_TE_RW1_p03s.Rdata")
RW1_N = Ns_sum

# Bias 
RW1N_bias = matrix(0, nsim, T)


# RMSE
RW1N_rmse = matrix(0, nsim, T)


# Coverage 
RW1_cov = matrix(0, nsim,T)

for (i in 1:nsim) {
  
  set.seed(numbers[i])
  data = simRW1(J,T, sigma, betas, ncov1)
  N = data$N
  RW1N_rmse[i,] = sqrt(mean((RW1_N[, 2,i]-N)^2))/N
  for (t in 1:T) {
    RW1N_bias[i,t] = (RW1_N[t,2,i] - N[t])/N[t]
    if( RW1_N[t,4, i] <= N[t] & RW1_N[t,5, i]>=N[t]) RW1_cov[i,t] = 1
    else RW1_cov[i,t] = 0
    
  }
  
  
}

# Absolute bias
y = abs( colMeans(RW1N_bias))
plot(y, col = "blue", pch = 20, ylim = c(min(y),max(y) ) )

#  Bias 
y =  colMeans(RW1N_bias)
plot(y, col = "black", pch = 20, ylim = c(min(y),max(y)), ylab= "Relative bias", xlab = "Sampling occasions" )
abline(h = 0)
median(y)

#  RMSE 
y =  colMeans(RW1N_rmse)
plot(y, col = "black", pch = 20, ylim = c(min(y),max(y)), ylab= "Relative bias", xlab = "Sampling occasions" )
abline(h = 0)
median(y)


# Coverage
b=  colMeans(RW1_cov)
plot(b, col = "black", pch = 20, ylim = c(min(b),max(b)), ylab= "Coverage", xlab = "Months" )
median(b)