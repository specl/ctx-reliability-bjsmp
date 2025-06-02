library(parallel)
library(mvtnorm)
library(R2jags)
library(abind)


mod = "
model{
  pTau2 ~ dgamma(.5, .5)
  for (n in 1:N){
    y[n] ~ dnorm(theta[sub[n],task[n]], pTau2)}
  for (i in 1:I){
      theta[i,1:J] ~ dmnorm(mu, omega)}
  omega~dscaled.wishart(S,2)
  for (j in 1:J) {
    mu[j] ~ dnorm(0, 1/100)
  }
} 
"


runJags = function(mod, dat, pars, nchains = 2, niter = 2000, nburnin = 500, inits = NULL) {
  if (is.null(inits)) {
    jags(model.file = textConnection(mod), data = dat, n.chains = nchains,
         n.iter = niter, n.burnin = nburnin, n.thin = 1, parameters.to.save = pars)
  } else {
    jags(model.file = textConnection(mod), data = dat, inits = inits,
         n.chains = nchains, n.iter = niter, n.burnin = nburnin, n.thin = 1,
         parameters.to.save = pars)
  }
}


makeTheta = function(I, J, t.Sig) {
  t.theta = mvtnorm::rmvnorm(I, mean = rep(0, J), sigma = t.Sig) 
  return(t.theta)
}


makeEmptyDat = function(I, J, L){
  sub = rep(1:I, each = J*L)
  task = rep(rep(1:J, each = L),I)
  trials = rep(rep(1:L, J), I)
  dat = data.frame("sub" = sub,
                   "task" = task,
                   "trial" = trials)
  return(dat)
}


addRT = function(dat, t.Theta, t.sig = .2){
  subtask = cbind(dat$sub, dat$task)
  mu = t.Theta[subtask]
  dat$y = rnorm(length(mu), mu, t.sig)
  return(dat)
}

simulateSNRdata = function(I, L, t.sig, t.rho, t.tau = 1, J = 2) {
  # Create correlation matrix
  t.cor = matrix(t.rho, nrow = J, ncol = J)
  diag(t.cor) = 1
  
  # Standard deviations for each task
  sigs = rep(t.sig, J)
  
  # Covariance matrix: scaled version of correlation matrix
  t.cov = diag(sigs) %*% t.cor %*% diag(sigs)
  
  # Simulate latent true abilities
  t.theta = makeTheta(I, J, t.cov)
  
  # Simulate observed data with noise
  dat = makeEmptyDat(I, J, L)
  dat = addRT(dat, t.theta, t.tau)
  
  # Compute observed and true correlations and Signal-to-noise ratio
  scores = tapply(dat$y, list(dat$sub, dat$task), mean)
  sample.cor = cor.test(scores[, 1], scores[, 2])
  samp.cor = sample.cor$estimate
  ci_dif = diff(sample.cor$conf.int)
  maifast.cor = cor(t.theta)
  SNR = t.sig/t.tau
  
  out = list(
    dat = dat,
    I = I,
    L = L,
    SNR = SNR,
    t.rho = t.rho,
    t.theta = t.theta,
    t.cov = t.cov,
    sample_cor = samp.cor,
    sample_cor_ci_diff = ci_dif,
    latent_cor = maifast.cor[1, 2]
  )
  
  return(out)
}

