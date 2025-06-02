require('R2jags')

mod = "
model {
  for (n in 1:N) {
    y[n] ~ dnorm(mu[n], tau2Inv[task[n]])
    mu[n] <- alpha[sub[n], task[n]] + tau[task[n]]*(cond[n]-1.5) * theta[sub[n], task[n]]
  }

  for (i in 1:I) {
    for (j in 1:J) {
      alpha[i,j] ~ dnorm(a1, 1/a2^2) 
      effect[i,j] = tau[j]*theta[i,j]
    }
  }

  for (j in 1:J) {
    tau2Inv[j] ~ dgamma(0.5, 0.5*t^2) 
    tau[j] = sqrt(1/tau2Inv[j])
  }

  # Multivariate normal prior for theta vectors
  for (i in 1:I) {
    theta[i,1:J] ~ dmnorm(nu[1:J], Sigma_inv[1:J,1:J])
  }

  # Priors for nu_j
  for (j in 1:J) {
    nu[j] ~ dnorm(v1, 1/v2^2) 
  }

  # Prior for inverse covariance matrix (Scaled Wishart)
  Sigma_inv[1:J,1:J] ~ dscaled.wishart(s,1)
}
"

runMod=function(dat,M=200){
  I=length(unique(dat$sub))
  J=length(unique(dat$task))
  setup=list(
    "y" = dat$y,
    "task" = dat$task,
    "sub" = dat$sub,
    "cond" = dat$cond,
    "I" = I,
    "J" = J,
    "N" = nrow(dat))
  prior=list(
    "a1" = 700,
    "a2" = 1000,
    "t" = 250,
    "v1" = 50,
    "v2" = 100,
    "s" = rep(30^2,J))
  pars = c("theta","Sigma_inv","effect")
  out=jags(data=c(setup,prior), 
           parameters=pars, 
           model.file = textConnection(mod), 
           n.chains=1,n.iter=M,n.burnin=M/10,n.thin=1)
  return(out)
}


