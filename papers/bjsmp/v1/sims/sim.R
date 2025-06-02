source("funs.R")

runSim = function(x){
  set.seed(x$seed)
  I = x$I
  J = x$J
  L = x$L
  t.rho = x$t.rho
  t.tau = x$t.tau
  t.sig = x$t.sig
  
  out = simulateSNRdata(I, L, t.sig, t.rho, t.tau, J)
  sample.cor = out$sample_cor
  sample.ci.diff = out$sample_cor_ci_diff
  
  dat = out$dat
  parms = c("mu", "omega", "pTau2")
  data_list = list(
    "y" = dat$y,
    "N" = nrow(dat),
    "I" = length(unique(dat$sub)),
    "sub" = dat$sub,
    "task" = dat$task,
    "S" = rep(t.sig, J),
    "J" = length(unique(dat$task))
  )
  
  fit = runJags(mod, data_list, parms)
  
  pTau2 = fit$BUGSoutput$sims.list$pTau2
  pSig2 = fit$BUGSoutput$sims.list$omega
  
  tau2 = sapply(1:length(pTau2), function(i) 1/pTau2[i])
  
  sig2 = lapply(1:nrow(pSig2), function(i) solve(pSig2[i,,]))
  sig2 = abind(sig2, along = 0)
  
  cor.mat = lapply(1:nrow(pSig2), function(i) cov2cor(solve(pSig2[i,,])))
  cor.mat = abind(cor.mat, along = 0)
  summary.cor = apply(cor.mat, 2:3, function(x) quantile(x, c(.025, .5, .975)))
  ci_width = summary.cor[3,2,1] - summary.cor[1,2,1]
  
  
  out = list(
    tau2 = tau2,
    sig2 = sig2,
    I = I,
    J = J,
    L = L,
    t.rho = t.rho,
    SNR = t.sig/t.tau,
    Seed = x$seed,
    summary.cor = summary.cor,
    ci_width = unname(ci_width),
    sample.ci.diff = sample.ci.diff,
    sample.cor = sample.cor
  )
  
  return(out)
}


J = 2
t.tau = 1
L_vals = rep(c(50, 100, 200, 400), each = 1000)

design_grid = expand.grid(
  I = c(100, 200),
  t.rho = c(0, .5),
  t.sig = c(0.1, 0.2),
  L = L_vals
)

set.seed(123)  
design_grid$seed = sample(1e6, nrow(design_grid), replace = FALSE)

design_grid$J = J
design_grid$t.tau = t.tau

design_list = apply(design_grid, 1, function(row) {
  list(
    I = as.integer(row["I"]),
    J = as.integer(row["J"]),
    L = as.integer(row["L"]),
    t.rho = as.numeric(row["t.rho"]),
    t.sig = as.numeric(row["t.sig"]),
    t.tau = as.numeric(row["t.tau"]),
    seed = as.integer(row["seed"])
  )
})

results = mclapply(seq_along(design_list), function(i) {
  res = runSim(design_list[[i]])
  saveRDS(res, file = sprintf("output2/output_%03d.rds", i))
  cat(sprintf("Saved: output_%03d.rds\n", i))
  return(res)
}, mc.cores = 18)

saveRDS(design_grid, "design_grid_full.rds")
saveRDS(results, "results_raw_full.rds")

