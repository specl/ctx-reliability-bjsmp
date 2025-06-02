proxI = function(SNR, L, t.rho, target_width = 0.2) {
  SNR2 = SNR^2
  A = SNR2 / (SNR2 + 1 / L)
  
  lossFun = function(logI) {
    I = exp(logI)
    if (I <= 3) return(1e6)
    
    SE_z = 1 / sqrt(I - 3)
    z_center = atanh(t.rho)
    z_lower = z_center - 1.96 * SE_z
    z_upper = z_center + 1.96 * SE_z
    unattenuated_width = tanh(z_upper) - tanh(z_lower)
    
    attenuated_width = unattenuated_width / A
    return((attenuated_width - target_width)^2)
  }
  
  opt = optimize(lossFun, interval = log(c(10, 10000000)))
  return(exp(opt$minimum))
}


proxCIwidth = function(SNR, L, I, t.rho){
  
  SNR2 = SNR^2
  A = SNR2 / (SNR2 + 1/L)
  SE_z = 1 / (sqrt((I  - 3))) 
  
  z_center = atanh(t.rho)
  z_lower = z_center - 1.96 * SE_z
  z_upper = z_center + 1.96 * SE_z
  r_width = (tanh(z_upper) - tanh(z_lower))/A
  
  return(r_width)
}

readSims = function(){
  abc = list.files("../../dev/extensions2/output2/", pattern = "\\.rds$", full.names = TRUE)
  all_results = lapply(abc, readRDS)
  
  summary_df = do.call(rbind, lapply(all_results, function(x) {
    if (is.null(x) || is.null(x$summary.cor)) return(NULL)
    
    data.frame(
      I = x$I,
      L = x$L,
      SNR = x$SNR,
      ci_width = x$ci_width,
      samp = x$sample.ci.diff,
      t.rho = x$t.rho,
      seeds = x$Seed,
      samp = x$sample.ci.diff
    )
  }))
  
  write.csv(summary_df, "sims/sim_summary.csv")
}



approxIandLs = function(target_width = .2){
  Ls = seq(1, 100000, by = 0.5)
  SNR = c(.1, .2, .5, 1, 2)
  
  Is = createLineValues(Ls, SNR, target_width)
  out = list(
    "Is" = Is,
    "Ls" = Ls,
    "SNR" = SNR
  )
  saveRDS(out, file = "sims/approxIandLs.rds")
}



drawLines = function(Ls, Is, cols){
  dwai = mapply(function(x){
    lines(Ls, Is[[x]], col = cols[x], lwd = 2)
  }, 1:length(Is))
}

plotIvsL = function(){
  if (file.exists("sims/approxIandLs.rds")){
    out = readRDS("sims/approxIandLs.rds")
  } else {
    approxIandLs()
    out = readRDS("approxIandLs.rds")
  }
  Is = out$Is
  Ls = out$Ls
  SNR = out$SNR
  
  
  plot(NA, NA, xlim = c(1, 10^3), ylim = c(10, 1100000), type = "n",
       xlab = "Trials", 
       ylab = "Individuals",
       main = "",
       log = "xy", axes = FALSE,
       cex.main = 1.5, cex.lab = 1.5)
  
  axis(1, at = 10^(1:4), labels = expression(10^1, 10^2, 10^3, 10^4))
  axis(1, at = c(1:10, as.vector(outer(1:9, 10^(1:2)))), labels = FALSE, tcl = -0.3)
  
  axis(2, at = 10^(1:6), labels = expression(10^1, 10^2, 10^3, 10^4, 10^5, 10^6), las = 2)
  axis(2, at = as.vector(outer(1:9, 10^(1:5))), labels = FALSE, tcl = -0.3)
  cols = c("red3", "#548B54", "#3A5FCD", "#FF69B4", "black" )
  drawLines(Ls, Is, cols)
  grid()
  box()
  labs = lapply(SNR, function(g) bquote(gamma == .(g)))
  legend("topright", 
         legend = labs,
         col = cols, 
         lwd = 2, 
         pt.cex = 1.2, cex = 1.1, 
         bty = "n", 
         inset = 0.01)
}



plotCIwidth= function(){
  if(file.exists("sims/sim_summary.csv")){
    summary_df = read.csv("sims/sim_summary.csv")
  } else {
    summary_df = readSims()
  }
  
  
  ind = which(
    summary_df$I == 100 &
      summary_df$t.rho == 0
  )
  
  res_1 = summary_df[ind,]
  
  res_1$width_ratio = res_1$ci_width / res_1$samp
  res_1$med_ci = sapply(1:nrow(res_1), function(x) median(res_1$ci_width[x], res_1$samp[x]))
  res_1$compress = res_1$SNR^2 / (res_1$SNR^2 + 1 / res_1$L)
  res_1$expansion = 1 / res_1$compress
  res_1$proxCIwidth = proxCIwidth(
    res_1$SNR,
    res_1$L,
    res_1$I,
    res_1$t.rho
  )
  
  line_med = as.data.frame.table(tapply(res_1$width_ratio, 
                                        list(res_1$I, res_1$L, res_1$SNR), 
                                        median))
  
  line_low = as.data.frame.table(tapply(res_1$width_ratio, 
                                        list(res_1$I, res_1$L, res_1$SNR), 
                                        quantile, probs = 0.10))
  
  line_high = as.data.frame.table(tapply(res_1$width_ratio, 
                                         list(res_1$I, res_1$L, res_1$SNR), 
                                         quantile, probs = 0.90))
  
  line_med$expansion = 1 / (as.numeric(as.character(line_med$Var3))^2 / 
                              (as.numeric(as.character(line_med$Var3))^2 + 1 / as.numeric(as.character(line_med$Var2))))
  line_high$expansion = 1 / (as.numeric(as.character(line_med$Var3))^2 / 
                               (as.numeric(as.character(line_med$Var3))^2 + 1 / as.numeric(as.character(line_med$Var2))))
  line_low$expansion = 1 / (as.numeric(as.character(line_med$Var3))^2 / 
                              (as.numeric(as.character(line_med$Var3))^2 + 1 / as.numeric(as.character(line_med$Var2))))
  
  par(mfrow = c(1,2))
  plot(NA, NA, 
       xlim = c(0.9, 3.2), ylim = c(1, 5),
       xlab = "Expansion Factor", ylab = "Width Ratio",
       main = "",
       cex.lab = 1.3, cex.main = 1.4, cex.axis = 1.1, axes = F)
  axis(1)
  axis(2, las = 2)
  box()
  
  jit = rnorm(length(res_1$width_ratio), 0 , .005)
  
  grid(col = "lightgray", lty = "dotted")
  
  points(res_1$expansion + jit, res_1$width_ratio,
         col = rgb(0.5, 0.5, 0.5, 0.4), pch = 19, cex = 0.6)
  
  points(line_low$expansion, line_low$Freq, col = "deepskyblue3", pch = 19)
  points(line_high$expansion, line_high$Freq, col = "indianred2", pch = 19)
  points(line_med$expansion, line_med$Freq, col = "black", pch = 19)
  
  lines(lowess(line_low$expansion, line_low$Freq, f = 1), col = "deepskyblue3", lwd = 2)
  lines(lowess(line_med$expansion, line_med$Freq, f = 1), col = "black", lwd = 2)
  lines(lowess(line_high$expansion, line_high$Freq, f = 1), col = "indianred2", lwd = 2)
  
  abline(0, 1, lty = 2, col = "darkgray")
  
  legend("topleft", 
         legend = c("10th percentile", "Median", "90th percentile"),
         col = c("deepskyblue3", "black", "indianred2"), 
         lwd = 2, pch = 19, 
         pt.cex = 1.2, cex = 1.1, 
         bty = "n", 
         inset = 0.01)
  
}

createLineValues = function(Ls, SNR, target_width) {
  
  if (length(SNR) == length(target_width) && length(SNR) > 1) {
    stop("SNR and target_width cannot be the same length unless one is scalar.")
  }
  if (length(SNR) == 1) {
    Is = lapply(target_width, function(w) 
      sapply(Ls, function(L) proxI(SNR, L = L, t.rho = 0, target_width = w)))
    
  } else if (length(target_width) == 1) {
    Is = lapply(SNR, function(s) 
      sapply(Ls, function(L) proxI(s, L = L, t.rho = 0, target_width = target_width)))
    
  } else {
    stop("Either SNR or target_width must be of length 1.")
  }
  
  return(Is)
}
