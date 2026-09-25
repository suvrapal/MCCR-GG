#  GENERALISED GAMMA MCCR 3 RISK SEER distant data analysis
library(numDeriv)
library(pracma)

upper_incgam_log <- function(x, a) {
  lgamma(a) + pgamma(x, shape = a, scale = 1, lower.tail = FALSE, log.p = TRUE)
}

lower_incgam_log <- function(x, a) {
  lgamma(a) +  pgamma(x, shape = a, scale = 1, lower.tail = TRUE, log.p = TRUE)
}

haz_fun <- function(x, q, lam, sig){
  
  if(abs(q)< 1e-4){
    
    f <- dlnorm(x,
                meanlog = -log(lam),
                sdlog   = sig)
    
    S <- plnorm(x,
                meanlog = -log(lam),
                sdlog   = sig,
                lower.tail = FALSE)
    
    return(f/S)
  }
  a <- q^(-2)
  
  z <- a * (lam*x)^(q/sig)
  
  
  log_num <- log(abs(q)) + a*log(a) + ((q^(-1))/sig)*log(lam*x) - a*(lam*x)^(q/sig)
  
    if(q > 0){
    
    log_denom <- log(sig) + log(x) + upper_incgam_log(z,a)
    
  }
    if(q < 0){
    
    log_denom <- log(sig) + log(x) + lower_incgam_log(z,a)
    
  }
  
  
  return(exp(log_num - log_denom))
}

cum_haz_fun <- function(x,q,lam,sig){
  
  if(abs(q)< 1e-4){
    
    S <- plnorm(x,
                meanlog=-log(lam),
                sdlog=sig,
                lower.tail=FALSE)
    
    return(-log(S))
  }
  a <- q^(-2)
  
  log_arg <- (-2*log(abs(q))) + ((q/sig)*log(lam*x))
  
  z <- exp(log_arg)
  
  if(q > 0){
    
    log_S <- upper_incgam_log(z,a) - lgamma(a)
    
  }
  
  if(q < 0){
    
    log_S <- lower_incgam_log(z,a) - lgamma(a)
    
  }
  
  return(-log_S)
  
}

log.lik.fn = function(param=c(cure,q1, lam1, sig1, q2, lam2, sig2, q3, lam3, sig3), y0, y1, y2,y3) {
  cure = pmin(pmax(param[1], 1e-6), 1-1e-6)
  q1   = param[2]; lam1 = param[3]; sig1 = param[4]
  q2   = param[5]; lam2 = param[6]; sig2 = param[7]
  q3   = param[8]; lam3 = param[9]; sig3 = param[10]
  
  h1.0 = haz_fun(y0, q1, lam1, sig1)
  h1.1 = haz_fun(y1, q1, lam1, sig1)
  h1.2 = haz_fun(y2, q1, lam1, sig1)
  h1.3 = haz_fun(y3, q1, lam1, sig1)
  
  h2.0 = haz_fun(y0, q2, lam2, sig2)
  h2.1 = haz_fun(y1, q2, lam2, sig2)
  h2.2 = haz_fun(y2, q2, lam2, sig2)
  h2.3 = haz_fun(y3, q2, lam2, sig2)
  
  h3.0 = haz_fun(y0, q3, lam3, sig3)
  h3.1 = haz_fun(y1, q3, lam3, sig3)
  h3.2 = haz_fun(y2, q3, lam3, sig3)
  h3.3 = haz_fun(y3, q3, lam3, sig3)
  
  H1.0 = cum_haz_fun(y0, q1, lam1, sig1)
  H1.1 = cum_haz_fun(y1, q1, lam1, sig1)
  H1.2 = cum_haz_fun(y2, q1, lam1, sig1)
  H1.3 = cum_haz_fun(y3, q1, lam1, sig1)
  
  H2.0 = cum_haz_fun(y0, q2, lam2, sig2)
  H2.1 = cum_haz_fun(y1, q2, lam2, sig2)
  H2.2 = cum_haz_fun(y2, q2, lam2, sig2)
  H2.3 = cum_haz_fun(y3, q2, lam2, sig2)
  
  H3.0 = cum_haz_fun(y0, q3, lam3, sig3)
  H3.1 = cum_haz_fun(y1, q3, lam3, sig3)
  H3.2 = cum_haz_fun(y2, q3, lam3, sig3)
  H3.3 = cum_haz_fun(y3, q3, lam3, sig3)
  
  l.obs = sum(log(1 - cure) + log(h1.1) - H1.1 - H2.1 - H3.1) +
    sum(log((cure*exp(-H2.0-H3.0)) + ((1 - cure)*exp(-H1.0 - H2.0 -H3.0)) )) +
    sum(log(h2.2)) +
    sum(log((cure*exp(-H2.2-H3.2)) + ((1 - cure)*exp(-H1.2 - H2.2 - H3.2)) )) +
    sum(log(h3.3)) +
    sum(log((cure*exp(-H2.3-H3.3)) + ((1 - cure)*exp(-H1.3 - H2.3 - H3.3)) ))
  
  return(l.obs)
}

log.lik.weibull = function(param_weibull, y0, y1, y2, y3) {
  cure = pmin(pmax(param_weibull[1], 1e-6), 1 - 1e-6)
  
  q1 = 1; lam1 = param_weibull[2]; sig1 = param_weibull[3]
  q2 = 1; lam2 = param_weibull[4]; sig2 = param_weibull[5]
  q3 = 1; lam3 = param_weibull[6]; sig3 = param_weibull[7]
  
  h1.0 = haz_fun(y0, q1, lam1, sig1)
  h1.1 = haz_fun(y1, q1, lam1, sig1)
  h1.2 = haz_fun(y2, q1, lam1, sig1)
  h1.3 = haz_fun(y3, q1, lam1, sig1)
  
  h2.0 = haz_fun(y0, q2, lam2, sig2)
  h2.1 = haz_fun(y1, q2, lam2, sig2)
  h2.2 = haz_fun(y2, q2, lam2, sig2)
  h2.3 = haz_fun(y3, q2, lam2, sig2)
  
  h3.0 = haz_fun(y0, q3, lam3, sig3)
  h3.1 = haz_fun(y1, q3, lam3, sig3)
  h3.2 = haz_fun(y2, q3, lam3, sig3)
  h3.3 = haz_fun(y3, q3, lam3, sig3)
  
  H1.0 = cum_haz_fun(y0, q1, lam1, sig1)
  H1.1 = cum_haz_fun(y1, q1, lam1, sig1)
  H1.2 = cum_haz_fun(y2, q1, lam1, sig1)
  H1.3 = cum_haz_fun(y3, q1, lam1, sig1)
  
  
  H2.0 = cum_haz_fun(y0, q2, lam2, sig2)
  H2.1 = cum_haz_fun(y1, q2, lam2, sig2)
  H2.2 = cum_haz_fun(y2, q2, lam2, sig2)
  H2.3 = cum_haz_fun(y3, q2, lam2, sig2)
  
  H3.0 = cum_haz_fun(y0, q3, lam3, sig3)
  H3.1 = cum_haz_fun(y1, q3, lam3, sig3)
  H3.2 = cum_haz_fun(y2, q3, lam3, sig3)
  H3.3 = cum_haz_fun(y3, q3, lam3, sig3)
  
  
  l.obs = sum(log(1-param_weibull[1])+log(h1.1)-H1.1-H2.1-H3.1)+sum(log((param_weibull[1]*exp(-H2.0-H3.0))+((1-param_weibull[1])*exp(-H1.0-H2.0-H3.0))))+
    sum(log(h2.2))+sum(log((param_weibull[1]*exp(-H2.2-H3.2))+((1-param_weibull[1])*exp(-H1.2-H2.2-H3.2))))+
    sum(log(h3.3))+sum(log((param_weibull[1]*exp(-H2.3-H3.3))+((1-param_weibull[1])*exp(-H1.3-H2.3-H3.3))))
  
  return(l.obs)
}

log.lik.gamma = function(param_gamma, y0, y1, y2, y3) {
  cure = pmin(pmax(param_gamma[1], 1e-6), 1 - 1e-6)
  lam1 = param_gamma[2]
  sig1 = param_gamma[3]
  q1   = sig1 
  lam2 = param_gamma[4]
  sig2 = param_gamma[5]
  q2   = sig2
  lam3 = param_gamma[6]
  sig3 = param_gamma[7]
  q3   = sig3
  
  h1.1 = haz_fun(y1, q1, lam1, sig1)
  h2.2 = haz_fun(y2, q2, lam2, sig2)
  h3.3 = haz_fun(y3, q3, lam3, sig3)
  
  H1.0 = cum_haz_fun(y0, q1, lam1, sig1)
  H1.1 = cum_haz_fun(y1, q1, lam1, sig1)
  H1.2 = cum_haz_fun(y2, q1, lam1, sig1)
  H1.3 = cum_haz_fun(y3, q1, lam1, sig1)
  
  
  H2.0 = cum_haz_fun(y0, q2, lam2, sig2)
  H2.1 = cum_haz_fun(y1, q2, lam2, sig2)
  H2.2 = cum_haz_fun(y2, q2, lam2, sig2)
  H2.3 = cum_haz_fun(y3, q2, lam2, sig2)
  
  H3.0 = cum_haz_fun(y0, q3, lam3, sig3)
  H3.1 = cum_haz_fun(y1, q3, lam3, sig3)
  H3.2 = cum_haz_fun(y2, q3, lam3, sig3)
  H3.3 = cum_haz_fun(y3, q3, lam3, sig3)
  
  
  l.obs = sum(log(1-param_gamma[1])+log(h1.1)-H1.1-H2.1-H3.1)+sum(log((param_gamma[1]*exp(-H2.0-H3.0))+((1-param_gamma[1])*exp(-H1.0-H2.0-H3.0))))+
    sum(log(h2.2))+sum(log((param_gamma[1]*exp(-H2.2-H3.2))+((1-param_gamma[1])*exp(-H1.2-H2.2-H3.2))))+
    sum(log(h3.3))+sum(log((param_gamma[1]*exp(-H2.3-H3.3))+((1-param_gamma[1])*exp(-H1.3-H2.3-H3.3))))
  
  return(l.obs)
}

log.lik.gg.fixed.q = function(param_gengamma,q, y0, y1, y2,y3) {
  cure = pmin(pmax(param_gengamma[1], 1e-6), 1 - 1e-6)
  
  q1 = q; lam1 =param_gengamma[2]; sig1 = param_gengamma[3]
  q2 = q; lam2 = param_gengamma[4]; sig2 = param_gengamma[5]
  q3 = q; lam3 = param_gengamma[6]; sig3 = param_gengamma[7]
  
  
  h1.1 = haz_fun(y1, q1, lam1, sig1)
  h2.2 = haz_fun(y2, q2, lam2, sig2)
  h3.3 = haz_fun(y3, q3, lam3, sig3)
  
  H1.0 = cum_haz_fun(y0, q1, lam1, sig1)
  H1.1 = cum_haz_fun(y1, q1, lam1, sig1)
  H1.2 = cum_haz_fun(y2, q1, lam1, sig1)
  H1.3 = cum_haz_fun(y3, q1, lam1, sig1)
  
  
  H2.0 = cum_haz_fun(y0, q2, lam2, sig2)
  H2.1 = cum_haz_fun(y1, q2, lam2, sig2)
  H2.2 = cum_haz_fun(y2, q2, lam2, sig2)
  H2.3 = cum_haz_fun(y3, q2, lam2, sig2)
  
  H3.0 = cum_haz_fun(y0, q3, lam3, sig3)
  H3.1 = cum_haz_fun(y1, q3, lam3, sig3)
  H3.2 = cum_haz_fun(y2, q3, lam3, sig3)
  H3.3 = cum_haz_fun(y3, q3, lam3, sig3)
  
  
  l.obs = sum(log(1-param_gengamma[1])+log(h1.1)-H1.1-H2.1-H3.1)+sum(log((param_gengamma[1]*exp(-H2.0-H3.0))+((1-param_gengamma[1])*exp(-H1.0-H2.0-H3.0))))+
    sum(log(h2.2))+sum(log((param_gengamma[1]*exp(-H2.2-H3.2))+((1-param_gengamma[1])*exp(-H1.2-H2.2-H3.2))))+
    sum(log(h3.3))+sum(log((param_gengamma[1]*exp(-H2.3-H3.3))+((1-param_gengamma[1])*exp(-H1.3-H2.3-H3.3))))
  return(l.obs)
}


log.lik.lognormal = function(param_lognormal, y0, y1, y2,y3) {
  
  cure = pmin(pmax(param_lognormal[1], 1e-6), 1 - 1e-6)
  
  q1 = 0; lam1 = param_lognormal[2]; sig1 = param_lognormal[3]
  q2 = 0; lam2 = param_lognormal[4]; sig2 = param_lognormal[5]
  q3 = 0; lam3 = param_lognormal[6]; sig3 = param_lognormal[7]
  
  h1.0 = haz_fun(y0, q1, lam1, sig1)
  h1.1 = haz_fun(y1, q1, lam1, sig1)
  h1.2 = haz_fun(y2, q1, lam1, sig1)
  h1.3 = haz_fun(y3, q1, lam1, sig1)
  
  h2.0 = haz_fun(y0, q2, lam2, sig2)
  h2.1 = haz_fun(y1, q2, lam2, sig2)
  h2.2 = haz_fun(y2, q2, lam2, sig2)
  h2.3 = haz_fun(y3, q2, lam2, sig2)
  
  h3.0 = haz_fun(y0, q3, lam3, sig3)
  h3.1 = haz_fun(y1, q3, lam3, sig3)
  h3.2 = haz_fun(y2, q3, lam3, sig3)
  h3.3 = haz_fun(y3, q3, lam3, sig3)
  
  H1.0 = cum_haz_fun(y0, q1, lam1, sig1)
  H1.1 = cum_haz_fun(y1, q1, lam1, sig1)
  H1.2 = cum_haz_fun(y2, q1, lam1, sig1)
  H1.3 = cum_haz_fun(y3, q1, lam1, sig1)
  
  
  H2.0 = cum_haz_fun(y0, q2, lam2, sig2)
  H2.1 = cum_haz_fun(y1, q2, lam2, sig2)
  H2.2 = cum_haz_fun(y2, q2, lam2, sig2)
  H2.3 = cum_haz_fun(y3, q2, lam2, sig2)
  
  H3.0 = cum_haz_fun(y0, q3, lam3, sig3)
  H3.1 = cum_haz_fun(y1, q3, lam3, sig3)
  H3.2 = cum_haz_fun(y2, q3, lam3, sig3)
  H3.3 = cum_haz_fun(y3, q3, lam3, sig3)
  
  
  l.obs = sum(log(1-param_lognormal[1])+log(h1.1)-H1.1-H2.1-H3.1)+sum(log((param_lognormal[1]*exp(-H2.0-H3.0))+((1-param_lognormal[1])*exp(-H1.0-H2.0-H3.0))))+
    sum(log(h2.2))+sum(log((param_lognormal[1]*exp(-H2.2-H3.2))+((1-param_lognormal[1])*exp(-H1.2-H2.2-H3.2))))+
    sum(log(h3.3))+sum(log((param_lognormal[1]*exp(-H2.3-H3.3))+((1-param_lognormal[1])*exp(-H1.3-H2.3-H3.3))))
  
  return(l.obs)
  
  
  
}

# EM Algorithm for 3 risks
EM.MCCR3_global_lognormal = function(data, tol, maxit, c, q_1, lambda1, sigma1, q_2, lambda2, sigma2, q_3, lambda3, sigma3) {
  
  data0 = data[data$d == 0, ]
  data1 = data[data$d == 1, ]
  data2 = data[data$d == 2, ]
  data3 = data[data$d == 3, ]
  
  y0 = data0$y; y1 = data1$y; y2 = data2$y; y3 = data3$y
  
  p.old = matrix(c(c, q_1, lambda1, sigma1, q_2, lambda2, sigma2, q_3, lambda3, sigma3), ncol = 1)
  p.new = matrix(0, ncol = 1, nrow = 10)
  
  continue = TRUE
  iter = 1
  
  global_lognormal_iter_count = 0
  
  while(continue) {
     all_lognormal_triggered = (abs(p.old[2,1]) < 1e-4) || (abs(p.old[5,1]) < 1e-4) || (abs(p.old[8,1]) < 1e-4)
    
    if (all_lognormal_triggered) {
      global_lognormal_iter_count = global_lognormal_iter_count + 1
    }
   
    h1.0 = haz_fun(y0, p.old[2,1], p.old[3,1], p.old[4,1])
    h1.1 = haz_fun(y1, p.old[2,1], p.old[3,1], p.old[4,1])
    h1.2 = haz_fun(y2, p.old[2,1], p.old[3,1], p.old[4,1])
    h1.3 = haz_fun(y3, p.old[2,1], p.old[3,1], p.old[4,1])
    
    h2.0 = haz_fun(y0, p.old[5,1], p.old[6,1], p.old[7,1])
    h2.1 = haz_fun(y1, p.old[5,1], p.old[6,1], p.old[7,1])
    h2.2 = haz_fun(y2, p.old[5,1], p.old[6,1], p.old[7,1])
    h2.3 = haz_fun(y3, p.old[5,1], p.old[6,1], p.old[7,1])
    
    h3.0 = haz_fun(y0, p.old[8,1], p.old[9,1], p.old[10,1])
    h3.1 = haz_fun(y1, p.old[8,1], p.old[9,1], p.old[10,1])
    h3.2 = haz_fun(y2, p.old[8,1], p.old[9,1], p.old[10,1])
    h3.3 = haz_fun(y3, p.old[8,1], p.old[9,1], p.old[10,1])
    
    H1.0 = cum_haz_fun(y0, p.old[2,1], p.old[3,1], p.old[4,1])
    H1.1 = cum_haz_fun(y1, p.old[2,1], p.old[3,1], p.old[4,1])
    H1.2 = cum_haz_fun(y2, p.old[2,1], p.old[3,1], p.old[4,1])
    H1.3 = cum_haz_fun(y3, p.old[2,1], p.old[3,1], p.old[4,1])
    
    H2.0 = cum_haz_fun(y0, p.old[5,1], p.old[6,1], p.old[7,1])
    H2.1 = cum_haz_fun(y1, p.old[5,1], p.old[6,1], p.old[7,1])
    H2.2 = cum_haz_fun(y2, p.old[5,1], p.old[6,1], p.old[7,1])
    H2.3 = cum_haz_fun(y3, p.old[5,1], p.old[6,1], p.old[7,1])
    
    H3.0 = cum_haz_fun(y0, p.old[8,1], p.old[9,1], p.old[10,1])
    H3.1 = cum_haz_fun(y1, p.old[8,1], p.old[9,1], p.old[10,1])
    H3.2 = cum_haz_fun(y2, p.old[8,1], p.old[9,1], p.old[10,1])
    H3.3 = cum_haz_fun(y3, p.old[8,1], p.old[9,1], p.old[10,1])
    
    w0 = ((1-p.old[1,1])*exp(-H1.0-H2.0-H3.0))/((p.old[1,1]*exp(-H2.0-H3.0)) + ((1-p.old[1,1])*exp(-H1.0-H2.0-H3.0))) 
    w2 = ((1-p.old[1,1])*exp(-H1.2-H2.2-H3.2))/((p.old[1,1]*exp(-H2.2-H3.2)) + ((1-p.old[1,1])*exp(-H1.2-H2.2-H3.2))) 
    w3 = ((1-p.old[1,1])*exp(-H1.3-H2.3-H3.3))/((p.old[1,1]*exp(-H2.3-H3.3)) + ((1-p.old[1,1])*exp(-H1.3-H2.3-H3.3))) 
    
  
    Q = function(par) {
      res = (length(y1)*log(1-par[1])) + 
        sum(((1-w0)*log(par[1]))+(w0*log(1-par[1]))) + 
        sum(((1-w2)*log(par[1]))+(w2*log(1-par[1]))) + 
        sum(((1-w3)*log(par[1]))+(w3*log(1-par[1]))) 
      return(res)
    }
    c.new = optimize(Q, c(1e-5, 1-1e-5), maximum = TRUE)$maximum
    
    Q1_gg = function(par1) {
      res1 = sum(log(haz_fun(y1, par1[1], par1[2], par1[3])) - cum_haz_fun(y1, par1[1], par1[2], par1[3])) - 
        sum(w0*cum_haz_fun(y0, par1[1], par1[2], par1[3])) - 
        sum(w2*cum_haz_fun(y2, par1[1], par1[2], par1[3])) - 
        sum(w3*cum_haz_fun(y3, par1[1], par1[2], par1[3]))
      return(-res1)
    }
    risk1.new = optim(par=p.old[2:4,1], fn=Q1_gg, method="Nelder-Mead")$par
    
    
    Q2_gg = function(par2) {
      res2 = sum(log(haz_fun(y2, par2[1], par2[2], par2[3]))) - 
        sum(cum_haz_fun(y1, par2[1], par2[2], par2[3])) - 
        sum(cum_haz_fun(y0, par2[1], par2[2], par2[3])) - 
        sum(cum_haz_fun(y2, par2[1], par2[2], par2[3])) - 
        sum(cum_haz_fun(y3, par2[1], par2[2], par2[3]))
      return(-res2)
    }
    risk2.new = optim(par=p.old[5:7,1], fn=Q2_gg, method="Nelder-Mead")$par
    
    Q3_gg = function(par3) {
      res3 = sum(log(haz_fun(y3, par3[1], par3[2], par3[3]))) - 
        sum(cum_haz_fun(y1, par3[1], par3[2], par3[3])) - 
        sum(cum_haz_fun(y0, par3[1], par3[2], par3[3])) - 
        sum(cum_haz_fun(y2, par3[1], par3[2], par3[3])) - 
        sum(cum_haz_fun(y3, par3[1], par3[2], par3[3]))
      return(-res3)
    }
    risk3.new = optim(par=p.old[8:10,1], fn=Q3_gg, method="Nelder-Mead")$par
    
    
    p.new = matrix(c(c.new, risk1.new, risk2.new, risk3.new), ncol = 1)
    

    continue = (any(abs((p.new - p.old) / (p.old )) > tol)) && (iter < maxit)
    
    
    p.old = p.new
    iter = iter + 1
  }
  
  print(paste("Converged in iterations:", iter))
  print(paste("Total iterations running systemic lognormal constraints:", global_lognormal_iter_count))
  
  
  final_logLik <- log.lik.fn(
    p.new,
    y0 = y0,
    y1 = y1,
    y2 = y2,
    y3 = y3
  )
  
  print(paste(
    "Maximized Log likelihood after EM algorithm =",
    final_logLik
  ))
  
  final_lognorm_flag = (abs(p.new[2,1]) < 1e-4) || (abs(p.new[5,1]) < 1e-4) || (abs(p.new[8,1]) < 1e-4)
  
  active_indices = c(1) 
  
  active_indices = c(active_indices, 2, 3, 4, 5, 6, 7, 8, 9, 10)
  
  
  p.active = p.new[active_indices, 1]
  
  log_lik_dynamic_profile = function(reduced_p, y0, y1, y2, y3, final_lognorm_flag) {
    full_p = numeric(10)
    full_p[1] = reduced_p[1]
    
    full_p[2:10] = reduced_p[2:10]                 
    result <- log.lik.fn(matrix(full_p, ncol=1), y0=y0, y1=y1, y2=y2, y3=y3)

    return(result)
  }
  
  hess.mat = hessian(log_lik_dynamic_profile, p.active, y0=y0, y1=y1, y2=y2, y3=y3, final_lognorm_flag=final_lognorm_flag)
  
  if (any(!is.finite(hess.mat))) {
    return(list(estimates = rep(0,20), global_lognormal_iter_count = global_lognormal_iter_count))
  }
  
  FI = try(solve(-1 * hess.mat), silent=TRUE)
  
  if (inherits(FI, "try-error") || any(!is.finite(FI)) || any(diag(FI) <= 0)) {
    return(list(estimates = rep(0,20), global_lognormal_iter_count = global_lognormal_iter_count))
  }
  
  se.se_est = numeric(10)
  se.se_est[active_indices] = sqrt(diag(FI))
  
  
  EM.out = if(sum(is.na(se.se_est))==0) c(p.new, se.se_est) else rep(0,20)
  return(list(estimates = EM.out, global_lognormal_iter_count = global_lognormal_iter_count))
}#end of the EM function


# Likelihood Ratio Test  Function 

compute_lrt_pvalue_weibull = function(data, em_estimates) {
  data0 = data[data$d == 0, ]; y0 = data0$y
  data1 = data[data$d == 1, ]; y1 = data1$y
  data2 = data[data$d == 2, ]; y2 = data2$y
  data3 = data[data$d == 3, ]; y3 = data3$y
  
  param_full_init <- em_estimates[1:10]
  
  L_full <- function(p){
    
    if(p[1] <= 0 || p[1] >= 1)
      return(1e10)
     if(any(p[c(3,4,6,7,9,10)] <= 0))
      return(1e10)
    val <- -log.lik.fn(p, y0, y1, y2, y3)
    if(is.na(val) || is.nan(val) || is.infinite(val))
      return(1e10)
    
    return(val)
  }
  
  fit_full <- try(optim(
    par = param_full_init, 
    fn = L_full, 
    method = "Nelder-Mead",
    control = list(maxit = 5000, reltol = 1e-12, parscale = abs(param_full_init))
  ), silent = TRUE)
  
  if (inherits(fit_full, "try-error")) {
    return(NA)
    cat("Full GG Optimization failed with an error.\n")
  } else if (fit_full$convergence == 0) {
    cat("Full GG Optimization converged successfully.\n")
  } else {
    cat("Full GG Optimization did NOT converge.\n")
    cat("Convergence code:", fit_full$convergence, "\n")
  }
  
  logL_full <- -fit_full$value

  cat("\n============================================\n")
  cat("FULL GENERALIZED GAMMA MODEL\n")
  cat("============================================\n")
  
  cat("Full GG log-likelihood:", logL_full, "\n\n")
  
  
  # Fit Restricted Model (Weibull: q=1) 
  r1 <- which.min(abs(profile_results$q - 1))
  init_weibull <- as.numeric(profile_results[r1, c("cure_hat", "lambda1_hat", "sigma1_hat",
                                                   "lambda2_hat", "sigma2_hat",
                                                   "lambda3_hat", "sigma3_hat")])
  
  L_W = function(p) {
    if (p[1] <= 0 || p[1] >= 1 || any(p[2:7] <= 0)) return(1e10)
    val <- -log.lik.weibull(p, y0, y1, y2, y3)
    if (is.na(val) || is.nan(val) || is.infinite(val)) return(1e10)
    return(val)
  }
  
  fit_weibull <- try(optim(
    par = init_weibull, 
    fn = L_W, 
    method = "Nelder-Mead",
    control = list(maxit = 7000, reltol = 1e-12)
  ), silent = TRUE)
  
  if (inherits(fit_weibull, "try-error")) {
    return(NA)
    cat("Reduced Weibull Optimization failed with an error.\n")
  } else if (fit_weibull$convergence == 0) {
    cat("Reduced Weibull Optimization converged successfully.\n")
  } else {
    cat("Reduced Weibull Optimization did NOT converge.\n")
    cat("Convergence code:", fit_weibull$convergence, "\n")
  }
  
  logL_reduced_weibull <- -fit_weibull$value
  estimates_reduced_weibull <- fit_weibull$par
  
  lrt_stat_weibull <- 2 * (logL_full - logL_reduced_weibull)
  p_val_weibull    <- pchisq(lrt_stat_weibull, df = 3, lower.tail = FALSE)
  return(list(
    lrt_stat_weibull = lrt_stat_weibull,
    p_val_weibull = p_val_weibull,
    logL_full = logL_full,
    logL_reduced_weibull = logL_reduced_weibull,
    estimates_reduced_weibull = estimates_reduced_weibull
  )) 
  
}

compute_lrt_pvalue_gamma = function(data, em_estimates) {
  data0 = data[data$d == 0, ]; y0 = data0$y
  data1 = data[data$d == 1, ]; y1 = data1$y
  data2 = data[data$d == 2, ]; y2 = data2$y
  data3 = data[data$d == 3, ]; y3 = data3$y
  
  param_full_init <- em_estimates[1:10]
  
  L_full <- function(p){
    
    if(p[1] <= 0 || p[1] >= 1)
      return(1e10)
    if(any(p[c(3,4,6,7,9,10)] <= 0))
      return(1e10)
    val <- -log.lik.fn(p, y0, y1, y2, y3)
    if(is.na(val) || is.nan(val) || is.infinite(val))
      return(1e10)
    return(val)
  }
  fit_full <- try(optim(
    par = param_full_init, 
    fn = L_full, 
    method = "Nelder-Mead",
    control = list(maxit = 5000, reltol = 1e-12, parscale = abs(param_full_init))
  ), silent = TRUE)
  
  if (inherits(fit_full, "try-error")) {
    return(NA)
    cat("Full GG Optimization failed with an error.\n")
  } else if (fit_full$convergence == 0) {
    cat("Full GG Optimization converged successfully.\n")
  } else {
    cat("Full GG Optimization did NOT converge.\n")
    cat("Convergence code:", fit_full$convergence, "\n")
  }
  
  logL_full <- -fit_full$value
  
  # logL_full <- -L_full(param_full_init)
  cat("\n============================================\n")
  cat("FULL GENERALIZED GAMMA MODEL\n")
  cat("============================================\n")
  
  cat("Full GG log-likelihood:", logL_full, "\n\n")
  
  # Fit Restricted Model (Gamma: q=sig)
  lambda1.init <- 1 / mean1
  sigma1.init  <- sqrt(var1) / mean1
  
  lambda2.init <- 1 / mean2
  sigma2.init  <- sqrt(var2) / mean2
  
  lambda3.init <- 1 / mean3
  sigma3.init  <- sqrt(var3) / mean3
   init_gamma <- c(
    cure.init,
    lambda1.init, sigma1.init,
    lambda2.init, sigma2.init,
    lambda3.init, sigma3.init
  )
  L_G = function(p) {
    if (p[1] <= 0 || p[1] >= 1 || any(p[2:7] <= 0)) {
      return(1e10) 
    }
    
    val <- -log.lik.gamma(p, y0, y1, y2, y3)
    
    if (is.na(val) || is.nan(val) || is.infinite(val)) {
      return(1e10)
    }
    
    return(val)
  }
  
  fit_gamma <- try(optim(
    par = init_gamma, 
    fn = L_G, 
    method = "Nelder-Mead",
    control = list(maxit = 7000, reltol = 1e-12)
  ), silent = TRUE)
  
  if (inherits(fit_gamma, "try-error")) {
    return(NA)
    cat("Reduced gamma Optimization failed with an error.\n")
  } else if (fit_gamma$convergence == 0) {
    cat("Reduced gamma Optimization converged successfully.\n")
  } else {
    cat("Reduced gamma Optimization did NOT converge.\n")
    cat("Convergence code:", fit_gamma$convergence, "\n")
  }
  
  
  logL_reduced_gamma <- -fit_gamma$value
  estimates_reduced_gamma <- fit_gamma$par
  
  lrt_stat_gamma <- 2 * (logL_full - logL_reduced_gamma)
  
  p_val_gamma <- pchisq(lrt_stat_gamma, df = 3, lower.tail = FALSE)
  return(list(
    lrt_stat_gamma = lrt_stat_gamma,
    p_val_gamma = p_val_gamma,
    logL_full = logL_full,
    logL_reduced_gamma = logL_reduced_gamma,
    estimates_reduced_gamma = estimates_reduced_gamma
  ))
  
}

compute_lrt_pvalue_lognormal = function(data, em_estimates) {
  data0 = data[data$d == 0, ]; y0 = data0$y
  data1 = data[data$d == 1, ]; y1 = data1$y
  data2 = data[data$d == 2, ]; y2 = data2$y
  data3 = data[data$d == 3, ]; y3 = data3$y
  
  param_full_init <- em_estimates[1:10]
  
  L_full <- function(p){
    
    if(p[1] <= 0 || p[1] >= 1)
      return(1e10)
    if(any(p[c(3,4,6,7,9,10)] <= 0))
      return(1e10)
     val <- -log.lik.fn(p, y0, y1, y2, y3)
     if(is.na(val) || is.nan(val) || is.infinite(val))
      return(1e10)
    
    return(val)
  }
  fit_full <- try(optim(
    par = param_full_init, 
    fn = L_full, 
    method = "Nelder-Mead",
    control = list(maxit = 5000, reltol = 1e-12, parscale = abs(param_full_init))
  ), silent = TRUE)
  
  if (inherits(fit_full, "try-error")) {
    return(NA)
    cat("Full GG Optimization failed with an error.\n")
  } else if (fit_full$convergence == 0) {
    cat("Full GG Optimization converged successfully.\n")
  } else {
    cat("Full GG Optimization did NOT converge.\n")
    cat("Convergence code:", fit_full$convergence, "\n")
  }
  
  
  
  logL_full <- -fit_full$value
  cat("\n============================================\n")
  cat("FULL GENERALIZED GAMMA MODEL\n")
  cat("============================================\n")
  
  cat("Full GG log-likelihood:", logL_full, "\n\n")
 
   n=1617
  AIC_GenGamma <- -2 * logL_full + 2 * 10
  BIC_GenGamma <- -2 * logL_full + 10 * log(n)
  AIC_GenGamma <- round(AIC_GenGamma, 6)
  BIC_GenGamma <- round(BIC_GenGamma, 6)
  
  cat("\n--------------------------------------------\n")
  cat("GENERALIZED GAMMA MODEL\n")
  cat("Log-Likelihood:", round(logL_full, 6), "\n")
  cat("Number of parameters:", 10, "\n")
  cat("AIC:", AIC_GenGamma, "\n")
  cat("BIC:", BIC_GenGamma, "\n")
  # cat("--------------------------------------------\n")
  
  r0 <- which.min(abs(profile_results$q - 0))
  init_lognormal <- as.numeric(profile_results[r0, c("cure_hat", "lambda1_hat", "sigma1_hat",
                                                     "lambda2_hat", "sigma2_hat",
                                                     "lambda3_hat", "sigma3_hat")])
  
  L_LG = function(p) {
    if (p[1] <= 0 || p[1] >= 1 || any(p[2:7] <= 0)) {
      return(1e10) 
    }
    
    val <- -log.lik.lognormal(p, y0, y1, y2, y3)
        if (is.na(val) || is.nan(val) || is.infinite(val)) {
      return(1e10)
    }
    
    return(val)
  }
  
  fit_lognormal <- try(optim(
    par = init_lognormal, 
    fn = L_LG, 
    method = "Nelder-Mead",
    control = list(maxit = 5000, reltol = 1e-12)
  ), silent = TRUE)
  
  if (inherits(fit_lognormal, "try-error")) {
    return(NA)
    cat("Reduced lognormal Optimization failed with an error.\n")
  } else if (fit_lognormal$convergence == 0) {
    cat("Reduced lognormal Optimization converged successfully.\n")
  } else {
    cat("Reduced lognormal Optimization did NOT converge.\n")
    cat("Convergence code:", fit_lognormal$convergence, "\n")
  }
  
  logL_reduced_lognormal <- -fit_lognormal$value
  estimates_reduced_lognormal <- fit_lognormal$par
  
  lrt_stat_lognormal <- 2 * (logL_full - logL_reduced_lognormal)
  
  
  
  p_val_lognormal <- 0.5*pchisq(lrt_stat_lognormal, df = 3, lower.tail = FALSE)
  return(list(
    lrt_stat_lognormal = lrt_stat_lognormal,
    p_val_lognormal = p_val_lognormal,
    logL_full = logL_full,
    logL_reduced_lognormal = logL_reduced_lognormal,
    estimates_reduced_lognormal = estimates_reduced_lognormal
  ))
  
  
}


# 1. Read real data

data_raw <- read.csv("File name")

data <- data.frame(
  y = data_raw[["Time in years"]],
  d = data_raw[["d"]]
)
data$y <- as.numeric(data$y)
data$d <- as.numeric(data$d)

data <- na.omit(data)

library(survival)
J <- ifelse(data$d != 0, 1, 0)
data.init <- data.frame(
  Y = data$y,
  J = J
)

# Kaplan-Meier estimate
KM.init <- survfit(Surv(Y, J) ~ 1, data = data.init)

plot(
  KM.init$time,
  KM.init$surv,
  type = "s",
   xlab = "Time (years)",
  ylab = "Survival Probability",
  ylim = c(0, 1),
  yaxt = "n"
)
axis(
  2,
  at = seq(0, 1, by = 0.2),
  labels = seq(0, 1, by = 0.2)
)

cure.init <- tail(KM.init$surv, 1)
cat("Initial cure-rate value =", cure.init, "\n")

mean_gg <- function(q, lambda, sigma) {
  
  if (lambda <= 0) {
    stop("lambda must be > 0")
  }
  
  if (sigma <= 0) {
    stop("sigma must be > 0")
  }
    if (abs(q) < 1e-4) {
    
    mu <- exp(sigma^2 / 2) / lambda
    
    return(mu)
  }
   a <- q^(-2)
  moment_argument <- a + (sigma / q)
  
  if (moment_argument <= 0) {
    
    return(NA_real_)
  }
  
  mu <- (abs(q)^(2 * sigma / q) / lambda) *
    (gamma(moment_argument) / gamma(a))
  
  return(mu)
}

var_gg <- function(q, lambda, sigma) {
  
  if (lambda <= 0) {
    stop("lambda must be > 0")
  }
  
  if (sigma <= 0) {
    stop("sigma must be > 0")
  }
  if (abs(q) < 1e-4) {
    
    variance <- (exp(sigma^2) / lambda^2) *
      (exp(sigma^2) - 1)
    
    return(variance)
  }
  a <- q^(-2)

  mean_argument <- a + (sigma / q)
  second_argument <- a + (2 * sigma / q)
  
  if (mean_argument <= 0) {
    
    return(NA_real_)
  }
  
  if (second_argument <= 0) {
    
    return(NA_real_)
  }
  
  mean_x <- (abs(q)^(2 * sigma / q) / lambda) *
    (gamma(mean_argument) / gamma(a))
 
  second_moment <- (abs(q)^(4 * sigma / q) / lambda^2) *
    (gamma(second_argument) / gamma(a))
  
  variance <- second_moment - mean_x^2
  
  return(variance)
}

data1 <- data$y[data$d == 1]
mean1 <- mean(data1)
var1  <- var(data1)
cat("Risk 1: Breast Cancer\n")
cat("Number of observations =", length(data1), "\n")
cat("Sample mean =", mean1, "\n")
cat("Sample variance =", var1, "\n\n")


data2 <- data$y[data$d == 2]
mean2 <- mean(data2)
var2  <- var(data2)
cat("Risk 2: Other Cancers\n")
cat("Number of observations =", length(data2), "\n")
cat("Sample mean =", mean2, "\n")
cat("Sample variance =", var2, "\n\n")


data3 <- data$y[data$d == 3]
mean3 <- mean(data3)
var3  <- var(data3)
cat("Risk 3: Other Causes\n")
cat("Number of observations =", length(data3), "\n")
cat("Sample mean =", mean3, "\n")
cat("Sample variance =", var3, "\n\n")

# profile likelihood on q
library(BB)
q_values <- seq(-10, 10, by = 0.1)

y0 <- data$y[data$d == 0]
y1 <- data$y[data$d == 1]
y2 <- data$y[data$d == 2]
y3 <- data$y[data$d == 3]

data1 <- data$y[data$d == 1]
mean1 <- mean(data1)
var1  <- var(data1)

data2 <- data$y[data$d == 2]
mean2 <- mean(data2)
var2  <- var(data2)

data3 <- data$y[data$d == 3]
mean3 <- mean(data3)
var3  <- var(data3)


cat("\n============================================================\n")
cat("SAMPLE MOMENTS\n")
cat("============================================================\n")

cat("Risk 1: mean =", mean1, " variance =", var1, "\n")
cat("Risk 2: mean =", mean2, " variance =", var2, "\n")
cat("Risk 3: mean =", mean3, " variance =", var3, "\n")

solve_moments_BB <- function(qi, sample_mean, sample_var) {
  moment_equations <- function(theta) {
    
    eta <- theta[1]
    xi  <- theta[2]
     if (!is.finite(eta) ||
        !is.finite(xi)) {
      
      return(c(1e10, 1e10))
    }
    
   if (eta < -50 || eta > 50 ||
        xi < -50 || xi > 50) {
      
      return(c(1e10, 1e10))
    }
   lambda <- exp(eta)
   if (abs(qi) < 1e-4) {
      sigma <- exp(xi)
      
    } else if (qi < 0) {
      sigma_max <- 1 / (2 * abs(qi))
      
      sigma <- sigma_max * plogis(xi)
      
    } else {
         sigma <- exp(xi)
    }
    if (!is.finite(lambda) ||
        !is.finite(sigma) ||
        lambda <= 0 ||
        sigma <= 0) {
      
      return(c(1e10, 1e10))
    }
    
    theoretical_mean <- tryCatch(
      
      mean_gg(
        q = qi,
        lambda = lambda,
        sigma = sigma
      ),
      
      error = function(e) NA_real_
    )
    
     theoretical_var <- tryCatch(
      
      var_gg(
        q = qi,
        lambda = lambda,
        sigma = sigma
      ),
      
      error = function(e) NA_real_
    )
    
    if (!is.finite(theoretical_mean) ||
        !is.finite(theoretical_var)) {
      
      return(c(1e10, 1e10))
    }
     return(
      c(
        theoretical_mean - sample_mean,
        theoretical_var  - sample_var
      )
    )
  }
  
 fit_BB <- tryCatch(
    
    BBsolve(
      par = c(0.5,0.5),
      fn = moment_equations,
      control = list(
        maxit = 5000,
        tol = 1e-8
      )
    ),
    
    error = function(e) NULL
  )
  
  if (is.null(fit_BB)) {
    return(NULL)
  }
  
 theta_hat <- fit_BB$par
  
  
  if (length(theta_hat) != 2 ||
      any(!is.finite(theta_hat))) {
    
    return(NULL)
  }
  
  eta_hat <- theta_hat[1]
  xi_hat  <- theta_hat[2]
  
  if (eta_hat < -50 || eta_hat > 50 ||
      xi_hat < -50 || xi_hat > 50) {
    
    return(NULL)
  }
  
  lambda_hat <- exp(eta_hat)
  
  
  if (abs(qi) < 1e-4) {
    
    sigma_hat <- exp(xi_hat)
    
  } else if (qi < 0) {
    
    sigma_max <- 1 / (2 * abs(qi))
    
    sigma_hat <- sigma_max * plogis(xi_hat)
    
  } else {
    
    sigma_hat <- exp(xi_hat)
  }
   if (!is.finite(lambda_hat) ||
      !is.finite(sigma_hat) ||
      lambda_hat <= 0 ||
      sigma_hat <= 0) {
    
    return(NULL)
  }
  
  residual <- moment_equations(theta_hat)
  
  residual_norm <- sqrt(
    sum(residual^2)
  )
  
  if (!is.finite(residual_norm)) {
    return(NULL)
  }
  return(
    list(
      
      lambda = lambda_hat,
      
      sigma = sigma_hat,
      
      residual_norm = residual_norm,
      
      convergence = fit_BB$convergence
      
    )
  )
}

fit_fixed_q <- function(
    qi,
    cure_start,
    lambda1_start,
    sigma1_start,
    lambda2_start,
    sigma2_start,
    lambda3_start,
    sigma3_start,
    y0, y1, y2, y3
) {
  
  start <- c(
    
    cure_start,
    
    lambda1_start,
    sigma1_start,
    
    lambda2_start,
    sigma2_start,
    
    lambda3_start,
    sigma3_start
  )
  
  objective <- function(p) {
    if (!is.finite(p[1]) ||
        p[1] <= 0 ||
        p[1] >= 1) {
      
      return(1e10)
    }
     if (any(!is.finite(p[2:7]))) {
      return(1e10)
    }
    
    
    if (any(p[2:7] <= 0)) {
      return(1e10)
    }
     if (qi < 0) {
      
      sigma_max <- 1 / (2 * abs(qi))
      
      if (p[3] >= sigma_max ||
          p[5] >= sigma_max ||
          p[7] >= sigma_max) {
        
        return(1e10)
      }
    }
    
    value <- tryCatch(
      
      -log.lik.gg.fixed.q(
        
        param_gengamma = p,
        
        q = qi,
        
        y0 = y0,
        y1 = y1,
        y2 = y2,
        y3 = y3
        
      ),
      
      error = function(e) NA_real_
    )
    if (!is.finite(value)) {
      return(1e10)
    }
    
    
    return(value)
  }
  
  fit <- tryCatch(
    
    optim(
      
      par = start,
      
      fn = objective,
      
      method = "Nelder-Mead",
      
      control = list(
        
        maxit = 10000,
        
        reltol = 1e-12,
        
        parscale = pmax(
          abs(start),
          1e-4
        )
      )
    ),
    
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(NULL)
  }
  
  
  if (!is.finite(fit$value)) {
    return(NULL)
  }
  
  return(
    list(
      
      par = fit$par,
      
      logLik = -fit$value,
      
      convergence = fit$convergence,
      
      counts = fit$counts,
      
      initial = start
      
    )
  )
}


profile_list <- list()
profile_counter <- 0


cat("\n\n============================================================\n")
cat("STARTING FIXED-q PROFILE\n")
cat("============================================================\n\n")


for (qi in q_values) {
  
  
  cat("\n------------------------------------------------------------\n")
  cat(
    "q =",
    sprintf("%.1f", qi),
    "\n"
  )
  cat("------------------------------------------------------------\n")
  
  
  cat("BBsolve Risk 1 ...\n")
  
  risk1_BB <- solve_moments_BB(
    
    qi = qi,
    
    sample_mean = mean1,
    
    sample_var = var1
  )
  
  
  if (is.null(risk1_BB)) {
    
    cat("Risk 1 BBsolve failed.\n")
    cat("Skipping q =", qi, "\n")
    
    next
  }
  
  
  cat(
    "Risk 1 BBsolve:",
    "lambda =",
    sprintf("%.10f", risk1_BB$lambda),
    "sigma =",
    sprintf("%.10f", risk1_BB$sigma),
    "\n"
  )
  
  cat("BBsolve Risk 2 ...\n")
  
  risk2_BB <- solve_moments_BB(
    
    qi = qi,
    
    sample_mean = mean2,
    
    sample_var = var2
  )
  
  
  if (is.null(risk2_BB)) {
    
    cat("Risk 2 BBsolve failed.\n")
    cat("Skipping q =", qi, "\n")
    
    next
  }
  
  
  cat(
    "Risk 2 BBsolve:",
    "lambda =",
    sprintf("%.10f", risk2_BB$lambda),
    "sigma =",
    sprintf("%.10f", risk2_BB$sigma),
    "\n"
  )
  
  cat("BBsolve Risk 3 ...\n")
  
  risk3_BB <- solve_moments_BB(
    
    qi = qi,
    
    sample_mean = mean3,
    
    sample_var = var3
  )
  
  
  if (is.null(risk3_BB)) {
    
    cat("Risk 3 BBsolve failed.\n")
    cat("Skipping q =", qi, "\n")
    
    next
  }
  
  
  cat(
    "Risk 3 BBsolve:",
    "lambda =",
    sprintf("%.10f", risk3_BB$lambda),
    "sigma =",
    sprintf("%.10f", risk3_BB$sigma),
    "\n"
  )
  
   fixed_q_start <- c(
    
    cure.init,
    
    risk1_BB$lambda,
    risk1_BB$sigma,
    
    risk2_BB$lambda,
    risk2_BB$sigma,
    
    risk3_BB$lambda,
    risk3_BB$sigma
  )
  
  
  cat("\nInitial values supplied to Nelder-Mead:\n")
  print(
    fixed_q_start,
    digits = 12
  )
  
 fixed_fit <- fit_fixed_q(
    
    qi = qi,
    
    cure_start = cure.init,
    
    lambda1_start = risk1_BB$lambda,
    sigma1_start  = risk1_BB$sigma,
    
    lambda2_start = risk2_BB$lambda,
    sigma2_start  = risk2_BB$sigma,
    
    lambda3_start = risk3_BB$lambda,
    sigma3_start  = risk3_BB$sigma,
    
    y0 = y0,
    y1 = y1,
    y2 = y2,
    y3 = y3
  )
  
  
  if (is.null(fixed_fit)) {
    
    cat("Fixed-q Nelder-Mead failed.\n")
    cat("Skipping q =", qi, "\n")
    
    next
  }
  
  p_hat <- fixed_fit$par
  
  
  cat("\nFinal Nelder-Mead estimates:\n")
  
  cat(
    "cure    =",
    sprintf("%.10f", p_hat[1]),
    "\n"
  )
  
  cat(
    "lambda1 =",
    sprintf("%.10f", p_hat[2]),
    "\n"
  )
  
  cat(
    "sigma1  =",
    sprintf("%.10f", p_hat[3]),
    "\n"
  )
  
  cat(
    "lambda2 =",
    sprintf("%.10f", p_hat[4]),
    "\n"
  )
  
  cat(
    "sigma2  =",
    sprintf("%.10f", p_hat[5]),
    "\n"
  )
  
  cat(
    "lambda3 =",
    sprintf("%.10f", p_hat[6]),
    "\n"
  )
  
  cat(
    "sigma3  =",
    sprintf("%.10f", p_hat[7]),
    "\n"
  )
  
  
  cat(
    "Log-likelihood =",
    sprintf("%.10f", fixed_fit$logLik),
    "\n"
  )
  
  
  cat(
    "Convergence code =",
    fixed_fit$convergence,
    "\n"
  )
  profile_counter <- profile_counter + 1
  
  
  profile_list[[profile_counter]] <- data.frame(
    
    q = qi,
    
    # BBsolve values
    lambda1_BB = risk1_BB$lambda,
    sigma1_BB  = risk1_BB$sigma,
    
    lambda2_BB = risk2_BB$lambda,
    sigma2_BB  = risk2_BB$sigma,
    
    lambda3_BB = risk3_BB$lambda,
    sigma3_BB  = risk3_BB$sigma,
    
    # Final Nelder-Mead values
    cure_hat = p_hat[1],
    
    lambda1_hat = p_hat[2],
    sigma1_hat  = p_hat[3],
    
    lambda2_hat = p_hat[4],
    sigma2_hat  = p_hat[5],
    
    lambda3_hat = p_hat[6],
    sigma3_hat  = p_hat[7],
    
    # Likelihood
    logLik = fixed_fit$logLik,
    
    # Optimization information
    convergence = fixed_fit$convergence,
    
    # BBsolve residuals
    risk1_moment_error = risk1_BB$residual_norm,
    risk2_moment_error = risk2_BB$residual_norm,
    risk3_moment_error = risk3_BB$residual_norm
  )
}
if (length(profile_list) == 0) {
  
  stop(
    "No q values produced a successful BBsolve + ",
    "fixed-q Nelder-Mead optimization."
  )
}


profile_results <- do.call(
  rbind,
  profile_list
)


rownames(profile_results) <- NULL


cat("\n\n============================================================\n")
cat("COMPLETE FIXED-q PROFILE RESULTS\n")
cat("============================================================\n\n")


print(
  profile_results,
  row.names = FALSE,
  digits = 10
)


best_row <- which.max(
  profile_results$logLik
)


best_q <- profile_results$q[best_row]


cat("\n\n============================================================\n")
cat("BEST q FROM FIXED-q PROFILE\n")
cat("============================================================\n\n")


cat(
  "Best q =",
  sprintf("%.10f", best_q),
  "\n"
)
cat(
  "Maximum fixed-q log-likelihood =",
  sprintf(
    "%.10f",
    profile_results$logLik[best_row]
  ),
  "\n\n"
)

cat("Corresponding estimates:\n")
cat(
  "cure    =",
  sprintf("%.10f", profile_results$cure_hat[best_row]),
  "\n"
)
cat(
  "lambda1 =",
  sprintf("%.10f", profile_results$lambda1_hat[best_row]),
  "\n"
)
cat(
  "sigma1  =",
  sprintf("%.10f", profile_results$sigma1_hat[best_row]),
  "\n"
)
cat(
  "lambda2 =",
  sprintf("%.10f", profile_results$lambda2_hat[best_row]),
  "\n"
)
cat(
  "sigma2  =",
  sprintf("%.10f", profile_results$sigma2_hat[best_row]),
  "\n"
)
cat(
  "lambda3 =",
  sprintf("%.10f", profile_results$lambda3_hat[best_row]),
  "\n"
)
cat(
  "sigma3  =",
  sprintf("%.10f", profile_results$sigma3_hat[best_row]),
  "\n"
)

profile_plot <- profile_results[
  profile_results$logLik > -10^10,
]

# Maximum log-likelihood
max_logLik <- max(profile_plot$logLik, na.rm = TRUE)

plot(
  profile_plot$q,
  profile_plot$logLik,
  type = "b",
  pch = 16,
  cex = 0.5,
  xlab = "q",
  ylab = "Profile log-likelihood",
  ylim = c(-5000, -3000)
)

abline(
  v = best_q,
  lty = 2
)

text(
  best_q,
  -4950,
  labels = paste0("q = ", round(best_q, 2)),
  pos = 3
)

abline(
  h = max_logLik,
  lty = 2
)

text(
  min(profile_plot$q, na.rm = TRUE),
  max_logLik,
  labels = paste0("                        logLik = ", round(max_logLik, 2)),
  pos = 3,
  cex = 0.8
)

# Use the initial values and corresponding q with maximum log likelihood value got from running the above code 
c.init <- 0.0248728748
q1.init <- 0.5
lam1.init <- 0.2340058151 
sig1.init <- 1.3255027620 
q2.init <- 0.5
lam2.init <- 0.0009520608 
sig2.init <- 2.3456369356 
q3.init <- 0.5
lam3.init <- 0.0112895948 
sig3.init <- 1.5751719889 

alpha = 0.05

EM.result <- EM.MCCR3_global_lognormal(
  data = data,
  tol = 0.001,
  maxit = 500,
  c = c.init,
  q_1 = q1.init,
  lambda1 = lam1.init,
  sigma1 = sig1.init,
  q_2 = q2.init,
  lambda2 = lam2.init,
  sigma2 = sig2.init,
  q_3 = q3.init,
  lambda3 = lam3.init,
  sigma3 = sig3.init
)

EM <- EM.result$estimates

if (all(EM == 0)) {
  
  stop("The EM algorithm did not converge. Try different initial values.")
  
}

estimates <- EM[1:10]
SE <- EM[11:20]

names(estimates) <- c(
  "c",
  "q1", "lambda1", "sigma1",
  "q2", "lambda2", "sigma2",
  "q3", "lambda3", "sigma3"
)

names(SE) <- c(
  "c",
  "q1", "lambda1", "sigma1",
  "q2", "lambda2", "sigma2",
  "q3", "lambda3", "sigma3"
)

parameter_table <- data.frame(
  Parameter = names(estimates),
  Estimate = as.numeric(estimates),
  SE = as.numeric(SE)
)

parameter_table$Estimate <- round(parameter_table$Estimate, 6)
parameter_table$SE <- round(parameter_table$SE, 6)

cat("\n==================================================\n")
cat("        GENERALIZED GAMMA MODEL FIT\n")
cat("==================================================\n")
#cat("Sample size:", n, "\n\n")

print(parameter_table)

weibull_LRT <- compute_lrt_pvalue_weibull(
  data,
  estimates
)

gamma_LRT <- compute_lrt_pvalue_gamma(
  data,
  estimates
)

lognormal_LRT <- compute_lrt_pvalue_lognormal(
  data,
  estimates
)

LRT_table <- data.frame(
  Model = c(
    "Weibull",
    "Gamma",
    "Lognormal"
  ),
  
  LRT = c(
    weibull_LRT$lrt_stat_weibull,
    gamma_LRT$lrt_stat_gamma,
    lognormal_LRT$lrt_stat_lognormal
  ),
  
  p_value = c(
    weibull_LRT$p_val_weibull,
    gamma_LRT$p_val_gamma,
    lognormal_LRT$p_val_lognormal
  )
)

LRT_table$LRT <- round(LRT_table$LRT, 6)
LRT_table$p_value <- round(LRT_table$p_value, 6)

LRT_table$Decision <- ifelse(
  LRT_table$p_value < alpha,
  "Reject H0",
  "Do not reject H0"
)

cat("\n==================================================\n")
cat("             LIKELIHOOD RATIO TESTS\n")
cat("==================================================\n")
cat("Significance level alpha =", alpha, "\n\n")

print(LRT_table)

logL_weibull <- weibull_LRT$logL_reduced_weibull
logL_gamma <- gamma_LRT$logL_reduced_gamma
logL_lognormal <- lognormal_LRT$logL_reduced_lognormal

estimates_reduced_weibull <- weibull_LRT$estimates_reduced_weibull
estimates_reduced_gamma <- gamma_LRT$estimates_reduced_gamma
estimates_reduced_lognormal <- lognormal_LRT$estimates_reduced_lognormal


hess.mat_reduced_weibull <- hessian( log.lik.weibull,estimates_reduced_weibull,y0 = data$y[data$d == 0], y1 = data$y[data$d == 1],y2 = data$y[data$d == 2], y3 = data$y[data$d == 3])
hess.mat_reduced_gamma <- hessian( log.lik.gamma,estimates_reduced_gamma,y0 = data$y[data$d == 0], y1 = data$y[data$d == 1],y2 = data$y[data$d == 2], y3 = data$y[data$d == 3])
hess.mat_reduced_lognormal <- hessian( log.lik.lognormal,estimates_reduced_lognormal,y0 = data$y[data$d == 0], y1 = data$y[data$d == 1],y2 = data$y[data$d == 2], y3 = data$y[data$d == 3])


if (any(!is.finite(hess.mat_reduced_weibull))) {
  stop("Hessian contains non-finite values")
}

FI_weibull <- try(solve(-hess.mat_reduced_weibull), silent = TRUE)

if (inherits(FI_weibull, "try-error") ||
    any(!is.finite(FI_weibull)) ||
    any(diag(FI_weibull) <= 0)) {
  stop("Hessian inversion failed for weibull or variance is not positive")
}



if (any(!is.finite(hess.mat_reduced_gamma))) {
  stop("Hessian contains non-finite values")
}

FI_gamma <- try(solve(-hess.mat_reduced_gamma), silent = TRUE)

if (inherits(FI_gamma, "try-error") ||
    any(!is.finite(FI_gamma)) ||
    any(diag(FI_gamma) <= 0)) {
  stop("Hessian inversion failed for gamma or variance is not positive")
}



if (any(!is.finite(hess.mat_reduced_lognormal))) {
  stop("Hessian contains non-finite values")
}

FI_lognormal <- try(solve(-hess.mat_reduced_lognormal), silent = TRUE)

if (inherits(FI_lognormal, "try-error") ||
    any(!is.finite(FI_lognormal)) ||
    any(diag(FI_lognormal) <= 0)) {
  stop("Hessian inversion failed for lognormal or variance is not positive")
}


SE_weibull_reduced <- sqrt(diag(FI_weibull))

SE_gamma_reduced <- sqrt(diag(FI_gamma))

SE_lognormal_reduced <- sqrt(diag(FI_lognormal))


red_names <- c("c", "lambda1", "sigma1", "lambda2", "sigma2", "lambda3", "sigma3")

parameter_table_weibull <- data.frame(
  Parameter = red_names,
  Estimate = as.numeric(weibull_LRT$estimates_reduced_weibull),
  SE = as.numeric(SE_weibull_reduced)
)

parameter_table_gamma <- data.frame(
  Parameter = red_names,
  Estimate = as.numeric(gamma_LRT$estimates_reduced_gamma),
  SE = as.numeric(SE_gamma_reduced)
)

parameter_table_lognormal <- data.frame(
  Parameter = red_names,
  Estimate = as.numeric(lognormal_LRT$estimates_reduced_lognormal),
  SE = as.numeric(SE_lognormal_reduced)
)
options(scipen = 999)

parameter_table_all <- rbind(
  cbind(Model = "Weibull", parameter_table_weibull),
  cbind(Model = "Gamma", parameter_table_gamma),
  cbind(Model = "Lognormal", parameter_table_lognormal)
)

rownames(parameter_table_all) <- NULL

parameter_table_all$Estimate <- sprintf("%.8f", parameter_table_all$Estimate)
parameter_table_all$SE <- sprintf("%.8f", parameter_table_all$SE)

print(parameter_table_all)

AIC_table <- data.frame(
  Model = c(
    "Lognormal",
    "Weibull",
    "Gamma"
  ),
  
  LogLikelihood = c(
    logL_lognormal,
    logL_weibull,
    logL_gamma
  ),
  
  Number_of_parameters = c(
    7,
    7,
    7
  )
)

AIC_table$AIC <- -2 * AIC_table$LogLikelihood +
  2 * AIC_table$Number_of_parameters

AIC_table$LogLikelihood <- round(
  AIC_table$LogLikelihood, 6
)

AIC_table$AIC <- round(
  AIC_table$AIC, 6
)

AIC_table$BIC <- -2 * AIC_table$LogLikelihood +
  AIC_table$Number_of_parameters * log(n)

AIC_table$BIC <- round(
  AIC_table$BIC, 6
)


cat("\n==================================================\n")
cat("                AIC / BIC COMPARISON\n")
cat("==================================================\n\n")

print(AIC_table)


best_AIC <- AIC_table$Model[
  which.min(AIC_table$AIC)
]

best_BIC <- AIC_table$Model[
  which.min(AIC_table$BIC)
]

cat("\n--------------------------------------------------\n")
cat("Best model according to AIC:", best_AIC, "\n")
cat("Best model according to BIC:", best_BIC, "\n")
cat("--------------------------------------------------\n")

write.csv(
  parameter_table,
  "parameter_estimates_real_data.csv",
  row.names = FALSE
)

write.csv(
  LRT_table,
  "LRT_results_real_data.csv",
  row.names = FALSE
)

write.csv(
  AIC_table,
  "AIC_BIC_results_real_data.csv",
  row.names = FALSE
)


# Cumulative Incidence Function and Overall Survival Function
# Replace these with the estimates from the reduced model (here, Weibull)
c = 0.20450838 
q1 = 1
lambda1 = 0.28289297
sigma1 = 0.96969222
q2 = 1
lambda2 = 0.00251791
sigma2 = 1.38200556
q3 = 1
lambda3 = 0.01455260
sigma3 = 1.05015427

f1 = function(t){
  h1 = haz_fun(t, q1, lambda1, sigma1)
  H1 = cum_haz_fun(t, q1, lambda1, sigma1) 
  H2 = cum_haz_fun(t, q2, lambda2, sigma2)
  H3 = cum_haz_fun(t, q3, lambda3, sigma3) 
  res1 = h1 * exp(-H1 - H2 - H3) 
  return(res1) }

f2 = function(t){ 
  h2 = haz_fun(t, q2, lambda2, sigma2) 
  H1 = cum_haz_fun(t, q1, lambda1, sigma1) 
  H2 = cum_haz_fun(t, q2, lambda2, sigma2)
  H3 = cum_haz_fun(t, q3, lambda3, sigma3)
  res2 = h2 * exp(-H1 - H2 - H3)
  return(res2) }

f22 = function(t){
  h2 = haz_fun(t, q2, lambda2, sigma2) 
  H2 = cum_haz_fun(t, q2, lambda2, sigma2) 
  H3 = cum_haz_fun(t, q3, lambda3, sigma3) 
  res22 = h2 * exp(-H2 - H3) 
  return(res22) }

f3 = function(t){ 
  h3 = haz_fun(t, q3, lambda3, sigma3) 
  H1 = cum_haz_fun(t, q1, lambda1, sigma1) 
  H2 = cum_haz_fun(t, q2, lambda2, sigma2) 
  H3 = cum_haz_fun(t, q3, lambda3, sigma3) 
  res3 = h3 * exp(-H1 - H2 - H3) 
  return(res3) }

f33 = function(t){ 
  h3 = haz_fun(t, q3, lambda3, sigma3) 
  H2 = cum_haz_fun(t, q2, lambda2, sigma2) 
  H3 = cum_haz_fun(t, q3, lambda3, sigma3) 
  res33 = h3 * exp(-H2 - H3) 
  return(res33) }

x = seq(0.1, 8, by = 0.1)

i1 = rep(NA, length(x)) 
i2 = rep(NA, length(x))
i3 = rep(NA, length(x)) 
i22 = rep(NA, length(x)) 
i33 = rep(NA, length(x)) 
for(i in 1:length(x)){ 
  i1[i] = integral(f1, xmin = 0, xmax = x[i])
  i2[i] = integral(f2, xmin = 0, xmax = x[i])
  i3[i] = integral(f3, xmin = 0, xmax = x[i]) 
  i22[i] = integral(f22, xmin = 0, xmax = x[i]) 
  i33[i] = integral(f33, xmin = 0, xmax = x[i]) } 

# CIF 
CIF1 = (1 - c) * i1
CIF2 = c * i22 + (1 - c) * i2 
CIF3 = c * i33 + (1 - c) * i3

plot(
  x, CIF1,
  xlab = "Time (in years)",
  ylab = "Cumulative Probability of Death",
  ylim = c(0, 1),
  yaxt = "n",
  type = "l",
  lty = 1,
  col = "red",
  lwd = 2
)

axis(2, at = seq(0.2, 1.2, by = 0.2))

lines(x, CIF2, lty = 2, col = "green", lwd = 2)
lines(x, CIF3, lty = 3, col = "blue", lwd = 2)

legend(
  "topleft",
  legend = c("Breast Cancer", "Other Cancers", "Other Causes"),
  cex = 0.6,
  lty = c(1, 2, 3),
  col = c("red", "green", "blue"),
  lwd = 2
)

# Overall Survival 
t = seq(0, 8, by = 0.1)

# Cumulative hazards
H1 = cum_haz_fun(t, q1, lambda1, sigma1)
H2 = cum_haz_fun(t, q2, lambda2, sigma2)
H3 = cum_haz_fun(t, q3, lambda3, sigma3)


S = c * exp(-H2 - H3) +
  (1 - c) * exp(-H1 - H2 - H3)

plot(
  t, S,
  xlab = "Time (in years)",
  ylab = "Overall survival probability",
  ylim = c(0, 1),
  type = "l",
  lty = 1,
  col = "black",
  lwd = 2
)
axis(2, at = seq(0.4, 2.2, by = 0.4))





