#  GENERALISED GAMMA MCCR 3 RISK Simulation Study
library(numDeriv)
library(pracma)

upper_incgam_log <- function(x, a) {
  lgamma(a) + pgamma(x, shape = a, scale = 1, lower.tail = FALSE, log.p = TRUE)
}

lower_incgam_log <- function(x, a) {
  lgamma(a) +  pgamma(x, shape = a, scale = 1, lower.tail = TRUE, log.p = TRUE)
}

# generating random variables from gen gamma
r_gen_gamma <- function(n, q, lam, sig){
  
  if(abs(q)< 1e-4){
  
    mu <- -log(lam)
    return(rlnorm(n, meanlog = mu, sdlog = sig))
  }
  
  a <- 1/(lam*(q^(-2*sig/q)))
  d <- (q^(-1))/sig
  p <- abs(q)/sig
  
  y <- rgamma(n, shape = d/p, scale = 1)
  x <- a*(y^(1/p))
  
  return(x)
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

# data generation 
MCCR.data.new3 <- function(n, c, q1, lam1, sig1, q2, lam2, sig2, q3, lam3, sig3, b) {
  
  y <- rep(NA, n)
  d <- rep(NA, n)
 
  U1 <- runif(n, 0, 1)
  cured <- (U1 <= c)
 
  R <- runif(n, min = 1e-6, max = b)
  
  T2_cure =  r_gen_gamma(n,q2,lam2,sig2)
  T3_cure = r_gen_gamma(n,q3,lam3,sig3)
  T_cure = pmin(T2_cure,T3_cure)
  
  T1_unc = r_gen_gamma(n,q1,lam1,sig1)
  T2_unc = r_gen_gamma(n,q2,lam2,sig2)
  T3_unc = r_gen_gamma(n,q3,lam3,sig3)
  T_unc = pmin(T1_unc,T2_unc,T3_unc)
  
  d = apply(cbind(T1_unc,T2_unc,T3_unc),1,which.min) 
  d_cure = apply(cbind(T2_cure,T3_cure),1,which.min) 
  
  T = ifelse(U1<=c,T_cure, T_unc)
  
  for(l in 1:n){
    if(U1[l]<=c){
      d[l] = d_cure[l]+1 # cure case
    }
  }
  
  y = pmin(T,R) # censoring
  d[T>R] = 0    # Cesored case
  return(data.frame(cbind(y,d)))
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

# log-likelihood function gen gamma with q = 0.5
log.lik.gengamma = function(param_gengamma, y0, y1, y2,y3) {
  cure = pmin(pmax(param_gengamma[1], 1e-6), 1 - 1e-6)
  
  q1 = 0.5; lam1 =param_gengamma[2]; sig1 = param_gengamma[3]
  q2 = 0.5; lam2 = param_gengamma[4]; sig2 = param_gengamma[5]
  q3 = 0.5; lam3 = param_gengamma[6]; sig3 = param_gengamma[7]
  
  
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



# EM algorithm with 3 competing risks

EM.MCCR3_global_lognormal = function(data, tol, maxit, c, q_1, lambda1, sigma1, q_2, lambda2, sigma2, q_3, lambda3, sigma3) {
  
  # Split data subsets
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
  
   final_lognorm_flag = (abs(p.new[2,1]) < 1e-4) || (abs(p.new[5,1]) < 1e-4) || (abs(p.new[8,1]) < 1e-4)
  
  active_indices = c(1) 
   active_indices = c(active_indices, 2, 3, 4, 5, 6, 7, 8, 9, 10)
 
  
  p.active = p.new[active_indices, 1]
  
 log_lik_dynamic_profile = function(reduced_p, y0, y1, y2, y3, final_lognorm_flag) {
    full_p = numeric(10)
    full_p[1] = reduced_p[1]
    
    if (final_lognorm_flag) {
      full_p[2]  = 0; full_p[3:4]   = reduced_p[2:3]  # Risk 1
      full_p[5]  = 0; full_p[6:7]   = reduced_p[4:5]  # Risk 2
      full_p[8]  = 0; full_p[9:10]  = reduced_p[6:7]  # Risk 3
    } else {
      full_p[2:10] = reduced_p[2:10]                 
    }
    
    return(log.lik.fn(matrix(full_p, ncol=1), y0=y0, y1=y1, y2=y2, y3=y3))
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
  
  logL_full <- -L_full(param_full_init)
  
  
  
  init_weibull <- c(
    em_estimates[1],  # c (cure fraction)
    em_estimates[3],  # lambda1
    em_estimates[4],  # sigma1
    em_estimates[6],  # lambda2
    em_estimates[7],  # sigma2
    em_estimates[9],  # lambda3
    em_estimates[10]  # sigma3
    
  )
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
    control = list(maxit = 5000, reltol = 1e-12, parscale = abs(init_weibull))
  ), silent = TRUE)
  
   if (inherits(fit_weibull, "try-error")) {
    return(NA)
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

compute_lrt_pvalue_gengamma = function(data, em_estimates) {
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
  
  logL_full <- -L_full(param_full_init)
  
  init_gengamma <- c(
    em_estimates[1],  # c (cure fraction)
    em_estimates[3],  # lambda1
    em_estimates[4],  # sigma1
    em_estimates[6], # lambda2
    em_estimates[7],  # sigma2
    em_estimates[9], # lambda3
    em_estimates[10]  # sigma3
  )
  
  L_GG = function(p) {
    if (p[1] <= 0 || p[1] >= 1 || any(p[2:7] <= 0)) return(1e10)
    val <- -log.lik.gengamma(p, y0, y1, y2, y3)
    if (is.na(val) || is.nan(val) || is.infinite(val)) return(1e10)
    return(val)
  }
  
  fit_gengamma <- try(optim(
    par = init_gengamma, 
    fn = L_GG, 
    method = "Nelder-Mead",
    control = list(maxit = 5000, reltol = 1e-12, parscale = abs(init_gengamma))
  ), silent = TRUE)
  
  if (inherits(fit_gengamma, "try-error")) {
    return(NA)
  } 
  logL_reduced_gengamma <- -fit_gengamma$value
  estimates_reduced_gengamma <- fit_gengamma$par
  
  lrt_stat_gengamma <- 2 * (logL_full - logL_reduced_gengamma)
  p_val_gengamma    <- pchisq(lrt_stat_gengamma, df = 3, lower.tail = FALSE)
  return(list(
    lrt_stat_gengamma = lrt_stat_gengamma,
    p_val_gengamma = p_val_gengamma,
    logL_full = logL_full,
    logL_reduced_gengamma = logL_reduced_gengamma,
    estimates_reduced_gengamma = estimates_reduced_gengamma
  )) }

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
  
  logL_full <- -L_full(param_full_init)
  
  init_gamma <- c(
    em_estimates[1],  # c (cure fraction)
    em_estimates[3],  # lambda1
    em_estimates[4],  # sigma1
    em_estimates[6],  # lambda2
    em_estimates[7],  # sigma2
    em_estimates[9],  # lambda3
    em_estimates[10]  # sigma3
    
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
    control = list(maxit = 5000, reltol = 1e-12, parscale = abs(init_gamma))
  ), silent = TRUE)
  
  if (inherits(fit_gamma, "try-error")) {
    return(NA) 
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
  logL_full <- -L_full(param_full_init)
  
  
  
  init_lognormal <- c(
    em_estimates[1],  # c (cure fraction)
    em_estimates[3],  # lambda1
    em_estimates[4],  # sigma1
    em_estimates[6],  # lambda2
    em_estimates[7],  # sigma2
    em_estimates[9],  # lambda3
    em_estimates[10]  # sigma3
    
  )
  
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
    control = list(maxit = 5000, reltol = 1e-12, parscale = abs(init_lognormal))
  ), silent = TRUE)
  
  if (inherits(fit_lognormal, "try-error")) {
    return(NA) 
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

# calling the functions

N = 500
n = 800
alpha = 0.05

c.true    = 0.3
q1.true   = 0
lam1.true = 1.8
sig1.true = 0.6
q2.true   = 0
lam2.true = 1.7
sig2.true = 0.7
q3.true   = 0
lam3.true = 2.1
sig3.true = 0.7

b = 12
incr = 0.2

c.hat = rep(NA,N)
q1.hat = rep(NA,N)
lam1.hat = rep(NA,N)
sig1.hat = rep(NA,N)
q2.hat = rep(NA,N)
lam2.hat = rep(NA,N)
sig2.hat = rep(NA,N)
q3.hat = rep(NA,N)
lam3.hat = rep(NA,N)
sig3.hat = rep(NA,N)

std.c = rep(NA,N)
std.q1 = rep(NA,N)
std.lam1 = rep(NA,N)
std.sig1 = rep(NA,N)
std.q2 = rep(NA,N)
std.lam2 = rep(NA,N)
std.sig2 = rep(NA,N)
std.q3 = rep(NA,N)
std.lam3 = rep(NA,N)
std.sig3 = rep(NA,N)

c.temp = rep(NA,N)
q1.temp = rep(NA,N)
lam1.temp = rep(NA,N)
sig1.temp = rep(NA,N)
q2.temp = rep(NA,N)
lam2.temp = rep(NA,N)
sig2.temp = rep(NA,N)
q3.temp = rep(NA,N)
lam3.temp = rep(NA,N)
sig3.temp = rep(NA,N)

lcl_c_95=rep(NA,N)
lcl_q1_95=rep(NA,N)
lcl_lam1_95=rep(NA,N)
lcl_sig1_95=rep(NA,N)
lcl_q2_95=rep(NA,N)
lcl_lam2_95=rep(NA,N)
lcl_sig2_95=rep(NA,N)
lcl_q3_95=rep(NA,N)
lcl_lam3_95=rep(NA,N)
lcl_sig3_95=rep(NA,N)

ucl_c_95=rep(NA,N)
ucl_q1_95=rep(NA,N)
ucl_lam1_95=rep(NA,N)
ucl_sig1_95=rep(NA,N)
ucl_q2_95=rep(NA,N)
ucl_lam2_95=rep(NA,N)
ucl_sig2_95=rep(NA,N)
ucl_q3_95=rep(NA,N)
ucl_lam3_95=rep(NA,N)
ucl_sig3_95=rep(NA,N)

t_c_95=0
t_q1_95=0
t_lam1_95=0
t_sig1_95=0
t_q2_95=0
t_lam2_95=0
t_sig2_95=0
t_q3_95=0
t_lam3_95=0
t_sig3_95=0

rejections_weibull <- 0
valid_runs_weibull <- 0
p_values_collected_weibull <- c()

rejections_gamma <- 0
valid_runs_gamma <- 0
p_values_collected_gamma <- c()

rejections_gengamma <- 0
valid_runs_gengamma <- 0
p_values_collected_gengamma <- c()


rejections_lognormal <- 0
valid_runs_lognormal <- 0
p_values_collected_lognormal <- c()



aic_winners   <- rep(NA, N)
bic_winners   <- rep(NA, N)

j = 1

count = 0

while(j <= N){
  print(j)
  data = MCCR.data.new3(n=n,c=c.true,q1=q1.true,lam1=lam1.true,sig1=sig1.true,q2=q2.true,lam2=lam2.true,sig2=sig2.true,q3=q3.true,lam3=lam3.true,sig3=sig3.true,b=b)
  
  summary(data$y)
  addmargins(table(data$d,useNA="ifany"))
  
  
  
  c.init = sample(seq((c.true-(incr*abs(c.true))),(c.true+(incr*abs(c.true))),by=0.01),1)
  q1.init = sample(seq((q1.true-(incr*abs(q1.true))),(q1.true+(incr*abs(q1.true))),by=0.01),1)
  lam1.init = sample(seq((lam1.true-(incr*abs(lam1.true))),(lam1.true+(incr*abs(lam1.true))),by=0.01),1)
  sig1.init = sample(seq((sig1.true-(incr*abs(sig1.true))),(sig1.true+(incr*abs(sig1.true))),by=0.01),1)
  
  q2.init = sample(seq((q2.true-(incr*abs(q2.true))),(q2.true+(incr*abs(q2.true))),by=0.01),1)
  lam2.init = sample(seq((lam2.true-(incr*abs(lam2.true))),(lam2.true+(incr*abs(lam2.true))),by=0.01),1)
  sig2.init = sample(seq((sig2.true-(incr*abs(sig2.true))),(sig2.true+(incr*abs(sig2.true))),by=0.01),1)
  
  q3.init = sample(seq((q3.true-(incr*abs(q3.true))),(q3.true+(incr*abs(q3.true))),by=0.01),1)
  lam3.init = sample(seq((lam3.true-(incr*abs(lam3.true))),(lam3.true+(incr*abs(lam3.true))),by=0.01),1)
  sig3.init = sample(seq((sig3.true-(incr*abs(sig3.true))),(sig3.true+(incr*abs(sig3.true))),by=0.01),1)
  
  
  EM.result = EM.MCCR3_global_lognormal(
    data=data,
    tol=0.001,
    maxit=500,
    c=c.init,
    q_1=q1.init,
    lambda1=lam1.init,
    sigma1=sig1.init,
    q_2=q2.init,
    lambda2=lam2.init,
    sigma2=sig2.init,
    q_3=q3.init,
    lambda3=lam3.init,
    sigma3=sig3.init
  )
  
  EM <- EM.result$estimates
  if(all(EM==0)){
    j = j + 0
    count=count+1
  }else{
    valid_runs_weibull <- valid_runs_weibull + 1
    valid_runs_gamma <- valid_runs_gamma + 1
    valid_runs_gengamma <- valid_runs_gengamma + 1
    valid_runs_lognormal <- valid_runs_lognormal + 1
    
    lrt_weibull <- compute_lrt_pvalue_weibull(data, EM[1:10])
    
    p_val_weibull <- lrt_weibull$p_val_weibull
    p_values_collected_weibull <- c(
      p_values_collected_weibull,
      p_val_weibull
    )
    
    lrt_lognormal <- compute_lrt_pvalue_lognormal(data, EM[1:10])
    
    p_val_lognormal <- lrt_lognormal$p_val_lognormal
    p_values_collected_lognormal <- c(
      p_values_collected_lognormal,
      p_val_lognormal
    )
    
    
    p_val_gamma <- compute_lrt_pvalue_gamma(data, EM[1:10])$p_val_gamma
    p_values_collected_gamma <- c(p_values_collected_gamma, p_val_gamma)
    
    p_val_gengamma <- compute_lrt_pvalue_gengamma(data, EM[1:10])$p_val_gengamma
    p_values_collected_gengamma <- c(p_values_collected_gengamma, p_val_gengamma)
   
       if(p_val_weibull < alpha) {
      rejections_weibull <- rejections_weibull + 1
    }
    if(p_val_gamma < alpha) {
      rejections_gamma <- rejections_gamma + 1
    }
    if(p_val_gengamma < alpha) {
      rejections_gengamma <- rejections_gengamma + 1
    }
    if(p_val_lognormal < alpha) {
      rejections_lognormal <- rejections_lognormal + 1
    }
    
    
   
    reject_lognormal <- !is.na(p_val_lognormal) && p_val_lognormal < alpha 
    reject_weibull <- !is.na(p_val_weibull) && p_val_weibull < alpha 
    
  
    
    logL_reduced_weibull <- compute_lrt_pvalue_weibull(data, EM[1:10])$logL_reduced_weibull
    logL_reduced_lognormal <- compute_lrt_pvalue_lognormal(data, EM[1:10])$logL_reduced_lognormal
    logL_reduced_gamma <- compute_lrt_pvalue_gamma(data, EM[1:10])$logL_reduced_gamma
    logL_reduced_gengamma <- compute_lrt_pvalue_gengamma(data, EM[1:10])$logL_reduced_gengamma
    
    aics <- c(Lognormal = -2 * logL_reduced_lognormal + 2 * 7,
              Weibull   = -2 * logL_reduced_weibull + 2 * 7,
              Gamma     = -2 * logL_reduced_gamma + 2 * 7,
              GenGamma  = -2 * logL_reduced_gengamma + 2 * 7)
    
    bics <- c(Lognormal = -2 * logL_reduced_lognormal + 7 * log(800),
              Weibull   = -2 * logL_reduced_weibull + 7 * log(800),
              Gamma     = -2 * logL_reduced_gamma + 7 * log(800),
              GenGamma  = -2 * logL_reduced_gengamma + 7 * log(800))
    
    aic_winners[j] <- names(which.min(aics))
    bic_winners[j] <- names(which.min(bics))
    
    
    
    c.hat[j] = EM[1]
    q1.hat[j] = EM[2]
    lam1.hat[j] = EM[3]
    sig1.hat[j] = EM[4]
    q2.hat[j] = EM[5]
    lam2.hat[j] = EM[6]
    sig2.hat[j] = EM[7]
    q3.hat[j] = EM[8]
    lam3.hat[j] = EM[9]
    sig3.hat[j] = EM[10]
    
    std.c[j] = EM[11]
    std.q1[j] = EM[12]
    std.lam1[j] = EM[13]
    std.sig1[j] = EM[14]
    std.q2[j] = EM[15]
    std.lam2[j] = EM[16]
    std.sig2[j] = EM[17]
    std.q3[j] = EM[18]
    std.lam3[j] = EM[19]
    std.sig3[j] = EM[20]
    
    c.temp[j] = c.hat[j] - c.true
    q1.temp[j] = q1.hat[j] - q1.true
    lam1.temp[j] = lam1.hat[j] - lam1.true
    sig1.temp[j] = sig1.hat[j] - sig1.true
    q2.temp[j] = q2.hat[j] - q2.true
    lam2.temp[j] = lam2.hat[j] - lam2.true
    sig2.temp[j] = sig2.hat[j] - sig2.true
    q3.temp[j] = q3.hat[j] - q3.true
    lam3.temp[j] = lam3.hat[j] - lam3.true
    sig3.temp[j] = sig3.hat[j] - sig3.true
    
    lcl_c_95[j] = c.hat[j] - (1.96*std.c[j])
    lcl_q1_95[j] = q1.hat[j] - (1.96*std.q1[j])
    lcl_lam1_95[j] = lam1.hat[j] - (1.96*std.lam1[j])
    lcl_sig1_95[j] = sig1.hat[j] - (1.96*std.sig1[j])
    lcl_q2_95[j] = q2.hat[j] - (1.96*std.q2[j])
    lcl_lam2_95[j] = lam2.hat[j] - (1.96*std.lam2[j])
    lcl_sig2_95[j] = sig2.hat[j] - (1.96*std.sig2[j])
    lcl_q3_95[j] = q3.hat[j] - (1.96*std.q3[j])
    lcl_lam3_95[j] = lam3.hat[j] - (1.96*std.lam3[j])
    lcl_sig3_95[j] = sig3.hat[j] - (1.96*std.sig3[j])
    
    ucl_c_95[j] = c.hat[j] + (1.96*std.c[j])
    ucl_q1_95[j] = q1.hat[j] + (1.96*std.q1[j])
    ucl_lam1_95[j] = lam1.hat[j] + (1.96*std.lam1[j])
    ucl_sig1_95[j] = sig1.hat[j] + (1.96*std.sig1[j])
    ucl_q2_95[j] = q2.hat[j] + (1.96*std.q2[j])
    ucl_lam2_95[j] = lam2.hat[j] + (1.96*std.lam2[j])
    ucl_sig2_95[j] = sig2.hat[j] + (1.96*std.sig2[j])
    ucl_q3_95[j] = q3.hat[j] + (1.96*std.q3[j])
    ucl_lam3_95[j] = lam3.hat[j] + (1.96*std.lam3[j])
    ucl_sig3_95[j] = sig3.hat[j] + (1.96*std.sig3[j])
    
    
    
    if(c.true > lcl_c_95[j] & c.true < ucl_c_95[j]){
      t_c_95 = t_c_95 + 1
    }
    if(q1.true > lcl_q1_95[j] & q1.true < ucl_q1_95[j]){
      t_q1_95 = t_q1_95 + 1
    }
    
    if(lam1.true > lcl_lam1_95[j] & lam1.true < ucl_lam1_95[j]){
      t_lam1_95 = t_lam1_95 + 1
    }
    if(sig1.true > lcl_sig1_95[j] & sig1.true < ucl_sig1_95[j]){
      t_sig1_95 = t_sig1_95 + 1
    }
    if(q2.true > lcl_q2_95[j] & q2.true < ucl_q2_95[j]){
      t_q2_95 = t_q2_95 + 1
    }
    
    if(lam2.true > lcl_lam2_95[j] & lam2.true < ucl_lam2_95[j]){
      t_lam2_95 = t_lam2_95 + 1
    }
    if(sig2.true > lcl_sig2_95[j] & sig2.true < ucl_sig2_95[j]){
      t_sig2_95 = t_sig2_95 + 1
    }
    if(q3.true > lcl_q3_95[j] & q3.true < ucl_q3_95[j]){
      t_q3_95 = t_q3_95 + 1
    }
    
    if(lam3.true > lcl_lam3_95[j] & lam3.true < ucl_lam3_95[j]){
      t_lam3_95 = t_lam3_95 + 1
    }
    if(sig3.true > lcl_sig3_95[j] & sig3.true < ucl_sig3_95[j]){
      t_sig3_95 = t_sig3_95 + 1
    }
    
    j = j+1
    
  }#end of else
}#end of while

print("count")
print(count)

avg.c = sum(c.hat)/N
avg.q1 = sum(q1.hat)/N
avg.lam1 = sum(lam1.hat)/N
avg.sig1 = sum(sig1.hat)/N
avg.q2 = sum(q2.hat)/N
avg.lam2 = sum(lam2.hat)/N
avg.sig2 = sum(sig2.hat)/N
avg.q3 = sum(q3.hat)/N
avg.lam3 = sum(lam3.hat)/N
avg.sig3 = sum(sig3.hat)/N

avg.std.c = sum(std.c)/N
avg.std.q1 = sum(std.q1)/N
avg.std.lam1 = sum(std.lam1)/N
avg.std.sig1 = sum(std.sig1)/N
avg.std.q2 = sum(std.q2)/N
avg.std.lam2 = sum(std.lam2)/N
avg.std.sig2 = sum(std.sig2)/N
avg.std.q3 = sum(std.q3)/N
avg.std.lam3 = sum(std.lam3)/N
avg.std.sig3 = sum(std.sig3)/N

bias.c = sum(c.temp)/N
bias.q1 = sum(q1.temp)/N
bias.lam1 = sum(lam1.temp)/N
bias.sig1 = sum(sig1.temp)/N
bias.q2 = sum(q2.temp)/N
bias.lam2 = sum(lam2.temp)/N
bias.sig2 = sum(sig2.temp)/N
bias.q3 = sum(q3.temp)/N
bias.lam3 = sum(lam3.temp)/N
bias.sig3 = sum(sig3.temp)/N

rmse.c = sqrt(sum(c.temp^2)/(N-1))
rmse.q1 = sqrt(sum(q1.temp^2)/(N-1))
rmse.lam1 = sqrt(sum(lam1.temp^2)/(N-1))
rmse.sig1 = sqrt(sum(sig1.temp^2)/(N-1))
rmse.q2 = sqrt(sum(q2.temp^2)/(N-1))
rmse.lam2 = sqrt(sum(lam2.temp^2)/(N-1))
rmse.sig2 = sqrt(sum(sig2.temp^2)/(N-1))
rmse.q3 = sqrt(sum(q3.temp^2)/(N-1))
rmse.lam3 = sqrt(sum(lam3.temp^2)/(N-1))
rmse.sig3 = sqrt(sum(sig3.temp^2)/(N-1))


cp_c_95 = t_c_95/N
cp_q1_95 = t_q1_95/N
cp_lam1_95 = t_lam1_95/N
cp_sig1_95 = t_sig1_95/N
cp_q2_95 = t_q2_95/N
cp_lam2_95 = t_lam2_95/N
cp_sig2_95 = t_sig2_95/N
cp_q3_95 = t_q3_95/N
cp_lam3_95 = t_lam3_95/N
cp_sig3_95 = t_sig3_95/N


avg.c
avg.q1
avg.lam1
avg.sig1
avg.q2
avg.lam2
avg.sig2
avg.q3
avg.lam3
avg.sig3

bias.c
bias.q1
bias.lam1
bias.sig1
bias.q2
bias.lam2
bias.sig2
bias.q3
bias.lam3
bias.sig3

rmse.c
rmse.q1
rmse.lam1
rmse.sig1
rmse.q2
rmse.lam2
rmse.sig2
rmse.q3
rmse.lam3
rmse.sig3

avg.std.c
avg.std.q1
avg.std.lam1
avg.std.sig1
avg.std.q2
avg.std.lam2
avg.std.sig2
avg.std.q3
avg.std.lam3
avg.std.sig3

cp_c_95 
cp_q1_95 
cp_lam1_95 
cp_sig1_95
cp_q2_95 
cp_lam2_95 
cp_sig2_95
cp_q3_95 
cp_lam3_95 
cp_sig3_95

par.est = round(c(avg.c,avg.q1,avg.lam1,avg.sig1,avg.q2,avg.lam2,avg.sig2,avg.q3,avg.lam3,avg.sig3),4)
par.bias = round(c(bias.c,bias.q1,bias.lam1,bias.sig1,bias.q2,bias.lam2,bias.sig2,bias.q3,bias.lam3,bias.sig3),4)
par.rmse = round(c(rmse.c,rmse.q1,rmse.lam1,rmse.sig1,rmse.q2,rmse.lam2,rmse.sig2,rmse.q3,rmse.lam3,rmse.sig3),4)
par.se = round(c(avg.std.c,avg.std.q1,avg.std.lam1,avg.std.sig1,avg.std.q2,avg.std.lam2,avg.std.sig2,avg.std.q3,avg.std.lam3,avg.std.sig3),4)
par.cp.95 = round(c(cp_c_95,cp_q1_95,cp_lam1_95,cp_sig1_95,cp_q2_95,cp_lam2_95,cp_sig2_95,cp_q3_95,cp_lam3_95,cp_sig3_95),4)

data.frame(par.est,par.bias,par.rmse,par.se,par.cp.95)


rejection_rate_weibull <- rejections_weibull / valid_runs_weibull
rejection_rate_gamma <- rejections_gamma / valid_runs_gamma
rejection_rate_gengamma <- rejections_gengamma / valid_runs_gengamma
rejection_rate_lognormal <- rejections_lognormal / valid_runs_lognormal



aic_selection_table    <- table(aic_winners) / N
bic_selection_table    <- table(bic_winners) / N

cat("\n==================================================\n")
cat("          SIMULATION REPORT: LRT FITTED = WEIBULL       \n")
cat("==================================================\n")
cat("Total Intended Runs (N):          ", N, "\n")
cat("Successfully Converged Runs:      ", valid_runs_weibull, "\n")
cat("Nominal Significance Level (Alpha):", alpha, "\n")
cat("Total Null Rejections Count:      ", rejections_weibull, "\n")
cat("--------------------------------------------------\n")
cat("Empirical Rejection Rate:         ", round(rejection_rate_weibull, 4), " (", rejection_rate_weibull * 100, "%)\n", sep="")

cat("\n==================================================\n")
cat("          SIMULATION REPORT: LRT FITTED = GAMMA       \n")
cat("==================================================\n")
cat("Total Intended Runs (N):          ", N, "\n")
cat("Successfully Converged Runs:      ", valid_runs_gamma, "\n")
cat("Nominal Significance Level (Alpha):", alpha, "\n")
cat("Total Null Rejections Count:      ", rejections_gamma, "\n")
cat("--------------------------------------------------\n")
cat("Empirical Rejection Rate:         ", round(rejection_rate_gamma, 4), " (", rejection_rate_gamma * 100, "%)\n", sep="")

cat("\n==================================================\n")
cat("          SIMULATION REPORT: LRT FITTED = GENGAMMA       \n")
cat("==================================================\n")
cat("Total Intended Runs (N):          ", N, "\n")
cat("Successfully Converged Runs:      ", valid_runs_gengamma, "\n")
cat("Nominal Significance Level (Alpha):", alpha, "\n")
cat("Total Null Rejections Count:      ", rejections_gengamma, "\n")
cat("--------------------------------------------------\n")
cat("Empirical Rejection Rate:         ", round(rejection_rate_gengamma, 4), " (", rejection_rate_gengamma * 100, "%)\n", sep="")

cat("\n==================================================\n")
cat("          SIMULATION REPORT: LRT FITTED = LOGNORMAL       \n")
cat("==================================================\n")
cat("Total Intended Runs (N):          ", N, "\n")
cat("Successfully Converged Runs:      ", valid_runs_lognormal, "\n")
cat("Nominal Significance Level (Alpha):", alpha, "\n")
cat("Total Null Rejections Count:      ", rejections_lognormal, "\n")
cat("--------------------------------------------------\n")
cat("Empirical Rejection Rate:         ", round(rejection_rate_lognormal, 4), " (", rejection_rate_lognormal * 100, "%)\n", sep="")



cat("\n==================================================\n")
cat("      SIMULATION REPORT: MODEL SELECTION RATES      \n")
cat("==================================================\n")
cat("Total Simulation Runs (N):         ", N, "\n")
cat("Data Generating Distribution:       Lognormal (q=0)\n")
cat("--------------------------------------------------\n\n")

cat("--- AIC SELECTION PROPORTIONS ---\n")
print(aic_selection_table)
cat("\n")

cat("--- BIC SELECTION PROPORTIONS ---\n")
print(bic_selection_table)
cat("==================================================\n")
