# Full
library("lamW")
library(nleqslv)
library(GenSA)
library(GA)

rm(list = ls())

alpha = 5
beta = -20
gamma = .5
nu = exp(-gamma)
# nu = 0.5
delta= 25
n = 10000
alpha_mme_vec <- c()
beta_mme_vec <- c()
beta_mle_vec <- c()
gamma_mme_vec <- c()
gamma_mle_vec <- c()
nu_mme_vec <- c()
delta_mme_vec <- c()
delta_mle_vec <- c()
rho_mme <- c()
rho_mle <- c()
pearson_corr_vec <- c()

solve_nu_safe <- function(alpha_mme, mu_mme, S12, S2, M2, z_start = 0) {
  rhs <- S12^2 / (S2 - M2)
  
  f_z <- function(z) {
    nu <- 1 / (1 + exp(-z))  # maps z ∈ ℝ to nu ∈ (0,1)
    num <- alpha_mme^2 * (nu - 1)^2
    exp_term <- mu_mme^((nu - 1)^2)
    denom <- exp_term - 1
    
    # if (!is.finite(denom) || denom == 0) return(Inf)
    num / denom - rhs
  }
  
  sol <- nleqslv(x = z_start, fn = f_z)
  
  if (sol$termcd != 1) warning("Solution may not have converged")
  
  # Map back to ν ∈ (0,1)
  nu <- 1 / (1 + exp(-sol$x))
  return(nu)
}

solve_nu_optim <- function(alpha_mme, mu_mme, S12, S2, M2, nu_start = nu) {
  rhs <- S12^2 / (S2 - M2)

  obj <- function(nu) {
    if (nu <= 0 || nu >= 1) return(1e10)  # enforce bounds
    num <- alpha_mme^2* (nu - 1)^2
    denom <- mu_mme^((nu - 1)^2) - 1
    if (!is.finite(denom) || denom == 0) return(1e10)
    (num / denom - rhs)^2
  }

  res <- optim(par = nu_start, fn = obj, method = "L-BFGS-B", lower = 1e-6, upper = 0.9999)

  if (res$convergence != 0) warning("Solution may not have converged")

  return(res$par)
}

solve_nu_trans <- function(alpha_mme, mu_mme, S12, S2, M2, trans_start = (nu-1)^2) {
  rhs <- S12^2 / (S2 - M2)

  obj <- function(trans) {
    # if (nu <= 0 || nu >= 1) return(1e10)  # enforce bounds
    num <- alpha_mme^2* trans
    denom <- mu_mme^trans - 1
    # if (!is.finite(denom) || denom == 0) return(1e10)
    (num / denom - rhs)^2
  }

  res <- optim(par = trans_start, fn = obj, method = "L-BFGS-B", lower = 1e-6, upper = 0.9999)

  if (res$convergence != 0) warning("Solution may not have converged")

  return(res$par)
}


# solve_nu_trans <- function(alpha_mme, mu_mme, S12, S2, M2, trans_start = (nu - 1)^2) {
#   rhs <- S12^2 / (S2 - M2)
#   obj <- function(trans) {
#     num <- alpha_mme^2 * trans
#     denom <- mu_mme^trans - 1
#     return(-(num / denom - rhs)^2)  # explicit return stops printing
#   }
#   
#   
#   # res <- suppressWarnings(GenSA(
#   #   par    = trans_start,
#   #   fn     = obj,
#   #   lower  = 0,
#   #   upper  = 1,
#   #   control = list(max.call = 5000, verbose = FALSE) # you can increase if needed
#   # ))
#   res = ga(type = "real-valued",  obj, lower = 0, upper = 1,popSize = 100,maxiter = 25,keepBest = TRUE, monitor = FALSE)
# 
# 
#   return( res@solution[1,])
# }



mle_solve = function(x1, x2, n, alpha_mme){
  obj = function(pars){ 
    beta  = pars[1]
    eta   = pars[2]
    zeta  = pars[3]
    
    gamma = exp(eta)             # ensures gamma >= 0
    delta = -beta + exp(zeta)    # ensures delta >= -beta
    nu = exp(-gamma)
    
    loglik = -n*(alpha_mme+beta+delta) + log(alpha_mme)*sum(x1) + beta*sum(nu^x1) + 
      sum(x2*log(delta + beta*(1 - nu^x1)))
    return(-loglik)  # optim minimizes
  }
  
  res <- optim(
    # par = c(beta, log(gamma), log(delta+beta)),       # starting values (beta, eta, zeta)
    par = c(-1, 0, 0),       # starting values (beta, eta, zeta)
    fn  = obj,
    method = "L-BFGS-B",
    lower = c(-50, -10, -10),  # bounds on (beta, eta, zeta)
    upper = c(0,   10,  10)
  )
  
  # Transform back to (beta, gamma, delta)
  beta  <- res$par[1]
  eta   <- res$par[2]
  zeta  <- res$par[3]
  gamma <- exp(eta)
  delta <- -beta + exp(zeta)
  return (c(beta, gamma, delta))
}



seed = 1
i = 1
runs = 1e4
while (length(gamma_mme_vec) < runs){
  set.seed(seed)
  x1 <- rpois(n = n, lambda = alpha)
  x2 <- rpois(n =n, lambda = (delta+ beta*(1 - exp(-gamma*x1))))
  # x2 <- rpois(n =n, lambda = (delta+ beta*(1 - nu^x1)))
  M1 <- sum(x1)/n
  M2<- sum(x2)/n
  S12 <- cov(x1, x2)*(n-1)/n
  S2 <- var(x2) *(n-1)/n
  if (S2 == M2){
    seed<- seed+1
    next
  }
  if ( M1^2/(exp(M1) - 1) > S12^2/(S2-M2)){
    seed<- seed+1
    next
  }
  if ( M1 < S12^2/(S2-M2)){
    seed<- seed+1
    next
  }

  

  
  
  alpha_mme <- M1
  mu_mme = exp(alpha_mme)
  # nu_mme = -sqrt(solve_nu_trans(alpha_mme, mu_mme, S12, S2, M2))+1
  nu_mme = -sqrt(solve_nu_trans(alpha_mme, mu_mme, S12, S2, M2))+1
    
  gamma_mme = -log(nu_mme)
  beta_mme = S12/(alpha_mme*(1 - nu_mme)*mu_mme^(nu_mme-1))
  delta_mme = M2- beta_mme*(1 - mu_mme^(nu_mme - 1))

  
  # print (c(S12^2/(S2- M2), alpha_mme^2*(nu_mme - 1)^2/(mu_mme^((nu_mme-1)^2) - 1)))
  
  
  alpha_mme_vec[i] <- alpha_mme
  beta_mme_vec[i] <- beta_mme
  gamma_mme_vec[i] <- gamma_mme
  nu_mme_vec[i] <- nu_mme
  delta_mme_vec[i] <- delta_mme
  
  num =  sqrt(alpha_mme)*beta_mme*(1- nu_mme)*mu_mme^(nu_mme -1)
  den = sqrt( beta_mme^2*(mu_mme^(nu_mme^2 - 1) - mu_mme^(2*nu_mme -2)) + beta_mme*(1 - mu_mme^(nu_mme - 1)) +delta_mme)
  rho_mme[i] = num/den
  
  #MLE 
  mles = mle_solve(x1, x2, n = n, alpha_mme)
  beta_mle = mles[1]
  gamma_mle = mles[2]
  delta_mle = mles[3]
  beta_mle_vec[i] <- beta_mle
  gamma_mle_vec[i] <- gamma_mle
  nu_mle = exp(-gamma_mle)
  delta_mle_vec[i] <- delta_mle
  
  num =  sqrt(alpha_mme)*beta_mle*(1- nu_mle)*mu_mme^(nu_mle -1)
  den = sqrt( beta_mle^2*(mu_mme^(nu_mle^2 - 1) - mu_mme^(2*nu_mle -2)) + beta_mle*(1 - mu_mme^(nu_mle - 1)) +delta_mle)
  rho_mle[i] = num/den
  
  
  pearson_corr_vec[i] <- cor(x= x1, y= x2, method = "pearson")
  i<- i+1
  seed = seed+1
}



remove_outliers_iqr <- function(x, coef = 1.5, na.rm = TRUE) {
  q <- quantile(x, probs = c(0.25, 0.75), na.rm = na.rm)
  iqr <- q[2] - q[1]
  lo <- q[1] - coef * iqr
  hi <- q[2] + coef * iqr
  x[x >= lo & x <= hi]
}


sim_func= function(mme, mle = mme){
  mme = remove_outliers_iqr(mme)
  mle = remove_outliers_iqr(mle)
  c(summary(mme)[3], summary(mle)[3], sd(mme), sd(mle), summary(mme)[3] - qnorm(0.975)*sd(mme),
    summary(mme)[3] + qnorm(0.975)*sd(mme),  summary(mle)[3] - qnorm(0.975)*sd(mle),  
    summary(mle)[3] + qnorm(0.975)*sd(mle))
}

sim_func(alpha_mme_vec)
# sim_func(beta_mme_vec, beta_mle_vec)
sim_func(gamma_mme_vec, gamma_mle_vec)
sim_func(delta_mme_vec, delta_mle_vec)
sim_func(rho_mme, rho_mle)



sqrt(var(alpha_mme_vec))
hist(beta_mme_vec, breaks = 1000)
summary(beta_mme_vec)
hist(gamma_mme_vec, breaks =100)
summary(gamma_mme_vec)
hist(delta_mme_vec, breaks = 100)
summary(delta_mme_vec)

hist(beta_mle_vec, breaks = 1000)
summary(beta_mle_vec)
hist(gamma_mle_vec, breaks =100)
summary(gamma_mle_vec)
hist(delta_mle_vec, breaks = 100)
summary(delta_mle_vec)


summary(rho_mme)
summary(pearson_corr_vec)






















# Case 1: beta  = -1

library("lamW")
library(nleqslv)
library(GA)

rm(list = ls())

alpha = 1
beta = -1
gamma = 0.3
nu = exp(-gamma)
delta= 25
n = 10000
alpha_mme_vec <- c()
gamma_mme_vec <- c()
gamma_mle_vec <- c()
nu_mme_vec <- c()
delta_mme_vec <- c()
delta_mle_vec <- c()
rho_mme <- c()
rho_mle <- c()
pearson_corr_vec <- c()

solve_nu_safe <- function(alpha_mme, mu_mme, S12, S2, M2, z_start = 0) {
  rhs <- S12^2 / (S2 - M2)
  
  f_z <- function(z) {
    nu <- 1 / (1 + exp(-z))  # maps z ∈ ℝ to nu ∈ (0,1)
    num <- alpha_mme^2 * (nu - 1)^2
    exp_term <- mu_mme^((nu - 1)^2)
    denom <- exp_term - 1
    
    # if (!is.finite(denom) || denom == 0) return(Inf)
    num / denom - rhs
  }
  
  sol <- nleqslv(x = z_start, fn = f_z)
  
  if (sol$termcd != 1) warning("Solution may not have converged")
  
  # Map back to ν ∈ (0,1)
  nu <- 1 / (1 + exp(-sol$x))
  return(nu)
}

solve_nu_optim <- function(alpha_mme, mu_mme, S12, S2, M2, nu_start = nu) {
  rhs <- S12^2 / (S2 - M2)
  
  obj <- function(nu) {
    if (nu <= 0 || nu >= 1) return(1e10)  # enforce bounds
    num <- alpha_mme^2* (nu - 1)^2
    denom <- mu_mme^((nu - 1)^2) - 1
    if (!is.finite(denom) || denom == 0) return(1e10)
    (num / denom - rhs)^2
  }
  
  res <- optim(par = nu_start, fn = obj, method = "L-BFGS-B", lower = 1e-6, upper = 0.9999)
  
  if (res$convergence != 0) warning("Solution may not have converged")
  
  return(res$par)
}

mle_solve = function(x1, x2, n, alpha_mme){
  obj = function(pars){ 
    eta   = pars[1]
    zeta  = pars[2]
    
    gamma = exp(eta)             # ensures gamma >= 0
    delta = -beta + exp(zeta)    # ensures delta >= -beta
    nu = exp(-gamma)
    
    loglik = -n*(alpha_mme+beta+delta) + log(alpha_mme)*sum(x1) + beta*sum(nu^x1) + 
      sum(x2*log(delta + beta*(1 - nu^x1)))
    return(loglik)  # optim minimizes
  }
  
  # res <- optim(
  #   # par = c(log(gamma), log(delta+beta)),       # starting values (beta, eta, zeta)
  #   par = c(-.1, .1),       # starting values (beta, eta, zeta)
  #   fn  = obj,
  #   method = "L-BFGS-B",
  #   lower = c(-Inf, -Inf),  # bounds on (beta, eta, zeta)
  #   upper = c( Inf,  Inf)
  # )
  res = ga(type = "real-valued",  obj, lower = c(-10, -10), upper = c(10, 10),popSize = 100,maxiter = 25,keepBest = TRUE, monitor = FALSE)
  
  # return( res@solution[1,])
  eta = res@solution[1,1]
  zeta = res@solution[1,2]
  # print (eta)
  gamma <- exp(eta)
  delta <- -beta + exp(zeta)
  return  (c(gamma, delta))
  # # Transform back to (beta, gamma, delta)
  # eta   <- res$par[1]
  # zeta  <- res$par[2]
  # gamma <- exp(eta)
  # delta <- -beta + exp(zeta)
  # # print (eta)
  # # print (zeta)
  # print (res)
  # if (res$convergence !=0){print ('help')}
  return (c(gamma, delta))
}


seed = 1
i = 1
runs = 1e4
while (length(gamma_mme_vec) < runs){
  set.seed(seed)
  x1 <- rpois(n = n, lambda = alpha)
  x2 <- rpois(n =n, lambda = (delta+ beta*(1 - exp(-gamma*x1))))
  # x2 <- rpois(n =n, lambda = (delta+ beta*(1 - nu^x1)))
  M1 <- sum(x1)/n
  M2<- sum(x2)/n
  S12 <- cov(x1, x2)*(n-1)/n
  S2 <- var(x2) *(n-1)/n
  # if (S2 == M2){
  #   seed<- seed+1
  #   next
  # }
  # if ( M1^2/(exp(M1) - 1) > S12^2/(S2-M2)){
  #   seed<- seed+1
  #   next
  # }
  # if ( M1 < S12^2/(S2-M2)){
  #   seed<- seed+1
  #   next
  # }
  
  alpha_mme <- M1
  mu_mme = exp(alpha_mme)
  # nu_mme = solve_nu_safe(alpha_mme, mu_mme, S12, S2, M2)
  # gamma_mme = -log(nu_mme)
  

  if (M1/S12*beta < exp(1)){
    seed = seed+1
    next
  }
  nu_mme = lambertW0(-S12/M1*beta)+1
  gamma_mme = -log(nu_mme)
  

  delta_mme = M2- beta*(1 - mu_mme^(nu_mme - 1))
  
  
  alpha_mme_vec[i] <- alpha_mme
  gamma_mme_vec[i] <- gamma_mme
  nu_mme_vec[i] = nu_mme
  delta_mme_vec[i] <- delta_mme
  
  # print (c(S12, beta*alpha_mme*(1 - nu_mme)*mu_mme^(nu_mme - 1), 1/S12/M1))

  num =  sqrt(alpha_mme)*beta*(1- nu_mme)*mu_mme^(nu_mme -1)
  den = sqrt( beta^2*(mu_mme^(nu_mme^2 - 1) - mu_mme^(2*nu_mme -2)) + beta*(1 - mu_mme^(nu_mme - 1)) +delta_mme)
  rho_mme[i] = num/den
  
  #MLE 
  mles = mle_solve(x1, x2, n = n, alpha_mme)
  gamma_mle = mles[1]
  delta_mle = mles[2]
  gamma_mle_vec[i] <- gamma_mle
  nu_mle = exp(-gamma_mle)
  delta_mle_vec[i] <- delta_mle
  
  num =  sqrt(alpha_mme)*beta*(1- nu_mle)*mu_mme^(nu_mle -1)
  den = sqrt( beta^2*(mu_mme^(nu_mle^2 - 1) - mu_mme^(2*nu_mle -2)) + beta*(1 - mu_mme^(nu_mle - 1)) +delta_mle)
  rho_mle[i] = num/den
  
  # print (c(S12, beta*alpha_mme*(1 - nu_mle)*mu_mme^(nu_mme - 1), 1/S12/M1))
  
  pearson_corr_vec[i] <- cor(x= x1, y= x2, method = "pearson")
  i<- i+1
  seed = seed+1
}

sim_func(alpha_mme_vec)
sim_func(gamma_mme_vec, gamma_mle_vec)
sim_func(delta_mme_vec, delta_mle_vec)
sim_func(rho_mme, rho_mle)



sqrt(var(alpha_mme_vec))
hist(gamma_mme_vec, breaks =100)
summary(gamma_mme_vec)
hist(delta_mme_vec, breaks = 100)
summary(delta_mme_vec)


hist(gamma_mle_vec, breaks =10000000, xlim = c(0, 1))
summary(gamma_mle_vec)
hist(delta_mle_vec, breaks = 100)
summary(delta_mle_vec)


summary(rho_mme)
summary(rho_mle)

summary(pearson_corr_vec)
























# Case 2: gamma  = 1

rm(list = ls())

alpha = 1
beta = -2
gamma = 1
nu = exp(-gamma)
delta= 25
n = 10000
alpha_mme_vec <- c()
beta_mme_vec <- c()
delta_mme_vec <- c()
rho_mme <- c()
alpha_mle_vec <- c()
beta_mle_vec <- c()
delta_mle_vec <- c()
rho_mle <- c()
pearson_corr_vec <- c()

mle_solve = function(x1, x2, n, alpha_mme){
  obj = function(pars){ 
    beta  = pars[1]
    zeta  = pars[2]
    
    delta = -beta + exp(zeta)    # ensures delta >= -beta
    nu = exp(-gamma)
    
    loglik = -n*(alpha_mme+beta+delta) + log(alpha_mme)*sum(x1) + beta*sum(nu^x1) + 
      sum(x2*log(delta + beta*(1 - nu^x1)))
    return(-loglik)  # optim minimizes
  }
  
  res <- optim(
    # par = c(beta, log(gamma), log(delta+beta)),       # starting values (beta, eta, zeta)
    par = c(-1,0),       # starting values (beta, eta, zeta)
    fn  = obj,
    method = "L-BFGS-B",
    lower = c(-50, -10),  # bounds on (beta, eta, zeta)
    upper = c(0,  10)
  )
  
  # Transform back to (beta, gamma, delta)
  beta  <- res$par[1]
  zeta  <- res$par[2]
  delta <- -beta + exp(zeta)
  return (c(beta, delta))
}


seed = 1
i = 1
runs = 1e4
while (length(beta_mme_vec) < runs){
  set.seed(seed)
  x1 <- rpois(n = n, lambda = alpha)
  x2 <- rpois(n =n, lambda = (delta+ beta*(1 - exp(-gamma*x1))))
  # x2 <- rpois(n =n, lambda = (delta+ beta*(1 - nu^x1)))
  M1 <- sum(x1)/n
  M2<- sum(x2)/n
  S12 <- cov(x1, x2)*(n-1)/n
  S2 <- var(x2) *(n-1)/n

  alpha_mme <- M1
  mu_mme = exp(alpha_mme)

  beta_mme = S12 / (alpha_mme*(1 - exp(-1))*mu_mme^(exp(-1)-1))
  
  delta_mme = M2- beta_mme*(1 - mu_mme^(exp(-1) - 1))
  
  
  alpha_mme_vec[i] <- alpha_mme
  beta_mme_vec[i] = beta_mme
  delta_mme_vec[i] <- delta_mme
  
  
  num =  sqrt(alpha_mme)*beta_mme*(1- exp(-1))*mu_mme^(exp(-1) -1)
  den = sqrt( beta_mme^2*(mu_mme^(exp(-1)^2 - 1) - mu_mme^(2*exp(-1) -2)) + beta_mme*(1 - mu_mme^(exp(-1) - 1)) +delta_mme)
  rho_mme[i] = num/den
  
  #MLE 
  mles = mle_solve(x1, x2, n = n, alpha_mme)
  beta_mle = mles[1]
  delta_mle = mles[2]
  beta_mle_vec[i] <- beta_mle
  delta_mle_vec[i] <- delta_mle
  
  num =  sqrt(alpha_mme)*beta_mle*(1- nu)*mu_mme^(nu -1)
  den = sqrt( beta_mle^2*(mu_mme^(nu^2 - 1) - mu_mme^(2*nu -2)) + beta_mle*(1 - mu_mme^(nu- 1)) +delta_mle)
  rho_mle[i] = num/den
  
  
  pearson_corr_vec[i] <- cor(x= x1, y= x2, method = "pearson")
  i<- i+1
  seed = seed+1
}

hist(alpha_mme_vec, breaks = 1000)
summary(alpha_mme_vec)
hist(beta_mme_vec, breaks =1000)
summary(beta_mme_vec)
hist(beta_mle_vec, breaks = 1000)
summary(beta_mle_vec)

hist(delta_mme_vec, breaks = 100)
summary(delta_mme_vec)



summary(rho_mme)
summary(pearson_corr_vec)





# Case 4: gamma  = 1, beta = 1

rm(list = ls())

alpha = 1
beta = 1
gamma = 1
nu = exp(-gamma)
delta= 4
n = 1000
alpha_mme_vec <- c()
delta_mme_vec <- c()
rho_mme <- c()
delta_mle_vec <- c()
rho_mle <- c()
pearson_corr_vec <- c()

mle_solve = function(x1, x2, n, alpha_mme){
  obj = function(pars){ 
    zeta  = pars[1]
    delta = -beta + exp(zeta)    # ensures delta >= -beta

    loglik = -n*(alpha_mme+beta+delta) + log(alpha_mme)*sum(x1) + beta*sum(nu^x1) + 
      sum(x2*log(delta + beta*(1 - nu^x1)))
    return(-loglik)  # optim minimizes
  }
  
  res <- optim(
    # par = c(beta, log(gamma), log(delta+beta)),       # starting values (beta, eta, zeta)
    par = c(0.1),       # starting values (beta, eta, zeta)
    fn  = obj,
    method = "L-BFGS-B",
    lower = c(-10),  # bounds on (beta, eta, zeta)
    upper = c(10)
  )
  
  # Transform back to (beta, gamma, delta)
  zeta  <- res$par[1]
  delta <- -beta + exp(zeta)
  return (c(delta))
}


seed = 1
i = 1
runs = 1e4
while (length(alpha_mme_vec) < runs){
  set.seed(seed)
  x1 <- rpois(n = n, lambda = alpha)
  x2 <- rpois(n =n, lambda = (delta+ beta*(1 - exp(-gamma*x1))))
  M1 <- sum(x1)/n
  M2<- sum(x2)/n
  S12 <- cov(x1, x2)*(n-1)/n
  S2 <- var(x2) *(n-1)/n
  
  alpha_mme <- M1
  mu_mme = exp(alpha_mme)
  
  delta_mme = M2- beta*(1 - mu_mme^(exp(-1) - 1))
  
  
  alpha_mme_vec[i] <- alpha_mme
  delta_mme_vec[i] = delta_mme

  num =  sqrt(alpha_mme)*beta*(1- exp(-1))*mu_mme^(exp(-1) -1)
  den = sqrt( beta^2*(mu_mme^(exp(-1)^2 - 1) - mu_mme^(2*exp(-1) -2)) + beta*(1 - mu_mme^(exp(-1) - 1)) +delta_mme)
  rho_mme[i] = num/den
  
  #MLE 
  mles = mle_solve(x1, x2, n = n, alpha_mme)
  delta_mle = mles[1]
  delta_mle_vec[i] <- delta_mle
  
  num =  sqrt(alpha_mme)*beta*(1- nu)*mu_mme^(nu -1)
  den = sqrt( beta^2*(mu_mme^(nu^2 - 1) - mu_mme^(2*nu -2)) + beta*(1 - mu_mme^(nu- 1)) +delta_mle)
  rho_mle[i] = num/den
  
  
  pearson_corr_vec[i] <- cor(x= x1, y= x2, method = "pearson")
  i<- i+1
  seed = seed+1
}

hist(alpha_mme_vec, breaks = 1000)
summary(alpha_mme_vec)

hist(delta_mme_vec, breaks = 100)
summary(delta_mme_vec)



summary(rho_mme)
summary(pearson_corr_vec)


#========================================================================
#                 Application
#========================================================================
rm(list = ls())
library(Brobdingnag)
library(readxl)

data<- read.csv("C:/Users/Jared/OneDrive - University of Cape Town/Documents/WST 795-20230920T125726Z-001/WST 795/Arr/Health_Retirement_Study.csv")
xvec <- rep(data$x, data$F)
yvec <- rep(data$y, data$F)
data <- data.frame(x = xvec, y  = yvec)



# data[is.na(data)] <- 0
# data <- data[, 2:9]
# rownames(data )<- c(0:8)
# data<- data.frame(data)
# colnames(data) <- c(0:7)
# xvec <-  c()
# yvec <- c()
# for (i in rownames(data)){
#   for (j in colnames(data)){
#     i <- as.numeric(i)
#     j<- as.numeric(j)
#     print (data[i+1, j+1])
#     xvec <- c(xvec,  rep ( i, data[i+1, j+1]))
#     yvec <- c(yvec, rep ( j, data[i+1, j+1] ))
#   }
# }
# data <- data.frame(x = xvec, y  = yvec)


table(data)
n<- nrow(data)
# x1<- data$x
# x2<- data$y
# Mirrored
x1<- data$y
x2<- data$x



M1 <- mean(x1)
M2<- mean(x2)
S12 <- cov(x1, x2)* (n-1)/n
M12 <- sum(x1*x2 )/ n
S2<-var(x2)*(n-1)/n 
cor(x= x1, y= x2, method = "pearson")


var(x1)/mean(x1)

M1^2/(exp(M1) - 1)
S12^2/(S2 - M2)
M1

M1/S12

term <- 1
for (i in 1:n){
  term <- as.brob(term* factorial(x1[i]) *factorial(x2[i]))
}
h <- log(1/term)





































rm(list = ls())

Accident_data <- read.csv('C:/Users/Jared/OneDrive - University of Cape Town/Documents/WST 795-20230920T125726Z-001/WST 795/Arr/Accident_data.csv')
n<- nrow(Accident_data)
x2<- Accident_data$x
x1<- Accident_data$y
# Mirrored
# x2<- Accident_data$y
# x1<- Accident_data$x
M1 <- mean(x1)
M2<- mean(x2)
S12 <- cov(x1, x2)* (n-1)/n
M12 <- sum(x1*x2 )/ n
S2<-var(x2)*(n-1)/n 
cor(x= x1, y= x2, method = "pearson")


var(x2)/mean(x2)
var(x1)/mean(x1)
cov(x1, x2)

M1^2/(exp(M1) - 1)
S12^2/(S2 - M2)
M1

M1/S12
M1
S12

term <- 1
for (i in 1:n){
  term <- term* factorial(x1[i]) *factorial(x2[i])
}
h = -log(term)





# Full
solve_nu_trans <- function(alpha_mme, mu_mme, S12, S2, M2, trans_start = 0.5) {
  rhs <- S12^2 / (S2 - M2)
  
  obj <- function(trans) {
    # if (nu <= 0 || nu >= 1) return(1e10)  # enforce bounds
    num <- alpha_mme^2* trans
    denom <- mu_mme^trans - 1
    # if (!is.finite(denom) || denom == 0) return(1e10)
    (num / denom - rhs)^2
  }
  
  res <- optim(par = trans_start, fn = obj, method = "L-BFGS-B", lower = 1e-6, upper = 0.9999)
  
  if (res$convergence != 0) warning("Solution may not have converged")
  print (res)
  return(res$par)
}




mle_solve = function(x1, x2, n, alpha_mme){
  obj = function(pars){ 
    dum1  = pars[1]
    eta   = pars[2]
    zeta  = pars[3]
    
    beta = exp(dum1)
    gamma = exp(eta)             # ensures gamma >= 0
    delta = exp(zeta)    # ensures delta >= -beta
    nu = exp(-gamma)
    
    loglik = -n*(alpha_mme+beta+delta) + log(alpha_mme)*sum(x1) + beta*sum(nu^x1) + 
      sum(x2*log(delta + beta*(1 - nu^x1)))
    return(-loglik)  # optim minimizes
  }
  # obj <- function(pars){
  #   dum1 = pars[1]
  #   eta   = pars[2]
  #   zeta  = pars[3]
  #   
  #   beta = exp(dum1)
  #   gamma = exp(eta)             # ensures gamma >= 0
  #   delta = exp(zeta)    # ensures delta >= -beta
  #   nu = exp(-gamma)                       # gamma >= 0
  #   
  #   n  <- length(x1)
  #   lam1 <- alpha_mme                           # X1 mean
  #   lam2 <- delta + beta * (1 - exp(-gamma * x1))       # delta = 0
  #   
  #   ll1 <- sum(dpois(x1, lam1, log = TRUE))
  #   ll2 <- sum(dpois(x2, lam2, log = TRUE))     # dpois(0,0) = 0; handles x1=0,x2=0
  #   ll  <- ll1 + ll2
  #   if (!is.finite(ll)) return(1e12)
  #   -ll
  # }
  
  res <- optim(
    # par = c(beta, log(gamma), log(delta+beta)),       # starting values (beta, eta, zeta)
    par = c(.1, .1, .1),       # starting values (beta, eta, zeta)
    fn  = obj,
    method = "L-BFGS-B",
    # lower = c(0, -10, -10),  # bounds on (beta, eta, zeta)
    # upper = c(10,   10,  10)
    lower = rep(-Inf, 3),
    upper = rep(Inf, 3)
  )
  
  print (res)
  # Transform back to (beta, gamma, delta)
  dum1  <- res$par[1]
  eta   <- res$par[2]
  zeta  <- res$par[3]
  beta =exp(dum1)
  gamma <- exp(eta)
  delta <- exp(zeta)
  return (c(beta, gamma, delta))
}


alpha_mme <- M1
mu_mme = exp(alpha_mme)
# nu_mme = -sqrt(solve_nu_trans(alpha_mme, mu_mme, S12, S2, M2))+1
nu_mme = -sqrt(solve_nu_trans(alpha_mme, mu_mme, S12, S2, M2))+1
gamma_mme = -log(nu_mme)
beta_mme = S12/(alpha_mme*(1 - nu_mme)*mu_mme^(nu_mme-1))
delta_mme = M2- beta_mme*(1 - mu_mme^(nu_mme - 1))
num =  sqrt(alpha_mme)*beta_mme*(1- nu_mme)*mu_mme^(nu_mme -1)
den = sqrt( beta_mme^2*(mu_mme^(nu_mme^2 - 1) - mu_mme^(2*nu_mme -2)) + beta_mme*(1 - mu_mme^(nu_mme - 1)) +delta_mme)
rho_mme = num/den
  
#MLE 
mles = mle_solve(x1, x2, n = n, alpha_mme)
mles
beta_mle = mles[1]
gamma_mle = mles[2]
delta_mle = mles[3]
nu_mle = exp(-gamma_mle)
num =  sqrt(alpha_mme)*beta_mle*(1- nu_mle)*mu_mme^(nu_mle -1)
den = sqrt( beta_mle^2*(mu_mme^(nu_mle^2 - 1) - mu_mme^(2*nu_mle -2)) + beta_mle*(1 - mu_mme^(nu_mle - 1)) +delta_mle)
rho_mle= num/den
  
#AIC
k<- 4
l_full = -n*(alpha_mme+beta_mle+delta_mle) + log(alpha_mme)*sum(x1) + beta_mle*sum(nu_mle^x1) + sum(x2*log(delta_mle + beta_mle*(1 - nu_mle^x1))) +h
aic <- 2*k - 2*l_full
 bic <- k*log(n) - 2*l_full




















# Case I
library("lamW")

solve_nu_optim <- function(alpha_mme, mu_mme, S12, S2, M2, nu_start = nu) {
  rhs <- S12^2 / (S2 - M2)
  
  obj <- function(nu) {
    if (nu <= 0 || nu >= 1) return(1e10)  # enforce bounds
    num <- alpha_mme^2* (nu - 1)^2
    denom <- mu_mme^((nu - 1)^2) - 1
    if (!is.finite(denom) || denom == 0) return(1e10)
    (num / denom - rhs)^2
  }
  
  res <- optim(par = nu_start, fn = obj, method = "L-BFGS-B", lower = 1e-6, upper = 0.9999)
  
  if (res$convergence != 0) warning("Solution may not have converged")
  
  return(res$par)
}

mle_solve = function(x1, x2, n, alpha_mme){
  beta = 1
  # obj = function(pars){
  #   eta   = pars[1]
  #   zeta  = pars[2]
  # 
  #   gamma = exp(eta)             # ensures gamma >= 0
  #   delta = exp(zeta)    # ensures delta >= -beta
  #   nu = exp(-gamma)
  # 
  #   loglik = -n*(alpha_mme+beta+delta) + log(alpha_mme)*sum(x1) + beta*sum(nu^x1) +
  #     sum(x2*log(delta + beta*(1 - nu^x1)))
  #   return(-loglik)  # optim minimizes
  # }
  
  obj <- function(pars){
    eta   = pars[1]
    zeta  = pars[2]

    gamma = exp(eta)             # ensures gamma >= 0
    delta = exp(zeta)    # ensures delta >= -beta
    nu = exp(-gamma)                       # gamma >= 0

    n  <- length(x1)
    lam1 <- alpha_mme                           # X1 mean
    lam2 <- delta + beta * (1 - exp(-gamma * x1))       # delta = 0

    ll1 <- sum(dpois(x1, lam1, log = TRUE))
    ll2 <- sum(dpois(x2, lam2, log = TRUE))     # dpois(0,0) = 0; handles x1=0,x2=0
    ll  <- ll1 + ll2
    if (!is.finite(ll)) return(1e12)
    -ll
  }
  
  res <- optim(
    # par = c(log(gamma_mme), log(delta_mme+beta)),       # starting values (beta, eta, zeta)
    par = c(.1, .1),       # starting values (beta, eta, zeta)
    fn  = obj,
    method = "L-BFGS-B",
    # lower = c(-10, -10),  # bounds on (beta, eta, zeta)
    # upper = c( 10,  10)
    lower = rep(-Inf, 2),
    upper = rep(Inf, 2)
  )
  # res = ga(type = "real-valued",  obj, lower = c(-10, -10), upper = c(10, 10),popSize = 100,maxiter = 25,keepBest = TRUE, monitor = FALSE)
  
  # return( res@solution[1,])
  # eta = res@solution[1,1]
  # zeta = res@solution[1,2]
  # print (eta)
  # gamma <- exp(eta)
  # delta <- -beta + exp(zeta)
  # return  (c(gamma, delta))
  # Transform back to (beta, gamma, delta)
  eta   <- res$par[1]
  zeta  <- res$par[2]
  gamma <- exp(eta)
  delta <- exp(zeta)
  # print (eta)
  # print (zeta)
  print (res)
  if (res$convergence !=0){print ('help')}
  return (c(gamma, delta))
}

alpha_mme <- M1
mu_mme = exp(alpha_mme)
nu_mme = lambertW0(-S12/M1*1)+1
gamma_mme = -log(nu_mme)
delta_mme = M2- 1*(1 - mu_mme^(nu_mme - 1))
num =  sqrt(alpha_mme)*1*(1- nu_mme)*mu_mme^(nu_mme -1)
den = sqrt( 1^2*(mu_mme^(nu_mme^2 - 1) - mu_mme^(2*nu_mme -2)) + 1*(1 - mu_mme^(nu_mme - 1)) +delta_mme)
rho_mme = num/den
  
#MLE 
mles = mle_solve(x1, x2, n = n, alpha_mme)
gamma_mle = mles[1]
delta_mle = mles[2]
nu_mle = exp(-gamma_mle)

num =  sqrt(alpha_mme)*1*(1- nu_mle)*mu_mme^(nu_mle -1)
den = sqrt( 1^2*(mu_mme^(nu_mle^2 - 1) - mu_mme^(2*nu_mle -2)) + 1*(1 - mu_mme^(nu_mle - 1)) +delta_mle)
rho_mle = num/den


#AIC
k<- 3
l = -n*(alpha_mme+1+delta_mle) + log(alpha_mme)*sum(x1) + 1*sum(nu_mle^x1) + sum(x2*log(delta_mle + 1*(1 - nu_mle^x1))) +h
aic <- 2*k - 2*l
bic <- k*log(n) - 2*l

neg2loglam <- -2*(l - l_full)
critical <- qchisq(p =0.95, df = 1)












# Independece Case
mle_solve_ind = function(x1, x2, n, alpha_mme){
  obj = function(pars){
    dum0 = pars[1]
    dum1   = pars[2]
    lam1 = exp(dum0)
    lam2  = exp(dum1)


    # loglik = -n*(alpha_mme+lam2) + log(alpha_mme)*sum(x1) + log(lam2)*sum(x2)
    loglik = -n*(lam1+lam2) + log(lam1)*sum(x1) + log(lam2)*sum(x2)
    
    return(-loglik)  # optim minimizes
  }

  res <- optim(
    par = c(.1, .1),      
    fn  = obj,
    method = "L-BFGS-B",
    lower = rep(-Inf, 2),
    upper = rep(Inf, 2)
  )
  dum0 =  res$par[1]
  dum1   <- res$par[2]
  lam1 = exp(dum0)
  lam2 <- exp(dum1)
  print (res)
  return (c(lam1, lam2))
}

mles = mle_solve_ind(x1, x2, n = n, alpha_mme)
lam1_mle = mles[1]
lam2_mle = mles[2]
#AIC
k<- 2
l_full = -n*(alpha_mme+beta_mle+delta_mle) + log(alpha_mme)*sum(x1) + beta_mle*sum(nu_mle^x1) + sum(x2*log(delta_mle + beta_mle*(1 - nu_mle^x1))) +h
l = -n*(lam1_mle+lam2_mle) + log(lam1_mle)*sum(x1) + log(lam2_mle)*sum(x2)+h
aic <- 2*k - 2*l
neg2loglam <- -2*(l - l_full)
critical <- qchisq(p =0.95, df = 2)






















# Case II


mle_solve = function(x1, x2, n, alpha_mme){
  obj = function(pars){ 
    dum1  = pars[1]
    beta = exp(dum1)
    zeta  = pars[2]
    
    delta = exp(zeta)  # ensures delta >= -beta
    nu = exp(-1)
    
    loglik = -n*(alpha_mme+beta+delta) + log(alpha_mme)*sum(x1) + beta*sum(nu^x1) + 
      sum(x2*log(delta + beta*(1 - nu^x1)))
    return(-loglik)  # optim minimizes
  }
  
  res <- optim(
    # par = c(beta_mme, log(delta_mme+beta_mme)),       # starting values (beta, eta, zeta)
    par = c(.1,.1),       # starting values (beta, eta, zeta)
    fn  = obj,
    method = "L-BFGS-B",
    # lower = c(0, -10),  # bounds on (beta, eta, zeta)
    # upper = c(10,  10)
    lower = rep(-Inf, 2),
    upper = rep(Inf, 2)
  )
  
  # Transform back to (beta, gamma, delta)
  dum1  <- res$par[1]
  beta = exp(dum1)
  zeta  <- res$par[2]
  delta <- exp(zeta)
  return (c(beta, delta))
}

alpha_mme <- M1
mu_mme = exp(alpha_mme)
beta_mme = S12 / (alpha_mme*(1 - exp(-1))*mu_mme^(exp(-1)-1))
delta_mme = M2- beta_mme*(1 - mu_mme^(exp(-1) - 1))


num =  sqrt(alpha_mme)*beta_mme*(1- exp(-1))*mu_mme^(exp(-1) -1)
den = sqrt( beta_mme^2*(mu_mme^(exp(-1)^2 - 1) - mu_mme^(2*exp(-1) -2)) + beta_mme*(1 - mu_mme^(exp(-1) - 1)) +delta_mme)
rho_mme = num/den

#MLE 
mles = mle_solve(x1, x2, n = n, alpha_mme)
mles
beta_mle = mles[1]
delta_mle = mles[2]

num =  sqrt(alpha_mme)*beta_mle*(1- exp(-1))*mu_mme^(exp(-1) -1)
den = sqrt( beta_mle^2*(mu_mme^(exp(-1)^2 - 1) - mu_mme^(2*exp(-1) -2)) + beta_mle*(1 - mu_mme^(exp(-1)- 1)) +delta_mle)
rho_mle= num/den

#AIC
k<- 3
l = -n*(alpha_mme+beta_mle+delta_mle) + log(alpha_mme)*sum(x1) + beta_mle*sum(exp(-1)^x1) + sum(x2*log(delta_mle + beta_mle*(1 - exp(-1)^x1))) +h
aic <- 2*k - 2*l
bic <- k*log(n) - 2*l

neg2loglam <- -2*(l - l_full)
critical <- qchisq(p =0.95, df = 1)





# Case III

solve_nu_trans <- function(alpha_mme, mu_mme, S12, S2, M2, trans_start = 0.5) {
  rhs <- S12^2 / (S2 - M2)
  
  obj <- function(trans) {
    # if (nu <= 0 || nu >= 1) return(1e10)  # enforce bounds
    num <- alpha_mme^2* trans
    denom <- mu_mme^trans - 1
    # if (!is.finite(denom) || denom == 0) return(1e10)
    (num / denom - rhs)^2
  }
  
  res <- optim(par = trans_start, fn = obj, method = "L-BFGS-B", lower = 1e-6, upper = 0.9999)
  
  if (res$convergence != 0) warning("Solution may not have converged")
  
  return(res$par)
}

 
mle_solve = function(x1, x2, n, alpha_mme){
  delta = 0
  # obj = function(pars){
  #   beta  = pars[1]
  #   eta   = pars[2]
  # 
  #   gamma = exp(eta)             # ensures gamma >= 0
  #   nu = exp(-gamma)
  # 
  #   loglik = -n*(alpha_mme+beta+delta) + log(alpha_mme)*sum(x1) + beta*sum(nu^x1)
  #   +sum((x2[(x1>0)])*log(delta + beta*(1 - nu^(x1[(x1>0)]))))
  #   return(-loglik)  # optim minimizes
  # }

  obj <- function(pars){
    dum1 <- pars[1]
    beta = exp(dum1)
    eta  <- pars[2]
    gamma <- exp(eta)                           # gamma >= 0

    n  <- length(x1)
    lam1 <- alpha_mme                           # X1 mean
    lam2 <- beta * (1 - exp(-gamma * x1))       # delta = 0

    ll1 <- sum(dpois(x1, lam1, log = TRUE))
    ll2 <- sum(dpois(x2, lam2, log = TRUE))     # dpois(0,0) = 0; handles x1=0,x2=0
    ll  <- ll1 + ll2
    if (!is.finite(ll)) return(1e12)
    -ll
  }

  
  res <- optim(
    # par = c(beta, log(gamma), log(delta+beta)),       # starting values (beta, eta, zeta)
    par = c(.1, 0.1),       # starting values (beta, eta, zeta)
    fn  = obj,
    method = "L-BFGS-B",
    lower = c(-Inf, Inf),  # bounds on (beta, eta, zeta)
    upper = c(Inf,   Inf)
  )
  
  # Transform back to (beta, gamma, delta)
  dum1  <- res$par[1]
  beta = exp(dum1)
  eta   <- res$par[2]
  gamma <- exp(eta)
  print (res)
  return (c(beta, gamma, res$value))
}


alpha_mme <- M1
mu_mme = exp(alpha_mme)
# nu_mme = -sqrt(solve_nu_trans(alpha_mme, mu_mme, S12, S2, M2))+1
nu_mme = -sqrt(solve_nu_trans(alpha_mme, mu_mme, S12, S2, M2))+1
gamma_mme = -log(nu_mme)
beta_mme = S12/(alpha_mme*(1 - nu_mme)*mu_mme^(nu_mme-1))
delta_mme = M2- beta_mme*(1 - mu_mme^(nu_mme - 1))
num =  sqrt(alpha_mme)*beta_mme*(1- nu_mme)*mu_mme^(nu_mme -1)
den = sqrt( beta_mme^2*(mu_mme^(nu_mme^2 - 1) - mu_mme^(2*nu_mme -2)) + beta_mme*(1 - mu_mme^(nu_mme - 1)) +delta_mme)
rho_mme = num/den

#MLE 
mles = mle_solve(x1, x2, n = n, alpha_mme)
mles
beta_mle = mles[1]
gamma_mle = mles[2]
nu_mle = exp(-gamma_mle)
num =  sqrt(alpha_mme)*beta_mle*(1- nu_mle)*mu_mme^(nu_mle -1)
den = sqrt( beta_mle^2*(mu_mme^(nu_mle^2 - 1) - mu_mme^(2*nu_mle -2)) + beta_mle*(1 - mu_mme^(nu_mle - 1)) + 0)
rho_mle= num/den

#AIC
k<- 3
l = -n*(alpha_mme+beta_mle+0) + log(alpha_mme)*sum(x1) + beta_mle*sum(nu_mle^x1) + sum(x2*log(0 + beta_mle*(1- nu_mle^x1)), na.rm = TRUE) +h
aic <- 2*k - 2*l
bic <- k*log(n) - 2*l

neg2loglam <- -2*(l - l_full)
critical <- qchisq(p =0.95, df = 1)



















































M2
1- (exp(M1)^(1/exp(1) - 1))

M2
1 - 1/M1*(1 - 1/exp(M1))


# Case IV

mle_solve = function(x1, x2, n, alpha_mme){
  beta = 1
  gamma = 1
  nu = exp(-gamma)
  obj = function(pars){ 
    zeta  = pars[1]
    delta = exp(zeta)    # ensures delta >= -beta
    
    loglik = -n*(alpha_mme+beta+delta) + log(alpha_mme)*sum(x1) + beta*sum(nu^x1) + 
      sum(x2*log(delta + beta*(1 - nu^x1)))
    return(-loglik)  # optim minimizes
  }
  
  res <- optim(
    # par = c(beta, log(gamma), log(delta+beta)),       # starting values (beta, eta, zeta)
    par = c(0.1),       # starting values (beta, eta, zeta)
    fn  = obj,
    method = "L-BFGS-B",
    lower = c(-Inf),  # bounds on (beta, eta, zeta)
    upper = c(Inf)
  )
  
  # Transform back to (beta, gamma, delta)
  zeta  <- res$par[1]
  delta <- exp(zeta)
  return (c(delta))
}

alpha_mme <- M1
mu_mme = exp(alpha_mme)
delta_mme = M2- 1*(1 - mu_mme^(exp(-1) - 1))
num =  sqrt(alpha_mme)*1*(1- exp(-1))*mu_mme^(exp(-1) -1)
den = sqrt( 1^2*(mu_mme^(exp(-1)^2 - 1) - mu_mme^(2*exp(-1) -2)) + 1*(1 - mu_mme^(exp(-1) - 1)) +delta_mme)
rho_mme = num/den
  
#MLE 
mles = mle_solve(x1, x2, n = n, alpha_mme)
delta_mle = mles[1]
num =  sqrt(alpha_mme)*1*(1- exp(-1))*mu_mme^(exp(-1) -1)
den = sqrt( 1^2*(mu_mme^(exp(-1)^2 - 1) - mu_mme^(2*exp(-1) -2)) + 1*(1 - mu_mme^(exp(-1)- 1)) +delta_mle)
rho_mle = num/den

#AIC
k<- 2
l = -n*(alpha_mme+1+delta_mle) + log(alpha_mme)*sum(x1) + 1*sum(exp(-1)^x1) + sum(x2*log(delta_mle + 1*(1 - exp(-1)^x1))) +h
aic <- 2*k - 2*l
bic <- k*log(n) - 2*l  


neg2loglam <- -2*(l - l_full)
critical <- qchisq(p =0.95, df = 2)










# Case V
library(GA)

mle_solve = function(x1, x2, n, alpha_mme){
  delta = 0
  gamma = 1
  # obj = function(pars){ 
  #   beta  = pars[1]
  #   nu = exp(-1)
  #   
  #   loglik = -n*(alpha_mme+beta+delta) + log(alpha_mme)*sum(x1) + beta*sum(nu^x1) + 
  #     sum(x2*log(delta + beta*(1 - nu^x1) + 1e-6))
  #   return(loglik)  # optim minimizes
  # }
  obj <- function(pars){
    beta <- pars[1]
                       # gamma >= 0
    
    n  <- length(x1)
    lam1 <- alpha_mme                           # X1 mean
    lam2 <- beta * (1 - exp(-gamma * x1))       # delta = 0
    
    ll1 <- sum(dpois(x1, lam1, log = TRUE))
    ll2 <- sum(dpois(x2, lam2, log = TRUE))     # dpois(0,0) = 0; handles x1=0,x2=0
    ll  <- ll1 + ll2
    if (!is.finite(ll)) return(1e12)
    -ll
  }
  
  res <- optim(
    # par = c(beta_mme),       # starting values (beta, eta, zeta)
    par = c(1),       # starting values (beta, eta, zeta)
    fn  = obj,
    method = "L-BFGS-B",
    lower = c(0),  # bounds on (beta, eta, zeta)
    upper = c(10)
  )
  
  # Transform back to (beta, gamma, delta)
  beta  <- res$par[1]
  return (c(beta))
}

alpha_mme <- M1
mu_mme = exp(alpha_mme)
# beta_mme = S12 / (alpha_mme*(1 - exp(-1))*mu_mme^(exp(-1)-1))
beta_mme = M2/(1 - mu_mme^(exp(-1) - 1))


num =  sqrt(alpha_mme)*beta_mme*(1- exp(-1))*mu_mme^(exp(-1) -1)
den = sqrt( beta_mme^2*(mu_mme^(exp(-1)^2 - 1) - mu_mme^(2*exp(-1) -2)) + beta_mme*(1 - mu_mme^(exp(-1) - 1)) +0)
rho_mme = num/den

#MLE 
mles = mle_solve(x1, x2, n = n, alpha_mme)
mles
beta_mle = mles[1]

num =  sqrt(alpha_mme)*beta_mle*(1- exp(-1))*mu_mme^(exp(-1) -1)
den = sqrt( beta_mle^2*(mu_mme^(exp(-1)^2 - 1) - mu_mme^(2*exp(-1) -2)) + beta_mle*(1 - mu_mme^(exp(-1)- 1)) +0)
rho_mle= num/den

#AIC
k<- 2
l = -n*(alpha_mme+beta_mle+0) + log(alpha_mme)*sum(x1) + beta_mle*sum(exp(-1)^x1) + sum(x2*log(0 + beta_mle*(1- exp(-1)^x1)), na.rm = TRUE) +h
# logterm = 0 + beta_mle*(1 - exp(-1)^x1)
# logterm[logterm>0]
# sum(x2*logterm)
# l = -n*(alpha_mme+beta_mle+0) + log(alpha_mme)*sum(x1) + beta_mle*sum(exp(-1)^x1) + sum(x2*logterm) +h
aic <- 2*k - 2*l
bic <- k*log(n) - 2*l
neg2loglam <- -2*(l - l_full)
critical <- qchisq(p =0.95, df = 2)




# Plots
## Plots Traffic


### Non-Mirrored

library(ggplot2)
library(dplyr)

freq_data <- Accident_data %>%
  transmute(x1 = y, x2 = x) %>%
  filter(!is.na(x1), !is.na(x2)) %>%
  count(x1, x2, name = "n") %>%
  arrange(x1, x2)

# Grid for smooth curves
x1_grid <- seq(min(x1, na.rm = TRUE), max(x1, na.rm = TRUE), length.out = 400)

# Three models
exp_case4   <- 0.821 + 1.0   * (1 - exp(-1 * x1_grid))
lomax_case2 <- 0.813 + 1.769 * (1 - (1 / (1 + x1_grid)))
bpp_sm1     <- 0.815 + 0.815 * x1_grid

# Full curves
curve_df <- tibble(
  x1 = rep(x1_grid, 3),
  x2_hat = c(exp_case4, lomax_case2, bpp_sm1),
  model = factor(rep(c("Exponential Case IV",
                       "Lomax Case II",
                       "A&M BPP SM-I"), each = length(x1_grid)))
)

# Model points at x1 = 0, 1, 2
x_points <- c(0, 1, 2)
points_df <- tibble(
  x1 = rep(x_points, 3),
  x2_hat = c(
    0.821 + 1.0   * (1 - exp(-1 * x_points)),            # Exponential
    0.813 + 1.769 * (1 - (1 / (1 + x_points))),          # Lomax
    0.815 + 0.815 * x_points                             # BPP SM1
  ),
  model = factor(rep(c("Exponential Case IV",
                       "Lomax Case II",
                       "A&M BPP SM-I"), each = length(x_points)))
)

# Colors + shapes
color_map <- c("Exponential Case IV" = "#191970",  # Midnight Blue
               "Lomax Case II"       = "#1C1F2A",  # Hermes Havane
               "A&M BPP SM-I"         = "dodgerblue")
shape_map <- c("Exponential Case IV" = 16,
               "Lomax Case II"       = 17,
               "A&M BPP SM-I"         = 15)

ggplot() +
  # Scatter points (black, no legend)
  geom_point(data = freq_data,
             aes(x = x1, y = x2, size = n),
             color = "black", alpha = 1, show.legend = FALSE) +
  
  # Frequency labels (black)
  geom_text(data = freq_data,
            aes(x = x1, y = x2, label = n),
            color = "black", size = 4,
            hjust = -0.5, vjust = -0.2) +
  
  # Model curves
  geom_line(data = curve_df,
            aes(x = x1, y = x2_hat, color = model),
            linewidth = 0.8) +
  
  # Model points (at x1 = 0,1,2)
  geom_point(data = points_df,
             aes(x = x1, y = x2_hat, color = model, shape = model),
             size = 3) +
  
  scale_color_manual(values = color_map) +
  scale_shape_manual(values = shape_map) +
  
  scale_x_continuous(breaks = seq(floor(min(x1, na.rm = TRUE)),
                                  ceiling(max(x1, na.rm = TRUE)), by = 1)) +
  
  labs(x = expression(x[1]),
       y = expression(E(X[2] ~ "|" ~ X[1] == x[1]))) +
  
  theme_minimal(base_size = 14) +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.95, 0.95),         # inside, top right
    legend.justification = c("right", "top"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text  = element_text(size = 14), 
    legend.background    = element_rect(fill = "white", color = "white"),
    legend.box.background = element_rect(fill = "white", color = "grey70", linewidth = 0.5),
    legend.box.margin     = margin(3, 3, 3, 3)
    
  )




### Mirrored

freq_data <- Accident_data %>%
  transmute(x1 = x, x2 = y) %>%
  filter(!is.na(x1), !is.na(x2)) %>%
  count(x1, x2, name = "n") %>%
  arrange(x1, x2)


library(ggplot2)
library(dplyr)

# Grid for smooth curves
x1_grid <- seq(min(x1, na.rm = TRUE), max(x1, na.rm = TRUE), length.out = 400)

# Three mirrored models
exp_case5   <- 0 + 0.143 * (1 - exp(-1 * x1_grid))
lomax_case5 <- 0 + 0.181 * (1 - (1 / (1 + x1_grid)))
bpp_msm2     <- 0 + 0.067 * x1_grid

# Full curves (use the mirrored series here)
curve_df <- tibble(
  x1 = rep(x1_grid, 3),
  x2_hat = c(exp_case5, lomax_case5, bpp_msm2),
  model = factor(rep(c("Mirrored Exponential Case V",
                       "Mirrored Lomax Case V",
                       "A&M BPP MSM-II"), each = length(x1_grid)))
)

# Model points at x1 = 0:5
x_points <- 0:5
points_df <- tibble(
  x1 = rep(x_points, 3),
  x2_hat = c(
    0 + 0.143 * (1 - exp(-1 * x_points)),           # Mirrored Exponential
    0 + 0.181 * (1 - (1 / (1 + x_points))),         # Mirrored Lomax
    0 + 0.067 * x_points                            # BPP SM-II
  ),
  model = factor(rep(c("Mirrored Exponential Case V",
                       "Mirrored Lomax Case V",
                       "A&M BPP MSM-II"), each = length(x_points)))
)

# Colors + shapes (same palette)
color_map <- c("Mirrored Exponential Case V" = "#191970",  # Midnight Blue
               "Mirrored Lomax Case V"       = "#1C1F2A",  # Hermes Havane
               "A&M BPP MSM-II"               = "dodgerblue")
shape_map <- c("Mirrored Exponential Case V" = 16,
               "Mirrored Lomax Case V"       = 17,
               "A&M BPP MSM-II"               = 15)

ggplot() +
  # Scatter points (black, no legend)
  geom_point(data = freq_data,
             aes(x = x1, y = x2, size = n),
             color = "black", alpha = 1, show.legend = FALSE) +
  # Frequency labels (black, top-right)
  geom_text(data = freq_data,
            aes(x = x1, y = x2, label = n),
            color = "black", size = 4,
            hjust = -0.5, vjust = 1.1) +
  # Model curves
  geom_line(data = curve_df,
            aes(x = x1, y = x2_hat, color = model),
            linewidth = 0.8) +
  # Model points only at x = 0:5
  geom_point(data = points_df,
             aes(x = x1, y = x2_hat, color = model, shape = model),
             size = 3) +
  scale_color_manual(values = color_map) +
  scale_shape_manual(values = shape_map) +
  # Integer ticks 0..5 (adjust if needed)
  scale_x_continuous(breaks = 0:5, limits = c(0, 5)) +
  labs(x = expression(x[1]),
       y = expression(E(X[2] ~ "|" ~ X[1] == x[1]))) +
  theme_minimal(base_size = 14) +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.95, 0.95),         # inside, top-right
    legend.justification = c("right", "top"),
    legend.background     = element_rect(fill = "white", color = "white"),
    legend.box.background = element_rect(fill = "white", color = "grey70", linewidth = 0.5),
    legend.box.margin     = margin(3, 3, 3, 3),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text  = element_text(size = 14)
  )



## Plots Health


### Non-Mirrored
library(ggplot2)
library(dplyr)

# Frequencies (x1 from x, x2 from y)
freq_data <- data %>%
  transmute(x1 = x, x2 = y) %>%
  filter(!is.na(x1), !is.na(x2)) %>%
  count(x1, x2, name = "n") %>%
  arrange(x1, x2)

# Grid and points over requested ranges
x1_grid   <- seq(0, 8, length.out = 9)
x_points  <- 0:8

# Three models
full_exp    <- 3.06066e-07 + 0.813 * (1 - exp(-2202170 * x1_grid))
lomax_case3 <- 0 + 0.813 * (1 - ( 7.312743e-06 / ( 7.312743e-06+ x1_grid)))
bpp_fm      <- 0.64 + 0.049 * x1_grid

# Curves
curve_df <- tibble(
  x1 = rep(x1_grid, 3),
  x2_hat = c(full_exp, lomax_case3, bpp_fm),
  model = factor(rep(c("Full Exponential","Lomax Case III","A&M BPP FM"), each = length(x1_grid)))
)

# Model points at x = 0:8
points_df <- tibble(
  x1 = rep(x_points, 3),
  x2_hat = c(
     3.06066e-07 + 0.813 * (1 - exp(-2202170 * x_points)),
     0 + 0.813 * (1 - ( 7.312743e-06 / ( 7.312743e-06+ x_points))),
     0.64 + 0.049 * x_points                   # A&M BPP FM
  ),
  model = factor(rep(c("Full Exponential","Lomax Case III","A&M BPP FM"), each = length(x_points)))
)

# Colors + shapes
color_map <- c("Full Exponential" = "#191970", "Lomax Case III" = "#1C1F2A", "A&M BPP FM" = "dodgerblue")
shape_map <- c("Full Exponential" = 16, "Lomax Case III" = 17, "A&M BPP FM" = 15)

ggplot() +
  geom_point(data = freq_data, aes(x = x1, y = x2, size = n), color = "black", alpha = 1, show.legend = FALSE) +
  geom_text(data = freq_data, aes(x = x1, y = x2, label = n), color = "black", size = 4, hjust = -0.3, vjust = 1.1) +
  geom_line(data = curve_df, aes(x = x1, y = x2_hat, color = model), linewidth = 0.8) +
  geom_point(data = points_df, aes(x = x1, y = x2_hat, color = model, shape = model), size = 3) +
  scale_color_manual(values = color_map) +
  scale_shape_manual(values = shape_map) +
  scale_x_continuous(breaks = 0:8, limits = c(0, 8)) +
  scale_y_continuous(breaks = 0:4, limits = c(0, 4)) +
  labs(x = expression(x[1]), y = expression(E(X[2] ~ "|" ~ X[1] == x[1]))) +
  theme_minimal(base_size = 14) +
  theme(
    legend.title         = element_blank(),
    legend.position      = c(0.85, 0.91),
    legend.justification = c("left", "top"),
    # legend.background     = element_rect(fill = "white", color = "white"),
    # legend.box.background = element_rect(fill = "white", color = "white", linewidth = 0.5),
    legend.box.margin     = margin(3, 3, 3, 3),
    panel.grid.major     = element_blank(),
    panel.grid.minor     = element_blank(),
    panel.border         = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.title           = element_text(size = 16, face = "bold"),
    axis.text            = element_text(size = 14)
  )



### Mirrored

library(ggplot2)
library(dplyr)

# Frequencies (x1 from y, x2 from x)
freq_data <- data %>%
  transmute(x1 = y, x2 = x) %>%
  filter(!is.na(x1), !is.na(x2)) %>%
  count(x1, x2, name = "n") %>%
  arrange(x1, x2)

# Grid and points over requested ranges
x1_grid <- seq(0, 4, length.out = 9)
x_points <- 0:4

# Three models
exp_case1    <- 2.559 + 1 * (1 - exp(-0.121 * x1_grid))
lomax_case1  <- 2.556 + 1 * (1 - (7.181 / (7.181 + x1_grid)))
bpp_mfm      <- 2.563 + 0.103 * x1_grid

# Curves
curve_df <- tibble(
  x1 = rep(x1_grid, 3),
  x2_hat = c(exp_case1, lomax_case1, bpp_mfm),
  model = factor(rep(c("Exponential Case I",
                       "Lomax Case I",
                       "A&M BPP MFM"), each = length(x1_grid)))
)

# Model points at x = 0:4
points_df <- tibble(
  x1 = rep(x_points, 3),
  x2_hat = c(
    2.559 + 1 * (1 - exp(-0.121 * x_points)),            # Exponential
    2.556 + 1 * (1 - (7.181 / (7.181 + x_points))),      # Lomax
    2.563 + 0.103 * x_points                              # BPP MFM
  ),
  model = factor(rep(c("Exponential Case I",
                       "Lomax Case I",
                       "A&M BPP MFM"), each = length(x_points)))
)

# Colors + shapes
color_map <- c("Exponential Case I" = "#191970",  # Midnight Blue
               "Lomax Case I"       = "#1C1F2A",  # Hermes Havane
               "A&M BPP MFM"        = "dodgerblue")
shape_map <- c("Exponential Case I" = 16,
               "Lomax Case I"       = 17,
               "A&M BPP MFM"        = 15)

ggplot() +
  # Scatter points (black, no legend)
  geom_point(data = freq_data,
             aes(x = x1, y = x2, size = n),
             color = "black", alpha = 1, show.legend = FALSE) +
  # Frequency labels (black, bottom-right of dots)
  geom_text(data = freq_data,
            aes(x = x1, y = x2, label = n),
            color = "black", size = 4,
            hjust = -0.3, vjust = -0.2) +
  # Model curves
  geom_line(data = curve_df,
            aes(x = x1, y = x2_hat, color = model),
            linewidth = 1.2) +
  # Model points only at x = 0:4
  geom_point(data = points_df,
             aes(x = x1, y = x2_hat, color = model, shape = model),
             size = 3) +
  scale_color_manual(values = color_map) +
  scale_shape_manual(values = shape_map) +
  # Integer ticks 0..4
  scale_x_continuous(breaks = 0:4, limits = c(0, 4)) +
  labs(x = expression(x[1]),
       y = expression(E(X[2] ~ "|" ~ X[1] == x[1]))) +
  theme_minimal(base_size = 14) +
  theme(
    legend.title         = element_blank(),
    legend.position      = c(0.93, 0.98),         # inside, top-right
    legend.justification = c("right", "top"),
    legend.background     = element_rect(fill = "white", color = "white"),
    legend.box.background = element_rect(fill = "white", color = "gray70", linewidth = 0.5),
    legend.box.margin     = margin(3, 3, 3, 3),
    panel.grid.major     = element_blank(),
    panel.grid.minor     = element_blank(),
    panel.border         = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.title           = element_text(size = 16, face = "bold"),
    axis.text            = element_text(size = 14)
  )


table(x1, x2)
