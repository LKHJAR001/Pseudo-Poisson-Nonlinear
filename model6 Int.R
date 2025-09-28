# Lomax(1, gamma)
rm(list = ls())

library("lamW")
library(nleqslv)
library(GenSA)
library(GA)
library(reticulate)
# py_install("mpmath")
mpmath <- import("mpmath")

alpha =5
beta = -20
gamma = .5
delta= 25
n = 100
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



solve_gamma = function(alpha_mme, mu_mme,  S2, M2, S12){
  func <- function(gamma){
    f1 <- as.numeric( as.character( mpmath$hyp1f1(gamma, gamma+1, alpha_mme) ) )
    f11 <-as.numeric( as.character( mpmath$hyp1f1(gamma+1, gamma+2, alpha_mme) )) 
    f2<- as.numeric( as.character( mpmath$hyp2f2(gamma, gamma, gamma+1, gamma+1, alpha_mme) ) ) 
    return ((alpha_mme^2*(f1- (gamma/(gamma+1)*f11))^2/(mu_mme*f2 - f1^2) - S12^2/(S2-M2))^2)
  }
  
  # options <- list("ftol" = 1e-20, "xtol" = 1e-20, "btol" = 1e-20)
  res <- optim(par = 0.1, fn = func, method = "L-BFGS-B", lower = 0, upper = gamma*5)
  return(res$par)
  
  # res = ga(type = "real-valued",  func, lower = 0, upper = gamma*2,popSize = 100,maxiter = 25,keepBest = TRUE, monitor = FALSE)
  # return( res@solution[1,])
  
}
# solve_gamma(alpha_mme, mu_mme,  S2, M2, S12)


solve_nu_trans <- function(alpha_mme, mu_mme, S12, S2, M2, trans_start) {
  rhs <- S12^2 / (S2 - M2)
  obj <- function(trans) {
    num <- alpha_mme^2 * trans
    denom <- mu_mme^trans - 1
    return(-(num / denom - rhs)^2)  # explicit return stops printing
  }
  
  # res <- optim(par = trans_start, fn = obj, method = "L-BFGS-B", lower = 0, upper = 1)
  # return(res$par)
  # res <- suppressWarnings(GenSA(
  #   par    = trans_start,
  #   fn     = obj,
  #   lower  = 0,
  #   upper  = 1,
  #   control = list(max.call = 5000, verbose = FALSE) # you can increase if needed
  # ))
  res = ga(type = "real-valued",  obj, lower = 0, upper = 1,popSize = 100,maxiter = 25,keepBest = TRUE, monitor= FALSE)
  
  
  return( res@solution[1,])
}




mle_solve = function(x1, x2, n, alpha_mme){
  obj = function(pars){ 
    beta  = pars[1]
    eta   = pars[2]
    zeta  = pars[3]
    
    gamma = exp(eta)             # ensures gamma >= 0
    delta = -beta + exp(zeta)    # ensures delta >= -beta
    
    loglik = -n*(alpha_mme+beta+delta) + 
      log(alpha_mme)*sum(x1) + 
      beta*sum(gamma/(x1+gamma)) + 
      sum(x2*log(delta + beta*(1 - (gamma/(x1+gamma)))))
    
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
  x2 <- rpois(n =n, lambda = (delta+ beta*(1 - gamma/(x1+gamma))))
  M1 <- sum(x1)/n
  M2<- sum(x2)/n
  S12 <- cov(x1, x2)*(n-1)/n
  S2 <- var(x2) *(n-1)/n

  if (S2 !=M2 & (M1^2/(exp(M1) - 1) < S12^2/(S2-M2)) & (M1 > S12^2/(S2-M2)) )
  {
  
    
  alpha_mme <- M1
  mu_mme = exp(alpha_mme)
  gamma_mme = solve_gamma(alpha_mme, mu_mme,  S2, M2, S12)
  f1 <- as.numeric( as.character( mpmath$hyp1f1(gamma_mme, gamma_mme+1, alpha_mme) ) )
  f11 <-as.numeric( as.character( mpmath$hyp1f1(gamma_mme+1, gamma_mme+2, alpha_mme) )) 
  f2<- as.numeric( as.character( mpmath$hyp2f2(gamma_mme, gamma_mme, gamma_mme+1, gamma_mme+1, alpha_mme) ) ) 
  
  
  beta_mme = S12/(alpha_mme/mu_mme*(f1 - (gamma_mme/(gamma_mme+1)*f11)))
  delta_mme = M2- beta_mme*(1 - 1/mu_mme*f1)
  
  

  alpha_mme_vec[i] <- alpha_mme
  beta_mme_vec[i] <- beta_mme
  gamma_mme_vec[i] <- gamma_mme
  delta_mme_vec[i] <- delta_mme
  

  num =  sqrt(alpha_mme)*beta_mme/mu_mme*(f1 - (gamma_mme/(1+gamma_mme)*f11))
  den = sqrt(beta_mme^2*(1/mu_mme*f2 - 1/mu_mme^2*f1^2) + beta_mme*(1 - 1/mu_mme*f1) + delta_mme)
  rho_mme[i] = num/den


  #MLE 
  mles = mle_solve(x1, x2, n = n, alpha_mme)
  beta_mle = mles[1]
  gamma_mle = mles[2]
  delta_mle = mles[3]
  beta_mle_vec[i] <- beta_mle
  gamma_mle_vec[i] <- gamma_mle
  delta_mle_vec[i] <- delta_mle
  f1 <- as.numeric( as.character( mpmath$hyp1f1(gamma_mle, gamma_mle+1, alpha_mme) ) )
  f11 <-as.numeric( as.character( mpmath$hyp1f1(gamma_mle+1, gamma_mle+2, alpha_mme) ))
  f2<- as.numeric( as.character( mpmath$hyp2f2(gamma_mle, gamma_mle, gamma_mle+1, gamma_mle+1, alpha_mme) ) )


  num =  sqrt(alpha_mme)*beta_mle/mu_mme*(f1 - (gamma_mle/(1+gamma_mle)*f11))
  den = sqrt(beta_mle^2*(1/mu_mme*f2 - 1/mu_mme^2*f1^2) + beta_mle*(1 - 1/mu_mme*f1) + delta_mle)
  rho_mle[i] = num/den
  pearson_corr_vec[i] <- cor(x= x1, y= x2, method = "pearson")

  i<- i+1
  seed = seed+1
  }else{seed = seed+1; next}
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
sim_func(beta_mme_vec, beta_mle_vec)
sim_func(gamma_mme_vec, gamma_mle_vec)
sim_func(delta_mme_vec, delta_mle_vec)
sim_func(rho_mme, rho_mle)





hist(beta_mme_vec_alt, breaks = 1000)

hist(beta_mme_vec, breaks = 1000)
summary(beta_mme_vec)
hist(gamma_mme_vec, breaks =1000)
hist(gamma_mme_vec_alt, breaks =1000)

summary(gamma_mme_vec)
summary(gamma_mle_vec)
dens = density(gamma_mme_vec, kernel = c('gaussian'))
dens$x[which.max(dens$y)]
hist(delta_mme_vec, breaks = 100)
summary(delta_mme_vec)
summary(delta_mle_vec)

hist(exprvec, breaks =100)

summary(rho_mme)
summary(pearson_corr_vec)

hist(rho_mme_alt, breaks=100)
hist(rho_mme, breaks = 100)
hist(pearson_corr_vec, breaks = 100)


# Case I
rm(list = ls())
library(reticulate)
# py_install("mpmath")
mpmath <- import("mpmath")

alpha = 0.7
beta = -1
gamma = 0.05
delta= 25
n = 10000
alpha_mme_vec <- c()
gamma_mme_vec <- c()
gamma_mle_vec <- c()
delta_mme_vec <- c()
delta_mle_vec <- c()
rho_mme <- c()
rho_mle <- c()
pearson_corr_vec <- c()


solve_gamma = function(alpha_mme, mu_mme, S12){
  func <- function(gamma){
    f1 <- as.numeric( as.character( mpmath$hyp1f1(gamma, gamma+1, alpha_mme) ) )
    f11 <-as.numeric( as.character( mpmath$hyp1f1(gamma+1, gamma+2, alpha_mme) )) 
    # f2<- as.numeric( as.character( mpmath$hyp2f2(gamma, gamma, gamma+1, gamma+1, alpha_mme) ) ) 
    return ((beta*alpha_mme/mu_mme*(f1 - (gamma/(1+gamma)*f11)) - S12)^2)
  }
  
  # options <- list("ftol" = 1e-20, "xtol" = 1e-20, "btol" = 1e-20)
  res <- optim(par = gamma, fn = func, method = "L-BFGS-B", lower = 0, upper = 10)
  return(res)
  
  # res = ga(type = "real-valued",  func, lower = 0, upper = gamma*2,popSize = 100,maxiter = 25,keepBest = TRUE, monitor = FALSE)
  # return( res@solution[1,])
  
}

mle_solve = function(x1, x2, n, alpha_mme){
  obj = function(pars){ 
    eta   = pars[1]
    zeta  = pars[2]
    
    gamma = exp(eta)             # ensures gamma >= 0
    delta = -beta + exp(zeta)    # ensures delta >= -beta
    
    loglik = -n*(alpha_mme+beta+delta) + 
      log(alpha_mme)*sum(x1) + 
      beta*sum(gamma/(x1+gamma)) + 
      sum(x2*log(delta + beta*(1 - (gamma/(x1+gamma)))))
    
    return(-loglik)  # optim minimizes
  }
  
  res <- optim(
    # par = c(log(gamma), log(delta+beta)),       # starting values (beta, eta, zeta)
    par = c(-.1, .1),       # starting values (beta, eta, zeta)
    fn  = obj,
    method = "L-BFGS-B",
    lower = c(-10, -10),  # bounds on (beta, eta, zeta)
    upper = c(10,  10)
  )
  
  # Transform back to (beta, gamma, delta)
  eta   <- res$par[1]
  zeta  <- res$par[2]
  gamma <- exp(eta)
  delta <- -beta + exp(zeta)
  return (c(gamma, delta))
}



seed = 1
i = 1
runs = 1e4
while (length(gamma_mme_vec) < runs){
  set.seed(seed)
  x1 <- rpois(n = n, lambda = alpha)
  x2 <- rpois(n =n, lambda = (delta+ beta*(1 - gamma/(x1+gamma))))
  M1 <- sum(x1)/n
  M2<- sum(x2)/n
  S12 <- cov(x1, x2)*(n-1)/n
  S2 <- var(x2) *(n-1)/n

  alpha_mme <- M1
  mu_mme = exp(alpha_mme)
  sol =  solve_gamma(alpha_mme, mu_mme, S12)
  if (sol$convergence !=0){seed = seed+1;next}
  gamma_mme =sol$par
  f1 <- as.numeric( as.character( mpmath$hyp1f1(gamma_mme, gamma_mme+1, alpha_mme) ) )
  f11 <-as.numeric( as.character( mpmath$hyp1f1(gamma_mme+1, gamma_mme+2, alpha_mme) ))
  f2<- as.numeric( as.character( mpmath$hyp2f2(gamma_mme, gamma_mme, gamma_mme+1, gamma_mme+1, alpha_mme) ) )
  delta_mme = M2- beta*(1 - 1/mu_mme*f1)
  
  print (c(S12, beta*alpha_mme/mu_mme*(f1 - (gamma_mme/(gamma_mme+1)*f11))))
  
  
  alpha_mme_vec[i] <- alpha_mme
  gamma_mme_vec[i] <- gamma_mme
  delta_mme_vec[i] <- delta_mme
  
  
  num =  sqrt(alpha_mme)*beta/mu_mme*(f1 - (gamma_mme/(1+gamma_mme)*f11))
  den = sqrt(beta^2*(1/mu_mme*f2 - 1/mu_mme^2*f1^2) + beta*(1 - 1/mu_mme*f1) + delta_mme)
  rho_mme[i] = num/den
  
  #MLE 
  mles = mle_solve(x1, x2, n = n, alpha_mme)
  gamma_mle = mles[1]
  delta_mle = mles[2]
  gamma_mle_vec[i] <- gamma_mle
  delta_mle_vec[i] <- delta_mle
  f1 <- as.numeric( as.character( mpmath$hyp1f1(gamma_mle, gamma_mle+1, alpha_mme) ) )
  f11 <-as.numeric( as.character( mpmath$hyp1f1(gamma_mle+1, gamma_mle+2, alpha_mme) ))
  f2<- as.numeric( as.character( mpmath$hyp2f2(gamma_mle, gamma_mle, gamma_mle+1, gamma_mle+1, alpha_mme) ) )
  
  
  num =  sqrt(alpha_mme)*beta/mu_mme*(f1 - (gamma_mle/(1+gamma_mle)*f11))
  den = sqrt(beta^2*(1/mu_mme*f2 - 1/mu_mme^2*f1^2) + beta*(1 - 1/mu_mme*f1) + delta_mle)
  rho_mle[i] = num/den
  pearson_corr_vec[i] <- cor(x= x1, y= x2, method = "pearson")
  
  
  # print (c(S12, beta*alpha_mme/mu_mme*(f1 - (gamma_mle/(gamma_mle+1)*f11))))
  
  

  
  
  i<- i+1
  seed = seed+1
}


hist(beta_mme_vec_alt, breaks = 1000)

hist(beta_mme_vec, breaks = 1000)
summary(beta_mme_vec)
hist(gamma_mme_vec, breaks =1000)
hist(gamma_mle_vec, breaks =1000)

summary(gamma_mme_vec)
dens = density(gamma_mme_vec, kernel = c('gaussian'))
dens$x[which.max(dens$y)]
hist(delta_mme_vec, breaks = 100)
summary(delta_mme_vec)
summary(gamma_mle_vec)

summary(rho_mme)
summary(pearson_corr_vec)
summary(gamma_mle_vec)

hist(rho_mme, breaks = 100)
hist(pearson_corr_vec, breaks = 100)





# Case II
rm(list = ls())
library(reticulate)
# py_install("mpmath")
mpmath <- import("mpmath")

alpha =1
beta = -2
gamma = 1
delta= 4
n = 100
alpha_mme_vec <- c()
alpha_mle_vec <- c()
beta_mme_vec <- c()
beta_mle_vec <- c()
delta_mme_vec <- c()
delta_mle_vec <- c()
rho_mme <- c()
rho_mle <- c()
pearson_corr_vec <- c()


mle_solve = function(x1, x2, n, alpha_mme){
  obj = function(pars){ 
    beta  = pars[1]
    delta  = pars[2]

    loglik = -n*(alpha_mme+beta+delta) + 
      log(alpha_mme)*sum(x1) + 
      beta*sum(gamma/(x1+gamma)) + 
      sum(x2*log(delta + beta*(1 - (gamma/(x1+gamma)))))
    
    return(-loglik)  # optim minimizes
  }
  
  res <- optim(
    # par = c(beta, log(gamma), log(delta+beta)),       # starting values (beta, eta, zeta)
    par = c(1, 1),       # starting values (beta, eta, zeta)
    fn  = obj,
    method = "L-BFGS-B",
    lower = c(1e-6, 1e-6),  # bounds on (beta, eta, zeta)
    upper = c(10,   10)
  )
  
  # Transform back to (beta, gamma, delta)
  beta  <- res$par[1]
  delta  <- res$par[2]
  return (c(beta, delta))
}


seed = 1
i = 1
runs = 1e4
while (length(beta_mme_vec) < runs){
  set.seed(seed)
  x1 <- rpois(n = n, lambda = alpha)
  x2 <- rpois(n =n, lambda = (delta+ beta*(1 - gamma/(x1+gamma))))
  M1 <- sum(x1)/n
  M2<- sum(x2)/n
  S12 <- cov(x1, x2)*(n-1)/n
  S2 <- var(x2) *(n-1)/n
  
  alpha_mme <- M1
  mu_mme = exp(alpha_mme)
  beta_mme = S12/(1/(alpha_mme*mu_mme)*(mu_mme - alpha_mme - 1))
  delta_mme = M2- beta_mme*(1 - 1/alpha_mme*(1 - 1/mu_mme))
  f1 <- as.numeric( as.character( mpmath$hyp1f1(gamma, gamma+1, alpha_mme) ) )
  f11 <-as.numeric( as.character( mpmath$hyp1f1(gamma+1, gamma+2, alpha_mme) ))
  f2<- as.numeric( as.character( mpmath$hyp2f2(gamma, gamma, gamma+1, gamma+1, alpha_mme) ) )
  alpha_mme_vec[i] <- alpha_mme
  beta_mme_vec[i] <- beta_mme
  delta_mme_vec[i] <- delta_mme
  
  num =  sqrt(alpha_mme)*beta_mme/mu_mme*(f1 - (1/(1+1)*f11))
  den = sqrt(beta_mme^2*(1/mu_mme*f2 - 1/mu_mme^2*f1^2) + beta_mme*(1 - 1/mu_mme*f1) + delta_mme)
  rho_mme[i] = num/den
  
  #MLE 
  mles = mle_solve(x1, x2, n = n, alpha_mme)
  beta_mle = mles[1]
  delta_mle = mles[2]
  beta_mle_vec[i] <- beta_mle
  delta_mle_vec[i] <- delta_mle
  f1 <- as.numeric( as.character( mpmath$hyp1f1(gamma, gamma+1, alpha_mme) ) )
  f11 <-as.numeric( as.character( mpmath$hyp1f1(gamma+1, gamma+2, alpha_mme) ))
  f2<- as.numeric( as.character( mpmath$hyp2f2(gamma, gamma, gamma+1, gamma+1, alpha_mme) ) )
  
  
  num =  sqrt(alpha_mme)*beta_mle/mu_mme*(f1 - (gamma/(1+gamma)*f11))
  den = sqrt(beta_mle^2*(1/mu_mme*f2 - 1/mu_mme^2*f1^2) + beta_mle*(1 - 1/mu_mme*f1) + delta_mle)
  rho_mle[i] = num/den
  
  pearson_corr_vec[i] <- cor(x= x1, y= x2, method = "pearson")
  i<- i+1
  seed = seed+1
}



hist(beta_mme_vec, breaks = 1000)
summary(beta_mme_vec)
hist(gamma_mme_vec, breaks =1000)

hist(delta_mme_vec, breaks = 100)
summary(delta_mme_vec)


summary(rho_mme)
summary(pearson_corr_vec)

hist(rho_mme, breaks = 100)
hist(pearson_corr_vec, breaks = 100)



# Case IV
rm(list = ls())
library(reticulate)
# py_install("mpmath")
mpmath <- import("mpmath")

alpha =2
beta = -1
gamma = 1
delta= 4
n = 100
alpha_mme_vec <- c()
delta_mme_vec <- c()
delta_mle_vec <- c()
rho_mme <- c()
rho_mle <- c()
pearson_corr_vec <- c()


mle_solve = function(x1, x2, n, alpha_mme){
  obj = function(pars){ 
    zeta  = pars[1]
    delta = -beta + exp(zeta)    # ensures delta >= -beta
    
    loglik = -n*(alpha_mme+beta+delta) + 
      log(alpha_mme)*sum(x1) + 
      beta*sum(gamma/(x1+gamma)) + 
      sum(x2*log(delta + beta*(1 - (gamma/(x1+gamma)))))
    
    return(-loglik)  # optim minimizes
  }
  
  res <- optim(
    # par = c(beta, log(gamma), log(delta+beta)),       # starting values (beta, eta, zeta)
    par = c(0),       # starting values (beta, eta, zeta)
    fn  = obj,
    method = "L-BFGS-B",
    lower = c(-10),  # bounds on (beta, eta, zeta)
    upper = c(10)
  )
  zeta  <- res$par[1]
  delta <- -beta + exp(zeta)
  return (c(delta))
}

seed = 1
i = 1
runs = 1e4
while (length(delta_mme_vec) < runs){
  set.seed(seed)
  x1 <- rpois(n = n, lambda = alpha)
  x2 <- rpois(n =n, lambda = (delta+ beta*(1 - gamma/(x1+gamma))))
  M1 <- sum(x1)/n
  M2<- sum(x2)/n
  S12 <- cov(x1, x2)*(n-1)/n
  S2 <- var(x2) *(n-1)/n
  
  alpha_mme <- M1
  mu_mme = exp(alpha_mme)
  delta_mme = M2- beta*(1 - 1/alpha_mme+ 1/(mu_mme*alpha_mme))
  f1 <- as.numeric( as.character( mpmath$hyp1f1(gamma, gamma+1, alpha_mme) ) )
  f11 <-as.numeric( as.character( mpmath$hyp1f1(gamma+1, gamma+2, alpha_mme) ))
  f2<- as.numeric( as.character( mpmath$hyp2f2(gamma, gamma, gamma+1, gamma+1, alpha_mme) ) )
  alpha_mme_vec[i] <- alpha_mme
  delta_mme_vec[i] <- delta_mme
  
  num =  sqrt(alpha_mme)*beta/mu_mme*(f1 - (1/(1+1)*f11))
  den = sqrt(beta^2*(1/mu_mme*f2 - 1/mu_mme^2*f1^2) + beta*(1 - 1/mu_mme*f1) + delta_mme)
  rho_mme[i] = num/den
  
  #MLE 
  mles = mle_solve(x1, x2, n = n, alpha_mme)
  delta_mle = mles[1]
  delta_mle_vec[i] <- delta_mle
  f1 <- as.numeric( as.character( mpmath$hyp1f1(gamma, gamma+1, alpha_mme) ) )
  f11 <-as.numeric( as.character( mpmath$hyp1f1(gamma+1, gamma+2, alpha_mme) ))
  f2<- as.numeric( as.character( mpmath$hyp2f2(gamma, gamma, gamma+1, gamma+1, alpha_mme) ) )
  
  
  num =  sqrt(alpha_mme)*beta/mu_mme*(f1 - (gamma/(1+gamma)*f11))
  den = sqrt(beta^2*(1/mu_mme*f2 - 1/mu_mme^2*f1^2) + beta*(1 - 1/mu_mme*f1) + delta_mle)
  rho_mle[i] = num/den
  
  pearson_corr_vec[i] <- cor(x= x1, y= x2, method = "pearson")
  i<- i+1
  seed = seed+1
  
}



hist(beta_mme_vec, breaks = 1000)
summary(beta_mme_vec)
hist(gamma_mme_vec, breaks =1000)

hist(delta_mme_vec, breaks = 100)
summary(delta_mme_vec)


summary(rho_mme)
summary(pearson_corr_vec)

hist(rho_mme, breaks = 100)
hist(pearson_corr_vec, breaks = 100)


#========================================================================
#                 Application
#========================================================================
rm(list = ls())
library(Brobdingnag)


data<- read.csv("C:/Users/Jared/OneDrive - University of Cape Town/Documents/WST 795-20230920T125726Z-001/WST 795/Arr/Health_Retirement_Study.csv")
xvec <- rep(data$x, data$F)
yvec <- rep(data$y, data$F)
data <- data.frame(x = xvec, y  = yvec)


# data <- read_excel("C:/Users/jared/Documents/WST 795/Arr/insurance.xlsx")
# 
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
x1<- data$x
x2<- data$y
# Mirrored
# x1<- data$y
# x2<- data$x
M1 <- mean(x1)
M2<- mean(x2)
S12 <- cov(x1, x2)* (n-1)/n
M12 <- sum(x1*x2 )/ n
S2<-var(x2)*(n-1)/n 
cor(x= x1, y= x2, method = "pearson")



M1^2/(exp(M1) - 1)
S12^2/(S2 - M2)
M1

M1/S12

term <- 1
for (i in 1:n){
  term <- as.brob(term* factorial(x1[i]) *factorial(x2[i]))
}
h <- log(1/term)






















Accident_data <- read.csv('C:/Users/Jared/OneDrive - University of Cape Town/Documents/WST 795-20230920T125726Z-001/WST 795/Arr/Accident_data.csv')
n<- nrow(Accident_data)
# x2<- Accident_data$x
# x1<- Accident_data$y
# Mirrored
x2<- Accident_data$y
x1<- Accident_data$x
M1 <- mean(x1)
M2<- mean(x2)
S12 <- cov(x1, x2)* (n-1)/n
M12 <- sum(x1*x2 )/ n
S2<-var(x2)*(n-1)/n 
cor(x= x1, y= x2, method = "pearson")

M1^2/(exp(M1) - 1)
S12^2/(S2 - M2)
M1

M1/S12


term <- 1
for (i in 1:n){
  term <- term* factorial(x1[i]) *factorial(x2[i])
}
h = -log(term)




# Full 

library(GenSA)
library(reticulate)
# py_install("mpmath")
mpmath <- import("mpmath")



mle_solve = function(x1, x2, n, alpha_mme){
  obj = function(pars){ 
    beta  = pars[1]
    dum1   = pars[2]
    dum2  = pars[3]
    dum3 = pars[4]
    
    gamma = exp(dum1)             # ensures gamma >= 0
    delta = exp(dum2)    
    eta = exp(dum3)
    
    
    loglik = -n*(alpha_mme+beta+delta) + 
      log(alpha_mme)*sum(x1) + 
      beta*sum((gamma/(x1+gamma))^eta) + 
      sum(x2*log(delta + beta*(1 - (gamma/(x1+gamma))^eta)))
    
    return(-loglik)  # optim minimizes
  }
  
  
  
  
  res <- optim(
    # par = c(beta, log(gamma), log(delta+beta)),       # starting values (beta, eta, zeta)
    par = c(.1, 0, 0, 0),       # starting values (beta, eta, zeta)
    fn  = obj,
    method = "L-BFGS-B",
    # lower = c(0, -10, -10, -10),  # bounds on (beta, eta, zeta)
    # upper = c(10,   10,  10, 10)
    lower = c(0, rep(-Inf, 3)),
    upper = rep(Inf, 4)
  )
  print (res)
  # pars <- GenSA(par = NULL, fn = obj, lower =c(-50, -10, -10), upper = c(0,   10,  10) )
  # print (pars)
  # Transform back to (beta, gamma, delta)
  beta  <- res$par[1]
  dum1   <- res$par[2]
  dum2  <- res$par[3]
  dum3 = res$par[4]
  gamma <- exp(dum1)
  delta <-  exp(dum2)
  eta = exp(dum3)
  return (c(beta, gamma, delta, eta))
}


alpha_mme <- M1
mu_mme = exp(alpha_mme)


#MLE 
mles = mle_solve(x1, x2, n = n, alpha_mme)
mles
beta_mle = mles[1]
gamma_mle = mles[2]
delta_mle = mles[3]
eta_mle = mles[4]

# f1 <- as.numeric( as.character( mpmath$hyp1f1(gamma_mle, gamma_mle+1, alpha_mme) ) )
# f11 <-as.numeric( as.character( mpmath$hyp1f1(gamma_mle+1, gamma_mle+2, alpha_mme) ))
# f2<- as.numeric( as.character( mpmath$hyp2f2(gamma_mle, gamma_mle, gamma_mle+1, gamma_mle+1, alpha_mme) ) )
# 
# num =  sqrt(alpha_mme)*beta_mle/mu_mme*(f1 - (gamma_mle/(1+gamma_mle)^eta_mle*f11))
# den = sqrt(beta_mle^2*(1/mu_mme*f2 - 1/mu_mme^2*f1^2) + beta_mle*(1 - 1/mu_mme*f1) + delta_mle)
# rho_mle = num/den
#AIC
k<- 5
l_full = -n*(alpha_mme+beta_mle+delta_mle) + log(alpha_mme)*sum(x1) + beta_mle*sum((gamma_mle/(x1+gamma_mle))^eta_mle) + sum(x2*log(delta_mle + beta_mle*(1 - (gamma_mle/(x1+gamma_mle))^eta_mle))) +h
aic <- 2*k - 2*l_full
bic <- k*log(n) - 2*l_full  































# Lomax(1, gamma)
library(GenSA)
library(reticulate)
# py_install("mpmath")
mpmath <- import("mpmath")


solve_gamma = function(alpha_mme, mu_mme,  S2, M2, S12){
  func <- function(gamma){
    f1 <- as.numeric( as.character( mpmath$hyp1f1(gamma, gamma+1, alpha_mme) ) )
    f11 <-as.numeric( as.character( mpmath$hyp1f1(gamma+1, gamma+2, alpha_mme) )) 
    f2<- as.numeric( as.character( mpmath$hyp2f2(gamma, gamma, gamma+1, gamma+1, alpha_mme) ) ) 
    return ((alpha_mme^2*(f1- (gamma/(gamma+1)*f11))^2/(mu_mme*f2 - f1^2) - S12^2/(S2-M2))^2)
  }
  
  # options <- list("ftol" = 1e-20, "xtol" = 1e-20, "btol" = 1e-20)
  res <- optim(par = 0.1, fn = func, method = "L-BFGS-B", lower = 0, upper = Inf)
  return(res$par)
  
  # res = ga(type = "real-valued",  func, lower = 0, upper = gamma*2,popSize = 100,maxiter = 25,keepBest = TRUE, monitor = FALSE)
  # return( res@solution[1,])
  
}
# solve_gamma(alpha_mme, mu_mme,  S2, M2, S12)


solve_nu_trans <- function(alpha_mme, mu_mme, S12, S2, M2, trans_start) {
  rhs <- S12^2 / (S2 - M2)
  obj <- function(trans) {
    num <- alpha_mme^2 * trans
    denom <- mu_mme^trans - 1
    return(-(num / denom - rhs)^2)  # explicit return stops printing
  }
  
  # res <- optim(par = trans_start, fn = obj, method = "L-BFGS-B", lower = 0, upper = 1)
  # return(res$par)
  # res <- suppressWarnings(GenSA(
  #   par    = trans_start,
  #   fn     = obj,
  #   lower  = 0,
  #   upper  = 1,
  #   control = list(max.call = 5000, verbose = FALSE) # you can increase if needed
  # ))
  res = ga(type = "real-valued",  obj, lower = 0, upper = 1,popSize = 100,maxiter = 25,keepBest = TRUE, monitor= FALSE)
  
  
  return( res@solution[1,])
}




mle_solve = function(x1, x2, n, alpha_mme){
  obj = function(pars){ 
    beta  = pars[1]
    eta   = pars[2]
    zeta  = pars[3]
    
    gamma = exp(eta)             # ensures gamma >= 0
    delta = exp(zeta)    # ensures delta >= -beta
    
    loglik = -n*(alpha_mme+beta+delta) + 
      log(alpha_mme)*sum(x1) + 
      beta*sum(gamma/(x1+gamma)) + 
      sum(x2*log(delta + beta*(1 - (gamma/(x1+gamma)))))
    
    return(-loglik)  # optim minimizes
  }
  
  res <- optim(
    # par = c(beta, log(gamma), log(delta+beta)),       # starting values (beta, eta, zeta)
    par = c(-1, 0, 0),       # starting values (beta, eta, zeta)
    fn  = obj,
    method = "L-BFGS-B",
    # lower = c(0, -10, -10),  # bounds on (beta, eta, zeta)
    # upper = c(10,   10,  10)
    lower = c(0, rep(-Inf, 3)),
    upper = c(Inf, 4)
  )
  
  # pars <- GenSA(par = NULL, fn = obj, lower =c(-50, -10, -10), upper = c(0,   10,  10) )
  # print (pars)
  # Transform back to (beta, gamma, delta)
  beta  <- res$par[1]
  eta   <- res$par[2]
  zeta  <- res$par[3]
  gamma <- exp(eta)
  delta <-  exp(zeta)
  return (c(beta, gamma, delta))
}


alpha_mme <- M1
mu_mme = exp(alpha_mme)
# gamma_mme = solve_gamma(alpha_mme, mu_mme,  S2, M2, S12)
# f1 <- as.numeric( as.character( mpmath$hyp1f1(gamma_mme, gamma_mme+1, alpha_mme) ) )
# f11 <-as.numeric( as.character( mpmath$hyp1f1(gamma_mme+1, gamma_mme+2, alpha_mme) ))
# f2<- as.numeric( as.character( mpmath$hyp2f2(gamma_mme, gamma_mme, gamma_mme+1, gamma_mme+1, alpha_mme) ) )
# 
# 
# beta_mme = S12/(alpha_mme/mu_mme*(f1 - (gamma_mme/(gamma_mme+1)*f11)))
# delta_mme = M2- beta_mme*(1 - 1/mu_mme*f1)
# 
# num =  sqrt(alpha_mme)*beta_mme/mu_mme*(f1 - (gamma_mme/(1+gamma_mme)*f11))
# den = sqrt(beta_mme^2*(1/mu_mme*f2 - 1/mu_mme^2*f1^2) + beta_mme*(1 - 1/mu_mme*f1) + delta_mme)
# rho_mme = num/den
    
    
#MLE 
mles = mle_solve(x1, x2, n = n, alpha_mme)
mles
beta_mle = mles[1]
gamma_mle = mles[2]
delta_mle = mles[3]

f1 <- as.numeric( as.character( mpmath$hyp1f1(gamma_mle, gamma_mle+1, alpha_mme) ) )
f11 <-as.numeric( as.character( mpmath$hyp1f1(gamma_mle+1, gamma_mle+2, alpha_mme) ))
f2<- as.numeric( as.character( mpmath$hyp2f2(gamma_mle, gamma_mle, gamma_mle+1, gamma_mle+1, alpha_mme) ) )

num =  sqrt(alpha_mme)*beta_mle/mu_mme*(f1 - (gamma_mle/(1+gamma_mle)*f11))
den = sqrt(beta_mle^2*(1/mu_mme*f2 - 1/mu_mme^2*f1^2) + beta_mle*(1 - 1/mu_mme*f1) + delta_mle)
rho_mle = num/den
#AIC
k<- 4
l = -n*(alpha_mme+beta_mle+delta_mle) + log(alpha_mme)*sum(x1) + beta_mle*sum(gamma_mle/(x1+gamma_mle)) + sum(x2*log(delta_mle + beta_mle*(1 - (gamma_mle/(x1+gamma_mle))))) +h
aic <- 2*k - 2*l
bic <- k*log(n) - 2*l  
neg2loglam <- -2*(l - l_full)
critical <- qchisq(p =0.95, df = 1)




# Case I
library(reticulate)
# py_install("mpmath")
mpmath <- import("mpmath")


solve_gamma = function(alpha_mme, mu_mme, S12){
  func <- function(dum1){
    gamma = exp(dum1)
    f1 <- as.numeric( as.character( mpmath$hyp1f1(gamma, gamma+1, alpha_mme) ) )
    f11 <-as.numeric( as.character( mpmath$hyp1f1(gamma+1, gamma+2, alpha_mme) )) 
    # f2<- as.numeric( as.character( mpmath$hyp2f2(gamma, gamma, gamma+1, gamma+1, alpha_mme) ) ) 
    return ((1*alpha_mme/mu_mme*(f1 - (gamma/(1+gamma)*f11)) - S12)^2)
  }
  
  # options <- list("ftol" = 1e-20, "xtol" = 1e-20, "btol" = 1e-20)
  res <- optim(par = 0.3, fn = func, method = "L-BFGS-B", lower = -Inf, upper = Inf)
  print (res)
  dum1 = res$par
  gamma = exp(dum1)
  return(gamma)
  
  # res = ga(type = "real-valued",  func, lower = 0, upper = gamma*2,popSize = 100,maxiter = 25,keepBest = TRUE, monitor = FALSE)
  # return( res@solution[1,])
}

mle_solve = function(x1, x2, n, alpha_mme){
  beta = 1
  obj = function(pars){ 
    eta   = pars[1]
    zeta  = pars[2]
    
    gamma = exp(eta)             # ensures gamma >= 0
    delta = exp(zeta)    # ensures delta >= -beta
    
    loglik = -n*(alpha_mme+beta+delta) + 
      log(alpha_mme)*sum(x1) + 
      beta*sum(gamma/(x1+gamma)) + 
      sum(x2*log(delta + beta*(1 - (gamma/(x1+gamma)))))
    
    return(-loglik)  # optim minimizes
  }
  
  res <- optim(
    # par = c(log(gamma), log(delta+beta)),       # starting values (beta, eta, zeta)
    par = c(-.1, .1),       # starting values (beta, eta, zeta)
    fn  = obj,
    method = "L-BFGS-B",
    # lower = c(-10, -10),  # bounds on (beta, eta, zeta)
    # upper = c(10,  10)
    lower = rep(-Inf, 2),
    upper = rep(Inf, 2)
  )
  res
  # Transform back to (beta, gamma, delta)
  eta   <- res$par[1]
  zeta  <- res$par[2]
  gamma <- exp(eta)
  delta <-exp(zeta)
  return (c(gamma, delta))
}

alpha_mme <- M1
mu_mme = exp(alpha_mme)
sol =  solve_gamma(alpha_mme, mu_mme, S12)
gamma_mme =sol
f1 <- as.numeric( as.character( mpmath$hyp1f1(gamma_mme, gamma_mme+1, alpha_mme) ) )
f11 <-as.numeric( as.character( mpmath$hyp1f1(gamma_mme+1, gamma_mme+2, alpha_mme) ))
f2<- as.numeric( as.character( mpmath$hyp2f2(gamma_mme, gamma_mme, gamma_mme+1, gamma_mme+1, alpha_mme) ) )
delta_mme = M2- 1*(1 - 1/mu_mme*f1)
num =  sqrt(alpha_mme)*1/mu_mme*(f1 - (gamma_mme/(1+gamma_mme)*f11))
den = sqrt(1^2*(1/mu_mme*f2 - 1/mu_mme^2*f1^2) + 1*(1 - 1/mu_mme*f1) + delta_mme)
rho_mme = num/den
  
#MLE 
mles = mle_solve(x1, x2, n = n, alpha_mme)
mles
gamma_mle = mles[1]
delta_mle = mles[2]
f1 <- as.numeric( as.character( mpmath$hyp1f1(gamma_mle, gamma_mle+1, alpha_mme) ) )
f11 <-as.numeric( as.character( mpmath$hyp1f1(gamma_mle+1, gamma_mle+2, alpha_mme) ))
f2<- as.numeric( as.character( mpmath$hyp2f2(gamma_mle, gamma_mle, gamma_mle+1, gamma_mle+1, alpha_mme) ) )
num =  sqrt(alpha_mme)*1/mu_mme*(f1 - (gamma_mle/(1+gamma_mle)*f11))
den = sqrt(1^2*(1/mu_mme*f2 - 1/mu_mme^2*f1^2) + 1*(1 - 1/mu_mme*f1) + delta_mle)
rho_mle = num/den

#AIC
k<- 3
l = -n*(alpha_mme+1+delta_mle) + log(alpha_mme)*sum(x1) + 1*sum(gamma_mle/(x1+gamma_mle)) + sum(x2*log(delta_mle + 1*(1 - (gamma_mle/(x1+gamma_mle))))) +h
aic <- 2*k - 2*l
bic <- k*log(n) - 2*l  

neg2loglam <- -2*(l - l_full)
critical <- qchisq(p =0.95, df = 2)





# Case II
library(reticulate)
# py_install("mpmath")
mpmath <- import("mpmath")



mle_solve = function(x1, x2, n, alpha_mme){
  gamma = 1
  obj = function(pars){ 
    dum1 = pars[1]
    dum2 = pars[2]
    beta  = exp(dum1)
    delta  =  exp(dum2)
    
    loglik = -n*(alpha_mme+beta+delta) + 
      log(alpha_mme)*sum(x1) + 
      beta*sum(gamma/(x1+gamma)) + 
      sum(x2*log(delta + beta*(1 - (gamma/(x1+gamma)))))
    
    return(-loglik)  # optim minimizes
  }
  
  res <- optim(
    # par = c(beta, log(gamma), log(delta+beta)),       # starting values (beta, eta, zeta)
    par = c(.1, .1),       # starting values (beta, eta, zeta)
    fn  = obj,
    method = "L-BFGS-B",
    # lower = c(1e-6, 1e-6),  # bounds on (beta, eta, zeta)
    # upper = c(10,   10)
    lower = rep(-Inf, 2),
    upper = rep(Inf, 2)
  )
  
  # Transform back to (beta, gamma, delta)
  dum1  <- res$par[1]
  dum2  <- res$par[2]
  beta  = exp(dum1)
  delta  =  exp(dum2)
  return (c(beta, delta))
}

alpha_mme <- M1
mu_mme = exp(alpha_mme)
beta_mme = S12/(1/(alpha_mme*mu_mme)*(mu_mme - alpha_mme - 1))
delta_mme = M2- beta_mme*(1 - 1/alpha_mme*(1 - 1/mu_mme))
f1 <- as.numeric( as.character( mpmath$hyp1f1(1, 1+1, alpha_mme) ) )
f11 <-as.numeric( as.character( mpmath$hyp1f1(1+1, 1+2, alpha_mme) ))
f2<- as.numeric( as.character( mpmath$hyp2f2(1, 1, 1+1, 1+1, alpha_mme) ) )
num =  sqrt(alpha_mme)*beta_mme/mu_mme*(f1 - (1/(1+1)*f11))
den = sqrt(beta_mme^2*(1/mu_mme*f2 - 1/mu_mme^2*f1^2) + beta_mme*(1 - 1/mu_mme*f1) + delta_mme)
rho_mme = num/den
#MLE 
mles = mle_solve(x1, x2, n = n, alpha_mme)
beta_mle = mles[1]
delta_mle = mles[2]
f1 <- as.numeric( as.character( mpmath$hyp1f1(1, 1+1, alpha_mme) ) )
f11 <-as.numeric( as.character( mpmath$hyp1f1(1+1, 1+2, alpha_mme) ))
f2<- as.numeric( as.character( mpmath$hyp2f2(1, 1, 1+1, 1+1, alpha_mme) ) )
num =  sqrt(alpha_mme)*beta_mle/mu_mme*(f1 - (1/(1+1)*f11))
den = sqrt(beta_mle^2*(1/mu_mme*f2 - 1/mu_mme^2*f1^2) + beta_mle*(1 - 1/mu_mme*f1) + delta_mle)
rho_mle = num/den

#AIC
k<- 3
l = -n*(alpha_mme+beta_mle+delta_mle) + log(alpha_mme)*sum(x1) + beta_mle*sum(1/(x1+1)) + sum(x2*log(delta_mle + beta_mle*(1 - (1/(x1+1))))) +h
aic <- 2*k - 2*l
bic <- k*log(n) - 2*l  

neg2loglam <- -2*(l - l_full)
critical <- qchisq(p =0.95, df = 2)




# Case III

library(GenSA)
library(GA)
library(reticulate)
# py_install("mpmath")
mpmath <- import("mpmath")


solve_gamma = function(alpha_mme, mu_mme,  S2, M2, S12){
  func <- function(gamma){
    f1 <- as.numeric( as.character( mpmath$hyp1f1(gamma, gamma+1, alpha_mme) ) )
    f11 <-as.numeric( as.character( mpmath$hyp1f1(gamma+1, gamma+2, alpha_mme) )) 
    f2<- as.numeric( as.character( mpmath$hyp2f2(gamma, gamma, gamma+1, gamma+1, alpha_mme) ) ) 
    return ((alpha_mme^2*(f1- (gamma/(gamma+1)*f11))^2/(mu_mme*f2 - f1^2) - S12^2/(S2-M2))^2)
  }
  
  # options <- list("ftol" = 1e-20, "xtol" = 1e-20, "btol" = 1e-20)
  res <- optim(par = 0.1, fn = func, method = "L-BFGS-B", lower = 0, upper = Inf)
  return(res$par)
  
  # res = ga(type = "real-valued",  func, lower = 0, upper = gamma*2,popSize = 100,maxiter = 25,keepBest = TRUE, monitor = FALSE)
  # return( res@solution[1,])
  
}
# solve_gamma(alpha_mme, mu_mme,  S2, M2, S12)


solve_nu_trans <- function(alpha_mme, mu_mme, S12, S2, M2, trans_start) {
  rhs <- S12^2 / (S2 - M2)
  obj <- function(trans) {
    num <- alpha_mme^2 * trans
    denom <- mu_mme^trans - 1
    return(-(num / denom - rhs)^2)  # explicit return stops printing
  }
  
  # res <- optim(par = trans_start, fn = obj, method = "L-BFGS-B", lower = 0, upper = 1)
  # return(res$par)
  # res <- suppressWarnings(GenSA(
  #   par    = trans_start,
  #   fn     = obj,
  #   lower  = 0,
  #   upper  = 1,
  #   control = list(max.call = 5000, verbose = FALSE) # you can increase if needed
  # ))
  res = ga(type = "real-valued",  obj, lower = 0, upper = 1,popSize = 100,maxiter = 25,keepBest = TRUE, monitor= FALSE)
  
  
  return( res@solution[1,])
}



mle_solve = function(x1, x2, n, alpha_mme){
  delta = 0
  # obj = function(pars){ 
  #   dum_beta  = pars[1]
  #   eta   = pars[2]
  #   
  #   beta = exp(dum_beta)
  #   gamma = exp(eta)             # ensures gamma >= 0
  # 
  #   loglik = -n*(alpha_mme+beta+delta) + 
  #     log(alpha_mme)*sum(x1) + 
  #     beta*sum(gamma/(x1+gamma)) + 
  #     sum(x2*log(delta + beta*(1 - (gamma/(x1+gamma)))), na.rm = TRUE)
  #   
  #   return(loglik)  # optim minimizes
  # }
  
  obj <- function(pars){
    dum1 <- pars[1]
    dum2  <- pars[2]
    beta = exp(dum1)
    gamma <- exp(dum2)                           # gamma >= 0
    
    n  <- length(x1)
    lam1 <- alpha_mme                           # X1 mean
    lam2 <- beta * (1 - (gamma/(x1+gamma)))       # delta = 0
    
    ll1 <- sum(dpois(x1, lam1, log = TRUE))
    ll2 <- sum(dpois(x2, lam2, log = TRUE))     # dpois(0,0) = 0; handles x1=0,x2=0
    ll  <- ll1 + ll2
    if (!is.finite(ll)) return(1e12)
    -ll
  }
  res <- optim(
    # par = c(beta, log(gamma), log(delta+beta)),       # starting values (beta, eta, zeta)
    par = c(0.5, 0.5),       # starting values (beta, eta, zeta)
    fn  = obj,
    method = "L-BFGS-B",
    # lower = c(1e-4, -10),  # bounds on (beta, eta, zeta)
    # upper = c(10,   10)
    lower = rep(-Inf, 2),
    upper = rep(Inf, 2)
  )
  # res = ga(type = "real-valued",  obj, lower = c(-10, -10), upper = c(10, 10),popSize = 100,maxiter = 100,keepBest = TRUE)
  # print (res@solution[1, ])
  # dum_beta = res@solution[1, 1]
  # eta = res@solution[1, 2]
  # beta = exp(dum_beta)
  # gamma = exp(eta)
  # pars <- GenSA(par = c(0.1, 0.1), fn = obj, lower =c(0, -10), upper = c(20,  10) )
  # print (pars)
  # Transform back to (beta, gamma, delta)
  dum1  <- res$par[1]
  dum2   <- res$par[2]
  beta = exp(dum1)
  gamma <- exp(dum2)
  return (c(beta, gamma))
}


alpha_mme <- M1
mu_mme = exp(alpha_mme)
# gamma_mme = solve_gamma(alpha_mme, mu_mme,  S2, M2, S12)
# f1 <- as.numeric( as.character( mpmath$hyp1f1(gamma_mme, gamma_mme+1, alpha_mme) ) )
# f11 <-as.numeric( as.character( mpmath$hyp1f1(gamma_mme+1, gamma_mme+2, alpha_mme) ))
# f2<- as.numeric( as.character( mpmath$hyp2f2(gamma_mme, gamma_mme, gamma_mme+1, gamma_mme+1, alpha_mme) ) )
# 
# 
# beta_mme = S12/(alpha_mme/mu_mme*(f1 - (gamma_mme/(gamma_mme+1)*f11)))
# delta_mme = M2- beta_mme*(1 - 1/mu_mme*f1)
# 
# num =  sqrt(alpha_mme)*beta_mme/mu_mme*(f1 - (gamma_mme/(1+gamma_mme)*f11))
# den = sqrt(beta_mme^2*(1/mu_mme*f2 - 1/mu_mme^2*f1^2) + beta_mme*(1 - 1/mu_mme*f1) + delta_mme)
# rho_mme = num/den


#MLE 
mles = mle_solve(x1, x2, n = n, alpha_mme)
mles
beta_mle = mles[1]
gamma_mle = mles[2]

f1 <- as.numeric( as.character( mpmath$hyp1f1(gamma_mle, gamma_mle+1, alpha_mme) ) )
f11 <-as.numeric( as.character( mpmath$hyp1f1(gamma_mle+1, gamma_mle+2, alpha_mme) ))
f2<- as.numeric( as.character( mpmath$hyp2f2(gamma_mle, gamma_mle, gamma_mle+1, gamma_mle+1, alpha_mme) ) )

num =  sqrt(alpha_mme)*beta_mle/mu_mme*(f1 - (gamma_mle/(1+gamma_mle)*f11))
den = sqrt(beta_mle^2*(1/mu_mme*f2 - 1/mu_mme^2*f1^2) + beta_mle*(1 - 1/mu_mme*f1) + 0)
rho_mle = num/den
#AIC
k<- 3
l= -n*(alpha_mme+beta_mle+0) + log(alpha_mme)*sum(x1) + beta_mle*sum(gamma_mle/(x1+gamma_mle)) + sum(x2*log(0 + beta_mle*(1 - (gamma_mle/(x1+gamma_mle)))), na.rm = TRUE) +h
aic <- 2*k - 2*l
bic <- k*log(n) - 2*l  
neg2loglam <- -2*(l - l_full)
critical <- qchisq(p =0.95, df = 2)



# Case IV
library(reticulate)
# py_install("mpmath")
mpmath <- import("mpmath")


mle_solve = function(x1, x2, n, alpha_mme){
  beta= 1
  gamma = 1
  obj = function(pars){ 
    dum1  = pars[1]
    delta =  exp(dum1)    # ensures delta >= -beta
    
    loglik = -n*(alpha_mme+beta+delta) + 
      log(alpha_mme)*sum(x1) + 
      beta*sum(gamma/(x1+gamma)) + 
      sum(x2*log(delta + beta*(1 - (gamma/(x1+gamma)))))
    
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
  dum1  <- res$par[1]
  delta <- exp(dum1)
  return (c(delta))
}


alpha_mme <- M1
mu_mme = exp(alpha_mme)
delta_mme = M2- 1*(1 - 1/alpha_mme+ 1/(alpha_mme*mu_mme))
f1 <- as.numeric( as.character( mpmath$hyp1f1(1, 1+1, alpha_mme) ) )
f11 <-as.numeric( as.character( mpmath$hyp1f1(1+1, 1+2, alpha_mme) ))
f2<- as.numeric( as.character( mpmath$hyp2f2(1, 1, 1+1, 1+1, alpha_mme) ) )
num =  sqrt(alpha_mme)*1/mu_mme*(f1 - (1/(1+1)*f11))
den = sqrt(1^2*(1/mu_mme*f2 - 1/mu_mme^2*f1^2) + 1*(1 - 1/mu_mme*f1) + delta_mme)
rho_mme = num/den
  
# MLE 
mles = mle_solve(x1, x2, n = n, alpha_mme)
delta_mle = mles[1]
f1 <- as.numeric( as.character( mpmath$hyp1f1(1, 1+1, alpha_mme) ) )
f11 <-as.numeric( as.character( mpmath$hyp1f1(1+1, 1+2, alpha_mme) ))
f2<- as.numeric( as.character( mpmath$hyp2f2(1, 1, 1+1, 1+1, alpha_mme) ) )
num =  sqrt(alpha_mme)*1/mu_mme*(f1 - (1/(1+1)*f11))
den = sqrt(1^2*(1/mu_mme*f2 - 1/mu_mme^2*f1^2) + 1*(1 - 1/mu_mme*f1) + delta_mle)
rho_mle = num/den
# AIC
k = 2
l = -n*(alpha_mme+1+delta_mle) + log(alpha_mme)*sum(x1) + 1*sum(1/(x1+1)) + sum(x2*log(delta_mle + 1*(1 - (1/(x1+1)))))+h
aic <- 2*k - 2*l
bic <- k*log(n) - 2*l  

neg2loglam <- -2*(l - l_full)
critical <- qchisq(p =0.95, df = 3)






# Case V

library(GenSA)
library(GA)
library(reticulate)
# py_install("mpmath")
mpmath <- import("mpmath")





mle_solve = function(x1, x2, n, alpha_mme){
  delta = 0
  gamma = 1
  # obj = function(pars){ 
  #   dum_beta  = pars[1]
  #   eta   = pars[2]
  #   
  #   beta = exp(dum_beta)
  #   gamma = exp(eta)             # ensures gamma >= 0
  # 
  #   loglik = -n*(alpha_mme+beta+delta) + 
  #     log(alpha_mme)*sum(x1) + 
  #     beta*sum(gamma/(x1+gamma)) + 
  #     sum(x2*log(delta + beta*(1 - (gamma/(x1+gamma)))), na.rm = TRUE)
  #   
  #   return(loglik)  # optim minimizes
  # }
  
  obj <- function(pars){
    dum1 <- pars[1]
    beta = exp(dum1)

    n  <- length(x1)
    lam1 <- alpha_mme                           # X1 mean
    lam2 <- beta * (1 - (gamma/(x1+gamma)))       # delta = 0
    
    ll1 <- sum(dpois(x1, lam1, log = TRUE))
    ll2 <- sum(dpois(x2, lam2, log = TRUE))     # dpois(0,0) = 0; handles x1=0,x2=0
    ll  <- ll1 + ll2
    if (!is.finite(ll)) return(1e12)
    -ll
  }
  res <- optim(
    # par = c(beta, log(gamma), log(delta+beta)),       # starting values (beta, eta, zeta)
    par = c(0.5),       # starting values (beta, eta, zeta)
    fn  = obj,
    method = "L-BFGS-B",
    # lower = c(1e-4, -10),  # bounds on (beta, eta, zeta)
    # upper = c(10,   10)
    lower = rep(-Inf, 1),
    upper = rep(Inf, 1)
  )
  # res = ga(type = "real-valued",  obj, lower = c(-10, -10), upper = c(10, 10),popSize = 100,maxiter = 100,keepBest = TRUE)
  # print (res@solution[1, ])
  # dum_beta = res@solution[1, 1]
  # eta = res@solution[1, 2]
  # beta = exp(dum_beta)
  # gamma = exp(eta)
  # pars <- GenSA(par = c(0.1, 0.1), fn = obj, lower =c(0, -10), upper = c(20,  10) )
  # print (pars)
  # Transform back to (beta, gamma, delta)
  dum1  <- res$par[1]
  beta = exp(dum1)
  return (c(beta))
}


alpha_mme <- M1
mu_mme = exp(alpha_mme)
f1 <- as.numeric( as.character( mpmath$hyp1f1(1, 1+1, alpha_mme) ) )
f11 <-as.numeric( as.character( mpmath$hyp1f1(1+1, 1+2, alpha_mme) ))
f2<- as.numeric( as.character( mpmath$hyp2f2(1, 1, 1+1, 1+1, alpha_mme) ) )


beta_mme = M2/ (1 - 1/mu_mme*f1)

num =  sqrt(alpha_mme)*beta_mme/mu_mme*(f1 - (1/(1+1)*f11))
den = sqrt(beta_mme^2*(1/mu_mme*f2 - 1/mu_mme^2*f1^2) + beta_mme*(1 - 1/mu_mme*f1) + 0)
rho_mme = num/den


#MLE 
mles = mle_solve(x1, x2, n = n, alpha_mme)
mles
beta_mle = mles[1]

f1 <- as.numeric( as.character( mpmath$hyp1f1(1, 1+1, alpha_mme) ) )
f11 <-as.numeric( as.character( mpmath$hyp1f1(1+1, 1+2, alpha_mme) ))
f2<- as.numeric( as.character( mpmath$hyp2f2(1, 1, 1+1, 1+1, alpha_mme) ) )

num =  sqrt(alpha_mme)*beta_mle/mu_mme*(f1 - (1/(1+1)*f11))
den = sqrt(beta_mle^2*(1/mu_mme*f2 - 1/mu_mme^2*f1^2) + beta_mle*(1 - 1/mu_mme*f1) + 0)
rho_mle = num/den
#AIC
k<- 2
l= -n*(alpha_mme+beta_mle+0) + log(alpha_mme)*sum(x1) + beta_mle*sum(1/(x1+1)) + sum(x2*log(0 + beta_mle*(1 - (1/(x1+1)))), na.rm = TRUE) +h
aic <- 2*k - 2*l
bic <- k*log(n) - 2*l  
neg2loglam <- -2*(l - l_full)
critical <- qchisq(p =0.95, df = 3)







