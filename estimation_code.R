rm(list=ls())
set.seed(2)
setwd("~/Desktop/paper/Lomax")
source(file = "./functions.R")
library(maxLik)
library(stats4)
library(numDeriv)

### Monte

LEORE <- function(p, alpha, beta, gamma, k){
  x_p <- log(1+(((1-p)^(-1/k)-1)/alpha)^(1/beta)) / gamma
  return(x_p)
}

### true parameters
alpha1 =1.5; beta1 =2.9; gamma1 = 1.3; k1=0.8
## fix alpha
init_cond = c(1.5, 3.0, 1.2, 0.7)


ns <- c(50, 100, 250, 500, 1000)
alphas <- c(); betas  <- c(); gammas <- c(); ks <- c()
error1 <- c(); error2 <- c(); error3 <- c(); error4 <- c()
parms1 <- c(); parms2 <- c(); parms3 <- c(); parms4 <- c()
mse1 <- c(); mse2 <- c(); mse3 <- c(); mse4 <- c()

NN = 500

for (j in 1:length(ns)){
  n <- ns[j]
  columns= c('param1','param2','param3')
  error = data.frame(matrix(nrow = 0, ncol = length(columns)))
  colnames(error) = columns
  
  for (k in 1:NN){
    Fx <- runif(n)
    x <- c(); ls_ins <- c(); wls_ins <- c(); ins <-c(); i21 <- c(); i2n1 <- c()
    for (i in 1:n) {
      ### simulation
      r <- LEORE(Fx[i], alpha1, beta1, gamma1, k1)
      x <- c(x, r)
      ins <- c(ins, (2*i-1)/(2*n))
      ls_ins <- c(ls_ins, i/(n+1))
      wls_ins <- c(wls_ins, (n+1)**2*(n+1)/(i*(n-i+1)))
      i21 <- c(i21, 2*i-1)
      i2n1 <- c(i2n1, 2*n+1-2*i)
    }
    x = sort(x)
    x_rev = sort(x, decreasing = TRUE)
    est_nlminb <- nlminb(init_cond, MPS_fun_param2)
    error[k,] = est_nlminb$par
  }
  
  error1 = c(error1, mean(error$param1)-beta1)
  error2 = c(error2, mean(error$param2)-gamma1)
  error3 = c(error3, mean(error$param3)-k1)

  parms1 = c(parms1, mean(error$param1))
  parms2 = c(parms2, mean(error$param2))
  parms3 = c(parms3, mean(error$param3))

  mse1 = c(mse1, sum((error$param1-beta1)**2)/NN)
  mse2 = c(mse2, sum((error$param2-gamma1)**2)/NN)
  mse3 = c(mse3, sum((error$param3-k1)**2)/NN)
}


error1
error2
error3

parms1
parms2
parms3

mse1
mse2
mse3



